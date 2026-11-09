//! Shared machinery for the cache-hosted k-mer panels ([crate::endpanel],
//! [crate::nucpanel]): parse a #k header + tab-separated kmer \t gene \t category
//! rows into a 2-bit lookup holding both strand orientations, and load/merge every
//! *.<ext> from the panel cache dir. Each panel module supplies only its category
//! enum, a row-category parser, and its own scoring.

use rustc_hash::FxHashMap;

use crate::sexpanel::encode_kmer_both_strands;

/// A panel lookup built on first use. On wasm it rebuilds when the front-end injects more
/// panel text, so a species panel loaded after the first check still takes effect.
pub struct PanelCache<V: 'static>(std::sync::Mutex<Option<(u64, &'static V)>>);

#[allow(clippy::new_without_default)] // const new backs the statics
impl<V> PanelCache<V> {
    pub const fn new() -> Self {
        Self(std::sync::Mutex::new(None))
    }

    // ponytail: each wasm rebuild leaks the old lookup; injections are one per species file.
    pub fn get(&self, build: impl FnOnce() -> V) -> &'static V {
        let gen = generation();
        let mut slot = self.0.lock().unwrap_or_else(|e| e.into_inner());
        match *slot {
            Some((g, v)) if g == gen => v,
            _ => {
                let v: &'static V = Box::leak(Box::new(build()));
                *slot = Some((gen, v));
                v
            }
        }
    }
}

#[cfg(not(target_arch = "wasm32"))]
fn generation() -> u64 {
    0
}

#[cfg(target_arch = "wasm32")]
fn generation() -> u64 {
    wasm_store::generation()
}

/// Fallback k when a panel omits its #k header (all shipped panels are k=50).
pub const DEFAULT_K: usize = 50;

/// A merged k-mer lookup (both orientations) keyed to a per-panel category T, plus the
/// panel's k.
pub struct KmerPanel<T> {
    pub map: FxHashMap<u128, T>,
    pub k: usize,
}

impl<T: Copy + PartialEq> KmerPanel<T> {
    /// Load and merge every *.{ext} panel in the cache dir. parse_cat maps the row's
    /// category column to T (returning None skips the row). Missing dir / no files /
    /// wasm → an empty panel (the caller's verdict degrades to "unknown").
    pub fn load(ext: &str, parse_cat: impl Fn(&str) -> Option<T>) -> KmerPanel<T> {
        let mut map = FxHashMap::default();
        let mut k = DEFAULT_K;
        for text in panel_texts(ext) {
            if let Some(hk) = header_k(&text) {
                k = hk;
            }
            merge(&mut map, &text, &parse_cat);
        }
        KmerPanel { map, k }
    }

    /// Parse a single panel's text (the non-cache seam the panel tests drive).
    #[cfg(test)]
    pub fn parse(text: &str, parse_cat: impl Fn(&str) -> Option<T>) -> KmerPanel<T> {
        let mut map = FxHashMap::default();
        merge(&mut map, text, &parse_cat);
        KmerPanel {
            map,
            k: header_k(text).unwrap_or(DEFAULT_K),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    /// Load panels grouped by their #species header (lowercased), merging same-species
    /// files. Panels lacking the header group under the empty key. Used where the caller
    /// must pick the panel matching a detected organism instead of merging all species.
    pub fn load_by_species(
        ext: &str,
        parse_cat: impl Fn(&str) -> Option<T>,
    ) -> std::collections::HashMap<String, KmerPanel<T>> {
        let mut out: std::collections::HashMap<String, KmerPanel<T>> =
            std::collections::HashMap::new();
        for text in panel_texts(ext) {
            let entry = out
                .entry(header_species(&text).unwrap_or_default())
                .or_insert_with(|| KmerPanel {
                    map: FxHashMap::default(),
                    k: DEFAULT_K,
                });
            if let Some(hk) = header_k(&text) {
                entry.k = hk;
            }
            merge(&mut entry.map, &text, &parse_cat);
        }
        out
    }
}

/// Text of every panel of this ext. Native: read every *.{ext} file in the panel
/// cache dir (sorted). Wasm (no filesystem): whatever the front-end injected via
/// [store_panel]; so the browser build must call load_panel(ext, bytes) before any
/// classify. Missing source yields an empty list (verdict degrades to "unknown").
#[cfg(not(target_arch = "wasm32"))]
fn panel_texts(ext: &str) -> Vec<String> {
    let Ok(dir) = crate::whitelists::panel_cache_dir() else {
        return Vec::new();
    };
    let Ok(entries) = std::fs::read_dir(&dir) else {
        return Vec::new();
    };
    let mut paths: Vec<_> = entries
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().and_then(|x| x.to_str()) == Some(ext))
        .collect();
    paths.sort();
    paths
        .iter()
        .filter_map(|p| std::fs::read_to_string(p).ok())
        .collect()
}

#[cfg(target_arch = "wasm32")]
fn panel_texts(ext: &str) -> Vec<String> {
    wasm_store::texts(ext)
}

/// Wasm-only in-memory panel store, populated by the browser front-end before it runs a
/// classify (there is no filesystem cache). Keyed by extension ("endpanel"/"nucpanel").
#[cfg(target_arch = "wasm32")]
pub mod wasm_store {
    use std::cell::RefCell;
    use std::collections::HashMap;

    thread_local! {
        static PANELS: RefCell<HashMap<String, Vec<String>>> = RefCell::new(HashMap::new());
        static GENERATION: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
    }

    /// Bumped on every injection so [super::PanelCache] knows to rebuild.
    pub fn generation() -> u64 {
        GENERATION.with(|g| g.get())
    }

    /// Append one panel's text under ext (multiple files of the same ext merge, as on
    /// native). Idempotent-ish: re-injecting the same ext accumulates, so the front-end
    /// should inject each panel once per page load.
    pub fn add(ext: &str, text: String) {
        PANELS.with(|p| {
            p.borrow_mut()
                .entry(ext.to_string())
                .or_default()
                .push(text)
        });
        GENERATION.with(|g| g.set(g.get() + 1));
    }

    pub fn texts(ext: &str) -> Vec<String> {
        PANELS.with(|p| p.borrow().get(ext).cloned().unwrap_or_default())
    }
}

/// Inject a panel's raw text for the wasm build (no-op signature kept simple). The
/// front-end fetches each *.endpanel/*.nucpanel and hands the bytes here.
#[cfg(target_arch = "wasm32")]
pub fn store_panel(ext: &str, text: String) {
    wasm_store::add(ext, text);
}

fn header_k(text: &str) -> Option<usize> {
    text.lines()
        .find_map(|l| l.strip_prefix("#k\t").and_then(|r| r.trim().parse().ok()))
}

fn header_species(text: &str) -> Option<String> {
    text.lines().find_map(|l| {
        l.strip_prefix("#species\t")
            .map(|r| r.trim().to_ascii_lowercase())
    })
}

/// Fold one panel's rows into map, storing both strand orientations so a raw read
/// window matches regardless of cDNA strand. A k-mer present with a *different*
/// category is a cross-panel conflict and is dropped (trustworthy as neither).
fn merge<T: Copy + PartialEq>(
    map: &mut FxHashMap<u128, T>,
    text: &str,
    parse_cat: &impl Fn(&str) -> Option<T>,
) {
    for line in text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut f = line.split('\t');
        let (kmer, cat) = match (f.next(), f.next(), f.next()) {
            (Some(km), Some(_gene), Some(c)) => (km, c),
            _ => continue,
        };
        let Some(t) = parse_cat(cat) else { continue };
        let Some((fwd, rc)) = encode_kmer_both_strands(kmer.as_bytes()) else {
            continue;
        };
        for enc in [fwd, rc] {
            match map.get(&enc) {
                None => {
                    map.insert(enc, t);
                }
                Some(existing) if *existing == t => {}
                Some(_) => {
                    map.remove(&enc);
                }
            }
        }
    }
}
