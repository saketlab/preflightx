//! seqspec ingestion: parse a seqspec YAML into a library structure and resolve
//! each sequencing read into an ordered list of segments (barcode / UMI / cDNA /
//! linker / ...). This is the generic substrate the detector scores against, so
//! any assay describable in seqspec can be reasoned about without a hand-coded
//! `AssaySpec`. See `docs/seqspec-support-plan.md`.

use anyhow::{Context, Result};
use serde::Deserialize;

/// A barcode/onlist reference inside a region. Only the identity fields matter to
/// preflightx; the actual list is fetched on demand, never bundled.
#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct Onlist {
    #[serde(default)]
    pub filename: String,
    #[serde(default)]
    pub md5: String,
    #[serde(default)]
    pub url: String,
}

/// One node in the seqspec `library_spec` region tree.
#[derive(Debug, Clone, Deserialize)]
pub struct Region {
    pub region_id: String,
    pub region_type: String,
    #[serde(default)]
    pub sequence_type: String,
    #[serde(default)]
    pub sequence: String,
    #[serde(default)]
    pub min_len: usize,
    #[serde(default)]
    pub max_len: usize,
    #[serde(default)]
    pub onlist: Option<Onlist>,
    #[serde(default)]
    pub regions: Vec<Region>,
}

impl Region {
    fn is_leaf(&self) -> bool {
        self.regions.is_empty()
    }
}

/// seqspec's `region_type` is free text, so the same concept appears under many names.
/// These normalize the in-read **cell/well-identifying barcode** synonyms and the UMI
/// synonyms seen across the catalogue, so geometry extraction isn't limited to specs
/// that happen to use the literal `barcode`/`umi` strings. Sample/lane sequencing
/// indices (`index7`, `i7`, `sample_index`, …) are deliberately **excluded** — they
/// identify the library, not the cell, and live on separate index reads.
fn is_barcode_type(region_type: &str) -> bool {
    matches!(
        region_type.to_ascii_lowercase().as_str(),
        "barcode"
            | "cell_bc"
            | "cell_barcode"
            | "cell_index"
            | "bead_bc"
            | "well_bc"
            | "rt_bc"
            | "plate_bc"
            | "sample_bc"
            | "subarray_bc"
            | "hairpin_barcode"
            | "bc1"
            | "bc2"
            | "bc3"
    )
}

fn is_umi_type(region_type: &str) -> bool {
    region_type.eq_ignore_ascii_case("umi")
}

/// Fixed inter-barcode sequences (and short variable spacers) that are *sequenced*
/// as part of the barcode read — they don't terminate it the way the cDNA insert does.
fn is_linker_type(region_type: &str) -> bool {
    let rt = region_type.to_ascii_lowercase();
    rt.starts_with("link") || rt.starts_with("spacer")
}

/// The biological insert (cDNA/gDNA/transcript/poly tail) that the barcode read runs
/// into. Some specs write these as fixed-length regions, so type — not just
/// variable-length — must be recognized.
fn is_insert_type(region_type: &str) -> bool {
    let rt = region_type.to_ascii_lowercase();
    matches!(
        rt.as_str(),
        "cdna" | "gdna" | "rna" | "dna" | "insert" | "hic" | "methyl" | "atac" | "protein"
    ) || rt.starts_with("poly") // poly_T / polyA / ploy_t (typo'd in some specs)
        || rt.starts_with("ploy")
}

/// A region that absorbs the rest of the read (the biological insert): an insert-type
/// region, or any *variable* region that is not a barcode, UMI, or linker. This — not
/// "any variable region" — is what terminates a read's structural layout.
fn is_read_terminator(s: &Segment) -> bool {
    is_insert_type(&s.region_type)
        || (s.max_len > s.len
            && !is_barcode_type(&s.region_type)
            && !is_umi_type(&s.region_type)
            && !is_linker_type(&s.region_type))
}

/// One sequencing read in the seqspec `sequence_spec`.
#[derive(Debug, Clone, Deserialize)]
pub struct Read {
    pub read_id: String,
    #[serde(default)]
    pub primer_id: String,
    #[serde(default)]
    pub modality: String,
    #[serde(default)]
    pub max_len: usize,
    #[serde(default)]
    pub strand: String,
}

/// A parsed seqspec assay.
#[derive(Debug, Clone, Deserialize)]
pub struct Assay {
    #[serde(default)]
    pub assay_id: String,
    /// Source spec filename (e.g. `split_seq.spec.yaml`); the stable catalogue key,
    /// since `assay_id` is frequently empty in upstream specs. Set by `embedded_assays`.
    #[serde(default)]
    pub spec_file: String,
    #[serde(default)]
    pub name: String,
    /// Publication DOI for the assay/protocol (often a `https://doi.org/...` URL;
    /// frequently empty upstream). Consumed by the doc generator to hyperlink each assay.
    #[serde(default)]
    #[allow(dead_code)]
    pub doi: String,
    #[serde(default)]
    pub modalities: Vec<String>,
    #[serde(default)]
    pub lib_struct: String,
    #[serde(default)]
    pub sequence_spec: Vec<Read>,
    #[serde(default)]
    pub library_spec: Vec<Region>,
}

/// One resolved segment a read covers, in read order.
#[derive(Debug, Clone, PartialEq)]
pub struct Segment {
    pub region_id: String,
    pub region_type: String,
    pub sequence_type: String,
    /// Expected length in this read (min_len of the region).
    pub len: usize,
    pub max_len: usize,
    /// Fixed sequence (oriented for this read's strand) for linker/adapter anchors.
    pub sequence: Option<String>,
    pub onlist: Option<Onlist>,
}

impl Assay {
    pub fn from_yaml(s: &str) -> Result<Assay> {
        serde_norway::from_str(s).context("failed to parse seqspec YAML")
    }

    /// Depth-first, in-order list of leaf regions (the physical layout of the molecule).
    pub fn leaf_regions(&self) -> Vec<&Region> {
        fn walk<'a>(regions: &'a [Region], out: &mut Vec<&'a Region>) {
            for r in regions {
                if r.is_leaf() {
                    out.push(r);
                } else {
                    walk(&r.regions, out);
                }
            }
        }
        let mut out = Vec::new();
        walk(&self.library_spec, &mut out);
        out
    }

    /// Resolve a read (by `read_id`) into the ordered segments it covers.
    ///
    /// A read primes from its `primer_id` region. A `pos`-strand read then covers
    /// the regions that follow the primer, in order; a `neg`-strand read covers the
    /// regions that precede the primer, in reverse (it sequences the other strand).
    /// We walk outward accumulating region `min_len` until the read's `max_len` is
    /// consumed, yielding the ordered barcode/UMI/cDNA/linker segments.
    pub fn read_layout(&self, read_id: &str) -> Option<Vec<Segment>> {
        let read = self.sequence_spec.iter().find(|r| r.read_id == read_id)?;
        let leaves = self.leaf_regions();
        let neg = read.strand == "neg";

        // The read covers the regions downstream of its primer in the strand direction.
        // A `generic_primer1`/`generic_primer2` (or any primer_id not present as a
        // region) is a placeholder for the molecule boundary: the read starts at the
        // 5' end (pos) or the 3' end (neg) and covers the whole molecule from there.
        let ordered: Vec<&Region> = match leaves.iter().position(|r| r.region_id == read.primer_id)
        {
            Some(primer_idx) if neg => leaves[..primer_idx].iter().rev().copied().collect(),
            Some(primer_idx) => leaves[primer_idx + 1..].to_vec(),
            None if neg => leaves.iter().rev().copied().collect(),
            None => leaves.clone(),
        };

        let cap = read.max_len.max(1);
        let mut segs = Vec::new();
        let mut acc = 0usize;
        for r in ordered {
            if acc >= cap {
                break;
            }
            let sequence = if matches!(r.sequence_type.as_str(), "fixed") && !r.sequence.is_empty()
            {
                Some(if neg {
                    revcomp(&r.sequence)
                } else {
                    r.sequence.clone()
                })
            } else {
                None
            };
            let seg = Segment {
                region_id: r.region_id.clone(),
                region_type: r.region_type.clone(),
                sequence_type: r.sequence_type.clone(),
                len: r.min_len,
                max_len: r.max_len,
                sequence,
                onlist: r.onlist.clone(),
            };
            // The cDNA/gDNA insert (a variable, non-barcode/UMI/linker region) absorbs
            // the rest of the read, so nothing structural can follow it. Short variable
            // spacers and variable-length barcodes do NOT terminate the layout.
            let terminator = is_read_terminator(&seg);
            segs.push(seg);
            if terminator {
                break;
            }
            acc += r.min_len.max(1);
        }
        Some(segs)
    }

    /// The read that carries the cell barcode(s) — the structural read the detector
    /// scores against. When an assay has several barcode-bearing reads (plate vs cell
    /// barcode, or RNA vs ATAC/protein modalities), prefer the one that best looks like
    /// the primary GEX cell-barcode read: most barcodes, barcodes backed by an onlist,
    /// the RNA modality, a UMI present, more fixed anchors, and no cDNA in the read.
    /// Returns `(read_id, layout)`.
    pub fn structural_read(&self) -> Option<(String, Vec<Segment>)> {
        let mut best: Option<(String, Vec<Segment>, i64)> = None;
        for read in &self.sequence_spec {
            let Some(layout) = self.read_layout(&read.read_id) else {
                continue;
            };
            let n_bc = layout
                .iter()
                .filter(|s| is_barcode_type(&s.region_type))
                .count();
            if n_bc == 0 {
                continue;
            }
            let onlist_bc = layout
                .iter()
                .filter(|s| is_barcode_type(&s.region_type) && s.onlist.is_some())
                .count();
            let has_umi = layout.iter().any(|s| is_umi_type(&s.region_type));
            let has_cdna = layout.iter().any(|s| s.max_len > s.len);
            let n_anchor = layout.iter().filter(|s| s.sequence.is_some()).count();
            let is_rna = read.modality.eq_ignore_ascii_case("rna");
            let score = n_bc as i64 * 100
                + onlist_bc as i64 * 30
                + if is_rna { 40 } else { 0 }
                + if has_umi { 20 } else { 0 }
                + n_anchor as i64 * 5
                + if has_cdna { 0 } else { 25 };
            if best.as_ref().map(|b| score > b.2).unwrap_or(true) {
                best = Some((read.read_id.clone(), layout, score));
            }
        }
        best.map(|(id, layout, _)| (id, layout))
    }
}

// Embedded seqspec catalogue, generated by build.rs from vendor/seqspec (submodule).
include!(concat!(env!("OUT_DIR"), "/specs.rs"));

/// Parse every embedded seqspec into an `Assay`. Specs that fail to parse are
/// skipped (returned in the error list) rather than aborting the catalogue.
pub fn embedded_assays() -> Vec<Assay> {
    EMBEDDED_SPECS
        .iter()
        .filter_map(|(name, yaml)| {
            let mut a = Assay::from_yaml(yaml).ok()?;
            a.spec_file = name.to_string();
            if a.assay_id.trim().is_empty() {
                a.assay_id = name.trim_end_matches(".spec.yaml").to_string();
            }
            Some(a)
        })
        .collect()
}

/// A single-cell barcode-read geometry derived from one catalogued assay's
/// structural (barcode-bearing) read — the seqspec-driven source of truth for "what
/// does a single-cell barcode read look like", used by `size`/`composition`/orientation.
#[derive(Debug, Clone)]
pub struct ScGeometry {
    pub spec_file: String,
    pub name: String,
    /// Expected length of the *observed* barcode read: the barcode/UMI/linker span
    /// after any leading primer/adapter and before the first variable (cDNA) region.
    pub read_len: usize,
}

impl ScGeometry {
    /// A human label: the assay `name`, falling back to the spec filename stem.
    pub fn label(&self) -> &str {
        if self.name.is_empty() {
            self.spec_file.trim_end_matches(".spec.yaml")
        } else {
            &self.name
        }
    }
}

/// Reduce a faithful structural-read layout to the *observed* barcode-read geometry.
///
/// seqspec's `map_read_id_to_regions` (which `read_layout` mirrors) yields every region
/// downstream of the read's primer — including any library handle/primer between the
/// primer and the barcode (e.g. the 39 bp Drop-seq Read1 custom primer), which is a
/// binding site, not sequenced. The sequenced barcode read therefore runs from the
/// first barcode/UMI region through the **last** barcode/UMI before the cDNA/gDNA
/// insert, *including* any adapters/linkers physically between them (they are
/// sequenced — e.g. scifi-RNA-seq's inter-barcode adapters). A trailing adapter after
/// the last barcode (a read-primer region a loose `primer_id` left in the read, e.g.
/// quartz-seq's read-2 adapter) is dropped. Variable-length barcodes (inDrop) are
/// kept. Returns the observed barcode-read length, or `None` if there is no barcode.
fn observed_read_len(layout: &[Segment]) -> Option<usize> {
    let is_bc_or_umi = |s: &Segment| is_barcode_type(&s.region_type) || is_umi_type(&s.region_type);
    let start = layout.iter().position(&is_bc_or_umi)?;
    // The window runs from the first barcode/UMI to the cDNA/gDNA insert (or the end).
    let end = layout[start..]
        .iter()
        .position(is_read_terminator)
        .map_or(layout.len(), |i| start + i);
    let window = &layout[start..end];
    // End at the last barcode/UMI; drop any trailing (post-barcode) adapter run.
    let last = window.iter().rposition(&is_bc_or_umi)?;
    let span = &window[..=last];
    span.iter()
        .any(|s| is_barcode_type(&s.region_type))
        .then(|| span.iter().map(|s| s.len).sum())
}

/// Every single-cell barcode-read geometry in the embedded catalogue (one per assay
/// with a barcode-bearing structural read). Parsed once and cached.
pub fn sc_geometries() -> &'static [ScGeometry] {
    static CACHE: std::sync::OnceLock<Vec<ScGeometry>> = std::sync::OnceLock::new();
    CACHE.get_or_init(|| {
        embedded_assays()
            .iter()
            .filter_map(|a| {
                let (_rid, layout) = a.structural_read()?;
                let read_len = observed_read_len(&layout)?;
                Some(ScGeometry {
                    spec_file: a.spec_file.clone(),
                    name: a.name.clone(),
                    read_len,
                })
            })
            .collect()
    })
}

/// Catalogued single-cell geometries whose observed barcode-read length is within
/// `tol` bp of `len` (empty for `len == 0`). The seqspec-driven source for "is this a
/// known single-cell barcode-read length", used by `size`/`composition`/orientation.
pub fn sc_geometries_matching(len: usize, tol: usize) -> Vec<&'static ScGeometry> {
    if len == 0 {
        return Vec::new();
    }
    sc_geometries()
        .iter()
        .filter(|g| g.read_len.abs_diff(len) <= tol)
        .collect()
}

/// Reverse-complement a fixed sequence (IUPAC-aware for the common bases).
fn revcomp(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            b'N' => 'N',
            other => other as char,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn geom_of(spec_file: &str) -> Option<usize> {
        sc_geometries()
            .iter()
            .find(|g| g.spec_file == spec_file)
            .map(|g| g.read_len)
    }

    /// The observed barcode-read span (first barcode/UMI through the last, before the
    /// insert) — the segments `observed_read_len` sums. Used by the doc generator.
    fn observed_span(layout: &[Segment]) -> Vec<&Segment> {
        let is_bcu = |s: &Segment| is_barcode_type(&s.region_type) || is_umi_type(&s.region_type);
        let Some(start) = layout.iter().position(&is_bcu) else {
            return vec![];
        };
        let end = layout[start..]
            .iter()
            .position(is_read_terminator)
            .map_or(layout.len(), |i| start + i);
        let window = &layout[start..end];
        let Some(last) = window.iter().rposition(&is_bcu) else {
            return vec![];
        };
        window[..=last].iter().collect()
    }

    /// Regenerate `docs/assay-geometries.md` from the catalogue. Run on demand:
    ///   cargo test gen_assay_docs -- --ignored
    #[test]
    #[ignore = "doc generator; run explicitly to refresh docs/assay-geometries.md"]
    fn gen_assay_docs() {
        fn short(region_type: &str) -> String {
            if is_barcode_type(region_type) {
                "CB".into()
            } else if is_umi_type(region_type) {
                "UMI".into()
            } else if is_linker_type(region_type) {
                "link".into()
            } else {
                region_type.to_ascii_lowercase()
            }
        }
        fn ascii(a: &Assay) -> Option<(String, usize, String)> {
            let (rid, layout) = a.structural_read()?;
            let span = observed_span(&layout);
            if span.is_empty() {
                return None;
            }
            let len: usize = span.iter().map(|s| s.len).sum();
            let boxes: Vec<String> = span
                .iter()
                .map(|s| format!("[ {} {} ]", short(&s.region_type), s.len))
                .collect();
            let end_idx = layout
                .iter()
                .position(|s| s.region_id == span.last().unwrap().region_id)
                .map(|i| i + 1)
                .unwrap_or(layout.len());
            let insert = layout
                .get(end_idx)
                .map(|s| format!("  ─▶ {}", s.region_type.to_ascii_lowercase()))
                .unwrap_or_default();
            Some((rid, len, format!("{}{}", boxes.join("─"), insert)))
        }

        let assays = embedded_assays();
        let mut out = String::new();
        out.push_str("# preflightx — resolved assay read geometries\n\n");
        out.push_str(
            "_Generated from the vendored seqspec catalogue (`vendor/seqspec`) by\n\
             `cargo test gen_assay_docs -- --ignored`. Do not edit by hand._\n\n\
             preflightx derives each assay's **observed barcode-read geometry** straight from\n\
             its seqspec description: it picks the barcode-bearing read (preferring the RNA\n\
             cell-barcode read for multimodal assays), strips leading library primers/adapters,\n\
             and spans the barcode / UMI / linker structure up to the cDNA insert. `[ CB n ]` is\n\
             a cell barcode, `[ UMI n ]` a UMI, `[ link n ]` a fixed inter-barcode linker, and\n\
             `─▶ cdna` marks where the biological insert begins.\n\n",
        );

        let featured = [
            "10x_rna_v3.spec.yaml",
            "10x_rna_5prime.spec.yaml",
            "drop_seq.spec.yaml",
            "split_seq.spec.yaml",
            "indropv2.spec.yaml",
            "bd_rhapsody_v1.spec.yaml",
            "scifi_rna_seq.spec.yaml",
            "10x_rna_atac_multi.spec.yaml",
        ];
        out.push_str("## Read structure (representative assays)\n\n```\n");
        for sf in featured {
            if let Some(a) = assays.iter().find(|a| a.spec_file == sf) {
                if let Some((rid, len, art)) = ascii(a) {
                    let name = if a.name.is_empty() {
                        a.spec_file.trim_end_matches(".spec.yaml")
                    } else {
                        &a.name
                    };
                    out.push_str(&format!("{name}   ({rid}, {len} bp)\n  {art}\n\n"));
                }
            }
        }
        out.push_str("```\n\n");

        let mut rows: Vec<(String, String)> = Vec::new();
        let (mut n_cov, mut n_plate) = (0, 0);
        for a in &assays {
            let name = if a.name.is_empty() {
                a.spec_file.trim_end_matches(".spec.yaml").to_string()
            } else {
                a.name.clone()
            };
            let geom = a.structural_read().and_then(|(rid, l)| {
                let span = observed_span(&l);
                if span.is_empty() {
                    return None;
                }
                let len: usize = span.iter().map(|s| s.len).sum();
                let cb: usize = span
                    .iter()
                    .filter(|s| is_barcode_type(&s.region_type))
                    .map(|s| s.len)
                    .sum();
                let umi: usize = span
                    .iter()
                    .filter(|s| is_umi_type(&s.region_type))
                    .map(|s| s.len)
                    .sum();
                Some(format!(
                    "{rid}, {len} bp (CB {cb}{})",
                    if umi > 0 {
                        format!(" + UMI {umi}")
                    } else {
                        String::new()
                    }
                ))
            });
            let cell = match geom {
                Some(g) => {
                    n_cov += 1;
                    g
                }
                None => {
                    n_plate += 1;
                    "_plate-based — no in-read cell barcode_".into()
                }
            };
            rows.push((
                a.spec_file.clone(),
                format!(
                    "| {} | `{}` | {} | {} |",
                    name,
                    a.spec_file.trim_end_matches(".spec.yaml"),
                    a.modalities.join(", "),
                    cell
                ),
            ));
        }
        rows.sort();
        out.push_str(&format!(
            "## All {} catalogued assays ({} with barcode geometry, {} plate-based)\n\n",
            assays.len(),
            n_cov,
            n_plate
        ));
        out.push_str("| Assay | spec | modalities | barcode read |\n");
        out.push_str("|-------|------|------------|--------------|\n");
        for (_, r) in &rows {
            out.push_str(r);
            out.push('\n');
        }
        out.push('\n');

        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/docs/assay-geometries.md");
        std::fs::create_dir_all(concat!(env!("CARGO_MANIFEST_DIR"), "/docs")).unwrap();
        std::fs::write(path, out).unwrap();
        eprintln!("wrote {path}");
    }

    /// Regenerate `docs/assay-structures.html` — a colour-coded visual of every
    /// assay's barcode-read structure. Run on demand:
    ///   cargo test gen_assay_html -- --ignored
    #[test]
    #[ignore = "doc generator; run explicitly to refresh docs/assay-structures.html"]
    fn gen_assay_html() {
        // Map a free-text region_type to one of the legend categories.
        fn category(rt: &str) -> &'static str {
            let l = rt.to_ascii_lowercase();
            if is_barcode_type(rt) {
                "cb"
            } else if is_umi_type(rt) {
                "umi"
            } else if is_linker_type(rt) {
                "linker"
            } else if l.starts_with("poly") || l.starts_with("ploy") {
                "poly"
            } else if matches!(
                l.as_str(),
                "cdna" | "gdna" | "rna" | "dna" | "insert" | "hic" | "methyl" | "atac" | "protein"
            ) {
                "insert"
            } else if l.starts_with("index")
                || l == "i5"
                || l == "i7"
                || l.ends_with("_index")
                || l == "sample_index"
            {
                "index"
            } else {
                "adapter"
            }
        }
        fn label(rt: &str) -> String {
            match category(rt) {
                "cb" => "BC".into(),
                "umi" => "UMI".into(),
                "linker" => "link".into(),
                "poly" => rt.to_ascii_lowercase().replace('_', ""),
                _ => rt.to_ascii_lowercase(),
            }
        }
        fn esc(s: &str) -> String {
            s.replace('&', "&amp;")
                .replace('<', "&lt;")
                .replace('>', "&gt;")
        }

        // Render one read's resolved layout into a track of segment blocks. Strips the
        // leading primer-binding handle (everything before the first sequenced
        // barcode/UMI/insert region) so each track shows only what is actually
        // sequenced. The variable cDNA/gDNA insert is drawn as a hatched terminal block
        // and absorbs the rest of the read. Returns (track_html, sequenced_bp), or None
        // if the read covers no sequenced content (a pure primer read).
        let read_track = |layout: &[Segment]| -> Option<(String, usize)> {
            let is_content = |s: &Segment| {
                is_barcode_type(&s.region_type)
                    || is_umi_type(&s.region_type)
                    || is_insert_type(&s.region_type)
            };
            let start = layout.iter().position(is_content)?;
            let mut track = String::new();
            let mut bp = 0usize;
            for s in &layout[start..] {
                let variable_insert =
                    is_insert_type(&s.region_type) && (s.max_len > s.len || s.len == 0);
                if variable_insert {
                    track.push_str(&format!(
                        "<div class=\"seg {cat} tail\" title=\"{rt} (insert)\">{lab} ▸</div>",
                        cat = category(&s.region_type),
                        rt = esc(&s.region_type),
                        lab = esc(&s.region_type.to_ascii_lowercase()),
                    ));
                    break;
                }
                bp += s.len;
                let px = (s.len.max(1) * 4).clamp(34, 200);
                track.push_str(&format!(
                    "<div class=\"seg {cat}\" style=\"flex:0 0 {px}px\" title=\"{rid_t} · {rt} · {n} bp\"><b>{lab}</b><i>{n}</i></div>",
                    cat = category(&s.region_type),
                    rid_t = esc(&s.region_id),
                    rt = esc(&s.region_type),
                    n = s.len,
                    lab = esc(&label(&s.region_type)),
                ));
            }
            Some((track, bp))
        };

        // Placeholder / template / mislabelled example specs upstream that aren't real
        // assays — excluded from the visual so the catalogue reads clean.
        const PLACEHOLDER_SPECS: &[&str] = &[
            "template.spec.yaml",
            "tang2009.spec.yaml", // upstream `name:` is the bogus "10x Genomics"
            "10xv3_scrnaseq_element_adept_truseq_dual.spec.yaml", // name: "Example assay"
            "10xv3_scrnaseq_illumina_truseq_dual.spec.yaml",      // name: "Example assay"
        ];
        let assays: Vec<Assay> = embedded_assays()
            .into_iter()
            .filter(|a| !PLACEHOLDER_SPECS.contains(&a.spec_file.as_str()))
            .collect();

        let card = |a: &Assay| -> Option<(usize, bool, String)> {
            let bc = a.structural_read();
            let bc_read_id = bc.as_ref().map(|(id, _)| id.as_str());
            let bc_len = bc.as_ref().and_then(|(_, l)| observed_read_len(l));

            let mut rows = String::new();
            let mut n_rows = 0;
            for read in &a.sequence_spec {
                let Some(layout) = a.read_layout(&read.read_id) else {
                    continue;
                };
                let Some((track, bp)) = read_track(&layout) else {
                    continue;
                };
                let rid = read
                    .read_id
                    .trim_end_matches(".gz")
                    .trim_end_matches(".fastq")
                    .trim_end_matches(".fq");
                let is_bc = bc_read_id == Some(read.read_id.as_str());
                let badge = if is_bc {
                    "<span class=\"bcflag\">barcode</span>".to_string()
                } else {
                    String::new()
                };
                let modtag = if read.modality.is_empty() {
                    String::new()
                } else {
                    format!("<span class=\"rmod\">{}</span>", esc(&read.modality))
                };
                let bptag = if bp > 0 {
                    format!("<span class=\"rbp\">{bp} bp</span>")
                } else {
                    String::new()
                };
                rows.push_str(&format!(
                    "<div class=\"readrow{bccls}\">\
                       <div class=\"rlabel\">{rid}{badge}{modtag}{bptag}</div>\
                       <div class=\"readline\"><span class=\"end\">5′</span><div class=\"track\">{track}</div><span class=\"end\">3′</span></div>\
                     </div>",
                    bccls = if is_bc { " is-bc" } else { "" },
                    rid = esc(rid),
                ));
                n_rows += 1;
            }
            if n_rows == 0 {
                return None;
            }

            let name = if a.name.is_empty() {
                a.spec_file.trim_end_matches(".spec.yaml")
            } else {
                &a.name
            };
            let mods = a
                .modalities
                .iter()
                .map(|m| format!("<span class=\"mod\">{}</span>", esc(m)))
                .collect::<String>();
            let doi = a.doi.trim();
            let doi_link = if doi.is_empty() {
                String::new()
            } else {
                let href = if doi.starts_with("http") {
                    doi.to_string()
                } else {
                    format!("https://doi.org/{}", doi.trim_start_matches("doi:").trim())
                };
                format!(
                    "<a class=\"doi\" href=\"{href}\" target=\"_blank\" rel=\"noopener\">DOI ↗</a>",
                    href = esc(&href)
                )
            };
            let is_plate = bc_read_id.is_none();
            let cap = match (bc_read_id, bc_len) {
                (Some(rid), Some(len)) => format!(
                    "{rid} · <strong>{len} bp</strong> barcode read",
                    rid = esc(rid)
                ),
                _ => "plate-based · cell identity is the well, no in-read barcode".to_string(),
            };
            let html = format!(
                "<article class=\"card{platecls}\">\
                   <header><h2>{name}</h2><div class=\"meta\">{mods}<code>{spec}</code>{doi_link}</div></header>\
                   <div class=\"reads\">{rows}</div>\
                   <div class=\"cap\">{cap}</div>\
                 </article>",
                platecls = if is_plate { " plate-card" } else { "" },
                name = esc(name),
                spec = esc(a.spec_file.trim_end_matches(".spec.yaml")),
            );
            Some((bc_len.unwrap_or(0), is_plate, html))
        };

        let all_cards: Vec<(usize, bool, String)> = assays.iter().filter_map(card).collect();
        let mut cards: Vec<(usize, String)> = all_cards
            .iter()
            .filter(|(_, is_plate, _)| !is_plate)
            .map(|(len, _, html)| (*len, html.clone()))
            .collect();
        cards.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
        let mut plate_cards: Vec<String> = all_cards
            .iter()
            .filter(|(_, is_plate, _)| *is_plate)
            .map(|(_, _, html)| html.clone())
            .collect();
        plate_cards.sort();
        // Any assay whose reads don't render at all (no barcode/UMI/insert content) —
        // listed by name so the catalogue stays complete.
        let plate: Vec<String> = assays
            .iter()
            .filter(|a| card(a).is_none())
            .map(|a| {
                esc(if a.name.is_empty() {
                    a.spec_file.trim_end_matches(".spec.yaml")
                } else {
                    &a.name
                })
            })
            .collect();

        // Clean, Apple-HIG-inspired styling: plain system sans (Arial), generous
        // whitespace, hairline separators, restrained colour used only to encode region
        // type. No webfonts, no texture — self-contained and offline.
        let css = r#"
*{box-sizing:border-box}
:root{
  --bg:#ffffff; --panel:#ffffff; --ink:#1d1d1f; --muted:#6e6e73; --faint:#86868b;
  --hair:#d2d2d7; --hair2:#e8e8ed; --tint:#f5f5f7;
  --cb:#0071e3; --umi:#34a853; --linker:#8e8e93; --insert:#ff6a00;
  --poly:#f5a623; --index:#a64ddb; --adapter:#c7c7cc;
}
html{scroll-behavior:smooth;-webkit-font-smoothing:antialiased;text-rendering:optimizeLegibility}
body{margin:0;color:var(--ink);background:var(--bg);
  font-family:Arial,-apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,sans-serif;
  font-size:15px;line-height:1.5}
.wrap{max-width:1000px;margin:0 auto;padding:72px 24px 110px}
.title{font-weight:700;font-size:clamp(32px,5.5vw,52px);line-height:1.05;
  letter-spacing:-.022em;margin:0 0 16px}
.title em{font-style:normal;color:var(--cb)}
.lede{max-width:68ch;color:var(--muted);font-size:16px;line-height:1.6;margin:0 0 34px}
.lede a{color:var(--cb);text-decoration:none}
.lede a:hover{text-decoration:underline}
.legend{display:flex;flex-wrap:wrap;align-items:center;gap:10px 18px;padding:14px 18px;margin:0 0 44px;
  background:rgba(255,255,255,.72);backdrop-filter:saturate(180%) blur(20px);
  border:1px solid var(--hair);border-radius:14px;position:sticky;top:12px;z-index:5}
.legend b{font-weight:600;margin-right:4px;font-size:13px;color:var(--faint);
  text-transform:uppercase;letter-spacing:.06em}
.key{display:inline-flex;align-items:center;gap:7px;font-size:13px;color:var(--ink)}
.sw{width:14px;height:14px;border-radius:4px}
.card{background:var(--panel);border:1px solid var(--hair2);border-radius:18px;
  padding:20px 24px;margin:0 0 16px;box-shadow:0 1px 2px rgba(0,0,0,.04);
  opacity:0;transform:translateY(8px);animation:rise .45s ease forwards}
@keyframes rise{to{opacity:1;transform:none}}
.card.plate-card{border-style:dashed;border-color:var(--hair);background:var(--tint)}
.card.plate-card .cap{color:var(--faint)}
.card header{display:flex;justify-content:space-between;align-items:baseline;
  gap:14px;flex-wrap:wrap;margin-bottom:16px}
.card h2{font-weight:600;font-size:19px;margin:0;letter-spacing:-.01em}
.meta{display:flex;align-items:center;gap:8px;flex-wrap:wrap}
.mod{font-size:11px;letter-spacing:.02em;color:var(--muted);
  background:var(--tint);border-radius:999px;padding:2px 10px}
.meta code{font-size:12px;color:var(--faint);font-family:ui-monospace,Menlo,Consolas,monospace}
.doi{font-size:11px;color:var(--cb);text-decoration:none;font-weight:600;
  border:1px solid var(--hair);border-radius:999px;padding:1px 9px}
.doi:hover{background:var(--cb);color:#fff;border-color:var(--cb)}
.reads{display:flex;flex-direction:column;gap:10px}
.readrow{display:grid;grid-template-columns:124px 1fr;align-items:center;gap:16px}
.rlabel{display:flex;align-items:center;flex-wrap:wrap;gap:6px;font-size:12px;color:var(--muted);
  font-weight:600;word-break:break-all}
.readrow.is-bc .rlabel{color:var(--ink)}
.bcflag{font-size:10px;letter-spacing:.01em;color:#fff;
  background:var(--cb);border-radius:999px;padding:1px 8px;font-weight:600}
.rmod{font-size:10px;letter-spacing:.02em;color:var(--muted);
  background:var(--tint);border-radius:999px;padding:1px 8px}
.rbp{font-size:11px;color:var(--faint);font-weight:500}
.readline{display:flex;align-items:center;gap:12px;min-width:0}
.end{color:var(--faint);font-size:12px;font-weight:600}
.track{display:flex;align-items:stretch;gap:3px;flex:1;overflow-x:auto;padding:1px 0}
.seg{flex:0 0 auto;min-height:30px;border-radius:7px;display:flex;flex-direction:row;
  align-items:center;justify-content:center;gap:5px;color:#fff;position:relative;padding:0 8px;
  transition:transform .15s ease, box-shadow .15s ease;cursor:default}
.seg:hover{transform:translateY(-1px);box-shadow:0 4px 12px rgba(0,0,0,.14);z-index:2}
.seg b{font-weight:600;font-size:11px;letter-spacing:.01em}
.seg i{font-style:normal;font-size:10px;opacity:.85}
.seg.cb{background:var(--cb)} .seg.umi{background:var(--umi)}
.seg.insert{background:var(--insert)} .seg.index{background:var(--index)}
.seg.poly{background:var(--poly);color:#3a2c00} .seg.poly i{opacity:.7}
.seg.linker{background:var(--linker)} .seg.linker i{opacity:.8}
.seg.adapter{background:var(--adapter);color:#3a3a3c} .seg.adapter i{opacity:.7}
.seg.tail{flex:0 0 auto;padding:0 14px;min-width:74px;color:var(--insert);font-weight:600;font-size:11px;
  background:transparent;border:1.5px dashed var(--insert)}
.cap{margin-top:14px;font-size:12px;color:var(--muted)}
.cap strong{color:var(--ink);font-weight:600}
.section{font-size:12px;letter-spacing:.06em;text-transform:uppercase;font-weight:600;
  color:var(--faint);margin:52px 0 18px;border-top:1px solid var(--hair2);padding-top:24px}
.plate{display:flex;flex-wrap:wrap;gap:8px}
.plate span{font-size:13px;color:var(--muted);background:var(--tint);border:1px solid var(--hair2);
  border-radius:8px;padding:6px 12px}
.foot{margin-top:60px;color:var(--faint);font-size:12px;line-height:1.7;
  border-top:1px solid var(--hair2);padding-top:20px}
"#;

        let legend = [
            ("cb", "cell barcode"),
            ("umi", "UMI"),
            ("linker", "linker / spacer"),
            ("poly", "polyA / polyT"),
            ("insert", "cDNA / gDNA insert"),
            ("index", "sample index"),
            ("adapter", "primer / adapter"),
        ]
        .iter()
        .map(|(c, name)| {
            format!(
                "<span class=\"key\"><span class=\"sw\" style=\"background:var(--{c})\"></span>{name}</span>"
            )
        })
        .collect::<String>();

        let mut h = String::new();
        h.push_str("<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">");
        h.push_str("<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">");
        h.push_str("<title>preflightx · single-cell read structures</title>");
        // No webfonts / CDN: the doc is fully self-contained and renders offline, in
        // keeping with preflightx's zero-network default. The CSS font stacks fall back
        // to system serif/mono (Georgia, system monospace).
        h.push_str(&format!(
            "<style>{css}</style></head><body><div class=\"wrap\">"
        ));
        h.push_str("<h1 class=\"title\">Single-cell<br>read <em>structures</em></h1>");
        h.push_str(&format!(
            "<p class=\"lede\">Every assay in the preflightx catalogue, drawn straight from its \
             <a href=\"https://github.com/pachterlab/seqspec\">seqspec</a> description. Each card shows \
             <strong>every sequencing read</strong> (R1, R2, …), left to right, 5′→3′, after stripping the \
             library primer that isn't sequenced; the read carrying the cell barcode is flagged. {} assays \
             with an in-read barcode; {} plate-based assays carry no in-read cell barcode (cell identity is \
             the well). Hover a block for its region id.</p>",
            cards.len(),
            plate_cards.len() + plate.len(),
        ));
        h.push_str(&format!(
            "<div class=\"legend\"><b>legend</b>{legend}</div>"
        ));

        for (i, (_, c)) in cards.iter().enumerate() {
            let delay = (i.min(24) * 40) as u32;
            h.push_str(&c.replacen(
                "class=\"card\"",
                &format!("class=\"card\" style=\"animation-delay:{delay}ms\""),
                1,
            ));
        }

        h.push_str("<div class=\"section\">Plate-based: no in-read cell barcode</div>");
        for (i, c) in plate_cards.iter().enumerate() {
            let delay = (i.min(12) * 40) as u32;
            h.push_str(&c.replacen(
                "class=\"card plate-card\"",
                &format!("class=\"card plate-card\" style=\"animation-delay:{delay}ms\""),
                1,
            ));
        }
        if !plate.is_empty() {
            h.push_str("<div class=\"plate\">");
            for p in &plate {
                h.push_str(&format!("<span>{p}</span>"));
            }
            h.push_str("</div>");
        }

        h.push_str(
            "<p class=\"foot\">Generated from the vendored seqspec catalogue by \
             <code>cargo test gen_assay_html -- --ignored</code>. Block widths are proportional to \
             length (clamped for very long/short segments); leading primer-binding regions are omitted \
             because they are not sequenced. Geometries follow each spec faithfully.</p>",
        );
        h.push_str("</div></body></html>");

        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/docs/assay-structures.html");
        std::fs::create_dir_all(concat!(env!("CARGO_MANIFEST_DIR"), "/docs")).unwrap();
        std::fs::write(path, h).unwrap();
        eprintln!("wrote {path}");
    }

    #[test]
    fn multimodal_prefers_rna_barcode_read() {
        // 10x multiome has an RNA read (CB16+UMI12=28) and an ATAC read (CB16 only);
        // structural_read must prefer the RNA cell-barcode read.
        assert_eq!(geom_of("10x_rna_atac_multi.spec.yaml"), Some(28));
        assert_eq!(geom_of("dogmaseq_lll.spec.yaml"), Some(28));
    }

    #[test]
    fn trailing_adapter_dropped_but_inter_barcode_adapter_kept() {
        // quartz-seq's read-2 adapter sits AFTER the cell barcode (loose primer_id) and
        // must be dropped -> 6 bp cell index. scifi-RNA-seq's adapters sit BETWEEN the
        // bead and RT barcodes and are sequenced -> kept in the 84 bp span.
        assert_eq!(geom_of("quartz_seq.spec.yaml"), Some(6));
        assert_eq!(geom_of("scifi_rna_seq.spec.yaml"), Some(84));
    }

    #[test]
    fn every_spec_is_covered_or_plate_based() {
        // The completeness guarantee: each catalogued assay either yields a barcode-read
        // geometry, or has no in-read cell barcode at all (plate-based). Nothing should
        // have a barcode region we fail to resolve into a geometry.
        let mut failed: Vec<String> = Vec::new();
        let mut covered = 0;
        let mut plate = 0;
        for a in embedded_assays() {
            let has_bc_region = a
                .leaf_regions()
                .iter()
                .any(|r| is_barcode_type(&r.region_type));
            let geom = a.structural_read().and_then(|(_, l)| observed_read_len(&l));
            match (geom, has_bc_region) {
                (Some(_), _) => covered += 1,
                (None, false) => plate += 1, // genuinely plate-based — correct
                (None, true) => failed.push(a.spec_file.clone()),
            }
        }
        assert!(
            failed.is_empty(),
            "{} spec(s) have a barcode region but no resolved geometry: {failed:?}",
            failed.len()
        );
        // Sanity: the catalogue is mostly droplet/combinatorial (in-read barcodes), with
        // a handful of plate-based assays.
        assert!(covered >= 45, "only {covered} specs covered");
        assert!(
            (1..=12).contains(&plate),
            "unexpected plate-based count {plate}"
        );
    }

    #[test]
    fn variable_and_synonym_barcodes_are_covered() {
        // inDrop's variable-length BC1 and BD Rhapsody's combinatorial CLS barcodes
        // (with a leading variable spacer) both resolve — exercising variable-length
        // barcodes, synonym region_types, and the linker-aware terminator.
        let g = sc_geometries();
        for spec in [
            "indropv2.spec.yaml",
            "bd_rhapsody_eb.spec.yaml",
            "split_seq.spec.yaml",
        ] {
            assert!(
                g.iter().any(|x| x.spec_file == spec),
                "{spec} should have a geometry"
            );
        }
    }

    #[test]
    fn observed_geometry_strips_leading_primer() {
        // Drop-seq R1's structural read includes the 39 bp Read1 custom primer before
        // the barcode (seqspec models it as a region after primer_id=illumina_p5);
        // the observed barcode read is the 20 bp (CB12+UMI8) after it.
        let drop = sc_geometries()
            .iter()
            .find(|g| g.spec_file == "drop_seq.spec.yaml")
            .expect("drop_seq present");
        assert_eq!(
            drop.read_len, 20,
            "drop_seq observed read should be CB12+UMI8"
        );
        // 10x v3: clean R1, no leading primer -> 16 + 12 = 28.
        let v3 = sc_geometries()
            .iter()
            .find(|g| g.spec_file == "10x_rna_v3.spec.yaml")
            .expect("v3 present");
        assert_eq!(v3.read_len, 28);
    }

    #[test]
    fn sc_geometries_matching_finds_known_lengths() {
        assert!(sc_geometries_matching(20, 1)
            .iter()
            .any(|g| g.spec_file == "drop_seq.spec.yaml"));
        assert!(sc_geometries_matching(28, 1)
            .iter()
            .any(|g| g.spec_file == "10x_rna_v3.spec.yaml"));
        assert!(sc_geometries_matching(0, 1).is_empty());
        assert!(sc_geometries_matching(1000, 1).is_empty());
    }

    #[test]
    fn parses_10x_v3_and_resolves_r1_layout() {
        let yaml = include_str!("../vendor/seqspec/docs/examples/assays/10x_rna_v3.spec.yaml");
        let assay = Assay::from_yaml(yaml).unwrap();
        assert_eq!(assay.assay_id, "10xRNAv3");
        assert_eq!(assay.modalities, vec!["rna"]);

        // R1 primes from truseq_read1 (pos strand), covering CB(16) then UMI(12).
        let layout = assay.read_layout("rna-R1.fastq.gz").unwrap();
        let kinds: Vec<(&str, usize)> = layout
            .iter()
            .map(|s| (s.region_type.as_str(), s.len))
            .collect();
        assert_eq!(kinds, vec![("barcode", 16), ("umi", 12)]);
    }

    #[test]
    fn variable_cdna_region_terminates_the_layout() {
        // SPLiT-seq R1 primes from s5 and covers ME1 then cDNA. cDNA is variable
        // (max_len > min_len) and absorbs the rest of the read, so NO barcode
        // segments may leak into the cDNA read's layout.
        let yaml = include_str!("../vendor/seqspec/docs/examples/assays/split_seq.spec.yaml");
        let assay = Assay::from_yaml(yaml).unwrap();
        let layout = assay.read_layout("R1.fastq.gz").unwrap();
        assert!(
            !layout.iter().any(|s| s.region_type == "barcode"),
            "cDNA read must not contain barcode segments, got {:?}",
            layout.iter().map(|s| &s.region_type).collect::<Vec<_>>()
        );
        // cDNA is the terminal segment.
        assert_eq!(layout.last().unwrap().region_type, "cDNA");
    }

    #[test]
    fn embedded_catalogue_loads_many_assays() {
        let assays = embedded_assays();
        // ~50 vendored specs; every one must parse cleanly.
        assert_eq!(assays.len(), EMBEDDED_SPECS.len());
        assert!(assays.len() >= 40, "only {} assays parsed", assays.len());
        assert!(assays
            .iter()
            .any(|a| a.assay_id.to_lowercase().contains("10x")));
        assert!(assays
            .iter()
            .any(|a| a.assay_id.to_lowercase().contains("split")
                || a.name.to_lowercase().contains("split")));
    }

    #[test]
    fn resolves_combinatorial_neg_strand_layout_splitseq() {
        // SPLiT-seq R2 primes from truseq_read2 on the neg strand and reads backward:
        // UMI(10) + barcode-3(8) + linker(30) + barcode-2(8) + linker(30) + barcode-1(8) = 94 bp.
        let yaml = include_str!("../vendor/seqspec/docs/examples/assays/split_seq.spec.yaml");
        let assay = Assay::from_yaml(yaml).unwrap();
        let layout = assay.read_layout("R2.fastq.gz").unwrap();
        let kinds: Vec<(&str, usize)> = layout
            .iter()
            .map(|s| (s.region_type.as_str(), s.len))
            .collect();
        assert_eq!(
            kinds,
            vec![
                ("umi", 10),
                ("barcode", 8),
                ("linker", 30),
                ("barcode", 8),
                ("linker", 30),
                ("barcode", 8),
            ]
        );
        // Linkers carry a fixed sequence (reverse-complemented for the neg strand).
        let linker = layout.iter().find(|s| s.region_type == "linker").unwrap();
        assert!(linker.sequence.is_some());
    }
}
