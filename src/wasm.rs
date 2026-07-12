//! Browser entry points (wasm-bindgen). Uploaded bytes go into [`fastq::vfs`]; the analysis
//! functions run unmodified against a synthetic path.
//!
//! No whitelist cache here. The front-end calls [`needed_whitelists`], fetches each, and
//! stashes it via [`load_whitelist`]; without them chemistry caps at medium confidence.

use std::path::Path;
use wasm_bindgen::prelude::*;

use crate::{assaydetect, chemistry, composition, fastq, report, sex, sexpanel, size};

const R1: &str = "r1";
const R2: &str = "r2";
pub(crate) const WL_DIR: &str = "wl";

/// Routes Rust panics to the JS console; runs inside `init()`.
#[wasm_bindgen(start)]
pub fn start() {
    console_error_panic_hook::set_once();
}

fn to_js<T: serde::Serialize>(value: &T) -> Result<String, JsValue> {
    serde_json::to_string(value).map_err(|e| JsValue::from_str(&e.to_string()))
}

fn err_to_js(e: anyhow::Error) -> JsValue {
    JsValue::from_str(&format!("{e:#}"))
}

/// Some(cache) once any onlist has been stashed via [`load_whitelist`].
fn wl_cache() -> Option<&'static Path> {
    if fastq::vfs::has_prefix(&format!("{WL_DIR}/")) {
        Some(Path::new(WL_DIR))
    } else {
        None
    }
}

/// Drops the read buffers after `f`; loaded onlists survive.
fn with_pair<T>(r1: Vec<u8>, r2: Vec<u8>, f: impl FnOnce() -> T) -> T {
    fastq::vfs::put(R1, r1);
    fastq::vfs::put(R2, r2);
    let out = f();
    fastq::vfs::remove(R1);
    fastq::vfs::remove(R2);
    out
}

/// Whitelist filenames matching the pair's R1 geometry. The front-end fetches each and
/// passes it back through [`load_whitelist`].
#[wasm_bindgen]
pub fn needed_whitelists(r1: Vec<u8>, r2: Vec<u8>, n_reads: usize) -> Result<String, JsValue> {
    let files = with_pair(r1, r2, || {
        chemistry::candidate_whitelist_files(Path::new(R1), Path::new(R2), n_reads)
    });
    to_js(&files.map_err(err_to_js)?)
}

/// Stash a fetched onlist (`.txt.gz` or plain). Cached across runs, so fetched once.
#[wasm_bindgen]
pub fn load_whitelist(name: &str, bytes: Vec<u8>) {
    fastq::vfs::put(&format!("{WL_DIR}/{name}"), bytes);
}

/// Inject a k-mer panel; no panel cache dir exists on wasm, so without this the sc-vs-sn
/// and 3'/5' verdicts stay `None`. `kind` is "endpanel" or "nucpanel". Call before [`check`].
#[wasm_bindgen]
pub fn load_panel(kind: &str, bytes: Vec<u8>) {
    crate::kmerpanel::store_panel(kind, String::from_utf8_lossy(&bytes).into_owned());
}

/// Onlists `detect_assay` can score. The front-end fetches each `url` and stashes it via
/// [`load_whitelist`] under `cache_name`.
#[wasm_bindgen]
pub fn onlist_sources() -> Result<String, JsValue> {
    to_js(&assaydetect::onlist_fetch_list())
}

/// Empty `r2` means single-end; `f` gets a flag for whether R2 is present.
fn with_reads<T>(r1: Vec<u8>, r2: Vec<u8>, f: impl FnOnce(bool) -> T) -> T {
    fastq::vfs::put(R1, r1);
    let has_r2 = !r2.is_empty();
    if has_r2 {
        fastq::vfs::put(R2, r2);
    }
    let out = f(has_r2);
    fastq::vfs::remove(R1);
    if has_r2 {
        fastq::vfs::remove(R2);
    }
    out
}

/// Full report. `r1`/`r2` are raw uploaded bytes (gzip auto-detected); empty `r2` is
/// single-end.
#[wasm_bindgen]
pub fn check(r1: Vec<u8>, r2: Vec<u8>, n_reads: usize) -> Result<String, JsValue> {
    with_reads(r1, r2, |has_r2| {
        if has_r2 {
            run_check(n_reads)
        } else {
            run_check_single(n_reads)
        }
    })
}

/// Chemistry and assay detection need both mates, so they are skipped here.
fn run_check_single(n: usize) -> Result<String, JsValue> {
    let r1 = Path::new(R1);
    let bulk = size::bulk_or_sc(r1, n).map_err(err_to_js)?;
    let sz = size::size_sanity(r1, None, n).map_err(err_to_js)?;
    let comp = composition::composition(r1, n).map_err(err_to_js)?;
    let (nucleus, end) = report::panel_verdicts(r1, None, n, None);
    let classification = report::Classification::derive(&None, &bulk, &nucleus);
    let full = report::FullReport {
        r1: R1.to_string(),
        r2: None,
        chemistry: None,
        assay_detection: None,
        bulk_vs_sc: bulk,
        size: sz,
        composition: comp,
        species: None,
        sex: None,
        nucleus,
        end,
        classification,
    };
    to_js(&full)
}

fn run_check(n: usize) -> Result<String, JsValue> {
    let r1 = Path::new(R1);
    let r2 = Path::new(R2);
    let (chem, assay_detection) =
        report::FullReport::chem_with_fallback(r1, r2, n, wl_cache(), None).map_err(err_to_js)?;
    let bulk = size::bulk_or_sc(r1, n).map_err(err_to_js)?;
    let sz = size::size_sanity(r1, Some(r2), n).map_err(err_to_js)?;
    let cdna = fastq::longer_read(r1, r2, 256).map_err(err_to_js)?;
    let comp = composition::composition(cdna, n).map_err(err_to_js)?;
    let (nucleus, end) = report::panel_verdicts(r1, Some(r2), n, None);
    let chemistry = Some(chem);
    let classification = report::Classification::derive(&chemistry, &bulk, &nucleus);
    let full = report::FullReport {
        r1: R1.to_string(),
        r2: Some(R2.to_string()),
        chemistry,
        assay_detection,
        bulk_vs_sc: bulk,
        size: sz,
        composition: comp,
        species: None,
        sex: None,
        nucleus,
        end,
        classification,
    };
    to_js(&full)
}

/// Generic seqspec-driven assay ranking (geometry + linker anchors).
#[wasm_bindgen]
pub fn detect_assay(r1: Vec<u8>, r2: Vec<u8>, n_reads: usize) -> Result<String, JsValue> {
    let result = with_pair(r1, r2, || {
        assaydetect::detect(Path::new(R1), Path::new(R2), n_reads).map_err(err_to_js)
    });
    to_js(&result?)
}

/// 10x chemistry detection only. Uses a stashed onlist when one is loaded.
#[wasm_bindgen]
pub fn chem(r1: Vec<u8>, r2: Vec<u8>, n_reads: usize) -> Result<String, JsValue> {
    let result = with_pair(r1, r2, || {
        chemistry::detect(Path::new(R1), Path::new(R2), n_reads, wl_cache(), None)
            .map_err(err_to_js)
    });
    to_js(&result?)
}

/// Species with a bundled sex panel; backs the front-end's selector. A non-XY_XIST
/// mechanism is not callable yet.
#[wasm_bindgen]
pub fn sex_species() -> Result<String, JsValue> {
    let panels: Vec<serde_json::Value> = sexpanel::species_list()
        .iter()
        .filter_map(|s| sexpanel::for_species(s))
        .map(|p| {
            serde_json::json!({
                "species": p.species,
                "status": p.status,
                "mechanism": p.mechanism,
            })
        })
        .collect();
    to_js(&panels)
}

/// Genetic-sex inference (Y-gene + XIST k-mers). Separate from [`check`] so the front-end
/// runs it on request. Y/XIST are sparse, so `sex_reads` samples deeper (native 200000).
#[wasm_bindgen]
pub fn sex(r1: Vec<u8>, r2: Vec<u8>, species: &str, sex_reads: usize) -> Result<String, JsValue> {
    let result = with_reads(r1, r2, |has_r2| {
        let r2_path = has_r2.then(|| Path::new(R2));
        sex::detect_sex_auto(Path::new(R1), r2_path, species, sex_reads).map_err(err_to_js)
    });
    to_js(&result?)
}

/// Sex against a front-end-supplied .spanel (see web/data/panels.json), which reaches any
/// catalogued species, not just the two embedded ones. Unusable panel gives a no-call.
#[wasm_bindgen]
pub fn sex_with_panel(
    r1: Vec<u8>,
    r2: Vec<u8>,
    species: &str,
    spanel_text: &str,
    sex_reads: usize,
) -> Result<String, JsValue> {
    let parsed = sexpanel::parse_spanel(spanel_text);
    let panel = parsed.as_ref().map(|(p, _)| p);
    let result: Result<_, JsValue> = with_reads(r1, r2, |has_r2| {
        let r2_path = has_r2.then(|| Path::new(R2));
        let reads = sex::sample_cdna(Path::new(R1), r2_path, sex_reads).map_err(err_to_js)?;
        Ok(sex::report_for_reads(panel, species, &reads))
    });
    to_js(&result?)
}

/// Same renderer as the CLI's `check --html`.
#[wasm_bindgen]
pub fn report_html_json(json: &str) -> Result<String, JsValue> {
    let v: serde_json::Value =
        serde_json::from_str(json).map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(crate::reporthtml::render(&v))
}

/// The embedded seqspec assay catalogue (no reads required).
#[wasm_bindgen]
pub fn assays() -> Result<String, JsValue> {
    to_js(&assaydetect::catalogue())
}

/// Library structure of a catalogued assay (by name or spec file).
#[wasm_bindgen]
pub fn describe(assay: &str) -> Result<String, JsValue> {
    match assaydetect::find_assay(assay) {
        Some(a) => Ok(assaydetect::describe(&a)),
        None => Err(JsValue::from_str(&format!(
            "no catalogued assay matching '{assay}'"
        ))),
    }
}

/// Structured library description matching CLI `describe --json`.
#[wasm_bindgen]
pub fn describe_json(assay: &str) -> Result<String, JsValue> {
    match assaydetect::find_assay(assay) {
        Some(a) => to_js(&assaydetect::describe_json(&a)),
        None => Err(JsValue::from_str(&format!(
            "no catalogued assay matching {assay}"
        ))),
    }
}
