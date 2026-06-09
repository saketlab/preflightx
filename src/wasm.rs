//! Browser entry points (wasm-bindgen).
//!
//! Each function mirrors a CLI subcommand: it loads the uploaded FASTQ bytes
//! into the in-memory [`fastq::vfs`], invokes the *unchanged* `&Path`-based
//! analysis function with a synthetic path, and returns the report as a JSON
//! string (the same shape the CLI emits under `--json`). Errors surface as a
//! JS exception carrying the anyhow message.
//!
//! 10x whitelist disambiguation works when the host serves the onlists: the
//! front-end calls [`needed_whitelists`] to learn which file matches the read
//! geometry, fetches it, hands the bytes to [`load_whitelist`], then calls
//! [`check`]. Without those files chemistry falls back to geometry + TSO and
//! confidence is capped at `medium`.

use std::path::Path;
use wasm_bindgen::prelude::*;

use crate::{assaydetect, chemistry, composition, fastq, report, size};

const R1: &str = "r1";
const R2: &str = "r2";
const WL_DIR: &str = "wl";

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

/// Stash R1/R2 in the vfs, run `f` against the synthetic paths, then drop the read
/// buffers (loaded onlists persist across runs).
fn with_pair<T>(r1: &[u8], r2: &[u8], f: impl FnOnce() -> T) -> T {
    fastq::vfs::put(R1, r1.to_vec());
    fastq::vfs::put(R2, r2.to_vec());
    let out = f();
    fastq::vfs::remove(R1);
    fastq::vfs::remove(R2);
    out
}

/// Whitelist filenames matching the pair's R1 geometry, as a JSON array. The
/// front-end fetches each and passes it back through [`load_whitelist`].
#[wasm_bindgen]
pub fn needed_whitelists(r1: &[u8], r2: &[u8], n_reads: usize) -> Result<String, JsValue> {
    let files = with_pair(r1, r2, || {
        chemistry::candidate_whitelist_files(Path::new(R1), Path::new(R2), n_reads)
    });
    to_js(&files.map_err(err_to_js)?)
}

/// Stash a fetched onlist (raw `.txt.gz` or plain bytes) for later [`check`]s.
/// Cached across runs so it is fetched once per session.
#[wasm_bindgen]
pub fn load_whitelist(name: &str, bytes: &[u8]) {
    fastq::vfs::put(&format!("{WL_DIR}/{name}"), bytes.to_vec());
}

/// Full report (chemistry + bulk-vs-sc + size + composition) for an R1/R2 pair.
/// `r1`/`r2` are the raw uploaded file bytes (gzip is auto-detected).
#[wasm_bindgen]
pub fn check(r1: &[u8], r2: &[u8], n_reads: usize) -> Result<String, JsValue> {
    with_pair(r1, r2, || run_check(n_reads))
}

fn run_check(n: usize) -> Result<String, JsValue> {
    let r1 = Path::new(R1);
    let r2 = Path::new(R2);
    let chem = chemistry::detect(r1, r2, n, wl_cache()).map_err(err_to_js)?;
    let bulk = size::bulk_or_sc(r1, n).map_err(err_to_js)?;
    let sz = size::size_sanity(r1, Some(r2), n).map_err(err_to_js)?;
    let comp = composition::composition(r2, n).map_err(err_to_js)?;
    let full = report::FullReport {
        r1: R1.to_string(),
        r2: Some(R2.to_string()),
        chemistry: Some(chem),
        bulk_vs_sc: bulk,
        size: sz,
        composition: comp,
    };
    to_js(&full)
}

/// Generic seqspec-driven assay ranking (geometry + linker anchors).
#[wasm_bindgen]
pub fn detect_assay(r1: &[u8], r2: &[u8], n_reads: usize) -> Result<String, JsValue> {
    let result = with_pair(r1, r2, || {
        assaydetect::detect(Path::new(R1), Path::new(R2), n_reads).map_err(err_to_js)
    });
    to_js(&result?)
}

/// 10x chemistry detection only. Uses a stashed onlist when one is loaded.
#[wasm_bindgen]
pub fn chem(r1: &[u8], r2: &[u8], n_reads: usize) -> Result<String, JsValue> {
    let result = with_pair(r1, r2, || {
        chemistry::detect(Path::new(R1), Path::new(R2), n_reads, wl_cache()).map_err(err_to_js)
    });
    to_js(&result?)
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
