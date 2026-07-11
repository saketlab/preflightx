use anyhow::Result;
use serde::Serialize;
use std::collections::{HashMap, HashSet};
use std::io::BufRead;
use std::path::Path;

use crate::assay::{AnchorMotif, AssaySpec, ASSAYS, TSO_PRESENT};
use crate::fastq;
use crate::stats::{bernoulli_llr, dominant_len};

// TSO background rate: fraction of R2 reads where the 13-mer appears by chance.
const P_TSO_BACKGROUND: f64 = 0.01;
// Expected CB match rates when the whitelist is correct vs wrong.
const P_WL_CORRECT: f64 = 0.70;
const P_WL_NULL: f64 = 0.01;

// Thresholds are per-read (LLR / n_reads) so they are sample-size-independent.
const HIGH_LL_PER_READ: f64 = 0.30; // strong total evidence
const HIGH_GAP_PER_READ: f64 = 0.10; // top assay uniquely ahead of second
const MEDIUM_LL_PER_READ: f64 = 0.001; // any positive discrimination

// SPRT early-stopping (no-whitelist path only; whitelist needs a large CB sample).
const SPRT_THRESHOLD: f64 = 200.0; // nats; |TSO LLR| above this decides the family
const SPRT_MIN_READS: usize = 5_000; // never stop before this many pairs
const SPRT_BATCH: usize = 2_000; // pairs per streaming batch

#[derive(Serialize, Debug)]
pub struct WhitelistScore {
    pub filename: String,
    pub hits: usize,
    pub pct: f64,
}

#[derive(Serialize, Debug)]
pub struct AnchorScore {
    pub motif: String,
    pub hit_pct: f64,
}

#[derive(Serialize, Debug)]
pub struct ToolHint {
    pub tool: String,
    pub hint: String,
}

/// One reasoning step in the verdict's evidence ledger: a signal, what was observed,
/// the nats of log-likelihood it contributed to the winner, and what it implies.
#[derive(Serialize, Debug)]
pub struct Evidence {
    pub signal: String,
    pub observation: String,
    /// Nats of log-likelihood this signal added to the winner (None for non-LLR steps).
    pub llr_contribution: Option<f64>,
    pub conclusion: String,
}

#[derive(Serialize, Debug)]
pub struct ChemReport {
    pub r1_length_dominant: usize,
    pub r1_length_distribution: Vec<(usize, usize)>,
    pub reads_sampled: usize,
    pub assay_id: String,
    pub assay_display: String,
    pub cb_read: String,
    pub cb_offset: usize,
    pub cb_len: usize,
    pub umi_offset: usize,
    pub umi_len: usize,
    pub n_cb_candidates: usize,
    /// False when the reads show no single-cell signal (bulk/genomic): no
    /// barcode-length R1, no TSO, no whitelist overlap. When false, `assay_id`
    /// is `"none"` and the call should be read as "not single-cell".
    pub is_single_cell: bool,
    /// True when sampling stopped early (SPRT decided), so `reads_sampled < n` is
    /// expected rather than a truncated input.
    pub early_stop: bool,
    pub whitelist: Option<String>,
    pub whitelist_scores: Vec<WhitelistScore>,
    pub anchor_scores: Vec<AnchorScore>,
    pub tool_hints: Vec<ToolHint>,
    pub confidence: String,
    /// ATAC-vs-RNA read-content modality: "atac" | "rna" | "unknown". "atac" overrides the
    /// RNA chemistry call (assay_id "atac", is_single_cell false).
    pub modality: String,
    /// promoter / (promoter + exon) k-mer hits; high -> ATAC. 0.0 when no panel/hits.
    pub atac_promoter_frac: f64,
    /// Transcript-end polarity from the 3'/5' k-mer panel, independent of R1 geometry
    /// (which cannot tell 3'v2 from 5'): "3'" | "5'" | "full-length" | "ambiguous" | "unknown".
    pub end: String,
    pub end_confidence: String,
    /// Reads hitting a 5'-end / 3'-end transcript marker, and r = n3/(n3+n5).
    pub n5: usize,
    pub n3: usize,
    pub end_ratio: f64,
    pub interpretation: Vec<Evidence>,
    pub notes: Vec<String>,
}

pub fn detect(
    r1: &Path,
    r2: &Path,
    n: usize,
    whitelist_cache: Option<&Path>,
    species: Option<&str>,
) -> Result<ChemReport> {
    // orient + sample + load whitelists (I/O) → analyze (pure) → build_report (pure).
    let (eff_r1, eff_r2, swapped) = orient(r1, r2)?;
    let wl_sets = match whitelist_cache {
        Some(cache) => {
            let (preview, _) = fastq::sample(eff_r1, n.min(2_000))?;
            let dominant = dominant_len(&preview);
            let candidates = candidates_for(dominant);
            load_candidate_whitelists(cache, &candidates)?
        }
        None => HashMap::new(),
    };
    // Full sample (no SPRT early-stop) when a whitelist cache or a k-mer panel is
    // available; all need depth (CB overlap, 3'/5' + ATAC housekeeping hits).
    let want_full =
        !wl_sets.is_empty() || crate::endpanel::have_panels() || crate::atacpanel::have_panels();
    let s = sample_reads(eff_r1, eff_r2, n, want_full)?;

    let analysis = analyze(&s.r1, &s.r2, &wl_sets, s.sprt_tso_hits);
    let end = crate::endpanel::classify(&s.r2, species);
    let atac = crate::atacpanel::classify_pair(&s.r1, &s.r2, species);
    Ok(build_report(analysis, swapped, n, s.early_stop, end, atac))
}

/// Distinct whitelist filenames relevant to the pair's observed R1 geometry.
/// The wasm front-end calls this to fetch only the matching onlist before
/// `detect`, instead of pulling every catalogued whitelist.
pub fn candidate_whitelist_files(r1: &Path, r2: &Path, n: usize) -> Result<Vec<&'static str>> {
    let (eff_r1, _, _) = orient(r1, r2)?;
    let (sample, _) = fastq::sample(eff_r1, n.min(2_000))?;
    let dominant = dominant_len(&sample);
    let mut files = Vec::new();
    for spec in candidates_for(dominant) {
        if let Some(wl) = &spec.whitelist {
            if !files.contains(&wl.filename) {
                files.push(wl.filename);
            }
        }
    }
    Ok(files)
}

/// Quick-sample both inputs to decide which file carries the barcode read. If R1
/// looks long and R2 looks like a barcode read, treat them as swapped. Returns the
/// effective (r1, r2) ordering and whether a swap was applied.
fn orient<'a>(r1: &'a Path, r2: &'a Path) -> Result<(&'a Path, &'a Path, bool)> {
    let (a_quick, b_quick) =
        fastq::par_pair(|| fastq::sample_seqs(r1, 1_000), || fastq::sample_seqs(r2, 1_000))?;
    let swapped =
        !fits_barcode_len(dominant_len(&a_quick)) && fits_barcode_len(dominant_len(&b_quick));
    Ok(if swapped {
        (r2, r1, true)
    } else {
        (r1, r2, false)
    })
}

/// Sampled reads plus how sampling ended: the running TSO count (SPRT path only)
/// and whether SPRT stopped early because the TSO signal was already decisive.
struct Sampled {
    r1: Vec<Vec<u8>>,
    r2: Vec<Vec<u8>>,
    sprt_tso_hits: Option<usize>,
    early_stop: bool,
}

/// Sample up to `n` read pairs. With a whitelist cache, collect the full `n` (the CB
/// set must be large to score overlap); without one, SPRT stops as soon as the TSO
/// signal is decisive.
fn sample_reads(eff_r1: &Path, eff_r2: &Path, n: usize, have_cache: bool) -> Result<Sampled> {
    if have_cache {
        let mut r1 = Vec::with_capacity(n);
        let mut r2 = Vec::with_capacity(n);
        fastq::stream_paired_batched(eff_r1, eff_r2, n, 10_000, |mut b1, mut b2| {
            r1.append(&mut b1);
            r2.append(&mut b2);
            true
        })?;
        Ok(Sampled {
            r1,
            r2,
            sprt_tso_hits: None,
            early_stop: false,
        })
    } else {
        let mut r1: Vec<Vec<u8>> = Vec::new();
        let mut r2: Vec<Vec<u8>> = Vec::new();
        let mut tso_hits: usize = 0;
        let mut early_stop = false;
        fastq::stream_paired_batched(eff_r1, eff_r2, n, SPRT_BATCH, |mut b1, mut b2| {
            tso_hits += b2.iter().filter(|s| anchor_hit(s, &TSO_PRESENT)).count();
            r1.append(&mut b1);
            r2.append(&mut b2);
            if r1.len() < SPRT_MIN_READS {
                return true;
            }
            let llr = bernoulli_llr(tso_hits, r1.len(), TSO_PRESENT.p_signal, P_TSO_BACKGROUND);
            let keep_going = llr.abs() <= SPRT_THRESHOLD;
            early_stop = !keep_going;
            keep_going
        })?;
        Ok(Sampled {
            r1,
            r2,
            sprt_tso_hits: Some(tso_hits),
            early_stop,
        })
    }
}

/// Load the whitelist barcode sets needed to score `candidates`, once per unique file.
/// Files that are absent from the cache are silently skipped (scored as no overlap).
fn load_candidate_whitelists(
    cache: &Path,
    candidates: &[&'static AssaySpec],
) -> Result<HashMap<&'static str, HashSet<Vec<u8>>>> {
    let mut sets: HashMap<&'static str, HashSet<Vec<u8>>> = HashMap::new();
    for spec in candidates {
        if let Some(wl) = &spec.whitelist {
            if sets.contains_key(wl.filename) {
                continue;
            }
            let path = cache.join(wl.filename);
            if let Some(set) = load_whitelist_opt(&path, wl.cb_len)? {
                sets.insert(wl.filename, set);
            }
        }
    }
    Ok(sets)
}

// pure scoring core (the testable seam)

/// The ranked outcome of scoring sampled reads against the assay registry. Holds
/// everything `build_report` needs; produced by `analyze` with no I/O.
pub(crate) struct Analysis {
    dominant: usize,
    len_dist: Vec<(usize, usize)>,
    candidates: Vec<&'static AssaySpec>,
    cb_sets: HashMap<(usize, usize), HashSet<Vec<u8>>>,
    anchor_pcts: HashMap<&'static [u8], f64>,
    wl_scores: HashMap<&'static str, (usize, f64)>,
    scored: Vec<(&'static AssaySpec, f64)>,
    confidence: &'static str,
    n_reads: usize,
    /// TSO hit rate (%) in R2 and the best whitelist overlap (%), kept so the gate
    /// and the report's notes share one computation.
    tso_pct: f64,
    max_wl_pct: f64,
    /// Whether the reads carry any single-cell signal at all (a barcode-length R1,
    /// a TSO, or whitelist overlap). False for bulk/genomic libraries, where the
    /// scored "winner" is a meaningless force-fit and must not be reported as a call.
    is_single_cell: bool,
}

/// Candidate assays for an observed R1 length: exact match, or +/-2 bp for
/// non-standard cycle counts. Falls back to the whole registry if none fit.
fn candidates_for(dominant: usize) -> Vec<&'static AssaySpec> {
    let fuzzy: Vec<_> = ASSAYS
        .iter()
        .filter(|s| s.r1_expected_len().abs_diff(dominant) <= 2)
        .collect();
    if fuzzy.is_empty() {
        ASSAYS.iter().collect()
    } else {
        fuzzy
    }
}

/// R1 length histogram: the dominant length and the full distribution (descending by count).
fn len_profile(seqs: &[Vec<u8>]) -> (usize, Vec<(usize, usize)>) {
    let mut hist: HashMap<usize, usize> = HashMap::new();
    for s in seqs {
        *hist.entry(s.len()).or_insert(0) += 1;
    }
    let mut dist: Vec<(usize, usize)> = hist.into_iter().collect();
    dist.sort_by(|a, b| b.1.cmp(&a.1));
    let dominant = dist.first().map(|(l, _)| *l).unwrap_or(0);
    (dominant, dist)
}

/// Pick the confidence bucket from the per-read evidence: `high` needs strong total
/// evidence, a clear lead over the runner-up, *and* a whitelist; `medium` needs any
/// positive discrimination (or a single unambiguous candidate); else `low`.
fn classify_confidence(
    scored: &[(&'static AssaySpec, f64)],
    candidates: &[&'static AssaySpec],
    dominant: usize,
    n_reads: usize,
    has_wl: bool,
) -> &'static str {
    let best = scored.first().map(|(_, s)| *s).unwrap_or(0.0);
    let second = scored.get(1).map(|(_, s)| *s).unwrap_or(0.0);
    // max(1) guards against 0/0 = NaN when no reads were parsed.
    let denom = n_reads.max(1) as f64;
    let ll_per_read = best / denom;
    let gap_per_read = (best - second) / denom;
    if ll_per_read > HIGH_LL_PER_READ && gap_per_read > HIGH_GAP_PER_READ && has_wl {
        "high"
    } else if ll_per_read > MEDIUM_LL_PER_READ || (candidates.len() == 1 && dominant > 0) {
        "medium"
    } else {
        "low"
    }
}

/// Score the sampled reads against every candidate assay. Pure: no file or network
/// access — whitelist barcode sets are passed in via `wl_sets`.
fn analyze(
    r1_seqs: &[Vec<u8>],
    r2_seqs: &[Vec<u8>],
    wl_sets: &HashMap<&'static str, HashSet<Vec<u8>>>,
    sprt_tso_hits: Option<usize>,
) -> Analysis {
    let n_reads = r1_seqs.len();
    let (dominant, len_dist) = len_profile(r1_seqs);
    let candidates = candidates_for(dominant);

    // CB candidates keyed by geometry, shared across assays with the same (offset, len).
    let mut cb_sets: HashMap<(usize, usize), HashSet<Vec<u8>>> = HashMap::new();
    for spec in &candidates {
        let key = (spec.r1_cb_offset, spec.r1_cb_len);
        cb_sets
            .entry(key)
            .or_insert_with(|| extract_cb_candidates(r1_seqs, spec.r1_cb_offset, spec.r1_cb_len));
    }

    // Whitelist overlap, computed once per unique whitelist file present in `wl_sets`.
    let mut wl_scores: HashMap<&'static str, (usize, f64)> = HashMap::new();
    for spec in &candidates {
        if let Some(wl) = &spec.whitelist {
            if wl_scores.contains_key(wl.filename) {
                continue;
            }
            if let Some(wl_set) = wl_sets.get(wl.filename) {
                let cbs = cb_sets
                    .get(&(spec.r1_cb_offset, spec.r1_cb_len))
                    .expect("cb_sets populated for all candidate geometries");
                let hits = cbs
                    .iter()
                    .filter(|cb| wl_set.contains(cb.as_slice()))
                    .count();
                let pct = if cbs.is_empty() {
                    0.0
                } else {
                    100.0 * hits as f64 / cbs.len() as f64
                };
                wl_scores.insert(wl.filename, (hits, pct));
            }
        }
    }

    // R2 anchor scan, once per unique motif. Seed the TSO motif from the SPRT count
    // (if taken) to avoid re-walking every R2 read.
    let mut anchor_pcts: HashMap<&'static [u8], f64> = HashMap::new();
    if let Some(hits) = sprt_tso_hits {
        let pct = if n_reads == 0 {
            0.0
        } else {
            100.0 * hits as f64 / n_reads as f64
        };
        anchor_pcts.insert(TSO_PRESENT.seq, pct);
    }
    for spec in &candidates {
        for anchor in spec.r2_anchors {
            anchor_pcts
                .entry(anchor.seq)
                .or_insert_with(|| scan_anchor(r2_seqs, anchor));
        }
    }

    let cb_sizes: HashMap<(usize, usize), usize> =
        cb_sets.iter().map(|(k, v)| (*k, v.len())).collect();
    let mut scored: Vec<(&'static AssaySpec, f64)> = candidates
        .iter()
        .map(|spec| {
            (
                *spec,
                log_likelihood_assay(spec, &anchor_pcts, n_reads, &wl_scores, &cb_sizes),
            )
        })
        .collect();
    scored.sort_by(|a, b| {
        b.1.partial_cmp(&a.1)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                // Tiebreak: prefer assays whose expected R1 length exactly matches observed.
                let da = a.0.r1_expected_len().abs_diff(dominant);
                let db = b.0.r1_expected_len().abs_diff(dominant);
                da.cmp(&db)
            })
    });

    let has_wl = !wl_scores.is_empty();
    let confidence = classify_confidence(&scored, &candidates, dominant, n_reads, has_wl);

    // Single-cell plausibility gate: real single-cell data has a barcode-length R1
    // (seqspec geometry), or a TSO, or whitelist overlap. Bulk/WGS has none, so this
    // stops a 100 bp genomic R1 from being force-fit to a 10x assay.
    let tso_pct = anchor_pcts.get(TSO_PRESENT.seq).copied().unwrap_or(0.0);
    let max_wl_pct = wl_scores.values().map(|(_, p)| *p).fold(0.0_f64, f64::max);
    let is_single_cell = fits_barcode_len(dominant) || tso_pct >= 10.0 || max_wl_pct >= 50.0;

    Analysis {
        dominant,
        len_dist,
        candidates,
        cb_sets,
        anchor_pcts,
        wl_scores,
        scored,
        confidence,
        n_reads,
        tso_pct,
        max_wl_pct,
        is_single_cell,
    }
}

/// Total log-likelihood ratio for an assay (TSO + whitelist contributions summed).
fn log_likelihood_assay(
    spec: &AssaySpec,
    anchor_pcts: &HashMap<&'static [u8], f64>,
    n_reads: usize,
    wl_scores: &HashMap<&str, (usize, f64)>,
    cb_sizes: &HashMap<(usize, usize), usize>,
) -> f64 {
    let (tso_ll, wl_ll) = llr_components(spec, anchor_pcts, n_reads, wl_scores, cb_sizes);
    tso_ll + wl_ll
}

/// Per-signal log-likelihood ratio for an assay, returned as `(tso_llr, whitelist_llr)`.
/// Summed for ranking; kept separate so the report can attribute the verdict per signal.
///
/// TSO contribution is **one-sided** (`.max(0.0)`): a rate above background rewards the
/// assays that expect it, but a rate at/below background adds 0 rather than penalizing.
/// Real 3' libraries often show ~background TSO-at-R2-start, so penalizing its absence
/// wrongly defaulted ambiguous 26 bp libraries (3'v2 vs 5'v1, same onlist) to 5'.
fn llr_components(
    spec: &AssaySpec,
    anchor_pcts: &HashMap<&'static [u8], f64>,
    n_reads: usize,
    wl_scores: &HashMap<&str, (usize, f64)>,
    cb_sizes: &HashMap<(usize, usize), usize>,
) -> (f64, f64) {
    let mut tso_ll = 0.0_f64;
    for anchor in spec.r2_anchors {
        let pct = anchor_pcts.get(anchor.seq).copied().unwrap_or(0.0);
        let k = (pct / 100.0 * n_reads as f64).round() as usize;
        tso_ll += bernoulli_llr(k, n_reads, anchor.p_signal, P_TSO_BACKGROUND).max(0.0);
    }

    let mut wl_ll = 0.0_f64;
    if let Some(wl) = &spec.whitelist {
        if let Some(&(k_match, _)) = wl_scores.get(wl.filename) {
            if let Some(&n_cands) = cb_sizes.get(&(spec.r1_cb_offset, spec.r1_cb_len)) {
                if n_cands > 0 {
                    wl_ll += bernoulli_llr(k_match, n_cands, P_WL_CORRECT, P_WL_NULL);
                }
            }
        }
    }

    (tso_ll, wl_ll)
}

/// Two assays share the same R1 geometry and whitelist, so they are indistinguishable by sequence.
fn same_geometry(a: &AssaySpec, b: &AssaySpec) -> bool {
    a.r1_cb_offset == b.r1_cb_offset
        && a.r1_cb_len == b.r1_cb_len
        && a.r1_umi_len == b.r1_umi_len
        && a.whitelist.as_ref().map(|w| w.filename) == b.whitelist.as_ref().map(|w| w.filename)
}

/// True if `seq` carries `anchor`'s motif (its first `min_match_len` bases) within
/// the `[scan_start, scan_end)` window. Single per-read predicate shared by the SPRT
/// TSO counter and `scan_anchor`.
#[inline]
fn anchor_hit(seq: &[u8], anchor: &AnchorMotif) -> bool {
    let end = seq.len().min(anchor.scan_end);
    let needle = &anchor.seq[..anchor.min_match_len.min(anchor.seq.len())];
    end > anchor.scan_start
        && end >= anchor.min_match_len
        && contains(&seq[anchor.scan_start..end], needle)
}

fn fits_barcode_len(len: usize) -> bool {
    !crate::libstruct::sc_geometries_matching(len, 2).is_empty()
}

fn uppercase_key(bytes: &[u8]) -> Vec<u8> {
    if bytes.iter().any(|b| b.is_ascii_lowercase()) {
        bytes.iter().map(|b| b.to_ascii_uppercase()).collect()
    } else {
        bytes.to_vec()
    }
}

fn extract_cb_candidates(seqs: &[Vec<u8>], offset: usize, len: usize) -> HashSet<Vec<u8>> {
    let mut set = HashSet::with_capacity(seqs.len());
    for s in seqs {
        if s.len() >= offset + len {
            set.insert(uppercase_key(&s[offset..offset + len]));
        }
    }
    set
}

fn scan_anchor(reads: &[Vec<u8>], anchor: &AnchorMotif) -> f64 {
    if reads.is_empty() {
        return 0.0;
    }
    let hits = reads.iter().filter(|s| anchor_hit(s, anchor)).count();
    100.0 * hits as f64 / reads.len() as f64
}

/// Load a whitelist if it is present in the cache, else `None`. On native this
/// reads the file from disk; on wasm it reads the bytes the front-end stashed in
/// [`fastq::vfs`] under the same `cache/filename` key.
#[cfg(not(target_arch = "wasm32"))]
fn load_whitelist_opt(path: &Path, cb_len: usize) -> Result<Option<HashSet<Vec<u8>>>> {
    if !path.exists() {
        return Ok(None);
    }
    let f = std::fs::File::open(path)?;
    let gz = path.extension().and_then(|e| e.to_str()) == Some("gz");
    let reader: Box<dyn std::io::Read> = if gz {
        Box::new(flate2::read::MultiGzDecoder::new(f))
    } else {
        Box::new(f)
    };
    Ok(Some(parse_whitelist(
        std::io::BufReader::new(reader),
        cb_len,
    )?))
}

#[cfg(target_arch = "wasm32")]
fn load_whitelist_opt(path: &Path, cb_len: usize) -> Result<Option<HashSet<Vec<u8>>>> {
    use std::io::Cursor;
    let key = path.to_string_lossy();
    let Some(bytes) = fastq::vfs::get(&key) else {
        return Ok(None);
    };
    let reader: Box<dyn std::io::Read> = if fastq::is_gzip(&bytes) {
        Box::new(flate2::read::MultiGzDecoder::new(Cursor::new(bytes)))
    } else {
        Box::new(Cursor::new(bytes))
    };
    Ok(Some(parse_whitelist(
        std::io::BufReader::new(reader),
        cb_len,
    )?))
}

fn trim_ascii(mut line: &[u8]) -> &[u8] {
    while line.first().is_some_and(|b| b.is_ascii_whitespace()) {
        line = &line[1..];
    }
    while line.last().is_some_and(|b| b.is_ascii_whitespace()) {
        line = &line[..line.len() - 1];
    }
    line
}

/// Read one barcode per line, keying on the uppercased first `cb_len` bytes.
fn parse_whitelist<R: BufRead>(mut reader: R, cb_len: usize) -> Result<HashSet<Vec<u8>>> {
    // 10x onlists run to several million barcodes; preallocate to avoid rehashing.
    let mut set = HashSet::with_capacity(4_000_000);
    let mut line = Vec::with_capacity(128);
    loop {
        line.clear();
        if reader.read_until(10, &mut line)? == 0 {
            break;
        }
        let line = trim_ascii(&line);
        if line.len() >= cb_len {
            set.insert(uppercase_key(&line[..cb_len]));
        }
    }
    Ok(set)
}

fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    if needle.is_empty() || haystack.len() < needle.len() {
        return false;
    }
    haystack.windows(needle.len()).any(|w| w == needle)
}

mod report;
pub(crate) use report::build_report;

#[cfg(test)]
mod tests;
