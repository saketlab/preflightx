use anyhow::Result;
use serde::Serialize;
use std::collections::{HashMap, HashSet};
use std::io::BufRead;
use std::path::Path;

use crate::assay::{AnchorMotif, AssaySpec, ASSAYS, TSO_PRESENT};
use crate::fastq;
use crate::stats::{bernoulli_llr, dominant_len};

// Rate at which the TSO 13-mer appears in R2 by chance.
const P_TSO_BACKGROUND: f64 = 0.01;
// CB match rates expected when the whitelist is correct vs wrong.
const P_WL_CORRECT: f64 = 0.70;
const P_WL_NULL: f64 = 0.01;
/// Onlist overlap (%) below which the whitelist contradicts the winner; real cell
/// barcodes hit their own onlist almost completely.
pub(crate) const WL_CONTRADICT_PCT: f64 = 10.0;

// Per-read (LLR / n_reads) so the thresholds are sample-size-independent.
const HIGH_LL_PER_READ: f64 = 0.30;
const HIGH_GAP_PER_READ: f64 = 0.10; // lead over the runner-up
const MEDIUM_LL_PER_READ: f64 = 0.001;

// SPRT early-stopping (no-whitelist path only; whitelist needs a large CB sample).
const SPRT_THRESHOLD: f64 = 200.0; // nats
const SPRT_MIN_READS: usize = 5_000;
const SPRT_BATCH: usize = 2_000;

/// Above this never-repeated-barcode fraction, a length-only match (see the
/// is_single_cell gate in analyze) reads as noise, not real cell multiplexing.
const MAX_NOVEL_BARCODE_RATIO: f64 = 0.5;

/// After CB+UMI in R1, 10x 5' has the TSO and 10x 3' has oligo-dT.
const R1_TSO_5P: &[u8] = b"TTTCTTATATGGG";
const R1_TSO_WINDOW: (usize, usize) = (24, 45);
const R1_POLYT: &[u8] = b"TTTTTTTTTTTT";
const R1_POLYT_WINDOW: (usize, usize) = (24, 60);
const R1_TSO_MIN_PCT: f64 = 30.0;
const R1_POLYT_MIN_PCT: f64 = 50.0;
const P_R1_MOTIF_SIGNAL: f64 = 0.50;
const P_R1_MOTIF_BACKGROUND: f64 = 0.01;

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

/// One reasoning step in the verdict's evidence ledger.
#[derive(Serialize, Debug)]
pub struct Evidence {
    pub signal: String,
    pub observation: String,
    /// Nats of log-likelihood this signal added to the winner; None for non-LLR steps.
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
    /// False when the reads show no single-cell signal at all; assay_id is then "none".
    pub is_single_cell: bool,
    /// SPRT stopped sampling early, so reads_sampled < n is expected, not a truncated input.
    pub early_stop: bool,
    pub whitelist: Option<String>,
    pub whitelist_scores: Vec<WhitelistScore>,
    pub anchor_scores: Vec<AnchorScore>,
    pub tool_hints: Vec<ToolHint>,
    pub confidence: String,
    /// "atac" | "rna" | "unknown"; "atac" overrides the RNA chemistry call.
    pub modality: String,
    /// promoter / (promoter + exon) k-mer hits; high -> ATAC. 0.0 when no panel/hits.
    pub atac_promoter_frac: f64,
    /// "3'" | "5'" | "full-length" | "ambiguous" | "unknown", from the k-mer panel;
    /// R1 geometry alone cannot tell 3'v2 from 5'.
    pub end: String,
    pub end_confidence: String,
    /// Reads hitting a 5'-end / 3'-end transcript marker, and r = n3/(n3+n5).
    pub n5: usize,
    pub n3: usize,
    pub end_ratio: f64,
    /// geometry scored on the R1 prefix
    pub long_r1: bool,
    /// 0 unless long_r1
    pub r1_tso_pct: f64,
    pub r1_polyt_pct: f64,
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
    // orient + sample + match whitelists (I/O) -> analyze (pure) -> build_report (pure).
    let (eff_r1, eff_r2, swapped) = orient(r1, r2)?;
    let wl_candidates: Vec<&'static AssaySpec> = match whitelist_cache {
        Some(cache) => {
            let (preview, _) = fastq::sample(eff_r1, n.min(2_000))?;
            candidates_for(dominant_len(&preview))
                .into_iter()
                .filter(|s| {
                    s.whitelist
                        .is_some_and(|wl| whitelist_present(&cache.join(wl.filename)))
                })
                .collect()
        }
        None => Vec::new(),
    };
    // Whitelist overlap and the k-mer panels all need depth, so skip the SPRT early-stop.
    let want_full = !wl_candidates.is_empty()
        || crate::endpanel::have_panels()
        || crate::atacpanel::have_panels();
    let s = sample_reads(eff_r1, eff_r2, n, want_full)?;
    let wl_sets = whitelist_cache
        .map(|c| load_candidate_whitelists(c, &wl_candidates, &s.r1))
        .transpose()?
        .unwrap_or_default();

    let analysis = analyze(&s.r1, &s.r2, &wl_sets, s.sprt_tso_hits);
    let end = crate::endpanel::classify(&s.r2, species);
    let atac = crate::atacpanel::classify_pair(&s.r1, &s.r2, species);
    Ok(build_report(analysis, swapped, n, s.early_stop, end, atac))
}

/// Distinct whitelist filenames matching the pair's R1 geometry; the wasm front-end
/// fetches only these instead of every catalogued whitelist.
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

/// Effective (r1, r2) ordering, plus whether a swap was applied: a long R1 with a
/// barcode-length R2 means the mates were handed over the wrong way round.
fn orient<'a>(r1: &'a Path, r2: &'a Path) -> Result<(&'a Path, &'a Path, bool)> {
    let (a_quick, b_quick) = fastq::par_pair(
        || fastq::sample_seqs(r1, 1_000),
        || fastq::sample_seqs(r2, 1_000),
    )?;
    let swapped =
        !fits_barcode_len(dominant_len(&a_quick)) && fits_barcode_len(dominant_len(&b_quick));
    Ok(if swapped {
        (r2, r1, true)
    } else {
        (r1, r2, false)
    })
}

/// Sampled reads plus how sampling ended; the TSO count is set on the SPRT path only.
struct Sampled {
    r1: Vec<Vec<u8>>,
    r2: Vec<Vec<u8>>,
    sprt_tso_hits: Option<usize>,
    early_stop: bool,
}

/// Sample up to n read pairs; without a whitelist cache, SPRT stops as soon as the
/// TSO signal is decisive.
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

/// Sampled R1 barcodes found on each cached candidate whitelist.
fn load_candidate_whitelists(
    cache: &Path,
    candidates: &[&'static AssaySpec],
    r1_seqs: &[Vec<u8>],
) -> Result<HashMap<&'static str, HashSet<Vec<u8>>>> {
    let mut wanted: HashMap<&'static str, (usize, HashSet<Vec<u8>>)> = HashMap::new();
    let mut done = HashSet::new();
    for spec in candidates {
        if let Some(wl) = &spec.whitelist {
            if !done.insert((wl.filename, spec.r1_cb_offset)) {
                continue;
            }
            wanted
                .entry(wl.filename)
                .or_insert_with(|| (wl.cb_len, HashSet::new()))
                .1
                .extend(extract_cb_candidates(r1_seqs, spec.r1_cb_offset, wl.cb_len));
        }
    }
    let mut sets: HashMap<&'static str, HashSet<Vec<u8>>> = HashMap::new();
    for (filename, (cb_len, keep)) in wanted {
        if let Some(reader) = open_whitelist(&cache.join(filename))? {
            sets.insert(filename, parse_whitelist(reader, cb_len, &keep)?);
        }
    }
    Ok(sets)
}

// pure scoring core (the testable seam)

/// The ranked outcome of scoring sampled reads against the assay registry. Holds
/// everything build_report needs; produced by analyze with no I/O.
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
    long_r1: bool,
    /// 0 unless long_r1
    r1_tso_pct: f64,
    r1_polyt_pct: f64,
}

/// Candidate assays for an observed R1 length: exact match, or +/-2 bp for
/// non-standard cycle counts.
///
fn candidates_for(dominant: usize) -> Vec<&'static AssaySpec> {
    if is_long_r1(dominant) {
        // inside a long read only an onlist can test the prefix
        return ASSAYS
            .iter()
            .filter(|s| s.whitelist.is_some() && s.r1_expected_len() <= dominant)
            .collect();
    }
    ASSAYS
        .iter()
        .filter(|s| s.r1_expected_len().abs_diff(dominant) <= 2)
        .collect()
}

pub(crate) fn is_long_r1(dominant: usize) -> bool {
    let longest = ASSAYS
        .iter()
        .map(|s| s.r1_expected_len())
        .max()
        .unwrap_or(0);
    dominant > longest + 2 && !fits_barcode_len(dominant)
}

fn motif_pct(reads: &[Vec<u8>], needle: &[u8], window: (usize, usize)) -> f64 {
    if reads.is_empty() {
        return 0.0;
    }
    let hits = reads
        .iter()
        .filter(|s| {
            let end = s.len().min(window.1);
            end > window.0 && contains(&s[window.0..end], needle)
        })
        .count();
    100.0 * hits as f64 / reads.len() as f64
}

fn r1_motif_llr(spec: &AssaySpec, r1_tso_pct: f64, r1_polyt_pct: f64, n_reads: usize) -> f64 {
    let pct = if spec.id.contains("5prime") {
        r1_tso_pct
    } else if spec.id.contains("3prime") {
        r1_polyt_pct
    } else {
        return 0.0;
    };
    let k = (pct / 100.0 * n_reads as f64).round() as usize;
    bernoulli_llr(k, n_reads, P_R1_MOTIF_SIGNAL, P_R1_MOTIF_BACKGROUND).max(0.0)
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

/// Pick the confidence bucket from the per-read evidence: high needs strong total
/// evidence, a clear lead over the runner-up, and a whitelist; medium needs any
/// positive discrimination (or a single unambiguous candidate); else low.
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
/// access; whitelist barcode sets are passed in via wl_sets.
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

    // Whitelist overlap, computed once per unique whitelist file present in wl_sets.
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

    let long_r1 = is_long_r1(dominant);
    let (r1_tso_pct, r1_polyt_pct) = if long_r1 {
        (
            motif_pct(r1_seqs, R1_TSO_5P, R1_TSO_WINDOW),
            motif_pct(r1_seqs, R1_POLYT, R1_POLYT_WINDOW),
        )
    } else {
        (0.0, 0.0)
    };

    let cb_sizes: HashMap<(usize, usize), usize> =
        cb_sets.iter().map(|(k, v)| (*k, v.len())).collect();
    let mut scored: Vec<(&'static AssaySpec, f64)> = candidates
        .iter()
        .map(|spec| {
            let mut ll = log_likelihood_assay(spec, &anchor_pcts, n_reads, &wl_scores, &cb_sizes);
            if long_r1 {
                ll += r1_motif_llr(spec, r1_tso_pct, r1_polyt_pct, n_reads);
            }
            (*spec, ll)
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
    //
    // Length alone is weak: a coincidental match (e.g. Smart-seq sharing a
    // no-whitelist assay's R1 length) rarely reuses the same barcode value, while a
    // real barcode pool does. A TSO or whitelist hit is direct evidence and skips
    // this reuse check.
    let tso_pct = anchor_pcts.get(TSO_PRESENT.seq).copied().unwrap_or(0.0);
    let max_wl_pct = wl_scores.values().map(|(_, p)| *p).fold(0.0_f64, f64::max);
    let denom = n_reads.max(1) as f64;
    let reused = |key: &(usize, usize)| {
        *cb_sizes.get(key).unwrap_or(&0) as f64 / denom <= MAX_NOVEL_BARCODE_RATIO
    };
    let barcode_reused = candidates
        .iter()
        .filter(|s| s.whitelist.is_none())
        .map(|s| (s.r1_cb_offset, s.r1_cb_len))
        .collect::<HashSet<_>>()
        .iter()
        .all(reused);
    // oligo-dT also occurs in 3'-tag bulk, so require repeating barcodes
    let prefix_reused = long_r1 && cb_sizes.contains_key(&(0, 16)) && reused(&(0, 16));
    let long_r1_sc = long_r1
        && (r1_tso_pct >= R1_TSO_MIN_PCT || (r1_polyt_pct >= R1_POLYT_MIN_PCT && prefix_reused));
    let is_single_cell = (fits_barcode_len(dominant) && barcode_reused)
        || tso_pct >= 10.0
        || max_wl_pct >= 50.0
        || long_r1_sc;

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
        long_r1,
        r1_tso_pct,
        r1_polyt_pct,
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

/// Per-signal log-likelihood ratio for an assay, returned as (tso_llr, whitelist_llr).
/// Summed for ranking; kept separate so the report can attribute the verdict per signal.
///
/// TSO contribution is **one-sided** (.max(0.0)): a rate above background rewards the
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

/// True if seq carries anchor's motif (its first min_match_len bases) within
/// the [scan_start, scan_end) window. Single per-read predicate shared by the SPRT
/// TSO counter and scan_anchor.
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

#[cfg(not(target_arch = "wasm32"))]
fn whitelist_present(path: &Path) -> bool {
    path.exists()
}

#[cfg(target_arch = "wasm32")]
fn whitelist_present(path: &Path) -> bool {
    fastq::vfs::get(&path.to_string_lossy()).is_some()
}

#[cfg(not(target_arch = "wasm32"))]
fn open_whitelist(path: &Path) -> Result<Option<Box<dyn BufRead>>> {
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
    Ok(Some(Box::new(std::io::BufReader::with_capacity(
        1 << 20,
        reader,
    ))))
}

#[cfg(target_arch = "wasm32")]
fn open_whitelist(path: &Path) -> Result<Option<Box<dyn BufRead>>> {
    use std::io::Cursor;
    let Some(bytes) = fastq::vfs::get(&path.to_string_lossy()) else {
        return Ok(None);
    };
    let reader: Box<dyn std::io::Read> = if fastq::is_gzip(&bytes) {
        Box::new(flate2::read::MultiGzDecoder::new(Cursor::new(bytes)))
    } else {
        Box::new(Cursor::new(bytes))
    };
    Ok(Some(Box::new(std::io::BufReader::new(reader))))
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

/// Whitelist barcodes that are also in keep.
fn parse_whitelist<R: BufRead>(
    mut reader: R,
    cb_len: usize,
    keep: &HashSet<Vec<u8>>,
) -> Result<HashSet<Vec<u8>>> {
    let mut set = HashSet::with_capacity(keep.len());
    let mut line = Vec::with_capacity(128);
    loop {
        line.clear();
        if reader.read_until(10, &mut line)? == 0 {
            break;
        }
        let line = trim_ascii(&line);
        if line.len() < cb_len {
            continue;
        }
        let bc = &line[..cb_len];
        if keep.contains(bc) {
            set.insert(bc.to_vec());
        } else if bc.iter().any(u8::is_ascii_lowercase) {
            let key = bc.to_ascii_uppercase();
            if keep.contains(&key) {
                set.insert(key);
            }
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
pub(crate) fn test_chem_report(is_single_cell: bool, confidence: &str) -> ChemReport {
    ChemReport {
        r1_length_dominant: 28,
        r1_length_distribution: vec![],
        reads_sampled: 0,
        assay_id: "none".into(),
        // matches what chemistry::report writes on a negative call
        assay_display: if is_single_cell {
            "10x 3' v3".into()
        } else {
            "Not single-cell".into()
        },
        cb_read: "R1".into(),
        cb_offset: 0,
        cb_len: 16,
        umi_offset: 16,
        umi_len: 12,
        n_cb_candidates: 0,
        is_single_cell,
        early_stop: false,
        whitelist: None,
        whitelist_scores: vec![],
        anchor_scores: vec![],
        tool_hints: vec![],
        confidence: confidence.into(),
        modality: "rna".into(),
        atac_promoter_frac: 0.0,
        end: "unknown".into(),
        end_confidence: "low".into(),
        n5: 0,
        n3: 0,
        end_ratio: 0.0,
        long_r1: false,
        r1_tso_pct: 0.0,
        r1_polyt_pct: 0.0,
        interpretation: vec![],
        notes: vec![],
    }
}

#[cfg(test)]
mod tests;
