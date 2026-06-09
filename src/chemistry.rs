use anyhow::Result;
use colored::{ColoredString, Colorize};
use serde::Serialize;
use std::collections::{HashMap, HashSet};
use std::path::Path;

use crate::assay::{AnchorMotif, AssaySpec, ASSAYS, TSO_PRESENT};
use crate::fastq;
use crate::stats::{bernoulli_llr, dominant_len};

// log-likelihood model parameters
// TSO background rate: fraction of R2 reads where the 13-mer appears by chance.
const P_TSO_BACKGROUND: f64 = 0.01;
// Whitelist: expected CB match rates when the whitelist is correct vs wrong.
const P_WL_CORRECT: f64 = 0.70;
const P_WL_NULL: f64 = 0.01;

// Dividing LLR by n_reads makes the threshold sample-size-independent.
const HIGH_LL_PER_READ: f64 = 0.30; // strong total evidence
const HIGH_GAP_PER_READ: f64 = 0.10; // top assay uniquely ahead of second
const MEDIUM_LL_PER_READ: f64 = 0.001; // any positive discrimination

// SPRT early-stopping parameters
// Used only when no whitelist cache is present (whitelist needs large CB sample).
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

/// One reasoning step in the verdict's evidence ledger. Together these form the
/// `interpretation` chain: each signal, what was observed, how many nats of
/// log-likelihood it contributed to the winning assay, and what it implies.
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
    pub interpretation: Vec<Evidence>,
    pub notes: Vec<String>,
}

pub fn detect(
    r1: &Path,
    r2: &Path,
    n: usize,
    whitelist_cache: Option<&Path>,
) -> Result<ChemReport> {
    // 1. Auto-orient (I/O), 2. sample reads (I/O), 3. load whitelists (I/O),
    // 4. score (pure), 5. assemble the report (pure). Steps 4-5 are unit-tested
    // directly via `analyze`; the I/O steps are thin and orientation-preserving.
    let (eff_r1, eff_r2, swapped) = orient(r1, r2)?;
    let s = sample_reads(eff_r1, eff_r2, n, whitelist_cache.is_some())?;

    let dominant = dominant_len(&s.r1);
    let candidates = candidates_for(dominant);
    let wl_sets = match whitelist_cache {
        Some(cache) => load_candidate_whitelists(cache, &candidates)?,
        None => HashMap::new(),
    };

    let analysis = analyze(&s.r1, &s.r2, &wl_sets, s.sprt_tso_hits);
    Ok(build_report(analysis, swapped, n, s.early_stop))
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
    let (a_quick, _) = fastq::sample(r1, 1_000)?;
    let (b_quick, _) = fastq::sample(r2, 1_000)?;
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
struct Analysis {
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
/// access — whitelist barcode sets are passed in via `wl_sets`. This is the seam the
/// chemistry tests drive directly.
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

/// Turn a scored `Analysis` into the serializable report: winner geometry, tool
/// hints, whitelist/anchor tables, the evidence ledger, and human-facing notes.
/// `n_requested` is the read budget asked for; `early_stop` is true when SPRT
/// decided before reaching it (so a short sample is explained, not mistaken for a
/// truncated input).
fn build_report(a: Analysis, swapped: bool, n_requested: usize, early_stop: bool) -> ChemReport {
    let Analysis {
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
    } = a;

    // A non-single-cell library has no real chemistry to report, and its scored
    // winner is a force-fit; never present it as a confident call.
    let confidence = if is_single_cell { confidence } else { "low" };
    let has_wl = !wl_scores.is_empty();
    let cb_sizes: HashMap<(usize, usize), usize> =
        cb_sets.iter().map(|(k, v)| (*k, v.len())).collect();
    let winner = scored.first().map(|(s, _)| *s);
    let runner = scored.get(1).map(|(s, _)| *s);

    // 3'v2 and 5'v1 share both geometry and onlist; once the TSO no longer penalizes
    // absence (one-sided LLR), the top two can tie and differ only in 3'/5' polarity.
    // That polarity isn't recoverable from the reads, so report both rather than guess.
    let gap = scored.first().map(|(_, s)| *s).unwrap_or(0.0)
        - scored.get(1).map(|(_, s)| *s).unwrap_or(0.0);
    let polarity_unresolved = match (winner, runner) {
        (Some(w), Some(r)) => {
            is_single_cell
                && gap.abs() < 1.0
                && same_geometry(w, r)
                && w.id.contains("3prime") != r.id.contains("3prime")
        }
        _ => false,
    };
    // Can't claim "high" for a specific polarity when the polarity is unresolved.
    let confidence = if polarity_unresolved && confidence == "high" {
        "medium"
    } else {
        confidence
    };

    let mut notes: Vec<String> = Vec::new();
    if !is_single_cell {
        notes.push(format!(
            "No single-cell signal detected: R1 is {dominant}bp (single-cell barcode reads are typically ~24-30bp), \
             TSO appears in {tso_pct:.0}% of R2 (background ~1%), and whitelist overlap is {max_wl_pct:.0}%. \
             This looks like bulk or genomic sequencing, not single-cell."
        ));
        if let Some(w) = winner {
            notes.push(format!(
                "If it were single-cell, the closest geometry would be {} ({}bp: CB{}+UMI{}) — but the evidence above does not support a single-cell call.",
                w.display_name, w.r1_expected_len(), w.r1_cb_len, w.r1_umi_len,
            ));
        }
    }
    if polarity_unresolved {
        if let (Some(w), Some(r)) = (winner, runner) {
            notes.insert(
                0,
                format!(
                    "3' vs 5' is unresolved: {} and {} share the same {dominant}bp barcode geometry and the same onlist, and the R2 TSO is at background, so the reads alone can't tell them apart. Check the library-prep records to confirm which.",
                    w.display_name, r.display_name,
                ),
            );
        }
    }
    if swapped {
        notes.push("Input files appeared swapped (R1 was long, R2 was short); auto-oriented to treat the shorter file as R1 (barcodes).".into());
    }
    if early_stop {
        notes.push(format!(
            "Stopped early after {n_reads} reads: the TSO signal was already decisive (sequential test), so more reads were not needed."
        ));
    } else if n_reads < n_requested {
        notes.push(format!(
            "Input ended after {n_reads} reads (fewer than the {n_requested} requested) — a short file or a sliced fetch."
        ));
    }
    // The assay-specific notes below only make sense for a real single-cell call;
    // skip them entirely when the reads are not single-cell.
    if is_single_cell {
        // Chemistries like Drop-seq use random (non-whitelisted) barcodes, so there is no
        // list to confirm them against and `fetch-whitelists` cannot help — say so instead
        // of suggesting a download that won't change the verdict.
        if winner.map(|w| w.whitelist.is_none()).unwrap_or(false) {
            notes.push(format!(
                "{} uses random (non-whitelisted) cell barcodes; no whitelist exists to confirm them, so confidence is capped at MEDIUM (geometry-only).",
                winner.map(|w| w.display_name).unwrap_or("This chemistry"),
            ));
        } else if !has_wl {
            notes.push(
                "No whitelist cache; run `preflightx fetch-whitelists` for full disambiguation."
                    .into(),
            );
        }
        if confidence != "high" && scored.len() > 1 {
            let alts: Vec<&str> = scored.iter().skip(1).take(3).map(|(s, _)| s.id).collect();
            notes.push(format!("Other candidates: {}.", alts.join(", ")));
        }
        if let Some(w) = winner {
            let expected = w.r1_expected_len();
            if expected != dominant && dominant > 0 {
                notes.push(format!(
                    "R1 reads are {}bp but {} expects {}bp; library may have been sequenced at non-standard cycles.",
                    dominant, w.id, expected,
                ));
            }
        }
    }

    let mut assay_id = winner.map(|w| w.id).unwrap_or("unknown").to_string();
    let mut assay_display = winner
        .map(|w| w.display_name)
        .unwrap_or("unknown")
        .to_string();
    if !is_single_cell {
        assay_id = "none".to_string();
        assay_display = "Not single-cell".to_string();
    } else if polarity_unresolved {
        if let (Some(w), Some(r)) = (winner, runner) {
            assay_display = format!("{} or {}", w.display_name, r.display_name);
        }
    }
    let cb_offset = winner.map(|w| w.r1_cb_offset).unwrap_or(0);
    let cb_len = winner.map(|w| w.r1_cb_len).unwrap_or(0);
    let umi_offset = winner.map(|w| w.umi_offset()).unwrap_or(0);
    let umi_len = winner.map(|w| w.r1_umi_len).unwrap_or(0);
    let whitelist_file =
        winner.and_then(|w| w.whitelist.as_ref().map(|wl| wl.filename.to_string()));
    let tool_hints = winner
        .map(|w| {
            let mut hints = vec![ToolHint {
                tool: "starsolo".into(),
                hint: w.starsolo_hint(),
            }];
            hints.extend(w.tool_hints.iter().map(|(t, h)| ToolHint {
                tool: t.to_string(),
                hint: h.to_string(),
            }));
            hints
        })
        .unwrap_or_default();

    let n_cb = winner
        .and_then(|w| cb_sets.get(&(w.r1_cb_offset, w.r1_cb_len)))
        .map(|s| s.len())
        .unwrap_or(0);

    let mut wl_scores_out: Vec<WhitelistScore> = wl_scores
        .iter()
        .map(|(f, (hits, pct))| WhitelistScore {
            filename: f.to_string(),
            hits: *hits,
            pct: fastq::round2(*pct),
        })
        .collect();
    wl_scores_out.sort_by(|a, b| {
        b.pct
            .partial_cmp(&a.pct)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let anchor_scores_out: Vec<AnchorScore> = anchor_pcts
        .iter()
        .map(|(seq, pct)| AnchorScore {
            motif: String::from_utf8_lossy(&seq[..seq.len().min(13)]).into_owned(),
            hit_pct: fastq::round2(*pct),
        })
        .collect();

    let mut interpretation = if is_single_cell {
        build_interpretation(
            winner,
            &scored,
            dominant,
            &candidates,
            &anchor_pcts,
            &wl_scores,
            &cb_sizes,
            n_reads,
            has_wl,
            confidence,
        )
    } else {
        build_not_single_cell_interpretation(dominant, tso_pct, max_wl_pct, n_reads, winner)
    };
    if polarity_unresolved {
        if let (Some(w), Some(r)) = (winner, runner) {
            interpretation.push(Evidence {
                signal: "3' vs 5' polarity".into(),
                observation: format!(
                    "{} and {} tie: same {dominant}bp geometry and onlist, TSO at background",
                    w.id, r.id,
                ),
                llr_contribution: None,
                conclusion:
                    "The reads carry no signal that separates 3' from 5' here, so both are reported. The R2 TSO would normally break this tie when present.".into(),
            });
        }
    }

    ChemReport {
        r1_length_dominant: dominant,
        r1_length_distribution: len_dist.into_iter().take(5).collect(),
        reads_sampled: n_reads,
        assay_id,
        assay_display,
        cb_read: "R1".into(),
        cb_offset,
        cb_len,
        umi_offset,
        umi_len,
        n_cb_candidates: n_cb,
        is_single_cell,
        early_stop,
        whitelist: whitelist_file,
        whitelist_scores: wl_scores_out,
        anchor_scores: anchor_scores_out,
        tool_hints,
        confidence: confidence.to_string(),
        interpretation,
        notes,
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
/// The two are summed for ranking; kept separate so the interpretation report can
/// attribute the verdict to each signal.
///
/// TSO contribution is **one-sided**: an observed TSO rate above background rewards
/// the assays that expect it (positive LLR), but a rate at/below background adds 0
/// rather than penalizing. Real 3' libraries often show ~background TSO-at-R2-start,
/// so penalizing its absence wrongly rejected 3' and defaulted ambiguous 26 bp
/// libraries (3'v2 vs 5'v1, same onlist) to 5'. `.max(0.0)` keeps presence as evidence
/// while making absence non-informative.
///
/// Whitelist contribution: LLR of observed match rate under P_WL_CORRECT vs P_WL_NULL.
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

/// Build the evidence ledger explaining *why* `winner` was chosen, signal by signal.
/// Mirrors the decision flow in `detect()`: length gives candidates, TSO gives
/// polarity, whitelist gives version, then the final scored verdict.
#[allow(clippy::too_many_arguments)]
fn build_interpretation(
    winner: Option<&AssaySpec>,
    scored: &[(&AssaySpec, f64)],
    dominant: usize,
    candidates: &[&AssaySpec],
    anchor_pcts: &HashMap<&'static [u8], f64>,
    wl_scores: &HashMap<&str, (usize, f64)>,
    cb_sizes: &HashMap<(usize, usize), usize>,
    n_reads: usize,
    has_wl: bool,
    confidence: &str,
) -> Vec<Evidence> {
    let mut ev = Vec::new();
    let winner = match winner {
        Some(w) => w,
        None => return ev,
    };

    let cand_ids: Vec<&str> = candidates.iter().map(|s| s.id).collect();
    ev.push(Evidence {
        signal: "R1 length".into(),
        observation: format!("{dominant} bp dominant across {n_reads} reads"),
        llr_contribution: None,
        conclusion: format!(
            "The barcode read length fixes the cell-barcode + UMI layout, so only {} \
             chemistries with a ~{dominant} bp R1 (\u{00b1}2 bp) are possible. The next \
             signals pick among them: {}",
            candidates.len(),
            cand_ids.join(", "),
        ),
    });

    let (tso_ll, wl_ll) = llr_components(winner, anchor_pcts, n_reads, wl_scores, cb_sizes);

    // TSO polarity, only for assays that model the motif.
    if !winner.r2_anchors.is_empty() {
        let pct = winner
            .r2_anchors
            .iter()
            .filter_map(|a| anchor_pcts.get(a.seq).copied())
            .fold(0.0_f64, f64::max);
        // Wording follows the observed level, not the winner's model: only call it a
        // TSO signal when it actually contributed evidence (one-sided, so absence = 0).
        let conclusion = if tso_ll > 1.0 {
            format!(
                "TSO present in {pct:.1}% of R2 (background ~{:.0}%): the template-switch oligo is detected, which marks 3' chemistry.",
                P_TSO_BACKGROUND * 100.0,
            )
        } else {
            format!(
                "TSO in only {pct:.1}% of R2, at the ~{:.0}% background: no template-switch signal, so the TSO can't tell 3' from 5' here.",
                P_TSO_BACKGROUND * 100.0,
            )
        };
        let needle = String::from_utf8_lossy(&TSO_PRESENT.seq[..TSO_PRESENT.min_match_len]);
        ev.push(Evidence {
            signal: "TSO motif (R2)".into(),
            observation: format!(
                "{}-mer {needle} in {pct:.1}% of R2[0..{}]",
                TSO_PRESENT.min_match_len, TSO_PRESENT.scan_end,
            ),
            llr_contribution: Some(fastq::round2(tso_ll)),
            conclusion,
        });
    }

    if let Some(wl) = &winner.whitelist {
        if has_wl {
            let win_pct = wl_scores.get(wl.filename).map(|(_, p)| *p).unwrap_or(0.0);
            let alt = wl_scores
                .iter()
                .filter(|(f, _)| **f != wl.filename)
                .max_by(|a, b| {
                    a.1 .1
                        .partial_cmp(&b.1 .1)
                        .unwrap_or(std::cmp::Ordering::Equal)
                });
            let alt_str = alt
                .map(|(f, (_, p))| format!("next pool {} at {:.1}%", f, p))
                .unwrap_or_else(|| "no other pool scored".into());
            let decisive = match alt {
                Some((_, (_, p))) => win_pct >= 50.0 && win_pct >= 3.0 * p.max(0.1),
                None => win_pct >= 50.0,
            };
            let conclusion = if decisive {
                format!(
                    "{:.1}% of barcodes in the {} pool ({}); decisive.",
                    win_pct, wl.filename, alt_str,
                )
            } else {
                format!(
                    "{} matched {:.1}% ({}); not decisive, version from geometry and registry order.",
                    wl.filename, win_pct, alt_str,
                )
            };
            ev.push(Evidence {
                signal: "Whitelist overlap".into(),
                observation: format!("{:.1}% of distinct CBs found in {}", win_pct, wl.filename),
                llr_contribution: Some(fastq::round2(wl_ll)),
                conclusion,
            });
        } else {
            ev.push(Evidence {
                signal: "Whitelist overlap".into(),
                observation: "no whitelist cache".into(),
                llr_contribution: None,
                conclusion: "version within family not checked against barcodes; resolved by length and registry order. Run `preflightx fetch-whitelists` to disambiguate.".into(),
            });
        }
    }

    let best = scored.first().map(|(_, s)| *s).unwrap_or(0.0);
    let second = scored.get(1).map(|(_, s)| *s).unwrap_or(0.0);
    let gap = best - second;
    let runner = scored.get(1).map(|(s, _)| *s);
    let tie = runner.is_some_and(|r| gap.abs() < 1.0 && same_geometry(winner, r));
    let conclusion = if tie {
        format!(
            "{} and {} are indistinguishable by sequence (same geometry and whitelist); \
             chose {} by registry order. Confidence {}.",
            winner.id,
            runner.map(|r| r.id).unwrap_or("?"),
            winner.id,
            confidence,
        )
    } else {
        format!(
            "{}: {:+.0} nats total ({:.3}/read), {:+.0} nats ahead of {}. Confidence {} \
             (high requires a whitelist and a clear lead).",
            winner.id,
            best,
            best / n_reads.max(1) as f64,
            gap,
            runner.map(|r| r.id).unwrap_or("(none)"),
            confidence,
        )
    };
    ev.push(Evidence {
        signal: "Verdict".into(),
        observation: format!("winner {} ({})", winner.id, winner.display_name),
        llr_contribution: Some(fastq::round2(best)),
        conclusion,
    });

    ev
}

/// Evidence ledger for a library with no single-cell signal. Walks the same three
/// signals the positive path uses, but each step says why it *fails* the test, so
/// the reader sees exactly how we concluded "not single-cell" rather than a verdict.
fn build_not_single_cell_interpretation(
    dominant: usize,
    tso_pct: f64,
    max_wl_pct: f64,
    n_reads: usize,
    winner: Option<&AssaySpec>,
) -> Vec<Evidence> {
    let closest = winner
        .map(|w| {
            format!(
                " The closest catalogued geometry is {} ({} bp).",
                w.id,
                w.r1_expected_len()
            )
        })
        .unwrap_or_default();
    vec![
        Evidence {
            signal: "Barcode geometry".into(),
            observation: format!("R1 is {dominant} bp dominant across {n_reads} reads"),
            llr_contribution: None,
            conclusion: format!(
                "Single-cell barcode reads are short (~24-30 bp: a cell barcode + UMI). \
                 {dominant} bp does not match any catalogued single-cell geometry, so no \
                 barcode read is present.{closest}"
            ),
        },
        Evidence {
            signal: "TSO (template-switch oligo)".into(),
            observation: format!("TSO 13-mer in {tso_pct:.1}% of R2 (chance background ~1%)"),
            llr_contribution: None,
            conclusion: "A 10x 3' library carries the TSO in ~20% of R2 reads. At background \
                         level here, there is no 3' single-cell signature."
                .into(),
        },
        Evidence {
            signal: "Barcode whitelist".into(),
            observation: format!("{max_wl_pct:.1}% of R1 barcodes match any 10x onlist"),
            llr_contribution: None,
            conclusion: "Real cell barcodes hit their onlist almost completely. Near-zero \
                         overlap means R1 is not a barcode read."
                .into(),
        },
        Evidence {
            signal: "Verdict".into(),
            observation: "no barcode geometry + no TSO + no whitelist overlap".into(),
            llr_contribution: None,
            conclusion: "All three single-cell signals are absent. This is bulk or genomic \
                         sequencing, not single-cell — no chemistry is reported."
                .into(),
        },
    ]
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

fn extract_cb_candidates(seqs: &[Vec<u8>], offset: usize, len: usize) -> HashSet<Vec<u8>> {
    let mut set = HashSet::new();
    for s in seqs {
        if s.len() >= offset + len {
            let mut cb = s[offset..offset + len].to_vec();
            for b in &mut cb {
                if b.is_ascii_lowercase() {
                    *b -= 32;
                }
            }
            set.insert(cb);
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

/// Read one barcode per line, keying on the uppercased first `cb_len` bytes.
fn parse_whitelist<R: std::io::BufRead>(reader: R, cb_len: usize) -> Result<HashSet<Vec<u8>>> {
    // 10x onlists run to several million barcodes; preallocate to avoid rehashing.
    let mut set = HashSet::with_capacity(4_000_000);
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.len() < cb_len {
            continue;
        }
        let mut k: Vec<u8> = line.bytes().take(cb_len).collect();
        for b in &mut k {
            if b.is_ascii_lowercase() {
                *b -= 32;
            }
        }
        set.insert(k);
    }
    Ok(set)
}

fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    if needle.is_empty() || haystack.len() < needle.len() {
        return false;
    }
    haystack.windows(needle.len()).any(|w| w == needle)
}

impl ChemReport {
    pub fn print(&self) {
        println!("{}", "Chemistry detection".bold().cyan());
        println!("  R1 dominant length:  {} bp", self.r1_length_dominant);
        if !self.r1_length_distribution.is_empty() {
            print!("  R1 length top-5:    ");
            for (l, c) in &self.r1_length_distribution {
                print!("{}bp={}  ", l, c);
            }
            println!();
        }
        println!("  Reads sampled:       {}", self.reads_sampled);

        let conf = self.confidence.as_str();
        let color = |s: &str| -> ColoredString {
            match conf {
                "high" => s.green(),
                "medium" => s.yellow(),
                _ => s.red(),
            }
        };
        println!(
            "  {} {}  ({})",
            "ASSAY:".bold(),
            color(&self.assay_display).bold(),
            color(conf)
        );
        println!(
            "  Geometry:  CB=R1[{}..{}]  UMI=R1[{}..{}]",
            self.cb_offset,
            self.cb_offset + self.cb_len,
            self.umi_offset,
            self.umi_offset + self.umi_len,
        );
        if let Some(wl) = &self.whitelist {
            println!("  Whitelist: {}", wl);
        }
        println!("  CB candidates (distinct): {}", self.n_cb_candidates);
        if !self.whitelist_scores.is_empty() {
            println!("  Whitelist overlap:");
            for ws in &self.whitelist_scores {
                println!(
                    "    {:<40}  {:>10} hits  {:>6.1}%",
                    ws.filename, ws.hits, ws.pct
                );
            }
        }
        for a in &self.anchor_scores {
            println!("  Anchor [{}...]: {:.1}%", a.motif, a.hit_pct);
        }
        if !self.interpretation.is_empty() {
            println!("  {}", "Why this call:".bold());
            for (i, e) in self.interpretation.iter().enumerate() {
                let llr = e
                    .llr_contribution
                    .map(|v| format!("  [LLR {:+.1}]", v).dimmed().to_string())
                    .unwrap_or_default();
                println!("    {}. {}{}", i + 1, e.signal.bold(), llr);
                println!("       observed: {}", e.observation);
                println!("       {} {}", "implies:".cyan(), e.conclusion);
            }
        }
        if !self.tool_hints.is_empty() {
            println!("  Tool hints:");
            for h in &self.tool_hints {
                println!("    {:12}  {}", h.tool, h.hint);
            }
        }
        for note in &self.notes {
            println!("  {} {}", "note:".dimmed(), note);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A distinct 16-bp cell barcode for index `i` (base-4 encoding of `i`).
    fn distinct_cb(i: usize) -> Vec<u8> {
        (0..16)
            .map(|j| b"ACGT"[(i / 4usize.pow(j as u32)) % 4] as char as u8)
            .collect()
    }

    /// `n` synthetic 10x R1 reads (16-bp CB + 12-bp UMI = 28 bp) and matching R2 reads
    /// (TSO at the start when `tso` is set). Returns `(r1, r2, cb_set)`.
    #[allow(clippy::type_complexity)]
    fn tenx_v3_reads(n: usize, tso: bool) -> (Vec<Vec<u8>>, Vec<Vec<u8>>, HashSet<Vec<u8>>) {
        let mut r1 = Vec::with_capacity(n);
        let mut cbs = HashSet::new();
        for i in 0..n {
            let cb = distinct_cb(i);
            let mut read = cb.clone();
            read.extend_from_slice(&[b'A'; 12]); // UMI
            r1.push(read);
            cbs.insert(cb);
        }
        let r2: Vec<Vec<u8>> = (0..n)
            .map(|_| {
                let mut read = Vec::new();
                if tso {
                    read.extend_from_slice(TSO_PRESENT.seq);
                }
                read.resize(90, b'C');
                read
            })
            .collect();
        (r1, r2, cbs)
    }

    #[test]
    fn tenx_v3_with_matching_whitelist_is_high_confidence() {
        let (r1, r2, cbs) = tenx_v3_reads(200, true);
        let mut wl_sets: HashMap<&'static str, HashSet<Vec<u8>>> = HashMap::new();
        wl_sets.insert("3M-february-2018.txt.gz", cbs);

        let a = analyze(&r1, &r2, &wl_sets, None);
        assert_eq!(
            a.scored[0].0.id, "10x_3prime_v3",
            "winner was {}",
            a.scored[0].0.id
        );
        assert_eq!(a.confidence, "high");
    }

    #[test]
    fn tenx_3prime_without_whitelist_picks_3prime_family_at_medium() {
        // TSO present but no whitelist cache: should land on a 3' assay, medium.
        let (r1, r2, _) = tenx_v3_reads(200, true);
        let a = analyze(&r1, &r2, &HashMap::new(), Some(200));
        assert!(
            a.scored[0].0.id.starts_with("10x_3prime"),
            "winner was {}",
            a.scored[0].0.id
        );
        assert_eq!(a.confidence, "medium");
    }

    #[test]
    fn dropseq_geometry_is_medium_and_correct() {
        // 20 bp R1 (12 CB + 8 UMI), no TSO, no whitelist -> Drop-seq, geometry-only.
        let r1: Vec<Vec<u8>> = (0..200)
            .map(|i| {
                let mut v = distinct_cb(i); // 16 distinct bases
                v.resize(20, b'A'); // pad to the 20 bp Drop-seq geometry
                v
            })
            .collect();
        let r2: Vec<Vec<u8>> = (0..200).map(|_| vec![b'C'; 90]).collect();
        let a = analyze(&r1, &r2, &HashMap::new(), Some(0));
        assert_eq!(
            a.scored[0].0.id, "drop_seq",
            "winner was {}",
            a.scored[0].0.id
        );
        assert_eq!(a.confidence, "medium");
    }

    #[test]
    fn bulk_wgs_reads_are_not_single_cell() {
        // 100 bp genomic R1, no TSO, no whitelist: must NOT be force-fit to a 10x assay.
        let r1: Vec<Vec<u8>> = (0..200).map(|_| vec![b'A'; 100]).collect();
        let r2: Vec<Vec<u8>> = (0..200).map(|_| vec![b'C'; 100]).collect();
        let a = analyze(&r1, &r2, &HashMap::new(), Some(0));
        assert!(
            !a.is_single_cell,
            "100bp genomic reads should not look single-cell"
        );
        let report = build_report(a, false, 200, false);
        assert_eq!(report.assay_id, "none");
        assert_eq!(report.assay_display, "Not single-cell");
        assert_eq!(report.confidence, "low");
        // The ledger walks all three signals and ends on a clear verdict.
        assert!(
            report
                .interpretation
                .iter()
                .any(|e| e.signal == "Barcode geometry"),
            "reasoning should explain the barcode-geometry check"
        );
        assert!(
            report
                .interpretation
                .last()
                .is_some_and(|e| e.signal == "Verdict" && e.conclusion.contains("not single-cell")),
            "reasoning should end on a not-single-cell verdict"
        );
        assert!(report.notes.iter().any(|n| n.contains("not single-cell")));
    }

    #[test]
    fn ambiguous_3p_5p_polarity_is_reported_as_both() {
        // 26 bp R1 (CB16+UMI10) whose barcodes match the 737K onlist, no TSO. 3'v2 and
        // 5'v1 share this geometry + onlist, so the polarity must be reported as
        // unresolved (both names), never a confident single pick.
        let n = 200;
        let mut r1 = Vec::new();
        let mut cbs = HashSet::new();
        for i in 0..n {
            let cb = distinct_cb(i);
            let mut read = cb.clone();
            read.extend_from_slice(&[b'A'; 10]); // UMI -> 26 bp total
            r1.push(read);
            cbs.insert(cb);
        }
        let r2: Vec<Vec<u8>> = (0..n).map(|_| vec![b'C'; 90]).collect();
        let mut wl_sets: HashMap<&'static str, HashSet<Vec<u8>>> = HashMap::new();
        wl_sets.insert("737K-august-2016.txt.gz", cbs);

        let a = analyze(&r1, &r2, &wl_sets, None);
        let report = build_report(a, false, n, false);
        assert!(
            report.assay_display.contains(" or "),
            "polarity should be reported as both, got {}",
            report.assay_display
        );
        assert_ne!(
            report.confidence, "high",
            "unresolved polarity can't be high"
        );
        assert!(report
            .interpretation
            .iter()
            .any(|e| e.signal == "3' vs 5' polarity"));
    }

    #[test]
    fn tenx_v3_reads_are_single_cell() {
        // Guard against the gate over-triggering on genuine single-cell geometry.
        let (r1, r2, _) = tenx_v3_reads(200, true);
        let a = analyze(&r1, &r2, &HashMap::new(), Some(180));
        assert!(a.is_single_cell);
    }

    #[test]
    fn candidates_for_includes_within_tolerance() {
        let ids: Vec<&str> = candidates_for(26).iter().map(|s| s.id).collect();
        assert!(ids.contains(&"10x_3prime_v2")); // exact 26 bp
        assert!(ids.contains(&"10x_3prime_v3")); // 28 bp, within +/-2
        assert!(!ids.contains(&"drop_seq")); // 20 bp, outside tolerance
    }

    #[test]
    fn len_profile_reports_the_mode() {
        let seqs = vec![vec![0u8; 28], vec![0u8; 28], vec![0u8; 26]];
        let (dom, dist) = len_profile(&seqs);
        assert_eq!(dom, 28);
        assert_eq!(dist[0], (28, 2));
    }
}
