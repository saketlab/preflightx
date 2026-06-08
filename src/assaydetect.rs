//! Generic, seqspec-driven assay detection.
//!
//! Unlike `chemistry.rs` (hand-tuned for the 10x family + Drop-seq), this scores
//! *any* assay describable in seqspec by combining two read-derived signals against
//! each catalogued assay's structural read:
//!
//! 1. **Geometry** — does the observed barcode-read length match the assay's
//!    expected structural length?
//! 2. **Linker anchors** — do the assay's fixed inter-barcode linker/adapter
//!    sequences appear at their expected offsets? This generalizes the 10x TSO
//!    anchor and is what makes combinatorial assays (SPLiT-seq, inDrop, PIP-seq,
//!    sci-RNA-seq, SHARE-seq) detectable *without* a barcode whitelist.
//!
//! Both signals are log-likelihood ratios (nats), summed and ranked, mirroring the
//! model in `chemistry.rs`.

use std::path::Path;

use anyhow::Result;
use colored::Colorize;
use serde::Serialize;

use crate::fastq;
use crate::libstruct::{embedded_assays, Assay, Segment};
use crate::stats::{bernoulli_llr, dominant_len};

// Anchor model: a real fixed linker is present in ~80% of eligible reads (allowing
// for sequencing error); a 6+ bp exact-ish motif at a fixed offset has a tiny
// chance rate. Same shape as the TSO model in chemistry.rs.
const ANCHOR_P_SIGNAL: f64 = 0.80;
const ANCHOR_P_NULL: f64 = 0.001;
const GEOM_BONUS: f64 = 10.0;
const GEOM_PENALTY: f64 = -50.0;
const MIN_ANCHOR_LEN: usize = 6;

/// One ranked detection result.
#[derive(Debug, Clone, Serialize)]
pub struct AssayHit {
    pub spec_file: String,
    pub assay_id: String,
    pub name: String,
    pub modalities: Vec<String>,
    pub lib_struct: String,
    pub structural_read_id: String,
    /// Which input file the barcode read was matched to.
    pub orientation: &'static str,
    pub expected_len: usize,
    pub observed_len: usize,
    pub anchors_found: usize,
    pub anchors_total: usize,
    /// Total log-likelihood ratio (nats); higher is a better match.
    pub score: f64,
}

/// Full generic-detection report for the `detect-assay` subcommand.
#[derive(Debug, Clone, Serialize)]
pub struct DetectReport {
    pub reads_sampled: usize,
    pub confidence: String,
    pub hits: Vec<AssayHit>,
}

/// Sample both files and rank the catalogued assays against them.
pub fn detect(r1: &Path, r2: &Path, n: usize) -> Result<DetectReport> {
    let (r1_seqs, _) = fastq::sample(r1, n)?;
    let (r2_seqs, _) = fastq::sample(r2, n)?;
    let assays = embedded_assays();
    let hits = detect_assays(&r1_seqs, &r2_seqs, &assays, 5);

    // Confidence from the per-read margin between the top two hits.
    let n_reads = r1_seqs.len().max(1) as f64;
    let confidence = match (hits.first(), hits.get(1)) {
        (Some(a), Some(b)) => {
            let gap = (a.score - b.score) / n_reads;
            if a.anchors_found >= 2 && gap > 0.30 {
                "high"
            } else if gap > 0.05 || a.anchors_found >= 1 {
                "medium"
            } else {
                "low"
            }
        }
        (Some(a), None) => {
            if a.anchors_found >= 1 {
                "medium"
            } else {
                "low"
            }
        }
        _ => "low",
    };

    Ok(DetectReport {
        reads_sampled: r1_seqs.len(),
        confidence: confidence.to_string(),
        hits,
    })
}

/// A one-line catalogue entry for the `assays` subcommand.
#[derive(Debug, Clone, Serialize)]
pub struct CatalogueEntry {
    pub spec_file: String,
    pub assay_id: String,
    pub name: String,
    pub modalities: Vec<String>,
    pub detectable: bool,
    pub anchors: usize,
    pub lib_struct: String,
}

/// List the embedded assay catalogue.
pub fn catalogue() -> Vec<CatalogueEntry> {
    let mut out: Vec<CatalogueEntry> = embedded_assays()
        .iter()
        .map(|a| {
            let (detectable, anchors) = match a.structural_read() {
                Some((_, layout)) => (true, anchors_of(&layout).len()),
                None => (false, 0),
            };
            CatalogueEntry {
                spec_file: a.spec_file.clone(),
                assay_id: a.assay_id.clone(),
                name: a.name.clone(),
                modalities: a.modalities.clone(),
                detectable,
                anchors,
                lib_struct: a.lib_struct.clone(),
            }
        })
        .collect();
    out.sort_by(|a, b| a.spec_file.cmp(&b.spec_file));
    out
}

/// Find a catalogued assay by spec filename, its stem, or assay_id (case-insensitive).
pub fn find_assay(query: &str) -> Option<Assay> {
    embedded_assays().into_iter().find(|a| {
        a.spec_file.eq_ignore_ascii_case(query)
            || a.spec_file
                .trim_end_matches(".spec.yaml")
                .eq_ignore_ascii_case(query)
            || a.assay_id.eq_ignore_ascii_case(query)
    })
}

/// A fixed motif expected at a known offset in the structural read.
struct Anchor {
    offset: usize,
    seq: Vec<u8>,
}

fn anchors_of(layout: &[Segment]) -> Vec<Anchor> {
    let mut out = Vec::new();
    let mut offset = 0usize;
    for s in layout {
        if let Some(seq) = &s.sequence {
            if seq.len() >= MIN_ANCHOR_LEN {
                out.push(Anchor {
                    offset,
                    seq: seq.as_bytes().to_vec(),
                });
            }
        }
        offset += s.len;
    }
    out
}

fn hamming_le(a: &[u8], b: &[u8], max_mismatch: usize) -> bool {
    let mut m = 0;
    for (x, y) in a.iter().zip(b.iter()) {
        if !x.eq_ignore_ascii_case(y) {
            m += 1;
            if m > max_mismatch {
                return false;
            }
        }
    }
    true
}

/// Score one observed read set against a structural layout. Returns
/// `(total_llr, anchors_found, dominant_len)`.
fn score_obs(obs: &[Vec<u8>], expected_len: usize, anchors: &[Anchor]) -> (f64, usize, usize) {
    if obs.is_empty() {
        return (f64::NEG_INFINITY, 0, 0);
    }
    let dominant = dominant_len(obs);
    let tol = (expected_len / 15).max(3);
    let diff = dominant.abs_diff(expected_len);
    let geometry_ok = expected_len > 0 && diff <= tol;
    // Graded: an exact length match scores the full bonus; each bp of deviation
    // within tolerance shaves it down, so the closest-fitting assay ranks highest.
    let mut total = if geometry_ok {
        GEOM_BONUS - diff as f64
    } else {
        GEOM_PENALTY
    };

    let mut found = 0;
    for a in anchors {
        let end = a.offset + a.seq.len();
        let max_mm = (a.seq.len() * 15 / 100).max(1);
        let mut eligible = 0usize;
        let mut hits = 0usize;
        for r in obs {
            if r.len() >= end {
                eligible += 1;
                if hamming_le(&r[a.offset..end], &a.seq, max_mm) {
                    hits += 1;
                }
            }
        }
        if eligible > 0 {
            total += bernoulli_llr(hits, eligible, ANCHOR_P_SIGNAL, ANCHOR_P_NULL);
            if hits * 100 >= eligible * 30 {
                found += 1;
            }
        }
    }
    (total, found, dominant)
}

/// Rank catalogued assays by how well the observed reads match each one. The two
/// input files are tried in both orientations (either can carry the barcode read).
pub fn detect_assays(
    r1: &[Vec<u8>],
    r2: &[Vec<u8>],
    assays: &[Assay],
    max_hits: usize,
) -> Vec<AssayHit> {
    let mut hits: Vec<AssayHit> = Vec::new();

    for a in assays {
        let Some((read_id, layout)) = a.structural_read() else {
            continue;
        };
        let expected: usize = layout.iter().map(|s| s.len).sum();
        let anchors = anchors_of(&layout);

        let mut best: Option<AssayHit> = None;
        for (orient, obs) in [("r1=barcode", r1), ("r2=barcode", r2)] {
            let (score, found, dom) = score_obs(obs, expected, &anchors);
            let candidate = AssayHit {
                spec_file: a.spec_file.clone(),
                assay_id: a.assay_id.clone(),
                name: a.name.clone(),
                modalities: a.modalities.clone(),
                lib_struct: a.lib_struct.clone(),
                structural_read_id: read_id.clone(),
                orientation: orient,
                expected_len: expected,
                observed_len: dom,
                anchors_found: found,
                anchors_total: anchors.len(),
                score,
            };
            if best
                .as_ref()
                .map(|b| candidate.score > b.score)
                .unwrap_or(true)
            {
                best = Some(candidate);
            }
        }
        if let Some(h) = best {
            hits.push(h);
        }
    }

    hits.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    hits.truncate(max_hits);
    hits
}

impl DetectReport {
    pub fn print(&self) {
        println!("{}", "Generic assay detection (seqspec)".bold().cyan());
        println!("  reads sampled: {}", self.reads_sampled);
        let conf = match self.confidence.as_str() {
            "high" => self.confidence.green(),
            "medium" => self.confidence.yellow(),
            _ => self.confidence.red(),
        };
        println!("  {} {}", "CONFIDENCE:".bold(), conf.bold());
        if self.hits.is_empty() {
            println!("  no catalogued assay matched the reads.");
            return;
        }
        println!("  {}", "ranked candidates:".bold());
        for (i, h) in self.hits.iter().enumerate() {
            let marker = if i == 0 { "→".green() } else { " ".normal() };
            let label = if h.name.is_empty() {
                &h.assay_id
            } else {
                &h.name
            };
            println!(
                "    {} {:<28} score={:>8.1}  anchors={}/{}  {} (exp {}bp, obs {}bp)",
                marker,
                label,
                h.score,
                h.anchors_found,
                h.anchors_total,
                h.orientation,
                h.expected_len,
                h.observed_len,
            );
        }
    }
}

/// Render an ASCII description of one assay: its structural read layout, the
/// barcode/UMI offsets, fixed linker anchors, and onlist references.
pub fn describe(a: &Assay) -> String {
    let mut s = String::new();
    s.push_str(&format!("assay: {} ({})\n", a.name, a.spec_file));
    if !a.modalities.is_empty() {
        s.push_str(&format!("modalities: {}\n", a.modalities.join(", ")));
    }
    if !a.lib_struct.is_empty() {
        s.push_str(&format!("lib_struct: {}\n", a.lib_struct));
    }
    match a.structural_read() {
        None => {
            s.push_str("structural read: none (no in-read barcode — plate-based or bulk-like)\n")
        }
        Some((rid, layout)) => {
            s.push_str(&format!("structural read: {rid}\n"));
            let mut offset = 0usize;
            for seg in &layout {
                let onlist = seg
                    .onlist
                    .as_ref()
                    .map(|o| format!("  onlist={}", o.filename))
                    .unwrap_or_default();
                let seq = seg
                    .sequence
                    .as_ref()
                    .map(|q| format!("  seq={q}"))
                    .unwrap_or_default();
                s.push_str(&format!(
                    "  [{:>3}..{:>3}] {:<10} {}{}{}\n",
                    offset,
                    offset + seg.len,
                    seg.region_type,
                    seg.region_id,
                    onlist,
                    seq,
                ));
                offset += seg.len;
            }
        }
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build one synthetic read matching a structural layout: fixed segments use
    /// their real sequence, variable barcode/UMI segments are filled with A's
    /// (content is irrelevant — only the linker anchors are scored).
    fn synth_read(layout: &[Segment]) -> Vec<u8> {
        let mut out = Vec::new();
        for s in layout {
            match &s.sequence {
                Some(seq) => out.extend_from_slice(seq.as_bytes()),
                None => out.extend(std::iter::repeat(b'A').take(s.len)),
            }
        }
        out
    }

    #[test]
    fn detects_combinatorial_assay_via_linkers_without_whitelist() {
        let assays = embedded_assays();
        let split = Assay::from_yaml(include_str!(
            "../vendor/seqspec/docs/examples/assays/split_seq.spec.yaml"
        ))
        .unwrap();
        let (_rid, layout) = split.structural_read().unwrap();

        // SPLiT-seq barcode read is R2; synthesize 500 of them, with random cDNA on R1.
        let r2: Vec<Vec<u8>> = (0..500).map(|_| synth_read(&layout)).collect();
        let r1: Vec<Vec<u8>> = (0..500).map(|_| vec![b'C'; 100]).collect();

        let hits = detect_assays(&r1, &r2, &assays, 5);
        assert!(!hits.is_empty());
        let top = &hits[0];
        assert_eq!(
            top.spec_file, "split_seq.spec.yaml",
            "top hit was {}",
            top.spec_file
        );
        assert_eq!(top.orientation, "r2=barcode");
        assert!(
            top.anchors_found >= 2,
            "expected >=2 linkers found, got {}",
            top.anchors_found
        );
    }

    #[test]
    fn exact_geometry_outranks_within_tolerance() {
        // A 28 bp barcode read (10x v3/v4 geometry) with no linkers should rank an
        // exact-length assay above one that is merely within the ±tolerance window
        // (e.g. 26 bp v2). Graded geometry, not a binary pass/fail.
        let assays = embedded_assays();
        let r1: Vec<Vec<u8>> = (0..300)
            .map(|i| (0..28).map(|j| b"ACGT"[(i + j) % 4]).collect())
            .collect();
        let r2: Vec<Vec<u8>> = (0..300).map(|_| vec![b'G'; 90]).collect();
        let hits = detect_assays(&r1, &r2, &assays, 10);
        assert_eq!(
            hits[0].expected_len, 28,
            "top hit should match the exact 28bp geometry, got {}bp ({})",
            hits[0].expected_len, hits[0].spec_file
        );
    }

    #[test]
    fn handcoded_10x_geometry_matches_seqspec() {
        // The 10x family is defined twice — as hand-coded `AssaySpec`s in assay.rs and
        // as vendored seqspec YAML — so this pins the two together: a future edit to
        // either side that changes a barcode/UMI length will fail here instead of
        // making `chem` and `detect-assay` silently disagree on the same input.
        use crate::assay::ASSAYS;
        let pairs = [
            ("10x_3prime_v1", "10x_rna_v1.spec.yaml"),
            ("10x_3prime_v2", "10x_rna_v2.spec.yaml"),
            ("10x_3prime_v3", "10x_rna_v3.spec.yaml"),
        ];
        for (id, spec_file) in pairs {
            let spec = ASSAYS
                .iter()
                .find(|s| s.id == id)
                .unwrap_or_else(|| panic!("hand-coded assay {id} missing"));
            let assay =
                find_assay(spec_file).unwrap_or_else(|| panic!("seqspec {spec_file} missing"));
            let (_rid, layout) = assay
                .structural_read()
                .unwrap_or_else(|| panic!("{spec_file} has no structural read"));
            let cb: usize = layout
                .iter()
                .filter(|s| s.region_type == "barcode")
                .map(|s| s.len)
                .sum();
            assert_eq!(
                cb, spec.r1_cb_len,
                "{id}: CB length disagrees with {spec_file}"
            );
            // UMI is co-located with the barcode only in some chemistries (v1 puts it on
            // the other read); cross-check it only when the structural read carries it.
            let umi: usize = layout
                .iter()
                .filter(|s| s.region_type == "umi")
                .map(|s| s.len)
                .sum();
            if umi > 0 {
                assert_eq!(
                    umi, spec.r1_umi_len,
                    "{id}: UMI length disagrees with {spec_file}"
                );
            }
        }
    }

    #[test]
    fn self_detection_ranks_every_assay_top_or_tied() {
        // For every catalogued assay that has a barcode-bearing read, reads
        // synthesized from its own structural layout must rank that assay first
        // (strictly, when it has distinguishing linkers; tied at most with other
        // assays of identical geometry and no anchors).
        let assays = embedded_assays();
        let mut checked = 0;
        for a in &assays {
            let Some((_rid, layout)) = a.structural_read() else {
                continue;
            };
            let expected: usize = layout.iter().map(|s| s.len).sum();
            if expected == 0 {
                continue;
            }
            let r1: Vec<Vec<u8>> = (0..300).map(|_| synth_read(&layout)).collect();
            let r2: Vec<Vec<u8>> = (0..300).map(|_| vec![b'G'; 90]).collect();
            let hits = detect_assays(&r1, &r2, &assays, assays.len());
            let top = hits[0].score;
            let ours = hits
                .iter()
                .find(|h| h.spec_file == a.spec_file)
                .expect("our assay must appear in hits");
            assert!(
                (ours.score - top).abs() < 1e-6,
                "{} ranked below top: {} vs {} (top = {})",
                a.spec_file,
                ours.score,
                top,
                hits[0].spec_file,
            );
            checked += 1;
        }
        assert!(checked >= 20, "only checked {checked} assays");
    }
}
