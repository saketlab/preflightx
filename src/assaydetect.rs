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
    /// Human-friendly name with 3'/5' orientation, e.g. "10x 3' v3", "Smart-seq2".
    pub display_name: String,
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
    /// UMI-based (digital counts) vs plate-based / full-length (no UMI, read counts).
    pub umi_based: bool,
    /// Per-barcode onlist overlap `(cache_name, hit_rate)`, empty when no onlist was
    /// loaded for this assay. Populated only when a barcode whitelist is cached.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub onlist_overlap: Vec<(String, f64)>,
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
    let sets = load_cached_onlists();
    let hits = detect_assays(&r1_seqs, &r2_seqs, assays, &sets, 5);

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
    /// Human-friendly name with 3'/5' orientation, e.g. "10x 3' v3", "CEL-seq2".
    pub display_name: String,
    pub name: String,
    pub modalities: Vec<String>,
    pub detectable: bool,
    pub anchors: usize,
    pub lib_struct: String,
    /// UMI-based (digital counts) vs plate-based / full-length (no UMI).
    pub umi_based: bool,
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
                display_name: pretty_name(&a.spec_file),
                name: a.name.clone(),
                modalities: a.modalities.clone(),
                detectable,
                anchors,
                lib_struct: a.lib_struct.clone(),
                umi_based: a.has_umi(),
            }
        })
        .collect();
    out.sort_by(|a, b| a.spec_file.cmp(&b.spec_file));
    out
}

/// Find a catalogued assay by spec filename, its stem, or assay_id (case-insensitive).
pub fn find_assay(query: &str) -> Option<Assay> {
    embedded_assays()
        .iter()
        .find(|a| {
            a.spec_file.eq_ignore_ascii_case(query)
                || a.spec_file
                    .trim_end_matches(".spec.yaml")
                    .eq_ignore_ascii_case(query)
                || a.assay_id.eq_ignore_ascii_case(query)
        })
        .cloned()
}

/// Human-friendly display name for a catalogued assay, from its spec filename — with the
/// 3'/5' orientation kept where it applies (e.g. `10x_rna_v3` → "10x 3' v3"). Falls back to
/// a dash-cased stem for anything not in the curated table.
pub fn pretty_name(spec_file: &str) -> String {
    let stem = spec_file.trim_end_matches(".spec.yaml");
    let known = match stem {
        "10x_rna_v1" => "10x 3' v1",
        "10x_rna_v2" => "10x 3' v2",
        "10x_rna_v3" => "10x 3' v3",
        "10x_rna_5prime" => "10x 5'",
        "10xv3_scrnaseq_element_adept_truseq_dual" => "10x 3' v3 (Element AdeptTruSeq)",
        "10xv3_scrnaseq_illumina_truseq_dual" => "10x 3' v3 (Illumina TruSeq)",
        "10x_atac" => "10x ATAC",
        "10xcrispr" => "10x CRISPR",
        "10xfb_3prime" => "10x Feature Barcode 3'",
        "10xfb_5prime" => "10x Feature Barcode 5'",
        "10xfb_vdj_5prime" => "10x VDJ 5'",
        "10x_rna_atac" => "10x RNA+ATAC",
        "10x_rna_atac_multi" => "10x Multiome (RNA+ATAC)",
        "bd_rhapsody_v1" => "BD Rhapsody",
        "bd_rhapsody_eb" => "BD Rhapsody (Enhanced Bead)",
        "cel_seq" => "CEL-seq",
        "cel_seq2" => "CEL-seq2",
        "dogmaseq_dig" => "DOGMA-seq (DIG)",
        "dogmaseq_lll" => "DOGMA-seq (LLL)",
        "drop_seq" => "Drop-seq",
        "indropv2" => "inDrop v2",
        "issaac_seq" => "ISSAAC-seq",
        "mars_seq" => "MARS-seq",
        "mcscrb_seq" => "mcSCRB-seq",
        "microwell_seq" => "Microwell-seq",
        "parse_wt_v2" => "Parse WT v2",
        "pi_atac_seq" => "Pi-ATAC-seq",
        "pipseqv2" => "PIP-seq v2",
        "pipseqv3" => "PIP-seq v3",
        "quartz_seq" => "Quartz-seq",
        "quartz_seq2" => "Quartz-seq2",
        "sccrispra" => "scCRISPR(a)",
        "scifi_rna_seq" => "scifi-RNA-seq",
        "sci_rna_seq" => "sci-RNA-seq",
        "sci_rna_seq3" => "sci-RNA-seq3",
        "seq_well" => "Seq-Well",
        "seq_well_s3" => "Seq-Well S3",
        "share_seq" => "SHARE-seq",
        "smart_seq2" => "Smart-seq2",
        "smart_seq3" => "Smart-seq3",
        "sn_m3c_seq" => "sn-m3C-seq",
        "snmctseq" => "snmC-seq",
        "split_seq" => "SPLiT-seq",
        "strt_seq" => "STRT-seq",
        "strt_seq_2i" => "STRT-seq-2i",
        "strt_seq_c1" => "STRT-seq-C1",
        "sugarseq" => "SUGAR-seq",
        "surecell" => "SureCell",
        "tang2009" => "Tang 2009",
        "vasa_seq_drop" => "VASA-seq (droplet)",
        "vasa_seq_plate" => "VASA-seq (plate)",
        _ => return stem.replace('_', "-"),
    };
    known.to_string()
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

// Real barcodes land in the onlist ~70% of the time (error + ambient), wrong ones at
// chance. Same model as the whitelist signal in chemistry.rs.
const P_ONLIST_CORRECT: f64 = 0.70;
const P_ONLIST_NULL: f64 = 0.01;

/// A barcode region in the structural read backed by a fetchable onlist.
#[derive(Debug, Clone)]
pub struct OnlistSeg {
    pub offset: usize,
    /// `len..=max_len` give the lengths matched against the onlist; for a variable
    /// barcode (inDrop bc1, 8-11 bp) they differ.
    pub len: usize,
    pub max_len: usize,
    /// Cache name of the resolved onlist, the key into the loaded `sets` map.
    pub filename: String,
}

/// A fetchable onlist: the mirror URL and the cache filename to store it under.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct OnlistSource {
    /// Globally-unique cache filename (the spec's onlist names collide across assays).
    pub cache_name: &'static str,
    pub url: &'static str,
}

/// Onlist sources keyed by `(spec_file, onlist_filename_in_spec)`. seqspec carries no
/// fetch URLs and its onlist filenames collide across assays, so each non-10x barcode
/// list is mapped here to a hosted copy on the scg_lib_structs mirror. Assays not in
/// this table stay geometry-only.
const ONLIST_SOURCES: &[(&str, &str, OnlistSource)] = &[
    (
        "indropv2.spec.yaml",
        "bc1_onlist.txt",
        OnlistSource {
            cache_name: "inDrop_barcode1_list.txt",
            url: "https://teichlab.github.io/scg_lib_structs/data/inDrop/inDrop_barcode1_list.txt",
        },
    ),
    // inDrop bc2 is dropped by the framing rule (downstream of variable bc1), so it is
    // intentionally not listed here — it would only be fetched and never scored.
    (
        "bd_rhapsody_v1.spec.yaml",
        "cls1_onlist.txt",
        OnlistSource {
            cache_name: "BD_CLS1.txt",
            url: "https://teichlab.github.io/scg_lib_structs/data/BD/BD_CLS1.txt",
        },
    ),
    (
        "bd_rhapsody_v1.spec.yaml",
        "cls2_onlist.txt",
        OnlistSource {
            cache_name: "BD_CLS2.txt",
            url: "https://teichlab.github.io/scg_lib_structs/data/BD/BD_CLS2.txt",
        },
    ),
    (
        "bd_rhapsody_v1.spec.yaml",
        "cls3_onlist.txt",
        OnlistSource {
            cache_name: "BD_CLS3.txt",
            url: "https://teichlab.github.io/scg_lib_structs/data/BD/BD_CLS3.txt",
        },
    ),
    (
        "pipseqv2.spec.yaml",
        "pip-seq_v2_bc1.tsv",
        OnlistSource {
            cache_name: "pip-seq_v2_bc1.tsv",
            url: "https://teichlab.github.io/scg_lib_structs/data/PIP-seq/pip-seq_v2_bc1.tsv",
        },
    ),
    (
        "pipseqv2.spec.yaml",
        "pip-seq_v2_bc2.tsv",
        OnlistSource {
            cache_name: "pip-seq_v2_bc2.tsv",
            url: "https://teichlab.github.io/scg_lib_structs/data/PIP-seq/pip-seq_v2_bc2.tsv",
        },
    ),
    (
        "pipseqv2.spec.yaml",
        "pip-seq_v2_bc3.tsv",
        OnlistSource {
            cache_name: "pip-seq_v2_bc3.tsv",
            url: "https://teichlab.github.io/scg_lib_structs/data/PIP-seq/pip-seq_v2_bc3.tsv",
        },
    ),
    (
        "pipseqv3.spec.yaml",
        "fb_v3_bc1.tsv",
        OnlistSource {
            cache_name: "fb_v3_bc1.tsv",
            url: "https://teichlab.github.io/scg_lib_structs/data/PIP-seq/fb_v3_bc1.tsv",
        },
    ),
    (
        "pipseqv3.spec.yaml",
        "fb_v3_bc2.tsv",
        OnlistSource {
            cache_name: "fb_v3_bc2.tsv",
            url: "https://teichlab.github.io/scg_lib_structs/data/PIP-seq/fb_v3_bc2.tsv",
        },
    ),
    (
        "pipseqv3.spec.yaml",
        "fb_v3_bc3.tsv",
        OnlistSource {
            cache_name: "fb_v3_bc3.tsv",
            url: "https://teichlab.github.io/scg_lib_structs/data/PIP-seq/fb_v3_bc3.tsv",
        },
    ),
    (
        "pipseqv3.spec.yaml",
        "fb_v3_bc4.tsv",
        OnlistSource {
            cache_name: "fb_v3_bc4.tsv",
            url: "https://teichlab.github.io/scg_lib_structs/data/PIP-seq/fb_v3_bc4.tsv",
        },
    ),
    (
        "quartz_seq2.spec.yaml",
        "cell_bc_onlist.txt",
        OnlistSource {
            cache_name: "Quartz-seq2_cell_barcodes.txt",
            url: "https://teichlab.github.io/scg_lib_structs/data/Quartz-seq2_cell_barcodes.txt",
        },
    ),
];

/// Resolve a spec-local onlist filename to a fetchable source, or `None` if preflightx
/// has no hosted copy for it (the assay stays geometry-only).
fn onlist_source(spec_file: &str, onlist_filename: &str) -> Option<OnlistSource> {
    ONLIST_SOURCES
        .iter()
        .find(|(s, f, _)| *s == spec_file && *f == onlist_filename)
        .map(|(_, _, src)| *src)
}

/// Every distinct onlist source preflightx can fetch (for `fetch-whitelists`).
pub fn onlist_fetch_list() -> Vec<OnlistSource> {
    let mut seen = std::collections::HashSet::new();
    ONLIST_SOURCES
        .iter()
        .filter(|(_, _, s)| seen.insert(s.cache_name))
        .map(|(_, _, s)| *s)
        .collect()
}

/// Load every cached onlist (keyed by cache name) so it can be scored. Files absent
/// from the cache are skipped, leaving those assays geometry-only.
#[cfg(not(target_arch = "wasm32"))]
fn load_cached_onlists() -> std::collections::HashMap<String, std::collections::HashSet<Vec<u8>>> {
    let mut sets = std::collections::HashMap::new();
    let Ok(dir) = crate::whitelists::cache_dir() else {
        return sets;
    };
    for src in onlist_fetch_list() {
        let path = dir.join(src.cache_name);
        if let Ok(text) = std::fs::read_to_string(&path) {
            sets.insert(src.cache_name.to_string(), parse_onlist(&text));
        }
    }
    sets
}

#[cfg(target_arch = "wasm32")]
fn load_cached_onlists() -> std::collections::HashMap<String, std::collections::HashSet<Vec<u8>>> {
    let mut sets = std::collections::HashMap::new();
    for src in onlist_fetch_list() {
        // Key must match where `wasm::load_whitelist` stashes onlists.
        let key = format!("{}/{}", crate::wasm::WL_DIR, src.cache_name);
        if let Some(bytes) = fastq::vfs::get(&key) {
            if let Ok(text) = std::str::from_utf8(bytes.as_ref()) {
                sets.insert(src.cache_name.to_string(), parse_onlist(&text));
            }
        }
    }
    sets
}

/// One barcode per line (first whitespace/comma field), uppercased. Handles `.txt`
/// and `.tsv`/`.csv` onlists; blank and comment lines are skipped.
fn parse_onlist(text: &str) -> std::collections::HashSet<Vec<u8>> {
    text.lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(|l| {
            l.split(['\t', ','])
                .next()
                .unwrap_or(l)
                .bytes()
                .map(|b| b.to_ascii_uppercase())
                .collect()
        })
        .collect()
}

fn set_contains_uppercase(
    set: &std::collections::HashSet<Vec<u8>>,
    bytes: &[u8],
    scratch: &mut Vec<u8>,
) -> bool {
    if bytes.iter().any(|b| b.is_ascii_lowercase()) {
        scratch.clear();
        scratch.extend(bytes.iter().map(|b| b.to_ascii_uppercase()));
        set.contains(scratch.as_slice())
    } else {
        set.contains(bytes)
    }
}

fn revcomp_upper_into(seq: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.reserve(seq.len());
    for &b in seq.iter().rev() {
        out.push(match b {
            65 | 97 => 84,
            67 | 99 => 71,
            71 | 103 => 67,
            84 | 116 => 65,
            _ => 78,
        });
    }
}

/// Summed onlist-overlap LLR (nats) plus a per-segment `(filename, hit_rate)`. Each
/// barcode is matched forward and reverse-complemented (sequencer-dependent), and at
/// every length for variable-length barcodes; the better orientation is scored.
/// Segments whose onlist is not loaded contribute nothing.
fn onlist_overlap_llr(
    obs: &[Vec<u8>],
    segs: &[OnlistSeg],
    sets: &std::collections::HashMap<String, std::collections::HashSet<Vec<u8>>>,
) -> (f64, Vec<(String, f64)>) {
    let mut total = 0.0_f64;
    let mut per_seg = Vec::new();
    for seg in segs {
        let Some(set) = sets.get(&seg.filename) else {
            continue;
        };
        let (mut eligible, mut fwd, mut rev) = (0usize, 0usize, 0usize);
        let mut fwd_scratch = Vec::new();
        let mut rev_scratch = Vec::new();
        for r in obs {
            if r.len() < seg.offset + seg.len {
                continue;
            }
            eligible += 1;
            // A variable-length barcode matches at whichever length is on-list.
            let mut hit_fwd = false;
            let mut hit_rev = false;
            for l in seg.len..=seg.max_len {
                if r.len() < seg.offset + l {
                    break;
                }
                let bc = &r[seg.offset..seg.offset + l];
                hit_fwd |= set_contains_uppercase(set, bc, &mut fwd_scratch);
                revcomp_upper_into(bc, &mut rev_scratch);
                hit_rev |= set.contains(rev_scratch.as_slice());
            }
            fwd += hit_fwd as usize;
            rev += hit_rev as usize;
        }
        if eligible == 0 {
            continue;
        }
        let hits = fwd.max(rev);
        total += bernoulli_llr(hits, eligible, P_ONLIST_CORRECT, P_ONLIST_NULL);
        per_seg.push((seg.filename.clone(), hits as f64 / eligible as f64));
    }
    (total, per_seg)
}

/// Onlist segments for `a`, offset to align with the reads `score_obs` sees.
/// `barcode_onlist_segments` counts from the first barcode/UMI, but the observed read
/// (like `anchors_of`) counts from the layout start, so shift by the leading length.
/// Keeps only segments with a hosted onlist, re-keyed to the cache name.
fn resolved_onlist_segs(a: &Assay, layout: &[Segment]) -> Vec<OnlistSeg> {
    let leading: usize = layout
        .iter()
        .take_while(|s| !crate::libstruct::is_barcode_or_umi(&s.region_type))
        .map(|s| s.len)
        .sum();
    crate::libstruct::barcode_onlist_segments_of(layout)
        .into_iter()
        .filter_map(|(off, len, max_len, local)| {
            onlist_source(&a.spec_file, &local).map(|src| OnlistSeg {
                offset: leading + off,
                len,
                max_len,
                filename: src.cache_name.to_string(),
            })
        })
        .collect()
}

/// Rank catalogued assays by how well the observed reads match each one. The two
/// input files are tried in both orientations (either can carry the barcode read).
pub fn detect_assays(
    r1: &[Vec<u8>],
    r2: &[Vec<u8>],
    assays: &[Assay],
    sets: &std::collections::HashMap<String, std::collections::HashSet<Vec<u8>>>,
    max_hits: usize,
) -> Vec<AssayHit> {
    let mut hits: Vec<AssayHit> = Vec::new();

    for a in assays {
        let Some((read_id, layout)) = a.structural_read() else {
            continue;
        };
        let expected: usize = layout.iter().map(|s| s.len).sum();
        let anchors = anchors_of(&layout);
        let onlist_segs = resolved_onlist_segs(a, &layout);

        let mut best: Option<AssayHit> = None;
        for (orient, obs) in [("r1=barcode", r1), ("r2=barcode", r2)] {
            let (mut score, found, dom) = score_obs(obs, expected, &anchors);
            let (ol_llr, overlap) = onlist_overlap_llr(obs, &onlist_segs, sets);
            score += ol_llr;
            let candidate = AssayHit {
                spec_file: a.spec_file.clone(),
                assay_id: a.assay_id.clone(),
                display_name: pretty_name(&a.spec_file),
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
                umi_based: a.has_umi(),
                onlist_overlap: overlap,
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
        let q = if self.hits[0].umi_based {
            "UMI-based (digital molecule counts)".green()
        } else {
            "plate-based / full-length — no UMI (read counts)".yellow()
        };
        println!("  {} {}", "QUANTIFICATION (top candidate):".bold(), q);
        println!("  {}", "ranked candidates:".bold());
        for (i, h) in self.hits.iter().enumerate() {
            let marker = if i == 0 { "→".green() } else { " ".normal() };
            println!(
                "    {} {:<26} score={:>8.1}  anchors={}/{}  {} (exp {}bp, obs {}bp)",
                marker,
                h.display_name,
                h.score,
                h.anchors_found,
                h.anchors_total,
                h.orientation,
                h.expected_len,
                h.observed_len,
            );
        }
        let top = &self.hits[0];
        if !top.onlist_overlap.is_empty() {
            let parts: Vec<String> = top
                .onlist_overlap
                .iter()
                .map(|(f, rate)| format!("{f} {:.0}%", rate * 100.0))
                .collect();
            println!(
                "  {} {}",
                "Onlist overlap (top candidate):".bold(),
                parts.join(", ")
            );
        }
        if let Some(a) = find_assay(&top.spec_file) {
            println!(
                "  {} {}",
                "STARsolo (top candidate):".bold(),
                crate::libstruct::starsolo(&a)
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
    s.push_str(&format!("STARsolo: {}\n", crate::libstruct::starsolo(a)));
    s
}

/// Structured form of [`describe`] for `--json`.
pub fn describe_json(a: &Assay) -> serde_json::Value {
    let structural = a.structural_read().map(|(rid, layout)| {
        let mut offset = 0usize;
        let segments: Vec<_> = layout
            .iter()
            .map(|seg| {
                let start = offset;
                offset += seg.len;
                serde_json::json!({
                    "start": start,
                    "end": offset,
                    "region_type": seg.region_type,
                    "region_id": seg.region_id,
                    "onlist": seg.onlist.as_ref().map(|o| &o.filename),
                    "sequence": seg.sequence,
                })
            })
            .collect();
        serde_json::json!({ "read_id": rid, "segments": segments })
    });
    serde_json::json!({
        "assay_id": a.assay_id,
        "spec_file": a.spec_file,
        "name": a.name,
        "modalities": a.modalities,
        "lib_struct": a.lib_struct,
        "structural_read": structural,
        "starsolo": crate::libstruct::starsolo(a),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pretty_names_with_orientation() {
        assert_eq!(pretty_name("10x_rna_v3.spec.yaml"), "10x 3' v3");
        assert_eq!(pretty_name("10x_rna_5prime.spec.yaml"), "10x 5'");
        assert_eq!(pretty_name("smart_seq2.spec.yaml"), "Smart-seq2");
        assert_eq!(pretty_name("cel_seq2.spec.yaml"), "CEL-seq2");
        assert_eq!(pretty_name("split_seq.spec.yaml"), "SPLiT-seq");
        // every catalogued spec yields a non-empty, non-underscored display name
        for e in catalogue() {
            assert!(!e.display_name.is_empty());
            assert!(!e.display_name.contains('_'), "{}", e.display_name);
        }
    }

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

        let hits = detect_assays(&r1, &r2, assays, &Default::default(), 5);
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
    fn onlist_overlap_rewards_in_whitelist_barcodes() {
        use std::collections::{HashMap, HashSet};
        // A small onlist of four 8 bp barcodes.
        let bcs = [
            b"AAAACCCC".to_vec(),
            b"GGGGTTTT".to_vec(),
            b"ACGTACGT".to_vec(),
            b"TGCATGCA".to_vec(),
        ];
        let set: HashSet<Vec<u8>> = bcs.iter().cloned().collect();
        let mut sets = HashMap::new();
        sets.insert("bc1.txt".to_string(), set);

        // 200 reads: an on-list barcode at offset 0, then 12 bp of filler.
        let obs: Vec<Vec<u8>> = (0..200)
            .map(|i| {
                let mut r = bcs[i % bcs.len()].clone();
                r.extend_from_slice(b"GGGGGGGGGGGG");
                r
            })
            .collect();
        let segs = vec![OnlistSeg {
            offset: 0,
            len: 8,
            max_len: 8,
            filename: "bc1.txt".into(),
        }];

        let (llr, per_seg) = onlist_overlap_llr(&obs, &segs, &sets);
        assert!(
            llr > 50.0,
            "on-list barcodes should give a strong positive LLR, got {llr}"
        );
        assert_eq!(per_seg.len(), 1);
        assert!(per_seg[0].1 > 0.9, "hit rate ~1.0, got {}", per_seg[0].1);
    }

    #[test]
    fn onlist_overlap_ignores_off_whitelist_and_missing_sets() {
        use std::collections::{HashMap, HashSet};
        let set: HashSet<Vec<u8>> = [b"AAAACCCC".to_vec(), b"GGGGTTTT".to_vec()]
            .into_iter()
            .collect();
        let mut sets = HashMap::new();
        sets.insert("bc1.txt".to_string(), set);

        // Barcodes absent from the onlist score a strongly negative LLR.
        let obs: Vec<Vec<u8>> = (0..200).map(|_| b"TTTTTTTTNNNNNNNNNNNN".to_vec()).collect();
        let segs = vec![OnlistSeg {
            offset: 0,
            len: 8,
            max_len: 8,
            filename: "bc1.txt".into(),
        }];
        let (llr, _) = onlist_overlap_llr(&obs, &segs, &sets);
        assert!(
            llr < -200.0,
            "200 off-list reads should be strongly penalized, got {llr}"
        );

        // A segment whose onlist was never fetched contributes nothing.
        let missing = vec![OnlistSeg {
            offset: 0,
            len: 8,
            max_len: 8,
            filename: "absent.txt".into(),
        }];
        let (llr0, per_seg0) = onlist_overlap_llr(&obs, &missing, &sets);
        assert_eq!(llr0, 0.0);
        assert!(per_seg0.is_empty());
    }

    #[test]
    fn onlist_overlap_matches_reverse_complement_orientation() {
        use std::collections::{HashMap, HashSet};
        // Onlist stores forward barcodes; reads carry their reverse complement
        // (some sequencers do this for inDrop). Overlap must still fire.
        let set: HashSet<Vec<u8>> = [b"AAAACCCC".to_vec(), b"ACGTACGT".to_vec()]
            .into_iter()
            .collect();
        let mut sets = HashMap::new();
        sets.insert("bc1.txt".to_string(), set);

        let rc = [b"GGGGTTTT".to_vec(), b"ACGTACGT".to_vec()];
        let obs: Vec<Vec<u8>> = (0..200)
            .map(|i| {
                let mut r = rc[i % rc.len()].clone();
                r.extend_from_slice(b"GGGGGGGGGGGG");
                r
            })
            .collect();
        let segs = vec![OnlistSeg {
            offset: 0,
            len: 8,
            max_len: 8,
            filename: "bc1.txt".into(),
        }];
        let (llr, _) = onlist_overlap_llr(&obs, &segs, &sets);
        assert!(llr > 50.0, "revcomp barcodes should still match, got {llr}");
    }

    #[test]
    fn parse_onlist_handles_columns_comments_and_case() {
        let text = "AAAA\nCCCC\tnote\nGGGG,desc\n# a comment\n\n   tttt  \n";
        let set = parse_onlist(text);
        let mut got: Vec<Vec<u8>> = set.into_iter().collect();
        got.sort();
        let mut want = vec![
            b"AAAA".to_vec(),
            b"CCCC".to_vec(),
            b"GGGG".to_vec(),
            b"TTTT".to_vec(),
        ];
        want.sort();
        assert_eq!(got, want);
    }

    #[test]
    fn onlist_overlap_matches_variable_length_barcodes() {
        use std::collections::{HashMap, HashSet};
        // inDrop bc1 is 8-11 bp: the onlist holds full-length barcodes, reads carry
        // whichever length, so a fixed 8 bp window would miss the 9-11 bp ones.
        let set: HashSet<Vec<u8>> = [
            b"AAAAAAAA".to_vec(),
            b"CCCCCCCCC".to_vec(),
            b"GGGGGGGGGG".to_vec(),
        ]
        .into_iter()
        .collect();
        let mut sets = HashMap::new();
        sets.insert("bc1.txt".to_string(), set);

        let bcs = [
            b"AAAAAAAA".to_vec(),
            b"CCCCCCCCC".to_vec(),
            b"GGGGGGGGGG".to_vec(),
        ];
        let obs: Vec<Vec<u8>> = (0..300)
            .map(|i| {
                let mut r = bcs[i % bcs.len()].clone();
                r.extend_from_slice(b"TTTTTTTTTTTT");
                r
            })
            .collect();
        let segs = vec![OnlistSeg {
            offset: 0,
            len: 8,
            max_len: 10,
            filename: "bc1.txt".into(),
        }];
        let (_, per_seg) = onlist_overlap_llr(&obs, &segs, &sets);
        assert!(
            per_seg[0].1 > 0.95,
            "variable-length barcodes should nearly all match, got {}",
            per_seg[0].1
        );
    }

    #[test]
    fn detect_assays_boosts_score_with_cached_onlist() {
        use std::collections::{HashMap, HashSet};
        let assays = embedded_assays();
        let indrop = assays
            .iter()
            .find(|a| a.spec_file == "indropv2.spec.yaml")
            .unwrap();
        let (_, layout) = indrop.structural_read().unwrap();

        // inDrop barcode read is R2; synth_read fills every barcode region with 'A'.
        let r2: Vec<Vec<u8>> = (0..500).map(|_| synth_read(&layout)).collect();
        let r1: Vec<Vec<u8>> = (0..500).map(|_| vec![b'C'; 100]).collect();

        let find = |hits: &[AssayHit]| {
            hits.iter()
                .find(|h| h.spec_file == "indropv2.spec.yaml")
                .cloned()
                .unwrap()
        };
        let without = find(&detect_assays(&r1, &r2, assays, &HashMap::new(), 60));

        // bc1 (the only framed barcode) is all-A in the synth reads; an onlist holding
        // that barcode boosts the score. bc2 is dropped (downstream of variable bc1).
        let a_set: HashSet<Vec<u8>> = [b"AAAAAAAA".to_vec()].into_iter().collect();
        let mut sets = HashMap::new();
        sets.insert("inDrop_barcode1_list.txt".to_string(), a_set);
        let with = find(&detect_assays(&r1, &r2, assays, &sets, 60));

        assert!(
            with.score > without.score + 100.0,
            "onlist overlap should raise the score: {} vs {}",
            with.score,
            without.score
        );
        assert_eq!(
            with.onlist_overlap.len(),
            1,
            "only the framed bc1 is scored"
        );
        assert!(with.onlist_overlap.iter().all(|(_, r)| *r > 0.9));
        assert!(without.onlist_overlap.is_empty());

        // A cached onlist of barcodes the reads carry in neither orientation (the synth
        // bc1 is AAAAAAAA, revcomp TTTTTTTT) must not raise the rank.
        let wrong: HashSet<Vec<u8>> = [b"GGGGGGGG".to_vec()].into_iter().collect();
        let mut wrong_sets = HashMap::new();
        wrong_sets.insert("inDrop_barcode1_list.txt".to_string(), wrong);
        let mismatch = find(&detect_assays(&r1, &r2, assays, &wrong_sets, 60));
        assert!(
            mismatch.score <= without.score,
            "off-list barcodes must not boost: {} vs {}",
            mismatch.score,
            without.score
        );
    }

    #[test]
    fn detect_assays_scores_all_framed_barcodes() {
        use std::collections::{HashMap, HashSet};
        // BD Rhapsody has three fixed 9 bp barcodes with no variable region among them,
        // so all three are framed and scored (unlike inDrop, where bc2 is dropped).
        let assays = embedded_assays();
        let bd = assays
            .iter()
            .find(|a| a.spec_file == "bd_rhapsody_v1.spec.yaml")
            .unwrap();
        let (_, layout) = bd.structural_read().unwrap();
        let r1: Vec<Vec<u8>> = (0..500).map(|_| synth_read(&layout)).collect();
        let r2: Vec<Vec<u8>> = (0..500).map(|_| vec![b'G'; 60]).collect();

        let a9: HashSet<Vec<u8>> = [b"AAAAAAAAA".to_vec()].into_iter().collect();
        let mut sets = HashMap::new();
        for f in ["BD_CLS1.txt", "BD_CLS2.txt", "BD_CLS3.txt"] {
            sets.insert(f.to_string(), a9.clone());
        }
        let hit = detect_assays(&r1, &r2, assays, &sets, 60)
            .into_iter()
            .find(|h| h.spec_file == "bd_rhapsody_v1.spec.yaml")
            .unwrap();
        assert_eq!(hit.onlist_overlap.len(), 3, "all three cls barcodes scored");
        assert!(hit.onlist_overlap.iter().all(|(_, r)| *r > 0.9));
    }

    #[test]
    fn onlist_sources_resolve_per_assay_and_are_unique() {
        // inDrop bc1 resolves; bc2 is intentionally absent (dropped by the framing rule).
        let bc1 = onlist_source("indropv2.spec.yaml", "bc1_onlist.txt").unwrap();
        assert!(bc1.url.starts_with("https://"));
        assert!(onlist_source("indropv2.spec.yaml", "bc2_onlist.txt").is_none());

        // The other wired assays resolve each barcode to a distinct hosted file.
        for (spec, local) in [
            ("bd_rhapsody_v1.spec.yaml", "cls1_onlist.txt"),
            ("bd_rhapsody_v1.spec.yaml", "cls2_onlist.txt"),
            ("bd_rhapsody_v1.spec.yaml", "cls3_onlist.txt"),
            ("pipseqv2.spec.yaml", "pip-seq_v2_bc1.tsv"),
            ("pipseqv2.spec.yaml", "pip-seq_v2_bc3.tsv"),
            ("pipseqv3.spec.yaml", "fb_v3_bc1.tsv"),
            ("pipseqv3.spec.yaml", "fb_v3_bc4.tsv"),
            ("quartz_seq2.spec.yaml", "cell_bc_onlist.txt"),
        ] {
            let src = onlist_source(spec, local)
                .unwrap_or_else(|| panic!("{spec}/{local} should resolve"));
            assert!(src.url.starts_with("https://"), "{spec}/{local}");
        }
        // The same spec-local onlist name (cell_bc_onlist.txt) is used by other specs,
        // so resolution must be keyed per assay, not by filename alone.
        assert!(onlist_source("drop_seq.spec.yaml", "cell_bc_onlist.txt").is_none());

        // Unframeable / non-barcode-list assays stay geometry-only.
        assert!(onlist_source("split_seq.spec.yaml", "barcode-1_onlist.txt").is_none());
        assert!(onlist_source("bd_rhapsody_eb.spec.yaml", "cls1_onlist.txt").is_none());
        assert!(onlist_source("indropv2.spec.yaml", "nope.txt").is_none());

        // Cache names are globally unique, so the loaded-set map never collides.
        let mut seen = std::collections::HashSet::new();
        for (_, _, src) in ONLIST_SOURCES {
            assert!(
                seen.insert(src.cache_name),
                "dup cache_name {}",
                src.cache_name
            );
        }
    }

    #[test]
    fn onlist_sources_reference_the_live_catalogue() {
        // Each entry keys on strings from the vendored seqspec YAML. A submodule bump
        // that renames a spec or onlist would otherwise silently degrade an assay to
        // geometry-only; this turns that into a red test.
        let assays = embedded_assays();
        for (spec, local, _) in ONLIST_SOURCES {
            let a = assays
                .iter()
                .find(|a| a.spec_file == *spec)
                .unwrap_or_else(|| panic!("{spec} not in embedded catalogue"));
            let found = a
                .leaf_regions()
                .iter()
                .any(|r| r.onlist.as_ref().map(|o| o.filename.as_str()) == Some(*local));
            assert!(found, "{spec} has no onlist named {local}");
        }
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
        let hits = detect_assays(&r1, &r2, assays, &Default::default(), 10);
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
        for a in assays {
            let Some((_rid, layout)) = a.structural_read() else {
                continue;
            };
            let expected: usize = layout.iter().map(|s| s.len).sum();
            if expected == 0 {
                continue;
            }
            let r1: Vec<Vec<u8>> = (0..300).map(|_| synth_read(&layout)).collect();
            let r2: Vec<Vec<u8>> = (0..300).map(|_| vec![b'G'; 90]).collect();
            let hits = detect_assays(&r1, &r2, assays, &Default::default(), assays.len());
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
