//! Transcript-end (3' vs 5') k-mer panels: an alignment-free discriminator between
//! 10x 3' and 5' chemistry (which share 16+10 R1 geometry) plus a full-length signature.
//!
//! Each panel k-mer is anchored to exactly one end of one ubiquitously-expressed
//! transcript, deduped so a hit is unambiguously the 5' or 3' end of a gene. Stream R2,
//! count `n5` vs `n3`, `r = n3/(n3+n5)`: high -> 3', low -> 5', roughly even ->
//! full-length. Panels are per-species (`#species` header) and selected by the detected
//! organism; a species with no cached panel abstains rather than score against another
//! organism's k-mers. Fetched by `fetch-panels` into
//! [`crate::whitelists::panel_cache_dir`]; no cache -> "unknown".
//!
//! Calibrated on whole-cell scRNA; single-nucleus reads spread into introns/gene bodies,
//! so the ratio is less reliable there (the geometry+whitelist `assay_id` is authoritative).

use rustc_hash::FxHashMap;
use serde::Serialize;
use std::collections::HashMap;
use std::sync::OnceLock;

use crate::kmerpanel::KmerPanel;
use crate::sexpanel::scan_encoded_windows;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum End {
    Five,  // 5' end of a mature transcript (cap-proximal)
    Three, // 3' end (polyA-proximal)
}

/// Map the `.endpanel` category column to an [`End`].
fn parse_end(cat: &str) -> Option<End> {
    match cat {
        "5p" => Some(End::Five),
        "3p" => Some(End::Three),
        _ => None,
    }
}

// Calibrated on 3'v2 (SRP254408, r=0.67) and 5'v3 (SRP329309, r=0.03).
const MIN_END_HITS: usize = 10; // floor on n3+n5 to call at all
const R_THREE: f64 = 0.60; // r >= this -> 3'
const R_FIVE: f64 = 0.40; // r <= this -> 5'
const FULL_LO: f64 = 0.45; // FULL_LO <= r <= FULL_HI -> full-length / plate-based
const FULL_HI: f64 = 0.55;

/// The transcript-end verdict for a sample of cDNA reads.
#[derive(Serialize, Debug)]
pub struct EndReport {
    /// "3'" | "5'" | "full-length" | "ambiguous" | "unknown".
    pub end: String,
    /// "high" | "medium" | "low".
    pub confidence: String,
    pub n5: usize,
    pub n3: usize,
    /// Reads scanned against the panel (length >= k).
    pub reads_scanned: usize,
    /// r = n3 / (n3 + n5); NaN-safe (0.0 when no end hits).
    pub ratio: f64,
    pub notes: Vec<String>,
}

impl EndReport {
    /// True when the panel gave a decisive polarity (3' or 5') — usable to correct a
    /// geometry-only chemistry call.
    pub fn is_decisive(&self) -> bool {
        (self.end == "3'" || self.end == "5'") && self.confidence != "low"
    }
}

/// The `.endpanel` lookups, one per species (`#species` header), loaded once.
fn species_panels() -> &'static HashMap<String, KmerPanel<End>> {
    static R: OnceLock<HashMap<String, KmerPanel<End>>> = OnceLock::new();
    R.get_or_init(|| KmerPanel::load_by_species("endpanel", parse_end))
}

/// Whether any end panel is cached (else `classify` returns a clean "unknown").
pub fn have_panels() -> bool {
    !species_panels().is_empty()
}

/// The panel matching `species` (canonical, case-folded). Falls back to the sole panel
/// when exactly one is cached; otherwise `None`, so a multi-species cache without a
/// matching hint abstains rather than score reads against another organism's k-mers.
fn pick_panel(species: Option<&str>) -> Option<&'static KmerPanel<End>> {
    let panels = species_panels();
    if let Some(sp) = species {
        if let Some(p) = panels.get(&crate::species::canonical_species(sp)) {
            return Some(p);
        }
    }
    if panels.len() == 1 {
        return panels.values().next();
    }
    None
}

/// Classify transcript-end bias from sampled cDNA (R2) `reads`, scoring against the panel
/// for `species`. Pure: no I/O.
pub fn classify(reads: &[Vec<u8>], species: Option<&str>) -> EndReport {
    match pick_panel(species) {
        Some(p) => classify_with(reads, &p.map, p.k),
        None => classify_with(reads, &FxHashMap::default(), crate::kmerpanel::DEFAULT_K),
    }
}

/// Testable core: classify against an explicit panel map + k.
fn classify_with(reads: &[Vec<u8>], panel: &FxHashMap<u128, End>, k: usize) -> EndReport {
    let mut n5 = 0usize;
    let mut n3 = 0usize;
    let mut scanned = 0usize;

    if !panel.is_empty() {
        for read in reads {
            let mut hit5 = false;
            let mut hit3 = false;
            let long_enough = scan_encoded_windows(read, k, |key| {
                if let Some(end) = panel.get(&key) {
                    match end {
                        End::Five => hit5 = true,
                        End::Three => hit3 = true,
                    }
                }
            });
            if long_enough {
                scanned += 1;
            }
            n5 += hit5 as usize;
            n3 += hit3 as usize;
        }
    }

    let total = n5 + n3;
    let ratio = if total == 0 {
        0.0
    } else {
        n3 as f64 / total as f64
    };
    let (end, confidence) = decide(n5, n3, panel.is_empty());

    let mut notes = Vec::new();
    if panel.is_empty() {
        notes.push("no transcript-end panel embedded; 3'/5' not assessed".into());
    } else if total < MIN_END_HITS {
        notes.push(format!(
            "only {total} transcript-end marker hits (< {MIN_END_HITS}); too few to place the reads on the transcript — sample more reads or the input isn't cDNA"
        ));
    } else if end == "full-length" {
        notes.push(
            "reads spread evenly across transcript ends — full-length / plate-based (Smart-seq), not a 3'/5'-primed droplet library".into(),
        );
    } else if end == "ambiguous" {
        notes.push(format!(
            "transcript-end ratio r={ratio:.2} is in the gray zone; not calling 3' or 5'"
        ));
    }

    EndReport {
        end,
        confidence: confidence.into(),
        n5,
        n3,
        reads_scanned: scanned,
        ratio: crate::fastq::round2(ratio),
        notes,
    }
}

/// The r = n3/(n3+n5) decision rule; gray zones return `ambiguous`, too-few-hits `unknown`.
fn decide(n5: usize, n3: usize, no_panel: bool) -> (String, &'static str) {
    if no_panel {
        return ("unknown".into(), "low");
    }
    let total = n5 + n3;
    if total < MIN_END_HITS {
        return ("unknown".into(), "low");
    }
    let r = n3 as f64 / total as f64;
    let strong = total >= 40;
    if r >= R_THREE {
        (
            "3'".into(),
            if strong && r >= 0.85 {
                "high"
            } else {
                "medium"
            },
        )
    } else if r <= R_FIVE {
        (
            "5'".into(),
            if strong && r <= 0.15 {
                "high"
            } else {
                "medium"
            },
        )
    } else if (FULL_LO..=FULL_HI).contains(&r) {
        ("full-length".into(), if strong { "medium" } else { "low" })
    } else {
        ("ambiguous".into(), "low")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sexpanel::encode_kmer;

    fn tiny_panel() -> (FxHashMap<u128, End>, usize, Vec<u8>, Vec<u8>) {
        let five = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC".to_vec();
        let three = b"GTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGT".to_vec();
        assert_eq!(five.len(), 50);
        assert_eq!(three.len(), 50);
        let mut map = FxHashMap::default();
        map.insert(encode_kmer(&five).unwrap(), End::Five);
        map.insert(encode_kmer(&three).unwrap(), End::Three);
        (map, 50, five, three)
    }

    fn reads_from(kmer: &[u8], n: usize) -> Vec<Vec<u8>> {
        (0..n)
            .map(|_| {
                let mut r = kmer.to_vec();
                r.extend(vec![b'A'; 40]);
                r
            })
            .collect()
    }

    #[test]
    fn three_prime_bias_calls_three() {
        let (panel, k, _five, three) = tiny_panel();
        let reads = reads_from(&three, 60);
        let r = classify_with(&reads, &panel, k);
        assert_eq!(r.end, "3'");
        assert_eq!(r.n3, 60);
        assert_eq!(r.n5, 0);
        assert_eq!(r.confidence, "high");
    }

    #[test]
    fn five_prime_bias_calls_five() {
        let (panel, k, five, _three) = tiny_panel();
        let reads = reads_from(&five, 60);
        let r = classify_with(&reads, &panel, k);
        assert_eq!(r.end, "5'");
        assert_eq!(r.confidence, "high");
    }

    #[test]
    fn even_spread_is_full_length() {
        let (panel, k, five, three) = tiny_panel();
        let mut reads = reads_from(&five, 30);
        reads.extend(reads_from(&three, 30));
        let r = classify_with(&reads, &panel, k);
        assert_eq!(r.end, "full-length", "r={} report={:?}", r.ratio, r);
    }

    #[test]
    fn too_few_hits_is_unknown() {
        let (panel, k, _five, three) = tiny_panel();
        let reads = reads_from(&three, 3); // below MIN_END_HITS
        let r = classify_with(&reads, &panel, k);
        assert_eq!(r.end, "unknown");
        assert_eq!(r.confidence, "low");
    }

    #[test]
    fn gray_zone_is_ambiguous() {
        let (panel, k, five, three) = tiny_panel();
        // r = 40/(40+30) ~= 0.57, between FULL_HI (0.55) and R_THREE (0.60) -> ambiguous
        let mut reads = reads_from(&three, 40);
        reads.extend(reads_from(&five, 30));
        let r = classify_with(&reads, &panel, k);
        assert_eq!(r.end, "ambiguous", "r={}", r.ratio);
    }

    #[test]
    fn committed_panel_parses_and_is_balanced() {
        let text = std::fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/data/endpanels/human.endpanel"
        ))
        .expect("committed human.endpanel present");
        let panel = KmerPanel::parse(&text, parse_end);
        let n5 = panel.map.values().filter(|e| **e == End::Five).count();
        let n3 = panel.map.values().filter(|e| **e == End::Three).count();
        // both orientations stored, so counts are ~2x the row count per end.
        assert!(n5 > 2000 && n3 > 2000, "5p={n5} 3p={n3}");

        let three: Vec<u8> = text
            .lines()
            .find(|l| l.ends_with("\t3p"))
            .map(|l| l.split('\t').next().unwrap().as_bytes().to_vec())
            .unwrap();
        let reads = reads_from(&three, 60);
        let r = classify_with(&reads, &panel.map, panel.k);
        assert_eq!(r.end, "3'", "report={r:?}");
    }

    const HUMAN_3P_R2: &str =
        "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/057/SRR11442257/SRR11442257_2.fastq.gz";
    const HUMAN_5P_R2: &str =
        "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR152/025/SRR15209425/SRR15209425_2.fastq.gz";

    // 20k reads dodges ENA mid-stream drops.
    fn end_of(url: &str) -> EndReport {
        let reads = crate::fastq::sample(std::path::Path::new(url), 20_000)
            .expect("stream R2 from ENA")
            .0;
        classify(&reads, Some("Homo_sapiens"))
    }

    #[test]
    #[ignore = "network: streams R2 from ENA"]
    fn real_3prime_srp254408_reads_three() {
        let r = end_of(HUMAN_3P_R2);
        assert_eq!(r.end, "3'", "ratio={} report={:?}", r.ratio, r);
    }

    #[test]
    #[ignore = "network: streams R2 from ENA"]
    fn real_5prime_srp329309_reads_five() {
        let r = end_of(HUMAN_5P_R2);
        assert_eq!(r.end, "5'", "ratio={} report={:?}", r.ratio, r);
    }
}
