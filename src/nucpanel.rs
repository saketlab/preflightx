//! Single-cell vs single-nucleus k-mer panels: an alignment-free read of whether a
//! library sequenced whole cells or isolated nuclei.
//!
//! Two signatures, both read off R2: `intron_frac = n_intron / (n_intron + n_exon)`
//! (nuclei carry unspliced pre-mRNA, high -> nucleus) and `mito_pct = n_mito /
//! reads_scanned` (nuclei are stripped of cytoplasmic mito, high -> cell). Each k-mer is
//! anchored to one category (`intron`/`exon`/`mito`) of a ubiquitous gene (or MT),
//! deduped so a hit is unambiguous. Fetched by `fetch-panels` into
//! [`crate::whitelists::panel_cache_dir`]; no cache -> `unknown`.

use rustc_hash::FxHashMap;
use serde::Serialize;
use std::sync::OnceLock;

use crate::kmerpanel::KmerPanel;
use crate::sexpanel::scan_encoded_windows;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Region {
    Intron,
    Exon,
    Mito,
}

/// Map the `.nucpanel` category column to a [`Region`].
fn parse_region(cat: &str) -> Option<Region> {
    match cat {
        "intron" => Some(Region::Intron),
        "exon" => Some(Region::Exon),
        "mito" => Some(Region::Mito),
        _ => None,
    }
}

const MIN_PANEL_HITS: usize = 20; // floor on total panel hits; below -> wrong species / not cDNA
const MITO_STRONG_CELL: f64 = 12.0; // mito this high -> whole cell
const MITO_CELL: f64 = 2.5; // mito_pct >= this -> whole-cell IF intron also low
const MITO_NUC: f64 = 1.5; // mito_pct <= this -> nuclei (mito-depleted)
const INTRON_NUC: f64 = 0.55; // intron_frac >= this -> nucleus-compatible
const INTRON_CELL: f64 = 0.45; // intron_frac < this -> cell-compatible
const INTRON_MIN_HITS: usize = 30; // need this many intron+exon hits to trust intron_frac

/// Single-cell vs single-nucleus verdict for a sample of cDNA (R2) reads.
#[derive(Serialize, Debug)]
pub struct NucReport {
    /// "single-cell" | "single-nucleus" | "ambiguous" | "unknown".
    pub verdict: String,
    /// "high" | "medium" | "low".
    pub confidence: String,
    pub n_intron: usize,
    pub n_exon: usize,
    pub n_mito: usize,
    pub reads_scanned: usize,
    /// n_intron / (n_intron + n_exon); high -> nucleus. 0.0 when no gene-body hits.
    pub intron_frac: f64,
    /// 100 * n_mito / reads_scanned; high -> whole cell. 0.0 when nothing scanned.
    pub mito_pct: f64,
    pub notes: Vec<String>,
}

/// The `.nucpanel` lookup, loaded from the panel cache dir once.
fn loaded() -> &'static KmerPanel<Region> {
    static R: OnceLock<KmerPanel<Region>> = OnceLock::new();
    R.get_or_init(|| KmerPanel::load("nucpanel", parse_region))
}

/// Whether any nucleus panel is loaded (else `classify` degrades to "unknown").
pub fn have_panels() -> bool {
    !loaded().is_empty()
}

/// Classify single-cell vs single-nucleus from sampled cDNA (R2) `reads`. Pure: no I/O.
pub fn classify(reads: &[Vec<u8>]) -> NucReport {
    let p = loaded();
    classify_with(reads, &p.map, p.k)
}

fn classify_with(reads: &[Vec<u8>], panel: &FxHashMap<u128, Region>, k: usize) -> NucReport {
    let mut n_intron = 0usize;
    let mut n_exon = 0usize;
    let mut n_mito = 0usize;
    let mut scanned = 0usize;

    if !panel.is_empty() {
        for read in reads {
            let mut hit_i = false;
            let mut hit_e = false;
            let mut hit_m = false;
            let long_enough = scan_encoded_windows(read, k, |key| match panel.get(&key) {
                Some(Region::Intron) => hit_i = true,
                Some(Region::Exon) => hit_e = true,
                Some(Region::Mito) => hit_m = true,
                None => {}
            });
            if long_enough {
                scanned += 1;
            }
            n_intron += hit_i as usize;
            n_exon += hit_e as usize;
            n_mito += hit_m as usize;
        }
    }

    let gene_hits = n_intron + n_exon;
    let intron_frac = if gene_hits == 0 {
        0.0
    } else {
        n_intron as f64 / gene_hits as f64
    };
    let mito_pct = if scanned == 0 {
        0.0
    } else {
        100.0 * n_mito as f64 / scanned as f64
    };
    let total_panel_hits = gene_hits + n_mito;
    let (verdict, confidence) = decide(
        gene_hits,
        total_panel_hits,
        mito_pct,
        intron_frac,
        panel.is_empty(),
    );

    let mut notes = Vec::new();
    if panel.is_empty() {
        notes.push("no nucleus panel cached; single-cell vs single-nucleus not assessed".into());
    } else if total_panel_hits < MIN_PANEL_HITS {
        notes.push(format!(
            "only {total_panel_hits} panel hits (< {MIN_PANEL_HITS}); the human/mouse panel doesn't match these reads — wrong species, or not cDNA"
        ));
    } else {
        match verdict.as_str() {
            "single-nucleus" => notes.push(format!(
                "single-nucleus: mito {mito_pct:.1}% (nuclei are stripped of cytoplasmic mito) / intron_frac {intron_frac:.2} (nuclear pre-mRNA)."
            )),
            "single-cell" => notes.push(format!(
                "single-cell: mito {mito_pct:.1}% (cytoplasmic mitochondria) with low intron_frac {intron_frac:.2}."
            )),
            _ => notes.push(format!(
                "ambiguous: mito {mito_pct:.1}% and intron_frac {intron_frac:.2} disagree or sit in the gray zone (mito-contaminated nuclei read like cells) — not calling."
            )),
        }
    }

    NucReport {
        verdict,
        confidence: confidence.into(),
        n_intron,
        n_exon,
        n_mito,
        reads_scanned: scanned,
        intron_frac: crate::fastq::round2(intron_frac),
        mito_pct: crate::fastq::round2(mito_pct),
        notes,
    }
}

/// Combine the two signals, requiring them to AGREE: mito-depleted -> nuclei; very high
/// mito -> cell; high mito with low intron -> cell; high mito with non-low intron ->
/// ambiguous (conflict); mid mito -> intron breaks the tie or abstains. The panel-hit
/// floor is checked first so a non-human/mouse sample reads `unknown`, not a false nucleus.
fn decide(
    gene_hits: usize,
    total_panel_hits: usize,
    mito_pct: f64,
    intron_frac: f64,
    no_panel: bool,
) -> (String, &'static str) {
    if no_panel || total_panel_hits < MIN_PANEL_HITS {
        return ("unknown".into(), "low");
    }
    let intron_ok = gene_hits >= INTRON_MIN_HITS; // else intron_frac is unreliable
    if mito_pct >= MITO_STRONG_CELL {
        return ("single-cell".into(), "high");
    }
    if mito_pct <= MITO_NUC {
        let conf = if mito_pct <= 0.5 {
            "high"
        } else if mito_pct <= 1.0 {
            "medium"
        } else {
            "low"
        };
        return ("single-nucleus".into(), conf);
    }
    if mito_pct >= MITO_CELL {
        if intron_ok && intron_frac < INTRON_CELL {
            let conf = if mito_pct >= 5.0 && intron_frac < 0.35 {
                "high"
            } else {
                "medium"
            };
            return ("single-cell".into(), conf);
        }
        return ("ambiguous".into(), "low");
    }
    if intron_ok && intron_frac >= INTRON_NUC {
        return ("single-nucleus".into(), "medium");
    }
    if intron_ok && intron_frac < INTRON_CELL {
        return ("single-cell".into(), "low");
    }
    ("ambiguous".into(), "low")
}

impl NucReport {
    pub fn print(&self) {
        use colored::Colorize;
        println!("{}", "Single-cell vs single-nucleus".bold().cyan());
        let v = match self.verdict.as_str() {
            "single-cell" => self.verdict.green(),
            "single-nucleus" => self.verdict.bright_blue(),
            "ambiguous" => self.verdict.yellow(),
            _ => self.verdict.dimmed(),
        };
        let c = match self.confidence.as_str() {
            "high" => self.confidence.green(),
            "medium" => self.confidence.yellow(),
            _ => self.confidence.red(),
        };
        println!("  VERDICT: {}   confidence: {}", v.bold(), c);
        println!(
            "  intron_frac: {:.2} (intron={} exon={})   mito: {:.1}% (n={})   reads: {}",
            self.intron_frac,
            self.n_intron,
            self.n_exon,
            self.mito_pct,
            self.n_mito,
            self.reads_scanned
        );
        for n in &self.notes {
            println!("    {} {}", "note:".dimmed(), n);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sexpanel::encode_kmer;

    #[allow(clippy::type_complexity)]
    fn tiny_panel() -> (FxHashMap<u128, Region>, usize, Vec<u8>, Vec<u8>, Vec<u8>) {
        let intron = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC".to_vec();
        let exon = b"GTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGT".to_vec();
        let mito = b"TTGCTTGCTTGCTTGCTTGCTTGCTTGCTTGCTTGCTTGCTTGCTTGCTT".to_vec();
        let mut map = FxHashMap::default();
        map.insert(encode_kmer(&intron).unwrap(), Region::Intron);
        map.insert(encode_kmer(&exon).unwrap(), Region::Exon);
        map.insert(encode_kmer(&mito).unwrap(), Region::Mito);
        (map, 50, intron, exon, mito)
    }

    fn reads_of(kmer: &[u8], n: usize) -> Vec<Vec<u8>> {
        (0..n)
            .map(|_| {
                let mut r = kmer.to_vec();
                r.extend(vec![b'N'; 40]); // pad; N never matches
                r
            })
            .collect()
    }

    #[test]
    fn mito_depleted_is_single_nucleus() {
        let (panel, k, _intron, exon, _mito) = tiny_panel();
        let reads = reads_of(&exon, 300);
        let r = classify_with(&reads, &panel, k);
        assert_eq!(r.verdict, "single-nucleus", "report={r:?}");
        assert_eq!(r.confidence, "high"); // 0% mito
    }

    #[test]
    fn high_mito_high_intron_conflict_is_ambiguous() {
        let (panel, k, intron, exon, mito) = tiny_panel();
        let mut reads = reads_of(&intron, 200); // intron_frac = 200/300 = 0.67
        reads.extend(reads_of(&exon, 100));
        reads.extend(reads_of(&mito, 40)); // ~12% mito, but < MITO_STRONG_CELL
        let r = classify_with(&reads, &panel, k);
        assert_eq!(r.verdict, "ambiguous", "report={r:?}");
        assert!(r.mito_pct >= 2.5 && r.intron_frac >= 0.55);
    }

    #[test]
    fn high_mito_is_single_cell() {
        let (panel, k, _intron, exon, mito) = tiny_panel();
        // ~11% mito over 900 reads
        let mut reads = reads_of(&mito, 100);
        reads.extend(reads_of(&exon, 800));
        let r = classify_with(&reads, &panel, k);
        assert_eq!(r.verdict, "single-cell", "report={r:?}");
        assert_eq!(r.confidence, "high");
        assert!(r.mito_pct >= 5.0);
    }

    #[test]
    fn middling_both_is_ambiguous() {
        let (panel, k, intron, exon, mito) = tiny_panel();
        let mut reads = reads_of(&intron, 40); // intron_frac = 40/80 = 0.5
        reads.extend(reads_of(&exon, 40));
        reads.extend(reads_of(&mito, 20));
        reads.extend((0..900).map(|_| vec![b'N'; 90])); // pad so mito ~= 2%
        let r = classify_with(&reads, &panel, k);
        assert_eq!(
            r.verdict, "ambiguous",
            "mito={} intron={}",
            r.mito_pct, r.intron_frac
        );
    }

    #[test]
    fn too_few_panel_hits_is_unknown() {
        let (panel, k, _intron, _exon, _mito) = tiny_panel();
        let reads: Vec<Vec<u8>> = (0..5000).map(|_| vec![b'N'; 90]).collect();
        let r = classify_with(&reads, &panel, k);
        assert_eq!(r.verdict, "unknown", "report={r:?}");
    }

    #[test]
    fn empty_panel_is_unknown() {
        let empty: FxHashMap<u128, Region> = FxHashMap::default();
        let (_p, k, intron, _e, _m) = tiny_panel();
        let r = classify_with(&reads_of(&intron, 100), &empty, k);
        assert_eq!(r.verdict, "unknown");
    }

    #[test]
    fn committed_panel_parses_and_has_all_categories() {
        let text = std::fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/data/nucpanels/human.nucpanel"
        ))
        .expect("committed human.nucpanel present");
        let map = KmerPanel::parse(&text, parse_region).map;
        let ni = map.values().filter(|r| **r == Region::Intron).count();
        let ne = map.values().filter(|r| **r == Region::Exon).count();
        let nm = map.values().filter(|r| **r == Region::Mito).count();
        assert!(
            ni > 5000 && ne > 2000 && nm > 5000,
            "intron={ni} exon={ne} mito={nm}"
        );
    }

    #[test]
    #[ignore = "network: streams R2 from ENA"]
    fn real_snrna_reads_single_nucleus() {
        let reads = crate::fastq::sample(
            std::path::Path::new(
                "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR263/005/SRR26313205/SRR26313205_2.fastq.gz",
            ),
            20_000,
        )
        .unwrap()
        .0;
        let r = classify(&reads);
        assert_eq!(r.verdict, "single-nucleus", "{r:?}");
    }
}
