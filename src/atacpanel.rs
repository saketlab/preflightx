//! ATAC-vs-RNA k-mer panel: an alignment-free read of whether a library is chromatin
//! accessibility (ATAC) or transcriptomic (RNA).
//!
//! Each k-mer is anchored to one of two categories of a ubiquitous gene: `exon` (mature
//! transcript; RNA reads land here) or `promoter` (genomic sequence around the TSS, open
//! in ATAC but never transcribed into mRNA). `promoter_frac = n_promoter / (n_promoter +
//! n_exon)`: RNA cDNA reads are exon-only (frac ~0), ATAC genomic reads also hit accessible
//! promoters (frac up). Fetched by `fetch-panels` into
//! [`crate::whitelists::panel_cache_dir`]; no cache -> `unknown`. This is a *modality*
//! signal folded into chemistry detection, not a standalone verdict.
//! Built by `experiments/atac_panel/build_atac_panel.py`.

use rustc_hash::FxHashMap;
use std::collections::HashMap;
use std::sync::OnceLock;

use crate::kmerpanel::KmerPanel;
use crate::sexpanel::scan_encoded_windows;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Cat {
    Exon,
    Promoter,
}

/// Map the `.atacpanel` category column to a [`Cat`].
fn parse_cat(cat: &str) -> Option<Cat> {
    match cat {
        "exon" => Some(Cat::Exon),
        "promoter" => Some(Cat::Promoter),
        _ => None,
    }
}

const MIN_HITS: usize = 10; // floor on exon+promoter hits to call at all
const FRAC_ATAC: f64 = 0.30; // promoter_frac >= this -> ATAC
const FRAC_RNA: f64 = 0.10; // promoter_frac <= this -> RNA

/// The ATAC-vs-RNA modality signal for a sample of reads.
#[derive(Debug, Clone)]
pub struct AtacSignal {
    /// "atac" | "rna" | "unknown".
    pub modality: &'static str,
    pub n_exon: usize,
    pub n_promoter: usize,
    /// n_promoter / (n_promoter + n_exon); high -> ATAC. 0.0 when no hits.
    pub promoter_frac: f64,
}

impl AtacSignal {
    /// True when the panel decisively called ATAC, usable to qualify a chemistry call.
    pub fn is_atac(&self) -> bool {
        self.modality == "atac"
    }
}

/// The `.atacpanel` lookups, one per species (`#species` header), loaded once.
fn species_panels() -> &'static HashMap<String, KmerPanel<Cat>> {
    static R: OnceLock<HashMap<String, KmerPanel<Cat>>> = OnceLock::new();
    R.get_or_init(|| KmerPanel::load_by_species("atacpanel", parse_cat))
}

/// Whether any ATAC panel is loaded (else `classify` returns "unknown").
pub fn have_panels() -> bool {
    !species_panels().is_empty()
}

/// The panel matching `species` (canonical, case-folded). Falls back to the sole panel
/// when exactly one is cached (single-species deployment); otherwise `None`, so a
/// multi-species cache without a species hint abstains rather than merge organisms.
fn select<'a>(
    panels: &'a HashMap<String, KmerPanel<Cat>>,
    species: Option<&str>,
) -> Option<&'a KmerPanel<Cat>> {
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

fn pick_panel(species: Option<&str>) -> Option<&'static KmerPanel<Cat>> {
    select(species_panels(), species)
}

fn unknown() -> AtacSignal {
    AtacSignal {
        modality: "unknown",
        n_exon: 0,
        n_promoter: 0,
        promoter_frac: 0.0,
    }
}

/// Classify ATAC vs RNA from sampled reads, scoring against the panel for `species`.
/// Pure: no I/O.
pub fn classify(reads: &[Vec<u8>], species: Option<&str>) -> AtacSignal {
    match pick_panel(species) {
        Some(p) => classify_with(reads, &p.map, p.k),
        None => unknown(),
    }
}

/// Classify over both mates. ATAC's two reads are both genomic, so scoring the pair
/// doubles the sparse hit count (RNA's short R1 contributes ~nothing). Chemistry's seam.
pub fn classify_pair(r1: &[Vec<u8>], r2: &[Vec<u8>], species: Option<&str>) -> AtacSignal {
    match pick_panel(species) {
        Some(p) => classify_with(r1.iter().chain(r2), &p.map, p.k),
        None => unknown(),
    }
}

/// Testable core: scan `reads` against an explicit panel map + k.
fn classify_with<'a>(
    reads: impl IntoIterator<Item = &'a Vec<u8>>,
    panel: &FxHashMap<u128, Cat>,
    k: usize,
) -> AtacSignal {
    let mut n_exon = 0usize;
    let mut n_promoter = 0usize;
    if !panel.is_empty() {
        for read in reads {
            let mut hit_e = false;
            let mut hit_p = false;
            scan_encoded_windows(read, k, |key| match panel.get(&key) {
                Some(Cat::Exon) => hit_e = true,
                Some(Cat::Promoter) => hit_p = true,
                None => {}
            });
            n_exon += hit_e as usize;
            n_promoter += hit_p as usize;
        }
    }

    let hits = n_exon + n_promoter;
    let promoter_frac = if hits == 0 {
        0.0
    } else {
        n_promoter as f64 / hits as f64
    };
    let modality = if hits < MIN_HITS {
        "unknown"
    } else if promoter_frac >= FRAC_ATAC {
        "atac"
    } else if promoter_frac <= FRAC_RNA {
        "rna"
    } else {
        "unknown"
    };

    AtacSignal {
        modality,
        n_exon,
        n_promoter,
        promoter_frac,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sexpanel::encode_kmer;

    fn tiny_panel() -> (FxHashMap<u128, Cat>, usize, Vec<u8>, Vec<u8>) {
        let exon = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC".to_vec();
        let prom = b"GTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGT".to_vec();
        let mut map = FxHashMap::default();
        map.insert(encode_kmer(&exon).unwrap(), Cat::Exon);
        map.insert(encode_kmer(&prom).unwrap(), Cat::Promoter);
        (map, 50, exon, prom)
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
    fn exon_only_is_rna() {
        let (panel, k, exon, _prom) = tiny_panel();
        let r = classify_with(&reads_from(&exon, 60), &panel, k);
        assert_eq!(r.modality, "rna", "{r:?}");
        assert_eq!(r.n_promoter, 0);
    }

    #[test]
    fn promoter_enriched_is_atac() {
        let (panel, k, exon, prom) = tiny_panel();
        let mut reads = reads_from(&prom, 40);
        reads.extend(reads_from(&exon, 60));
        let r = classify_with(&reads, &panel, k);
        assert_eq!(r.modality, "atac", "{r:?}");
        assert!(r.is_atac());
    }

    #[test]
    fn select_routes_by_species_and_falls_back() {
        let (map, k, _, _) = tiny_panel();
        let mut panels = HashMap::new();
        panels.insert("homo_sapiens".to_string(), KmerPanel { map, k });
        assert!(select(&panels, Some("human")).is_some());
        assert!(select(&panels, Some("mus_musculus")).is_some());
        assert!(select(&panels, None).is_some());

        let (map2, k2, _, _) = tiny_panel();
        panels.insert("mus_musculus".to_string(), KmerPanel { map: map2, k: k2 });
        assert!(select(&panels, None).is_none());
        assert!(select(&panels, Some("danio_rerio")).is_none());
        assert!(select(&panels, Some("Mus musculus")).is_some());
    }

    #[test]
    fn too_few_hits_is_unknown() {
        let (panel, k, exon, _prom) = tiny_panel();
        let r = classify_with(&reads_from(&exon, 5), &panel, k); // below MIN_HITS
        assert_eq!(r.modality, "unknown");
    }

    fn ena(u: &str, n: usize) -> Vec<Vec<u8>> {
        crate::fastq::sample(std::path::Path::new(u), n).unwrap().0
    }
    fn ena_url(srr: &str, tail: &str, mate: u8) -> String {
        format!("https://ftp.sra.ebi.ac.uk/vol1/fastq/{}/{tail}/{srr}/{srr}_{mate}.fastq.gz", &srr[..6])
    }

    const ATAC_R1: &str =
        "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR263/034/SRR26313334/SRR26313334_1.fastq.gz";
    const ATAC_R2: &str =
        "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR263/034/SRR26313334/SRR26313334_2.fastq.gz";
    const RNA_R2: &str =
        "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR263/005/SRR26313205/SRR26313205_2.fastq.gz";

    #[test]
    #[ignore = "network: streams from ENA"]
    fn real_snatac_reads_atac() {
        let r = classify_pair(&ena(ATAC_R1, 60_000), &ena(ATAC_R2, 60_000), Some("Homo_sapiens"));
        assert_eq!(r.modality, "atac", "{r:?}");
    }

    #[test]
    #[ignore = "network: streams from ENA"]
    fn real_snrna_reads_rna() {
        let r = classify(&ena(RNA_R2, 60_000), Some("Homo_sapiens"));
        assert_eq!(r.modality, "rna", "{r:?}");
    }

    const RNA_3P_R2: &str =
        "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/057/SRR11442257/SRR11442257_2.fastq.gz";

    #[test]
    #[ignore = "network: streams from ENA"]
    fn real_3prime_rna_reads_rna() {
        let r = classify(&ena(RNA_3P_R2, 30_000), Some("Homo_sapiens"));
        assert_eq!(r.modality, "rna", "{r:?}");
    }

    fn assert_snatac(srr: &str, tail: &str, species: &str) {
        let r1 = ena(&ena_url(srr, tail, 1), 40_000);
        let r2 = ena(&ena_url(srr, tail, 2), 40_000);
        let r = classify_pair(&r1, &r2, Some(species));
        assert_eq!(r.modality, "atac", "{species} {srr}: {r:?}");
    }

    #[test]
    #[ignore = "network: streams from ENA"]
    fn real_cat_snatac_is_atac() {
        assert_snatac("SRR34653007", "007", "Felis_catus");
    }

    #[test]
    #[ignore = "network: streams from ENA"]
    fn real_cow_snatac_is_atac() {
        assert_snatac("SRR34653001", "001", "Bos_taurus");
    }

    #[test]
    #[ignore = "network: streams from ENA"]
    fn real_horse_snatac_is_atac() {
        assert_snatac("SRR34652998", "098", "Equus_caballus");
    }

    #[test]
    #[ignore = "network: streams from ENA"]
    fn real_rabbit_snatac_is_atac() {
        assert_snatac("SRR34652995", "095", "Oryctolagus_cuniculus");
    }

    #[test]
    #[ignore = "network: streams from ENA"]
    fn real_cat_snrna_not_atac() {
        let r1 = ena(&ena_url("SRR34653040", "040", 1), 40_000);
        let r2 = ena(&ena_url("SRR34653040", "040", 2), 40_000);
        let r = classify_pair(&r1, &r2, Some("Felis_catus"));
        assert_ne!(r.modality, "atac", "{r:?}");
    }
}
