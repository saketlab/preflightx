//! Sex-marker k-mer panels: the build.rs-embedded data/sexmarkers/*.panel tables of
//! canonical k-mers, parsed into the lookup used by [crate::sex].
//!
//! A k-mer is in a panel only if it uniquely identifies its marker gene across the genome
//! and the transcriptome (see experiments/sex_panel/), so a read carrying it can only have
//! come from that gene. Both orientations are stored, so a read's forward k-mer matches
//! regardless of cDNA strand.

use std::collections::HashMap;
use std::sync::OnceLock;

use rustc_hash::FxHashMap;

include!(concat!(env!("OUT_DIR"), "/sexpanels.rs")); // pub static EMBEDDED_PANELS

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Category {
    Female, // XIST / Xist
    Male,   // a Y-linked gene
}

/// One species' parsed panel: metadata header + canonical-k-mer lookup.
pub struct SexPanel {
    pub species: String,
    pub mechanism: String, // e.g. "XY_XIST"; selects the decision rule
    pub status: String,    // "validated" | "experimental"
    pub k: usize,
    pub genes: Vec<String>,
    /// 2-bit encoded k-mer (both orientations) -> (category, index into genes)
    kmers: FxHashMap<u128, (Category, u16)>,
    #[cfg(test)]
    sample_kmers: Vec<(Vec<u8>, Category)>,
}

impl SexPanel {
    /// Look up a read k-mer (either orientation is present in the map).
    pub fn get(&self, kmer: &[u8]) -> Option<(Category, u16)> {
        encode_kmer(kmer).and_then(|key| self.get_encoded(key))
    }

    /// Look up an already 2-bit encoded read k-mer.
    pub(crate) fn get_encoded(&self, key: u128) -> Option<(Category, u16)> {
        self.kmers.get(&key).copied()
    }
    pub fn len(&self) -> usize {
        self.kmers.len()
    }
    pub fn is_empty(&self) -> bool {
        self.kmers.is_empty()
    }
    /// Whether the panel carries Y-linked markers; XIST-only panels can't call males.
    pub fn has_y(&self) -> bool {
        self.kmers.values().any(|(c, _)| *c == Category::Male)
    }

    /// One stored k-mer of the given category (for tests that synthesise reads).
    #[cfg(test)]
    pub(crate) fn sample(&self, cat: Category) -> Option<&[u8]> {
        self.sample_kmers
            .iter()
            .find(|(_, c)| *c == cat)
            .map(|(k, _)| k.as_slice())
    }
}

/// Reverse-complement (uppercase ACGTN; other bytes complemented to N).
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

#[inline]
fn base_bits(b: u8) -> Option<u8> {
    match b {
        65 | 97 => Some(0),
        67 | 99 => Some(1),
        71 | 103 => Some(2),
        84 | 116 => Some(3),
        _ => None,
    }
}

/// 2-bit encode one A/C/G/T k-mer. Ambiguous bases and k > 64 return None.
pub(crate) fn encode_kmer(seq: &[u8]) -> Option<u128> {
    if seq.len() > 64 {
        return None;
    }
    let mut key = 0u128;
    for &b in seq {
        key = (key << 2) | base_bits(b)? as u128;
    }
    Some(key)
}

/// Encode a k-mer and its reverse complement in one pass, without allocating the revcomp:
/// the revcomp key is sum_i comp(base_i) << 2i, where comp = 3 - bits (A<->T, C<->G).
pub(crate) fn encode_kmer_both_strands(seq: &[u8]) -> Option<(u128, u128)> {
    if seq.len() > 64 {
        return None;
    }
    let mut fwd = 0u128;
    let mut rc = 0u128;
    for (i, &b) in seq.iter().enumerate() {
        let bits = base_bits(b)?;
        fwd = (fwd << 2) | bits as u128;
        rc |= ((3 - bits) as u128) << (2 * i);
    }
    Some((fwd, rc))
}

/// Slide a k-window over seq, calling hit with the 2-bit key of each A/C/G/T-only window.
/// Returns false only when the read is shorter than k; a long all-N read counts as scanned.
pub(crate) fn scan_encoded_windows(seq: &[u8], k: usize, mut hit: impl FnMut(u128)) -> bool {
    if seq.len() < k {
        return false;
    }
    if k == 0 || k > 64 {
        return true;
    }
    let mask = if k == 64 {
        u128::MAX
    } else {
        (1u128 << (2 * k)) - 1
    };
    let mut key = 0u128;
    let mut valid = 0usize;
    for &b in seq {
        if let Some(bits) = base_bits(b) {
            key = ((key << 2) | bits as u128) & mask;
            valid += 1;
            if valid >= k {
                hit(key);
            }
        } else {
            key = 0;
            valid = 0;
        }
    }
    true
}

fn parse_panel(text: &str) -> Option<SexPanel> {
    let mut meta: HashMap<&str, &str> = HashMap::new();
    let mut genes: Vec<String> = Vec::new();
    let mut gene_idx: HashMap<String, u16> = HashMap::new();
    let mut kmers: FxHashMap<u128, (Category, u16)> = FxHashMap::default();
    #[cfg(test)]
    let mut sample_kmers: Vec<(Vec<u8>, Category)> = Vec::new();

    for line in text.lines() {
        if line.is_empty() {
            continue;
        }
        if let Some(rest) = line.strip_prefix('#') {
            // metadata "key\tvalue"; ignore the "#columns" descriptor row
            if let Some((k, v)) = rest.split_once('\t') {
                meta.insert(k, v);
            }
            continue;
        }
        // data row: canonical_kmer \t gene \t category  (extra cols ignored)
        let mut f = line.split('\t');
        let (kmer, gene, cat) = match (f.next(), f.next(), f.next()) {
            (Some(a), Some(b), Some(c)) => (a, b, c),
            _ => continue,
        };
        let category = match cat {
            "female" => Category::Female,
            "male" => Category::Male,
            _ => continue,
        };
        let idx = *gene_idx.entry(gene.to_string()).or_insert_with(|| {
            genes.push(gene.to_string());
            (genes.len() - 1) as u16
        });
        // both orientations, so scan-time lookup needs no per-k-mer canonicalisation
        if let Some((fwd, rc)) = encode_kmer_both_strands(kmer.as_bytes()) {
            kmers.insert(fwd, (category, idx));
            kmers.entry(rc).or_insert((category, idx));
            #[cfg(test)]
            if !sample_kmers.iter().any(|(_, c)| *c == category) {
                sample_kmers.push((kmer.as_bytes().to_vec(), category));
            }
        }
    }

    if kmers.is_empty() {
        return None;
    }
    Some(SexPanel {
        species: meta.get("species").unwrap_or(&"unknown").to_string(),
        mechanism: meta.get("mechanism").unwrap_or(&"XY_XIST").to_string(),
        status: meta.get("status").unwrap_or(&"experimental").to_string(),
        k: meta.get("k").and_then(|s| s.parse().ok()).unwrap_or(50),
        genes,
        kmers,
        #[cfg(test)]
        sample_kmers,
    })
}

/// Parse the 3-pool .spanel format (k=50, canonical) into the female/male
/// [SexPanel] plus the list of canonical autosomal (species-ID) k-mers.
///
/// Row layout is robust to the column quirk in the catalogue: autosomal rows are
/// kmer\tautosomal (2 cols) while female/male rows are count\tkmer\tcategory
/// (3 cols). The **category is always the last field** and the **k-mer is the
/// longest field**, so we pick those two regardless of column count. .spanel
/// rows carry no gene name, so female k-mers are tagged XIST and male Y
/// (one synthetic gene each; the panels guarantee a single Y/XIST locus).
pub(crate) fn parse_spanel(text: &str) -> Option<(SexPanel, Vec<Vec<u8>>)> {
    let mut meta: HashMap<&str, &str> = HashMap::new();
    let mut genes: Vec<String> = Vec::new();
    let mut gene_idx: HashMap<&str, u16> = HashMap::new();
    let mut kmers: FxHashMap<u128, (Category, u16)> = FxHashMap::default();
    let mut autosomal: Vec<Vec<u8>> = Vec::new();
    #[cfg(test)]
    let mut sample_kmers: Vec<(Vec<u8>, Category)> = Vec::new();

    for line in text.lines() {
        if line.is_empty() {
            continue;
        }
        if let Some(rest) = line.strip_prefix('#') {
            if let Some((k, v)) = rest.split_once('\t') {
                meta.insert(k, v);
            }
            continue;
        }
        // category = last field; k-mer = longest field (handles 2- and 3-col rows).
        let mut cat = "";
        let mut kmer = "";
        for f in line.split('\t') {
            cat = f;
            if f.len() > kmer.len() {
                kmer = f;
            }
        }
        match cat {
            "autosomal" => autosomal.push(kmer.as_bytes().to_vec()),
            "female" | "male" => {
                let (category, gene) = if cat == "female" {
                    (Category::Female, "XIST")
                } else {
                    (Category::Male, "Y")
                };
                let idx = *gene_idx.entry(gene).or_insert_with(|| {
                    genes.push(gene.to_string());
                    (genes.len() - 1) as u16
                });
                if let Some((fwd, rc)) = encode_kmer_both_strands(kmer.as_bytes()) {
                    kmers.insert(fwd, (category, idx));
                    kmers.entry(rc).or_insert((category, idx));
                    #[cfg(test)]
                    if !sample_kmers.iter().any(|(_, c)| *c == category) {
                        sample_kmers.push((kmer.as_bytes().to_vec(), category));
                    }
                }
            }
            _ => continue,
        }
    }

    if autosomal.is_empty() && kmers.is_empty() {
        return None;
    }
    let panel = SexPanel {
        species: meta.get("species").unwrap_or(&"unknown").to_string(),
        // Header-declared if present; current catalogue is all XY/XIST, auto-derived.
        mechanism: meta.get("mechanism").unwrap_or(&"XY_XIST").to_string(),
        status: meta.get("status").unwrap_or(&"experimental").to_string(),
        k: meta.get("k").and_then(|s| s.parse().ok()).unwrap_or(50),
        genes,
        kmers,
        #[cfg(test)]
        sample_kmers,
    };
    Some((panel, autosomal))
}

fn registry() -> &'static Vec<SexPanel> {
    static R: OnceLock<Vec<SexPanel>> = OnceLock::new();
    R.get_or_init(|| {
        EMBEDDED_PANELS
            .iter()
            .filter_map(|(_, t)| parse_panel(t))
            .collect()
    })
}

/// The panel for species (case-insensitive), if one is embedded.
pub fn for_species(species: &str) -> Option<&'static SexPanel> {
    registry()
        .iter()
        .find(|p| p.species.eq_ignore_ascii_case(species))
}

/// Names of the species with an embedded panel (for error messages / --help).
pub fn species_list() -> Vec<&'static str> {
    registry().iter().map(|p| p.species.as_str()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn revcomp_basic() {
        assert_eq!(revcomp(b"AAAACCCGGT"), b"ACCGGGTTTT");
    }

    #[test]
    fn both_strands_matches_reference() {
        // The one-pass bit-trick must equal encode(fwd) + encode(revcomp).
        for kmer in [b"A".as_slice(), b"ACGT", b"AAAACCCGGT", b"GTCAGTCAGT"] {
            let (fwd, rc) = encode_kmer_both_strands(kmer).unwrap();
            assert_eq!(Some(fwd), encode_kmer(kmer), "fwd {kmer:?}");
            assert_eq!(Some(rc), encode_kmer(&revcomp(kmer)), "rc {kmer:?}");
        }
        assert_eq!(encode_kmer_both_strands(b"ACGTN"), None); // ambiguous base
    }

    #[test]
    fn parse_minimal_panel() {
        let txt = "#species\thuman\n#k\t6\n#columns\tkmer\tgene\tcategory\n\
                   AAAAAA\tXIST\tfemale\nCCCCCC\tDDX3Y\tmale\n";
        let p = parse_panel(txt).unwrap();
        assert_eq!(p.species, "human");
        assert_eq!(p.k, 6);
        // both orientations stored: AAAAAA+TTTTTT (XIST), CCCCCC+GGGGGG (DDX3Y)
        assert_eq!(p.len(), 4);
        let xist = p.genes.iter().position(|g| g == "XIST").unwrap() as u16;
        assert_eq!(p.get(b"AAAAAA"), Some((Category::Female, xist)));
        assert_eq!(p.get(b"TTTTTT"), Some((Category::Female, xist))); // revcomp
        assert_eq!(p.get(b"CCCCCC").unwrap().0, Category::Male);
        assert_eq!(p.get(b"GGGGGG").unwrap().0, Category::Male); // revcomp of CCCCCC
        assert_eq!(p.get(b"ACGTAC"), None);
    }

    #[test]
    fn embedded_panels_load_and_are_specific() {
        // The real embedded panels parse, are non-trivial, and carry both categories.
        let human = for_species("human").expect("human panel embedded");
        assert!(
            human.len() > 10_000,
            "human panel too small: {}",
            human.len()
        );
        assert_eq!(human.k, 50);
        let has_female = human.genes.iter().any(|g| g.eq_ignore_ascii_case("xist"));
        assert!(has_female, "human panel missing XIST");
        assert!(species_list().contains(&"human"));
    }
}
