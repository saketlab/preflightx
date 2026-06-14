//! Sex-marker k-mer panels: parse the `build.rs`-embedded `data/sexmarkers/*.panel`
//! tables into a fast lookup used by [`crate::sex`].
//!
//! A panel is a tiny, public-genome-derived table of **canonical** k-mers (one strand of
//! each, `min(seq, revcomp)`), each tagged with a marker gene and a category
//! (`female` = the X-inactivation lncRNA XIST/Xist; `male` = a Y-linked gene). A k-mer is
//! in the panel only if it uniquely identifies its gene across the genome *and* the
//! transcriptome (see `experiments/sex_panel/`), so a read containing it can only have
//! come from that gene. Both orientations are stored, so a read's forward k-mer matches
//! regardless of cDNA strand.

use std::collections::HashMap;
use std::sync::OnceLock;

include!(concat!(env!("OUT_DIR"), "/sexpanels.rs")); // pub static EMBEDDED_PANELS

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Category {
    Female, // XIST / Xist
    Male,   // a Y-linked gene
}

/// One species' parsed panel: metadata header + canonical-k-mer lookup.
pub struct SexPanel {
    pub species: String,
    pub mechanism: String, // e.g. "XY_XIST"; selects the decision interpretation
    pub status: String,    // "validated" | "experimental"
    pub k: usize,
    pub genes: Vec<String>,
    /// k-mer (both orientations stored) -> (category, index into `genes`)
    kmers: HashMap<Vec<u8>, (Category, u16)>,
}

impl SexPanel {
    /// Look up a read k-mer (either orientation is present in the map).
    pub fn get(&self, kmer: &[u8]) -> Option<(Category, u16)> {
        self.kmers.get(kmer).copied()
    }
    pub fn len(&self) -> usize {
        self.kmers.len()
    }
    pub fn is_empty(&self) -> bool {
        self.kmers.is_empty()
    }

    /// One stored k-mer of the given category (for tests that synthesise reads).
    #[cfg(test)]
    pub(crate) fn sample(&self, cat: Category) -> Option<&[u8]> {
        self.kmers
            .iter()
            .find(|(_, (c, _))| *c == cat)
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

fn parse_panel(text: &str) -> Option<SexPanel> {
    let mut meta: HashMap<&str, &str> = HashMap::new();
    let mut genes: Vec<String> = Vec::new();
    let mut gene_idx: HashMap<String, u16> = HashMap::new();
    let mut kmers: HashMap<Vec<u8>, (Category, u16)> = HashMap::new();

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
        // Store both orientations so a read's forward k-mer slice matches directly,
        // without canonicalising (and allocating) per k-mer at scan time.
        let fwd = kmer.as_bytes().to_vec();
        let rc = revcomp(fwd.as_slice());
        kmers.insert(fwd, (category, idx));
        kmers.entry(rc).or_insert((category, idx));
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
    })
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

/// The panel for `species` (case-insensitive), if one is embedded.
pub fn for_species(species: &str) -> Option<&'static SexPanel> {
    registry()
        .iter()
        .find(|p| p.species.eq_ignore_ascii_case(species))
}

/// Names of the species with an embedded panel (for error messages / `--help`).
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
