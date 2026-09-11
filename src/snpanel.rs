//! Single-nucleus species-ID intron panel: the complement to the exon .spanel catalogue.
//!
//! snRNA reads are intron-dominated (unspliced pre-mRNA) and miss the exon-biased species
//! panels. Each .snpanel holds species-unique intronic k-mers. Loaded from .snpanel files
//! in the panel cache dir (fetched by fetch-panels) and scored by [crate::species] as a
//! second channel when the exon call is weak.

use rustc_hash::FxHashMap;

use crate::sexpanel::encode_kmer_both_strands;

/// Intron k-mers of every cached species, keyed to a species index (both strand orientations).
pub struct IntronCatalogue {
    species: Vec<String>,
    kmers: FxHashMap<u128, u16>,
    k: usize,
}

impl IntronCatalogue {
    #[cfg(not(target_arch = "wasm32"))]
    pub fn load() -> anyhow::Result<IntronCatalogue> {
        let dir = crate::whitelists::panel_cache_dir()?;
        let mut texts = Vec::new();
        if dir.is_dir() {
            let mut files: Vec<_> = std::fs::read_dir(&dir)?
                .filter_map(|e| e.ok())
                .map(|e| e.path())
                .filter(|p| p.extension().and_then(|x| x.to_str()) == Some("snpanel"))
                .collect();
            files.sort();
            for path in files {
                if let Ok(t) = std::fs::read_to_string(&path) {
                    texts.push(t);
                }
            }
        }
        Ok(IntronCatalogue::from_texts(
            &texts.iter().map(String::as_str).collect::<Vec<_>>(),
        ))
    }

    /// [load], or an empty catalogue when no .snpanel is cached / on error.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn load_or_empty() -> IntronCatalogue {
        Self::load().unwrap_or_else(|_| Self::from_texts(&[]))
    }

    pub fn from_texts(texts: &[&str]) -> IntronCatalogue {
        let mut species = Vec::new();
        let mut kmers = FxHashMap::default();
        let mut k = 50;
        for text in texts {
            let Some(name) = header(text, "species") else {
                continue;
            };
            if let Some(hk) = header(text, "k").and_then(|s| s.parse().ok()) {
                k = hk;
            }
            let idx = species.len() as u16;
            let mut any = false;
            for line in text.lines() {
                if line.is_empty() || line.starts_with('#') {
                    continue;
                }
                let kmer = line.split('\t').next().unwrap_or("");
                if let Some((fwd, rc)) = encode_kmer_both_strands(kmer.as_bytes()) {
                    kmers.entry(fwd).or_insert(idx);
                    kmers.entry(rc).or_insert(idx);
                    any = true;
                }
            }
            if any {
                species.push(name);
            }
        }
        IntronCatalogue { species, kmers, k }
    }

    pub fn is_empty(&self) -> bool {
        self.species.is_empty()
    }
    pub fn len(&self) -> usize {
        self.species.len()
    }
    pub fn k(&self) -> usize {
        self.k
    }
    pub fn species_of(&self, idx: u16) -> &str {
        &self.species[idx as usize]
    }
    pub fn hit_encoded(&self, key: u128) -> Option<u16> {
        self.kmers.get(&key).copied()
    }
    #[cfg(test)]
    pub fn hit(&self, kmer: &[u8]) -> Option<u16> {
        crate::sexpanel::encode_kmer(kmer).and_then(|k| self.hit_encoded(k))
    }
}

/// Value of a #key\t<value> header line.
fn header(text: &str, key: &str) -> Option<String> {
    let pfx = format!("#{key}\t");
    text.lines()
        .find_map(|l| l.strip_prefix(&pfx).map(|r| r.trim().to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loads_and_indexes_by_species() {
        let cat = IntronCatalogue::from_texts(&[
            "#species\tFelis_catus\n#k\t10\nAAAAAAAAAA\tACTB\tdeep_intron\n",
            "#species\tBos_taurus\n#k\t10\nGGGAGGGAGG\tGAPDH\tdeep_intron\n",
        ]);
        assert_eq!(cat.len(), 2);
        assert_eq!(cat.k(), 10);
        assert_eq!(
            cat.species_of(cat.hit(b"AAAAAAAAAA").unwrap()),
            "Felis_catus"
        );
        assert_eq!(
            cat.species_of(cat.hit(b"TTTTTTTTTT").unwrap()),
            "Felis_catus"
        ); // revcomp
        assert_eq!(
            cat.species_of(cat.hit(b"GGGAGGGAGG").unwrap()),
            "Bos_taurus"
        );
        assert!(cat.hit(b"ACACACACAC").is_none());
    }
}
