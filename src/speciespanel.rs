use anyhow::Result;
use std::collections::HashMap;
use std::fs;

use crate::sexpanel::{self, encode_kmer, revcomp, SexPanel};
use crate::species::same_species;
use crate::whitelists;

pub struct Catalogue {
    panels: Vec<SexPanel>,
    autosomal: HashMap<u128, u16>,
}

impl Catalogue {
    pub fn load() -> Result<Catalogue> {
        let dir = whitelists::panel_cache_dir()?;
        let mut panels = Vec::new();
        let mut autosomal = HashMap::new();
        if dir.is_dir() {
            let mut files: Vec<_> = fs::read_dir(&dir)?
                .filter_map(|e| e.ok())
                .map(|e| e.path())
                .filter(|p| p.extension().and_then(|x| x.to_str()) == Some("spanel"))
                .collect();
            files.sort();
            for path in files {
                let text = fs::read_to_string(&path)?;
                if let Some((panel, auto)) = sexpanel::parse_spanel(&text) {
                    push_panel(&mut panels, &mut autosomal, panel, auto);
                }
            }
        }
        Ok(Catalogue { panels, autosomal })
    }

    pub fn from_texts(texts: &[&str]) -> Catalogue {
        let mut panels = Vec::new();
        let mut autosomal = HashMap::new();
        for text in texts {
            if let Some((panel, auto)) = sexpanel::parse_spanel(text) {
                push_panel(&mut panels, &mut autosomal, panel, auto);
            }
        }
        Catalogue { panels, autosomal }
    }

    pub fn is_empty(&self) -> bool {
        self.panels.is_empty()
    }
    pub fn len(&self) -> usize {
        self.panels.len()
    }
    pub fn k(&self) -> usize {
        self.panels.first().map(|p| p.k).unwrap_or(50)
    }
    pub fn species_of(&self, idx: u16) -> &str {
        &self.panels[idx as usize].species
    }
    pub fn autosomal_hit(&self, kmer: &[u8]) -> Option<u16> {
        encode_kmer(kmer).and_then(|key| self.autosomal_hit_encoded(key))
    }
    pub(crate) fn autosomal_hit_encoded(&self, key: u128) -> Option<u16> {
        self.autosomal.get(&key).copied()
    }
    pub fn for_species(&self, name: &str) -> Option<&SexPanel> {
        self.panels.iter().find(|p| same_species(&p.species, name))
    }
}

fn push_panel(
    panels: &mut Vec<SexPanel>,
    autosomal: &mut HashMap<u128, u16>,
    panel: SexPanel,
    auto: Vec<Vec<u8>>,
) {
    let idx = panels.len() as u16;
    for kmer in auto {
        let rc = revcomp(&kmer);
        if let Some(key) = encode_kmer(&kmer) {
            autosomal.entry(key).or_insert(idx);
        }
        if let Some(key) = encode_kmer(&rc) {
            autosomal.entry(key).or_insert(idx);
        }
    }
    panels.push(panel);
}

#[cfg(test)]
mod tests {
    use crate::sexpanel::{parse_spanel, Category};

    #[test]
    fn parse_spanel_splits_three_pools() {
        let txt = "#species\tHomo_sapiens\n#k\t6\n#n_autosomal\t2\n\
                   AAAAAA\tautosomal\n\
                   CCCCCC\tautosomal\n\
                   1\tACGTAA\tfemale\n\
                   1\tTTTGGG\tmale\n";
        let (panel, auto) = parse_spanel(txt).unwrap();
        assert_eq!(panel.species, "Homo_sapiens");
        assert_eq!(panel.k, 6);
        assert_eq!(auto.len(), 2); // two autosomal k-mers
                                   // female/male land in the SexPanel, both orientations (4 entries).
        assert_eq!(panel.len(), 4);
        assert_eq!(panel.get(b"ACGTAA").unwrap().0, Category::Female);
        assert_eq!(panel.get(b"TTACGT").unwrap().0, Category::Female); // revcomp
        assert_eq!(panel.get(b"CCCCCC"), None); // autosomal is NOT a sex k-mer
        assert_eq!(panel.get(b"TTTGGG").unwrap().0, Category::Male);
        assert_eq!(panel.get(b"CCCAAA").unwrap().0, Category::Male); // revcomp of TTTGGG
    }
}
