//! seqspec ingestion: parse a seqspec YAML and resolve each read into ordered
//! segments (barcode / UMI / cDNA / linker). The generic substrate the detector
//! scores against, so any seqspec-describable assay needs no hand-coded AssaySpec.
//! See docs/seqspec-support-plan.md.

use anyhow::{Context, Result};
use serde::Deserialize;

/// A barcode/onlist reference inside a region; the list itself is fetched on demand, never bundled.
#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct Onlist {
    #[serde(default)]
    pub filename: String,
    #[serde(default)]
    pub md5: String,
    #[serde(default)]
    pub url: String,
}

/// One node in the seqspec library_spec region tree.
#[derive(Debug, Clone, Deserialize)]
pub struct Region {
    pub region_id: String,
    pub region_type: String,
    #[serde(default)]
    pub sequence_type: String,
    #[serde(default)]
    pub sequence: String,
    #[serde(default)]
    pub min_len: usize,
    #[serde(default)]
    pub max_len: usize,
    #[serde(default)]
    pub onlist: Option<Onlist>,
    #[serde(default)]
    pub regions: Vec<Region>,
}

impl Region {
    fn is_leaf(&self) -> bool {
        self.regions.is_empty()
    }
}

/// Cell/well barcode synonyms; seqspec region_type is free text.
/// Sample/lane indices (index7, i7, sample_index) are excluded: they identify the library, not the cell.
fn is_barcode_type(region_type: &str) -> bool {
    matches!(
        region_type.to_ascii_lowercase().as_str(),
        "barcode"
            | "cell_bc"
            | "cell_barcode"
            | "cell_index"
            | "bead_bc"
            | "well_bc"
            | "rt_bc"
            | "plate_bc"
            | "sample_bc"
            | "subarray_bc"
            | "hairpin_barcode"
            | "bc1"
            | "bc2"
            | "bc3"
    )
}

fn is_umi_type(region_type: &str) -> bool {
    region_type.eq_ignore_ascii_case("umi")
}

/// A barcode or UMI region; the first one marks the start of the observed read.
fn is_barcode_or_umi(region_type: &str) -> bool {
    is_barcode_type(region_type) || is_umi_type(region_type)
}

/// Inter-barcode sequences and short spacers; sequenced as part of the barcode read, so they do not terminate it.
fn is_linker_type(region_type: &str) -> bool {
    let rt = region_type.to_ascii_lowercase();
    rt.starts_with("link") || rt.starts_with("spacer")
}

/// The biological insert the barcode read runs into; matched by type, since some specs give it a fixed length.
fn is_insert_type(region_type: &str) -> bool {
    let rt = region_type.to_ascii_lowercase();
    matches!(
        rt.as_str(),
        "cdna" | "gdna" | "rna" | "dna" | "insert" | "hic" | "methyl" | "atac" | "protein"
    ) || rt.starts_with("poly") // ploy_t: typo'd upstream in some specs
        || rt.starts_with("ploy")
}

/// A region that absorbs the rest of the read: insert-type, or variable and not a barcode/UMI/linker.
/// Terminating on "any variable region" would truncate variable-length barcodes and spacers.
fn is_read_terminator(s: &Segment) -> bool {
    is_insert_type(&s.region_type)
        || (s.max_len > s.len
            && !is_barcode_type(&s.region_type)
            && !is_umi_type(&s.region_type)
            && !is_linker_type(&s.region_type))
}

/// One sequencing read in the seqspec sequence_spec.
#[derive(Debug, Clone, Deserialize)]
pub struct Read {
    pub read_id: String,
    #[serde(default)]
    pub primer_id: String,
    #[serde(default)]
    pub modality: String,
    #[serde(default)]
    pub max_len: usize,
    #[serde(default)]
    pub strand: String,
}

/// A parsed seqspec assay.
#[derive(Debug, Clone, Deserialize)]
pub struct Assay {
    #[serde(default)]
    pub assay_id: String,
    /// Source spec filename; the stable catalogue key, since assay_id is frequently empty upstream.
    #[serde(default)]
    pub spec_file: String,
    #[serde(default)]
    pub name: String,
    /// Publication DOI, often empty upstream; used by the doc generator.
    #[serde(default)]
    #[allow(dead_code)]
    pub doi: String,
    /// Free-text publication date; formats vary upstream and it is often empty.
    #[serde(default)]
    #[allow(dead_code)]
    pub date: String,
    #[serde(default)]
    pub modalities: Vec<String>,
    #[serde(default)]
    pub lib_struct: String,
    #[serde(default)]
    pub sequence_spec: Vec<Read>,
    #[serde(default)]
    pub library_spec: Vec<Region>,
}

/// One resolved segment a read covers, in read order.
#[derive(Debug, Clone, PartialEq)]
pub struct Segment {
    pub region_id: String,
    pub region_type: String,
    pub sequence_type: String,
    /// Expected length in this read (region min_len).
    pub len: usize,
    pub max_len: usize,
    /// Fixed sequence, oriented for this read's strand; the linker/adapter anchor.
    pub sequence: Option<String>,
    pub onlist: Option<Onlist>,
}

impl Assay {
    pub fn from_yaml(s: &str) -> Result<Assay> {
        serde_norway::from_str(s).context("failed to parse seqspec YAML")
    }

    /// Depth-first, in-order list of leaf regions (the physical layout of the molecule).
    pub fn leaf_regions(&self) -> Vec<&Region> {
        fn walk<'a>(regions: &'a [Region], out: &mut Vec<&'a Region>) {
            for r in regions {
                if r.is_leaf() {
                    out.push(r);
                } else {
                    walk(&r.regions, out);
                }
            }
        }
        let mut out = Vec::new();
        walk(&self.library_spec, &mut out);
        out
    }

    /// Whether the assay captures a UMI anywhere; false for full-length plate protocols
    /// (Smart-seq2/3), which quantify by read counts rather than deduplicated counts.
    pub fn has_umi(&self) -> bool {
        self.leaf_regions()
            .iter()
            .any(|r| is_umi_type(&r.region_type))
    }

    /// Resolve a read into the ordered segments it covers.
    ///
    /// Mirrors seqspec's map_read_id_to_regions: a pos read covers the regions after its
    /// primer; a neg read covers those before it, reversed and revcomp'd.
    pub fn read_layout(&self, read_id: &str) -> Option<Vec<Segment>> {
        let read = self.sequence_spec.iter().find(|r| r.read_id == read_id)?;
        let leaves = self.leaf_regions();
        let neg = read.strand == "neg";

        // A primer_id with no matching region (generic_primer1/2) means the molecule
        // boundary: read from the 5' end (pos) or the 3' end (neg).
        let ordered: Vec<&Region> = match leaves.iter().position(|r| r.region_id == read.primer_id)
        {
            Some(primer_idx) if neg => leaves[..primer_idx].iter().rev().copied().collect(),
            Some(primer_idx) => leaves[primer_idx + 1..].to_vec(),
            None if neg => leaves.iter().rev().copied().collect(),
            None => leaves.clone(),
        };

        let cap = read.max_len.max(1);
        let mut segs = Vec::new();
        let mut acc = 0usize;
        for r in ordered {
            if acc >= cap {
                break;
            }
            let sequence = if matches!(r.sequence_type.as_str(), "fixed") && !r.sequence.is_empty()
            {
                Some(if neg {
                    revcomp(&r.sequence)
                } else {
                    r.sequence.clone()
                })
            } else {
                None
            };
            let seg = Segment {
                region_id: r.region_id.clone(),
                region_type: r.region_type.clone(),
                sequence_type: r.sequence_type.clone(),
                len: r.min_len,
                max_len: r.max_len,
                sequence,
                onlist: r.onlist.clone(),
            };
            // The cDNA/gDNA insert (a variable, non-barcode/UMI/linker region) absorbs
            // the rest of the read, so nothing structural can follow it. Short variable
            // spacers and variable-length barcodes do NOT terminate the layout.
            let terminator = is_read_terminator(&seg);
            segs.push(seg);
            if terminator {
                break;
            }
            acc += r.min_len.max(1);
        }
        Some(segs)
    }

    /// The read that carries the cell barcode(s); the structural read the detector
    /// scores against. When an assay has several barcode-bearing reads (plate vs cell
    /// barcode, or RNA vs ATAC/protein modalities), prefer the one that best looks like
    /// the primary GEX cell-barcode read: most barcodes, barcodes backed by an onlist,
    /// the RNA modality, a UMI present, more fixed anchors, and no cDNA in the read.
    /// Returns (read_id, layout).
    pub fn structural_read(&self) -> Option<(String, Vec<Segment>)> {
        let mut best: Option<(String, Vec<Segment>, i64)> = None;
        for read in &self.sequence_spec {
            let Some(layout) = self.read_layout(&read.read_id) else {
                continue;
            };
            let n_bc = layout
                .iter()
                .filter(|s| is_barcode_type(&s.region_type))
                .count();
            if n_bc == 0 {
                continue;
            }
            let onlist_bc = layout
                .iter()
                .filter(|s| is_barcode_type(&s.region_type) && s.onlist.is_some())
                .count();
            let has_umi = layout.iter().any(|s| is_umi_type(&s.region_type));
            let has_cdna = layout.iter().any(|s| s.max_len > s.len);
            let n_anchor = layout.iter().filter(|s| s.sequence.is_some()).count();
            let is_rna = read.modality.eq_ignore_ascii_case("rna");
            let score = n_bc as i64 * 100
                + onlist_bc as i64 * 30
                + if is_rna { 40 } else { 0 }
                + if has_umi { 20 } else { 0 }
                + n_anchor as i64 * 5
                + if has_cdna { 0 } else { 25 };
            if best.as_ref().map(|b| score > b.2).unwrap_or(true) {
                best = Some((read.read_id.clone(), layout, score));
            }
        }
        best.map(|(id, layout, _)| (id, layout))
    }

    /// Onlist-backed barcode regions of the structural read. See
    /// [barcode_onlist_segments_of] for the offset convention.
    pub fn barcode_onlist_segments(&self) -> Vec<(usize, usize, usize, String)> {
        match self.structural_read() {
            Some((_, layout)) => barcode_onlist_segments_of(&layout),
            None => Vec::new(),
        }
    }
}

/// Onlist-backed barcode regions of a resolved structural-read layout, as
/// (offset, min_len, max_len, onlist_filename). Offsets are relative to the observed
/// read (the first barcode/UMI), matching cb_umi_offsets/starsolo. A variable-length
/// barcode at an exact offset is emitted (inDrop bc1); once any region is variable, every
/// downstream offset shifts with the real read, so those barcodes (inDrop bc2) are
/// unframed and dropped rather than scored out of frame.
pub fn barcode_onlist_segments_of(layout: &[Segment]) -> Vec<(usize, usize, usize, String)> {
    let Some(first) = layout
        .iter()
        .position(|s| is_barcode_or_umi(&s.region_type))
    else {
        return Vec::new();
    };
    let mut off = 0usize;
    let mut unframed = false;
    let mut out = Vec::new();
    for s in &layout[first..] {
        if is_insert_type(&s.region_type) {
            break;
        }
        if is_barcode_type(&s.region_type) && !unframed {
            if let Some(ol) = &s.onlist {
                if !ol.filename.is_empty() {
                    out.push((off, s.len, s.max_len.max(s.len), ol.filename.clone()));
                }
            }
        }
        if s.max_len > s.len {
            unframed = true;
        }
        off += s.len;
    }
    out
}

// Embedded seqspec catalogue, generated by build.rs from vendor/seqspec (submodule).
include!(concat!(env!("OUT_DIR"), "/specs.rs"));

/// Parse every embedded seqspec into an Assay. Specs that fail to parse are
/// skipped (returned in the error list) rather than aborting the catalogue.
pub fn embedded_assays() -> &'static [Assay] {
    static CACHE: std::sync::OnceLock<Vec<Assay>> = std::sync::OnceLock::new();
    CACHE.get_or_init(|| {
        EMBEDDED_SPECS
            .iter()
            .filter_map(|(name, yaml)| {
                let mut a = Assay::from_yaml(yaml).ok()?;
                a.spec_file = name.to_string();
                if a.assay_id.trim().is_empty() {
                    a.assay_id = name.trim_end_matches(".spec.yaml").to_string();
                }
                Some(a)
            })
            .collect()
    })
}

/// A single-cell barcode-read geometry derived from one catalogued assay's
/// structural (barcode-bearing) read; the seqspec-driven source of truth for "what
/// does a single-cell barcode read look like", used by size/composition/orientation.
#[derive(Debug, Clone)]
pub struct ScGeometry {
    pub spec_file: String,
    pub name: String,
    /// Expected length of the *observed* barcode read: the barcode/UMI/linker span
    /// after any leading primer/adapter and before the first variable (cDNA) region.
    pub read_len: usize,
}

impl ScGeometry {
    /// A human label: the assay name, falling back to the spec filename stem.
    pub fn label(&self) -> &str {
        if self.name.is_empty() {
            self.spec_file.trim_end_matches(".spec.yaml")
        } else {
            &self.name
        }
    }
}

/// Reduce a faithful structural-read layout to the *observed* barcode-read geometry.
///
/// seqspec's map_read_id_to_regions (which read_layout mirrors) yields every region
/// downstream of the read's primer; including any library handle/primer between the
/// primer and the barcode (e.g. the 39 bp Drop-seq Read1 custom primer), which is a
/// binding site, not sequenced. The sequenced barcode read therefore runs from the
/// first barcode/UMI region through the **last** barcode/UMI before the cDNA/gDNA
/// insert, *including* any adapters/linkers physically between them (they are
/// sequenced; e.g. scifi-RNA-seq's inter-barcode adapters). A trailing adapter after
/// the last barcode (a read-primer region a loose primer_id left in the read, e.g.
/// quartz-seq's read-2 adapter) is dropped. Variable-length barcodes (inDrop) are
/// kept. Returns those segments, or None if the read carries no barcode.
pub fn observed_span(layout: &[Segment]) -> Option<&[Segment]> {
    let is_bc_or_umi = |s: &Segment| is_barcode_type(&s.region_type) || is_umi_type(&s.region_type);
    let start = layout.iter().position(&is_bc_or_umi)?;
    // The window runs from the first barcode/UMI to the cDNA/gDNA insert (or the end).
    let end = layout[start..]
        .iter()
        .position(is_read_terminator)
        .map_or(layout.len(), |i| start + i);
    let window = &layout[start..end];
    // End at the last barcode/UMI; drop any trailing (post-barcode) adapter run.
    let last = window.iter().rposition(&is_bc_or_umi)?;
    let span = &window[..=last];
    span.iter()
        .any(|s| is_barcode_type(&s.region_type))
        .then_some(span)
}

fn observed_read_len(layout: &[Segment]) -> Option<usize> {
    Some(observed_span(layout)?.iter().map(|s| s.len).sum())
}

/// Every single-cell barcode-read geometry in the embedded catalogue (one per assay
/// with a barcode-bearing structural read). Parsed once and cached.
pub fn sc_geometries() -> &'static [ScGeometry] {
    static CACHE: std::sync::OnceLock<Vec<ScGeometry>> = std::sync::OnceLock::new();
    CACHE.get_or_init(|| {
        embedded_assays()
            .iter()
            .filter_map(|a| {
                let (_rid, layout) = a.structural_read()?;
                let read_len = observed_read_len(&layout)?;
                Some(ScGeometry {
                    spec_file: a.spec_file.clone(),
                    name: a.name.clone(),
                    read_len,
                })
            })
            .collect()
    })
}

/// Catalogued single-cell geometries whose observed barcode-read length is within
/// tol bp of len (empty for len == 0). The seqspec-driven source for "is this a
/// known single-cell barcode-read length", used by size/composition/orientation.
pub fn sc_geometries_matching(len: usize, tol: usize) -> Vec<&'static ScGeometry> {
    if len == 0 {
        return Vec::new();
    }
    sc_geometries()
        .iter()
        .filter(|g| g.read_len.abs_diff(len) <= tol)
        .collect()
}

/// True when a read's length profile looks like a fixed-length single-cell
/// barcode/UMI read (near-constant length matching a known barcode geometry),
/// as opposed to variable-length cDNA. Shared by composition and sex so the
/// "is this the barcode read?" test can't drift between them.
pub fn is_fixed_barcode_read(mean_len: f64, stddev: f64) -> bool {
    stddev < 1.0 && !sc_geometries_matching(mean_len.round() as usize, 1).is_empty()
}

/// STARsolo flags derived from an assay's seqspec geometry: SmartSeq for plate
/// assays with no in-read barcode, CB_UMI_Simple for a single contiguous cell
/// barcode + UMI (10x, Drop-seq), and CB_UMI_Complex (0-based read-anchored
/// soloCBposition/soloUMIposition) for split/combinatorial barcodes.
///
/// Leading sample-index/primer/adapter regions are stripped, so positions are
/// relative to the first sequenced barcode/UMI base. A variable-length region
/// upstream would shift the fixed offsets, so it is flagged rather than guessed.
pub fn starsolo(a: &Assay) -> String {
    const SMARTSEQ: &str = "--soloType SmartSeq";
    let Some((_, layout)) = a.structural_read() else {
        return SMARTSEQ.to_string();
    };
    let Some(first) = layout
        .iter()
        .position(|s| is_barcode_type(&s.region_type) || is_umi_type(&s.region_type))
    else {
        return SMARTSEQ.to_string();
    };

    let mut off = 0usize;
    let mut cbs: Vec<(usize, usize)> = Vec::new();
    let mut umi: Option<(usize, usize)> = None;
    let mut seen_variable = false;
    let mut unreliable = false;
    for s in &layout[first..] {
        if is_insert_type(&s.region_type) {
            break;
        }
        if is_barcode_type(&s.region_type) {
            unreliable |= seen_variable;
            cbs.push((off, s.len));
        } else if is_umi_type(&s.region_type) && umi.is_none() {
            unreliable |= seen_variable;
            umi = Some((off, s.len));
        }
        // A variable-length region shifts every offset after it; only barcode/UMI
        // blocks that follow one have unreliable fixed positions (a trailing polyA does not).
        seen_variable |= s.max_len > s.len;
        off += s.len;
    }
    if cbs.is_empty() {
        return SMARTSEQ.to_string();
    }
    let warn = if unreliable {
        "  # note: variable-length region upstream, verify positions".to_string()
    } else {
        String::new()
    };

    if cbs.len() == 1 {
        let (cs, cl) = cbs[0];
        return match umi {
            Some((us, ul)) => format!(
                "--soloType CB_UMI_Simple --soloCBstart {} --soloCBlen {} --soloUMIstart {} --soloUMIlen {}{warn}",
                cs + 1, cl, us + 1, ul
            ),
            // CB_UMI_Simple requires a UMI; a CB-only protocol (e.g. Quartz-seq v1) has none.
            None => format!(
                "--soloType CB_UMI_Simple --soloCBstart {} --soloCBlen {}  # note: no UMI in geometry; read-counted, not UMI-deduplicated{warn}",
                cs + 1, cl
            ),
        };
    }

    let cb_pos = cbs
        .iter()
        .map(|&(s, l)| format!("0_{}_0_{}", s, s + l - 1))
        .collect::<Vec<_>>()
        .join(" ");
    let umi_pos = umi
        .map(|(s, l)| format!(" --soloUMIposition 0_{}_0_{}", s, s + l - 1))
        .unwrap_or_default();
    format!("--soloType CB_UMI_Complex --soloCBposition {cb_pos}{umi_pos}{warn}")
}

/// Reverse-complement a fixed sequence (IUPAC-aware for the common bases).
fn revcomp(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            b'N' => 'N',
            other => other as char,
        })
        .collect()
}

#[cfg(test)]
mod tests;
