#[derive(Clone, Copy)]
pub struct WhitelistDef {
    pub filename: &'static str,
    pub url: &'static str,
    pub cb_len: usize,
}

#[derive(Clone, Copy)]
pub struct AnchorMotif {
    pub seq: &'static [u8],
    pub scan_start: usize,
    pub scan_end: usize,
    pub min_match_len: usize,
    /// Expected hit rate (0 to 1) for reads from an assay where this anchor is canonical.
    /// Use 0.20 when the motif should be present (3' GEX TSO); 0.01 when absent (5' GEX,
    /// background noise). 0.01 is the baseline, so TSO_ABSENT contributes LLR = 0.
    pub p_signal: f64,
}

pub struct AssaySpec {
    pub id: &'static str,
    pub display_name: &'static str,
    pub r1_cb_offset: usize,
    pub r1_cb_len: usize,
    pub r1_umi_len: usize,
    pub whitelist: Option<WhitelistDef>,
    pub r2_anchors: &'static [AnchorMotif],
    /// Per-tool hints that cannot be derived from geometry (e.g. kb protocol names,
    /// cellranger chemistry strings). STARsolo flags are generated from geometry via
    /// `starsolo_hint()` and are not stored here.
    pub tool_hints: &'static [(&'static str, &'static str)],
}

impl AssaySpec {
    pub fn r1_expected_len(&self) -> usize {
        self.r1_cb_offset + self.r1_cb_len + self.r1_umi_len
    }
    pub fn umi_offset(&self) -> usize {
        self.r1_cb_offset + self.r1_cb_len
    }
    /// Generate STARsolo CB_UMI_Simple flags from this assay's geometry.
    /// Keeps the hint in sync with the spec automatically.
    pub fn starsolo_hint(&self) -> String {
        let wl = self
            .whitelist
            .as_ref()
            .map(|w| format!("--soloCBwhitelist {} ", w.filename))
            .unwrap_or_else(|| "--soloCBwhitelist None ".to_string());
        format!(
            "--soloType CB_UMI_Simple {}--soloCBstart {} --soloCBlen {} --soloUMIstart {} --soloUMIlen {}",
            wl,
            self.r1_cb_offset + 1,
            self.r1_cb_len,
            self.umi_offset() + 1,
            self.r1_umi_len,
        )
    }
}

const WL_737K_2014: WhitelistDef = WhitelistDef {
    filename: "737K-april-2014_rc.txt.gz",
    url: "https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/737K-april-2014_rc.txt.gz",
    cb_len: 14,
};
const WL_737K_2016: WhitelistDef = WhitelistDef {
    filename: "737K-august-2016.txt.gz",
    url: "https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/737K-august-2016.txt.gz",
    cb_len: 16,
};
const WL_3M_FEB2018: WhitelistDef = WhitelistDef {
    filename: "3M-february-2018.txt.gz",
    url: "https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/3M-february-2018.txt.gz",
    cb_len: 16,
};
const WL_3M_3PGEX2023: WhitelistDef = WhitelistDef {
    filename: "3M-3pgex-may-2023.txt.gz",
    url: "https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/3M-3pgex-may-2023.txt.gz",
    cb_len: 16,
};
const WL_3M_5PGEX2023: WhitelistDef = WhitelistDef {
    filename: "3M-5pgex-jan-2023.txt.gz",
    url: "https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/3M-5pgex-jan-2023.txt.gz",
    cb_len: 16,
};

const TSO: &[u8] = b"AAGCAGTGGTATCAACGCAGAGT";

// In 3' GEX, R2 reads from the P7-to-TSO end of the cDNA: TSO appears at the START of R2.
// In 5' GEX, R2 reads from the P7-to-polyA end: TSO is absent from R2[0..30].
// p_signal for PRESENT = 0.20 (typical 3' TSO rate); for ABSENT = 0.01 (background, near baseline).
pub const TSO_PRESENT: AnchorMotif = AnchorMotif {
    seq: TSO,
    scan_start: 0,
    scan_end: 30,
    min_match_len: 13,
    p_signal: 0.20,
};
const TSO_ABSENT: AnchorMotif = AnchorMotif {
    seq: TSO,
    scan_start: 0,
    scan_end: 30,
    min_match_len: 13,
    p_signal: 0.01, // same as P_TSO_BACKGROUND, so LLR contribution is zero
};

static A_TSO_PRESENT: [AnchorMotif; 1] = [TSO_PRESENT];
static A_TSO_ABSENT: [AnchorMotif; 1] = [TSO_ABSENT];
static A_NONE: [AnchorMotif; 0] = [];

pub static ASSAYS: &[AssaySpec] = &[
    AssaySpec {
        id: "10x_3prime_v1",
        display_name: "10x Chromium 3' v1",
        r1_cb_offset: 0,
        r1_cb_len: 14,
        r1_umi_len: 10,
        whitelist: Some(WL_737K_2014),
        r2_anchors: &A_NONE,
        tool_hints: &[("kb", "-x 10xv1")],
    },
    AssaySpec {
        id: "10x_3prime_v2",
        display_name: "10x Chromium 3' v2",
        r1_cb_offset: 0,
        r1_cb_len: 16,
        r1_umi_len: 10,
        whitelist: Some(WL_737K_2016),
        r2_anchors: &A_TSO_PRESENT,
        tool_hints: &[
            ("kb", "-x 10xv2"),
            ("alevin", "--chromiumV2"),
            ("cellranger", "--chemistry SC3Pv2"),
        ],
    },
    AssaySpec {
        id: "10x_5prime_v1",
        display_name: "10x Chromium 5' v1",
        r1_cb_offset: 0,
        r1_cb_len: 16,
        r1_umi_len: 10,
        whitelist: Some(WL_737K_2016),
        r2_anchors: &A_TSO_ABSENT,
        tool_hints: &[("cellranger", "--chemistry SC5P-PE")],
    },
    AssaySpec {
        id: "10x_5prime_v2",
        display_name: "10x Chromium 5' v2",
        r1_cb_offset: 0,
        r1_cb_len: 16,
        r1_umi_len: 10,
        whitelist: Some(WL_3M_5PGEX2023),
        r2_anchors: &A_TSO_ABSENT,
        tool_hints: &[("cellranger", "--chemistry SC5P-PE-v2")],
    },
    AssaySpec {
        id: "10x_3prime_v3",
        display_name: "10x Chromium 3' v3",
        r1_cb_offset: 0,
        r1_cb_len: 16,
        r1_umi_len: 12,
        whitelist: Some(WL_3M_FEB2018),
        r2_anchors: &A_TSO_PRESENT,
        tool_hints: &[
            ("kb", "-x 10xv3"),
            ("alevin", "--chromiumV3"),
            ("cellranger", "--chemistry SC3Pv3"),
        ],
    },
    AssaySpec {
        id: "10x_3prime_v4",
        display_name: "10x Chromium 3' v4 (GEM-X)",
        r1_cb_offset: 0,
        r1_cb_len: 16,
        r1_umi_len: 12,
        whitelist: Some(WL_3M_3PGEX2023),
        r2_anchors: &A_TSO_PRESENT,
        tool_hints: &[("cellranger", "--chemistry SC3Pv4")],
    },
    AssaySpec {
        id: "10x_5prime_v3",
        display_name: "10x Chromium 5' v3",
        r1_cb_offset: 0,
        r1_cb_len: 16,
        r1_umi_len: 12,
        whitelist: Some(WL_3M_5PGEX2023),
        r2_anchors: &A_TSO_ABSENT,
        tool_hints: &[("cellranger", "--chemistry SC5P-PE-v3")],
    },
    AssaySpec {
        id: "drop_seq",
        display_name: "Drop-seq",
        r1_cb_offset: 0,
        r1_cb_len: 12,
        r1_umi_len: 8,
        whitelist: None,
        r2_anchors: &A_NONE,
        tool_hints: &[("kb", "-x DROPSEQ")],
    },
];

#[cfg(test)]
mod tests {
    use super::*;

    fn spec(id: &str) -> &'static AssaySpec {
        ASSAYS.iter().find(|s| s.id == id).expect("assay present")
    }

    #[test]
    fn r1_expected_len_sums_geometry() {
        // 16 CB + 12 UMI = 28 (3' v3); 12 CB + 8 UMI = 20 (Drop-seq).
        assert_eq!(spec("10x_3prime_v3").r1_expected_len(), 28);
        assert_eq!(spec("drop_seq").r1_expected_len(), 20);
    }

    #[test]
    fn starsolo_hint_tracks_geometry() {
        // 1-based CB/UMI start + lengths derived from the spec.
        let hint = spec("10x_3prime_v3").starsolo_hint();
        assert!(hint.contains("--soloCBstart 1 --soloCBlen 16"));
        assert!(hint.contains("--soloUMIstart 17 --soloUMIlen 12"));
        // Whitelisted assays name their whitelist; Drop-seq has none.
        assert!(hint.contains("--soloCBwhitelist 3M-february-2018.txt.gz"));
        assert!(spec("drop_seq")
            .starsolo_hint()
            .contains("--soloCBwhitelist None"));
    }
}
