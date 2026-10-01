use super::*;
use crate::testutil::write_fastq;

/// A distinct 16-bp cell barcode for index i (base-4 encoding of i).
fn distinct_cb(i: usize) -> Vec<u8> {
    (0..16)
        .map(|j| b"ACGT"[(i / 4usize.pow(j as u32)) % 4] as char as u8)
        .collect()
}

/// n synthetic 10x R1 reads (16-bp CB + 12-bp UMI = 28 bp) and matching R2 reads
/// (TSO at the start when tso is set). Returns (r1, r2, cb_set).
#[allow(clippy::type_complexity)]
fn tenx_v3_reads(n: usize, tso: bool) -> (Vec<Vec<u8>>, Vec<Vec<u8>>, HashSet<Vec<u8>>) {
    let mut r1 = Vec::with_capacity(n);
    let mut cbs = HashSet::new();
    for i in 0..n {
        let cb = distinct_cb(i);
        let mut read = cb.clone();
        read.extend_from_slice(&[b'A'; 12]); // UMI
        r1.push(read);
        cbs.insert(cb);
    }
    let r2: Vec<Vec<u8>> = (0..n)
        .map(|_| {
            let mut read = Vec::new();
            if tso {
                read.extend_from_slice(TSO_PRESENT.seq);
            }
            read.resize(90, b'C');
            read
        })
        .collect();
    (r1, r2, cbs)
}

#[test]
fn tenx_v3_with_matching_whitelist_is_high_confidence() {
    let (r1, r2, cbs) = tenx_v3_reads(200, true);
    let mut wl_sets: HashMap<&'static str, HashSet<Vec<u8>>> = HashMap::new();
    wl_sets.insert("3M-february-2018.txt.gz", cbs);

    let a = analyze(&r1, &r2, &wl_sets, None);
    assert_eq!(
        a.scored[0].0.id, "10x_3prime_v3",
        "winner was {}",
        a.scored[0].0.id
    );
    assert_eq!(a.confidence, "high");
}

#[test]
fn tenx_3prime_without_whitelist_picks_3prime_family_at_medium() {
    // TSO present but no whitelist cache: should land on a 3' assay, medium.
    let (r1, r2, _) = tenx_v3_reads(200, true);
    let a = analyze(&r1, &r2, &HashMap::new(), Some(200));
    assert!(
        a.scored[0].0.id.starts_with("10x_3prime"),
        "winner was {}",
        a.scored[0].0.id
    );
    assert_eq!(a.confidence, "medium");
}

#[test]
fn dropseq_geometry_is_medium_and_correct() {
    // 20 bp R1 (12 CB + 8 UMI), no TSO, no whitelist -> Drop-seq, geometry-only.
    let r1: Vec<Vec<u8>> = (0..200)
        .map(|i| {
            let mut v = distinct_cb(i); // 16 distinct bases
            v.resize(20, b'A'); // pad to the 20 bp Drop-seq geometry
            v
        })
        .collect();
    let r2: Vec<Vec<u8>> = (0..200).map(|_| vec![b'C'; 90]).collect();
    let a = analyze(&r1, &r2, &HashMap::new(), Some(0));
    assert_eq!(
        a.scored[0].0.id, "drop_seq",
        "winner was {}",
        a.scored[0].0.id
    );
    assert_eq!(a.confidence, "medium");
}

#[test]
fn bulk_wgs_reads_are_not_single_cell() {
    // 100 bp genomic R1, no TSO, no whitelist: must NOT be force-fit to a 10x assay.
    let r1: Vec<Vec<u8>> = (0..200).map(|_| vec![b'A'; 100]).collect();
    let r2: Vec<Vec<u8>> = (0..200).map(|_| vec![b'C'; 100]).collect();
    let a = analyze(&r1, &r2, &HashMap::new(), Some(0));
    assert!(
        !a.is_single_cell,
        "100bp genomic reads should not look single-cell"
    );
    let report = build_report(
        a,
        false,
        200,
        false,
        crate::endpanel::classify(&[], None),
        crate::atacpanel::classify(&[], None),
    );
    assert_eq!(report.assay_id, "none");
    assert_eq!(report.assay_display, "Not single-cell");
    assert_eq!(report.confidence, "low");
    // The ledger walks all three signals and ends on a clear verdict.
    assert!(
        report
            .interpretation
            .iter()
            .any(|e| e.signal == "Barcode geometry"),
        "reasoning should explain the barcode-geometry check"
    );
    assert!(
        report
            .interpretation
            .last()
            .is_some_and(|e| e.signal == "Verdict" && e.conclusion.contains("detect-assay")),
        "reasoning should end on a not-10x verdict that points at detect-assay"
    );
    assert!(report.notes.iter().any(|n| n.contains("detect-assay")));
}

#[test]
fn unmatched_barcode_length_abstains_instead_of_force_fitting() {
    // 33 bp (sci-RNA-seq3): single-cell, but no 10x or Drop-seq geometry within +/-2 bp.
    let r1: Vec<Vec<u8>> = (0..200)
        .map(|i| {
            let mut v = distinct_cb(i);
            v.resize(33, b'A');
            v
        })
        .collect();
    let r2: Vec<Vec<u8>> = (0..200).map(|_| vec![b'C'; 90]).collect();
    let a = analyze(&r1, &r2, &HashMap::new(), Some(0));
    assert!(a.is_single_cell, "33bp is a catalogued barcode-read length");
    assert!(
        a.scored.is_empty(),
        "no 10x/Drop-seq candidate should score"
    );
    let report = build_report(
        a,
        false,
        200,
        false,
        crate::endpanel::classify(&[], None),
        crate::atacpanel::classify(&[], None),
    );
    assert_eq!(report.assay_id, "unknown");
    assert_eq!(report.confidence, "low");
    assert!(
        report.notes.iter().any(|n| n.contains("detect-assay")),
        "abstention should point at the seqspec tier: {:?}",
        report.notes
    );
    assert!(
        report.notes.iter().any(|n| n.contains("sci-RNA-seq3")),
        "abstention should name the matching seqspec geometries: {:?}",
        report.notes
    );
    assert!(report.tool_hints.is_empty(), "no aligner recipe to offer");
}

#[test]
fn onlist_miss_overrides_a_coincidental_length_match() {
    // 30 bp (mCEL-Seq2) shortlists 10x by length, but onlist overlap is near zero.
    let r1: Vec<Vec<u8>> = (0..200)
        .map(|i| {
            let mut v = distinct_cb(i);
            v.resize(30, b'A');
            v
        })
        .collect();
    let r2: Vec<Vec<u8>> = (0..200).map(|_| vec![b'C'; 90]).collect();
    // onlist sharing none of the observed barcodes
    let mut wl_sets: HashMap<&'static str, HashSet<Vec<u8>>> = HashMap::new();
    for f in [
        "3M-february-2018.txt.gz",
        "3M-3pgex-may-2023.txt.gz",
        "3M-5pgex-jan-2023.txt.gz",
    ] {
        wl_sets.insert(f, HashSet::from([vec![b'T'; 16]]));
    }
    let a = analyze(&r1, &r2, &wl_sets, None);
    assert!(!a.scored.is_empty(), "10x is shortlisted by length");
    let report = build_report(
        a,
        false,
        200,
        false,
        crate::endpanel::classify(&[], None),
        crate::atacpanel::classify(&[], None),
    );
    assert_eq!(report.assay_id, "unknown", "onlist miss must veto the call");
    assert_eq!(report.confidence, "low");
    assert!(
        report.notes.iter().any(|n| n.contains("contradicts")),
        "note should say the onlist contradicts the geometry: {:?}",
        report.notes
    );
}

#[test]
fn ambiguous_3p_5p_polarity_is_reported_as_both() {
    let n = 200;
    let mut r1 = Vec::new();
    let mut cbs = HashSet::new();
    for i in 0..n {
        let cb = distinct_cb(i);
        let mut read = cb.clone();
        read.extend_from_slice(&[b'A'; 10]); // UMI -> 26 bp total
        r1.push(read);
        cbs.insert(cb);
    }
    let r2: Vec<Vec<u8>> = (0..n).map(|_| vec![b'C'; 90]).collect();
    let mut wl_sets: HashMap<&'static str, HashSet<Vec<u8>>> = HashMap::new();
    wl_sets.insert("737K-august-2016.txt.gz", cbs);

    let a = analyze(&r1, &r2, &wl_sets, None);
    let report = build_report(
        a,
        false,
        n,
        false,
        crate::endpanel::classify(&[], None),
        crate::atacpanel::classify(&[], None),
    );
    assert!(
        report.assay_display.contains(" or "),
        "polarity should be reported as both, got {}",
        report.assay_display
    );
    assert_ne!(
        report.confidence, "high",
        "unresolved polarity can't be high"
    );
    assert!(report
        .interpretation
        .iter()
        .any(|e| e.signal == "3' vs 5' polarity"));
}

#[test]
fn tenx_v3_reads_are_single_cell() {
    // Guard against the gate over-triggering on genuine single-cell geometry.
    let (r1, r2, _) = tenx_v3_reads(200, true);
    let a = analyze(&r1, &r2, &HashMap::new(), Some(180));
    assert!(a.is_single_cell);
}

#[test]
fn candidates_for_includes_within_tolerance() {
    let ids: Vec<&str> = candidates_for(26).iter().map(|s| s.id).collect();
    assert!(ids.contains(&"10x_3prime_v2")); // exact 26 bp
    assert!(ids.contains(&"10x_3prime_v3")); // 28 bp, within +/-2
    assert!(!ids.contains(&"drop_seq")); // 20 bp, outside tolerance
}

#[test]
fn len_profile_reports_the_mode() {
    let seqs = vec![vec![0u8; 28], vec![0u8; 28], vec![0u8; 26]];
    let (dom, dist) = len_profile(&seqs);
    assert_eq!(dom, 28);
    assert_eq!(dist[0], (28, 2));
}

#[test]
fn swapped_inputs_report_barcode_read_as_r2() {
    let (barcode_r1, cdna_r2, _) = tenx_v3_reads(200, true);
    let long_read: Vec<Vec<u8>> = cdna_r2;
    let short_read: Vec<Vec<u8>> = barcode_r1;
    let a = analyze(&short_read, &long_read, &HashMap::new(), Some(180));
    let report = build_report(
        a,
        true,
        200,
        false,
        crate::endpanel::classify(&[], None),
        crate::atacpanel::classify(&[], None),
    );
    assert_eq!(report.cb_read, "R2");
}

/// 16CB+10UMI = 26 bp R1, no TSO: geometry alone can't tell 3'v2 from 5'v1.
fn ambiguous_26bp_reads(n: usize) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let mut r1 = Vec::with_capacity(n);
    for i in 0..n {
        let mut read = distinct_cb(i);
        read.extend_from_slice(&[b'A'; 10]); // 10 bp UMI -> 26 bp
        r1.push(read);
    }
    let r2: Vec<Vec<u8>> = (0..n).map(|_| vec![b'C'; 90]).collect();
    (r1, r2)
}

fn end_report(end: &str, n3: usize, n5: usize) -> crate::endpanel::EndReport {
    crate::endpanel::EndReport {
        end: end.into(),
        confidence: "medium".into(),
        n5,
        n3,
        reads_scanned: n3 + n5,
        ratio: n3 as f64 / (n3 + n5) as f64,
        notes: vec![],
    }
}

#[test]
fn decisive_five_prime_end_corrects_registry_default() {
    let (r1, r2) = ambiguous_26bp_reads(500);
    let a = analyze(&r1, &r2, &HashMap::new(), Some(0));
    assert!(a.is_single_cell);
    let report = build_report(
        a,
        false,
        500,
        false,
        end_report("5'", 100, 300),
        crate::atacpanel::classify(&[], None),
    );
    assert_eq!(
        report.assay_id, "10x_5prime_v1",
        "display={}",
        report.assay_display
    );
    assert_eq!(report.end, "5'");
    assert!(report.notes.iter().any(|n| n.contains("corrected")));
}

#[test]
fn decisive_three_prime_end_keeps_three_prime() {
    let (r1, r2) = ambiguous_26bp_reads(500);
    let a = analyze(&r1, &r2, &HashMap::new(), Some(0));
    let report = build_report(
        a,
        false,
        500,
        false,
        end_report("3'", 300, 100),
        crate::atacpanel::classify(&[], None),
    );
    assert_eq!(report.assay_id, "10x_3prime_v2");
    assert_eq!(report.end, "3'");
    assert_eq!(report.confidence, "medium");
    assert!(report.notes.iter().any(|n| n.contains("confirms")));
}

#[test]
fn empty_whitelist_cache_keeps_sprt_fast_path() {
    let (r1_reads, r2_reads, _) = tenx_v3_reads(8_000, true);
    let r1 = write_fastq("chem_detect_r1", &r1_reads);
    let r2 = write_fastq("chem_detect_r2", &r2_reads);
    let cache_dir = std::env::temp_dir().join(format!(
        "preflightx-empty-cache-{}",
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&cache_dir).unwrap();

    let report = detect(&r1, &r2, 8_000, Some(&cache_dir), None).unwrap();

    std::fs::remove_file(&r1).ok();
    std::fs::remove_file(&r2).ok();
    std::fs::remove_dir_all(&cache_dir).ok();

    // A cached end/atac panel forces a full sample, so only assert early-stop without one.
    if crate::endpanel::have_panels() || crate::atacpanel::have_panels() {
        assert!(
            !report.early_stop,
            "a cached k-mer panel forces a full sample, not SPRT early-stop"
        );
        return;
    }
    assert!(
        report.early_stop,
        "expected SPRT early stop with no whitelist data"
    );
    assert!(
        report.reads_sampled < 8_000,
        "expected fewer than the full read budget, got {}",
        report.reads_sampled
    );
}

/// deterministic motif-free bases
fn pseudo_random(seed: usize, len: usize) -> Vec<u8> {
    let mut x = (seed as u64)
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (0..len)
        .map(|_| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(x >> 33) as usize % 4]
        })
        .collect()
}

/// 150 bp R1: CB16 + UMI + after + filler, over 40 cells
fn long_r1_reads(n: usize, umi_len: usize, after: &[u8]) -> (Vec<Vec<u8>>, HashSet<Vec<u8>>) {
    let mut cbs = HashSet::new();
    let reads = (0..n)
        .map(|i| {
            let cb = distinct_cb(i % 40);
            cbs.insert(cb.clone());
            let mut r = cb;
            r.extend_from_slice(&pseudo_random(i, umi_len));
            r.extend_from_slice(after);
            let rest = 150 - r.len();
            r.extend_from_slice(&pseudo_random(i + 1_000_000, rest));
            r
        })
        .collect();
    (reads, cbs)
}

#[test]
fn long_r1_is_long() {
    assert!(is_long_r1(150));
    assert!(is_long_r1(151));
    assert!(is_long_r1(100));
    // seqspec barcode read lengths
    assert!(!is_long_r1(28));
    assert!(!is_long_r1(33));
    assert!(!is_long_r1(50));
}

#[test]
fn long_r1_10x_5prime_tso_is_single_cell_without_cache() {
    // TSO in R1 only
    let (r1, _) = long_r1_reads(400, 10, R1_TSO_5P);
    let r2: Vec<Vec<u8>> = (0..400).map(|i| pseudo_random(i + 7, 150)).collect();
    let a = analyze(&r1, &r2, &HashMap::new(), Some(0));
    assert!(a.long_r1);
    assert!(a.r1_tso_pct > 90.0, "r1_tso_pct={}", a.r1_tso_pct);
    assert!(a.is_single_cell, "long-R1 5' library must be single-cell");
    assert!(
        a.scored[0].0.id.starts_with("10x_5prime"),
        "winner {}",
        a.scored[0].0.id
    );
}

#[test]
fn long_r1_10x_3prime_polyt_with_whitelist_is_v3() {
    let (r1, cbs) = long_r1_reads(400, 12, &[b'T'; 30]);
    let r2: Vec<Vec<u8>> = (0..400).map(|i| pseudo_random(i + 9, 150)).collect();
    let mut wl_sets: HashMap<&'static str, HashSet<Vec<u8>>> = HashMap::new();
    wl_sets.insert("3M-february-2018.txt.gz", cbs);
    let a = analyze(&r1, &r2, &wl_sets, Some(0));
    assert!(a.is_single_cell);
    assert!(a.max_wl_pct >= 99.0);
    assert!(
        a.scored[0].0.id.starts_with("10x_3prime"),
        "winner {}",
        a.scored[0].0.id
    );
    let report = build_report(
        a,
        false,
        400,
        false,
        crate::endpanel::classify(&[], None),
        crate::atacpanel::classify(&[], None),
    );
    assert!(report.is_single_cell && report.long_r1);
    assert!(
        report.notes.iter().any(|n| n.contains("cDNA run-on")),
        "notes={:?}",
        report.notes
    );
}

#[test]
fn long_r1_random_bulk_stays_bulk() {
    let r1: Vec<Vec<u8>> = (0..400).map(|i| pseudo_random(i, 150)).collect();
    let r2: Vec<Vec<u8>> = (0..400).map(|i| pseudo_random(i + 5_000, 150)).collect();
    let a = analyze(&r1, &r2, &HashMap::new(), Some(0));
    assert!(a.long_r1 && !a.is_single_cell);
}

#[test]
fn long_r1_oligo_dt_without_barcode_reuse_is_not_single_cell() {
    // 3'-tag bulk: fixed-offset oligo-dT, unique prefixes
    let r1: Vec<Vec<u8>> = (0..400)
        .map(|i| {
            let mut r = pseudo_random(i, 28);
            r.extend_from_slice(&[b'T'; 30]);
            r.extend_from_slice(&pseudo_random(i + 9_999, 92));
            r
        })
        .collect();
    let r2: Vec<Vec<u8>> = (0..400).map(|i| pseudo_random(i + 3, 150)).collect();
    let a = analyze(&r1, &r2, &HashMap::new(), Some(0));
    assert!(a.r1_polyt_pct > 90.0);
    assert!(!a.is_single_cell);
}

#[test]
fn parse_whitelist_keeps_only_sampled_barcodes() {
    let onlist = b"AAAACCCCGGGGTTTT\nacgtacgtacgtacgt-1\nTTTTGGGGCCCCAAAA\nshort\n";
    let keep: HashSet<Vec<u8>> = [b"ACGTACGTACGTACGT".to_vec(), b"GGGGGGGGGGGGGGGG".to_vec()]
        .into_iter()
        .collect();
    let got = parse_whitelist(&onlist[..], 16, &keep).unwrap();
    let want: HashSet<Vec<u8>> = [b"ACGTACGTACGTACGT".to_vec()].into_iter().collect();
    assert_eq!(
        got, want,
        "lowercase line uppercased; unsampled and short lines dropped"
    );
}

#[test]
fn filtered_whitelist_gives_same_overlap_as_full_onlist() {
    // pre-filtering the onlist to the sample must not change the score
    let (r1, r2, cbs) = tenx_v3_reads(200, true);
    let mut full = cbs.clone();
    for i in 1_000..5_000 {
        full.insert(distinct_cb(i));
    }
    let keep: HashSet<Vec<u8>> = extract_cb_candidates(&r1, 0, 16);
    let text: Vec<u8> = full
        .iter()
        .flat_map(|b| [b.as_slice(), b"\n"].concat())
        .collect();
    let filtered = parse_whitelist(&text[..], 16, &keep).unwrap();
    let mut a_sets: HashMap<&'static str, HashSet<Vec<u8>>> = HashMap::new();
    a_sets.insert("3M-february-2018.txt.gz", full);
    let mut b_sets: HashMap<&'static str, HashSet<Vec<u8>>> = HashMap::new();
    b_sets.insert("3M-february-2018.txt.gz", filtered);
    let a = analyze(&r1, &r2, &a_sets, None);
    let b = analyze(&r1, &r2, &b_sets, None);
    assert_eq!(a.wl_scores, b.wl_scores);
    assert_eq!(a.scored[0].0.id, b.scored[0].0.id);
    assert_eq!(a.confidence, b.confidence);
}
