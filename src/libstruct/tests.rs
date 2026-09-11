use super::*;

// Re-derive the expected read-relative CB/UMI offsets straight from the vendored
// seqspec layout, so the test tracks the spec files rather than hand-written numbers.
#[allow(clippy::type_complexity)]
fn cb_umi_offsets(a: &Assay) -> (Vec<(usize, usize)>, Option<(usize, usize)>) {
    let (_, layout) = a.structural_read().unwrap();
    let first = layout
        .iter()
        .position(|s| is_barcode_type(&s.region_type) || is_umi_type(&s.region_type))
        .unwrap();
    let (mut off, mut cbs, mut umi) = (0usize, Vec::new(), None);
    for s in &layout[first..] {
        if is_insert_type(&s.region_type) {
            break;
        }
        if is_barcode_type(&s.region_type) {
            cbs.push((off, s.len));
        } else if is_umi_type(&s.region_type) && umi.is_none() {
            umi = Some((off, s.len));
        }
        off += s.len;
    }
    (cbs, umi)
}

#[test]
fn barcode_onlist_segments_from_layout() {
    let assays = embedded_assays();
    let by = |f: &str| assays.iter().find(|a| a.spec_file == f).unwrap();

    // SPLiT-seq is all fixed-length, so every barcode is framed and its offsets must
    // equal the validated CB geometry from cb_umi_offsets, in read order (bc3..bc1).
    let split = by("split_seq.spec.yaml").barcode_onlist_segments();
    let geom: Vec<(usize, usize)> = split.iter().map(|(o, l, _, _)| (*o, *l)).collect();
    let (split_cbs, _) = cb_umi_offsets(by("split_seq.spec.yaml"));
    assert_eq!(geom, split_cbs);
    assert_eq!(
        split
            .iter()
            .map(|(_, _, _, f)| f.as_str())
            .collect::<Vec<_>>(),
        [
            "barcode-3_onlist.txt",
            "barcode-2_onlist.txt",
            "barcode-1_onlist.txt"
        ]
    );

    // inDrop bc1 is variable (8-11 bp); bc2 sits downstream of it, so its offset can't
    // be framed and it is dropped. Only the exactly-framed bc1 is emitted.
    let ind = by("indropv2.spec.yaml").barcode_onlist_segments();
    assert_eq!(
        ind.iter()
            .map(|(_, _, _, f)| f.as_str())
            .collect::<Vec<_>>(),
        ["bc1_onlist.txt"]
    );
    assert_eq!(ind[0].0, 0, "inDrop bc1 sits at the read start");
    assert_eq!((ind[0].1, ind[0].2), (8, 11), "inDrop bc1 is 8-11 bp");

    // Plate / no-onlist-barcode assays yield nothing.
    assert!(by("smart_seq2.spec.yaml")
        .barcode_onlist_segments()
        .is_empty());
}

#[test]
fn starsolo_flags_derive_from_geometry() {
    let assays = embedded_assays();
    let by = |f: &str| assays.iter().find(|a| a.spec_file == f).unwrap();

    // Single contiguous CB + UMI -> CB_UMI_Simple (1-based start, length).
    assert_eq!(
        starsolo(by("10x_rna_v3.spec.yaml")),
        "--soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12"
    );
    // No in-read barcode (plate) -> SmartSeq.
    assert_eq!(starsolo(by("smart_seq2.spec.yaml")), "--soloType SmartSeq");
    // CB-only protocol (no UMI) flags that STARsolo can't UMI-dedup it.
    assert!(starsolo(by("quartz_seq.spec.yaml")).contains("no UMI"));

    // Complex assays: every position emitted must match the offsets the seqspec
    // layout actually reports (0-based inclusive), for each barcode block and the UMI.
    for spec in [
        "split_seq.spec.yaml",
        "indropv2.spec.yaml",
        "sci_rna_seq3.spec.yaml",
    ] {
        let a = by(spec);
        let out = starsolo(a);
        assert!(out.contains("--soloType CB_UMI_Complex"), "{spec}: {out}");
        let (cbs, umi) = cb_umi_offsets(a);
        for (s, l) in cbs {
            assert!(
                out.contains(&format!("0_{s}_0_{}", s + l - 1)),
                "{spec}: {out}"
            );
        }
        let (us, ul) = umi.unwrap();
        assert!(
            out.contains(&format!("--soloUMIposition 0_{us}_0_{}", us + ul - 1)),
            "{spec}: {out}"
        );
    }
}

#[test]
fn has_umi_distinguishes_droplet_from_plate() {
    let assays = embedded_assays();
    let find = |f: &str| {
        assays
            .iter()
            .find(|a| a.spec_file == f)
            .unwrap_or_else(|| panic!("spec {f} not embedded"))
    };
    assert!(
        find("10x_rna_v3.spec.yaml").has_umi(),
        "10x v3 is UMI-based"
    );
    assert!(
        !find("smart_seq2.spec.yaml").has_umi(),
        "Smart-seq2 is full-length, no UMI"
    );
    assert!(assays.iter().any(|a| a.has_umi()));
    assert!(assays.iter().any(|a| !a.has_umi()));
}

fn geom_of(spec_file: &str) -> Option<usize> {
    sc_geometries()
        .iter()
        .find(|g| g.spec_file == spec_file)
        .map(|g| g.read_len)
}

/// Observed barcode-read span as Vec, for the doc generators.
fn observed_span(layout: &[Segment]) -> Vec<&Segment> {
    super::observed_span(layout).map_or_else(Vec::new, |s| s.iter().collect())
}

/// Regenerate docs/assay-geometries.md from the catalogue. Run on demand:
///   cargo test gen_assay_docs -- --ignored
#[test]
#[ignore = "doc generator; run explicitly to refresh docs/assay-geometries.md"]
fn gen_assay_docs() {
    fn short(region_type: &str) -> String {
        if is_barcode_type(region_type) {
            "CB".into()
        } else if is_umi_type(region_type) {
            "UMI".into()
        } else if is_linker_type(region_type) {
            "link".into()
        } else {
            region_type.to_ascii_lowercase()
        }
    }
    fn ascii(a: &Assay) -> Option<(String, usize, String)> {
        let (rid, layout) = a.structural_read()?;
        let span = observed_span(&layout);
        if span.is_empty() {
            return None;
        }
        let len: usize = span.iter().map(|s| s.len).sum();
        let boxes: Vec<String> = span
            .iter()
            .map(|s| format!("[ {} {} ]", short(&s.region_type), s.len))
            .collect();
        let end_idx = layout
            .iter()
            .position(|s| s.region_id == span.last().unwrap().region_id)
            .map(|i| i + 1)
            .unwrap_or(layout.len());
        let insert = layout
            .get(end_idx)
            .map(|s| format!("  ─▶ {}", s.region_type.to_ascii_lowercase()))
            .unwrap_or_default();
        Some((rid, len, format!("{}{}", boxes.join("─"), insert)))
    }

    let assays = embedded_assays();
    let mut out = String::new();
    out.push_str("# preflightx: resolved assay read geometries\n\n");
    out.push_str(
        "_Generated from the vendored seqspec catalogue (`vendor/seqspec`) by\n\
             `cargo test gen_assay_docs -- --ignored`. Do not edit by hand._\n\n\
             preflightx derives each assay's **observed barcode-read geometry** straight from\n\
             its seqspec description: it picks the barcode-bearing read (preferring the RNA\n\
             cell-barcode read for multimodal assays), strips leading library primers/adapters,\n\
             and spans the barcode / UMI / linker structure up to the cDNA insert. `[ CB n ]` is\n\
             a cell barcode, `[ UMI n ]` a UMI, `[ link n ]` a fixed inter-barcode linker, and\n\
             `─▶ cdna` marks where the biological insert begins.\n\n",
    );

    let featured = [
        "10x_rna_v3.spec.yaml",
        "10x_rna_5prime.spec.yaml",
        "drop_seq.spec.yaml",
        "split_seq.spec.yaml",
        "indropv2.spec.yaml",
        "bd_rhapsody_v1.spec.yaml",
        "scifi_rna_seq.spec.yaml",
        "10x_rna_atac_multi.spec.yaml",
    ];
    out.push_str("## Read structure (representative assays)\n\n```\n");
    for sf in featured {
        if let Some(a) = assays.iter().find(|a| a.spec_file == sf) {
            if let Some((rid, len, art)) = ascii(a) {
                let name = if a.name.is_empty() {
                    a.spec_file.trim_end_matches(".spec.yaml")
                } else {
                    &a.name
                };
                out.push_str(&format!("{name}   ({rid}, {len} bp)\n  {art}\n\n"));
            }
        }
    }
    out.push_str("```\n\n");

    let mut rows: Vec<(String, String)> = Vec::new();
    let (mut n_cov, mut n_plate) = (0, 0);
    for a in assays {
        let name = if a.name.is_empty() {
            a.spec_file.trim_end_matches(".spec.yaml").to_string()
        } else {
            a.name.clone()
        };
        let geom = a.structural_read().and_then(|(rid, l)| {
            let span = observed_span(&l);
            if span.is_empty() {
                return None;
            }
            let len: usize = span.iter().map(|s| s.len).sum();
            let cb: usize = span
                .iter()
                .filter(|s| is_barcode_type(&s.region_type))
                .map(|s| s.len)
                .sum();
            let umi: usize = span
                .iter()
                .filter(|s| is_umi_type(&s.region_type))
                .map(|s| s.len)
                .sum();
            Some(format!(
                "{rid}, {len} bp (CB {cb}{})",
                if umi > 0 {
                    format!(" + UMI {umi}")
                } else {
                    String::new()
                }
            ))
        });
        let cell = match geom {
            Some(g) => {
                n_cov += 1;
                g
            }
            None => {
                n_plate += 1;
                "_plate-based, no in-read cell barcode_".into()
            }
        };
        rows.push((
            a.spec_file.clone(),
            format!(
                "| {} | `{}` | {} | {} |",
                name,
                a.spec_file.trim_end_matches(".spec.yaml"),
                a.modalities.join(", "),
                cell
            ),
        ));
    }
    rows.sort();
    out.push_str(&format!(
        "## All {} catalogued assays ({} with barcode geometry, {} plate-based)\n\n",
        assays.len(),
        n_cov,
        n_plate
    ));
    out.push_str("| Assay | spec | modalities | barcode read |\n");
    out.push_str("|-------|------|------------|--------------|\n");
    for (_, r) in &rows {
        out.push_str(r);
        out.push('\n');
    }
    out.push('\n');

    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/docs/assay-geometries.md");
    std::fs::create_dir_all(concat!(env!("CARGO_MANIFEST_DIR"), "/docs")).unwrap();
    std::fs::write(path, out).unwrap();
    eprintln!("wrote {path}");
}

/// Regenerate web/assay-structures.html, a colour-coded visual of every
/// assay's barcode-read structure (also deployed at preflightx.saketlab.org).
/// Run on demand:
///   cargo test gen_assay_html -- --ignored
#[test]
#[ignore = "doc generator; run explicitly to refresh web/assay-structures.html"]
fn gen_assay_html() {
    fn category(rt: &str) -> &'static str {
        let l = rt.to_ascii_lowercase();
        if is_barcode_type(rt) {
            "cb"
        } else if is_umi_type(rt) {
            "umi"
        } else if is_linker_type(rt) {
            "linker"
        } else if l.starts_with("poly") || l.starts_with("ploy") {
            "poly"
        } else if matches!(
            l.as_str(),
            "cdna" | "gdna" | "rna" | "dna" | "insert" | "hic" | "methyl" | "atac" | "protein"
        ) {
            "insert"
        } else if l.starts_with("index")
            || l == "i5"
            || l == "i7"
            || l.ends_with("_index")
            || l == "sample_index"
        {
            "index"
        } else {
            "adapter"
        }
    }
    fn label(rt: &str) -> String {
        match category(rt) {
            "cb" => "BC".into(),
            "umi" => "UMI".into(),
            "linker" => "link".into(),
            "poly" => rt.to_ascii_lowercase().replace('_', ""),
            _ => rt.to_ascii_lowercase(),
        }
    }
    fn esc(s: &str) -> String {
        s.replace('&', "&amp;")
            .replace('<', "&lt;")
            .replace('>', "&gt;")
    }

    // Barcoding mechanism. Geometry alone can't separate droplet from single-barcode
    // well-plates, or microwell from split-pool, so it's curated; unmapped specs fall to Other.
    #[derive(PartialEq, Eq, Clone, Copy)]
    enum Mechanism {
        Droplet,
        Microwell,
        SplitPool,
        PlateBarcode,
        PlateNoBarcode,
        Other,
    }
    // Keyed by spec filename without the .spec.yaml suffix.
    fn mechanism_of(spec: &str) -> Mechanism {
        use Mechanism::*;
        match spec {
            "10x_atac" | "10xcrispr" | "10xfb_3prime" | "10xfb_5prime" | "10xfb_vdj_5prime"
            | "10x_rna_5prime" | "10x_rna_atac_multi" | "10x_rna_atac" | "10x_rna_v1"
            | "10x_rna_v2" | "10x_rna_v3" | "dogmaseq_dig" | "dogmaseq_lll" | "drop_seq"
            | "indropv2" | "pipseqv2" | "pipseqv3" | "sccrispra" | "sugarseq" | "surecell"
            | "vasa_seq_drop" => Droplet,
            "bd_rhapsody_eb" | "bd_rhapsody_v1" | "microwell_seq" | "seq_well" | "seq_well_s3" => {
                Microwell
            }
            "parse_wt_v2" | "scifi_rna_seq" | "sci_rna_seq" | "sci_rna_seq3" | "share_seq"
            | "split_seq" => SplitPool,
            "cel_seq" | "cel_seq2" | "mars_seq" | "mcscrb_seq" | "quartz_seq" | "quartz_seq2"
            | "strt_seq" | "strt_seq_2i" | "strt_seq_c1" | "sn_m3c_seq" | "snmctseq"
            | "vasa_seq_plate" | "issaac_seq" => PlateBarcode,
            "smart_seq2" | "smart_seq3" | "pi_atac_seq" => PlateNoBarcode,
            _ => Other,
        }
    }
    struct CardInfo {
        mechanism: Mechanism,
        year: Option<u16>,
        sort_len: usize,
        html: String,
    }

    // First plausible 4-digit year in a free-text date (formats vary upstream).
    fn year(date: &str) -> Option<u16> {
        let b = date.as_bytes();
        for w in b.windows(4) {
            if w.iter().all(u8::is_ascii_digit) {
                let y: u16 = std::str::from_utf8(w).unwrap().parse().unwrap();
                if (1990..=2100).contains(&y) {
                    return Some(y);
                }
            }
        }
        None
    }

    // Render one read's resolved layout into a track of segment blocks. Strips the
    // leading primer-binding handle (everything before the first sequenced
    // barcode/UMI/insert region) so each track shows only what is actually
    // sequenced. The variable cDNA/gDNA insert is drawn as a hatched terminal block
    // and absorbs the rest of the read. Returns (track_html, sequenced_bp), or None
    // if the read covers no sequenced content (a pure primer read).
    let read_track = |layout: &[Segment]| -> Option<(String, usize)> {
        let is_content = |s: &Segment| {
            is_barcode_type(&s.region_type)
                || is_umi_type(&s.region_type)
                || is_insert_type(&s.region_type)
        };
        let start = layout.iter().position(is_content)?;
        let mut track = String::new();
        let mut bp = 0usize;
        for s in &layout[start..] {
            let variable_insert =
                is_insert_type(&s.region_type) && (s.max_len > s.len || s.len == 0);
            if variable_insert {
                track.push_str(&format!(
                    "<div class=\"seg {cat} tail\" title=\"{rt} (insert)\">{lab} ▸</div>",
                    cat = category(&s.region_type),
                    rt = esc(&s.region_type),
                    lab = esc(&s.region_type.to_ascii_lowercase()),
                ));
                break;
            }
            bp += s.len;
            let px = (s.len.max(1) * 4).clamp(34, 200);
            track.push_str(&format!(
                    "<div class=\"seg {cat}\" style=\"flex:0 0 {px}px\" title=\"{rid_t} · {rt} · {n} bp\"><b>{lab}</b><i>{n}</i></div>",
                    cat = category(&s.region_type),
                    rid_t = esc(&s.region_id),
                    rt = esc(&s.region_type),
                    n = s.len,
                    lab = esc(&label(&s.region_type)),
                ));
        }
        Some((track, bp))
    };

    // Placeholder / template / mislabelled example specs upstream that aren't real
    // assays, excluded from the visual so the catalogue reads clean.
    const PLACEHOLDER_SPECS: &[&str] = &[
        "template.spec.yaml",
        "tang2009.spec.yaml", // upstream name: is the bogus "10x Genomics"
        "10xv3_scrnaseq_element_adept_truseq_dual.spec.yaml", // name: "Example assay"
        "10xv3_scrnaseq_illumina_truseq_dual.spec.yaml", // name: "Example assay"
    ];
    let assays: Vec<Assay> = embedded_assays()
        .iter()
        .filter(|a| !PLACEHOLDER_SPECS.contains(&a.spec_file.as_str()))
        .cloned()
        .collect();

    let card = |a: &Assay| -> Option<CardInfo> {
        let bc = a.structural_read();
        let bc_read_id = bc.as_ref().map(|(id, _)| id.as_str());
        let bc_len = bc.as_ref().and_then(|(_, l)| observed_read_len(l));

        let mut rows = String::new();
        let mut n_rows = 0;
        for read in &a.sequence_spec {
            let Some(layout) = a.read_layout(&read.read_id) else {
                continue;
            };
            let Some((track, bp)) = read_track(&layout) else {
                continue;
            };
            let rid = read
                .read_id
                .trim_end_matches(".gz")
                .trim_end_matches(".fastq")
                .trim_end_matches(".fq");
            let is_bc = bc_read_id == Some(read.read_id.as_str());
            let badge = if is_bc {
                "<span class=\"bcflag\">barcode</span>".to_string()
            } else {
                String::new()
            };
            let modtag = if read.modality.is_empty() {
                String::new()
            } else {
                format!("<span class=\"rmod\">{}</span>", esc(&read.modality))
            };
            let bptag = if bp > 0 {
                format!("<span class=\"rbp\">{bp} bp</span>")
            } else {
                String::new()
            };
            rows.push_str(&format!(
                    "<div class=\"readrow{bccls}\">\
                       <div class=\"rlabel\">{rid}{badge}{modtag}{bptag}</div>\
                       <div class=\"readline\"><span class=\"end\">5′</span><div class=\"track\">{track}</div><span class=\"end\">3′</span></div>\
                     </div>",
                    bccls = if is_bc { " is-bc" } else { "" },
                    rid = esc(rid),
                ));
            n_rows += 1;
        }
        if n_rows == 0 {
            return None;
        }

        let name = if a.name.is_empty() {
            a.spec_file.trim_end_matches(".spec.yaml")
        } else {
            &a.name
        };
        let mods = a
            .modalities
            .iter()
            .map(|m| format!("<span class=\"mod\">{}</span>", esc(m)))
            .collect::<String>();
        let doi = a.doi.trim();
        let doi_link = if doi.is_empty() {
            String::new()
        } else {
            let href = if doi.starts_with("http") {
                doi.to_string()
            } else {
                format!("https://doi.org/{}", doi.trim_start_matches("doi:").trim())
            };
            format!(
                "<a class=\"doi\" href=\"{href}\" target=\"_blank\" rel=\"noopener\">DOI ↗</a>",
                href = esc(&href)
            )
        };
        let mechanism = mechanism_of(a.spec_file.trim_end_matches(".spec.yaml"));
        let yr = year(&a.date);
        let yr_badge = yr
            .map(|y| format!("<span class=\"yr\">{y}</span>"))
            .unwrap_or_default();
        let cap = match (bc_read_id, bc_len) {
            (Some(rid), Some(len)) => format!(
                "{rid} · <strong>{len} bp</strong> barcode read",
                rid = esc(rid)
            ),
            _ => "plate-based · cell identity is the well, no in-read barcode".to_string(),
        };
        let html = format!(
                "<article class=\"card{platecls}\">\
                   <header><h3>{name}</h3><div class=\"meta\">{yr_badge}{mods}<code>{spec}</code>{doi_link}</div></header>\
                   <div class=\"reads\">{rows}</div>\
                   <div class=\"cap\">{cap}</div>\
                 </article>",
                platecls = if bc_read_id.is_none() { " plate-card" } else { "" },
                name = esc(name),
                spec = esc(a.spec_file.trim_end_matches(".spec.yaml")),
            );
        Some(CardInfo {
            mechanism,
            year: yr,
            sort_len: bc_len.unwrap_or(0),
            html,
        })
    };

    let all_cards: Vec<CardInfo> = assays.iter().filter_map(card).collect();
    // Cards of one mechanism, ordered chronologically (undated last), then by barcode length.
    let section = |m: Mechanism| -> Vec<String> {
        let mut v: Vec<&CardInfo> = all_cards.iter().filter(|c| c.mechanism == m).collect();
        v.sort_by(|a, b| {
            a.year
                .unwrap_or(u16::MAX)
                .cmp(&b.year.unwrap_or(u16::MAX))
                .then_with(|| a.sort_len.cmp(&b.sort_len))
                .then_with(|| a.html.cmp(&b.html))
        });
        v.into_iter().map(|c| c.html.clone()).collect()
    };
    let droplet_cards = section(Mechanism::Droplet);
    let microwell_cards = section(Mechanism::Microwell);
    let splitpool_cards = section(Mechanism::SplitPool);
    let platebc_cards = section(Mechanism::PlateBarcode);
    let plate_cards = section(Mechanism::PlateNoBarcode);
    let other_cards = section(Mechanism::Other);
    // Any assay whose reads don't render at all (no barcode/UMI/insert content),
    // listed by name so the catalogue stays complete.
    let plate: Vec<String> = assays
        .iter()
        .filter(|a| card(a).is_none())
        .map(|a| {
            esc(if a.name.is_empty() {
                a.spec_file.trim_end_matches(".spec.yaml")
            } else {
                &a.name
            })
        })
        .collect();

    // Palette/fonts mirror web/index.html; self-contained, light + dark.
    let css = r#"
*{box-sizing:border-box}
:root{
  --bg:#f4f1e8; --panel:#fcfbf6; --panel2:#efebdf; --ink:#191a26; --muted:#6b6a5f; --faint:#6b6a5f;
  --hair:#ddd8c8; --hair2:#e6e1d2;
  --cb:#3a52c6; --umi:#0a7d5f; --linker:#5b6472; --insert:#c4420e;
  --poly:#e0a124; --index:#8e44ad; --adapter:#d3ccb9;
  --on-accent:#fcfbf6; --shadow:0 1px 2px rgba(25,26,38,.05);
  --ui:-apple-system,BlinkMacSystemFont,'SF Pro Text','SF Pro Display',system-ui,sans-serif;
  --mono:ui-monospace,'SF Mono',SFMono-Regular,Menlo,monospace}
@media (prefers-color-scheme:dark){:root{
  --bg:#16150f; --panel:#1e1d15; --panel2:#25241a; --ink:#ece9dd; --muted:#9a978a; --faint:#8f8b7d;
  --hair:#34322a; --hair2:#2b2a22;
  --cb:#6b80ff; --umi:#2bb98f; --linker:#8b93a3; --insert:#ff7a45;
  --poly:#f0bd55; --index:#b974d6; --adapter:#6a6553;
  --on-accent:#16150f; --shadow:0 1px 2px rgba(0,0,0,.35)}}
html{-webkit-font-smoothing:antialiased;text-rendering:optimizeLegibility}
body{margin:0;color:var(--ink);background:var(--bg);font-family:var(--ui);font-size:15px;line-height:1.5}
.wrap{max-width:1000px;margin:0 auto;padding:64px 24px 96px}
.title{font-weight:600;font-size:clamp(30px,5vw,46px);line-height:1.06;letter-spacing:-.02em;margin:0 0 16px}
.title em{font-style:normal;color:var(--cb)}
.lede{max-width:68ch;color:var(--muted);font-size:15px;line-height:1.6;margin:0 0 34px}
.lede a{color:var(--cb)}
.legend{display:flex;flex-wrap:wrap;align-items:center;gap:10px 18px;padding:12px 16px;margin:0 0 40px;
  background:var(--panel);border:1px solid var(--hair);border-radius:8px;position:sticky;top:12px;z-index:5}
.legend b{font-weight:600;margin-right:4px;font-size:12px;color:var(--muted);
  text-transform:uppercase;letter-spacing:.06em}
.key{display:inline-flex;align-items:center;gap:7px;font-size:13px;color:var(--ink)}
.sw{width:14px;height:14px;border-radius:3px}
.card{background:var(--panel);border:1px solid var(--hair);border-radius:8px;
  padding:18px 22px;margin:0 0 14px;box-shadow:var(--shadow);
  opacity:0;transform:translateY(8px);animation:rise .4s ease forwards}
@keyframes rise{to{opacity:1;transform:none}}
.card.plate-card{border-style:dashed;background:var(--panel2)}
.card.plate-card .cap{color:var(--faint)}
.card header{display:flex;justify-content:space-between;align-items:baseline;
  gap:14px;flex-wrap:wrap;margin-bottom:16px}
.card h3{font-weight:600;font-size:18px;margin:0;letter-spacing:-.01em}
.meta{display:flex;align-items:center;gap:8px;flex-wrap:wrap}
.mod{font-size:11px;letter-spacing:.02em;color:var(--muted);
  background:var(--panel2);border-radius:999px;padding:2px 10px}
.yr{font-size:11px;color:var(--muted);background:var(--panel2);border-radius:999px;
  padding:2px 9px;font-weight:600;font-variant-numeric:tabular-nums}
.meta code{font-size:12px;color:var(--muted);font-family:var(--mono)}
.doi{font-size:11px;color:var(--cb);text-decoration:none;font-weight:600;
  border:1px solid var(--hair);border-radius:999px;padding:4px 10px;
  display:inline-flex;align-items:center;min-height:24px}
.doi:hover{background:var(--cb);color:var(--on-accent);border-color:var(--cb)}
.reads{display:flex;flex-direction:column;gap:10px}
.readrow{display:grid;grid-template-columns:124px 1fr;align-items:center;gap:16px}
.rlabel{display:flex;align-items:center;flex-wrap:wrap;gap:6px;font-size:12px;color:var(--muted);
  font-weight:600;word-break:break-all}
.readrow.is-bc .rlabel{color:var(--ink)}
.bcflag{font-size:10px;letter-spacing:.01em;color:var(--on-accent);
  background:var(--cb);border-radius:999px;padding:1px 8px;font-weight:600}
.rmod{font-size:10px;letter-spacing:.02em;color:var(--muted);
  background:var(--panel2);border-radius:999px;padding:1px 8px}
.rbp{font-size:11px;color:var(--muted);font-weight:500}
.readline{display:flex;align-items:center;gap:12px;min-width:0}
.end{color:var(--muted);font-size:12px;font-weight:600}
.track{display:flex;align-items:stretch;gap:3px;flex:1;overflow-x:auto;padding:1px 0}
.seg{flex:0 0 auto;min-height:30px;border-radius:4px;display:flex;flex-direction:row;
  align-items:center;justify-content:center;gap:5px;color:var(--on-accent);position:relative;padding:0 8px;cursor:default}
.seg b{font-weight:600;font-size:11px;letter-spacing:.01em}
.seg i{font-style:normal;font-size:10px;opacity:.85}
.seg.cb{background:var(--cb)} .seg.umi{background:var(--umi)}
.seg.insert{background:var(--insert)} .seg.index{background:var(--index)}
.seg.poly{background:var(--poly);color:#3a2c00} .seg.poly i{opacity:.7}
.seg.linker{background:var(--linker)} .seg.linker i{opacity:.8}
.seg.adapter{background:var(--adapter);color:var(--ink)} .seg.adapter i{opacity:.7}
.seg.tail{flex:0 0 auto;padding:0 14px;min-width:74px;color:var(--insert);font-weight:600;font-size:11px;
  background:transparent;border:1.5px dashed var(--insert)}
.cap{margin-top:14px;font-size:12px;color:var(--muted)}
.cap strong{color:var(--ink);font-weight:600}
.section{font-size:12px;letter-spacing:.06em;text-transform:uppercase;font-weight:600;
  color:var(--muted);margin:48px 0 18px;border-top:1px solid var(--hair);padding-top:24px}
.shint{text-transform:none;letter-spacing:0;font-weight:400;color:var(--faint);font-size:12px;margin-left:8px}
.plate{display:flex;flex-wrap:wrap;gap:8px}
.plate span{font-size:13px;color:var(--muted);background:var(--panel2);border:1px solid var(--hair);
  border-radius:6px;padding:6px 12px}
.foot{margin-top:56px;color:var(--muted);font-size:12px;line-height:1.7;
  border-top:1px solid var(--hair);padding-top:20px}
@media (max-width:480px){.readrow{grid-template-columns:1fr;gap:6px}.wrap{padding:48px 16px 72px}}
@media (prefers-reduced-motion:reduce){*{animation:none!important;transition:none!important}}
"#;

    let legend = [
            ("cb", "cell barcode"),
            ("umi", "UMI"),
            ("linker", "linker / spacer"),
            ("poly", "polyA / polyT"),
            ("insert", "cDNA / gDNA insert"),
            ("index", "sample index"),
            ("adapter", "primer / adapter"),
        ]
        .iter()
        .map(|(c, name)| {
            format!(
                "<span class=\"key\"><span class=\"sw\" style=\"background:var(--{c})\"></span>{name}</span>"
            )
        })
        .collect::<String>();

    let mut h = String::new();
    h.push_str("<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">");
    h.push_str("<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">");
    h.push_str("<title>preflightx · single-cell read structures</title>");
    h.push_str(&format!(
        "<style>{css}</style></head><body><div class=\"wrap\">"
    ));
    h.push_str("<h1 class=\"title\">Single-cell<br>read <em>structures</em></h1>");
    h.push_str(&format!(
            "<p class=\"lede\">Every assay in the preflightx catalogue, drawn straight from its \
             <a href=\"https://github.com/pachterlab/seqspec\">seqspec</a> description, grouped by how each \
             isolates a cell. Each card shows <strong>every sequencing read</strong> (R1, R2, …), left to \
             right, 5′→3′, after stripping the library primer that isn't sequenced; the read carrying the \
             cell barcode is flagged. {} assays across {} mechanisms, oldest first within each. Hover a block \
             for its region id.</p>",
            droplet_cards.len()
                + microwell_cards.len()
                + splitpool_cards.len()
                + platebc_cards.len()
                + plate_cards.len()
                + plate.len()
                + other_cards.len(),
            5,
        ));
    h.push_str(&format!(
        "<div class=\"legend\"><b>legend</b>{legend}</div>"
    ));

    // Stagger the entrance animation; the style is injected on the <article> tag so it
    // is independent of the card's class. Empty sections print no header.
    fn emit(h: &mut String, title: &str, sub: &str, cards: &[String]) {
        if cards.is_empty() {
            return;
        }
        h.push_str(&format!(
            "<h2 class=\"section\">{title} <span class=\"shint\">{sub}</span></h2>"
        ));
        for (i, c) in cards.iter().enumerate() {
            let delay = (i.min(24) * 40) as u32;
            h.push_str(&c.replacen(
                "<article ",
                &format!("<article style=\"animation-delay:{delay}ms\" "),
                1,
            ));
        }
    }
    emit(
        &mut h,
        "Droplet",
        "cells in droplets · one barcode block",
        &droplet_cards,
    );
    emit(
        &mut h,
        "Microwell",
        "cells in nanowells · segmented barcode",
        &microwell_cards,
    );
    emit(
        &mut h,
        "Split-pool",
        "combinatorial indexing · barcode in blocks",
        &splitpool_cards,
    );
    emit(
        &mut h,
        "Plate / well",
        "one cell per well · in-read well barcode",
        &platebc_cards,
    );
    emit(
        &mut h,
        "Plate",
        "one cell per well · no in-read barcode",
        &plate_cards,
    );
    emit(
        &mut h,
        "Other",
        "unmapped, add to the mechanism table",
        &other_cards,
    );
    if !plate.is_empty() {
        h.push_str("<div class=\"plate\">");
        for p in &plate {
            h.push_str(&format!("<span>{p}</span>"));
        }
        h.push_str("</div>");
    }

    h.push_str(
            "<p class=\"foot\">Generated from the vendored seqspec catalogue by \
             <code>cargo test gen_assay_html -- --ignored</code>. Block widths are proportional to \
             length (clamped for very long/short segments); leading primer-binding regions are omitted \
             because they are not sequenced. Geometries follow each spec faithfully.</p>",
        );
    h.push_str("</div></body></html>");

    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/web/assay-structures.html");
    std::fs::write(path, h).unwrap();
    eprintln!("wrote {path}");
}

#[test]
fn multimodal_prefers_rna_barcode_read() {
    // 10x multiome has an RNA read (CB16+UMI12=28) and an ATAC read (CB16 only);
    // structural_read must prefer the RNA cell-barcode read.
    assert_eq!(geom_of("10x_rna_atac_multi.spec.yaml"), Some(28));
    assert_eq!(geom_of("dogmaseq_lll.spec.yaml"), Some(28));
}

#[test]
fn trailing_adapter_dropped_but_inter_barcode_adapter_kept() {
    // quartz-seq's read-2 adapter sits AFTER the cell barcode (loose primer_id) and
    // must be dropped -> 6 bp cell index. scifi-RNA-seq's adapters sit BETWEEN the
    // bead and RT barcodes and are sequenced -> kept in the 84 bp span.
    assert_eq!(geom_of("quartz_seq.spec.yaml"), Some(6));
    assert_eq!(geom_of("scifi_rna_seq.spec.yaml"), Some(84));
}

#[test]
fn every_spec_is_covered_or_plate_based() {
    // The completeness guarantee: each catalogued assay either yields a barcode-read
    // geometry, or has no in-read cell barcode at all (plate-based). Nothing should
    // have a barcode region we fail to resolve into a geometry.
    let mut failed: Vec<String> = Vec::new();
    let mut covered = 0;
    let mut plate = 0;
    for a in embedded_assays() {
        let has_bc_region = a
            .leaf_regions()
            .iter()
            .any(|r| is_barcode_type(&r.region_type));
        let geom = a.structural_read().and_then(|(_, l)| observed_read_len(&l));
        match (geom, has_bc_region) {
            (Some(_), _) => covered += 1,
            (None, false) => plate += 1, // genuinely plate-based, correct
            (None, true) => failed.push(a.spec_file.clone()),
        }
    }
    assert!(
        failed.is_empty(),
        "{} spec(s) have a barcode region but no resolved geometry: {failed:?}",
        failed.len()
    );
    // Sanity: the catalogue is mostly droplet/combinatorial (in-read barcodes), with
    // a handful of plate-based assays.
    assert!(covered >= 45, "only {covered} specs covered");
    assert!(
        (1..=12).contains(&plate),
        "unexpected plate-based count {plate}"
    );
}

#[test]
fn variable_and_synonym_barcodes_are_covered() {
    // inDrop's variable-length BC1 and BD Rhapsody's combinatorial CLS barcodes
    // (with a leading variable spacer) both resolve, exercising variable-length
    // barcodes, synonym region_types, and the linker-aware terminator.
    let g = sc_geometries();
    for spec in [
        "indropv2.spec.yaml",
        "bd_rhapsody_eb.spec.yaml",
        "split_seq.spec.yaml",
    ] {
        assert!(
            g.iter().any(|x| x.spec_file == spec),
            "{spec} should have a geometry"
        );
    }
}

#[test]
fn observed_geometry_strips_leading_primer() {
    // Drop-seq R1's structural read includes the 39 bp Read1 custom primer before
    // the barcode (seqspec models it as a region after primer_id=illumina_p5);
    // the observed barcode read is the 20 bp (CB12+UMI8) after it.
    let drop = sc_geometries()
        .iter()
        .find(|g| g.spec_file == "drop_seq.spec.yaml")
        .expect("drop_seq present");
    assert_eq!(
        drop.read_len, 20,
        "drop_seq observed read should be CB12+UMI8"
    );
    // 10x v3: clean R1, no leading primer -> 16 + 12 = 28.
    let v3 = sc_geometries()
        .iter()
        .find(|g| g.spec_file == "10x_rna_v3.spec.yaml")
        .expect("v3 present");
    assert_eq!(v3.read_len, 28);
}

#[test]
fn sc_geometries_matching_finds_known_lengths() {
    assert!(sc_geometries_matching(20, 1)
        .iter()
        .any(|g| g.spec_file == "drop_seq.spec.yaml"));
    assert!(sc_geometries_matching(28, 1)
        .iter()
        .any(|g| g.spec_file == "10x_rna_v3.spec.yaml"));
    assert!(sc_geometries_matching(0, 1).is_empty());
    assert!(sc_geometries_matching(1000, 1).is_empty());
}

#[test]
fn parses_10x_v3_and_resolves_r1_layout() {
    let yaml = include_str!("../../vendor/seqspec/docs/examples/assays/10x_rna_v3.spec.yaml");
    let assay = Assay::from_yaml(yaml).unwrap();
    assert_eq!(assay.assay_id, "10xRNAv3");
    assert_eq!(assay.modalities, vec!["rna"]);

    // R1 primes from truseq_read1 (pos strand), covering CB(16) then UMI(12).
    let layout = assay.read_layout("rna-R1.fastq.gz").unwrap();
    let kinds: Vec<(&str, usize)> = layout
        .iter()
        .map(|s| (s.region_type.as_str(), s.len))
        .collect();
    assert_eq!(kinds, vec![("barcode", 16), ("umi", 12)]);
}

#[test]
fn variable_cdna_region_terminates_the_layout() {
    // SPLiT-seq R1 primes from s5 and covers ME1 then cDNA. cDNA is variable
    // (max_len > min_len) and absorbs the rest of the read, so NO barcode
    // segments may leak into the cDNA read's layout.
    let yaml = include_str!("../../vendor/seqspec/docs/examples/assays/split_seq.spec.yaml");
    let assay = Assay::from_yaml(yaml).unwrap();
    let layout = assay.read_layout("R1.fastq.gz").unwrap();
    assert!(
        !layout.iter().any(|s| s.region_type == "barcode"),
        "cDNA read must not contain barcode segments, got {:?}",
        layout.iter().map(|s| &s.region_type).collect::<Vec<_>>()
    );
    // cDNA is the terminal segment.
    assert_eq!(layout.last().unwrap().region_type, "cDNA");
}

#[test]
fn embedded_catalogue_loads_many_assays() {
    let assays = embedded_assays();
    // ~50 vendored specs; every one must parse cleanly.
    assert_eq!(assays.len(), EMBEDDED_SPECS.len());
    assert!(assays.len() >= 40, "only {} assays parsed", assays.len());
    assert!(assays
        .iter()
        .any(|a| a.assay_id.to_lowercase().contains("10x")));
    assert!(assays
        .iter()
        .any(|a| a.assay_id.to_lowercase().contains("split")
            || a.name.to_lowercase().contains("split")));
}

#[test]
fn resolves_combinatorial_neg_strand_layout_splitseq() {
    // SPLiT-seq R2 primes from truseq_read2 on the neg strand and reads backward:
    // UMI(10) + barcode-3(8) + linker(30) + barcode-2(8) + linker(30) + barcode-1(8) = 94 bp.
    let yaml = include_str!("../../vendor/seqspec/docs/examples/assays/split_seq.spec.yaml");
    let assay = Assay::from_yaml(yaml).unwrap();
    let layout = assay.read_layout("R2.fastq.gz").unwrap();
    let kinds: Vec<(&str, usize)> = layout
        .iter()
        .map(|s| (s.region_type.as_str(), s.len))
        .collect();
    assert_eq!(
        kinds,
        vec![
            ("umi", 10),
            ("barcode", 8),
            ("linker", 30),
            ("barcode", 8),
            ("linker", 30),
            ("barcode", 8),
        ]
    );
    // Linkers carry a fixed sequence (reverse-complemented for the neg strand).
    let linker = layout.iter().find(|s| s.region_type == "linker").unwrap();
    assert!(linker.sequence.is_some());
}

#[test]
fn scatac_structural_reads_are_modality_aware() {
    let assays = embedded_assays();
    let by = |f: &str| assays.iter().find(|a| a.spec_file == f).unwrap();

    let (_rid, atac) = by("10x_atac.spec.yaml")
        .structural_read()
        .expect("10x_atac resolves to its barcode read");
    let bc: Vec<usize> = atac
        .iter()
        .filter(|s| is_barcode_type(&s.region_type))
        .map(|s| s.len)
        .collect();
    assert_eq!(
        bc,
        [16],
        "10x scATAC has one 16 bp cell barcode, got {bc:?}"
    );
    assert!(
        !atac.iter().any(|s| is_umi_type(&s.region_type)),
        "scATAC has no UMI"
    );

    let (rid, rna) = by("10x_rna_atac_multi.spec.yaml")
        .structural_read()
        .expect("multiome resolves to its RNA read");
    assert!(
        rid.starts_with("rna"),
        "multiome must resolve to the RNA read, got {rid}"
    );
    let rbc: Vec<usize> = rna
        .iter()
        .filter(|s| is_barcode_type(&s.region_type))
        .map(|s| s.len)
        .collect();
    let umi: Vec<usize> = rna
        .iter()
        .filter(|s| is_umi_type(&s.region_type))
        .map(|s| s.len)
        .collect();
    assert_eq!(rbc, [16], "10x RNA cell barcode is 16 bp, got {rbc:?}");
    assert_eq!(umi, [12], "10x multiome RNA UMI is 12 bp, got {umi:?}");
}
