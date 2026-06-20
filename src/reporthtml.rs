//! Self-contained HTML report renderer, shared by the CLI (`check --html`) and the
//! browser (`wasm::report_html_json`, the page's "Export report" button).
//!
//! Renders from a `serde_json::Value`, not the typed [`crate::report::FullReport`], so the
//! same renderer serves both the CLI's typed report and the browser's client-merged object.
//! A bare object renders one report; an array renders one `<section>` per element.

use serde_json::Value;

/// Render one report object, or an array of them, into a complete HTML document.
pub fn render(v: &Value) -> String {
    let (body, sub) = match v {
        Value::Array(items) => (
            items.iter().map(section).collect::<String>(),
            format!("{} report{}", items.len(), if items.len() == 1 { "" } else { "s" }),
        ),
        _ => (section(v), label(v).unwrap_or_else(|| "single report".into())),
    };
    document(&body, &sub)
}

/// One `<section>` per report. Fields are pulled by name, so this is a third surface to keep
/// in sync with the report structs (alongside `--json` and each `print()`): a missing or
/// typo'd key renders nothing, with no compile error.
fn section(r: &Value) -> String {
    let chem = r.get("chemistry").filter(|c| !c.is_null());
    let single_end = chem.is_none();
    let not_sc = chem.and_then(|c| c.get("is_single_cell")) == Some(&Value::Bool(false));

    let (head, conf, warn) = if single_end {
        ("Single-end / bulk input".to_string(), String::new(), false)
    } else if not_sc {
        ("Not single-cell".to_string(), String::new(), true)
    } else {
        let c = chem.unwrap();
        let head = str_of(c, "assay_display")
            .or_else(|| str_of(c, "assay_id"))
            .unwrap_or_else(|| "unknown".into());
        (head, str_of(c, "confidence").as_deref().map(pill).unwrap_or_default(), false)
    };

    let label_html = label(r)
        .map(|l| format!("<div class=\"rlabel\">{}</div>", esc(&l)))
        .unwrap_or_default();

    format!(
        "<section class=\"report\">{label_html}\
         <div class=\"vmain{}\">{}</div>{conf}\
         <div class=\"grid\">{}</div>{}{}{}</section>",
        if warn { " warn" } else { "" },
        esc(&head),
        stat_grid(r, chem, not_sc),
        track_block(if not_sc { None } else { chem }),
        sex_block(r.get("sex")),
        assay_block(r.get("assay_detection")),
    )
}

/// The R1 read-structure diagram (adapter / cell-barcode / UMI segments), mirroring the
/// page's renderTrack. Only the barcode read is drawn; the cDNA mate is read 2.
fn track_block(chem: Option<&Value>) -> String {
    let Some(c) = chem else { return String::new() };
    let Some(cb) = num_of(c, "cb_len").filter(|&n| n > 0.0) else {
        return String::new();
    };
    let mut segs = vec![];
    let off = num_of(c, "cb_offset").unwrap_or(0.0);
    if off > 0.0 {
        segs.push((off, "adapter", "adapter".to_string()));
    }
    segs.push((cb, "cb", format!("BC {}", fnum(cb))));
    if let Some(umi) = num_of(c, "umi_len").filter(|&n| n > 0.0) {
        segs.push((umi, "umi", format!("UMI {}", fnum(umi))));
    }
    let total: f64 = segs.iter().map(|s| s.0).sum();
    let blocks: String = segs
        .iter()
        .map(|(len, cls, lbl)| {
            let n = fnum(*len);
            format!(
                "<div class=\"seg seg-{cls}\" style=\"flex:{n}\" title=\"{n} bp\">{}</div>",
                esc(lbl),
            )
        })
        .collect();
    format!(
        "<div class=\"track-lbl\">R1 read structure ({} bp barcode read) — cell barcode (CB) \
         identifies the cell, UMI tags the molecule; cDNA is read 2</div>\
         <div class=\"track\">{blocks}</div>",
        fnum(total),
    )
}

/// The report's stat grid; each scalar is emitted only when present. Add a report field
/// here to surface it (see the sync note on `section`).
fn stat_grid(r: &Value, chem: Option<&Value>, not_sc: bool) -> String {
    let mut g = Vec::new();
    if !not_sc {
        if let Some(c) = chem {
            if let Some(cb) = num_of(c, "cb_len") {
                g.push(stat(
                    "barcode geometry",
                    &format!(
                        "CB {} + UMI {} ({})",
                        fnum(cb),
                        num_of(c, "umi_len").map(fnum).unwrap_or_default(),
                        str_of(c, "cb_read").unwrap_or_default()
                    ),
                ));
            }
            if let Some(d) = num_of(c, "r1_length_dominant") {
                g.push(stat("R1 dominant len", &format!("{} bp", fnum(d))));
            }
            if let Some(n) = num_of(c, "reads_sampled") {
                g.push(stat("reads used", &group(n)));
            }
        }
    }
    if let Some(v) = r.get("bulk_vs_sc").and_then(|b| str_of(b, "verdict")) {
        g.push(stat("library", &v));
    }
    if let Some(p) = r.get("size").and_then(|s| s.get("r1")).and_then(|f| num_of(f, "mean_phred")) {
        g.push(stat("mean PHRED", &fnum(p)));
    }
    if let Some(c) = r.get("composition") {
        if let Some(gc) = num_of(c, "gc_pct") {
            g.push(stat("GC %", &fnum(gc)));
        }
        if let Some(o) = str_of(c, "guess_organism_family") {
            g.push(stat("organism hint", &o));
        }
    }
    if let Some(s) = r.get("size") {
        if let Some(a) = s.get("adequate_for_sc").and_then(Value::as_bool) {
            g.push(stat("adequate for SC", if a { "yes" } else { "no" }));
        }
    }
    if let Some(sp) = r.get("species").filter(|s| !s.is_null()) {
        if let Some(c) = str_of(sp, "called") {
            g.push(stat("species (reads)", &c));
        }
    }
    g.join("")
}

/// The opt-in genetic-sex panel, when the report carries one.
fn sex_block(sex: Option<&Value>) -> String {
    let Some(s) = sex.filter(|s| !s.is_null()) else {
        return String::new();
    };
    let verdict = str_of(s, "verdict").unwrap_or_else(|| "unknown".into());
    let cls = match verdict.as_str() {
        "male" | "female" | "uncertain" => verdict.as_str(),
        _ => "unknown",
    };
    let mut g = Vec::new();
    if let Some(sp) = str_of(s, "species") {
        let st = str_of(s, "panel_status").map(|p| format!(" ({p})")).unwrap_or_default();
        g.push(stat("species", &format!("{}{}", esc(&sp), esc(&st))));
    }
    for (key, lbl) in [
        ("y_hit_reads", "Y-panel hits"),
        ("xist_hit_reads", "XIST hits"),
        ("y_xist_ratio", "Y/XIST ratio"),
        ("reads_scanned", "reads scanned"),
    ] {
        if let Some(n) = num_of(s, key) {
            let val = if key == "reads_scanned" { group(n) } else { fnum(n) };
            g.push(stat(lbl, &val));
        }
    }
    let lead = format!(
        "<div><span class=\"sexverdict {cls}\">{}</span> {}</div>",
        esc(&verdict),
        str_of(s, "confidence").as_deref().map(pill).unwrap_or_default(),
    );
    detail_block("Genetic sex", "· opt-in, sample-level signal", &lead, &g.join(""))
}

/// Shared shell for the genetic-sex and seqspec-assay sub-sections.
fn detail_block(title: &str, sub: &str, lead: &str, stats: &str) -> String {
    format!(
        "<div class=\"sexblock\"><h3>{title} <span class=\"muted\">{sub}</span></h3>\
         {lead}<div class=\"grid\">{stats}</div></div>"
    )
}

/// The seqspec assay fallback's top hit, shown only when chemistry found no 10x match.
fn assay_block(assay: Option<&Value>) -> String {
    let Some(a) = assay.filter(|a| !a.is_null()) else {
        return String::new();
    };
    let top = a.get("hits").and_then(|h| h.get(0));
    let name = top
        .and_then(|t| str_of(t, "name").or_else(|| str_of(t, "assay_id")).or_else(|| str_of(t, "spec_file")))
        .unwrap_or_else(|| "no match".into());
    let mut g = Vec::new();
    if let Some(t) = top {
        if let Some(o) = str_of(t, "orientation") {
            g.push(stat("orientation", &o));
        }
        if let (Some(e), Some(ob)) = (num_of(t, "expected_len"), num_of(t, "observed_len")) {
            g.push(stat("expected / observed", &format!("{} / {} bp", fnum(e), fnum(ob))));
        }
        if let (Some(f), Some(tot)) = (num_of(t, "anchors_found"), num_of(t, "anchors_total")) {
            g.push(stat("anchors", &format!("{}/{}", fnum(f), fnum(tot))));
        }
    }
    let lead = format!(
        "<div class=\"vmain\">{}</div>{}",
        esc(&name),
        a.get("confidence").and_then(Value::as_str).map(pill).unwrap_or_default(),
    );
    detail_block("Best seqspec assay", "· geometry + linker anchors", &lead, &g.join(""))
}

fn label(r: &Value) -> Option<String> {
    str_of(r, "_label").or_else(|| str_of(r, "r1"))
}

fn stat(k: &str, v: &str) -> String {
    format!("<div class=\"stat\"><div class=\"k\">{}</div><div class=\"v\">{}</div></div>", esc(k), v)
}

fn pill(conf: &str) -> String {
    let cls = match conf {
        "high" | "medium" | "low" => conf,
        _ => "low",
    };
    format!("<span class=\"pill {cls}\">{}</span>", esc(conf))
}

/// A field's string value, or a number stringified.
fn str_of(v: &Value, key: &str) -> Option<String> {
    match v.get(key)? {
        Value::String(s) if !s.is_empty() => Some(s.clone()),
        Value::Number(n) => Some(n.to_string()),
        _ => None,
    }
}

fn num_of(v: &Value, key: &str) -> Option<f64> {
    v.get(key)?.as_f64()
}

/// Trim a float for display: integers print bare, others to one decimal.
fn fnum(n: f64) -> String {
    if n.fract().abs() < 1e-9 {
        format!("{}", n as i64)
    } else {
        format!("{n:.2}").trim_end_matches('0').trim_end_matches('.').to_string()
    }
}

/// Thousands-separated integer (reads counts).
fn group(n: f64) -> String {
    let n = n as i64;
    let s = n.abs().to_string();
    let mut out = String::new();
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(c);
    }
    if n < 0 {
        format!("-{out}")
    } else {
        out
    }
}

fn esc(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn document(body: &str, sub: &str) -> String {
    format!(
        "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">\
<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\
<title>preflightx report</title><style>{CSS}</style></head><body>\
<header><h1>preflight<span class=\"x\">x</span> report</h1>\
<div class=\"sub\">{}</div></header><main>{body}</main>\
<footer>Generated by <a href=\"https://github.com/saketlab/preflightx\">preflightx</a> \
— sampled the first reads only.</footer></body></html>",
        esc(sub)
    )
}

const CSS: &str = "\
:root{--bg:#f4f1e8;--panel:#fcfbf6;--panel2:#efebdf;--line:#ddd8c8;--fg:#191a26;--muted:#6b6a5f;\
--paper-on-accent:#fcfbf6;--accent:#3a52c6;--good:#1f7a4d;--warn:#9a6b00;--bad:#c0392b;--umi:#0a7d5f;\
--shadow:0 1px 2px rgba(25,26,38,.05);--ui:-apple-system,BlinkMacSystemFont,'SF Pro Text','SF Pro Display',system-ui,sans-serif;\
--mono:ui-monospace,'SF Mono',SFMono-Regular,Menlo,monospace}\
@media(prefers-color-scheme:dark){:root{--bg:#16150f;--panel:#1e1d15;--panel2:#25241a;--line:#34322a;\
--fg:#ece9dd;--muted:#9a978a;--paper-on-accent:#16150f;--accent:#6b80ff;--good:#45c98a;--warn:#f0bd55;\
--bad:#ff7a6b;--umi:#2bb98f;--shadow:0 1px 2px rgba(0,0,0,.35)}}\
*{box-sizing:border-box}\
body{margin:0;background:var(--bg);color:var(--fg);font:14.5px/1.55 var(--ui);\
padding:40px 24px;max-width:860px;margin:0 auto;-webkit-font-smoothing:antialiased}\
header{margin-bottom:24px}\
h1{margin:0;font-weight:600;font-size:30px;letter-spacing:-.01em}h1 .x{color:var(--accent)}\
.sub{color:var(--muted);margin-top:6px;font-size:13.5px;word-break:break-all}\
.report{background:var(--panel);border:1px solid var(--line);border-radius:14px;\
padding:22px 24px;margin:16px 0;box-shadow:var(--shadow)}\
.rlabel{font-family:var(--mono);font-size:12.5px;color:var(--muted);margin-bottom:10px;word-break:break-all}\
.vmain{font-size:22px;font-weight:600;letter-spacing:-.01em}\
.vmain.warn{color:var(--warn)}\
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(150px,1fr));gap:12px;margin-top:16px}\
.stat{background:var(--panel2);border:1px solid var(--line);border-radius:9px;padding:10px 12px}\
.stat .k{color:var(--muted);font-size:11.5px;text-transform:uppercase;letter-spacing:.03em}\
.stat .v{font-size:15px;font-weight:550;margin-top:3px;word-break:break-word}\
.pill{display:inline-block;margin-left:8px;padding:2px 10px;border-radius:20px;font-size:12px;\
font-weight:600;vertical-align:middle;color:var(--paper-on-accent)}\
.pill.high{background:var(--good)}.pill.medium{background:var(--warn)}.pill.low{background:var(--bad)}\
.track-lbl{color:var(--muted);font-size:12px;margin:18px 0 6px}\
.track{display:flex;gap:2px;height:30px;font-family:var(--mono);font-size:11px;font-weight:500}\
.seg{display:flex;align-items:center;justify-content:center;color:var(--paper-on-accent);\
border-radius:5px;overflow:hidden;white-space:nowrap;min-width:0;padding:0 4px}\
.seg-cb{background:var(--accent)}.seg-umi{background:var(--umi)}\
.seg-adapter{background:var(--panel2);color:var(--muted);border:1px solid var(--line)}\
.sexblock{margin-top:20px;padding-top:18px;border-top:1px solid var(--line)}\
.sexblock h3{margin:0 0 10px;font-size:15px;font-weight:600}.muted{color:var(--muted);font-weight:400;font-size:13px}\
.sexverdict{font-size:18px;font-weight:650;text-transform:capitalize}\
.sexverdict.male{color:var(--accent)}.sexverdict.female{color:var(--umi)}\
.sexverdict.uncertain,.sexverdict.unknown{color:var(--muted)}\
footer{margin-top:28px;color:var(--muted);font-size:12.5px}a{color:var(--accent)}";

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn renders_single_cell_report_with_geometry_and_sex() {
        let r = json!({
            "r1": "SRR123_1.fastq.gz",
            "chemistry": {"assay_display": "10x 3' v3", "confidence": "high",
                "is_single_cell": true, "cb_len": 16, "umi_len": 12, "cb_read": "R1",
                "r1_length_dominant": 28, "reads_sampled": 50000},
            "bulk_vs_sc": {"verdict": "single-cell"},
            "size": {"r1": {"mean_phred": 35.4}, "adequate_for_sc": true},
            "composition": {"gc_pct": 47.2, "guess_organism_family": "mammal"},
            "sex": {"verdict": "female", "confidence": "high", "species": "human",
                "panel_status": "validated", "y_hit_reads": 0, "xist_hit_reads": 41,
                "y_xist_ratio": 0.0, "reads_scanned": 200000}
        });
        let html = render(&r);
        assert!(html.starts_with("<!doctype html>"));
        assert!(html.contains("10x 3' v3"));
        assert!(html.contains("CB 16 + UMI 12 (R1)"));
        assert!(html.contains("R1 read structure")); // barcode-structure diagram
        assert!(html.contains("seg-cb") && html.contains("seg-umi"));
        assert!(html.contains("200,000")); // grouped reads_scanned
        assert!(html.contains("sexverdict female"));
        assert!(html.contains("pill high"));
    }

    #[test]
    fn report_palette_matches_app_tokens() {
        // The export report re-declares the palette as literals (it can't reach index.html's
        // :root tokens). Guard against drift: every hex it uses must still exist in the app CSS.
        let page = include_str!("../web/index.html");
        let mut hexes: Vec<&str> = CSS
            .split(|c: char| !c.is_ascii_hexdigit() && c != '#')
            .filter(|t| t.starts_with('#') && t.len() == 7)
            .collect();
        hexes.sort();
        hexes.dedup();
        assert!(hexes.len() >= 10, "expected the report's palette hexes, got {hexes:?}");
        for h in hexes {
            assert!(
                page.contains(h),
                "report CSS uses {h} but it is absent from web/index.html — palettes have \
                 drifted; update src/reporthtml.rs CSS to match the app tokens"
            );
        }
    }

    #[test]
    fn array_renders_one_section_per_run_with_labels() {
        let runs = json!([
            {"_label": "SRR1 · study", "chemistry": {"assay_id": "10x_3prime_v2", "is_single_cell": true}},
            {"_label": "SRR2 · study", "chemistry": {"is_single_cell": false, "r1_length_dominant": 150}}
        ]);
        let html = render(&runs);
        assert_eq!(html.matches("<section class=\"report\">").count(), 2);
        assert!(html.contains("SRR1 · study"));
        assert!(html.contains("Not single-cell"));
    }

    #[test]
    fn escapes_html_in_fields() {
        let r = json!({"r1": "<script>x</script>", "chemistry": null,
            "bulk_vs_sc": {"verdict": "bulk"}, "size": {}, "composition": {}});
        let html = render(&r);
        assert!(!html.contains("<script>x"));
        assert!(html.contains("&lt;script&gt;"));
    }
}
