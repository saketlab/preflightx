//! Embed compile-time data tables into the binary as `(filename, contents)` slices, so
//! the crate ships with no runtime file access:
//!   * the seqspec assay catalogue (`vendor/seqspec`, a git submodule), and
//!   * the per-species sex-marker k-mer panels (`data/sexmarkers/*.panel`).
//!
//! Only tiny descriptions / public-genome-derived k-mers are embedded — never barcode data.

use std::{env, fs, path::Path, path::PathBuf};

fn main() {
    let manifest = env::var("CARGO_MANIFEST_DIR").unwrap();
    let out_dir = env::var("OUT_DIR").unwrap();

    // seqspec catalogue (a git submodule, so a missing dir is a setup error worth flagging).
    let specs_dir = Path::new(&manifest).join("vendor/seqspec/docs/examples/assays");
    let specs = read_sorted(&specs_dir, "yaml").unwrap_or_else(|e| {
        panic!(
            "cannot read {} ({e}).\n\
             The seqspec catalogue is a git submodule — run:\n    \
             git submodule update --init vendor/seqspec",
            specs_dir.display()
        )
    });
    write_table(&out_dir, "specs.rs", "EMBEDDED_SPECS", &specs);
    println!("cargo:rerun-if-changed=vendor/seqspec/docs/examples/assays");

    // sex-marker panels (built/validated in experiments/sex_panel; absence is non-fatal).
    let panels =
        read_sorted(&Path::new(&manifest).join("data/sexmarkers"), "panel").unwrap_or_default();
    write_table(&out_dir, "sexpanels.rs", "EMBEDDED_PANELS", &panels);
    println!("cargo:rerun-if-changed=data/sexmarkers");
}

/// Sorted paths of files in `dir` with the given extension.
fn read_sorted(dir: &Path, ext: &str) -> std::io::Result<Vec<PathBuf>> {
    let mut entries: Vec<_> = fs::read_dir(dir)?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().and_then(|x| x.to_str()) == Some(ext))
        .collect();
    entries.sort();
    Ok(entries)
}

/// Generate `pub static <name>: &[(&str, &str)]` of `(file_name, include_str!(path))`.
fn write_table(out_dir: &str, out_file: &str, name: &str, entries: &[PathBuf]) {
    let mut s = format!("pub static {name}: &[(&str, &str)] = &[\n");
    for p in entries {
        s.push_str(&format!(
            "    ({:?}, include_str!({:?})),\n",
            p.file_name().unwrap().to_str().unwrap(),
            p.to_str().unwrap()
        ));
    }
    s.push_str("];\n");
    fs::write(Path::new(out_dir).join(out_file), s).unwrap();
}
