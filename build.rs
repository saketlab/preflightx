//! Embed every seqspec assay YAML into the binary as a (filename, contents) table,
//! so the catalogue ships with no runtime file access. Only the tiny YAML
//! descriptions are embedded — never barcode data.
//!
//! The specs come from the upstream `pachterlab/seqspec` repo, pinned as a git
//! submodule under `vendor/seqspec` (not hand-copied), so the catalogue stays in
//! sync with upstream by bumping the submodule rev.

use std::{env, fs, path::Path};

fn main() {
    let manifest = env::var("CARGO_MANIFEST_DIR").unwrap();
    let dir = Path::new(&manifest).join("vendor/seqspec/docs/examples/assays");
    let out_dir = env::var("OUT_DIR").unwrap();

    let mut entries: Vec<_> = fs::read_dir(&dir)
        .unwrap_or_else(|e| {
            panic!(
                "cannot read {} ({e}).\n\
                 The seqspec catalogue is a git submodule — run:\n    \
                 git submodule update --init vendor/seqspec",
                dir.display()
            )
        })
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.to_str().map(|s| s.ends_with(".yaml")).unwrap_or(false))
        .collect();
    entries.sort();

    let mut s = String::from("pub static EMBEDDED_SPECS: &[(&str, &str)] = &[\n");
    for p in &entries {
        let name = p.file_name().unwrap().to_str().unwrap();
        s.push_str(&format!(
            "    ({:?}, include_str!({:?})),\n",
            name,
            p.to_str().unwrap()
        ));
    }
    s.push_str("];\n");

    fs::write(Path::new(&out_dir).join("specs.rs"), s).unwrap();
    println!("cargo:rerun-if-changed=vendor/seqspec/docs/examples/assays");
}
