//! preflightx — pre-flight sanity checks for single-cell RNA-seq FASTQ pairs.
//!
//! Consumed two ways: the native `preflightx` binary (`main.rs` → [`cli::run`])
//! and a WebAssembly library ([`wasm`]). Analysis modules are platform-agnostic
//! and funnel I/O through [`fastq`], which on wasm reads from an in-memory file
//! map instead of the filesystem, so analysis code is shared across both targets.

pub mod assay;
pub mod assaydetect;
pub mod atacpanel;
pub mod chemistry;
pub mod composition;
pub mod endpanel;
pub mod fastq;
pub mod kmerpanel;
pub mod libstruct;
pub mod nucpanel;
pub mod report;
pub mod reporthtml;
pub mod sex;
pub mod sexpanel;
pub mod size;
pub mod snpanel;
pub mod species;
pub mod stats;

// Native-only: the cached `.spanel` catalogue (filesystem cache + ~130 MB bundle).
#[cfg(not(target_arch = "wasm32"))]
pub mod speciespanel;

// Native-only: the clap CLI front-end and whitelist download/cache management
// (reqwest + dirs). Neither is reachable from the wasm build.
#[cfg(not(target_arch = "wasm32"))]
pub mod cli;
#[cfg(not(target_arch = "wasm32"))]
pub mod whitelists;

// Browser entry points (wasm-bindgen).
#[cfg(target_arch = "wasm32")]
pub mod wasm;

#[cfg(test)]
pub mod testutil;
