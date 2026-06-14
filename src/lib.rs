//! preflightx — pre-flight sanity checks for single-cell RNA-seq FASTQ pairs.
//!
//! The crate is consumed two ways:
//!   * the native `preflightx` binary (`main.rs` → [`cli::run`]), and
//!   * a WebAssembly library ([`wasm`]) exposing the same checks to the browser.
//!
//! Every analysis module is platform-agnostic and funnels its I/O through
//! [`fastq`], which on wasm reads from an in-memory virtual file map instead of
//! the filesystem — so the analysis code is shared verbatim across both targets.

pub mod assay;
pub mod assaydetect;
pub mod chemistry;
pub mod composition;
pub mod fastq;
pub mod libstruct;
pub mod report;
pub mod sex;
pub mod sexpanel;
pub mod size;
pub mod stats;

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
