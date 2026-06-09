//! Native binary entry point. All logic lives in the library crate
//! (`preflightx::cli`); the wasm build produces a cdylib instead and never
//! compiles this binary's real body.

#[cfg(not(target_arch = "wasm32"))]
fn main() -> anyhow::Result<()> {
    preflightx::cli::run()
}

#[cfg(target_arch = "wasm32")]
fn main() {}
