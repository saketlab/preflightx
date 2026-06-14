use anyhow::{anyhow, Context, Result};
use colored::Colorize;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use crate::assay::ASSAYS;

pub fn cache_dir() -> Result<PathBuf> {
    if let Ok(d) = std::env::var("PREFLIGHTX_CACHE") {
        return Ok(PathBuf::from(d));
    }
    let base = dirs::cache_dir().ok_or_else(|| anyhow!("could not resolve a cache directory"))?;
    Ok(base.join("preflightx").join("whitelists"))
}

pub fn fetch_all() -> Result<()> {
    let dir = cache_dir()?;
    fs::create_dir_all(&dir)?;
    println!("Fetching whitelists into {}", dir.display());
    println!(
        "{}",
        "Source: scg_lib_structs mirror (teichlab.github.io). Subject to each vendor's EULA."
            .dimmed()
    );
    println!();

    let mut seen = std::collections::HashSet::new();
    let mut whitelists = Vec::new();
    for spec in ASSAYS {
        if let Some(wl) = &spec.whitelist {
            if seen.insert(wl.filename) {
                whitelists.push(wl);
            }
        }
    }

    let client = reqwest::blocking::ClientBuilder::new()
        .timeout(std::time::Duration::from_secs(900))
        .build()?;
    let mut any_failed = false;
    // 10x whitelists are large; a partial file under 1 MB is treated as not cached.
    for wl in &whitelists {
        let cached = dir.join(wl.filename).exists()
            && fs::metadata(dir.join(wl.filename))?.len() > 1_000_000;
        any_failed |= !fetch_into(&client, &dir, wl.filename, wl.url, cached)?;
    }
    // Non-10x barcode onlists for `detect-assay`: small plain-text lists.
    for src in crate::assaydetect::onlist_fetch_list() {
        let cached = dir.join(src.cache_name).exists();
        any_failed |= !fetch_into(&client, &dir, src.cache_name, src.url, cached)?;
    }
    println!();
    if any_failed {
        println!("{}", "Some whitelists could not be downloaded.".yellow());
        println!(
            "{}",
            "Manual fallback: copy whitelist .txt.gz files into the cache dir shown below."
                .dimmed()
        );
        println!();
    }
    println!("Cache directory: {}", dir.display());
    Ok(())
}

/// Download `name` into `dir` unless already `cached`, printing one status line.
/// Returns `false` only on a download failure (so callers can track `any_failed`).
fn fetch_into(
    client: &reqwest::blocking::Client,
    dir: &Path,
    name: &str,
    url: &str,
    cached: bool,
) -> Result<bool> {
    if cached {
        println!("  [{}] {} -- already cached", name, "skip".dimmed());
        return Ok(true);
    }
    print!("  [{name}] downloading ... ");
    std::io::stdout().flush()?;
    match download_one(client, url, &dir.join(name)) {
        Ok(n) if n >= 1024 * 1024 => {
            println!("{} ({} MB)", "ok".green(), n / (1024 * 1024));
            Ok(true)
        }
        Ok(n) => {
            println!("{} ({} KB)", "ok".green(), (n / 1024).max(1));
            Ok(true)
        }
        Err(e) => {
            println!("{} ({})", "failed".red(), e);
            Ok(false)
        }
    }
}

fn download_one(client: &reqwest::blocking::Client, url: &str, target: &Path) -> Result<usize> {
    let resp = client
        .get(url)
        .send()
        .with_context(|| format!("GET {url}"))?;
    if !resp.status().is_success() {
        return Err(anyhow!("HTTP {}", resp.status()));
    }
    let bytes = resp.bytes()?;
    fs::write(target, &bytes).with_context(|| format!("writing {}", target.display()))?;
    Ok(bytes.len())
}
