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
    for wl in &whitelists {
        let target = dir.join(wl.filename);
        if target.exists() && fs::metadata(&target)?.len() > 1_000_000 {
            println!("  [{}] {} -- already cached", wl.filename, "skip".dimmed());
            continue;
        }
        print!("  [{}] downloading ... ", wl.filename);
        std::io::stdout().flush()?;
        match download_one(&client, wl.url, &target) {
            Ok(n) => println!("{} ({} MB)", "ok".green(), n / (1024 * 1024)),
            Err(e) => {
                any_failed = true;
                println!("{} ({})", "failed".red(), e);
            }
        }
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
