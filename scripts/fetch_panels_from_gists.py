#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

API = "https://api.github.com/gists/{gist_id}"


def _req(url: str) -> urllib.request.Request:
    headers = {"User-Agent": "preflightx-fetch-panels"}
    token = os.environ.get("GITHUB_TOKEN")
    if token:
        headers["Authorization"] = f"Bearer {token}"
    return urllib.request.Request(url, headers=headers)


def fetch(panel: dict, out: Path) -> tuple[str, bool, str]:
    name = panel["spanel"]
    dst = out / name
    if dst.exists() and dst.stat().st_size > 1024:
        return name, True, "cached"
    try:
        with urllib.request.urlopen(
            _req(API.format(gist_id=panel["gist_id"])), timeout=60
        ) as r:
            meta = json.load(r)
        f = meta["files"].get(name) or next(iter(meta["files"].values()))
        with urllib.request.urlopen(_req(f["raw_url"]), timeout=120) as r:
            data = r.read()
        if len(data) < 1024:
            return name, False, "too small"
        dst.write_bytes(data)
        return name, True, "fetched"
    except Exception as e:  # noqa: BLE001
        return name, False, repr(e)[:60]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prov", type=Path, default="web/data/panels.json")
    ap.add_argument("--out", type=Path, default="web/data/panels")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    panels = [
        p for p in json.loads(args.prov.read_text())["panels"] if p.get("gist_id")
    ]
    args.out.mkdir(parents=True, exist_ok=True)
    ok = fail = 0
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        for fut in as_completed([ex.submit(fetch, p, args.out) for p in panels]):
            name, good, msg = fut.result()
            ok += good
            fail += not good
            if not good:
                print(f"  FAIL {name}: {msg}", file=sys.stderr)
    print(
        f"fetched {ok}/{len(panels)} panels -> {args.out} ({fail} failed)",
        file=sys.stderr,
    )
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
