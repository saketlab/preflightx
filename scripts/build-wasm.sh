#!/usr/bin/env bash
# Build the preflightx WebAssembly module + JS bindings into web/pkg/.
set -euo pipefail
cd "$(dirname "$0")/.."

TOOLCHAIN="${PREFLIGHTX_WASM_TOOLCHAIN:-stable}"
OUT_DIR="web/pkg"
PROFILE="release"

# Resolve the rustup toolchain's bin dir (works even when conda shadows cargo).
RUSTC_PATH="$(rustup which --toolchain "$TOOLCHAIN" rustc)"
TC_BIN="$(dirname "$RUSTC_PATH")"
export RUSTC="$TC_BIN/rustc"
export PATH="$TC_BIN:$PATH"
CARGO="$TC_BIN/cargo"

echo ">> toolchain: $("$RUSTC" --version)"

# Ensure the wasm target is present for this toolchain.
if ! "$RUSTC" --print target-list | grep -qx wasm32-unknown-unknown \
   || ! ls "$("$RUSTC" --print sysroot)/lib/rustlib/wasm32-unknown-unknown/lib/" >/dev/null 2>&1; then
  echo ">> installing wasm32-unknown-unknown for $TOOLCHAIN"
  rustup target add --toolchain "$TOOLCHAIN" wasm32-unknown-unknown
fi

echo ">> cargo build (wasm, $PROFILE)"
"$CARGO" build --lib --target wasm32-unknown-unknown --no-default-features --"$PROFILE"

WASM="target/wasm32-unknown-unknown/$PROFILE/preflightx.wasm"

echo ">> wasm-bindgen ($(wasm-bindgen --version))"
wasm-bindgen --target web --no-typescript --out-dir "$OUT_DIR" "$WASM"

# Size optimization (best-effort). Modern rustc emits bulk-memory/sign-ext/etc., which
# older wasm-opt rejects unless enabled; write to a temp file and only swap on success
# so a failed/old wasm-opt never corrupts or blocks the build.
if command -v wasm-opt >/dev/null 2>&1; then
  echo ">> wasm-opt -Oz ($(wasm-opt --version 2>/dev/null | head -1))"
  if wasm-opt -Oz \
       --enable-bulk-memory --enable-nontrapping-float-to-int --enable-sign-ext \
       --enable-mutable-globals --enable-multivalue --enable-reference-types \
       "$OUT_DIR/preflightx_bg.wasm" -o "$OUT_DIR/preflightx_bg.wasm.opt"; then
    mv "$OUT_DIR/preflightx_bg.wasm.opt" "$OUT_DIR/preflightx_bg.wasm"
  else
    echo ">> wasm-opt failed; shipping the unoptimized wasm"
    rm -f "$OUT_DIR/preflightx_bg.wasm.opt"
  fi
fi

# Precompress for nginx gzip_static/brotli_static.
echo ">> precompressing bundle (gzip$(command -v brotli >/dev/null 2>&1 && echo ' + brotli'))"
for f in "$OUT_DIR/preflightx_bg.wasm" "$OUT_DIR/preflightx.js"; do
  gzip -9 -k -f "$f"
  command -v brotli >/dev/null 2>&1 && brotli -q 11 -f -o "$f.br" "$f"
done

echo ">> done -> $OUT_DIR/"
ls -lh "$OUT_DIR"
