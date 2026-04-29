#!/usr/bin/env bash
# Build the synthetic 1 Mbp fixture three times with different -t values.
# libsais's OpenMP induced-sorting is deterministic: output .bwt.2bit.64
# must be byte-identical across all thread counts.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
BWAMEM3="$ROOT/bwa-mem3"
FA_SRC="$HERE/fixtures/synthetic_1mb.fa"
[[ -s "$FA_SRC" ]] || { echo "FAIL: $FA_SRC missing"; exit 1; }

RUN_DIR="$(mktemp -d)"
trap 'rm -rf "$RUN_DIR"' EXIT

for t in 1 4 8; do
    d="$RUN_DIR/t$t"
    mkdir -p "$d"
    cp "$FA_SRC" "$d/t.fa"
    "$BWAMEM3" index -t "$t" "$d/t.fa" >/dev/null
done

for ext in bwt.2bit.64 0123; do
    cmp -s "$RUN_DIR/t1/t.fa.$ext" "$RUN_DIR/t4/t.fa.$ext" \
        || { echo "FAIL: t=1 vs t=4 diverge on .$ext"; exit 1; }
    cmp -s "$RUN_DIR/t1/t.fa.$ext" "$RUN_DIR/t8/t.fa.$ext" \
        || { echo "FAIL: t=1 vs t=8 diverge on .$ext"; exit 1; }
done
echo "OK: libsais_determinism (t=1 == t=4 == t=8 on .bwt.2bit.64 + .0123)"
