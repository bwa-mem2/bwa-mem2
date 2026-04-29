#!/usr/bin/env bash
# 10 Mbp hg38 slice byte-diff. Extracts the first 10 Mbp of hg38 chr1 once
# into a cache and cmps subsequent libsais builds against its frozen
# baseline. Skips if hg38 or samtools is unavailable.
#
# This test is intentionally opt-in: the hg38 FASTA is too large to check in,
# and CI does not export BWA_TEST_HG38_FASTA, so the test silently skips on
# every CI platform. Run it manually on a developer host that has hg38 staged.
#
# Required env vars (no defaults — the test SKIPs if unset):
#   BWA_TEST_HG38_FASTA      - path to hg38 FASTA (samtools faidx-able)
#
# Optional env vars:
#   BWA_TEST_HG38_SLICE_DIR  - cache directory for the extracted slice and
#                              frozen baseline (default: $TMPDIR/hg38-slice)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
BWAMEM3="$ROOT/bwa-mem3"

HG38="${BWA_TEST_HG38_FASTA:-}"
SLICE_DIR="${BWA_TEST_HG38_SLICE_DIR:-${TMPDIR:-/tmp}/hg38-slice}"

if [[ -z "$HG38" || ! -s "$HG38" ]]; then
    echo "SKIP: BWA_TEST_HG38_FASTA unset or fixture missing"
    exit 0
fi
if ! command -v samtools >/dev/null; then
    echo "SKIP: samtools not in PATH"
    exit 0
fi

mkdir -p "$SLICE_DIR"
if [[ ! -s "$SLICE_DIR/slice.fa" ]]; then
    echo "INFO: extracting 10 Mbp slice (chr1:1-10000000)"
    # Write to a temp sibling and rename on success: a partial faidx
    # write (signal, network FS hiccup, disk full) would otherwise leave
    # a truncated slice.fa in place, and the next run's `[[ -s ... ]]`
    # guard would happily reuse it and silently build the baseline
    # against truncated input.
    samtools faidx "$HG38" chr1:1-10000000 > "$SLICE_DIR/slice.fa.tmp"
    mv "$SLICE_DIR/slice.fa.tmp" "$SLICE_DIR/slice.fa"
fi

# Baseline must be a known-good snapshot from a trusted commit (e.g.
# legacy/sais-lite). Auto-bootstrapping with the CURRENT $BWAMEM3 would
# defeat the test: the freshly captured baseline and the just-built
# index are then both produced by the same backend at the same SHA, so
# cmp trivially passes regardless of correctness. The one-time bootstrap
# path is gated behind BWA_TEST_HG38_SLICE_BASELINE_BOOTSTRAP=1 so it
# can't fire by accident.
if [[ ! -s "$SLICE_DIR/baseline/slice.fa.bwt.2bit.64" ]]; then
    if [[ "${BWA_TEST_HG38_SLICE_BASELINE_BOOTSTRAP:-0}" != "1" ]]; then
        echo "FAIL: baseline at $SLICE_DIR/baseline is missing or incomplete." >&2
        echo "      Seed it from a known-good commit (e.g. checkout legacy/sais-lite," >&2
        echo "      build, run \`bwa-mem3 index\` on the slice FASTA, copy the" >&2
        echo "      .pac/.ann/.amb/.0123/.bwt.2bit.64 outputs into \$SLICE_DIR/baseline)," >&2
        echo "      or one-shot bootstrap with BWA_TEST_HG38_SLICE_BASELINE_BOOTSTRAP=1." >&2
        exit 1
    fi
    echo "INFO: BWA_TEST_HG38_SLICE_BASELINE_BOOTSTRAP=1 -- bootstrapping baseline at $SLICE_DIR/baseline"
    mkdir -p "$SLICE_DIR/baseline"
    cp "$SLICE_DIR/slice.fa" "$SLICE_DIR/baseline/slice.fa"
    "$BWAMEM3" index "$SLICE_DIR/baseline/slice.fa" >/dev/null
fi

TD="$(mktemp -d)"
trap 'rm -rf "$TD"' EXIT

cp "$SLICE_DIR/slice.fa" "$TD/slice.fa"
"$BWAMEM3" index "$TD/slice.fa" >/dev/null

FAIL=0
for ext in pac ann amb 0123 bwt.2bit.64; do
    if ! cmp -s "$SLICE_DIR/baseline/slice.fa.$ext" "$TD/slice.fa.$ext"; then
        echo "FAIL: slice.fa.$ext diverges from baseline"
        cmp -l "$SLICE_DIR/baseline/slice.fa.$ext" "$TD/slice.fa.$ext" | head -3
        FAIL=1
    fi
done

if [[ $FAIL -ne 0 ]]; then
    echo "libsais_hg38_slice_diff_test: FAILED"
    exit 1
fi
echo "OK: libsais_hg38_slice_diff (10 Mbp slice all 5 artifacts match baseline)"
