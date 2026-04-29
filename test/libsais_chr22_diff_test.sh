#!/usr/bin/env bash
# Scale-up byte-diff test: index chr22 with libsais at HEAD and cmp against
# a reference index pre-built from the same FASTA (captured at legacy/sais-lite
# or a prior libsais run).
#
# This test is intentionally opt-in: the chr22 FASTA and baseline index are
# too large to check in, and CI does not export the env vars below, so the
# test silently skips on every CI platform. Run it manually on a developer
# host that has chr22 staged locally.
#
# Required env vars (no defaults — the test SKIPs if unset):
#   BWA_TEST_CHR22_FASTA    - path to chr22 FASTA
#   BWA_TEST_CHR22_BASELINE - directory holding the trusted baseline
#                             chr22.fa.{pac,ann,amb,0123,bwt.2bit.64}
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
BWAMEM3="$ROOT/bwa-mem3"

FIXTURE="${BWA_TEST_CHR22_FASTA:-}"
BASELINE_DIR="${BWA_TEST_CHR22_BASELINE:-}"

if [[ -z "$FIXTURE" || ! -s "$FIXTURE" ]]; then
    echo "SKIP: BWA_TEST_CHR22_FASTA unset or fixture missing"
    exit 0
fi
if [[ -z "$BASELINE_DIR" ]]; then
    echo "SKIP: BWA_TEST_CHR22_BASELINE unset"
    exit 0
fi

# Baseline must be a known-good snapshot from a trusted commit (e.g.
# legacy/sais-lite). Auto-bootstrapping with the CURRENT $BWAMEM3 would
# defeat the test: the freshly captured baseline and the just-built
# index are then both produced by the same backend at the same SHA, so
# cmp trivially passes regardless of correctness. Operators must seed
# $BASELINE_DIR explicitly. The one-time bootstrap path is gated behind
# BWA_TEST_CHR22_BASELINE_BOOTSTRAP=1 so it can't fire by accident.
# Check ALL five artifacts so a partial earlier build can't trick the
# cmp loop into reporting a spurious divergence on the missing extension.
EXTS=(pac ann amb 0123 bwt.2bit.64)
baseline_complete=1
for ext in "${EXTS[@]}"; do
    if [[ ! -s "$BASELINE_DIR/chr22.fa.$ext" ]]; then
        baseline_complete=0
        break
    fi
done
if [[ $baseline_complete -eq 0 ]]; then
    if [[ "${BWA_TEST_CHR22_BASELINE_BOOTSTRAP:-0}" != "1" ]]; then
        echo "FAIL: baseline at $BASELINE_DIR is missing or incomplete." >&2
        echo "      Seed it from a known-good commit (e.g. checkout legacy/sais-lite," >&2
        echo "      build, run \`bwa-mem3 index\` on the chr22 FASTA, copy the" >&2
        echo "      .pac/.ann/.amb/.0123/.bwt.2bit.64 outputs into \$BASELINE_DIR)," >&2
        echo "      or one-shot bootstrap with BWA_TEST_CHR22_BASELINE_BOOTSTRAP=1." >&2
        exit 1
    fi
    echo "INFO: BWA_TEST_CHR22_BASELINE_BOOTSTRAP=1 -- bootstrapping baseline at $BASELINE_DIR"
    mkdir -p "$BASELINE_DIR"
    cp "$FIXTURE" "$BASELINE_DIR/chr22.fa"
    "$BWAMEM3" index "$BASELINE_DIR/chr22.fa" >/dev/null
fi

TD="$(mktemp -d)"
trap 'rm -rf "$TD"' EXIT

cp "$FIXTURE" "$TD/chr22.fa"
"$BWAMEM3" index "$TD/chr22.fa" >/dev/null

FAIL=0
for ext in "${EXTS[@]}"; do
    if [[ ! -s "$BASELINE_DIR/chr22.fa.$ext" ]]; then
        echo "FAIL: baseline missing chr22.fa.$ext"
        FAIL=1
        continue
    fi
    if ! cmp -s "$BASELINE_DIR/chr22.fa.$ext" "$TD/chr22.fa.$ext"; then
        echo "FAIL: chr22.fa.$ext diverges from baseline"
        cmp -l "$BASELINE_DIR/chr22.fa.$ext" "$TD/chr22.fa.$ext" | head -3
        FAIL=1
    fi
done

if [[ $FAIL -ne 0 ]]; then
    echo "libsais_chr22_diff_test: FAILED"
    exit 1
fi
echo "OK: libsais_chr22_diff (all 5 artifacts match baseline)"
