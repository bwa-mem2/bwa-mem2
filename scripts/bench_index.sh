#!/usr/bin/env bash
# Run the bwa-mem3 index build across a few reference sizes, capture wall
# time and peak RSS via /usr/bin/time, emit a TSV and a human-readable
# summary. Outputs go under BENCH_DIR (default: $TMPDIR/bwa-index-bench;
# override with BENCH_DIR=... for a larger scratch volume).
#
# Usage:
#   scripts/bench_index.sh [targets...]
# Targets: chr1_4_6mb synthetic_1mb chr22 hg38_slice hg38 hg38_meth
# (default: all that have their source FASTA available)
#
# Note: chr1_4_6mb is a 4.6 Mbp slice of hg38 chr1, NOT a real E. coli
# reference. The slice has hg38 base composition / repeat structure, so
# its FM-index construction characteristics differ from authentic E. coli
# K-12. Treat its numbers as a "small-reference" data point, not as
# E. coli-comparable.
#
# Environment overrides:
#   BENCH_DIR        - where to stage work and write results
#   BWAMEM3          - path to bwa-mem3 binary (default: $ROOT/bwa-mem3)
#   BWA_TEST_HG38_FASTA / BWA_TEST_HG38_MET_C2T - source FASTA paths
#   BWA_INDEX_THREADS    - passed to `bwa-mem3 index -t`. Unset = auto-detect.
#   BWA_INDEX_MAX_MEMORY - passed to `bwa-mem3 index --max-memory`. Unset =
#       auto-detect (min(50% RAM, 32G)). The 32G cap rejects hg38 (~58 GiB
#       peak) and hg38_meth, so override here to reproduce the PR's numbers
#       (e.g. BWA_INDEX_MAX_MEMORY=128G on a 256 GiB host).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
BWAMEM3="${BWAMEM3:-$ROOT/bwa-mem3}"

# Build the index-builder argv once. Empty defaults mean "let bwa-mem3 pick",
# which is the documented auto-detect path. Setting an explicit budget makes
# wall-time / peak-RSS comparable across hosts.
INDEX_ARGS=()
[[ -n "${BWA_INDEX_THREADS:-}"    ]] && INDEX_ARGS+=(-t "$BWA_INDEX_THREADS")
[[ -n "${BWA_INDEX_MAX_MEMORY:-}" ]] && INDEX_ARGS+=(--max-memory "$BWA_INDEX_MAX_MEMORY")
BENCH_DIR="${BENCH_DIR:-${TMPDIR:-/tmp}/bwa-index-bench}"
# hg38 / hg38_meth FASTAs are too large to ship; if unset, the run_one
# helper SKIPs cleanly via its `[[ -s "$fasta" ]]` guard.
HG38_FASTA="${BWA_TEST_HG38_FASTA:-}"
HG38_MET_C2T="${BWA_TEST_HG38_MET_C2T:-}"

mkdir -p "$BENCH_DIR"
RESULTS_TSV="$BENCH_DIR/results.tsv"
SUMMARY="$BENCH_DIR/summary.txt"

UNAME_S="$(uname -s)"

parse_peak() {
    # $1: timing file. emits peak RSS in bytes.
    if [[ "$UNAME_S" == "Darwin" ]]; then
        grep -E '^ *[0-9]+ +maximum resident set size' "$1" | awk '{print $1}'
    else
        grep -E 'Maximum resident set size' "$1" | awk -F': +' '{print $2 * 1024}'
    fi
}

parse_wall_sec() {
    # Emits seconds (float). macOS: "N.NN real". Linux: "Elapsed ... m:ss".
    if [[ "$UNAME_S" == "Darwin" ]]; then
        grep -E '^ +[0-9.]+ +real' "$1" | awk '{print $1}'
    else
        # "Elapsed (wall clock) time (h:mm:ss or m:ss): X:YY.ZZ"
        grep -E 'Elapsed' "$1" | awk -F': +' '{print $2}' | awk -F: '
            NF==3 {print ($1*3600)+($2*60)+$3}
            NF==2 {print ($1*60)+$2}
            NF==1 {print $1}'
    fi
}

run_one() {
    local label="$1" fasta="$2"
    [[ -s "$fasta" ]] || { echo "SKIP: $label ($fasta missing)"; return 0; }
    echo "=== $label: $fasta ==="
    local wd="$BENCH_DIR/$label"
    rm -rf "$wd"
    mkdir -p "$wd"
    # Symlink the source FASTA: bwa-mem3 index reads it once, and copying
    # multi-GiB references (hg38 ≈ 3 GiB, hg38_meth doubled) wastes disk
    # and adds wall-clock noise before the timed step we actually care about.
    ln -sf "$fasta" "$wd/ref.fa"
    local timing="$wd/time.out"
    local logfile="$wd/index.log"

    # ${arr[@]+"${arr[@]}"} is the bash 3.2-safe idiom for "expand array if
    # set, else nothing" under set -u. Plain "${INDEX_ARGS[@]}" trips
    # `unbound variable` on an empty array under macOS's system bash.
    if [[ "$UNAME_S" == "Darwin" ]]; then
        /usr/bin/time -l "$BWAMEM3" index ${INDEX_ARGS[@]+"${INDEX_ARGS[@]}"} "$wd/ref.fa" >"$logfile" 2>"$timing" || {
            echo "FAIL: $label build failed"; cat "$timing"; return 1;
        }
    else
        /usr/bin/time -v "$BWAMEM3" index ${INDEX_ARGS[@]+"${INDEX_ARGS[@]}"} "$wd/ref.fa" >"$logfile" 2>"$timing" || {
            echo "FAIL: $label build failed"; cat "$timing"; return 1;
        }
    fi

    local peak wall size
    peak="$(parse_peak "$timing")"
    wall="$(parse_wall_sec "$timing")"
    if [[ -z "$peak" || -z "$wall" ]]; then
        echo "WARN: $label: failed to parse $timing (peak='$peak' wall='$wall')" >&2
        cat "$timing" >&2
        return 1
    fi
    size="$(stat -f%z "$fasta" 2>/dev/null || stat -c%s "$fasta")"
    printf "%s\t%s\t%.2f\t%s\n" "$label" "$size" "$wall" "$peak" >> "$RESULTS_TSV"
    printf "%-18s wall %7.2fs  peak %5.2f GiB  fasta %5.2f MiB\n" \
        "$label" "$wall" "$(awk -v b="$peak" 'BEGIN{printf "%.2f", b/1024/1024/1024}')" \
        "$(awk -v s="$size" 'BEGIN{printf "%.2f", s/1024/1024}')" | tee -a "$SUMMARY"
}

# Reset outputs for this run.
: > "$RESULTS_TSV"
: > "$SUMMARY"
printf "label\tfasta_bytes\twall_sec\tpeak_rss_bytes\n" >> "$RESULTS_TSV"

TARGETS=("$@")
if [[ ${#TARGETS[@]} -eq 0 ]]; then
    TARGETS=(chr1_4_6mb synthetic_1mb chr22 hg38_slice hg38 hg38_meth)
fi

# Source fixtures -- generated on first run and cached.
# chr1_4_6mb: 4.6 Mbp slice of hg38 chr1. Acts as a small-reference data
# point in the sweep; deliberately NOT named "ecoli" because its base
# composition / repeat structure are hg38-derived and diverge from real
# E. coli K-12.
CHR1_4_6MB="$BENCH_DIR/chr1_4_6mb.fa"
if [[ " ${TARGETS[*]} " == *" chr1_4_6mb "* && ! -s "$CHR1_4_6MB" ]]; then
    if command -v samtools >/dev/null 2>&1 && [[ -s "$HG38_FASTA" ]]; then
        echo "INFO: no chr1_4_6mb fixture; slicing hg38 chr1:1-4600000"
        samtools faidx "$HG38_FASTA" chr1:1-4600000 > "$CHR1_4_6MB"
    fi
fi

SYNTH="$ROOT/test/fixtures/synthetic_1mb.fa"
CHR22="$BENCH_DIR/chr22.fa"
if [[ " ${TARGETS[*]} " == *" chr22 "* && ! -s "$CHR22" ]]; then
    if command -v samtools >/dev/null 2>&1 && [[ -s "$HG38_FASTA" ]]; then
        samtools faidx "$HG38_FASTA" chr22 > "$CHR22"
    fi
fi

HG38_SLICE_DIR="$BENCH_DIR/hg38-slice"
HG38_SLICE="$HG38_SLICE_DIR/slice.fa"
if [[ " ${TARGETS[*]} " == *" hg38_slice "* && ! -s "$HG38_SLICE" ]]; then
    if command -v samtools >/dev/null 2>&1 && [[ -s "$HG38_FASTA" ]]; then
        mkdir -p "$HG38_SLICE_DIR"
        samtools faidx "$HG38_FASTA" chr1:1-10000000 > "$HG38_SLICE"
    fi
fi

for t in "${TARGETS[@]}"; do
    case "$t" in
        chr1_4_6mb)     run_one "chr1_4_6mb"     "$CHR1_4_6MB"  || true ;;
        synthetic_1mb)  run_one "synthetic_1mb"  "$SYNTH"       || true ;;
        chr22)          run_one "chr22"          "$CHR22"       || true ;;
        hg38_slice)     run_one "hg38_slice"     "$HG38_SLICE"  || true ;;
        hg38)           run_one "hg38_plain"     "$HG38_FASTA"  || true ;;
        hg38_meth)      run_one "hg38_meth_c2t"  "$HG38_MET_C2T" || true ;;
        *) echo "WARN: unknown target '$t'" ;;
    esac
done

echo
echo "--- results written to $RESULTS_TSV ---"
echo "--- summary in $SUMMARY ---"
