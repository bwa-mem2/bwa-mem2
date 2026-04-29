#!/usr/bin/env bash
# Memory-budget-honored test. Builds the 1 Mbp synthetic fixture under
# several --max-memory settings, measures actual peak RSS via /usr/bin/time,
# and fails if peak exceeds the budget by more than 10% slack. The preflight
# should gate any case where the libsais estimate would exceed the budget.
#
# Parses /usr/bin/time -l output on macOS and /usr/bin/time -v on Linux.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
BWAMEM3="$ROOT/bwa-mem3"
FA_SRC="$HERE/fixtures/synthetic_1mb.fa"
[[ -s "$FA_SRC" ]] || { echo "FAIL: $FA_SRC missing"; exit 1; }

UNAME_S="$(uname -s)"
if [[ "$UNAME_S" == "Darwin" ]]; then
    TIME_CMD=(/usr/bin/time -l)
else
    TIME_CMD=(/usr/bin/time -v)
fi

parse_peak() {
    # $1: timing output path. echoes peak RSS in bytes.
    local out="$1"
    if [[ "$UNAME_S" == "Darwin" ]]; then
        grep -E '^ *[0-9]+ +maximum resident set size' "$out" \
            | awk '{print $1}'
    else
        grep -E 'Maximum resident set size' "$out" \
            | awk -F': +' '{print $2 * 1024}'
    fi
}

# Budgets chosen (128 / 512 / 2048 MiB) exercise the preflight + actual
# peak-RSS contract above the fixed thread-pool / libomp / mimalloc-arena
# infrastructure cost (~15 MiB on this host that's independent of the
# budget). Tighter budgets below ~64 MiB would fail the +10% slack check
# even on a trivially small input, so we lean on the preflight test
# (libsais_index_diff) to cover the "does the budget get respected"
# contract at larger scales.
FAIL=0
for budget_mib in 128 512 2048; do
    budget_bytes=$(( budget_mib * 1024 * 1024 ))
    slack_bytes=$(( budget_bytes * 11 / 10 ))
    TD="$(mktemp -d)"
    # Cleanup unconditionally on exit (parse-failure / set -e / signal)
    # so $TD doesn't leak between iterations or on the early-exit paths.
    trap 'rm -rf "$TD"' EXIT
    cp "$FA_SRC" "$TD/t.fa"
    TIMING="$TD/time.out"
    "${TIME_CMD[@]}" "$BWAMEM3" index --max-memory "${budget_mib}M" "$TD/t.fa" >"$TD/stdout" 2>"$TIMING" || {
        echo "FAIL: build failed at --max-memory ${budget_mib}M"
        cat "$TIMING" | tail -20
        exit 1
    }
    peak="$(parse_peak "$TIMING")"
    [[ -n "$peak" ]] || { echo "FAIL: could not parse peak RSS"; cat "$TIMING"; exit 1; }
    if [[ "$peak" -gt "$slack_bytes" ]]; then
        echo "FAIL: --max-memory ${budget_mib}M -> peak $(( peak / 1024 / 1024 )) MiB (budget ${budget_mib}M + 10% = $(( slack_bytes / 1024 / 1024 )) MiB)"
        FAIL=1
    else
        printf "OK:   --max-memory %dM -> peak %d MiB (under budget+10%%)\n" \
            "$budget_mib" "$(( peak / 1024 / 1024 ))"
    fi
    rm -rf "$TD"
    trap - EXIT
done

if [[ $FAIL -ne 0 ]]; then
    echo "libsais_memory_budget_test: FAILED"
    exit 1
fi
echo "OK: libsais_memory_budget (all budgets respected)"
