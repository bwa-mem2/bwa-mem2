#!/usr/bin/env bash
# Byte-diff bwa-mem3 index output against the baselines committed under
# test/fixtures/baselines/. The baselines were produced at the legacy/sais-lite
# tag and capture the byte-exact reference artifacts we expect.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
BWAMEM3="$ROOT/bwa-mem3"
BASELINES="$HERE/fixtures/baselines"

FAIL=0
for fa in phix.fa synthetic_1mb.fa; do
    for ext in 0123 bwt.2bit.64; do
        [[ -s "$BASELINES/$fa.$ext" ]] || { echo "FAIL: baseline $BASELINES/$fa.$ext missing"; exit 1; }
    done
done

TD="$(mktemp -d)"
trap 'rm -rf "$TD"' EXIT

for fa in phix.fa synthetic_1mb.fa; do
    [[ -s "$HERE/fixtures/$fa" ]] || { echo "FAIL: source fixture $HERE/fixtures/$fa missing"; exit 1; }
    cp "$HERE/fixtures/$fa" "$TD/$fa"
    "$BWAMEM3" index "$TD/$fa" >/dev/null 2>&1
    for ext in 0123 bwt.2bit.64; do
        if ! cmp -s "$BASELINES/$fa.$ext" "$TD/$fa.$ext"; then
            echo "FAIL: $fa.$ext diverges from baseline"
            cmp -l "$BASELINES/$fa.$ext" "$TD/$fa.$ext" | head -5
            FAIL=1
        fi
    done
done

if [[ $FAIL -ne 0 ]]; then
    echo "libsais_index_diff_test: FAILED"
    exit 1
fi
echo "OK: libsais_index_diff (phix + synthetic_1mb match baselines)"
