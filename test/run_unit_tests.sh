#!/usr/bin/env bash
# Local unit-test harness for the five C++ binaries under test/.
# Builds, indexes the phiX fixture, runs each binary, asserts exit 0
# and non-empty output. No hash pinning — that lives in the CI workflow.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
BWAMEM3="$ROOT/bwa-mem3"
FIXTURES="$HERE/fixtures"

fail() { echo "FAIL: $*" >&2; exit 1; }
ok()   { echo "OK:   $*"; }

[[ -x "$BWAMEM3" ]] || fail "bwa-mem3 not built at $BWAMEM3 (run 'make' or 'make arm64' first)"

# Build the five unit binaries.
(cd "$HERE" && make) || fail "test/ make failed"

# Build any test binaries that live outside test/Makefile.
( cd "$ROOT" && make -j4 shm_section_find_test shm_pack_round_trip_test ) >/dev/null

# synthetic_1mb.fa is checked in alongside its committed baselines so the
# byte-diff test stays reproducible across Python versions and platforms.
# A previous version of this script regenerated the FASTA via random.choice;
# random.choice's bytes are deterministic per Python version, but the
# baselines were built from a different generator, so CI ran the diff
# against a mismatched source and flaked.
[[ -s "$FIXTURES/synthetic_1mb.fa" ]] || fail "synthetic_1mb.fa fixture missing at $FIXTURES/synthetic_1mb.fa"

# Build the phiX FMI index if not already present. Check all five artifacts
# that `bwa-mem3 index` produces so a corrupt/partial prior run is re-indexed.
if [[ ! -s "$FIXTURES/phix.fa.bwt.2bit.64" || \
      ! -s "$FIXTURES/phix.fa.0123"        || \
      ! -s "$FIXTURES/phix.fa.amb"         || \
      ! -s "$FIXTURES/phix.fa.ann"         || \
      ! -s "$FIXTURES/phix.fa.pac" ]]; then
    "$BWAMEM3" index "$FIXTURES/phix.fa" >/dev/null 2>&1 || fail "bwa-mem3 index on phix.fa failed"
fi

# --- fmi_test --------------------------------------------------------------
OUT="$(cd "$HERE" && ./fmi_test "$FIXTURES/phix.fa" "$FIXTURES/reads.fa" 10 19 1 2>&1)"
echo "$OUT" | grep -q 'numReads = 10'       || fail "fmi_test: numReads line missing"
echo "$OUT" | grep -q 'totalSmems ='        || fail "fmi_test: totalSmems line missing"
echo "$OUT" | grep -qE '\[[0-9]+,[0-9]+\]'  || fail "fmi_test: no SMEM output (PRINT_OUTPUT patch not applied?)"
ok "fmi_test"

# --- smem2_test ------------------------------------------------------------
OUT="$(cd "$HERE" && ./smem2_test "$FIXTURES/phix.fa" "$FIXTURES/smem2_input.txt" 5 50 19 1 2>&1)"
echo "$OUT" | grep -q 'numReads = 5'        || fail "smem2_test: numReads line missing"
echo "$OUT" | grep -q 'totalSmems ='        || fail "smem2_test: totalSmems line missing"
ok "smem2_test"

# --- bwt_seed_strategy_test ------------------------------------------------
# Note: unlike fmi_test and smem2_test, this binary doesn't emit a
# "numReads = N" line. Only assert the post-SW totalSmems marker.
OUT="$(cd "$HERE" && ./bwt_seed_strategy_test "$FIXTURES/phix.fa" "$FIXTURES/bwt_seed_input.fa" 5 50 19 20 2>&1)"
echo "$OUT" | grep -q 'minSeedLen ='        || fail "bwt_seed_strategy_test: minSeedLen line missing (read-parsing broken?)"
echo "$OUT" | grep -q 'totalSmems ='        || fail "bwt_seed_strategy_test: totalSmems line missing"
ok "bwt_seed_strategy_test"

# --- sa2ref_test -----------------------------------------------------------
OUT_FILE="$(mktemp)"
FKSW="$HERE/fksw.txt"
HEADER_OUT="$(mktemp)"
trap 'rm -f "$OUT_FILE" "$FKSW" "$HEADER_OUT"' EXIT
(cd "$HERE" && ./sa2ref_test "$FIXTURES/phix.fa" "$FIXTURES/sa2ref_input.txt" "$OUT_FILE" 20 >/dev/null 2>&1) || fail "sa2ref_test crashed"
[[ -s "$OUT_FILE" ]] || fail "sa2ref_test: output file empty"
ok "sa2ref_test (wrote $(wc -l < "$OUT_FILE" | tr -d ' ') coords)"

# --- header_insert_test ----------------------------------------------------
(cd "$HERE" && ./header_insert_test 2>&1) | tee "$HEADER_OUT" >/dev/null || fail "header_insert_test crashed"
grep -q 'ALL HEADER INSERT TESTS PASSED' "$HEADER_OUT" \
    || fail "header_insert_test: final banner missing ($(tail -n1 "$HEADER_OUT"))"
ok "header_insert_test"

# --- xeonbsw (main_banded) -------------------------------------------------
(cd "$HERE" && ./xeonbsw -pairs "$FIXTURES/pairs.txt" >/dev/null 2>&1) || fail "xeonbsw crashed"
[[ -s "$FKSW" ]] || fail "xeonbsw: fksw.txt empty (main_banded fksw-write patch not applied?)"
LINES="$(wc -l < "$FKSW" | tr -d ' ')"
[[ "$LINES" -eq 3 ]] || fail "xeonbsw: expected 3 lines in fksw.txt, got $LINES"
ok "xeonbsw ($LINES pair scores emitted)"

# --- pg_cl_escape_test ----------------------------------------------------
# Regression for issue #45 / upstream #293: tabs inside `-R` must not
# bleed into the @PG CL: value. Uses the same phiX fixture as the rest
# of the harness (indexed above if needed).
"$HERE/pg_cl_escape_test.sh" "$BWAMEM3" "$FIXTURES" || fail "pg_cl_escape_test failed"
ok "pg_cl_escape_test"

# --- help_prescan_test ----------------------------------------------------
# `mem --help` pre-scan must not match `--help` when it is the value of
# an option that takes an argument (-R, -o, --set-as-failed, ...).
"$HERE/help_prescan_test.sh" "$BWAMEM3" "$FIXTURES" || fail "help_prescan_test failed"
ok "help_prescan_test"

# --- smem_lockstep_parity_test --------------------------------------------
OUT="$(cd "$HERE" && ./smem_lockstep_parity_test "$FIXTURES/phix.fa" 2>&1)"
CASES_PASSED="$(echo "$OUT" | sed -nE 's/^([0-9]+) \/ ([0-9]+) cases passed$/\1/p')"
CASES_TOTAL="$(echo "$OUT"  | sed -nE 's/^([0-9]+) \/ ([0-9]+) cases passed$/\2/p')"
[[ -n "$CASES_TOTAL" ]]              || fail "smem_lockstep_parity_test: no summary line"
[[ "$CASES_PASSED" == "$CASES_TOTAL" ]] || fail "smem_lockstep_parity_test: $CASES_PASSED / $CASES_TOTAL cases passed"
ok "smem_lockstep_parity_test ($CASES_PASSED / $CASES_TOTAL)"

# --- long-read end-to-end (issue 44) --------------------------------------
# Pre-fix, reads > 151 bp overran the per-thread SMEM buffer (segfault) and
# reads > 512 bp tripped the MAX_READ_LEN_FOR_LOCKSTEP assert. Post-fix the
# SMEM and lockstep-slot buffers are sized per-batch, so these must align
# cleanly and produce SAM lines.
for LEN_FQ in long_read_300bp long_read_1kbp long_read_3kbp; do
    RAW="$("$BWAMEM3" mem "$FIXTURES/phix.fa" "$FIXTURES/${LEN_FQ}.fq" 2>/dev/null)" \
        || fail "bwa-mem3 mem ${LEN_FQ}.fq: non-zero exit (crash regression)"
    OUT="$(printf '%s\n' "$RAW" | grep -v '^@' || true)"
    [[ -n "$OUT" ]] || fail "bwa-mem3 mem ${LEN_FQ}.fq: no SAM records emitted"
    # Every record should be a successful alignment (not flag 4 = unmapped).
    UNMAPPED="$(echo "$OUT" | awk '$2 == 4 || $2 == 2052' | wc -l | tr -d ' ')"
    [[ "$UNMAPPED" == "0" ]] || fail "bwa-mem3 mem ${LEN_FQ}.fq: $UNMAPPED unmapped record(s)"
    ok "bwa-mem3 mem ${LEN_FQ}.fq (all records mapped)"
done


# --- interleaved -p mode regression --------------------------------------
# Bugs covered:
#  - Empty/zero-length-read batches (which a malformed FASTQ that bseq parses
#    as a single 0-bp record would trigger) used to fall through the lazy-init
#    grow checks in mem_collect_smem with NULL per-thread buffers and SEGV in
#    the populate loop. mem_collect_smem now early-returns on zero-base batches.
#  - The OLD `tot_len >= mmc->wsize_mem[tid]` grow block in mem_kernel1_core
#    used to fire when both sides were 0, calling _mm_realloc(NULL, 0, 0, ...)
#    which returns NULL and then exit(1) on the post-collect check. That block
#    is gone; mem_collect_smem owns all SMEM-family sizing.
PE_TMP="$(mktemp -d)"
# Append to (don't replace) the EXIT trap installed earlier for OUT_FILE /
# FKSW / HEADER_OUT — overwriting it would leak those paths once this block
# runs.
trap 'rm -f "$OUT_FILE" "$FKSW" "$HEADER_OUT"; rm -rf "$PE_TMP"' EXIT

# Build a small properly-interleaved FASTQ from phix. Use python so the
# fixture is generated programmatically (per repo convention) and so the
# interleaving is correct at FASTQ-record granularity.
python3 - "$FIXTURES/phix.fa" "$PE_TMP/interleaved.fq" <<'PY'
import random, sys
ref_path, out_path = sys.argv[1], sys.argv[2]
ref = ''
with open(ref_path) as f:
    next(f)
    for line in f:
        ref += line.strip()
random.seed(42)
N, RL, INS = 50, 120, 280
comp = str.maketrans('ACGTN', 'TGCAN')
with open(out_path, 'w') as out:
    for i in range(N):
        pos = random.randint(0, len(ref) - INS)
        r1 = ref[pos:pos+RL]
        r2 = ref[pos+INS-RL:pos+INS][::-1].translate(comp)
        q  = 'I' * RL
        out.write(f'@pe{i}/1\n{r1}\n+\n{q}\n')
        out.write(f'@pe{i}/2\n{r2}\n+\n{q}\n')
PY

OUT="$("$BWAMEM3" mem -p "$FIXTURES/phix.fa" "$PE_TMP/interleaved.fq" 2>/dev/null)" \
    || fail "bwa-mem3 mem -p interleaved.fq: non-zero exit (crash regression)"
RECORDS="$(printf '%s\n' "$OUT" | grep -cv '^@' || true)"
[[ "$RECORDS" -ge 100 ]] || fail "bwa-mem3 mem -p: expected >=100 records (50 pairs × 2), got $RECORDS"
ok "bwa-mem3 mem -p interleaved.fq ($RECORDS records)"

# Empty FASTQ + zero-length-read FASTQ regressions for the empty-batch
# defensive path in mem_collect_smem / mem_kernel1_core.
: > "$PE_TMP/empty.fq"
"$BWAMEM3" mem "$FIXTURES/phix.fa" "$PE_TMP/empty.fq" >/dev/null 2>&1 \
    || fail "bwa-mem3 mem on empty FASTQ: non-zero exit"
ok "bwa-mem3 mem on empty FASTQ (no crash)"

printf '@zero\n\n+\n\n' > "$PE_TMP/zero.fq"
"$BWAMEM3" mem "$FIXTURES/phix.fa" "$PE_TMP/zero.fq" >/dev/null 2>&1 \
    || fail "bwa-mem3 mem on zero-length read: non-zero exit"
ok "bwa-mem3 mem on zero-length read (no crash)"

# --- packed_text_test ------------------------------------------------------
(cd "$HERE" && ./packed_text_test) || fail "packed_text_test"
ok "packed_text_test"

# --- system_test -----------------------------------------------------------
(cd "$HERE" && ./system_test) || fail "system_test"
ok "system_test"

# --- libsais byte-diff against committed baselines -------------------------
(cd "$HERE" && ./libsais_index_diff_test.sh) || fail "libsais_index_diff_test"
ok "libsais_index_diff_test"

# --- libsais determinism across thread counts ------------------------------
(cd "$HERE" && ./libsais_determinism_test.sh) || fail "libsais_determinism_test"
ok "libsais_determinism_test"

# --- libsais chr22 byte-diff (skips if scratch fixture absent) -------------
(cd "$HERE" && ./libsais_chr22_diff_test.sh) || fail "libsais_chr22_diff_test"
ok "libsais_chr22_diff_test"

# --- libsais 10 Mbp hg38 slice byte-diff (skips if scratch fixture absent) -
(cd "$HERE" && ./libsais_hg38_slice_diff_test.sh) || fail "libsais_hg38_slice_diff_test"
ok "libsais_hg38_slice_diff_test"

# --- libsais --max-memory peak-RSS budget test -----------------------------
(cd "$HERE" && ./libsais_memory_budget_test.sh) || fail "libsais_memory_budget_test"
ok "libsais_memory_budget_test"

# --- shm_section_find_test (unit: section-find helper) -----------------------
echo "==> shm_section_find_test"
( cd "$ROOT" && ./shm_section_find_test )
echo "OK:   shm_section_find_test"

# --- shm pack round-trip (unit: pack/unpack a phiX segment) ------------------
(cd "$HERE" && ./shm_pack_round_trip_test.sh) || fail "shm_pack_round_trip_test"
ok "shm_pack_round_trip_test"

# --- shm end-to-end round-trip (byte-identical SAM with and without shm) -----
(cd "$HERE" && ./shm_round_trip_test.sh) || fail "shm_round_trip_test"
ok "shm_round_trip_test"

# --- shm --meth flag (parity with mem --meth's c2t-suffix resolution) --------
(cd "$HERE" && ./shm_meth_test.sh) || fail "shm_meth_test"
ok "shm_meth_test"

echo "ALL UNIT TESTS PASSED"
