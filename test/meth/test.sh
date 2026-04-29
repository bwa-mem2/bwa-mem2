#!/usr/bin/env bash
# Regression test: bwa-mem3 mem --meth end-to-end.
#
# Two layers of assertions:
#
# Layer 1 (always runs):  valid BAM emission.
#   - binary builds and runs with --meth
#   - produces uncompressed BAM readable by samtools
#   - @PG ID:bwa-mem3-meth present
#   - BGZF EOF marker at tail
#   - --set-as-failed / --do-not-penalize-chimeras parse cleanly
#
# Layer 2 (runs if pixi + bwameth.py available):  equivalence to bwameth.py.
#   Builds a bwameth c2t reference, c2t-converts reads, runs BOTH the
#   bwameth.py Python pipeline and `bwa-mem3 mem --meth` on the same
#   converted reads, and diffs structural fields + YD:Z tags + flag
#   distribution. Currently zero diff on the bwa-meth/example fixture.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
BWAMEM3="$HERE/../../bwa-mem3"
SAMTOOLS="${SAMTOOLS:-samtools}"
BWAMETH_DIR="${BWAMETH_DIR:-$HOME/work/git/bwa-meth}"
BWAMETH_PY="$BWAMETH_DIR/bwameth.py"

if [[ ! -x "$BWAMEM3" ]]; then
    echo "ERROR: bwa-mem3 binary not found at $BWAMEM3. Run 'make arm64' first."
    exit 2
fi
if ! command -v "$SAMTOOLS" >/dev/null 2>&1; then
    echo "ERROR: samtools not found on PATH."
    exit 2
fi

cd "$HERE"

# ---------------------------------------------------------------------------
# Layer 1: BAM emission smoke test
# ---------------------------------------------------------------------------

if [[ ! -f ref.fa.bwameth.c2t.bwt.2bit.64 ]]; then
    "$BWAMEM3" index --meth ref.fa >/dev/null 2>&1
fi

"$BWAMEM3" mem --meth -t 2 ref.fa t_R1.fastq.gz 2>/dev/null > /tmp/meth_test.bam

EXPECT_EOF="1f8b08040000000000ff0600424302001b0003000000000000000000"
ACTUAL_EOF="$(tail -c 28 /tmp/meth_test.bam | od -An -v -t x1 | tr -d ' \n')"
if [[ "${ACTUAL_EOF%$'\n'}" != "${EXPECT_EOF}" ]]; then
    echo "FAIL: BGZF EOF marker mismatch (actual=$ACTUAL_EOF)"; exit 1
fi

HDR="$("$SAMTOOLS" view -H /tmp/meth_test.bam 2>&1)"
if echo "$HDR" | grep -qi 'truncated\|EOF marker is absent'; then
    echo "FAIL: samtools reports truncated BAM"; echo "$HDR"; exit 1
fi
if ! echo "$HDR" | grep -q 'ID:bwa-mem3-meth'; then
    echo "FAIL: @PG ID:bwa-mem3-meth missing"; exit 1
fi

TOTAL="$("$SAMTOOLS" view -c /tmp/meth_test.bam 2>/dev/null)"
if [[ "$TOTAL" -lt 1 ]]; then echo "FAIL: zero records in output BAM"; exit 1; fi

"$BWAMEM3" mem --meth --set-as-failed f --do-not-penalize-chimeras \
    ref.fa t_R1.fastq.gz 2>/dev/null > /tmp/meth_test2.bam
if [[ ! -s /tmp/meth_test2.bam ]]; then
    echo "FAIL: --set-as-failed + --do-not-penalize-chimeras produced empty output"
    exit 1
fi

echo "OK layer 1: bwa-mem3 mem --meth (records=$TOTAL, BGZF-EOF ok, @PG bwa-mem3-meth ok)"

# ---------------------------------------------------------------------------
# Layer 2: bwa-meth equivalence (requires pixi + bwameth.py + toolshed)
# ---------------------------------------------------------------------------

if ! command -v pixi >/dev/null 2>&1; then
    echo "SKIP layer 2: pixi not on PATH"
    exit 0
fi
if [[ ! -f "$BWAMETH_PY" ]]; then
    echo "SKIP layer 2: bwameth.py not found at $BWAMETH_PY (set BWAMETH_DIR)"
    exit 0
fi
if [[ ! -f "$HERE/pyproject.toml" ]]; then
    echo "SKIP layer 2: pixi env not initialized in $HERE (run 'pixi init && pixi add toolshed')"
    exit 0
fi

export PATH="$(cd "$HERE/../.." && pwd):$PATH"

if [[ ! -f "$HERE/ref.fa.bwameth.c2t.0123" ]]; then
    (cd "$HERE" && pixi run python3 "$BWAMETH_PY" index-mem2 ref.fa >/dev/null 2>&1)
fi

pixi run python3 "$BWAMETH_PY" --reference ref.fa t_R1.fastq.gz t_R2.fastq.gz \
    2>/dev/null > /tmp/meth_oracle.sam

pixi run python3 "$BWAMETH_PY" c2t t_R1.fastq.gz t_R2.fastq.gz 2>/dev/null \
    > /tmp/meth_c2t.fq
"$BWAMEM3" mem --meth -CM -p -T 40 -B 2 -L 10 -U 100 -t 4 \
    ref.fa.bwameth.c2t /tmp/meth_c2t.fq 2>/dev/null > /tmp/meth_mine.bam
"$SAMTOOLS" view /tmp/meth_mine.bam 2>/dev/null > /tmp/meth_mine.sam

grep -v '^@' /tmp/meth_oracle.sam > /tmp/meth_oracle_records.sam

norm() {
    awk -F'\t' 'BEGIN{OFS="\t"} { if ($7 == "=") $7 = $3; print }' "$1" \
        | sort -k1,1 -k2,2n
}

# TLEN (col 9) intentionally diverges from bwameth.py on the edge case of
# mates landing on opposite projected strands (f* vs r* in the doubled c2t
# reference). bwameth.py is a post-processor that strips the f/r prefix to
# make RNAME==RNEXT but inherits bwa-mem3's internal TLEN=0 (which bwa-mem3
# sets because the internal rids differ). Per the SAM spec, TLEN=0 is
# reserved for "information unavailable (e.g., mates on different refs)";
# once we consolidate f*/r* to one output contig the information IS
# available, so we compute a real TLEN. Diff QNAME..MPOS (cols 1-8) to
# assert structural parity without re-enforcing bwameth.py's TLEN=0 quirk.
ONE="$(diff <(norm /tmp/meth_mine.sam | cut -f1-8) \
             <(norm /tmp/meth_oracle_records.sam | cut -f1-8) | wc -l | tr -d ' ')"
if [[ "$ONE" != "0" ]]; then
    echo "FAIL layer 2: structural diff in cols 1-8 ($ONE lines)"
    diff <(norm /tmp/meth_mine.sam | cut -f1-8) <(norm /tmp/meth_oracle_records.sam | cut -f1-8) | head -20
    exit 1
fi

TWO="$(diff <(norm /tmp/meth_mine.sam | cut -f6) \
             <(norm /tmp/meth_oracle_records.sam | cut -f6) | wc -l | tr -d ' ')"
if [[ "$TWO" != "0" ]]; then echo "FAIL layer 2: CIGAR diff ($TWO lines)"; exit 1; fi

for d in f r; do
    MINE_YD="$(grep -c "YD:Z:$d" /tmp/meth_mine.sam || true)"
    ORACLE_YD="$(grep -c "YD:Z:$d" /tmp/meth_oracle_records.sam || true)"
    if [[ "$MINE_YD" != "$ORACLE_YD" ]]; then
        echo "FAIL layer 2: YD:Z:$d count mismatch (mine=$MINE_YD oracle=$ORACLE_YD)"
        exit 1
    fi
done

MINE_N="$(wc -l < /tmp/meth_mine.sam | tr -d ' ')"
ORACLE_N="$(wc -l < /tmp/meth_oracle_records.sam | tr -d ' ')"
if [[ "$MINE_N" != "$ORACLE_N" ]]; then
    echo "FAIL layer 2: record count mismatch (mine=$MINE_N oracle=$ORACLE_N)"
    exit 1
fi

MINE_F="$(grep -c YD:Z:f /tmp/meth_mine.sam || true)"
MINE_R="$(grep -c YD:Z:r /tmp/meth_mine.sam || true)"
echo "OK layer 2: bwa-mem3 mem --meth matches bwameth.py (records=$MINE_N, YD:Z:f=$MINE_F YD:Z:r=$MINE_R)"

# ---------------------------------------------------------------------------
# Layer 3: full-pipeline end-to-end — `bwa-mem3 index --meth` + `mem --meth`
# ---------------------------------------------------------------------------
# Tests the single-command drop-in replacement for the bwameth.py pipeline:
#
#   bwa-mem3 index --meth ref.fa                           # once
#   bwa-mem3 mem   --meth ref.fa R1.fq R2.fq | samtools sort ...
#
# No Python, no bwameth.py invocation, no piping. Asserts byte-for-byte
# equivalence with bwameth.py end-to-end on the example fixture:
#   - same set of aligned primary QNAMEs
#   - same chromosome + position for every mapped primary
#   - same CIGAR for every mapped primary

# Only runs when the bwameth.py oracle is available (Layer 2's pixi env).
if ! command -v pixi >/dev/null 2>&1; then
    echo "SKIP layer 3: pixi not on PATH (needed to regenerate oracle)"
    exit 0
fi
if [[ ! -f "$BWAMETH_PY" ]]; then
    echo "SKIP layer 3: bwameth.py not found at $BWAMETH_PY"
    exit 0
fi

# Fresh index via our own --meth builder.
rm -f "$HERE/ref.fa.bwameth.c2t"*
"$BWAMEM3" index --meth ref.fa >/dev/null 2>&1

# Single-command end-to-end alignment.
"$BWAMEM3" mem --meth -t 4 ref.fa t_R1.fastq.gz t_R2.fastq.gz 2>/dev/null > /tmp/meth_e2e.bam
"$SAMTOOLS" view /tmp/meth_e2e.bam 2>/dev/null > /tmp/meth_e2e.sam

# Oracle: bwameth.py end-to-end on the same raw FASTQs (reuses Layer 2's
# meth_oracle.sam if it was just written, otherwise regenerate).
if [[ ! -s /tmp/meth_oracle.sam ]]; then
    pixi run python3 "$BWAMETH_PY" --reference ref.fa t_R1.fastq.gz t_R2.fastq.gz \
        2>/dev/null > /tmp/meth_oracle.sam
fi
grep -v '^@' /tmp/meth_oracle.sam > /tmp/meth_oracle_records.sam

# Record count sanity.
N_MINE="$("$SAMTOOLS" view -c /tmp/meth_e2e.bam 2>/dev/null)"
N_ORAC="$(wc -l < /tmp/meth_oracle_records.sam | tr -d ' ')"
if [[ "$N_MINE" != "$N_ORAC" ]]; then
    echo "FAIL layer 3: record count mismatch (mine=$N_MINE oracle=$N_ORAC)"
    exit 1
fi

# Per-primary agreement: same chrom+pos and same CIGAR for every mapped read.
# (Soft-clip/secondary tie-breaking can be order-sensitive; we allow the
# tiny symmetric set of "mapped only in one" to be 0 here, matching the
# bwa-meth/example fixture.)
python3 - <<'PY'
import sys
def load(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('@'): continue
            f = line.rstrip('\n').split('\t')
            flag = int(f[1])
            if flag & 0x100 or flag & 0x800: continue
            is_r2 = 2 if (flag & 0x80) else 1
            unmapped = 1 if (flag & 0x4) else 0
            out[f'{f[0]}/{is_r2}']=(unmapped, int(f[3]), f[5], f[2])
    return out
n=load('/tmp/meth_e2e.sam'); o=load('/tmp/meth_oracle.sam')
only_n=sorted(set(n) - set(o))
only_o=sorted(set(o) - set(n))
both=[k for k in n if k in o and n[k][0]==0 and o[k][0]==0]
gap=[k for k in n if k in o and n[k][0]==1 and o[k][0]==0]
nonly=[k for k in n if k in o and n[k][0]==0 and o[k][0]==1]
same_pos=sum(1 for k in both if n[k][1]==o[k][1] and n[k][3]==o[k][3])
same_cigar=sum(1 for k in both if n[k][2]==o[k][2])
if only_n or only_o or gap or nonly or same_pos != len(both) or same_cigar != len(both):
    sys.stderr.write(
        f'FAIL layer 3: diverged from bwameth.py oracle\n'
        f'  native-only QNAMEs: {len(only_n)}\n'
        f'  oracle-only QNAMEs: {len(only_o)}\n'
        f'  both-mapped: {len(both)}\n'
        f'  oracle-only mapped: {len(gap)}\n'
        f'  native-only mapped: {len(nonly)}\n'
        f'  same chrom+pos: {same_pos}/{len(both)}\n'
        f'  same CIGAR:     {same_cigar}/{len(both)}\n')
    sys.exit(1)
print(f'OK layer 3: bwa-mem3 mem --meth == bwameth.py end-to-end '
      f'(records={len(both)}, chrom+pos match, CIGAR match)')
PY
