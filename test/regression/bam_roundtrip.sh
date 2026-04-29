#!/usr/bin/env bash
# test/regression/bam_roundtrip.sh
#
# Regression: bwa-mem3 mem --bam=6 produces a BAM that decodes cleanly
# and has the same record count as the SAM path.
#
# Was: the "--bam=6 roundtrip smoke (phiX)" step inline in ci.yml.
#
# Inputs:
#   BWA_MEM3    — path to bwa-mem3 binary
#   CI_TEST_DIR — directory containing bwamem2_phix174.fa (pre-indexed),
#                 bwamem2.sam (from the phix_parity step), and FASTQs
set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${CI_TEST_DIR:?CI_TEST_DIR must be set}"

cd "$CI_TEST_DIR"
"$BWA_MEM3" mem --bam=6 bwamem2_phix174.fa \
    reads.bwa.read1.fastq.gz reads.bwa.read2.fastq.gz > bwamem2.bam
samtools quickcheck bwamem2.bam
# grep -c exits 1 on zero matches, which would abort the script under
# `set -euo pipefail` before we can report a real "0 vs 0" result —
# match the `|| true` pattern used elsewhere (thread_determinism.sh).
sam_records=$(grep -cv '^@' bwamem2.sam || true)
bam_records=$(samtools view -c bwamem2.bam)
if [ "$sam_records" != "$bam_records" ]; then
    echo "FAIL: SAM ($sam_records) vs --bam=6 BAM ($bam_records) record count mismatch"
    exit 1
fi
echo "PASS: --bam=6 roundtrip ($bam_records records)"
