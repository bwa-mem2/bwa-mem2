#!/usr/bin/env bash
# test/regression/phix_parity.sh
#
# Regression: bwa and bwa-mem3 produce identical SAM output on phiX174
# (500 PE reads, simulated with dwgsim -z 42 for reproducibility).
#
# Was: the "Compare bwa and bwa-mem3 output" step inline in ci.yml.
#
# Inputs (read from environment):
#   BWA_MEM3    — path to the bwa-mem3 binary under test
#   CI_TEST_DIR — working directory containing phix174.fa and dwgsim-simulated
#                 reads.bwa.read[12].fastq.gz (produced by upstream CI steps)
#
# Exits 0 on parity, 1 on diff.
set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${CI_TEST_DIR:?CI_TEST_DIR must be set}"

cd "$CI_TEST_DIR"

cp phix174.fa bwa_phix174.fa
bwa index bwa_phix174.fa
bwa mem bwa_phix174.fa reads.bwa.read1.fastq.gz reads.bwa.read2.fastq.gz \
    > bwa.sam 2>bwa.log

cp phix174.fa bwamem2_phix174.fa
"$BWA_MEM3" index bwamem2_phix174.fa
"$BWA_MEM3" mem bwamem2_phix174.fa \
    reads.bwa.read1.fastq.gz reads.bwa.read2.fastq.gz \
    > bwamem2.sam 2>bwamem2.log

normalize() {
    grep -v '^@PG' "$1" | grep -v '^@HD' \
        | sed 's/\tMQ:i:[0-9]*//' \
        | sed 's/\tHN:i:[0-9]*//'
}
normalize bwa.sam     | sort > bwa.normalized.sam
normalize bwamem2.sam | sort > bwamem2.normalized.sam

echo "bwa aligned records:      $(grep -cv '^@' bwa.normalized.sam)"
echo "bwa-mem3 aligned records: $(grep -cv '^@' bwamem2.normalized.sam)"

if diff bwa.normalized.sam bwamem2.normalized.sam > /dev/null 2>&1; then
    echo "PASS: bwa and bwa-mem3 output is identical (after normalizing headers and MQ tag)"
else
    echo "FAIL: bwa and bwa-mem3 output differs"
    diff bwa.normalized.sam bwamem2.normalized.sam | head -40
    exit 1
fi
