#!/usr/bin/env bash
# Reproducible PGO+LTO build for bwa-mem2 (gcc). ~11-13% faster than stock gcc -O3,
# output bit-identical. Run from a clean bwa-mem2 source tree that already has the
# two code fixes applied (bwa-mem2-code-fixes.patch).
#
# Usage: build_pgo.sh <ref.fa> <train_1.fq> <train_2.fq>   (run inside the source tree)
set -euo pipefail
REF="$1"; T1="$2"; T2="$3"
PROF="$(pwd)/.pgoprof"; mkdir -p "$PROF"
GEN="-fprofile-generate=$PROF -fprofile-update=prefer-atomic"
USE="-fprofile-use=$PROF -fprofile-correction -flto=auto -funroll-loops"

echo "[1/3] instrumented build"
make clean >/dev/null 2>&1; (cd ext/safestringlib && make clean >/dev/null 2>&1) || true
make arch=native CXX=g++ -j"$(nproc)" \
     CXXFLAGS="-g -O3 $GEN -fpermissive -march=native" LDFLAGS="$GEN -O3"

echo "[2/3] training run (index must already exist for \$REF)"
[ -f "$REF.bwt.2bit.64" ] || ./bwa-mem2 index "$REF"
./bwa-mem2 mem -t "$(nproc)" "$REF" "$T1" "$T2" > /dev/null

echo "[3/3] optimized rebuild with profile"
make clean >/dev/null 2>&1; (cd ext/safestringlib && make clean >/dev/null 2>&1) || true
make arch=native CXX=g++ -j"$(nproc)" \
     CXXFLAGS="-g -O3 $USE -fpermissive -march=native" LDFLAGS="$USE -O3"
echo "done -> ./bwa-mem2  (PGO+LTO+native)"
