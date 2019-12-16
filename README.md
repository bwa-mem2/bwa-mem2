### BWA-MEM-M

## Getting Started
```sh
# Compile from source
git clone https://github.com/arun-sub/bwa-mem2.git ert
cd ert

# To find out vectorization features supported in your machine
cat /proc/cpuinfo

# If AVX512BW (512-bit SIMD) is supported
make clean
make -j<num_threads> arch=avx512bw

# If AVX2 (256-bit SIMD) is supported
make clean
make -j<num_threads> arch=avx2

# If SSE4.1 (128-bit SIMD) is supported (default)
make -j<num_threads>

# Build index (Takes ~2 hr for human genome with 56 threads. 1 hr for BWT, 1 hr for ERT)
./bwa-mem2 index -a ert -t <num_threads> -p <index prefix> <input.fasta>

# Perform alignment
./bwa-mem2 mem -Y -K 100000000 -t <num_threads> -Z <index prefix> <input_1.fastq> <input_2.fastq> -o <output_ert.sam>

# To verify output with BWA-MEM
git clone https://github.com/lh3/bwa.git
cd bwa
make

# Perform alignment
./bwa mem -Y -K 100000000 -t <num_threads> <index prefix> <input_1.fastq> <input_2.fastq> -o <output_mem.sam>

# Compare output SAM files
diff <output_mem.sam> <output_ert.sam>

# To diff large SAM files use https://github.com/unhammer/diff-large-files

```

## Notes

* BWA-MEM-M has been tested for read lengths up to 151 bp. For larger read lengths, please update READ_LEN in src/macro.h
* BWA-MEM-M requires atleast 64 GB RAM. For WGS runs on human genome, it is recommended to have 128-192 GB RAM.


### BWA-MEM2 (old)

## Getting Started
```sh
# Use precompiled binaries (recommended)
curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.0pre1/bwa-mem2-2.0pre1_x64-linux.tar.bz2 \
  | tar jxf -
bwa-mem2-2.0pre1_x64-linux/bwa-mem2 index ref.fa
bwa-mem2-2.0pre1_x64-linux/bwa-mem2 mem ref.fa read1.fq read2.fq > out.sam

# Compile from source (not recommended for general users)
git clone https://github.com/bwa-mem2/bwa-mem2
cd bwa-mem2
make
./bwa-mem2
```

## Introduction

Bwa-mem2 is the next version of the bwa-mem algorithm in [bwa][bwa]. It
produces alignment identical to bwa and is ~80% faster.

The original bwa was developed by Heng Li (@lh3). Performance enhancement in
bwa-mem2 was primarily done by Vasimuddin Md (@yuk12) and Sanchit Misra (@sanchit-misra)
from Parallel Computing Lab, Intel.
Bwa-mem2 is distributed under the MIT license.

## Installation

For general users, it is recommended to use the precompiled binaries from the
[release page][rel]. These binaries were compiled with the Intel compiler and
runs faster than gcc-compiled binaries. The precompiled binaries also
indirectly support CPU dispatch. The `bwa-mem2` binary can automatically choose
the most efficient implementation based on the SIMD instruction set available
on the running machine. Precompiled binaries were generated on a CentOS6
machine using the following command line:
```sh
make CXX=icpc multi
```

[bwa]: https://github.com/lh3/bwa
[rel]: https://github.com/bwa-mem2/bwa-mem2/releases

## Citation

Vasimuddin Md, Sanchit Misra, Heng Li, Srinivas Aluru.
<b> Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. </b>
<i> IEEE Parallel and Distributed Processing Symposium (IPDPS), 2019. </i>
