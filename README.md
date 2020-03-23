### BWA-Mich

BWA-Mich builds upon BWA-MEM2 and includes performance improvements to the seeding and mate-rescue steps. 
It uses the Enumerated Radix Tree (ERT) index which is ~60 GB for the human genome.
BWA-Mich produces identical results as BWA-MEM2 and is 1.2-1.4x faster. 

## Getting Started
```sh
# Compile from source
git clone https://github.com/arun-sub/bwa-mem2.git ert
cd ert

# To find out vectorization features supported in your machine
cat /proc/cpuinfo

# If AVX512BW (512-bit SIMD) is supported
make clean
make -j<num_threads> arch=avx512

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

* BWA-Mich has been tested for read lengths up to 251 bp. For larger read lengths, please update READ_LEN (default = 151) in src/macro.h and rebuild the index. Note that variable read-lengths
  in the input data are supported, so READ_LEN can be conservatively set higher, and the index need not be rebuilt every time read length is changed.
* BWA-Mich requires atleast 70 GB RAM. For WGS runs on human genome (>32 threads), it is recommended to have 128-192 GB RAM.

## Performance Results

* Evaluation performed on 25 publicly available whole human genome paired-end datasets from Illumina Platinum Genomes, Illumina BaseSpaceHub and 1000 Genomes Project-Phase3.
* Human reference genome was downloaded from [here](https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta).

<p align="center">
<img src="https://github.com/arun-sub/bwa-mem2/blob/master/images/BWA-MEM2-ERT-Performance.png" height="400"/a></br>
</p>

## Citation

If you use BWA-Mich, please cite the following [paper](https://github.com/arun-sub/bwa-mem2/blob/master/doc/PaperDraft-Enumerated-Radix-Trees.pdf):

>Arun Subramaniyan, Jack Wadden, Kush Goliya, Nathan Ozog, Xiao Wu, Satish Narayanasamy, David Blaauw, Reetuparna Das. Accelerating Maximal-Exact-Match Seeding with Enumerated Radix Trees. BIORXIV/2020/003897


### BWA-MEM2 (old)
## Important Information

***Index structure has changed since commit 6743183. Rebuild the Index if you are using a later commit.***

## Getting Started
```sh
# Use precompiled binaries (recommended)
curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.0pre2/bwa-mem2-2.0pre2_x64-linux.tar.bz2 \
  | tar jxf -
bwa-mem2-2.0pre2_x64-linux/bwa-mem2 index ref.fa
bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem ref.fa read1.fq read2.fq > out.sam

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

## Usage

The usage is exactly same as the original BWA MEM tool. Here is a brief synopsys. Run ./bwa-mem2 for available commands.

```sh
# Indexing the reference sequence (Requires 28N GB memory where N is the size of the reference sequence).
./bwa-mem2 index [-p prefix] <in.fasta>
Where 
<in.fasta> is the path to reference sequence fasta file and 
<prefix> is the prefix of the names of the files that store the resultant index. Default is in.fasta.

# Mapping 
# Run "./bwa-mem2 mem" to get all options
./bwa-mem2 mem -t <num_threads> <prefix> <reads.fq/fa> > out.sam
Where <prefix> is the prefix specified when creating the index or the path to the reference fasta file in case no prefix was provided.
```

## Performance

Datasets:


 Alias	    |  Dataset source				|  No. of reads	| Read length 
 --------- | --------- | --------- | --------- 
 D1	|  Broad Institute				|  2 x 2.5M	bp	|	151bp
 D2	|  SRA: SRR7733443				|  2 x 2.5M	bp	|	151bp  
 D3	|  SRA: SRR9932168				|  2 x 2.5M	bp	|	151bp  
 D4	|  SRA: SRX6999918				|  2 x 2.5M	bp	|	151bp  



Machine details:  
Processor: Intel(R) Xeon(R) 8280 CPU @ 2.70GHz  
OS: CentOS Linux release 7.6.1810  
Memory: 100GB  


<p align="center">
<img src="https://github.com/bwa-mem2/bwa-mem2/blob/master/images/bwa-mem2-1.png" height="400"/a></br>
<img src="https://github.com/bwa-mem2/bwa-mem2/blob/master/images/bwa-mem2-2.png" height="400"/a></br>
<img src="https://github.com/bwa-mem2/bwa-mem2/blob/master/images/bwa-mem2-3.png" height="400"/a></br>
<img src="https://github.com/bwa-mem2/bwa-mem2/blob/master/images/bwa-mem2-4.png" height="400"/a></br>
</p>


## Citation

Vasimuddin Md, Sanchit Misra, Heng Li, Srinivas Aluru.
<b> Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. </b>
<i> IEEE Parallel and Distributed Processing Symposium (IPDPS), 2019. </i>
