### BWA-Mich

BWA-Mich builds upon BWA-MEM2 and includes performance improvements to the seeding and mate-rescue steps. 
It uses the Enumerated Radix Tree (ERT) index which is ~60 GB for the human genome.
BWA-Mich produces identical results as BWA-MEM2 and is 1.2-1.4x faster. 

## Getting Started
```sh
# Compile from source
git clone --recursive https://github.com/arun-sub/bwa-mem2.git ert
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

# Perform alignment (please specify the path prefix to the ERT index with -Z option)
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

* BWA-Mich requires atleast 70 GB RAM. For WGS runs on human genome (>32 threads), it is recommended to have 128-192 GB RAM.

## Performance Results

* Evaluation performed on 25 publicly available whole human genome paired-end datasets from Illumina Platinum Genomes, Illumina BaseSpaceHub and 1000 Genomes Project-Phase3.
* Human reference genome was downloaded from [here](https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta).

<p align="center">
<img src="https://github.com/arun-sub/bwa-mem2/blob/master/images/BWA-MEM2-ERT-Performance.png" height="400"/a></br>
</p>

## Citation

If you use BWA-Mich, please cite the following [paper](https://biorxiv.org/cgi/content/short/2020.03.23.003897v1):

> **Arun Subramaniyan, Jack Wadden, Kush Goliya, Nathan Ozog, Xiao Wu, Satish Narayanasamy, David Blaauw, Reetuparna Das. *Accelerating Maximal-Exact-Match Seeding with Enumerated Radix Trees.* https://doi.org/10.1101/2020.03.23.003897 (biorxiv). ISCA'2021 (to appear)**

