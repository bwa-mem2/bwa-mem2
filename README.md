[![GitHub Downloads](https://img.shields.io/github/downloads/bwa-mem2/bwa-mem2/total?label=GitHub%20Downloads)](https://github.com/bwa-mem2/bwa-mem2/releases)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/bwa-mem2?label=BioConda%20Installs)](https://anaconda.org/bioconda/bwa-mem2)
[![European Galaxy server](https://img.shields.io/badge/usegalaxy-.eu-brightgreen?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAACXBIWXMAAAsTAAALEwEAmpwYAAACC2lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjE8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlBob3RvbWV0cmljSW50ZXJwcmV0YXRpb24+MjwvdGlmZjpQaG90b21ldHJpY0ludGVycHJldGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KD0UqkwAAAn9JREFUOBGlVEuLE0EQruqZiftwDz4QYT1IYM8eFkHFw/4HYX+GB3/B4l/YP+CP8OBNTwpCwFMQXAQPKtnsg5nJZpKdni6/6kzHvAYDFtRUT71f3UwAEbkLch9ogQxcBwRKMfAnM1/CBwgrbxkgPAYqlBOy1jfovlaPsEiWPROZmqmZKKzOYCJb/AbdYLso9/9B6GppBRqCrjSYYaquZq20EUKAzVpjo1FzWRDVrNay6C/HDxT92wXrAVCH3ASqq5VqEtv1WZ13Mdwf8LFyyKECNbgHHAObWhScf4Wnj9CbQpPzWYU3UFoX3qkhlG8AY2BTQt5/EA7qaEPQsgGLWied0A8VKrHAsCC1eJ6EFoUd1v6GoPOaRAtDPViUr/wPzkIFV9AaAZGtYB568VyJfijV+ZBzlVZJ3W7XHB2RESGe4opXIGzRTdjcAupOK09RA6kzr1NTrTj7V1ugM4VgPGWEw+e39CxO6JUw5XhhKihmaDacU2GiR0Ohcc4cZ+Kq3AjlEnEeRSazLs6/9b/kh4eTC+hngE3QQD7Yyclxsrf3cpxsPXn+cFdenF9aqlBXMXaDiEyfyfawBz2RqC/O9WF1ysacOpytlUSoqNrtfbS642+4D4CS9V3xb4u8P/ACI4O810efRu6KsC0QnjHJGaq4IOGUjWTo/YDZDB3xSIxcGyNlWcTucb4T3in/3IaueNrZyX0lGOrWndstOr+w21UlVFokILjJLFhPukbVY8OmwNQ3nZgNJNmKDccusSb4UIe+gtkI+9/bSLJDjqn763f5CQ5TLApmICkqwR0QnUPKZFIUnoozWcQuRbC0Km02knj0tPYx63furGs3x/iPnz83zJDVNtdP3QAAAABJRU5ErkJggg==)](https://usegalaxy.eu/root?tool_id=bwa_mem2)
[![US Galaxy server](https://img.shields.io/badge/usegalaxy-.org-orange?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAACXBIWXMAAAsTAAALEwEAmpwYAAACC2lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjE8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlBob3RvbWV0cmljSW50ZXJwcmV0YXRpb24+MjwvdGlmZjpQaG90b21ldHJpY0ludGVycHJldGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KD0UqkwAAAn9JREFUOBGlVEuLE0EQruqZiftwDz4QYT1IYM8eFkHFw/4HYX+GB3/B4l/YP+CP8OBNTwpCwFMQXAQPKtnsg5nJZpKdni6/6kzHvAYDFtRUT71f3UwAEbkLch9ogQxcBwRKMfAnM1/CBwgrbxkgPAYqlBOy1jfovlaPsEiWPROZmqmZKKzOYCJb/AbdYLso9/9B6GppBRqCrjSYYaquZq20EUKAzVpjo1FzWRDVrNay6C/HDxT92wXrAVCH3ASqq5VqEtv1WZ13Mdwf8LFyyKECNbgHHAObWhScf4Wnj9CbQpPzWYU3UFoX3qkhlG8AY2BTQt5/EA7qaEPQsgGLWied0A8VKrHAsCC1eJ6EFoUd1v6GoPOaRAtDPViUr/wPzkIFV9AaAZGtYB568VyJfijV+ZBzlVZJ3W7XHB2RESGe4opXIGzRTdjcAupOK09RA6kzr1NTrTj7V1ugM4VgPGWEw+e39CxO6JUw5XhhKihmaDacU2GiR0Ohcc4cZ+Kq3AjlEnEeRSazLs6/9b/kh4eTC+hngE3QQD7Yyclxsrf3cpxsPXn+cFdenF9aqlBXMXaDiEyfyfawBz2RqC/O9WF1ysacOpytlUSoqNrtfbS642+4D4CS9V3xb4u8P/ACI4O810efRu6KsC0QnjHJGaq4IOGUjWTo/YDZDB3xSIxcGyNlWcTucb4T3in/3IaueNrZyX0lGOrWndstOr+w21UlVFokILjJLFhPukbVY8OmwNQ3nZgNJNmKDccusSb4UIe+gtkI+9/bSLJDjqn763f5CQ5TLApmICkqwR0QnUPKZFIUnoozWcQuRbC0Km02knj0tPYx63furGs3x/iPnz83zJDVNtdP3QAAAABJRU5ErkJggg==)](https://usegalaxy.org/root?tool_id=bwa_mem2)

## Important Information

***We are happy to announce that the index size on disk is down by 8 times and in memory by 4 times due to moving to only one type of FM-index (2bit.64 instead of 2bit.64 and 8bit.32) and 8x compression of suffix array. For example, for human genome, index size on disk is down to ~10GB from ~80GB and memory footprint is down to ~10GB from ~40GB.***
***There is a substantial reduction in index IO time due to the reduction and hardly any performance impact on read mapping.***
***Due to this change in index structure (in commit #4b59796, 10th October 2020), you will need to rebuild the index.***

***Added MC flag in the output sam file in commit a591e22. Output should match original bwa-mem version 0.7.17.***

***As of commit e0ac59e, we have a git submodule safestringlib. To get it, use --recursive while cloning or use "git submodule init" and "git submodule update" in an already cloned repository (See below for more details).***


## Getting Started
```sh
# Use precompiled binaries (recommended)
curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 \
  | tar jxf -
bwa-mem2-2.2.1_x64-linux/bwa-mem2 index ref.fa
bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem ref.fa read1.fq read2.fq > out.sam

# Compile from source (not recommended for general users)
# Get the source
git clone --recursive https://github.com/bwa-mem2/bwa-mem2
cd bwa-mem2
# Or
git clone https://github.com/bwa-mem2/bwa-mem2
cd bwa-mem2
git submodule init
git submodule update
# Compile and run
make
./bwa-mem2
```

## Introduction

The tool bwa-mem2 is the next version of the bwa-mem algorithm in [bwa][bwa]. It
produces alignment identical to bwa and is ~1.3-3.1x faster depending on the use-case, dataset and the running machine.

The original bwa was developed by Heng Li (@lh3). Performance enhancement in
bwa-mem2 was primarily done by Vasimuddin Md (@yuk12) and Sanchit Misra (@sanchit-misra)
from Parallel Computing Lab, Intel.
bwa-mem2 is distributed under the MIT license.

## Installation

For general users, it is recommended to use the precompiled binaries from the
[release page][rel]. These binaries were compiled with the Intel compiler and
runs faster than gcc-compiled binaries. The precompiled binaries also
indirectly support CPU dispatch. The `bwa-mem2` binary can automatically choose
the most efficient implementation based on the SIMD instruction set available
on the running machine. Precompiled binaries were generated on a CentOS7
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
Reference Genome: human_g1k_v37.fasta

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


We followed the steps below to collect the performance results:  
A. Data download steps:
1. Download SRA toolkit from https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software#header-global    
2. tar xfzv sratoolkit.2.10.5-centos_linux64.tar.gz  
3. Download D2: sratoolkit.2.10.5-centos_linux64/bin/fastq-dump --split-files SRR7733443   
4. Download D3: sratoolkit.2.10.5-centos_linux64/bin/fastq-dump --split-files SRR9932168   
5. Download D4: sratoolkit.2.10.5-centos_linux64/bin/fastq-dump --split-files SRX6999918   



B. Alignment steps:   
1. git clone https://github.com/bwa-mem2/bwa-mem2.git   
2. cd bwa-mem2   
3. ```make CXX=icpc``` (using intel C/C++ compiler)   
or   ```make``` (using gcc compiler)   
4. ./bwa-mem2 index <ref.fa>   
5. ./bwa-mem2 mem [-t <#threads>] <ref.fa> <in_1.fastq> [<in_2.fastq>]  >  <output.sam>   

For example,  in our double socket (56 threads each) and double numa compute node, we used the following command line to align D2 to human_g1k_v37.fasta reference genome.  
```
numactl -m 0 -C 0-27,56-83 ./bwa-mem2 index human_g1k_v37.fasta  
numactl -m 0 -C 0-27,56-83 ./bwa-mem2 mem -t 56 human_g1k_v37.fasta SRR7733443_1.fastq SRR7733443_2.fastq > d2_align.sam
```

<p align="center">
<img src="https://github.com/bwa-mem2/bwa-mem2/blob/master/images/bwa-mem2-1.png" height="400"/a></br>
<img src="https://github.com/bwa-mem2/bwa-mem2/blob/master/images/bwa-mem2-2.png" height="400"/a></br>
<img src="https://github.com/bwa-mem2/bwa-mem2/blob/master/images/bwa-mem2-3.png" height="400"/a></br>
<img src="https://github.com/bwa-mem2/bwa-mem2/blob/master/images/bwa-mem2-4.png" height="400"/a></br>
</p> 

## bwa-mem2 seeding phase accelerated using LISA (Learned-Indexes for Sequence Analysis)

bwa-mem2-lisa is an accelerated version of bwa-mem2 where we apply learned-indexes to the seeding phase. bwa-mem2-lisa branch contains the source code of the implementation. Following are the features of bwa-mem2-lisa:
1. Exact same output as bwa-mem2.
2. All command-lines for creating an index and the read mapping are exactly same as bwa-mem2.
3. bwa-mem2-lisa accelerates seeding phase (one of the major bottlenecks in bwa-mem2) by up to 4.5x compared to bwa-mem2.
4. The memory footprint of bwa-mem2-lisa index is ~120GB for human genome.
5. The code is present in bwa-mem2-lisa branch: https://github.com/bwa-mem2/bwa-mem2/tree/bwa-mem2-lisa


## bwa-mem2 seeding speedup with Enumerated Radix Trees (Code in ert branch)

The ert branch of bwa-mem2 repository contains codebase of enuerated radix tree based acceleration of bwa-mem2. The ert code is built on the top of bwa-mem2 (thanks to the hard work by @arun-sub). 
The following are the highlights of the ert based bwa-mem2 tool: 
1. Exact same output as bwa-mem(2) 
2. The tool has two additional flags to enable the use of ert solution (for index creation and mapping), else it runs in vanilla bwa-mem2 mode 
3. It uses 1 additional flag to create ert index (different from bwa-mem2 index) and 1 additional flag for using that ert index (please see the readme of ert branch) 
4. The ert solution is 10% - 30% faster (tested on above machine configuration) in comparison to vanilla bwa-mem2 -- users are adviced to use option `-K 1000000` to see the speedups 
5. The memory foot print of the ert index is ~60GB 
6. The code is present in ert branch: https://github.com/bwa-mem2/bwa-mem2/tree/ert


## Citation

Vasimuddin Md, Sanchit Misra, Heng Li, Srinivas Aluru.
<b> Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. </b>
<i> IEEE Parallel and Distributed Processing Symposium (IPDPS), 2019. [10.1109/IPDPS.2019.00041](https://doi.org/10.1109/IPDPS.2019.00041) </i>
