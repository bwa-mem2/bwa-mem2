
# bwa-mem2
----------------------------------------------------------------------------------
IMPORTANT:  
1.'bwa-mem2 index' requires ~48GB of space for full human genome plus index.  
2. We recommend to run the code in AVX512 mode to get maximum performance.  
3. The command-line of bwa-mem2 is exactly same as bwa-mem.  
----------------------------------------------------------------------------------


Compile:  
$ ```make```  
'make' compiles using SSE4.1 vector mode by default.

Compile options:  
```make CXX=<compiler>```, e.g compiler can 'g++' or 'icpc'  
```make arch=<mode>```  

mode:  
a. native - generates 'bwa-mem2' binary for the architechture on which it is compiled. It detects the underlying harware flags for AVX512/AVX2/SSE2 vector modes and compiles accordingly. 
For example, Intel Xeon Haswell supports AVX2 vector mode; Intel Xeon Skylake supports AVX512 vector mode. AVX vector support by the processor can be checked in /proc/cpuinfo file.  
b. sse - generates 'bwa-mem2' SSE4.1 vector code binary  
c. avx2 - generates 'bwa-mem2' AVX2 vector code binary  
d. avx512bw - generates 'bwa-mem2' AVX512BW vector code binary  


 ```make multi```  
Generates binaries for all the three vector modes: bwa-mem2.sse41, bwa-mem2.avx2, bwa-mem2.avx512bw and also generates a binary called bwa-mem2 that, if run, detects the underlying architecture and runs the corresponding binary out of bwa-mem2.sse41, bwa-mem2.avx2, bwa-mem2.avx512bw.  


Run:   
   - Create the index: 'bwa-mem2 index <reference file>'  
   - Mapping: 'bwa-mem2 mem' can be run with appropriate paramters  
   e.g. bwa-mem2 mem -t 1 <index> <read1.fq> [<read2.fq>] > <out.sam>  

Notes:  
- In SSE4.1 vector mode, `bwa-mem2` uses bwt index with 2-bit representation  
- In AVX512/AVX2 vector mode, `bwa-mem2` uses bwt index with 8-bit representation  


Example run:  
```sh bwa-mem2 index datasets/hg38.fa```   
```bwa-mem2 mem -t 1 datasets/hg38.fa datasets/SRR099967.filt.fastq > datasets/aln-se.sa```



NUMA:  
    In multi-socket system, it is beneficial to bind the process memory to a socket  
    Linux command `numactl` can be used to bind process memory to a given numa domain  
    For exmaple, the following command binds the memory to socket-0  
``` $ numactl -m 0 bwa-mem2 mem -t 10 -o datasets/aln-se.sa datasets/hg38.fa datasets/SRR099967.filt.fastq```  

Similarly, the process threads can be bound to cores.  
```$ numactl -m 0 -C 0-9 bwa-mem2 mem -t 10 -o datasets/aln-se.sa datasets/hg38.fa datasets/SRR099967.filt.fastq```  

```$ lscpu``` - Linux command provides the architectural details such as list of NUMA domains and their CPU lists  


Directory structure (bwa-mem2):  
- src/ - contains source code files  
- test/ - contains test files to run individual benchmarks for smem, sal, bsw (disabled for now)  
- ext/ - contains external Intel safestring library for secure string operations  
- images/ - contains bwa-mem2 perfomance charts, demonstrating its speedup against original bwa-mem  


Notes:  
- As mentioned previously, the code supports four exceution modes:  
  1. Vector SSE2: It runs bwa-mem2 with the kernels running SSE2 vector instructions (128 bits vector-width)  
  2. Vector AVX2: It runs bwa-mem2 with the kernels running AVX2 vector instructions (256 bits vector-width)  
  3. Vector AVX512: It runs bwa-mem2 with the kernels running AVX512 vector instructions (512 bits vector-width)  

- BWA-MEM2 displays run-time performance profiling of the code on standard output. (Hopefully, we managed it to be self-explanatory)  

----------------------------------------------------------------------------------
