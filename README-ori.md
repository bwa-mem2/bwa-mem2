
# bwa-mem2
----------------------------------------------------------------------------------
IMPORTANT:  
1.'bwa-mem2 index' requires ~48GB of space for full genome plus index.  
2. We recommend to run the code in AVX512 mode to get maximum performance.  
3. The command-line of bwa-mem2 is exactly same as bwa-mem.  
----------------------------------------------------------------------------------


Compile:  
$ ```make```  
'make' detects the underlying harware flags for AVX512/AVX2/SSE2 vector modes and compiles accordingly.
If the platform does not support any AVX/SSE vector mode then it compiles the code in fully scalar mode.
For example, Intel Xeon Haswell supports AVX2 vector mode; Intel Xeon Skylake supports AVX512
vector mode. AVX vector support by the processor can be check in /proc/cpuinfo file.  

Compile options:  
```make CXX=<compiler>```, e.g compiler can 'g++' or 'icpc'  
```make arch=<mode>```  

mode:  
a. native - generates 'bwa-mem2' binary for platform native architechture  
b. sse - generates 'bwa-mem2' SSE vector code binary  
c. avx2 - generates 'bwa-mem2' avx2 vector code binary  
d. avx512bw - generates 'bwa-mem2' avx512bw vector code binary  

 ```make multi```  
Generates binaries for all the three vector modes: bwa-mem2.sse41, bwa-mem2.avx2, bwa-mem2.avx512bw.  


Run:   
   - Create the index: 'bwa-mem2 index <reference file>'  
   - 'bwa-mem2 mem' can be run with appropriate paramters  
   e.g. bwa-mem2 mem -t 1 <idx> <read1.fq> [<read2.fq>] > <out.sam>  

Notes:  
- In SSE2 vector mode, `bwa-mem2 index` creates bwt index with 2-bit representation  
- In AVX512/AVX2 vector mode, 'bwa-mem2 index' creates bwt index with 8-bit representation  


Example run:  
```sh bwa-mem2 index datasets/hg38.fa```   
```bwa-mem2 mem -t 1 datasets/hg38.fa datasets/SRR099967.filt.fastq > datasets/aln-se.sa```



NUMA:  
    In multi-socket system, it is beneficial to bind the process memory to a socket  
    Linux command `numactl` can be used to bind process memory to a given numa domain  
    For exmaple, the following command binds the memory to socket-0  

``` $ numactl -m 0 bwa-mem2 mem -t 10 -o datasets/aln-se.sa datasets/hg38.fa datasets/SRR099967.filt.fastq```  
Similarly, the process threads can be bind to cores.  
```$ numactl -m 0 -C 0-10 bwa-mem2 mem -t 10 -o datasets/aln-se.sa datasets/hg38.fa datasets/SRR099967.filt.fastq```  

```$ lscpu``` - Linux command provides the architectural details such as list of NUMAs domains and their CPU lists  


Directory structure (bwa-mem2):  
- src/ - contains source code files  
- test/ - contains test files to run individual benchmarks for smem, sal, bsw (disabled for now)  
- ext/ - contains external Intel safestring library for secure string operations  
- images/ - contains bwa-mem2 perofmrnace charts, demonstrating its speedup against original bwa-mem  


Notes:  
- As mentioned previously, the code supports four exceution modes:  
  1. Scalar: It runs bwa-mem2 with the kernels running with scalar instructions  
  2. Vector SSE2: It runs bwa-mem2 with the kernels running SSE2 vector instructions (128 bits vector-width)  
  3. Vector AVX2: It runs bwa-mem2 with the kernels running AVX2 vector instructions (256 bits vector-width)  
  4. Vector AVX512: It runs bwa-mem2 with the kernels running AVX512 vector instructions (512 bits vector-width)  

- BWA-MEM2 displays run-time performance profiling of the code on standard output. (Hopefully, we managed it to be self-explanatory)  

----------------------------------------------------------------------------------
