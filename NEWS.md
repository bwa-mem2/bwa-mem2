Release 2.2.1 (17 March 2021)
---------------------------------
Hotfix for v2.2: Fixed the bug mentioned in #135.


Release 2.2 (8 March 2021)
---------------------------------
Changes since the last release (2.1):

* Passed the validation test on ~88 billions reads (Credits: Keiran Raine, CASM division, Sanger Institute)
* Fixed bugs reported in #109 causing mismatch between bwa-mem and bwa-mem2
* Fixed the issue (# 112) causing crash due to corrupted thread id 
* Using all the SSE flags to create optimized SSE41 and SSE42 binaries


Release 2.1 (16 October 2020)
---------------------------------
Release 2.1 of BWA-MEM2.

Changes since the last release (2.0):
* *Smaller index*: the index size on disk is down by 8 times and in memory by 4 times due to moving to only one type of FM-index (2bit.64 instead of 2bit.64 and 8bit.32) and 8x compression of suffix array. For example, for human genome, index size on disk is down to ~10GB from ~80GB and memory footprint is down to ~10GB from ~40GB. There is a substantial decrease in index IO time due to the reduction and hardly any performance impact on read mapping.

* Added support for 2 more execution modes: sse4.2 and avx.

* Fixed multiple bugs including those reported in Issues #71, #80 and #85.

* Merged multiple pull requests.


Release 2.0 (9 July 2020)
---------------------------------
This is the first production release of BWA-MEM2.

Changes since the last release:
* Made the source code more secure with more than 300 changes all across it.

* Added support for memory re-allocations in case the pre-allocated fixed memory is insufficient.

* Added support for MC flag in the sam file and support for -5, -q flags in the command line.

* The output is now identical to the output of bwa-mem-0.7.17.

* Merged index building code with FMI_Search class.

* Added support for different ways to input read files, now, it is same as bwa-mem.

* Fixed a bug in AVX512 sam processing part, which was leading to incorrect output.


Release 2.0pre2 (4 February 2020)
---------------------------------

Miscellaneous changes:

 * Changed the license from GPL to MIT.

 * IMPORTANT: the index structure has changed since commit 6743183. Please
   rebuild the index if you are using a later commit or the new release.

 * Added charts in README.md comparing the performance of bwa-mem2 with bwa-mem.

Major code changes:

 * Fixed working for variable length reads.

 * Fixed a bug involving reads of length greater than 250bp.

 * Added support for allocation of more memory in small chunks if large
   pre-allocated fixed memory is insufficient. This is needed very rarely
   (thus, having no impact on performance) but prevents asserts from failing
   (code from crashing) in that scenario. 

 * Fixed a memory leak due to not releasing the memory allocated for seeds
   after smem.

 * Fixed a segfault due to non-alignment of small allocated memory in the
   optimized banded Smith-Waterman.

 * Enabled working with genomes larger than 7-8 billion nucleotides (e.g. Wheat
   genome).

 * Fixed a segfault occuring (with gcc compiler) while reading the index.
