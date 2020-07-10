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
