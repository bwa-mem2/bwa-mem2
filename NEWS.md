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
