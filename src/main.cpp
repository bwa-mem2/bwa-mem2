/*************************************************************************************
                    GNU GENERAL PUBLIC LICENSE
           		      Version 3, 29 June 2007

BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
Copyright (C) 2019  Vasimuddin Md, Sanchit Misra, Intel Corporation, Heng Li.
    
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License at https://www.gnu.org/licenses/ for more details.


TERMS AND CONDITIONS FOR DISTRIBUTION OF THE CODE
                                             
1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 
3. Neither the name of Intel Corporation nor the names of its contributors may
   be used to endorse or promote products derived from this software without
   specific prior written permission.

Authors: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>.
*****************************************************************************************/

// ----------------------------------
#include "main.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "2.0pre1"
#endif


// ----------------------------------
uint64_t proc_freq, tprof[LIM_R][LIM_C], prof[LIM_R];
int nthreads;
int num_ranks = 1, myrank = 0;
int64_t reference_seq_len;
// ----------------------------------


int usage()
{
	fprintf(stderr, "Usage: bwa-mem2 <command> <arguments>\n");
	fprintf(stderr, "Commands:\n");
	fprintf(stderr, "  index         create index\n");
	fprintf(stderr, "  mem           alignment\n");
	fprintf(stderr, "  version       print version number\n");
	return 1;
}

int main(int argc, char* argv[])
{
		
	// ---------------------------------	
	uint64_t tim = __rdtsc();
	sleep(1);
	proc_freq = __rdtsc() - tim;

	int ret = -1;
	if (argc < 2) return usage();

	if (strcmp(argv[1], "index") == 0)
	{
		 uint64_t tim = __rdtsc();
		 ret = bwa_index(argc-1, argv+1);
		 tprof[INDEX][0] += __rdtsc() - tim;
		 return ret;
	}
	else if (strcmp(argv[1], "mem") == 0)
	{
		uint64_t tim = __rdtsc();
		kstring_t pg = {0,0,0};
		extern char *bwa_pg;

		fprintf(stderr, "-----------------------------\n");
#if __AVX512BW__
		fprintf(stderr, "Executing in AVX512 mode!!\n");
#endif
#if ((!__AVX512BW__) & (__AVX2__))
		fprintf(stderr, "Executing in AVX2 mode!!\n");
#endif
#if ((!__AVX512BW__) && (!__AVX2__) && (__SSE2__))
		fprintf(stderr, "Executing in SSE4.1 mode!!\n");
#endif
#if ((!__AVX512BW__) && (!__AVX2__) && (!__SSE2__))
		fprintf(stderr, "Executing in Scalar mode!!\n");
#endif
		fprintf(stderr, "-----------------------------\n");

		ksprintf(&pg, "@PG\tID:bwa\tPN:bwa\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);
		for (int i = 1; i < argc; ++i) ksprintf(&pg, " %s", argv[i]);
		ksprintf(&pg, "\n");
		bwa_pg = pg.s;
		ret = main_mem(argc-1, argv+1);
		free(bwa_pg);
		tprof[MEM][0] = __rdtsc() - tim;
		// return ret;
	}
	else if (strcmp(argv[1], "version") == 0)
	{
		puts(PACKAGE_VERSION);
		return 0;
	}
	
	/* Display runtime profiling stats */
	display_stats();
	
	fprintf(stderr, "\nImportant parameter settings: \n");
	fprintf(stderr, "\tBATCH_SIZE: %d\n", BATCH_SIZE);
	fprintf(stderr, "\tMAX_SEQ_LEN_REF: %d\n", MAX_SEQ_LEN_REF);
	fprintf(stderr, "\tMAX_SEQ_LEN_QER: %d\n", MAX_SEQ_LEN_QER);
	fprintf(stderr, "\tMAX_SEQ_LEN8: %d\n", MAX_SEQ_LEN8);
	fprintf(stderr, "\tSEEDS_PER_READ: %d\n", SEEDS_PER_READ);
	fprintf(stderr, "\tSIMD_WIDTH8 X: %d\n", SIMD_WIDTH8);
	fprintf(stderr, "\tSIMD_WIDTH16 X: %d\n", SIMD_WIDTH16);
	fprintf(stderr, "\tAVG_SEEDS_PER_READ: %d\n", AVG_SEEDS_PER_READ);

	return 0;
}
