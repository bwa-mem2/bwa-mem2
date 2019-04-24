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
#define PACKAGE_VERSION "0.1.0"  // alpha release
#endif


// ----------------------------------
uint64_t proc_freq, tprof[LIM_R][LIM_C], prof[LIM_R];
int nthreads;
int num_ranks = 1, myrank = 0;
int64_t reference_seq_len;
// ----------------------------------


int usage()
{
	printf("\nProgram: BWA-MEM2 (Sequence alignment using Burrows-Wheeeler Transform)\n");
	printf("Version: %s\n", PACKAGE_VERSION);
	printf("Contacts: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;\n\t  Heng Li <hli@jimmy.harvard.edu>\n\n");
	printf("This program comes with ABSOLUTELY NO WARRANTY.\n"
		   "This program is free software: you can redistribute it and/or modify\n"
		   "it under the terms of the GNU General Public License as published by\n"
		   "the Free Software Foundation, either version 3 of the License, or\n"
		   "(at your option) any later version.\n\n"
		   "This program is distributed in the hope that it will be useful,\n"
		   "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
		   "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
		   "GNU General Public License for more details.\n\n"
		   "For more information, see <https://www.gnu.org/licenses/>.\n\n");
	printf("This is free software, and you are welcome to redistribute it\n"
		   "under certain conditions.\n"
		   "1. Redistributions of source code must retain the above copyright notice,\n"
		   "this list of conditions and the following disclaimer.\n" 
		   "2. Redistributions in binary form must reproduce the above copyright notice,\n"
		   "this list of conditions and the following disclaimer in the documentation\n"
		   "and/or other materials provided with the distribution. \n"
		   "3. Neither the name of Intel Corporation nor the names of its contributors may\n"
		   "be used to endorse or promote products derived from this software without\n"
		   "specific prior written permission.\n");
	
	printf("\n\tusage: bwa <index | mem> [options]\n\n");
	return 1;
}

int main(int argc, char* argv[])
{
		
	// ---------------------------------	
	uint64_t tim = _rdtsc();
	sleep(1);
	proc_freq = _rdtsc() - tim;

	extern char *bwa_pg;
	kstring_t pg = {0,0,0};
	ksprintf(&pg, "@PG\tID:bwa\tPN:bwa\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);
	for (int i = 1; i < argc; ++i) ksprintf(&pg, " %s", argv[i]); ksprintf(&pg, "\n");
	bwa_pg = pg.s;
	
	int ret = -1;
	if (argc < 2) return usage();

	fprintf(stderr, "\nBWA-MEM2 v0.1.0\n"
		   "This program is under GPLv3. "
		   "For more details, refer to the LICENSE file or run 'bwa-mem2'.\n\n");

	
	fprintf(stderr, "-----------------------------\n");
#if __AVX512BW__
	fprintf(stderr, "Executing in AVX512 mode!!\n");
#endif
#if ((!__AVX512BW__) & (__AVX2__))
	fprintf(stderr, "Executing in AVX2 mode!!\n");
#endif

#if ((!__AVX512BW__) && (!__AVX2__))
	fprintf(stderr, "Executing in Scalar mode!!\n");
#endif
	printf("-----------------------------\n");

	if (strcmp(argv[1], "index") == 0)
	{
		 uint64_t tim = _rdtsc();
		 ret = bwa_index(argc-1, argv+1);
		 fprintf(stderr, "Index creation completed!!!\n");
		 tprof[INDEX][0] += _rdtsc() - tim;

		 fprintf(stderr, "Index creation time: %0.2lf\n\n", tprof[INDEX][0]*1.0/proc_freq);
		 free(bwa_pg);
		 return 1;
	}
	else if (strcmp(argv[1], "mem") == 0)
	{
		uint64_t tim = _rdtsc();
		ret = main_mem(argc-1, argv+1);
		tprof[MEM][0] = _rdtsc() - tim;
		
		if (ret == 1) {
			free(bwa_pg);
			return 1;
		}
	}
	
	/* Display runtime profiling stats */
	if (myrank == 0) {
		display_stats();
	}
	
	if (myrank == 0) {
		fprintf(stderr, "\nImportant parameter settings: \n");
		fprintf(stderr, "\tBATCH_SIZE: %d\n", BATCH_SIZE);
		fprintf(stderr, "\tMAX_SEQ_LEN_REF: %d\n", MAX_SEQ_LEN_REF);
		fprintf(stderr, "\tMAX_SEQ_LEN_QER: %d\n", MAX_SEQ_LEN_QER);
		fprintf(stderr, "\tMAX_SEQ_LEN8: %d\n", MAX_SEQ_LEN8);
		fprintf(stderr, "\tSEEDS_PER_READ: %d\n", SEEDS_PER_READ);
		fprintf(stderr, "\tSIMD_WIDTH8 X: %d\n", SIMD_WIDTH8);
		fprintf(stderr, "\tSIMD_WIDTH16 X: %d\n", SIMD_WIDTH16);
		fprintf(stderr, "\tAVG_SEEDS_PER_READ: %d\n", AVG_SEEDS_PER_READ);
		fprintf(stderr, "\tAVG_AUX_SEEDS_PER_READ: %d\n", AVG_AUX_SEEDS_PER_READ);
	}
	free(bwa_pg);	
	return 1;
}



