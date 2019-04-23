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

//#include "sampling.h"
#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include "fasta_file.h"
#include "FMI_search.h"
#include <omp.h>
#include <string.h>

#ifdef VTUNE_ANALYSIS
#include <ittnotify.h>
#endif

#define MAX_NUM_SMEM 20000000
#define MAX_NUM_OFFSET 400000000

int myrank, num_ranks;

int64_t loadData(FILE *infp, SMEM *smemArray) {

	int64_t numSmems = 0;
    while(numSmems < MAX_NUM_SMEM) {
		int32_t val = 0;
        SMEM smem;
        smem.rid = 0;
        smem.m = 0;
        smem.n = 0;
        smem.l = 0;
        smem.k = 0;
        smem.s = 0;
		if ((val = fscanf(infp, "%ld, %ld", &(smem.k), &(smem.s))) == -1)
			break;
		smemArray[numSmems++] = smem;
	}
	
	return numSmems;
}


int main(int argc, char **argv) {
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
    printf("argc = %d\n", argc);
    if(argc!=5)
    {
        printf("Need five arguments : ref_file sa_coord_file out_file max_occ\n");
        return 1;
    }

	// Open file
	FILE *infp = NULL;
	infp = fopen(argv[2], "r");
	if (infp == NULL) {
		printf("Unable to open sa2ref input file!!"); exit(0);
	}
	FILE *outfp = NULL;
	outfp = fopen(argv[3], "w");
	if (outfp == NULL) {
		printf("Unable to open sa2ref output file!!"); exit(0);
	}
	
    FMI_search *fmiSearch = new FMI_search(argv[1]);

	SMEM *smemArray = (SMEM *) _mm_malloc (MAX_NUM_SMEM *
												   sizeof(SMEM), 64);
	int64_t numSmems = loadData(infp, smemArray);
	int64_t *refCoordArray = (int64_t *) _mm_malloc (MAX_NUM_OFFSET *
												   sizeof(int64_t), 64);
    
    int max_occ = atoi(argv[4]);

    int64_t startTick, endTick;
#ifdef VTUNE_ANALYSIS
    __itt_resume();
#endif
    startTick = __rdtsc();

    int64_t i;
    int64_t numOffs = 0;
    for(i = 0; i < numSmems; i++)
    {
        int32_t count = 0;
        fmiSearch->get_sa_entries(smemArray + i, refCoordArray + numOffs, &count, 1, max_occ);
        numOffs += count;
    }
    endTick = __rdtsc();
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
    printf("Consumed: %ld cycles\n", endTick - startTick);

    for(i = 0; i < numOffs; i++)
    {
        int64_t coord = refCoordArray[i];
        fprintf(outfp, "%ld\n", coord);
    }

    fclose(infp);
    fclose(outfp);
    _mm_free(smemArray);
    _mm_free(refCoordArray);
    delete fmiSearch;
    return 0;
}

