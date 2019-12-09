/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Vasimuddin Md, Sanchit Misra, Intel Corporation, Heng Li.

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

Authors: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>.
*****************************************************************************************/

//#include "sampling.h"
#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
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

