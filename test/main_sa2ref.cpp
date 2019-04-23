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


#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <omp.h>
#include "sa2ref.hpp"

#ifdef VTUNE_ANALYSIS
#include <ittnotify.h>
#endif

#define MAX_NUM_OFFSET 40000000
#define MAX_NUM_OFFSET_ALLOC 40000000

double clock_freq;
FILE *infp;

int loadData(dataPair *dataPairArray) {

	int32_t numOffs = 0;
    while(numOffs < MAX_NUM_OFFSET_ALLOC) {
		uint64_t offset = 0;
		int val = 0;
		if ((val = fscanf(infp, "%ld", &offset)) == -1)
			break;
		//int val = fscanf(infp, "%d", &offset);
		//fprintf(stderr, "%ld\n", offset);
		//printf("%ld, val: %d\n", offset, val);
		dataPairArray[numOffs++].bwt_offset = offset;
	}
	
	
	return numOffs;
}


int main(int argc, char *argv[])
{
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
	if (argc != 4) {
		printf("usage: <exec> <Ref. genome file> <sa2ref infile> <sa2ref outfile>!!\n");
		exit(0);
	}
	
	clock_freq = _rdtsc();
	sleep(1);
	clock_freq = _rdtsc() - clock_freq;

    int32_t numThreads = 1;
#pragma omp parallel
    {
        int32_t tid = omp_get_thread_num();
        int32_t nt = omp_get_num_threads();
        if(tid == (nt - 1))
        {
            numThreads = nt;
        }
    }

	// Open file
	infp = NULL;
	infp = fopen(argv[2], "r");
	if (infp == NULL) {
		printf("Unable to open sa2ref input file!!"); exit(0);
	}
	FILE *outfp = NULL;
	outfp = fopen(argv[3], "w");
	if (outfp == NULL) {
		printf("Unable to open sa2ref output file!!"); exit(0);
	}
	
	// alloc memory
	dataPair *dataPairArray = (dataPair *) malloc (MAX_NUM_OFFSET *
												   sizeof(dataPair));
	
	mysa2ref *sa2ref = new mysa2ref(argv[1]);
	
    int64_t startTick, endTick, totalTicks = 0, pTotalTicks = 0;
	clock_t clk_tim;
	int32_t numOffs = 0;
	int i=0;

	numOffs = loadData(dataPairArray);

	fprintf(stderr, "Processing BWT offsets...\n");
#ifdef VTUNE_ANALYSIS
    __itt_resume();
#endif
	startTick = __rdtsc();
    for(i = 0; i < numOffs; i++)
    {
        dataPair *d = dataPairArray + i;
        sa2ref->run_sa2ref(d);
    }
	endTick = __rdtsc();
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
	totalTicks += endTick - startTick;

	fprintf(stderr, "Total SA2Refs performed: %d\n", numOffs);
	fprintf(stderr, "Total time taken: %0.2lf\n", totalTicks*1.0/clock_freq);
    for(i = 0; i < numOffs; i++)
    {
        dataPair *d = dataPairArray + i;
        fprintf(outfp, "%ld\n", d->ref_coord);
    }

	fclose(infp);
	fclose(outfp);
	free(dataPairArray);
}
