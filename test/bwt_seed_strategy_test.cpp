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

#define QUERY_DB_SIZE 12500000000L

int myrank, num_ranks;

int main(int argc, char **argv) {
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif

    if(argc!=7)
    {
        printf("Need seven arguments : ref_file query_set batch_size readlength minSeedLen max_mem_intv\n");
        return 1;
    }

    char *query_seq=(char *)malloc(QUERY_DB_SIZE*sizeof(char));
    int64_t numReads;
    numReads=read_multi_fasta_file(argv[2],query_seq);

    if(numReads==-1)
    {
        printf("Error opening query'%s'. Bailing out.",argv[2]);
        free(query_seq);
        return -1;
    }

    FMI_search *fmiSearch = new FMI_search(argv[1]);

    int readlength=atoi(argv[4]);
    assert(readlength > 0);
    assert(readlength < 10000);
    assert(numReads > 0);
    assert(numReads * readlength < QUERY_DB_SIZE);

    uint8_t *enc_qdb=(uint8_t *)malloc(numReads*readlength*sizeof(uint8_t));

    int64_t cind,st;
    uint64_t r;
    for (st=0; st < numReads; st++) {
        cind=st*readlength;
        for(r = 0; r < readlength; ++r) {
            switch(query_seq[r+cind])
            {
                case '0': enc_qdb[r+cind]=0;
                          break;
                case '1': enc_qdb[r+cind]=1;
                          break;
                case '2': enc_qdb[r+cind]=2;
                          break;
                case '3': enc_qdb[r+cind]=3;
                          break;
                default: enc_qdb[r+cind]=4;
            }
            //printf("%c %d\n", query_seq[r+cind], enc_qdb[r + cind]);
            //printf("%d", enc_qdb[r + cind]);
        }
    }

    int batch_size=0;
    batch_size=atoi(argv[3]);
    assert(batch_size > 0);
    assert(batch_size <= numReads);

    SMEM *matchArray = (SMEM *)_mm_malloc(numReads * readlength * sizeof(SMEM), 64);

    int32_t minSeedLen = atoi(argv[5]);
    int32_t max_mem_intv = atoi(argv[6]);
    printf("minSeedLen = %d, max_mem_intv = %d\n", minSeedLen, max_mem_intv);
    int numthreads=1;
    int64_t numTotalSmem[numthreads];
#pragma omp parallel num_threads(numthreads)
    {
        int tid = omp_get_thread_num();

        if(tid == 0)
            printf("Running %d threads\n", omp_get_num_threads());
    }

    int32_t *max_intv_array = (int32_t *)_mm_malloc(numReads * sizeof(int32_t), 64);
    int32_t i;
    for(i = 0; i < numReads; i++)
    {
        max_intv_array[i] = max_mem_intv;
    }

    int64_t startTick, endTick;
#ifdef VTUNE_ANALYSIS
    __itt_resume();
#endif
    startTick = __rdtsc();
    numthreads = 1;
    numTotalSmem[0] = fmiSearch->bwtSeedStrategyAllPosOneThread(enc_qdb,
            max_intv_array,
            numReads,
            readlength,
            minSeedLen,
            matchArray);
    endTick = __rdtsc();
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
    printf("Consumed: %ld cycles\n", endTick - startTick);

    int64_t totalSmem = 0;
    int tid;
    for(tid = 0; tid < numthreads; tid++)
    {
        totalSmem += numTotalSmem[tid];
    }
    printf("totalSmems = %ld\n", totalSmem);

    fmiSearch->sortSMEMs(matchArray,
            numTotalSmem,
            numReads,
            readlength,
            numthreads);

    int32_t perThreadQuota = (numReads + (numthreads - 1)) / numthreads;
    for(tid = 0; tid < numthreads; tid++)
    {
        int32_t first = tid * perThreadQuota;
        SMEM *myMatchArray = matchArray + first * readlength;
        int64_t i = 0;
        int64_t rid;
        for(rid = 0; rid < numReads; rid++)
        {
            while((i < numTotalSmem[tid]) && (myMatchArray[i].rid == rid))
            {
                SMEM smem = myMatchArray[i];
                //printf("\n%u: ", smem.rid);
                printf("[ %ld %ld %ld %u %u ] ", smem.k, smem.l, smem.s, smem.m, smem.n + 1);
                i++;
            }
            printf("\n");
        }
    }
    
    free(query_seq);
    free(enc_qdb);
    _mm_free(matchArray);
    _mm_free(max_intv_array);
    delete fmiSearch;
    return 0;
}

