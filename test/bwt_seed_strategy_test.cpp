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

#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include "FMI_search.h"
#include <omp.h>
#include <string.h>

#ifdef VTUNE_ANALYSIS
#include <ittnotify.h>
#endif

#define QUERY_DB_SIZE 512000000

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

    int32_t numReads = 0;
    int64_t total_size = 0;
    gzFile fp = gzopen(argv[2], "r");
	if (fp == 0)
	{
		fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[2]);
        exit(EXIT_FAILURE);
	}
    
    printf("before reading sequences\n");
    bseq1_t *seqs = bseq_read_one_fasta_file(QUERY_DB_SIZE, &numReads, fp, &total_size);

    if(seqs == NULL)
    {
        printf("ERROR! seqs = NULL\n");
        exit(0);
    }
    int32_t *query_cum_len_ar = (int32_t *)_mm_malloc(numReads * sizeof(int32_t), 64);

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
        query_cum_len_ar[st] = st * readlength;
        cind=st*readlength;
        for(r = 0; r < readlength; ++r) {
            switch(seqs[st].seq[r])
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
            seqs,
            query_cum_len_ar,
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
    
    free(enc_qdb);
    _mm_free(matchArray);
    _mm_free(max_intv_array);
    delete fmiSearch;
    return 0;
}

