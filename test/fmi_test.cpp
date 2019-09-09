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
        printf("Need seven arguments : ref_file query_set batch_size readlength minSeedLen n_threads\n");
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
    int numthreads=atoi(argv[6]);
    assert(numthreads > 0);
    assert(numthreads <= omp_get_max_threads());

    int64_t num_batches = (numReads + batch_size - 1 ) / batch_size;
    int64_t numTotalSmem[num_batches];
#pragma omp parallel num_threads(numthreads)
    {
        int tid = omp_get_thread_num();

        if(tid == 0)
            printf("Running %d threads\n", omp_get_num_threads());
    }

    int32_t *min_intv_array = (int32_t *)_mm_malloc(numReads * sizeof(int32_t), 64);
    int64_t i;
    for(i = 0; i < numReads; i++)
    {
        min_intv_array[i] = 1;
    }

    int64_t startTick, endTick;
#ifdef VTUNE_ANALYSIS
    __itt_resume();
#endif
    startTick = __rdtsc();
#if 0
    fmiSearch->getSMEMs(enc_qdb,
            numReads,
            batch_size,
            readlength,
            minIntv,
            numthreads,
            matchArray,
            numTotalSmem);
#else
    memset(numTotalSmem, 0, num_batches * sizeof(int64_t));
    int64_t workTicks[numthreads];
    memset(workTicks, 0, numthreads * sizeof(int64_t));

#pragma omp parallel num_threads(numthreads)
    {
        int32_t *rid_array = (int32_t *)_mm_malloc(batch_size * sizeof(int32_t), 64);
        int32_t tid = omp_get_thread_num();
        int64_t startTick = __rdtsc();

#pragma omp for schedule(dynamic)
        for(i = 0; i < numReads; i += batch_size)
        {
            int64_t st1 = __rdtsc();
            int32_t batch_count = batch_size;
            if((i + batch_count) > numReads) batch_count = numReads - i;
            int32_t j;
            for(j = 0; j < batch_count; j++)
            {
                rid_array[j] = j;
            }
            int32_t batch_id = i/batch_size;
            //printf("%d] i = %d, batch_count = %d, batch_size = %d\n", tid, i, batch_count, batch_size);
            //fflush(stdout);
            fmiSearch->getSMEMsAllPosOneThread(enc_qdb + i * readlength,
                    min_intv_array + i,
                    rid_array,
                    batch_count,
                    batch_size,
                    readlength,
                    minSeedLen,
                    matchArray + i * readlength,
                    numTotalSmem + batch_id);
            //printf("numTotalSmem = %d\n", numTotalSmem[0]);
            //fflush(stdout);
            fmiSearch->sortSMEMs(matchArray + i * readlength,
                    numTotalSmem + batch_id,
                    batch_count,
                    readlength,
                    1);
            for(j = 0; j < numTotalSmem[batch_id]; j++)
            {
                matchArray[i * readlength + j].rid += i;
            }
            int64_t et1 = __rdtsc();
            workTicks[tid] += (et1 - st1);
        }

        int64_t endTick = __rdtsc();
        printf("%d] %ld ticks, workTicks = %ld\n", tid, endTick - startTick, workTicks[tid]);
        _mm_free(rid_array);
    }

#endif
    endTick = __rdtsc();
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
    int64_t sumTicks = 0;
    int64_t maxTicks = 0;
    for(i = 0; i < numthreads; i++)
    {
        sumTicks += workTicks[i];
        if(workTicks[i] > maxTicks) maxTicks = workTicks[i];
    }
    double avgTicks = (sumTicks * 1.0) / numthreads;
    printf("avgTicks = %lf, maxTicks = %ld, load imbalance = %lf\n", avgTicks, maxTicks, maxTicks/avgTicks);

    printf("Consumed: %ld cycles\n", endTick - startTick);
    //printf("beCalls = %lld\n", fmiSearch->beCalls);

    int64_t totalSmem = 0;
    int32_t batch_id = 0;
    for(batch_id = 0; batch_id < num_batches; batch_id++)
    {
        totalSmem += numTotalSmem[batch_id];
    }
    printf("totalSmems = %ld\n", totalSmem);

    int32_t prevRid = -1;
    for(batch_id = 0; batch_id < num_batches; batch_id++)
    {
        int32_t first = batch_id * batch_size;;
        SMEM *myMatchArray = matchArray + first * readlength;
        int64_t i;
        for(i = 0; i < numTotalSmem[batch_id]; i++)
        {
            SMEM smem = myMatchArray[i];
            if(smem.rid != prevRid)
            {
                int32_t j;
                for(j = prevRid + 1; j <= smem.rid; j++)
                    printf("%u:\n", j);
            }
            prevRid = smem.rid;
            printf("[%u,%u]", smem.m, smem.n + 1);
            // printf("%u, %u]", smem.k, smem.s);
#if 1
            printf(" ["); 
            int64_t u1, u2, u3;
            u1 = smem.k;
            u2 = smem.k + smem.s;
            for(u3 = u1; u3 < u2; u3++)
            {
                printf("%ld,", fmiSearch->get_sa_entry(u3));
            }
            printf("]"); 
#endif
            printf("\n");
        }
    }

    free(query_seq);
    free(enc_qdb);
    _mm_free(matchArray);
    _mm_free(min_intv_array);
    delete fmiSearch;
    return 0;
}

