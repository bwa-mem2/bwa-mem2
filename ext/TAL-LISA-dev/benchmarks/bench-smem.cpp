/*************************************************************************************
MIT License

Copyright (c) 2020 Intel Labs

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Authors: Sanchit Misra <sanchit.misra@intel.com>; Vasimuddin Md <vasimuddin.md@intel.com>
*****************************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include <omp.h>
#include <string.h>
// #include "bwa.h"
#include "FMI_search.h"

#ifdef VTUNE_ANALYSIS
#include <ittnotify.h>
#endif

#define QUERY_DB_SIZE 1280000000
int myrank, num_ranks;

void bseq_destroy(bseq1_t *s)
{
    if(s)
    {
        if(s->name) free(s->name);
        if(s->comment) free(s->comment);
        if(s->seq) free(s->seq);
        if(s->qual) free(s->qual);
        if(s->sam) free(s->sam);
        free(s);
    }
}

int main(int argc, char **argv) {
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
    {
        printf("Running:\n");
        int i;
        for(i = 0; i < argc; i++)
        {
            printf("%s ", argv[i]);
        }
        printf("\n");
    }
    if(argc!=6)
    {
        printf("Need five arguments : ref_file query_set batch_size minSeedLen n_threads\n");
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
    
    FMI_search *fmiSearch = new FMI_search(argv[1]);
    fmiSearch->load_index_with_rev_complement();

    printf("before reading sequences\n");
    bseq1_t *seqs = bseq_read_one_fasta_file(QUERY_DB_SIZE, &numReads, fp, &total_size);

    if(seqs == NULL)
    {
        printf("ERROR! seqs = NULL\n");
        exit(EXIT_FAILURE);
    }
    int32_t *query_cum_len_ar = (int32_t *)_mm_malloc(numReads * sizeof(int32_t), 64);

    int max_readlength = seqs[0].l_seq;
    int min_readlength = seqs[0].l_seq;
    for(int i = 1; i < numReads; i++)
    {
        if(max_readlength < seqs[i].l_seq)
            max_readlength = seqs[i].l_seq;
        if(min_readlength > seqs[i].l_seq)
            min_readlength = seqs[i].l_seq;
    }
    assert(max_readlength > 0);
    assert(max_readlength < 10000);
    assert(numReads > 0);
    assert(numReads * max_readlength < QUERY_DB_SIZE);
    printf("numReads = %d, max_readlength = %d, min_readlength = %d\n", numReads, max_readlength, min_readlength);
    uint8_t *enc_qdb=(uint8_t *)malloc(numReads * max_readlength * sizeof(uint8_t));
    assert(enc_qdb != NULL);

    int64_t cind,st;
#if 0
    printf("Priting query\n");
    for(st = 0; st < max_readlength; st++)
    {
        printf("%c", seqs[0].seq[st]);
    }
    printf("\n");
#endif
    uint64_t r;
    for (st=0; st < numReads; st++) {
        query_cum_len_ar[st] = st * max_readlength;
        cind=st*max_readlength;
        for(r = 0; r < max_readlength; ++r) {
            switch(seqs[st].seq[r])
            {
                case 'A': enc_qdb[r+cind]=0;
                          break;
                case 'C': enc_qdb[r+cind]=1;
                          break;
                case 'G': enc_qdb[r+cind]=2;
                          break;
                case 'T': enc_qdb[r+cind]=3;
                          break;
                default: enc_qdb[r+cind]=4;
            }
            //printf("%c %d\n", seqs[st].seq[r], enc_qdb[r + cind]);
        }
    }

    int batch_size=0;
    batch_size=atoi(argv[3]);
    assert(batch_size > 0);
    assert(batch_size <= numReads);


    int32_t minSeedLen = atoi(argv[4]);
    int numthreads=atoi(argv[5]);
    assert(numthreads > 0);
    assert(numthreads <= omp_get_max_threads());
    SMEM *matchArray[numthreads];

    int64_t num_batches = (numReads + batch_size - 1 ) / batch_size;
    int64_t *numTotalSmem = (int64_t *)_mm_malloc(num_batches * sizeof(int64_t), 64);;
    SMEM **batchStart = (SMEM **)_mm_malloc(num_batches * sizeof(SMEM *), 64);;
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
    memset(numTotalSmem, 0, num_batches * sizeof(int64_t));
    memset(batchStart, 0, num_batches * sizeof(int64_t));
    int64_t workTicks[numthreads];
    memset(workTicks, 0, numthreads * sizeof(int64_t));
    int64_t perThreadQuota = numReads / numthreads;

#pragma omp parallel num_threads(numthreads)
    {
        int32_t *rid_array = (int32_t *)_mm_malloc(batch_size * sizeof(int32_t), 64);
        int32_t tid = omp_get_thread_num();
        int64_t matchArrayAlloc = perThreadQuota * 20;
        matchArray[tid] = (SMEM *)malloc(matchArrayAlloc * sizeof(SMEM));
        int64_t myTotalSmems = 0;
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
            if((matchArrayAlloc - myTotalSmems) < (batch_size * max_readlength))
            {
                matchArrayAlloc *= 2;
                matchArray[tid] = (SMEM *)realloc(matchArray[tid], matchArrayAlloc * sizeof(SMEM)); 
            }
            fmiSearch->getSMEMsAllPosOneThread(enc_qdb + i * max_readlength,
                    min_intv_array + i,
                    rid_array,
                    batch_count,
                    batch_size,
                    seqs + i,
                    query_cum_len_ar,
                    max_readlength,
                    minSeedLen,
                    matchArray[tid] + myTotalSmems,
                    numTotalSmem + batch_id);
            batchStart[batch_id] = matchArray[tid] + myTotalSmems;
            fmiSearch->sortSMEMs(matchArray[tid] + myTotalSmems,
                    numTotalSmem + batch_id,
                    batch_count,
                    max_readlength,
                    1);
            for(j = 0; j < numTotalSmem[batch_id]; j++)
            {
                matchArray[tid][myTotalSmems + j].rid += i;
            }
            myTotalSmems += numTotalSmem[batch_id];
            int64_t et1 = __rdtsc();
            workTicks[tid] += (et1 - st1);
        }

        int64_t endTick = __rdtsc();
        printf("%d] %ld ticks, workTicks = %ld\n", tid, endTick - startTick, workTicks[tid]);
        _mm_free(rid_array);
    }

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

    int64_t totalSmem = 0;
    int32_t batch_id = 0;
    for(batch_id = 0; batch_id < num_batches; batch_id++)
    {
        totalSmem += numTotalSmem[batch_id];
    }
    printf("totalSmems = %ld\n", totalSmem);

#ifdef PRINT_OUTPUT
    int32_t prevRid = -1;
    for(batch_id = 0; batch_id < num_batches; batch_id++)
    {
        SMEM *myMatchArray = batchStart[batch_id];
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
#endif
    _mm_free(query_cum_len_ar);
    free(enc_qdb);
    for(int tid = 0; tid < numthreads; tid++)
    {
        free(matchArray[tid]);
    }
    _mm_free(min_intv_array);
    _mm_free(numTotalSmem);
    _mm_free(batchStart);
    delete fmiSearch;
    bseq_destroy(seqs);
    return 0;
}

