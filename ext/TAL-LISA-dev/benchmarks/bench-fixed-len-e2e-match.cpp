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

Authors: Sanchit Misra <sanchit.misra@intel.com>; Vasimuddin Md <vasimuddin.md@intel.com>; Kanak Mahadik
*****************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <omp.h>
#include <string.h>
#include <immintrin.h>

#include "FMI_search.h"

#ifdef VTUNE_ANALYSIS
#include <ittnotify.h>
#endif

#ifdef IACA_ANALYSIS
#include "iacaMarks.h"
#endif


#define QUERY_DB_SIZE 12500000000L

#define DUMMY_CHAR 6

#ifdef __AVX512PF__
#include <hbwmalloc.h>

#define TID_SECOND_NUMA 999999

#else

#ifdef __AVX512BW__
#define TID_SECOND_NUMA 56
#else
#define TID_SECOND_NUMA 18
#endif

#endif

#define BEST_BATCH_SIZE 24
#define MAX_BATCH_SIZE 128

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


int compare_uint(const void *a, const void *b)
{
    uint32_t *pa = (uint32_t *)a;
    uint32_t *pb = (uint32_t *)b;

    uint32_t va = *pa;
    uint32_t vb = *pb;

    if(va < vb) return -1;
    if(va > vb) return 1;
    return 0;
}


int64_t loopTicks;
int64_t loopCount;

int64_t get_freq()
{
    int64_t startTick = __rdtsc();
    sleep(1);
    int64_t endTick = __rdtsc();
    return (endTick - startTick);
}

void exact_match(FMI_search *fmiSearch,
                 uint8_t *enc_qdb,
                 int64_t numQueries, uint32_t query_length,
                 int64_t *k_l_range,
                 int32_t batch_size, int32_t numthreads)
{
#pragma omp parallel num_threads(numthreads)
    {
        int32_t tid = omp_get_thread_num();
        //int64_t start_thread = __rdtsc();
        uint32_t perThreadQuota = (numQueries + numthreads - 1) / numthreads;
        uint32_t first = tid * perThreadQuota;
        uint32_t last  = (tid + 1) * perThreadQuota;
        if(last > numQueries) last = numQueries;

        fmiSearch->exact_search(enc_qdb + first * query_length, last - first, query_length, k_l_range + 2 * first, batch_size);
        //int64_t end_thread = __rdtsc();
        //printf("%d] %ld ticks\n", tid, end_thread - start_thread);
    }

}

void inexact_match(FMI_search *fmiSearch,
                 uint8_t *enc_qdb,
                 int32_t z,
                 int64_t numQueries, uint32_t query_length,
                 int64_t *k_l_range, int32_t batch_size,
                 int32_t numthreads)
{

#pragma omp parallel num_threads(numthreads)
    {
        int32_t tid = omp_get_thread_num();
        uint32_t perThreadQuota = (numQueries + numthreads - 1) / numthreads;
        uint32_t first = tid * perThreadQuota;
        uint32_t last  = (tid + 1) * perThreadQuota;
        if(last > numQueries) last = numQueries;

        fmiSearch->inexact_search(enc_qdb + first * query_length, z,
                                  last - first, query_length, k_l_range + first * (query_length * 3 * z + 1), batch_size);
    }
}

int main(int argc, char **argv) {
    loopTicks = 0;
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
    if(argc < 5)
    {
        printf("Need following arguments : reference_file query_set z n_threads [batch_size]\n");
        return 1;
    }
    int32_t numQueries = 0;
    int64_t total_size = 0;
    gzFile fp = gzopen(argv[2], "r");
	if (fp == 0)
	{
		fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[2]);
        exit(EXIT_FAILURE);
	}
    
    FMI_search *fmiSearch = new FMI_search(argv[1]);
    fmiSearch->load_index_forward_only();

    printf("before reading queries\n");
    bseq1_t *seqs = bseq_read_one_fasta_file(QUERY_DB_SIZE, &numQueries, fp, &total_size);

    if(seqs == NULL)
    {
        printf("ERROR! seqs = NULL\n");
        exit(EXIT_FAILURE);
    }

    int max_query_length = seqs[0].l_seq;
    int min_query_length = seqs[0].l_seq;
    int max_len_id = 0;
    int min_len_id = 0;
    for(int i = 1; i < numQueries; i++)
    {
        if(max_query_length < seqs[i].l_seq)
        {
            max_query_length = seqs[i].l_seq;
            max_len_id = i;
        }
        if(min_query_length > seqs[i].l_seq)
        {
            min_query_length = seqs[i].l_seq;
            min_len_id = i;
        }
    }
    printf("max_query_length = %d, min_query_length = %d\n", max_query_length, min_query_length);
    printf("max_len_id = %d, min_len_id = %d\n", max_len_id, min_len_id);
    assert(max_query_length > 0);
    assert(max_query_length < 10000);
    assert(max_query_length == min_query_length);
    assert(numQueries > 0);
    assert(numQueries * (int64_t)max_query_length < QUERY_DB_SIZE);
    int32_t query_length = max_query_length;
    printf("numQueries = %d, query_length = %d\n", numQueries, query_length);
    uint8_t *enc_qdb = (uint8_t *)_mm_malloc(numQueries * query_length * sizeof(uint8_t), 64);
    assert(enc_qdb != NULL);

    int64_t cind,st;
#if 0
    printf("Priting query\n");
    for(st = 0; st < query_length; st++)
    {
        printf("%c", seqs[0].seq[st]);
    }
    printf("\n");
#endif
    uint64_t r;
    int32_t orig_query_count = numQueries;
    numQueries = 0;
    for (st=0; st < orig_query_count; st++) {
        cind = numQueries * query_length;
        int nflag = 0;
        for(r = 0; r < query_length; ++r) {
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
                default: nflag = 1;
            }
        }
        if(nflag == 0)
            numQueries++;
    }
    int32_t z = atoi(argv[3]);
    int numthreads = atoi(argv[4]);
    assert(numthreads > 0);
    assert(numthreads <= omp_get_max_threads());

    if((z != 0) && (z != 1))
    {
        printf("ERROR! z can only be zero or one.\n");
        exit(-1);
    }

    int64_t numSeq = numQueries * (query_length * 3 * z + 1);
    printf("numSeq = %d, k_l size = %ld\n", numSeq, numQueries * (query_length * 3 * z + 1) * 2 * sizeof(int64_t));
    int64_t *k_l_range = (int64_t *)_mm_malloc(numSeq * 2 * sizeof(int64_t), 64);
    memset(k_l_range, 0, numSeq * 2 * sizeof(int64_t));
    int64_t startTick, endTick;

    int32_t batch_size = BEST_BATCH_SIZE;
    if(argc > 5)
        batch_size = atoi(argv[5]);
    assert(batch_size < MAX_BATCH_SIZE);
    if(numQueries < (batch_size * numthreads))
    {
        printf("very few reads for %d threads\n", numthreads);
        printf("numQueries = %d, batch_size = %d, numthreads = %d\n", numQueries, batch_size, numthreads);
        exit(0);
    }

#ifdef VTUNE_ANALYSIS
    __itt_resume();
#endif
    startTick = __rdtsc();

    if(z > 0)
    {
        inexact_match(fmiSearch,
                enc_qdb, z,
                numQueries, query_length,
                k_l_range, batch_size,
                numthreads);
    }
    else
    {
        exact_match(fmiSearch,
                 enc_qdb,
                 numQueries, query_length,
                 k_l_range,
                 batch_size, numthreads);
    }
    endTick = __rdtsc();
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
    int64_t totalTicks = endTick - startTick;
   
    int64_t numMatches = 0;
    int64_t numMapped  = 0;
    //uint32_t buf[1048576];

    _mm_free(enc_qdb);    

    //uint32_t *sa_bwt = NULL;
    //sa_bwt = (uint32_t *)_mm_malloc(ref_seq_len * sizeof(uint32_t), 64);
    //fread(sa_bwt, sizeof(uint32_t), ref_seq_len, cpstream);

    int64_t i = 0;
    for(i = 0; i < numQueries; i++)
    {
        int32_t j;
        int32_t max_map_per_read = z * query_length * 3 + 1;
        int64_t read_num_matches = 0;
        for(j = 0; j < max_map_per_read; j++)
        {
            int64_t u1, u2;
            u1 = k_l_range[(i * max_map_per_read + j) * 2];
            u2 = k_l_range[(i * max_map_per_read + j) * 2 + 1];
            printf("%ld: ", i);
            if(u2 - u1 > 0)
                printf("%ld, %ld", u1, u2);
            printf("\n");
#if 0
#ifdef PRINT_OUTPUT
            memcpy(buf + read_num_matches, sa_bwt + u1, (u2 - u1) * sizeof(uint32_t));
#endif
#endif
            read_num_matches += (u2 - u1);
        }
#if 0
#ifdef PRINT_OUTPUT
            qsort(buf, read_num_matches, sizeof(uint32_t), compare_uint);
            for(j = 0; j < read_num_matches; j++)
            {
                printf("%u ", buf[j]);
            }
            printf("\n");
#endif
#endif
        if(read_num_matches > 0) numMapped++;
        numMatches += read_num_matches;
    }
    printf("numQueries = %d\n",   numQueries);
    printf("numMatches = %ld\n", numMatches);
    printf("numMapped = %ld\n",  numMapped);
    printf("loopTicks = %ld\n", loopTicks);

    printf("Consumed %ld cycles\n", totalTicks);
    printf("Consumed %0.3lf secs\n", totalTicks * 1.0 / get_freq());
    printf("log, %d, %u, %d, %ld, %0.3lf\n", query_length, orig_query_count, numthreads, totalTicks, totalTicks * 1.0 / get_freq());

    //fclose(cpstream);

    //_mm_free(sa_bwt);
    _mm_free(k_l_range);    
    delete fmiSearch;
    bseq_destroy(seqs);
    return 0;
}

