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

Authors: Sanchit Misra <sanchit.misra@intel.com>; Vasimuddin Md <vasimuddin.md@intel.com>;
*****************************************************************************************/

#include <stdio.h>
#include "FMI_search.h"

extern int myrank, num_ranks;

#define \
GET_OCC_SCALAR(pp, c, occ_id_pp, y_pp, occ_pp, bwt_str_pp) \
                int64_t occ_id_pp = pp >> CP_SHIFT; \
                int64_t y_pp = pp & CP_MASK; \
                int64_t occ_pp = cp_occ[occ_id_pp].cp_count[c]; \
                uint8_t *bwt_str_pp = cp_occ[occ_id_pp].bwt_str; \
                for(bwt_id = 0; bwt_id < y_pp; bwt_id++) \
                {\
                    occ_pp += (c == bwt_str_pp[bwt_id]); \
                }

#define \
GET_OCC(pp, c, c256, occ_id_pp, y_pp, occ_pp, bwt_str_pp, bwt_pp_vec, mask_pp_vec, mask_pp) \
                int64_t occ_id_pp = pp >> CP_SHIFT; \
                int64_t y_pp = pp & CP_MASK; \
                int64_t occ_pp = cp_occ[occ_id_pp].cp_count[c]; \
                uint8_t *bwt_str_pp = cp_occ[occ_id_pp].bwt_str; \
                __m256i bwt_pp_vec = _mm256_load_si256((const __m256i *)(bwt_str_pp)); \
                __m256i mask_pp_vec = _mm256_cmpeq_epi8(bwt_pp_vec, c256); \
                uint64_t mask_pp = _mm256_movemask_epi8(mask_pp_vec); \
                mask_pp = mask_pp << (32 - y_pp); \
                occ_pp += _mm_countbits_32(mask_pp);


FMI_search::FMI_search(char *ref_file_name)
{
    //beCalls = 0;
    char cp_file_name[1000];
    assert(strnlen(ref_file_name, 1000) + 12 < 1000);
    sprintf(cp_file_name, "%s.bwt.8bit.%d", ref_file_name, CP_BLOCK_SIZE);
    // Read the BWT and FM index of the reference sequence
    FILE *cpstream = NULL;
    cpstream = fopen(cp_file_name,"rb");
    if (cpstream == NULL)
    {
      fprintf(stderr, "ERROR! Unable to open the file: %s\n", cp_file_name);
      exit(0);
    }

    fread(&reference_seq_len, sizeof(int64_t), 1, cpstream);
    assert(reference_seq_len > 0);
    assert(reference_seq_len <= (UINT32_MAX * (int64_t)CP_BLOCK_SIZE));
	if(myrank == 0)
		fprintf(stderr, "reference seq len = %ld\n", reference_seq_len);

    // create checkpointed occ
    uint32_t cp_occ_size = (reference_seq_len >> CP_SHIFT) + 1;
    cp_occ = NULL;

    fread(&count[0], sizeof(int64_t), 5, cpstream);
    cp_occ = (CP_OCC *)_mm_malloc(cp_occ_size * sizeof(CP_OCC), 64);

    fread(cp_occ, sizeof(CP_OCC), cp_occ_size, cpstream);
    int64_t ii = 0;
    for(ii = 0; ii < 5; ii++)// update read count structure
    {
        count[ii] = count[ii] + 1;
    }
    sa = (int64_t *)_mm_malloc(reference_seq_len * sizeof(int64_t), 64);
    fread(sa, sizeof(int64_t), reference_seq_len, cpstream);
    fclose(cpstream);

    sentinel_index = -1;
    int64_t x;
    for(x = 0; x < reference_seq_len; x++)
    {
        if(sa[x] == 0)
            sentinel_index = x;
    }
	if(myrank == 0) {
		fprintf(stderr, "count\n");
		for(x = 0; x < 5; x++)
		{
			fprintf(stderr, "%ld,\t%u\n", x, count[x]);
		}
		fprintf(stderr, "\n");
	}

#if BWA_OTHER_ELE
	if(myrank == 0)
        fprintf(stderr, "Reading other elements of the index from files %s\n",
			   ref_file_name);
        bwa_idx_load_ele(ref_file_name, BWA_IDX_ALL);
#endif
    for(ii = 0; ii < 4; ii++)
    {
        int32_t j;
        for(j = 0; j < 64; j++)
        {
            c_bcast_array[ii * 64 + j] = ii;
        }
    }
    fprintf(stderr, "Done reading Index!!\n");
}

FMI_search::~FMI_search()
{
    _mm_free(sa);
    _mm_free(cp_occ);
}

#if 0
int32_t FMI_search::getSMEMsOnePosOneThread2(uint8_t *enc_qdb,
											 int16_t *query_pos_array,
											 int32_t *min_intv_array,
											 int32_t *rid,
											 int32_t numReads,
											 int32_t numActive,
											 int32_t batch_size,
											 int32_t readlength,
											 int32_t minSeedLen,
											 SMEM *matchArray,
											 int64_t *__numTotalSmem)
{
    int64_t numTotalSmem = *__numTotalSmem;
    SMEM *prevArray = (SMEM *)_mm_malloc(readlength * sizeof(SMEM), 64);
    SMEM *currArray = (SMEM *)_mm_malloc(readlength * sizeof(SMEM), 64);

    numActive = 0;
    uint32_t i;
    // Perform SMEM for original reads
    for(i = 0; i < numReads; i++)
    {
        //printf("i = %d\n", i);
        int x = query_pos_array[i];
        if(x == readlength) continue;
        numActive++;
        int numPrev = 0;
        int numSmem = 0;
        int next_x = x + 1;

        //printf("x = %d\n", x);
        //printf("numPrev = %d\n", numPrev);
        //printf("Forward search\n");
        // Forward search
        SMEM smem;
        smem.rid = rid[i];
        smem.m = x;
        smem.n = x;
                //printf("i: %d, readlen: %d, x: %d\n", i, readlength, x);
        uint8_t a = enc_qdb[i * readlength + x];
        //printf("a = %d\n", a);
        if(a < 4)
        {
            smem.k = count[a];
            smem.l = count[3 - a];
            smem.s = count[a+1] - count[a];
            //printf("[k,l,s] = %d,%d,%d\n", smem.k, smem.l, smem.s);

            int j;
            for(j = x + 1; j < readlength; j++)
            {
                a = enc_qdb[i * readlength + j];
                next_x = j + 1;
                if(a < 4)
                {
                    SMEM smem_ = smem;

                    // Forward extension is backward extension with the BWT of reverse complement
                    smem_.k = smem.l;
                    smem_.l = smem.k;
					SMEM newSmem_ = backwardExt(smem_, 3 - a);
					// SMEM newSmem_ = forwardExt(smem_, 3 - a);
                    SMEM newSmem = newSmem_;
                    newSmem.k = newSmem_.l;
                    newSmem.l = newSmem_.k;
                    newSmem.n = j;

                    //printf("New smem: %u, %u, %u, %u, %u\n", newSmem.m, newSmem.n, newSmem.k, newSmem.l, newSmem.s);
                    if(newSmem.s != smem.s)
                    {
                        //printf("Add to prev: %u, %u, %u, %u, %u\n", smem.m, smem.n, smem.k, smem.l, smem.s);
                        prevArray[numPrev] = smem;
                        numPrev++;
                    }
                    smem = newSmem;
                    _mm_prefetch((const char *)(&cp_occ[(smem.k) >> CP_SHIFT]), _MM_HINT_T0);
                    _mm_prefetch((const char *)(&cp_occ[(smem.l) >> CP_SHIFT]), _MM_HINT_T0);
                    if(newSmem.s < min_intv_array[i])
                    {
                        next_x = j;
                        break;
                    }
                }
                else
                {
                    //printf("else");
                    prevArray[numPrev] = smem;
                    numPrev++;
                    break;
                }
            }
            if(smem.s >= min_intv_array[i])
            {
                //printf("Add to prev1: %u, %u, %u, %u, %u\n", smem.m, smem.n, smem.k, smem.l, smem.s);
                prevArray[numPrev++] = smem;
            }

            SMEM *curr, *prev;
            prev = prevArray;
            curr = currArray;

            int p;
            for(p = 0; p < (numPrev/2); p++)
            {
                SMEM temp = prev[p];
                prev[p] = prev[numPrev - p - 1];
                prev[numPrev - p - 1] = temp;
            }

            //printf("Backward search\n");
            // Backward search
            int cur_j = readlength;
            for(j = x - 1; j >= 0; j--)
            {
                //printf("i = %d, j = %d\n", i, j);
                int numCurr = 0;
                int curr_s = -1;
                a = enc_qdb[i * readlength + j];
                //printf("a = %d\n", a);
                if(a > 3)
                {
                    break;
                }
                for(p = 0; p < numPrev; p++)
                {
                    SMEM smem = prev[p];
                    //printf("smem: %u, %u, %u, %u, %u\n", smem.m, smem.n, smem.k, smem.l, smem.s);
                    SMEM newSmem = backwardExt(smem, a);
                    newSmem.m = j;
                    //printf("newSmem: %u, %u, %u, %u, %u\n", newSmem.m, newSmem.n, newSmem.k, newSmem.l, newSmem.s);

                    if(newSmem.s < min_intv_array[i])
                    {
                        if((numCurr == 0) && (j < cur_j))
                        {
                            cur_j = j;
                            if((smem.s >= min_intv_array[i]) && ((smem.n - smem.m + 1) >= minSeedLen))
                            {
                                //printf("Add to match: %u, %u, %u, %u, %u\n", smem.m, smem.n, smem.k, smem.l, smem.s);
                                matchArray[numTotalSmem + numSmem++] = smem;
                            }
                        }
                    }
                    if((newSmem.s >= min_intv_array[i]) && (newSmem.s != curr_s))
                    {
                        curr_s = newSmem.s;
                        //printf("Add to curr: %u, %u, %u, %u, %u\n", newSmem.m, newSmem.n, newSmem.k, newSmem.l, newSmem.s);
                        curr[numCurr++] = newSmem;
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k) >> CP_SHIFT]), _MM_HINT_T0);
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.l) >> CP_SHIFT]), _MM_HINT_T0);
                    }
                }
                SMEM *temp = prev;
                prev = curr;
                curr = temp;
                numPrev = numCurr;
                if(numCurr == 0)
                {
                    break;
                }
            }
            if(numPrev != 0)
            {
                SMEM smem = prev[0];
                if((smem.s >= min_intv_array[i]) && ((smem.n - smem.m + 1) >= minSeedLen))
                {
                    //printf("Add to match1: %u, %u, %u, %u, %u\n", prev[0].m, prev[0].n, prev[0].k, prev[0].l, prev[0].s);
                    matchArray[numTotalSmem + numSmem++] = smem;
                }
                numPrev = 0;
            }
        }
        query_pos_array[i] = next_x;
        numTotalSmem += numSmem;
    }
    //printf("numTotalSmem = %ld\n", numTotalSmem);
    (*__numTotalSmem) = numTotalSmem;
    _mm_free(prevArray);
    _mm_free(currArray);
    return numActive;
}
#endif

void FMI_search::getSMEMsOnePosOneThread(uint8_t *enc_qdb,
        int16_t *query_pos_array,
        int32_t *min_intv_array,
        int32_t *rid_array,
        int32_t numReads,
        int32_t batch_size,
        int32_t readlength,
        int32_t minSeedLen,
        SMEM *matchArray,
        int64_t *__numTotalSmem)
{
    int64_t numTotalSmem = *__numTotalSmem;
    SMEM prevArray[readlength];

    uint32_t i;
    // Perform SMEM for original reads
    for(i = 0; i < numReads; i++)
    {
        int x = query_pos_array[i];
        int32_t rid = rid_array[i];
        int next_x = x + 1;

        uint8_t a = enc_qdb[rid * readlength + x];

        if(a < 4)
        {
            SMEM smem;
            smem.rid = rid;
            smem.m = x;
            smem.n = x;
            smem.k = count[a];
            smem.l = count[3 - a];
            smem.s = count[a+1] - count[a];
            //printf("[k,l,s] = %d,%d,%d\n", smem.k, smem.l, smem.s);
            int numPrev = 0;

            int j;
            for(j = x + 1; j < readlength; j++)
            {
                a = enc_qdb[rid * readlength + j];
                next_x = j + 1;
                if(a < 4)
                {
                    SMEM smem_ = smem;

                    // Forward extension is backward extension with the BWT of reverse complement
                    smem_.k = smem.l;
                    smem_.l = smem.k;
					SMEM newSmem_ = backwardExt(smem_, 3 - a);
					//SMEM newSmem_ = forwardExt(smem_, 3 - a);
                    SMEM newSmem = newSmem_;
                    newSmem.k = newSmem_.l;
                    newSmem.l = newSmem_.k;
                    newSmem.n = j;

                    int32_t s_neq_mask = newSmem.s != smem.s;

                    //printf("New smem: %u, %u, %u, %u, %u\n", newSmem.m, newSmem.n, newSmem.k, newSmem.l, newSmem.s);
                    //printf("Add to prev: %u, %u, %u, %u, %u\n", smem.m, smem.n, smem.k, smem.l, smem.s);
                    prevArray[numPrev] = smem;
                    numPrev += s_neq_mask;
                    if(newSmem.s < min_intv_array[i])
                    {
                        next_x = j;
                        break;
                    }
                    smem = newSmem;
#ifdef ENABLE_PREFETCH
                    _mm_prefetch((const char *)(&cp_occ[(smem.k) >> CP_SHIFT]), _MM_HINT_T0);
                    _mm_prefetch((const char *)(&cp_occ[(smem.l) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                }
                else
                {
                    break;
                }
            }
            if(smem.s >= min_intv_array[i])
            {
                //printf("Add to prev1: %u, %u, %u, %u, %u\n", smem.m, smem.n, smem.k, smem.l, smem.s);
                prevArray[numPrev] = smem;
                numPrev++;
            }

            SMEM *prev;
            prev = prevArray;

            int p;
            for(p = 0; p < (numPrev/2); p++)
            {
                SMEM temp = prev[p];
                prev[p] = prev[numPrev - p - 1];
                prev[numPrev - p - 1] = temp;
            }

            // Backward search
            int cur_j = readlength;
            for(j = x - 1; j >= 0; j--)
            {
                int numCurr = 0;
                int curr_s = -1;
                a = enc_qdb[rid * readlength + j];

                if(a > 3)
                {
                    break;
                }
                for(p = 0; p < numPrev; p++)
                {
                    SMEM smem = prev[p];
                    SMEM newSmem = backwardExt(smem, a);
                    newSmem.m = j;

                    if((newSmem.s < min_intv_array[i]) && ((smem.n - smem.m + 1) >= minSeedLen))
                    {
                        cur_j = j;
                        //printf("Add to match: %u, %u, %u, %u, %u\n", smem.m, smem.n, smem.k, smem.l, smem.s);
                        matchArray[numTotalSmem++] = smem;
                        break;
                    }
                    if((newSmem.s >= min_intv_array[i]) && (newSmem.s != curr_s))
                    {
                        curr_s = newSmem.s;
                        prev[numCurr++] = newSmem;
#ifdef ENABLE_PREFETCH
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k) >> CP_SHIFT]), _MM_HINT_T0);
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k + newSmem.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                        break;
                    }
                }
                p++;
                for(; p < numPrev; p++)
                {
                    SMEM smem = prev[p];
                    //printf("smem: %u, %u, %u, %u, %u\n", smem.m, smem.n, smem.k, smem.l, smem.s);
                    SMEM newSmem = backwardExt(smem, a);
                    newSmem.m = j;
                    //printf("newSmem: %u, %u, %u, %u, %u\n", newSmem.m, newSmem.n, newSmem.k, newSmem.l, newSmem.s);

                    if((newSmem.s >= min_intv_array[i]) && (newSmem.s != curr_s))
                    {
                        curr_s = newSmem.s;
                        prev[numCurr++] = newSmem;
#ifdef ENABLE_PREFETCH
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k) >> CP_SHIFT]), _MM_HINT_T0);
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k + newSmem.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                    }
                }
                numPrev = numCurr;
                if(numCurr == 0)
                {
                    break;
                }
            }
            if(numPrev != 0)
            {
                SMEM smem = prev[0];
                if(((smem.n - smem.m + 1) >= minSeedLen))
                {
                    //printf("Add to match1: %u, %u, %u, %u, %u\n", prev[0].m, prev[0].n, prev[0].k, prev[0].l, prev[0].s);
                    matchArray[numTotalSmem++] = smem;
                }
                numPrev = 0;
            }
        }
        query_pos_array[i] = next_x;
    }
    (*__numTotalSmem) = numTotalSmem;
}

void FMI_search::getSMEMsAllPosOneThread(uint8_t *enc_qdb,
        int32_t *min_intv_array,
        int32_t *rid_array,
        int32_t numReads,
        int32_t batch_size,
        int32_t readlength,
        int32_t minSeedLen,
        SMEM *matchArray,
        int64_t *__numTotalSmem)
{
    int16_t *query_pos_array = (int16_t *)_mm_malloc(numReads * sizeof(int16_t), 64);

    int32_t i;
    for(i = 0; i < numReads; i++)
        query_pos_array[i] = 0;

    int32_t numActive = numReads;
    (*__numTotalSmem) = 0;

    do
    {
        int32_t head = 0;
        int32_t tail = 0;
        for(head = 0; head < numActive; head++)
        {
            if(query_pos_array[head] < readlength)
            {
                rid_array[tail] = rid_array[head];
                query_pos_array[tail] = query_pos_array[head];
                min_intv_array[tail] = min_intv_array[head];
                tail++;
            }
        }
        getSMEMsOnePosOneThread(enc_qdb,
                    query_pos_array,
                    min_intv_array,
                    rid_array,
                    tail,
                    batch_size,
                    readlength,
                    minSeedLen,
                    matchArray,
                    __numTotalSmem);
        numActive = tail;
    } while(numActive > 0);

    _mm_free(query_pos_array);
}

#if 0
int32_t FMI_search::bwtSeedStrategyOnePosOneThread(uint8_t *enc_qdb,
        int16_t *query_pos_array,
        int32_t *max_intv_array,
        int32_t numReads,
        int32_t rid,
        int32_t readlength,
        int32_t minSeedLen,
        SMEM *matchArray,
        int64_t numTotalSeed)
{
    //numActive = 0;
    uint32_t i;
    // Perform SMEM for original reads
    for(i = 0; i < numReads; i++)
    {
        //printf("i = %d\n", i);
        int x = query_pos_array[i];
        //printf("i = %d, x = %d\n", i, x);
        if(x == readlength) continue;
        //numActive++;
        int next_x = x + 1;

        //printf("x = %d\n", x);
        //printf("Forward search\n");
        // Forward search
        SMEM smem;
        smem.rid = rid;
        smem.m = x;
        smem.n = x;
        uint8_t a = enc_qdb[i * readlength + x];
        //printf("a = %d\n", a);
        if(a < 4)
        {
            smem.k = count[a];
            smem.l = count[3 - a];
            smem.s = count[a+1] - count[a];
            //printf("%d) [a,k,l,s,m,n] = %d,%d,%d,%d,%d,%d\n", i, a, smem.k, smem.l, smem.s, smem.m, smem.n);

            int j;
            for(j = x + 1; j < readlength; j++)
            {
                next_x = j + 1;
                a = enc_qdb[i * readlength + j];
                if(a < 4)
                {
                    SMEM smem_ = smem;

                    // Forward extension is backward extension with the BWT of reverse complement
                    smem_.k = smem.l;
                    smem_.l = smem.k;
                    SMEM newSmem_ = backwardExt(smem_, 3 - a);
                    //SMEM smem = backwardExt(smem, 3 - a);
                    //smem.n = j;
                    SMEM newSmem = newSmem_;
                    newSmem.k = newSmem_.l;
                    newSmem.l = newSmem_.k;
                    newSmem.n = j;
                    smem = newSmem;
#ifdef ENABLE_PREFETCH
                    _mm_prefetch((const char *)(&cp_occ[(smem.k) >> CP_SHIFT]), _MM_HINT_T0);
                    _mm_prefetch((const char *)(&cp_occ[(smem.l) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                    //printf("%d) [a,k,l,s,m,n] = %d,%d,%d,%d,%d,%d\n", i, a, smem.k, smem.l, smem.s, smem.m, smem.n);

                    if((smem.s < max_intv_array[i]) && ((smem.n - smem.m + 1) >= minSeedLen))
                    {
                        //printf("%d) match found\n", i);
                        if(smem.s > 0)
                        {
                            matchArray[numTotalSeed++] = smem;
                        }
                        break;
                    }
                }
                else
                {
                    //printf("else");
                    break;
                }
            }

        }
        query_pos_array[i] = next_x;
    }

    return numTotalSeed;
}
#endif

int64_t FMI_search::bwtSeedStrategyAllPosOneThread(uint8_t *enc_qdb,
        int32_t *max_intv_array,
        int32_t numReads,
        int32_t readlength,
        int32_t minSeedLen,
        SMEM *matchArray)
{
    int32_t i;

    int64_t numTotalSeed = 0;

    for(i = 0; i < numReads; i++)
    {
        int16_t x = 0;
        while(x < readlength)
        {
            //printf("i = %d\n", i);
            //printf("i = %d, x = %d\n", i, x);
            int next_x = x + 1;

            //printf("x = %d\n", x);
            //printf("Forward search\n");
            // Forward search
            SMEM smem;
            smem.rid = i;
            smem.m = x;
            smem.n = x;
            uint8_t a = enc_qdb[i * readlength + x];
            //printf("a = %d\n", a);
            if(a < 4)
            {
                smem.k = count[a];
                smem.l = count[3 - a];
                smem.s = count[a+1] - count[a];
                //printf("%d) [a,k,l,s,m,n] = %d,%d,%d,%d,%d,%d\n", i, a, smem.k, smem.l, smem.s, smem.m, smem.n);

                int j;
                for(j = x + 1; j < readlength; j++)
                {
                    next_x = j + 1;
                    a = enc_qdb[i * readlength + j];
                    if(a < 4)
                    {
                        SMEM smem_ = smem;

                        // Forward extension is backward extension with the BWT of reverse complement
                        smem_.k = smem.l;
                        smem_.l = smem.k;
                        SMEM newSmem_ = backwardExt(smem_, 3 - a);
                        //SMEM smem = backwardExt(smem, 3 - a);
                        //smem.n = j;
                        SMEM newSmem = newSmem_;
                        newSmem.k = newSmem_.l;
                        newSmem.l = newSmem_.k;
                        newSmem.n = j;
                        smem = newSmem;
#ifdef ENABLE_PREFETCH
                        _mm_prefetch((const char *)(&cp_occ[(smem.k) >> CP_SHIFT]), _MM_HINT_T0);
                        _mm_prefetch((const char *)(&cp_occ[(smem.l) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                        //printf("%d) [a,k,l,s,m,n] = %d,%d,%d,%d,%d,%d\n", i, a, smem.k, smem.l, smem.s, smem.m, smem.n);

                        if((smem.s < max_intv_array[i]) && ((smem.n - smem.m + 1) >= minSeedLen))
                        {
                            //printf("%d) match found\n", i);
                            if(smem.s > 0)
                            {
                                matchArray[numTotalSeed++] = smem;
                            }
                            break;
                        }
                    }
                    else
                    {
                        //printf("else");
                        break;
                    }
                }

            }
            x = next_x;
        }
    }
    return numTotalSeed;
}



#if 0
void FMI_search::getSMEMsOnePos(uint8_t *enc_qdb,
        int16_t *query_pos_array,
        int16_t *min_seed_len_array,
        int32_t numReads,
        int32_t batch_size,
        int32_t readlength,
        int32_t nthreads,
        SMEM *matchArray,
        int64_t *numTotalSmem)
{
    SMEM *prevArray = (SMEM *)_mm_malloc(nthreads * readlength * sizeof(SMEM), 64);
    SMEM *currArray = (SMEM *)_mm_malloc(nthreads * readlength * sizeof(SMEM), 64);


#pragma omp parallel num_threads(nthreads)
    {
        int tid = omp_get_thread_num();
        numTotalSmem[tid] = 0;
        SMEM *myPrevArray = prevArray + tid * readlength;
        SMEM *myCurrArray = prevArray + tid * readlength;

        int32_t perThreadQuota = (numReads + (nthreads - 1)) / nthreads;
        int32_t first = tid * perThreadQuota;
        int32_t last  = (tid + 1) * perThreadQuota;
        if(last > numReads) last = numReads;
        SMEM *myMatchArray = matchArray + first * readlength;

        uint32_t i;
        // Perform SMEM for original reads
        for(i = first; i < last; i++)
        {
            //printf("i = %d\n", i);
            int x = query_pos_array[i];
            int numPrev = 0;
            int numSmem = 0;

            //while (x >= 0)
            {
                //printf("x = %d\n", x);
                //printf("numPrev = %d\n", numPrev);
                //printf("Forward search\n");
                // Forward search
                SMEM smem;
                smem.rid = i;
                smem.m = x;
                smem.n = x;
                uint8_t a = enc_qdb[i * readlength + x];
                //printf("a = %d\n", a);
                if(a > 3)
                {
                    x--;
                    continue;
                }
                smem.k = count[a];
                smem.l = count[3 - a];
                smem.s = count[a+1] - count[a];
                //printf("[k,l,s] = %d,%d,%d\n", smem.k, smem.l, smem.s);

                int j;
                for(j = x + 1; j < readlength; j++)
                {
                    a = enc_qdb[i * readlength + j];
                    if(a < 4)
                    {
                        SMEM smem_ = smem;

                        // Forward extension is backward extension with the BWT of reverse complement
                        smem_.k = smem.l;
                        smem_.l = smem.k;
                        SMEM newSmem_ = backwardExt(smem_, 3 - a);
                        SMEM newSmem = newSmem_;
                        newSmem.k = newSmem_.l;
                        newSmem.l = newSmem_.k;
                        newSmem.n = j;

                        //printf("New smem: %u, %u, %u, %u, %u\n", newSmem.m, newSmem.n, newSmem.k, newSmem.l, newSmem.s);
                        if(newSmem.s != smem.s)
                        {
                            //printf("Add to prev: %u, %u, %u, %u, %u\n", smem.m, smem.n, smem.k, smem.l, smem.s);
                            myPrevArray[numPrev] = smem;
                            numPrev++;
                        }
                        smem = newSmem;
                        if(newSmem.s == 0)
                        {
                            break;
                        }
                    }
                    else
                    {
                        //printf("else");
                        myPrevArray[numPrev] = smem;
                        numPrev++;
                        break;
                    }
                }
                if(smem.s != 0)
                {
                    //printf("Add to prev1: %u, %u, %u, %u, %u\n", smem.m, smem.n, smem.k, smem.l, smem.s);
                    myPrevArray[numPrev++] = smem;
                }

                SMEM *curr, *prev;
                prev = myPrevArray;
                curr = myCurrArray;

                int p;
                for(p = 0; p < (numPrev/2); p++)
                {
                    SMEM temp = prev[p];
                    prev[p] = prev[numPrev - p - 1];
                    prev[numPrev - p - 1] = temp;
                }

                int next_x = x - 1;
                //printf("Backward search\n");
                // Backward search
                int cur_j = readlength;
                for(j = x - 1; j >= 0; j--)
                {
                    //printf("j = %d\n", j);
                    int numCurr = 0;
                    int curr_s = -1;
                    a = enc_qdb[i * readlength + j];
                    //printf("a = %d\n", a);
                    if(a > 3)
                    {
                        next_x = j - 1;
                        break;
                    }
                    for(p = 0; p < numPrev; p++)
                    {
                        SMEM smem = prev[p];
                        //printf("smem: %u, %u, %u, %u, %u\n", smem.m, smem.n, smem.k, smem.l, smem.s);
                        SMEM newSmem = backwardExt(smem, a);
                        newSmem.m = j;
                        //printf("newSmem: %u, %u, %u, %u, %u\n", newSmem.m, newSmem.n, newSmem.k, newSmem.l, newSmem.s);

                        if(newSmem.s == 0)
                        {
                            if((numCurr == 0) && (j < cur_j))
                            {
                                cur_j = j;
                                //printf("Add to match: %u, %u, %u, %u, %u\n", smem.m, smem.n, smem.k, smem.l, smem.s);
                                if((smem.n - smem.m + 1) >= minSeedLen)
                                    myMatchArray[numTotalSmem[tid] + numSmem++] = smem;
                            }
                        }
                        if((newSmem.s != 0) && (newSmem.s != curr_s))
                        {
                            curr_s = newSmem.s;
                            //printf("Add to curr: %u, %u, %u, %u, %u\n", newSmem.m, newSmem.n, newSmem.k, newSmem.l, newSmem.s);
                            curr[numCurr++] = newSmem;
                        }
                    }
                    SMEM *temp = prev;
                    prev = curr;
                    curr = temp;
                    numPrev = numCurr;
                    if(numCurr == 0)
                    {
                        next_x = j;
                        break;
                    }
                    else
                    {
                        next_x = j - 1;
                    }
                }
                if(numPrev != 0)
                {
                    //printf("Add to match1: %u, %u, %u, %u, %u\n", prev[0].m, prev[0].n, prev[0].k, prev[0].l, prev[0].s);
                    SMEM smem = prev[0];
                    if((smem.n - smem.m + 1) >= minSeedLen)
                        myMatchArray[numTotalSmem[tid] + numSmem++] = smem;
                    numPrev = 0;
                }
                x = next_x;
            }
            numTotalSmem[tid] += numSmem;
        }
    }

    _mm_free(prevArray);
    _mm_free(currArray);
}
#endif

void FMI_search::getSMEMs(uint8_t *enc_qdb,
        int32_t numReads,
        int32_t batch_size,
        int32_t readlength,
        int32_t minSeedLen,
        int32_t nthreads,
        SMEM *matchArray,
        int64_t *numTotalSmem)
{
    SMEM *prevArray = (SMEM *)_mm_malloc(nthreads * readlength * sizeof(SMEM), 64);
    SMEM *currArray = (SMEM *)_mm_malloc(nthreads * readlength * sizeof(SMEM), 64);


#pragma omp parallel num_threads(nthreads)
    {
        int tid = omp_get_thread_num();
        numTotalSmem[tid] = 0;
        SMEM *myPrevArray = prevArray + tid * readlength;
        SMEM *myCurrArray = prevArray + tid * readlength;

        int32_t perThreadQuota = (numReads + (nthreads - 1)) / nthreads;
        int32_t first = tid * perThreadQuota;
        int32_t last  = (tid + 1) * perThreadQuota;
        if(last > numReads) last = numReads;
        SMEM *myMatchArray = matchArray + first * readlength;

        uint32_t i;
        // Perform SMEM for original reads
        for(i = first; i < last; i++)
        {
            //printf("i = %d\n", i);
            int x = readlength - 1;
            int numPrev = 0;
            int numSmem = 0;

            while (x >= 0)
            {
                //printf("x = %d\n", x);
                //printf("numPrev = %d\n", numPrev);
                //printf("Forward search\n");
                // Forward search
                SMEM smem;
                smem.rid = i;
                smem.m = x;
                smem.n = x;
                uint8_t a = enc_qdb[i * readlength + x];
                //printf("a = %d\n", a);
                if(a > 3)
                {
                    x--;
                    continue;
                }
                smem.k = count[a];
                smem.l = count[3 - a];
                smem.s = count[a+1] - count[a];
                //printf("[k,l,s] = %d,%d,%d\n", smem.k, smem.l, smem.s);

                int j;
                for(j = x + 1; j < readlength; j++)
                {
                    a = enc_qdb[i * readlength + j];
                    if(a < 4)
                    {
                        SMEM smem_ = smem;

                        // Forward extension is backward extension with the BWT of reverse complement
                        smem_.k = smem.l;
                        smem_.l = smem.k;
                        SMEM newSmem_ = backwardExt(smem_, 3 - a);
                        SMEM newSmem = newSmem_;
                        newSmem.k = newSmem_.l;
                        newSmem.l = newSmem_.k;
                        newSmem.n = j;

                        //printf("New smem: %u, %u, %u, %u, %u\n", newSmem.m, newSmem.n, newSmem.k, newSmem.l, newSmem.s);
                        if(newSmem.s != smem.s)
                        {
                            //printf("Add to prev: %u, %u, %u, %u, %u\n", smem.m, smem.n, smem.k, smem.l, smem.s);
                            myPrevArray[numPrev] = smem;
                            numPrev++;
                        }
                        smem = newSmem;
                        if(newSmem.s == 0)
                        {
                            break;
                        }
                    }
                    else
                    {
                        //printf("else");
                        myPrevArray[numPrev] = smem;
                        numPrev++;
                        break;
                    }
                }
                if(smem.s != 0)
                {
                    //printf("Add to prev1: %u, %u, %u, %u, %u\n", smem.m, smem.n, smem.k, smem.l, smem.s);
                    myPrevArray[numPrev++] = smem;
                }

                SMEM *curr, *prev;
                prev = myPrevArray;
                curr = myCurrArray;

                int p;
                for(p = 0; p < (numPrev/2); p++)
                {
                    SMEM temp = prev[p];
                    prev[p] = prev[numPrev - p - 1];
                    prev[numPrev - p - 1] = temp;
                }

                int next_x = x - 1;
                //printf("Backward search\n");
                // Backward search
                int cur_j = readlength;
                for(j = x - 1; j >= 0; j--)
                {
                    //printf("j = %d\n", j);
                    int numCurr = 0;
                    int curr_s = -1;
                    a = enc_qdb[i * readlength + j];
                    //printf("a = %d\n", a);
                    if(a > 3)
                    {
                        next_x = j - 1;
                        break;
                    }
                    for(p = 0; p < numPrev; p++)
                    {
                        SMEM smem = prev[p];
                        //printf("smem: %u, %u, %u, %u, %u\n", smem.m, smem.n, smem.k, smem.l, smem.s);
                        SMEM newSmem = backwardExt(smem, a);
                        newSmem.m = j;
                        //printf("newSmem: %u, %u, %u, %u, %u\n", newSmem.m, newSmem.n, newSmem.k, newSmem.l, newSmem.s);

                        if(newSmem.s == 0)
                        {
                            if((numCurr == 0) && (j < cur_j))
                            {
                                cur_j = j;
                                //printf("Add to match: %u, %u, %u, %u, %u\n", smem.m, smem.n, smem.k, smem.l, smem.s);
                                if((smem.n - smem.m + 1) >= minSeedLen)
                                    myMatchArray[numTotalSmem[tid] + numSmem++] = smem;
                            }
                        }
                        if((newSmem.s != 0) && (newSmem.s != curr_s))
                        {
                            curr_s = newSmem.s;
                            //printf("Add to curr: %u, %u, %u, %u, %u\n", newSmem.m, newSmem.n, newSmem.k, newSmem.l, newSmem.s);
                            curr[numCurr++] = newSmem;
                        }
                    }
                    SMEM *temp = prev;
                    prev = curr;
                    curr = temp;
                    numPrev = numCurr;
                    if(numCurr == 0)
                    {
                        next_x = j;
                        break;
                    }
                    else
                    {
                        next_x = j - 1;
                    }
                }
                if(numPrev != 0)
                {
                    //printf("Add to match1: %u, %u, %u, %u, %u\n", prev[0].m, prev[0].n, prev[0].k, prev[0].l, prev[0].s);
                    SMEM smem = prev[0];
                    if((smem.n - smem.m + 1) >= minSeedLen)
                        myMatchArray[numTotalSmem[tid] + numSmem++] = smem;
                    numPrev = 0;
                }
                x = next_x;
            }
            numTotalSmem[tid] += numSmem;
        }
    }

    _mm_free(prevArray);
    _mm_free(currArray);
}


int compare_smem(const void *a, const void *b)
{
    SMEM *pa = (SMEM *)a;
    SMEM *pb = (SMEM *)b;

    if(pa->rid < pb->rid)
        return -1;
    if(pa->rid > pb->rid)
        return 1;

    if(pa->m < pb->m)
        return -1;
    if(pa->m > pb->m)
        return 1;
    if(pa->n > pb->n)
        return -1;
    if(pa->n < pb->n)
        return 1;
    return 0;
}

void FMI_search::sortSMEMs(SMEM *matchArray,
        int64_t numTotalSmem[],
        int32_t numReads,
        int32_t readlength,
        int nthreads)
{
    int tid;
    int32_t perThreadQuota = (numReads + (nthreads - 1)) / nthreads;
    for(tid = 0; tid < nthreads; tid++)
    {
        int32_t first = tid * perThreadQuota;
        SMEM *myMatchArray = matchArray + first * readlength;
        qsort(myMatchArray, numTotalSmem[tid], sizeof(SMEM), compare_smem);
    }
}


SMEM FMI_search::backwardExt(SMEM smem, uint8_t a)
{
    //beCalls++;
    uint8_t b;

    int64_t k[4], l[4], s[4];
    for(b = 0; b < 4; b++)
    {
        __m256i b256;
        int64_t sp = (int64_t)(smem.k);
        int64_t ep = (int64_t)(smem.k) + (int64_t)(smem.s);
#if ((!__AVX512BW__) && (!__AVX2__))
        int32_t bwt_id;
        GET_OCC_SCALAR(sp, b, occ_id_sp, y_sp, occ_sp, bwt_str_sp);
        GET_OCC_SCALAR(ep, b, occ_id_ep, y_ep, occ_ep, bwt_str_ep);
#else
        b256 = _mm256_load_si256((const __m256i *)(c_bcast_array + b * 64));
        GET_OCC(sp, b, b256, occ_id_sp, y_sp, occ_sp, bwt_str_sp, bwt_sp_vec, mask_sp_vec, mask_sp);
        GET_OCC(ep, b, b256, occ_id_ep, y_ep, occ_ep, bwt_str_ep, bwt_ep_vec, mask_ep_vec, mask_ep);
#endif
        k[b] = count[b] + occ_sp;
        s[b] = occ_ep - occ_sp;
    }

    int64_t sentinel_offset = 0;
    // if((smem.k < sentinel_index) && ((smem.k + smem.s) >= sentinel_index)) sentinel_offset = 1;
	if((smem.k <= sentinel_index) && ((smem.k + smem.s) > sentinel_index)) sentinel_offset = 1;
    l[3] = smem.l + sentinel_offset;
    l[2] = l[3] + s[3];
    l[1] = l[2] + s[2];
    l[0] = l[1] + s[1];

    smem.k = k[a];
    smem.l = l[a];
    smem.s = s[a];
    return smem;
}

int64_t FMI_search::get_sa_entry(int64_t pos)
{
    return sa[pos];
}

void FMI_search::get_sa_entries(int64_t *posArray, int64_t *coordArray, uint32_t count, int32_t nthreads)
{
    uint32_t i;
    //int PFD = 100;
    #pragma omp parallel for num_threads(nthreads)
    for(i = 0; i < count; i++)
    {
	  // _mm_prefetch((const char *)(sa + posArray[i + PFD]), _MM_HINT_T0);
	  coordArray[i] = sa[posArray[i]];
    }
}

void FMI_search::get_sa_entries(SMEM *smemArray, int64_t *coordArray, int32_t *coordCountArray, uint32_t count, int32_t max_occ)
{
    uint32_t i;
    //int PFD = 100;
    int32_t totalCoordCount = 0;
    for(i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        int64_t hi = smem.k + smem.s;
        int64_t step = (smem.s > max_occ) ? smem.s / max_occ : 1;
        int64_t j;
        for(j = smem.k; (j < hi) && (c < max_occ); j+=step, c++)
        {
            // _mm_prefetch((const char *)(sa + posArray[i + PFD]), _MM_HINT_T0);
            coordArray[totalCoordCount + c] = sa[j];
        }
        coordCountArray[i] = c;
        totalCoordCount += c;
    }
}
