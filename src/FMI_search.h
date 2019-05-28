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

#ifndef _FMI_SEARCH_H
#define _FMI_SEARCH_H

#include<stdlib.h>
#include<stdint.h>
#include <string.h>
//#include <xmmintrin.h>
#include <immintrin.h>
// #include <omp.h>
#include "bntseq.h"
#include "read_index_ele.h"
#include "bwa.h"

#if ((!__AVX2__))

#define CP_BLOCK_SIZE 64
#define CP_MASK 63
#define CP_SHIFT 6
#define BIT_DATA_TYPE uint64_t
#define PADDING 24

#if defined(__clang__) || defined(__GNUC__)
static inline int _mm_countbits_64(unsigned long x) {
	return __builtin_popcountl(x);
}
#endif

#define \
GET_OCC(pp, c, occ_id_pp, y_pp, occ_pp, bwt_str_bit0_pp, bwt_str_bit1_pp, bit0_cmp_pp, bit1_cmp_pp, mismatch_mask_pp) \
                int64_t occ_id_pp = pp >> CP_SHIFT; \
                int64_t y_pp = pp & CP_MASK; \
                int64_t occ_pp = cp_occ[occ_id_pp].cp_count[c]; \
                if(y_pp > 0) \
                { \
                BIT_DATA_TYPE bwt_str_bit0_pp = cp_occ[occ_id_pp].bwt_str_bit0; \
                BIT_DATA_TYPE bwt_str_bit1_pp = cp_occ[occ_id_pp].bwt_str_bit1; \
                BIT_DATA_TYPE bit0_cmp_pp = bwt_str_bit0_pp ^ base_mask[c][0]; \
                BIT_DATA_TYPE bit1_cmp_pp = bwt_str_bit1_pp ^ base_mask[c][1]; \
                uint64_t mismatch_mask_pp = bit0_cmp_pp | bit1_cmp_pp | cp_occ[occ_id_pp].dollar_mask; \
                mismatch_mask_pp = mismatch_mask_pp >> (CP_BLOCK_SIZE - y_pp); \
                occ_pp += y_pp - _mm_countbits_64(mismatch_mask_pp); \
                }

typedef struct checkpoint_occ
{
    BIT_DATA_TYPE bwt_str_bit0;
    BIT_DATA_TYPE bwt_str_bit1;
    BIT_DATA_TYPE dollar_mask;
    uint32_t cp_count[4];
    uint8_t  pad[PADDING];
}CP_OCC;

#else

#define CP_BLOCK_SIZE 32
#define CP_MASK 31
#define CP_SHIFT 5

#if defined(__clang__) || defined(__GNUC__)
static inline int _mm_countbits_32(unsigned x) {
	return __builtin_popcount(x);
}
#endif

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

typedef struct checkpoint_occ
{
    uint8_t  bwt_str[CP_BLOCK_SIZE];
    uint32_t cp_count[4];
    uint8_t  pad[16];
}CP_OCC;

#endif

typedef struct smem_struct
{
#ifdef DEBUG
	uint64_t info; // for debug
#endif
    uint32_t rid;
    uint32_t m, n;
    int64_t k, l, s;
}SMEM;

#define SAL_PFD 16

class FMI_search: public indexEle
{
public:
	FMI_search(char *fmi_file_name);
	~FMI_search();
	//int64_t beCalls;
	
	void getSMEMs(uint8_t *enc_qdb,
				  int32_t numReads,
				  int32_t batch_size,
				  int32_t readlength,
				  int32_t minSeedLengh,
				  int32_t numthreads,
				  SMEM *matchArray,
				  int64_t *numTotalSmem);
	
	void getSMEMsOnePosOneThread(uint8_t *enc_qdb,
								 int16_t *query_pos_array,
								 int32_t *min_intv_array,
								 int32_t *rid_array,
								 int32_t numReads,
								 int32_t batch_size,
								 const bseq1_t *seq_,
								 int32_t *query_cum_len_ar,
								 int32_t  max_readlength,
								 int32_t minSeedLen,
								 SMEM *matchArray,
								 int64_t *__numTotalSmem);
	
	void getSMEMsAllPosOneThread(uint8_t *enc_qdb,
								 int32_t *min_intv_array,
								 int32_t *rid_array,
								 int32_t numReads,
								 int32_t batch_size,
								 const bseq1_t *seq_,
								 int32_t *query_cum_len_ar,
								 int32_t max_readlength,
								 int32_t minSeedLen,
								 SMEM *matchArray,
								 int64_t *__numTotalSmem);
		
	
	int64_t bwtSeedStrategyAllPosOneThread(uint8_t *enc_qdb,
										   int32_t *max_intv_array,
										   int32_t numReads,
										   const bseq1_t *seq_,
										   int32_t *query_cum_len_ar,
										   int32_t minSeedLen,
										   SMEM *matchArray);
		
	void sortSMEMs(SMEM *matchArray,
				   int64_t numTotalSmem[],
				   int32_t numReads,
				   int32_t readlength,
				   int nthreads);
	int64_t get_sa_entry(int64_t pos);
	void get_sa_entries(int64_t *posArray, int64_t *coordArray, uint32_t count, int32_t nthreads);
	void get_sa_entries(SMEM *smemArray, int64_t *coordArray, int32_t *coordCountArray, uint32_t count, int32_t max_occ);

	int64_t reference_seq_len;
	uint32_t sentinel_index;
private:
        int64_t count[5];
        uint32_t *sa_ls_word;
        int8_t *sa_ms_byte;
        CP_OCC *cp_occ;

#if ((!__AVX2__))
        BIT_DATA_TYPE base_mask[4][2];
#else
        uint8_t c_bcast_array[256] __attribute__((aligned(64)));
#endif

        SMEM backwardExt(SMEM smem, uint8_t a);
};

#endif
