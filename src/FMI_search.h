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
#include <omp.h>
#include "bntseq.h"
#include "read_index_ele.h"

#define CP_BLOCK_SIZE 32
#define CP_MASK 31
#define CP_SHIFT 5

typedef struct smem_struct
{
#ifdef DEBUG
	uint64_t info; // for debug
#endif
    uint32_t rid;
    uint32_t m, n;
    int64_t k, l, s;
}SMEM;

typedef struct checkpoint_occ
{
    uint8_t  bwt_str[CP_BLOCK_SIZE];
    uint32_t cp_count[4];
    uint8_t  pad[16];
}CP_OCC;

#if BWA_OTHER_ELE
class FMI_search: public indexEle
#else
class FMI_search
#endif
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

#if 0
        int32_t getSMEMsOnePosOneThread2(
                uint8_t *enc_qdb,
                int16_t *query_pos_array,
                int32_t *min_intv_array,
				int32_t *rid,
                int32_t numReads,
                int32_t numActive,
                int32_t batch_size,
                int32_t readlength,
                int32_t minSeedLen,
                SMEM *matchArray,
                int64_t *__numTotalSmem);
#endif
        	   
        void getSMEMsOnePosOneThread(
                uint8_t *enc_qdb,
                int16_t *query_pos_array,
                int32_t *min_intv_array,
                int32_t *rid_array,
                int32_t numReads,
                int32_t batch_size,
                int32_t readlength,
                int32_t minSeedLen,
                SMEM *matchArray,
                int64_t *__numTotalSmem);

        void getSMEMsAllPosOneThread(
                uint8_t *enc_qdb,
                int32_t *min_intv_array,
                int32_t *rid_array,
                int32_t numReads,
                int32_t batch_size,
                int32_t readlength,
                int32_t minSeedLen,
                SMEM *matchArray,
                int64_t *__numTotalSmem);

#if 0
        int32_t bwtSeedStrategyOnePosOneThread(
                uint8_t *enc_qdb,
                int16_t *query_pos_array,
                int32_t *max_intv_array,
                int32_t numReads,
                int32_t rid,
                int32_t readlength,
                int32_t minSeedLen,
                SMEM *matchArray,
                int64_t numTotalSeed);
#endif

        int64_t bwtSeedStrategyAllPosOneThread(
                uint8_t *enc_qdb,
                int32_t *max_intv_array,
                int32_t numReads,
                int32_t readlength,
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
        int64_t *sa;
        CP_OCC *cp_occ;

        uint8_t c_bcast_array[256] __attribute__((aligned(64)));

        SMEM backwardExt(SMEM smem, uint8_t a);
};

#endif
