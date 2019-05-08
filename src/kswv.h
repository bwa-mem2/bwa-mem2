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

Authors: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
*****************************************************************************************/

#ifndef _KSWV_H_
#define  _KSWV_H_

#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
// #include <emmintrin.h>
#include <omp.h>
#include "macro.h"
#if !MAINY
#include "ksw.h"
// #include "bwamem.h"
#include "bandedSWA.h"
#else
#include <immintrin.h>
#endif

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif


#define MAX_SEQ_LEN_REF_SAM 2048
#define MAX_SEQ_LEN_QER_SAM 256

#if MAINY
#define KSW_XBYTE  0x10000
#define KSW_XSTOP  0x20000
#define KSW_XSUBO  0x40000
#define KSW_XSTART 0x80000

// extern uint64_t proc_freq, tprof[LIM_R][LIM_C];

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif


// #define MAX_SEQ_LEN8 256

#define MAX_SEQ_LEN_EXT 256

#if __AVX512BW__
#define SIMD_WIDTH8 64
#define SIMD_WIDTH16 32
#endif

#define max(x, y) ((x)>(y)?(x):(y))
#define min(x, y) ((x)>(y)?(y):(x))

#define MAX_NUM_PAIRS 1000000
#define MAX_NUM_PAIRS_ALLOC 20000

#define DEFAULT_AMBIG -1

typedef struct dnaSeqPair
{
	int32_t idr, idq, id;
	int32_t len1, len2;
	int32_t h0;
	int seqid, regid;
	// int32_t score, tle, gtle, qle;
	// int32_t gscore, max_off;
	int score; // best score
	int te, qe; // target end and query end
	int score2, te2; // second best score and ending position on the target
	int tb, qb; // target start and query start
}SeqPair;

typedef struct {
	int qlen, slen;
	uint8_t shift, mdiff, max, size;
	__m128i *qp, *H0, *H1, *E, *Hmax;
} kswq_t;

typedef struct {
	int score; // best score
	int te, qe; // target end and query end
	int score2, te2; // second best score and ending position on the target
	int tb, qb; // target start and query start
} kswr_t;


const kswr_t g_defr = { 0, -1, -1, -1, -1, -1, -1 };

#define __max_16(ret, xx) do { \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 8)); \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 4)); \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 2)); \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 1)); \
    	(ret) = _mm_extract_epi16((xx), 0) & 0x00ff; \
	} while (0)

#define DP  6
#define DP1 7
#define DP2 8
#define DP3 9

#else

//extern struct _kswq_t;
extern const kswr_t g_defr;

#endif

// typedef struct {
// 	int qlen, slen;
// 	uint8_t shift, mdiff, max, size;
// 	__m512i *qp, *H0, *H1, *E, *Hmax;
// } kswqi_t;


class kswv {
public:
	uint64_t SW_cells;

	kswv(const int o_del, const int e_del, const int o_ins,
		 const int e_ins, int8_t w_match, int8_t w_mismatch,
		 int numThreads);
	
	~kswv();

	void getScores8(SeqPair *pairArray,
					uint8_t *seqBufRef,
					uint8_t *seqBufQer,
					kswr_t* aln,
					int32_t numPairs,
					uint16_t numThreads,
					int phase);

	void getScores16(SeqPair *pairArray,
					 uint8_t *seqBufRef,
					 uint8_t *seqBufQer,
					 kswr_t* aln,
					 int32_t numPairs,
					 uint16_t numThreads,
					 int phase);

	void kswvScalaWrapper(SeqPair *seqPairArray,
						  uint8_t *seqBufRef,
						  uint8_t *seqBufQer,
						  kswr_t* aln,
						  int numPairs,
						  int nthreads);

	kswq_t* ksw_qinit(int size, int qlen, uint8_t *query, int m, const int8_t *mat);
	// kswqi_t* ksw_qinit_intra(int size, int qlen, uint8_t *query, int m, const int8_t *mat);
	// void kswvBatchWrapper16_intra(SeqPair *pairArray,
	// 							  uint8_t *seqBufRef,
	// 							  uint8_t *seqBufQer,
	// 							  int32_t numPairs,
	// 							  uint16_t numThreads);
	// 
	// kswr_t kswv512_16_intra(uint8_t seq1SoA[],
	// 						kswqi_t *q,
	// 						int16_t nrow,
	// 						int16_t ncol,
	// 						SeqPair *p,
	// 						uint16_t tid,
	// 						int32_t numPairs,
	// 						int &maxi);
	
private:
#if __AVX512BW__
	void kswvBatchWrapper8(SeqPair *pairArray,
						   uint8_t *seqBufRef,
						   uint8_t *seqBufQer,
						   kswr_t* aln,
						   int32_t numPairs,
						   uint16_t numThreads,
						   int phase);

	int kswv512_u8(uint8_t seq1SoA[],
				   uint8_t seq2SoA[],
				   int16_t nrow,
				   int16_t ncol,
				   SeqPair *p,
				   kswr_t *aln,
				   int po_ind,
				   uint16_t tid,
				   int32_t numPairs,
				   int phase);
	
	void kswvBatchWrapper16(SeqPair *pairArray,
							uint8_t *seqBufRef,
							uint8_t *seqBufQer,
							kswr_t* aln,
							int32_t numPairs,
							uint16_t numThreads,
							int phase);
	
	int kswv512_16_exp(int16_t seq1SoA[],
					   int16_t seq2SoA[],
					   int16_t nrow,
					   int16_t ncol,
					   SeqPair *p,
					   kswr_t* aln,
					   int po_ind,
					   uint16_t tid,
					   int32_t numPairs,
					   int phase);

	void kswv512_16(int16_t seq1SoA[],
				   int16_t seq2SoA[],
				   int16_t nrow,
				   int16_t ncol,
				   SeqPair *p,
				   kswr_t* aln,
				   int po_ind,
				   uint16_t tid,
				   int32_t numPairs);
#endif
	
	kswr_t kswvScalar_u8(kswq_t *q, int tlen, const uint8_t *target,
						int _o_del, int _e_del, int _o_ins, int _e_ins,
						int xtra);  // the first gap costs -(_o+_e)

	kswr_t kswvScalar_u8_exp(kswq_t *q, int tlen, const uint8_t *target,
							 int _o_del, int _e_del, int _o_ins, int _e_ins,
							 int xtra);  // the first gap costs -(_o+_e)
	
	kswr_t kswvScalar_i16(kswq_t *q, int tlen, const uint8_t *target,
						  int _o_del, int _e_del, int _o_ins, int _e_ins,
						  int xtra); // the first gap costs -(_o+_e)
	
	kswr_t kswvScalar_i16_exp(kswq_t *q, int tlen, const uint8_t *target,
							  int _o_del, int _e_del, int _o_ins, int _e_ins,
							  int xtra); // the first gap costs -(_o+_e)
		

	void bwa_fill_scmat(int8_t mat[25]);
		
	int m;
	int o_del, o_ins, e_del, e_ins;
	// const int8_t *mat;

	int8_t w_match;
	int8_t w_mismatch;
	int8_t w_open;
	int8_t w_extend;
	int8_t w_ambig;
	uint8_t *F8;
	uint8_t *H8_0, *H8_max, *H8_1;
	uint8_t *rowMax8;
	
	int16_t *F16;
	int16_t *H16_0, *H16_max, *H16_1;
	// int16_t *qp16;
	int16_t *rowMax16;

	int g_qmax;
	int64_t sort1Ticks;
	int64_t setupTicks;
	int64_t swTicks;
	int64_t sort2Ticks;
};

#endif
