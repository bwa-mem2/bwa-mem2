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

#ifndef SCALAR_BANDEDSWA_HPP
#define SCALAR_BANDEDSWA_HPP

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <immintrin.h>
#include <assert.h>
#include "macro.h"

#define MAX_SEQ_LEN_REF 256
#define MAX_SEQ_LEN_QER 128
#define MAX_SEQ_LEN_EXT 256
#define MAX_NUM_PAIRS 10000000
#define MAX_NUM_PAIRS_ALLOC 20000

// used in BSW and SAM-SW
#define DEFAULT_AMBIG -1


// SIMD_WIDTH in bits
// AVX2
#if ((!__AVX512BW__) & (__AVX2__))
#define SIMD_WIDTH8 32
#define SIMD_WIDTH16 16
#endif

// AVX512
#if __AVX512BW__
#define SIMD_WIDTH8 64
#define SIMD_WIDTH16 32
#endif

#if ((!__AVX512BW__) & (!__AVX2__) & (__SSE2__))
#define SIMD_WIDTH8 16
#define SIMD_WIDTH16 8
#endif

// Scalar
#if ((!__AVX512BW__) & (!__AVX2__) & (!__SSE2__))
#define SIMD_WIDTH8 1
#define SIMD_WIDTH16 1
#endif

#define MAX_LINE_LEN 256
#define MAX_SEQ_LEN8 128
#define MAX_SEQ_LEN16 32768
#define MATRIX_MIN_CUTOFF 0
#define LOW_INIT_VALUE -128
#define SORT_BLOCK_SIZE 16384
#define max(x, y) ((x)>(y)?(x):(y))
#define min(x, y) ((x)>(y)?(y):(x))

typedef struct dnaSeqPair
{
	int32_t idr, idq, id;
	int32_t len1, len2;
	int32_t h0;
	int seqid, regid;
	int32_t score, tle, gtle, qle;
	int32_t gscore, max_off;

}SeqPair;


typedef struct dnaOutScore
{
	int32_t score, tle, gtle, qle;
	int32_t gscore, max_off;
} OutScore;

typedef struct {
	int32_t h, e;
} eh_t;


class BandedPairWiseSW {
	
public:
	uint64_t SW_cells;

	BandedPairWiseSW(const int o_del, const int e_del, const int o_ins,
					 const int e_ins, const int zdrop,
					 const int end_bonus, const int8_t *mat_,
					 int8_t w_match, int8_t w_mismatch, int numThreads);
	~BandedPairWiseSW();
	// Scalar code section
	int scalarBandedSWA(int qlen, const uint8_t *query, int tlen,
						const uint8_t *target, int w,
						int h0, int *_qle, int *_tle,
						int *_gtle, int *_gscore,
						int *_max_off);

	void scalarBandedSWAWrapper(SeqPair *seqPairArray,
								uint8_t *seqBufRef,
								uint8_t *seqBufQer,
								int numPairs,
								int nthreads,
								int w);

#if ((!__AVX512BW__) & (!__AVX2__) & (__SSE2__))
	// AVX256 is not updated for banding and separate ins/del in the inner loop.
	// 8 bit vector code section	
	void getScores8(SeqPair *pairArray,
					uint8_t *seqBufRef,
					uint8_t *seqBufQer,
					int32_t numPairs,
					uint16_t numThreads,
					int8_t w);

	void smithWatermanBatchWrapper8(SeqPair *pairArray,
								   uint8_t *seqBufRef,
								   uint8_t *seqBufQer,
								   int32_t numPairs,
								   uint16_t numThreads,
								   int8_t w);

	void smithWaterman128_8(uint8_t seq1SoA[],
							uint8_t seq2SoA[],
							uint8_t nrow,
							uint8_t ncol,
							SeqPair *p,
							uint8_t h0[],
							uint16_t tid,
							int32_t numPairs,
							int zdrop,
							uint8_t w,
							uint8_t qlen[],
							uint8_t myband[]);
	// 16 bit vector code section
	void getScores16(SeqPair *pairArray,
					 uint8_t *seqBufRef,
					 uint8_t *seqBufQer,
					 int32_t numPairs,
					 uint16_t numThreads,
					 int8_t w);

	void smithWatermanBatchWrapper16(SeqPair *pairArray,
									 uint8_t *seqBufRef,
									 uint8_t *seqBufQer,
									 int32_t numPairs,
									 uint16_t numThreads,
									 int8_t w);
	
	void smithWaterman128_16(uint16_t seq1SoA[],
							 uint16_t seq2SoA[],
							 uint16_t nrow,
							 uint16_t ncol,
							 SeqPair *p,
							 uint16_t h0[],
							 uint16_t tid,
							 int32_t numPairs,
							 int zdrop,
							 uint16_t w,
							 uint16_t qlen[],
							 uint16_t myband[]);
	
#endif  //SSE2

#if ((!__AVX512BW__) & (__AVX2__))
	// AVX256 is not updated for banding and separate ins/del in the inner loop.
	// 8 bit vector code section	
	void getScores8(SeqPair *pairArray,
					uint8_t *seqBufRef,
					uint8_t *seqBufQer,
					int32_t numPairs,
					uint16_t numThreads,
					int8_t w);

	void smithWatermanBatchWrapper8(SeqPair *pairArray,
								   uint8_t *seqBufRef,
								   uint8_t *seqBufQer,
								   int32_t numPairs,
								   uint16_t numThreads,
								   int8_t w);

	void smithWaterman256_8(uint8_t seq1SoA[],
							uint8_t seq2SoA[],
							uint8_t nrow,
							uint8_t ncol,
							SeqPair *p,
							uint8_t h0[],
							uint16_t tid,
							int32_t numPairs,
							int zdrop,
							uint8_t w,
							uint8_t qlen[],
							uint8_t myband[]);
	// 16 bit vector code section
	void getScores16(SeqPair *pairArray,
					 uint8_t *seqBufRef,
					 uint8_t *seqBufQer,
					 int32_t numPairs,
					 uint16_t numThreads,
					 int8_t w);

	void smithWatermanBatchWrapper16(SeqPair *pairArray,
									 uint8_t *seqBufRef,
									 uint8_t *seqBufQer,
									 int32_t numPairs,
									 uint16_t numThreads,
									 int8_t w);
	
	void smithWaterman256_16(uint16_t seq1SoA[],
							 uint16_t seq2SoA[],
							 uint16_t nrow,
							 uint16_t ncol,
							 SeqPair *p,
							 uint16_t h0[],
							 uint16_t tid,
							 int32_t numPairs,
							 int zdrop,
							 uint16_t w,
							 uint16_t qlen[],
							 uint16_t myband[]);
	
#endif  //axv2

#if __AVX512BW__
	// 8 bit vector code section	
	void getScores8(SeqPair *pairArray,
					uint8_t *seqBufRef,
					uint8_t *seqBufQer,
					int32_t numPairs,
					uint16_t numThreads,
					int8_t w);

	void smithWatermanBatchWrapper8(SeqPair *pairArray,
								   uint8_t *seqBufRef,
								   uint8_t *seqBufQer,
								   int32_t numPairs,
								   uint16_t numThreads,
								   int8_t w);

	void smithWaterman512_8(uint8_t seq1SoA[],
							uint8_t seq2SoA[],
							uint8_t nrow,
							uint8_t ncol,
							SeqPair *p,
							uint8_t h0[],
							uint16_t tid,
							int32_t numPairs,
							int zdrop,
							uint8_t w,
							uint8_t qlen[],
							uint8_t myband[]);

	// 16 bit vector code section
	void getScores16(SeqPair *pairArray,
					 uint8_t *seqBufRef,
					 uint8_t *seqBufQer,
					 int32_t numPairs,
					 uint16_t numThreads,
					 int8_t w);

	void smithWatermanBatchWrapper16(SeqPair *pairArray,
									 uint8_t *seqBufRef,
									 uint8_t *seqBufQer,
									 int32_t numPairs,
									 uint16_t numThreads,
									 int8_t w);
	
	void smithWaterman512_16(uint16_t seq1SoA[],
							 uint16_t seq2SoA[],
							 uint16_t nrow,
							 uint16_t ncol,
							 SeqPair *p,
							 uint16_t h0[],
							 uint16_t tid,
							 int32_t numPairs,
							 int zdrop,
							 uint16_t w,
							 uint16_t qlen[],
							 uint16_t myband[]);
#endif

	int64_t getTicks();
	
private:
	int m;
	int end_bonus, zdrop;
	int o_del, o_ins, e_del, e_ins;
	const int8_t *mat;

	int8_t w;  //band size
	int8_t w_match;
	int8_t w_mismatch;
	int8_t w_open;
	int8_t w_extend;
	int8_t w_ambig;
	int8_t *F8_;
	int8_t *H8_, *H8__;
	
	int16_t *F16_;
	int16_t *H16_, *H16__;

	int64_t sort1Ticks;
	int64_t setupTicks;
	int64_t swTicks;
	int64_t sort2Ticks;
};


#define DP  4
#define DP1 5
#define DP2 6
#define DP3 7


// Some Debug macros

#define DEB1									\
	SW_cells2 += nrow * ncol * SIMD_WIDTH;

#define DEB2										\
	prof[2][tid] += SIMD_WIDTH*4;


#define DEB3													\
	uint32_t ex = _mm256_movemask_epi8(exit0), msk=0x1;			\
	for (int l=0;l<32;l++) {									\
		if ((ex & (msk << l)) == 1) {									\
			if (p[l].len2 < end) prof[3][tid] += (end - p[l].len2) * 4;	\
		}																\
	}																	\
	msk = 0x1;															\
	uint64_t gval1 = 0, gval2 = 0, gval3 = 0, gval4 = 0;				\
	for (int l=0;l<32;l++) {											\
		if ((ex & (msk << l)) == 0) gval1 ++;							\
	}																	\
	prof[3][tid] += (end-beg)*4*gval1;				

#define DEB4										\
	ex = _mm256_movemask_epi8(exit0);				\
	msk=0x1;										\
	for (int l=0;l<32;l++) {						\
		if ((ex & (msk << l)) == 0) gval2 ++;		\
	}												\
	prof[3][tid] += (end-beg)*3*(gval2 - gval1);	\

#define DEB5									\
	ex = _mm256_movemask_epi8(exit0);			\
	msk=0x1;									\
	for (int l=0;l<32;l++) {						\
		if ((ex & (msk << l)) == 0) gval3 ++;		\
	}												\
	prof[3][tid] += (end-beg)*2*(gval3-gval2);		\
	
#define DEB6									\
	ex = _mm256_movemask_epi8(exit0);			\
	msk=0x1;									\
	for (int l=0;l<32;l++) {						\
		if ((ex & (msk << l)) == 0) gval4 ++;		\
	}												\
	prof[3][tid] += (end-beg)*(gval4-gval3);		\
	

#define STAT1									\
	int8_t mat[200][200]= {0};

#define STAT2											\
	mat[i][0] = *(H2 + i * SIMD_WIDTH + lane);				\
	/*mat[i+1][0] = *(H2 + (i+1) * SIMD_WIDTH + lane);		\
	mat[i+2][0] = *(H2 + (i+2) * SIMD_WIDTH + lane);		\
	mat[i+3][0] = *(H2 + (i+3) * SIMD_WIDTH + lane);*/

#define STAT3												\
	int8_t temp[SIMD_WIDTH]  __attribute((aligned(64)));	\
	mat[i][j+1] = *(H + j * SIMD_WIDTH + lane);				\
	/*_mm256_store_si256((__m256i *)(temp), h10);			\
	mat[i+1][j+1] = temp[lane];								\
	_mm256_store_si256((__m256i *)(temp), h20);				\
	mat[i+2][j+1] = temp[lane];								\
	_mm256_store_si256((__m256i *)(temp), h30);				\
	mat[i+3][j+1] = temp[lane];*/

#define STAT4									\
	for (int l=0; l<p[lane].len1; l++)						\
		printf("%d ", *(seq1SoA + l*SIMD_WIDTH + lane));	\
	printf("\n");											\
	for (int l=0; l<p[lane].len2; l++)						\
		printf("%d ", *(seq2SoA + l*SIMD_WIDTH + lane));	\
	printf("\n");											\
	printf("nrow: %d, ncol: %d\n", nrow, ncol);				\
	for(i = 0; i <=60; i++) {								\
		printf("%2d ", i);									\
		for(j = 0; j <=50 ; j++)							\
			printf("%2d ", mat[i][j]);						\
		printf("\n");										\
	}														\
	printf("\n");											\
	exit(0);

#endif
