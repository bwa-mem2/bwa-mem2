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

Authors: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>
*****************************************************************************************/

#include <immintrin.h>
#include <string.h>
#include <unistd.h>
#include "kswv.h"
//#include "bwamem.h"
#include "limits.h"


// ------------------------------------------------------------------------------------
// MACROs for vector code
#if MAINY
uint64_t prof[10][112], data, SW_cells2;
FILE *fsam;
#else
extern uint64_t prof[10][112], data, SW_cells2;
#endif

#define AMBIG_ 4  // ambiguous base

// for 16 bit
#define DUMMY1_ 4
#define DUMMY2_ 5
#define DUMMY3 26

#define AMBR16 15
#define AMBQ16 16

// for 8-bit
#define DUMMY8 8
#define DUMMY5 5
#define AMBRQ 0xFF

#define AMBR 4
#define AMBQ 8

#ifndef max_
#define max_(x, y) ((x)>(y)?(x):(y))
#endif
#ifndef min_
#define min_(x, y) ((x)<(y)?(x):(y))
#endif

int spot = 42419;
int lim = 20;
uint8_t *query, *target, *query_;
int len1, len2;
int iid = 0;
FILE *fp;
// -----------------------------------------------------------------------------------


#define MAIN_SAM_CODE8_OPT(s1, s2, h00, h11, e11, f11, f21, max512, sft512) \
	{																	\
		__m512i sbt11, xor11, or11;										\
		xor11 = _mm512_xor_si512(s1, s2);								\
		sbt11 = _mm512_shuffle_epi8(permSft512, xor11);					\
		__mmask64 cmpq = _mm512_cmpeq_epu8_mask(s2, five512);			\
		sbt11 = _mm512_mask_blend_epi8(cmpq, sbt11, sft512);			\
		or11 =  _mm512_or_si512(s1, s2);								\
		__mmask64 cmp = _mm512_movepi8_mask(or11);						\
		sbt11 = _mm512_mask_blend_epi8(cmp, sbt11, mismatch512);		\
		__m512i m11 = _mm512_adds_epu8(h00, sbt11);						\
		m11 = _mm512_subs_epu8(m11, sft512);							\
		h11 = _mm512_max_epu8(m11, e11);								\
		h11 = _mm512_max_epu8(h11, f11);								\
		__mmask64 cmp0 = _mm512_cmpgt_epu8_mask(h11, imax512);			\
		imax512 = _mm512_max_epu8(imax512, h11);						\
		iqe512 = _mm512_mask_blend_epi8(cmp0, iqe512, l512);			\
		__m512i gapE512 = _mm512_subs_epu8(h11, oe_ins512);				\
		e11 = _mm512_subs_epu8(e11, e_ins512);							\
		e11 = _mm512_max_epu8(gapE512, e11);							\
		__m512i gapD512 = _mm512_subs_epu8(h11, oe_del512);				\
		f21 = _mm512_subs_epu8(f11, e_del512);							\
		f21 = _mm512_max_epu8(gapD512, f21);							\
	}

#define MAIN_SAM_CODE16_OPT(s1, s2, h00, h11, e11, f11, f21, max512)	\
	{																	\
		__m512i sbt11, xor11;											\
		xor11 = _mm512_xor_si512(s1, s2);								\
		sbt11 = _mm512_permutexvar_epi16(xor11, perm512);				\
		__m512i m11 = _mm512_add_epi16(h00, sbt11);						\
		h11 = _mm512_max_epi16(m11, e11);								\
		h11 = _mm512_max_epi16(h11, f11);								\
		h11 = _mm512_max_epi16(h11, zero512);							\
		__mmask32 cmp0 = _mm512_cmpgt_epi16_mask(h11, imax512);			\
		imax512 = _mm512_max_epi16(imax512, h11);						\
		iqe512 = _mm512_mask_blend_epi16(cmp0, iqe512, l512);			\
		__m512i gapE512 = _mm512_sub_epi16(h11, oe_ins512);				\
		e11 = _mm512_sub_epi16(e11, e_ins512);							\
		e11 = _mm512_max_epi16(gapE512, e11);							\
		__m512i gapD512 = _mm512_sub_epi16(h11, oe_del512);				\
		f21 = _mm512_sub_epi16(f11, e_del512);						    \
		f21 = _mm512_max_epi16(gapD512, f21);							\
	}


#define MAIN_SAM_CODE16_DEB(s1, s2, h00, h11, e11, f11, f21, max512)		\
	{																	\
		__mmask32 dum11 = _mm512_cmpeq_epi16_mask(s2, dum512);			\
		__mmask32 cmp11 = _mm512_cmpeq_epi16_mask(s1, s2);				\
		__m512i sbt11 = _mm512_mask_blend_epi16(cmp11, mismatch512, match512); \
		__m512i tmp512 = _mm512_max_epu16(s1, s2);						\
		cmp11 = _mm512_movepi16_mask(tmp512);							\
		sbt11 = _mm512_mask_blend_epi16(cmp11, sbt11, w_ambig_512);		\
		sbt11 = _mm512_mask_blend_epi16(dum11, sbt11, zero512);			\
		__m512i m11 = _mm512_add_epi16(h00, sbt11);						\
		h11 = _mm512_max_epi16(m11, e11);								\
		h11 = _mm512_max_epi16(h11, f11);								\
		if (i == 3 && j == 150) {										\
			_mm512_store_si512((__m512i *)(temp), h00);					\
			_mm512_store_si512((__m512i *)(temp1), sbt11);				\
			_mm512_store_si512((__m512i *)(temp2), s1);					\
			_mm512_store_si512((__m512i *)(temp3), s2);					\
			printf("%d %d %d %d\n", temp[lane], temp1[lane], temp2[lane], temp3[lane]);	\
		}																\
		max512 = _mm512_max_epi16(max512, h11);							\
		__m512i gapE512 = _mm512_subs_epu16(h11, oe_ins512);			\
		e11 = _mm512_subs_epu16(e11, e_ins512);							\
		e11 = _mm512_max_epi16(gapE512, e11);							\
		__m512i gapD512 = _mm512_subs_epu16(h11, oe_del512);			\
		f21 = _mm512_subs_epu16(f11, e_del512);							\
		f21 = _mm512_max_epi16(gapD512, f21);							\
	}


// constructor
kswv::kswv(const int o_del, const int e_del, const int o_ins,
		   const int e_ins, int8_t w_match, int8_t w_mismatch,
		   int numThreads)
{
	SW_cells = 0;
	// this->mat = mat_;
	this->m = 5;
	this->o_del = o_del;
	this->o_ins = o_ins;
	this->e_del = e_del;
	this->e_ins = e_ins;
	
	this->w_match	 = w_match;
	this->w_mismatch = w_mismatch;
	this->w_open	 = o_del;  // redundant, used in vector code.
	this->w_extend	 = e_del;  // redundant, used in vector code.
	this->w_ambig	 = DEFAULT_AMBIG;
	this->g_qmax = max_(w_match, w_mismatch);
	this->g_qmax = max_(this->g_qmax, w_ambig);
	
	this->swTicks = 0;
	setupTicks = 0;
	sort1Ticks = 0;
	swTicks = 0;
	sort2Ticks = 0;
	//printf("match: %d, mismatch: %d, open: %d, extend: %d, ambig: %d\n",
	//	   w_match, this->w_mismatch, w_open, w_extend, w_ambig);
    F16	    = (int16_t *)_mm_malloc(MAX_SEQ_LEN_QER_SAM * SIMD_WIDTH16 * numThreads * sizeof(int16_t), 64);
    H16_0   = (int16_t *)_mm_malloc(MAX_SEQ_LEN_QER_SAM * SIMD_WIDTH16 * numThreads * sizeof(int16_t), 64);
	H16_1   = (int16_t *)_mm_malloc(MAX_SEQ_LEN_QER_SAM * SIMD_WIDTH16 * numThreads * sizeof(int16_t), 64);
	H16_max = (int16_t *)_mm_malloc(MAX_SEQ_LEN_QER_SAM * SIMD_WIDTH16 * numThreads * sizeof(int16_t), 64);
	rowMax16 = (int16_t *)_mm_malloc(MAX_SEQ_LEN_REF_SAM * SIMD_WIDTH16 * numThreads * sizeof(int16_t), 64);
	// 8-bit
	//F8	   = (uint8_t *)_mm_malloc(MAX_SEQ_LEN_QER_SAM * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 64);
    //H8_0   = (uint8_t *)_mm_malloc(MAX_SEQ_LEN_QER_SAM * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 64);
	//H8_1   = (uint8_t *)_mm_malloc(MAX_SEQ_LEN_QER_SAM * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 64);
	//H8_max = (uint8_t *)_mm_malloc(MAX_SEQ_LEN_QER_SAM * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 64);
	//rowMax8 = (uint8_t *)_mm_malloc(MAX_SEQ_LEN_REF_SAM * SIMD_WIDTH8 * numThreads* sizeof(uint8_t), 64);
	F8 = (uint8_t*) F16;
	H8_0 = (uint8_t*) H16_0;
	H8_1 = (uint8_t*) H16_1;
	H8_max = (uint8_t*) H16_max;
	rowMax8 = (uint8_t*) rowMax16;
}

// destructor 
kswv::~kswv() {
	_mm_free(F16); _mm_free(H16_0); _mm_free(H16_max); _mm_free(H16_1);
	_mm_free(rowMax16);
	
	//_mm_free(F8); _mm_free(H8_0); _mm_free(H8_max); _mm_free(H8_1);
	//_mm_free(rowMax8);	
}

//void kswv::bwa_fill_scmat(int a, int b, int ambig, int8_t mat[25]) {
void kswv::bwa_fill_scmat(int8_t mat[25]) {
	int a = this->w_match;
	int b = this->w_mismatch;
	int ambig = this->w_ambig;
	
	int i, j, k;
	for (i = k = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j)
			mat[k++] = i == j? a : b;
		mat[k++] = ambig; // ambiguous base
	}
	for (j = 0; j < 5; ++j) mat[k++] = ambig;
}

void kswv::getScores8(SeqPair *pairArray,
					  uint8_t *seqBufRef,
					  uint8_t *seqBufQer,
					  kswr_t* aln,					  
					  int32_t numPairs,
					  uint16_t numThreads,
					  int phase)
{
    int64_t startTick, endTick;
    //F8	   = (uint8_t *)_mm_malloc(MAX_SEQ_LEN_QER_SAM * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 64);
    //H8_0   = (uint8_t *)_mm_malloc(MAX_SEQ_LEN_QER_SAM * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 64);
	//H8_1   = (uint8_t *)_mm_malloc(MAX_SEQ_LEN_QER_SAM * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 64);
	//H8_max = (uint8_t *)_mm_malloc(MAX_SEQ_LEN_QER_SAM * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 64);
	//rowMax8 = (uint8_t *)_mm_malloc(MAX_SEQ_LEN_REF_SAM * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 64);
	//rowMax16 = (int16_t *)_mm_malloc(MAX_SEQ_LEN_REF_SAM * SIMD_WIDTH16 * numThreads * sizeof(int16_t), 64);
	
    // startTick = __rdtsc();
	kswvBatchWrapper8(pairArray, seqBufRef, seqBufQer, aln, numPairs, numThreads, phase);
    // endTick = __rdtsc();

	// _mm_free(F8); _mm_free(H8_0); _mm_free(H8_max); _mm_free(H8_1);
	// _mm_free(rowMax8);
}

#define PFD_ 2
void kswv::kswvBatchWrapper8(SeqPair *pairArray,
							 uint8_t *seqBufRef,
							 uint8_t *seqBufQer,
							 kswr_t* aln,
							 int32_t numPairs,
							 uint16_t numThreads,
							 int phase)
{
	// printf("numThreads: %d %d\n", numThreads, omp_get_thread_num());

	int64_t st1, st2, st3, st4, st5;
    // st1 = __rdtsc();
    uint8_t *seq1SoA = NULL;
	seq1SoA = (uint8_t *)_mm_malloc(MAX_SEQ_LEN_REF_SAM * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 64);
	
    uint8_t *seq2SoA = NULL;
	seq2SoA = (uint8_t *)_mm_malloc(MAX_SEQ_LEN_QER_SAM * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 64);
	
    assert(seq1SoA != NULL);
    assert(seq2SoA != NULL);

    int32_t ii;
    int32_t roundNumPairs = ((numPairs + SIMD_WIDTH8 - 1) / SIMD_WIDTH8 ) * SIMD_WIDTH8;
	// assert(roundNumPairs < BATCH_SIZE * SEEDS_PER_READ);
    for(ii = numPairs; ii < roundNumPairs; ii++)
    {
		pairArray[ii].regid = ii;
        pairArray[ii].id = ii;
        pairArray[ii].len1 = 0;
        pairArray[ii].len2 = 0;
    }
		
    // st2 = __rdtsc();	
#if SORT_PAIRS     // disbaled in bwa-mem2 (only used in separate benchmark sw code)
	{
    // Sort the sequences according to decreasing order of lengths
    SeqPair *tempArray = (SeqPair *)_mm_malloc(SORT_BLOCK_SIZE * numThreads *
											   sizeof(SeqPair), 64);
    int16_t *hist = (int16_t *)_mm_malloc((MAX_SEQ_LEN_QER_SAM + 32) * numThreads *
										  sizeof(int16_t), 64);

#pragma omp parallel num_threads(numThreads)
    {
        int32_t tid = omp_get_thread_num();
        SeqPair *myTempArray = tempArray + tid * SORT_BLOCK_SIZE;
        int16_t *myHist = hist + tid * (MAX_SEQ_LEN_QER_SAM + 32);

#pragma omp for
        for(ii = 0; ii < roundNumPairs; ii+=SORT_BLOCK_SIZE)
        {
            int32_t first, last;
            first = ii;
            last  = ii + SORT_BLOCK_SIZE;
            if(last > roundNumPairs) last = roundNumPairs;
            sortPairsLen(pairArray + first, last - first, myTempArray, myHist);
        }
    }
    _mm_free(hist);
	}
#endif
	
    // st3 = __rdtsc();

//#pragma omp parallel num_threads(numThreads)
    {
        int32_t i;
        // uint16_t tid = omp_get_thread_num();
		uint16_t tid = 0;
        uint8_t *mySeq1SoA = seq1SoA + tid * MAX_SEQ_LEN_REF_SAM * SIMD_WIDTH8;
        uint8_t *mySeq2SoA = seq2SoA + tid * MAX_SEQ_LEN_QER_SAM * SIMD_WIDTH8;
        uint8_t *seq1;
        uint8_t *seq2;
				
		int nstart = 0, nend = numPairs;

		
//#pragma omp for schedule(dynamic, 128)
		for(i = nstart; i < nend; i+=SIMD_WIDTH8)
		{
			// prof[4][0]++;
            int32_t j, k;
            int maxLen1 = 0;
            int maxLen2 = 0;
			
            for(j = 0; j < SIMD_WIDTH8; j++)
            {
                SeqPair sp = pairArray[i + j];
#if MAINY				
				seq1 = seqBufRef + (int64_t)sp.id * MAX_SEQ_LEN_REF_SAM;
#else
				seq1 = seqBufRef + sp.idr;
#endif
                for(k = 0; k < sp.len1; k++)
                {
                    // mySeq1SoA[k * SIMD_WIDTH8 + j] = (seq1[k] == AMBIG_?0xFF:seq1[k]);
					mySeq1SoA[k * SIMD_WIDTH8 + j] = (seq1[k] == AMBIG_ ? AMBR:seq1[k]);
                }
                if(maxLen1 < sp.len1) maxLen1 = sp.len1;
            }
            for(j = 0; j < SIMD_WIDTH8; j++)
            {
                SeqPair sp = pairArray[i + j];
                for(k = sp.len1; k <= maxLen1; k++) //removed "="
                {
					//mySeq1SoA[k * SIMD_WIDTH8 + j] = DUMMY1_;
					mySeq1SoA[k * SIMD_WIDTH8 + j] = 0xFF;
                }
            }

            for(j = 0; j < SIMD_WIDTH8; j++)
            {				
                SeqPair sp = pairArray[i + j];
#if MAINY
				seq2 = seqBufQer + (int64_t)sp.id * MAX_SEQ_LEN_QER_SAM;
#else
				seq2 = seqBufQer + sp.idq;
#endif
				assert(sp.len2 < MAX_SEQ_LEN_QER_SAM);
				int quanta = 16 - sp.len2 % 16;  // based on SSE2-8 bit lane
                for(k = 0; k < sp.len2; k++)
                {
					// mySeq2SoA[k * SIMD_WIDTH8 + j] = (seq2[k]==AMBIG_?0xFF:seq2[k]);
					mySeq2SoA[k * SIMD_WIDTH8 + j] = (seq2[k]==AMBIG_? AMBQ:seq2[k]);
                }
				for(k = sp.len2; k < sp.len2 + quanta; k++) {
					// mySeq2SoA[k * SIMD_WIDTH16 + j] = DUMMY3;
					mySeq2SoA[k * SIMD_WIDTH8 + j] = DUMMY5;  // SSE2 qunata
				}
                if(maxLen2 < (sp.len2 + quanta)) maxLen2 = sp.len2 + quanta;
            }
			
            for(j = 0; j < SIMD_WIDTH8; j++)
            {
                SeqPair sp = pairArray[i + j];
				int quanta = 16 - sp.len2 % 16;  // based on SSE2-8 bit lane				
                for(k = sp.len2 + quanta; k <= maxLen2; k++)
                {
					// mySeq2SoA[k * SIMD_WIDTH8 + j] = DUMMY2_;
					mySeq2SoA[k * SIMD_WIDTH8 + j] = 0xFF;
                }
            }
			
            kswv512_u8(mySeq1SoA, mySeq2SoA,
					   maxLen1, maxLen2,
					   pairArray + i,
					   aln, i,
					   tid,
					   numPairs,
					   phase);
        }
    }
	
    // st4 = __rdtsc();
	
#if SORT_PAIRS     // disbaled in bwa-mem2 (only used in separate benchmark sw code)
	{
    // Sort the sequences according to increasing order of id
#pragma omp parallel num_threads(numThreads)
    {
        int32_t tid = omp_get_thread_num();
        SeqPair *myTempArray = tempArray + tid * SORT_BLOCK_SIZE;

#pragma omp for
        for(ii = 0; ii < roundNumPairs; ii+=SORT_BLOCK_SIZE)
        {
            int32_t first, last;
            first = ii;
            last  = ii + SORT_BLOCK_SIZE;
            if(last > roundNumPairs) last = roundNumPairs;
            sortPairsId(pairArray + first, first, last - first, myTempArray);
        }
    }
    _mm_free(tempArray);
	}
#endif
	
    // st5 = __rdtsc();
    setupTicks = st2 - st1;
    sort1Ticks = st3 - st2;
    swTicks = st4 - st3;
    sort2Ticks = st5 - st4;

	// free mem
	_mm_free(seq1SoA);
	_mm_free(seq2SoA);
	
    return;
}
int kswv::kswv512_u8(uint8_t seq1SoA[],
					 uint8_t seq2SoA[],
					 int16_t nrow,
					 int16_t ncol,
					 SeqPair *p,
					 kswr_t *aln,
					 int po_ind,
					 uint16_t tid,
					 int32_t numPairs,
					 int phase)
{
#ifdef __AVX512BW__
	int m_b, n_b;
	uint8_t minsc[SIMD_WIDTH8] = {0}, endsc[SIMD_WIDTH8] = {0};
	uint64_t *b;

	__m512i zero512 = _mm512_setzero_si512();
	__m512i one512  = _mm512_set1_epi8(1);
	// kswr_t r;  //output struct
	
	// initialization
	// r = g_defr;   // check for its role later.
	// int xtra = 0;
	// minsc = (xtra & KSW_XSUBO)? xtra & 0xffff : 0x10000;
	// endsc = (xtra & KSW_XSTOP)? xtra & 0xffff : 0x10000;
	m_b = n_b = 0; b = 0;

	// calculate sft

	int8_t temp[SIMD_WIDTH8] __attribute((aligned(64))) = {0};

	uint8_t shift = 127, mdiff = 0, qmax_;
	mdiff = max_(this->w_match, (int8_t) this->w_mismatch);
	mdiff = max_(mdiff, (int8_t) this->w_ambig);
	shift = min_(this->w_match, (int8_t) this->w_mismatch);
	shift = min_((int8_t) shift, this->w_ambig);

	qmax_ = mdiff;
	shift = 256 - (uint8_t) shift;
	mdiff += shift;
	
	temp[0] = this->w_match;                                   // states: 1. matches
	temp[1] = temp[2] = temp[3] =  this->w_mismatch;           // 2. mis-matches
	temp[4] = temp[5] = temp[6] = temp[7] =  this->w_ambig;    // 3. beyond boundary
	temp[8] = temp[9] = temp[10] = temp[11] = this->w_ambig;   // 4. 0 - sse2 region
	temp[12] = this->w_ambig;                                  // 5. ambig

	// printf("shift: %d\n", shift);
	for (int i=0; i<16; i++) // for shuffle_epi8
		temp[i] += shift;

	int pos = 0;
	for (int i=16; i<SIMD_WIDTH8; i++) {
		temp[i] = temp[pos++];
		if (pos % 16 == 0) pos = 0;
	}
	//for (int i=0; i<SIMD_WIDTH8; i++)
	//	printf("%d ", temp[i]);
	//printf("\n\n");
	//exit(0);
	
	__m512i permSft512 = _mm512_load_si512(temp);
	__m512i sft512 = _mm512_set1_epi8(shift);
	__m512i cmax512 = _mm512_set1_epi8(255);
	
	// __m512i minsc512, endsc512;
	__mmask64 minsc_msk_a = 0x0000, endsc_msk_a = 0x0000;
	int val = 0;
	for (int i=0; i<SIMD_WIDTH8; i++) {
		int xtra = p[i].h0;
		val = (xtra & KSW_XSUBO)? xtra & 0xffff : 0x10000;
		if (val <= 255) {
			minsc[i] = val;
			//if (val < 0)
			//	minsc[i] = 0;
			minsc_msk_a |= (0x1L << i);
		}
		// msc_mask;
		val = (xtra & KSW_XSTOP)? xtra & 0xffff : 0x10000;
		if (val <= 255) {
			endsc[i] = val;
			//if (val < 0)
			//	endsc[i] = 0;
			endsc_msk_a |= (0x1L << i);
		}
	}
	// printf();
	__m512i minsc512 = _mm512_load_si512((__m512i*) minsc);
	__m512i endsc512 = _mm512_load_si512((__m512i*) endsc);
	   
    __m512i mismatch512	= _mm512_set1_epi8(this->w_mismatch + shift);
	__m512i e_del512	= _mm512_set1_epi8(this->e_del);
	__m512i oe_del512	= _mm512_set1_epi8(this->o_del + this->e_del);
	__m512i e_ins512	= _mm512_set1_epi8(this->e_ins);
	__m512i oe_ins512	= _mm512_set1_epi8(this->o_ins + this->e_ins);
	__m512i five512	    = _mm512_set1_epi8(DUMMY5);	// ambig mapping element
	__m512i gmax512		= zero512; // exit1 = zero512;
	__m512i te512		= zero512;	// change to -1
	__m512i te512_		= zero512;	// change to -1
	__mmask64 exit0 = 0xFFFFFFFFFFFFFFFF;

	tid = 0;  // no threading for now !!
	uint8_t	*H0		= H8_0 + tid * SIMD_WIDTH8 * MAX_SEQ_LEN_QER_SAM;
	uint8_t	*H1		= H8_1 + tid * SIMD_WIDTH8 * MAX_SEQ_LEN_QER_SAM;
	uint8_t	*Hmax	= H8_max + tid * SIMD_WIDTH8 * MAX_SEQ_LEN_QER_SAM;
	uint8_t	*F		= F8 + tid * SIMD_WIDTH8 * MAX_SEQ_LEN_QER_SAM;
	uint8_t	*rowMax	= rowMax8 + tid * SIMD_WIDTH8 * MAX_SEQ_LEN_REF_SAM;
	
	
	for (int i=0; i <=ncol; i++) {
		_mm512_store_si512((__m512*) (H0 + i * SIMD_WIDTH8), zero512);
		_mm512_store_si512((__m512*) (Hmax + i * SIMD_WIDTH8), zero512);
		_mm512_store_si512((__m512*) (F + i * SIMD_WIDTH8), zero512);
	}

	__m512i max512 = zero512, imax512, pimax512 = zero512;
	__mmask64 mask512 = 0x0000;
	__mmask64 minsc_msk = 0x0000;

	__m512i qe512 = _mm512_set1_epi8(0);
	_mm512_store_si512((__m512i *)(H0), zero512);
	_mm512_store_si512((__m512i *)(H1), zero512);

	int i, limit = nrow;
	for (i=0; i < nrow; i++)
	{
        __m512i e11 = zero512;
        __m512i h00, h11, h10, s1;
		__m512i i512 = _mm512_set1_epi16(i);
		int j ;
		
		s1 = _mm512_load_si512((__m512i *)(seq1SoA + (i + 0) * SIMD_WIDTH8));
		h10 = zero512;
		imax512 = zero512;
		__m512i iqe512 = _mm512_set1_epi8(-1);

		__m512i l512 = zero512;
		for (j=0; j<ncol; j++)
		{
			__m512i f11, s2, f21;
			h00 = _mm512_load_si512((__m512i *)(H0 + j * SIMD_WIDTH8));  // check for col "0"
			s2  = _mm512_load_si512((__m512i *)(seq2SoA + (j) * SIMD_WIDTH8));
			f11 = _mm512_load_si512((__m512i *)(F + (j+1) * SIMD_WIDTH8));

			MAIN_SAM_CODE8_OPT(s1, s2, h00, h11, e11, f11, f21, max512, sft512);

			_mm512_store_si512((__m512i *)(H1 + (j + 1) * SIMD_WIDTH8), h11);  // check for col "0"
			_mm512_store_si512((__m512i *)(F + (j + 1)* SIMD_WIDTH8), f21);
			l512 = _mm512_add_epi8(l512, one512);
		}

		// TODO: Block I
		if (i > 0) {
			__mmask64 msk64 = _mm512_cmpgt_epu8_mask(imax512, pimax512);
			msk64 |= mask512;
			pimax512 = _mm512_mask_blend_epi8(msk64, pimax512, zero512);
			pimax512 = _mm512_mask_blend_epi8(minsc_msk, zero512, pimax512);
			pimax512 = _mm512_mask_blend_epi8(exit0, zero512, pimax512);
			
			_mm512_store_si512((__m512i *) (rowMax + (i-1)*SIMD_WIDTH8), pimax512);
			mask512 = ~msk64;
		}
		pimax512 = imax512;
		minsc_msk = _mm512_cmpge_epu8_mask(imax512, minsc512);
		minsc_msk &= minsc_msk_a;

		// Block II: gmax, te
		__mmask64 cmp0 = _mm512_cmpgt_epu8_mask(imax512, gmax512);
		cmp0 &= exit0;
		gmax512 = _mm512_mask_blend_epi8(cmp0, gmax512, imax512);
		te512 = _mm512_mask_blend_epi16((__mmask32)cmp0, te512,  i512);
		te512_ = _mm512_mask_blend_epi16((__mmask32) (cmp0 >> SIMD_WIDTH16), te512_,  i512);
		qe512 = _mm512_mask_blend_epi8(cmp0, qe512,  iqe512);
		
		cmp0 = _mm512_cmpge_epu8_mask(gmax512, endsc512);
		cmp0 &= endsc_msk_a;
		
		__m512i left512 = _mm512_adds_epu8(gmax512, sft512);
		__mmask64 cmp2 = _mm512_cmpge_epu8_mask(left512, cmax512);

		exit0 = (~(cmp0 | cmp2)) & exit0;
		if (exit0 == 0) {
			limit = i++;
			// printf("Breaking, limit: %d\n", limit);
			break;
		}		

		uint8_t *S = H1; H1 = H0; H0 = S;
		i512 = _mm512_add_epi16(i512, one512);
	} // for nrow

	pimax512 = _mm512_mask_blend_epi8(mask512, pimax512, zero512);
	pimax512 = _mm512_mask_blend_epi8(minsc_msk, zero512, pimax512);
	pimax512 = _mm512_mask_blend_epi8(exit0, zero512, pimax512);
	_mm512_store_si512((__m512i *) (rowMax + (i-1) * SIMD_WIDTH8), pimax512);
	
	/******************* DP loop over *****************************/   
    /*************** Partial output setting ***************/
	uint8_t score[SIMD_WIDTH8] __attribute((aligned(64)));

	int16_t te1[SIMD_WIDTH8] __attribute((aligned(64)));
	
	uint8_t qe[SIMD_WIDTH8] __attribute((aligned(64)));
	
	int16_t low[SIMD_WIDTH8] __attribute((aligned(64)));
	int16_t high[SIMD_WIDTH8] __attribute((aligned(64)));
	
	_mm512_store_si512((__m512i *) score, gmax512);	
	_mm512_store_si512((__m512i *) te1, te512);
	_mm512_store_si512((__m512i *) (te1 + SIMD_WIDTH16), te512_);
	_mm512_store_si512((__m512i *) qe, qe512);

	int live = 0;
	for (int l=0; l<SIMD_WIDTH8; l++) {
		int ind = po_ind + l;
		int16_t *te;
		if (i < SIMD_WIDTH16) te = te1;
		else te = te1;
#if !MAINY
		ind = p[l].regid;    // index of corr. aln
		// if (ind != po_ind + l) printf("ind: %d, po_ind: %d, l: %d\n", ind, po_ind, l);
		// assert(ind == po_ind + l);
		if (phase) {
			if (aln[ind].score == score[l]) {
				aln[ind].tb = aln[ind].te - te[l];
				aln[ind].qb = aln[ind].qe - qe[l];
			}
		} else {
			aln[ind].score = score[l] + shift < 255? score[l] : 255;
			aln[ind].te = te[l];
			aln[ind].qe = qe[l];
			if (aln[ind].score != 255) {
				qe[l] = 1;
				live ++;				
			}
			else qe[l] = 0;
		}
#else
		aln[ind].score = score[l] + shift < 255? score[l] : 255;
		aln[ind].te = te[l];
		aln[ind].qe = qe[l];
		if (aln[ind].score != 255) {
			qe[l] = 1;
			live ++;				
		}
		else qe[l] = 0;		
#endif
	}
	
#if !MAINY
	if (phase) return 1;
#endif

	if (live == 0) return 1;

	/*************** Score2 and te2 *******************/
	int qmax = this->g_qmax;
	// assert(qmax == 1);
	// int qmax = 1;
	int maxl = 0 , minh = nrow;
	for (int i=0; i<SIMD_WIDTH8; i++)
	{
		int val = (score[i] + qmax - 1) / qmax;

		int16_t *te;
		if (i < SIMD_WIDTH16) te = te1;
		else te = te1;

		low[i] = te[i] - val;
		high[i] = te[i] + val;
		if (qe[i]) {
			maxl = maxl < low[i] ? low[i] : maxl;
			minh = minh > high[i] ? high[i] : minh;
		}
	}

	max512 = zero512;
	te512 = _mm512_set1_epi16(-1);
	te512_ = _mm512_set1_epi16(-1);
	__m512i low512 = _mm512_load_si512((__m512i*) low);    // make it int16
	__m512i high512 = _mm512_load_si512((__m512i*) high);  // int16
	__m512i low512_ = _mm512_load_si512((__m512i*) (low + SIMD_WIDTH16));    // make it int16
	__m512i high512_ = _mm512_load_si512((__m512i*) (high + SIMD_WIDTH16));  // int16

	
	__m512i rmax512;
	for (int i=0; i< maxl; i++)
	{
		__m512i i512 = _mm512_set1_epi16(i);
		rmax512 = _mm512_load_si512((__m512i*) (rowMax + i*SIMD_WIDTH8));

		__mmask64 mask11 = _mm512_cmpgt_epi16_mask(low512, i512);
		__mmask64 mask12 = _mm512_cmpgt_epi16_mask(low512_, i512);
		__mmask64 mask2 = _mm512_cmpgt_epu8_mask(rmax512, max512);
		__mmask64 mask1 = mask11 | (mask12 << SIMD_WIDTH16);		
		mask2 &= mask1;
		max512 = _mm512_mask_blend_epi8(mask2, max512, rmax512);
		te512  = _mm512_mask_blend_epi16(mask2, te512, i512);
		te512_ = _mm512_mask_blend_epi16(mask2 >> SIMD_WIDTH16, te512_, i512);
	}	
#if 1
	for (int i=minh+1; i<limit; i++)
	{
		__m512i i512 = _mm512_set1_epi16(i);
		rmax512 = _mm512_load_si512((__m512i*) (rowMax + i*SIMD_WIDTH8));
		__mmask64 mask11 = _mm512_cmpgt_epi16_mask(i512, high512);
		__mmask64 mask12 = _mm512_cmpgt_epi16_mask(i512, high512_);
		__mmask64 mask2 = _mm512_cmpgt_epu8_mask(rmax512, max512);
		__mmask64 mask1 = mask11 | (mask12 << SIMD_WIDTH16);
		mask2 &= mask1;
		max512 = _mm512_mask_blend_epi8(mask2, max512, rmax512);
		te512 = _mm512_mask_blend_epi16(mask2, te512, i512);
		te512_ = _mm512_mask_blend_epi16(mask2 >> SIMD_WIDTH16, te512_, i512);	
	}
#endif
	
	int16_t temp4[SIMD_WIDTH8] __attribute((aligned(64)));
	// _mm512_store_si512((__m512i *) temp, max512_);
	_mm512_store_si512((__m512i *) temp, max512);
	_mm512_store_si512((__m512i *) temp4, te512);
	_mm512_store_si512((__m512i *) (temp4 + SIMD_WIDTH16), te512_);		
	for (int i=0; i<SIMD_WIDTH8; i++) {
		int ind = po_ind + i;
		int16_t *te2;
		if (i < SIMD_WIDTH16) te2 = temp4;
		else te2 = temp4;
#if !MAINY
		ind = p[i].regid;    // index of corr. aln
		assert(ind == po_ind + i);
#endif
		if (qe[i]) {
			aln[ind].score2 = (temp[i] == 0? (int)-1: (uint8_t) temp[i]);
			aln[ind].te2 = te2[i];
		} else {
			aln[ind].score2 = -1;
			aln[ind].te2 = -1;
		}
		
#if OUT
		fprintf(stderr, "score: %d, te: %d, qe: %d, score2: %d, te2: %d\n",
				aln[ind].score, aln[ind].te, aln[ind].qe, aln[ind].score2, aln[ind].te2);
#endif
	}

	// printf("Check5..\n");	
#endif // ~__AVX512BW__
	return 1;	
}

/**************** Scalar code *************************/
/**
 * Initialize the query data structure
 *
 * @param size   Number of bytes used to store a score; valid valures are 1 or 2
 * @param qlen   Length of the query sequence
 * @param query  Query sequence
 * @param m      Size of the alphabet
 * @param mat    Scoring matrix in a one-dimension array
 *
 * @return       Query data structure
 */
kswq_t* kswv::ksw_qinit(int size, int qlen, uint8_t *query, int m, const int8_t *mat)
{
	kswq_t *q;
	int slen, a, tmp, p;

	size = size > 1? 2 : 1;
	p = 8 * (3 - size); // # values per __m128i
	slen = (qlen + p - 1) / p; // segmented length
	{ // fillers
		// printf("qlen: %d %d\n", qlen, p*slen);
		for (int j=qlen; j<p*slen; j++)
			query[j] = DUMMY1_;
	}   
	q = (kswq_t*)malloc(sizeof(kswq_t) + 256 + 16 * slen * (m + 4)); // a single block of memory
	q->qp = (__m128i*)(((size_t)q + sizeof(kswq_t) + 15) >> 4 << 4); // align memory
	q->H0 = q->qp + slen * m;
	q->H1 = q->H0 + slen;
	q->E  = q->H1 + slen;
	q->Hmax = q->E + slen;
	q->slen = slen; q->qlen = qlen; q->size = size;
	// compute shift
	tmp = m * m;
	for (a = 0, q->shift = 127, q->mdiff = 0; a < tmp; ++a) { // find the minimum and maximum score
		if (mat[a] < (int8_t)q->shift) q->shift = mat[a];
		if (mat[a] > (int8_t)q->mdiff) q->mdiff = mat[a];
	}
	//for (int  l=0; l<25; l++)
	//	printf("%d ", mat[l]); printf("\n");
	//for (int  l=0; l<25; l++)
	//	printf("%d ", (uint8_t)mat[l]); printf("\n");	
	//printf("Pre %d %d %d\n", q->max, q->shift, q->mdiff);
	
	q->max = q->mdiff;
	q->shift = 256 - q->shift; // NB: q->shift is uint8_t
	q->mdiff += q->shift; // this is the difference between the min and max scores
	// printf("POST %d %d %d\n", q->max, q->shift, q->mdiff);
	// exit(0);
	// An example: p=8, qlen=19, slen=3 and segmentation:
	//  {{0,3,6,9,12,15,18,-1},{1,4,7,10,13,16,-1,-1},{2,5,8,11,14,17,-1,-1}}
	if (size == 1) {
		int8_t *t = (int8_t*)q->qp;
		for (a = 0; a < m; ++a) {
			int i, k, nlen = slen * p;
			const int8_t *ma = mat + a * m;
			for (i = 0; i < slen; ++i)
				for (k = i; k < nlen; k += slen) // p iterations
					*t++ = (k >= qlen? 0 : ma[query[k]]) + q->shift;
		}
	} else {
		int16_t *t = (int16_t*)q->qp;
		for (a = 0; a < m; ++a) {
			int i, k, nlen = slen * p;
			const int8_t *ma = mat + a * m;
			for (i = 0; i < slen; ++i)
				for (k = i; k < nlen; k += slen) // p iterations
					*t++ = (k >= qlen? 0 : ma[query[k]]);
		}
	}
	return q;
}

kswr_t kswv::kswvScalar_u8(kswq_t *q, int tlen, const uint8_t *target,
						  int _o_del, int _e_del, int _o_ins, int _e_ins,
						  int xtra) // the first gap costs -(_o+_e)
{
	int slen, i, m_b, n_b, te = -1, gmax = 0, minsc, endsc;
	uint64_t *b;
	__m128i zero, oe_del, e_del, oe_ins, e_ins, shift, *H0, *H1, *E, *Hmax;
	kswr_t r;

#define __max_16(ret, xx) do { \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 8)); \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 4)); \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 2)); \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 1)); \
    	(ret) = _mm_extract_epi16((xx), 0) & 0x00ff;		\
	} while (0)
	
	// initialization
	r = g_defr;
	minsc = (xtra & KSW_XSUBO)? xtra & 0xffff : 0x10000;
	endsc = (xtra & KSW_XSTOP)? xtra & 0xffff : 0x10000;
	m_b = n_b = 0; b = 0;	zero = _mm_set1_epi32(0);
	oe_del = _mm_set1_epi8(_o_del + _e_del);
	e_del = _mm_set1_epi8(_e_del);
	oe_ins = _mm_set1_epi8(_o_ins + _e_ins);
	e_ins = _mm_set1_epi8(_e_ins);
	shift = _mm_set1_epi8(q->shift);
	H0 = q->H0; H1 = q->H1; E = q->E; Hmax = q->Hmax;
	slen = q->slen;
	for (i = 0; i < slen; ++i) {
		_mm_store_si128(E + i, zero);
		_mm_store_si128(H0 + i, zero);
		_mm_store_si128(Hmax + i, zero);
	}
	// the core loop
	for (i = 0; i < tlen; ++i) {
		int j, k, cmp, imax;
		__m128i e, h, t, f = zero, max = zero, *S = q->qp + target[i] * slen; // s is the 1st score vector
		h = _mm_load_si128(H0 + slen - 1); // h={2,5,8,11,14,17,-1,-1} in the above example
		h = _mm_slli_si128(h, 1); // h=H(i-1,-1); << instead of >> because x64 is little-endian
		for (j = 0; LIKELY(j < slen); ++j) {
			/* SW cells are computed in the following order:
			 *   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
			 *   E(i+1,j) = max{H(i,j)-q, E(i,j)-r}
			 *   F(i,j+1) = max{H(i,j)-q, F(i,j)-r}
			 */
			// compute H'(i,j); note that at the beginning, h=H'(i-1,j-1)
			h = _mm_adds_epu8(h, _mm_load_si128(S + j));
			h = _mm_subs_epu8(h, shift); // h=H'(i-1,j-1)+S(i,j)
			e = _mm_load_si128(E + j); // e=E'(i,j)
			h = _mm_max_epu8(h, e);
			h = _mm_max_epu8(h, f); // h=H'(i,j)
			max = _mm_max_epu8(max, h); // set max
			_mm_store_si128(H1 + j, h); // save to H'(i,j)
			// now compute E'(i+1,j)
			e = _mm_subs_epu8(e, e_del); // e=E'(i,j) - e_del
			t = _mm_subs_epu8(h, oe_del); // h=H'(i,j) - o_del - e_del
			e = _mm_max_epu8(e, t); // e=E'(i+1,j)
			_mm_store_si128(E + j, e); // save to E'(i+1,j)
			// now compute F'(i,j+1)
			f = _mm_subs_epu8(f, e_ins);
			t = _mm_subs_epu8(h, oe_ins); // h=H'(i,j) - o_ins - e_ins
			f = _mm_max_epu8(f, t);
			// get H'(i-1,j) and prepare for the next j
			h = _mm_load_si128(H0 + j); // h=H'(i-1,j)
		}
		// NB: we do not need to set E(i,j) as we disallow adjecent insertion and then deletion
		for (k = 0; LIKELY(k < 16); ++k) { // this block mimics SWPS3; NB: H(i,j) updated in the lazy-F loop cannot exceed max
			f = _mm_slli_si128(f, 1);
			for (j = 0; LIKELY(j < slen); ++j) {
				h = _mm_load_si128(H1 + j);
				h = _mm_max_epu8(h, f); // h=H'(i,j)
				_mm_store_si128(H1 + j, h);
				h = _mm_subs_epu8(h, oe_ins);
				f = _mm_subs_epu8(f, e_ins);
				cmp = _mm_movemask_epi8(_mm_cmpeq_epi8(_mm_subs_epu8(f, h), zero));
				if (UNLIKELY(cmp == 0xffff)) goto end_loop16;
			}
		}
end_loop16:
		//int k;for (k=0;k<16;++k)printf("%d ", ((uint8_t*)&max)[k]);printf("\n");
		__max_16(imax, max); // imax is the maximum number in max
		if (imax >= minsc) { // write the b array; this condition adds branching unfornately
			if (n_b == 0 || (int32_t)b[n_b-1] + 1 != i) { // then append
				if (n_b == m_b) {
					m_b = m_b? m_b<<1 : 8;
					b = (uint64_t*) realloc (b, 8 * m_b);
				}
				b[n_b++] = (uint64_t)imax<<32 | i;
			} else if ((int)(b[n_b-1]>>32) < imax) b[n_b-1] = (uint64_t)imax<<32 | i; // modify the last
		}
		if (imax > gmax) {
			gmax = imax; te = i; // te is the end position on the target
			for (j = 0; LIKELY(j < slen); ++j) // keep the H1 vector
				_mm_store_si128(Hmax + j, _mm_load_si128(H1 + j));
			if (gmax + q->shift >= 255 || gmax >= endsc) break;
		}
		S = H1; H1 = H0; H0 = S; // swap H0 and H1
	}
	r.score = gmax + q->shift < 255? gmax : 255;
	r.te = te;
	if (r.score != 255) { // get a->qe, the end of query match; find the 2nd best score
		int max = -1, tmp, low, high, qlen = slen * 16;
		uint8_t *t = (uint8_t*) Hmax;
		for (i = 0; i < qlen; ++i, ++t)
			if ((int)*t > max) max = *t, r.qe = i / 16 + i % 16 * slen;
			else if ((int)*t == max && (tmp = i / 16 + i % 16 * slen) < r.qe) r.qe = tmp; 
		//printf("%d,%d\n", max, gmax);
		if (b) {
			assert(q->max != 0);
			i = (r.score + q->max - 1) / q->max;
			low = te - i; high = te + i;
			for (i = 0; i < n_b; ++i) {
				int e = (int32_t)b[i];
				if ((e < low || e > high) && (int)(b[i]>>32) > r.score2)
					r.score2 = b[i]>>32, r.te2 = e;
			}
		}
	}
	
#if OUT
	fprintf(stderr, "score: %d, te: %d, qe: %d, score2: %d, te2: %d\n",
			r.score, r.te, r.qe, r.score2, r.te2);
#endif
	
	free(b);
	return r;
}

// -------------------------------------------------------------
// kswc scalar, wrapper function
//-------------------------------------------------------------
void kswv::kswvScalaWrapper(SeqPair *seqPairArray,
							uint8_t *seqBufRef,
							uint8_t *seqBufQer,
							kswr_t* aln,
							int numPairs,
							int nthreads) {

	int8_t mat[25];
	// bwa_fill_scmat(1, 4, mat);
	bwa_fill_scmat(mat);

	int st = 0, ed = numPairs;
	
	for (int i=st; i<ed; i++) {
		SeqPair *p = seqPairArray + i;
		kswr_t *myaln = aln + i;
			
		//uint8_t *seq1 = seqBuf + p->id * 2 * MAX_SEQ_LEN;
        //uint8_t *seq2 = seqBuf + (p->id * 2 + 1) * MAX_SEQ_LEN;
		uint8_t *target = seqBufRef + p->id * MAX_SEQ_LEN_REF_SAM;
		uint8_t *query = seqBufQer + p->id * MAX_SEQ_LEN_QER_SAM;
		int tlen = p->len1;
		int qlen = p->len2;
		int xtra = p->h0;

#if 0
		kswq_t *q = ksw_qinit((xtra & KSW_XBYTE)? 1 : 2, qlen, query, this->m, mat);
		kswr_t ks = kswvScalar_u8_exp(q, tlen, target,
									  o_del, e_del,
									  o_ins, e_ins,
									  xtra);
#else
		kswq_t *q = ksw_qinit(2, qlen, query, this->m, mat);
		kswr_t ks = kswvScalar_i16_exp(q, tlen, target,
									   o_del, e_del,
									   o_ins, e_ins,
									   xtra);
#endif


		myaln->score = ks.score;
		myaln->tb = ks.tb;
		myaln->te = ks.te;
		myaln->qb = ks.qb;
		myaln->qe = ks.qe;
		myaln->score2 = ks.score2;
		myaln->te2 = ks.te2;

		free(q->qp);
		free(q);
	}

}

kswr_t kswv::kswvScalar_u8_exp(kswq_t *q, int tlen, const uint8_t *target,
							   int _o_del, int _e_del, int _o_ins, int _e_ins,
							   int xtra) // the first gap costs -(_o+_e)
{
	int slen, i, m_b, n_b, te = -1, gmax = 0, minsc, endsc;
	uint64_t *b;
	__m128i zero, oe_del, e_del, oe_ins, e_ins, shift, *H0, *H1, *E, *Hmax;
	kswr_t r;
#define SIMD 16
	
#define __max_16(ret, xx) do { \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 8)); \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 4)); \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 2)); \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 1)); \
    	(ret) = _mm_extract_epi16((xx), 0) & 0x00ff;		\
	} while (0)
	
	// initialization
	r = g_defr;
	minsc = (xtra & KSW_XSUBO)? xtra & 0xffff : 0x10000;
	endsc = (xtra & KSW_XSTOP)? xtra & 0xffff : 0x10000;
	assert(minsc < 256);
	// if (endsc > 256)
	//	printf("endsc: %x\n", endsc);
	
	m_b = n_b = 0; b = 0;
	zero = _mm_set1_epi32(0);
	oe_del = _mm_set1_epi8(_o_del + _e_del);
	e_del = _mm_set1_epi8(_e_del);
	oe_ins = _mm_set1_epi8(_o_ins + _e_ins);
	e_ins = _mm_set1_epi8(_e_ins);
	shift = _mm_set1_epi8(q->shift);
	H0 = q->H0; H1 = q->H1; E = q->E; Hmax = q->Hmax;
	slen = q->slen;
	for (i = 0; i < slen; ++i) {
		_mm_store_si128(E + i, zero);
		_mm_store_si128(H0 + i, zero);
		_mm_store_si128(Hmax + i, zero);
	}
	// printf("tlen: %d, slen: %d\n", tlen, slen);
	// the core loop
	for (i = 0; i < tlen; ++i) {
		int j, k, cmp, imax;
		__m128i e, h, t, f = zero, max = zero, *S = q->qp + target[i] * slen; // s is the 1st score vector
		h = _mm_load_si128(H0 + slen - 1); // h={2,5,8,11,14,17,-1,-1} in the above example
		h = _mm_slli_si128(h, 1); // h=H(i-1,-1); << instead of >> because x64 is little-endian
		for (j = 0; LIKELY(j < slen); ++j) {
			/* SW cells are computed in the following order:
			 *   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
			 *   E(i+1,j) = max{H(i,j)-q, E(i,j)-r}
			 *   F(i,j+1) = max{H(i,j)-q, F(i,j)-r}
			 */
			// compute H'(i,j); note that at the beginning, h=H'(i-1,j-1)
			h = _mm_adds_epu8(h, _mm_load_si128(S + j));
			h = _mm_subs_epu8(h, shift); // h=H'(i-1,j-1)+S(i,j)
			e = _mm_load_si128(E + j); // e=E'(i,j)
			h = _mm_max_epu8(h, e);
			h = _mm_max_epu8(h, f); // h=H'(i,j)
			max = _mm_max_epu8(max, h); // set max
			_mm_store_si128(H1 + j, h); // save to H'(i,j)
			// now compute E'(i+1,j)
			e = _mm_subs_epu8(e, e_del); // e=E'(i,j) - e_del
			t = _mm_subs_epu8(h, oe_del); // h=H'(i,j) - o_del - e_del
			e = _mm_max_epu8(e, t); // e=E'(i+1,j)
			_mm_store_si128(E + j, e); // save to E'(i+1,j)
			// now compute F'(i,j+1)
			f = _mm_subs_epu8(f, e_ins);
			t = _mm_subs_epu8(h, oe_ins); // h=H'(i,j) - o_ins - e_ins
			f = _mm_max_epu8(f, t);
			// get H'(i-1,j) and prepare for the next j
			h = _mm_load_si128(H0 + j); // h=H'(i-1,j)
		}
		
		// NB: we do not need to set E(i,j) as we disallow adjecent insertion and then deletion
		for (k = 0; LIKELY(k < 16); ++k) { // this block mimics SWPS3; NB: H(i,j) updated in the lazy-F loop cannot exceed max
			f = _mm_slli_si128(f, 1);
			for (j = 0; LIKELY(j < slen); ++j) {
				h = _mm_load_si128(H1 + j);
				h = _mm_max_epu8(h, f); // h=H'(i,j)
				_mm_store_si128(H1 + j, h);
				h = _mm_subs_epu8(h, oe_ins);
				f = _mm_subs_epu8(f, e_ins);
				cmp = _mm_movemask_epi8(_mm_cmpeq_epi8(_mm_subs_epu8(f, h), zero));
				if (UNLIKELY(cmp == 0xffff)) goto end_loop16;
			}
		}		
	end_loop16:
		// int k;for (k=0;k<16;++k)printf("%d ", ((uint8_t*)&max)[k]);printf("\n");
		__max_16(imax, max); // imax is the maximum number in max

		if (imax >= minsc) { // write the b array; this condition adds branching unfornately
			if (n_b == 0 || (int32_t)b[n_b-1] + 1 != i) { // then append
				if (n_b == m_b) {
					m_b = m_b? m_b<<1 : 8;
					b = (uint64_t*) realloc (b, 8 * m_b);
				}
				b[n_b++] = (uint64_t)imax<<32 | i;
			} else if ((int)(b[n_b-1]>>32) < imax) b[n_b-1] = (uint64_t)imax<<32 | i; // modify the last
		}
		if (imax > gmax) {
			gmax = imax; te = i; // te is the end position on the target
			assert(gmax < 256);
			for (j = 0; LIKELY(j < slen); ++j) // keep the H1 vector
				_mm_store_si128(Hmax + j, _mm_load_si128(H1 + j));
			if (gmax + q->shift >= 255) 
				printf("gmax: %d\n", gmax);
			
			if (gmax + q->shift >= 255 || gmax >= endsc) break;			
		}
		S = H1; H1 = H0; H0 = S; // swap H0 and H1
	}
	r.score = gmax + q->shift < 255? gmax : 255;
	r.te = te;
	// assert(r.score != 255);
	
	
	if (r.score != 255) { // get a->qe, the end of query match; find the 2nd best score
		int max = -1, tmp, low, high, qlen = slen * 16;
		uint8_t *t = (uint8_t*) Hmax;
		for (i = 0; i < qlen; ++i, ++t)
			if ((int)*t > max) max = *t, r.qe = i / 16 + i % 16 * slen;
			else if ((int)*t == max && (tmp = i / 16 + i % 16 * slen) < r.qe) r.qe = tmp; 
		if (b) {
			assert(q->max != 0);
			i = (r.score + q->max - 1) / q->max;
#if 0
			int range = i;
#endif
			low = te - i; high = te + i;
			// printf("n_b: %d, low: %d\n", n_b, low);
			for (i = 0; i < n_b; ++i) {
				int e = (int32_t)b[i];
				// printf("%d %d %d\n", i, e, (int)(b[i]>>32));
				if ((e < low || e > high) && (int)(b[i]>>32) > r.score2)
					r.score2 = b[i]>>32, r.te2 = e;
			}
			{
#if 0			
			int min = -1, t = -1;
			for (i = 0; i < n_b; ++i) {
				int e = (int) (b[i] >> 32);
				if (r.score != e)
				{
					min = e > min? e : min;
					t = (int32_t) b[i];
				}
				printf("i: %d, score: %d\n", (int32_t) b[i], (int) (b[i] >> 32));
			}
			if (r.score2 != min) {
				printf("i (range): %d, low: %d, high: %d, te: %d\n", range, low, high, te);
				printf("score: %d, te: %d, score2: %d, te2: %d, min: %d, i: %d\n",
					   r.score, r.te, r.score2, r.te2, min, t);
			}
			exit(0);
#endif
			}
		}
	}

#if OUT
	fprintf(stderr, "score: %d, te: %d, qe: %d, score2: %d, te2: %d\n",
			r.score, r.te, r.qe, r.score2, r.te2);
#endif

	free(b);
	return r;
}



kswr_t kswv::kswvScalar_i16_exp(kswq_t *q, int tlen, const uint8_t *target,
								int _o_del, int _e_del, int _o_ins, int _e_ins,
								int xtra) // the first gap costs -(_o+_e)
{
	static int itr = 0;
	
	int slen, i, m_b, n_b, te = -1, gmax = 0, minsc, endsc;
	uint64_t *b;
	__m128i zero, oe_del, e_del, oe_ins, e_ins, *H0, *H1, *E, *Hmax;
	kswr_t r;
	
#define SIMD16 8

#define __max_8(ret, xx) do {								 \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 8)); \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 4)); \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 2)); \
    	(ret) = _mm_extract_epi16((xx), 0);					 \
	} while (0)

	// initialization
	r = g_defr;
	minsc = (xtra&KSW_XSUBO)? xtra&0xffff : 0x10000;
	endsc = (xtra&KSW_XSTOP)? xtra&0xffff : 0x10000;
	m_b = n_b = 0; b = 0;
	zero = _mm_set1_epi32(0);
	oe_del = _mm_set1_epi16(_o_del + _e_del);
	e_del = _mm_set1_epi16(_e_del);
	oe_ins = _mm_set1_epi16(_o_ins + _e_ins);
	e_ins = _mm_set1_epi16(_e_ins);
	H0 = q->H0; H1 = q->H1; E = q->E; Hmax = q->Hmax;
	slen = q->slen;
	for (i = 0; i < slen; ++i) {
		_mm_store_si128(E + i, zero);
		_mm_store_si128(H0 + i, zero);
		_mm_store_si128(Hmax + i, zero);
	}


#ifdef VTUNE_ANALYSIS
    __itt_resume();
#endif
	
	// the core loop
	for (i = 0; i < tlen; ++i) {
		int j, k, imax;
		__m128i e, t, h, f = zero, max = zero, *S = q->qp + target[i] * slen; // s is the 1st score vector
		h = _mm_load_si128(H0 + slen - 1); // h={2,5,8,11,14,17,-1,-1} in the above example
		h = _mm_slli_si128(h, 2);

		for (j = 0; LIKELY(j < slen); ++j)
		{
			
			h = _mm_adds_epi16(h, *S++);			
			e = _mm_load_si128(E + j);			
			h = _mm_max_epi16(h, e);
			h = _mm_max_epi16(h, f);
			max = _mm_max_epi16(max, h);
			_mm_store_si128(H1 + j, h);
			e = _mm_subs_epu16(e, e_del);
			t = _mm_subs_epu16(h, oe_del);
			e = _mm_max_epi16(e, t);
			_mm_store_si128(E + j, e);
			f = _mm_subs_epu16(f, e_ins);
			t = _mm_subs_epu16(h, oe_ins);
			f = _mm_max_epi16(f, t);
			h = _mm_load_si128(H0 + j);
			// prof[DP2][0] ++;
			prof[DP2][0] += 14;
		}
				
		for (k = 0; LIKELY(k < 16); ++k) {
			f = _mm_slli_si128(f, 2);
			for (j = 0; LIKELY(j < slen); ++j) {
				h = _mm_load_si128(H1 + j);
				h = _mm_max_epi16(h, f);
				_mm_store_si128(H1 + j, h);
				h = _mm_subs_epu16(h, oe_ins);
				f = _mm_subs_epu16(f, e_ins);
				if(UNLIKELY(!_mm_movemask_epi8(_mm_cmpgt_epi16(f, h)))) goto end_loop8;
			}
		}
end_loop8:
		__max_8(imax, max);
		
		if (imax >= minsc) {
			if (n_b == 0 || (int32_t)b[n_b-1] + 1 != i) {
				if (n_b == m_b) {
					m_b = m_b? m_b<<1 : 8;
					b = (uint64_t*)realloc(b, 8 * m_b);
				}
				b[n_b++] = (uint64_t)imax<<32 | i;
			} else if ((int)(b[n_b-1]>>32) < imax) b[n_b-1] = (uint64_t)imax<<32 | i; // modify the last
		}
		if (imax > gmax) {
			gmax = imax; te = i;
			for (j = 0; LIKELY(j < slen); ++j)
				_mm_store_si128(Hmax + j, _mm_load_si128(H1 + j));
			if (gmax >= endsc) {
				// printf("Breaking, gmax: %d, endsc: %d\n", gmax, endsc);
				break;
			}
		}
		S = H1; H1 = H0; H0 = S;
	}


#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif

	int max_ = -1;
	r.score = gmax; r.te = te;
	int max = -1, tmp, low, high, qlen = slen * 8;
	{
		uint16_t *t = (uint16_t*)Hmax;
		for (i = 0, r.qe = -1; i < qlen; ++i, ++t)
			if ((int)*t > max) max = *t, r.qe = i / 8 + i % 8 * slen;
			else if ((int)*t == max && (tmp = i / 8 + i % 8 * slen) < r.qe) r.qe = tmp;
		max_ = max;
		if (b) {
			assert(q->max != 0);
			i = (r.score + q->max - 1) / q->max;
			low = te - i; high = te + i;
			for (i = 0; i < n_b; ++i) {
				int e = (int32_t)b[i];
				if ((e < low || e > high) && (int)(b[i]>>32) > r.score2)
					r.score2 = b[i]>>32, r.te2 = e;
			}
		}
	}

	// fprintf(stderr, "gmax: %d, te: %d, qe: %d, max: %d, max2: %d, te2: %d\n",
	// gmax, te, r.qe, max_, r.score2, r.te2);
#if OUT
	fprintf(stderr, "score: %d, te: %d, qe: %d, score2: %d, te2: %d\n",
			gmax, te, r.qe, r.score2, r.te2);
#endif

	itr ++;	
	free(b);
	return r;
}

#if 1   //original
kswr_t kswv::kswvScalar_i16(kswq_t *q, int tlen, const uint8_t *target,
							int _o_del, int _e_del, int _o_ins, int _e_ins,
							int xtra) // the first gap costs -(_o+_e)
{
	int slen, i, m_b, n_b, te = -1, gmax = 0, minsc, endsc;
	uint64_t *b;
	__m128i zero, oe_del, e_del, oe_ins, e_ins, *H0, *H1, *E, *Hmax;
	kswr_t r;
#define SIMD16 8

#define __max_8(ret, xx) do { \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 8)); \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 4)); \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 2)); \
    	(ret) = _mm_extract_epi16((xx), 0); \
	} while (0)

	// initialization
	r = g_defr;
	minsc = (xtra&KSW_XSUBO)? xtra&0xffff : 0x10000;
	endsc = (xtra&KSW_XSTOP)? xtra&0xffff : 0x10000;
	m_b = n_b = 0; b = 0;
	zero = _mm_set1_epi32(0);
	oe_del = _mm_set1_epi16(_o_del + _e_del);
	e_del = _mm_set1_epi16(_e_del);
	oe_ins = _mm_set1_epi16(_o_ins + _e_ins);
	e_ins = _mm_set1_epi16(_e_ins);
	H0 = q->H0; H1 = q->H1; E = q->E; Hmax = q->Hmax;
	slen = q->slen;
	for (i = 0; i < slen; ++i) {
		_mm_store_si128(E + i, zero);
		_mm_store_si128(H0 + i, zero);
		_mm_store_si128(Hmax + i, zero);
	}
	// the core loop
	for (i = 0; i < tlen; ++i) {
		int j, k, imax;
		__m128i e, t, h, f = zero, max = zero, *S = q->qp + target[i] * slen; // s is the 1st score vector
		h = _mm_load_si128(H0 + slen - 1); // h={2,5,8,11,14,17,-1,-1} in the above example
		h = _mm_slli_si128(h, 2);
		for (j = 0; LIKELY(j < slen); ++j) {

			h = _mm_adds_epi16(h, *S++);
			e = _mm_load_si128(E + j);
			h = _mm_max_epi16(h, e);
			h = _mm_max_epi16(h, f);
			max = _mm_max_epi16(max, h);
			_mm_store_si128(H1 + j, h);
			e = _mm_subs_epu16(e, e_del);
			t = _mm_subs_epu16(h, oe_del);
			e = _mm_max_epi16(e, t);
			_mm_store_si128(E + j, e);
			f = _mm_subs_epu16(f, e_ins);
			t = _mm_subs_epu16(h, oe_ins);
			f = _mm_max_epi16(f, t);
			h = _mm_load_si128(H0 + j);
		}
		
		for (k = 0; LIKELY(k < 16); ++k) {
			f = _mm_slli_si128(f, 2);
			for (j = 0; LIKELY(j < slen); ++j) {
				h = _mm_load_si128(H1 + j);
				h = _mm_max_epi16(h, f);
				_mm_store_si128(H1 + j, h);
				h = _mm_subs_epu16(h, oe_ins);
				f = _mm_subs_epu16(f, e_ins);
				if(UNLIKELY(!_mm_movemask_epi8(_mm_cmpgt_epi16(f, h)))) goto end_loop8;
			}
		}
end_loop8:
		__max_8(imax, max);
		if (imax >= minsc) {
			if (n_b == 0 || (int32_t)b[n_b-1] + 1 != i) {
				if (n_b == m_b) {
					m_b = m_b? m_b<<1 : 8;
					b = (uint64_t*)realloc(b, 8 * m_b);
				}
				b[n_b++] = (uint64_t)imax<<32 | i;
			} else if ((int)(b[n_b-1]>>32) < imax) b[n_b-1] = (uint64_t)imax<<32 | i; // modify the last
		}
		if (imax > gmax) {
			gmax = imax; te = i;
			for (j = 0; LIKELY(j < slen); ++j)
				_mm_store_si128(Hmax + j, _mm_load_si128(H1 + j));
			if (gmax >= endsc) break;
		}
		S = H1; H1 = H0; H0 = S;
	}
	
	r.score = gmax; r.te = te;
	{
		int max = -1, tmp, low, high, qlen = slen * 8;
		uint16_t *t = (uint16_t*)Hmax;
		for (i = 0, r.qe = -1; i < qlen; ++i, ++t)
			if ((int)*t > max) max = *t, r.qe = i / 8 + i % 8 * slen;
			else if ((int)*t == max && (tmp = i / 8 + i % 8 * slen) < r.qe) r.qe = tmp; 
		if (b) {
			assert(q->max != 0);
			i = (r.score + q->max - 1) / q->max;
			low = te - i; high = te + i;
			for (i = 0; i < n_b; ++i) {
				int e = (int32_t)b[i];
				if ((e < low || e > high) && (int)(b[i]>>32) > r.score2)
					r.score2 = b[i]>>32, r.te2 = e;
			}
		}
	}

	free(b);
	return r;
}
#endif

#if 0
kswr_t kswv::kswvScalarPure(kswq_t *q, int tlen, const uint8_t *target,
							int _o_del, int _e_del, int _o_ins, int _e_ins,
							int xtra) // the first gap costs -(_o+_e)
{
	static int itr = 0;
	int lane_ = 0;
	bool gate = 0;
	
	int slen, i, m_b, n_b, te = -1, gmax = 0, minsc, endsc;
	uint64_t *b;
	int zero, oe_del, e_del, oe_ins, e_ins, *H0, *H1, *E, *Hmax;
	kswr_t r;
	
#define SIMD16 8
	int16_t dp[1048][200], rMax[tlen+1];
	int16_t temp[SIMD16];


	// initialization
	r = g_defr;
	minsc = (xtra&KSW_XSUBO)? xtra&0xffff : 0x10000;
	endsc = (xtra&KSW_XSTOP)? xtra&0xffff : 0x10000;
	m_b = n_b = 0; b = 0;
	zero = 0;
	oe_del = _o_del + _e_del;
	e_del = _e_del;
	oe_ins = o_ins + _e_ins;
	e_ins = _e_ins;
	H0 = q->H0; H1 = q->H1; E = q->E; Hmax = q->Hmax;  // fix
	slen = q->slen;

	for (i = 0; i < slen*SIMD16; ++i) {
		E[i] = H0[i] = Hmax[i] = 0;
	}

	// the core loop
	for (i = 0; i < tlen; ++i) {
		int j, k, imax;
		int e, t, h, f = zero, max = zero, *S = q->qp + target[i] * slen; // s is the 1st score vector

		h = H0[0];
		for (j = 0; j < slen*SIMD16; ++j)
		{

			h += *S++;
			e = E[j+1];
			h = max_(h, e);
			h = max_(h, f);
			H1[j+1] = h;
			max = max_(max, h);
			e -= e_del;
			t = h - oe_del;
			e = max_(e, t);
			E[j+1] = e;
			f -= e_ins;
			t = h - oe_ins;
			f = max_(f, t);
			h = H0[j+1];
			prof[DP2][0] ++;
		}	
		imax = max;
		if (imax >= minsc) {
			if (n_b == 0 || (int32_t)b[n_b-1] + 1 != i) {
				if (n_b == m_b) {
					m_b = m_b? m_b<<1 : 8;
					b = (uint64_t*)realloc(b, 8 * m_b);
				}
				b[n_b++] = (uint64_t)imax<<32 | i;
			} else if ((int)(b[n_b-1]>>32) < imax) b[n_b-1] = (uint64_t)imax<<32 | i; // modify the last
		}
		if (imax > gmax) {
			gmax = imax; te = i;
			for (j = 0; LIKELY(j < slen*SIMD16); ++j)
				//_mm_store_si128(Hmax + j, _mm_load_si128(H1 + j));
				Hmax[j] = H1[j];
			if (gmax >= endsc) break;
		}
		S = H1; H1 = H0; H0 = S;		
	}
	r.score = gmax; r.te = te;
	int max = -1, tmp, low, high, qlen = slen * 8;
	{
		int *t = (int*) Hmax;
		for (i = 0, r.qe = -1; i < qlen; ++i, ++t)
			if ((int)*t > max) max = *t, r.qe = i;
			else if ((int)*t == max && (tmp = i) < r.qe) r.qe = tmp;
		if (b) {
			assert(q->max != 0);
			i = (r.score + q->max - 1) / q->max;
			low = te - i; high = te + i;
			for (i = 0; i < n_b; ++i) {
				int e = (int32_t)b[i];
				if ((e < low || e > high) && (int)(b[i]>>32) > r.score2)
					r.score2 = b[i]>>32, r.te2 = e;
			}
		}
	}
	
}
#endif
/******************************************************************************/
// main function:
#if MAINY
#define DEFAULT_MATCH 1
#define DEFAULT_MISMATCH -4
#define DEFAULT_OPEN 6
#define DEFAULT_EXTEND 1
#define DEFAULT_AMBIG -1

//#define MAX_NUM_PAIRS 1000
//#define MATRIX_MIN_CUTOFF -100000000
//#define LOW_INIT_VALUE (INT32_MIN/2)
// #define AMBIG 52

double freq = 2.3*1e9;
int32_t w_match, w_mismatch, w_open, w_extend, w_ambig;
uint64_t SW_cells;
char *pairFileName;
FILE *pairFile;
int8_t h0 = 0;
double clock_freq;
//uint64_t prof[10][112], data, SW_cells2;


void parseCmdLine(int argc, char *argv[])
{
    int i;
    w_match = DEFAULT_MATCH;
    w_mismatch = DEFAULT_MISMATCH;
    w_open = DEFAULT_OPEN;
    w_extend = DEFAULT_EXTEND;
	w_ambig = DEFAULT_AMBIG;
	
    int pairFlag = 0;
    for(i = 1; i < argc; i+=2)
    {
        if(strcmp(argv[i], "-match") == 0)
        {
            w_match = atoi(argv[i + 1]);
        }
        if(strcmp(argv[i], "-mismatch") == 0) //penalty, +ve number
        {
            w_mismatch = atoi(argv[i + 1]);
        }
		if(strcmp(argv[i], "-ambig") == 0)
        {
            w_ambig = atoi(argv[i + 1]);
        }

        if(strcmp(argv[i], "-gapo") == 0)
        {
            w_open = atoi(argv[i + 1]);
        }
        if(strcmp(argv[i], "-gape") == 0)
        {
            w_extend = atoi(argv[i + 1]);
        }
        if(strcmp(argv[i], "-pairs") == 0)
        {
            pairFileName = argv[i + 1];
            pairFlag = 1;
        }
        if(strcmp(argv[i], "-h0") == 0)
        {
            h0 = atoi(argv[i + 1]);
        }
    }
    if(pairFlag == 0)
    {
        printf("ERROR! pairFileName not specified.\n");
        exit(0);
    }
}

int loadPairs(SeqPair *seqPairArray, uint8_t *seqBufRef, uint8_t* seqBufQer, FILE *pairFile)
{
	static int32_t cnt = 0;
    int32_t numPairs = 0;
    while(numPairs < MAX_NUM_PAIRS_ALLOC)
    {
		int32_t xtra = 0;
		char temp[10];
		fgets(temp, 10, pairFile);
		sscanf(temp, "%x", &xtra);
		//printf("xtra: %d, %x, %s\n", xtra, xtra, temp);

        //if(!fgets((char *)(seqBuf + numPairs * 2 * MAX_SEQ_LEN), MAX_SEQ_LEN, pairFile))
		if(!fgets((char *)(seqBufRef + numPairs * MAX_SEQ_LEN_REF_SAM), MAX_SEQ_LEN_REF_SAM, pairFile))
        {
            break;
        }
        //if(!fgets((char *)(seqBuf + (numPairs * 2 + 1) * MAX_SEQ_LEN), MAX_SEQ_LEN, pairFile))
		if(!fgets((char *)(seqBufQer + numPairs * MAX_SEQ_LEN_QER_SAM), MAX_SEQ_LEN_QER_SAM, pairFile))	
        {
            printf("ERROR! Odd number of sequences in %s\n", pairFileName);
            break;
        }

        SeqPair sp;
        sp.id = numPairs;
        // sp.seq1 = seqBuf + numPairs * 2 * MAX_SEQ_LEN;
        // sp.seq2 = seqBuf + (numPairs * 2 + 1) * MAX_SEQ_LEN;
        sp.len1 = strnlen((char *)(seqBufRef + numPairs * MAX_SEQ_LEN_REF_SAM), MAX_SEQ_LEN_REF_SAM) - 1;
        sp.len2 = strnlen((char *)(seqBufQer + numPairs * MAX_SEQ_LEN_QER_SAM), MAX_SEQ_LEN_QER_SAM) - 1;
		sp.h0 = xtra;
        // sp.score = 0;

		//uint8_t *seq1 = seqBuf + numPairs * 2 * MAX_SEQ_LEN;
        //uint8_t *seq2 = seqBuf + (numPairs * 2 + 1) * MAX_SEQ_LEN;
		uint8_t *seq1 = seqBufRef + numPairs * MAX_SEQ_LEN_REF_SAM;
        uint8_t *seq2 = seqBufQer + numPairs * MAX_SEQ_LEN_QER_SAM;

		for (int l=0; l<sp.len1; l++)
			seq1[l] -= 48;
		for (int l=0; l<sp.len2; l++)
			seq2[l] -= 48;

        seqPairArray[numPairs] = sp;
        numPairs++;
        // SW_cells += (sp.len1 * sp.len2);
    }
    // fclose(pairFile);
    return numPairs;
}

// profiling stats
uint64_t find_stats(uint64_t *val, int nt, double &min, double &max, double &avg) {
	min = 1e10;
	max = 0;
	avg = 0;
	for (int i=0; i<nt; i++) {
		avg += val[i];
		if (max < val[i]) max = val[i];
		if (min > val[i]) min = val[i];		
	}
	avg /= nt;
}

int main(int argc, char *argv[]) {

#ifdef VTUNE_ANALYSIS
	printf("Vtune analysis enabled....\n");
    __itt_pause();
#endif
	fsam = fopen("fsam.txt", "w");
	
	kswv *mate;

    parseCmdLine(argc, argv);

	SeqPair *seqPairArray = (SeqPair *)_mm_malloc((MAX_NUM_PAIRS + SIMD_WIDTH8) * sizeof(SeqPair), 64);
	uint8_t *seqBufRef = NULL, *seqBufQer = NULL;
	seqBufRef = (uint8_t *)_mm_malloc((MAX_SEQ_LEN_REF_SAM * MAX_NUM_PAIRS + MAX_LINE_LEN)
									  * sizeof(int8_t), 64);
	seqBufQer = (uint8_t *)_mm_malloc((MAX_SEQ_LEN_QER_SAM * MAX_NUM_PAIRS + MAX_LINE_LEN)
									  * sizeof(int8_t), 64);

	kswr_t *aln = NULL;
	aln = (kswr_t *) _mm_malloc ((MAX_NUM_PAIRS + SIMD_WIDTH8) * sizeof(kswr_t), 64);
	
	if (seqBufRef == NULL || seqBufQer == NULL || aln == NULL)
	{
		printf("Memory not allocated\nExiting...\n");
		exit(0);
	} else {
		printf("Memory allocated: %0.2lf MB\n",
			   ((int64_t)(MAX_SEQ_LEN_REF_SAM + MAX_SEQ_LEN_QER_SAM ) * MAX_NUM_PAIRS + MAX_LINE_LEN)/1e6);
	}
	uint64_t tim = __rdtsc(), readTim = 0;
	
    int32_t numThreads = 1;
#pragma omp parallel
    {
        int32_t tid = omp_get_thread_num();
        int32_t nt = omp_get_num_threads();
        if(tid == (nt - 1))
        {
            numThreads = nt;
        }
    }
	numThreads =1 ;
	//printf("Done reading input file!!, numPairs: %d, nt: %d\n",
	//	   numPairs, numThreads);

#if SORT_PAIRS     // disbaled in bwa-mem2 (only used in separate benchmark sw code)
	printf("\tSorting is enabled !!\n");
#endif

#if OUT
	printf("\tResults printing (on console) is enabled.....\n");
#endif

	tim = __rdtsc();
	sleep(1);
	freq = __rdtsc() - tim;

	int numPairs = 0, totNumPairs = 0;
    FILE *pairFile = fopen(pairFileName, "r");
    if(pairFile == NULL)
    {
        fprintf(stderr, "Could not open file: %s\n", pairFileName);
        exit(0);
    }

	// BandedPairWiseSW *pwsw = new BandedPairWiseSW(w_match, w_mismatch, w_open,
	//											  w_extend, w_ambig, end_bonus);
	kswv *pwsw = new kswv(w_open, w_extend, w_open, w_extend, w_match, w_mismatch, numThreads);
	

	int64_t myTicks = 0;
	printf("Processor freq: %0.2lf MHz\n", freq/1e6);
	printf("Executing int8 code....\n");
#if VEC
	printf("Executing AVX512 vectorized code!!\n");
	
	while (1) {
		uint64_t tim = __rdtsc();

		//printf("Loading current batch of pairs..\n");
		numPairs = loadPairs(seqPairArray, seqBufRef, seqBufQer, pairFile);
		totNumPairs += numPairs;
#if STAT
		printf("Loading done, numPairs: %d\n", numPairs);	   
		if (totNumPairs > spot) spot -= totNumPairs - numPairs;
		else continue;
#endif
		readTim += __rdtsc() - tim;
		if (numPairs == 0) break;

		tim = __rdtsc();
		int phase = 0;
		pwsw->getScores8(seqPairArray, seqBufRef, seqBufQer, aln, numPairs, numThreads, phase);
		// pwsw->getScores16(seqPairArray, seqBufRef, seqBufQer, aln, numPairs, numThreads, phase);
		myTicks += __rdtsc() - tim;

#if STAT
		for (int l=0; l<10; l++)
		{
			// SeqPair r = seqPairArray[l];
			kswr_t r = aln[l];
			fprintf(stderr, "%d %d %d %d %d %d %d\n", r.score, r.tb, r.te, r.qb, r.qe, r.score2, r.te2);
		}		
		break;
#endif
		
#if OUT
		// printf("Execution complete!!, writing output!\n");
		//for (int l=0; l<numPairs; l++)
		//{
		//	// SeqPair r = seqPairArray[l];
		//	kswr_t r = aln[l];
		//	fprintf(stderr, "%d %d %d %d %d %d %d\n", r.score, r.tb, r.te, r.qb, r.qe, r.score2, r.te2);
		//}
		// printf("Vector code: Writing output completed!!!\n\n");
		// printf("Vector code -- Wrote output to the file\n");
#endif
	}
	
	// int64_t myTicks = 0;//pwsw->getTicks();
	printf("Read time  = %0.2lf\n", readTim/freq);
    printf("Overall SW cycles = %ld, %0.2lf\n", myTicks, myTicks*1.0/freq);
    printf("SW cells(T)  = %ld\n", SW_cells);
	printf("SW cells(||)  = %ld\n", SW_cells2);
    printf("SW GCUPS  = %lf\n", SW_cells * freq / myTicks);
	
	{
		printf("More stats:\n");
		// double freq = 2.3*1e9;
		double min, max, avg;
		find_stats(prof[1], numThreads, min, max, avg);
		printf("Time in pre-processing: %0.2lf (%0.2lf, %0.2lf)\n",
			   avg*1.0/freq, min*1.0/freq, max*1.0/freq);
		find_stats(prof[0], numThreads, min, max, avg);
		printf("Time spent in smithWaterman(): %0.2lf (%0.2lf, %0.2lf)\n",
			   avg*1.0/freq, min*1.0/freq, max*1.0/freq);
		printf("\nTotal cells computed: %ld, %0.2lf\n",
			   prof[2][0], prof[2][0]*1.0/(totNumPairs));
		printf("\tTotal useful cells computed: %ld, %0.2lf\n",
			   prof[3][0], prof[3][0]*1.0/(totNumPairs));
		// printf("Total bp read from memory: %ld\n", data);
		printf("Computations after exit: %ld\n", prof[4][0]);
		printf("Cumulative seq1/2 len: %ld, %ld\n", prof[4][1], prof[4][2]);
	
#if 1
		printf("\nDebugging info:\n");
		printf("Time taken for DP loop: %0.2lf\n", prof[DP][0]*1.0/freq);
		printf("Time taken for DP loop upper part: %0.2lf\n", prof[DP3][0]*1.0/freq);	
		printf("Time taken for DP inner loop: %0.2lf\n", prof[DP1][0]*1.0/freq);
		printf("Time taken for DP loop lower part: %ld\n", prof[DP2][0]);
#endif
	}

	
#else
	printf("Executing scalar code!!\n");
	while (1) {
		uint64_t tim = __rdtsc();
		numPairs = loadPairs(seqPairArray, seqBufRef, seqBufQer, pairFile);
		totNumPairs += numPairs;
#if STAT
		printf("Loading done, numPairs: %d\n", numPairs);
		if (totNumPairs > spot) spot -= totNumPairs - numPairs;
		else continue;
#endif
		readTim += __rdtsc() - tim;
		if (numPairs == 0) break;

		tim = __rdtsc();
		pwsw->kswvScalaWrapper(seqPairArray,
							   seqBufRef,
							   seqBufQer,
							   aln,
							   numPairs,
							   numThreads);
		myTicks += __rdtsc() - tim;
#if STAT
		for (int l=0; l<10; l++)
		{
			// SeqPair r = seqPairArray[l];
			kswr_t r = aln[l];
			fprintf(stderr, "%d %d %d %d %d %d %d\n", r.score, r.tb, r.te, r.qb, r.qe, r.score2, r.te2);
		}		
		break;
#endif
		
	} // while
	
	printf("Read time  = %0.2lf\n", readTim/freq);
	printf("Overall SW cycles = %ld, %0.2lf\n", myTicks, myTicks*1.0/freq);
	printf("Time taken for DP loop lower part: %ld \n", prof[DP2][0]);

#if 1
		printf("\nDebugging info:\n");
		printf("Time taken for DP loop: %0.2lf\n", prof[DP][0]*1.0/freq);
		printf("Time taken for DP loop upper part: %0.2lf\n", prof[DP3][0]*1.0/freq);	
		printf("Time taken for DP inner loop: %0.2lf\n", prof[DP1][0]*1.0/freq);
		printf("Time taken for DP loop lower part: %ld\n", prof[DP2][0]);
#endif

#endif


#ifdef VTUNE_ANALYSIS
	printf("Vtune analysis enabled....\n");
#endif
		
	// free memory
	_mm_free(seqPairArray);
	_mm_free(seqBufRef);
	_mm_free(seqBufQer);
	_mm_free(aln);
	
	fclose(pairFile);
	fclose(fsam);
	return 1;
}
#endif

/// 16 bit lanes
void kswv::getScores16(SeqPair *pairArray,
					   uint8_t *seqBufRef,
					   uint8_t *seqBufQer,
					   kswr_t* aln,
					   int32_t numPairs,
					   uint16_t numThreads,
					   int phase)
{
    int64_t startTick, endTick;

    // startTick = __rdtsc();
	kswvBatchWrapper16(pairArray, seqBufRef, seqBufQer, aln, numPairs, numThreads, phase);
	// kswvBatchWrapper16_intra(pairArray, seqBufRef, seqBufQer, numPairs, numThreads);
    // endTick = __rdtsc();

}

// #define PFD 2
void kswv::kswvBatchWrapper16(SeqPair *pairArray,
							  uint8_t *seqBufRef,
							  uint8_t *seqBufQer,
							  kswr_t* aln,
							  int32_t numPairs,							  
							  uint16_t numThreads,
							  int phase)
{
	// printf("numThreads: %d %d\n", numThreads, omp_get_thread_num());

	int64_t st1, st2, st3, st4, st5;
    // st1 = __rdtsc();
    int16_t *seq1SoA = NULL;
	seq1SoA = (int16_t *)_mm_malloc(MAX_SEQ_LEN_REF_SAM * SIMD_WIDTH16 * numThreads * sizeof(int16_t), 64);
	
    int16_t *seq2SoA = NULL;
	seq2SoA = (int16_t *)_mm_malloc(MAX_SEQ_LEN_QER_SAM * SIMD_WIDTH16 * numThreads * sizeof(int16_t), 64);

    assert(seq1SoA != NULL);
    assert(seq2SoA != NULL);	
	
	// int16_t	*qp	= qp16 + tid * SIMD_WIDTH16 * MAX_SEQ_LEN_QER_SAM * this->m;
	
    int32_t ii;
    int32_t roundNumPairs = ((numPairs + SIMD_WIDTH16 - 1) / SIMD_WIDTH16 ) * SIMD_WIDTH16;
	// assert(roundNumPairs < BATCH_SIZE * SEEDS_PER_READ);
    for(ii = numPairs; ii < roundNumPairs; ii++)
    {
		pairArray[ii].regid = ii;
        pairArray[ii].id = ii;
        pairArray[ii].len1 = 0;
        pairArray[ii].len2 = 0;
    }
		
    // st2 = __rdtsc();	
#if SORT_PAIRS     // disbaled in bwa-mem2 (only used in separate benchmark sw code)
	{
    // Sort the sequences according to decreasing order of lengths
    SeqPair *tempArray = (SeqPair *)_mm_malloc(SORT_BLOCK_SIZE * numThreads *
											   sizeof(SeqPair), 64);
    int16_t *hist = (int16_t *)_mm_malloc((MAX_SEQ_LEN_QER_SAM + 32) * numThreads *
										  sizeof(int16_t), 64);

#pragma omp parallel num_threads(numThreads)
    {
        int32_t tid = omp_get_thread_num();
        SeqPair *myTempArray = tempArray + tid * SORT_BLOCK_SIZE;
        int16_t *myHist = hist + tid * (MAX_SEQ_LEN_QER_SAM + 32);

#pragma omp for
        for(ii = 0; ii < roundNumPairs; ii+=SORT_BLOCK_SIZE)
        {
            int32_t first, last;
            first = ii;
            last  = ii + SORT_BLOCK_SIZE;
            if(last > roundNumPairs) last = roundNumPairs;
            sortPairsLen(pairArray + first, last - first, myTempArray, myHist);
        }
    }
    _mm_free(hist);
	}
#endif
	
    // st3 = __rdtsc();

//#pragma omp parallel num_threads(numThreads)
    {
        int32_t i;
        // uint16_t tid = omp_get_thread_num();
		uint16_t tid = 0;
        int16_t *mySeq1SoA = NULL;
		mySeq1SoA = seq1SoA + tid * MAX_SEQ_LEN_REF_SAM * SIMD_WIDTH16;
		assert(mySeq1SoA != NULL);
			
		int16_t *mySeq2SoA = NULL;
		mySeq2SoA = seq2SoA + tid * MAX_SEQ_LEN_QER_SAM * SIMD_WIDTH16;
		assert(mySeq1SoA != NULL);
		
        uint8_t *seq1;
        uint8_t *seq2;
		
		int nstart = 0, nend = numPairs;


//#pragma omp for schedule(dynamic, 128)
		for(i = nstart; i < nend; i+=SIMD_WIDTH16)
		{
            int32_t j, k;
            int maxLen1 = 0;
            int maxLen2 = 0;

            for(j = 0; j < SIMD_WIDTH16; j++)
            {
                SeqPair sp = pairArray[i + j];
#if MAINY
				seq1 = seqBufRef + (int64_t)sp.id * MAX_SEQ_LEN_REF_SAM;
#else
				seq1 = seqBufRef + sp.idr;
#endif
                for(k = 0; k < sp.len1; k++)
                {
                    // mySeq1SoA[k * SIMD_WIDTH16 + j] = (seq1[k] == AMBIG?0xFFFF:seq1[k]);
					mySeq1SoA[k * SIMD_WIDTH16 + j] = (seq1[k] == AMBIG_ ? AMBR16:seq1[k]);
                }
                if(maxLen1 < sp.len1) maxLen1 = sp.len1;
            }
            for(j = 0; j < SIMD_WIDTH16; j++)
            {
                SeqPair sp = pairArray[i + j];
                for(k = sp.len1; k <= maxLen1; k++) //removed "="
                {
					mySeq1SoA[k * SIMD_WIDTH16 + j] = DUMMY1_;
                }
            }

			// int16_t offset1 = (-1 * this->w_mismatch + this->w_ambig) << 4;
			// int16_t offset2 = (-1 * this->w_mismatch) << 4;
            for(j = 0; j < SIMD_WIDTH16; j++)
            {
                SeqPair sp = pairArray[i + j];
#if MAINY
				seq2 = seqBufQer + (int64_t)sp.id * MAX_SEQ_LEN_QER_SAM;
#else
				seq2 = seqBufQer + sp.idq;
#endif
				assert(sp.len2 < MAX_SEQ_LEN_QER_SAM);

#if 1
				int quanta = 8 - sp.len2 % 8;  // based on SSE2-16 bit lane
#else
				int quanta = 16 - sp.len2 % 16;  // based on SSE2-8 bit lane
#endif
				// printf("quanta: %d %d %d\n", quanta, sp.len2, sp.len2%8);
                for(k = 0; k < sp.len2; k++)
                {
					// mySeq2SoA[k * SIMD_WIDTH16 + j] = (seq2[k]==AMBIG?0xFFFF:seq2[k]);
					mySeq2SoA[k * SIMD_WIDTH16 + j] = (seq2[k]==AMBIG_? AMBQ16:seq2[k]);
                }
				for(k = sp.len2; k < sp.len2 + quanta; k++) {
					mySeq2SoA[k * SIMD_WIDTH16 + j] = DUMMY3;
				}
                if(maxLen2 < (sp.len2 + quanta)) maxLen2 = sp.len2 + quanta;
            }
			
            for(j = 0; j < SIMD_WIDTH16; j++)
            {
                SeqPair sp = pairArray[i + j];
#if 1
				int quanta = 8 - sp.len2 % 8;  // based on SSE2-16 bit lane
#else
				int quanta = 16 - sp.len2 % 16;  // based on SSE2-8 bit lane
#endif
                for(k = sp.len2 + quanta; k <= maxLen2; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH16 + j] = DUMMY2_;
                }
            }

            kswv512_16_exp(mySeq1SoA, mySeq2SoA,
						   maxLen1, maxLen2,
						   pairArray + i,
						   aln, i,
						   tid,
						   numPairs,
						   phase);
        }
    }
	
    // st4 = __rdtsc();
	
#if SORT_PAIRS     // disbaled in bwa-mem2 (only used in separate benchmark sw code)
	{
    // Sort the sequences according to increasing order of id
#pragma omp parallel num_threads(numThreads)
    {
        int32_t tid = omp_get_thread_num();
        SeqPair *myTempArray = tempArray + tid * SORT_BLOCK_SIZE;

#pragma omp for
        for(ii = 0; ii < roundNumPairs; ii+=SORT_BLOCK_SIZE)
        {
            int32_t first, last;
            first = ii;
            last  = ii + SORT_BLOCK_SIZE;
            if(last > roundNumPairs) last = roundNumPairs;
            sortPairsId(pairArray + first, first, last - first, myTempArray);
        }
    }
    _mm_free(tempArray);
	}
#endif
	
    // st5 = __rdtsc();
    setupTicks = st2 - st1;
    sort1Ticks = st3 - st2;
    swTicks = st4 - st3;
    sort2Ticks = st5 - st4;

	// free mem
	_mm_free(seq1SoA);
	_mm_free(seq2SoA);
    return;
}

/********************** Inter-Task Execution ****************************/
int kswv::kswv512_16_exp(int16_t seq1SoA[],
						 int16_t seq2SoA[],
						 int16_t nrow,
						 int16_t ncol,
						 SeqPair *p,
						 kswr_t *aln,
						 int po_ind,
						 uint16_t tid,
						 int32_t numPairs,
						 int phase)
{
#ifdef __AVX512BW__
	int m_b, n_b;
	// int16_t minsc[SIMD_WIDTH16], endsc[SIMD_WIDTH16];
	int16_t minsc[SIMD_WIDTH16] = {0}, endsc[SIMD_WIDTH16] = {0};
	uint64_t *b;
	int limit = nrow;
	
	__m512i zero512 = _mm512_setzero_si512();
	__m512i one512  = _mm512_set1_epi16(1);
	__m512i minus1  = _mm512_set1_epi16(-1);

	//int16_t temp[SIMD_WIDTH16] __attribute((aligned(64))) = {0};
	int16_t temp[SIMD_WIDTH16] __attribute((aligned(64))) = {0};
	int16_t temp1[SIMD_WIDTH16] __attribute((aligned(64)));
	int16_t temp2[SIMD_WIDTH16] __attribute((aligned(64)));

	//temp[0] = this->w_match;
	//temp[1] = temp[2] = temp[3] =  this->w_mismatch;
	//temp[4] = temp[5] = temp[6] = temp[7] =  this->w_mismatch;
	//temp[16] = temp[17] = temp[18] = temp[19] = this->w_ambig;
	//temp[20] = temp[31] =  this->w_ambig;
	//temp[10] = temp[12] = temp[13] = temp[14] = temp[15] =  this->w_ambig;
	
	temp[0] = this->w_match;
	temp[1] = temp[2] = temp[3] =  this->w_mismatch;
	temp[4] = temp[5] = temp[6] = temp[7] =  this->w_mismatch;
	temp[16] = temp[17] = temp[18] = temp[19] = this->w_ambig;
	temp[31] =  this->w_ambig;
	temp[12] = temp[13] = temp[14] = temp[15] =  this->w_ambig;
	temp[10] = temp[20] = temp[30] = this->w_mismatch;
		
	__m512i perm512 = _mm512_load_si512(temp);
	

	
	m_b = n_b = 0; b = 0;

	__mmask32 minsc_msk_a = 0x0000, endsc_msk_a = 0x0000;
	int val = 0;
	for (int i=0; i<SIMD_WIDTH16; i++) {
		int xtra = p[i].h0;
		val = (xtra & KSW_XSUBO)? xtra & 0xffff : 0x10000;
		if (val <= SHRT_MAX) {
			minsc[i] = val;
			//if (val < SHRT_MIN)
			//	minsc[i] = -1;
			minsc_msk_a |= (0x1 << i);
		}
		// msc_mask;
		val = (xtra & KSW_XSTOP)? xtra & 0xffff : 0x10000;
		if (val <= SHRT_MAX) {
			endsc[i] = val;
			//if (val < SHRT_MIN)
			//	endsc[i] = -1;
			endsc_msk_a |= (0x1 << i);
		}
		// esc_mask;
		//if (po_ind + i == 4385) {
		//	printf("xtra: %x, minsc: %d %x %d %x %d, nrow: %d %d\n",
		//		   xtra, minsc[i], minsc_msk_a, endsc[i], endsc_msk_a, SIMD_WIDTH16,
		//		   nrow, ncol);
		//}
	}

	__m512i minsc512 = _mm512_load_si512((__m512i*) minsc);
	__m512i endsc512 = _mm512_load_si512((__m512i*) endsc);
	
	__m512i e_del512	= _mm512_set1_epi16(this->e_del);
	__m512i oe_del512	= _mm512_set1_epi16(this->o_del + this->e_del);
	__m512i e_ins512	= _mm512_set1_epi16(this->e_ins);
	__m512i oe_ins512	= _mm512_set1_epi16(this->o_ins + this->e_ins);
	__m512i gmax512		= zero512; // exit1 = zero512;
	__m512i te512		= zero512;	// change to -1
	__mmask32 exit0 = 0xFFFFFFFF;

	
	tid = 0;  // no threading for now !!
	int16_t	*H0		= H16_0 + tid * SIMD_WIDTH16 * MAX_SEQ_LEN_QER_SAM;
	int16_t	*H1		= H16_1 + tid * SIMD_WIDTH16 * MAX_SEQ_LEN_QER_SAM;
	int16_t	*Hmax	= H16_max + tid * SIMD_WIDTH16 * MAX_SEQ_LEN_QER_SAM;
	int16_t	*F		= F16 + tid * SIMD_WIDTH16 * MAX_SEQ_LEN_QER_SAM;
	int16_t	*rowMax	= rowMax16 + tid * SIMD_WIDTH16 * MAX_SEQ_LEN_REF_SAM;
	
	_mm_prefetch((const char*) (F + SIMD_WIDTH16), 0);
	_mm_prefetch((const char*) seq2SoA, 0);
	_mm_prefetch((const char*) seq1SoA, 0);
	_mm_prefetch((const char*) (H1 + SIMD_WIDTH16), 0);
	_mm_prefetch((const char*) (F + SIMD_WIDTH16), 0);

	for (int i=ncol; i >= 0; i--) {
		_mm512_store_si512((__m512*) (H0 + i * SIMD_WIDTH16), zero512);
		_mm512_store_si512((__m512*) (Hmax + i * SIMD_WIDTH16), zero512);
		_mm512_store_si512((__m512*) (F + i * SIMD_WIDTH16), zero512);
	}

	__m512i max512 = zero512, imax512, pimax512 = zero512;
	__mmask32 mask512 = 0x0000;
	__mmask32 minsc_msk = 0x0000;

	__m512i qe512 = _mm512_set1_epi16(-1);
	_mm512_store_si512((__m512i *)(H0), zero512);
	_mm512_store_si512((__m512i *)(H1), zero512);
	__m512i i512 = zero512;
	

#ifdef VTUNE_ANALYSIS
    __itt_resume();
#endif
	int i;
	for (i=0; i < nrow; i++)
	{
        __m512i e11 = zero512;
        __m512i h00, h11, h10, s1;
		int j;
		
		s1 = _mm512_load_si512((__m512i *)(seq1SoA + (i + 0) * SIMD_WIDTH16));
		h10 = zero512;
		imax512 = zero512;
		__m512i iqe512 = _mm512_set1_epi16(-1);
		
		__m512i l512 = zero512;
		for (j=0; j<ncol; j++)
 		{
			__m512i f11, s2, f21;
			h00 = _mm512_load_si512((__m512i *)(H0 + j * SIMD_WIDTH16));
			s2  = _mm512_load_si512((__m512i *)(seq2SoA + (j) * SIMD_WIDTH16));
			f11 = _mm512_load_si512((__m512i *)(F + (j+1) * SIMD_WIDTH16));
			
			MAIN_SAM_CODE16_OPT(s1, s2, h00, h11, e11, f11, f21, max512);

			_mm512_store_si512((__m512i *)(H1 + (j+1) * SIMD_WIDTH16), h11);
			_mm512_store_si512((__m512i *)(F + (j+1) * SIMD_WIDTH16), f21);
			l512 = _mm512_add_epi16(l512, one512);
			// prof[DP2][0] += 22;
			
		}   /* Inner DP loop */
		
		// Block I
		if (i > 0) {
			__mmask32 msk32 = _mm512_cmpgt_epi16_mask(imax512, pimax512);
			msk32 |= mask512;
			// msk32 &= ~minsc_msk;
			pimax512 = _mm512_mask_blend_epi16(msk32, pimax512, minus1);			
			pimax512 = _mm512_mask_blend_epi16(minsc_msk, minus1, pimax512);
			//pimax512 = _mm512_mask_mov_epi16(minus1, minsc_msk, pimax512);

			// new
			// pimax512 = _mm512_mask_blend_epi16(exit0, zero512, pimax512);
			pimax512 = _mm512_mask_blend_epi16(exit0, minus1, pimax512);
			
			_mm512_store_si512((__m512i *) (rowMax + (i-1)*SIMD_WIDTH16), pimax512);
			mask512 = ~msk32;
		}
		pimax512 = imax512;
		minsc_msk = _mm512_cmpge_epi16_mask(imax512, minsc512);
		minsc_msk &= minsc_msk_a;
		
		// Block II: gmax, te
		__mmask32 cmp0 = _mm512_cmpgt_epi16_mask(imax512, gmax512);
		cmp0 &= exit0;
		// gmax512 = _mm512_max_epi16(gmax512, imax512);
		gmax512 = _mm512_mask_blend_epi16(cmp0, gmax512, imax512);
		// te512 = _mm512_mask_mov_epi16(te512, cmp0, i512);
		te512 = _mm512_mask_blend_epi16(cmp0, te512,  i512);
		qe512 = _mm512_mask_blend_epi16(cmp0, qe512,  iqe512);
		
		cmp0 = _mm512_cmpge_epi16_mask(gmax512, endsc512);
		cmp0 &= endsc_msk_a;
		
		exit0 = (~cmp0) & exit0;
		if (exit0 == 0) {
			limit = i++;
			break;
		}
		
		int16_t *S = H1; H1 = H0; H0 = S;
		i512 = _mm512_add_epi16(i512, one512);
	} // for nrow
	pimax512 = _mm512_mask_blend_epi16(mask512, pimax512, minus1);
	pimax512 = _mm512_mask_blend_epi16(minsc_msk, minus1, pimax512);
	pimax512 = _mm512_mask_blend_epi16(exit0, minus1, pimax512);
	_mm512_store_si512((__m512i *) (rowMax + (i-1) * SIMD_WIDTH16), pimax512);
	// __m512i max512_ = max512;
	
	/******************* DP loop over *****************************/
    /*************** Partial output setting ***************/
	int16_t score[SIMD_WIDTH16] __attribute((aligned(64)));
	int16_t te[SIMD_WIDTH16] __attribute((aligned(64)));
	int16_t qe[SIMD_WIDTH16] __attribute((aligned(64)));
	int16_t low[SIMD_WIDTH16] __attribute((aligned(64)));
	int16_t high[SIMD_WIDTH16] __attribute((aligned(64)));	
	_mm512_store_si512((__m512i *) score, gmax512);	
	_mm512_store_si512((__m512i *) te, te512);
	_mm512_store_si512((__m512i *) qe, qe512);


	// int8_t sadj[SIMD_WIDTH16] = {0};
	for (int l=0; l<SIMD_WIDTH16; l++) {
		int ind = po_ind + l;
#if !MAINY
		ind = p[l].regid;    // index of corr. aln
		// if (ind != po_ind + l) printf("ind: %d, po_ind: %d, l: %d\n", ind, po_ind, l);
		// assert(ind == po_ind + l);
		if (phase) {
			if (aln[ind].score == score[l]) {
				aln[ind].tb = aln[ind].te - te[l];
				aln[ind].qb = aln[ind].qe - qe[l];
			}
		} else {
			aln[ind].score = score[l];
			aln[ind].te = te[l];
			aln[ind].qe = qe[l];
		}
#else
		aln[ind].score = score[l];
		aln[ind].te = te[l];
		aln[ind].qe = qe[l];
#endif
	}
	
#if !MAINY
	if (phase) return 1;
#endif


	/*************** Score2 and te2 *******************/
	int qmax = this->g_qmax;	
	// int qmax = 1;
	int maxl = 0 , minh = nrow;
	for (int i=0; i<SIMD_WIDTH16; i++)
	{
		int val = (score[i] + qmax - 1) / qmax;
		low[i] = te[i] - val;
		high[i] = te[i] + val;
		maxl = maxl < low[i] ? low[i] : maxl;
		minh = minh > high[i] ? high[i] : minh;
	}
	max512 = _mm512_set1_epi16(-1);
	te512 = _mm512_set1_epi16(-1);
	__m512i low512 = _mm512_load_si512((__m512i*) low);
	__m512i high512 = _mm512_load_si512((__m512i*) high);

	__m512i rmax512;
	for (int i=0; i< maxl; i++)
	{
		__m512i i512 = _mm512_set1_epi16(i);
		rmax512 = _mm512_load_si512((__m512i*) (rowMax + i*SIMD_WIDTH16));
		__mmask32 mask1 = _mm512_cmpgt_epi16_mask(low512, i512);
		__mmask32 mask2 = _mm512_cmpgt_epi16_mask(rmax512, max512);
		mask2 &= mask1;
		max512 = _mm512_mask_blend_epi16(mask2, max512, rmax512);
		te512 = _mm512_mask_blend_epi16(mask2, te512, i512);
	}	

	for (int i=minh+1; i<limit; i++)
	{
		__m512i i512 = _mm512_set1_epi16(i);
		rmax512 = _mm512_load_si512((__m512i*) (rowMax + i*SIMD_WIDTH16));
		__mmask32 mask1 = _mm512_cmpgt_epi16_mask(i512, high512);
		__mmask32 mask2 = _mm512_cmpgt_epi16_mask(rmax512, max512);
		mask2 &= mask1;
		max512 = _mm512_mask_blend_epi16(mask2, max512, rmax512);
		te512 = _mm512_mask_blend_epi16(mask2, te512, i512);	
	}

	// _mm512_store_si512((__m512i *) temp, max512_);
	_mm512_store_si512((__m512i *) temp1, max512);
	_mm512_store_si512((__m512i *) temp2, te512);		
	for (int i=0; i<SIMD_WIDTH16; i++) {
		int ind = po_ind + i;
#if !MAINY
		ind = p[i].regid;    // index of corr. aln
		assert(ind == po_ind + i);
#endif		
		aln[ind].score2 = temp1[i];
		aln[ind].te2 = temp2[i];
		//if (ind == 4385) {
		//	printf("4385, %d %d %d\n", temp1[i], temp2[i], n_b);
		//	for (int l=0; l<lim; l++)
		//		printf("%d %d (%d %d)\n", l, rowMax[l*SIMD_WIDTH16 + i], po_ind, i);
		//}
#if OUT
		fprintf(stderr, "score: %d, te: %d, qe: %d, score2: %d, te2: %d\n",
				aln[ind].score, aln[ind].te, aln[ind].qe, temp1[i], temp2[i]);
#endif
	}
#endif // ~__AVX512BW__
	return 1;
}

void kswv::kswv512_16(int16_t seq1SoA[],
					 int16_t seq2SoA[],
					 int16_t nrow,
					 int16_t ncol,
					 SeqPair *p,
					 kswr_t *aln,
					 int po_ind,
					 uint16_t tid,
					 int32_t numPairs)
{
#ifdef __AVX512BW__
	int m_b, n_b;
	int16_t minsc[SIMD_WIDTH16], endsc[SIMD_WIDTH16];
	uint64_t *b;

	__m512i zero512 = _mm512_setzero_si512();
	__m512i one512  = _mm512_set1_epi16(1);
	__m512i minus1  = _mm512_set1_epi16(-1);

	int16_t temp[SIMD_WIDTH16] __attribute((aligned(64))) = {0};
	int16_t temp1[SIMD_WIDTH16] __attribute((aligned(64)));
	int16_t temp2[SIMD_WIDTH16] __attribute((aligned(64)));

	temp[0] = this->w_match;
	temp[1] = temp[2] = temp[3] =  this->w_mismatch;
	temp[4] = temp[5] = temp[6] = temp[7] =  this->w_mismatch;
	temp[16] = temp[17] = temp[18] = temp[19] = this->w_ambig;
	temp[20] = temp[31] =  this->w_ambig;
	temp[10] = temp[12] = temp[13] = temp[14] = temp[15] =  this->w_ambig;

	__m512i perm512 = _mm512_load_si512(temp);	
	m_b = n_b = 0; b = 0;

	int val = 0;
	for (int i=0; i<SIMD_WIDTH16; i++) {
		int xtra = p[i].h0;
		val = (xtra & KSW_XSUBO)? xtra & 0xffff : 0x10000;
		minsc[i] = val;
		// msc_mask;
		val = (xtra & KSW_XSTOP)? xtra & 0xffff : 0x10000;
		endsc[i] = val;
		// esc_mask;
	}

	__m512i minsc512 = _mm512_load_si512((__m512i*) minsc);
	
	__m512i e_del512	= _mm512_set1_epi16(this->e_del);
	__m512i oe_del512	= _mm512_set1_epi16(this->o_del + this->e_del);
	__m512i e_ins512	= _mm512_set1_epi16(this->e_ins);
	__m512i oe_ins512	= _mm512_set1_epi16(this->o_ins + this->e_ins);
	__m512i gmax512		= zero512; // exit1 = zero512;
	__m512i te512		= zero512;	// change to -1
	__mmask32 exit0 = 0xFFFFFFFF;

	
	tid = 0;  // no threading for now !!
	int16_t	*H0		= H16_0 + tid * SIMD_WIDTH16 * MAX_SEQ_LEN_QER_SAM;
	int16_t	*H1		= H16_1 + tid * SIMD_WIDTH16 * MAX_SEQ_LEN_QER_SAM;
	int16_t	*Hmax	= H16_max + tid * SIMD_WIDTH16 * MAX_SEQ_LEN_QER_SAM;
	int16_t	*F		= F16 + tid * SIMD_WIDTH16 * MAX_SEQ_LEN_QER_SAM;
	// int16_t	*qp		= qp16 + tid * SIMD_WIDTH16 * MAX_SEQ_LEN_QER_SAM;

	int16_t	*rowMax	= rowMax16 + tid * SIMD_WIDTH16 * MAX_SEQ_LEN_REF_SAM;
	
	_mm_prefetch((const char*) (F + SIMD_WIDTH16), 0);
	_mm_prefetch((const char*) seq2SoA, 0);
	_mm_prefetch((const char*) seq1SoA, 0);
	_mm_prefetch((const char*) (H1 + SIMD_WIDTH16), 0);
	_mm_prefetch((const char*) (F + SIMD_WIDTH16), 0);

	for (int i=ncol; i >= 0; i--) {
		_mm512_store_si512((__m512*) (H0 + i * SIMD_WIDTH16), zero512);
		_mm512_store_si512((__m512*) (Hmax + i * SIMD_WIDTH16), zero512);
		_mm512_store_si512((__m512*) (F + i * SIMD_WIDTH16), zero512);
	}

	__m512i max512 = zero512, imax512, pimax512 = zero512;
	__mmask32 mask512 = 0x0000;
	__mmask32 minsc_msk = 0x0000;

	__m512i qe512 = _mm512_set1_epi16(-1);
	_mm512_store_si512((__m512i *)(H0), zero512);
	_mm512_store_si512((__m512i *)(H1), zero512);
	__m512i i512 = zero512;;

	// uint64_t tim = __rdtsc();	
	int i;
	for (i=0; i < nrow; i++)
	{
        __m512i e11 = zero512;
        __m512i h00, h11, h10, s1;
		
		s1 = _mm512_load_si512((__m512i *)(seq1SoA + (i + 0) * SIMD_WIDTH16));
		h10 = zero512;
		imax512 = zero512;
		__m512i iqe512 = _mm512_set1_epi16(-1);
		
		// int64_t tim_ = __rdtsc();		
		__m512i l512 = zero512;
		int j;
		for (j=0; j<ncol; j++)
		{
			__m512i f11, s2, f21;
			h00 = _mm512_load_si512((__m512i *)(H0 + j * SIMD_WIDTH16));
			s2  = _mm512_load_si512((__m512i *)(seq2SoA + (j) * SIMD_WIDTH16));
			f11 = _mm512_load_si512((__m512i *)(F + (j+1) * SIMD_WIDTH16));

			// my main code
			MAIN_SAM_CODE16_OPT(s1, s2, h00, h11, e11, f11, f21, max512);
			
			_mm512_store_si512((__m512i *)(H1 + (j+1) * SIMD_WIDTH16), h11);  
			_mm512_store_si512((__m512i *)(F + (j+1) * SIMD_WIDTH16), f21);
			l512 = _mm512_add_epi16(l512, one512);
			// prof[DP2][0] += 22;
			
		} // core DP loop
		// prof[DP1][0] += __rdtsc() - tim_;
				
		// TODO: Block I
		if (i > 0) {
			__mmask32 msk32 = _mm512_cmpgt_epi16_mask(imax512, pimax512);
			msk32 |= mask512;
			// msk32 &= ~minsc_msk;
			pimax512 = _mm512_mask_blend_epi16(msk32, pimax512, minus1);			
			pimax512 = _mm512_mask_blend_epi16(minsc_msk, minus1, pimax512);
			//pimax512 = _mm512_mask_mov_epi16(minus1, minsc_msk, pimax512);
			
			_mm512_store_si512((__m512i *) (rowMax + (i-1)*SIMD_WIDTH16), pimax512);
			
			mask512 = ~msk32;
		}
		pimax512 = imax512;
		minsc_msk = _mm512_cmpge_epi16_mask(imax512, minsc512);
		
		// Block II: gmax, te
		__mmask32 cmp0 = _mm512_cmpgt_epi16_mask(imax512, gmax512);
		gmax512 = _mm512_max_epi16(gmax512, imax512);
		// te512 = _mm512_mask_mov_epi16(te512, cmp0, i512);
		te512 = _mm512_mask_blend_epi16(cmp0, te512,  i512);
		qe512 = _mm512_mask_blend_epi16(cmp0, qe512,  iqe512);
		
		cmp0 &= exit0; // check exit code
		
		int16_t *S = H1; H1 = H0; H0 = S;
		i512 = _mm512_add_epi16(i512, one512);
	} // for nrow
	
	pimax512 = _mm512_mask_blend_epi16(mask512, pimax512, minus1);
	pimax512 = _mm512_mask_blend_epi16(minsc_msk, minus1, pimax512);
	_mm512_store_si512((__m512i *) (rowMax + (i-1) * SIMD_WIDTH16), pimax512);
	
	// prof[DP][0] += __rdtsc() - tim;

    /*************** Partial output setting ***************/
	int16_t score[SIMD_WIDTH16] __attribute((aligned(64)));
	int16_t te[SIMD_WIDTH16] __attribute((aligned(64)));
	int16_t qe[SIMD_WIDTH16] __attribute((aligned(64)));
	int16_t low[SIMD_WIDTH16] __attribute((aligned(64)));
	int16_t high[SIMD_WIDTH16] __attribute((aligned(64)));	
	_mm512_store_si512((__m512i *) score, gmax512);	
	_mm512_store_si512((__m512i *) te, te512);
	_mm512_store_si512((__m512i *) qe, qe512);
	
	for (int l=0; l<SIMD_WIDTH16; l++) {
		int ind = po_ind + l;
#if !MAINY
		ind = p[l].regid;    // index of corr. aln
		assert(ind == po_ind + l);
#endif		
		aln[ind].score = score[l];
		aln[ind].te = te[l];
		aln[ind].qe = qe[l];
	}

	/*************** Score2 and te2 *******************/
	int qmax = this->g_qmax;	
	// int qmax = 1;
	int maxl = 0 , minh = nrow;
	for (int i=0; i<SIMD_WIDTH16; i++)
	{
		int val = (score[i] + qmax - 1) / qmax;
		low[i] = te[i] - val;
		high[i] = te[i] + val;
		maxl = maxl < low[i] ? low[i] : maxl;
		minh = minh > high[i] ? high[i] : minh;
	}
	max512 = _mm512_set1_epi16(-1);
	te512 = _mm512_set1_epi16(-1);
	__m512i low512 = _mm512_load_si512((__m512i*) low);
	__m512i high512 = _mm512_load_si512((__m512i*) high);

	
	__m512i rmax512;
	for (int i=0; i< maxl; i++)
	{
		__m512i i512 = _mm512_set1_epi16(i);
		rmax512 = _mm512_load_si512((__m512i*) (rowMax + i*SIMD_WIDTH16));
		__mmask32 mask1 = _mm512_cmpgt_epi16_mask(low512, i512);
		__mmask32 mask2 = _mm512_cmpgt_epi16_mask(rmax512, max512);
		mask2 &= mask1;
		max512 = _mm512_mask_blend_epi16(mask2, max512, rmax512);
		te512 = _mm512_mask_blend_epi16(mask2, te512, i512);
	}	
	for (int i=minh+1; i< nrow; i++)
	{
		__m512i i512 = _mm512_set1_epi16(i);
		rmax512 = _mm512_load_si512((__m512i*) (rowMax + i*SIMD_WIDTH16));
		__mmask32 mask1 = _mm512_cmpgt_epi16_mask(i512, high512);
		__mmask32 mask2 = _mm512_cmpgt_epi16_mask(rmax512, max512);
		mask2 &= mask1;
		max512 = _mm512_mask_blend_epi16(mask2, max512, rmax512);
		te512 = _mm512_mask_blend_epi16(mask2, te512, i512);	
	}

	//_mm512_store_si512((__m512i *) temp, max512_);
	_mm512_store_si512((__m512i *) temp1, max512);
	_mm512_store_si512((__m512i *) temp2, te512);		
	for (int i=0; i<SIMD_WIDTH16; i++) {
		int ind = po_ind + i;
#if !MAINY
		ind = p[i].regid;    // index of corr. aln
		assert(ind == po_ind + i);
#endif		
		aln[ind].score2 = temp1[i];
		aln[ind].te2 = temp2[i];
	//	// fprintf(stderr, "gmax: %d, te: %d, qe: %d, max: %d, max2: %d, te2: %d\n",
	//	// aln[ind].score, aln[ind].te, aln[ind].qe, temp[i], temp1[i], temp2[i]);		
	}
#endif
}

/********************** Intra-Task stuff ************************/
kswqi_t* kswv::ksw_qinit_intra(int size, int qlen, uint8_t *query, int m, const int8_t *mat)
{
	kswqi_t *q;
	int slen, a, tmp, p;

	size = size > 1? 2 : 1;
	p = 32 * (3 - size); // # values per __m128i
	slen = (qlen + p - 1) / p; // segmented length

	//q = (kswqi_t*)_mm_malloc(sizeof(kswqi_t) + 256 + 16 * slen * (m + 4), 64); // a single block of memory
	q = (kswqi_t*)_mm_malloc(sizeof(kswqi_t) + 256 + 64 * slen * (m + 4), 64);
	q->qp = (__m512i*)(((size_t)q + sizeof(kswqi_t) + 64) >> 4 << 4); // align memory
	//q->qp = (__m512i*)((size_t)q); // align memory
	
	q->H0 = q->qp + slen * m;
	q->H1 = q->H0 + slen;
	q->E  = q->H1 + slen;
	q->Hmax = q->E + slen;
	q->slen = slen; q->qlen = qlen; q->size = size;
	// compute shift
	tmp = m * m;
	for (a = 0, q->shift = 127, q->mdiff = 0; a < tmp; ++a) { // find the minimum and maximum score
		if (mat[a] < (int8_t)q->shift) q->shift = mat[a];
		if (mat[a] > (int8_t)q->mdiff) q->mdiff = mat[a];
	}
	q->max = q->mdiff;
	q->shift = 256 - q->shift; // NB: q->shift is uint8_t
	q->mdiff += q->shift; // this is the difference between the min and max scores
	// An example: p=8, qlen=19, slen=3 and segmentation:
	//  {{0,3,6,9,12,15,18,-1},{1,4,7,10,13,16,-1,-1},{2,5,8,11,14,17,-1,-1}}
	if (size == 1) {
		int8_t *t = (int8_t*)q->qp;
		assert(slen * p < MAX_SEQ_LEN_QER_SAM);
		for (a = 0; a < m; ++a) {
			int i, k, nlen = slen * p;
			const int8_t *ma = mat + a * m;
			for (i = 0; i < slen; ++i)
				for (k = i; k < nlen; k += slen) // p iterations
					*t++ = (k >= qlen? 0 : ma[query[k]]) + q->shift;
		}
	} else {
		int16_t *t = (int16_t*)q->qp;
		assert(slen * p < MAX_SEQ_LEN_QER_SAM);
		for (a = 0; a < m; ++a) {
			int i, k, nlen = slen * p;
			const int8_t *ma = mat + a * m;
			for (i = 0; i < slen; ++i)
				for (k = i; k < nlen; k += slen) // p iterations
					*t++ = (k >= qlen? 0 : ma[query[k]]);
		}
	}
	return q;
}

void kswv::kswvBatchWrapper16_intra(SeqPair *pairArray,
									uint8_t *seqBufRef,
									uint8_t *seqBufQer,
									int32_t numPairs,
									uint16_t numThreads)
{

	//int16_t *seq1SoA = NULL;
	//int16_t *seq2SoA = NULL;		
	//seq1SoA = (int16_t *)_mm_malloc(MAX_SEQ_LEN_REF_SAM * numThreads * sizeof(int16_t), 64);
	//seq2SoA = (int16_t *)_mm_malloc(MAX_SEQ_LEN_QER_SAM * numThreads * sizeof(int16_t), 64);
	//if (seq1SoA == NULL || seq2SoA == NULL)	printf("Mem not allocated!!!\n");

	int i = 0;
	int8_t mat[25];
	// bwa_fill_scmat(1, 4, mat);
	bwa_fill_scmat(mat);

//#pragma omp parallel num_threads(numThreads)
    {
		uint16_t tid = 0;
		uint8_t *seq1;
		uint8_t *seq2;

		int nstart = 0, nend = numPairs;
		
//#pragma omp for schedule(dynamic, 128)
		for(i = nstart; i < nend; i++)
		{
			SeqPair sp = pairArray[i];
			seq1 = seqBufRef + (int64_t)sp.id * MAX_SEQ_LEN_REF_SAM;
			//for(int k = 0; k < sp.len1; k++)
			//	seq1SoA[k] = seq1[k];

			target = seq1;
			seq2 = seqBufQer + (int64_t)sp.id * MAX_SEQ_LEN_QER_SAM;
			query = seq2;
			//for(k = 0; k < sp.len2; k++)
			//	seq2SoA[k] = seq2[k];

			// printf("Check 1\n");
			// query profile
			// kswq_t *q = ksw_qinit(2, qlen, query, this->m, mat);
			kswqi_t *q = ksw_qinit_intra(2, sp.len2, seq2, this->m, mat);
			
			int maxi = -1;
			uint64_t tim_mate = __rdtsc();
            kswv512_16_intra(seq1, q,
                             sp.len1, sp.len2,
                             pairArray + i,
                             tid,
                             numPairs,
                             maxi);
			prof[0][tid] += __rdtsc() - tim_mate;
			_mm_free(q);
		}
	}

	//_mm_free(seq1SoA);
	//_mm_free(seq2SoA);
	return;
}

kswr_t kswv::kswv512_16_intra(uint8_t seq1SoA[],
							  kswqi_t *q,
							  int16_t nrow,
							  int16_t ncol,
							  SeqPair *p,
							  uint16_t tid,
							  int32_t numPairs,
							  int &maxi)
{
	kswr_t r = g_defr;
#ifdef __AVX512BW__
	int slen, i, m_b, n_b, te = -1, gmax = 0, minsc, endsc;
	uint64_t *b;
	__m512i zero, oe_del, e_del, oe_ins, e_ins, *H0, *H1, *E, *Hmax;

	// printf ("nrow: %d, ncol: %d\n", nrow, ncol);
    // exit(0);
	
#define SIMD16I 32
	int16_t temp[SIMD16I] __attribute((aligned(64)));

#define __max_8I(ret, xx) do {								  \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 16)); \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 8));  \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 4));  \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 2));  \
    	(ret) = _mm_extract_epi16((xx), 0);					  \
	} while (0)
	
	// initialization
	int xtra = p->h0;

	r = g_defr;
	minsc = (xtra&KSW_XSUBO)? xtra&0xffff : 0x10000;
	endsc = (xtra&KSW_XSTOP)? xtra&0xffff : 0x10000;
	m_b = n_b = 0; b = 0;
	zero = _mm512_setzero_si512();
	oe_del = _mm512_set1_epi16(this->o_del + this->e_del);
	e_del  = _mm512_set1_epi16(this->e_del);
	oe_ins = _mm512_set1_epi16(this->o_ins + this->e_ins);
	e_ins  = _mm512_set1_epi16(this->e_ins);
	H0 = q->H0; H1 = q->H1; E = q->E; Hmax = q->Hmax;
	slen = q->slen;
	for (i = 0; i < slen; ++i) {
		_mm512_store_si512(E + i, zero);
		_mm512_store_si512(H0 + i, zero);
		_mm512_store_si512(Hmax + i, zero);
	}
	__m512i idx = _mm512_set_epi16(30, 29, 28, 27, 26, 25, 24, 23, 22,
								   21, 20, 19, 18, 17, 16, 15, 14, 13,
								   12, 11, 10, 9 , 8, 7, 6, 5, 4, 3,
								   2, 1 ,0 ,31);

	__m512i des = _mm512_set_epi16(0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
								   0xFFFF, 0xFFFF,
								   0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
								   0xFFFF, 0xFFFF,
								   0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
								   0xFFFF, 0xFFFF,
								   0xFFFF, 0xFFFF, 0xFFFF, 0);
	__mmask32 lmask[slen];
	__mmask32 msk = 0xFFFFFFFF;

	int val = (ncol + 1) % 8;
	int lim = ncol + 0 + (val > 0 ? 8 - val: 0);
	// printf ("ncol: %d, lim: %d, slen: %d, len: %d\n", ncol, lim, slen, slen * SIMD16I);
	for (int i=0; i<slen; i++) {
		int val = slen * (SIMD16I - 1) + i - 1;
		// printf("i: %d, val: %d, lim: %d, div: %d\n", i, val, lim, (val - lim)/slen);
		lmask[i] = msk >> (val - lim < 0? 0 : (val - lim)/slen+1) ;
	}
	
	// the core loop
	// uint64_t tim = __rdtsc();
	for (i = 0; i < nrow; ++i) {
		int j, k, imax;
		__m512i e, t, h, f = zero, max = zero, *S = q->qp + seq1SoA[i] * slen; // s
		h = _mm512_load_si512(H0 + slen - 1); // h={2,5,8,11,14,17,-1,-1} in the above example
		// h = _mm512_slli_si512(h, 2);

		h = _mm512_permutexvar_epi16(idx, h);
		h = _mm512_and_si512(h, des);

		// uint64_t tim_ = __rdtsc();
		for (j = 0; LIKELY(j < slen); ++j)
		{
			h = _mm512_mask_blend_epi16(lmask[j], zero, h);
			h = _mm512_adds_epi16(h, *S++);						
			e = _mm512_load_si512(E + j);
			e = _mm512_mask_blend_epi16(lmask[j], zero, e);  // new
			h = _mm512_max_epi16(h, e);
			h = _mm512_max_epi16(h, f);
			max = _mm512_max_epi16(max, h);
			_mm512_store_si512(H1 + j, h);			
			e = _mm512_subs_epu16(e, e_del);
			t = _mm512_subs_epu16(h, oe_del);
			e = _mm512_max_epi16(e, t);
			_mm512_store_si512(E + j, e);
			f = _mm512_subs_epu16(f, e_ins);
			t = _mm512_subs_epu16(h, oe_ins);
			f = _mm512_max_epi16(f, t);
			h = _mm512_load_si512(H0 + j);
			prof[DP2][0] += 16;
		}
		// prof[DP1][0] += __rdtsc() - tim_;
		
		for (k = 0; LIKELY(k < SIMD16I); ++k) {
			//f = _mm512_slli_si512(f, 2);
			f = _mm512_permutexvar_epi16(idx, f);
			f = _mm512_and_si512(f, des);
			for (j = 0; LIKELY(j < slen); ++j) {
				h = _mm512_load_si512(H1 + j);
				h = _mm512_max_epi16(h, f);
				_mm512_store_si512(H1 + j, h);
				h = _mm512_subs_epu16(h, oe_ins);
				f = _mm512_subs_epu16(f, e_ins);
				//if(UNLIKELY(!_mm512_movepi16_mask(_mm512_cmpgt_epi16(f, h)))) goto end_loop8;
				if((!_mm512_cmpgt_epi16_mask(f, h))) goto end_loop8;
			}
		}
end_loop8:
		// __max_8(imax, max);
		_mm512_store_si512(temp, max);
		imax = temp[0];		
		for (int l=1; l<SIMD16I; l++)
			if (imax < temp[l]) imax = temp[l];
		// maxi = imax;
		
		if (imax >= minsc) {
			if (n_b == 0 || (int32_t)b[n_b-1] + 1 != i) {
				if (n_b == m_b) {
					m_b = m_b? m_b<<1 : 8;
					b = (uint64_t*)realloc(b, 8 * m_b);
				}
				b[n_b++] = (uint64_t)imax<<32 | i;
			} else if ((int)(b[n_b-1]>>32) < imax) b[n_b-1] = (uint64_t)imax<<32 | i; // modify the last
		}
		if (imax > gmax) {
			maxi = imax;
			gmax = imax; te = i;
			for (j = 0; LIKELY(j < slen); ++j)
				_mm512_store_si512(Hmax + j, _mm512_load_si512(H1 + j));
			if (gmax >= endsc) break;
		}
		S = H1; H1 = H0; H0 = S;
	}
	// prof[DP][0] += __rdtsc() - tim;

	int max_ = -1;
	r.score = gmax; r.te = te;
	int max = -1, tmp, low, high, qlen = slen * 32;
	{
		uint16_t *t = (uint16_t*)Hmax;
		for (i = 0, r.qe = -1; i < qlen; ++i, ++t)
			if ((int)*t > max) max = *t, r.qe = i / 32 + i % 32 * slen;
			else if ((int)*t == max && (tmp = i / 32 + i % 32 * slen) < r.qe) r.qe = tmp;
		max_ = max;
		if (b) {
			assert(q->max != 0);
			i = (r.score + q->max - 1) / q->max;
			low = te - i; high = te + i;
			for (i = 0; i < n_b; ++i) {
				int e = (int32_t)b[i];
				if ((e < low || e > high) && (int)(b[i]>>32) > r.score2)
					r.score2 = b[i]>>32, r.te2 = e;
			}
		}
	}

	free(b);
#endif
	return r;
}
