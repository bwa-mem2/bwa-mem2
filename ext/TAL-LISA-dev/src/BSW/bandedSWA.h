/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

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

Authors: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
*****************************************************************************************/

#ifndef SCALAR_BANDEDSWA_HPP
#define SCALAR_BANDEDSWA_HPP

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include "macro.h"

#if (__AVX512BW__ || __AVX2__)
#include <immintrin.h>
#else
#include <smmintrin.h>  // for SSE4.1
#define __mmask8 uint8_t
#define __mmask16 uint16_t
#endif

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
#define min_(x, y) ((x)>(y)?(y):(x))
#define max_(x, y) ((x)>(y)?(x):(y))

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
                     const int8_t w_match, const int8_t w_mismatch, int numThreads);
    ~BandedPairWiseSW();
    // Scalar code section
    int scalarBandedSWA(int qlen, const uint8_t *query, int tlen,
                        const uint8_t *target, int32_t w,
                        int h0, int *_qle, int *_tle,
                        int *_gtle, int *_gscore,
                        int *_max_off);

    void scalarBandedSWAWrapper(SeqPair *seqPairArray,
                                uint8_t *seqBufRef,
                                uint8_t *seqBufQer,
                                int numPairs,
                                int nthreads,
                                int32_t w);

#if ((!__AVX512BW__) & (!__AVX2__) & (__SSE2__))
    // AVX256 is not updated for banding and separate ins/del in the inner loop.
    // 8 bit vector code section    
    void getScores8(SeqPair *pairArray,
                    uint8_t *seqBufRef,
                    uint8_t *seqBufQer,
                    int32_t numPairs,
                    uint16_t numThreads,
                    int32_t w);

    void smithWatermanBatchWrapper8(SeqPair *pairArray,
                                   uint8_t *seqBufRef,
                                   uint8_t *seqBufQer,
                                   int32_t numPairs,
                                   uint16_t numThreads,
                                   int32_t w);

    void smithWaterman128_8(uint8_t seq1SoA[],
                            uint8_t seq2SoA[],
                            uint8_t nrow,
                            uint8_t ncol,
                            SeqPair *p,
                            uint8_t h0[],
                            uint16_t tid,
                            int32_t numPairs,
                            int zdrop,
                            int32_t w,
                            uint8_t qlen[],
                            uint8_t myband[]);
    // 16 bit vector code section
    void getScores16(SeqPair *pairArray,
                     uint8_t *seqBufRef,
                     uint8_t *seqBufQer,
                     int32_t numPairs,
                     uint16_t numThreads,
                     int32_t w);

    void smithWatermanBatchWrapper16(SeqPair *pairArray,
                                     uint8_t *seqBufRef,
                                     uint8_t *seqBufQer,
                                     int32_t numPairs,
                                     uint16_t numThreads,
                                     int32_t w);
    
    void smithWaterman128_16(uint16_t seq1SoA[],
                             uint16_t seq2SoA[],
                             uint16_t nrow,
                             uint16_t ncol,
                             SeqPair *p,
                             uint16_t h0[],
                             uint16_t tid,
                             int32_t numPairs,
                             int zdrop,
                             int32_t w,
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
                    int32_t w);

    void smithWatermanBatchWrapper8(SeqPair *pairArray,
                                   uint8_t *seqBufRef,
                                   uint8_t *seqBufQer,
                                   int32_t numPairs,
                                   uint16_t numThreads,
                                   int32_t w);

    void smithWaterman256_8(uint8_t seq1SoA[],
                            uint8_t seq2SoA[],
                            uint8_t nrow,
                            uint8_t ncol,
                            SeqPair *p,
                            uint8_t h0[],
                            uint16_t tid,
                            int32_t numPairs,
                            int zdrop,
                            int32_t w,
                            uint8_t qlen[],
                            uint8_t myband[]);
    // 16 bit vector code section
    void getScores16(SeqPair *pairArray,
                     uint8_t *seqBufRef,
                     uint8_t *seqBufQer,
                     int32_t numPairs,
                     uint16_t numThreads,
                     int32_t w);

    void smithWatermanBatchWrapper16(SeqPair *pairArray,
                                     uint8_t *seqBufRef,
                                     uint8_t *seqBufQer,
                                     int32_t numPairs,
                                     uint16_t numThreads,
                                     int32_t w);
    
    void smithWaterman256_16(uint16_t seq1SoA[],
                             uint16_t seq2SoA[],
                             uint16_t nrow,
                             uint16_t ncol,
                             SeqPair *p,
                             uint16_t h0[],
                             uint16_t tid,
                             int32_t numPairs,
                             int zdrop,
                             int32_t w,
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
                    int32_t w);

    void smithWatermanBatchWrapper8(SeqPair *pairArray,
                                   uint8_t *seqBufRef,
                                   uint8_t *seqBufQer,
                                   int32_t numPairs,
                                   uint16_t numThreads,
                                   int32_t w);

    void smithWaterman512_8(uint8_t seq1SoA[],
                            uint8_t seq2SoA[],
                            uint8_t nrow,
                            uint8_t ncol,
                            SeqPair *p,
                            uint8_t h0[],
                            uint16_t tid,
                            int32_t numPairs,
                            int zdrop,
                            int32_t w,
                            uint8_t qlen[],
                            uint8_t myband[]);

    // 16 bit vector code section
    void getScores16(SeqPair *pairArray,
                     uint8_t *seqBufRef,
                     uint8_t *seqBufQer,
                     int32_t numPairs,
                     uint16_t numThreads,
                     int32_t w);

    void smithWatermanBatchWrapper16(SeqPair *pairArray,
                                     uint8_t *seqBufRef,
                                     uint8_t *seqBufQer,
                                     int32_t numPairs,
                                     uint16_t numThreads,
                                     int32_t w);
    
    void smithWaterman512_16(uint16_t seq1SoA[],
                             uint16_t seq2SoA[],
                             uint16_t nrow,
                             uint16_t ncol,
                             SeqPair *p,
                             uint16_t h0[],
                             uint16_t tid,
                             int32_t numPairs,
                             int zdrop,
                             int32_t w,
                             uint16_t qlen[],
                             uint16_t myband[]);
#endif

    int64_t getTicks();
    
private:
    int m;
    int end_bonus, zdrop;
    int o_del, o_ins, e_del, e_ins;
    const int8_t *mat;

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

#endif
