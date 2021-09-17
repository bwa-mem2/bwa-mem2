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

#ifndef _FMI_SEARCH_H
#define _FMI_SEARCH_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <immintrin.h>
#include <limits.h>
#include <fstream>

#include "bntseq.h"
#include "bwt.h"
#include "bseq.h"
#include "utils.h"

#define DUMMY_CHAR 6

#define SA_COMPRESSION 1
#define SA_COMPX 03 // (= power of 2)
#define SA_COMPX_MASK 0x7    // 0x7 or 0x3 or 0x1


#define assert_not_null(x, size, cur_alloc) \
        if (x == NULL) { fprintf(stderr, "Allocation of %0.2lf GB for " #x " failed.\nCurrent Allocation = %0.2lf GB\n", size * 1.0 /(1024*1024*1024), cur_alloc * 1.0 /(1024*1024*1024)); exit(EXIT_FAILURE); }

#define MAX_FIFO_SIZE 1048576

#define CP_BLOCK_SIZE 64
#define CP_FILENAME_SUFFIX ".bwt.2bit.64"
#define CP_MASK 63
#define CP_SHIFT 6

typedef struct checkpoint_occ_scalar
{
    int64_t cp_count[4];
    uint64_t one_hot_bwt_str[4];
}CP_OCC;

#if defined(__clang__) || defined(__GNUC__)
static inline int _mm_countbits_64(unsigned long x) {
    return __builtin_popcountl(x);
}
#endif

#define \
GET_OCC(pp, c, occ_id_pp, y_pp, occ_pp, one_hot_bwt_str_c_pp, match_mask_pp) \
                int64_t occ_id_pp = pp >> CP_SHIFT; \
                int64_t y_pp = pp & CP_MASK; \
                int64_t occ_pp = cp_occ[occ_id_pp].cp_count[c]; \
                uint64_t one_hot_bwt_str_c_pp = cp_occ[occ_id_pp].one_hot_bwt_str[c]; \
                uint64_t match_mask_pp = one_hot_bwt_str_c_pp & one_hot_mask_array[y_pp]; \
                occ_pp += _mm_countbits_64(match_mask_pp);



#if OLD
#define CP_BLOCK_SIZE_SCALAR 64
#define CP_FILENAME_SUFFIX_SCALAR ".bwt.2bit.64"
#define CP_MASK_SCALAR 63
#define CP_SHIFT_SCALAR 6
#define BIT_DATA_TYPE uint64_t
#define PADDING_SCALAR 8

#define CP_BLOCK_SIZE_AVX 32
#define CP_FILENAME_SUFFIX_AVX ".bwt.8bit.32"
#define CP_MASK_AVX 31
#define CP_SHIFT_AVX 5

typedef struct checkpoint_occ_scalar
{
    BIT_DATA_TYPE bwt_str_bit0;
    BIT_DATA_TYPE bwt_str_bit1;
    BIT_DATA_TYPE dollar_mask;
    int64_t cp_count[4];
    uint8_t  pad[PADDING_SCALAR];
}CP_OCC_SCALAR;

typedef struct checkpoint_occ_avx
{
    uint8_t  bwt_str[CP_BLOCK_SIZE_AVX];
    int64_t cp_count[4];
}CP_OCC_AVX;

#if ((!__AVX2__))

typedef CP_OCC_SCALAR CP_OCC;
#define CP_SHIFT CP_SHIFT_SCALAR
#define CP_FILENAME_SUFFIX CP_FILENAME_SUFFIX_SCALAR
#define CP_MASK CP_MASK_SCALAR
#define CP_BLOCK_SIZE CP_BLOCK_SIZE_SCALAR

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

#else

typedef CP_OCC_AVX CP_OCC;
#define CP_SHIFT CP_SHIFT_AVX
#define CP_FILENAME_SUFFIX CP_FILENAME_SUFFIX_AVX
#define CP_MASK CP_MASK_AVX
#define CP_BLOCK_SIZE CP_BLOCK_SIZE_AVX

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

#endif
#endif


#if 0
#if defined (__2BIT_LEVEL__)

#if (CP_BLOCK_SIZE == 32) || (CP_BLOCK_SIZE == 64)

#if (CP_BLOCK_SIZE == 64)

#define CP_MASK 63
#define CP_SHIFT 6
#define BIT_DATA_TYPE uint64_t
#define PADDING 24
#define _MM_COUNTBITS _mm_countbits_64

#else
#define CP_MASK 31
#define CP_SHIFT 5
#define BIT_DATA_TYPE uint32_t
#define PADDING 4
#define _MM_COUNTBITS _mm_countbits_32

#endif

#define \
GET_OCC(pp, c, occ_id_pp, y_pp, occ_pp, bwt_str_bit0_pp, bwt_str_bit1_pp, bit0_cmp_pp, bit1_cmp_pp, mismatch_mask_pp) \
                uint32_t occ_id_pp = pp >> CP_SHIFT; \
                uint32_t y_pp = pp & CP_MASK; \
                uint32_t occ_pp = my_cp_occ[occ_id_pp].cp_count[c]; \
                if(y_pp > 0) \
                { \
                BIT_DATA_TYPE bwt_str_bit0_pp = my_cp_occ[occ_id_pp].bwt_str_bit0; \
                BIT_DATA_TYPE bwt_str_bit1_pp = my_cp_occ[occ_id_pp].bwt_str_bit1; \
                BIT_DATA_TYPE bit0_cmp_pp = bwt_str_bit0_pp ^ base_mask[c][0]; \
                BIT_DATA_TYPE bit1_cmp_pp = bwt_str_bit1_pp ^ base_mask[c][1]; \
                uint64_t mismatch_mask_pp = bit0_cmp_pp | bit1_cmp_pp | my_cp_occ[occ_id_pp].dollar_mask; \
                mismatch_mask_pp = mismatch_mask_pp >> (CP_BLOCK_SIZE - y_pp); \
                occ_pp += y_pp - _MM_COUNTBITS(mismatch_mask_pp); \
                }

typedef struct checkpoint_occ
{
    BIT_DATA_TYPE bwt_str_bit0;
    BIT_DATA_TYPE bwt_str_bit1;
    BIT_DATA_TYPE dollar_mask;
    uint32_t cp_count[4];
    uint8_t  pad[PADDING];
}CP_OCC;

#elif (CP_BLOCK_SIZE == 128)

#define CP_MASK 127
#define CP_SHIFT 7
#define BIT_DATA_TYPE uint64_t
#define PADDING 0
#define ITER 2
#define _MM_COUNTBITS _mm_countbits_64


#define \
GET_OCC(pp, c, occ_id_pp, y_pp, y0_pp, y1_pp, occ_pp, bwt_str_bit0_pp, bwt_str_bit1_pp, bit0_cmp_pp, bit1_cmp_pp, mismatch_mask_pp, num_match_pp) \
                uint32_t occ_id_pp = pp >> CP_SHIFT; \
                uint32_t y_pp = pp & CP_MASK; \
                uint32_t occ_pp = my_cp_occ[occ_id_pp].cp_count[c]; \
                \
                uint64_t bwt_str_bit0_pp, bwt_str_bit1_pp; \
                uint64_t bit0_cmp_pp, bit1_cmp_pp; \
                uint64_t mismatch_mask_pp; \
                uint32_t y0_pp = y_pp; \
                if(y0_pp > 64) y0_pp = 64;\
                if(y0_pp > 0) \
                { \
                bwt_str_bit0_pp = my_cp_occ[occ_id_pp].bwt_str_bit0[0]; \
                bwt_str_bit1_pp = my_cp_occ[occ_id_pp].bwt_str_bit1[0]; \
                bit0_cmp_pp = bwt_str_bit0_pp ^ base_mask[c][0]; \
                bit1_cmp_pp = bwt_str_bit1_pp ^ base_mask[c][1]; \
                mismatch_mask_pp = bit0_cmp_pp | bit1_cmp_pp | my_cp_occ[occ_id_pp].dollar_mask[0]; \
                mismatch_mask_pp = mismatch_mask_pp >> (64 - y0_pp); \
                occ_pp += y0_pp - _MM_COUNTBITS(mismatch_mask_pp); \
                \
                uint32_t y1_pp = y_pp - y0_pp; \
                if(y1_pp > 0) \
                {\
                bwt_str_bit0_pp = my_cp_occ[occ_id_pp].bwt_str_bit0[1]; \
                bwt_str_bit1_pp = my_cp_occ[occ_id_pp].bwt_str_bit1[1]; \
                bit0_cmp_pp = bwt_str_bit0_pp ^ base_mask[c][0]; \
                bit1_cmp_pp = bwt_str_bit1_pp ^ base_mask[c][1]; \
                mismatch_mask_pp = bit0_cmp_pp | bit1_cmp_pp | my_cp_occ[occ_id_pp].dollar_mask[1]; \
                mismatch_mask_pp = mismatch_mask_pp >> (64 - y1_pp); \
                occ_pp += (y1_pp - _MM_COUNTBITS(mismatch_mask_pp)); \
                }\
                }\

typedef struct checkpoint_occ
{
    uint64_t bwt_str_bit0[ITER];
    uint64_t bwt_str_bit1[ITER];
    uint64_t dollar_mask[ITER];
    uint32_t cp_count[4];
    //uint8_t  pad[0];
}CP_OCC;

#elif (CP_BLOCK_SIZE == 256)

#define CP_MASK 255
#define CP_SHIFT 8
#define BIT_DATA_TYPE uint64_t
#define PADDING 16
#define ITER 4
#define _MM_COUNTBITS _mm_countbits_64


#define \
GET_OCC(pp, c, occ_id_pp, y_pp, y0_pp, y1_pp, occ_pp, bwt_str_bit0_pp, bwt_str_bit1_pp, bit0_cmp_pp, bit1_cmp_pp, mismatch_mask_pp, num_match_pp) \
                uint32_t occ_id_pp = pp >> CP_SHIFT; \
                uint32_t y_pp = pp & CP_MASK; \
                uint32_t occ_pp = my_cp_occ[occ_id_pp].cp_count[c]; \
                \
                uint64_t bwt_str_bit0_pp, bwt_str_bit1_pp; \
                uint64_t bit0_cmp_pp, bit1_cmp_pp; \
                uint64_t mismatch_mask_pp; \
                x = 0; \
                uint32_t y0_pp = y_pp; \
                while(y0_pp > 64) \
                {\
                bwt_str_bit0_pp = my_cp_occ[occ_id_pp].bwt_str_bit0[x]; \
                bwt_str_bit1_pp = my_cp_occ[occ_id_pp].bwt_str_bit1[x]; \
                bit0_cmp_pp = bwt_str_bit0_pp ^ base_mask[c][0]; \
                bit1_cmp_pp = bwt_str_bit1_pp ^ base_mask[c][1]; \
                mismatch_mask_pp = bit0_cmp_pp | bit1_cmp_pp | my_cp_occ[occ_id_pp].dollar_mask[x]; \
                occ_pp += 64 - _MM_COUNTBITS(mismatch_mask_pp); \
                x++;\
                y0_pp -= 64;\
                }\
                \
                if(y0_pp > 0) \
                {\
                bwt_str_bit0_pp = my_cp_occ[occ_id_pp].bwt_str_bit0[x]; \
                bwt_str_bit1_pp = my_cp_occ[occ_id_pp].bwt_str_bit1[x]; \
                bit0_cmp_pp = bwt_str_bit0_pp ^ base_mask[c][0]; \
                bit1_cmp_pp = bwt_str_bit1_pp ^ base_mask[c][1]; \
                mismatch_mask_pp = bit0_cmp_pp | bit1_cmp_pp | my_cp_occ[occ_id_pp].dollar_mask[x]; \
                mismatch_mask_pp = mismatch_mask_pp >> (64 - y0_pp); \
                occ_pp += (y0_pp - _MM_COUNTBITS(mismatch_mask_pp)); \
                }\

typedef struct checkpoint_occ
{
    uint64_t bwt_str_bit0[ITER];
    uint64_t bwt_str_bit1[ITER];
    uint64_t dollar_mask[ITER];
    uint32_t cp_count[4];
    uint8_t  pad[PADDING];
}CP_OCC;

#endif
#else

#define CP_BLOCK_SIZE 32
#define CP_MASK 31
#define CP_SHIFT 5

#define \
GET_OCC(pp, c, c256, occ_id_pp, y_pp, occ_pp, bwt_str_pp, bwt_pp_vec, mask_pp_vec, mask_pp) \
                uint32_t occ_id_pp = pp >> CP_SHIFT; \
                uint32_t y_pp = pp & CP_MASK; \
                uint32_t occ_pp = my_cp_occ[occ_id_pp].cp_count[c]; \
                uint8_t *bwt_str_pp = my_cp_occ[occ_id_pp].bwt_str; \
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
#endif



// from index ele
#define BWA_IDX_BWT 0x1
#define BWA_IDX_BNS 0x2
#define BWA_IDX_PAC 0x4
#define BWA_IDX_ALL 0x7

typedef struct {
	bntseq_t *bns; // information on the reference sequences
	uint8_t  *pac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base
	int    is_shm;
	int64_t l_mem;
	uint8_t  *mem;
} bwaidx_fm_t;

typedef struct batch_meta
{
    int64_t k, l;
    int32_t rid;
    int32_t cur_ind;
}BatchMetadata;

typedef struct inexact_batch_meta
{
    int32_t rid;
    int32_t cur_fifo_id;
    int32_t fifo_size;
    int32_t num_matches;
}InexactBatchMetadata;

typedef struct align_state
{
    uint32_t k, l;
    int32_t i;
    int32_t z;
}AlignState;


typedef struct smem_struct
{
    uint32_t rid;
    uint32_t m, n;
    int64_t k, l, s;
    //int smem_id;
}SMEM;

#define SAL_PFD 16

// class FMI_search: public indexEle
class FMI_search
{
    public:
    FMI_search(const char *fname);
    ~FMI_search();
    //int64_t beCalls;
    
    int build_index(int for_only);
    int build_index_forward_only();
    int build_index_with_rev_complement();
    void load_index(int for_only);
    void load_index_forward_only();
    void load_index_with_rev_complement();

    void exact_search(uint8_t *enc_qdb,
                     int64_t num_queries,
                     uint32_t query_length,
                     int64_t *k_l_range,
                     int32_t batch_size);

    void inexact_search(uint8_t *enc_qdb,
                        int32_t z,
                        int64_t num_queries,
                        uint32_t query_length,
                        int64_t *k_l_range,
                        int32_t batch_size);

    void getSMEMs(uint8_t *enc_qdb,
                  int32_t numReads,
                  int32_t batch_size,
                  int32_t readlength,
                  int32_t minSeedLengh,
                  int32_t numthreads,
                  SMEM *matchArray,
                  int64_t *numTotalSmem);

     void getSMEMsSeedForwardShrink(uint8_t *enc_qdb,
                                         int16_t *query_pos_array,
                                         int32_t *min_intv_array,
                                         int32_t *rid_array,
                                         int32_t numReads,
                                         int32_t batch_size,
                                         const bseq1_t *seq_,
                                         int32_t *query_cum_len_ar,
                                         int32_t max_readlength,
                                         int32_t minSeedLen,
                                         SMEM *matchArray,
                                         int64_t *__numTotalSmem,
					 int *lisa_min_intv);
    
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
    void get_sa_entries(int64_t *posArray,
                        int64_t *coordArray,
                        uint32_t count,
                        int32_t nthreads);
    void get_sa_entries(SMEM *smemArray,
                        int64_t *coordArray,
                        int32_t *coordCountArray,
                        uint32_t count,
                        int32_t max_occ);
    // SA Compression
    int64_t get_sa_entry_compressed(int64_t pos, int tid);
    void get_sa_entries(SMEM *smemArray,
                        int64_t *coordArray,
                        int32_t *coordCountArray,
                        uint32_t count,
                        int32_t max_occ,
                        int tid);
    int64_t call_one_step(int64_t pos, int64_t &sa_entry, int64_t &offset);

    void get_sa_entries_prefetch(SMEM *smemArray, int64_t *coordArray,
                                 int64_t *coordCountArray, int64_t count,
                                 const int32_t max_occ, int tid, int64_t &id_);

    void bwa_idx_load_ele(const char *hint, int which);
    
    int64_t reference_seq_len;
    int64_t sentinel_index;
    bwaidx_fm_t *idx;    
        
    CP_OCC *cp_occ;
    int64_t count[5];
    SMEM backwardExt(SMEM smem, uint8_t a);
private:
        char file_name[PATH_MAX];
        int64_t index_alloc;
        uint32_t *sa_ls_word;
        int8_t *sa_ms_byte;

    uint64_t *one_hot_mask_array;

    // #if ((!__AVX2__))
    // BIT_DATA_TYPE base_mask[4][2];
    // #else
    // uint8_t *c_bcast_array;
    // #endif

        int64_t pac_seq_len(const char *fn_pac);
        void pac2nt(const char *fn_pac,
                    std::string &reference_seq,
                    int for_only);
        int build_fm_index_avx(const char *ref_file_name,
                               char *binary_seq,
                               int64_t ref_seq_len,
                               int64_t *sa_bwt,
                               int64_t *count);
        int build_fm_index_scalar(const char *ref_file_name,
                               char *binary_seq,
                               int64_t ref_seq_len,
                               int64_t *sa_bwt,
                               int64_t *count);
        int build_fm_index(const char *ref_file_name,
                           char *binary_seq,
                           int64_t ref_seq_len,
                           int64_t *sa_bwt,
                           int64_t *count);    
};

#endif
