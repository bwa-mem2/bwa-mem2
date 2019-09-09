/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Vasimuddin Md, Sanchit Misra, Intel Corporation, Heng Li.

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

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include<cstring>
#include<vector>
#include<set>
#include <ctime>
#include<fstream>
#include <emmintrin.h>

#include "sais.h"

#include "utils.h"
#include "bntseq.h"

#if !SAIS
#include<seqan/index.h>
using namespace seqan;
#endif

using namespace std;

#define DUMMY_CHAR 6

// #if ((!__AVX2__))
// SSE stuff
#define CP_BLOCK_SIZE_SSE 64
#define CP_MASK_SSE 63
#define CP_SHIFT_SSE 6
#define BIT_DATA_TYPE uint64_t
#define PADDING_SSE 24

typedef struct checkpoint_occ_sse
{
    BIT_DATA_TYPE bwt_str_bit0;
    BIT_DATA_TYPE bwt_str_bit1;
    BIT_DATA_TYPE dollar_mask;
    uint32_t cp_count[4];
    uint8_t  pad[PADDING_SSE];
}CP_OCC_SSE;

// #else
// AVX stuff
#define CP_BLOCK_SIZE_AVX 32
#define CP_MASK_AVX 31
#define CP_SHIFT_AVX 5
#define PADDING_AVX 16

typedef struct checkpoint_occ_avx
{
    uint8_t  bwt_str[CP_BLOCK_SIZE_AVX];
    uint32_t cp_count[4];
    uint8_t  pad[PADDING_AVX];
}CP_OCC_AVX;

// #endif

int64_t pac_seq_len(const char *fn_pac)
{
	FILE *fp;
	int64_t pac_len;
	uint8_t c;
	fp = xopen(fn_pac, "rb");
	err_fseek(fp, -1, SEEK_END);
	pac_len = err_ftell(fp);
	err_fread_noeof(&c, 1, 1, fp);
	err_fclose(fp);
	return (pac_len - 1) * 4 + (int)c;
}

void pac2nt(const char *fn_pac, std::string &reference_seq)
{
	uint8_t *buf2;
	int64_t i, pac_size, seq_len;
	FILE *fp;

	// initialization
	seq_len = pac_seq_len(fn_pac);
	fp = xopen(fn_pac, "rb");

	// prepare sequence
	pac_size = (seq_len>>2) + ((seq_len&3) == 0? 0 : 1);
	buf2 = (uint8_t*)calloc(pac_size, 1);
	err_fread_noeof(buf2, 1, pac_size, fp);
	err_fclose(fp);
	for (i = 0; i < seq_len; ++i) {
		int nt = buf2[i>>2] >> ((3 - (i&3)) << 1) & 3;
        switch(nt)
        {
            case 0:
                reference_seq += "A";
            break;
            case 1:
                reference_seq += "C";
            break;
            case 2:
                reference_seq += "G";
            break;
            case 3:
                reference_seq += "T";
            break;
            default:
                printf("ERROR! Value of nt is not in 0,1,2,3!");
                exit(0);
        }
	}
    for(i = seq_len - 1; i >= 0; i--)
    {
        char c = reference_seq[i];
        switch(c)
        {
            case 'A':
                reference_seq += "T";
            break;
            case 'C':
                reference_seq += "G";
            break;
            case 'G':
                reference_seq += "C";
            break;
            case 'T':
                reference_seq += "A";
            break;
        }
    }
	free(buf2);
}

#if LEGACY
int build_fm_index_generic(const char *ref_file_name, char *binary_seq, int64_t ref_seq_len, int64_t *sa_bwt, int64_t *count) {
    printf("ref_seq_len = %ld\n", ref_seq_len);
    fflush(stdout);

    char outname[200];
#if ((!__AVX2__))
    sprintf(outname, "%s.bwt.2bit.%d", ref_file_name, CP_BLOCK_SIZE);
#else
    sprintf(outname, "%s.bwt.8bit.%d", ref_file_name, CP_BLOCK_SIZE);
#endif
    std::fstream outstream (outname, ios::out | ios::binary);
    outstream.seekg(0);	

    printf("count = %ld, %ld, %ld, %ld, %ld\n", count[0], count[1], count[2], count[3], count[4]);
    fflush(stdout);

    uint8_t *bwt;

    ref_seq_len++;
    outstream.write((char *)(&ref_seq_len), 1 * sizeof(int64_t));
    outstream.write((char*)count, 5 * sizeof(int64_t));

    int64_t i;
    int64_t ref_seq_len_aligned = ((ref_seq_len + CP_BLOCK_SIZE - 1) / CP_BLOCK_SIZE) * CP_BLOCK_SIZE;
    bwt = (uint8_t *)_mm_malloc(ref_seq_len_aligned * sizeof(uint8_t), 64);

// #pragma omp parallel for
    for(i=0; i< ref_seq_len; i++)
    {
        if(sa_bwt[i] == 0)
        {
            bwt[i] = 4;
            printf("BWT[%ld] = 4\n", i);
        }
        else
        {
            char c = binary_seq[sa_bwt[i]-1];
            switch(c)
            {
                case 0: bwt[i] = 0;
                          break;
                case 1: bwt[i] = 1;
                          break;
                case 2: bwt[i] = 2;
                          break;
                case 3: bwt[i] = 3;
                          break;
                default:
                        printf("ERROR! i = %ld, c = %c\n", i, c);
                        exit(1);
            }
        }
    }
    for(i = ref_seq_len; i < ref_seq_len_aligned; i++)
        bwt[i] = DUMMY_CHAR;


    printf("CP_SHIFT = %d, CP_MASK = %d\n", CP_SHIFT, CP_MASK);
    printf("sizeof CP_OCC = %ld\n", sizeof(CP_OCC));
    fflush(stdout);
    // create checkpointed occ
    int64_t cp_occ_size = (ref_seq_len >> CP_SHIFT) + 1;
    CP_OCC *cp_occ = NULL;

    cp_occ = (CP_OCC *)_mm_malloc(cp_occ_size * sizeof(CP_OCC), 64);
    memset(cp_occ, 0, cp_occ_size * sizeof(CP_OCC));
    uint32_t cp_count[16];

    memset(cp_count, 0, 16 * sizeof(uint32_t));
    for(i = 0; i < ref_seq_len; i++)
    {
        if((i & CP_MASK) == 0)
        {
            CP_OCC cpo;
            cpo.cp_count[0] = cp_count[0];
            cpo.cp_count[1] = cp_count[1];
            cpo.cp_count[2] = cp_count[2];
            cpo.cp_count[3] = cp_count[3];
#if ((!__AVX2__))
			BIT_DATA_TYPE bwt_str_bit0 = 0;
			BIT_DATA_TYPE bwt_str_bit1 = 0;
			BIT_DATA_TYPE dollar_mask = 0;
			int32_t j;
			for(j = 0; j < CP_BLOCK_SIZE; j++)
			{
				uint8_t c = bwt[i + j];
				if((c == 4) || (c == DUMMY_CHAR))
				{
					dollar_mask <<= 1;
					dollar_mask += 1;
					c = 0;
				}
				else if(c > 3)
				{
					printf("ERROR! [%ld, %d] c = %u\n", (long)i, j, c);
					exit(0);
				}
				else
				{
					dollar_mask <<= 1;
					dollar_mask += 0;
				}
				bwt_str_bit0 = bwt_str_bit0 << 1;
				bwt_str_bit0 += (c & 1);
				bwt_str_bit1 = bwt_str_bit1 << 1;
				bwt_str_bit1 += ((c >> 1) & 1);
			}
			cpo.bwt_str_bit0 = bwt_str_bit0;
			cpo.bwt_str_bit1 = bwt_str_bit1;
			cpo.dollar_mask  = dollar_mask;
#else
			memcpy(cpo.bwt_str, bwt + i, CP_BLOCK_SIZE * sizeof(uint8_t));
#endif
            memset(cpo.pad, 0, PADDING);
            cp_occ[i >> CP_SHIFT] = cpo;
        }
        cp_count[bwt[i]]++;
    }
    outstream.write((char*)cp_occ, cp_occ_size * sizeof(CP_OCC));


    uint32_t *sa_ls_word = (uint32_t *)_mm_malloc(ref_seq_len * sizeof(uint32_t), 64);
    int8_t *sa_ms_byte = (int8_t *)_mm_malloc(ref_seq_len * sizeof(int8_t), 64);
    for(i = 0; i < ref_seq_len; i++)
    {
        sa_ls_word[i] = sa_bwt[i] & 0xffffffff;
        sa_ms_byte[i] = (sa_bwt[i] >> 32) & 0xff;
    }
    outstream.write((char*)sa_ms_byte, ref_seq_len * sizeof(int8_t));
    outstream.write((char*)sa_ls_word, ref_seq_len * sizeof(uint32_t));
    outstream.close();
    printf("max_occ_ind = %ld\n", i >> CP_SHIFT);    
    fflush(stdout);

    _mm_free(cp_occ);
    _mm_free(bwt);
    _mm_free(sa_ms_byte);
    _mm_free(sa_ls_word);
    return 0;
}
#endif

int build_fm_index_avx(const char *ref_file_name, char *binary_seq, int64_t ref_seq_len, int64_t *sa_bwt, int64_t *count) {
    printf("ref_seq_len = %ld\n", ref_seq_len);
    fflush(stdout);

    char outname[200];
    sprintf(outname, "%s.bwt.8bit.%d", ref_file_name, CP_BLOCK_SIZE_AVX);

    std::fstream outstream (outname, ios::out | ios::binary);
    outstream.seekg(0);	

    printf("count = %ld, %ld, %ld, %ld, %ld\n", count[0], count[1], count[2], count[3], count[4]);
    fflush(stdout);

    uint8_t *bwt;

    ref_seq_len++;
    outstream.write((char *)(&ref_seq_len), 1 * sizeof(int64_t));
    outstream.write((char*)count, 5 * sizeof(int64_t));

    int64_t i;
    int64_t ref_seq_len_aligned = ((ref_seq_len + CP_BLOCK_SIZE_AVX - 1) / CP_BLOCK_SIZE_AVX) * CP_BLOCK_SIZE_AVX;
    bwt = (uint8_t *)_mm_malloc(ref_seq_len_aligned * sizeof(uint8_t), 64);

// #pragma omp parallel for
    for(i=0; i< ref_seq_len; i++)
    {
        if(sa_bwt[i] == 0)
        {
            bwt[i] = 4;
            printf("BWT[%ld] = 4\n", i);
        }
        else
        {
            char c = binary_seq[sa_bwt[i]-1];
            switch(c)
            {
                case 0: bwt[i] = 0;
                          break;
                case 1: bwt[i] = 1;
                          break;
                case 2: bwt[i] = 2;
                          break;
                case 3: bwt[i] = 3;
                          break;
                default:
                        printf("ERROR! i = %ld, c = %c\n", i, c);
                        exit(1);
            }
        }
    }
    for(i = ref_seq_len; i < ref_seq_len_aligned; i++)
        bwt[i] = DUMMY_CHAR;


    printf("CP_SHIFT = %d, CP_MASK = %d\n", CP_SHIFT_AVX, CP_MASK_AVX);
    printf("sizeof CP_OCC = %ld\n", sizeof(CP_OCC_AVX));
    fflush(stdout);
    // create checkpointed occ
    int64_t cp_occ_size = (ref_seq_len >> CP_SHIFT_AVX) + 1;
    CP_OCC_AVX *cp_occ = NULL;

    cp_occ = (CP_OCC_AVX *)_mm_malloc(cp_occ_size * sizeof(CP_OCC_AVX), 64);
    memset(cp_occ, 0, cp_occ_size * sizeof(CP_OCC_AVX));
    uint32_t cp_count[16];

    memset(cp_count, 0, 16 * sizeof(uint32_t));
    for(i = 0; i < ref_seq_len; i++)
    {
        if((i & CP_MASK_AVX) == 0)
        {
            CP_OCC_AVX cpo;
            cpo.cp_count[0] = cp_count[0];
            cpo.cp_count[1] = cp_count[1];
            cpo.cp_count[2] = cp_count[2];
            cpo.cp_count[3] = cp_count[3];
			memcpy(cpo.bwt_str, bwt + i, CP_BLOCK_SIZE_AVX * sizeof(uint8_t));

            memset(cpo.pad, 0, PADDING_AVX);
            cp_occ[i >> CP_SHIFT_AVX] = cpo;
        }
        cp_count[bwt[i]]++;
    }
    outstream.write((char*)cp_occ, cp_occ_size * sizeof(CP_OCC_AVX));


    uint32_t *sa_ls_word = (uint32_t *)_mm_malloc(ref_seq_len * sizeof(uint32_t), 64);
    int8_t *sa_ms_byte = (int8_t *)_mm_malloc(ref_seq_len * sizeof(int8_t), 64);
    for(i = 0; i < ref_seq_len; i++)
    {
        sa_ls_word[i] = sa_bwt[i] & 0xffffffff;
        sa_ms_byte[i] = (sa_bwt[i] >> 32) & 0xff;
    }
    outstream.write((char*)sa_ms_byte, ref_seq_len * sizeof(int8_t));
    outstream.write((char*)sa_ls_word, ref_seq_len * sizeof(uint32_t));
    outstream.close();
    printf("max_occ_ind = %ld\n", i >> CP_SHIFT_AVX);    
    fflush(stdout);

    _mm_free(cp_occ);
    _mm_free(bwt);
    _mm_free(sa_ms_byte);
    _mm_free(sa_ls_word);
    return 0;
}

int build_fm_index_sse(const char *ref_file_name, char *binary_seq, int64_t ref_seq_len, int64_t *sa_bwt, int64_t *count) {
    printf("ref_seq_len = %ld\n", ref_seq_len);
    fflush(stdout);

    char outname[200];

    sprintf(outname, "%s.bwt.2bit.%d", ref_file_name, CP_BLOCK_SIZE_SSE);

    std::fstream outstream (outname, ios::out | ios::binary);
    outstream.seekg(0);	

    printf("count = %ld, %ld, %ld, %ld, %ld\n", count[0], count[1], count[2], count[3], count[4]);
    fflush(stdout);

    uint8_t *bwt;

    ref_seq_len++;
    outstream.write((char *)(&ref_seq_len), 1 * sizeof(int64_t));
    outstream.write((char*)count, 5 * sizeof(int64_t));

    int64_t i;
    int64_t ref_seq_len_aligned = ((ref_seq_len + CP_BLOCK_SIZE_SSE - 1) / CP_BLOCK_SIZE_SSE) * CP_BLOCK_SIZE_SSE;
    bwt = (uint8_t *)_mm_malloc(ref_seq_len_aligned * sizeof(uint8_t), 64);

// #pragma omp parallel for
    for(i=0; i< ref_seq_len; i++)
    {
        if(sa_bwt[i] == 0)
        {
            bwt[i] = 4;
            printf("BWT[%ld] = 4\n", i);
        }
        else
        {
            char c = binary_seq[sa_bwt[i]-1];
            switch(c)
            {
                case 0: bwt[i] = 0;
                          break;
                case 1: bwt[i] = 1;
                          break;
                case 2: bwt[i] = 2;
                          break;
                case 3: bwt[i] = 3;
                          break;
                default:
                        printf("ERROR! i = %ld, c = %c\n", i, c);
                        exit(1);
            }
        }
    }
    for(i = ref_seq_len; i < ref_seq_len_aligned; i++)
        bwt[i] = DUMMY_CHAR;


    printf("CP_SHIFT = %d, CP_MASK = %d\n", CP_SHIFT_SSE, CP_MASK_SSE);
    printf("sizeof CP_OCC = %ld\n", sizeof(CP_OCC_SSE));
    fflush(stdout);
    // create checkpointed occ
    int64_t cp_occ_size = (ref_seq_len >> CP_SHIFT_SSE) + 1;
    CP_OCC_SSE *cp_occ = NULL;

    cp_occ = (CP_OCC_SSE *)_mm_malloc(cp_occ_size * sizeof(CP_OCC_SSE), 64);
    memset(cp_occ, 0, cp_occ_size * sizeof(CP_OCC_SSE));
    uint32_t cp_count[16];

    memset(cp_count, 0, 16 * sizeof(uint32_t));
    for(i = 0; i < ref_seq_len; i++)
    {
        if((i & CP_MASK_SSE) == 0)
        {
            CP_OCC_SSE cpo;
            cpo.cp_count[0] = cp_count[0];
            cpo.cp_count[1] = cp_count[1];
            cpo.cp_count[2] = cp_count[2];
            cpo.cp_count[3] = cp_count[3];

			BIT_DATA_TYPE bwt_str_bit0 = 0;
			BIT_DATA_TYPE bwt_str_bit1 = 0;
			BIT_DATA_TYPE dollar_mask = 0;
			int32_t j;
			for(j = 0; j < CP_BLOCK_SIZE_SSE; j++)
			{
				uint8_t c = bwt[i + j];
				if((c == 4) || (c == DUMMY_CHAR))
				{
					dollar_mask <<= 1;
					dollar_mask += 1;
					c = 0;
				}
				else if(c > 3)
				{
					printf("ERROR! [%ld, %d] c = %u\n", (long)i, j, c);
					exit(0);
				}
				else
				{
					dollar_mask <<= 1;
					dollar_mask += 0;
				}
				bwt_str_bit0 = bwt_str_bit0 << 1;
				bwt_str_bit0 += (c & 1);
				bwt_str_bit1 = bwt_str_bit1 << 1;
				bwt_str_bit1 += ((c >> 1) & 1);
			}
			cpo.bwt_str_bit0 = bwt_str_bit0;
			cpo.bwt_str_bit1 = bwt_str_bit1;
			cpo.dollar_mask  = dollar_mask;

            memset(cpo.pad, 0, PADDING_SSE);
            cp_occ[i >> CP_SHIFT_SSE] = cpo;
        }
        cp_count[bwt[i]]++;
    }
    outstream.write((char*)cp_occ, cp_occ_size * sizeof(CP_OCC_SSE));


    uint32_t *sa_ls_word = (uint32_t *)_mm_malloc(ref_seq_len * sizeof(uint32_t), 64);
    int8_t *sa_ms_byte = (int8_t *)_mm_malloc(ref_seq_len * sizeof(int8_t), 64);
    for(i = 0; i < ref_seq_len; i++)
    {
        sa_ls_word[i] = sa_bwt[i] & 0xffffffff;
        sa_ms_byte[i] = (sa_bwt[i] >> 32) & 0xff;
    }
    outstream.write((char*)sa_ms_byte, ref_seq_len * sizeof(int8_t));
    outstream.write((char*)sa_ls_word, ref_seq_len * sizeof(uint32_t));
    outstream.close();
    printf("max_occ_ind = %ld\n", i >> CP_SHIFT_SSE);    
    fflush(stdout);

    _mm_free(cp_occ);
    _mm_free(bwt);
    _mm_free(sa_ms_byte);
    _mm_free(sa_ls_word);
    return 0;
}

#if !SAIS
// native ubild index routines, needs seqAn library
int build_sa(std::string ref_db, String<int64_t> sa, int64_t *count, int64_t *suffix_array)
{
    int64_t startTick;
    startTick = __rdtsc();
    String<char> text;
    text=ref_db;
    int64_t j;
    int64_t text_length = length(text);
    resize(sa, text_length);
    printf("SA init ticks = %ld\n", __rdtsc() - startTick);
    startTick = __rdtsc();
    createSuffixArray(sa, text, Skew7());
    printf("SA ticks = %ld\n", __rdtsc() - startTick);
    startTick = __rdtsc();
    for(int64_t i = 0; i < text_length; i++)
    {
        suffix_array[i]=sa[i];

    }
    printf("SA copy ticks = %ld\n", __rdtsc() - startTick);
    startTick = __rdtsc();

    int64_t ref_seq_len = length(text) - 1;
    printf("ref_seq_len = %ld\n", ref_seq_len);
    memset(count, 0, 16 * sizeof(int64_t));
    for (int64_t i=0; i < ref_seq_len; i++)  
    {
        switch(text[i])
        {
            case 'A': ++count[0];
                      break;
            case 'C': ++count[1];
                      break;
            case 'G': ++count[2];
                      break;
            case 'T': ++count[3];
                      break;
            default: break;
        }
    }	
    count[4]=count[0]+count[1]+count[2]+count[3];
    count[3]=count[0]+count[1]+count[2];
    count[2]=count[0]+count[1];
    count[1]=count[0];
    count[0]=0;	
    printf("count ticks = %ld\n", __rdtsc() - startTick);
    startTick = __rdtsc();
    printf("count = %ld, %ld, %ld, %ld, %ld\n", count[0], count[1], count[2], count[3], count[4]);

    return 0;
}

int build_index(const char *prefix) {

    int64_t startTick;
    startTick = __rdtsc();
       
    std::string reference_seq;
    char pac_file_name[200];
    sprintf(pac_file_name, "%s.pac", prefix);
    pac2nt(pac_file_name, reference_seq);
    int status;
    char *binary_ref_seq = (char *)_mm_malloc((reference_seq.length())*sizeof(char), 64);
    char binary_ref_name[200];
    sprintf(binary_ref_name, "%s.0123", prefix);
    std::fstream binary_ref_stream (binary_ref_name, ios::out | ios::binary);
    binary_ref_stream.seekg(0);	
    printf("init ticks = %ld\n", __rdtsc() - startTick);
    startTick = __rdtsc();
    int64_t i;
    for(i = 0; i < reference_seq.length(); i++)
    {
        switch(reference_seq[i])
        {
            case 'A':
            binary_ref_seq[i] = 0;
            break;
            case 'C':
            binary_ref_seq[i] = 1;
            break;
            case 'G':
            binary_ref_seq[i] = 2;
            break;
            case 'T':
            binary_ref_seq[i] = 3;
            break;
            default:
            binary_ref_seq[i] = 4;

        }
    }
    printf("ref seq len = %ld\n", reference_seq.length());
    binary_ref_stream.write(binary_ref_seq, reference_seq.length() * sizeof(char));
    printf("binary seq ticks = %ld\n", __rdtsc() - startTick);
    startTick = __rdtsc();

    reference_seq+='$';
    int64_t count[16];
    String<int64_t> sa;
    String<char> bwt;

    int64_t *suffix_array=(int64_t *)_mm_malloc(length(reference_seq) * sizeof(int64_t), 64);
    startTick = __rdtsc();
    status=build_sa(reference_seq, sa, count, suffix_array);
    printf("build index ticks = %ld\n", __rdtsc() - startTick);
    startTick = __rdtsc();

    build_fm_index_avx(prefix, binary_ref_seq, length(reference_seq) - 1, suffix_array, count);
	build_fm_index_sse(prefix, binary_ref_seq, length(reference_seq) - 1, suffix_array, count);
    _mm_free(binary_ref_seq);
    _mm_free(suffix_array);
    return 0;
}

#else

int build_index(const char *prefix) {

    int64_t startTick;
    startTick = __rdtsc();

    std::string reference_seq;
    char pac_file_name[200];
    sprintf(pac_file_name, "%s.pac", prefix);
    pac2nt(pac_file_name, reference_seq);
	int64_t pac_len = reference_seq.length();
    int status;
    char *binary_ref_seq = (char *)_mm_malloc(pac_len * sizeof(char), 64);
    char binary_ref_name[200];
    sprintf(binary_ref_name, "%s.0123", prefix);
    std::fstream binary_ref_stream (binary_ref_name, ios::out | ios::binary);
    binary_ref_stream.seekg(0);
    fprintf(stderr, "init ticks = %ld\n", __rdtsc() - startTick);
    startTick = __rdtsc();
    int64_t i, count[16];
	memset(count, 0, sizeof(int64_t) * 16);
    for(i = 0; i < pac_len; i++)
    {
        switch(reference_seq[i])
        {
            case 'A':
            binary_ref_seq[i] = 0, ++count[0];
            break;
            case 'C':
            binary_ref_seq[i] = 1, ++count[1];
            break;
            case 'G':
            binary_ref_seq[i] = 2, ++count[2];
            break;
            case 'T':
            binary_ref_seq[i] = 3, ++count[3];
            break;
            default:
            binary_ref_seq[i] = 4;

        }
    }
    count[4]=count[0]+count[1]+count[2]+count[3];
    count[3]=count[0]+count[1]+count[2];
    count[2]=count[0]+count[1];
    count[1]=count[0];
    count[0]=0;
    fprintf(stderr, "ref seq len = %ld\n", pac_len);
    binary_ref_stream.write(binary_ref_seq, pac_len * sizeof(char));
    fprintf(stderr, "binary seq ticks = %ld\n", __rdtsc() - startTick);
    startTick = __rdtsc();

    int64_t *suffix_array=(int64_t *)_mm_malloc((pac_len + 2) * sizeof(int64_t), 64);
    startTick = __rdtsc();
	status = saisxx(reference_seq.c_str(), suffix_array + 1, pac_len);
	suffix_array[0] = pac_len;
    fprintf(stderr, "build index ticks = %ld\n", __rdtsc() - startTick);
    startTick = __rdtsc();

    build_fm_index_avx(prefix, binary_ref_seq, pac_len, suffix_array, count);
	build_fm_index_sse(prefix, binary_ref_seq, pac_len, suffix_array, count);
    _mm_free(binary_ref_seq);
    _mm_free(suffix_array);
    return 0;
}

#endif
