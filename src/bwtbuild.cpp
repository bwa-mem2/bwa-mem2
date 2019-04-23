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

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include<cstring>
#include<vector>
#include<set>
#include <ctime>
#include<fstream>

#include "sais.h"

#include "utils.h"
#include "bntseq.h"

using namespace std;

#define DUMMY_CHAR 6

#if defined (__2BIT_LEVEL__)

#if (CP_BLOCK_SIZE == 32) || (CP_BLOCK_SIZE == 64)

#if (CP_BLOCK_SIZE == 64)

#define CP_MASK 63
#define CP_SHIFT 6
#define BIT_DATA_TYPE uint64_t
#define PADDING 24

#else
#define CP_MASK 31
#define CP_SHIFT 5
#define BIT_DATA_TYPE uint32_t
#define PADDING 4

#endif


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



typedef struct checkpoint_occ
{
    uint64_t bwt_str_bit0[ITER];
    uint64_t bwt_str_bit1[ITER];
    uint64_t dollar_mask[ITER];
    uint32_t cp_count[4];
}CP_OCC;

#elif (CP_BLOCK_SIZE == 256)

#define CP_MASK 255
#define CP_SHIFT 8
#define BIT_DATA_TYPE uint64_t
#define PADDING 16
#define ITER 4
#define _MM_COUNTBITS _mm_countbits_64

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
#define PADDING 16

typedef struct checkpoint_occ
{
    uint8_t  bwt_str[CP_BLOCK_SIZE];
    uint32_t cp_count[4];
    uint8_t  pad[PADDING];
}CP_OCC;

#endif

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

int build_fm_index(const char *ref_file_name, char *binary_seq, int64_t ref_seq_len, int64_t *sa_bwt, int64_t *count) {
    printf("ref_seq_len = %ld\n", ref_seq_len);
    fflush(stdout);

    char outname[200];
#ifdef __2BIT_LEVEL__
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

#pragma omp parallel for
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
            memcpy(cpo.bwt_str, bwt + i, CP_BLOCK_SIZE * sizeof(uint8_t));
            memset(cpo.pad, 0, 16);
            cp_occ[i >> CP_SHIFT] = cpo;
        }
        cp_count[bwt[i]]++;
    }
    outstream.write((char*)cp_occ, cp_occ_size * sizeof(CP_OCC));
    outstream.write((char*)sa_bwt, ref_seq_len * sizeof(int64_t));
    outstream.close();
    printf("max_occ_ind = %ld\n", i >> CP_SHIFT);    
    fflush(stdout);

    _mm_free(cp_occ);
    _mm_free(bwt);
    return 0;
}

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

    build_fm_index(prefix, binary_ref_seq, pac_len, suffix_array, count);
    _mm_free(binary_ref_seq);
    _mm_free(suffix_array);
    return 0;
}

