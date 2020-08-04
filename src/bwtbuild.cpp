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
#include <inttypes.h>
#include <limits.h>

#include "sais.h"

#include "utils.h"
#include "bntseq.h"

#if !SAIS
#include<seqan/index.h>
using namespace seqan;
#endif

#include "memcpy_bwamem.h"

#ifdef __cplusplus
extern "C" {
#endif
#include "safe_str_lib.h"
#ifdef __cplusplus
}
#endif


#define assert_not_null(x, size) \
        if (x == NULL) { fprintf(stderr, "Allocation of %0.2lf GB for " #x " failed.\nCurrent Allocation = %0.2lf GB\n", size * 1.0 /(1024*1024*1024), cur_alloc * 1.0 /(1024*1024*1024)); exit(EXIT_FAILURE); }

int64_t cur_alloc = 0;
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
                fprintf(stderr, "ERROR! Value of nt is not in 0,1,2,3!");
                exit(EXIT_FAILURE);
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

void build_binaryRef(const char* prefix) {

    std::string reference_seq;
    char pac_file_name[PATH_MAX];
    strcpy_s(pac_file_name, PATH_MAX, prefix);
    strcat_s(pac_file_name, PATH_MAX, ".pac"); 
    pac2nt(pac_file_name, reference_seq);
    int64_t size = reference_seq.length() * sizeof(char);
    char *binary_ref_seq = (char *)_mm_malloc(size, 64);
    assert_not_null(binary_ref_seq, size);
    char binary_ref_name[PATH_MAX];
    strcpy_s(binary_ref_name, PATH_MAX, prefix);
    strcat_s(binary_ref_name, PATH_MAX, ".0123");
    std::fstream binary_ref_stream (binary_ref_name, std::ios::out | std::ios::binary);
    binary_ref_stream.seekg(0);	
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
    _mm_free(binary_ref_seq);

}
