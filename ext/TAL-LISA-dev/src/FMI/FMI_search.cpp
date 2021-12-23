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
#include <bits/stdc++.h>
using namespace std;
#include <stdio.h>
#include <cstring>
#include "sais.h"
#include "FMI_search.h"
#ifdef __cplusplus
extern "C" {
#endif
#include "safe_str_lib.h"
#ifdef __cplusplus
}
#endif

FMI_search::FMI_search(const char *fname)
{
    fprintf(stderr, "* Entering FMI_search\n");
    // assert(strlen(fname) < PATH_MAX);
    strcpy_s(file_name, PATH_MAX, fname);
    reference_seq_len = 0;
    sentinel_index = 0;
    index_alloc = 0;
    memset(count, 0, 5 * sizeof(int64_t));
    sa_ls_word = NULL;
    sa_ms_byte = NULL;
    cp_occ = NULL;
    one_hot_mask_array = NULL;

    // read_index_ele
    idx = (bwaidx_fm_t*) calloc(1, sizeof(bwaidx_fm_t));
    assert(idx != NULL);    
}

FMI_search::~FMI_search()
{
    if(sa_ms_byte)
        _mm_free(sa_ms_byte);
    if(sa_ls_word)
        _mm_free(sa_ls_word);
    if(cp_occ)
        _mm_free(cp_occ);    
    if(one_hot_mask_array)
        _mm_free(one_hot_mask_array);

    // read_index_ele
    if (idx == 0) return;
    if (idx->mem == 0)
    {
        if (idx->bns) bns_destroy(idx->bns);
        if (idx->pac) free(idx->pac);
    } else {
        free(idx->bns->anns); free(idx->bns);
        if (!idx->is_shm) free(idx->mem);
    }
    free(idx);    
}

void FMI_search::bwa_idx_load_ele(const char *hint, int which)
{
    char *prefix;
    int l_hint = strlen(hint);
    prefix = (char *) malloc(l_hint + 3 + 4 + 1);
    assert(prefix != NULL);
    strcpy_s(prefix, l_hint + 3 + 4 + 1, hint);

    fprintf(stderr, "* Index prefix: %s\n", prefix);
    
    // idx = (bwaidx_fm_t*) calloc(1, sizeof(bwaidx_fm_t));
    if (which & BWA_IDX_BNS) {
        int i, c;
        idx->bns = bns_restore(prefix);
        if (idx->bns == 0) {
            printf("Error!! : [%s] bns is NULL!!\n", __func__);
            exit(EXIT_FAILURE);
        }
        for (i = c = 0; i < idx->bns->n_seqs; ++i)
            if (idx->bns->anns[i].is_alt) ++c;
        
        fprintf(stderr, "* Read %d ALT contigs\n", c);
        
        if (which & BWA_IDX_PAC)
        {
            idx->pac = (uint8_t*) calloc(idx->bns->l_pac/4+1, 1);
            assert(idx->pac != NULL);
            err_fread_noeof(idx->pac, 1, idx->bns->l_pac/4+1, idx->bns->fp_pac); // concatenated 2-bit encoded sequence
            err_fclose(idx->bns->fp_pac);
            idx->bns->fp_pac = 0;
        }
    }
    free(prefix);
}

int64_t FMI_search::pac_seq_len(const char *fn_pac)
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

void FMI_search::pac2nt(const char *fn_pac, std::string &reference_seq, int for_only)
{
    uint8_t *buf2;
    int64_t i, pac_size, seq_len;
    FILE *fp;

    // initialization
    seq_len = pac_seq_len(fn_pac);
    assert(seq_len > 0);
    assert(seq_len <= 0x7fffffffffL);
    fp = xopen(fn_pac, "rb");

    // prepare sequence
    pac_size = (seq_len>>2) + ((seq_len&3) == 0? 0 : 1);
    buf2 = (uint8_t*)calloc(pac_size, 1);
    assert(buf2 != NULL);
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
    if(!for_only)
    {
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
    }
    free(buf2);
}


#if OLD
int FMI_search::build_fm_index_avx(const char *ref_file_name, char *binary_seq, int64_t ref_seq_len, int64_t *sa_bwt, int64_t *count) {
    printf("ref_seq_len = %ld\n", ref_seq_len);
    fflush(stdout);

    char outname[PATH_MAX];
    strcpy_s(outname, PATH_MAX, ref_file_name);
    strcat_s(outname, PATH_MAX, CP_FILENAME_SUFFIX_AVX);
    //sprintf(outname, "%s.bwt.8bit.%d", ref_file_name, CP_BLOCK_SIZE_AVX);

    std::fstream outstream (outname, std::ios::out | std::ios::binary);
    outstream.seekg(0); 

    printf("count = %ld, %ld, %ld, %ld, %ld\n", count[0], count[1], count[2], count[3], count[4]);
    fflush(stdout);

    uint8_t *bwt;

    ref_seq_len++;
    outstream.write((char *)(&ref_seq_len), 1 * sizeof(int64_t));
    outstream.write((char*)count, 5 * sizeof(int64_t));

    int64_t i;
    int64_t ref_seq_len_aligned = ((ref_seq_len + CP_BLOCK_SIZE_AVX - 1) / CP_BLOCK_SIZE_AVX) * CP_BLOCK_SIZE_AVX;
    int64_t size = ref_seq_len_aligned * sizeof(uint8_t);
    bwt = (uint8_t *)_mm_malloc(size, 64);
    index_alloc += size;
    assert_not_null(bwt, size, index_alloc);

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
                        fprintf(stderr, "ERROR! i = %ld, c = %c\n", i, c);
                        exit(EXIT_FAILURE);
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

    size = cp_occ_size * sizeof(CP_OCC_AVX);
    cp_occ = (CP_OCC_AVX *)_mm_malloc(size, 64);
    index_alloc += size;
    assert_not_null(cp_occ, size, index_alloc);
    memset(cp_occ, 0, cp_occ_size * sizeof(CP_OCC_AVX));
    int64_t cp_count[16];

    memset(cp_count, 0, 16 * sizeof(int64_t));
    for(i = 0; i < ref_seq_len; i++)
    {
        if((i & CP_MASK_AVX) == 0)
        {
            CP_OCC_AVX cpo;
            cpo.cp_count[0] = cp_count[0];
            cpo.cp_count[1] = cp_count[1];
            cpo.cp_count[2] = cp_count[2];
            cpo.cp_count[3] = cp_count[3];
            memcpy_s(cpo.bwt_str, CP_BLOCK_SIZE_AVX * sizeof(uint8_t), bwt + i, CP_BLOCK_SIZE_AVX * sizeof(uint8_t));

            cp_occ[i >> CP_SHIFT_AVX] = cpo;
        }
        cp_count[bwt[i]]++;
    }
    outstream.write((char*)cp_occ, cp_occ_size * sizeof(CP_OCC_AVX));
    _mm_free(cp_occ);
    _mm_free(bwt);
    index_alloc -= (ref_seq_len_aligned * sizeof(uint8_t) + cp_occ_size * sizeof(CP_OCC_AVX));

    size = ref_seq_len * sizeof(uint32_t);
    uint32_t *sa_ls_word = (uint32_t *)_mm_malloc(size, 64);
    index_alloc += size;
    assert_not_null(sa_ls_word, size, index_alloc);
    size = ref_seq_len * sizeof(int8_t);
    int8_t *sa_ms_byte = (int8_t *)_mm_malloc(size, 64);
    index_alloc += size;
    assert_not_null(sa_ms_byte, size, index_alloc);
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

    _mm_free(sa_ms_byte);
    _mm_free(sa_ls_word);
    return 0;
}

int FMI_search::build_fm_index_scalar(const char *ref_file_name, char *binary_seq, int64_t ref_seq_len, int64_t *sa_bwt, int64_t *count) {
    printf("ref_seq_len = %ld\n", ref_seq_len);
    fflush(stdout);

    char outname[PATH_MAX];

    strcpy_s(outname, PATH_MAX, ref_file_name);
    strcat_s(outname, PATH_MAX, CP_FILENAME_SUFFIX_SCALAR);
    //sprintf(outname, "%s.bwt.2bit.%d", ref_file_name, CP_BLOCK_SIZE_SCALAR);

    std::fstream outstream (outname, std::ios::out | std::ios::binary);
    outstream.seekg(0); 

    printf("count = %ld, %ld, %ld, %ld, %ld\n", count[0], count[1], count[2], count[3], count[4]);
    fflush(stdout);

    uint8_t *bwt;

    ref_seq_len++;
    outstream.write((char *)(&ref_seq_len), 1 * sizeof(int64_t));
    outstream.write((char*)count, 5 * sizeof(int64_t));

    int64_t i;
    int64_t ref_seq_len_aligned = ((ref_seq_len + CP_BLOCK_SIZE_SCALAR - 1) / CP_BLOCK_SIZE_SCALAR) * CP_BLOCK_SIZE_SCALAR;
    int64_t size = ref_seq_len_aligned * sizeof(uint8_t);
    bwt = (uint8_t *)_mm_malloc(size, 64);
    assert_not_null(bwt, size, index_alloc);

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
                        fprintf(stderr, "ERROR! i = %ld, c = %c\n", i, c);
                        exit(EXIT_FAILURE);
            }
        }
    }
    for(i = ref_seq_len; i < ref_seq_len_aligned; i++)
        bwt[i] = DUMMY_CHAR;


    printf("CP_SHIFT = %d, CP_MASK = %d\n", CP_SHIFT_SCALAR, CP_MASK_SCALAR);
    printf("sizeof CP_OCC = %ld\n", sizeof(CP_OCC_SCALAR));
    fflush(stdout);
    // create checkpointed occ
    int64_t cp_occ_size = (ref_seq_len >> CP_SHIFT_SCALAR) + 1;
    CP_OCC_SCALAR *cp_occ = NULL;

    size = cp_occ_size * sizeof(CP_OCC_SCALAR);
    cp_occ = (CP_OCC_SCALAR *)_mm_malloc(size, 64);
    assert_not_null(cp_occ, size, index_alloc);
    memset(cp_occ, 0, cp_occ_size * sizeof(CP_OCC_SCALAR));
    int64_t cp_count[16];

    memset(cp_count, 0, 16 * sizeof(int64_t));
    for(i = 0; i < ref_seq_len; i++)
    {
        if((i & CP_MASK_SCALAR) == 0)
        {
            CP_OCC_SCALAR cpo;
            cpo.cp_count[0] = cp_count[0];
            cpo.cp_count[1] = cp_count[1];
            cpo.cp_count[2] = cp_count[2];
            cpo.cp_count[3] = cp_count[3];

            BIT_DATA_TYPE bwt_str_bit0 = 0;
            BIT_DATA_TYPE bwt_str_bit1 = 0;
            BIT_DATA_TYPE dollar_mask = 0;
            int32_t j;
            for(j = 0; j < CP_BLOCK_SIZE_SCALAR; j++)
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
                    fprintf(stderr, "ERROR! [%ld, %d] c = %u\n", (long)i, j, c);
                    exit(EXIT_FAILURE);
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

            memset(cpo.pad, 0, PADDING_SCALAR);
            cp_occ[i >> CP_SHIFT_SCALAR] = cpo;
        }
        cp_count[bwt[i]]++;
    }
    outstream.write((char*)cp_occ, cp_occ_size * sizeof(CP_OCC_SCALAR));
    _mm_free(cp_occ);
    _mm_free(bwt);

    size = ref_seq_len * sizeof(uint32_t);
    uint32_t *sa_ls_word = (uint32_t *)_mm_malloc(size, 64);
    assert_not_null(sa_ls_word, size, index_alloc);
    size = ref_seq_len * sizeof(int8_t);
    int8_t *sa_ms_byte = (int8_t *)_mm_malloc(size, 64);
    assert_not_null(sa_ms_byte, size, index_alloc);
    for(i = 0; i < ref_seq_len; i++)
    {
        sa_ls_word[i] = sa_bwt[i] & 0xffffffff;
        sa_ms_byte[i] = (sa_bwt[i] >> 32) & 0xff;
    }
    outstream.write((char*)sa_ms_byte, ref_seq_len * sizeof(int8_t));
    outstream.write((char*)sa_ls_word, ref_seq_len * sizeof(uint32_t));
    outstream.close();
    printf("max_occ_ind = %ld\n", i >> CP_SHIFT_SCALAR);    
    fflush(stdout);

    _mm_free(sa_ms_byte);
    _mm_free(sa_ls_word);
    return 0;
}

int FMI_search::build_index(int for_only) {

    char *prefix = file_name;
    uint64_t startTick, endTick;
    startTick = __rdtsc();
    index_alloc = 0;

    std::string reference_seq;
    char pac_file_name[PATH_MAX];
    assert(strlen(prefix) < PATH_MAX - 32);
    strcpy_s(pac_file_name, PATH_MAX, prefix);
    strcat_s(pac_file_name, PATH_MAX, ".pac");
    //sprintf(pac_file_name, "%s.pac", prefix);
    if(for_only)
    {
        strcat_s(prefix, PATH_MAX, ".for_only");
        //sprintf(prefix, "%s.for_only", prefix);
    }
    pac2nt(pac_file_name, reference_seq, for_only);
    int64_t pac_len = reference_seq.length();
    int status;
    int64_t size = pac_len * sizeof(char);
    char *binary_ref_seq = (char *)_mm_malloc(size, 64);
    index_alloc += size;
    assert_not_null(binary_ref_seq, size, index_alloc);
    char binary_ref_name[PATH_MAX];
    strcpy_s(binary_ref_name, PATH_MAX, prefix);
    strcat_s(binary_ref_name, PATH_MAX, ".0123");
    //sprintf(binary_ref_name, "%s.0123", prefix);
    std::fstream binary_ref_stream (binary_ref_name, std::ios::out | std::ios::binary);
    binary_ref_stream.seekg(0);
    endTick = __rdtsc();
    fprintf(stderr, "init ticks = %ld\n", endTick - startTick);
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
    endTick = __rdtsc();
    fprintf(stderr, "binary seq ticks = %ld\n", endTick - startTick);
    startTick = __rdtsc();

    size = (pac_len + 2) * sizeof(int64_t);
    int64_t *suffix_array=(int64_t *)_mm_malloc(size, 64);
    index_alloc += size;
    assert_not_null(suffix_array, size, index_alloc);
    startTick = __rdtsc();
    //status = saisxx<const char *, int64_t *, int64_t>(reference_seq.c_str(), suffix_array + 1, pac_len, 4);
    status = saisxx(reference_seq.c_str(), suffix_array + 1, pac_len);
    suffix_array[0] = pac_len;
    endTick = __rdtsc();
    fprintf(stderr, "build index ticks = %ld\n", endTick - startTick);
    startTick = __rdtsc();

    build_fm_index_avx(prefix, binary_ref_seq, pac_len, suffix_array, count);
    build_fm_index_scalar(prefix, binary_ref_seq, pac_len, suffix_array, count);
    _mm_free(binary_ref_seq);
    _mm_free(suffix_array);
    return 0;
}

#else

int FMI_search::build_fm_index(const char *ref_file_name,
                               char *binary_seq,
                               int64_t ref_seq_len,
                               int64_t *sa_bwt,
                               int64_t *count) {
    
    printf("ref_seq_len = %ld\n", ref_seq_len);
    fflush(stdout);

    char outname[PATH_MAX];

    strcpy_s(outname, PATH_MAX, ref_file_name);
    strcat_s(outname, PATH_MAX, CP_FILENAME_SUFFIX);
    //sprintf(outname, "%s.bwt.2bit.%d", ref_file_name, CP_BLOCK_SIZE);

    std::fstream outstream (outname, std::ios::out | std::ios::binary);
    outstream.seekg(0); 

    printf("count = %ld, %ld, %ld, %ld, %ld\n", count[0], count[1], count[2], count[3], count[4]);
    fflush(stdout);

    uint8_t *bwt;

    ref_seq_len++;
    outstream.write((char *)(&ref_seq_len), 1 * sizeof(int64_t));
    outstream.write((char*)count, 5 * sizeof(int64_t));

    int64_t i;
    int64_t ref_seq_len_aligned = ((ref_seq_len + CP_BLOCK_SIZE - 1) / CP_BLOCK_SIZE) * CP_BLOCK_SIZE;
    int64_t size = ref_seq_len_aligned * sizeof(uint8_t);
    bwt = (uint8_t *)_mm_malloc(size, 64);
    assert_not_null(bwt, size, index_alloc);

    int64_t sentinel_index = -1;
    for(i=0; i< ref_seq_len; i++)
    {
        if(sa_bwt[i] == 0)
        {
            bwt[i] = 4;
            printf("BWT[%ld] = 4\n", i);
            sentinel_index = i;
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
                        fprintf(stderr, "ERROR! i = %ld, c = %c\n", i, c);
                        exit(EXIT_FAILURE);
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

    size = cp_occ_size * sizeof(CP_OCC);
    cp_occ = (CP_OCC *)_mm_malloc(size, 64);
    assert_not_null(cp_occ, size, index_alloc);
    memset(cp_occ, 0, cp_occ_size * sizeof(CP_OCC));
    int64_t cp_count[16];

    memset(cp_count, 0, 16 * sizeof(int64_t));
    for(i = 0; i < ref_seq_len; i++)
    {
        if((i & CP_MASK) == 0)
        {
            CP_OCC cpo;
            cpo.cp_count[0] = cp_count[0];
            cpo.cp_count[1] = cp_count[1];
            cpo.cp_count[2] = cp_count[2];
            cpo.cp_count[3] = cp_count[3];

            int32_t j;
            cpo.one_hot_bwt_str[0] = 0;
            cpo.one_hot_bwt_str[1] = 0;
            cpo.one_hot_bwt_str[2] = 0;
            cpo.one_hot_bwt_str[3] = 0;

            for(j = 0; j < CP_BLOCK_SIZE; j++)
            {
                cpo.one_hot_bwt_str[0] = cpo.one_hot_bwt_str[0] << 1;
                cpo.one_hot_bwt_str[1] = cpo.one_hot_bwt_str[1] << 1;
                cpo.one_hot_bwt_str[2] = cpo.one_hot_bwt_str[2] << 1;
                cpo.one_hot_bwt_str[3] = cpo.one_hot_bwt_str[3] << 1;
                uint8_t c = bwt[i + j];
                //printf("c = %d\n", c);
                if(c < 4)
                {
                    cpo.one_hot_bwt_str[c] += 1;
                }
            }

            cp_occ[i >> CP_SHIFT] = cpo;
        }
        cp_count[bwt[i]]++;
    }
    outstream.write((char*)cp_occ, cp_occ_size * sizeof(CP_OCC));
    _mm_free(cp_occ);
    _mm_free(bwt);

    #if SA_COMPRESSION  

    size = ((ref_seq_len >> SA_COMPX)+ 1)  * sizeof(uint32_t);
    uint32_t *sa_ls_word = (uint32_t *)_mm_malloc(size, 64);
    assert_not_null(sa_ls_word, size, index_alloc);
    size = ((ref_seq_len >> SA_COMPX) + 1) * sizeof(int8_t);
    int8_t *sa_ms_byte = (int8_t *)_mm_malloc(size, 64);
    assert_not_null(sa_ms_byte, size, index_alloc);
    int64_t pos = 0;
    for(i = 0; i < ref_seq_len; i++)
    {
        if ((i & SA_COMPX_MASK) == 0)
        {
            sa_ls_word[pos] = sa_bwt[i] & 0xffffffff;
            sa_ms_byte[pos] = (sa_bwt[i] >> 32) & 0xff;
            pos++;
        }
    }
    fprintf(stderr, "pos: %d, ref_seq_len__: %ld\n", pos, ref_seq_len >> SA_COMPX);
    outstream.write((char*)sa_ms_byte, ((ref_seq_len >> SA_COMPX) + 1) * sizeof(int8_t));
    outstream.write((char*)sa_ls_word, ((ref_seq_len >> SA_COMPX) + 1) * sizeof(uint32_t));
    
    #else
    
    size = ref_seq_len * sizeof(uint32_t);
    uint32_t *sa_ls_word = (uint32_t *)_mm_malloc(size, 64);
    assert_not_null(sa_ls_word, size, index_alloc);
    size = ref_seq_len * sizeof(int8_t);
    int8_t *sa_ms_byte = (int8_t *)_mm_malloc(size, 64);
    assert_not_null(sa_ms_byte, size, index_alloc);
    for(i = 0; i < ref_seq_len; i++)
    {
        sa_ls_word[i] = sa_bwt[i] & 0xffffffff;
        sa_ms_byte[i] = (sa_bwt[i] >> 32) & 0xff;
    }
    outstream.write((char*)sa_ms_byte, ref_seq_len * sizeof(int8_t));
    outstream.write((char*)sa_ls_word, ref_seq_len * sizeof(uint32_t));
    
    #endif

    outstream.write((char *)(&sentinel_index), 1 * sizeof(int64_t));
    outstream.close();
    printf("max_occ_ind = %ld\n", i >> CP_SHIFT);    
    fflush(stdout);

    _mm_free(sa_ms_byte);
    _mm_free(sa_ls_word);
    return 0;
}

int FMI_search::build_index(int for_only) {

    char *prefix = file_name;
    uint64_t startTick;
    startTick = __rdtsc();
    index_alloc = 0;

    std::string reference_seq;
    char pac_file_name[PATH_MAX];
    strcpy_s(pac_file_name, PATH_MAX, prefix);
    strcat_s(pac_file_name, PATH_MAX, ".pac");
    //sprintf(pac_file_name, "%s.pac", prefix);
    if(for_only)
    {
        strcat_s(prefix, PATH_MAX, ".for_only");
        //sprintf(prefix, "%s.for_only", prefix);
    }
    // pac2nt(pac_file_name, reference_seq);
    pac2nt(pac_file_name, reference_seq, for_only);
    int64_t pac_len = reference_seq.length();
    int status;
    int64_t size = pac_len * sizeof(char);
    char *binary_ref_seq = (char *)_mm_malloc(size, 64);
    index_alloc += size;
    assert_not_null(binary_ref_seq, size, index_alloc);
    char binary_ref_name[PATH_MAX];
    strcpy_s(binary_ref_name, PATH_MAX, prefix);
    strcat_s(binary_ref_name, PATH_MAX, ".0123");
    //sprintf(binary_ref_name, "%s.0123", prefix);
    std::fstream binary_ref_stream (binary_ref_name, std::ios::out | std::ios::binary);
    binary_ref_stream.seekg(0);
    fprintf(stderr, "init ticks = %llu\n", __rdtsc() - startTick);
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
    fprintf(stderr, "binary seq ticks = %llu\n", __rdtsc() - startTick);
    startTick = __rdtsc();

    size = (pac_len + 2) * sizeof(int64_t);
    int64_t *suffix_array=(int64_t *)_mm_malloc(size, 64);
    index_alloc += size;
    assert_not_null(suffix_array, size, index_alloc);
    startTick = __rdtsc();
    //status = saisxx<const char *, int64_t *, int64_t>(reference_seq.c_str(), suffix_array + 1, pac_len, 4);
    status = saisxx(reference_seq.c_str(), suffix_array + 1, pac_len);
    suffix_array[0] = pac_len;
    fprintf(stderr, "build suffix-array ticks = %llu\n", __rdtsc() - startTick);
    startTick = __rdtsc();

    build_fm_index(prefix, binary_ref_seq, pac_len, suffix_array, count);
    fprintf(stderr, "build fm-index ticks = %llu\n", __rdtsc() - startTick);
    _mm_free(binary_ref_seq);
    _mm_free(suffix_array);
    return 0;
}

#endif


int FMI_search::build_index_forward_only()
{
    return build_index(1);
}
    
int FMI_search::build_index_with_rev_complement()
{
    return build_index(0);
}


#if OLD
void FMI_search::load_index(int for_only)
{
    char *prefix = file_name;
    assert(strlen(prefix) < PATH_MAX - 32);
    if(for_only)
    {
        strcat_s(prefix, PATH_MAX, ".for_only");
        //sprintf(prefix, "%s.for_only", prefix);
    }
    //beCalls = 0;
    char cp_file_name[PATH_MAX];
    assert(strnlen(prefix, 1000) + 12 < 1000);
    strcpy_s(cp_file_name, PATH_MAX, prefix);
    strcat_s(cp_file_name, PATH_MAX, CP_FILENAME_SUFFIX);
#if 0
#if ((!__AVX2__))
    sprintf(cp_file_name, "%s.bwt.2bit.%d", prefix, CP_BLOCK_SIZE_SCALAR);
#else
    sprintf(cp_file_name, "%s.bwt.8bit.%d", prefix, CP_BLOCK_SIZE_AVX);
#endif
#endif
    // Read the BWT and FM index of the reference sequence
    FILE *cpstream = NULL;
    cpstream = fopen(cp_file_name,"rb");
    if (cpstream == NULL)
    {
        fprintf(stderr, "ERROR! Unable to open the file: %s\n", cp_file_name);
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stderr, "Index file found. Loading index from %s\n", cp_file_name);
    }

    err_fread_noeof(&reference_seq_len, sizeof(int64_t), 1, cpstream);
    assert(reference_seq_len > 0);
    assert(reference_seq_len <= 0x7fffffffffL);

    fprintf(stderr, "* Reference seq len for bi-index = %ld\n", reference_seq_len);

    // create checkpointed occ
    int64_t cp_occ_size = (reference_seq_len >> CP_SHIFT) + 1;
    cp_occ = NULL;

    err_fread_noeof(&count[0], sizeof(int64_t), 5, cpstream);
    if ((cp_occ = (CP_OCC *)_mm_malloc(cp_occ_size * sizeof(CP_OCC), 64)) == NULL) {
        fprintf(stderr, "ERROR! unable to allocated cp_occ memory\n");
        exit(EXIT_FAILURE);
    }

    err_fread_noeof(cp_occ, sizeof(CP_OCC), cp_occ_size, cpstream);
    int64_t ii = 0;
    for(ii = 0; ii < 5; ii++)// update read count structure
    {
        count[ii] = count[ii] + 1;
    }
    sa_ms_byte = (int8_t *)_mm_malloc(reference_seq_len * sizeof(int8_t), 64);
    sa_ls_word = (uint32_t *)_mm_malloc(reference_seq_len * sizeof(uint32_t), 64);
    err_fread_noeof(sa_ms_byte, sizeof(int8_t), reference_seq_len, cpstream);
    err_fread_noeof(sa_ls_word, sizeof(uint32_t), reference_seq_len, cpstream);
    fclose(cpstream);

    sentinel_index = -1;
    int64_t x;
    for(x = 0; x < reference_seq_len; x++)
    {
        if(get_sa_entry(x) == 0)
            sentinel_index = x;
    }

    fprintf(stderr, "* Count:\n");
    for(x = 0; x < 5; x++)
    {
        fprintf(stderr, "%ld,\t%lu\n", x, (unsigned long)count[x]);
    }
    fprintf(stderr, "\n");  

#if ((!__AVX2__))
    base_mask[0][0] = 0;
    base_mask[0][1] = 0;
    base_mask[1][0] = 0xffffffffffffffffL;
    base_mask[1][1] = 0;
    base_mask[2][0] = 0;
    base_mask[2][1] = 0xffffffffffffffffL;
    base_mask[3][0] = 0xffffffffffffffffL;
    base_mask[3][1] = 0xffffffffffffffffL;
#else
    c_bcast_array = (uint8_t *)_mm_malloc(256 * sizeof(uint8_t), 64);
    for(ii = 0; ii < 4; ii++)
    {
        int32_t j;
        for(j = 0; j < 64; j++)
        {
            c_bcast_array[ii * 64 + j] = ii;
        }
    }
#endif
    fprintf(stderr, "* Done reading Index!!\n");
}
#else

void FMI_search::load_index(int for_only)
{
    one_hot_mask_array = (uint64_t *)_mm_malloc(64 * sizeof(uint64_t), 64);
    one_hot_mask_array[0] = 0;
    uint64_t base = 0x8000000000000000L;
    one_hot_mask_array[1] = base;
    int64_t i = 0;
    for(i = 2; i < 64; i++)
    {
        one_hot_mask_array[i] = (one_hot_mask_array[i - 1] >> 1) | base;
    }

    char *ref_file_name = file_name;
    if(for_only)
    {
        strcat_s(ref_file_name, PATH_MAX, ".for_only");
        //sprintf(prefix, "%s.for_only", prefix);
    }

    //beCalls = 0;
    char cp_file_name[PATH_MAX];
    strcpy_s(cp_file_name, PATH_MAX, ref_file_name);
    strcat_s(cp_file_name, PATH_MAX, CP_FILENAME_SUFFIX);

    // Read the BWT and FM index of the reference sequence
    FILE *cpstream = NULL;
    cpstream = fopen(cp_file_name,"rb");
    if (cpstream == NULL)
    {
        fprintf(stderr, "ERROR! Unable to open the file: %s\n", cp_file_name);
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stderr, "* Index file found. Loading index from %s\n", cp_file_name);
    }

    err_fread_noeof(&reference_seq_len, sizeof(int64_t), 1, cpstream);
    assert(reference_seq_len > 0);
    assert(reference_seq_len <= 0x7fffffffffL);

    fprintf(stderr, "* Reference seq len for bi-index = %ld\n", reference_seq_len);

    // create checkpointed occ
    int64_t cp_occ_size = (reference_seq_len >> CP_SHIFT) + 1;
    cp_occ = NULL;

    err_fread_noeof(&count[0], sizeof(int64_t), 5, cpstream);
    if ((cp_occ = (CP_OCC *)_mm_malloc(cp_occ_size * sizeof(CP_OCC), 64)) == NULL) {
        fprintf(stderr, "ERROR! unable to allocated cp_occ memory\n");
        exit(EXIT_FAILURE);
    }

    err_fread_noeof(cp_occ, sizeof(CP_OCC), cp_occ_size, cpstream);
    int64_t ii = 0;
    for(ii = 0; ii < 5; ii++)// update read count structure
    {
        count[ii] = count[ii] + 1;
    }

    #if SA_COMPRESSION

    int64_t reference_seq_len_ = (reference_seq_len >> SA_COMPX) + 1;
    sa_ms_byte = (int8_t *)_mm_malloc(reference_seq_len_ * sizeof(int8_t), 64);
    sa_ls_word = (uint32_t *)_mm_malloc(reference_seq_len_ * sizeof(uint32_t), 64);
    err_fread_noeof(sa_ms_byte, sizeof(int8_t), reference_seq_len_, cpstream);
    err_fread_noeof(sa_ls_word, sizeof(uint32_t), reference_seq_len_, cpstream);
    
    #else
    
    sa_ms_byte = (int8_t *)_mm_malloc(reference_seq_len * sizeof(int8_t), 64);
    sa_ls_word = (uint32_t *)_mm_malloc(reference_seq_len * sizeof(uint32_t), 64);
    err_fread_noeof(sa_ms_byte, sizeof(int8_t), reference_seq_len, cpstream);
    err_fread_noeof(sa_ls_word, sizeof(uint32_t), reference_seq_len, cpstream);

    #endif

    sentinel_index = -1;
    #if SA_COMPRESSION
    err_fread_noeof(&sentinel_index, sizeof(int64_t), 1, cpstream);
    fprintf(stderr, "* sentinel-index: %ld\n", sentinel_index);
    #endif
    fclose(cpstream);

    int64_t x;
    #if !SA_COMPRESSION
    for(x = 0; x < reference_seq_len; x++)
    {
        // fprintf(stderr, "x: %ld\n", x);
        #if SA_COMPRESSION
        if(get_sa_entry_compressed(x) == 0) {
            sentinel_index = x;
            break;
        }
        #else
        if(get_sa_entry(x) == 0) {
            sentinel_index = x;
            break;
        }
        #endif
    }
    fprintf(stderr, "\nsentinel_index: %ld\n", x);    
    #endif

    fprintf(stderr, "* Count:\n");
    for(x = 0; x < 5; x++)
    {
        fprintf(stderr, "%ld,\t%lu\n", x, (unsigned long)count[x]);
    }
    fprintf(stderr, "\n");  

    fprintf(stderr, "* Reading other elements of the index from files %s\n",
            ref_file_name);
    bwa_idx_load_ele(ref_file_name, BWA_IDX_ALL);

    fprintf(stderr, "* Done reading Index!!\n");
}

#endif


void FMI_search::load_index_forward_only()
{
    load_index(1);
}

void FMI_search::load_index_with_rev_complement()
{
    load_index(0);
}

void FMI_search::getSMEMsSeedForwardShrink(uint8_t *enc_qdb,
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
					 int *lisa_min_intv)
{
    int64_t numTotalSmem = *__numTotalSmem;

    uint32_t i;
    // Perform SMEM for original reads
    int numPrev = 0;
    for(i = 0; i < numReads; i++)
    {
        int x = query_pos_array[i];
        int32_t rid = rid_array[i];

        int readlength = seq_[rid].l_seq;
        int offset = query_cum_len_ar[rid];
        // uint8_t a = enc_qdb[rid * readlength + x];
        uint8_t a = enc_qdb[offset + x];

        if(a < 4)
        {
            SMEM smem;
            smem.rid = rid;
            smem.m = x;
            smem.n = x;
            smem.k = count[a];
            smem.l = count[3 - a];
            smem.s = count[a+1] - count[a];
            
            int j;
            for(j = x + 1; j < readlength; j++)
            {
                // a = enc_qdb[rid * readlength + j];
                a = enc_qdb[offset + j];
                if(a < 4)
                {
                    SMEM smem_ = smem;

                    // Forward extension is backward extension with the BWT of reverse complement
                    smem_.k = smem.l;
                    smem_.l = smem.k;
                    SMEM newSmem_ = backwardExt(smem_, 3 - a);
                    //SMEM newSmem_ = forwardExt(smem_, 3 - a);
                    SMEM newSmem = newSmem_;
                    newSmem.k = newSmem_.l;
                    newSmem.l = newSmem_.k;
                    newSmem.n = j;

             //       int32_t s_neq_mask = newSmem.s != smem.s;

            //      prevArray[numPrev] = smem;
            //      numPrev += s_neq_mask;
                    if(newSmem.s < min_intv_array[i])
                    {
                        break;
                    }
                    smem = newSmem;
#ifdef ENABLE_PREFETCH
                    _mm_prefetch((const char *)(&cp_occ[(smem.k) >> CP_SHIFT]), _MM_HINT_T0);
                    _mm_prefetch((const char *)(&cp_occ[(smem.l) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                }
                else
                {
                    break;
                }
            }
            if(smem.s >= min_intv_array[i])
            {
//                prevArray[numPrev] = smem;
		  matchArray[numPrev] = smem;		
		  lisa_min_intv[numPrev++] = min_intv_array[i] - 1;		
            }

        }
    }
*__numTotalSmem = numPrev;
}



void FMI_search::getSMEMsOnePosOneThread(uint8_t *enc_qdb,
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
                                         int64_t *__numTotalSmem)
{
    int64_t numTotalSmem = *__numTotalSmem;
    SMEM prevArray[max_readlength];

    uint32_t i;
    // Perform SMEM for original reads
    for(i = 0; i < numReads; i++)
    {
        int x = query_pos_array[i];
        int32_t rid = rid_array[i];
        int next_x = x + 1;

        int readlength = seq_[rid].l_seq;
        int offset = query_cum_len_ar[rid];
        // uint8_t a = enc_qdb[rid * readlength + x];
        uint8_t a = enc_qdb[offset + x];

        if(a < 4)
        {
            SMEM smem;
            smem.rid = rid;
            smem.m = x;
            smem.n = x;
            smem.k = count[a];
            smem.l = count[3 - a];
            smem.s = count[a+1] - count[a];
            int numPrev = 0;
            
            int j;
            for(j = x + 1; j < readlength; j++)
            {
                // a = enc_qdb[rid * readlength + j];
                a = enc_qdb[offset + j];
                next_x = j + 1;
                if(a < 4)
                {
                    SMEM smem_ = smem;

                    // Forward extension is backward extension with the BWT of reverse complement
                    smem_.k = smem.l;
                    smem_.l = smem.k;
                    SMEM newSmem_ = backwardExt(smem_, 3 - a);
                    //SMEM newSmem_ = forwardExt(smem_, 3 - a);
                    SMEM newSmem = newSmem_;
                    newSmem.k = newSmem_.l;
                    newSmem.l = newSmem_.k;
                    newSmem.n = j;

                    int32_t s_neq_mask = newSmem.s != smem.s;

                    prevArray[numPrev] = smem;
                    numPrev += s_neq_mask;
                    if(newSmem.s < min_intv_array[i])
                    {
                        next_x = j;
                        break;
                    }
                    smem = newSmem;
#ifdef ENABLE_PREFETCH
                    _mm_prefetch((const char *)(&cp_occ[(smem.k) >> CP_SHIFT]), _MM_HINT_T0);
                    _mm_prefetch((const char *)(&cp_occ[(smem.l) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                }
                else
                {
                    break;
                }
            }
            if(smem.s >= min_intv_array[i])
            {

                prevArray[numPrev] = smem;
                numPrev++;
		//fprintf(stderr, "FMI: backward staart %ld %ld %ld %ld\n", smem.m, smem.n, min_intv_array[i], smem.s);
            }

            SMEM *prev;
            prev = prevArray;

            int p;
            for(p = 0; p < (numPrev/2); p++)
            {
                SMEM temp = prev[p];
                prev[p] = prev[numPrev - p - 1];
                prev[numPrev - p - 1] = temp;
            }

            // Backward search
            int cur_j = readlength;
            for(j = x - 1; j >= 0; j--)
            {
                int numCurr = 0;
                int curr_s = -1;
                // a = enc_qdb[rid * readlength + j];
                a = enc_qdb[offset + j];

                if(a > 3)
                {
                    break;
                }
                for(p = 0; p < numPrev; p++)
                {
                    SMEM smem = prev[p];
                    SMEM newSmem = backwardExt(smem, a);
                    newSmem.m = j;

                    if((newSmem.s < min_intv_array[i]) && ((smem.n - smem.m + 1) >= minSeedLen))
                    {
                        cur_j = j;

                        matchArray[numTotalSmem++] = smem;
                        break;
                    }
                    if((newSmem.s >= min_intv_array[i]) && (newSmem.s != curr_s))
                    {
                        curr_s = newSmem.s;
                        prev[numCurr++] = newSmem;
#ifdef ENABLE_PREFETCH
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k) >> CP_SHIFT]), _MM_HINT_T0);
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k + newSmem.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                        break;
                    }
                }
                p++;
                for(; p < numPrev; p++)
                {
                    SMEM smem = prev[p];

                    SMEM newSmem = backwardExt(smem, a);
                    newSmem.m = j;


                    if((newSmem.s >= min_intv_array[i]) && (newSmem.s != curr_s))
                    {
                        curr_s = newSmem.s;
                        prev[numCurr++] = newSmem;
#ifdef ENABLE_PREFETCH
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k) >> CP_SHIFT]), _MM_HINT_T0);
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k + newSmem.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                    }
                }
                numPrev = numCurr;
                if(numCurr == 0)
                {
                    break;
                }
            }
            if(numPrev != 0)
            {
                SMEM smem = prev[0];
                if(((smem.n - smem.m + 1) >= minSeedLen))
                {

                    matchArray[numTotalSmem++] = smem;
                }
                numPrev = 0;
            }
        }
        query_pos_array[i] = next_x;
    }
    (*__numTotalSmem) = numTotalSmem;
}

void FMI_search::getSMEMsAllPosOneThread(uint8_t *enc_qdb,
                                         int32_t *min_intv_array,
                                         int32_t *rid_array,
                                         int32_t numReads,
                                         int32_t batch_size,
                                         const bseq1_t *seq_,
                                         int32_t *query_cum_len_ar,
                                         int32_t max_readlength,
                                         int32_t minSeedLen,
                                         SMEM *matchArray,
                                         int64_t *__numTotalSmem)
{
    int16_t *query_pos_array = (int16_t *)_mm_malloc(numReads * sizeof(int16_t), 64);
    
    int32_t i;
    for(i = 0; i < numReads; i++)
        query_pos_array[i] = 0;

    int32_t numActive = numReads;
    (*__numTotalSmem) = 0;

    do
    {
        int32_t head = 0;
        int32_t tail = 0;
        for(head = 0; head < numActive; head++)
        {
            int readlength = seq_[rid_array[head]].l_seq;
            if(query_pos_array[head] < readlength)
            {
                rid_array[tail] = rid_array[head];
                query_pos_array[tail] = query_pos_array[head];
                min_intv_array[tail] = min_intv_array[head];
                tail++;             
            }               
        }
        getSMEMsOnePosOneThread(enc_qdb,
                                query_pos_array,
                                min_intv_array,
                                rid_array,
                                tail,
                                batch_size,
                                seq_,
                                query_cum_len_ar,
                                max_readlength,
                                minSeedLen,
                                matchArray,
                                __numTotalSmem);
        numActive = tail;
    } while(numActive > 0);

    _mm_free(query_pos_array);
}



int64_t FMI_search::bwtSeedStrategyAllPosOneThread(uint8_t *enc_qdb,
                                                   int32_t *max_intv_array,
                                                   int32_t numReads,
                                                   const bseq1_t *seq_,
                                                   int32_t *query_cum_len_ar,
                                                   int32_t minSeedLen,
                                                   SMEM *matchArray)
{
    int32_t i;

    int64_t numTotalSeed = 0;

    for(i = 0; i < numReads; i++)
    {
        int readlength = seq_[i].l_seq;
        int16_t x = 0;
        while(x < readlength)
        {
            int next_x = x + 1;

            // Forward search
            SMEM smem;
            smem.rid = i;
            smem.m = x;
            smem.n = x;
            
            int offset = query_cum_len_ar[i];
            uint8_t a = enc_qdb[offset + x];
            // uint8_t a = enc_qdb[i * readlength + x];

            if(a < 4)
            {
                smem.k = count[a];
                smem.l = count[3 - a];
                smem.s = count[a+1] - count[a];


                int j;
                for(j = x + 1; j < readlength; j++)
                {
                    next_x = j + 1;
                    // a = enc_qdb[i * readlength + j];
                    a = enc_qdb[offset + j];
                    if(a < 4)
                    {
                        SMEM smem_ = smem;

                        // Forward extension is backward extension with the BWT of reverse complement
                        smem_.k = smem.l;
                        smem_.l = smem.k;
                        SMEM newSmem_ = backwardExt(smem_, 3 - a);
                        //SMEM smem = backwardExt(smem, 3 - a);
                        //smem.n = j;
                        SMEM newSmem = newSmem_;
                        newSmem.k = newSmem_.l;
                        newSmem.l = newSmem_.k;
                        newSmem.n = j;
                        smem = newSmem;
#ifdef ENABLE_PREFETCH
                        _mm_prefetch((const char *)(&cp_occ[(smem.k) >> CP_SHIFT]), _MM_HINT_T0);
                        _mm_prefetch((const char *)(&cp_occ[(smem.l) >> CP_SHIFT]), _MM_HINT_T0);
#endif


		//	fprintf(stderr, "fmi-tal-log: %ld %ld %ld %ld %ld %ld %ld\n", 3 - a, smem.m, smem.n, smem.k, smem.k + smem.s, smem.s, (long) max_intv_array[i]);
                        if((smem.s < max_intv_array[i]) && ((smem.n - smem.m + 1) >= minSeedLen)) // acgtcg
                        {

                            if(smem.s > 0)
                            {
                                matchArray[numTotalSeed++] = smem;
                            }
                            break;
                        }
                    }
                    else
                    {

                        break;
                    }
                }

            }
            x = next_x;
        }
    }
    return numTotalSeed;
}


void FMI_search::getSMEMs(uint8_t *enc_qdb,
        int32_t numReads,
        int32_t batch_size,
        int32_t readlength,
        int32_t minSeedLen,
        int32_t nthreads,
        SMEM *matchArray,
        int64_t *numTotalSmem)
{
    SMEM *prevArray = (SMEM *)_mm_malloc(nthreads * readlength * sizeof(SMEM), 64);
    SMEM *currArray = (SMEM *)_mm_malloc(nthreads * readlength * sizeof(SMEM), 64);


// #pragma omp parallel num_threads(nthreads)
    {
        int tid = 0; //omp_get_thread_num();   // removed omp
        numTotalSmem[tid] = 0;
        SMEM *myPrevArray = prevArray + tid * readlength;
        SMEM *myCurrArray = prevArray + tid * readlength;

        int32_t perThreadQuota = (numReads + (nthreads - 1)) / nthreads;
        int32_t first = tid * perThreadQuota;
        int32_t last  = (tid + 1) * perThreadQuota;
        if(last > numReads) last = numReads;
        SMEM *myMatchArray = matchArray + first * readlength;

        uint32_t i;
        // Perform SMEM for original reads
        for(i = first; i < last; i++)
        {
            int x = readlength - 1;
            int numPrev = 0;
            int numSmem = 0;

            while (x >= 0)
            {
                // Forward search
                SMEM smem;
                smem.rid = i;
                smem.m = x;
                smem.n = x;
                uint8_t a = enc_qdb[i * readlength + x];

                if(a > 3)
                {
                    x--;
                    continue;
                }
                smem.k = count[a];
                smem.l = count[3 - a];
                smem.s = count[a+1] - count[a];

                int j;
                for(j = x + 1; j < readlength; j++)
                {
                    a = enc_qdb[i * readlength + j];
                    if(a < 4)
                    {
                        SMEM smem_ = smem;

                        // Forward extension is backward extension with the BWT of reverse complement
                        smem_.k = smem.l;
                        smem_.l = smem.k;
                        SMEM newSmem_ = backwardExt(smem_, 3 - a);
                        SMEM newSmem = newSmem_;
                        newSmem.k = newSmem_.l;
                        newSmem.l = newSmem_.k;
                        newSmem.n = j;

                        if(newSmem.s != smem.s)
                        {
                            myPrevArray[numPrev] = smem;
                            numPrev++;
                        }
                        smem = newSmem;
                        if(newSmem.s == 0)
                        {
                            break;
                        }
                    }
                    else
                    {
                        myPrevArray[numPrev] = smem;
                        numPrev++;
                        break;
                    }
                }
                if(smem.s != 0)
                {
                    myPrevArray[numPrev++] = smem;
                }

                SMEM *curr, *prev;
                prev = myPrevArray;
                curr = myCurrArray;

                int p;
                for(p = 0; p < (numPrev/2); p++)
                {
                    SMEM temp = prev[p];
                    prev[p] = prev[numPrev - p - 1];
                    prev[numPrev - p - 1] = temp;
                }

                int next_x = x - 1;

                // Backward search
                int cur_j = readlength;
                for(j = x - 1; j >= 0; j--)
                {
                    int numCurr = 0;
                    int curr_s = -1;
                    a = enc_qdb[i * readlength + j];
                    //printf("a = %d\n", a);
                    if(a > 3)
                    {
                        next_x = j - 1;
                        break;
                    }
                    for(p = 0; p < numPrev; p++)
                    {
                        SMEM smem = prev[p];
                        SMEM newSmem = backwardExt(smem, a);
                        newSmem.m = j;

                        if(newSmem.s == 0)
                        {
                            if((numCurr == 0) && (j < cur_j))
                            {
                                cur_j = j;
                                if((smem.n - smem.m + 1) >= minSeedLen)
                                    myMatchArray[numTotalSmem[tid] + numSmem++] = smem;
                            }
                        }
                        if((newSmem.s != 0) && (newSmem.s != curr_s))
                        {
                            curr_s = newSmem.s;
                            curr[numCurr++] = newSmem;
                        }
                    }
                    SMEM *temp = prev;
                    prev = curr;
                    curr = temp;
                    numPrev = numCurr;
                    if(numCurr == 0)
                    {
                        next_x = j;
                        break;
                    }
                    else
                    {
                        next_x = j - 1;
                    }
                }
                if(numPrev != 0)
                {
                    SMEM smem = prev[0];
                    if((smem.n - smem.m + 1) >= minSeedLen)
                        myMatchArray[numTotalSmem[tid] + numSmem++] = smem;
                    numPrev = 0;
                }
                x = next_x;
            }
            numTotalSmem[tid] += numSmem;
        }
    }

    _mm_free(prevArray);
    _mm_free(currArray);
}


int compare_smem(const void *a, const void *b)
{
    SMEM *pa = (SMEM *)a;
    SMEM *pb = (SMEM *)b;

    if(pa->rid < pb->rid)
        return -1;
    if(pa->rid > pb->rid)
        return 1;

    if(pa->m < pb->m)
        return -1;
    if(pa->m > pb->m)
        return 1;
    if(pa->n > pb->n)
        return -1;
    if(pa->n < pb->n)
        return 1;
    return 0;
}

void FMI_search::sortSMEMs(SMEM *matchArray,
        int64_t numTotalSmem[],
        int32_t numReads,
        int32_t readlength,
        int nthreads)
{
    int tid;
    int32_t perThreadQuota = (numReads + (nthreads - 1)) / nthreads;
    for(tid = 0; tid < nthreads; tid++)
    {
        int32_t first = tid * perThreadQuota;
        SMEM *myMatchArray = matchArray + first * readlength;
        qsort(myMatchArray, numTotalSmem[tid], sizeof(SMEM), compare_smem);
    }
}

#if 0
SMEM FMI_search::backwardExt_back(SMEM smem, uint8_t a)
{
    //beCalls++;
    uint8_t b;

    int64_t k[4], l[4], s[4];
    for(b = 0; b < 4; b++)
    {
        int64_t sp = (int64_t)(smem.k);
        int64_t ep = (int64_t)(smem.k) + (int64_t)(smem.s);

        // #if ((!__AVX2__))
        // GET_OCC(sp, b, occ_id_sp, y_sp, occ_sp, bwt_str_bit0_sp, bwt_str_bit1_sp, bit0_cmp_sp, bit1_cmp_sp, mismatch_mask_sp);
        // GET_OCC(ep, b, occ_id_ep, y_ep, occ_ep, bwt_str_bit0_ep, bwt_str_bit1_ep, bit0_cmp_ep, bit1_cmp_ep, mismatch_mask_ep);
        GET_OCC(sp, b, occ_id_sp, y_sp, occ_sp, one_hot_bwt_str_c_sp, match_mask_sp);
        GET_OCC(ep, b, occ_id_ep, y_ep, occ_ep, one_hot_bwt_str_c_ep, match_mask_ep);
        // #else
        // __m256i b256;
        // b256 = _mm256_load_si256((const __m256i *)(c_bcast_array + b * 64));
        // GET_OCC(sp, b, b256, occ_id_sp, y_sp, occ_sp, bwt_str_sp, bwt_sp_vec, mask_sp_vec, mask_sp);
        // GET_OCC(ep, b, b256, occ_id_ep, y_ep, occ_ep, bwt_str_ep, bwt_ep_vec, mask_ep_vec, mask_ep);
        // #endif

        k[b] = count[b] + occ_sp;
        s[b] = occ_ep - occ_sp;
    }

    int64_t sentinel_offset = 0;
    if((smem.k <= sentinel_index) && ((smem.k + smem.s) > sentinel_index)) sentinel_offset = 1;
    l[3] = smem.l + sentinel_offset;
    l[2] = l[3] + s[3];
    l[1] = l[2] + s[2];
    l[0] = l[1] + s[1];

    smem.k = k[a];
    smem.l = l[a];
    smem.s = s[a];
    return smem;
}
#endif

SMEM FMI_search::backwardExt(SMEM smem, uint8_t a)
{
    //beCalls++;
    int8_t b;
    int64_t sentinel_offset = 0;
    if((smem.k <= sentinel_index) && ((smem.k + smem.s) > sentinel_index)) sentinel_offset = 1;

    int64_t k[4], l[4], s[4];
    l[3] = smem.l + sentinel_offset;
    for(b = 3; b >= a; b--)
    {
        int64_t sp = (int64_t)(smem.k);
        int64_t ep = (int64_t)(smem.k) + (int64_t)(smem.s);

        // #if ((!__AVX2__))
        // GET_OCC(sp, b, occ_id_sp, y_sp, occ_sp, bwt_str_bit0_sp, bwt_str_bit1_sp, bit0_cmp_sp, bit1_cmp_sp, mismatch_mask_sp);
        // GET_OCC(ep, b, occ_id_ep, y_ep, occ_ep, bwt_str_bit0_ep, bwt_str_bit1_ep, bit0_cmp_ep, bit1_cmp_ep, mismatch_mask_ep);
        GET_OCC(sp, b, occ_id_sp, y_sp, occ_sp, one_hot_bwt_str_c_sp, match_mask_sp);
        GET_OCC(ep, b, occ_id_ep, y_ep, occ_ep, one_hot_bwt_str_c_ep, match_mask_ep);
        // #else
        // __m256i b256;
        // b256 = _mm256_load_si256((const __m256i *)(c_bcast_array + b * 64));
        // GET_OCC(sp, b, b256, occ_id_sp, y_sp, occ_sp, bwt_str_sp, bwt_sp_vec, mask_sp_vec, mask_sp);
        // GET_OCC(ep, b, b256, occ_id_ep, y_ep, occ_ep, bwt_str_ep, bwt_ep_vec, mask_ep_vec, mask_ep);
        // #endif

        k[b] = count[b] + occ_sp;
        s[b] = occ_ep - occ_sp;
    }

    for(b = 2; b >= a; b--)
  	l[b] = l[b+1] + s[b+1];
  //  l[2] = l[3] + s[3];
  //  l[1] = l[2] + s[2];
  //  l[0] = l[1] + s[1];

    smem.k = k[a];
    smem.l = l[a];
    smem.s = s[a];
    return smem;
}

int64_t FMI_search::get_sa_entry(int64_t pos)
{
    int64_t sa_entry = sa_ms_byte[pos];
    sa_entry = sa_entry << 32;
    sa_entry = sa_entry + sa_ls_word[pos];
    return sa_entry;
}

void FMI_search::get_sa_entries(int64_t *posArray, int64_t *coordArray, uint32_t count, int32_t nthreads)
{
    uint32_t i;
// #pragma omp parallel for num_threads(nthreads)
    for(i = 0; i < count; i++)
    {
        int64_t pos = posArray[i];
        int64_t sa_entry = sa_ms_byte[pos];
        sa_entry = sa_entry << 32;
        sa_entry = sa_entry + sa_ls_word[pos];
        //_mm_prefetch((const char *)(sa_ms_byte + pos + SAL_PFD), _MM_HINT_T0);
        coordArray[i] = sa_entry;
    }
}

void FMI_search::get_sa_entries(SMEM *smemArray, int64_t *coordArray, int32_t *coordCountArray, uint32_t count, int32_t max_occ)
{
    uint32_t i;
    int32_t totalCoordCount = 0;
    for(i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        int64_t hi = smem.k + smem.s;
        int64_t step = (smem.s > max_occ) ? smem.s / max_occ : 1;
        int64_t j;
        for(j = smem.k; (j < hi) && (c < max_occ); j+=step, c++)
        {
            int64_t pos = j;
            int64_t sa_entry = sa_ms_byte[pos];
            sa_entry = sa_entry << 32;
            sa_entry = sa_entry + sa_ls_word[pos];
            //_mm_prefetch((const char *)(sa_ms_byte + pos + SAL_PFD * step), _MM_HINT_T0);
            coordArray[totalCoordCount + c] = sa_entry;
        }
        coordCountArray[i] = c;
        totalCoordCount += c;
    }
}

/****************************** sa_compression ***************************/
int64_t FMI_search::get_sa_entry_compressed(int64_t pos, int tid)
{
    if ((pos & SA_COMPX_MASK) == 0) {
        
        #if  SA_COMPRESSION
        int64_t sa_entry = sa_ms_byte[pos >> SA_COMPX];
        #else
        int64_t sa_entry = sa_ms_byte[pos];     // simulation
        #endif
        
        sa_entry = sa_entry << 32;
        
        #if  SA_COMPRESSION
        sa_entry = sa_entry + sa_ls_word[pos >> SA_COMPX];
        #else
        sa_entry = sa_entry + sa_ls_word[pos];   // simulation
        #endif
        
        return sa_entry;        
    }
    else {
        // tprof[MEM_CHAIN][tid] ++;
        int64_t offset = 0; 
        int64_t sp = pos;
        while(true)
        {
            int64_t occ_id_pp_ = sp >> CP_SHIFT;
            int64_t y_pp_ = CP_BLOCK_SIZE - (sp & CP_MASK) - 1; 
            uint64_t *one_hot_bwt_str = cp_occ[occ_id_pp_].one_hot_bwt_str;
            uint8_t b;

            if((one_hot_bwt_str[0] >> y_pp_) & 1)
                b = 0;
            else if((one_hot_bwt_str[1] >> y_pp_) & 1)
                b = 1;
            else if((one_hot_bwt_str[2] >> y_pp_) & 1)
                b = 2;
            else if((one_hot_bwt_str[3] >> y_pp_) & 1)
                b = 3;
            else
                b = 4;

            if (b == 4) {
                return offset;
            }

            GET_OCC(sp, b, occ_id_sp, y_sp, occ_sp, one_hot_bwt_str_c_sp, match_mask_sp);

            sp = count[b] + occ_sp;
            
            offset ++;
            if ((sp & SA_COMPX_MASK) == 0) break;
        }
        // assert((reference_seq_len >> SA_COMPX) - 1 >= (sp >> SA_COMPX));
        #if  SA_COMPRESSION
        int64_t sa_entry = sa_ms_byte[sp >> SA_COMPX];
        #else
        int64_t sa_entry = sa_ms_byte[sp];      // simultion
        #endif
        
        sa_entry = sa_entry << 32;

        #if  SA_COMPRESSION
        sa_entry = sa_entry + sa_ls_word[sp >> SA_COMPX];
        #else
        sa_entry = sa_entry + sa_ls_word[sp];      // simulation
        #endif
        
        sa_entry += offset;
        return sa_entry;
    }
}

void FMI_search::get_sa_entries(SMEM *smemArray, int64_t *coordArray, int32_t *coordCountArray, uint32_t count, int32_t max_occ, int tid)
{
    
    uint32_t i;
    int32_t totalCoordCount = 0;
    for(i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        int64_t hi = smem.k + smem.s;
        int64_t step = (smem.s > max_occ) ? smem.s / max_occ : 1;
        int64_t j;
        for(j = smem.k; (j < hi) && (c < max_occ); j+=step, c++)
        {
            int64_t pos = j;
            int64_t sa_entry = get_sa_entry_compressed(pos, tid);
            coordArray[totalCoordCount + c] = sa_entry;
        }
        // coordCountArray[i] = c;
        *coordCountArray += c;
        totalCoordCount += c;
    }
}

// SA_COPMRESSION w/ PREFETCH
int64_t FMI_search::call_one_step(int64_t pos, int64_t &sa_entry, int64_t &offset)
{
    if ((pos & SA_COMPX_MASK) == 0) {        
        sa_entry = sa_ms_byte[pos >> SA_COMPX];        
        sa_entry = sa_entry << 32;        
        sa_entry = sa_entry + sa_ls_word[pos >> SA_COMPX];        
        // return sa_entry;
        return 1;
    }
    else {
        // int64_t offset = 0; 
        int64_t sp = pos;

        int64_t occ_id_pp_ = sp >> CP_SHIFT;
        int64_t y_pp_ = CP_BLOCK_SIZE - (sp & CP_MASK) - 1; 
        uint64_t *one_hot_bwt_str = cp_occ[occ_id_pp_].one_hot_bwt_str;
        uint8_t b;

        if((one_hot_bwt_str[0] >> y_pp_) & 1)
            b = 0;
        else if((one_hot_bwt_str[1] >> y_pp_) & 1)
            b = 1;
        else if((one_hot_bwt_str[2] >> y_pp_) & 1)
            b = 2;
        else if((one_hot_bwt_str[3] >> y_pp_) & 1)
            b = 3;
        else
            b = 4;
        if (b == 4) {
            sa_entry = 0;
            return 1;
        }
        
        GET_OCC(sp, b, occ_id_sp, y_sp, occ_sp, one_hot_bwt_str_c_sp, match_mask_sp);
        
        sp = count[b] + occ_sp;
        
        offset ++;
        if ((sp & SA_COMPX_MASK) == 0) {
    
            sa_entry = sa_ms_byte[sp >> SA_COMPX];        
            sa_entry = sa_entry << 32;
            sa_entry = sa_entry + sa_ls_word[sp >> SA_COMPX];
            
            sa_entry += offset;
            // return sa_entry;
            return 1;
        }
        else {
            sa_entry = sp;
            return 0;
        }
    } // else
}

void FMI_search::get_sa_entries_prefetch(SMEM *smemArray, int64_t *coordArray,
                                         int64_t *coordCountArray, int64_t count,
                                         const int32_t max_occ, int tid, int64_t &id_)
{
    
    // uint32_t i;
    int32_t totalCoordCount = 0;
    int32_t mem_lim = 0, id = 0;
    
    for(int i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        mem_lim += smem.s;
    }

    int64_t *pos_ar = (int64_t *) _mm_malloc( mem_lim * sizeof(int64_t), 64);
    int64_t *map_ar = (int64_t *) _mm_malloc( mem_lim * sizeof(int64_t), 64);
//LISA

#if REV_CMP_SEARCH
//    bool *lisa_rev_comp_ar = (bool *) _mm_malloc( mem_lim * sizeof(bool), 64);
    int lisa_rev_comp_ar[mem_lim];// * sizeof(bool), 64);
#endif
    for(int i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
#if REV_CMP_SEARCH
	bool rev_comp_flag = smem.l == -1000 ? true : false;
#endif     
   	int64_t hi = smem.k + smem.s;
        int64_t step = (smem.s > max_occ) ? smem.s / max_occ : 1;
        int64_t j;
        for(j = smem.k; (j < hi) && (c < max_occ); j+=step, c++)
        {
            int64_t pos = j;
             pos_ar[id]  = pos;
//LISA -----------------------------
#if REV_CMP_SEARCH
	    lisa_rev_comp_ar[id] = rev_comp_flag;
	     if (rev_comp_flag)
	     	coordArray[id] = (int64_t)774967607 - (smem.n - smem.m + 2);
// ------------------------------
#endif
             map_ar[id++] = totalCoordCount + c;
            // int64_t sa_entry = get_sa_entry_compressed(pos, tid);
            // coordArray[totalCoordCount + c] = sa_entry;
        }
        //coordCountArray[i] = c;
        *coordCountArray += c;
        totalCoordCount += c;
    }
    
    id_ += id;
    
    const int32_t sa_batch_size = 20;
    int64_t working_set[sa_batch_size], map_pos[sa_batch_size];;
// LISA
#if REV_CMP_SEARCH
    bool lisa_rev_comp_batch[sa_batch_size];
#endif
    int64_t offset[sa_batch_size] = {-1};
    
    int i = 0, j = 0;    
    while(i<id && j<sa_batch_size)
    {
        int64_t pos =  pos_ar[i];
        working_set[j] = pos;
        map_pos[j] = map_ar[i];
//LISA
#if REV_CMP_SEARCH
	lisa_rev_comp_batch[j] = lisa_rev_comp_ar[i];		
#endif
        offset[j] = 0;
        
        if (pos & SA_COMPX_MASK == 0) {
            _mm_prefetch(&sa_ms_byte[pos >> SA_COMPX], _MM_HINT_T0);
            _mm_prefetch(&sa_ls_word[pos >> SA_COMPX], _MM_HINT_T0);
        }
        else {
            int64_t occ_id_pp_ = pos >> CP_SHIFT;
            _mm_prefetch(&cp_occ[occ_id_pp_], _MM_HINT_T0);
        }
        i++;
        j++;
    }
        
    int lim = j, all_quit = 0;
    while (all_quit < id)
    {
        
        for (int k=0; k<lim; k++)
        {
            int64_t sp = 0, pos = 0;
            bool quit;
            if (offset[k] >= 0) {
                quit = call_one_step(working_set[k], sp, offset[k]);
            }
            else
                continue;
            
            if (quit) {
#if REV_CMP_SEARCH
		if(lisa_rev_comp_batch[k])
		{
			coordArray[map_pos[k]] -= sp;
	//		fprintf(stderr, " lisa:%lld\n", coordArray[map_pos[k]]);
		}
                else 
#endif
		{
			coordArray[map_pos[k]] = sp;
	//		fprintf(stderr, "%lld\n", sp);
                }
		all_quit ++;
                
                if (i < id)
                {
                    pos = pos_ar[i];
                    working_set[k] = pos;
//LISA
#if REV_CMP_SEARCH
		    lisa_rev_comp_batch[k] = lisa_rev_comp_ar[i];		
#endif
                    map_pos[k] = map_ar[i++];
                    offset[k] = 0;
                    
                    if (pos & SA_COMPX_MASK == 0) {
                        _mm_prefetch(&sa_ms_byte[pos >> SA_COMPX], _MM_HINT_T0);
                        _mm_prefetch(&sa_ls_word[pos >> SA_COMPX], _MM_HINT_T0);
                    }
                    else {
                        int64_t occ_id_pp_ = pos >> CP_SHIFT;
                        _mm_prefetch(&cp_occ[occ_id_pp_], _MM_HINT_T0);
                    }
                }
                else
                    offset[k] = -1;
            }
            else {
                working_set[k] = sp;
                if (sp & SA_COMPX_MASK == 0) {
                    _mm_prefetch(&sa_ms_byte[sp >> SA_COMPX], _MM_HINT_T0);
                    _mm_prefetch(&sa_ls_word[sp >> SA_COMPX], _MM_HINT_T0);
                }
                else {
                    int64_t occ_id_pp_ = sp >> CP_SHIFT;
                    _mm_prefetch(&cp_occ[occ_id_pp_], _MM_HINT_T0);
                }                
            }
        }
    }
  
	//sort(&coordArray[0], &coordArray[0] + id_);
//	for(int i = 0; i < id_; i++){
//		fprintf(stderr, "%lld\n", coordArray[i] );
//	}
 
    _mm_free(pos_ar);
    _mm_free(map_ar);
}



/***************************** Exact Match ******************************/
void FMI_search::exact_search(uint8_t *enc_qdb,
                 int64_t num_queries, uint32_t query_length,
                 int64_t *k_l_range,
                 int32_t batch_size)
{
    BatchMetadata metadata[batch_size];

    int64_t ii;
    for(ii = 0; ii < batch_size; ii++)
    {
        BatchMetadata m;
        int64_t rid = ii;
        int64_t ind = enc_qdb[rid * query_length + query_length - 1];
        m.k = count[ind];
        m.l = count[ind + 1];
        m.rid = rid;
        m.cur_ind = query_length - 2;
        metadata[ii] = m;
    }
    int64_t next_query = batch_size;
    int32_t num_active = batch_size;

    //loopCount = 0;

    while(num_active > 0)
    {
        int32_t new_num_active = num_active;
        for(ii = 0; ii < num_active; ii++)
        {
            int64_t sp, ep;
            int64_t rid;
            int j;
            sp = metadata[ii].k;
            ep = metadata[ii].l;
            rid = metadata[ii].rid;
            j = metadata[ii].cur_ind;

            //loopCount++;
            int64_t ind = enc_qdb[rid * query_length + j];
            // #if ((!__AVX2__))
            // #if (CP_BLOCK_SIZE <= 64)
            // GET_OCC(sp, ind, occ_id_sp, y_sp, occ_sp, bwt_str_bit0_sp, bwt_str_bit1_sp, bit0_cmp_sp, bit1_cmp_sp, mismatch_mask_sp);
            GET_OCC(sp, ind, occ_id_sp, y_sp, occ_sp, one_hot_bwt_str_c_sp, match_mask_sp);
            // #else
            // GET_OCC(sp, ind, occ_id_sp, y_sp, y0_sp, y1_sp, occ_sp, bwt_str_bit0_sp, bwt_str_bit1_sp, bit0_cmp_sp, bit1_cmp_sp, mismatch_mask_sp);
            // #endif
            // #else
            // __m256i c256;
            // //c256 = _mm256_set1_epi8(ind);
            // c256 = _mm256_load_si256((const __m256i *)(c_bcast_array + ind * 64));
            // GET_OCC(sp, ind, c256, occ_id_sp, y_sp, occ_sp, bwt_str_sp, bwt_sp_vec, mask_sp_vec, mask_sp);
            // #endif
            
            sp  = count[ind] + occ_sp;

            metadata[ii].k = sp;
            
            // #if ((!__AVX2__))
            // #if (CP_BLOCK_SIZE <= 64)
            // GET_OCC(ep, ind, occ_id_ep, y_ep, occ_ep, bwt_str_bit0_ep, bwt_str_bit1_ep, bit0_cmp_ep, bit1_cmp_ep, mismatch_mask_ep);
            GET_OCC(ep, ind, occ_id_ep, y_ep, occ_ep, one_hot_bwt_str_c_ep, match_mask_ep);
            // #else
            // GET_OCC(ep, ind, occ_id_ep, y_ep, y0_ep, y1_ep, occ_ep, bwt_str_bit0_ep, bwt_str_bit1_ep, bit0_cmp_ep, bit1_cmp_ep, mismatch_mask_ep);
            // #endif
            // #else
            // GET_OCC(ep, ind, c256, occ_id_ep, y_ep, occ_ep, bwt_str_ep, bwt_ep_vec, mask_ep_vec, mask_ep);
            // #endif
            
            ep  = count[ind] + occ_ep;

            metadata[ii].l = ep;
            metadata[ii].cur_ind--;
            if((ep <= sp) || (metadata[ii].cur_ind < 0))
            {
                k_l_range[2 * rid] = metadata[ii].k;
                k_l_range[2 * rid + 1] = metadata[ii].l;
                if(next_query < num_queries)
                {
                    BatchMetadata m;
                    int64_t ind = enc_qdb[next_query * query_length + query_length - 1];
                    m.k = count[ind];
                    m.l = count[ind + 1];
                    m.rid = next_query;
                    m.cur_ind = query_length - 2;
                    metadata[ii] = m;
                    next_query++;
                }
                else
                {
                    metadata[ii].rid = -1;
                    new_num_active--;
                }
            }
#ifdef ENABLE_PREFETCH
            _mm_prefetch((const char *)(&cp_occ[(metadata[ii].k) >> CP_SHIFT]), _MM_HINT_T0);
            _mm_prefetch((const char *)(&cp_occ[(metadata[ii].l) >> CP_SHIFT]),_MM_HINT_T0);
#if ((!__AVX2__)) && (CP_BLOCK_SIZE == 256)
            _mm_prefetch((const char *)(&cp_occ[(metadata[ii].k) >> CP_SHIFT]) + 64, _MM_HINT_T0);
            _mm_prefetch((const char *)(&cp_occ[(metadata[ii].l) >> CP_SHIFT]) + 64,_MM_HINT_T0);
#endif
#endif
        }
        if(num_active > new_num_active)
        {
            int p1 = 0;
            for(ii = 0; ii < num_active; ii++)
            {
                if(metadata[ii].rid >= 0)
                {
                    metadata[p1] = metadata[ii];
                    p1++;
                }
            }
        }
        num_active = new_num_active;
    }
    //int64_t end_thread = __rdtsc();
    //printf("%d] %ld ticks\n", tid, end_thread - start_thread);
}

void FMI_search::inexact_search(uint8_t *enc_qdb, int32_t z,
                                int64_t num_queries, uint32_t query_length,
                                int64_t *k_l_range, int32_t batch_size)
{
    AlignState *FIFO = (AlignState *)_mm_malloc(batch_size * MAX_FIFO_SIZE * sizeof(AlignState), 64);
    InexactBatchMetadata metadata[batch_size];

    int64_t ii;
    for(ii = 0; ii < batch_size; ii++)
    {
        InexactBatchMetadata m;
        m.rid = ii;

        m.cur_fifo_id = 0;
        AlignState as;
        as.k = 0;
        as.l = reference_seq_len;
        as.i = query_length - 1;
        as.z = z; 
        FIFO[ii * MAX_FIFO_SIZE] = as;
        m.fifo_size = 1;
        m.num_matches = 0;
        metadata[ii] = m;
    }
    int64_t next_read = batch_size;
    int32_t num_active = batch_size;
    while(num_active > 0)
    {
        int32_t new_num_active = num_active;
        for(ii = 0; ii < batch_size; ii++)
        {
            int64_t rid = metadata[ii].rid;
            if(rid == -1) continue;
            int32_t cur_fifo_id = metadata[ii].cur_fifo_id;
            int32_t fifo_size = metadata[ii].fifo_size;
            int32_t num_matches = metadata[ii].num_matches;
            int64_t *my_k_l_range = k_l_range +  rid * 2 * (query_length * 3 * z + 1);  

            AlignState as = FIFO[ii * MAX_FIFO_SIZE + cur_fifo_id];

            if(as.i < 0)
            {
                my_k_l_range[num_matches * 2] = as.k;
                my_k_l_range[num_matches * 2 + 1] = as.l;
                metadata[ii].num_matches++;
            }
            else
            {
                if(as.z < 0)
                {
                    printf("ERROR! rid = %ld, z = %d\n", rid, as.z);
                    exit(-1);
                }
                if(as.z == 0)
                {
                    uint8_t c = enc_qdb[rid * query_length + as.i];

                    // #if ((!__AVX2__))
                    // #if (CP_BLOCK_SIZE <= 64)
                    // GET_OCC(as.k, c, occ_id_sp, y_sp, occ_sp, bwt_str_bit0_sp, bwt_str_bit1_sp, bit0_cmp_sp, bit1_cmp_sp, mismatch_mask_sp);
                    GET_OCC(as.k, c, occ_id_sp, y_sp, occ_sp, one_hot_bwt_str_c_sp, match_mask_sp);
                    // #else
                    // GET_OCC(as.k, c, occ_id_sp, y_sp, y0_sp, y1_sp, occ_sp, bwt_str_bit0_sp, bwt_str_bit1_sp, bit0_cmp_sp, bit1_cmp_sp, mismatch_mask_sp);
                    // #endif
                    // #else
                    // __m256i c256;
                    // c256 = _mm256_set1_epi8(c);
                    // GET_OCC(as.k, c, c256, occ_id_sp, y_sp, occ_sp, bwt_str_sp, bwt_sp_vec, mask_sp_vec, mask_sp);
                    // #endif
                    
                    uint32_t sp  = count[c] + occ_sp;

                    // #if ((!__AVX2__))
                    // #if (CP_BLOCK_SIZE <= 64)
                    // GET_OCC(as.l, c, occ_id_ep, y_ep, occ_ep, bwt_str_bit0_ep, bwt_str_bit1_ep, bit0_cmp_ep, bit1_cmp_ep, mismatch_mask_ep);
                    GET_OCC(as.l, c, occ_id_ep, y_ep, occ_ep, one_hot_bwt_str_c_ep, match_mask_ep);
                    // #else
                    // GET_OCC(as.l, c, occ_id_ep, y_ep, y0_ep, y1_ep, occ_ep, bwt_str_bit0_ep, bwt_str_bit1_ep, bit0_cmp_ep, bit1_cmp_ep, mismatch_mask_ep);
                    // #endif
                    // #else
                    // GET_OCC(as.l, c, c256, occ_id_ep, y_ep, occ_ep, bwt_str_ep, bwt_ep_vec, mask_ep_vec, mask_ep);
                    // #endif
                    
                    uint32_t ep  = count[c] + occ_ep;
                    
                    if(sp < ep)
                    {
                        AlignState new_as;
                        new_as.k = sp;
                        new_as.l = ep;
                        new_as.i = as.i - 1;
                        new_as.z = as.z;
                        FIFO[ii * MAX_FIFO_SIZE + fifo_size] = new_as;
                        fifo_size++;
                        if(fifo_size > (MAX_FIFO_SIZE - 10))
                        {
                            printf("ERROR! fifo_size = %d\n", fifo_size);
                            exit(-1);
                        }
                    }
                }
                else
                {
                    uint8_t c;
                    for(c = 0; c < 4; c++)
                    {
                        // #if ((!__AVX2__))
                        // #if (CP_BLOCK_SIZE <= 64)
                        // GET_OCC(as.k, c, occ_id_sp, y_sp, occ_sp, bwt_str_bit0_sp, bwt_str_bit1_sp, bit0_cmp_sp, bit1_cmp_sp, mismatch_mask_sp);
                        GET_OCC(as.k, c, occ_id_sp, y_sp, occ_sp, one_hot_bwt_str_c_sp, match_mask_sp);
                        // #else
                        // GET_OCC(as.k, c, occ_id_sp, y_sp, y0_sp, y1_sp, occ_sp, bwt_str_bit0_sp, bwt_str_bit1_sp, bit0_cmp_sp, bit1_cmp_sp, mismatch_mask_sp);
                        // #endif
                        // #else
                        // __m256i c256;
                        // c256 = _mm256_set1_epi8(c);
                        // GET_OCC(as.k, c, c256, occ_id_sp, y_sp, occ_sp, bwt_str_sp, bwt_sp_vec, mask_sp_vec, mask_sp);
                        // #endif
                        
                        uint32_t sp  = count[c] + occ_sp;

                        // #if ((!__AVX2__))
                        // #if (CP_BLOCK_SIZE <= 64)
                        // GET_OCC(as.l, c, occ_id_ep, y_ep, occ_ep, bwt_str_bit0_ep, bwt_str_bit1_ep, bit0_cmp_ep, bit1_cmp_ep, mismatch_mask_ep);
                        GET_OCC(as.l, c, occ_id_ep, y_ep, occ_ep, one_hot_bwt_str_c_ep, match_mask_ep);
                        // #else
                        // GET_OCC(as.l, c, occ_id_ep, y_ep, y0_ep, y1_ep, occ_ep, bwt_str_bit0_ep, bwt_str_bit1_ep, bit0_cmp_ep, bit1_cmp_ep, mismatch_mask_ep);
                        // #endif
                        // #else
                        // GET_OCC(as.l, c, c256, occ_id_ep, y_ep, occ_ep, bwt_str_ep, bwt_ep_vec, mask_ep_vec, mask_ep);
                        // #endif
                        
                        uint32_t ep  = count[c] + occ_ep;

                        if(sp < ep)
                        {
                            AlignState new_as;
                            new_as.k = sp;
                            new_as.l = ep;
                            new_as.i = as.i - 1;
                            if(c == enc_qdb[rid * query_length + as.i])
                            {
                                new_as.z = as.z;
                            }
                            else
                            {
                                new_as.z = as.z - 1;
                            }
                            FIFO[ii * MAX_FIFO_SIZE + fifo_size] = new_as;
                            fifo_size++;
                            if(fifo_size > (MAX_FIFO_SIZE - 10))
                            {
                                printf("ERROR! fifo_size = %d\n", fifo_size);
                                exit(-1);
                            }
                        }

                    }
                }
            }
            metadata[ii].fifo_size = fifo_size;
            metadata[ii].cur_fifo_id++;
            if(metadata[ii].cur_fifo_id == metadata[ii].fifo_size)
            {
                if(next_read < num_queries)
                {
                    //Replace with another
                    InexactBatchMetadata m;
                    m.rid = next_read;
                    m.cur_fifo_id = 0;
                    AlignState as;
                    as.k = 0;
                    as.l = reference_seq_len;
                    as.i = query_length - 1;
                    as.z = z; 
                    FIFO[ii * MAX_FIFO_SIZE] = as;
                    m.fifo_size = 1;
                    m.num_matches = 0;
                    metadata[ii] = m;
                    next_read++;
                }
                else
                {
                    metadata[ii].rid = -1;
                    new_num_active--;
                }
            }
#ifdef ENABLE_PREFETCH
            _mm_prefetch((const char *)(&cp_occ[(FIFO[ii * MAX_FIFO_SIZE + metadata[ii].cur_fifo_id].k) >> CP_SHIFT]), _MM_HINT_T0);
            _mm_prefetch((const char *)(&cp_occ[(FIFO[ii * MAX_FIFO_SIZE + metadata[ii].cur_fifo_id].l) >> CP_SHIFT]), _MM_HINT_T0);
#if ((!__AVX2__)) && (CP_BLOCK_SIZE == 256)
            _mm_prefetch((const char *)(&cp_occ[(FIFO[ii * MAX_FIFO_SIZE + metadata[ii].cur_fifo_id].k) >> CP_SHIFT]) + 64, _MM_HINT_T0);
            _mm_prefetch((const char *)(&cp_occ[(FIFO[ii * MAX_FIFO_SIZE + metadata[ii].cur_fifo_id].l) >> CP_SHIFT]) + 64, _MM_HINT_T0);
#endif
#endif
        }

        num_active = new_num_active;
    }
}
