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
         Heng Li <hli@jimmy.harvard.edu>
*****************************************************************************************/

#ifndef BWA_H_
#define BWA_H_

#include <stdint.h>
#include "bntseq.h"
#include "bwt.h"
#include "macro.h"

#define BWA_IDX_BWT 0x1
#define BWA_IDX_BNS 0x2
#define BWA_IDX_PAC 0x4
#define BWA_IDX_ALL 0x7

#define BWA_CTL_SIZE 0x10000

#define BWTALGO_AUTO  0
#define BWTALGO_RB2   1
#define BWTALGO_BWTSW 2
#define BWTALGO_IS    3

typedef struct {
	// bwt2_t   *bwt2;
	bwt_t    *bwt; // FM-index
	bntseq_t *bns; // information on the reference sequences
	uint8_t  *pac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base

	int    is_shm;
	int64_t l_mem;
	uint8_t  *mem;
} bwaidx_t;


typedef struct {
	int l_seq, id;
	char *name, *comment, *seq, *qual, *sam;
} bseq1_t;

extern int bwa_verbose;
extern char bwa_rg_id[256];

bseq1_t *bseq_read_orig(int chunk_size, int *n_, void *ks1_, void *ks2_, int64_t *s);

bseq1_t *bseq_read(int chunk_size, int *n_, void *ks1_,
				   void *ks2_, FILE* fpp, int len,
				   int64_t *sz);

void bseq_classify(int n, bseq1_t *seqs, int m[2], bseq1_t *sep[2]);

	void bwa_fill_scmat(int a, int b, int8_t mat[25]);

	uint32_t *bwa_gen_cigar(const int8_t mat[25], int q, int r, int w_,
							int64_t l_pac, const uint8_t *pac, int l_query,
							uint8_t *query, int64_t rb, int64_t re,
							int *score, int *n_cigar, int *NM);
	
	uint32_t *bwa_gen_cigar2(const int8_t mat[25], int o_del, int e_del,
							 int o_ins, int e_ins, int w_, int64_t l_pac,
							 const uint8_t *pac, int l_query, uint8_t *query,
							 int64_t rb, int64_t re, int *score,
							 int *n_cigar, int *NM);

	int bwa_idx_build(const char *fa, const char *prefix,
					  int algo_type, int block_size);

	char *bwa_idx_infer_prefix(const char *hint);
	bwt_t *bwa_idx_load_bwt(const char *hint);
	bwt2_t *bwa_idx_load_bwt2(const char *hint);
	
	bwaidx_t *bwa_idx_load_from_shm(const char *hint);
	bwaidx_t *bwa_idx_load_from_disk(const char *hint, int which);
	bwaidx_t *bwa_idx_load(const char *hint, int which);
	
	void bwa_idx_destroy(bwaidx_t *idx);
	void bwa_print_sam_hdr(const bntseq_t *bns, const char *hdr_line, FILE *fp);
	char *bwa_set_rg(const char *s);
	char *bwa_insert_header(const char *s, char *hdr);

#endif
