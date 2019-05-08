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

#ifndef FASTMAP_HPP
#define FASTMAP_HPP

#include <ctype.h>
#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>
#include <fstream>
#include "bwa.h"
#include "bwamem.h"
#include "kthread.h"
#include "kvec.h"
#include "utils.h"
#include "bntseq.h"
#include "kseq.h"
// #if MPI_ENABLED
// #include <mpi.h>
// #endif

KSEQ_DECLARE(gzFile)

typedef struct {
	kseq_t *ks, *ks2;
	mem_opt_t *opt;
	mem_pestat_t *pes0;
	int64_t n_processed;
	int copy_comment, actual_chunk_size;
	int64_t my_ntasks, ntasks, task_size;
	// bwaidx_t *idx;
	int64_t totEl;
	FILE *fp;
	// MPI_File mfp;
} ktp_aux_t;

typedef struct {
	ktp_aux_t *aux;
	int n_seqs;
	bseq1_t *seqs;
} ktp_data_t;



#define OPT_ARGS								\
	{																	\
	while ((c = getopt(argc, argv, "1paMCSPVYjk:c:v:s:r:t:R:A:B:O:E:U:w:L:d:T:Q:D:m:I:N:W:x:G:h:y:K:X:H:o:q:")) >= 0) \
	{\
		if (c == 'k') opt->min_seed_len = atoi(optarg), opt0.min_seed_len = 1;\
		else if (c == '1') no_mt_io = 1;\
		else if (c == 'x') mode = optarg;\
		else if (c == 'w') opt->w = atoi(optarg), opt0.w = 1;\
		else if (c == 'A') opt->a = atoi(optarg), opt0.a = 1, assert(opt->a >= INT_MIN && opt->a <= INT_MAX); \
		else if (c == 'B') opt->b = atoi(optarg), opt0.b = 1, assert(opt->b >= INT_MIN && opt->b <= INT_MAX); \
		else if (c == 'T') opt->T = atoi(optarg), opt0.T = 1, assert(opt->T >= INT_MIN && opt->T <= INT_MAX);\
		else if (c == 'U')\
			opt->pen_unpaired = atoi(optarg), opt0.pen_unpaired = 1, assert(opt->pen_unpaired >= INT_MIN && opt->pen_unpaired <= INT_MAX);\
		else if (c == 't')\
			opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1, assert(opt->n_threads >= INT_MIN && opt->n_threads <= INT_MAX);\
		else if (c == 'o')\
		{\
			is_o = 1;\
			aux.fp = fopen(optarg, "w");								\
			if (aux.fp == NULL) {										\
				fprintf(stderr, "Error: can't open %s input file\n", optarg);	\
				exit(0);												\
			}															\
			/*fclose(aux.fp);*/											\
			/*MPI_File_open(MPI_COMM_WORLD, optarg, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &aux.mfp);*/ \
			aux.totEl = 0;\
		}\
		else if (c == 'q') {\
			nread_lim = atoi(optarg);\
		}\
		else if (c == 'P') opt->flag |= MEM_F_NOPAIRING;	\
		else if (c == 'a') opt->flag |= MEM_F_ALL;\
		else if (c == 'p') opt->flag |= MEM_F_PE | MEM_F_SMARTPE;\
		else if (c == 'M') opt->flag |= MEM_F_NO_MULTI;\
		else if (c == 'S') opt->flag |= MEM_F_NO_RESCUE;\
		else if (c == 'Y') opt->flag |= MEM_F_SOFTCLIP;\
		else if (c == 'V') opt->flag |= MEM_F_REF_HDR;\
		else if (c == 'c') opt->max_occ = atoi(optarg), opt0.max_occ = 1;\
		else if (c == 'd') opt->zdrop = atoi(optarg), opt0.zdrop = 1;\
		else if (c == 'v') bwa_verbose = atoi(optarg);\
		else if (c == 'j') ignore_alt = 1;\
		else if (c == 'r')\
			opt->split_factor = atof(optarg), opt0.split_factor = 1.;\
		else if (c == 'D') opt->drop_ratio = atof(optarg), opt0.drop_ratio = 1.;\
		else if (c == 'm') opt->max_matesw = atoi(optarg), opt0.max_matesw = 1;\
		else if (c == 's') opt->split_width = atoi(optarg), opt0.split_width = 1;\
		else if (c == 'G')\
			opt->max_chain_gap = atoi(optarg), opt0.max_chain_gap = 1;\
		else if (c == 'N')\
			opt->max_chain_extend = atoi(optarg), opt0.max_chain_extend = 1;\
		else if (c == 'W')												\
			opt->min_chain_weight = atoi(optarg), opt0.min_chain_weight = 1; \
		else if (c == 'y')												\
			opt->max_mem_intv = atol(optarg), opt0.max_mem_intv = 1;	\
		else if (c == 'C') aux.copy_comment = 1;						\
		else if (c == 'K') fixed_chunk_size = atoi(optarg);				\
		else if (c == 'X') opt->mask_level = atof(optarg);				\
		else if (c == 'h')												\
		{																\
			opt0.max_XA_hits = opt0.max_XA_hits_alt = 1;				\
			opt->max_XA_hits = opt->max_XA_hits_alt = strtol(optarg, &p, 10);\
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))\
				opt->max_XA_hits_alt = strtol(p+1, &p, 10);\
		}\
		else if (c == 'Q')\
		{\
			opt0.mapQ_coef_len = 1;\
			opt->mapQ_coef_len = atoi(optarg);\
			opt->mapQ_coef_fac = opt->mapQ_coef_len > 0? log(opt->mapQ_coef_len) : 0;\
		}\
		else if (c == 'O')\
		{\
			opt0.o_del = opt0.o_ins = 1;\
			opt->o_del = opt->o_ins = strtol(optarg, &p, 10);\
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))\
				opt->o_ins = strtol(p+1, &p, 10);\
		}\
		else if (c == 'E')\
		{\
			opt0.e_del = opt0.e_ins = 1;\
			opt->e_del = opt->e_ins = strtol(optarg, &p, 10);\
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))\
				opt->e_ins = strtol(p+1, &p, 10);\
		}\
		else if (c == 'L')\
		{\
			opt0.pen_clip5 = opt0.pen_clip3 = 1;\
			opt->pen_clip5 = opt->pen_clip3 = strtol(optarg, &p, 10);\
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))\
				opt->pen_clip3 = strtol(p+1, &p, 10);\
		}\
		else if (c == 'R')\
		{															   \
		if ((rg_line = bwa_set_rg(optarg)) == 0) {					   \
			free(opt);												   \
			if (is_o)												   \
				fclose(aux.fp); \
				return 1;\
            }\
		}\
		else if (c == 'H')\
		{\
			if (optarg[0] != '@')\
			{\
				FILE *fp;\
				if ((fp = fopen(optarg, "r")) != 0)\
				{\
					char *buf;\
					buf = (char *) calloc(1, 0x10000);\
					while (fgets(buf, 0xffff, fp))\
					{\
						i = strlen(buf);\
						assert(buf[i-1] == '\n'); \
						buf[i-1] = 0;\
						hdr_line = bwa_insert_header(buf, hdr_line);\
					}\
					free(buf);\
					fclose(fp);\
				}\
			} else hdr_line = bwa_insert_header(optarg, hdr_line);\
		}\
		else if (c == 'I')\
		{\
			aux.pes0 = pes;\
			pes[1].failed = 0;\
			pes[1].avg = strtod(optarg, &p);\
			pes[1].std = pes[1].avg * .1;\
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))\
				pes[1].std = strtod(p+1, &p);\
			pes[1].high = (int)(pes[1].avg + 4. * pes[1].std + .499);\
			pes[1].low  = (int)(pes[1].avg - 4. * pes[1].std + .499);\
			if (pes[1].low < 1) pes[1].low = 1;\
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))\
				pes[1].high = (int)(strtod(p+1, &p) + .499);\
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))\
				pes[1].low  = (int)(strtod(p+1, &p) + .499);\
		}\
		else {									\
			free(opt);												   \
			if (is_o)												   \
				fclose(aux.fp); \
			return 1;							\
		}\
	}\
	}

#define MODE									\
	if (mode)\
	{\
		if (strcmp(mode, "intractg") == 0)\
		{\
			if (!opt0.o_del) opt->o_del = 16;\
			if (!opt0.o_ins) opt->o_ins = 16;\
			if (!opt0.b) opt->b = 9;\
			if (!opt0.pen_clip5) opt->pen_clip5 = 5;\
			if (!opt0.pen_clip3) opt->pen_clip3 = 5;\
		}\
		else if (strcmp(mode, "pacbio") == 0 || strcmp(mode, "pbref") == 0 || strcmp(mode, "ont2d") == 0)\
		{\
			if (!opt0.o_del) opt->o_del = 1;\
			if (!opt0.e_del) opt->e_del = 1;\
			if (!opt0.o_ins) opt->o_ins = 1;\
			if (!opt0.e_ins) opt->e_ins = 1;\
			if (!opt0.b) opt->b = 1;\
			if (opt0.split_factor == 0.) opt->split_factor = 10.;\
			if (strcmp(mode, "ont2d") == 0)\
			{\
				if (!opt0.min_chain_weight) opt->min_chain_weight = 20;\
				if (!opt0.min_seed_len) opt->min_seed_len = 14;\
				if (!opt0.pen_clip5) opt->pen_clip5 = 0;\
				if (!opt0.pen_clip3) opt->pen_clip3 = 0;\
			}\
			else\
			{\
				if (!opt0.min_chain_weight) opt->min_chain_weight = 40;\
				if (!opt0.min_seed_len) opt->min_seed_len = 17;\
				if (!opt0.pen_clip5) opt->pen_clip5 = 0;\
				if (!opt0.pen_clip3) opt->pen_clip3 = 0;\
			}\
		}\
		else\
		{\
			fprintf(stderr, "[E::%s] unknown read type '%s'\n", __func__, mode);\
			free(opt);												   \
			if (is_o)												   \
				fclose(aux.fp); \
			return 1;						\
        }\
	} else update_a(opt, &opt0);

#endif
