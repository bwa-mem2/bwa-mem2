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

#include <stdio.h>
#include "fastmap.h"
#include "FMI_search.h"
#include "fasta_file.h"

// --------------
// global vars
FMI_search *fmi;
extern uint64_t proc_freq, tprof[LIM_R][LIM_C];
extern unsigned char nst_nt4_table[256];
extern int num_ranks, myrank;
extern int nthreads;
uint8_t *ref_string;
int readLen;
int64_t nreads;
// ---------------
int64_t get_limit_fsize(FILE *fpp, int64_t nread_lim,
						char buf[], char buf1[]) {
	
	int64_t val = 0, len = 10000;
	fseek(fpp, 0, SEEK_END);
	val = ftell(fpp);
	if (nread_lim >= val)
		return val;

	fseek(fpp, nread_lim, SEEK_SET);
	int64_t position = nread_lim;
	while(true) {
		if (fgets((char*) buf, len, fpp) != NULL) {
			if (buf[0] == '@') {
				int64_t pos = ftell(fpp);
				fgets((char*) buf1, len, fpp);
				fgets((char*) buf1, len, fpp);
				if (buf1[0] == '+')
					break;
				fseek(fpp, pos, SEEK_SET);
			}
			position = ftell(fpp);
		}
		else break;
	}
	return position;	
}

void memoryAlloc(ktp_aux_t *aux, worker_t &w, int64_t nreads)
{
	mem_opt_t	*opt			  = aux->opt;	
	int memSize = nreads;

	/* Mem allocation section for core kernels */
	w.regs = NULL; w.chain_ar = NULL; w.seedBuf = NULL;
	w.regs = (mem_alnreg_v *) calloc(memSize, sizeof(mem_alnreg_v));
	w.chain_ar = (mem_chain_v*) malloc (memSize * sizeof(mem_chain_v));
    w.seedBuf = (mem_seed_t *)_mm_malloc(sizeof(mem_seed_t) * memSize * AVG_SEEDS_PER_READ, 64);
	assert(w.seedBuf != NULL);
    w.seedBufSize = BATCH_SIZE * AVG_SEEDS_PER_READ;
    w.auxSeedBuf = (mem_seed_t *)_mm_malloc(sizeof(mem_seed_t) * memSize * AVG_AUX_SEEDS_PER_READ, 64);
    w.auxSeedBufSize = BATCH_SIZE * AVG_AUX_SEEDS_PER_READ;
	
	if (w.regs == NULL || w.chain_ar == NULL) {
		fprintf(stderr, "Memory not allocated!!\nExiting...\n");
		exit(0);
	}

	int64_t allocMem = memSize * sizeof(mem_alnreg_v) +
		memSize * sizeof(mem_chain_v) +
		sizeof(mem_seed_t) * memSize * AVG_SEEDS_PER_READ +
		sizeof(mem_seed_t) * memSize * AVG_AUX_SEEDS_PER_READ;
	fprintf(stderr, "------------------------------------------\n");
	fprintf(stderr, "Memory pre-allocation for chaining: %ld\n", allocMem);

	
	/* SWA mem allocation */
	// int avg_seed_per_read = 35;
	w.size = BATCH_SIZE * SEEDS_PER_READ;
	
	w.mmc.seqBufLeftRef = NULL;	w.mmc.seqBufLeftQer = NULL;
	w.mmc.seqBufRightRef = NULL; w.mmc.seqBufRightQer = NULL;
	w.mmc.seqPairArrayLeft = NULL; 	w.mmc.seqPairArrayAux = NULL;
	w.mmc.seqPairArrayLeft128 = NULL; w.mmc.seqPairArrayRight128 = NULL;
	
	w.mmc.seqBufLeftRef = (uint8_t *)_mm_malloc((w.size * MAX_SEQ_LEN_REF * sizeof(int8_t)
												 + MAX_LINE_LEN) * opt->n_threads, 64);
	assert(w.mmc.seqBufLeftRef != NULL);

	w.mmc.seqBufLeftQer = (uint8_t *)_mm_malloc((w.size * MAX_SEQ_LEN_QER * sizeof(int8_t)
												 + MAX_LINE_LEN) * opt->n_threads, 64);	   
	assert(w.mmc.seqBufLeftQer != NULL);
	

	w.mmc.seqBufRightRef = (uint8_t *)_mm_malloc((w.size * MAX_SEQ_LEN_REF * sizeof(int8_t)
												  + MAX_LINE_LEN) * opt->n_threads, 64);
	w.mmc.seqBufRightQer = (uint8_t *)_mm_malloc((w.size * MAX_SEQ_LEN_QER * sizeof(int8_t)
												  + MAX_LINE_LEN) * opt->n_threads, 64);
	
	//w.mmc.seqPairArrayLeft = (SeqPair *)_mm_malloc((w.size * sizeof(SeqPair)
	//												   ) * opt->n_threads, 64);
	//w.mmc.seqPairArrayRight = (SeqPair *)_mm_malloc((w.size * sizeof(SeqPair)
	//						 ) * opt->n_threads, 64);
	w.mmc.seqPairArrayAux = (SeqPair *)_mm_malloc((w.size * sizeof(SeqPair)
						       ) * opt->n_threads, 64);
	w.mmc.seqPairArrayLeft128 = (SeqPair *)_mm_malloc((w.size * sizeof(SeqPair)
							   ) * opt->n_threads, 64);
	w.mmc.seqPairArrayRight128 = (SeqPair *)_mm_malloc((w.size * sizeof(SeqPair)
							    ) * opt->n_threads, 64);

	assert(w.mmc.seqPairArrayRight128 != NULL);

	allocMem = (w.size * MAX_SEQ_LEN_REF * sizeof(int8_t) + MAX_LINE_LEN) * opt->n_threads * 2+
		(w.size * MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_LINE_LEN) * opt->n_threads	* 2 +		
		w.size * sizeof(SeqPair) * opt->n_threads * 3;
	
	fprintf(stderr, "Memory pre-allocation for BSW: %ld\n", allocMem);
	
	w.mmc.matchArray = (SMEM *)_mm_malloc
		(nthreads * BATCH_MUL * BATCH_SIZE * readLen * sizeof(SMEM), 64);
	w.mmc.min_intv_ar = (int32_t *) malloc
		(nthreads * BATCH_MUL * BATCH_SIZE * readLen *sizeof(int32_t));
	w.mmc.query_pos_ar = (int16_t *) malloc
		(nthreads * BATCH_MUL * BATCH_SIZE * readLen *sizeof(int16_t));
	w.mmc.enc_qdb = (uint8_t *) malloc
		(nthreads * BATCH_MUL * BATCH_SIZE * readLen *sizeof(uint8_t));
	w.mmc.rid = (int32_t *) malloc
		(nthreads * BATCH_MUL * BATCH_SIZE * readLen *sizeof(int32_t));
	w.mmc.lim = (int32_t *) _mm_malloc
		(nthreads * (BATCH_SIZE + 32) * sizeof(int32_t), 64);

	allocMem = nthreads * BATCH_MUL * BATCH_SIZE * readLen * sizeof(SMEM) +
		nthreads * BATCH_MUL * BATCH_SIZE * readLen *sizeof(int32_t) +
		nthreads * BATCH_MUL * BATCH_SIZE * readLen *sizeof(int16_t) +
		nthreads * BATCH_MUL * BATCH_SIZE * readLen *sizeof(int32_t) +
		nthreads * (BATCH_SIZE + 32) * sizeof(int32_t);
	fprintf(stderr, "Memory pre-allocation for BWT: %ld\n", allocMem);
	fprintf(stderr, "------------------------------------------\n");
}

ktp_data_t *kt_pipeline(void *shared, int step, void *data, mem_opt_t *opt, worker_t &w)
{
	ktp_aux_t *aux = (ktp_aux_t*) shared;
	ktp_data_t *ret = (ktp_data_t*) data;

	if (step == 0)
	{
		// printf("Thread entering step 0, with task_size: %d, CPU: %d\n", aux->task_size, sched_getcpu());
		ktp_data_t *ret = (ktp_data_t *) calloc(1, sizeof(ktp_data_t));
		uint64_t tim = __rdtsc();
		/* Read "reads" from input file (fread) */
		int64_t sz = 0;
		//if (temp++ == 0) 
		ret->seqs = bseq_read_orig(aux->task_size,
								   &ret->n_seqs,
								   aux->ks,
								   aux->ks2,
								   &sz);

		// *ret = *ret2;
		tprof[READ_IO][0] += __rdtsc() - tim;
		tprof[0][0] += sz;  // debug info, for accuracy checks!!
		tprof[0][2] += ret->n_seqs;
		// iteration ++;
		// fprintf(stderr, "Read %d seqs, task_size: %ld\n", ret->n_seqs, aux->task_size);
		assert(ret->n_seqs <= nreads);
		
		fprintf(stderr, "[%0.4d] read_chunk: %ld, work_chunk_size: %ld, nseq: %d\n",
				myrank, aux->task_size, sz, ret->n_seqs);   

		if (ret->seqs == 0) {
			free(ret);
			// printf("Threads leaving step 0 (Break), CPU: %d\n", sched_getcpu());
			return 0;
		}
		if (!aux->copy_comment){
			for (int i = 0; i < ret->n_seqs; ++i){
				free(ret->seqs[i].comment);
				ret->seqs[i].comment = 0;
			}
		}
		{
			int64_t size = 0;
			for (int i = 0; i < ret->n_seqs; ++i) size += ret->seqs[i].l_seq;

			fprintf(stderr, "\t[%0.4d][ M::%s] read %d sequences (%ld bp)...\n",
					myrank, __func__, ret->n_seqs, (long)size);
		}
				
		// printf("Threads leaving step 0, CPU: %d\n", sched_getcpu());
		return ret;
	} // Step 0			
	else if (step == 1)  /* Step 2: Main processing-engine */
	{
		static int task = 0;
		// printf("Thread entering step 1, CPU: %d\n", sched_getcpu());
								
		fprintf(stderr, "[%0.4d] 2. Calling mem_process_seqs.., task: %d\n", myrank, task++);

		uint64_t tim = __rdtsc();
		if (opt->flag & MEM_F_SMARTPE)
		{
			bseq1_t *sep[2];
			int n_sep[2];
			mem_opt_t tmp_opt = *opt;

			bseq_classify(ret->n_seqs, ret->seqs, n_sep, sep);

			fprintf(stderr, "[M::%s] %d single-end sequences; %d paired-end sequences.....\n",
					__func__, n_sep[0], n_sep[1]);
			
			if (n_sep[0]) {
				tmp_opt.flag &= ~MEM_F_PE;
				/* single-end sequences, in the mixture */
				mem_process_seqs(&tmp_opt,
								 fmi->idx->bns,
								 fmi->idx->pac,
								 aux->n_processed,
								 n_sep[0],
								 sep[0],
								 0,
								 w);
				// mem_process_seqs(&tmp_opt, idx->bwt, idx->bns, idx->pac,
				//				 aux->n_processed, n_sep[0], sep[0], 0);
				
				for (int i = 0; i < n_sep[0]; ++i)
					ret->seqs[sep[0][i].id].sam = sep[0][i].sam;
			}
			if (n_sep[1]) {
				tmp_opt.flag |= MEM_F_PE;
				/* paired-end sequences, in the mixture */
				mem_process_seqs(&tmp_opt,
								 fmi->idx->bns,
								 fmi->idx->pac,
								 aux->n_processed + n_sep[0],
								 n_sep[1],
								 sep[1],
								 aux->pes0,
								 w);
				
				// mem_process_seqs(&tmp_opt, idx->bwt, idx->bns, idx->pac,
				// aux->n_processed + n_sep[0], n_sep[1], sep[1], aux->pes0);
				
				for (int i = 0; i < n_sep[1]; ++i)
					ret->seqs[sep[1][i].id].sam = sep[1][i].sam;
			}
			free(sep[0]); free(sep[1]);
		}
		else {
			/* pure (single/paired-end), reads processing */
			mem_process_seqs(opt,
							 fmi->idx->bns,
							 fmi->idx->pac,
							 aux->n_processed,
							 ret->n_seqs,
							 ret->seqs,
							 aux->pes0,
							 w);
		}				
		tprof[MEM_PROCESS2][0] += __rdtsc() - tim;
				
		// printf("Thread leaving step 1, CPU: %d\n", sched_getcpu());
		return ret;
	} 			
	/* Step 3: Write output */
	else if (step == 2) {
		// printf("Thread entering step 2, CPU: %d\n", sched_getcpu());
				
		aux->n_processed += ret->n_seqs;   //modified!!
		uint64_t tim = __rdtsc();

		for (int i = 0; i < ret->n_seqs; ++i)
		{
			if (ret->seqs[i].sam) {
				// err_fputs(ret->seqs[i].sam, stderr);
				fputs(ret->seqs[i].sam, aux->fp);
			}
			free(ret->seqs[i].name); free(ret->seqs[i].comment);
		    free(ret->seqs[i].seq); free(ret->seqs[i].qual);
			free(ret->seqs[i].sam);
		}
		free(ret->seqs);
		tprof[SAM_IO][0] += __rdtsc() - tim;

		// printf("Thread leaving step 2, CPU: %d\n", sched_getcpu());
		// printf("[%0.4d] mem usage: %ld %ld\n", myrank, getCurrentRSS(), getPeakRSS());
		return 0;
	} // step 2

	return 0;
}

static void *ktp_worker(void *data) {
	
	// printf("ptid: %d\n", *(int*)ptid);
	ktp_worker_t *w = (ktp_worker_t*) data;
	ktp_t *p = w->pl;
	// args_t args = w->args;
	
	while (w->step < p->n_steps) {
		// test whether we can kick off the job with this worker
		pthread_mutex_lock(&p->mutex);
		for (;;) {
			int i;
			// test whether another worker is doing the same step
			for (i = 0; i < p->n_workers; ++i) {
				if (w == &p->workers[i]) continue; // ignore itself
				if (p->workers[i].step <= w->step && p->workers[i].index < w->index)
					break;
			}
			if (i == p->n_workers) break; // no workers with smaller indices are doing w->step or the previous steps
			pthread_cond_wait(&p->cv, &p->mutex);
		}
		pthread_mutex_unlock(&p->mutex);

		// working on w->step
		// printf("[%0.4d] Calling pipeline.., step: %d\n", w->i, w->step);

		w->data = kt_pipeline(p->shared, w->step, w->step? w->data : 0, w->opt, *(w->w)); // for the first step, input is NULL
		// if (w->data == 0) printf("NULL data recv from step %d by thread %d\n", w->step, w->i);

		// update step and let other workers know
		pthread_mutex_lock(&p->mutex);
		w->step = w->step == p->n_steps - 1 || w->data? (w->step + 1) % p->n_steps : p->n_steps;

		// printf("Step: %d, %d\n", w->step, w->data? (w->step + 1) % p->n_steps : p->n_steps);

		if (w->step == 0) w->index = p->index++;
		pthread_cond_broadcast(&p->cv);
		pthread_mutex_unlock(&p->mutex);
	}
	pthread_exit(0);
}

static int process(void *shared, gzFile gfp, gzFile gfp2, int pipe_threads)
{
	ktp_aux_t	*aux			  = (ktp_aux_t*) shared;
	ktp_data_t	*ret;
	int			 i, nt, iteration = 0;
	int64_t		 sz				  = 0;
	worker_t	 w;
	mem_opt_t	*opt			  = aux->opt;

	nthreads = opt->n_threads; // global variable for profiling!
	fprintf(stderr, "Threads used (compute): %d\n",
			nthreads);
	
	nreads = aux->actual_chunk_size/ (readLen) + 10;
	// fprintf(stderr, "Input read length: %d\n", readLen);
	fprintf(stderr, "Projected #read in a task: %d\n", nreads);
	
	/* All memory allocation */
	memoryAlloc(aux, w, nreads);
	
	int task = 0;
	/* pipeline using pthreads */
	ktp_t aux_;
	int p_nt = pipe_threads; // 2;
	int n_steps = 3;
	
	//if (n_threads < 1) n_threads = 1;
	aux_.n_workers = p_nt;
	aux_.n_steps = n_steps;
	// aux_.func = process;
	aux_.shared = aux;
	aux_.index = 0;
	pthread_mutex_init(&aux_.mutex, 0);
	pthread_cond_init(&aux_.cv, 0);
	fprintf(stderr, "No. of pipeline threads: %d\n", p_nt);
	aux_.workers = (ktp_worker_t*) malloc(p_nt * sizeof(ktp_worker_t));
	
	for (int i = 0; i < p_nt; ++i) {
		ktp_worker_t *wr = &aux_.workers[i];
		wr->step = 0; wr->pl = &aux_; wr->data = 0;
		wr->index = aux_.index++;
		wr->i = i;
		wr->opt = opt;
		wr->w = &w;
	}
	
	pthread_t *ptid = (pthread_t *) calloc(p_nt, sizeof(pthread_t));
	
	for (int i = 0; i < p_nt; ++i)
		pthread_create(&ptid[i], 0, ktp_worker, (void*) &aux_.workers[i]);
	
	for (int i = 0; i < p_nt; ++i)
		pthread_join(ptid[i], 0);
	
	free(ptid);
	free(aux_.workers);
	/***** pipeline ends ******/
	
	fprintf(stderr, "[%0.4d] Computation ends..\n", myrank);
	
	/* Dealloc memory allcoated in the header section */	
	free(w.chain_ar);
	free(w.regs);
    _mm_free(w.seedBuf);
    _mm_free(w.auxSeedBuf);
	_mm_free(w.mmc.seqBufLeftRef);
	_mm_free(w.mmc.seqBufRightRef);
	_mm_free(w.mmc.seqBufLeftQer);
	_mm_free(w.mmc.seqBufRightQer);
	
	_mm_free(w.mmc.seqPairArrayAux);
	_mm_free(w.mmc.seqPairArrayLeft128);
	_mm_free(w.mmc.seqPairArrayRight128);
	
	_mm_free(w.mmc.matchArray);
	free(w.mmc.min_intv_ar);
	free(w.mmc.query_pos_ar);
	free(w.mmc.enc_qdb);
	free(w.mmc.rid);
	_mm_free(w.mmc.lim);
	return 0;
}

static void update_a(mem_opt_t *opt, const mem_opt_t *opt0)
{
	if (opt0->a) { // matching score is changed
		if (!opt0->b) opt->b *= opt->a;
		if (!opt0->T) opt->T *= opt->a;
		if (!opt0->o_del) opt->o_del *= opt->a;
		if (!opt0->e_del) opt->e_del *= opt->a;
		if (!opt0->o_ins) opt->o_ins *= opt->a;
		if (!opt0->e_ins) opt->e_ins *= opt->a;
		if (!opt0->zdrop) opt->zdrop *= opt->a;
		if (!opt0->pen_clip5) opt->pen_clip5 *= opt->a;
		if (!opt0->pen_clip3) opt->pen_clip3 *= opt->a;
		if (!opt0->pen_unpaired) opt->pen_unpaired *= opt->a;
	}
}

int main_mem(int argc, char *argv[])
{
	int			 i, c, ignore_alt = 0, no_mt_io = 0;
	int			 fixed_chunk_size		   = -1;
	char		*p, *rg_line			   = 0, *hdr_line = 0;
	const char	*mode					   = 0;
	
	mem_opt_t		*opt, opt0;
	gzFile			 fp, fp2   = 0;
	mem_pestat_t	 pes[4];
	ktp_aux_t		 aux;
	bool			 is_o	   = 0;
	int64_t			 nread_lim = 0;

	// uint64_t tim = __rdtsc();
	memset(&aux, 0, sizeof(ktp_aux_t));
	memset(pes, 0, 4 * sizeof(mem_pestat_t));
	for (i = 0; i < 4; ++i) pes[i].failed = 1;
	
	// opterr = 0;
	aux.fp = stdout;
	aux.opt = opt = mem_opt_init();
	memset(&opt0, 0, sizeof(mem_opt_t));
	
	/* Parse input arguments */
	OPT_ARGS;
	
	/* Check output file name */
	if (rg_line)
	{
		hdr_line = bwa_insert_header(rg_line, hdr_line);
		free(rg_line);
	}

	if (opt->n_threads < 1) opt->n_threads = 1;
	if (optind + 2 != argc && optind + 3 != argc) {
		fprintf(stderr, "optind: %d, argc: %d\n", optind, argc);
		usage();
		free(opt);
		if (is_o) 
			fclose(aux.fp);
		return 1;
	}

	/* Further input parsing */
	MODE;
	
	/* Matrix for SWA */
	bwa_fill_scmat(opt->a, opt->b, opt->mat);
	// tprof[PREPROCESS][0] += __rdtsc() - tim;
	
	/* Load bwt2/FMI index */
	{
		uint64_t tim = __rdtsc();

		fprintf(stderr, "Ref file: %s\n", argv[optind]);			
		fmi = new FMI_search(argv[optind]);
		tprof[FMI][0] += __rdtsc() - tim;
		
		// reading ref string from the file
		tim = __rdtsc();
		fprintf(stderr, "Reading reference genome..\n");
		
        char binary_seq_file[200];
        sprintf(binary_seq_file, "%s.0123", argv[optind]);
		
		fprintf(stderr, "Binary seq file = %s\n", binary_seq_file);
		FILE *fr = fopen(binary_seq_file, "r");
		
		if (fr == NULL) {
			fprintf(stderr, "Error: can't open %s input file\n", binary_seq_file);
			exit(0);
		}
		
		int64_t rlen = 0;
		fseek(fr, 0, SEEK_END); 
		rlen = ftell(fr);
		ref_string = (uint8_t*) _mm_malloc(rlen, 64);
		rewind(fr);

		/* Reading ref. sequence */
		fread(ref_string, 1, rlen, fr);

		uint64_t timer  = __rdtsc();
		tprof[REF_IO][0] += timer - tim;
		
		fclose(fr);
		fprintf(stderr, "Reference genome size: %ld bp\n", rlen);
		fprintf(stderr, "Done readng reference genome !!\n\n");
	}

	if (ignore_alt)
		for (i = 0; i < fmi->idx->bns->n_seqs; ++i)
			fmi->idx->bns->anns[i].is_alt = 0;

	/* READS file operations */
	fp = gzopen(argv[optind + 1], "r");
	if (fp == 0)
	{
		fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 1]);
		free(opt);
		if (is_o) 
			fclose(aux.fp);
		return 1;
	}
	aux.ks = kseq_init(fp);
	
#if PAIRED_END
	/* Handling Paired-end reads */
	aux.ks2 = 0;
	if (optind + 2 < argc) {
		if (opt->flag & MEM_F_PE) {
			fprintf(stderr, "[W::%s] when '-p' is in use, the second query file is ignored.\n", __func__);
		}
		else
		{
			fp2 = gzopen(argv[optind + 2], "r");
			if (fp2 == 0) {
				fprintf(stderr, "[E::%s] failed to open file `%s'.\n", __func__, argv[optind + 2]);
				free(opt);
				err_gzclose(fp);
				kseq_destroy(aux.ks);
				if (is_o) 
					fclose(aux.fp);				
				return 1;
			}
			aux.ks2 = kseq_init(fp2);
			opt->flag |= MEM_F_PE;
			assert(aux.ks2 != 0);
		}
	}
#endif

	bwa_print_sam_hdr(fmi->idx->bns, hdr_line, aux.fp);

	// aux.totEl += ftell(aux.fp);
	// aux.actual_chunk_size = fixed_chunk_size > 0? fixed_chunk_size : opt->chunk_size * opt->n_threads; // IMP modification

	if (fixed_chunk_size > 0)
		aux.task_size = fixed_chunk_size;
	else {
		//aux.task_size = 10000000 * opt->n_threads; //aux.actual_chunk_size;
		aux.task_size = opt->chunk_size * opt->n_threads; //aux.actual_chunk_size;
		//aux.task_size = 11387500 * opt->n_threads; //aux.actual_chunk_size;
		//aux.task_size =   50000000;  //50 MB task_size
	}
	tprof[MISC][1] = opt->chunk_size = aux.actual_chunk_size = aux.task_size;

	{ // to find read length
		int dummy_task_size = 4000000, ns = 0;
		int64_t sz = 0;
		bseq1_t *seqs = bseq_read_orig(dummy_task_size,
									   &ns,
									   aux.ks,
									   NULL,
									   &sz);
		assert(seqs != NULL);
		readLen = seqs[0].l_seq;
		fprintf(stderr, "readLen: %d\n\n", readLen);
		for (int i = 0; i < ns; ++i) {
			free(seqs[i].name); free(seqs[i].comment);
			free(seqs[i].seq); free(seqs[i].qual);
		}
		//kseq_rewind(aux.ks);
		gzrewind(fp);
		kseq_destroy(aux.ks);
		aux.ks = kseq_init(fp);
	}
	
	// Major function
	fprintf(stderr, "[%0.4d] 1: Calling process()\n", myrank);

	uint64_t tim = __rdtsc();
	/* Relay process function */
	process(&aux, fp, fp2, no_mt_io? 1:2);
	
	tprof[PROCESS][0] += __rdtsc() - tim;

	// free memory
	_mm_free(ref_string);
	free(hdr_line);
	free(opt);
	kseq_destroy(aux.ks);	
	err_gzclose(fp);

#if PAIRED_END
	if (aux.ks2) {
		kseq_destroy(aux.ks2);
		err_gzclose(fp2);
	}
#endif
	
	if (is_o) {
		fclose(aux.fp);
// #if MPI_ENABLED
// 		MPI_File_close(&(aux.mfp));
// #endif
	}

	// new bwt/FMI
	delete(fmi);	
	
	return 0;
}

