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
         Heng Li <hli@jimmy.harvard.edu>
*****************************************************************************************/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if NUMA_ENABLED
#include <numa.h>
#endif
#include <sstream>
#include "fastmap.h"
#include "FMI_search.h"
// #include "fasta_file.h"
#if 0 && (__linux__)
#include <sys/sysinfo.h>
#endif

// --------------
// global vars
FMI_search *fmi;
extern uint64_t proc_freq, tprof[LIM_R][LIM_C];
extern unsigned char nst_nt4_table[256];
extern int num_ranks, myrank;
extern int nthreads;
uint8_t *ref_string;
int readLen, affy[256];
int64_t nreads, memSize;
// ---------------
void __cpuid(unsigned int i, unsigned int cpuid[4]) {
#ifdef _WIN32
  __cpuid((int *) cpuid, (int)i);

#else
  asm volatile
    ("cpuid" : "=a" (cpuid[0]), "=b" (cpuid[1]), "=c" (cpuid[2]), "=d" (cpuid[3])
     : "0" (i), "2" (0));
#endif
}


int HTStatus()
{
  unsigned int cpuid[4];
  char platform_vendor[12];
  __cpuid(0, cpuid);
  ((unsigned int *)platform_vendor)[0] = cpuid[1]; // B
  ((unsigned int *)platform_vendor)[1] = cpuid[3]; // D
  ((unsigned int *)platform_vendor)[2] = cpuid[2]; // C
  std::string platform = std::string(platform_vendor, 12);

  __cpuid(1, cpuid);
  unsigned int platform_features = cpuid[3]; //D

  // __cpuid(1, cpuid);
  unsigned int num_logical_cpus = (cpuid[1] >> 16) & 0xFF; // B[23:16]
  // fprintf(stderr, "#logical cpus: ", num_logical_cpus);
  
  unsigned int num_cores = -1;
  if (platform == "GenuineIntel") {
	  __cpuid(4, cpuid);
	  num_cores = ((cpuid[0] >> 26) & 0x3f) + 1; //A[31:26] + 1
	  fprintf(stderr, "Platform vendor: Intel.\n");
  } else  {
	  fprintf(stderr, "Platform vendor unknown.\n");
  }

  // fprintf(stderr, "#physical cpus: ", num_cores);

  int ht = platform_features & (1 << 28) && num_cores < num_logical_cpus;
  if (ht)
	  fprintf(stderr, "CPUs support hyperThreading !!\n");

  return ht;
}

//---------------
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

void memoryAllocErt(ktp_aux_t *aux, worker_t &w, int ntid, char* idx_prefix) {
	mem_opt_t	*opt			  = aux->opt;	
	memSize = nreads;

	/* Mem allocation section for core kernels */
	w.regs = NULL; w.chain_ar = NULL; w.hits_ar = NULL; w.seedBuf = NULL;
	w.regs = (mem_alnreg_v *) calloc(memSize, sizeof(mem_alnreg_v));
	w.chain_ar = (mem_chain_v*) malloc (memSize * sizeof(mem_chain_v));
	w.seedBuf = (mem_seed_t *) calloc(memSize * AVG_SEEDS_PER_READ, sizeof(mem_seed_t));
	assert(w.seedBuf != NULL);
    w.seedBufSize = BATCH_SIZE * AVG_SEEDS_PER_READ;
	
	if (w.regs == NULL || w.chain_ar == NULL || w.seedBuf == NULL) {
		fprintf(stderr, "Memory not allocated!!\nExiting...\n");
		exit(0);
	}

	int64_t allocMem = memSize * sizeof(mem_alnreg_v) +
		memSize * sizeof(mem_chain_v) +
		sizeof(mem_seed_t) * memSize * AVG_SEEDS_PER_READ;
	fprintf(stderr, "------------------------------------------\n");
	fprintf(stderr, "Memory pre-allocation for chaining: %0.4lf MB\n", allocMem/1e6);

	
	/* SWA mem allocation */
	// int avg_seed_per_read = 35;
	w.size = BATCH_SIZE * SEEDS_PER_READ;
	for(int l=0; l<ntid; l++) {
        w.mmc.seqBufLeftRef[l * MAX_LINE_LEN] = (uint8_t *)malloc(w.size * MAX_SEQ_LEN_REF * sizeof(int8_t) + MAX_LINE_LEN);
        w.mmc.seqBufLeftQer[l * MAX_LINE_LEN] = (uint8_t *)malloc(w.size * MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_LINE_LEN);
        w.mmc.seqBufRightRef[l * MAX_LINE_LEN] = (uint8_t *)malloc(w.size * MAX_SEQ_LEN_REF * sizeof(int8_t) + MAX_LINE_LEN);
        w.mmc.seqBufRightQer[l * MAX_LINE_LEN] = (uint8_t *)malloc(w.size * MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_LINE_LEN);
        w.mmc.seqPairArrayAux[l * MAX_LINE_LEN]      = (SeqPair *) malloc(w.size * sizeof(SeqPair) + MAX_LINE_LEN);
		w.mmc.seqPairArrayLeft128[l * MAX_LINE_LEN]  = (SeqPair *) malloc(w.size * sizeof(SeqPair) + MAX_LINE_LEN);
		w.mmc.seqPairArrayRight128[l * MAX_LINE_LEN] = (SeqPair *) malloc(w.size * sizeof(SeqPair) + MAX_LINE_LEN);
		w.mmc.wsize[l * MAX_LINE_LEN] = w.size;
	}

	allocMem = (w.size * MAX_SEQ_LEN_REF * sizeof(int8_t) + MAX_LINE_LEN) * opt->n_threads * 2 +
		(w.size * MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_LINE_LEN) * opt->n_threads * 2 +		
		(w.size * sizeof(SeqPair) + MAX_LINE_LEN) * opt->n_threads * 3;
	
	fprintf(stderr, "Memory pre-allocation for BSW: %0.4lf MB\n", allocMem/1e6);
    
    w.mmc.lim = (int32_t *) _mm_malloc
		(nthreads * (BATCH_SIZE + 32) * sizeof(int32_t), 64);

	allocMem = nthreads * (BATCH_SIZE + 32) * sizeof(int32_t);
	fprintf(stderr, "Memory pre-allocation for BWT: %0.4lf MB\n", allocMem/1e6);
	fprintf(stderr, "------------------------------------------\n");
   
    char* kmer_tbl_file_name = (char*) malloc(strlen(idx_prefix) + 12);
    strcpy(kmer_tbl_file_name, idx_prefix); 
    strcat(kmer_tbl_file_name, ".kmer_table");
    char* ml_tbl_file_name = (char*) malloc(strlen(idx_prefix) + 12);
    strcpy(ml_tbl_file_name, idx_prefix);
    strcat(ml_tbl_file_name, ".mlt_table");

    FILE *kmer_tbl_fd, *ml_tbl_fd;

    kmer_tbl_fd = fopen(kmer_tbl_file_name, "rb");
    if (kmer_tbl_fd == NULL) {
        fprintf(stderr, "[M::%s::ERT] Can't open k-mer index\n.", __func__);
        exit(1);
    }
    ml_tbl_fd = fopen(ml_tbl_file_name, "rb");
    if (ml_tbl_fd == NULL) {
        fprintf(stderr, "[M::%s::ERT] Can't open multi-level tree index\n.", __func__);
        exit(1);
    }
        
	free(kmer_tbl_file_name);
	free(ml_tbl_file_name);

    double ctime, rtime;
    ctime = cputime(); rtime = realtime();
    allocMem = numKmers * 8L; 
    //
    // Read k-mer index
    //
    w.kmer_offsets = (uint64_t*) malloc(numKmers * sizeof(uint64_t));
    if (bwa_verbose >= 3) {
        fprintf(stderr, "[M::%s::ERT] Reading kmer index to memory\n", __func__);
    }
    fread(w.kmer_offsets, sizeof(uint64_t), numKmers, kmer_tbl_fd);
    // 
    // Read multi-level tree index
    //
    fseek(ml_tbl_fd, 0L, SEEK_END);
    long size = ftell(ml_tbl_fd);
    allocMem += size;
    w.mlt_table = (uint8_t*) malloc(size * sizeof(uint8_t));
    fseek(ml_tbl_fd, 0L, SEEK_SET);        
    if (bwa_verbose >= 3) {
        fprintf(stderr, "[M::%s::ERT] Reading multi-level tree index to memory\n", __func__);
    }
    fread(w.mlt_table, sizeof(uint8_t), size, ml_tbl_fd); 
    
    fclose(kmer_tbl_fd);
    fclose(ml_tbl_fd);

    if (bwa_verbose >= 3) {        
        fprintf(stderr, "[M::%s::ERT] Index tables loaded in %.3f CPU sec, %.3f real sec...\n", __func__, cputime() - ctime, realtime() - rtime);
    }

    allocMem += ((nthreads * MAX_LINE_LEN * sizeof(mem_v)) + (nthreads * MAX_LINE_LEN * sizeof(u64v)));
    allocMem += ((nthreads * BATCH_MUL * READ_LEN * sizeof(mem_t)) + (nthreads * MAX_HITS_PER_READ * sizeof(uint64_t)));
    w.smemBufSize = MAX_LINE_LEN * sizeof(mem_v);
    w.smems = (mem_v*) malloc(nthreads * w.smemBufSize);
    w.hitBufSize = MAX_LINE_LEN * sizeof(u64v);
    w.hits_ar = (u64v*) malloc(nthreads * w.hitBufSize);
    for (int i = 0 ; i < nthreads; ++i) {
        kv_init_base(mem_t, w.smems[i * MAX_LINE_LEN], BATCH_MUL * READ_LEN);
        kv_init_base(uint64_t, w.hits_ar[i * MAX_LINE_LEN], MAX_HITS_PER_READ);
    }
    w.useErt = 1;

    fprintf(stderr, "Memory pre-allocation for ERT: %0.4lf GB\n", allocMem/1e9);
	fprintf(stderr, "------------------------------------------\n");

}

void memoryAlloc(ktp_aux_t *aux, worker_t &w, int ntid)
{
	mem_opt_t	*opt			  = aux->opt;	
	memSize = nreads;

	/* Mem allocation section for core kernels */
	w.regs = NULL; w.chain_ar = NULL; w.seedBuf = NULL;
	w.regs = (mem_alnreg_v *) calloc(memSize, sizeof(mem_alnreg_v));
	w.chain_ar = (mem_chain_v*) malloc (memSize * sizeof(mem_chain_v));
    w.seedBuf = (mem_seed_t *) calloc(sizeof(mem_seed_t),  memSize * AVG_SEEDS_PER_READ);
	assert(w.seedBuf != NULL);
    w.seedBufSize = BATCH_SIZE * AVG_SEEDS_PER_READ;
	
	if (w.regs == NULL || w.chain_ar == NULL) {
		fprintf(stderr, "Memory not allocated!!\nExiting...\n");
		exit(EXIT_FAILURE);
	}

	int64_t allocMem = memSize * sizeof(mem_alnreg_v) +
		memSize * sizeof(mem_chain_v) +
		sizeof(mem_seed_t) * memSize * AVG_SEEDS_PER_READ;
	fprintf(stderr, "------------------------------------------\n");
	fprintf(stderr, "Memory pre-allocation for chaining: %0.4lf MB\n", allocMem/1e6);

	
	/* SWA mem allocation */
	// int avg_seed_per_read = 35;
	w.size = BATCH_SIZE * SEEDS_PER_READ;		
#if 0	
	w.mmc.seqPairArrayAux	   = (SeqPair *)_mm_malloc((w.size + MAX_LINE_LEN) * sizeof(SeqPair)
													   * opt->n_threads, 64);
	w.mmc.seqPairArrayLeft128  = (SeqPair *)_mm_malloc((w.size + MAX_LINE_LEN) * sizeof(SeqPair)
													   * opt->n_threads, 64);
	w.mmc.seqPairArrayRight128 = (SeqPair *)_mm_malloc((w.size + MAX_LINE_LEN) * sizeof(SeqPair)
													   * opt->n_threads, 64);
#else
	
	for(int l=0; l<ntid; l++) {
        w.mmc.seqBufLeftRef[l * MAX_LINE_LEN] = (uint8_t *)malloc(w.size * MAX_SEQ_LEN_REF * sizeof(int8_t) + MAX_LINE_LEN);
        w.mmc.seqBufLeftQer[l * MAX_LINE_LEN] = (uint8_t *)malloc(w.size * MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_LINE_LEN);
        w.mmc.seqBufRightRef[l * MAX_LINE_LEN] = (uint8_t *)malloc(w.size * MAX_SEQ_LEN_REF * sizeof(int8_t) + MAX_LINE_LEN);
        w.mmc.seqBufRightQer[l * MAX_LINE_LEN] = (uint8_t *)malloc(w.size * MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_LINE_LEN);
        w.mmc.seqPairArrayAux[l * MAX_LINE_LEN]      = (SeqPair *) malloc(w.size * sizeof(SeqPair) + MAX_LINE_LEN);
		w.mmc.seqPairArrayLeft128[l * MAX_LINE_LEN]  = (SeqPair *) malloc(w.size * sizeof(SeqPair) + MAX_LINE_LEN);
		w.mmc.seqPairArrayRight128[l * MAX_LINE_LEN] = (SeqPair *) malloc(w.size * sizeof(SeqPair) + MAX_LINE_LEN);
		w.mmc.wsize[l * MAX_LINE_LEN] = w.size;
	}
#endif
	
	allocMem = (w.size * MAX_SEQ_LEN_REF * sizeof(int8_t) + MAX_LINE_LEN) * opt->n_threads * 2 +
		(w.size * MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_LINE_LEN) * opt->n_threads * 2 +		
		(w.size * sizeof(SeqPair) + MAX_LINE_LEN) * opt->n_threads * 3;
	
	fprintf(stderr, "Memory pre-allocation for BSW: %0.4lf MB\n", allocMem/1e6);
	
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
	fprintf(stderr, "Memory pre-allocation for BWT: %0.4lf MB\n", allocMem/1e6);
	fprintf(stderr, "------------------------------------------\n");
    
    w.useErt = 0;
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
		ret->seqs = bseq_read_orig(aux->task_size,
								   &ret->n_seqs,
								   aux->ks, aux->ks2,
								   &sz);

		// *ret = *ret2;
		tprof[READ_IO][0] += __rdtsc() - tim;
		tprof[0][0] += sz;  // debug info, for accuracy checks!!
		tprof[0][2] += ret->n_seqs;
		// assert(ret->n_seqs <= nreads);
		
		fprintf(stderr, "[%.4d] read_chunk: %ld, work_chunk_size: %ld, nseq: %d\n",
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

			fprintf(stderr, "\t[%.4d][ M::%s] read %d sequences (%ld bp)...\n",
					myrank, __func__, ret->n_seqs, (long)size);
		}
				
		// printf("Threads leaving step 0, CPU: %d\n", sched_getcpu());
		return ret;
	} // Step 0			
	else if (step == 1)  /* Step 2: Main processing-engine */
	{
		static int task = 0;
		// printf("Thread entering step 1, CPU: %d\n", sched_getcpu());
		if (nreads < ret->n_seqs) {
			fprintf(stderr, "Reallocating initial memory allocations!!\n");
			free(w.regs); free(w.chain_ar); free(w.seedBuf);
			nreads = ret->n_seqs;
			w.regs = (mem_alnreg_v *) calloc(nreads, sizeof(mem_alnreg_v));
			w.chain_ar = (mem_chain_v*) malloc (nreads * sizeof(mem_chain_v));
			w.seedBuf = (mem_seed_t *) calloc(sizeof(mem_seed_t), nreads * AVG_SEEDS_PER_READ);
		}		
								
		fprintf(stderr, "[%.4d] 2. Calling mem_process_seqs.., task: %d\n", myrank, task++);

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

static int process(void *shared, gzFile gfp, gzFile gfp2, int pipe_threads, char* ert_idx_prefix)
{
	ktp_aux_t	*aux			  = (ktp_aux_t*) shared;
	worker_t	 w;
	mem_opt_t	*opt			  = aux->opt;

	nthreads = opt->n_threads; // global variable for profiling!
	int  deno = 1;
#if NUMA_ENABLED
	int tc = numa_num_task_cpus();
	int tn = numa_num_task_nodes();
	int tcc = numa_num_configured_cpus();
	fprintf(stderr, "num_cpus: %d, num_numas: %d, configured cpus: %d\n", tc, tn, tcc);
	int ht = HTStatus();
	if (ht) deno = 2;
	
	if (nthreads < tcc/tn/deno) {
		fprintf(stderr, "Enabling single numa domain...\n\n");
		// numa_set_preferred(0);
		// bitmask mask(0);
		struct bitmask *mask = numa_bitmask_alloc(numa_num_possible_nodes());
		numa_bitmask_clearall(mask);
		numa_bitmask_setbit(mask, 0);
		numa_bind(mask);
		numa_bitmask_free(mask);
	}
#endif
#if 0 && (__linux__)
	{ // Affinity/HT stuff
		unsigned int cpuid[4];
		asm volatile
			("cpuid" : "=a" (cpuid[0]), "=b" (cpuid[1]), "=c" (cpuid[2]), "=d" (cpuid[3])
			 : "0" (0xB), "2" (1));
		int num_logical_cpus = cpuid[1] & 0xFFFF;

		asm volatile
			("cpuid" : "=a" (cpuid[0]), "=b" (cpuid[1]), "=c" (cpuid[2]), "=d" (cpuid[3])
			 : "0" (0xB), "2" (0));
		int num_ht = cpuid[1] & 0xFFFF;
		int num_total_logical_cpus = get_nprocs_conf();
		int num_sockets = num_total_logical_cpus / num_logical_cpus;
		fprintf(stderr, "#sockets: %d, #cores/socket: %d, #logical_cpus: %d, #ht/core: %d\n",
				num_sockets, num_logical_cpus/num_ht, num_total_logical_cpus, num_ht);
		
		for (int i=0; i<num_total_logical_cpus; i++) affy[i] = i;
		int slookup[256] = {-1};

		if (num_ht == 2 && num_sockets == 2)  // generalize it for n sockets
		{
			for (int i=0; i<num_total_logical_cpus; i++) {
				std::ostringstream ss;
				ss << i;
				std::string str = "/sys/devices/system/cpu/cpu"+ ss.str();
				str = str +"/topology/thread_siblings_list";
				// std::cout << str << std::endl;
				// std::string str = "cpu.txt";
				FILE *fp = fopen(str.c_str(), "r");
				if (fp == NULL) {
					printf("Cant open the file..\n");
					break;
				}
				else {
					int a, b, v;
					char ch[10] = {'\0'};
					fgets(ch, 10, fp);
					// std::cout << ch;
					v = sscanf(ch, "%u,%u",&a,&b);
					if (v == 1) v = sscanf(ch, "%u-%u",&a,&b);
					// printf("2v: %d, %d %d\n", v, a, b);
					if (v == 1) {
						fprintf(stderr, "Mis-match between HT and threads_sibling_list...%s\n", ch);
						fprintf(stderr, "Continuing with default affinity settings..\n");
						break;
					}
					// affy[i] = a;
					// affy[i + lim] = b;
					slookup[a] = 1;
					slookup[b] = 2;
					fclose(fp);
				}
			}
			int a = 0, b = num_total_logical_cpus / num_ht;
			for (int i=0; i<num_total_logical_cpus; i++) {
				if (slookup[i] == -1) {
					fprintf(stderr, "Unseen cpu topology..\n");
					break;
				}
				if (slookup[i] == 1) affy[a++] = i;
				else affy[b++] = i;
			}
		}
	}
#endif
	
	fprintf(stderr, "\nThreads used (compute): %d\n",
			nthreads);
	nreads = aux->actual_chunk_size/ (readLen) + 10;
	fprintf(stderr, "Info: projected #read in a task: %ld\n", (long)nreads);
	
    
    /* All memory allocation */
    if (ert_idx_prefix) {
        memoryAllocErt(aux, w, nthreads, ert_idx_prefix);
    }
    else {
        memoryAlloc(aux, w, nthreads);
    }
	
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

	pthread_mutex_destroy(&aux_.mutex);
	pthread_cond_destroy(&aux_.cv);

	free(ptid);
	free(aux_.workers);
	/***** pipeline ends ******/
	
	fprintf(stderr, "[%.4d] Computation ends..\n", myrank);
	
	/* Dealloc memory allcoated in the header section */	
	free(w.chain_ar);
	free(w.regs);
    free(w.seedBuf);

#if 0
	_mm_free(w.mmc.seqPairArrayAux);
	_mm_free(w.mmc.seqPairArrayLeft128);
	_mm_free(w.mmc.seqPairArrayRight128);
#else
	for(int l=0; l<nthreads; l++) {
        free(w.mmc.seqBufLeftRef[l * MAX_LINE_LEN]);
        free(w.mmc.seqBufLeftQer[l * MAX_LINE_LEN]);
        free(w.mmc.seqBufRightRef[l * MAX_LINE_LEN]);
        free(w.mmc.seqBufRightQer[l * MAX_LINE_LEN]);
		free(w.mmc.seqPairArrayAux[l * MAX_LINE_LEN]);
		free(w.mmc.seqPairArrayLeft128[l * MAX_LINE_LEN]);
		free(w.mmc.seqPairArrayRight128[l * MAX_LINE_LEN]);
	}
#endif

    if (ert_idx_prefix) {
        free(w.kmer_offsets);
        free(w.mlt_table);
        for (int i = 0 ; i < nthreads; ++i) {
            kv_destroy(w.smems[i * MAX_LINE_LEN]);
            kv_destroy(w.hits_ar[i * MAX_LINE_LEN]);
        }
        free(w.smems);
        free(w.hits_ar);
    }
    else {    
        _mm_free(w.mmc.matchArray);
        free(w.mmc.min_intv_ar);
        free(w.mmc.query_pos_ar);
        free(w.mmc.enc_qdb);
        free(w.mmc.rid);
    }
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

static void usage(const mem_opt_t *opt)
{
	fprintf(stderr, "Usage: bwa2 mem [options] <idxbase> <in1.fq> [in2.fq]\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "  Algorithm options:\n");
	fprintf(stderr, "    -o STR        Output SAM file name\n");
	fprintf(stderr, "    -t INT        number of threads [%d]\n", opt->n_threads);
	fprintf(stderr, "    -k INT        minimum seed length [%d]\n", opt->min_seed_len);
	fprintf(stderr, "    -w INT        band width for banded alignment [%d]\n", opt->w);
	fprintf(stderr, "    -d INT        off-diagonal X-dropoff [%d]\n", opt->zdrop);
	fprintf(stderr, "    -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [%g]\n", opt->split_factor);
	fprintf(stderr, "    -y INT        seed occurrence for the 3rd round seeding [%ld]\n", (long)opt->max_mem_intv);
	fprintf(stderr, "    -c INT        skip seeds with more than INT occurrences [%d]\n", opt->max_occ);
	fprintf(stderr, "    -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [%.2f]\n", opt->drop_ratio);
	fprintf(stderr, "    -W INT        discard a chain if seeded bases shorter than INT [0]\n");
	fprintf(stderr, "    -m INT        perform at most INT rounds of mate rescues for each read [%d]\n", opt->max_matesw);
	fprintf(stderr, "    -S            skip mate rescue\n");
	fprintf(stderr, "    -o            output file name missing\n");
	fprintf(stderr, "    -P            skip pairing; mate rescue performed unless -S also in use\n");
	fprintf(stderr, "Scoring options:\n");
	fprintf(stderr, "   -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [%d]\n", opt->a);
	fprintf(stderr, "   -B INT        penalty for a mismatch [%d]\n", opt->b);
	fprintf(stderr, "   -O INT[,INT]  gap open penalties for deletions and insertions [%d,%d]\n", opt->o_del, opt->o_ins);
	fprintf(stderr, "   -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [%d,%d]\n", opt->e_del, opt->e_ins);
	fprintf(stderr, "   -L INT[,INT]  penalty for 5'- and 3'-end clipping [%d,%d]\n", opt->pen_clip5, opt->pen_clip3);
	fprintf(stderr, "   -U INT        penalty for an unpaired read pair [%d]\n", opt->pen_unpaired);
//	fprintf(stderr, "   -x STR        read type. Setting -x changes multiple parameters unless overriden [null]\n");
//	fprintf(stderr, "                 pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)\n");
//	fprintf(stderr, "                 ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)\n");
//	fprintf(stderr, "                 intractg: -B9 -O16 -L5  (intra-species contigs to ref)\n");
	fprintf(stderr, "Input/output options:\n");
	fprintf(stderr, "   -p            smart pairing (ignoring in2.fq)\n");
	fprintf(stderr, "   -R STR        read group header line such as '@RG\\tID:foo\\tSM:bar' [null]\n");
	fprintf(stderr, "   -H STR/FILE   insert STR to header if it starts with @; or insert lines in FILE [null]\n");
	fprintf(stderr, "   -j            treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)\n");
	fprintf(stderr, "   -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
	fprintf(stderr, "   -T INT        minimum score to output [%d]\n", opt->T);
	fprintf(stderr, "   -h INT[,INT]  if there are <INT hits with score >80%% of the max score, output all in XA [%d,%d]\n", opt->max_XA_hits, opt->max_XA_hits_alt);
	fprintf(stderr, "   -a            output all alignments for SE or unpaired PE\n");
	fprintf(stderr, "   -C            append FASTA/FASTQ comment to SAM output\n");
	fprintf(stderr, "   -V            output the reference FASTA header in the XR tag\n");
	fprintf(stderr, "   -Y            use soft clipping for supplementary alignments\n");
	fprintf(stderr, "   -M            mark shorter split hits as secondary\n");
	fprintf(stderr, "   -I FLOAT[,FLOAT[,INT[,INT]]]\n");
	fprintf(stderr, "                 specify the mean, standard deviation (10%% of the mean if absent), max\n");
	fprintf(stderr, "                 (4 sigma from the mean if absent) and min of the insert size distribution.\n");
	fprintf(stderr, "                 FR orientation only. [inferred]\n");
	fprintf(stderr, "   -Z            Use ERT index for seeding\n");
	fprintf(stderr, "Note: Please read the man page for detailed description of the command line and options.\n");
}

int main_mem(int argc, char *argv[])
{
	int			 i, c, ignore_alt = 0, no_mt_io = 0;
	int			 fixed_chunk_size		   = -1;
	char		*p, *rg_line			   = 0, *hdr_line = 0;
	const char	*mode					   = 0;
    int useErt = 0;
    char  *idx_prefix                = 0;

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
	while ((c = getopt(argc, argv, "1paMCSPVYjk:c:v:s:r:t:R:A:B:O:E:U:w:L:d:T:Q:D:m:I:N:W:x:G:h:y:K:X:H:o:q:Z")) >= 0)
	{
		if (c == 'k') opt->min_seed_len = atoi(optarg), opt0.min_seed_len = 1;
		else if (c == '1') no_mt_io = 1;
		else if (c == 'x') mode = optarg;
		else if (c == 'w') opt->w = atoi(optarg), opt0.w = 1;
		else if (c == 'A') opt->a = atoi(optarg), opt0.a = 1, assert(opt->a >= INT_MIN && opt->a <= INT_MAX);
		else if (c == 'B') opt->b = atoi(optarg), opt0.b = 1, assert(opt->b >= INT_MIN && opt->b <= INT_MAX);
		else if (c == 'T') opt->T = atoi(optarg), opt0.T = 1, assert(opt->T >= INT_MIN && opt->T <= INT_MAX);
		else if (c == 'U')
			opt->pen_unpaired = atoi(optarg), opt0.pen_unpaired = 1, assert(opt->pen_unpaired >= INT_MIN && opt->pen_unpaired <= INT_MAX);
		else if (c == 't')
			opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1, assert(opt->n_threads >= INT_MIN && opt->n_threads <= INT_MAX);
		else if (c == 'o')
		{
			is_o = 1;
			aux.fp = fopen(optarg, "w");
			if (aux.fp == NULL) {
				fprintf(stderr, "Error: can't open %s input file\n", optarg);
				exit(EXIT_FAILURE);
			}
			/*fclose(aux.fp);*/
			/*MPI_File_open(MPI_COMM_WORLD, optarg, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &aux.mfp);*/
			aux.totEl = 0;
		}
		else if (c == 'q') {
			nread_lim = atoi(optarg);
		}
		else if (c == 'P') opt->flag |= MEM_F_NOPAIRING;
		else if (c == 'a') opt->flag |= MEM_F_ALL;
		else if (c == 'p') opt->flag |= MEM_F_PE | MEM_F_SMARTPE;
		else if (c == 'M') opt->flag |= MEM_F_NO_MULTI;
		else if (c == 'S') opt->flag |= MEM_F_NO_RESCUE;
		else if (c == 'Y') opt->flag |= MEM_F_SOFTCLIP;
		else if (c == 'V') opt->flag |= MEM_F_REF_HDR;
		else if (c == 'c') opt->max_occ = atoi(optarg), opt0.max_occ = 1;
		else if (c == 'd') opt->zdrop = atoi(optarg), opt0.zdrop = 1;
		else if (c == 'v') bwa_verbose = atoi(optarg);
		else if (c == 'j') ignore_alt = 1;
		else if (c == 'r')
			opt->split_factor = atof(optarg), opt0.split_factor = 1.;
		else if (c == 'D') opt->drop_ratio = atof(optarg), opt0.drop_ratio = 1.;
		else if (c == 'm') opt->max_matesw = atoi(optarg), opt0.max_matesw = 1;
		else if (c == 's') opt->split_width = atoi(optarg), opt0.split_width = 1;
		else if (c == 'G')
			opt->max_chain_gap = atoi(optarg), opt0.max_chain_gap = 1;
		else if (c == 'N')
			opt->max_chain_extend = atoi(optarg), opt0.max_chain_extend = 1;
		else if (c == 'W')
			opt->min_chain_weight = atoi(optarg), opt0.min_chain_weight = 1;
		else if (c == 'y')
			opt->max_mem_intv = atol(optarg), opt0.max_mem_intv = 1;
		else if (c == 'C') aux.copy_comment = 1;
		else if (c == 'K') fixed_chunk_size = atoi(optarg);
		else if (c == 'X') opt->mask_level = atof(optarg);
		else if (c == 'h')
		{
			opt0.max_XA_hits = opt0.max_XA_hits_alt = 1;
			opt->max_XA_hits = opt->max_XA_hits_alt = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->max_XA_hits_alt = strtol(p+1, &p, 10);
		}
		else if (c == 'Q')
		{
			opt0.mapQ_coef_len = 1;
			opt->mapQ_coef_len = atoi(optarg);
			opt->mapQ_coef_fac = opt->mapQ_coef_len > 0? log(opt->mapQ_coef_len) : 0;
		}
		else if (c == 'O')
		{
			opt0.o_del = opt0.o_ins = 1;
			opt->o_del = opt->o_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->o_ins = strtol(p+1, &p, 10);
		}
		else if (c == 'E')
		{
			opt0.e_del = opt0.e_ins = 1;
			opt->e_del = opt->e_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->e_ins = strtol(p+1, &p, 10);
		}
		else if (c == 'L')
		{
			opt0.pen_clip5 = opt0.pen_clip3 = 1;
			opt->pen_clip5 = opt->pen_clip3 = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->pen_clip3 = strtol(p+1, &p, 10);
		}
		else if (c == 'R')
		{
		if ((rg_line = bwa_set_rg(optarg)) == 0) {
			free(opt);
			if (is_o)
				fclose(aux.fp);
				return 1;
            }
		}
		else if (c == 'H')
		{
			if (optarg[0] != '@')
			{
				FILE *fp;
				if ((fp = fopen(optarg, "r")) != 0)
				{
					char *buf;
					buf = (char *) calloc(1, 0x10000);
					while (fgets(buf, 0xffff, fp))
					{
						i = strlen(buf);
						assert(buf[i-1] == '\n');
						buf[i-1] = 0;
						hdr_line = bwa_insert_header(buf, hdr_line);
					}
					free(buf);
					fclose(fp);
				}
			} else hdr_line = bwa_insert_header(optarg, hdr_line);
		}
		else if (c == 'I')
		{
			aux.pes0 = pes;
			pes[1].failed = 0;
			pes[1].avg = strtod(optarg, &p);
			pes[1].std = pes[1].avg * .1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].std = strtod(p+1, &p);
			pes[1].high = (int)(pes[1].avg + 4. * pes[1].std + .499);
			pes[1].low  = (int)(pes[1].avg - 4. * pes[1].std + .499);
			if (pes[1].low < 1) pes[1].low = 1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].high = (int)(strtod(p+1, &p) + .499);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].low  = (int)(strtod(p+1, &p) + .499);
		}
		else if (c == 'Z') {
            useErt = 1;
        }
		else {
			free(opt);
			if (is_o)
				fclose(aux.fp);
			return 1;
		}
	}
	
	/* Check output file name */
	if (rg_line)
	{
		hdr_line = bwa_insert_header(rg_line, hdr_line);
		free(rg_line);
	}

	if (opt->n_threads < 1) opt->n_threads = 1;
	if (optind + 2 != argc && optind + 3 != argc) {
		usage(opt);
		free(opt);
		if (is_o) 
			fclose(aux.fp);
		return 1;
	}

	/* Further input parsing */
	if (mode)
	{
		fprintf(stderr, "WARNING: bwa-mem2 doesn't work well with long reads or contigs; please use minimap2 instead.\n");
		if (strcmp(mode, "intractg") == 0)
		{
			if (!opt0.o_del) opt->o_del = 16;
			if (!opt0.o_ins) opt->o_ins = 16;
			if (!opt0.b) opt->b = 9;
			if (!opt0.pen_clip5) opt->pen_clip5 = 5;
			if (!opt0.pen_clip3) opt->pen_clip3 = 5;
		}
		else if (strcmp(mode, "pacbio") == 0 || strcmp(mode, "pbref") == 0 || strcmp(mode, "ont2d") == 0)
		{
			if (!opt0.o_del) opt->o_del = 1;
			if (!opt0.e_del) opt->e_del = 1;
			if (!opt0.o_ins) opt->o_ins = 1;
			if (!opt0.e_ins) opt->e_ins = 1;
			if (!opt0.b) opt->b = 1;
			if (opt0.split_factor == 0.) opt->split_factor = 10.;
			if (strcmp(mode, "ont2d") == 0)
			{
				if (!opt0.min_chain_weight) opt->min_chain_weight = 20;
				if (!opt0.min_seed_len) opt->min_seed_len = 14;
				if (!opt0.pen_clip5) opt->pen_clip5 = 0;
				if (!opt0.pen_clip3) opt->pen_clip3 = 0;
			}
			else
			{
				if (!opt0.min_chain_weight) opt->min_chain_weight = 40;
				if (!opt0.min_seed_len) opt->min_seed_len = 17;
				if (!opt0.pen_clip5) opt->pen_clip5 = 0;
				if (!opt0.pen_clip3) opt->pen_clip3 = 0;
			}
		}
		else
		{
			fprintf(stderr, "[E::%s] unknown read type '%s'\n", __func__, mode);
			free(opt);
			if (is_o)
				fclose(aux.fp);
			return 1;
        }
	} else update_a(opt, &opt0);
	
	/* Matrix for SWA */
	bwa_fill_scmat(opt->a, opt->b, opt->mat);
	// tprof[PREPROCESS][0] += __rdtsc() - tim;

	/* Load bwt2/FMI index */
	{
		uint64_t tim = __rdtsc();

		fprintf(stderr, "Reference index file: %s\n", argv[optind]);			
		
        if (!useErt) {
            fmi = new FMI_search(argv[optind]);
        }
        else {
            idx_prefix = argv[optind];
            fmi = new FMI_search(argv[optind], BWA_IDX_BNS | BWA_IDX_PAC);
        }
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
			exit(EXIT_FAILURE);
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
		fprintf(stderr, "Done reading reference genome !!\n\n");

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

	readLen = READ_LEN;  // for memory pre-allocation // we realloc when required
	
	// Major function
	fprintf(stderr, "[%.4d] 1: Calling process()\n", myrank);

	uint64_t tim = __rdtsc();

    /* Relay process function */
    process(&aux, fp, fp2, no_mt_io? 1:2, idx_prefix);
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
	}

	// new bwt/FMI
    delete(fmi);
	
	return 0;
}

