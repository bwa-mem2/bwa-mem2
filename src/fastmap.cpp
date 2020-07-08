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

#if AFF && (__linux__)
#include <sys/sysinfo.h>
int affy[256];
#endif

// --------------
extern uint64_t tprof[LIM_R][LIM_C];
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


/*** Memory pre-allocations ***/
void memoryAlloc(ktp_aux_t *aux, worker_t &w, int32_t nreads, int32_t nthreads)
{
    mem_opt_t *opt = aux->opt;  
    int32_t memSize = nreads;
    int32_t readLen = READ_LEN;

    /* Mem allocation section for core kernels */
    w.regs = NULL; w.chain_ar = NULL; w.seedBuf = NULL;

    w.regs = (mem_alnreg_v *) calloc(memSize, sizeof(mem_alnreg_v));
    w.chain_ar = (mem_chain_v*) malloc (memSize * sizeof(mem_chain_v));
    w.seedBuf = (mem_seed_t *) calloc(sizeof(mem_seed_t),  memSize * AVG_SEEDS_PER_READ);

    assert(w.seedBuf  != NULL);
    assert(w.regs     != NULL);
    assert(w.chain_ar != NULL);

    w.seedBufSize = BATCH_SIZE * AVG_SEEDS_PER_READ;

    /*** printing ***/
    int64_t allocMem = memSize * sizeof(mem_alnreg_v) +
        memSize * sizeof(mem_chain_v) +
        sizeof(mem_seed_t) * memSize * AVG_SEEDS_PER_READ;
    fprintf(stderr, "------------------------------------------\n");
    fprintf(stderr, "1. Memory pre-allocation for Chaining: %0.4lf MB\n", allocMem/1e6);

    
    /* SWA mem allocation */
    int64_t wsize = BATCH_SIZE * SEEDS_PER_READ;
    for(int l=0; l<nthreads; l++)
    {
        w.mmc.seqBufLeftRef[l*CACHE_LINE]  = (uint8_t *)
            _mm_malloc(wsize * MAX_SEQ_LEN_REF * sizeof(int8_t) + MAX_LINE_LEN, 64);
        w.mmc.seqBufLeftQer[l*CACHE_LINE]  = (uint8_t *)
            _mm_malloc(wsize * MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_LINE_LEN, 64);
        w.mmc.seqBufRightRef[l*CACHE_LINE] = (uint8_t *)
            _mm_malloc(wsize * MAX_SEQ_LEN_REF * sizeof(int8_t) + MAX_LINE_LEN, 64);
        w.mmc.seqBufRightQer[l*CACHE_LINE] = (uint8_t *)
            _mm_malloc(wsize * MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_LINE_LEN, 64);
        
        w.mmc.wsize_buf_ref[l*CACHE_LINE] = wsize * MAX_SEQ_LEN_REF;
        w.mmc.wsize_buf_qer[l*CACHE_LINE] = wsize * MAX_SEQ_LEN_QER;
        
        assert(w.mmc.seqBufLeftRef[l*CACHE_LINE]  != NULL);
        assert(w.mmc.seqBufLeftQer[l*CACHE_LINE]  != NULL);
        assert(w.mmc.seqBufRightRef[l*CACHE_LINE] != NULL);
        assert(w.mmc.seqBufRightQer[l*CACHE_LINE] != NULL);
    }
    
    for(int l=0; l<nthreads; l++) {
        w.mmc.seqPairArrayAux[l]      = (SeqPair *) malloc((wsize + MAX_LINE_LEN)* sizeof(SeqPair));
        w.mmc.seqPairArrayLeft128[l]  = (SeqPair *) malloc((wsize + MAX_LINE_LEN)* sizeof(SeqPair));
        w.mmc.seqPairArrayRight128[l] = (SeqPair *) malloc((wsize + MAX_LINE_LEN)* sizeof(SeqPair));
        w.mmc.wsize[l] = wsize;

        assert(w.mmc.seqPairArrayAux[l] != NULL);
        assert(w.mmc.seqPairArrayLeft128[l] != NULL);
        assert(w.mmc.seqPairArrayRight128[l] != NULL);
    }   


    allocMem = (wsize * MAX_SEQ_LEN_REF * sizeof(int8_t) + MAX_LINE_LEN) * opt->n_threads * 2+
        (wsize * MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_LINE_LEN) * opt->n_threads  * 2 +       
        wsize * sizeof(SeqPair) * opt->n_threads * 3;   
    fprintf(stderr, "2. Memory pre-allocation for BSW: %0.4lf MB\n", allocMem/1e6);

    for (int l=0; l<nthreads; l++)
    {
        w.mmc.wsize_mem[l]     = BATCH_MUL * BATCH_SIZE *               readLen;
        w.mmc.matchArray[l]    = (SMEM *) _mm_malloc(w.mmc.wsize_mem[l] * sizeof(SMEM), 64);
        w.mmc.min_intv_ar[l]   = (int32_t *) malloc(w.mmc.wsize_mem[l] * sizeof(int32_t));
        w.mmc.query_pos_ar[l]  = (int16_t *) malloc(w.mmc.wsize_mem[l] * sizeof(int16_t));
        w.mmc.enc_qdb[l]       = (uint8_t *) malloc(w.mmc.wsize_mem[l] * sizeof(uint8_t));
        w.mmc.rid[l]           = (int32_t *) malloc(w.mmc.wsize_mem[l] * sizeof(int32_t));
        w.mmc.lim[l]           = (int32_t *) _mm_malloc((BATCH_SIZE + 32) * sizeof(int32_t), 64); // candidate not for reallocation, deferred for next round of changes.
    }

    allocMem = nthreads * BATCH_MUL * BATCH_SIZE * readLen * sizeof(SMEM) +
        nthreads * BATCH_MUL * BATCH_SIZE * readLen *sizeof(int32_t) +
        nthreads * BATCH_MUL * BATCH_SIZE * readLen *sizeof(int16_t) +
        nthreads * BATCH_MUL * BATCH_SIZE * readLen *sizeof(int32_t) +
        nthreads * (BATCH_SIZE + 32) * sizeof(int32_t);
    fprintf(stderr, "3. Memory pre-allocation for BWT: %0.4lf MB\n", allocMem/1e6);
    fprintf(stderr, "------------------------------------------\n");
}

ktp_data_t *kt_pipeline(void *shared, int step, void *data, mem_opt_t *opt, worker_t &w)
{
    ktp_aux_t *aux = (ktp_aux_t*) shared;
    ktp_data_t *ret = (ktp_data_t*) data;

    if (step == 0)
    {
        ktp_data_t *ret = (ktp_data_t *) calloc(1, sizeof(ktp_data_t));
        assert(ret != NULL);
        uint64_t tim = __rdtsc();

        /* Read "reads" from input file (fread) */
        int64_t sz = 0;
        ret->seqs = bseq_read_orig(aux->task_size,
                                   &ret->n_seqs,
                                   aux->ks, aux->ks2,
                                   &sz);

        tprof[READ_IO][0] += __rdtsc() - tim;
        
        fprintf(stderr, "[0000] read_chunk: %ld, work_chunk_size: %ld, nseq: %d\n",
                aux->task_size, sz, ret->n_seqs);   

        if (ret->seqs == 0) {
            free(ret);
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

            fprintf(stderr, "\t[0000][ M::%s] read %d sequences (%ld bp)...\n",
                    __func__, ret->n_seqs, (long)size);
        }
                
        return ret;
    } // Step 0         
    else if (step == 1)  /* Step 2: Main processing-engine */
    {
        static int task = 0;
        if (w.nreads < ret->n_seqs)
        {
            fprintf(stderr, "[0000] Reallocating initial memory allocations!!\n");
            free(w.regs); free(w.chain_ar); free(w.seedBuf);
            w.nreads = ret->n_seqs;
            w.regs = (mem_alnreg_v *) calloc(w.nreads, sizeof(mem_alnreg_v));
            w.chain_ar = (mem_chain_v*) malloc (w.nreads * sizeof(mem_chain_v));
            w.seedBuf = (mem_seed_t *) calloc(sizeof(mem_seed_t), w.nreads * AVG_SEEDS_PER_READ);
            assert(w.regs != NULL); assert(w.chain_ar != NULL); assert(w.seedBuf != NULL);
        }       
                                
        fprintf(stderr, "[0000] Calling mem_process_seqs.., task: %d\n", task++);

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
                                 aux->n_processed,
                                 n_sep[0],
                                 sep[0],
                                 0,
                                 w);
                
                for (int i = 0; i < n_sep[0]; ++i)
                    ret->seqs[sep[0][i].id].sam = sep[0][i].sam;
            }
            if (n_sep[1]) {
                tmp_opt.flag |= MEM_F_PE;
                /* paired-end sequences, in the mixture */
                mem_process_seqs(&tmp_opt,
                                 aux->n_processed + n_sep[0],
                                 n_sep[1],
                                 sep[1],
                                 aux->pes0,
                                 w);
                                
                for (int i = 0; i < n_sep[1]; ++i)
                    ret->seqs[sep[1][i].id].sam = sep[1][i].sam;
            }
            free(sep[0]); free(sep[1]);
        }
        else {
            /* pure (single/paired-end), reads processing */
            mem_process_seqs(opt,
                             aux->n_processed,
                             ret->n_seqs,
                             ret->seqs,
                             aux->pes0,
                             w);
        }               
        tprof[MEM_PROCESS2][0] += __rdtsc() - tim;
                
        return ret;
    }           
    /* Step 3: Write output */
    else if (step == 2)
    {
        aux->n_processed += ret->n_seqs;
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
        free(ret);
        tprof[SAM_IO][0] += __rdtsc() - tim;

        return 0;
    } // step 2
    
    return 0;
}

static void *ktp_worker(void *data)
{   
    ktp_worker_t *w = (ktp_worker_t*) data;
    ktp_t *p = w->pl;
    
    while (w->step < p->n_steps) {
        // test whether we can kick off the job with this worker
        int pthread_ret = pthread_mutex_lock(&p->mutex);
        assert(pthread_ret == 0);
        for (;;) {
            int i;
            // test whether another worker is doing the same step
            for (i = 0; i < p->n_workers; ++i) {
                if (w == &p->workers[i]) continue; // ignore itself
                if (p->workers[i].step <= w->step && p->workers[i].index < w->index)
                    break;
            }
            if (i == p->n_workers) break; // no workers with smaller indices are doing w->step or the previous steps
            pthread_ret = pthread_cond_wait(&p->cv, &p->mutex);
            assert(pthread_ret == 0);
        }
        pthread_ret = pthread_mutex_unlock(&p->mutex);
        assert(pthread_ret == 0);

        // working on w->step
        w->data = kt_pipeline(p->shared, w->step, w->step? w->data : 0, w->opt, *(w->w)); // for the first step, input is NULL

        // update step and let other workers know
        pthread_ret = pthread_mutex_lock(&p->mutex);
        assert(pthread_ret == 0);
        w->step = w->step == p->n_steps - 1 || w->data? (w->step + 1) % p->n_steps : p->n_steps;

        if (w->step == 0) w->index = p->index++;
        pthread_ret = pthread_cond_broadcast(&p->cv);
        assert(pthread_ret == 0);
        pthread_ret = pthread_mutex_unlock(&p->mutex);
        assert(pthread_ret == 0);
    }
    pthread_exit(0);
}

static int process(void *shared, gzFile gfp, gzFile gfp2, int pipe_threads)
{
    ktp_aux_t   *aux = (ktp_aux_t*) shared;
    worker_t     w;
    mem_opt_t   *opt = aux->opt;
    int32_t nthreads = opt->n_threads; // global variable for profiling!
    w.nthreads = opt->n_threads;
    
#if NUMA_ENABLED
    int  deno = 1;
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
#if AFF && (__linux__)
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
                    fprintf("Error: Cant open the file..\n");
                    break;
                }
                else {
                    int a, b, v;
                    char ch[10] = {'\0'};
                    fgets(ch, 10, fp);
                    v = sscanf(ch, "%u,%u",&a,&b);
                    if (v == 1) v = sscanf(ch, "%u-%u",&a,&b);
                    if (v == 1) {
                        fprintf(stderr, "Mis-match between HT and threads_sibling_list...%s\n", ch);
                        fprintf(stderr, "Continuing with default affinity settings..\n");
                        break;
                    }
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
    
    int32_t nreads = aux->actual_chunk_size/ READ_LEN + 10;
    
    /* All memory allocation */
    memoryAlloc(aux, w, nreads, nthreads);
    fprintf(stderr, "* Threads used (compute): %d\n", nthreads);
    
    /* pipeline using pthreads */
    ktp_t aux_;
    int p_nt = pipe_threads; // 2;
    int n_steps = 3;
    
    w.ref_string = aux->ref_string;
    w.fmi = aux->fmi;
    w.nreads  = nreads;
    // w.memSize = nreads;
    
    aux_.n_workers = p_nt;
    aux_.n_steps = n_steps;
    aux_.shared = aux;
    aux_.index = 0;
    int pthread_ret = pthread_mutex_init(&aux_.mutex, 0);
    assert(pthread_ret == 0);
    pthread_ret = pthread_cond_init(&aux_.cv, 0);
    assert(pthread_ret == 0);

    fprintf(stderr, "* No. of pipeline threads: %d\n\n", p_nt);
    aux_.workers = (ktp_worker_t*) malloc(p_nt * sizeof(ktp_worker_t));
    assert(aux_.workers != NULL);
    
    for (int i = 0; i < p_nt; ++i) {
        ktp_worker_t *wr = &aux_.workers[i];
        wr->step = 0; wr->pl = &aux_; wr->data = 0;
        wr->index = aux_.index++;
        wr->i = i;
        wr->opt = opt;
        wr->w = &w;
    }
    
    pthread_t *ptid = (pthread_t *) calloc(p_nt, sizeof(pthread_t));
    assert(ptid != NULL);
    
    for (int i = 0; i < p_nt; ++i)
        pthread_create(&ptid[i], 0, ktp_worker, (void*) &aux_.workers[i]);
    
    for (int i = 0; i < p_nt; ++i)
        pthread_join(ptid[i], 0);

    pthread_ret = pthread_mutex_destroy(&aux_.mutex);
    assert(pthread_ret == 0);
    pthread_ret = pthread_cond_destroy(&aux_.cv);
    assert(pthread_ret == 0);

    free(ptid);
    free(aux_.workers);
    /***** pipeline ends ******/
    
    fprintf(stderr, "[0000] Computation ends..\n");
    
    /* Dealloc memory allcoated in the header section */    
    free(w.chain_ar);
    free(w.regs);
    free(w.seedBuf);
    
    for(int l=0; l<nthreads; l++) {
        _mm_free(w.mmc.seqBufLeftRef[l*CACHE_LINE]);
        _mm_free(w.mmc.seqBufRightRef[l*CACHE_LINE]);
        _mm_free(w.mmc.seqBufLeftQer[l*CACHE_LINE]);
        _mm_free(w.mmc.seqBufRightQer[l*CACHE_LINE]);
    }

    for(int l=0; l<nthreads; l++) {
        free(w.mmc.seqPairArrayAux[l]);
        free(w.mmc.seqPairArrayLeft128[l]);
        free(w.mmc.seqPairArrayRight128[l]);
    }

    for(int l=0; l<nthreads; l++) {
        _mm_free(w.mmc.matchArray[l]);
        free(w.mmc.min_intv_ar[l]);
        free(w.mmc.query_pos_ar[l]);
        free(w.mmc.enc_qdb[l]);
        free(w.mmc.rid[l]);
        _mm_free(w.mmc.lim[l]);
    }

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
    fprintf(stderr, "Usage: bwa-mem2 mem [options] <idxbase> <in1.fq> [in2.fq]\n");
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
//  fprintf(stderr, "   -x STR        read type. Setting -x changes multiple parameters unless overriden [null]\n");
//  fprintf(stderr, "                 pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)\n");
//  fprintf(stderr, "                 ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)\n");
//  fprintf(stderr, "                 intractg: -B9 -O16 -L5  (intra-species contigs to ref)\n");
    fprintf(stderr, "Input/output options:\n");
    fprintf(stderr, "   -p            smart pairing (ignoring in2.fq)\n");
    fprintf(stderr, "   -R STR        read group header line such as '@RG\\tID:foo\\tSM:bar' [null]\n");
    fprintf(stderr, "   -H STR/FILE   insert STR to header if it starts with @; or insert lines in FILE [null]\n");
    fprintf(stderr, "   -j            treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)\n");
    fprintf(stderr, "   -5            for split alignment, take the alignment with the smallest coordinate as primary\n");
    fprintf(stderr, "   -q            don't modify mapQ of supplementary alignments\n");
    fprintf(stderr, "   -K INT        process INT input bases in each batch regardless of nThreads (for reproducibility) []\n");    
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
    fprintf(stderr, "Note: Please read the man page for detailed description of the command line and options.\n");
}

int main_mem(int argc, char *argv[])
{
    int          i, c, ignore_alt = 0, no_mt_io = 0;
    int          fixed_chunk_size          = -1;
    char        *p, *rg_line               = 0, *hdr_line = 0;
    const char  *mode                      = 0;
    
    mem_opt_t    *opt, opt0;
    gzFile        fp, fp2 = 0;
    void         *ko = 0, *ko2 = 0;
    int           fd, fd2;
    mem_pestat_t  pes[4];
    ktp_aux_t     aux;
    bool          is_o    = 0;
    uint8_t      *ref_string;
    
    memset(&aux, 0, sizeof(ktp_aux_t));
    memset(pes, 0, 4 * sizeof(mem_pestat_t));
    for (i = 0; i < 4; ++i) pes[i].failed = 1;
    
    // opterr = 0;
    aux.fp = stdout;
    aux.opt = opt = mem_opt_init();
    memset(&opt0, 0, sizeof(mem_opt_t));
    
    /* Parse input arguments */
    // comment: added option '5' in the list
    while ((c = getopt(argc, argv, "51qpaMCSPVYjk:c:v:s:r:t:R:A:B:O:E:U:w:L:d:T:Q:D:m:I:N:W:x:G:h:y:K:X:H:o:f:")) >= 0)
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
        else if (c == 'o' || c == 'f')
        {
            is_o = 1;
            aux.fp = fopen(optarg, "w");
            if (aux.fp == NULL) {
                fprintf(stderr, "Error: can't open %s input file\n", optarg);
                exit(EXIT_FAILURE);
            }
            /*fclose(aux.fp);*/
        }
        else if (c == 'P') opt->flag |= MEM_F_NOPAIRING;
        else if (c == 'a') opt->flag |= MEM_F_ALL;
        else if (c == 'p') opt->flag |= MEM_F_PE | MEM_F_SMARTPE;
        else if (c == 'M') opt->flag |= MEM_F_NO_MULTI;
        else if (c == 'S') opt->flag |= MEM_F_NO_RESCUE;
        else if (c == 'Y') opt->flag |= MEM_F_SOFTCLIP;
        else if (c == 'V') opt->flag |= MEM_F_REF_HDR;
        else if (c == '5') opt->flag |= MEM_F_PRIMARY5 | MEM_F_KEEP_SUPP_MAPQ; // always apply MEM_F_KEEP_SUPP_MAPQ with -5
        else if (c == 'q') opt->flag |= MEM_F_KEEP_SUPP_MAPQ;
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
                    assert(buf != NULL);
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
    
    /* Load bwt2/FMI index */
    uint64_t tim = __rdtsc();
    
    fprintf(stderr, "* Ref file: %s\n", argv[optind]);          
    aux.fmi = new FMI_search(argv[optind]);
    aux.fmi->load_index();
    tprof[FMI][0] += __rdtsc() - tim;
    
    // reading ref string from the file
    tim = __rdtsc();
    fprintf(stderr, "* Reading reference genome..\n");
    
    char binary_seq_file[PATH_MAX];
    strcpy_s(binary_seq_file, PATH_MAX, argv[optind]);
    strcat_s(binary_seq_file, PATH_MAX, ".0123");
    //sprintf(binary_seq_file, "%s.0123", argv[optind]);
    
    fprintf(stderr, "* Binary seq file = %s\n", binary_seq_file);
    FILE *fr = fopen(binary_seq_file, "r");
    
    if (fr == NULL) {
        fprintf(stderr, "Error: can't open %s input file\n", binary_seq_file);
        exit(EXIT_FAILURE);
    }
    
    int64_t rlen = 0;
    fseek(fr, 0, SEEK_END); 
    rlen = ftell(fr);
    ref_string = (uint8_t*) _mm_malloc(rlen, 64);
    aux.ref_string = ref_string;
    rewind(fr);
    
    /* Reading ref. sequence */
    err_fread_noeof(ref_string, 1, rlen, fr);
    
    uint64_t timer  = __rdtsc();
    tprof[REF_IO][0] += timer - tim;
    
    fclose(fr);
    fprintf(stderr, "* Reference genome size: %ld bp\n", rlen);
    fprintf(stderr, "* Done reading reference genome !!\n\n");
    
    if (ignore_alt)
        for (i = 0; i < aux.fmi->idx->bns->n_seqs; ++i)
            aux.fmi->idx->bns->anns[i].is_alt = 0;

    /* READS file operations */
    ko = kopen(argv[optind + 1], &fd);
	if (ko == 0) {
		fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 1]);
        free(opt);
        if (is_o) 
            fclose(aux.fp);
        delete aux.fmi;
        kclose(ko);
        return 1;
    }
    // fp = gzopen(argv[optind + 1], "r");
    fp = gzdopen(fd, "r");
    aux.ks = kseq_init(fp);
    
    // PAIRED_END
    /* Handling Paired-end reads */
    aux.ks2 = 0;
    if (optind + 2 < argc) {
        if (opt->flag & MEM_F_PE) {
            fprintf(stderr, "[W::%s] when '-p' is in use, the second query file is ignored.\n",
                    __func__);
        }
        else
        {
            ko2 = kopen(argv[optind + 2], &fd2);
            if (ko2 == 0) {
                fprintf(stderr, "[E::%s] failed to open file `%s'.\n", __func__, argv[optind + 2]);
                free(opt);
                free(ko);
                err_gzclose(fp);
                kseq_destroy(aux.ks);
                if (is_o) 
                    fclose(aux.fp);             
                delete aux.fmi;
                kclose(ko);
                kclose(ko2);
                return 1;
            }            
            // fp2 = gzopen(argv[optind + 2], "r");
            fp2 = gzdopen(fd2, "r");
            aux.ks2 = kseq_init(fp2);
            opt->flag |= MEM_F_PE;
            assert(aux.ks2 != 0);
        }
    }

    bwa_print_sam_hdr(aux.fmi->idx->bns, hdr_line, aux.fp);

    if (fixed_chunk_size > 0)
        aux.task_size = fixed_chunk_size;
    else {
        //aux.task_size = 10000000 * opt->n_threads; //aux.actual_chunk_size;
        aux.task_size = opt->chunk_size * opt->n_threads; //aux.actual_chunk_size;
    }
    tprof[MISC][1] = opt->chunk_size = aux.actual_chunk_size = aux.task_size;

    tim = __rdtsc();

    /* Relay process function */
    process(&aux, fp, fp2, no_mt_io? 1:2);
    
    tprof[PROCESS][0] += __rdtsc() - tim;

    // free memory
    _mm_free(ref_string);
    free(hdr_line);
    free(opt);
    kseq_destroy(aux.ks);   
    err_gzclose(fp); kclose(ko);

    // PAIRED_END
    if (aux.ks2) {
        kseq_destroy(aux.ks2);
        err_gzclose(fp2); kclose(ko2);
    }
    
    if (is_o) {
        fclose(aux.fp);
    }

    // new bwt/FMI
    delete(aux.fmi);    

    /* Display runtime profiling stats */
    tprof[MEM][0] = __rdtsc() - tprof[MEM][0];
    display_stats(aux.opt->n_threads);
    
    return 0;
}

