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

#ifndef BWAMEM_HPP
#define BWAMEM_HPP
#include "sais.h"
#include "bwt.h"
#include "bntseq.h"
#include "bwa.h"
#include "kthread.h"
#include "macro.h"
#include "bandedSWA.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <vector>
#include "kstring.h"
#include "ksw.h"
#include "kvec.h"
#include "ksort.h"
#include "utils.h"
#include "macro.h"
#include "profiling.h"
#include "FMI_search.h"
#include "LISA_search.h"
#define MEM_MAPQ_COEF 30.0
#define MEM_MAPQ_MAX  60

struct __smem_i;
typedef struct __smem_i smem_i;

#define MEM_F_PE        0x2
#define MEM_F_NOPAIRING 0x4
#define MEM_F_ALL       0x8
#define MEM_F_NO_MULTI  0x10
#define MEM_F_NO_RESCUE 0x20
#define MEM_F_REF_HDR   0x100
#define MEM_F_SOFTCLIP  0x200
#define MEM_F_SMARTPE   0x400

// V17
#define MEM_F_PRIMARY5  0x800
#define MEM_F_KEEP_SUPP_MAPQ 0x1000


typedef struct mem_opt_t {
    int a, b;               // match score and mismatch penalty
    int o_del, e_del;
    int o_ins, e_ins;
    int pen_unpaired;       // phred-scaled penalty for unpaired reads
    int pen_clip5,pen_clip3;// clipping penalty. This score is not deducted from the DP score.
    int w;                  // band width
    int zdrop;              // Z-dropoff

    uint64_t max_mem_intv;

    int T;                  // output score threshold; only affecting output
    int flag;               // see MEM_F_* macros
    int min_seed_len;       // minimum seed length
    int min_chain_weight;
    int max_chain_extend;
    float split_factor;     // split into a seed if MEM is longer than min_seed_len*split_factor
    int split_width;        // split into a seed if its occurence is smaller than this value
    int max_occ;            // skip a seed if its occurence is larger than this value
    int max_chain_gap;      // do not chain seed if it is max_chain_gap-bp away from the closest seed
    int n_threads;          // number of threads
    int64_t chunk_size;         // process chunk_size-bp sequences in a batch
    float mask_level;       // regard a hit as redundant if the overlap with another better hit is over mask_level times the min length of the two hits
    float drop_ratio;       // drop a chain if its seed coverage is below drop_ratio times the seed coverage of a better chain overlapping with the small chain
    float XA_drop_ratio;    // when counting hits for the XA tag, ignore alignments with score < XA_drop_ratio * max_score; only effective for the XA tag
    float mask_level_redun;
    float mapQ_coef_len;
    int mapQ_coef_fac;
    int max_ins;            // when estimating insert size distribution, skip pairs with insert longer than this value
    int max_matesw;         // perform maximally max_matesw rounds of mate-SW for each end
    int max_XA_hits, max_XA_hits_alt; // if there are max_hits or fewer, output them all
    int8_t mat[25];         // scoring matrix; mat[0] == 0 if unset
} mem_opt_t;


struct mem_alnreg_t;
// * Chaining *
typedef struct abc {
    abc() {
        done = 0;
        rbeg = qbeg = len = score = aln = 0;
    }
    int64_t rbeg;
    int32_t qbeg;
    int32_t len;
    int32_t score;
    int8_t done;
    int aln;
} mem_seed_t; // unaligned memory

typedef struct {
    int32_t seqid, cseed;
    int32_t n, m, first, rid;
    uint32_t w:29, kept:2, is_alt:1;
    float frac_rep;
    int64_t pos;
    mem_seed_t *seeds;
} mem_chain_t;

typedef struct { size_t n, m, cc; mem_chain_t *a;  } mem_chain_v;

typedef struct mem_alnreg_t {
    // mem_alnreg_t() {c=NULL;}
    int64_t rb, re; // [rb,re): reference sequence in the alignment
    int qb, qe;     // [qb,qe): query sequence in the alignment
    int rid;        // reference seq ID
    mem_chain_t *c;
    int score;      // best local SW score
    int truesc;     // actual score corresponding to the aligned region; possibly smaller than $score
    int sub;        // 2nd best SW score
    int alt_sc;
    int csub;       // SW score of a tandem hit
    int sub_n;      // approximate number of suboptimal hits
    int w;          // actual band width used in extension
    int seedcov;    // length of regions coverged by seeds
    int secondary;  // index of the parent hit shadowing the current hit; <0 if primary
    int secondary_all;
    int seedlen0;   // length of the starting seed
    int n_comp:30, is_alt:2; // number of sub-alignments chained together
    float frac_rep;
    uint64_t hash;
    int flg;
} mem_alnreg_t;

typedef struct { size_t n, m; mem_alnreg_t *a; } mem_alnreg_v;

typedef struct {
    int low, high;   // lower and upper bounds within which a read pair is considered to be properly paired
    int failed;      // non-zero if the orientation is not supported by sufficient data
    double avg, std; // mean and stddev of the insert size distribution
} mem_pestat_t;

typedef struct { // This struct is only used for the convenience of API.
    int64_t pos;     // forward strand 5'-end mapping position
    int rid;         // reference sequence index in bntseq_t; <0 for unmapped
    int flag;        // extra flag
    uint32_t is_rev:1, is_alt:1, mapq:8, NM:22; // is_rev: whether on the reverse strand; mapq: mapping quality; NM: edit distance
    int n_cigar;     // number of CIGAR operations
    uint32_t *cigar; // CIGAR in the BAM encoding: opLen<<4|op; op to integer mapping: MIDSH=>01234
    char *XA;        // alternative mappings

    int score, sub, alt_sc;
} mem_aln_t;

// struct
typedef struct {
    bwtintv_v mem, mem1, *tmpv[2];
} smem_aux_t;

typedef struct
{
    SeqPair *seqPairArrayAux[MAX_THREADS];
    SeqPair *seqPairArrayLeft128[MAX_THREADS];
    SeqPair *seqPairArrayRight128[MAX_THREADS];
    
    int64_t wsize[MAX_THREADS];

    int64_t wsize_buf_ref[MAX_THREADS*CACHE_LINE]; 
    int64_t wsize_buf_qer[MAX_THREADS*CACHE_LINE];

    uint8_t *seqBufLeftRef[MAX_THREADS*CACHE_LINE];
    uint8_t *seqBufRightRef[MAX_THREADS*CACHE_LINE];
    uint8_t *seqBufLeftQer[MAX_THREADS*CACHE_LINE];
    uint8_t *seqBufRightQer[MAX_THREADS*CACHE_LINE];    

    SMEM *matchArray[MAX_THREADS];
    int32_t *min_intv_ar[MAX_THREADS];
    int32_t *rid[MAX_THREADS];
    int32_t *lim[MAX_THREADS];
    int16_t *query_pos_ar[MAX_THREADS];
    uint8_t *enc_qdb[MAX_THREADS];
    
    int64_t wsize_mem[MAX_THREADS];

    threadData *td[MAX_THREADS];
} mem_cache;

// chain moved to 
typedef struct worker_t {
    const mem_opt_t      *opt;
    //const bntseq_t         *bns;
    // const uint8_t         *pac;
    const mem_pestat_t   *pes;
    smem_aux_t      **aux;
    bseq1_t          *seqs;
    mem_alnreg_v     *regs;
    int64_t           n_processed;
    mem_chain_v      *chain_ar;
    mem_cache         mmc;
    mem_seed_t       *seedBuf;
    int64_t           seedBufSize;
    mem_seed_t       *auxSeedBuf;
    int64_t           auxSeedBufSize;
    uint8_t          *ref_string;
    int16_t           nthreads;
    int32_t           nreads;
    FMI_search       *fmi; 
    //QBWT_HYBRID<index_t> *lisa; 
    LISA_search<index_t> *lisa; 
} worker_t;


typedef kvec_t(int) int_v;

smem_i *smem_itr_init(const bwt_t *bwt);
void smem_itr_destroy(smem_i *itr);
void smem_set_query(smem_i *itr, int len, const uint8_t *query);
void smem_config(smem_i *itr, int min_intv, int max_len, uint64_t max_intv);
const bwtintv_v *smem_next(smem_i *itr);

mem_opt_t *mem_opt_init(void);
void mem_fill_scmat(int a, int b, int8_t mat[25]);

void mem_reg2sam(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
                 bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m);

int mem_approx_mapq_se(const mem_opt_t *opt, const mem_alnreg_t *a) ;

int mem_mark_primary_se(const mem_opt_t *opt, int n, mem_alnreg_t *a, int64_t id);

static void mem_mark_primary_se_core(const mem_opt_t *opt, int n, mem_alnreg_t *a, int_v *z);

char **mem_gen_alt(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
                   const mem_alnreg_v *a, int l_query, const char *query); // ONLY work after mem_mark_primary_se()
void mem_aln2sam(const mem_opt_t *opt, const bntseq_t *bns, kstring_t *str, bseq1_t *s,
                 int n, const mem_aln_t *list, int which, const mem_aln_t *m_);

static inline int get_rlen(int n_cigar, const uint32_t *cigar);
static inline int infer_bw(int l1, int l2, int score, int a, int q, int r);

int mem_kernel1_core(FMI_search *fmi, 
                     //QBWT_HYBRID<index_t> *lisa, 
                     LISA_search<index_t> *lisa, 
		     const mem_opt_t *opt,
                     bseq1_t *seq_,
                     int nseq,
                     mem_chain_v *chain_ar,
                     mem_cache *mmc,
                     int tid);

void* _mm_realloc(void *ptr, int64_t csize, int64_t nsize, int16_t dsize);

void mem_chain2aln_across_reads_V2(const mem_opt_t *opt, const bntseq_t *bns,
                                   const uint8_t *pac, bseq1_t *seq_, int nseq,
                                   mem_chain_v* chain_ar, mem_alnreg_v *av_v,
                                   mem_cache *mmc, uint8_t *ref_string, int tid);

int mem_sam_pe_batch_pre(const mem_opt_t *opt, const bntseq_t *bns,
                         const uint8_t *pac, const mem_pestat_t pes[4],
                         uint64_t id, bseq1_t s[2], mem_alnreg_v a[2],
                         mem_cache *mmc,  int64_t &pcnt, int32_t &gcnt,
                         int32_t&, int32_t&, int tid);

int mem_matesw_batch_pre(const mem_opt_t *opt, const bntseq_t *bns,
                         const uint8_t *pac, const mem_pestat_t pes[4],
                         const mem_alnreg_t *a, int l_ms, const uint8_t *ms,
                         mem_alnreg_v *ma, mem_cache *mmc, int pcnt, int32_t gcnt,
                         int32_t &maxRefLen, int32_t &maxQerLen, int32_t tid);

int mem_sam_pe_batch(const mem_opt_t *opt, mem_cache *mmc,
                     int64_t &pcnt, int64_t &pcnt8, kswr_t *aln,
                     int32_t, int32_t, int tid);

int mem_sam_pe_batch_post(const mem_opt_t *opt, const bntseq_t *bns,
                          const uint8_t *pac, const mem_pestat_t pes[4],
                          uint64_t id, bseq1_t s[2], mem_alnreg_v a[2],
                          kswr_t **myaln, mem_cache *mmc,
                          int32_t &gcnt, int tid);

int mem_matesw_batch_post(const mem_opt_t *opt, const bntseq_t *bns,
                          const uint8_t *pac, const mem_pestat_t pes[4],
                          const mem_alnreg_t *a, int l_ms, const uint8_t *ms,
                          mem_alnreg_v *ma, kswr_t **myaln, int32_t gcnt,
                          int32_t *gar, mem_cache *mmc);

int mem_sam_pe(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
               const mem_pestat_t pes[4], uint64_t id, bseq1_t s[2],
               mem_alnreg_v a[2]);
/**
 * Align a batch of sequences and generate the alignments in the SAM format
 *
 * This routine requires $seqs[i].{l_seq,seq,name} and write $seqs[i].sam.
 * Note that $seqs[i].sam may consist of several SAM lines if the
 * corresponding sequence has multiple primary hits.
 *
 * In the paired-end mode (i.e. MEM_F_PE is set in $opt->flag), query
 * sequences must be interleaved: $n must be an even number and the 2i-th
 * sequence and the (2i+1)-th sequence constitute a read pair. In this
 * mode, there should be enough (typically >50) unique pairs for the
 * routine to infer the orientation and insert size.
 *
 * @param opt    alignment parameters
 * @param bwt    FM-index of the reference sequence
 * @param bns    Information of the reference
 * @param pac    2-bit encoded reference
 * @param n      number of query sequences
 * @param seqs   query sequences; $seqs[i].seq/sam to be modified after the call
 * @param pes0   insert-size info; if NULL, infer from data; if not NULL, it should be an array with 4 elements,
 *               corresponding to each FF, FR, RF and RR orientation. See mem_pestat() for more info.
 */
void mem_process_seqs(mem_opt_t *opt, int64_t n_processed,
                      int n, bseq1_t *seqs, const mem_pestat_t *pes0,
                      worker_t &w);


/**
 * Generate CIGAR and forward-strand position from alignment region
 *
 * @param opt    alignment parameters
 * @param bns    Information of the reference
 * @param pac    2-bit encoded reference
 * @param l_seq  length of query sequence
 * @param seq    query sequence
 * @param ar     one alignment region
 *
 * @return       CIGAR, strand, mapping quality and forward-strand position
 */
mem_aln_t mem_reg2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
                      int l_seq, const char *seq, const mem_alnreg_t *ar);


/**
 * Infer the insert size distribution from interleaved alignment regions
 *
 * This function can be called after mem_align1(), as long as paired-end
 * reads are properly interleaved.
 *
 * @param opt    alignment parameters
 * @param l_pac  length of concatenated reference sequence
 * @param n      number of query sequences; must be an even number
 * @param regs   region array of size $n; 2i-th and (2i+1)-th elements constitute a pair
 * @param pes    inferred insert size distribution (output)
 */
void mem_pestat(const mem_opt_t *opt, int64_t l_pac, int n, const mem_alnreg_v *regs,
                mem_pestat_t pes[4]);

void mem_reorder_primary5(int T, mem_alnreg_v *a);

#endif
