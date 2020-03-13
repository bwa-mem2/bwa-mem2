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

#include "bwamem.h"
#include "FMI_search.h"

//----------------
extern uint64_t proc_freq, tprof[LIM_R][LIM_C];
extern int nthreads;
extern int num_ranks, myrank;
extern FMI_search *fmi;
extern uint8_t *ref_string;       //defined in fastmap.cpp
extern int readLen;
extern int64_t nreads, memSize;
//----------------
#include "kbtree.h"

#define chain_cmp(a, b) (((b).pos < (a).pos) - ((a).pos < (b).pos))
KBTREE_INIT(chn, mem_chain_t, chain_cmp)
			
#define intv_lt(a, b) ((a).info < (b).info)
KSORT_INIT(mem_intv, bwtintv_t, intv_lt)
#define intv_lt1(a, b) ((((uint64_t)(a).m) <<32 | ((uint64_t)(a).n)) < (((uint64_t)(b).m) <<32 | ((uint64_t)(b).n)))  // trial
KSORT_INIT(mem_intv1, SMEM, intv_lt1)  // debug
			
#define max_(x, y) ((x)>(y)?(x):(y))
#define min_(x, y) ((x)>(y)?(y):(x))
			
#define MAX_BAND_TRY  2
			
			int tcnt = 0;
/********************
 * Filtering chains *
 ********************/

#define chn_beg(ch) ((ch).seeds->qbeg)
#define chn_end(ch) ((ch).seeds[(ch).n-1].qbeg + (ch).seeds[(ch).n-1].len)

#define flt_lt(a, b) ((a).w > (b).w)
KSORT_INIT(mem_flt, mem_chain_t, flt_lt)
//------------------------------------------------------------------
// Alignment: Construct the alignment from a chain *

static inline int cal_max_gap(const mem_opt_t *opt, int qlen)
{
	int l_del = (int)((double)(qlen * opt->a - opt->o_del) / opt->e_del + 1.);
	int l_ins = (int)((double)(qlen * opt->a - opt->o_ins) / opt->e_ins + 1.);
	//int l_del = (int)((double)(qlen * opt->a - opt->o_del) + 1.);
	//int l_ins = (int)((double)(qlen * opt->a - opt->o_ins) + 1.);

	int l = l_del > l_ins? l_del : l_ins;
	l = l > 1? l : 1;
	return l < opt->w<<1? l : opt->w<<1;
}

//------------------------------------------------------------------
// SMEMs
static smem_aux_t *smem_aux_init()
{
	smem_aux_t *a;
	if ((a = (smem_aux_t *) calloc(BATCH_SIZE, sizeof(smem_aux_t))) == NULL) { fprintf(stderr, "ERROR: out of memory %s\n", __func__); exit(1); }
	for (int i=0; i<BATCH_SIZE; i++)
	{
		a[i].tmpv[0] = (bwtintv_v *) calloc(1, sizeof(bwtintv_v));
		a[i].tmpv[1] = (bwtintv_v *) calloc(1, sizeof(bwtintv_v));
		if (!a[i].tmpv[0] || !a[i].tmpv[1]) { fprintf(stderr, "ERROR: out of memory %s\n", __func__); exit(1); }
	}
	return a;
}

static void smem_aux_destroy(smem_aux_t *a)
{
	for (int i=0; i<BATCH_SIZE; i++)
	{
		free(a[i].tmpv[0]->a);
		free(a[i].tmpv[0]);
		free(a[i].tmpv[1]->a);
		free(a[i].tmpv[1]);
		free(a[i].mem.a);
		free(a[i].mem1.a);
	}
	free(a);
}

mem_opt_t *mem_opt_init()
{
	mem_opt_t *o;
	if ((o = (mem_opt_t *) calloc(1, sizeof(mem_opt_t))) == NULL)  { fprintf(stderr, "ERROR: out of memory\n"); exit(1); }
	o->flag = 0;
	o->a = 1; o->b = 4;
	o->o_del = o->o_ins = 6;
	o->e_del = o->e_ins = 1;
	o->w = 100;
	o->T = 30;
	o->zdrop = 100;
	o->pen_unpaired = 17;
	o->pen_clip5 = o->pen_clip3 = 5;

	o->max_mem_intv = 20;

	o->min_seed_len = 19;
	o->split_width = 10;
	o->max_occ = 500;
	o->max_chain_gap = 10000;
	o->max_ins = 10000;
	o->mask_level = 0.50;
	o->drop_ratio = 0.50;
	o->XA_drop_ratio = 0.80;
	o->split_factor = 1.5;
	o->chunk_size = 10000000;
	o->n_threads = 1;
	o->max_XA_hits = 5;
	o->max_XA_hits_alt = 200;
	o->max_matesw = 50;
	o->mask_level_redun = 0.95;
	o->min_chain_weight = 0;
	o->max_chain_extend = 1<<30;
	o->mapQ_coef_len = 50; o->mapQ_coef_fac = log(o->mapQ_coef_len);
	bwa_fill_scmat(o->a, o->b, o->mat);
	return o;
}

/******************************
 * De-overlap single-end hits *
 ******************************/

#define alnreg_slt2(a, b) ((a).re < (b).re)
KSORT_INIT(mem_ars2, mem_alnreg_t, alnreg_slt2)

#define alnreg_slt(a, b) ((a).score > (b).score || ((a).score == (b).score && ((a).rb < (b).rb || ((a).rb == (b).rb && (a).qb < (b).qb))))
KSORT_INIT(mem_ars, mem_alnreg_t, alnreg_slt)

#define alnreg_hlt(a, b)  ((a).score > (b).score || ((a).score == (b).score && ((a).is_alt < (b).is_alt || ((a).is_alt == (b).is_alt && (a).hash < (b).hash))))
KSORT_INIT(mem_ars_hash, mem_alnreg_t, alnreg_hlt)

#define alnreg_hlt2(a, b) ((a).is_alt < (b).is_alt || ((a).is_alt == (b).is_alt && ((a).score > (b).score || ((a).score == (b).score && (a).hash < (b).hash))))
KSORT_INIT(mem_ars_hash2, mem_alnreg_t, alnreg_hlt2)

#define PATCH_MAX_R_BW 0.05f
#define PATCH_MIN_SC_RATIO 0.90f

int mem_patch_reg(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
				  uint8_t *query, const mem_alnreg_t *a, const mem_alnreg_t *b,
				  int *_w)
{
	int w, score, q_s, r_s;
	double r;
	if (bns == 0 || pac == 0 || query == 0) return 0;
	assert(a->rid == b->rid && a->rb <= b->rb);
	
	if (a->rb < bns->l_pac && b->rb >= bns->l_pac) return 0; // on different strands
	if (a->qb >= b->qb || a->qe >= b->qe || a->re >= b->re) return 0; // not colinear
	w = (a->re - b->rb) - (a->qe - b->qb); // required bandwidth
	w = w > 0? w : -w; // l = abs(l)
	r = (double)(a->re - b->rb) / (b->re - a->rb) - (double)(a->qe - b->qb) / (b->qe - a->qb); // relative bandwidth
	r = r > 0.? r : -r; // r = fabs(r)

	if (bwa_verbose >= 4)
		fprintf(stderr, "* potential hit merge between [%d,%d)<=>[%ld,%ld) and "
			   "[%d,%d)<=>[%ld,%ld), @ %s; w=%d, r=%.4g\n",
			   a->qb, a->qe, (long)a->rb, (long)a->re, b->qb, b->qe,
			   (long)b->rb, (long)b->re, bns->anns[a->rid].name, w, r);
	
	if (a->re < b->rb || a->qe < b->qb) // no overlap on query or on ref
	{
		if (w > opt->w<<1 || r >= PATCH_MAX_R_BW) return 0; // the bandwidth or the relative bandwidth is too large
	} else if (w > opt->w<<2 || r >= PATCH_MAX_R_BW*2) return 0; // more permissive if overlapping on both ref and query
	
	// global alignment
	w += a->w + b->w;
	w = w < opt->w<<2? w : opt->w<<2;

	if (bwa_verbose >= 4)
		fprintf(stderr, "* test potential hit merge with global alignment; w=%d\n", w);
	
	bwa_gen_cigar2(opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, w,
				   bns->l_pac, pac, b->qe - a->qb, query + a->qb, a->rb, b->re,
				   &score, 0, 0);
	
	q_s = (int)((double)(b->qe - a->qb) / ((b->qe - b->qb) + (a->qe - a->qb)) *
				(b->score + a->score) + .499); // predicted score from query
	
	r_s = (int)((double)(b->re - a->rb) / ((b->re - b->rb) + (a->re - a->rb)) *
				(b->score + a->score) + .499); // predicted score from ref
	
	if (bwa_verbose >= 4)
		fprintf(stderr, "* score=%d;(%d,%d)\n", score, q_s, r_s);
	
	if ((double)score / (q_s > r_s? q_s : r_s) < PATCH_MIN_SC_RATIO) return 0;
	*_w = w;
	return score;
}
/*********************************
 * Test if a seed is good enough *
 *********************************/

#define MEM_SHORT_EXT 50
#define MEM_SHORT_LEN 200

#define MEM_HSP_COEF 1.1f
#define MEM_MINSC_COEF 5.5f
#define MEM_SEEDSW_COEF 0.05f
int stat;

int mem_sort_dedup_patch(const mem_opt_t *opt, const bntseq_t *bns,
						 const uint8_t *pac, uint8_t *query, int n,
						 mem_alnreg_t *a)
{
	int m, i, j;
	if (n <= 1) return n;
	ks_introsort(mem_ars2, n, a); // sort by the END position, not START!
	
	for (i = 0; i < n; ++i) a[i].n_comp = 1;
	for (i = 1; i < n; ++i)
	{
		mem_alnreg_t *p = &a[i];
		if (p->rid != a[i-1].rid || p->rb >= a[i-1].re + opt->max_chain_gap)
			continue; // then no need to go into the loop below
		
		for (j = i - 1; j >= 0 && p->rid == a[j].rid && p->rb < a[j].re + opt->max_chain_gap; --j) {
			mem_alnreg_t *q = &a[j];
			int64_t or_, oq, mr, mq;
			int score, w;
			if (q->qe == q->qb) continue; // a[j] has been excluded
			or_ = q->re - p->rb; // overlap length on the reference
			oq = q->qb < p->qb? q->qe - p->qb : p->qe - q->qb; // overlap length on the query
			mr = q->re - q->rb < p->re - p->rb? q->re - q->rb : p->re - p->rb; // min ref len in alignment
			mq = q->qe - q->qb < p->qe - p->qb? q->qe - q->qb : p->qe - p->qb; // min qry len in alignment
			if (or_ > opt->mask_level_redun * mr && oq > opt->mask_level_redun * mq) { // one of the hits is redundant
				if (p->score < q->score)
				{
					p->qe = p->qb;
					break;
				}
				else q->qe = q->qb;
			}
			else if (q->rb < p->rb && (score = mem_patch_reg(opt, bns, pac, query, q, p, &w)) > 0) { // then merge q into p
				p->n_comp += q->n_comp + 1;
				p->seedcov = p->seedcov > q->seedcov? p->seedcov : q->seedcov;
				p->sub = p->sub > q->sub? p->sub : q->sub;
				p->csub = p->csub > q->csub? p->csub : q->csub;
				p->qb = q->qb, p->rb = q->rb;
				p->truesc = p->score = score;
				p->w = w;
				q->qb = q->qe;
			}
		}
	}
	for (i = 0, m = 0; i < n; ++i) // exclude identical hits
		if (a[i].qe > a[i].qb) {
			if (m != i) a[m++] = a[i];
			else ++m;
		}
	n = m;
	ks_introsort(mem_ars, n, a);
	for (i = 1; i < n; ++i) { // mark identical hits
		if (a[i].score == a[i-1].score && a[i].rb == a[i-1].rb && a[i].qb == a[i-1].qb)
			a[i].qe = a[i].qb;
	}
	for (i = 1, m = 1; i < n; ++i) // exclude identical hits
		if (a[i].qe > a[i].qb) {
			if (m != i) a[m++] = a[i];
			else ++m;
		}
	return m;
}

int mem_sort_dedup_patch_rev(const mem_opt_t *opt, const bntseq_t *bns,
							 const uint8_t *pac, uint8_t *query, int n,
							 mem_alnreg_t *a)
{
	int m, i, j;
	if (n <= 1) return n;
	ks_introsort(mem_ars2, n, a); // sort by the END position, not START!
	
	for (i = 0; i < n; ++i) a[i].n_comp = 1;
	for (i = 1; i < n; ++i)
	{
		mem_alnreg_t *p = &a[i];
		if (p->rid != a[i-1].rid || p->rb >= a[i-1].re + opt->max_chain_gap)
			continue; // then no need to go into the loop below
		
		for (j = i - 1; j >= 0 && p->rid == a[j].rid && p->rb < a[j].re + opt->max_chain_gap; --j) {
			mem_alnreg_t *q = &a[j];
			int64_t or_, oq, mr, mq;
			int score, w;
			if (q->qe == q->qb) continue; // a[j] has been excluded
			or_ = q->re - p->rb; // overlap length on the reference
			oq = q->qb < p->qb? q->qe - p->qb : p->qe - q->qb; // overlap length on the query
			mr = q->re - q->rb < p->re - p->rb? q->re - q->rb : p->re - p->rb; // min ref len in alignment
			mq = q->qe - q->qb < p->qe - p->qb? q->qe - q->qb : p->qe - p->qb; // min qry len in alignment
			if (or_ > opt->mask_level_redun * mr && oq > opt->mask_level_redun * mq) { // one of the hits is redundant
				if (p->score < q->score)
				{
					p->qe = p->qb;
					break;
				}
				else q->qe = q->qb;
			}
			else if (q->rb < p->rb && (score = mem_patch_reg(opt, bns, pac, query, q, p, &w)) > 0) { // then merge q into p
				p->n_comp += q->n_comp + 1;
				p->seedcov = p->seedcov > q->seedcov? p->seedcov : q->seedcov;
				p->sub = p->sub > q->sub? p->sub : q->sub;
				p->csub = p->csub > q->csub? p->csub : q->csub;
				p->qb = q->qb, p->rb = q->rb;
				p->truesc = p->score = score;
				p->w = w;
				q->qb = q->qe;
			}
		}
	}
	for (i = 0, m = 0; i < n; ++i) // exclude identical hits
		if (a[i].qe > a[i].qb) {
			if (m != i) a[m++] = a[i];
			else ++m;
		}
	n = m;
	ks_introsort(mem_ars, n, a);
	for (i = 1; i < n; ++i) { // mark identical hits
		if (a[i].score == a[i-1].score && a[i].rb == a[i-1].rb && a[i].qb == a[i-1].qb) {
			if (a[i].flg  != a[i-1].flg)
			{
				if (a[i].flg == 1)
					a[i].qe = a[i].qb;
				else
					a[i-1].qe = a[i-1].qb;
			}
			else
				a[i].qe = a[i].qb;
		}
	}
	for (i = 1, m = 1; i < n; ++i) // exclude identical hits
		if (a[i].qe > a[i].qb) {
			if (m != i) a[m++] = a[i];
			else ++m;
		}
	return m;
}

// return 1 if the seed is merged into the chain
static int test_and_merge(const mem_opt_t *opt, int64_t l_pac, mem_chain_t *c,
						  const mem_seed_t *p, int seed_rid, int tid)
{
	int64_t qend, rend, x, y;
	const mem_seed_t *last = &c->seeds[c->n-1];
	qend = last->qbeg + last->len;
	rend = last->rbeg + last->len;

	if (seed_rid != c->rid) return 0; // different chr; request a new chain
	if (p->qbeg >= c->seeds[0].qbeg && p->qbeg + p->len <= qend &&
		p->rbeg >= c->seeds[0].rbeg && p->rbeg + p->len <= rend)
		return 1; // contained seed; do nothing
	
	if ((last->rbeg < l_pac || c->seeds[0].rbeg < l_pac) &&
		p->rbeg >= l_pac) return 0; // don't chain if on different strand
	
	x = p->qbeg - last->qbeg; // always non-negtive
	y = p->rbeg - last->rbeg;
	if (y >= 0 && x - y <= opt->w && y - x <= opt->w &&
		x - last->len < opt->max_chain_gap &&
		y - last->len < opt->max_chain_gap) { // grow the chain
		if (c->n == c->m)
		{
			mem_seed_t *auxSeedBuf = NULL;
			int pm = c->m;			
			c->m <<= 1;
			if (pm == SEEDS_PER_CHAIN) {  // re-new memory
				if ((auxSeedBuf = (mem_seed_t *) calloc(c->m, sizeof(mem_seed_t))) == NULL) { fprintf(stderr, "ERROR: out of memory auxSeedBuf\n"); exit(1); }
				memcpy((char*) (auxSeedBuf), c->seeds, c->n * sizeof(mem_seed_t));
				c->seeds = auxSeedBuf;
				tprof[PE13][tid]++;
			} else {  // new memory
				// fprintf(stderr, "[%0.4d] re-allocing old seed, m: %d\n", tid, c->m);
				if ((auxSeedBuf = (mem_seed_t *) realloc(c->seeds, c->m * sizeof(mem_seed_t))) == NULL) { fprintf(stderr, "ERROR: out of memory auxSeedBuf\n"); exit(1); }
				c->seeds = auxSeedBuf;
			}
            memset((char*) (c->seeds + c->n), 0, (c->m - c->n) * sizeof(mem_seed_t));			
		}
		c->seeds[c->n++] = *p;
		return 1;
	}
	return 0; // request to add a new chain
}

int mem_seed_sw(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
				int l_query, const uint8_t *query, const mem_seed_t *s)
{
	int qb, qe, rid;
	int64_t rb, re, mid, l_pac = bns->l_pac;
	uint8_t *rseq = 0;
	kswr_t x;

	if (s->len >= MEM_SHORT_LEN) return -1; // the seed is longer than the max-extend; no need to do SW
	qb = s->qbeg, qe = s->qbeg + s->len;
	rb = s->rbeg, re = s->rbeg + s->len;
	mid = (rb + re) >> 1;
	qb -= MEM_SHORT_EXT; qb = qb > 0? qb : 0;
	qe += MEM_SHORT_EXT; qe = qe < l_query? qe : l_query;
	rb -= MEM_SHORT_EXT; rb = rb > 0? rb : 0;
	re += MEM_SHORT_EXT; re = re < l_pac<<1? re : l_pac<<1;
	if (rb < l_pac && l_pac < re) {
		if (mid < l_pac) re = l_pac;
		else rb = l_pac;
	}
	if (qe - qb >= MEM_SHORT_LEN || re - rb >= MEM_SHORT_LEN) return -1; // the seed seems good enough; no need to do SW

	rseq = bns_fetch_seq(bns, pac, &rb, mid, &re, &rid);
	x = ksw_align2(qe - qb, (uint8_t*)query + qb, re - rb, rseq, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, KSW_XSTART, 0);
	free(rseq);
	return x.score;
}

int mem_chain_weight(const mem_chain_t *c)
{
	int64_t end;
	int j, w = 0, tmp;
	for (j = 0, end = 0; j < c->n; ++j) {
		const mem_seed_t *s = &c->seeds[j];
		if (s->qbeg >= end) w += s->len;
		else if (s->qbeg + s->len > end) w += s->qbeg + s->len - end;
		end = end > s->qbeg + s->len? end : s->qbeg + s->len;
	}
	tmp = w; w = 0;
	for (j = 0, end = 0; j < c->n; ++j) {
		const mem_seed_t *s = &c->seeds[j];
		if (s->rbeg >= end) w += s->len;
		else if (s->rbeg + s->len > end) w += s->rbeg + s->len - end;
		end = end > s->rbeg + s->len? end : s->rbeg + s->len;
	}
	w = w < tmp? w : tmp;
	return w < 1<<30? w : (1<<30)-1;
}

void mem_print_chain(const bntseq_t *bns, mem_chain_v *chn)
{
	int i, j;
	for (i = 0; i < chn->n; ++i)
	{
		mem_chain_t *p = &chn->a[i];
		fprintf(stderr, "* Found CHAIN(%d): n=%d; weight=%d", i, p->n, mem_chain_weight(p));
		for (j = 0; j < p->n; ++j)
		{
			bwtint_t pos;
			int is_rev;
			pos = bns_depos(bns, p->seeds[j].rbeg, &is_rev);
			if (is_rev) pos -= p->seeds[j].len - 1;
			fprintf(stderr, "\t%d;%d;%d,%ld(%s:%c%ld)",
					   p->seeds[j].score, p->seeds[j].len, p->seeds[j].qbeg,
					   (long)p->seeds[j].rbeg, bns->anns[p->rid].name,
					   "+-"[is_rev], (long)(pos - bns->anns[p->rid].offset) + 1);
		}
		fputc('\n', stderr);
	}
}

void mem_flt_chained_seeds(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
						   bseq1_t *seq_, int n_chn, mem_chain_t *a)
{
	//double min_l = opt->min_chain_weight?
	// MEM_HSP_COEF * opt->min_chain_weight : MEM_MINSC_COEF * log(l_query);
	
	int i, j, k;// min_HSP_score = (int)(opt->a * min_l + .499);
	//if (min_l > MEM_SEEDSW_COEF * l_query) return; // don't run the following for short reads

	for (i = 0; i < n_chn; ++i)
	{
		mem_chain_t *c = &a[i];
		const uint8_t *query = (uint8_t*) seq_[c->seqid].seq;
		int l_query = seq_[c->seqid].l_seq;

		double min_l = opt->min_chain_weight?
		MEM_HSP_COEF * opt->min_chain_weight : MEM_MINSC_COEF * log(l_query);
		int min_HSP_score = (int)(opt->a * min_l + .499);
		if (min_l > MEM_SEEDSW_COEF * l_query) continue;
			
		for (j = k = 0; j < c->n; ++j)
		{
			mem_seed_t *s = &c->seeds[j];
			s->score = mem_seed_sw(opt, bns, pac, l_query, query, s);
			if (s->score < 0 || s->score >= min_HSP_score)
			{
				s->score = s->score < 0? s->len * opt->a : s->score;
				c->seeds[k++] = *s;
			}
		}
		c->n = k;
	}
}

int mem_chain_flt(const mem_opt_t *opt, int n_chn_, mem_chain_t *a_, int tid)
{
	int i, k, n_numc = 0;
	if (n_chn_ == 0) return 0; // no need to filter
	// compute the weight of each chain and drop chains with small weight
	for (i = k = 0; i < n_chn_; ++i)
	{
		mem_chain_t *c = &a_[i];
		c->first = -1; c->kept = 0;
		c->w = mem_chain_weight(c);
		if (c->w < opt->min_chain_weight)
        {
			if (c->m > SEEDS_PER_CHAIN) {
				tprof[PE11][tid] ++;
				free(c->seeds);
			}
			//free(c->seeds);
        }
		else a_[k++] = *c;
	}
	n_chn_ = k;
	std::vector<std::pair<int, int> > range;
	std::pair<int, int> pr;
	int pseqid = a_[0].seqid;
	pr.first = 0;
	for (i=1; i<n_chn_; i++)
	{
		mem_chain_t *c =&a_[i];
		if (c->seqid != pseqid) {
			//if (flag == -1) {
			//	pr.first = i;
			//flag = 1;
			//}
			//else
			{
				pr.second = i;
				range.push_back(pr);
				pr.first = i;
				// flag = -1;
			}
		}
		pseqid = c->seqid;
	}
	pr.second = i;
	range.push_back(pr);

	int ilag = 0;
	for (int l=0; l<range.size(); l++)
	{
		// this keeps int indices of the non-overlapping chains
		kvec_t(int) chains = {0,0,0};
		mem_chain_t *a =&a_[range[l].first];
		int n_chn = range[l].second - range[l].first;
		// original code block starts
		ks_introsort(mem_flt, n_chn, a);

        // pairwise chain comparisons
		a[0].kept = 3;
		kv_push(int, chains, 0);
		for (i = 1; i < n_chn; ++i)
		{
			int large_ovlp = 0;
			for (k = 0; k < chains.n; ++k)
			{
				int j = chains.a[k];
				int b_max = chn_beg(a[j]) > chn_beg(a[i])? chn_beg(a[j]) : chn_beg(a[i]);
				int e_min = chn_end(a[j]) < chn_end(a[i])? chn_end(a[j]) : chn_end(a[i]);
				if (e_min > b_max && (!a[j].is_alt || a[i].is_alt)) { // have overlap; don't consider ovlp where the kept chain is ALT while the current chain is primary
					int li = chn_end(a[i]) - chn_beg(a[i]);
					int lj = chn_end(a[j]) - chn_beg(a[j]);
					int min_l = li < lj? li : lj;
					if (e_min - b_max >= min_l * opt->mask_level && min_l < opt->max_chain_gap) { // significant overlap
						large_ovlp = 1;
						if (a[j].first < 0) a[j].first = i; // keep the first shadowed hit s.t. mapq can be more accurate
						if (a[i].w < a[j].w * opt->drop_ratio && a[j].w - a[i].w >= opt->min_seed_len<<1)
							break;
					}
				}
			}
			if (k == chains.n)
			{
				kv_push(int, chains, i);
				a[i].kept = large_ovlp? 2 : 3;
			}
		}
		for (i = 0; i < chains.n; ++i)
		{
			mem_chain_t *c = &a[chains.a[i]];
			if (c->first >= 0) a[c->first].kept = 1;
		}
		free(chains.a);
		for (i = k = 0; i < n_chn; ++i) { // don't extend more than opt->max_chain_extend .kept=1/2 chains
			if (a[i].kept == 0 || a[i].kept == 3) continue;
			if (++k >= opt->max_chain_extend) break;
		}

		for (; i < n_chn; ++i)
			if (a[i].kept < 3) a[i].kept = 0;
	
		for (i = k = 0; i < n_chn; ++i)  // free discarded chains
		{
			mem_chain_t *c = &a[i];
			if (c->kept == 0)
            {
				if (c->m > SEEDS_PER_CHAIN) {
					tprof[PE11][tid] ++;
					free(c->seeds);
				}
                //free(c->seeds);
            }
			else a[k++ - ilag] = a[i];
		}
		// original code block ends
		ilag += n_chn - k;
		n_numc += k;
	}

	return n_numc;
}

SMEM *mem_collect_smem(const mem_opt_t *opt,
					   const bseq1_t *seq_,
					   int nseq,
					   SMEM *matchArray,
					   int32_t *min_intv_ar,
					   int16_t *query_pos_ar,
					   uint8_t *enc_qdb,
					   int32_t *rid,
					   int64_t &tot_smem)
{
	int64_t pos = 0;
	int split_len = (int)(opt->min_seed_len * opt->split_factor + .499);
	int64_t num_smem1 = 0, num_smem2 = 0, num_smem3 = 0;
	int max_readlength = -1;
	
	int32_t *query_cum_len_ar = (int32_t *)_mm_malloc(nseq * sizeof(int32_t), 64);
	
	int offset = 0;
	for (int l=0; l<nseq; l++)
	{
		min_intv_ar[l] = 1;
		for (int j=0; j<seq_[l].l_seq; j++) 
			// enc_qdb[l*seq_[l].l_seq + j] = seq_[l].seq[j];
			enc_qdb[offset + j] = seq_[l].seq[j];
		
		offset += seq_[l].l_seq;		
        rid[l] = l;
	}

	max_readlength = seq_[0].l_seq;
	query_cum_len_ar[0] = 0;
	for(int i = 1; i < nseq; i++) {
        query_cum_len_ar[i] = query_cum_len_ar[i - 1] + seq_[i-1].l_seq;
		if (max_readlength < seq_[i].l_seq)
			max_readlength = seq_[i].l_seq;
	}

	fmi->getSMEMsAllPosOneThread(enc_qdb, min_intv_ar, rid, nseq, nseq,
								 seq_, query_cum_len_ar, max_readlength, opt->min_seed_len,
								 matchArray, &num_smem1);


	for (int64_t i=0; i<num_smem1; i++)
	{
		SMEM *p = &matchArray[i];
		int start = p->m, end = p->n +1;
		if (end - start < split_len || p->s > opt->split_width)
			continue;
		
		int len = seq_[p->rid].l_seq;

		rid[pos] = p->rid;
		
		query_pos_ar[pos] = (end + start)>>1;

		assert(query_pos_ar[pos] < len);

		min_intv_ar[pos] = p->s + 1;
		pos ++;
	}

	fmi->getSMEMsOnePosOneThread(enc_qdb,
								 query_pos_ar,
								 min_intv_ar,
								 rid,
								 pos,
								 pos,
								 seq_,
								 query_cum_len_ar,
								 max_readlength,
								 opt->min_seed_len,
								 matchArray + num_smem1,
								 &num_smem2);

	if (opt->max_mem_intv > 0)
	{
		for (int l=0; l<nseq; l++)
			min_intv_ar[l] = opt->max_mem_intv;

		num_smem3 = fmi->bwtSeedStrategyAllPosOneThread(enc_qdb, min_intv_ar,
														nseq, seq_, query_cum_len_ar, 
														opt->min_seed_len + 1,
														matchArray + num_smem1 + num_smem2);		
	}
	tot_smem = num_smem1 + num_smem2 + num_smem3;

	fmi->sortSMEMs(matchArray, &tot_smem, nseq, seq_[0].l_seq, 1); // seq_[0].l_seq - only used for blocking when using nthreads

	pos = 0;
	int64_t smem_ptr = 0;
	for (int l=0; l<nseq && pos < tot_smem - 1; l++) {
		pos = smem_ptr - 1;
		do {
			pos++;
		} while (pos < tot_smem - 1 && matchArray[pos].rid == matchArray[pos + 1].rid);
		int64_t n = pos + 1 - smem_ptr;
		
		if (n > 0)
			ks_introsort(mem_intv1, n, &matchArray[smem_ptr]);
		smem_ptr = pos + 1;
	}

	_mm_free(query_cum_len_ar);
	return matchArray;
}

/** NEW ONE **/
void mem_chain_seeds(const mem_opt_t *opt,
					 const bntseq_t *bns,
					 const bseq1_t *seq_,
					 int nseq,
					 int tid,
					 mem_chain_v *chain_ar,
                     mem_seed_t *seedBuf,
                     int64_t seedBufSize,
					 SMEM *matchArray,
					 int64_t num_smem)
{
	int b, e, l_rep, size = 0;
    int64_t i, pos = 0;
	int64_t smem_ptr = 0;
	int64_t l_pac = bns->l_pac;

	int num[nseq];
	memset(num, 0, nseq*sizeof(int));
	int64_t *sa_coord = (int64_t *) _mm_malloc(sizeof(int64_t) * opt->max_occ, 64);
    int64_t seedBufCount = 0;
    int64_t auxSeedBufCount = 0;
	
	for (int l=0; l<nseq; l++)
		kv_init(chain_ar[l]);
	
	// filter seq at early stage than this!, shifted to collect!!!
	// if (len < opt->min_seed_len) return chain; // if the query is shorter than the seed length, no match
	
	uint64_t tim = __rdtsc();
	for (int l=0; l<nseq && pos < num_smem - 1; l++)
	{
		// addition, FIX FIX FIX!!!!! THIS!!!!
		//if (aux_->mem.n == 0) continue;

		if (matchArray[smem_ptr].rid > l) continue;
		if (seq_[l].l_seq < opt->min_seed_len) continue;
		assert(matchArray[smem_ptr].rid == l);

		kbtree_t(chn) *tree;
		tree = kb_init(chn, KB_DEFAULT_SIZE + 8); // +8, due to addition of counters in chain
		mem_chain_v *chain = &chain_ar[l];
		size = 0;
			
		b = e = l_rep = 0;
		pos = smem_ptr - 1;
		//for (i = 0, b = e = l_rep = 0; i < aux_->mem.n; ++i) // compute frac_rep
		do
		{
			pos ++;
			SMEM *p = &matchArray[pos];
			int sb = p->m, se = p->n + 1;
			if (p->s <= opt->max_occ) continue;
			if (sb > e) l_rep += e - b, b = sb, e = se;
			else e = e > se? e : se;
		} while (pos < num_smem - 1 && matchArray[pos].rid == matchArray[pos + 1].rid);
		l_rep += e - b;

		for (i = smem_ptr; i <= pos; i++)
		{
			SMEM *p = &matchArray[i];
			int64_t step;
            int32_t count, slen = p->n + 1 - p->m; // seed length
			int64_t k;
			// if (slen < opt->min_seed_len) continue; // ignore if too short or too repetitive
			step = p->s > opt->max_occ? p->s / opt->max_occ : 1;

			// uint64_t tim = __rdtsc();
			int cnt = 0;			
            fmi->get_sa_entries(p, sa_coord, &cnt, 1, opt->max_occ);
			cnt = 0;
			// tprof[MEM_SA][tid] += __rdtsc() - tim;
			
			for (k = count = 0; k < p->s && count < opt->max_occ; k += step, ++count)
			{
				mem_chain_t tmp, *lower, *upper;
				mem_seed_t s;
				int rid, to_add = 0;
								
				s.rbeg = tmp.pos = sa_coord[cnt++];
				s.qbeg = p->m;
				s.score= s.len = slen;
				if (s.rbeg < 0 || s.len < 0) 
					fprintf(stderr, "rbeg: %ld, slen: %d, cnt: %d, n: %d, m: %d, num_smem: %ld\n",
						   s.rbeg, s.len, cnt-1, p->n, p->m, num_smem);
				
				rid = bns_intv2rid(bns, s.rbeg, s.rbeg + s.len);
				// bridging multiple reference sequences or the
				// forward-reverse boundary; TODO: split the seed;
				// don't discard it!!!
				if (rid < 0) continue; 
				if (kb_size(tree))
				{
					kb_intervalp(chn, tree, &tmp, &lower, &upper); // find the closest chain

					if (!lower || !test_and_merge(opt, l_pac, lower, &s, rid, tid))
						to_add = 1;
				}
				else to_add = 1;

				//uint64_t tim = __rdtsc();
				if (to_add) // add the seed as a new chain
				{
					tmp.n = 1; tmp.m = SEEDS_PER_CHAIN;
                    if((seedBufCount + tmp.m) > seedBufSize)
                    {
                        // fprintf(stderr, "ERROR! seedBuf count exceeds size. count = %ld, "
						// 	   "tmp.m = %d, size = %ld\n",
						// 	   seedBufCount, tmp.m, seedBufSize);
						// fprintf(stderr, "Reserved memory allocated for storing seeds has "
						// 		"falling short!!!\n");
						// fprintf(stderr, "Please increase the value of macros: "
						// 	   "AVG_SEEDS_PER_READ & AVG_AUX_SEEDS_PER_READ, "
						// 	   "in src/macro.h, re-compile and re-run.\n");
                        // exit(0);
						// tmp.m <<= 1;
						tmp.m += 1;
						// fprintf(stderr, "[%0.4d] Callocing new seed..\n", tid);
						tmp.seeds = (mem_seed_t *)calloc (tmp.m, sizeof(mem_seed_t));
						tprof[PE13][tid]++;
                    }
					else {
						tmp.seeds = seedBuf + seedBufCount;
						seedBufCount += tmp.m;
					}
                    memset((char*) (tmp.seeds), 0, tmp.m * sizeof(mem_seed_t));
					tmp.seeds[0] = s;
					tmp.rid = rid;
					tmp.seqid = l;
					tmp.is_alt = !!bns->anns[rid].is_alt;
					kb_putp(chn, tree, &tmp);
					num[l]++;
				}
			}
		} // seeds

		smem_ptr = pos + 1;		
		size = kb_size(tree);
		// tprof[PE21][0] += kb_size(tree) * sizeof(mem_chain_t);
		
		//kv_resize(mem_chain_t, *chain, kb_size(tree));
		kv_resize(mem_chain_t, *chain, size);

#define traverse_func(p_) (chain->a[chain->n++] = *(p_))
		__kb_traverse(mem_chain_t, tree, traverse_func);
#undef traverse_func

		for (i = 0; i < chain->n; ++i)
			chain->a[i].frac_rep = (float)l_rep / seq_[l].l_seq;

		kb_destroy(chn, tree);		
		
	} // iterations over input reads
	tprof[MEM_SA_BLOCK][tid] += __rdtsc() - tim;

	_mm_free(sa_coord);
}

int mem_kernel1_core(const mem_opt_t *opt,
					 const bntseq_t *bns,
					 const uint8_t *pac,
					 bseq1_t *seq_,
					 int nseq,
					 mem_chain_v *chain_ar,
                     mem_seed_t *seedBuf,
                     int64_t seedBufSize,
					 mem_cache *mmc,
					 int tid)
{
	SMEM *matchArray = mmc->matchArray[tid];
	int32_t *min_intv_ar = mmc->min_intv_ar[tid];
	int16_t *query_pos_ar = mmc->query_pos_ar[tid];
	uint8_t *enc_qdb = mmc->enc_qdb[tid];
	int32_t *rid = mmc->rid[tid];
	
	int i;
	int64_t num_smem = 0;
 	mem_chain_v *chn;
	
	uint64_t tim;
	// tim = __rdtsc();
	/* convert to 2-bit encoding if we have not done so */
	for (int l=0; l<nseq; l++)
	{
		char *seq = seq_[l].seq;
		int len = seq_[l].l_seq;
		// kv_init(regs[l]);
		
		for (i = 0; i < len; ++i)
			seq[i] = seq[i] < 4? seq[i] : nst_nt4_table[(int)seq[i]]; //nst_nt4??		
	}
	// tprof[MEM_ALN_M2][tid] += __rdtsc() - tim;
	tim = __rdtsc();	
	/********************** Kernel 1: FM+SMEMs *************************/
	printf_(VER, "6. Calling mem_collect_smem.., tid: %d\n", tid);
	//SMEM *matchArray = mem_collect_smem(opt, seq_, nseq, num_smem);
	mem_collect_smem(opt,
					 seq_,
					 nseq,
					 matchArray,
					 min_intv_ar,
					 query_pos_ar,
					 enc_qdb,
					 rid,
					 num_smem);

	if (num_smem >= BATCH_MUL * BATCH_SIZE * readLen){
		fprintf(stderr, "num_smem: %ld\n", num_smem);
		assert(num_smem < BATCH_MUL * BATCH_SIZE * readLen);
	}
	printf_(VER, "6. Done! mem_collect_smem, num_smem: %ld\n", num_smem);
	tprof[MEM_COLLECT][tid] += __rdtsc() - tim;	


	/********************* Kernel 1.1: SA2REF **********************/
	printf_(VER, "6.1. Calling mem_chain..\n");
	//mem_chain_seeds(opt, bwt, bns, seq_, neseq, aux, tid, chain_ar);
	mem_chain_seeds(opt, fmi->idx->bns,
					seq_, nseq, tid,
					chain_ar,
                    seedBuf,
                    seedBufSize,
                    matchArray,
					num_smem);
	
	printf_(VER, "5. Done mem_chain..\n");
	// tprof[MEM_CHAIN][tid] += __rdtsc() - tim;

	/************** Post-processing of collected smems/chains ************/
	// tim = __rdtsc();
	printf_(VER, "6.1. Calling mem_chain_flt..\n");
	for (int l=0; l<nseq; l++)
	{
		chn = &chain_ar[l];
		chn->n = mem_chain_flt(opt, chn->n, chn->a, tid);
	}
	printf_(VER, "7. Done mem_chain_flt..\n");
	// tprof[MEM_ALN_M1][tid] += __rdtsc() - tim;

	
	printf_(VER, "8. Calling mem_flt_chained_seeds..\n");
	for (int l=0; l<nseq; l++) {
		chn = &chain_ar[l];
		//mem_flt_chained_seeds(opt, bns, pac, seq_, chn->n, chn->a);
		mem_flt_chained_seeds(opt, fmi->idx->bns, fmi->idx->pac, seq_, chn->n, chn->a);
	}
	printf_(VER, "8. Done mem_flt_chained_seeds..\n");
	// tprof[MEM_ALN_M2][tid] += __rdtsc() - tim;


	return 1;
}

int mem_kernel2_core(const mem_opt_t *opt,
					 const bntseq_t *bns,
					 const uint8_t *pac,
					 bseq1_t *seq_,
					 mem_alnreg_v *regs,
					 int nseq,
					 mem_chain_v *chain_ar,
					 mem_cache *mmc,
					 int64_t offset1,
					 int64_t offset2,
					 int64_t offset3,
					 int tid)
{
	int i;
	for (int l=0; l<nseq; l++)
	{
		kv_init(regs[l]);
	}
	/****************** Kernel 2: B-SWA *********************/
	uint64_t tim = __rdtsc();
	printf_(VER, "9. Calling mem_chain2aln...\n");
	mem_chain2aln_across_reads_V2(opt,
								  fmi->idx->bns,
								  fmi->idx->pac,
								  seq_,
								  nseq,
								  chain_ar,
								  regs,
								  mmc,
								  offset1,
								  offset2,
								  offset3,
								  tid);

	printf_(VER, "9. Done mem_chain2aln...\n\n");
	tprof[MEM_ALN2][tid] += __rdtsc() - tim;

	// tim = __rdtsc();
	for (int l=0; l<nseq; l++) {
		mem_chain_v *chain = &chain_ar[l];
		for (int i = 0; i < chain->n; ++i) {
			mem_chain_t chn = chain->a[i];
			if (chn.m > SEEDS_PER_CHAIN) {
				tprof[PE11][tid] ++;
				free(chn.seeds);
			}
			tprof[PE12][tid]++;
		}
		free(chain_ar[l].a);
	}
	
	int m = 0;
	for (int l=0; l<nseq; l++)
	{
		mem_alnreg_t *a = regs[l].a;
		int n = regs[l].n;
		for (i = 0, m = 0; i < n; ++i) // exclude identical hits
			if (a[i].qe > a[i].qb) {
				if (m != i) a[m++] = a[i];
				else ++m;
			}
		regs[l].n = m;
	}

	for (int l=0; l<nseq; l++) {		
		regs[l].n = mem_sort_dedup_patch(opt, fmi->idx->bns,
										 fmi->idx->pac,
										 (uint8_t*) seq_[l].seq,
										 regs[l].n, regs[l].a);
	}

	for (int l=0; l<nseq; l++)
	{
		for (i = 0; i < regs[l].n; ++i)
		{
			mem_alnreg_t *p = &regs[l].a[i];
			if (p->rid >= 0 && fmi->idx->bns->anns[p->rid].is_alt)
				p->is_alt = 1;
		}
	}
	// tprof[POST_SWA][tid] += __rdtsc() - tim;
	
	return 1;
}

static void worker_aln(void *data, int seq_id, int batch_size, int tid)
{
	worker_t *w = (worker_t*) data;
	
	printf_(VER, "11. Calling mem_kernel2_core..\n");
	int64_t sizeA = (w->size) *tid;
	int64_t sizeB = (w->size * MAX_SEQ_LEN_REF + MAX_LINE_LEN) * tid;
	int64_t sizeC = (w->size * MAX_SEQ_LEN_QER + MAX_LINE_LEN) * tid;
	
	mem_kernel2_core(w->opt, w->bns, w->pac,
							   w->seqs + seq_id,
							   w->regs + seq_id,
							   batch_size,
							   w->chain_ar + seq_id,
							   &w->mmc, sizeA, sizeB, 
							   sizeC, tid);
	printf_(VER, "11. Done mem_kernel2_core....\n");

}

/* Kernel, called by threads */
static void worker_bwt(void *data, int seq_id, int batch_size, int tid)
{
	worker_t *w = (worker_t*) data;
	printf_(VER, "4. Calling mem_kernel1_core..%d %d\n", seq_id, tid);
	// uint64_t offset = BATCH_SIZE * readLen * tid * BATCH_MUL;
	int seedBufSz = w->seedBufSize;

	if (batch_size < BATCH_SIZE) {
		seedBufSz = memSize - seq_id;
		// fprintf(stderr, "[%0.4d] Info: adjusted seedBufSz %d\n", tid, seedBufSz);
	}

	mem_kernel1_core(w->opt, w->bns, w->pac,
					 w->seqs + seq_id,
					 batch_size,
					 w->chain_ar + seq_id,
					 w->seedBuf + seq_id * AVG_SEEDS_PER_READ,
					 seedBufSz,
					 &(w->mmc),
					 tid);
	printf_(VER, "4. Done mem_kernel1_core....\n");
}

int64_t sort_classify(mem_cache *mmc, int offset1, int64_t pcnt, int tid)
{

#if 0
	SeqPair *seqPairArray = mmc->seqPairArrayLeft128 + offset1;
	// SeqPair *seqPairArrayAux = mmc->seqPairArrayAux + offset1;
#else
	SeqPair *seqPairArray = mmc->seqPairArrayLeft128[tid];
	// SeqPair *seqPairArrayAux = mmc->seqPairArrayAux[tid];
	SeqPair *seqPairArrayAux = mmc->seqPairArrayRight128[tid];
#endif

	int64_t pos8 = 0, pos16 = 0;
	for (int i=0; i<pcnt; i++) {
		SeqPair *s = seqPairArray + i;
		int xtra = s->h0;
		int size = (xtra & KSW_XBYTE)? 1 : 2;
		if (size == 1) // 8
		{
			seqPairArray[pos8++] = seqPairArray[i];
		} else { // 16
			seqPairArrayAux[pos16++] = seqPairArray[i];
		}
	}
	assert(pos8 + pos16 == pcnt);
	
	for (int i=pos8; i<pcnt; i++) {
		seqPairArray[i] = seqPairArrayAux[i-pos8];
	}

	return pos8;
}

static void worker_sam(void *data, int seqid, int batch_size, int tid)
{
	worker_t *w = (worker_t*) data;
	int64_t sizeA = (w->size) *tid;
	int64_t sizeB = (w->size * MAX_SEQ_LEN_REF + MAX_LINE_LEN) * tid;
	int64_t sizeC = (w->size * MAX_SEQ_LEN_QER + MAX_LINE_LEN) * tid;
	
	if (w->opt->flag & MEM_F_PE)
	{
		int64_t pcnt = 0;
		int start = seqid;
		int end = seqid + batch_size;
		int pos = start >> 1;
		
#if (((!__AVX512BW__) && (!__AVX2__)) || ((!__AVX512BW__) && (__AVX2__))) 
		for (int i=start; i< end; i+=2)
		{
			// orig mem_sam_pe() function
			mem_sam_pe(w->opt, w->bns,
					   w->pac, w->pes,
					   (w->n_processed >> 1) + pos++,   // check!
					   &w->seqs[i],
					   &w->regs[i]);
			
			free(w->regs[i].a);
			free(w->regs[i+1].a);
		}
#else   // re-structured
		// pre-processing
		// uint64_t tim = __rdtsc();
		int32_t maxRefLen = 0, maxQerLen = 0;
		int32_t gcnt = 0;
		for (int i=start; i< end; i+=2)
		{
			mem_sam_pe_batch_pre(w->opt, w->bns,
								 w->pac, w->pes,
								 (w->n_processed >> 1) + pos++,   // check!
								 &w->seqs[i],
								 &w->regs[i],
								 &w->mmc, sizeA, sizeB, sizeC,
								 pcnt, gcnt,
								 maxRefLen,
								 maxQerLen,
								 tid);
		}
		
		// tprof[SAM1][tid] += __rdtsc() - tim;
		int64_t pcnt8 = sort_classify(&w->mmc, sizeA, pcnt, tid);

		kswr_t *aln = (kswr_t *) _mm_malloc ((pcnt + SIMD_WIDTH8) * sizeof(kswr_t), 64);
		assert(aln != NULL);

		// processing
		mem_sam_pe_batch(w->opt, &w->mmc, sizeA, sizeB, sizeC,
						 pcnt, pcnt8, aln, maxRefLen, maxQerLen,
						 tid);		

		// post-processing
		// tim = __rdtsc();
		gcnt = 0;
		pos = start >> 1;
		kswr_t *myaln = aln;
		for (int i=start; i< end; i+=2)
		{
			mem_sam_pe_batch_post(w->opt, w->bns,
								  w->pac, w->pes,
								  (w->n_processed >> 1) + pos++,   // check!
								  &w->seqs[i],
								  &w->regs[i],
								  &myaln,
								  &w->mmc, sizeA, sizeB, sizeC,
								  gcnt,
								  tid);

			free(w->regs[i].a);
			free(w->regs[i+1].a);
		}
		//tprof[SAM3][tid] += __rdtsc() - tim;		
		_mm_free(aln);  // kswr_t
#endif
	}
	else
	{
		for (int i=seqid; i<seqid + batch_size; i++)
		{
			mem_mark_primary_se(w->opt, w->regs[i].n,
								w->regs[i].a,
								w->n_processed + i);
			
			mem_reg2sam(w->opt, w->bns, w->pac, &w->seqs[i],
						&w->regs[i], 0, 0);
			free(w->regs[i].a);
		}
	}
}

void mem_process_seqs(mem_opt_t *opt,
					  const bntseq_t *bns,
					  const uint8_t *pac,
					  int64_t n_processed,
					  int n,
					  bseq1_t *seqs,
					  const mem_pestat_t *pes0,
					  worker_t &w)
{
	// worker_t w;
	mem_pestat_t pes[4];
	double ctime, rtime;
	
	ctime = cputime(); rtime = realtime();
	// global_bns = bns;
	w.opt = opt;
	//w.bwt = bwt;
	w.bns = bns; w.pac = pac;
	w.seqs = seqs; w.n_processed = n_processed;
	w.pes = &pes[0];

	//int n_ = (opt->flag & MEM_F_PE) ? n : n;   // this requires n%2==0
	int n_ = n;
	
	uint64_t tim = __rdtsc();	
	fprintf(stderr, "[%0.4d] 3. Calling kt_for - worker_bwt\n", myrank);
	
	kt_for(worker_bwt, &w, n_); // SMEMs (+SAL)

	fprintf(stderr, "[%0.4d] 3. Calling kt_for - worker_aln\n", myrank);
	
	kt_for(worker_aln, &w, n_); // BSW
	tprof[WORKER10][0] += __rdtsc() - tim;		


#if PAIRED_END
	if (opt->flag & MEM_F_PE) { // infer insert sizes if not provided
		if (pes0)
			memcpy(pes, pes0, 4 * sizeof(mem_pestat_t)); // if pes0 != NULL, set the insert-size
		                                                 // distribution as pes0
		else {
			fprintf(stderr, "Inferring insert size distribution from data, l_pac: %ld, n: %d\n", bns->l_pac, n);
			mem_pestat(opt, bns->l_pac, n, w.regs, pes); // otherwise, infer the insert size
		                                                 // distribution from data
		}
	}
#endif
	
	tim = __rdtsc();
	fprintf(stderr, "[%0.4d] 10. Calling kt_for - worker_sam\n", myrank);
	
	kt_for(worker_sam, &w,  n_);   // SAM	
  	tprof[WORKER20][0] += __rdtsc() - tim;

	fprintf(stderr, "\t[%0.4d][ M::%s] Processed %d reads in %.3f "
			"CPU sec, %.3f real sec\n", myrank,
			__func__, n, cputime() - ctime, realtime() - rtime);

}

static void mem_mark_primary_se_core(const mem_opt_t *opt, int n, mem_alnreg_t *a, int_v *z)
{ // similar to the loop in mem_chain_flt()
	int i, k, tmp;
	tmp = opt->a + opt->b;
	tmp = opt->o_del + opt->e_del > tmp? opt->o_del + opt->e_del : tmp;
	tmp = opt->o_ins + opt->e_ins > tmp? opt->o_ins + opt->e_ins : tmp;
	z->n = 0;
	kv_push(int, *z, 0);
	for (i = 1; i < n; ++i) {
		for (k = 0; k < z->n; ++k) {
			int j = z->a[k];
			int b_max = a[j].qb > a[i].qb? a[j].qb : a[i].qb;
			int e_min = a[j].qe < a[i].qe? a[j].qe : a[i].qe;
			if (e_min > b_max) { // have overlap
				int min_l = a[i].qe - a[i].qb < a[j].qe - a[j].qb? a[i].qe - a[i].qb : a[j].qe - a[j].qb;
				if (e_min - b_max >= min_l * opt->mask_level) { // significant overlap
					if (a[j].sub == 0) a[j].sub = a[i].score;
					if (a[j].score - a[i].score <= tmp && (a[j].is_alt || !a[i].is_alt))
						++a[j].sub_n;
					break;
				}
			}
		}
		if (k == z->n) kv_push(int, *z, i);
		else a[i].secondary = z->a[k];
	}
}

int mem_mark_primary_se(const mem_opt_t *opt, int n, mem_alnreg_t *a, int64_t id)
{
	int i, n_pri;
	int_v z = {0,0,0};
	if (n == 0) return 0;
	
	for (i = n_pri = 0; i < n; ++i)
	{
		a[i].sub = a[i].alt_sc = 0, a[i].secondary = a[i].secondary_all = -1, a[i].hash = hash_64(id+i);
		if (!a[i].is_alt) ++n_pri;
	}
	ks_introsort(mem_ars_hash, n, a);
	mem_mark_primary_se_core(opt, n, a, &z);
	for (i = 0; i < n; ++i)
	{
		mem_alnreg_t *p = &a[i];
		p->secondary_all = i; // keep the rank in the first round
		if (!p->is_alt && p->secondary >= 0 && a[p->secondary].is_alt)
			p->alt_sc = a[p->secondary].score;
	}
	if (n_pri >= 0 && n_pri < n)
	{
		kv_resize(int, z, n);
		if (n_pri > 0) ks_introsort(mem_ars_hash2, n, a);
		for (i = 0; i < n; ++i) z.a[a[i].secondary_all] = i;
		for (i = 0; i < n; ++i)
		{
			if (a[i].secondary >= 0)
			{
				a[i].secondary_all = z.a[a[i].secondary];
				if (a[i].is_alt) a[i].secondary = INT_MAX;
			} else a[i].secondary_all = -1;
		}
		if (n_pri > 0) { // mark primary for hits to the primary assembly only
			for (i = 0; i < n_pri; ++i) a[i].sub = 0, a[i].secondary = -1;
			mem_mark_primary_se_core(opt, n_pri, a, &z);
		}
	}
	else {
		for (i = 0; i < n; ++i)
			a[i].secondary_all = a[i].secondary;
	}
	free(z.a);
	return n_pri;
}

/************************
 * Integrated interface *
 ************************/

int mem_approx_mapq_se(const mem_opt_t *opt, const mem_alnreg_t *a)
{
	int mapq, l, sub = a->sub? a->sub : opt->min_seed_len * opt->a;
	double identity;
	sub = a->csub > sub? a->csub : sub;
	if (sub >= a->score) return 0;
	l = a->qe - a->qb > a->re - a->rb? a->qe - a->qb : a->re - a->rb;
	identity = 1. - (double)(l * opt->a - a->score) / (opt->a + opt->b) / l;
	if (a->score == 0) {
		mapq = 0;
	} else if (opt->mapQ_coef_len > 0) {
		double tmp;
		tmp = l < opt->mapQ_coef_len? 1. : opt->mapQ_coef_fac / log(l);
		tmp *= identity * identity;
		mapq = (int)(6.02 * (a->score - sub) / opt->a * tmp * tmp + .499);
	} else {
		mapq = (int)(MEM_MAPQ_COEF * (1. - (double)sub / a->score) * log(a->seedcov) + .499);
		mapq = identity < 0.95? (int)(mapq * identity * identity + .499) : mapq;
	}
	if (a->sub_n > 0) mapq -= (int)(4.343 * log(a->sub_n+1) + .499);
	if (mapq > 60) mapq = 60;
	if (mapq < 0) mapq = 0;
	mapq = (int)(mapq * (1. - a->frac_rep) + .499);
	return mapq;
}

// TODO (future plan): group hits into a uint64_t[] array. This will be cleaner and more flexible
void mem_reg2sam(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
				 bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m)
{
	kstring_t str;
	kvec_t(mem_aln_t) aa;
	int k, l;
	char **XA = 0;

	if (!(opt->flag & MEM_F_ALL))
		XA = mem_gen_alt(opt, bns, pac, a, s->l_seq, s->seq);
	
	kv_init(aa);
	str.l = str.m = 0; str.s = 0;
	for (k = l = 0; k < a->n; ++k)
	{
		mem_alnreg_t *p = &a->a[k];
		mem_aln_t *q;
		if (p->score < opt->T) continue;
		//fprintf(stderr, "%d %d %d %d\n", p->secondary, p->is_alt, opt->flag&MEM_F_ALL,
		//		(p->secondary >= 0 && (p->is_alt || !(opt->flag&MEM_F_ALL))));
		if (p->secondary >= 0 && (p->is_alt || !(opt->flag&MEM_F_ALL))) continue;
	    // assert(p->secondary < INT_MAX);
		if (p->secondary >= 0 && p->secondary < INT_MAX && p->score < a->a[p->secondary].score * opt->drop_ratio) continue;
		q = kv_pushp(mem_aln_t, aa);

		*q = mem_reg2aln(opt, bns, pac, s->l_seq, s->seq, p);
		assert(q->rid >= 0); // this should not happen with the new code
		q->XA = XA? XA[k] : 0;		
		q->flag |= extra_flag; // flag secondary
		if (p->secondary >= 0) q->sub = -1; // don't output sub-optimal score
		if (l && p->secondary < 0) // if supplementary
			q->flag |= (opt->flag&MEM_F_NO_MULTI)? 0x10000 : 0x800;
		if (l && !p->is_alt && q->mapq > aa.a[0].mapq) q->mapq = aa.a[0].mapq;
		++l;
	}
	if (aa.n == 0) { // no alignments good enough; then write an unaligned record
		mem_aln_t t;
		t = mem_reg2aln(opt, bns, pac, s->l_seq, s->seq, 0);
		t.flag |= extra_flag;
		mem_aln2sam(opt, bns, &str, s, 1, &t, 0, m);		
	} else {
		for (k = 0; k < aa.n; ++k)
			mem_aln2sam(opt, bns, &str, s, aa.n, aa.a, k, m);
		for (k = 0; k < aa.n; ++k) free(aa.a[k].cigar);
		free(aa.a);
	}
	s->sam = str.s;
	if (XA) {
		for (k = 0; k < a->n; ++k) free(XA[k]);
		free(XA);
	}
}

void mem_aln2sam(const mem_opt_t *opt, const bntseq_t *bns, kstring_t *str,
				 bseq1_t *s, int n, const mem_aln_t *list, int which, const mem_aln_t *m_)
{	
	int i, l_name;
	mem_aln_t ptmp = list[which], *p = &ptmp, mtmp, *m = 0; // make a copy of the alignment to convert

	if (m_) mtmp = *m_, m = &mtmp;
	// set flag
	p->flag |= m? 0x1 : 0; // is paired in sequencing
	p->flag |= p->rid < 0? 0x4 : 0; // is mapped
	p->flag |= m && m->rid < 0? 0x8 : 0; // is mate mapped
	if (p->rid < 0 && m && m->rid >= 0) // copy mate to alignment
		p->rid = m->rid, p->pos = m->pos, p->is_rev = m->is_rev, p->n_cigar = 0;
	if (m && m->rid < 0 && p->rid >= 0) // copy alignment to mate
		m->rid = p->rid, m->pos = p->pos, m->is_rev = p->is_rev, m->n_cigar = 0;
	p->flag |= p->is_rev? 0x10 : 0; // is on the reverse strand
	p->flag |= m && m->is_rev? 0x20 : 0; // is mate on the reverse strand

	// print up to CIGAR
	l_name = strlen(s->name);
	ks_resize(str, str->l + s->l_seq + l_name + (s->qual? s->l_seq : 0) + 20);
	kputsn(s->name, l_name, str); kputc('\t', str); // QNAME
	kputw((p->flag&0xffff) | (p->flag&0x10000? 0x100 : 0), str); kputc('\t', str); // FLAG
	if (p->rid >= 0) { // with coordinate
		kputs(bns->anns[p->rid].name, str); kputc('\t', str); // RNAME
		kputl(p->pos + 1, str); kputc('\t', str); // POS
		kputw(p->mapq, str); kputc('\t', str); // MAPQ
		if (p->n_cigar) { // aligned
			for (i = 0; i < p->n_cigar; ++i) {
				int c = p->cigar[i]&0xf;
				if (!(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt && (c == 3 || c == 4))
					c = which? 4 : 3; // use hard clipping for supplementary alignments
				kputw(p->cigar[i]>>4, str); kputc("MIDSH"[c], str);
			}
		} else kputc('*', str); // having a coordinate but unaligned (e.g. when copy_mate is true)
	} else kputsn("*\t0\t0\t*", 7, str); // without coordinte
	kputc('\t', str);

	// print the mate position if applicable
	if (m && m->rid >= 0) {
		if (p->rid == m->rid) kputc('=', str);
		else kputs(bns->anns[m->rid].name, str);
		kputc('\t', str);
		kputl(m->pos + 1, str); kputc('\t', str);
		if (p->rid == m->rid) {
			int64_t p0 = p->pos + (p->is_rev? get_rlen(p->n_cigar, p->cigar) - 1 : 0);
			int64_t p1 = m->pos + (m->is_rev? get_rlen(m->n_cigar, m->cigar) - 1 : 0);
			if (m->n_cigar == 0 || p->n_cigar == 0) kputc('0', str);
			else kputl(-(p0 - p1 + (p0 > p1? 1 : p0 < p1? -1 : 0)), str);
		} else kputc('0', str);
	} else kputsn("*\t0\t0", 5, str);
	kputc('\t', str);

	// print SEQ and QUAL
	if (p->flag & 0x100) { // for secondary alignments, don't write SEQ and QUAL
		kputsn("*\t*", 3, str);
	} else if (!p->is_rev) { // the forward strand
		int i, qb = 0, qe = s->l_seq;
		if (p->n_cigar && which && !(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt) { // have cigar && not the primary alignment && not softclip all
			if ((p->cigar[0]&0xf) == 4 || (p->cigar[0]&0xf) == 3) qb += p->cigar[0]>>4;
			if ((p->cigar[p->n_cigar-1]&0xf) == 4 || (p->cigar[p->n_cigar-1]&0xf) == 3) qe -= p->cigar[p->n_cigar-1]>>4;
		}
		ks_resize(str, str->l + (qe - qb) + 1);
		for (i = qb; i < qe; ++i) str->s[str->l++] = "ACGTN"[(int)s->seq[i]];
		kputc('\t', str);
		if (s->qual) { // printf qual
			ks_resize(str, str->l + (qe - qb) + 1);
			for (i = qb; i < qe; ++i) str->s[str->l++] = s->qual[i];
			str->s[str->l] = 0;
		} else kputc('*', str);
	} else { // the reverse strand
		int i, qb = 0, qe = s->l_seq;
		if (p->n_cigar && which && !(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt) {
			if ((p->cigar[0]&0xf) == 4 || (p->cigar[0]&0xf) == 3) qe -= p->cigar[0]>>4;
			if ((p->cigar[p->n_cigar-1]&0xf) == 4 || (p->cigar[p->n_cigar-1]&0xf) == 3) qb += p->cigar[p->n_cigar-1]>>4;
		}
		ks_resize(str, str->l + (qe - qb) + 1);
		for (i = qe-1; i >= qb; --i) str->s[str->l++] = "TGCAN"[(int)s->seq[i]];
		kputc('\t', str);
		if (s->qual) { // printf qual
			ks_resize(str, str->l + (qe - qb) + 1);
			for (i = qe-1; i >= qb; --i) str->s[str->l++] = s->qual[i];
			str->s[str->l] = 0;
		} else kputc('*', str);
	}

	// print optional tags
	if (p->n_cigar) {
		kputsn("\tNM:i:", 6, str); kputw(p->NM, str);
		kputsn("\tMD:Z:", 6, str); kputs((char*)(p->cigar + p->n_cigar), str);
	}
	if (p->score >= 0) { kputsn("\tAS:i:", 6, str); kputw(p->score, str); }
	if (p->sub >= 0) { kputsn("\tXS:i:", 6, str); kputw(p->sub, str); }
	if (bwa_rg_id[0]) { kputsn("\tRG:Z:", 6, str); kputs(bwa_rg_id, str); }
	if (!(p->flag & 0x100)) { // not multi-hit
		for (i = 0; i < n; ++i)
			if (i != which && !(list[i].flag&0x100)) break;
		if (i < n) { // there are other primary hits; output them
			kputsn("\tSA:Z:", 6, str);
			for (i = 0; i < n; ++i) {
				const mem_aln_t *r = &list[i];
				int k;
				if (i == which || (r->flag&0x100)) continue; // proceed if: 1) different from the current; 2) not shadowed multi hit
				kputs(bns->anns[r->rid].name, str); kputc(',', str);
				kputl(r->pos+1, str); kputc(',', str);
				kputc("+-"[r->is_rev], str); kputc(',', str);
				for (k = 0; k < r->n_cigar; ++k) {
					kputw(r->cigar[k]>>4, str); kputc("MIDSH"[r->cigar[k]&0xf], str);
				}
				kputc(',', str); kputw(r->mapq, str);
				kputc(',', str); kputw(r->NM, str);
				kputc(';', str);
			}
		}
		if (p->alt_sc > 0)
			ksprintf(str, "\tpa:f:%.3f", (double)p->score / p->alt_sc);
	}
#if 1
	if (p->XA) { kputsn("\tXA:Z:", 6, str); kputs(p->XA, str); }
#else
	//if (p->XA) { kputsn("\tXA:Z:", 6, str); kputs(p->XA, str);
	//	// printf("%d %s\n", strlen(p->XA), p->XA);
	//	char *t = strstr(p->XA, "-161728");
	//	if (t) {
	//		printf(">> %s\n\n", p->XA);
	//		exit(0);
	//	}
	//}
#endif
	
	if (s->comment) { kputc('\t', str); kputs(s->comment, str); }
	if ((opt->flag&MEM_F_REF_HDR) && p->rid >= 0 && bns->anns[p->rid].anno != 0 && bns->anns[p->rid].anno[0] != 0) {
		int tmp;
		kputsn("\tXR:Z:", 6, str);
		tmp = str->l;
		kputs(bns->anns[p->rid].anno, str);
		for (i = tmp; i < str->l; ++i) // replace TAB in the comment to SPACE
			if (str->s[i] == '\t') str->s[i] = ' ';
	}
	kputc('\n', str);
}

mem_aln_t mem_reg2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const char *query_, const mem_alnreg_t *ar)
{
	mem_aln_t a;
	int i, w2, tmp, qb, qe, NM, score, is_rev, last_sc = -(1<<30), l_MD;
	int64_t pos, rb, re;
	uint8_t *query;

	memset(&a, 0, sizeof(mem_aln_t));
	if (ar == 0 || ar->rb < 0 || ar->re < 0) { // generate an unmapped record
		a.rid = -1; a.pos = -1; a.flag |= 0x4;
		return a;
	}
	qb = ar->qb, qe = ar->qe;
	rb = ar->rb, re = ar->re;
	query = (uint8_t*) malloc(l_query);
	for (i = 0; i < l_query; ++i) // convert to the nt4 encoding
		query[i] = query_[i] < 5? query_[i] : nst_nt4_table[(int)query_[i]];
	a.mapq = ar->secondary < 0? mem_approx_mapq_se(opt, ar) : 0;
	if (ar->secondary >= 0) a.flag |= 0x100; // secondary alignment
	tmp = infer_bw(qe - qb, re - rb, ar->truesc, opt->a, opt->o_del, opt->e_del);
	w2  = infer_bw(qe - qb, re - rb, ar->truesc, opt->a, opt->o_ins, opt->e_ins);
	w2 = w2 > tmp? w2 : tmp;
	if (bwa_verbose >= 4) fprintf(stderr, "* Band width: inferred=%d, cmd_opt=%d, alnreg=%d\n", w2, opt->w, ar->w);
	if (w2 > opt->w) w2 = w2 < ar->w? w2 : ar->w;
	i = 0; a.cigar = 0;
	do {
		free(a.cigar);
		w2 = w2 < opt->w<<2? w2 : opt->w<<2;
		a.cigar = bwa_gen_cigar2(opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, w2, bns->l_pac, pac, qe - qb, (uint8_t*)&query[qb], rb, re, &score, &a.n_cigar, &NM);
		if (bwa_verbose >= 4) fprintf(stderr, "* Final alignment: w2=%d, global_sc=%d, local_sc=%d\n", w2, score, ar->truesc);
		if (score == last_sc || w2 == opt->w<<2) break; // it is possible that global alignment and local alignment give different scores
		last_sc = score;
		w2 <<= 1;
	} while (++i < 3 && score < ar->truesc - opt->a);
	l_MD = strlen((char*)(a.cigar + a.n_cigar)) + 1;
	a.NM = NM;
	pos = bns_depos(bns, rb < bns->l_pac? rb : re - 1, &is_rev);
	a.is_rev = is_rev;
	if (a.n_cigar > 0) { // squeeze out leading or trailing deletions
		assert(a.cigar != NULL);
		if ((a.cigar[0]&0xf) == 2) {
			pos += a.cigar[0]>>4;
			--a.n_cigar;
			memmove(a.cigar, a.cigar + 1, a.n_cigar * 4 + l_MD);
		} else if ((a.cigar[a.n_cigar-1]&0xf) == 2) {
			--a.n_cigar;
			memmove(a.cigar + a.n_cigar, a.cigar + a.n_cigar + 1, l_MD); // MD needs to be moved accordingly
		}
	}
	if (qb != 0 || qe != l_query) { // add clipping to CIGAR
		int clip5, clip3;
		clip5 = is_rev? l_query - qe : qb;
		clip3 = is_rev? qb : l_query - qe;
		a.cigar = (uint32_t*) realloc(a.cigar, 4 * (a.n_cigar + 2) + l_MD);
		if (clip5) {
			memmove(a.cigar+1, a.cigar, a.n_cigar * 4 + l_MD); // make room for 5'-end clipping
			a.cigar[0] = clip5<<4 | 3;
			++a.n_cigar;
		}
		if (clip3) {
			memmove(a.cigar + a.n_cigar + 1, a.cigar + a.n_cigar, l_MD); // make room for 3'-end clipping
			a.cigar[a.n_cigar++] = clip3<<4 | 3;
		}
	}
	a.rid = bns_pos2rid(bns, pos);
	assert(a.rid == ar->rid);
	a.pos = pos - bns->anns[a.rid].offset;
	a.score = ar->score; a.sub = ar->sub > ar->csub? ar->sub : ar->csub;
	a.is_alt = ar->is_alt; a.alt_sc = ar->alt_sc;
	free(query);
	return a;
}

/*****************************
 * Basic hit->SAM conversion *
 *****************************/

static inline int infer_bw(int l1, int l2, int score, int a, int q, int r)
{
	int w;
	if (l1 == l2 && l1 * a - score < (q + r - a)<<1) return 0; // to get equal alignment length, we need at least two gaps
	w = ((double)((l1 < l2? l1 : l2) * a - score - q) / r + 2.);
	if (w < abs(l1 - l2)) w = abs(l1 - l2);
	return w;
}

static inline int get_rlen(int n_cigar, const uint32_t *cigar)
{
	int k, l;
	for (k = l = 0; k < n_cigar; ++k) {
		int op = cigar[k]&0xf;
		if (op == 0 || op == 2)
			l += cigar[k]>>4;
	}
	return l;
}

/************************ New functions, version2*****************************************/
#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)
//#define _get_pac(pac, l) 0

// NOTE: shift these new version of functions from bntseq.cpp to bntseq.cpp,
// once they are incorporated in the code.

uint8_t *bns_get_seq_v2(int64_t l_pac, const uint8_t *pac, int64_t beg, int64_t end, int64_t *len, uint8_t *seqb)
{
	uint8_t *seq = 0;
	if (end < beg) end ^= beg, beg ^= end, end ^= beg; // if end is smaller, swap
	if (end > l_pac<<1) end = l_pac<<1;
	if (beg < 0) beg = 0;
	if (beg >= l_pac || end <= l_pac) {
		int64_t k, l = 0;
		*len = end - beg;
		assert(end-beg < BATCH_SIZE * SEEDS_PER_READ * sizeof(SeqPair));
		
		//seq = (uint8_t*) malloc(end - beg);
		// seq = seqb;
		if (beg >= l_pac) { // reverse strand
#if 0   // orig			
			int64_t beg_f = (l_pac<<1) - 1 - end;
			int64_t end_f = (l_pac<<1) - 1 - beg;
			for (k = end_f; k > beg_f; --k) {
				seq[l++] = 3 - _get_pac(pac, k);
				assert(seq[l-1] == ref_string[beg + l - 1]);
			}
#else
			seq = ref_string + beg;
#endif
		} else { // forward strand
#if 0
			for (k = beg; k < end; ++k) {
				seq[l++] = _get_pac(pac, k);
				assert(seq[l-1] == ref_string[k]);
			}
#else
			seq = ref_string + beg;
#endif			
		}

	} else *len = 0; // if bridging the forward-reverse boundary, return nothing
	return seq;
}

uint8_t *bns_fetch_seq_v2(const bntseq_t *bns, const uint8_t *pac,
						  int64_t *beg, int64_t mid, int64_t *end, int *rid,
						  uint8_t *seqb)
{
	int64_t far_beg, far_end, len;
	int is_rev;
	uint8_t *seq;

	if (*end < *beg) *end ^= *beg, *beg ^= *end, *end ^= *beg; // if end is smaller, swap
	// if (*beg > mid || mid >= *end)
	//	fprintf(stderr, "%ld %ld %ld\n", *beg, mid, *end);
	assert(*beg <= mid && mid < *end);
	
	*rid = bns_pos2rid(bns, bns_depos(bns, mid, &is_rev));
	far_beg = bns->anns[*rid].offset;
	far_end = far_beg + bns->anns[*rid].len;
	if (is_rev) { // flip to the reverse strand
		int64_t tmp = far_beg;
		far_beg = (bns->l_pac<<1) - far_end;
		far_end = (bns->l_pac<<1) - tmp;
	}
	*beg = *beg > far_beg? *beg : far_beg;
	*end = *end < far_end? *end : far_end;

	seq = bns_get_seq_v2(bns->l_pac, pac, *beg, *end, &len, seqb);
	
	if (seq == 0 || *end - *beg != len) {
		fprintf(stderr, "[E::%s] begin=%ld, mid=%ld, end=%ld, len=%ld, seq=%p, rid=%d, far_beg=%ld, far_end=%ld\n",
				__func__, (long)*beg, (long)mid, (long)*end, (long)len, seq, *rid, (long)far_beg, (long)far_end);
	}
	assert(seq && *end - *beg == len); // assertion failure should never happen

	return seq;
}


inline void sortPairsLenExt(SeqPair *pairArray, int32_t count, SeqPair *tempArray,
							int32_t *hist, int &numPairs128, int &numPairs16,
							int &numPairs1)
{
	int32_t i;
	numPairs128 = numPairs16 = numPairs1 = 0;
    // __m256i zero256 = _mm256_setzero_si256();

	int32_t *hist2 = hist + MAX_SEQ_LEN8;
	int32_t *hist3 = hist + MAX_SEQ_LEN8 + MAX_SEQ_LEN16;
	
    for(i = 0; i <= MAX_SEQ_LEN8 + MAX_SEQ_LEN16; i+=1)
        //_mm256_store_si256((__m256i *)(hist + i), zero256);
		hist[i] = 0;
	
	int *arr = (int*) calloc (count, sizeof(int));
	
    for(i = 0; i < count; i++)
    {
		SeqPair sp = pairArray[i];
		int val = max_(sp.len1, sp.len2);
		int minval = sp.h0 + min_(sp.len1, sp.len2);
		//int maxval = sp.len1;
		//if (maxval < sp.len2) maxval = sp.len2;
		//int minval = sp.h0 + maxval;
		
		if (val < MAX_SEQ_LEN8 && minval < MAX_SEQ_LEN8)
			hist[minval]++;
		else if (val < MAX_SEQ_LEN16 && minval < MAX_SEQ_LEN16)
			hist2[minval] ++;
		else
			hist3[0] ++;

		arr[i] = 0;
    }
	
    int32_t cumulSum = 0;
    for(i = 0; i < MAX_SEQ_LEN8; i++)
    {
        int32_t cur = hist[i];
        hist[i] = cumulSum;
        cumulSum += cur;
    }
    for(i = 0; i < MAX_SEQ_LEN16; i++)
    {
        int32_t cur = hist2[i];
        hist2[i] = cumulSum;
        cumulSum += cur;
    }
	hist3[0] = cumulSum;


	for(i = 0; i < count; i++)
    {
		SeqPair sp = pairArray[i];
		int val = max_(sp.len1, sp.len2);
		int minval = sp.h0 + min_(sp.len1, sp.len2);
		//int maxval = sp.len1;
		//if (maxval < sp.len2) maxval = sp.len2;
		//int minval = sp.h0 + maxval;
		
		if (val < MAX_SEQ_LEN8 && minval < MAX_SEQ_LEN8)
		{
			int32_t pos = hist[minval];
			tempArray[pos] = sp;
			hist[minval]++;
			numPairs128 ++;
			if (arr[pos] != 0)
			{
				fprintf(stderr, "Error 1: repreat, pos: %d, arr: %d, minval: %d, (%d %d)\n",
					   pos, arr[pos], minval, sp.len1, sp.len2);
			}
			arr[pos] = 1;
		}
		else if (val < MAX_SEQ_LEN16 && minval < MAX_SEQ_LEN16) {
			int32_t pos = hist2[minval];
			tempArray[pos] = sp;
			hist2[minval]++;
			numPairs16 ++;
			if (arr[pos] != 0)
			{
				SeqPair spt = pairArray[arr[pos]-1];
				fprintf(stderr, "Error XX: repeat, "
					   "i: %d, pos: %d, arr: %d, hist2: %d, minval: %d, (%d %d %d) (%d %d %d)\n",
					   i, pos, arr[pos], hist2[minval],  minval, sp.h0, sp.len1, sp.len2,
					   spt.h0, spt.len1, spt.len2);
			}			
			arr[pos] = i + 1;
		}
		else {
			int32_t pos = hist3[0];
			tempArray[pos] = sp;
			hist3[0]++;
			arr[pos] = i + 1;
			numPairs1 ++;
		}
    }
	
    for(i = 0; i < count; i++) {
		pairArray[i] = tempArray[i];
	}

	free(arr);

#if 1  // DEB
	for(i = 0; i < numPairs128; i++) {
		if (pairArray[i].len1 >= 128 || pairArray[i].len2 >= 128 || pairArray[i].h0 >= 128)
			fprintf(stderr, "Error: Not matching..1 %d %d %d\n",
					pairArray[i].len1, pairArray[i].len2, pairArray[i].h0);
	}
	for(i = numPairs128; i < numPairs16; i++) {
		if (pairArray[i].len1 >= MAX_SEQ_LEN16 || pairArray[i].len2 >= MAX_SEQ_LEN16 || pairArray[i].h0 >= MAX_SEQ_LEN16)
			fprintf(stderr, "Error: Not matching..2\n");
	}
#endif
}

inline void sortPairsLen(SeqPair *pairArray, int32_t count, SeqPair *tempArray, int32_t *hist)
{

    int32_t i;
#if ((!__AVX512BW__) & (__AVX2__ | __SSE2__))
    for(i = 0; i <= MAX_SEQ_LEN16; i++) hist[i] = 0;
#else	
    __m512i zero512 = _mm512_setzero_si512();
    for(i = 0; i <= MAX_SEQ_LEN16; i+=16)
    {
        _mm512_store_si512((__m512i *)(hist + i), zero512);
    }
#endif
    
    for(i = 0; i < count; i++)
    {
        SeqPair sp = pairArray[i];
        hist[sp.len1]++;
    }
    int32_t prev = 0;
    int32_t cumulSum = 0;
    for(i = 0; i <= MAX_SEQ_LEN16; i++)
    {
        int32_t cur = hist[i];
        hist[i] = cumulSum;
        cumulSum += cur;
    }

    for(i = 0; i < count; i++)
    {
        SeqPair sp = pairArray[i];
        int32_t pos = hist[sp.len1];

        tempArray[pos] = sp;
        hist[sp.len1]++;
    }

    for(i = 0; i < count; i++) {
        pairArray[i] = tempArray[i];
	}
}

/* Restructured BSW parent function */
#define FAC 8
#define PFD 2
void mem_chain2aln_across_reads_V2(const mem_opt_t *opt, const bntseq_t *bns,
								   const uint8_t *pac, bseq1_t *seq_, int nseq,
								   mem_chain_v* chain_ar, mem_alnreg_v *av_v,
								   mem_cache *mmc, int64_t offset1, int64_t offset2,
								   int64_t offset3, int tid)
{
#if 0
	SeqPair *seqPairArrayAux	  = mmc->seqPairArrayAux + offset1;
	SeqPair *seqPairArrayLeft128  = mmc->seqPairArrayLeft128 + offset1;
	SeqPair *seqPairArrayRight128 = mmc->seqPairArrayRight128 + offset1;
#else
	SeqPair *seqPairArrayAux	  = mmc->seqPairArrayAux[tid];
	SeqPair *seqPairArrayLeft128  = mmc->seqPairArrayLeft128[tid];
	SeqPair *seqPairArrayRight128 = mmc->seqPairArrayRight128[tid];
	int64_t *wsize_pair = &(mmc->wsize[tid]);
	// fprintf(stderr, "wsize: %d\n", *wsize);
#endif

#if 0
	uint8_t *seqBufLeftRef	= mmc->seqBufLeftRef + offset2;
	uint8_t *seqBufRightRef = mmc->seqBufRightRef + offset2;
	uint8_t *seqBufLeftQer	= mmc->seqBufLeftQer + offset3;
	uint8_t *seqBufRightQer = mmc->seqBufRightQer + offset3;
#else
	uint8_t *seqBufLeftRef	= mmc->seqBufLeftRef[tid]; 
	uint8_t *seqBufRightRef = mmc->seqBufRightRef[tid];
	uint8_t *seqBufLeftQer	= mmc->seqBufLeftQer[tid]; 
	uint8_t *seqBufRightQer = mmc->seqBufRightQer[tid];
	int64_t *wsize_buf_ref = &(mmc->wsize_buf_ref[tid]);
	int64_t *wsize_buf_qer = &(mmc->wsize_buf_qer[tid]);
#endif
	
 	// int32_t *lim_g = mmc->lim + (BATCH_SIZE + 32) * tid;
	int32_t *lim_g = mmc->lim[tid];
	
	mem_seed_t *s;
	int64_t l_pac = bns->l_pac, rmax[8] __attribute__((aligned(64)));
	// std::vector<int8_t> nexitv(nseq, 0);
	
	int numPairsLeft = 0, numPairsRight = 0;
	int numPairsLeft1 = 0, numPairsRight1 = 0;
	int numPairsLeft128 = 0, numPairsRight128 = 0;
	int numPairsLeft16 = 0, numPairsRight16 = 0;

	int64_t leftRefOffset = 0, rightRefOffset = 0;
	int64_t leftQerOffset = 0, rightQerOffset = 0;

	int srt_size = MAX_SEEDS_PER_READ, fac = FAC;
	uint64_t *srt = (uint64_t *) malloc(srt_size * 8);
	uint32_t *srtgg = (uint32_t*) malloc(nseq * SEEDS_PER_READ * fac * sizeof(uint32_t));

	int spos = 0;
	
	// uint64_t timUP = __rdtsc();
	for (int l=0; l<nseq; l++)
	{
		int max = 0;
		uint8_t *rseq = 0;
		
		uint32_t *srtg = srtgg;
		lim_g[l+1] = 0;
		
		const uint8_t *query = (uint8_t *) seq_[l].seq;
		int l_query = seq_[l].l_seq;
	
		mem_chain_v *chn = &chain_ar[l];		
		mem_alnreg_v *av = &av_v[l];  // alignment
		mem_chain_t *c;

		_mm_prefetch((const char*) query, _MM_HINT_NTA);
		
		// aln mem allocation
		av->m = 0;
		for (int j=0; j<chn->n; j++) {
			c = &chn->a[j];	av->m += c->n;
		}
		av->a = (mem_alnreg_t*)calloc(av->m, sizeof(mem_alnreg_t));

		// aln mem allocation ends
		for (int j=0; j<chn->n; j++)
		{
			c = &chn->a[j];
			assert(c->seqid == l);
			
			int64_t tmp = 0;
			if (c->n == 0) continue;
			
			_mm_prefetch((const char*) (srtg + spos + 64), _MM_HINT_NTA);
			_mm_prefetch((const char*) (lim_g), _MM_HINT_NTA);
			
			// get the max possible span
			rmax[0] = l_pac<<1; rmax[1] = 0;

			for (int i = 0; i < c->n; ++i) {
				int64_t b, e;
				const mem_seed_t *t = &c->seeds[i];
				b = t->rbeg - (t->qbeg + cal_max_gap(opt, t->qbeg));
				e = t->rbeg + t->len + ((l_query - t->qbeg - t->len) +
										cal_max_gap(opt, l_query - t->qbeg - t->len));

				tmp = rmax[0];
				rmax[0] = tmp < b? rmax[0] : b;
				rmax[1] = (rmax[1] > e)? rmax[1] : e;
				if (t->len > max) max = t->len;
			}
			
			rmax[0] = rmax[0] > 0? rmax[0] : 0;
			rmax[1] = rmax[1] < l_pac<<1? rmax[1] : l_pac<<1;
			if (rmax[0] < l_pac && l_pac < rmax[1])
			{
				if (c->seeds[0].rbeg < l_pac) rmax[1] = l_pac;
				else rmax[0] = l_pac;
			}

			/* retrieve the reference sequence */
			{
				int rid = 0;
				// free rseq
				rseq = bns_fetch_seq_v2(bns, pac, &rmax[0],
										c->seeds[0].rbeg,
										&rmax[1], &rid,
										(uint8_t*) seqPairArrayAux);
				assert(c->rid == rid);
			}

			_mm_prefetch((const char*) rseq, _MM_HINT_NTA);
			// _mm_prefetch((const char*) rseq + 64, _MM_HINT_NTA);
			
			// assert(c->n < MAX_SEEDS_PER_READ);  // temp
			if (c->n > srt_size) {
				srt_size = c->n + 10;
				srt = (uint64_t *) realloc(srt, srt_size * 8);
			}
			
			for (int i = 0; i < c->n; ++i) 
				srt[i] = (uint64_t)c->seeds[i].score<<32 | i;

			if (c->n > 1) 
				ks_introsort_64(c->n, srt);
			
			// assert((spos + c->n) < SEEDS_PER_READ * FAC * nseq);
			if ((spos + c->n) > SEEDS_PER_READ * fac * nseq) {
				fac <<= 1;
				srtgg = (uint32_t *) realloc(srtgg, nseq * SEEDS_PER_READ * fac * sizeof(uint32_t));
			}
			
			for (int i = 0; i < c->n; ++i)
				srtg[spos++] = srt[i];

			lim_g[l+1] += c->n;
			
			// uint64_t tim = __rdtsc();
			for (int k=c->n-1; k >= 0; k--)
			{
				s = &c->seeds[(uint32_t)srt[k]];

				mem_alnreg_t *a;
				// a = kv_pushp(mem_alnreg_t, *av);
				a = &av->a[av->n++];
				memset(a, 0, sizeof(mem_alnreg_t));
				
				s->aln = av->n-1;
				
				a->w = opt->w;
				a->score = a->truesc = -1;
				a->rid = c->rid;
				a->frac_rep = c->frac_rep;
				a->seedlen0 = s->len;
				a->c = c; //ptr
				a->rb = a->qb = a->re = a->qe = H0_;
				// a->secondary_all = 0;
				
				tprof[PE19][tid] ++;
				
				int flag = 0;
				std::pair<int, int> pr;
				if (s->qbeg)  // left extension
				{
					SeqPair sp;
					sp.h0 = s->len * opt->a;
					sp.seqid = c->seqid;
					sp.regid = av->n - 1;
						
					// assert(numPairsLeft < BATCH_SIZE * SEEDS_PER_READ);
					if (numPairsLeft >= *wsize_pair) {
						fprintf(stderr, "[LOG][%0.4d] Re-allocating seqPairArrays, in Left\n", tid);
						*wsize_pair += 1000;
						seqPairArrayAux = (SeqPair *) realloc(seqPairArrayAux,
															  (*wsize_pair + MAX_LINE_LEN)
															  * sizeof(SeqPair));
						mmc->seqPairArrayAux[tid] = seqPairArrayAux;
						seqPairArrayLeft128 = (SeqPair *) realloc(seqPairArrayLeft128,
																  (*wsize_pair + MAX_LINE_LEN)
																  * sizeof(SeqPair));
						mmc->seqPairArrayLeft128[tid] = seqPairArrayLeft128;
						seqPairArrayRight128 = (SeqPair *) realloc(seqPairArrayRight128,
																   (*wsize_pair + MAX_LINE_LEN)
																   * sizeof(SeqPair));
						mmc->seqPairArrayRight128[tid] = seqPairArrayRight128;
					}

					
					sp.idq = leftQerOffset;
					sp.idr = leftRefOffset;
					uint8_t *qs = seqBufLeftQer + sp.idq;
					uint8_t *rs = seqBufLeftRef + sp.idr;
					
					leftQerOffset += s->qbeg;
					// assert(leftQerOffset < MAX_SEQ_LEN_QER * BATCH_SIZE * SEEDS_PER_READ);
					if (leftQerOffset >= *wsize_buf_qer) {
						fprintf(stderr, "[%0.4d] Re-allocating (doubling) seqBufQers in %s (left)\n",
								tid, __func__);
						*wsize_buf_qer *= 2;
						uint8_t *seqBufQer_ = realloc(mmc->seqBufLeftQer[tid], *wsize_buf_qer);
						mmc->seqBufLeftQer[tid] = seqBufLeftQer = seqBufQer_;
						
						seqBufQer_ = realloc(mmc->seqBufRightQer[tid], *wsize_buf_qer);
						mmc->seqBufRightQer[tid] = seqBufRightQer = seqBufQer_;		
					}
					

					for (int i = 0; i < s->qbeg; ++i) qs[i] = query[s->qbeg - 1 - i];  // vectorize it!!
					
					tmp = s->rbeg - rmax[0];
					leftRefOffset += tmp;
					// assert(leftRefOffset < MAX_SEQ_LEN_REF * BATCH_SIZE * SEEDS_PER_READ);
					if (leftRefOffset >= *wsize_buf_ref) {
						fprintf(stderr, "[%0.4d] Re-allocating (doubling) seqBufRefs in %s (left)\n",
								tid, __func__);
						*wsize_buf_ref *= 2;
						uint8_t *seqBufRef_ = realloc(mmc->seqBufLeftRef[tid], *wsize_buf_ref);
						mmc->seqBufLeftRef[tid] = seqBufLeftRef = seqBufRef_;
						
						seqBufRef_ = realloc(mmc->seqBufRightRef[tid], *wsize_buf_ref);
						mmc->seqBufRightRef[tid] = seqBufRightRef = seqBufRef_;				
					}
					

					for (int64_t i = 0; i < tmp; ++i) rs[i] = rseq[tmp - 1 - i]; //seq1
					
					sp.len2 = s->qbeg;
					sp.len1 = tmp;
					int minval = sp.h0 + max_(sp.len1, sp.len2);
					// int minval = sp.h0 +  min_(sp.len1, sp.len2);
					//fprintf(fsam, "%d %d %d %d n128: %d n16: %d\n",
					//		sp.len1, sp.len2, minval, MAX_SEQ_LEN8, numPairsLeft128, numPairsLeft16);

					//if (max_(sp.len1, sp.len2) < MAX_SEQ_LEN8 && minval < MAX_SEQ_LEN8) {
					if (sp.len1 < MAX_SEQ_LEN8 && sp.len2 < MAX_SEQ_LEN8 && minval < MAX_SEQ_LEN8) {
						numPairsLeft128++;
					}
					else if (minval < MAX_SEQ_LEN16) {
						numPairsLeft16++;
					}
					else {							
						numPairsLeft1++;
					}
					
					seqPairArrayLeft128[numPairsLeft] = sp;
					numPairsLeft ++;
					a->qb = s->qbeg; a->rb = s->rbeg;
				}
				else
				{
					flag = 1;
					a->score = a->truesc = s->len * opt->a, a->qb = 0, a->rb = s->rbeg;
				}

				if (s->qbeg + s->len != l_query)  // right extension
				{
					int64_t qe = s->qbeg + s->len;
					int64_t re = s->rbeg + s->len - rmax[0];
					assert(re >= 0);					
					SeqPair sp;

					sp.h0 = H0_; //random number
					sp.seqid = c->seqid;
					sp.regid = av->n - 1;

					// assert(numPairsRight < BATCH_SIZE * SEEDS_PER_READ);
					if (numPairsRight >= *wsize_pair) {
						fprintf(stderr, "[LOG] [%0.4d] Re-allocating seqPairArrays Right\n", tid);
						*wsize_pair += 1000;
						seqPairArrayAux = (SeqPair *) realloc(seqPairArrayAux,
															  (*wsize_pair + MAX_LINE_LEN)
															  * sizeof(SeqPair));
						mmc->seqPairArrayAux[tid] = seqPairArrayAux;
						seqPairArrayLeft128 = (SeqPair *) realloc(seqPairArrayLeft128,
																  (*wsize_pair + MAX_LINE_LEN)
																  * sizeof(SeqPair));
						mmc->seqPairArrayLeft128[tid] = seqPairArrayLeft128;
						seqPairArrayRight128 = (SeqPair *) realloc(seqPairArrayRight128,
																   (*wsize_pair + MAX_LINE_LEN)
																   * sizeof(SeqPair));
						mmc->seqPairArrayRight128[tid] = seqPairArrayRight128;
					}

					
					sp.len2 = l_query - qe;
					sp.len1 = rmax[1] - rmax[0] - re;

					sp.idq = rightQerOffset;
					sp.idr = rightRefOffset;
					uint8_t *qs = seqBufRightQer + sp.idq;
					uint8_t *rs = seqBufRightRef + sp.idr;
					
					rightQerOffset += sp.len2;
					// assert(rightQerOffset < MAX_SEQ_LEN_QER * BATCH_SIZE * SEEDS_PER_READ);
					if (rightQerOffset >= *wsize_buf_qer) {
						fprintf(stderr, "[%0.4d] Re-allocating (doubling) seqBufQers in %s (right)\n",
								tid, __func__);
						*wsize_buf_qer *= 2;
						uint8_t *seqBufQer_ = realloc(mmc->seqBufLeftQer[tid], *wsize_buf_qer);
						mmc->seqBufLeftQer[tid] = seqBufLeftQer = seqBufQer_;
						
						seqBufQer_ = realloc(mmc->seqBufRightQer[tid], *wsize_buf_qer);
						mmc->seqBufRightQer[tid] = seqBufRightQer = seqBufQer_;		
					}

					rightRefOffset += sp.len1;
					// assert(rightRefOffset < MAX_SEQ_LEN_REF * BATCH_SIZE * SEEDS_PER_READ);
					if (rightRefOffset >= *wsize_buf_ref) {
						fprintf(stderr, "[%0.4d] Re-allocating (doubling) seqBufRefs in %s (right)\n",
								tid, __func__);
						*wsize_buf_ref *= 2;
						uint8_t *seqBufRef_ = realloc(mmc->seqBufLeftRef[tid], *wsize_buf_ref);
						mmc->seqBufLeftRef[tid] = seqBufLeftRef = seqBufRef_;
						
						seqBufRef_ = realloc(mmc->seqBufRightRef[tid], *wsize_buf_ref);
						mmc->seqBufRightRef[tid] = seqBufRightRef = seqBufRef_;				
					}
					
					tprof[PE23][tid] += sp.len1 + sp.len2;

					
					for (int i = 0; i < sp.len2; ++i) qs[i] = query[qe + i];

					for (int i = 0; i < sp.len1; ++i) rs[i] = rseq[re + i]; //seq1


					int minval = sp.h0 + max_(sp.len1, sp.len2);
					// int minval = sp.h0 +  min_(sp.len1, sp.len2);					
					// if (max_(sp.len1, sp.len2) < MAX_SEQ_LEN8 && minval < MAX_SEQ_LEN8) {
					if (minval < MAX_SEQ_LEN8) {
						numPairsRight128++;
					}
					else if(minval < MAX_SEQ_LEN16) {
						numPairsRight16++;
					}
					else {
						numPairsRight1++;
					}
					seqPairArrayRight128[numPairsRight] = sp;
					numPairsRight ++;
					a->qe = qe; a->re = rmax[0] + re;
				}
				else
				{
					a->qe = l_query, a->re = s->rbeg + s->len;
					// seedcov business, this "if" block should be redundant, check and remove.
					if (a->rb != H0_ && a->qb != H0_)
					{
						int i;
						for (i = 0, a->seedcov = 0; i < c->n; ++i)
						{
							const mem_seed_t *t = &c->seeds[i];
							if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe &&
								t->rbeg >= a->rb && t->rbeg + t->len <= a->re) // seed fully contained
								a->seedcov += t->len;
						}
					}
				}
			}
			// free(rseq);
			// tprof[MEM_ALN2_DOWN1][tid] += __rdtsc() - tim;
		}
	}
	// tprof[MEM_ALN2_UP][tid] += __rdtsc() - timUP;
	

	int32_t *hist = (int32_t *)_mm_malloc((MAX_SEQ_LEN8 + MAX_SEQ_LEN16 + 32) *
										  sizeof(int32_t), 64);
	
	/* Sorting based score is required as that affects the use of SIMD lanes */
	sortPairsLenExt(seqPairArrayLeft128, numPairsLeft, seqPairArrayAux, hist,
				  numPairsLeft128, numPairsLeft16, numPairsLeft1);
	assert(numPairsLeft == (numPairsLeft128 + numPairsLeft16 + numPairsLeft1));
	

	// SWA
	// uint64_t timL = __rdtsc();
	int nthreads = 1;

	// Now, process all the collected seq-pairs
	// First, left alignment, move out these calls
	BandedPairWiseSW bswLeft(opt->o_del, opt->e_del, opt->o_ins,
							 opt->e_ins, opt->zdrop, opt->pen_clip5,
							 opt->mat, opt->a, opt->b, nthreads);
	
	BandedPairWiseSW bswRight(opt->o_del, opt->e_del, opt->o_ins,
							  opt->e_ins, opt->zdrop, opt->pen_clip3,
							  opt->mat, opt->a, opt->b, nthreads);
	
	int i;
	// Left
	SeqPair *pair_ar = seqPairArrayLeft128 + numPairsLeft128 + numPairsLeft16;
	SeqPair *pair_ar_aux = seqPairArrayAux;
	int nump = numPairsLeft1;
	
#if 1
	// scalar
	for ( i=0; i<MAX_BAND_TRY; i++)
	{
		int w = opt->w << i;
		// uint64_t tim = __rdtsc();
		bswLeft.scalarBandedSWAWrapper(pair_ar,
									   seqBufLeftRef,
									   seqBufLeftQer,
									   nump,
									   nthreads,
									   w);
		//tprof[PE5][0] += nump;
		//tprof[PE6][0] ++;
		// tprof[MEM_ALN2_B][tid] += __rdtsc() - tim;
			
		int num = 0;
		for (int l=0; l<nump; l++)
		{
			mem_alnreg_t *a;		
			SeqPair *sp = &pair_ar[l];
			a = &(av_v[sp->seqid].a[sp->regid]); // prev
			//OutScore *o = outScoreArray + l;
			int prev = a->score;
			a->score = sp->score;

			if (a->score == prev || sp->max_off < (w >> 1) + (w >> 2) ||
				i+1 == MAX_BAND_TRY)
			{
				if (sp->gscore <= 0 || sp->gscore <= a->score - opt->pen_clip5) {
					a->qb -= sp->qle; a->rb -= sp->tle;
					a->truesc = a->score;
				} else {
					a->qb = 0; a->rb -= sp->gtle;
					a->truesc = sp->gscore;
				}

				a->w = max_(a->w, w);
				if (a->rb != H0_ && a->qb != H0_ && a->qe != H0_ && a->re != H0_)
				{
					int i = 0;
					for (i = 0, a->seedcov = 0; i < a->c->n; ++i){
						const mem_seed_t *t = &(a->c->seeds[i]);
						if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe &&
							t->rbeg >= a->rb && t->rbeg + t->len <= a->re)
							a->seedcov += t->len;
					}
				}
			
			} else {
				pair_ar_aux[num++] = *sp;
			}
		}
		nump = num;
		SeqPair *tmp = pair_ar;
		pair_ar = pair_ar_aux;
		pair_ar_aux = tmp;
	}


	//****************** Left - vector int16 ***********************
#if 1
	assert(numPairsLeft == (numPairsLeft128 + numPairsLeft16 + numPairsLeft1));

	pair_ar = seqPairArrayLeft128 + numPairsLeft128;
	pair_ar_aux = seqPairArrayAux;
	
	nump = numPairsLeft16;
	for ( i=0; i<MAX_BAND_TRY; i++)
	{
		int w = opt->w << i;
		// int64_t tim = __rdtsc();
#if ((!__AVX512BW__) && (!__AVX2__) && (!__SSE2__))
		bswLeft.scalarBandedSWAWrapper(pair_ar, seqBufLeftRef, seqBufLeftQer, nump, nthreads, w);
#else
		sortPairsLen(pair_ar, nump, seqPairArrayAux, hist);
		bswLeft.getScores16(pair_ar,
							seqBufLeftRef,
							seqBufLeftQer,
							nump,
							nthreads,
							w);
#endif
		
		tprof[PE5][0] += nump;
		tprof[PE6][0] ++;				
		// tprof[MEM_ALN2_B][tid] += __rdtsc() - tim;

		int num = 0;
		for (int l=0; l<nump; l++)
		{			
			mem_alnreg_t *a;
			SeqPair *sp = &pair_ar[l];
			a = &(av_v[sp->seqid].a[sp->regid]); // prev

			int prev = a->score;
			a->score = sp->score;

			
			if (a->score == prev || sp->max_off < (w >> 1) + (w >> 2) ||
				i+1 == MAX_BAND_TRY)
			{
				if (sp->gscore <= 0 || sp->gscore <= a->score - opt->pen_clip5) {
					a->qb -= sp->qle; a->rb -= sp->tle;
					a->truesc = a->score;
				} else {
					a->qb = 0; a->rb -= sp->gtle;
					a->truesc = sp->gscore;
				}

				a->w = max_(a->w, w);
				if (a->rb != H0_ && a->qb != H0_ && a->qe != H0_ && a->re != H0_)
				{
					int i = 0;
					for (i = 0, a->seedcov = 0; i < a->c->n; ++i){
						const mem_seed_t *t = &(a->c->seeds[i]);
						if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe &&
							t->rbeg >= a->rb && t->rbeg + t->len <= a->re)
							a->seedcov += t->len;
					}
				}
			} else {
				pair_ar_aux[num++] = *sp;
			}
		}
		nump = num;
		SeqPair *tmp = pair_ar;
		pair_ar = pair_ar_aux;
		pair_ar_aux = tmp;
	}
#endif

	//****************** Left - vector int8 ***********************
#if 1
	pair_ar = seqPairArrayLeft128;
	pair_ar_aux = seqPairArrayAux;
	
	nump = numPairsLeft128;
	for ( i=0; i<MAX_BAND_TRY; i++)
	{
		int w = opt->w << i;
		// int64_t tim = __rdtsc();
		
#if ((!__AVX512BW__) && (!__AVX2__) && (!__SSE2__))
		bswLeft.scalarBandedSWAWrapper(pair_ar, seqBufLeftRef, seqBufLeftQer, nump, nthreads, w);
#else
		sortPairsLen(pair_ar, nump, seqPairArrayAux, hist);
		bswLeft.getScores8(pair_ar,
						   seqBufLeftRef,
						   seqBufLeftQer,
						   nump,
						   nthreads,
						   w);
#endif	
		
		tprof[PE1][0] += nump;
		tprof[PE2][0] ++;
		// tprof[MEM_ALN2_D][tid] += __rdtsc() - tim;

		int num = 0;
		for (int l=0; l<nump; l++)
		{			
			mem_alnreg_t *a;
			SeqPair *sp = &pair_ar[l];
			a = &(av_v[sp->seqid].a[sp->regid]); // prev

			int prev = a->score;
			a->score = sp->score;
			
			if (a->score == prev || sp->max_off < (w >> 1) + (w >> 2) ||
				i+1 == MAX_BAND_TRY)
			{
				if (sp->gscore <= 0 || sp->gscore <= a->score - opt->pen_clip5) {
					a->qb -= sp->qle; a->rb -= sp->tle;
					a->truesc = a->score;
				} else {
					a->qb = 0; a->rb -= sp->gtle;
					a->truesc = sp->gscore;
				}

				a->w = max_(a->w, w);
				if (a->rb != H0_ && a->qb != H0_ && a->qe != H0_ && a->re != H0_)
				{
					int i = 0;
					for (i = 0, a->seedcov = 0; i < a->c->n; ++i){
						const mem_seed_t *t = &(a->c->seeds[i]);
						if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe &&
							t->rbeg >= a->rb && t->rbeg + t->len <= a->re)
							a->seedcov += t->len;
					}
				}
			} else {
				pair_ar_aux[num++] = *sp;
			}
		}
		nump = num;
		SeqPair *tmp = pair_ar;
		pair_ar = pair_ar_aux;
		pair_ar_aux = tmp;
	}
#endif

	// tprof[CLEFT][tid] += __rdtsc() - timL;
	
	// uint64_t timR = __rdtsc();
	// **********************************************************
	// Right, scalar
	for (int l=0; l<numPairsRight; l++) {
		mem_alnreg_t *a;		
		SeqPair *sp = &seqPairArrayRight128[l];
		a = &(av_v[sp->seqid].a[sp->regid]); // prev
		sp->h0 = a->score;
	}

	sortPairsLenExt(seqPairArrayRight128, numPairsRight, seqPairArrayAux,
					hist, numPairsRight128, numPairsRight16, numPairsRight1);

	assert(numPairsRight == (numPairsRight128 + numPairsRight16 + numPairsRight1));

	pair_ar = seqPairArrayRight128 + numPairsRight128 + numPairsRight16;
	pair_ar_aux = seqPairArrayAux;
	nump = numPairsRight1;

	for ( i=0; i<MAX_BAND_TRY; i++)
	{
		int w = opt->w << i;
		//tim = __rdtsc();		
		bswRight.scalarBandedSWAWrapper(pair_ar,
						seqBufRightRef,
						seqBufRightQer,
						nump,
						nthreads,
						w);
		//tprof[PE7][0] += nump;
		//tprof[PE8][0] ++;
		//tprof[MEM_ALN2_C][tid] += __rdtsc() - tim;
		int num = 0;

		for (int l=0; l<nump; l++)
		{
			mem_alnreg_t *a;		
			SeqPair *sp = &pair_ar[l];
			a = &(av_v[sp->seqid].a[sp->regid]); // prev
			//OutScore *o = outScoreArray + l;
			int prev = a->score;
			a->score = sp->score;
			
			// no further banding
			if (a->score == prev || sp->max_off < (w >> 1) + (w >> 2) ||
				i+1 == MAX_BAND_TRY)
			{
				if (sp->gscore <= 0 || sp->gscore <= a->score - opt->pen_clip3) {
					a->qe += sp->qle, a->re += sp->tle;
					a->truesc += a->score - sp->h0;
				} else {
					int l_query = seq_[sp->seqid].l_seq;
					a->qe = l_query, a->re += sp->gtle;
					a->truesc += sp->gscore - sp->h0;
				}
				a->w = max_(a->w, w);
				if (a->rb != H0_ && a->qb != H0_ && a->qe != H0_ && a->re != H0_)
				{
					int i = 0;
					for (i = 0, a->seedcov = 0; i < a->c->n; ++i) {
						const mem_seed_t *t = &a->c->seeds[i];
						if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe &&
							t->rbeg >= a->rb && t->rbeg + t->len <= a->re)
							a->seedcov += t->len;
					}
				}
			} else {
				pair_ar_aux[num++] = *sp;
			}
		}
		nump = num;
		SeqPair *tmp = pair_ar;
		pair_ar = pair_ar_aux;
		pair_ar_aux = tmp;
	}

	// ************************* Right - vector int16 **********************
	pair_ar = seqPairArrayRight128 + numPairsRight128;
	pair_ar_aux = seqPairArrayAux;
	nump = numPairsRight16;
		
	for ( i=0; i<MAX_BAND_TRY; i++)
	{
		int w = opt->w << i;
		// uint64_t tim = __rdtsc();
#if ((!__AVX512BW__) && (!__AVX2__) && (!__SSE2__))
		bswRight.scalarBandedSWAWrapper(pair_ar, seqBufRightRef, seqBufRightQer, nump, nthreads, w);
#else
		sortPairsLen(pair_ar, nump, seqPairArrayAux, hist);
		bswRight.getScores16(pair_ar,
							 seqBufRightRef,
							 seqBufRightQer,
							 nump,
							 nthreads,
							 w);
#endif

		tprof[PE7][0] += nump;
		tprof[PE8][0] ++;
		// tprof[MEM_ALN2_C][tid] += __rdtsc() - tim;
		
		int num = 0;

		for (int l=0; l<nump; l++)
		{
			mem_alnreg_t *a;		
			SeqPair *sp = &pair_ar[l];
			a = &(av_v[sp->seqid].a[sp->regid]); // prev
			//OutScore *o = outScoreArray + l;
			int prev = a->score;
			a->score = sp->score;
			
			// no further banding
			if (a->score == prev || sp->max_off < (w >> 1) + (w >> 2) ||
				i+1 == MAX_BAND_TRY)
			{
				if (sp->gscore <= 0 || sp->gscore <= a->score - opt->pen_clip3) {
					a->qe += sp->qle, a->re += sp->tle;
					a->truesc += a->score - sp->h0;
				} else {
					int l_query = seq_[sp->seqid].l_seq;
					a->qe = l_query, a->re += sp->gtle;
					a->truesc += sp->gscore - sp->h0;
				}
				a->w = max_(a->w, w);
				if (a->rb != H0_ && a->qb != H0_ && a->qe != H0_ && a->re != H0_)
				{
					int i = 0;
					for (i = 0, a->seedcov = 0; i < a->c->n; ++i) {
						const mem_seed_t *t = &a->c->seeds[i];
						if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe &&
							t->rbeg >= a->rb && t->rbeg + t->len <= a->re)
							a->seedcov += t->len;
					}
				}
			} else {
				pair_ar_aux[num++] = *sp;
			}
		}
		nump = num;
		SeqPair *tmp = pair_ar;
		pair_ar = pair_ar_aux;
		pair_ar_aux = tmp;
	}


	// ************************* Right, vector int8 **********************
	pair_ar = seqPairArrayRight128;
	pair_ar_aux = seqPairArrayAux;
	nump = numPairsRight128;
		
	for ( i=0; i<MAX_BAND_TRY; i++)
	{
		int w = opt->w << i;
		// uint64_t tim = __rdtsc();
		
#if ((!__AVX512BW__) && (!__AVX2__) && (!__SSE2__))
		bswRight.scalarBandedSWAWrapper(pair_ar, seqBufRightRef, seqBufRightQer, nump, nthreads, w); 
#else
		sortPairsLen(pair_ar, nump, seqPairArrayAux, hist);
		bswRight.getScores8(pair_ar,
							seqBufRightRef,
							seqBufRightQer,
							nump,
							nthreads,
							w);
#endif	
			
		tprof[PE3][0] += nump;
		tprof[PE4][0] ++;
		// tprof[MEM_ALN2_E][tid] += __rdtsc() - tim;
		int num = 0;

		for (int l=0; l<nump; l++)
		{
			mem_alnreg_t *a;		
			SeqPair *sp = &pair_ar[l];
			a = &(av_v[sp->seqid].a[sp->regid]); // prev
			//OutScore *o = outScoreArray + l;
			int prev = a->score;
			a->score = sp->score;
			// no further banding
			if (a->score == prev || sp->max_off < (w >> 1) + (w >> 2) ||
				i+1 == MAX_BAND_TRY)
			{
				if (sp->gscore <= 0 || sp->gscore <= a->score - opt->pen_clip3) {
					a->qe += sp->qle, a->re += sp->tle;
					a->truesc += a->score - sp->h0;
				} else {
					int l_query = seq_[sp->seqid].l_seq;
					a->qe = l_query, a->re += sp->gtle;
					a->truesc += sp->gscore - sp->h0;
				}
				a->w = max_(a->w, w);
				if (a->rb != H0_ && a->qb != H0_ && a->qe != H0_ && a->re != H0_)
				{
					int i = 0;
					for (i = 0, a->seedcov = 0; i < a->c->n; ++i) {
						const mem_seed_t *t = &a->c->seeds[i];
						if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe &&
							t->rbeg >= a->rb && t->rbeg + t->len <= a->re)
							a->seedcov += t->len;
					}
				}
			} else {
				pair_ar_aux[num++] = *sp;
			}
		}
		nump = num;
		SeqPair *tmp = pair_ar;
		pair_ar = pair_ar_aux;
		pair_ar_aux = tmp;
	}
#endif
	_mm_free(hist);
	// tprof[CRIGHT][tid] += __rdtsc() - timR;
	
	if (numPairsLeft >= *wsize_pair || numPairsRight >= *wsize_pair)
	{   // refine it!
		fprintf(stderr, "Error: This should not have happened!!!\n");
		fprintf(stderr, "Error: assert failed for seqPair size, "
				"numPairsLeft: %d, numPairsRight %d\nExiting.\n",
				numPairsLeft, numPairsRight);
		exit(EXIT_FAILURE);
	}
	/* Discard seeds and hence their alignemnts */
	
	lim_g[0] = 0;
	for (int l=1; l<nseq; l++)
		lim_g[l] += lim_g[l-1];
			
	// uint64_t tim = __rdtsc();			
	int *lim = (int *) calloc(BATCH_SIZE, sizeof(int));


	for (int l=0; l<nseq; l++)
	{
		int s_start = 0, s_end = 0;
		uint32_t *srtg = srtgg + lim_g[l];
		
		int l_query = seq_[l].l_seq;
		mem_chain_v *chn = &chain_ar[l];		
		mem_alnreg_v *av = &av_v[l];  // alignment
		mem_chain_t *c;
		
		for (int j=0; j<chn->n; j++)
		{
			c = &chn->a[j];
			assert(c->seqid == l);

			s_end = s_start + c->n;

			uint32_t *srt2 = srtg + s_start;
			s_start += c->n;
			
			int k = 0;
			for (k = c->n-1; k >= 0; k--)
			{
				s = &c->seeds[srt2[k]];
				int i = 0;
				int v = 0;
				for (i = 0; i < av->n && v < lim[l]; ++i)  // test whether extension has been made before
				{
					mem_alnreg_t *p = &av->a[i];
					//mem_alnreg_t *p = &av->a[v];
					if (p->qb == -1 && p->qe == -1) {
						continue;
					}

					int64_t rd;
					int qd, w, max_gap;
					if (s->rbeg < p->rb || s->rbeg + s->len > p->re || s->qbeg < p->qb
						|| s->qbeg + s->len > p->qe) {
						v++; continue; // not fully contained
					}
				
					if (s->len - p->seedlen0 > .1 * l_query) { v++; continue;}
					// qd: distance ahead of the seed on query; rd: on reference
					qd = s->qbeg - p->qb; rd = s->rbeg - p->rb;
					// the maximal gap allowed in regions ahead of the seed
					max_gap = cal_max_gap(opt, qd < rd? qd : rd); 
					w = max_gap < p->w? max_gap : p->w; // bounded by the band width
					if (qd - rd < w && rd - qd < w) break; // the seed is "around" a previous hit
					// similar to the previous four lines, but this time we look at the region behind
					qd = p->qe - (s->qbeg + s->len); rd = p->re - (s->rbeg + s->len);
					max_gap = cal_max_gap(opt, qd < rd? qd : rd);
					w = max_gap < p->w? max_gap : p->w;
					if (qd - rd < w && rd - qd < w) break;
					
					v++;
				}
				
				// the seed is (almost) contained in an existing alignment;
				// further testing is needed to confirm it is not leading
				// to a different aln
				//if (i < av->n)
				if (v < lim[l])
				{
					for (v = k + 1; v < c->n; ++v)
					{
						const mem_seed_t *t;
						if (srt2[v] == UINT_MAX) continue;
						//t = &c->seeds[(uint32_t)srt[v]];
						t = &c->seeds[srt2[v]];
						//if (t->done == H0_) continue;  //check for interferences!!!
						// only check overlapping if t is long enough;
						// TODO: more efficient by early stopping
						if (t->len < s->len * .95) continue; 
						if (s->qbeg <= t->qbeg && s->qbeg + s->len - t->qbeg >= s->len>>2 &&
							t->qbeg - s->qbeg != t->rbeg - s->rbeg) break;
						if (t->qbeg <= s->qbeg && t->qbeg + t->len - s->qbeg >= s->len>>2 &&
							s->qbeg - t->qbeg != s->rbeg - t->rbeg) break;
					}
					if (v == c->n) {                  // no overlapping seeds; then skip extension
						mem_alnreg_t *ar = &(av_v[l].a[s->aln]);
						ar->qb = ar->qe = -1;         // purge the alingment
						srt2[k] = UINT_MAX;
						tprof[PE18][tid]++;
						continue;
					}
				}
				
				lim[l]++;
			}
		}
	}
	free(srtgg);
	free(srt);
	free(lim);
	// tprof[MEM_ALN2_DOWN][tid] += __rdtsc() - tim;	
}
