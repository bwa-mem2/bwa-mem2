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

#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include <assert.h>
#include "bntseq.h"
#include "bwa.h"
#include "ksw.h"
#include "utils.h"
#include "kstring.h"
#include "kvec.h"
#include <string>

int bwa_verbose = 3;
char bwa_rg_id[256];
char *bwa_pg;

long long int count_fasta_file(char *file_path, char *str_db);
/************************
 * Batch FASTA/Q reader *
 ************************/

#include "kseq.h"
KSEQ_DECLARE(gzFile)

static inline void trim_readno(kstring_t *s)
{
	if (s->l > 2 && s->s[s->l-2] == '/' && isdigit(s->s[s->l-1]))
		s->l -= 2, s->s[s->l] = 0;
}

static inline void kseq2bseq1(const kseq_t *ks, bseq1_t *s)
{ // TODO: it would be better to allocate one chunk of memory, but probably it does not matter in practice
	s->name = strdup(ks->name.s);
	s->comment = ks->comment.l? strdup(ks->comment.s) : 0;
	s->seq = strdup(ks->seq.s);
	s->qual = ks->qual.l? strdup(ks->qual.s) : 0;
	s->l_seq = strlen(s->seq);
}

/* Customized for MPI processing */
//bseq1_t *bseq_read(int chunk_size, int *n_, void *ks1_, void *ks2_)
bseq1_t *bseq_read(int64_t chunk_size, int *n_, void *ks1_, void *ks2_,
				   FILE* fpp, int len, int64_t *s)
{
	kseq_t *ks = (kseq_t*)ks1_, *ks2 = (kseq_t*)ks2_;
	int64_t size = 0, m, n, size2 = 0;
	bseq1_t *seqs;
	m = n = 0; seqs = 0;
	char buf[len];
	
	while (kseq_read(ks) >= 0)
	{
		if (ks2 && kseq_read(ks2) < 0) { // the 2nd file has fewer reads
			fprintf(stderr, "[W::%s] the 2nd file has fewer sequences.\n", __func__);
			break;
		}
		if (n >= m) {
			m = m? m<<1 : 256;
			seqs = (bseq1_t*) realloc(seqs, m * sizeof(bseq1_t));
		}
		trim_readno(&ks->name);
		kseq2bseq1(ks, &seqs[n]);
		seqs[n].id = n;
		{
			//kseq_t *ksd = ks;
			//kstream_t *kst = ksd->f;
#if 0
			//printf("Check D..\n%s\n%s\n%s\n%s\n",
			//	   seqs[n].name, seqs[n].seq,
			//	   seqs[n].comment, seqs[n].qual);
			
			if (seqs[n].name != NULL)
				size += strlen(seqs[n].name);
			//printf("%d ", strlen(seqs[n].name)+strlen(seqs[n].comment)+1);
			if (seqs[n].comment != NULL) {
				size += strlen(seqs[n].comment);
				std::string str = seqs[n].comment;
				std::size_t found = str.find("length");
				if (found != std::string::npos) {
					size += strlen(seqs[n].comment) + strlen(seqs[n].name) + 1;
				}
			}
			else
				size += 1;
			
			if (seqs[n].qual != NULL)
				size += strlen(seqs[n].qual);

			size += 7; // non accounted chars
#else
			//kstring_t kstr;
			//printf("%d\n", ks_getuntil2(kst, KS_SEP_LINE, &kstr, 0, 0));
			fgets((char*) buf, len, fpp);
			size2 += strlen((char*) buf);
			// printf("First line: %d, %s\n", strlen(buf), buf);
			fgets((char*) buf, len, fpp);
			size2 += strlen((char*) buf);
			if (seqs[n].qual != NULL) {
				fgets((char*) buf, len, fpp);
				size2 += strlen((char*) buf);
				fgets((char*) buf, len, fpp);
				size2 += strlen((char*) buf);
			}
#endif
		}
		//size += seqs[n++].l_seq;
		size = size2;		n++;
		
		//printf("size: %d, size2: %d\n", size, size2);
		//static int cnt = 0;
		//if (cnt++ == 4)exit(0);
		
		if (ks2) {
			trim_readno(&ks2->name);
			kseq2bseq1(ks2, &seqs[n]);
			seqs[n].id = n;
			n++;
			// size += seqs[n++].l_seq;
		}
		//if (size >= chunk_size && (n&1) == 0) break;
		if (size >= chunk_size) {
			break;
		}
	}
	if (size == 0) { // test if the 2nd file is finished
		if (ks2 && kseq_read(ks2) >= 0)
			fprintf(stderr, "[W::%s] the 1st file has fewer sequences.\n", __func__);
	}
	*n_ = n;
	*s = size;
	return seqs;
}

bseq1_t *bseq_read_orig(int64_t chunk_size, int *n_, void *ks1_, void *ks2_, int64_t *s)
{
	kseq_t *ks = (kseq_t*)ks1_, *ks2 = (kseq_t*)ks2_;
	int64_t size = 0, m, n;
	bseq1_t *seqs;
	m = n = 0; seqs = 0;
	while (kseq_read(ks) >= 0)
	{
		if (ks2 && kseq_read(ks2) < 0) { // the 2nd file has fewer reads
			fprintf(stderr, "[W::%s] the 2nd file has fewer sequences.\n", __func__);
			break;
		}
		if (n >= m) {
			m = m? m<<1 : 256;
			seqs = (bseq1_t*) realloc(seqs, m * sizeof(bseq1_t));
		}
		trim_readno(&ks->name);
		kseq2bseq1(ks, &seqs[n]);
		seqs[n].id = n;
		//{
		//	size += strlen(seqs[n].name);
		//	size += strlen(seqs[n].comment);
		//	size += strlen(seqs[n].qual);
		//	// fprintf(stderr, "qual len: %d %d\n", strlen(seqs[n].qual), seqs[n].l_seq);
		//	size += 7; // non accounted chars
		//}
		size += seqs[n++].l_seq;

		if (ks2) {
			trim_readno(&ks2->name);
			kseq2bseq1(ks2, &seqs[n]);
			seqs[n].id = n;
			size += seqs[n++].l_seq;
		}
		if (size >= chunk_size && (n&1) == 0) break;
		// if (size >= chunk_size) {
		//	break;
		// }
	}
	if (size == 0) { // test if the 2nd file is finished
		if (ks2 && kseq_read(ks2) >= 0)
			fprintf(stderr, "[W::%s] the 1st file has fewer sequences.\n", __func__);
	}
	*n_ = n;
	*s = size;
	return seqs;
}

bseq1_t *bseq_read_one_fasta_file(int64_t chunk_size, int *n_, gzFile fp, int64_t *s)
{
    kseq_t *ks = kseq_init(fp);
    return bseq_read_orig(chunk_size, n_, ks, NULL, s);
}

void bseq_classify(int n, bseq1_t *seqs, int m[2], bseq1_t *sep[2])
{
	int i, has_last;
	kvec_t(bseq1_t) a[2] = {{0,0,0}, {0,0,0}};
	for (i = 1, has_last = 1; i < n; ++i) {
		if (has_last) {
			if (strcmp(seqs[i].name, seqs[i-1].name) == 0) {
				kv_push(bseq1_t, a[1], seqs[i-1]);
				kv_push(bseq1_t, a[1], seqs[i]);
				has_last = 0;
			} else kv_push(bseq1_t, a[0], seqs[i-1]);
		} else has_last = 1;
	}
	if (has_last) kv_push(bseq1_t, a[0], seqs[i-1]);
	sep[0] = a[0].a, m[0] = a[0].n;
	sep[1] = a[1].a, m[1] = a[1].n;
}

/*****************
 * CIGAR related *
 *****************/

void bwa_fill_scmat(int a, int b, int8_t mat[25])
{
	int i, j, k;
	for (i = k = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j)
			mat[k++] = i == j? a : -b;
		mat[k++] = -1; // ambiguous base
	}
	for (j = 0; j < 5; ++j) mat[k++] = -1;   // DEFAULT AMBIG
}

// Generate CIGAR when the alignment end points are known
uint32_t *bwa_gen_cigar2(const int8_t mat[25], int o_del, int e_del, int o_ins, int e_ins, int w_, int64_t l_pac, const uint8_t *pac, int l_query, uint8_t *query, int64_t rb, int64_t re, int *score, int *n_cigar, int *NM)
{
	uint32_t *cigar = 0;
	uint8_t tmp, *rseq;
	int i;
	int64_t rlen;
	kstring_t str;
	const char *int2base;

	if (n_cigar) *n_cigar = 0;
	if (NM) *NM = -1;
	if (l_query <= 0 || rb >= re || (rb < l_pac && re > l_pac)) return 0; // reject if negative length or bridging the forward and reverse strand
	rseq = bns_get_seq(l_pac, pac, rb, re, &rlen);
	if (re - rb != rlen) goto ret_gen_cigar; // possible if out of range
	if (rb >= l_pac) { // then reverse both query and rseq; this is to ensure indels to be placed at the leftmost position
		for (i = 0; i < l_query>>1; ++i)
			tmp = query[i], query[i] = query[l_query - 1 - i], query[l_query - 1 - i] = tmp;
		for (i = 0; i < rlen>>1; ++i)
			tmp = rseq[i], rseq[i] = rseq[rlen - 1 - i], rseq[rlen - 1 - i] = tmp;
	}
	if (l_query == re - rb && w_ == 0) { // no gap; no need to do DP
		// UPDATE: we come to this block now... FIXME: due to an issue in mem_reg2aln(), we never come to this block. This does not affect accuracy, but it hurts performance.
		if (n_cigar) {
			cigar = (uint32_t*) malloc(4);
			cigar[0] = l_query<<4 | 0;
			*n_cigar = 1;
		}
		for (i = 0, *score = 0; i < l_query; ++i)
			*score += mat[rseq[i]*5 + query[i]];
	} else {
		int w, max_gap, max_ins, max_del, min_w;
		// set the band-width
		max_ins = (int)((double)(((l_query+1)>>1) * mat[0] - o_ins) / e_ins + 1.);
		max_del = (int)((double)(((l_query+1)>>1) * mat[0] - o_del) / e_del + 1.);
		max_gap = max_ins > max_del? max_ins : max_del;
		max_gap = max_gap > 1? max_gap : 1;
		w = (max_gap + abs(rlen - l_query) + 1) >> 1;
		w = w < w_? w : w_;
		min_w = abs(rlen - l_query) + 3;
		w = w > min_w? w : min_w;
		// NW alignment
		if (bwa_verbose >= 4) {
			fprintf(stderr, "* Global bandwidth: %d\n", w);
			fprintf(stderr, "* Global ref:   "); for (i = 0; i < rlen; ++i) fputc("ACGTN"[(int)rseq[i]], stderr); fputc('\n', stderr);
			fprintf(stderr, "* Global query: "); for (i = 0; i < l_query; ++i) fputc("ACGTN"[(int)query[i]], stderr); fputc('\n', stderr);
		}
		*score = ksw_global2(l_query, query, rlen, rseq, 5, mat, o_del, e_del, o_ins, e_ins, w, n_cigar, &cigar);
	}
	if (NM && n_cigar) {// compute NM and MD
		int k, x, y, u, n_mm = 0, n_gap = 0;
		str.l = str.m = *n_cigar * 4; str.s = (char*)cigar; // append MD to CIGAR
		int2base = rb < l_pac? "ACGTN" : "TGCAN";
		for (k = 0, x = y = u = 0; k < *n_cigar; ++k) {
			int op, len;
			cigar = (uint32_t*)str.s;
			op  = cigar[k]&0xf, len = cigar[k]>>4;
			if (op == 0) { // match
				for (i = 0; i < len; ++i) {
					if (query[x + i] != rseq[y + i]) {
						kputw(u, &str);
						kputc(int2base[rseq[y+i]], &str);
						++n_mm; u = 0;
					} else ++u;
				}
				x += len; y += len;
			} else if (op == 2) { // deletion
				if (k > 0 && k < *n_cigar - 1) { // don't do the following if D is the first or the last CIGAR
					kputw(u, &str); kputc('^', &str);
					for (i = 0; i < len; ++i)
						kputc(int2base[rseq[y+i]], &str);
					u = 0; n_gap += len;
				}
				y += len;
			} else if (op == 1) x += len, n_gap += len; // insertion
		}
		kputw(u, &str); kputc(0, &str);
		*NM = n_mm + n_gap;
		cigar = (uint32_t*)str.s;
	}
	if (rb >= l_pac) // reverse back query
		for (i = 0; i < l_query>>1; ++i)
			tmp = query[i], query[i] = query[l_query - 1 - i], query[l_query - 1 - i] = tmp;

ret_gen_cigar:
	free(rseq);
	return cigar;
}

uint32_t *bwa_gen_cigar(const int8_t mat[25], int q, int r, int w_, int64_t l_pac, const uint8_t *pac, int l_query, uint8_t *query, int64_t rb, int64_t re, int *score, int *n_cigar, int *NM)
{
	return bwa_gen_cigar2(mat, q, r, q, r, w_, l_pac, pac, l_query, query, rb, re, score, n_cigar, NM);
}

/*********************
 * Full index reader *
 *********************/
// #if 0
char *bwa_idx_infer_prefix(const char *hint)
{
	char *prefix;
	int l_hint;
	FILE *fp;
	l_hint = strlen(hint);
	prefix = (char *) malloc(l_hint + 3 + 4 + 1);
	strcpy(prefix, hint);
	strcpy(prefix + l_hint, ".64.bwt");
	if ((fp = fopen(prefix, "rb")) != 0) {
		fclose(fp);
		prefix[l_hint + 3] = 0;
		return prefix;
	} else {
		strcpy(prefix + l_hint, ".bwt");
		if ((fp = fopen(prefix, "rb")) == 0) {
			free(prefix);
			return 0;
		} else {
			fclose(fp);
			prefix[l_hint] = 0;
			return prefix;
		}
	}
}

bwt_t *bwa_idx_load_bwt(const char *hint)
{
	char *tmp, *prefix;
	bwt_t *bwt;
	prefix = bwa_idx_infer_prefix(hint);
	if (prefix == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to locate the index files\n", __func__);
		return 0;
	}
	tmp = (char*) calloc(strlen(prefix) + 5, 1);
	strcat(strcpy(tmp, prefix), ".bwt"); // FM-index
	bwt = bwt_restore_bwt(tmp);
	strcat(strcpy(tmp, prefix), ".sa");  // partial suffix array (SA)
	bwt_restore_sa(tmp, bwt);
	free(tmp); free(prefix);
	return bwt;
}

bwaidx_t *bwa_idx_load_from_disk(const char *hint, int which)
{
	bwaidx_t *idx;
	char *prefix;
	prefix = bwa_idx_infer_prefix(hint);
	
	if (prefix == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to locate the index files\n", __func__);
		return 0;
	}
	idx = (bwaidx_t*) calloc(1, sizeof(bwaidx_t));
	if (which & BWA_IDX_BWT) idx->bwt = bwa_idx_load_bwt(hint);
	if (which & BWA_IDX_BNS) {
		int i, c;
		idx->bns = bns_restore(prefix);
		assert(idx->bns != 0);
		for (i = c = 0; i < idx->bns->n_seqs; ++i)
			if (idx->bns->anns[i].is_alt) ++c;
		if (bwa_verbose >= 3)
			fprintf(stderr, "[M::%s] read %d ALT contigs\n", __func__, c);
		if (which & BWA_IDX_PAC) {
			idx->pac = (uint8_t*) calloc(idx->bns->l_pac/4+1, 1);
			err_fread_noeof(idx->pac, 1, idx->bns->l_pac/4+1, idx->bns->fp_pac); // concatenated 2-bit encoded sequence
			err_fclose(idx->bns->fp_pac);
			idx->bns->fp_pac = 0;
		}
	}
	free(prefix);
	return idx;
}

bwaidx_t *bwa_idx_load(const char *hint, int which)
{
	return bwa_idx_load_from_disk(hint, which);
}

void bwa_idx_destroy(bwaidx_t *idx)
{
	if (idx == 0) return;
	if (idx->mem == 0) {
		if (idx->bwt) bwt_destroy(idx->bwt);
		if (idx->bns) bns_destroy(idx->bns);
		if (idx->pac) free(idx->pac);
	} else {
		free(idx->bwt); free(idx->bns->anns); free(idx->bns);
		if (!idx->is_shm) free(idx->mem);
	}
	free(idx);
}

int bwa_mem2idx(int64_t l_mem, uint8_t *mem, bwaidx_t *idx)
{
	int64_t k = 0, x;
	int i;

	// generate idx->bwt
	x = sizeof(bwt_t); idx->bwt = (bwt_t*) malloc(x); memcpy(idx->bwt, mem + k, x); k += x;
	x = idx->bwt->bwt_size * 4; idx->bwt->bwt = (uint32_t*)(mem + k); k += x;
	x = idx->bwt->n_sa * sizeof(bwtint_t); idx->bwt->sa = (bwtint_t*)(mem + k); k += x;

	// generate idx->bns and idx->pac
	x = sizeof(bntseq_t); idx->bns = (bntseq_t*) malloc(x); memcpy(idx->bns, mem + k, x); k += x;
	x = idx->bns->n_holes * sizeof(bntamb1_t); idx->bns->ambs = (bntamb1_t*)(mem + k); k += x;
	x = idx->bns->n_seqs  * sizeof(bntann1_t); idx->bns->anns = (bntann1_t*) malloc(x); memcpy(idx->bns->anns, mem + k, x); k += x;
	for (i = 0; i < idx->bns->n_seqs; ++i) {
		idx->bns->anns[i].name = (char*)(mem + k); k += strlen(idx->bns->anns[i].name) + 1;
		idx->bns->anns[i].anno = (char*)(mem + k); k += strlen(idx->bns->anns[i].anno) + 1;
	}
	idx->pac = (uint8_t*)(mem + k); k += idx->bns->l_pac/4+1;
	assert(k == l_mem);

	idx->l_mem = k; idx->mem = mem;
	return 0;
}

int bwa_idx2mem(bwaidx_t *idx)
{
	int i;
	int64_t k, x, tmp;
	uint8_t *mem;

	// copy idx->bwt
	x = idx->bwt->bwt_size * 4;
	mem = (uint8_t*) realloc(idx->bwt->bwt, sizeof(bwt_t) + x); idx->bwt->bwt = 0;
	memmove(mem + sizeof(bwt_t), mem, x);
	memcpy(mem, idx->bwt, sizeof(bwt_t)); k = sizeof(bwt_t) + x;
	x = idx->bwt->n_sa * sizeof(bwtint_t); mem = (uint8_t*) realloc(mem, k + x); memcpy(mem + k, idx->bwt->sa, x); k += x;
	free(idx->bwt->sa);
	free(idx->bwt); idx->bwt = 0;

	// copy idx->bns
	tmp = idx->bns->n_seqs * sizeof(bntann1_t) + idx->bns->n_holes * sizeof(bntamb1_t);
	for (i = 0; i < idx->bns->n_seqs; ++i) // compute the size of heap-allocated memory
		tmp += strlen(idx->bns->anns[i].name) + strlen(idx->bns->anns[i].anno) + 2;
	mem = (uint8_t*) realloc(mem, k + sizeof(bntseq_t) + tmp);
	x = sizeof(bntseq_t); memcpy(mem + k, idx->bns, x); k += x;
	x = idx->bns->n_holes * sizeof(bntamb1_t); memcpy(mem + k, idx->bns->ambs, x); k += x;
	free(idx->bns->ambs);
	x = idx->bns->n_seqs * sizeof(bntann1_t); memcpy(mem + k, idx->bns->anns, x); k += x;
	for (i = 0; i < idx->bns->n_seqs; ++i) {
		x = strlen(idx->bns->anns[i].name) + 1; memcpy(mem + k, idx->bns->anns[i].name, x); k += x;
		x = strlen(idx->bns->anns[i].anno) + 1; memcpy(mem + k, idx->bns->anns[i].anno, x); k += x;
		free(idx->bns->anns[i].name); free(idx->bns->anns[i].anno);
	}
	free(idx->bns->anns);

	// copy idx->pac
	x = idx->bns->l_pac/4+1;
	mem = (uint8_t*) realloc(mem, k + x);
	memcpy(mem + k, idx->pac, x); k += x;
	free(idx->bns); idx->bns = 0;
	free(idx->pac); idx->pac = 0;

	return bwa_mem2idx(k, mem, idx);
}
// #endif

/***********************
 * SAM header routines *
 ***********************/

void bwa_print_sam_hdr(const bntseq_t *bns, const char *hdr_line, FILE *fp)
{
	int i, n_SQ = 0;
	extern char *bwa_pg;
	if (hdr_line) {
		const char *p = hdr_line;
		while ((p = strstr(p, "@SQ\t")) != 0) {
			if (p == hdr_line || *(p-1) == '\n') ++n_SQ;
			p += 4;
		}
	}
	if (n_SQ == 0) {
		for (i = 0; i < bns->n_seqs; ++i) {
#if ORIG
			err_printf("@SQ\tSN:%s\tLN:%d", bns->anns[i].name, bns->anns[i].len);
			if (bns->anns[i].is_alt) err_printf("\tAH:*\n");
			else err_fputc('\n', stdout);
#else
			char buf[500];
			sprintf(buf, "@SQ\tSN:%s\tLN:%d", bns->anns[i].name, bns->anns[i].len);
			err_fputs(buf, fp);
			if (bns->anns[i].is_alt) {
				sprintf(buf, "\tAH:*\n");
				err_fputs(buf, fp);
			} else
				err_fputc('\n', fp);
#endif
		}
	} else if (n_SQ != bns->n_seqs && bwa_verbose >= 2)
		fprintf(stderr, "[W::%s] %d @SQ lines provided with -H; %d sequences in the index. "
			   "Continue anyway.\n", __func__, n_SQ, bns->n_seqs);
	
#if ORIG
	if (hdr_line) err_printf("%s\n", hdr_line);
	if (bwa_pg) err_printf("%s\n", bwa_pg);
#else
	if (hdr_line) {
		err_fputs(hdr_line, fp);
		err_fputs("\n", fp);
	}
	if (bwa_pg) err_fputs(bwa_pg, fp);
#endif
}

static char *bwa_escape(char *s)
{
	char *p, *q;
	for (p = q = s; *p; ++p) {
		if (*p == '\\') {
			++p;
			if (*p == 't') *q++ = '\t';
			else if (*p == 'n') *q++ = '\n';
			else if (*p == 'r') *q++ = '\r';
			else if (*p == '\\') *q++ = '\\';
		} else *q++ = *p;
	}
	*q = '\0';
	return s;
}

char *bwa_set_rg(const char *s)
{
	char *p, *q, *r, *rg_line = 0;
	memset(bwa_rg_id, 0, 256);
	if (strstr(s, "@RG") != s) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] the read group line is not started with @RG\n", __func__);
		goto err_set_rg;
	}
	rg_line = strdup(s);
	bwa_escape(rg_line);
	if ((p = strstr(rg_line, "\tID:")) == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] no ID at the read group line\n", __func__);
		goto err_set_rg;
	}
	p += 4;
	for (q = p; *q && *q != '\t' && *q != '\n'; ++q);
	if (q - p + 1 > 256) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] @RG:ID is longer than 255 characters\n", __func__);
		goto err_set_rg;
	}
	for (q = p, r = bwa_rg_id; *q && *q != '\t' && *q != '\n'; ++q)
		*r++ = *q;
	return rg_line;

err_set_rg:
	free(rg_line);
	return 0;
}

char *bwa_insert_header(const char *s, char *hdr)
{
	int len = 0;
	if (s == 0 || s[0] != '@') return hdr;
	if (hdr) {
		len = strlen(hdr);
		hdr = (char*) realloc(hdr, len + strlen(s) + 2);
		hdr[len++] = '\n';
		strcpy(hdr + len, s);
	} else hdr = strdup(s);
	bwa_escape(hdr + len);
	return hdr;
}
