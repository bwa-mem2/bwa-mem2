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

Authors: Saurabh Kalikar <saurabh.kalikar@intel.com>; Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>; Chirag Jain
*****************************************************************************************/

#include <immintrin.h>
#include <assert.h>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <climits>
#include <vector>
#include <map>
using namespace std;


static const char LogTable256_dp_lib[256] = {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
	-1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
	LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
	LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
};

static inline int ilog2_32_dp_lib(uint32_t v)
{
	uint32_t t, tt;
	if ((tt = v>>16)) return (t = tt>>8) ? 24 + LogTable256_dp_lib[t] : 16 + LogTable256_dp_lib[tt];
	return (t = v>>8) ? 8 + LogTable256_dp_lib[t] : LogTable256_dp_lib[v];
}


class anchor_t {
	public:
		uint64_t r;
		int32_t q;
		int32_t l;
		anchor_t(){}

		anchor_t(uint64_t x, int32_t y, int32_t length){
			r = x; q = y; l = length;
		}	
};

typedef uint32_t num_bits_t;
void create_SoA_Anchors_32_bit(anchor_t* anc, num_bits_t anc_size, num_bits_t* &anchor_r, num_bits_t* &anchor_q, num_bits_t* &anchor_l){
	//num_bits_t anc_size = anc.size();

	anchor_r = (num_bits_t*) malloc((16+anc_size)*sizeof(num_bits_t));
	anchor_q = (num_bits_t*) malloc((16+anc_size)*sizeof(num_bits_t));
	anchor_l = (num_bits_t*) malloc((16+anc_size)*sizeof(num_bits_t));

	anchor_r = &anchor_r[16];
	anchor_q = &anchor_q[16];
	anchor_l = &anchor_l[16];

	for(uint32_t i = 0; i < anc_size; i++){
		anchor_r[i] = (num_bits_t)anc[i].r;
		anchor_q[i] = (num_bits_t)anc[i].q;
		anchor_l[i] = (num_bits_t)anc[i].l;
	}


}
class dp_chain {
	public:
	//Tunable parameters 
	int max_dist_x, max_dist_y, bw, max_skip, max_iter, is_cdna, n_segs;
	float gap_scale;
#ifdef __AVX512BW__
	__m512i zero_v;// = _mm512_setzero_si512();
#elif __AVX2__
	__m256i zero_avx2_v;// = _mm512_setzero_si512();
#endif

	dp_chain(){}

	void test(){
		printf("hyper-parameters: %d %d %d %d %d %f %d %d\n",max_dist_x, max_dist_y, bw, max_skip, max_iter, gap_scale, is_cdna, n_segs);
	}

	dp_chain(int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter, float gap_scale, int is_cdna, int n_segs){
 		this->max_dist_x = max_dist_x;
		this->max_dist_y = max_dist_y;
		this->bw = bw;
		this->max_skip = max_skip;
		this->max_iter = max_iter;
		this->gap_scale = gap_scale;
		this->is_cdna = is_cdna;
		this->n_segs = n_segs;
#ifdef __AVX512BW__
		zero_v = _mm512_setzero_si512();
#elif __AVX2__
		zero_avx2_v = _mm256_setzero_si256();
#endif

	}
	int32_t get_overlap_cost(anchor_t b, anchor_t a){
		int32_t ref_overlap = a.r + a.l - b.r;
		int32_t query_overlap = a.q + a.l - b.q;

		return max(0, max(ref_overlap, query_overlap));
	}


#ifdef __AVX512BW__
	inline __m512i get_gap_cost_vectorized_int32(__m512i dd_v, float avg_qspan, float gap_scale){
		//Vectorized log2
		uint32_t base = 31;
		__m512i vbase = _mm512_set1_epi32(base);
		__mmask16 msk = 0xFFFF;
        	__m512i vout  = _mm512_mask_lzcnt_epi32(dd_v, msk, dd_v);
        	//int res = base - temp[0];
    		__m512i r_v = _mm512_sub_epi32(vbase, vout);

		// log_dd = dd?ilog2:0; log_dd>>1
		__mmask16 neg_mask = _mm512_cmpneq_epi32_mask(dd_v, zero_v); 
		__m512i log_dd_v = _mm512_srli_epi32(_mm512_maskz_or_epi32(neg_mask, r_v, zero_v), 1);

		//dd * 0.01*avg_qspan 
		float avg_qspan_val = 0.01*avg_qspan;
		__m512 avg_qspan_v = _mm512_set1_ps(avg_qspan_val);
		__m512i cost = _mm512_cvt_roundps_epi32(_mm512_mul_ps(_mm512_cvtepi32_ps(dd_v), avg_qspan_v),_MM_FROUND_TO_NEG_INF |_MM_FROUND_NO_EXC);

		// gap_cost = (dd * 0.01*avg_qspan) + (log_dd>>1 + 0.499)
		__m512i gap_cost_v = (_mm512_add_epi32(cost, log_dd_v));

		return gap_cost_v;
	}
#elif __AVX2__
//#endif
	inline __m256i get_gap_cost_vectorized_int32(__m256i dd_v, float avg_qspan, float gap_scale){

		/*		
		//Vectorized log2
		uint32_t base = 31;
		__m256i vbase = _mm256_set1_epi32(base);
		__mmask8 msk = 0xFF;
        	__m256i vout  = _mm256_mask_lzcnt_epi32(dd_v, msk, dd_v);
        	//int res = base - temp[0];
    		__m256i r_v = _mm256_sub_epi32(vbase, vout);
		*/
		uint32_t lg[8], dd_array[8];
		_mm256_storeu_si256((__m256i *)dd_array, dd_v);
		for(int i = 0; i < 8; i++){
			lg[i] = ilog2_32_dp_lib(dd_array[i]);
		}
		__m256i r_v = _mm256_loadu_si256((__m256i*)lg);
		
	

		// log_dd = dd?ilog2:0; log_dd>>1
		//__mmask8 neg_mask = _mm256_cmpneq_epi32_mask(dd_v, zero_avx2_v); 
		//__m256i log_dd_v = _mm256_srli_epi32(_mm256_maskz_or_epi32(neg_mask, r_v, zero_avx2_v), 1);

		__m256i neg_mask = _mm256_cmpeq_epi32(dd_v, zero_avx2_v);
		__m256i log_dd_v = _mm256_srli_epi32(_mm256_andnot_si256 (neg_mask, r_v), 1);

		//dd * 0.01*avg_qspan 
		float avg_qspan_val = 0.01*avg_qspan;
		__m256 avg_qspan_v = _mm256_set1_ps(avg_qspan_val);
		__m256i cost = _mm256_cvtps_epi32(_mm256_round_ps((_mm256_mul_ps(_mm256_cvtepi32_ps(dd_v), avg_qspan_v)),_MM_FROUND_TO_NEG_INF |_MM_FROUND_NO_EXC ));

		// gap_cost = (dd * 0.01*avg_qspan) + (log_dd>>1 + 0.499)
		__m256i gap_cost_v = (_mm256_add_epi32(cost, log_dd_v));

		return gap_cost_v;
	}
#endif


	
	int32_t get_gap_cost(anchor_t a, anchor_t b, float avg_qspan, float gap_scale){
		int32_t dr = a.r - b.r;
		int32_t dq = a.q - b.q;
		//bool flag = true;
		//bool is_cdna = false;	
		int32_t dd = dr > dq? dr - dq : dq - dr; //dd = |dr-dq|;
		int32_t log_dd = dd? ilog2_32_dp_lib(dd) : 0;
		int32_t gap_cost = 0;
		
		gap_cost = (int)(dd * .01 * avg_qspan) + (log_dd>>1) + 0.00; //TODO: Only multiplication should be casted to (int)

		return gap_cost;
	}

	int32_t get_gap_cost(anchor_t a, anchor_t b){
		return 0;
	}


#ifdef __AVX512BW__
	//Vector code with SoA function parameters 32-bit number representation - avx512
	void mm_dp_vectorized(int64_t n, anchor_t *anchors, uint32_t* anchor_r, uint32_t* anchor_q, uint32_t* anchor_l, uint32_t* &f, int32_t* &p, int32_t* &v, int32_t dr, int32_t dq, int (*gap_cost)(anchor_t a, anchor_t b, void* meta_data), void* meta_data)
	{ 
		uint64_t sum_qspan = 0;	
		float avg_qspan;
		for (int i = 0; i < n; ++i) sum_qspan += anchor_l[i];
		avg_qspan = (float)sum_qspan / n;
		int st = 0;	
		
		__m512i dr_v = _mm512_set1_epi32((int64_t)dr);
		__m512i dq_v = _mm512_set1_epi32((int64_t)dq);
		__m512i j_idx_base = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);		
		int32_t maxfVector_v[16];
		int32_t maxjVector_v[16];
		__m512i neg_one_v = _mm512_set1_epi32((int32_t)-1);

		for(int i = 0; i < n; i++){


			int32_t max_j = -1;
			int32_t max_f = anchor_l[i];

			

			//uint64_t ri = anchors[i].r;
			while (st < i && !(anchors[i].r - anchors[st].r <= (uint32_t)dr)) ++st; //predecessor's position is too far

			// TODO: Minimap specific max_iter parameter
			if (i - st > max_iter) st = i - max_iter; //predecessor's index is too far

			
			int j = i - 1;
			if(!(j - st <= 5))
			{
				//broadcast ri and qi
				__m512i ri_v = _mm512_set1_epi32(anchor_r[i]); 
				__m512i qi_v = _mm512_set1_epi32(anchor_q[i]);
			

				__m512i maxj_v = neg_one_v;
				__m512i maxf_v = _mm512_set1_epi32((int32_t)anchor_l[i]);
				__m512i li_v = maxf_v;

			// 16-way vectorized
			//_mm_prefetch((const char *)(&anchor_r[j - 30]), _MM_HINT_T2);	
			//_mm_prefetch((const char *)(&anchor_q[j - 30]), _MM_HINT_T2);	

			for(j = i - 1; (j - 15) >= st; j = j - 16){
		
				_mm_prefetch((const char *)(&anchor_r[j-60]), _MM_HINT_T0);	
				_mm_prefetch((const char *)(&anchor_q[j-60]), _MM_HINT_T0);	
	
				uint32_t *rj, *qj; 
				rj = &anchor_r[j-15];	qj = &anchor_q[j-15];

 
				// Load rj and qj
    				__m512i rj_v = _mm512_loadu_si512(rj);
    				__m512i qj_v = _mm512_loadu_si512(qj);
			 	
					
    				__m512i ddr_v = _mm512_sub_epi32(ri_v, rj_v);
    				__m512i ddq_v = _mm512_sub_epi32(qi_v, qj_v);
			
				//TODO: Minimap2 specific continue condition
				__m512i dd_v = _mm512_abs_epi32(_mm512_sub_epi32(ddr_v, ddq_v)); 
				__m512i bw_v = _mm512_set1_epi32((int32_t) bw);
				__mmask16 bw_gt = _mm512_cmpgt_epi32_mask(dd_v, bw_v);						
				__mmask16 mask_eq = _mm512_cmpeq_epi32_mask(ddr_v, zero_v);						
				__mmask16 mask_leq = _mm512_cmple_epi32_mask(ddq_v, zero_v);						
				__mmask16 mask_gt1 = _mm512_cmpgt_epi32_mask(ddq_v, dq_v);						
				__mmask16 mask_gt2 = _mm512_cmpgt_epi32_mask(ddq_v, dr_v);						
				
				__mmask16 loopContinueMask = ~(bw_gt | mask_eq | mask_leq | mask_gt1 | mask_gt2);

				// Load f[j-8, j]
    				__m512i fj_v = _mm512_loadu_si512(&f[j-15]);
				
				//Vectorized gap cost function
				__m512i gc_v = get_gap_cost_vectorized_int32(dd_v, avg_qspan, gap_scale);

		//---------------- Inline get_overlap_cost function -------------------
				__m512i min1 = _mm512_min_epi32(ddr_v, ddq_v);	
				__m512i oc_v = _mm512_min_epi32(li_v, min1);	
		//----------------------------------------------------------------------		
				__m512i f_plus_oc_v = _mm512_add_epi32(fj_v, oc_v);
				__m512i sc_v = _mm512_maskz_sub_epi32(loopContinueMask,f_plus_oc_v, gc_v);

				// Update Maxf and Maxj
				__mmask16 mask_max_sc = _mm512_cmpgt_epi32_mask(sc_v, maxf_v);
				__m512i j_idx_v = _mm512_add_epi32(j_idx_base, _mm512_set1_epi32(j - 15));

				maxf_v = _mm512_max_epi32(sc_v, maxf_v);
				maxj_v = _mm512_mask_blend_epi32(mask_max_sc, maxj_v, j_idx_v);
			}
			

			if(j >= st)
			{
				uint32_t *rj, *qj; 
				rj = &anchor_r[j-15];	qj = &anchor_q[j-15];

 
				// Load rj and qj
    				__m512i rj_v = _mm512_loadu_si512(rj);
    				__m512i qj_v = _mm512_loadu_si512(qj);
			 	
					
    				__m512i ddr_v = _mm512_sub_epi32(ri_v, rj_v);
    				__m512i ddq_v = _mm512_sub_epi32(qi_v, qj_v);
			
				//TODO: Minimap2 specific continue condition
				__m512i dd_v = _mm512_abs_epi32(_mm512_sub_epi32(ddr_v, ddq_v)); 
				__m512i bw_v = _mm512_set1_epi32((int32_t) bw);
				__mmask16 bw_gt = _mm512_cmpgt_epi32_mask(dd_v, bw_v);						
				__mmask16 mask_eq = _mm512_cmpeq_epi32_mask(ddr_v, zero_v);						
				__mmask16 mask_leq = _mm512_cmple_epi32_mask(ddq_v, zero_v);						
				__mmask16 mask_gt1 = _mm512_cmpgt_epi32_mask(ddq_v, dq_v);						
				__mmask16 mask_gt2 = _mm512_cmpgt_epi32_mask(ddq_v, dr_v);						

				
				__mmask16 loopContinueMask = ~(bw_gt | mask_eq | mask_leq | mask_gt1 | mask_gt2);


				//Last partial vector processing mask - To enable, change loop condition to  -> j >= st
				
  				int shift = st - (j-15);
					
				loopContinueMask = loopContinueMask>>(shift);
				loopContinueMask = loopContinueMask<<(shift);
				if(loopContinueMask != 0x00){



				// Load f[j-8, j]
    				__m512i fj_v = _mm512_loadu_si512(&f[j-15]);
				
				//Vectorized gap cost function
				__m512i gc_v = get_gap_cost_vectorized_int32(dd_v, avg_qspan, gap_scale);

		//---------------- Inline get_overlap_cost function -------------------
				__m512i min1 = _mm512_min_epi32(ddr_v, ddq_v);	
				__m512i oc_v = _mm512_min_epi32(li_v, min1);	
		//----------------------------------------------------------------------		
				__m512i f_plus_oc_v = _mm512_add_epi32(fj_v, oc_v);
				__m512i sc_v = _mm512_maskz_sub_epi32(loopContinueMask,f_plus_oc_v, gc_v);


				// Update Maxf and Maxj
				__mmask16 mask_max_sc = _mm512_cmpgt_epi32_mask(sc_v, maxf_v);
				__m512i j_idx_v = _mm512_add_epi32(j_idx_base, _mm512_set1_epi32(j - 15));

				maxf_v = _mm512_max_epi32(sc_v, maxf_v);
				maxj_v = _mm512_mask_blend_epi32(mask_max_sc, maxj_v, j_idx_v);
		
				}
			}


				_mm512_store_epi32(maxfVector_v, maxf_v);
				_mm512_store_epi32(maxjVector_v, maxj_v);
			
				for(int iter = 15; iter >=0; iter--){
					if(maxfVector_v[iter] > max_f) {
						max_f = maxfVector_v[iter];
						max_j = maxjVector_v[iter];
					}
					else if (maxfVector_v[iter] == max_f){
						max_j = max (max_j, maxjVector_v[iter]);
						if((uint32_t)max_f == anchor_l[i]) max_j = -1;
					}
				}
		
			
			}
			else{
			int32_t ri = anchor_r[i], qi = anchor_q[i];
			for(; j >=st; j--){

				int32_t ddr, ddq;

				int32_t rj = anchor_r[j], 
					qj = anchor_q[j]; 
				ddr = ri - rj;
				ddq = qi - qj;

				if(abs(ddr - ddq) > bw) continue;
				if(ddr == 0 || ddq <= 0) continue;

				if(ddq > dq || ddq > dr) continue;

				
				int32_t oc = 0;

				int32_t score = f[j];//anchor_l[i]; 
				oc = ddr < ddq? ddr: ddq;
				oc = oc < (int32_t)anchor_l[i]? oc : anchor_l[i];
				score += oc;

				int32_t dr = ddr;
				int32_t dq = ddq;
				int32_t dd = abs(dr - dq);//dr > dq? dr - dq : dq - dr; //dd = |dr-dq|;
				int32_t log_dd = dd? ilog2_32_dp_lib(dd) : 0;
				int32_t gap_cost = 0;
				
				gap_cost = (int)(dd * .01 * avg_qspan) + (log_dd>>1) + 0.00; //TODO: Only multiplication should be casted to (int)
				score -= gap_cost;

				if(score > max_f){
					max_f = score;
					max_j = j;
				}

			}

			}

			f[i] = max_f; 
			p[i] = max_j;
			v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f; // v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
		}

	}
#elif __AVX2__
	void mm_dp_vectorized(int64_t n, anchor_t *anchors, uint32_t* anchor_r, uint32_t* anchor_q, uint32_t* anchor_l, uint32_t* &f, int32_t* &p, int32_t* &v, int32_t dr, int32_t dq, int (*gap_cost)(anchor_t a, anchor_t b, void* meta_data), void* meta_data)
	{ 
		uint64_t sum_qspan = 0;	
		float avg_qspan;
		for (int i = 0; i < n; ++i) sum_qspan += anchor_l[i];
		avg_qspan = (float)sum_qspan / n;
		int st = 0;	
		
		__m256i dr_v = _mm256_set1_epi32(dr);
		__m256i dq_v = _mm256_set1_epi32(dq);
		__m256i j_idx_base = _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);		
		int32_t maxfVector_v[8];
		int32_t maxjVector_v[8];
		__m256i neg_one_v = _mm256_set1_epi32((int32_t)-1);

		for(int i = 0; i < n; i++){


			int32_t max_j = -1;
			int32_t max_f = anchor_l[i];

			

			//uint64_t ri = anchors[i].r;
			while (st < i && !(anchors[i].r - anchors[st].r <= dr)) ++st; //predecessor's position is too far

			// TODO: Minimap specific max_iter parameter
			if (i - st > max_iter) st = i - max_iter; //predecessor's index is too far

			
			int j = i - 1;
			if(!(j - st <= 5))
			{
				//broadcast ri and qi
				__m256i ri_v = _mm256_set1_epi32(anchor_r[i]); 
				__m256i qi_v = _mm256_set1_epi32(anchor_q[i]);
			

				__m256i maxj_v = neg_one_v;
				__m256i maxf_v = _mm256_set1_epi32((int32_t)anchor_l[i]);
				__m256i li_v = maxf_v;

			// 8-way vectorized
			for(j = i - 1; (j - 7) >= st; j = j - 8){
		
				_mm_prefetch((const char *)(&anchor_r[j-60]), _MM_HINT_T0);	
				_mm_prefetch((const char *)(&anchor_q[j-60]), _MM_HINT_T0);	
	
				uint32_t *rj, *qj; 
				int j_stride = j - 7;
				rj = &anchor_r[j_stride];	qj = &anchor_q[j_stride];

 
				// Load rj and qj
    				__m256i rj_v = _mm256_loadu_si256((__m256i*) rj);
    				__m256i qj_v = _mm256_loadu_si256((__m256i*) qj);
			 	
					
    				__m256i ddr_v = _mm256_sub_epi32(ri_v, rj_v);
    				__m256i ddq_v = _mm256_sub_epi32(qi_v, qj_v);
			
				//TODO: Minimap2 specific continue condition
				__m256i dd_v = _mm256_abs_epi32(_mm256_sub_epi32(ddr_v, ddq_v)); 
				__m256i bw_v = _mm256_set1_epi32((int32_t) bw);
				
			/*
				__mmask8 bw_gt = _mm256_cmpgt_epi32_mask(dd_v, bw_v);						
				__mmask8 mask_eq = _mm256_cmpeq_epi32_mask(ddr_v, zero_avx2_v);						
				__mmask8 mask_leq1 = _mm256_cmpgt_epi32_mask(zero_avx2_v, ddq_v);						
				__mmask8 mask_leq2 = _mm256_cmpeq_epi32_mask(ddq_v, zero_avx2_v);						
				__mmask8 mask_gt1 = _mm256_cmpgt_epi32_mask(ddq_v, dq_v);						
				__mmask8 mask_gt2 = _mm256_cmpgt_epi32_mask(ddq_v, dr_v);						

				
				__mmask8 loopContinueMask = ~(bw_gt | mask_eq | (mask_leq1 | mask_leq2 ) | mask_gt1 | mask_gt2);
			*/	

				__m256i bw_gt = _mm256_cmpgt_epi32(dd_v, bw_v);						
				__m256i mask_eq = _mm256_cmpeq_epi32(ddr_v, zero_avx2_v);						
				__m256i mask_leq1 = _mm256_cmpgt_epi32(zero_avx2_v, ddq_v);						
				__m256i mask_leq2 = _mm256_cmpeq_epi32(ddq_v, zero_avx2_v);						
				__m256i mask_gt1 = _mm256_cmpgt_epi32(ddq_v, dq_v);						
				__m256i mask_gt2 = _mm256_cmpgt_epi32(ddq_v, dr_v);						
				
				__m256i tmp1 = _mm256_or_si256(bw_gt, mask_eq); 
				__m256i tmp2 = _mm256_or_si256(mask_leq1, mask_leq2); 
				__m256i tmp3 = _mm256_or_si256(mask_gt1, mask_gt2); 

				__m256i loopContinueMask = _mm256_or_si256(_mm256_or_si256 (tmp1, tmp2), tmp3);	

				// Load f[j-8, j]
    				__m256i fj_v = _mm256_loadu_si256((__m256i*) &f[j_stride]);
				
				//Vectorized gap cost function
				__m256i gc_v = get_gap_cost_vectorized_int32(dd_v, avg_qspan, gap_scale);

		//---------------- Inline get_overlap_cost function -------------------
				__m256i min1 = _mm256_min_epi32(ddr_v, ddq_v);	
				__m256i oc_v = _mm256_min_epi32(li_v, min1);	
		//----------------------------------------------------------------------		
				__m256i f_plus_oc_v = _mm256_add_epi32(fj_v, oc_v);
				//__m256i sc_v = _mm256_maskz_sub_epi32(loopContinueMask,f_plus_oc_v, gc_v);
				__m256i sc_v = _mm256_andnot_si256(loopContinueMask, _mm256_sub_epi32(f_plus_oc_v, gc_v));

				// Update Maxf and Maxj
				__m256i mask_max_sc = _mm256_cmpgt_epi32(sc_v, maxf_v);
				__m256i j_idx_v = _mm256_add_epi32(j_idx_base, _mm256_set1_epi32(j_stride));

				maxf_v = _mm256_max_epi32(sc_v, maxf_v);
				maxj_v = _mm256_or_si256(_mm256_andnot_si256(mask_max_sc, maxj_v), _mm256_and_si256(mask_max_sc, j_idx_v));
			}
			

			if(j >= st)
			{
				uint32_t *rj, *qj; 
				int j_stride = j - 7;
				rj = &anchor_r[j_stride];	qj = &anchor_q[j_stride];

 
				// Load rj and qj
    				__m256i rj_v = _mm256_loadu_si256((__m256i*) rj);
    				__m256i qj_v = _mm256_loadu_si256((__m256i*) qj);
			 	
					
    				__m256i ddr_v = _mm256_sub_epi32(ri_v, rj_v);
    				__m256i ddq_v = _mm256_sub_epi32(qi_v, qj_v);
			
				//TODO: Minimap2 specific continue condition
				__m256i dd_v = _mm256_abs_epi32(_mm256_sub_epi32(ddr_v, ddq_v)); 
				__m256i bw_v = _mm256_set1_epi32((int32_t) bw);
/*
				__mmask8 bw_gt = _mm256_cmpgt_epi32_mask(dd_v, bw_v);						
				__mmask8 mask_eq = _mm256_cmpeq_epi32_mask(ddr_v, zero_avx2_v);						
				__mmask8 mask_leq = _mm256_cmple_epi32_mask(ddq_v, zero_avx2_v);						
				__mmask8 mask_gt1 = _mm256_cmpgt_epi32_mask(ddq_v, dq_v);						
				__mmask8 mask_gt2 = _mm256_cmpgt_epi32_mask(ddq_v, dr_v);						

				
				__mmask8 loopContinueMask = ~(bw_gt | mask_eq | mask_leq | mask_gt1 | mask_gt2);
*/
				__m256i bw_gt = _mm256_cmpgt_epi32(dd_v, bw_v);						
				__m256i mask_eq = _mm256_cmpeq_epi32(ddr_v, zero_avx2_v);						
				__m256i mask_leq1 = _mm256_cmpgt_epi32(zero_avx2_v, ddq_v);						
				__m256i mask_leq2 = _mm256_cmpeq_epi32(ddq_v, zero_avx2_v);						
				__m256i mask_gt1 = _mm256_cmpgt_epi32(ddq_v, dq_v);						
				__m256i mask_gt2 = _mm256_cmpgt_epi32(ddq_v, dr_v);						
				
				__m256i tmp1 = _mm256_or_si256(bw_gt, mask_eq); 
				__m256i tmp2 = _mm256_or_si256(mask_leq1, mask_leq2); 
				__m256i tmp3 = _mm256_or_si256(mask_gt1, mask_gt2); 

				__m256i loopContinueMask = _mm256_or_si256(_mm256_or_si256 (tmp1, tmp2), tmp3);	

				//Last partial vector processing mask - To enable, change loop condition to  -> j >= st
				
  				int shift = st - (j_stride);
				
				int32_t msk_ar[8];
				for(int it = 0; it < 8; it++){
					msk_ar[it] = (it < (shift))?0xFFFFFFFF:0;
				}
				loopContinueMask = _mm256_or_si256(loopContinueMask, _mm256_loadu_si256((__m256i*)msk_ar));
				//loopContinueMask = loopContinueMask>>(shift);
				//loopContinueMask = loopContinueMask<<(shift);
				//if(loopContinueMask != 0x0)
				{



				// Load f[j-8, j]
    				__m256i fj_v = _mm256_loadu_si256((__m256i*) &f[j_stride]);
				
				//Vectorized gap cost function
				__m256i gc_v = get_gap_cost_vectorized_int32(dd_v, avg_qspan, gap_scale);

		//---------------- Inline get_overlap_cost function -------------------
				__m256i min1 = _mm256_min_epi32(ddr_v, ddq_v);	
				__m256i oc_v = _mm256_min_epi32(li_v, min1);	
		//----------------------------------------------------------------------		
				__m256i f_plus_oc_v = _mm256_add_epi32(fj_v, oc_v);
				//__m256i sc_v = _mm256_maskz_sub_epi32(loopContinueMask,f_plus_oc_v, gc_v);
				__m256i sc_v = _mm256_andnot_si256(loopContinueMask, _mm256_sub_epi32(f_plus_oc_v, gc_v));


				// Update Maxf and Maxj
				__m256i mask_max_sc = _mm256_cmpgt_epi32(sc_v, maxf_v);
				__m256i j_idx_v = _mm256_add_epi32(j_idx_base, _mm256_set1_epi32(j_stride));

				maxf_v = _mm256_max_epi32(sc_v, maxf_v);
				maxj_v = _mm256_or_si256(_mm256_andnot_si256(mask_max_sc, maxj_v), _mm256_and_si256(mask_max_sc, j_idx_v));
		
				}
			}


				//_mm256_store_epi32(maxfVector_v, maxf_v);
				//_mm256_store_epi32(maxjVector_v, maxj_v);
				_mm256_store_si256((__m256i*) maxfVector_v, maxf_v);
				_mm256_store_si256((__m256i*) maxjVector_v, maxj_v);
			
				for(int iter = 7; iter >=0; iter--){
					if(maxfVector_v[iter] > max_f) {
						max_f = maxfVector_v[iter];
						max_j = maxjVector_v[iter];
					}
					else if (maxfVector_v[iter] == max_f){
						max_j = max (max_j, maxjVector_v[iter]);
						if(max_f == anchor_l[i]) max_j = -1;
					}
				}
		
			
			}
			else{
			int32_t ri = anchor_r[i], qi = anchor_q[i];
			for(; j >=st; j--){

				int32_t ddr, ddq;

				int32_t rj = anchor_r[j], 
					qj = anchor_q[j]; 
				ddr = ri - rj;
				ddq = qi - qj;

				if(abs(ddr - ddq) > bw) continue;
				if(ddr == 0 || ddq <= 0) continue;

				if(ddq > dq || ddq > dr) continue;

				
				int32_t oc = 0;
				int32_t lj = anchor_l[j];
				int32_t ref_overlap = rj + lj - ri;
				int32_t query_overlap = qj + lj - qi;

				int32_t score = f[j];//anchor_l[i]; 
				oc = ddr < ddq? ddr: ddq;
				oc = oc < anchor_l[i]? oc : anchor_l[i];
				score += oc;

				int32_t dr = ddr;
				int32_t dq = ddq;
				int32_t dd = abs(dr - dq);//dr > dq? dr - dq : dq - dr; //dd = |dr-dq|;
				int32_t log_dd = dd? ilog2_32_dp_lib(dd) : 0;
				int32_t gap_cost = 0;
				
				gap_cost = (int)(dd * .01 * avg_qspan) + (log_dd>>1) + 0.00; //TODO: Only multiplication should be casted to (int)
				score -= gap_cost;

				bool check = score > max_f;
				if(score > max_f){
					max_f = score;
					max_j = j;
				}

			}

			}

			f[i] = max_f; 
			p[i] = max_j;
			v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f; // v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
		}

	}
#endif

};
