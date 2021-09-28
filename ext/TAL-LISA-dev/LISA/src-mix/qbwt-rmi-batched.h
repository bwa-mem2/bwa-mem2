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

Authors: Saurabh Kalikar <saurabh.kalikar@intel.com>; Sanchit Misra <sanchit.misra@intel.com>
*****************************************************************************************/
#ifndef QBWT_RMI_BATCHED_H
#define QBWT_RMI_BATCHED_H
#ifdef VTUNE_ANALYSIS
#include <ittnotify.h>
#endif

#include "qbwt-ipbwt-rmi.h"
#include "chunkEncode.h"
class threadData {

	public:

		Info *chunk_pool;
		int chunk_cnt ;
		uint64_t *str_enc;
		int64_t *intv_all;
		
    		int64_t numSMEMs;
		Info *fmi_pool;
		int fmi_cnt;

		Info *tree_pool;
		int tree_cnt;
		index_t **s_siz;
		QBWT_HYBRID<index_t>::LcpInfo **s_info;
		uint8_t *s_msk;
		threadData(int64_t pool_size);
		void dealloc_td();
};

#define S_SWP_END do{\
        int c; \
	if(cnt > pref_dist)\
		c = --cnt;\
	else\
		c = --shrink_batch_size;\
        q = tree_pool[c]; \
        siz[0] = s_siz[c][0]; siz[1] = s_siz[c][1]; \
        info[0] = s_info[c][0]; info[1] = s_info[c][1]; \
        msk = s_msk[c]; \
}while(0)

#define S_RUN do{\
    bool small = siz[1] < siz[0]; \
    if(info[small].bw_ext_msk & msk) { \
        if(small) { \
            q.r = q.l + qbwt.get_lcp(q.intv.second); \
            q.intv.second = q.intv.first + siz[1]; \
        } else { \
            q.r = q.l + qbwt.get_lcp(q.intv.first); \
            q.intv.first = q.intv.second - siz[0]; \
        } \
	if(q.r > q.mid)\
		fmi_pool[fmi_cnt++] = q;\
        S_SWP_END;\
    } else { \
        if(small) { \
            q.intv.second = q.intv.first + siz[1]; \
            info[1] = qbwt.lcpi[q.intv.second]; \
            siz[1] = GET_WIDTH(info[1].s_width, q.intv.second); \
        } else { \
            q.intv.first = q.intv.second - siz[0]; \
            info[0] = qbwt.lcpi[q.intv.first]; \
            siz[0] = GET_WIDTH(info[0].s_width, q.intv.first); \
        } \
    }\
}while(0)


#define S_PREFETCH do{\
            bool small = siz[1] < siz[0]; \
                my_prefetch((const char*)(qbwt.lcpp1 + (small ? q.intv.second : q.intv.first)), _MM_HINT_T0); \
                my_prefetch((const char*)(qbwt.lcpi + (small ?  \
                                q.intv.first + siz[1] : q.intv.second - siz[0])), _MM_HINT_T0); \
}while(0)

#define GET_WIDTH(A, B) (A == qbwt.WID_MAX ? \
        (index_t)lower_bound(qbwt.b_width.begin(), qbwt.b_width.end(), B,  \
            [&](pair<index_t, index_t> p, index_t qq){return p.first < qq;})->second: \
            (index_t)A)


#define S_LOAD(i) \
            Info &q = tree_pool[i]; \
            index_t* siz = s_siz[i]; \
            QBWT_HYBRID<index_t>::LcpInfo* info = s_info[i]; \
            uint8_t &msk = s_msk[i] 




inline void s_pb(QBWT_HYBRID<index_t> &qbwt, Info &_q, int cnt, threadData &td);

void tree_shrink_batched( QBWT_HYBRID<index_t> &qbwt, int cnt, threadData &td);

void fmi_extend_batched( QBWT_HYBRID<index_t> &qbwt, int cnt, Info* q_batch, threadData &td, Output* output, int min_seed_len);

void fmi_extend_batched_exact_search( QBWT_HYBRID<index_t> &qbwt, int cnt, Info* q_batch, threadData &td, Output* output, int min_seed_len);

void fmi_shrink_batched( QBWT_HYBRID<index_t> &qbwt, int cnt, Info* q_batch, threadData &td, Info* output, int min_seed_len);

void smem_rmi_batched(Info *qs, int64_t qs_size, int64_t batch_size, QBWT_HYBRID<index_t> &qbwt, threadData &td, Output* output, int min_seed_len, bool apply_lisa = true);

void exact_search_rmi_batched(Info *qs, int64_t qs_size, int64_t batch_size, QBWT_HYBRID<index_t> &qbwt, threadData &td, Output* output, int min_seed_len, bool apply_lisa = true);

void exact_search_rmi_batched_k3(Info *qs, int64_t qs_size, int64_t batch_size, QBWT_HYBRID<index_t> &qbwt, threadData &td, Output* output, int min_seed_len, FMI_search* tal_fmi);

int64_t bwtSeedStrategyAllPosOneThread_with_info(
                                                   int32_t numReads,
                                                   int32_t minSeedLen,
                                                   SMEM *matchArray,
						FMI_search* tal_fmi,
						Info* qs, threadData & td, uint64_t qbwt_n);

#endif
