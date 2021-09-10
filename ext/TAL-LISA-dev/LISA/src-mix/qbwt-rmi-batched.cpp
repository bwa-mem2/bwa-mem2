#include "qbwt-rmi-batched.h"

threadData::threadData(int64_t pool_size){

	numSMEMs = 0;
	chunk_pool = (Info*) aligned_alloc(64, sizeof(Info)*pool_size);
	chunk_cnt = 0;
	str_enc = (uint64_t *)aligned_alloc(64, pool_size * sizeof(uint64_t));
	intv_all = (int64_t *)aligned_alloc(64, pool_size * 2 * sizeof(int64_t));
	fmi_pool = (Info*) aligned_alloc(64, sizeof(Info)*pool_size);
	fmi_cnt = 0;

	tree_pool = (Info*) aligned_alloc(64, sizeof(Info)*pool_size);
	tree_cnt = 0;
	index_t* s_siz_one_d = (index_t*) aligned_alloc(64, sizeof(index_t) * 2 * pool_size);
	s_siz = (index_t**) aligned_alloc(64, sizeof(index_t*) * pool_size);
	for(int64_t i = 0; i < 2*pool_size; i = i + 2)
		s_siz[i/2] = &s_siz_one_d[i];


	QBWT_HYBRID<index_t>::LcpInfo* s_info_one_d = (QBWT_HYBRID<index_t>::LcpInfo*) aligned_alloc(64, sizeof(QBWT_HYBRID<index_t>::LcpInfo) * 2 * pool_size);

	s_info = (QBWT_HYBRID<index_t>::LcpInfo**) aligned_alloc(64, sizeof(QBWT_HYBRID<index_t>::LcpInfo*)*pool_size);
	for(int64_t i = 0; i < 2*pool_size; i = i + 2)
		s_info[i/2] = &s_info_one_d[i];

	s_msk = (uint8_t*) aligned_alloc(64, sizeof(uint8_t) * pool_size);

}

void threadData::dealloc_td(){

	free(chunk_pool);// = (Info*) aligned_alloc(64, sizeof(Info)*pool_size);
	free(str_enc);// = (uint64_t *)aligned_alloc(64, pool_size * sizeof(uint64_t));
	free(intv_all);// = (int64_t *)aligned_alloc(64, pool_size * 2 * sizeof(int64_t));
	free(fmi_pool);// = (Info*) aligned_alloc(64, sizeof(Info)*pool_size);

	free(tree_pool);// = (Info*) aligned_alloc(64, sizeof(Info)*pool_size);
	free(s_siz[0]);// = (index_t*) aligned_alloc(64, sizeof(index_t) * 2 * pool_size);

	free(s_siz);

	free(s_info[0]);// = (QBWT_HYBRID<index_t>::LcpInfo*) aligned_alloc(64, sizeof(QBWT_HYBRID<index_t>::LcpInfo) * 2 * pool_size);

	free(s_info);// = (QBWT_HYBRID<index_t>::LcpInfo**) aligned_alloc(64, sizeof(QBWT_HYBRID<index_t>::LcpInfo*)*pool_size);

	free(s_msk);// = (uint8_t*) aligned_alloc(64, sizeof(uint8_t) * pool_size);

}



inline void s_pb(QBWT_HYBRID<index_t> &qbwt, Info &_q, int cnt, threadData &td) {

    Info &q = td.tree_pool[cnt]; 
    index_t* siz = td.s_siz[cnt]; 
    QBWT_HYBRID<index_t>::LcpInfo* info = td.s_info[cnt]; 
    uint8_t &msk = td.s_msk[cnt]; 
    q = _q;
    info[0] = qbwt.lcpi[q.intv.first]; info[1] = qbwt.lcpi[q.intv.second];
    siz[0] = GET_WIDTH(info[0].s_width, q.intv.first); siz[1] = GET_WIDTH(info[1].s_width, q.intv.second);
#ifndef NO_DNA_ORD
    msk = 1<<dna_ord(q.p[q.l-1]);
#else 
    msk = 1<<q.p[q.l-1];
#endif 
}



void tree_shrink_batched( QBWT_HYBRID<index_t> &qbwt, int cnt, threadData &td){

	Info* tree_pool = td.tree_pool;
	index_t** s_siz = td.s_siz;
	QBWT_HYBRID<index_t>::LcpInfo** s_info = td.s_info;
	uint8_t* s_msk = td.s_msk;

	Info* fmi_pool = td.fmi_pool;
	int &fmi_cnt = td.fmi_cnt;

	int pref_dist = 50;
	int shrink_batch_size = pref_dist = min(pref_dist, cnt);
	pref_dist = shrink_batch_size;

	while(shrink_batch_size > 0) {
		for(int i = 0; i < shrink_batch_size; i++){
			S_LOAD(i);
			S_RUN;
			S_PREFETCH;
		}
	}
}


void fmi_shrink_batched( QBWT_HYBRID<index_t> &qbwt, int cnt, Info* q_batch, threadData &td, Output* output, int min_seed_len){

	Info* tree_pool = td.tree_pool;
	int &tree_cnt = td.tree_cnt;
	
	int pref_dist = 30;
	int fmi_batch_size = pref_dist = min(pref_dist, cnt);
	pref_dist = fmi_batch_size;

	Info pf_batch[fmi_batch_size];
	
	auto cnt1 = fmi_batch_size;

	// prepare first batch
	for(int i = 0; i < fmi_batch_size; i++){
		Info &q = q_batch[i];	
		static constexpr int INDEX_T_BITS = sizeof(index_t)*__CHAR_BIT__;
		static constexpr int shift = __lg(INDEX_T_BITS);
		auto ls = ((q.intv.first>>shift)<<3), hs = ((q.intv.second>>shift)<<3); 
		my_prefetch((const char*)(qbwt.fmi->occb + ls), _MM_HINT_T0); 
		my_prefetch((const char*)(qbwt.fmi->occb + hs), _MM_HINT_T0); 
                my_prefetch((const char*)(q.p + q.l) , _MM_HINT_T0); 
	}
	
	while(fmi_batch_size > 0) {

		for(int i = 0; i < fmi_batch_size; i++){
			Info &q = q_batch[i];	
			//process one step
			int it = q.l;

			if(it < q.r){			
				auto next = qbwt.fmi->backward_extend({q.intv.first, q.intv.second}, 3 - q.p[it]);
				if(next.low >= next.high) { 
				
					tree_pool[tree_cnt++] = q;  //State change

					if(cnt1< cnt) //More queries to be processed?
						q = q_batch[cnt1++]; //direction +
					else
						q = q_batch[--fmi_batch_size];
          			      	my_prefetch((const char*)(q.p + q.l) , _MM_HINT_T0);

				}
				else {
					q.l = it + 1, q.intv={next.low, next.high};  //fmi-continue
				}
			}	
			else{
					//query finished
					if(cnt1 < cnt) //More queries to be processed?
						q = q_batch[cnt1++];
					else
						q = q_batch[--fmi_batch_size];
                			my_prefetch((const char*)(q.p + q.l) , _MM_HINT_T0);
			}
			static constexpr int INDEX_T_BITS = sizeof(index_t)*__CHAR_BIT__;
			static constexpr int shift = __lg(INDEX_T_BITS);
			auto ls = ((q.intv.first>>shift)<<3), hs = ((q.intv.second>>shift)<<3); 
			my_prefetch((const char*)(qbwt.fmi->occb + ls + 4), _MM_HINT_T0); 
			my_prefetch((const char*)(qbwt.fmi->occb + hs + 4), _MM_HINT_T0); 
		}
	}
}

void fmi_extend_batched_exact_search( QBWT_HYBRID<index_t> &qbwt, int cnt, Info* q_batch, threadData &td, Output* output, int min_seed_len){

	Info* tree_pool = td.tree_pool;
	int &tree_cnt = td.tree_cnt;
	
	int pref_dist = 30;
	int fmi_batch_size = pref_dist = min(pref_dist, cnt);
	pref_dist = fmi_batch_size;

	Info pf_batch[fmi_batch_size];
	
	auto cnt1 = fmi_batch_size;

	// prepare first batch
	for(int i = 0; i < fmi_batch_size; i++){
		Info &q = q_batch[i];	
		static constexpr int INDEX_T_BITS = sizeof(index_t)*__CHAR_BIT__;
		static constexpr int shift = __lg(INDEX_T_BITS);
		auto ls = ((q.intv.first>>shift)<<3), hs = ((q.intv.second>>shift)<<3); 
		my_prefetch((const char*)(qbwt.fmi->occb + ls), _MM_HINT_T0); 
		my_prefetch((const char*)(qbwt.fmi->occb + hs), _MM_HINT_T0); 
                my_prefetch((const char*)(q.p + q.l -  1) , _MM_HINT_T0); 
	}
	
	while(fmi_batch_size > 0) {

		for(int i = 0; i < fmi_batch_size; i++){
			Info &q = q_batch[i];	
			//process one step
			int it = q.l - 1;

			if(it >=0 && (int)q.p[it] < 4){ // considering encoding acgt -> 0123			
				auto next = qbwt.fmi->backward_extend({q.intv.first, q.intv.second}, q.p[it]);
				//fprintf(stderr, "fmi log: %d %ld %ld %ld %ld %ld\n",q.p[it], q.mid - q.r, q.mid - q.l, next.high, next.low , next.high - next.low);
				if((next.high - next.low) < q.min_intv && (q.r - (q.l-1) >= min_seed_len)) { 
				
					q.l--;
					q.intv={next.low, next.high};
					if(next.high - next.low > 0) { 
						if(q.r - q.l >= min_seed_len){
	#ifdef OUTPUT
							//[0, 2] == [98, 101]
				//			fprintf(stderr, "**** fmi log output ******: %ld %ld %ld %ld %ld\n", q.mid - q.r, q.mid - q.l - 1, next.low, next.high , next.high - next.low);
							output->tal_smem[td.numSMEMs].rid = q.id;                                                                                               			
							output->tal_smem[td.numSMEMs].n = q.l;                                                                                               			
							output->tal_smem[td.numSMEMs].m = (q.r - 1);                                                                                               			
							output->tal_smem[td.numSMEMs].k = q.intv.first;                                                                                               			
							output->tal_smem[td.numSMEMs].l = q.mid; // Here length of the read                                                                                               			
							output->tal_smem[td.numSMEMs].s = q.intv.second - q.intv.first;
							//output->tal_smem[td.numSMEMs].smem_id = q.smem_id;                                                                                               			
							td.numSMEMs++;                                                                                               			
	#endif
						}

					}
					//TODO: update l and r pointer 
					// 2. conditional insertion to rmi pool or continue in fmi pool
					q.r = q.l;
					
					if(q.r >= qbwt.rmi->K){
						q.intv = {0, qbwt.n};
						td.chunk_pool[td.chunk_cnt++] = q;

						if(cnt1< cnt) //More queries to be processed?
							q = q_batch[cnt1++]; //direction +
						else
							q = q_batch[--fmi_batch_size];
          			      	}
					else {
						q.intv = {0, qbwt.n};
					}
					my_prefetch((const char*)(q.p + q.l -  1) , _MM_HINT_T0);

				}
				else {
					q.l = it, q.intv={next.low, next.high};  //fmi-continue
				}
			}	
			else{
					//query finished
					if((q.intv.second - q.intv.first) < q.min_intv && q.r - q.l >= min_seed_len){
#ifdef OUTPUT
						
							output->tal_smem[td.numSMEMs].rid = q.id;                                                                                               			
							output->tal_smem[td.numSMEMs].n = q.l;                                                                                               			
							output->tal_smem[td.numSMEMs].m = (q.r - 1);                                                                                               			
							output->tal_smem[td.numSMEMs].k = q.intv.first;                                                                                               			
							output->tal_smem[td.numSMEMs].l = q.mid; // Here length of the read                                                                                               			
							output->tal_smem[td.numSMEMs].s = q.intv.second - q.intv.first;
							//output->tal_smem[td.numSMEMs].smem_id = q.smem_id;                                                                                               			
							td.numSMEMs++;                                                                                               			
#endif
					}
					if(cnt1 < cnt) //More queries to be processed?
						q = q_batch[cnt1++];
					else
						q = q_batch[--fmi_batch_size];
                			my_prefetch((const char*)(q.p + q.l -  1) , _MM_HINT_T0);
			}
			static constexpr int INDEX_T_BITS = sizeof(index_t)*__CHAR_BIT__;
			static constexpr int shift = __lg(INDEX_T_BITS);
			auto ls = ((q.intv.first>>shift)<<3), hs = ((q.intv.second>>shift)<<3); 
			my_prefetch((const char*)(qbwt.fmi->occb + ls + 4), _MM_HINT_T0); 
			my_prefetch((const char*)(qbwt.fmi->occb + hs + 4), _MM_HINT_T0); 
		}
	}
}

void fmi_extend_batched( QBWT_HYBRID<index_t> &qbwt, int cnt, Info* q_batch, threadData &td, Output* output, int min_seed_len){

	Info* tree_pool = td.tree_pool;
	int &tree_cnt = td.tree_cnt;
	
	int pref_dist = 30;
	int fmi_batch_size = pref_dist = min(pref_dist, cnt);
	pref_dist = fmi_batch_size;

	Info pf_batch[fmi_batch_size];
	
	auto cnt1 = fmi_batch_size;

	// prepare first batch
	for(int i = 0; i < fmi_batch_size; i++){
		Info &q = q_batch[i];	
		static constexpr int INDEX_T_BITS = sizeof(index_t)*__CHAR_BIT__;
		static constexpr int shift = __lg(INDEX_T_BITS);
		auto ls = ((q.intv.first>>shift)<<3), hs = ((q.intv.second>>shift)<<3); 
		my_prefetch((const char*)(qbwt.fmi->occb + ls), _MM_HINT_T0); 
		my_prefetch((const char*)(qbwt.fmi->occb + hs), _MM_HINT_T0); 
                my_prefetch((const char*)(q.p + q.l -  1) , _MM_HINT_T0); 
	}
	
	while(fmi_batch_size > 0) {

		for(int i = 0; i < fmi_batch_size; i++){
			Info &q = q_batch[i];	
			//process one step
			int it = q.l -1;

			if(it >=0 && (int)q.p[it] < 4){ // considering encoding acgt -> 0123			
				auto next = qbwt.fmi->backward_extend({q.intv.first, q.intv.second}, q.p[it]);
				//if(next.low >= next.high) { //TODO: min_intv_size
				if(!((next.high - next.low) > q.min_intv)) { 
				

					if(q.r - q.l >= min_seed_len && q.l != q.prev_l){
#ifdef OUTPUT
/*						output[q.id].qPos.push_back({q.l, q.r});
						output[q.id].refPos.push_back(q.intv);
						output[q.id].smem.push_back(SMEM_out(q.id, q.l, q.r, q.intv.first, q.intv.second));
*/
//						output->smem[td.numSMEMs] = SMEM_out(q.id, q.l, q.r, q.intv.first, q.intv.second);             
						
						output->tal_smem[td.numSMEMs].rid = q.id;                                                                                               			
						output->tal_smem[td.numSMEMs].m = q.l;                                                                                               			
						output->tal_smem[td.numSMEMs].n = q.r - 1;                                                                                               			
						output->tal_smem[td.numSMEMs].k = q.intv.first;                                                                                               			
						output->tal_smem[td.numSMEMs].l = 0;                                                                                               			
						output->tal_smem[td.numSMEMs].s = q.intv.second - q.intv.first;
						//output->tal_smem[td.numSMEMs].smem_id = q.smem_id;                                                                                               			
						td.numSMEMs++;                       
						q.prev_l = q.l;                                                                        			
#endif
					}
					tree_pool[tree_cnt++] = q;  //State change
					if(cnt1< cnt) //More queries to be processed?
						q = q_batch[cnt1++]; //direction +
					else
						q = q_batch[--fmi_batch_size];
          			      	my_prefetch((const char*)(q.p + q.l -  1) , _MM_HINT_T0);

				}
				else {
					q.l = it, q.intv={next.low, next.high};  //fmi-continue
				}
			}	
			else{
					//query finished
					if(q.r - q.l >= min_seed_len && q.l != q.prev_l){
#ifdef OUTPUT
//						output[q.id].qPos.push_back({q.l, q.r});
//						output[q.id].refPos.push_back(q.intv);
//						output[q.id].smem.push_back(SMEM_out(q.id, q.l, q.r, q.intv.first, q.intv.second));

//						output->smem[td.numSMEMs++] = SMEM_out(q.id, q.l, q.r, q.intv.first, q.intv.second);                                                                                                            			
//						output->smem[td.numSMEMs] = SMEM_out(q.id, q.l, q.r, q.intv.first, q.intv.second);             
						
						output->tal_smem[td.numSMEMs].rid = q.id;                                                                                               			
						output->tal_smem[td.numSMEMs].m = q.l;                                                                                               			
						output->tal_smem[td.numSMEMs].n = q.r - 1;                                                                                               			
						output->tal_smem[td.numSMEMs].k = q.intv.first;                                                                                               			
						output->tal_smem[td.numSMEMs].l = 0;                                                                                               			
						output->tal_smem[td.numSMEMs].s = q.intv.second - q.intv.first;
						//output->tal_smem[td.numSMEMs].smem_id = q.smem_id;                                                                                               			
						td.numSMEMs++;                                                                                               			
						q.prev_l = q.l;                                                                        			
#endif
					}
					if(cnt1 < cnt) //More queries to be processed?
						q = q_batch[cnt1++];
					else
						q = q_batch[--fmi_batch_size];
                			my_prefetch((const char*)(q.p + q.l -  1) , _MM_HINT_T0);
			}
			static constexpr int INDEX_T_BITS = sizeof(index_t)*__CHAR_BIT__;
			static constexpr int shift = __lg(INDEX_T_BITS);
			auto ls = ((q.intv.first>>shift)<<3), hs = ((q.intv.second>>shift)<<3); 
			my_prefetch((const char*)(qbwt.fmi->occb + ls + 4), _MM_HINT_T0); 
			my_prefetch((const char*)(qbwt.fmi->occb + hs + 4), _MM_HINT_T0); 
		}
	}
}

void smem_rmi_batched(Info *qs, int64_t qs_size, int64_t batch_size, QBWT_HYBRID<index_t> &qbwt, threadData &td, Output* output, int min_seed_len, bool apply_lisa){
	Info *chunk_pool = td.chunk_pool;
	int &chunk_cnt = td.chunk_cnt;

	uint64_t *str_enc = td.str_enc;
	int64_t *intv_all = td.intv_all;

	Info* fmi_pool = td.fmi_pool;
	int &fmi_cnt = td.fmi_cnt;

	Info* tree_pool = td.tree_pool;
	int &tree_cnt = td.tree_cnt;

	int K = qbwt.rmi->K;		
	
	
	int64_t next_q = 0;
	while(next_q < qs_size || (chunk_cnt + fmi_cnt ) > 0){
		while(next_q < qs_size && chunk_cnt < batch_size && fmi_cnt < batch_size){
				if(apply_lisa == true && qs[next_q].r >= K){
					chunk_pool[chunk_cnt++] = qs[next_q++];
				}
				else
					fmi_pool[fmi_cnt++] = qs[next_q++];
		}
		
		// process chunk batch
		if(next_q >= qs_size || !(chunk_cnt < batch_size)){
			// prepareChunkBatchVectorized_v1(chunk_pool, chunk_cnt, str_enc, intv_all);	
			prepareChunkBatch(chunk_pool, chunk_cnt, str_enc, intv_all, K);

			qbwt.rmi->backward_extend_chunk_batched(&str_enc[0], chunk_cnt, intv_all);
			auto cnt = chunk_cnt;
			chunk_cnt = 0;
			for(int64_t j = 0; j < cnt; j++) {
				Info &q = chunk_pool[j];
				auto next_intv = make_pair(intv_all[2*j],intv_all[2*j + 1]);

				// next state: chunk batching
				//if(next_intv.first < next_intv.second){//TODO: min_intv_size
				if((next_intv.second - next_intv.first) > q.min_intv){//TODO: min_intv_size
					q.intv = next_intv;
					q.l -= K;

					if(q.l >= K)// heuristics: && !(q.intv.second - q.intv.first > 1 && q.intv.second - q.intv.first < 100))	
						chunk_pool[chunk_cnt++] = q;
					else if(q.l > 0){
						fmi_pool[fmi_cnt++] = q;
					}
					else if(q.r - q.l >= min_seed_len && q.l != q.prev_l){

#ifdef OUTPUT
						// output[q.id].qPos.push_back({q.l, q.r});
						// output[q.id].refPos.push_back(q.intv);
						// output[q.id].smem.push_back(SMEM_out(q.id, q.l, q.r, q.intv.first, q.intv.second));
						//output->smem[td.numSMEMs++] = SMEM_out(q.id, q.l, q.r, q.intv.first, q.intv.second);                                                                                                            			
//						output->smem[td.numSMEMs] = SMEM_out(q.id, q.l, q.r, q.intv.first, q.intv.second);             
						
						output->tal_smem[td.numSMEMs].rid = q.id;                                                                                               			
						output->tal_smem[td.numSMEMs].m = q.l;                                                                                               			
						output->tal_smem[td.numSMEMs].n = q.r - 1;                                                                                               			
						output->tal_smem[td.numSMEMs].k = q.intv.first;                                                                                               			
						output->tal_smem[td.numSMEMs].l = 0;                                                                                               			
						output->tal_smem[td.numSMEMs].s = q.intv.second - q.intv.first;
						//output->tal_smem[td.numSMEMs].smem_id = q.smem_id;                                                                                               			
						td.numSMEMs++;                                                                                               			
						q.prev_l = q.l;                                                                        			
#endif
					}
				}
				else {
					//Next state: fmi procssing
					fmi_pool[fmi_cnt++] = q;
				}

			}
		}
		
		// fmi processing
		if(next_q >= qs_size || !(fmi_cnt < batch_size)){
			auto cnt = fmi_cnt;
			fmi_cnt = 0;	

			fmi_extend_batched(qbwt, cnt, &fmi_pool[0], td, output, min_seed_len);	

			cnt = tree_cnt;
			tree_cnt = 0;
			for(int i = 0; i < cnt; i++) {
				auto &q = tree_pool[i];
				s_pb(qbwt, q, i, td);
				my_prefetch((const char*)(qbwt.lcpi + tree_pool[i+50].intv.first) , _MM_HINT_T0);
				my_prefetch((const char*)(qbwt.lcpi + tree_pool[i+50].intv.second) , _MM_HINT_T0);
				my_prefetch((const char*)(tree_pool[i + 50].p + tree_pool[i + 50].l - 1) , _MM_HINT_T0);
			}
			tree_shrink_batched(qbwt, cnt, td);
		}
	}


}

void exact_search_rmi_batched(Info *qs, int64_t qs_size, int64_t batch_size, QBWT_HYBRID<index_t> &qbwt, threadData &td, Output* output, int min_seed_len, bool apply_lisa){
	Info *chunk_pool = td.chunk_pool;
	int &chunk_cnt = td.chunk_cnt;

	uint64_t *str_enc = td.str_enc;
	int64_t *intv_all = td.intv_all;

	Info* fmi_pool = td.fmi_pool;
	int &fmi_cnt = td.fmi_cnt;


	int K = qbwt.rmi->K;		
	
	
	int64_t next_q = 0;
	while(next_q < qs_size || (chunk_cnt + fmi_cnt ) > 0){
		while(next_q < qs_size && chunk_cnt < batch_size && fmi_cnt < batch_size){
				if(apply_lisa == true && qs[next_q].r >= K){
					chunk_pool[chunk_cnt++] = qs[next_q++];
				}
				else
					fmi_pool[fmi_cnt++] = qs[next_q++];
		}
		
		// process chunk batch
		if(next_q >= qs_size || !(chunk_cnt < batch_size)){
			// prepareChunkBatchVectorized_v1(chunk_pool, chunk_cnt, str_enc, intv_all);	
			prepareChunkBatch(chunk_pool, chunk_cnt, str_enc, intv_all, K);

			qbwt.rmi->backward_extend_chunk_batched(&str_enc[0], chunk_cnt, intv_all);
			auto cnt = chunk_cnt;
			chunk_cnt = 0;
			for(int64_t j = 0; j < cnt; j++) {
				Info &q = chunk_pool[j];
				auto next_intv = make_pair(intv_all[2*j],intv_all[2*j + 1]);

				// next state: chunk batching
				//if(next_intv.first < next_intv.second){//TODO: min_intv_size
				
				//fprintf(stderr, "rmi log: %ld %ld %ld %ld %ld %ld\n", j, q.l, q.r, next_intv.second , next_intv.first , q.min_intv);
				if((next_intv.second - next_intv.first) > q.min_intv){//TODO: min_intv_size
					q.intv = next_intv;
					q.l -= K;

					if(q.l >= K)// heuristics: && !(q.intv.second - q.intv.first > 1 && q.intv.second - q.intv.first < 100))	
						chunk_pool[chunk_cnt++] = q;
					else if(q.l > 0){
						fmi_pool[fmi_cnt++] = q;
					}
					else if(false && q.r - q.l >= min_seed_len){

#ifdef OUTPUT
						// output[q.id].qPos.push_back({q.l, q.r});
						// output[q.id].refPos.push_back(q.intv);
						// output[q.id].smem.push_back(SMEM_out(q.id, q.l, q.r, q.intv.first, q.intv.second));
						//output->smem[td.numSMEMs++] = SMEM_out(q.id, q.l, q.r, q.intv.first, q.intv.second);                                                                                                            			
//						output->smem[td.numSMEMs] = SMEM_out(q.id, q.l, q.r, q.intv.first, q.intv.second);             
						
						output->tal_smem[td.numSMEMs].rid = q.id;                                                                                               			
						output->tal_smem[td.numSMEMs].m = q.l;                                                                                               			
						output->tal_smem[td.numSMEMs].n = q.r - 1;                                                                                               			
						output->tal_smem[td.numSMEMs].k = q.intv.first;                                                                                               			
						output->tal_smem[td.numSMEMs].l = 0;                                                                                               			
						output->tal_smem[td.numSMEMs].s = q.intv.second - q.intv.first;
						//output->tal_smem[td.numSMEMs].smem_id = q.smem_id;                                                                                               			
						td.numSMEMs++;                                                                                               			
#endif
					}
				}
				else {
					//Next state: fmi procssing
					fmi_pool[fmi_cnt++] = q;
				}

			}
		}
		
		// fmi processing
		if(next_q >= qs_size || !(fmi_cnt < batch_size)){
			auto cnt = fmi_cnt;
			fmi_cnt = 0;	

			fmi_extend_batched_exact_search(qbwt, cnt, &fmi_pool[0], td, output, min_seed_len);	

		}
	}


}

