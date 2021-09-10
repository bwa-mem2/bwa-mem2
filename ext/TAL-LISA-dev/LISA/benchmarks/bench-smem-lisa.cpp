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
#include "qbwt-rmi-batched.h"
#include "read.h"
#include <fstream>
// Batch pools
#define pool_size 30000
int batch_size = 10000;


bool smem_sort(SMEM_out a, SMEM_out b) { 
	return a.id < b.id || (a.id == b.id && a.q_l < b.q_l);
}

int smem_qsort(const void* a_p, const void* b_p) { 

	SMEM_out a = *((SMEM_out*) a_p);
	SMEM_out b = *((SMEM_out*) b_p);

	if(a.id < b.id || (a.id == b.id && a.q_l < b.q_l)) return -1;
	else return 1;
}

int main(int argc, char** argv) {
#ifdef VTUNE_ANALYSIS
	__itt_pause();
#endif
	if(argc != 7) {
		error_quit("Need 6 args: ref_file query_set k num_rmi_leaf_nodes number_of_threads min_seed_length");
	}
	uint64_t totalTicks = 0;
	int K = atoi(argv[3]);
	eprintln("using K = %d", K);

	int64_t num_rmi_leaf_nodes = atol(argv[4]);
	eprintln("using num_rmi_leaf_nodes = %ld", num_rmi_leaf_nodes);



    	string seq;
//  	read_seq_lisa(argv[1], seq);
	eprintln("Read ref file done.");
	eprintln("seq.size() = %lu", seq.size());
	string queries; int max_query_len = 0;
	tie(queries, max_query_len) = read_query_separated_with_dot(argv[2]);
	eprintln("Read query file done.");


	int numThreads = atoi(argv[5]);
	int min_seed_len = atoi(argv[6]);
	if(numThreads < 0 || numThreads > 128) return 0;
	if(min_seed_len < 0 || min_seed_len > 151) return 0;
	assert(numThreads > 0 && 128);

	threadData* v_td = (threadData*) malloc(numThreads * sizeof(threadData));
	assert(v_td != NULL);
	for(int i = 0; i < numThreads; i++){
		threadData td(pool_size);
		v_td[i] = td;
	}

	int seq_size = seq.size();
	string size_file_name = (string) argv[1] + "_SIZE";
	ifstream fi(size_file_name.c_str());
	int64_t size_file;
	fi>>size_file;
	eprintln("Read ref file done. %lld",  size_file);
	 

	QBWT_HYBRID<index_t> qbwt(seq, size_file, argv[1], K, num_rmi_leaf_nodes);


	int64_t qs_size ;
	vector<Info> qs;
	int idCount = 0;
	{
		for(int64_t i=0, j; i<(int64_t)queries.size(); i = j + 1) {
			for(j=i; queries[j] != '.' && queries[j] != ';'; j++);
			if(j == i) continue;
			Info q;
			q.p = queries.c_str() + i;
			q.l = q.r = j-i;
			q.intv = {0, qbwt.n};

#ifdef NO_DNA_ORD 
			for(int k=0; k<q.r; k++) {
				queries[i+k] = __lg(queries[i+k]-'A'+2)-1; // "ACGT" -> 0123
			}
#endif 
			if(j - i >= min_seed_len)
			{
				q.id = idCount++;
				qs.push_back(q);
			}
		}
	}


	qs_size = qs.size();

#pragma omp parallel num_threads(numThreads)
	{
		eprintln("Thread created");
	}

#ifdef VTUNE_ANALYSIS
	__itt_resume();
#endif
	totalTicks -= __rdtsc();

	int64_t TotalSMEM = 0;
	int64_t parallel_batch_size  = ceil((qs_size/numThreads + 1)/80);

	Output *output = (Output*) malloc(numThreads * sizeof(Output));
	int64_t perThreadQuota = qs.size()/numThreads + 1;

	int64_t num_batches = ceil((double)qs_size/parallel_batch_size);	
	SMEM_out** batch_start = (SMEM_out**)malloc(num_batches * sizeof(SMEM_out*));
	int64_t* num_smem_per_batch = (int64_t*)malloc(num_batches * sizeof(int64_t)); 
	
	
#pragma omp parallel num_threads(numThreads)
	{
		int64_t workTicks = 0;
		int64_t q_processed = 0;
		int64_t vAnsAllocation = perThreadQuota * 20;
		int tid = omp_get_thread_num();	

		output[tid] = Output(tid);
		output[tid].smem = (SMEM_out*) malloc(vAnsAllocation * sizeof(SMEM_out));

#pragma omp for schedule(dynamic, 1) 
		for(int64_t i = 0; i < qs_size; i = i + parallel_batch_size){

			int64_t qs_sz = ((i + parallel_batch_size) <= qs_size)? parallel_batch_size: qs_size - i;
			int64_t batch_id = i/parallel_batch_size;
			int64_t prev_smem_count = v_td[tid].numSMEMs;	
			batch_start[batch_id] = output[tid].smem + prev_smem_count;
		
			// Ensures thread local memory is sufficient
			if((vAnsAllocation -  v_td[tid].numSMEMs < max_query_len * qs_sz)){
				eprintln("Insufficient memory!! Allocating more memory for stroing SMEMs");
				vAnsAllocation *= 2;
                
				output[tid].smem = (SMEM_out *)realloc(output[tid].smem, vAnsAllocation * sizeof(SMEM_out)); 
			}
			// SMEM search
			smem_rmi_batched(&qs[i], qs_sz, batch_size, qbwt, v_td[tid], &output[tid], min_seed_len);
		
			num_smem_per_batch[batch_id] = v_td[tid].numSMEMs - prev_smem_count;
			
			//Sort SMEM to rearrange the SMEMs within a batch			
			sort(batch_start[batch_id], batch_start[batch_id] + num_smem_per_batch[batch_id], smem_sort);
			//qsort(batch_start[batch_id], num_smem_per_batch[batch_id], sizeof(SMEM_out), smem_qsort);
		}
	}
	totalTicks += __rdtsc();

#ifdef VTUNE_ANALYSIS
	__itt_pause();
#endif


	for(int i = 0; i < numThreads; i++)
		TotalSMEM += v_td[i].numSMEMs;

	int64_t num_queries = count(queries.begin(), queries.end(), ';');
	assert(num_queries > 0);
	eprintln("Search Done.");
	eprintln("numSMEMs = %lld", (long long)TotalSMEM);
	eprintln("SMEMs per query = %.3f",  TotalSMEM * 1.0 / num_queries);
	eprintln("totalTicks = %lld", (long long)totalTicks);



#ifdef PRINT_OUTPUT	


	int64_t prev_qid = 0;
	for(int i = 0; i< num_batches; i++)
	{
		int64_t num_smem =  num_smem_per_batch[i];

		for(int j = 0; j < num_smem; j++){
			while (prev_qid <= batch_start[i][j].id){
				printf("%lld:\n", prev_qid);
				prev_qid++;
			}
			printf("[%d,%d][%lld, %lld]\n", batch_start[i][j].q_l, batch_start[i][j].q_r, batch_start[i][j].ref_l, batch_start[i][j].ref_r - batch_start[i][j].ref_l);		
		}
	}

#endif


#if 0
	for (int i = 0; i < output.size(); i++){
		eprintln("%d:", output[i].id);
		for(int j = output[i].qPos.size() - 1; j >= 0; j--)
			eprintln("[%d,%d][%lld, %lld]", output[i].qPos[j].first, output[i].qPos[j].second, (long long)output[i].refPos[j].first, (long long)(output[i].refPos[j].second - output[i].refPos[j].first));
	}

#endif
	free(v_td);
	return 0;
}
#undef flip
#undef rev_comp
