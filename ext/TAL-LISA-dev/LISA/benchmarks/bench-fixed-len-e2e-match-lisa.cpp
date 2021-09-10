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
#ifdef VTUNE_ANALYSIS
#include <ittnotify.h>
#endif


#include<fstream>
#include "common.h"
#include "read.h"
#include <immintrin.h>
#include "sais.h"
#include "ipbwt_rmi.h"
#include <omp.h>
#ifdef _64BIT 
    typedef int64_t index_t;
#else 
    typedef uint32_t index_t;
#endif 
/*
struct Info {
    const char* p;
    int q_len;
    int l, r; //[l, r)
    pair<index_t, index_t> intv, next_intv;
    uint64_t nxt_ext;
};
*/
int64_t one_calls = 0;

int64_t bin_search_walk = 0;





// TODO try only do chunk-based on the first SMEM?

int main(int argc, char** argv) {
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
    

    if(!(argc == 5 || argc == 6)) {
        error_quit("Need 5 args: ref_file query_set K num_rmi_leaf_nodes num_threads");
    }

    int K = atoi(argv[3]);
    eprintln("using K = %d", K);
	
    if(K < 0 || K > 30) return 0;

    int64_t num_rmi_leaf_nodes = atol(argv[4]);
    eprintln("using num_rmi_leaf_nodes = %ld", num_rmi_leaf_nodes);

    string seq;// = read_seq(argv[1]);
    read_seq_lisa(argv[1], seq);
    eprintln("Read ref file done.");
    eprintln("seq.size() = %lu", seq.size());
    string queries; int max_query_len = 0;
    tie(queries, max_query_len) = read_query_separated_with_dot(argv[2]);


    eprintln("Read query file done.");



#ifdef REV_COMP
    eprintln("No char placed between ref seq and reverse complement, to replicate BWA-MEM bug.");
    // appending reverse complement
    for(int64_t i=(index_t)seq.size()-1-(seq.back()=='@');i>=0;i--) {
#ifndef NO_DNA_ORD 
        seq.push_back(dna[3-dna_ord(seq[i])]);
#else 
        seq.push_back(dna[3-(__lg(seq[i]-'A'+2)-1)]);
#endif 
    }
#endif 

    seq.push_back('$');



    


string ref_seq_filename = argv[1];
#ifdef REV_COMP
    string rmi_filename = ref_seq_filename + ".qbwt4.walg.rev_comp";
#else
    string rmi_filename = ref_seq_filename + ".qbwt4.walg";
#endif


    //IPBWT_RMI<index_t, uint64_t> rmi(seq, rmi_filename, sa.data());
 


   IPBWT_RMI<index_t, uint64_t> rmi(seq, seq.size(), rmi_filename, K, num_rmi_leaf_nodes, NULL);
//  IPBWT_RMI<index_t, uint64_t> rmi(seq, 3101804740, rmi_filename, K, num_rmi_leaf_nodes, NULL);

    int64_t numMatches = 0;
    //int64_t sorting_ticks = 0;
    //int64_t one_ticks = 0;
    //int64_t chunk_ticks = 0;
    int64_t totalTicks = 0;
        vector<Info> qs;
        {
            for(int64_t i=0, j; i<(int64_t)queries.size(); i = j + 1) {
                for(j=i; queries[j] != '.' && queries[j] != ';'; j++){}
                if(j == i) continue;
                Info q;
                q.p = queries.c_str() + i;
 //             q.q_len = j-i;
                q.l = q.r = j - i;//q.q_len;
                q.intv = {0, rmi.n};
#ifdef NO_DNA_ORD 
                for(int k=0; k<q.r; k++) {
                    queries[i+k] = __lg(queries[i+k]-'A'+2)-1; // "ACGT" -> 0123
                }
#endif 
                qs.push_back(q);
            }
        }

    uint64_t *str_enc = (uint64_t *)malloc(qs.size() * sizeof(uint64_t));
    
    int64_t q_size = qs.size();
    int64_t *intv_all = (int64_t *)malloc(q_size * 2 * sizeof(int64_t));
    //int64_t readIntv;

    assert(str_enc != NULL && intv_all != NULL);

	int numThreads = atoi(argv[5]);//56;
//	int numThreads = 1;

  #pragma omp parallel num_threads(numThreads)
  {

	int id = omp_get_thread_num();
	if(id == 0)
		eprintln("Thread created");

  }



    for(int64_t j = 0; j < q_size; j++)
    {
        intv_all[2 * j] = 0;
        intv_all[2 * j + 1] = rmi.n;
    }




#ifdef VTUNE_ANALYSIS
        __itt_resume();
#endif

totalTicks -= __rdtsc();
uint64_t matchCount = 0;
int optimal_num_threads = min(34, numThreads);


int num_iter = max_query_len;
while(num_iter)
{
	num_iter -= K;
#pragma omp parallel num_threads(optimal_num_threads)
//#pragma omp parallel num_threads(1)
{
	int64_t workTicks = 0;
	int64_t q_processed = 0;	

    #pragma omp for schedule (static)
    for(int64_t j = 0; j < q_size; j++)
    {
	    
//	int64_t start =  __rdtsc();	
    	uint64_t nxt_ext = 0;
        Info &q = qs[j];


#ifndef NO_DNA_ORD
	for(int it = q.l-K; it!=q.l; it++) {
	    nxt_ext = (nxt_ext<<2) | dna_ord(q.p[it]); 
	}
#else			
	for(int it = q.l-K; it != q.l; it++) {
	    nxt_ext = (nxt_ext<<2) | (q.p[it]); 
	}
#endif
        str_enc[j] = nxt_ext;
	q.l = q.l - K;			
	q_processed ++;
    }
            
}



 

#if ENABLE_PREFETCH_OPT

	int64_t parallel_batch_size  = ceil((q_size/numThreads + 1)/80);
	#pragma omp parallel num_threads(numThreads)
	{
		int64_t workTicks = 0;
		int64_t q_processed = 0;	
		#pragma omp for schedule(dynamic, 1) 
		for(int64_t i = 0; i < q_size; i = i + parallel_batch_size){
			//int tid = omp_get_thread_num();			
//				
			int64_t qs_sz = ((i + parallel_batch_size) <= q_size)? parallel_batch_size: q_size - i;


//			rmi.backward_extend_chunk_batched(&str_enc[0], q_size, intv_all); 
			rmi.backward_extend_chunk_batched(&str_enc[i], qs_sz, &intv_all[i*2]); 
			q_processed += qs_sz;
		}
	}



#else
    	#pragma omp parallel for num_threads(numThreads)
        for(int64_t i = 0; i < q_size; i++)
        {
            auto q = qs[i];
            auto q_intv = rmi.backward_extend_chunk(str_enc[i], {intv_all[2 * i], intv_all[2 * i + 1]});//}q.intv);
    
	    intv_all[2 * i] = q_intv.first;
	    intv_all[2*i + 1] = q_intv.second;
        }
#endif


}

totalTicks += __rdtsc();

#ifdef VTUNE_ANALYSIS
        __itt_pause();
#endif



        for(int64_t i = 0; i < q_size; i++)
        {
            int64_t nm = intv_all[2 * i + 1] - intv_all[2 * i];
            if(nm > 0)
            {
                numMatches += nm;
		matchCount++;
            }
        }



    int64_t num_queries = count(queries.begin(), queries.end(), ';');
    assert(num_queries > 0);
    eprintln("Search Done.");
    eprintln("Number of exact matchs = %lld, match counnt %lld num queries %lld qs_size %lld", (long long)numMatches, (long long)matchCount, (long long)num_queries, (long long)q_size);
    eprintln("totalTicks = %lld", (long long)totalTicks);
    eprintln("Ticks per query = %.3f", (double)(totalTicks * 1.0 / num_queries));
    eprintln("%lld: Binary search per query = %.3f", num_rmi_leaf_nodes,(double)(bin_search_walk * 1.0 / num_queries));


#ifdef PRINT_OUTPUT
        for(int64_t i = 0; i < q_size; i++)
        {
            int64_t nm = intv_all[2 * i + 1] - intv_all[2 * i];
            printf("%ld: ",i);
            if(nm > 0)
            {
            	printf("%ld, %ld", intv_all[2 * i], intv_all[2 * i + 1]);
            }
          printf("\n");
        }
#endif
    free(str_enc); free(intv_all);
    return 0;
}
#undef flip
#undef rev_comp

