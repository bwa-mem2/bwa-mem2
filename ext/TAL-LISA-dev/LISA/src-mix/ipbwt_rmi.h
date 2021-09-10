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

#ifndef IPBWT_RMI_H 
#define IPBWT_RMI_H
#ifdef VECTORIZE
#include <immintrin.h>
#endif 
#include "common.h"
#include <omp.h>
#include <unistd.h>
#include <sys/mman.h>

#ifdef ENABLE_PREFETCH
enum query_state
{
    GUESS_RMI_ROOT,
    GUESS_RMI_LEAF,
    LAST_MILE
};
#endif

        

template<typename index_t, typename kenc_t>
class IPBWT_RMI {
    public:
        IPBWT_RMI(){}
        IPBWT_RMI(const string &t, string ref_seq_filename, int K_, int64_t num_rmi_leaf_nodes, index_t *__sa = NULL);
        IPBWT_RMI(const string &t, index_t t_size, string ref_seq_filename, int K_, int64_t num_rmi_leaf_nodes, index_t *__sa = NULL);
        ~IPBWT_RMI();

        struct ref_pos_t;
#ifdef ENABLE_PREFETCH
        struct BatchMetadata;
#endif
        int K;
        double inv_second_size;
        int NUM_POS_BITS = 38;
        uint64_t POS_MASK = 0x3fffffffffL;
        int NUM_CHUNK_BITS = 42;
        uint64_t CHUNK_MASK = 0x3ffffffffffL;
        int NUM_IPBWT_BITS = NUM_CHUNK_BITS + NUM_POS_BITS;
        uint64_t NUM_IPBWT_BYTES = (NUM_IPBWT_BITS + 7 ) / 8;
        __mmask64 mask_ipbwt_load = (1 << NUM_IPBWT_BYTES) - 1;
        double L0_PARAMETER0;
        double L0_PARAMETER1;
        int64_t L1_SIZE;
        double *L1_PARAMETERS;
        int16_t m0[32];
        int16_t m1[32];
        int32_t m_one_bits[9];
        pair<index_t, index_t> backward_extend_chunk(kenc_t str_enc, pair<index_t, index_t> intv) const;
        kenc_t gen_str_enc(const char *a, int len) const;
        int64_t bs_ticks;
        int64_t vbs_ticks;
        int64_t bs_calls;
        int64_t vbs_calls;

    // private:
        static constexpr int kenc_bits = sizeof(kenc_t)*__CHAR_BIT__;
        index_t n, second_size;

        uint8_t *ipbwt_array;
#ifdef ENABLE_PREFETCH
        static constexpr int64_t UNROLL = 28;

#if 0
        inline bool process_partial_query_one_step(BatchMetadataPartial &meta, index_t &p_low, index_t &p_high, bool is_partial);
        void backward_extend_chunk_batched_info (Info* str_enc_list, int64_t qs_size, index_t *intv_all, bool is_partial);
        void backward_partial_extend_chunk_batched (kenc_t* str_enc_list, int64_t qs_size, index_t *intv_all, bool is_partial);
        void backward_partial_extend_chunk_batched_v1 (kenc_t* str_enc_list, int64_t qs_size, index_t *intv_all, bool is_partial);
#endif
        inline int64_t get_guess_root_step(double key);
        inline int64_t get_guess_leaf_step(double key, int64_t modelIndex, int64_t *err);
        inline void last_mile_binary_search_one_step(ipbwt_t ipb_x, int64_t &first, int64_t &m);
        inline void last_mile_vectorized_search_final_step(ipbwt_t ipb_x, int64_t &first, int64_t &m);
        inline bool process_query_one_step(BatchMetadata &meta, index_t &p_low, index_t &p_high);
        void backward_extend_chunk_batched (kenc_t* str_enc_list, int64_t qs_size, index_t *intv_all);
        void backward_extend_multi_chunk_batched (Info* str_enc_list, int64_t qs_size, index_t *intv_all);
#endif

        double to_floating_point(pair<uint64_t, uint64_t> p) const;
        index_t get_guess_from_rmi(ipbwt_t ipb_x, index_t *err) const;
        index_t last_mile_from_guess(ipbwt_t ipb_x, index_t guess, index_t err) const;
/*
	
inline ipbwt_t ipbwt(uint64_t chunk, uint64_t pos) const 
{ 
    uint64_t second = (pos & POS_MASK) + (chunk << NUM_POS_BITS);
    uint64_t first = chunk >> (64 - NUM_POS_BITS);
    return {first, second}; 
}

inline ipbwt_t ipbwt(const index_t &i) const 
{ 
    uint8_t *base_addr = ipbwt_array + i * NUM_IPBWT_BYTES;
    uint64_t second = ((uint64_t *)base_addr)[0];
    uint64_t first = ((uint16_t *)(base_addr + 8))[0];
    return {first, second}; 
}
*/
	ipbwt_t ipbwt(uint64_t chunk, uint64_t pos) const;
	ipbwt_t ipbwt(const index_t &i) const;
	uint64_t get_ref_pos(ref_pos_t A, index_t i);
	kenc_t get_enc_rev_comp(kenc_t nxt_ext, int K = 21);
	kenc_t get_enc_suffix(kenc_t nxt_ext, int suffix_length, int K = 21);
};

template<typename index_t, typename kenc_t> 
struct IPBWT_RMI<index_t, kenc_t>::ref_pos_t
{
    uint32_t *ls_word;
    uint8_t *ms_byte;
};

#ifdef ENABLE_PREFETCH
template<typename index_t, typename kenc_t> 
struct IPBWT_RMI<index_t, kenc_t>::BatchMetadata
{
    ipbwt_t ipb_x[2];
    double key[2];
    index_t modelIndex[2];
    index_t first[2];
    index_t m[2];
    int64_t qid;
    query_state state;
};

#endif

template<typename index_t, typename kenc_t> 
inline ipbwt_t IPBWT_RMI<index_t, kenc_t>::ipbwt(uint64_t chunk, uint64_t pos) const 
{ 
    uint64_t second = (pos & POS_MASK) + (chunk << NUM_POS_BITS);
    uint64_t first = chunk >> (64 - NUM_POS_BITS);
    return {first, second}; 
}
template<typename index_t, typename kenc_t> 
inline ipbwt_t IPBWT_RMI<index_t, kenc_t>::ipbwt(const index_t &i) const 
{ 
    uint8_t *base_addr = ipbwt_array + i * NUM_IPBWT_BYTES;
    uint64_t second = ((uint64_t *)base_addr)[0];
    uint64_t first = ((uint16_t *)(base_addr + 8))[0];
    return {first, second}; 
}
template<typename index_t, typename kenc_t> 
inline uint64_t IPBWT_RMI<index_t, kenc_t>::get_ref_pos(ref_pos_t A, index_t i)
{
    uint64_t sa_ms_byte = A.ms_byte[i];
    uint64_t sa_ls_word = A.ls_word[i];
    uint64_t pos = (sa_ms_byte << 32) + sa_ls_word;
    return pos;
}

template<typename index_t, typename kenc_t> 
inline kenc_t IPBWT_RMI<index_t, kenc_t>::get_enc_rev_comp(kenc_t nxt_ext, int K) {
		
	auto temp = nxt_ext;
	temp = ~temp;
	kenc_t new_ext = 0, mask = 3;

	for (int i = 0; i < K; i++){
		new_ext = new_ext << 2 | temp & mask;
		temp = temp >>2;
	}
	return new_ext;
}

template<typename index_t, typename kenc_t> 
inline kenc_t IPBWT_RMI<index_t, kenc_t>::get_enc_suffix(kenc_t nxt_ext, int suffix_length, int K){

	kenc_t new_ext;
	nxt_ext = nxt_ext << ((64 - 2*K) + 2*(K - (suffix_length + 1)));
	new_ext = nxt_ext >> (64 - 2*K);
	return new_ext;
}


template<typename index_t, typename kenc_t>
inline kenc_t IPBWT_RMI<index_t, kenc_t>::gen_str_enc(const char *a, int len) const {
    kenc_t str_enc = 0;
    int i;
    for(i=0; i<len && a[i]!='$'; i++) {
#ifndef NO_DNA_ORD 
        str_enc = (str_enc<<2) | dna_ord(a[i]);
#else
        str_enc = (str_enc<<2) | (a[i]);
#endif 
    }
    return str_enc << (2*(len-i)); // need this due to t.c_str() format
}


template<typename index_t, typename kenc_t>
double IPBWT_RMI<index_t, kenc_t>::to_floating_point(pair<uint64_t, uint64_t> p) const {
    return p.first + p.second * inv_second_size;
}



template<typename index_t, typename kenc_t>
IPBWT_RMI<index_t, kenc_t>::IPBWT_RMI(const string &t, index_t t_size, string ref_seq_filename, int K_, int64_t num_rmi_leaf_nodes, index_t *__sa):
    n(t_size), K(K_), second_size(n+K_) {

    assert(K <= 21);
    eprintln("n = %ld, K = %d", (long)n, K);
    eprintln("NUM_IPBWT_BYTES = %ld", NUM_IPBWT_BYTES);
    inv_second_size = ((double)1.0)/second_size;
    auto build_ipbwt = [&]() { 
        int64_t startTick, endTick;
         
        startTick = __rdtsc();
        if(__sa == NULL) {
            __sa = (int64_t *)malloc(n * sizeof(int64_t));
            assert(__sa != NULL);
            saisxx(t.c_str(), __sa, (int64_t)n);
         }
        ref_pos_t sa;
       // sa.ls_word = (uint32_t *)malloc(n * sizeof(uint32_t));
       // sa.ms_byte = (uint8_t *)malloc(n * sizeof(uint8_t));
       sa.ls_word = (uint32_t *)calloc(n, sizeof(uint32_t));
       sa.ms_byte = (uint8_t *)calloc(n, sizeof(uint8_t));
        assert(sa.ls_word != NULL);
        assert(sa.ms_byte != NULL);
        for(index_t i = 0; i < n; i++)
        {
            sa.ls_word[i] = __sa[i] & 0xffffffff;
            sa.ms_byte[i] = (__sa[i] >> 32) & 0xff;
        }
        for(index_t i = 0; i < n; i++)
        {
            uint64_t sa_ms_byte = sa.ms_byte[i];
            uint64_t sa_ls_word = sa.ls_word[i];
            uint64_t pos = (sa_ms_byte << 32) + sa_ls_word;
            if(pos != __sa[i])
            {
                eprintln("%ld] pos = %ld, __sa[i] = %ld", (long)i, (long)pos, (long)__sa[i]);
                exit(0);
            }
        }
        free(__sa);
        endTick = __rdtsc();
        eprintln("SA done in %ld cycles.", endTick - startTick);

        {
            startTick = __rdtsc();
            ref_pos_t rnk;
            rnk.ls_word = (uint32_t *)malloc(n * sizeof(uint32_t));
            rnk.ms_byte = (uint8_t *)malloc(n * sizeof(uint8_t));
            assert(rnk.ls_word != NULL);
            assert(rnk.ms_byte != NULL);
            eprintln("rnk allocated, n = %ld", (long)n);
//#pragma omp parallel for
            for(index_t i=0; i<n; i++)
            {
                uint64_t sa_ms_byte = sa.ms_byte[i];
                uint64_t sa_ls_word = sa.ls_word[i];
                uint64_t pos = (sa_ms_byte << 32) + sa_ls_word;
                rnk.ls_word[pos] = i & 0xffffffff;
                rnk.ms_byte[pos] = (i >> 32) & 0xff;
            }
            endTick = __rdtsc();
            eprintln("RNK done in %ld cycles.", endTick - startTick);
            startTick = __rdtsc();
            ipbwt_array = (uint8_t *)aligned_alloc(64, n * NUM_IPBWT_BYTES);
            assert(ipbwt_array != NULL);
        eprintln("ipbwt_array allocated with size %lu", (unsigned long)n * NUM_IPBWT_BYTES);
            assert(K<n);
//#pragma omp parallel for
            for(index_t i=0; i<n; i++) {
                uint64_t pos;
                index_t p = get_ref_pos(sa, i) + K;
                if(p>=n) {
                    p-=n;
                    pos = K-1-p; // dollar_sign_pos
                    //((int64_t *)ipbwt_array)[2 * i] = K-1-p; // dollar_sign_pos
                } else {
                    pos = get_ref_pos(rnk, p) + K;
                    //((int64_t *)ipbwt_array)[2 * i] = get_ref_pos(rnk, p) + K;
                }
                uint64_t str_enc = gen_str_enc(t.c_str() + get_ref_pos(sa, i), K);
                //eprintln("%ld] pos = %lu, str_enc = %lu", i, pos, str_enc);
                //((int64_t *)ipbwt_array)[2 * i + 1] = str_enc;
                pos = (pos & POS_MASK) + (str_enc << NUM_POS_BITS);
                uint16_t chunk = str_enc >> (64 - NUM_POS_BITS);
                //eprintln("%ld] pos = %lu, chunk = %u", i, pos, chunk);
                //fflush(stderr);
                uint64_t *pos_ptr = (uint64_t *)(&(ipbwt_array[i * NUM_IPBWT_BYTES]));
                pos_ptr[0] = pos;
                uint16_t *chunk_ptr = (uint16_t *)(&(ipbwt_array[i * NUM_IPBWT_BYTES])) + 4;
                chunk_ptr[0] = chunk;
            }
            endTick = __rdtsc();
            eprintln("ipbwt_array done in %ld cycles.", endTick - startTick);
            for(index_t i = 1; i < n; i++)
            {
                if(ipbwt(i - 1) >= ipbwt(i))
                {
                    //eprintln("%ld] {%lu, %lu} {%lu, %lu}", i, ((int64_t *)ipbwt_array)[2 * (i - 1) + 1], ((int64_t *)ipbwt_array)[2 * (i - 1)], ((int64_t *)ipbwt_array)[2 * i + 1], ((int64_t *)ipbwt_array)[2 * i]);
                    eprintln("%ld] {%lu, %lu} {%lu, %lu}", (long)i, (unsigned long)ipbwt(i - 1).first, (unsigned long)ipbwt(i - 1). second, (unsigned long)ipbwt(i).first, (unsigned long)ipbwt(i).second);
                    exit(0);
                }
            }
            free(sa.ls_word); free(sa.ms_byte);
            free(rnk.ls_word); free(rnk.ms_byte);
        }
    };


    string ipbwt_filename = ref_seq_filename + ".ipbwt";
    for(const auto &s:{sizeof(index_t), sizeof(kenc_t), (size_t)K, (size_t)NUM_POS_BITS, (size_t)NUM_CHUNK_BITS}) {
        ipbwt_filename += string(".") + to_string(s);
    }
    if(ifstream(ipbwt_filename.c_str()).good()) {
        eprintln("Found existing %s!!", (char *)ipbwt_filename.c_str());

#ifndef HUGE_PAGE
        ipbwt_array = (uint8_t *)aligned_alloc(64, n * NUM_IPBWT_BYTES);
#else
	eprintln("Using Huge-Page 2MB");
	ipbwt_array = (uint8_t *) mmap(NULL, n * NUM_IPBWT_BYTES, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB, -1, 0);
#endif
	eprintln("ipbwt Array size: %lld", (long long)n * NUM_IPBWT_BYTES);

        load(ipbwt_filename, {((char *)ipbwt_array)}, {(n * NUM_IPBWT_BYTES)}); 
        eprintln("Load successful.");
        for(index_t i = 0; i < n; i++)
        {
            if(i > 0)
            {
                if(ipbwt(i - 1) >= ipbwt(i))
                {
                    //eprintln("%ld] {%lu, %lu} {%lu, %lu}", i, ((int64_t *)ipbwt_array)[2 * (i - 1) + 1], ((int64_t *)ipbwt_array)[2 * (i - 1)], ((int64_t *)ipbwt_array)[2 * i + 1], ((int64_t *)ipbwt_array)[2 * i]);
                    eprintln("%ld] {%lu, %lu} {%lu, %lu}", (long)i, (unsigned long)ipbwt(i - 1).first, (unsigned long)ipbwt(i - 1). second, (unsigned long)ipbwt(i).first, (unsigned long)ipbwt(i).second);
                    exit(0);
                }
                assert(ipbwt(i-1) < ipbwt(i));
            }
        }
    } else {
        eprintln("No existing %s. Building...", (char*)ipbwt_filename.c_str());
        build_ipbwt();
        eprintln("%s build done.", (char *)ipbwt_filename.c_str());
        save(ipbwt_filename, {((char *)ipbwt_array)}, {(n * NUM_IPBWT_BYTES)}); 
        eprintln("%s saved.", (char*)ipbwt_filename.c_str());
    }
    eprintln("ipbwt correctness ok.");

    string rmi_filename = ipbwt_filename + "." + to_string(num_rmi_leaf_nodes) + ".rmi_PARAMETERS";
    if(ifstream(rmi_filename.c_str()).good()) {
        eprintln("Found existing %s!!", (char*)rmi_filename.c_str());
        vector<char*> ptrs;vector<size_t> sizes;
        ptrs.push_back((char *)&(L0_PARAMETER0));
        sizes.push_back(sizeof(double));
        ptrs.push_back((char *)&(L0_PARAMETER1));
        sizes.push_back(sizeof(double));
        ptrs.push_back((char *)&(L1_SIZE));
        sizes.push_back(sizeof(int64_t));
        load(rmi_filename, ptrs, sizes);

#ifndef HUGE_PAGE
  	L1_PARAMETERS = (double*) _mm_malloc(L1_SIZE * 3 * sizeof(double), 64);
#else      
	L1_PARAMETERS = (double*) mmap(NULL, L1_SIZE * 3 * sizeof(double), PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB, -1, 0);
#endif


	eprintln("L1 size in Bytes: %lld", (long long)L1_SIZE * 3 * sizeof(double));

        if (L1_PARAMETERS == NULL)
        {
            eprintln("L1_PARAMETERS == NULL");
            exit(0);;
        }
        ptrs.push_back((char *)L1_PARAMETERS);
        sizes.push_back(L1_SIZE * 3 * sizeof(double));
        load(rmi_filename, ptrs, sizes);
        eprintln("Load successful.");
        eprintln("L1_SIZE = %ld", L1_SIZE);
    } else {
        eprintln("No existing %s. Building train data...", (char*)rmi_filename.c_str());
        double *train_data = (double *)_mm_malloc(n * sizeof(double), 64);
        for(index_t i = 0; i < n; i++)
        {
            ipbwt_t p = ipbwt(i);
            uint64_t second = p.second & POS_MASK;
            uint64_t first = (p.first << (64 - NUM_POS_BITS)) + (p.second >> NUM_POS_BITS);
            train_data[i] = to_floating_point({first, second});
        }
        vector<char*> ptrs;vector<size_t> sizes;
        ptrs.push_back((char *)&n);
        sizes.push_back(sizeof(int64_t));
        ptrs.push_back((char *)train_data);
        sizes.push_back(n * sizeof(double));
        string ipbwt_f64_filename = ipbwt_filename + ".f64";
        save(ipbwt_f64_filename, ptrs, sizes);
        eprintln("Training data saved at %s. Building RMI using Ryan's code. Run this code again after the RMI is built.", (char*)ipbwt_f64_filename.c_str());
        eprintln("executing ./build-rmi.linear_spline.linear.sh");
        //const char *execv_argv[] = {"./build-rmi.linear_spline.linear.sh",
        //                            ipbwt_f64_filename.c_str(),
        //                            rmi_filename.c_str(),
        //                            to_string(num_rmi_leaf_nodes).c_str()};
        //execv_argv[0] = ipbwt_f64_filename.c_str();
        //execv_argv[1] = rmi_filename.c_str();
        //execv_argv[2] = to_string(num_rmi_leaf_nodes);
        int status = execl("./scripts/build-rmi.linear_spline.linear.sh",
              "./scripts/build-rmi.linear_spline.linear.sh",
              ipbwt_f64_filename.c_str(),
              rmi_filename.c_str(),
              to_string(num_rmi_leaf_nodes).c_str(), (char *)NULL);
        eprintln("ERROR! execl failed with err %d, errno = %d: %s", status, errno, strerror( errno));
        exit(0);
    }

    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            m0[i * 4 + j] = i * 5 + j;
            m1[i * 4 + j] = 4 + i * 5 + j;
        }
    }
    for(int i = 4; i < 8; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            m0[i * 4 + j] = i * 5 + j + 24;
            m1[i * 4 + j] = 4 + i * 5 + j + 24;
        }
    }
#if 0
    for(int i = 0; i < 32; i++)
    {
        eprintln("%d] %d %d", i, m0[i], m1[i]);
    }
#endif
    for(int i = 0; i < 9; i++)
    {
        m_one_bits[i] = (1 << i) - 1;
    }
    m_one_bits[8] = 0xff;
    bs_calls = vbs_calls = 0;
    bs_ticks = vbs_ticks = 0;
}




template<typename index_t, typename kenc_t>
IPBWT_RMI<index_t, kenc_t>::IPBWT_RMI(const string &t, string ref_seq_filename, int K_, int64_t num_rmi_leaf_nodes, index_t *__sa):
    n((index_t)t.size()), K(K_), second_size(n+K_) {

    assert(K <= 21);
    eprintln("n = %ld, K = %d", (long)n, K);
    eprintln("NUM_IPBWT_BYTES = %ld", NUM_IPBWT_BYTES);
    inv_second_size = ((double)1.0)/second_size;
    auto build_ipbwt = [&]() { 
        int64_t startTick, endTick;
         
        startTick = __rdtsc();
        if(__sa == NULL) {
            __sa = (int64_t *)malloc(n * sizeof(int64_t));
            assert(__sa != NULL);
            saisxx(t.c_str(), __sa, (int64_t)n);
         }
        ref_pos_t sa;
        //sa.ls_word = (uint32_t *)malloc(n * sizeof(uint32_t));
        //sa.ms_byte = (uint8_t *)malloc(n * sizeof(uint8_t));
        sa.ls_word = (uint32_t *)calloc(n, sizeof(uint32_t));
        sa.ms_byte = (uint8_t *)calloc(n, sizeof(uint8_t));
        assert(sa.ls_word != NULL);
        assert(sa.ms_byte != NULL);
        for(index_t i = 0; i < n; i++)
        {
            sa.ls_word[i] = __sa[i] & 0xffffffff;
            sa.ms_byte[i] = (__sa[i] >> 32) & 0xff;
        }
        for(index_t i = 0; i < n; i++)
        {
            uint64_t sa_ms_byte = sa.ms_byte[i];
            uint64_t sa_ls_word = sa.ls_word[i];
            uint64_t pos = (sa_ms_byte << 32) + sa_ls_word;
            if(pos != __sa[i])
            {
                eprintln("%ld] pos = %ld, __sa[i] = %ld", (long)i, (long)pos, (long)__sa[i]);
                exit(0);
            }
        }
        free(__sa);
        endTick = __rdtsc();
        eprintln("SA done in %ld cycles.", endTick - startTick);

        {
            startTick = __rdtsc();
            ref_pos_t rnk;
           //rnk.ls_word = (uint32_t *)malloc(n * sizeof(uint32_t));
           //rnk.ms_byte = (uint8_t *)malloc(n * sizeof(uint8_t));
            rnk.ls_word = (uint32_t *)calloc(n, sizeof(uint32_t));
            rnk.ms_byte = (uint8_t *)calloc(n, sizeof(uint8_t));
            assert(rnk.ls_word != NULL);
            assert(rnk.ms_byte != NULL);
            eprintln("rnk allocated, n = %ld", (long)n);
//#pragma omp parallel for
            for(index_t i=0; i<n; i++)
            {
                uint64_t sa_ms_byte = sa.ms_byte[i];
                uint64_t sa_ls_word = sa.ls_word[i];
                uint64_t pos = (sa_ms_byte << 32) + sa_ls_word;
                rnk.ls_word[pos] = i & 0xffffffff;
                rnk.ms_byte[pos] = (i >> 32) & 0xff;
            }
            endTick = __rdtsc();
            eprintln("RNK done in %ld cycles.", endTick - startTick);
            startTick = __rdtsc();
            ipbwt_array = (uint8_t *)aligned_alloc(64, n * NUM_IPBWT_BYTES);
            assert(ipbwt_array != NULL);
        eprintln("ipbwt_array allocated with size %lu", (unsigned long)n * NUM_IPBWT_BYTES);
            assert(K<n);
//#pragma omp parallel for
            for(index_t i=0; i<n; i++) {
                uint64_t pos;
                index_t p = get_ref_pos(sa, i) + K;
                if(p>=n) {
                    p-=n;
                    pos = K-1-p; // dollar_sign_pos
                    //((int64_t *)ipbwt_array)[2 * i] = K-1-p; // dollar_sign_pos
                } else {
                    pos = get_ref_pos(rnk, p) + K;
                    //((int64_t *)ipbwt_array)[2 * i] = get_ref_pos(rnk, p) + K;
                }
                uint64_t str_enc = gen_str_enc(t.c_str() + get_ref_pos(sa, i), K);
                //eprintln("%ld] pos = %lu, str_enc = %lu", i, pos, str_enc);
                //((int64_t *)ipbwt_array)[2 * i + 1] = str_enc;
                pos = (pos & POS_MASK) + (str_enc << NUM_POS_BITS);
                uint16_t chunk = str_enc >> (64 - NUM_POS_BITS);
                //eprintln("%ld] pos = %lu, chunk = %u", i, pos, chunk);
                //fflush(stderr);
                uint64_t *pos_ptr = (uint64_t *)(&(ipbwt_array[i * NUM_IPBWT_BYTES]));
                pos_ptr[0] = pos;
                uint16_t *chunk_ptr = (uint16_t *)(&(ipbwt_array[i * NUM_IPBWT_BYTES])) + 4;
                chunk_ptr[0] = chunk;
            }
            endTick = __rdtsc();
            eprintln("ipbwt_array done in %ld cycles.", endTick - startTick);
            for(index_t i = 1; i < n; i++)
            {
                if(ipbwt(i - 1) >= ipbwt(i))
                {
                    //eprintln("%ld] {%lu, %lu} {%lu, %lu}", i, ((int64_t *)ipbwt_array)[2 * (i - 1) + 1], ((int64_t *)ipbwt_array)[2 * (i - 1)], ((int64_t *)ipbwt_array)[2 * i + 1], ((int64_t *)ipbwt_array)[2 * i]);
                    eprintln("%ld] {%lu, %lu} {%lu, %lu}", (long)i, (unsigned long)ipbwt(i - 1).first, (unsigned long)ipbwt(i - 1). second, (unsigned long)ipbwt(i).first, (unsigned long)ipbwt(i).second);
                    exit(0);
                }
            }
            free(sa.ls_word); free(sa.ms_byte);
            free(rnk.ls_word); free(rnk.ms_byte);
        }
    };


    string ipbwt_filename = ref_seq_filename + ".ipbwt";
    for(const auto &s:{sizeof(index_t), sizeof(kenc_t), (size_t)K, (size_t)NUM_POS_BITS, (size_t)NUM_CHUNK_BITS}) {
        ipbwt_filename += string(".") + to_string(s);
    }
    if(ifstream(ipbwt_filename.c_str()).good()) {
        eprintln("Found existing %s!!", (char*)ipbwt_filename.c_str());
        ipbwt_array = (uint8_t *)aligned_alloc(64, n * NUM_IPBWT_BYTES);
        load(ipbwt_filename, {((char *)ipbwt_array)}, {(n * NUM_IPBWT_BYTES)}); 
        eprintln("Load successful.");
        for(index_t i = 0; i < n; i++)
        {
            if(i > 0)
            {
                if(ipbwt(i - 1) >= ipbwt(i))
                {
                    //eprintln("%ld] {%lu, %lu} {%lu, %lu}", i, ((int64_t *)ipbwt_array)[2 * (i - 1) + 1], ((int64_t *)ipbwt_array)[2 * (i - 1)], ((int64_t *)ipbwt_array)[2 * i + 1], ((int64_t *)ipbwt_array)[2 * i]);
                    eprintln("%ld] {%lu, %lu} {%lu, %lu}", (long)i, (unsigned long)ipbwt(i - 1).first, (unsigned long)ipbwt(i - 1). second, (unsigned long)ipbwt(i).first, (unsigned long)ipbwt(i).second);
                    exit(0);
                }
                assert(ipbwt(i-1) < ipbwt(i));
            }
        }
    } else {
        eprintln("No existing %s. Building...", (char*)ipbwt_filename.c_str());
        build_ipbwt();
        eprintln("%s build done.", (char*)ipbwt_filename.c_str());
        save(ipbwt_filename, {((char *)ipbwt_array)}, {(n * NUM_IPBWT_BYTES)}); 
        eprintln("%s saved.", (char*)ipbwt_filename.c_str());
    }
    eprintln("ipbwt correctness ok.");

    string rmi_filename = ipbwt_filename + "." + to_string(num_rmi_leaf_nodes) + ".rmi_PARAMETERS";
    if(ifstream(rmi_filename.c_str()).good()) {
        eprintln("Found existing %s!!", (char*)rmi_filename.c_str());
        vector<char*> ptrs;vector<size_t> sizes;
        ptrs.push_back((char *)&(L0_PARAMETER0));
        sizes.push_back(sizeof(double));
        ptrs.push_back((char *)&(L0_PARAMETER1));
        sizes.push_back(sizeof(double));
        ptrs.push_back((char *)&(L1_SIZE));
        sizes.push_back(sizeof(int64_t));
        load(rmi_filename, ptrs, sizes);
        L1_PARAMETERS = (double*) _mm_malloc(L1_SIZE * 3 * sizeof(double), 64);
        if (L1_PARAMETERS == NULL)
        {
            eprintln("L1_PARAMETERS == NULL");
            exit(0);;
        }
        ptrs.push_back((char *)L1_PARAMETERS);
        sizes.push_back(L1_SIZE * 3 * sizeof(double));
        load(rmi_filename, ptrs, sizes);
        eprintln("Load successful.");
        eprintln("L1_SIZE = %ld", L1_SIZE);
    } else {
        eprintln("No existing %s. Building train data...", (char*)rmi_filename.c_str());
        double *train_data = (double *)_mm_malloc(n * sizeof(double), 64);
        for(index_t i = 0; i < n; i++)
        {
            ipbwt_t p = ipbwt(i);
            uint64_t second = p.second & POS_MASK;
            uint64_t first = (p.first << (64 - NUM_POS_BITS)) + (p.second >> NUM_POS_BITS);
            train_data[i] = to_floating_point({first, second});
        }
        vector<char*> ptrs;vector<size_t> sizes;
        ptrs.push_back((char *)&n);
        sizes.push_back(sizeof(int64_t));
        ptrs.push_back((char *)train_data);
        sizes.push_back(n * sizeof(double));
        string ipbwt_f64_filename = ipbwt_filename + ".f64";
        save(ipbwt_f64_filename, ptrs, sizes);
        eprintln("Training data saved at %s. Building RMI using Ryan's code. Run this code again after the RMI is built.", (char*)ipbwt_f64_filename.c_str());
        eprintln("executing ./build-rmi.linear_spline.linear.sh");
        //const char *execv_argv[] = {"./build-rmi.linear_spline.linear.sh",
        //                            ipbwt_f64_filename.c_str(),
        //                            rmi_filename.c_str(),
        //                            to_string(num_rmi_leaf_nodes).c_str()};
        //execv_argv[0] = ipbwt_f64_filename.c_str();
        //execv_argv[1] = rmi_filename.c_str();
        //execv_argv[2] = to_string(num_rmi_leaf_nodes);
        int status = execl("./scripts/build-rmi.linear_spline.linear.sh",
              "./scripts/build-rmi.linear_spline.linear.sh",
              ipbwt_f64_filename.c_str(),
              rmi_filename.c_str(),
              to_string(num_rmi_leaf_nodes).c_str(), (char *)NULL);
        eprintln("ERROR! execl failed with err %d, errno = %d: %s", status, errno, strerror( errno));
        exit(0);
    }

    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            m0[i * 4 + j] = i * 5 + j;
            m1[i * 4 + j] = 4 + i * 5 + j;
        }
    }
    for(int i = 4; i < 8; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            m0[i * 4 + j] = i * 5 + j + 24;
            m1[i * 4 + j] = 4 + i * 5 + j + 24;
        }
    }
#if 0
    for(int i = 0; i < 32; i++)
    {
        eprintln("%d] %d %d", i, m0[i], m1[i]);
    }
#endif
    for(int i = 0; i < 9; i++)
    {
        m_one_bits[i] = (1 << i) - 1;
    }
    m_one_bits[8] = 0xff;
    bs_calls = vbs_calls = 0;
    bs_ticks = vbs_ticks = 0;
}

template<typename index_t, typename kenc_t>
IPBWT_RMI<index_t, kenc_t>::~IPBWT_RMI()
{
#ifndef HUGE_PAGE
    free(ipbwt_array);
    _mm_free(L1_PARAMETERS);
#else
    munmap(ipbwt_array, n * NUM_IPBWT_BYTES);
    munmap(L1_PARAMETERS, L1_SIZE * 3 * sizeof(double));
#endif

    //eprintln("bs_calls = %ld, vbs_calls = %ld", bs_calls, vbs_calls);
    //eprintln("bs_ticks = %ld, vbs_ticks = %ld", bs_ticks, vbs_ticks);
}

template<typename index_t, typename kenc_t>
inline index_t IPBWT_RMI<index_t, kenc_t>::get_guess_from_rmi(ipbwt_t ipb_x, index_t *err) const {

    double key = to_floating_point(ipb_x);
    int64_t modelIndex;
    double fpred;
    fpred = std::fma(L0_PARAMETER1, key, L0_PARAMETER0);
    modelIndex = FCLAMP(fpred, L1_SIZE - 1.0);
    fpred = std::fma(L1_PARAMETERS[modelIndex * 3 + 1], key, L1_PARAMETERS[modelIndex * 3]);
    *err = *((uint64_t*) (L1_PARAMETERS + (modelIndex * 3 + 2)));

    //int64_t pos = FCLAMP(fpred, n - 1.0);
    //printf("pos = %ld\n", pos);
    return FCLAMP(fpred, n - 1.0);
}

template<typename index_t, typename kenc_t>
inline index_t IPBWT_RMI<index_t, kenc_t>::last_mile_from_guess
    (ipbwt_t ipb_x, index_t guess, index_t err) const {
    //index_t p = guess;
    index_t low = guess - err;
    if(low < 0) low = 0;
    index_t high = guess + err + 1;
    if(high > n) high = n;

    index_t m = high - low;
    index_t first = low;

    while(m > 0)
    {
        index_t half = m / 2;
        index_t middle = first + half;

        if(ipb_x >= ipbwt(middle))
        {
            first = middle + 1;
            m -= half + 1;
        }
        else
        {
            m = half;
        }
    }
    first--;
    return first;

    // TODO: if the final pos is not in the err range, look for it
}

template<typename index_t, typename kenc_t>
inline pair<index_t, index_t> IPBWT_RMI<index_t, kenc_t>::backward_extend_chunk
(kenc_t str_enc, pair<index_t, index_t> intv) const {
    index_t low = intv.first, high = intv.second;
    ipbwt_t ipb_x_low = {str_enc, low+K};
    ipbwt_t ipb_x_high = {str_enc, high+K};

    int64_t err_low = 0, err_high = 0;
    auto low_guess = get_guess_from_rmi(ipb_x_low, &err_low);

    auto high_guess = get_guess_from_rmi(ipb_x_high, &err_high);

    low = last_mile_from_guess(ipb_x_low, low_guess, err_low);
    high = last_mile_from_guess(ipb_x_high, high_guess, err_high);

    return {low, high};
}

template<typename index_t, typename kenc_t>
inline int64_t IPBWT_RMI<index_t, kenc_t>::get_guess_root_step(double key)
{
    int64_t modelIndex;
    double fpred = std::fma(L0_PARAMETER1, key, L0_PARAMETER0);
    modelIndex = FCLAMP(fpred, L1_SIZE - 1.0);
    return modelIndex;
}

template<typename index_t, typename kenc_t>
inline int64_t IPBWT_RMI<index_t, kenc_t>::get_guess_leaf_step(double key, int64_t modelIndex, int64_t *err)
{
    double fpred = std::fma(L1_PARAMETERS[modelIndex * 3 + 1], key, L1_PARAMETERS[modelIndex * 3]);
    *err = *((uint64_t*) (L1_PARAMETERS + (modelIndex * 3 + 2)));

    int64_t guess = FCLAMP(fpred, n - 1.0);
    return guess;
}

template<typename index_t, typename kenc_t>
inline void IPBWT_RMI<index_t, kenc_t>::last_mile_binary_search_one_step(ipbwt_t ipb_x, int64_t &first, int64_t &m)
{
    //bs_calls++;
    //int64_t startTick = __rdtsc();
    int64_t half = m >> 1;
    int64_t middle = first + half;
    ipbwt_t ipb_middle = ipbwt(middle);
    int64_t cond = (ipb_x.first > ipb_middle.first) | ((ipb_x.first == ipb_middle.first) & (ipb_x.second > ipb_middle.second));
    //bool cond = ipb_x >= ipb_middle;
    int64_t middle_plus_one = middle + 1;
    int64_t m_minus_half_plus_one = m - (half + 1);
    first = cond * middle_plus_one + (1 - cond) * first;
    m = cond * m_minus_half_plus_one + (1 - cond) * half;
    //int64_t endTick = __rdtsc();
    //bs_ticks += (endTick - startTick);
}

template<typename index_t, typename kenc_t>
inline void IPBWT_RMI<index_t, kenc_t>::last_mile_vectorized_search_final_step(ipbwt_t ipb_x, int64_t &first, int64_t &m)
{
    //vbs_calls++;
    __m512i y0 = _mm512_loadu_si512(&(ipbwt_array[NUM_IPBWT_BYTES * first])); 
    __m512i y1 = _mm512_loadu_si512(&(ipbwt_array[NUM_IPBWT_BYTES * (first)]) + 16);
    
    __m512i x_first = _mm512_set1_epi64(ipb_x.first); 
    __m512i x_second = _mm512_set1_epi64(ipb_x.second); 

    __m512i idx0 = _mm512_loadu_si512(m0);
    __m512i idx1 = _mm512_loadu_si512(m1);
    __m512i y_first = _mm512_maskz_permutex2var_epi16(0x11111111, y0, idx1, y1);
    __m512i y_second = _mm512_permutex2var_epi16(y0, idx0, y1);

    __mmask8 mask_gt_first = _mm512_cmpgt_epu64_mask(x_first, y_first);
    __mmask8 mask_eq_first = _mm512_cmpeq_epu64_mask(x_first, y_first);
    __mmask8 mask_gt_second = _mm512_cmpgt_epu64_mask(x_second, y_second);
   
    
    __mmask8 mask_gt = mask_gt_first | (mask_eq_first & mask_gt_second);
    int32_t numgt = _mm_countbits_32(mask_gt & m_one_bits[m]);
    first = first + numgt;
    m = 0;
}

template<typename index_t, typename kenc_t>
inline bool IPBWT_RMI<index_t, kenc_t>::process_query_one_step(BatchMetadata &meta, index_t &p_low, index_t &p_high)
{
    if(meta.state == GUESS_RMI_ROOT)
    {
        meta.modelIndex[0] = get_guess_root_step(meta.key[0]);
        meta.modelIndex[1] = get_guess_root_step(meta.key[1]);
        meta.state = GUESS_RMI_LEAF;
        _mm_prefetch((const char *)(&L1_PARAMETERS[meta.modelIndex[0] * 3]), _MM_HINT_T0);
        _mm_prefetch((const char *)(&L1_PARAMETERS[meta.modelIndex[0] * 3 + 2]), _MM_HINT_T0);
        _mm_prefetch((const char *)(&L1_PARAMETERS[meta.modelIndex[1] * 3]), _MM_HINT_T0);
        _mm_prefetch((const char *)(&L1_PARAMETERS[meta.modelIndex[1] * 3 + 2]), _MM_HINT_T0);
    }
    else if(meta.state == GUESS_RMI_LEAF)
    {
        int64_t err[2];
        int64_t guess[2];
        int64_t last[2];
        err[0] = err[1] = 0;
        guess[0] = get_guess_leaf_step(meta.key[0], meta.modelIndex[0], &err[0]);
        guess[1] = get_guess_leaf_step(meta.key[1], meta.modelIndex[1], &err[1]);
        meta.first[0] = guess[0] - err[0];
        if(meta.first[0] < 0) meta.first[0] = 0;
        last[0] = guess[0] + err[0] + 1;
        if(last[0] > n) last[0] = n;
        meta.m[0] = last[0] - meta.first[0];
        meta.first[1] = guess[1] - err[1];
        if(meta.first[1] < 0) meta.first[1] = 0;
        last[1] = guess[1] + err[1] + 1;
        if(last[1] > n) last[1] = n;
        meta.m[1] = last[1] - meta.first[1];
        meta.state = LAST_MILE;
        if(meta.m[0] > 8)
        {
            _mm_prefetch((const char *)(&ipbwt_array[(meta.first[0] + (meta.m[0] >> 1)) * NUM_IPBWT_BYTES]), _MM_HINT_T0);
            _mm_prefetch((const char *)(&ipbwt_array[(meta.first[0] + (meta.m[0] >> 1)) * NUM_IPBWT_BYTES]) + 8, _MM_HINT_T0);
        }
        else
        {
            _mm_prefetch((const char *)(&ipbwt_array[meta.first[0] * NUM_IPBWT_BYTES]), _MM_HINT_T0);
            _mm_prefetch((const char *)(&ipbwt_array[meta.first[0] * NUM_IPBWT_BYTES]) + 40, _MM_HINT_T0);
            _mm_prefetch((const char *)(&ipbwt_array[meta.first[0] * NUM_IPBWT_BYTES]) + 79, _MM_HINT_T0);
        }
        if(meta.m[1] > 8)
        {
            _mm_prefetch((const char *)(&ipbwt_array[(meta.first[1] + (meta.m[1] >> 1)) * NUM_IPBWT_BYTES]), _MM_HINT_T0);
            _mm_prefetch((const char *)(&ipbwt_array[(meta.first[1] + (meta.m[1] >> 1)) * NUM_IPBWT_BYTES]) + 8, _MM_HINT_T0);
        }
        else
        {
            _mm_prefetch((const char *)(&ipbwt_array[meta.first[1] * NUM_IPBWT_BYTES]), _MM_HINT_T0);
            _mm_prefetch((const char *)(&ipbwt_array[meta.first[1] * NUM_IPBWT_BYTES]) + 40, _MM_HINT_T0);
            _mm_prefetch((const char *)(&ipbwt_array[meta.first[1] * NUM_IPBWT_BYTES]) + 79, _MM_HINT_T0);
        }
    }
    else
    {
        if(meta.m[0] > 8)
        {
            last_mile_binary_search_one_step(meta.ipb_x[0], meta.first[0], meta.m[0]);
            if(meta.m[0] > 8)
            {
                _mm_prefetch((const char *)(&ipbwt_array[(meta.first[0] + (meta.m[0] >> 1)) * NUM_IPBWT_BYTES]), _MM_HINT_T0);
                _mm_prefetch((const char *)(&ipbwt_array[(meta.first[0] + (meta.m[0] >> 1)) * NUM_IPBWT_BYTES]) + 8, _MM_HINT_T0);
            }
            else
            {
                _mm_prefetch((const char *)(&ipbwt_array[meta.first[0] * NUM_IPBWT_BYTES]), _MM_HINT_T0);
                _mm_prefetch((const char *)(&ipbwt_array[meta.first[0] * NUM_IPBWT_BYTES]) + 40, _MM_HINT_T0);
                _mm_prefetch((const char *)(&ipbwt_array[meta.first[0] * NUM_IPBWT_BYTES]) + 79, _MM_HINT_T0);
            }
        }
        else
        {
            last_mile_vectorized_search_final_step(meta.ipb_x[0], meta.first[0], meta.m[0]);
        }
        if(meta.m[1] > 8)
        {
            last_mile_binary_search_one_step(meta.ipb_x[1], meta.first[1], meta.m[1]);
            if(meta.m[1] > 8)
            {
                _mm_prefetch((const char *)(&ipbwt_array[(meta.first[1] + (meta.m[1] >> 1)) * NUM_IPBWT_BYTES]), _MM_HINT_T0);
                _mm_prefetch((const char *)(&ipbwt_array[(meta.first[1] + (meta.m[1] >> 1)) * NUM_IPBWT_BYTES]) + 8, _MM_HINT_T0);
            }
            else
            {
                _mm_prefetch((const char *)(&ipbwt_array[meta.first[1] * NUM_IPBWT_BYTES]), _MM_HINT_T0);
                _mm_prefetch((const char *)(&ipbwt_array[meta.first[1] * NUM_IPBWT_BYTES]) + 40, _MM_HINT_T0);
                _mm_prefetch((const char *)(&ipbwt_array[meta.first[1] * NUM_IPBWT_BYTES]) + 79, _MM_HINT_T0);
            }
        }
        else
        {
            last_mile_vectorized_search_final_step(meta.ipb_x[1], meta.first[1], meta.m[1]);
        }
        if(meta.m[0] == 0 && meta.m[1] == 0)
        {
            p_low = meta.first[0];
            p_high = meta.first[1];
            return 0;
        }
    }
    return 1;
}


template<typename index_t, typename kenc_t>
void IPBWT_RMI<index_t, kenc_t>::backward_extend_chunk_batched
    (kenc_t* str_enc_list, int64_t qs_size, index_t *intv_all) {

    //eprintln("backward_extend_chunk_batched");
    BatchMetadata metadata[UNROLL];
    {
        int64_t init_batch_size = min(UNROLL, qs_size);
        int64_t next_query;
        for(next_query = 0; next_query < init_batch_size; next_query++)
        {
            BatchMetadata m;
            m.qid = next_query;
            m.state = GUESS_RMI_ROOT;
            m.ipb_x[0] = ipbwt(str_enc_list[next_query], intv_all[next_query * 2]+K); 
            m.ipb_x[1] = ipbwt(str_enc_list[next_query], intv_all[next_query * 2 + 1]+K); 
            m.key[0] = to_floating_point({str_enc_list[next_query], intv_all[next_query * 2]+K});
            m.key[1] = to_floating_point({str_enc_list[next_query], intv_all[next_query * 2 + 1]+K});
            metadata[next_query] = m;
        }
        int64_t num_active = next_query;

        while(num_active > 0)
        {
            for(int64_t i = 0; i < num_active; i++)
            {
                index_t pos_low, pos_high;
                int status = process_query_one_step(metadata[i], pos_low, pos_high);
                if(!status)
                {
                    intv_all[2 * metadata[i].qid] = pos_low;
                    intv_all[2 * metadata[i].qid + 1] = pos_high;
                    if(next_query < qs_size)
                    {
                        // Get next query
                        BatchMetadata m;
                        m.qid = next_query;
                        m.state = GUESS_RMI_ROOT;
                        m.ipb_x[0] = ipbwt(str_enc_list[next_query], intv_all[next_query * 2]+K); 
                        m.ipb_x[1] = ipbwt(str_enc_list[next_query], intv_all[next_query * 2 + 1]+K); 
                        m.key[0] = to_floating_point({str_enc_list[next_query], intv_all[next_query * 2]+K});
                        m.key[1] = to_floating_point({str_enc_list[next_query], intv_all[next_query * 2 + 1]+K});
                        metadata[i] = m;
                        next_query++;
                    }
                    else
                    {
                        // No more new queries
                        metadata[i] = metadata[num_active - 1];
                        num_active--;
                    }
                }
            }
        }
    }
}

template<typename index_t, typename kenc_t>
void IPBWT_RMI<index_t, kenc_t>::backward_extend_multi_chunk_batched
    (Info* str_enc_list, int64_t qs_size, index_t *intv_all) {

    //eprintln("backward_extend_chunk_batched");
    BatchMetadata metadata[UNROLL];
    {
        int64_t init_batch_size = min(UNROLL, qs_size);
        int64_t next_query;
        for(next_query = 0; next_query < init_batch_size; next_query++)
        {
            BatchMetadata m;
            m.qid = next_query;
            m.state = GUESS_RMI_ROOT;
	    uint64_t enc_str = str_enc_list[next_query].get_enc_str();
            m.ipb_x[0] = ipbwt(enc_str, intv_all[next_query * 2]+K); 
            m.ipb_x[1] = ipbwt(enc_str, intv_all[next_query * 2 + 1]+K); 
            m.key[0] = to_floating_point({enc_str, intv_all[next_query * 2]+K});
            m.key[1] = to_floating_point({enc_str, intv_all[next_query * 2 + 1]+K});
            metadata[next_query] = m;
        }
        int64_t num_active = next_query;

        while(num_active > 0)
        {
            for(int64_t i = 0; i < num_active; i++)
            {
                index_t pos_low, pos_high;
                int status = process_query_one_step(metadata[i], pos_low, pos_high);
                if(!status)
                {
                    intv_all[2 * metadata[i].qid] = pos_low;
                    intv_all[2 * metadata[i].qid + 1] = pos_high;
                    if(next_query < qs_size)
                    {
                        // Get next query
                        BatchMetadata m;
                        m.qid = next_query;
                        m.state = GUESS_RMI_ROOT;
			uint64_t enc_str = str_enc_list[next_query].get_enc_str();
                        m.ipb_x[0] = ipbwt(enc_str, intv_all[next_query * 2]+K); 
                        m.ipb_x[1] = ipbwt(enc_str, intv_all[next_query * 2 + 1]+K); 
                        m.key[0] = to_floating_point({enc_str, intv_all[next_query * 2]+K});
                        m.key[1] = to_floating_point({enc_str, intv_all[next_query * 2 + 1]+K});
                        metadata[i] = m;
                        next_query++;
                    }
                    else
                    {
                        // No more new queries
                        metadata[i] = metadata[num_active - 1];
                        num_active--;
                    }
                }
            }
        }
    }
}






#endif
