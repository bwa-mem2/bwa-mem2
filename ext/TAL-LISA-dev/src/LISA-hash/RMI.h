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

Authors: Sanchit Misra <sanchit.misra@intel.com>; Saurabh Kalikar <saurabh.kalikar@intel.com>; 
*****************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdint.h>
#include <limits.h>
#include <immintrin.h>
#include <fstream>
#include <cmath>
#include <sys/mman.h>

#define BATCH_SIZE 32

//typedef double rmi_key_t;
//typedef uint64_t rmi_key_t;

enum query_state
{
    GUESS_RMI_ROOT,
    GUESS_RMI_LEAF,
    LAST_MILE
};



template<typename rmi_key_t>
class RMI
{
    public:
    RMI(char *prefix);
    ~RMI();


	typedef struct batchMetadata
	{
	    int64_t qid;
	    query_state state;
	    rmi_key_t key;
	    int64_t modelIndex;
	    int64_t first;
	    int64_t m;
	}BatchMetadata;

    void get_random_keys(int64_t nq, rmi_key_t *key_array, int64_t *orig_pos_array);
    rmi_key_t get_element(int64_t index);
    uint64_t get_guess(rmi_key_t key, int64_t *err);
    uint64_t last_mile_search(rmi_key_t key, int64_t guess, size_t err);
    int64_t lookup(rmi_key_t key);
    void lookup_batched(rmi_key_t *key_array, int64_t num_queries, int64_t *pos_array);
    void print_stats();

    private:
    int64_t n;
    rmi_key_t *sorted_array;
    //double L0_PARAMETER0 = 0.0;
    //double L0_PARAMETER1 = 0.0000009536740890325746;
    //int64_t L1_SIZE = 4194304;
    double L0_PARAMETER0 = 0.0;
    double L0_PARAMETER1 = 0.0;
    int64_t L1_SIZE = 0;
    double *L1_PARAMETERS;
    bool load_sorted_array(char *prefix);
    bool load_rmi(char *prefix);
    inline int64_t FCLAMP(double inp, double bound);
    inline int64_t get_guess_root_step(rmi_key_t key);
    inline int64_t get_guess_leaf_step(rmi_key_t key, int64_t modelIndex, int64_t *err);
    inline void last_mile_search_one_step(rmi_key_t key, int64_t &first, int64_t &m);
    inline void last_mile_search_vectorized_step(rmi_key_t key, int64_t &first, int64_t &m);
    int process_query_one_step(BatchMetadata &meta, int64_t &pos);
#ifdef STATS
    double total_log_err = 0;
    double total_log_gap = 0;
    int64_t max_err = 0;
    int64_t max_gap = 0;
    int64_t err_hist[20];
    int64_t nq = 0;
#endif
};


template<typename rmi_key_t>
RMI<rmi_key_t>::RMI(char *prefix)
{
    load_sorted_array(prefix);
    load_rmi(prefix);
#ifdef STATS
    memset(err_hist, 0, 20*8);
#endif
}

template<typename rmi_key_t>
RMI<rmi_key_t>::~RMI()
{
    free(sorted_array);
    _mm_free(L1_PARAMETERS);
}

template<typename rmi_key_t>
void RMI<rmi_key_t>::get_random_keys(int64_t nq, rmi_key_t *key_array, int64_t *orig_pos_array)
{
    srand(0);
    for(int64_t i = 0; i < nq; i++)
    {
        double rand_num = rand() * 1.0 /(1L + RAND_MAX);
        int64_t query_id = rand_num * nq;
        key_array[i] = sorted_array[query_id];
        orig_pos_array[i] = query_id;
    }
}

template<typename rmi_key_t>
rmi_key_t RMI<rmi_key_t>::get_element(int64_t index)
{
    if(index < 0 || index >= n)
    {
        printf("index is out of bounds, index = %ld\n", index);
        exit(0);
    }
    return sorted_array[index];
}

template<typename rmi_key_t>
bool RMI<rmi_key_t>::load_sorted_array(char *prefix)
{
    char filename[PATH_MAX];
    strcpy(filename, prefix);

	//suffix

#ifdef UINT64
    strcat(filename, ".uint64");
#endif
#ifdef F64
    strcat(filename, ".f64");
#endif
    std::ifstream infile(filename, std::ios::in | std::ios::binary);
    if (!infile.good())
    {
        printf("%s file not found\n", filename);
        exit(0);
    }
    infile.read((char *)&(this->n), sizeof(uint64_t));
    fprintf(stderr, "n = %ld\n", n);
#ifdef MMAP
    sorted_array = (rmi_key_t*) mmap(NULL, n * sizeof(rmi_key_t), PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB, -1, 0);
#else
    sorted_array = (rmi_key_t*) malloc(n * sizeof(rmi_key_t));
#endif
    if (sorted_array == NULL) return false;

    infile.read((char*)sorted_array, n * sizeof(rmi_key_t));



    if (!infile.good()) return false;
    return true;
}

template<typename rmi_key_t>
bool RMI<rmi_key_t>::load_rmi(char *prefix)
{
    char filename[PATH_MAX];
    strcpy(filename, prefix);
    strcat(filename, ".rmi_PARAMETERS");
    std::ifstream infile(filename, std::ios::in | std::ios::binary);
    if (!infile.good())
    {
        printf("%s file not found\n", filename);
        exit(0);
    }
    infile.read((char *)&(this->L0_PARAMETER0), sizeof(double));
    infile.read((char *)&(this->L0_PARAMETER1), sizeof(double));
    infile.read((char *)&(this->L1_SIZE), sizeof(int64_t));
    fprintf(stderr, "L0_PARAMETER0 = %E, L0_PARAMETER1 = %E, L1_SIZE = %ld\n", 
            L0_PARAMETER0, L0_PARAMETER1, L1_SIZE);
    L1_PARAMETERS = (double*) _mm_malloc(L1_SIZE * 3 * sizeof(double), 64);
    if (L1_PARAMETERS == NULL) return false;
    infile.read((char*)L1_PARAMETERS, L1_SIZE * 3 * sizeof(double));
    if (!infile.good()) return false;
    return true;
}

template<typename rmi_key_t>
uint64_t RMI<rmi_key_t>::get_guess(rmi_key_t key, int64_t *err)
{
    //printf("entering get_guess n = %ld\n", n);
    size_t modelIndex;
    double fpred;
    fpred = std::fma(L0_PARAMETER1, key, L0_PARAMETER0);
    modelIndex = FCLAMP(fpred, L1_SIZE - 1.0);
    //printf("modelIndex = %lu, n = %ld\n", modelIndex, n);
    fpred = std::fma(L1_PARAMETERS[modelIndex * 3 + 1], key, L1_PARAMETERS[modelIndex * 3]);
    //printf("%E, %E, %E\n", L1_PARAMETERS[modelIndex * 3], L1_PARAMETERS[modelIndex * 3 + 1], key);
    //printf("fpred = %lf, n = %ld\n", fpred, n);
    *err = *((uint64_t*) (L1_PARAMETERS + (modelIndex * 3 + 2)));

    int64_t pos = FCLAMP(fpred, n - 1.0);
    //printf("pos = %ld\n", pos);
    return FCLAMP(fpred, n - 1.0);
}


template<typename rmi_key_t>
uint64_t RMI<rmi_key_t>::last_mile_search(rmi_key_t key, int64_t guess, size_t err)
{
    int64_t first = guess - err;
    if(first < 0) first = 0;
    int64_t last = guess + err + 1;
    if(last > n) last = n;
    int64_t m = last - first;

    while(m > 1)
    {
        int64_t half = m / 2;
        int64_t middle = first + half;
        if(key >= sorted_array[middle])
        {
            first = middle;
            m -= half;
        }
        else
        {
            m = half;
        }
    }
    return first;
}

template<typename rmi_key_t>
int64_t RMI<rmi_key_t>::lookup(rmi_key_t key)
{
    int64_t err;
    //printf("entering lookup n = %ld\n", n);
    int64_t guess = get_guess(key, &err);
    //printf("guess = %ld, err = %ld\n", guess, err);
    //if(guess != 0) printf("guess non zero\n");
    int64_t pos = last_mile_search(key, guess, err);
    if(key != sorted_array[pos])
    {
        printf("key = %lld, sorted_array[%lld] = %lld, err = %ld, gap = %ld\n", key, pos, sorted_array[pos], err, labs(pos - guess));
        exit(0);
    }
#ifdef STATS
    nq++;
    int64_t gap = labs(guess-pos);
    int64_t log2_err = log2(1.0 * err);
    err_hist[log2_err]++;
    if(max_err < err) max_err = err;
    if(max_gap < gap) max_gap = gap;
    total_log_err += log2(1.0 + err);
    total_log_gap += log2(1.0 + gap);
#endif
    return pos;
}

template<typename rmi_key_t>
inline int64_t RMI<rmi_key_t>::get_guess_root_step(rmi_key_t key)
{
    int64_t modelIndex;
    double fpred = std::fma(L0_PARAMETER1, key, L0_PARAMETER0);
    modelIndex = FCLAMP(fpred, L1_SIZE - 1.0);
    return modelIndex;
}

template<typename rmi_key_t>
inline int64_t RMI<rmi_key_t>::get_guess_leaf_step(rmi_key_t key, int64_t modelIndex, int64_t *err)
{
    double fpred = std::fma(L1_PARAMETERS[modelIndex * 3 + 1], key, L1_PARAMETERS[modelIndex * 3]);
    *err = *((uint64_t*) (L1_PARAMETERS + (modelIndex * 3 + 2)));

    int64_t guess = FCLAMP(fpred, n - 1.0);
    return guess;
}

template<typename rmi_key_t>
inline void RMI<rmi_key_t>::last_mile_search_one_step(rmi_key_t key, int64_t &first, int64_t &m)
{
    int64_t half = m >> 1;
    int64_t middle = first + half;
    int64_t cond = (key >= sorted_array[middle]);
    first = middle * cond + first * (1 - cond);
    m = (m - half) * cond + half * (1 - cond);
}

#if VECTORIZE && __AVX512BW__
template<typename rmi_key_t>
inline void RMI<rmi_key_t>::last_mile_search_vectorized_step(rmi_key_t key, int64_t &first, int64_t &m)
{
    //    printf("Vectorized\n");

#ifdef F64
    __m512d key_vec = _mm512_set1_pd(key);
    __m512d element_vec = _mm512_loadu_pd(sorted_array + first);
    __mmask8 mask_lt = _mm512_cmplt_pd_mask(element_vec, key_vec);
#endif
#ifdef UINT64
    __m512i key_vec = _mm512_set1_epi64(key);
//    __m512i element_vec = _mm512_loadu_epi64(sorted_array + first);
    __m512i element_vec = _mm512_loadu_si512(sorted_array + first);
    __mmask8 mask_lt = _mm512_cmplt_epi64_mask(element_vec, key_vec);
#endif
    //int32_t numlt = _mm_countbits_32(mask_lt & ((1 << m) - 1));
    int32_t numlt = _mm_popcnt_u32(mask_lt & ((1 << m) - 1));
    first = first + numlt;
    m = 0;
}

template<typename rmi_key_t>
int RMI<rmi_key_t>::process_query_one_step(BatchMetadata &bm, int64_t &pos)
{
    if(bm.state == GUESS_RMI_ROOT)
    {
        bm.modelIndex = get_guess_root_step(bm.key);
        bm.state = GUESS_RMI_LEAF;
        _mm_prefetch((const char *)(&L1_PARAMETERS[bm.modelIndex * 3]), _MM_HINT_T0);
        _mm_prefetch((const char *)(&L1_PARAMETERS[bm.modelIndex * 3 + 2]), _MM_HINT_T0);
    }
    else if(bm.state == GUESS_RMI_LEAF)
    {
        int64_t err;
        int64_t guess = get_guess_leaf_step(bm.key, bm.modelIndex, &err);
        bm.first = guess - err;
        if(bm.first < 0) bm.first = 0;
        int64_t last = guess + err + 1;
        if(last > n) last = n;
        bm.m = last - bm.first;
        bm.state = LAST_MILE;
        int64_t middle = bm.m >> 1;
        int64_t cond = (bm.m > 8);
        int64_t pf0 = bm.first + cond * middle;
        int64_t pf1 = bm.first + (1 - cond) * 7;
        _mm_prefetch((const char *)(&sorted_array[pf0]), _MM_HINT_T0);
        _mm_prefetch((const char *)(&sorted_array[pf1]), _MM_HINT_T0);
    }
    else
    {
        if(bm.m > 8)
        {
            last_mile_search_one_step(bm.key, bm.first, bm.m);
            int64_t middle = bm.m >> 1;
            int64_t cond = (bm.m > 8);
            int64_t pf0 = bm.first + cond * middle;
            int64_t pf1 = bm.first + (1 - cond) * 7;
            _mm_prefetch((const char *)(&sorted_array[pf0]), _MM_HINT_T0);
            _mm_prefetch((const char *)(&sorted_array[pf1]), _MM_HINT_T0);
        }
#if 1
        else
        {
            last_mile_search_vectorized_step(bm.key, bm.first, bm.m);
            pos = bm.first;
	    if(sorted_array[pos] != bm.key)
		pos = -1;
            return 0;
        }
#else
        if(bm.m == 1)
        {
            pos = bm.first;
            return 0;
        }
#endif
    }
    return 1;
}

#else

template<typename rmi_key_t>
int RMI<rmi_key_t>::process_query_one_step(BatchMetadata &bm, int64_t &pos)
{
    if(bm.state == GUESS_RMI_ROOT)
    {
        bm.modelIndex = get_guess_root_step(bm.key);
        bm.state = GUESS_RMI_LEAF;
        _mm_prefetch((const char *)(&L1_PARAMETERS[bm.modelIndex * 3]), _MM_HINT_T0);
        _mm_prefetch((const char *)(&L1_PARAMETERS[bm.modelIndex * 3 + 2]), _MM_HINT_T0);
    }
    else if(bm.state == GUESS_RMI_LEAF)
    {
        int64_t err;
        int64_t guess = get_guess_leaf_step(bm.key, bm.modelIndex, &err);
        bm.first = guess - err;
        if(bm.first < 0) bm.first = 0;
        int64_t last = guess + err + 1;
        if(last > n) last = n;
        bm.m = last - bm.first;
        bm.state = LAST_MILE;
        int64_t middle = bm.m >> 1;
        _mm_prefetch((const char *)(&sorted_array[bm.first + middle]), _MM_HINT_T0);
    }
    else
    {
        if(bm.m > 1)
        {
            last_mile_search_one_step(bm.key, bm.first, bm.m);
            int64_t middle = bm.m >> 1;
            _mm_prefetch((const char *)(&sorted_array[bm.first + middle]), _MM_HINT_T0);
        }
        if(bm.m == 1)
        {
            pos = bm.first;
	    if(sorted_array[pos] != bm.key)
		pos = -1;
            return 0;
        }
    }
    return 1;
}
#endif

#if 0
int RMI<rmi_key_t>::process_query_one_step(BatchMetadata &bm, int64_t &pos)
{
    if(bm.state == GUESS_RMI_ROOT)
    {
        bm.modelIndex = get_guess_root_step(bm.key);
        bm.state = GUESS_RMI_LEAF;
        _mm_prefetch((const char *)(&L1_PARAMETERS[bm.modelIndex * 3]), _MM_HINT_T0);
        _mm_prefetch((const char *)(&L1_PARAMETERS[bm.modelIndex * 3 + 2]), _MM_HINT_T0);
    }
    else if(bm.state == GUESS_RMI_LEAF)
    {
        int64_t err;
        int64_t guess = get_guess_leaf_step(bm.key, bm.modelIndex, &err);
        bm.first = guess - err;
        if(bm.first < 0) bm.first = 0;
        int64_t last = guess + err + 1;
        if(last > n) last = n;
        bm.m = last - bm.first;
        bm.state = LAST_MILE;
        if(bm.m > 8)
        {
            _mm_prefetch((const char *)(&sorted_array[bm.first + (bm.m >> 1)]), _MM_HINT_T0);
        }
        else
        {
            _mm_prefetch((const char *)(&sorted_array[bm.first]), _MM_HINT_T0);
            _mm_prefetch((const char *)(&sorted_array[bm.first + 7]), _MM_HINT_T0);
        }
    }
    else
    {
        if(bm.m > 8)
        {
            last_mile_search_one_step(bm.key, bm.first, bm.m);
            if(bm.m > 8)
            {
                _mm_prefetch((const char *)(&sorted_array[bm.first + (bm.m >> 1)]), _MM_HINT_T0);
            }
            else
            {
                _mm_prefetch((const char *)(&sorted_array[bm.first]), _MM_HINT_T0);
                _mm_prefetch((const char *)(&sorted_array[bm.first + 7]), _MM_HINT_T0);
            }
        }
#if 1
        else
        {
            last_mile_search_vectorized_step(bm.key, bm.first, bm.m);
            pos = bm.first;
            return 0;
        }
#else
        if(bm.m == 1)
        {
            pos = bm.first;
            return 0;
        }
#endif
    }
    return 1;
}
#endif

template<typename rmi_key_t>
void RMI<rmi_key_t>::lookup_batched(rmi_key_t *key_array, int64_t num_queries, int64_t *pos_array)
{
    int64_t next_query;

    BatchMetadata metadata[BATCH_SIZE];

    for(next_query = 0; next_query < BATCH_SIZE && next_query < num_queries; next_query++)
    {
        BatchMetadata bm;
        bm.qid = next_query;
        bm.state = GUESS_RMI_ROOT;
        bm.key = key_array[next_query];
        metadata[next_query] = bm;
    }
    int64_t num_active = next_query;

    while(num_active > 0)
    {
        for(int64_t i = 0; i < num_active; i++)
        {
            int64_t pos;
            int status = process_query_one_step(metadata[i], pos);
            //printf("%ld] qid = %ld, status = %d, pos = %ld\n", i, bm.qid, status, pos);
            if(!status)
            {
                BatchMetadata bm = metadata[i];
                //if(bm.key != sorted_array[pos])
                //{
                //    printf("%ld] key = %E, sorted_array[%d] = %E\n", bm.qid, bm.key, pos, sorted_array[pos]);
                //    exit(0);
                //}
                //printf("correct\n");
                pos_array[bm.qid] = pos;
                //printf("next_query = %ld, num_queries = %ld\n", next_query, num_queries);
                if(next_query < num_queries)
                {
                    // Get next query
                    BatchMetadata bm;
                    bm.qid = next_query;
                    bm.state = GUESS_RMI_ROOT;
                    bm.key = key_array[next_query];
                    metadata[i] = bm;
                    next_query++;
                }
                else
                {
                    // No more new queries
                    metadata[i] = metadata[num_active - 1];
                    num_active--;
                }
                //printf("After getting next_query\n");
            }
        }
    }
}



template<typename rmi_key_t>
inline int64_t RMI<rmi_key_t>::FCLAMP(double inp, double bound) {
  if (inp < 0.0) return 0;
  return (inp > bound ? bound : (size_t)inp);
}

template<typename rmi_key_t>
void RMI<rmi_key_t>::print_stats()
{
#ifdef STATS
    printf("avg_log2_err = %lf, avg_log2_gap = %lf\n", total_log_err / nq, total_log_gap / nq);
    printf("max_err = %ld, max_gap = %ld\n", max_err, max_gap);
    for(int i = 0; i < 20; i++)
    {
        printf("%ld] log2_err freq = %ld\n", i, err_hist[i]);
    }
#endif
}
