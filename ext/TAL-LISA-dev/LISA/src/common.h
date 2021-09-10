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
#ifndef COMMON_H
#define COMMON_H

#include <vector>
#include <cstring>
#include <string>
#include <tuple>
#include <algorithm>
#include <cassert>
#include <numeric>
#include <utility>
#include <cstdio>
#include <type_traits>
#include <fstream>
#include <array>
#include <initializer_list>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <climits>
#include <map>
#include <future>
#include <omp.h>
#include <array>
using namespace std;

#ifndef __rdtsc 
#ifdef _rdtsc
#define __rdtsc _rdtsc
#else
#define __rdtsc __builtin_ia32_rdtsc 
#endif
#endif

#ifdef __lg
#undef __lg
#endif

#ifndef _MM_HINT_NT
#define _MM_HINT_NT _MM_HINT_NTA
#endif 

template<typename T>
inline constexpr unsigned long __lg(T n) {
    return sizeof(uint64_t) * __CHAR_BIT__  - 1 - __builtin_clzll(n);
}

#define eprintln(...) do{\
    fprintf(stderr,__VA_ARGS__);\
    fprintf(stderr,"\n");\
}while(0)

#define error_quit(...) do{\
    eprintln(__VA_ARGS__);\
    exit(1);\
}while(0)

const string dna = "ACGT";
constexpr int dna_ord(const char &a) {
#ifdef NO_DNA_ORD 
    __builtin_unreachable();
    // assert(0 && "dna_ord is not supported");
#else
    return __lg(a-'A'+2)-1; // "ACGT" -> 0123
#endif 
}

#ifdef _64BIT 
    typedef int64_t index_t;
#else 
    typedef uint32_t index_t;
#endif


class SMEM_out {
    public:
	    int id, q_l, q_r; 
	    index_t ref_l, ref_r;
	    SMEM_out(int _id, int _q_l, int _q_r, index_t _ref_l, index_t _ref_r){
		    id = _id;
		    q_l = _q_l;
		    q_r = _q_r;
		    ref_l = _ref_l;
		    ref_r = _ref_r;	
	    }

};


class vector_based_output {
	public:
		int id;
		vector<pair<int, int>> qPos;
		vector<pair<index_t, index_t>> refPos;
		vector_based_output(int a){ id  = a;

		}
};


class Output {
	public:
		int id;
		SMEM_out* smem;
		Output(int a){ id  = a;
		}
};

 
struct Info {
    // TODO use big int?
    const char* p;
    int l, r; 
    uint64_t id;
    pair<index_t, index_t> intv;




	void set(int a, int b, index_t c, index_t d){
		l = a; r = b; intv = make_pair(c,d);
	}
	void print(){
//		printf(" %d %d %lld %lld %lld %d %d ", l, r, intv.first, intv.second, intv.second - intv.first, numPrevSuccChk, treeShrinkLength);
	}


    uint64_t get_enc_str(){
	uint64_t nxt_ext = 0;
	int i = l - 21; //K;
	//TODO: hard coded K
	nxt_ext = (nxt_ext<<2) | dna_ord(p[i++]);
	nxt_ext = (nxt_ext<<2) | dna_ord(p[i++]);
	nxt_ext = (nxt_ext<<2) | dna_ord(p[i++]);
	nxt_ext = (nxt_ext<<2) | dna_ord(p[i++]);
	nxt_ext = (nxt_ext<<2) | dna_ord(p[i++]);
	
	nxt_ext = (nxt_ext<<2) | dna_ord(p[i++]);
	nxt_ext = (nxt_ext<<2) | dna_ord(p[i++]);
	nxt_ext = (nxt_ext<<2) | dna_ord(p[i++]);
	nxt_ext = (nxt_ext<<2) | dna_ord(p[i++]);
	nxt_ext = (nxt_ext<<2) | dna_ord(p[i++]);
	
	nxt_ext = (nxt_ext<<2) | dna_ord(p[i++]);
	nxt_ext = (nxt_ext<<2) | dna_ord(p[i++]);
	nxt_ext = (nxt_ext<<2) | dna_ord(p[i++]);
	nxt_ext = (nxt_ext<<2) | dna_ord(p[i++]);
	nxt_ext = (nxt_ext<<2) | dna_ord(p[i++]);
	
	nxt_ext = (nxt_ext<<2) | dna_ord(p[i++]);
	nxt_ext = (nxt_ext<<2) | dna_ord(p[i++]);
	nxt_ext = (nxt_ext<<2) | dna_ord(p[i++]);
	nxt_ext = (nxt_ext<<2) | dna_ord(p[i++]);
	nxt_ext = (nxt_ext<<2) | dna_ord(p[i++]);
	nxt_ext = (nxt_ext<<2) | dna_ord(p[i++]);

	return nxt_ext;
    }

    void printRead(){
	printf(">\n");
	for(int i = 0; i < r; i++){
		printf("%c",p[i]);
	}
	printf("\n");
	
	
    }   

};
#endif
