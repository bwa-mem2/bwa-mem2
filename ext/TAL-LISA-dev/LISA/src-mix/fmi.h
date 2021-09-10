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
#ifndef FMI_H
#define FMI_H

#include "sais.h"
#include <sys/mman.h>


template<typename index_t>
class FMI {
    public: 
        FMI(){}
        FMI(const string &t, index_t *sa, string _bases, string parent_filename); 
        FMI(const string &t, index_t t_size, index_t *sa, string _bases, string parent_filename); 
        struct Interval{ index_t low, high; }; // left-inclusive
        Interval backward_extend(Interval intv, char a) const;
	~FMI();
    // private:
        static constexpr int INDEX_T_BITS = sizeof(index_t)*__CHAR_BIT__;
        string bases;
        int ALP_SIZE;
        index_t n;
        int64_t m;
        index_t *cnt;
        index_t *occb;
        inline index_t occ(int i, index_t j) const {
            static constexpr int shift = __lg(INDEX_T_BITS), msk = INDEX_T_BITS-1;
            const index_t* d = occb + ((j>>shift)<<3) + ((i-1)<<1);
            return d[1] + __builtin_popcountll(((typename make_unsigned<index_t>::type)d[0])<<(msk&~j));
        }
        void load(string filename);
        void save(string filename) const;
};

template<typename index_t>
void FMI<index_t>::load(string filename) {
    ifstream instream(filename.c_str(), ifstream::binary);
    instream.seekg(0);

    cnt = new index_t[ALP_SIZE]();
    eprintln("MEM-SIZE: %lld", (long long)sizeof(index_t)*ALP_SIZE);
    instream.read((char*)cnt, ALP_SIZE*sizeof(cnt[0]));


#ifndef HUGE_PAGE
    eprintln("FMI - NO -- Huge Page allocation..");
    occb = (index_t*)aligned_alloc(64, (sizeof(occb[0]) * 4 * m + 63) / 64 * 64);
#else
    eprintln("FMI - Huge Page allocation..");
    occb = (index_t*) mmap(NULL, (sizeof(occb[0]) * 4 * m + 63) / 64 * 64, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB, -1, 0);
#endif

    eprintln("MEM-SIZE: FMI OCCB: %lld",(long long)((sizeof(occb[0]) * 4 * m + 63) / 64 * 64));
    instream.read((char*)occb, 4*m*sizeof(occb[0]));
    instream.close();
}

template<typename index_t>
void FMI<index_t>::save(string filename) const {
    ofstream outstream(filename.c_str(), ofstream::binary);
    outstream.seekp(0);
    outstream.write((char*)cnt, ALP_SIZE*sizeof(cnt[0]));
    outstream.write((char*)occb, 4*m*sizeof(occb[0]));
    outstream.close();
}

template<typename index_t>
FMI<index_t>::FMI(const string &t, index_t *sa, string _bases, string parent_filename): 
    bases(_bases), 
    ALP_SIZE((int)_bases.size()), 
    n((index_t)t.size()),
    m(((n+1+INDEX_T_BITS-1) / INDEX_T_BITS * 2LL)) {
    assert(bases.find('$') == string::npos);
    assert(bases=="@ACGT");
    assert(is_sorted(bases.begin(), bases.end()));
    assert(t[0] == '$' || t.back() == '$');

    string bin_filename = parent_filename + ".fmi4.compressed.binned.u";
    for(const auto &s:{sizeof(index_t)}) {
        bin_filename += string(".") + to_string(s);
    }
    if(ifstream(bin_filename.c_str()).good()) {
        eprintln("Found existing %s!!", (char*)bin_filename.c_str());
        load(bin_filename); 
        eprintln("Load successful.");
        return;
    } else {
        eprintln("No existing %s. Building...", (char*)bin_filename.c_str());
    }
    
    vector<index_t> true_sa;
    if(sa == NULL) {
        vector<int64_t> _sa(n);
        saisxx(t.c_str(), _sa.data(), (int64_t)n);
        true_sa.resize(n);
        for(index_t i=0;i<n;i++) true_sa[i]=_sa[i];
        sa = true_sa.data();
        eprintln("SA done (in fmi4).");
    }

    cnt = new index_t[ALP_SIZE]();
    for(char c:t) {
        for(int i=0; i<ALP_SIZE; i++) {
            cnt[i] += (c < bases[i]);
        }
    }

    auto bwt = [&](int64_t i){ 
        return t[sa[i]==0 ? n-1 : sa[i]-1]; 
    };


    occb = (index_t*)aligned_alloc(64, (sizeof(occb[0]) * 4 * m + 63) / 64 * 64);
    memset(occb, 0, sizeof(occb[0]) * 4 * m);
    for(int i=1; i<ALP_SIZE; i++) {
        const char c = bases[i];
        for(int64_t j=0; j<n; j++) {
            occb[((j+1) / INDEX_T_BITS * 2) * 4 + (i-1) * 2] |= (bwt(j)==c?1ULL:0ULL)<<((j+1)%INDEX_T_BITS);
        }
    }
    for(int64_t j=9; j<4*m; j+=2) {
        occb[j] = occb[j-8] + __builtin_popcountll((typename make_unsigned<index_t>::type)occb[j-9]);
    }
    eprintln("%s build done.", (char*)bin_filename.c_str());
    save(bin_filename);
    eprintln("save done.");
}


template<typename index_t>
FMI<index_t>::FMI(const string &t, index_t t_size, index_t *sa, string _bases, string parent_filename): 
    bases(_bases), 
    ALP_SIZE((int)_bases.size()), 
    n((index_t)t_size),
    m(((n+1+INDEX_T_BITS-1) / INDEX_T_BITS * 2LL)) {
    assert(bases.find('$') == string::npos);
    assert(bases=="@ACGT");
    assert(is_sorted(bases.begin(), bases.end()));
//    assert(t[0] == '$' || t.back() == '$');

    string bin_filename = parent_filename + ".fmi4.compressed.binned.u";
    for(const auto &s:{sizeof(index_t)}) {
        bin_filename += string(".") + to_string(s);
    }
    if(ifstream(bin_filename.c_str()).good()) {
        eprintln("Found existing %s!!", (char*)bin_filename.c_str());
        load(bin_filename); 
        eprintln("Load successful.");
        return;
    } else {
        eprintln("No existing %s. Building...", (char*)bin_filename.c_str());
    }
    
    vector<index_t> true_sa;
    if(sa == NULL) {
        vector<int64_t> _sa(n);
        saisxx(t.c_str(), _sa.data(), (int64_t)n);
        true_sa.resize(n);
        for(index_t i=0;i<n;i++) true_sa[i]=_sa[i];
        sa = true_sa.data();
        eprintln("SA done (in fmi4).");
    }

    cnt = new index_t[ALP_SIZE]();
    for(char c:t) {
        for(int i=0; i<ALP_SIZE; i++) {
            cnt[i] += (c < bases[i]);
        }
    }

    auto bwt = [&](int64_t i){ 
        return t[sa[i]==0 ? n-1 : sa[i]-1]; 
    };


    occb = (index_t*)aligned_alloc(64, (sizeof(occb[0]) * 4 * m + 63) / 64 * 64);
    memset(occb, 0, sizeof(occb[0]) * 4 * m);
    for(int i=1; i<ALP_SIZE; i++) {
        const char c = bases[i];
        for(int64_t j=0; j<n; j++) {
            occb[((j+1) / INDEX_T_BITS * 2) * 4 + (i-1) * 2] |= (bwt(j)==c?1ULL:0ULL)<<((j+1)%INDEX_T_BITS);
        }
    }
    for(int64_t j=9; j<4*m; j+=2) {
        occb[j] = occb[j-8] + __builtin_popcountll((typename make_unsigned<index_t>::type)occb[j-9]);
    }
    eprintln("%s build done.", (char*)bin_filename.c_str());
    save(bin_filename);
    eprintln("save done.");
}



template<typename index_t>
typename FMI<index_t>::Interval FMI<index_t>::backward_extend(Interval intv, char a) const {
#ifndef NO_DNA_ORD 
    const int i = __lg(a-'@'+1); // @ACGT -> 01234
#else
    const int i = a+1; // ACGT -> 1234
#endif 
    Interval ret = {
        cnt[i] + occ(i, intv.low/*-1*/),
        cnt[i] + occ(i, intv.high/*-1*/)
    };
    return ret;
}

template<typename index_t>
FMI<index_t>::~FMI(){
    eprintln("FMI memory deallocated\n");
#ifndef HUGE_PAGE
    free(occb);
#else
    munmap(occb, (sizeof(occb[0]) * 4 * m + 63) / 64 * 64);
#endif
}
#endif
