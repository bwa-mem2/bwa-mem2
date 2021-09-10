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
#ifndef QBWT_IPBWT_RMI_H
#define QBWT_IPBWT_RMI_H
#include "common.h"
#include "fmi.h"
#include "ipbwt_rmi.h"



template<typename index_t>
class QBWT_HYBRID {
    public:
	QBWT_HYBRID(){};
        QBWT_HYBRID(string t, index_t t_size, string ref_seq_filename, int K, int64_t num_rmi_leaf_nodes);
        ~QBWT_HYBRID();
        pair<int,int>* all_SMEMs(const char* p, const int p_len, pair<int,int>* ans_ptr, const int min_seed_length) const;
        pair<int,int>* print_all_SMEMs(const char* p, const int p_len, pair<int,int>* ans_ptr, const int min_seed_length, const int &shift) const;
        index_t n;
        typedef typename FMI<index_t>::Interval Interval;
        Interval init_intv;
        void forward_step(const char* p, Interval &intv, int &l, int &r) const;
        void backward_step(const char* p, Interval &intv, int &l, int &r) const;

    // private:
        FMI<index_t> *fmi;
        typedef float linreg_t;
        typedef uint64_t kenc_t; 
        IPBWT_RMI<index_t, kenc_t> *rmi;

        struct LcpInfo {
            uint16_t s_width: 12;
            uint8_t bw_ext_msk: 4;
        };
        static constexpr uint16_t WID_MAX = (1<<12)-1;
        LcpInfo *lcpi;
        vector<pair<index_t, index_t>> b_width;
        uint8_t *lcpp1;
        static constexpr uint8_t LCPP1_MAX = (1<<8)-1;
        vector<pair<index_t, index_t>> large_lcpp1;
        int64_t get_lcp(index_t i) const;
        pair<Interval, index_t> forward_shrink_phase(Interval intv, char a) const;
        pair<index_t, index_t> advance_chunk(kenc_t first, pair<index_t, index_t> intv) const;
        void load(string filename);
        void save(string filename) const;

	void test();
};
template<typename index_t>
int64_t QBWT_HYBRID<index_t>::get_lcp(index_t i) const {
    return (lcpp1[i] == LCPP1_MAX ?
	    (int64_t)lower_bound(large_lcpp1.begin(), large_lcpp1.end(), i,
		[&](pair<index_t, index_t> p, index_t q){return p.first < q;})->second:
	    (int64_t)lcpp1[i]) - 1;
}

template<typename index_t>
pair<index_t, index_t> QBWT_HYBRID<index_t>::advance_chunk(kenc_t first, pair<index_t, index_t> intv) const {
    return rmi->backward_extend_chunk(first, {intv.first, intv.second});
}
template<typename index_t>
void QBWT_HYBRID<index_t>::load(string filename) {
    ifstream instream(filename.c_str(), ifstream::binary);
    instream.seekg(0);

//#ifndef HUGE_PAGE  
    lcpp1 = new uint8_t[n+1]();
//#else
//    lcpp1 = (uint8_t*) mmap(NULL, sizeof(lcpp1[0]) * (n+1), PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB, -1, 0);
//#endif


    eprintln("MEM-SIZE: lcp1 %lld", (long long)sizeof(lcpp1[0])*(n+1));
    instream.read((char*)lcpp1, (n+1)*sizeof(lcpp1[0]));
  
#ifndef HUGE_PAGE  
    lcpi = new LcpInfo[n+1]();
#else
    lcpi = (LcpInfo*) mmap(NULL, sizeof(lcpi[0])*(n+1), PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB, -1, 0);
#endif
    eprintln("MEM-SIZE: lcpi %lld", (long long)sizeof(lcpi[0])*(n+1));
    instream.read((char*)lcpi, (n+1)*sizeof(lcpi[0]));

    int64_t siz;
    instream.read((char*)&siz, sizeof(siz));
    b_width.resize(siz);
    eprintln("MEM-SIZE: b_width %lld", (long long)sizeof(b_width[0])*siz);

#define LOADV(v) do{\
    for(auto itl=v.begin(), itr=itl; itl<v.end(); itl=itr) {\
        static constexpr int Z = 1<<16;\
        itr = min(itl+Z, v.end());\
        static typename std::remove_const<typename std::remove_reference<decltype(*itl)>::type>::type tmp[Z];\
        int len = itr-itl;\
        instream.read((char*)tmp, len*sizeof(tmp[0]));\
        for(int i=0; i<len; i++) itl[i] = tmp[i];\
    }\
}while(0)

    LOADV(b_width);

    instream.read((char*)&siz, sizeof(siz));
//    large_lcpp1.resize(siz);
//    LOADV(large_lcpp1);

    instream.close();
#undef LOADV
}

template<typename index_t>
void QBWT_HYBRID<index_t>::save(string filename) const {
    ofstream outstream(filename.c_str(), ofstream::binary);
    outstream.seekp(0);

    outstream.write((char*)lcpp1, (n+1)*sizeof(lcpp1[0]));
    
    outstream.write((char*)lcpi, (n+1)*sizeof(lcpi[0]));

    int64_t siz = (int64_t)b_width.size();
    outstream.write((char*)&siz, sizeof(siz));
#define SAVEV(v) do{\
    for(auto itl=v.begin(), itr=itl; itl<v.end(); itl=itr) {\
        static constexpr int Z = 1<<16;\
        itr = min(itl+Z, v.end());\
        static typename std::remove_const<typename std::remove_reference<decltype(*itl)>::type>::type tmp[Z];\
        int len = itr-itl;\
        for(int i=0; i<len; i++) tmp[i] = itl[i];\
        outstream.write((char*)tmp, len*sizeof(tmp[0]));\
    }\
}while(0)

    SAVEV(b_width);

    siz = (int64_t)large_lcpp1.size();
    outstream.write((char*)&siz, sizeof(siz));
    SAVEV(large_lcpp1);

    outstream.close();
#undef SAVEV
}

template<typename index_t>
pair<typename FMI<index_t>::Interval, index_t> QBWT_HYBRID<index_t>::forward_shrink_phase(Interval intv, char a) const {

    index_t e[2] = {intv.low, intv.high};
    LcpInfo info[2] = {lcpi[e[0]], lcpi[e[1]]};
#define GET_WIDTH(i) (info[i].s_width == WID_MAX ? \
        (index_t)lower_bound(b_width.begin(), b_width.end(), e[i],  \
            [&](pair<index_t, index_t> p, index_t q){return p.first < q;})->second: \
        (index_t)info[i].s_width)
    index_t siz[2] = {GET_WIDTH(0), GET_WIDTH(1)};

#ifndef NO_DNA_ORD 
    const uint8_t msk = 1<<dna_ord(a);
#else 
    const uint8_t msk = 1<<(a);
#endif 
#define CAN_BW_EXTEND(i) (info[i].bw_ext_msk & msk)
    
    bool small;
#define big (!small)
#define FORWARD_SHRINK do{\
        e[small] = (e[small] > e[big] ? e[big] + siz[small] : e[big] - siz[small]);\
        info[small] = lcpi[e[small]];\
        siz[small] = GET_WIDTH(small);\
    }while(0)

    for(small = siz[1] < siz[0]; !CAN_BW_EXTEND(small); small = siz[1] < siz[0]) {
        FORWARD_SHRINK;
    }

    if(e[small] > e[big]) {
        return {{e[big], e[big] + siz[small]}, get_lcp(e[small])};
    } else {
        return {{e[big]-siz[small], e[big]}, get_lcp(e[small])};
    }
#undef FORWARD_SHRINK
#undef big
#undef CAN_BW_EXTEND
#undef GET_WIDTH
}

template<typename index_t>
inline void QBWT_HYBRID<index_t>::forward_step(const char *p, Interval &intv, int &l, int &r) const
{
    index_t siz;
    tie(intv, siz) = forward_shrink_phase(intv, p[l-1]);
    r = l + (int)siz;
}

template<typename index_t>
inline void QBWT_HYBRID<index_t>::backward_step(const char *p, Interval &intv, int &l, int &r) const
{
    for(int i=l-1; i>=0; i--) {
        auto next = fmi->backward_extend(intv, p[i]);
        if(next.low >= next.high) break;
        else l=i, intv=next;
    }
}

template<typename index_t>
pair<int,int>* QBWT_HYBRID<index_t>::all_SMEMs(const char* p, const int p_len, pair<int,int>* ans_ptr, const int min_seed_length) const {
    Interval intv = init_intv;
    auto end = ans_ptr;
    int l = p_len, r = p_len; // [l,r)
    backward_step(p, intv, l, r);
    if(r-l>=min_seed_length) {
        *(end++) = make_pair(l, r);
    }

    while(l != 0) {
        forward_step(p, intv, l, r);
        backward_step(p, intv, l, r);
        if(r-l>=min_seed_length) {
            *(end++) = make_pair(l, r);
        }
    }

    reverse(ans_ptr, end);
    return end;
}

template<typename index_t>
pair<int,int>* QBWT_HYBRID<index_t>::print_all_SMEMs(const char* p, const int p_len, pair<int,int>* ans_ptr, const int min_seed_length, const int &shift) const {
    Interval intv = init_intv;
    vector<Interval> vs;
    auto end = ans_ptr;
    int l = p_len, r = p_len; // [l,r)
    backward_step(p, intv, l, r);
    if(r-l>=min_seed_length) {
        *(end++) = make_pair(l, r);
        vs.push_back(intv);
    }

    while(l != 0) {
        forward_step(p, intv, l, r);
        backward_step(p, intv, l, r);
        if(r-l>=min_seed_length) {
            *(end++) = make_pair(l, r);
            vs.push_back(intv);
        }
    }

    reverse(ans_ptr, end);
    reverse(vs.begin(), vs.end());
    for(int i=0; i<(int)vs.size(); i++) {
        printf("+ %ld %ld %lld %lld\n", 
                (long)ans_ptr[i].first+shift, (long)ans_ptr[i].second+shift,(long long) 1LL*vs[i].low, (long long)1LL*vs[i].high);
    }
    return end;
}


template<typename index_t>
void QBWT_HYBRID<index_t>::test(){
	std::cout<<"Test called";
}
template<typename index_t>
QBWT_HYBRID<index_t>::~QBWT_HYBRID(){
	eprintln("qbwt rmi deallocated");
	delete fmi;
	delete rmi;

	delete lcpp1;

#ifndef HUGE_PAGE
	delete lcpi;
#else
	munmap(lcpi, sizeof(lcpi[0])*(n+1));
#endif

}


template<typename index_t>
QBWT_HYBRID<index_t>::QBWT_HYBRID(string t, index_t t_size, string ref_seq_filename, int K, int64_t num_rmi_leaf_nodes):
#ifdef REV_COMP 
#ifndef BWA_MEM_BUG 
    n(2*(index_t)t_size+2),
#else
    n(2*(index_t)t_size+1),
#endif
#else 
    n(t_size+1),
#endif 
    init_intv({0, n}) {

#ifdef REV_COMP
    string bin_filename = ref_seq_filename + ".qbwt4.walg.rev_comp";
#else
    string bin_filename = ref_seq_filename + ".qbwt4.walg";
#endif
    string rmi_filename = bin_filename;
    for(const auto &s:{sizeof(index_t)}) {
        bin_filename += string(".") + to_string(s);
    }
    if(ifstream(bin_filename.c_str()).good()) {
        eprintln("Found existing %s!!", (char*)bin_filename.c_str());
        load(bin_filename);
        fmi = new FMI<index_t>(t.c_str(), n, NULL, "@"+dna, bin_filename);
        rmi = new IPBWT_RMI<index_t, uint64_t>(t, n, rmi_filename, K, num_rmi_leaf_nodes, NULL);
        eprintln("Load successful.");
        eprintln("large lcp size = %lu", large_lcpp1.size());
        eprintln("large lcp space usage = %.6fN", (double)(large_lcpp1.size()*sizeof(large_lcpp1[0])*1.0/n));
        eprintln("large width size = %lu", b_width.size());
        eprintln("large width space usage = %.6fN", (double)(b_width.size()*sizeof(b_width[0])*1.0/n));
        return;
    } else {
        eprintln("No existing %s. Building...", (char*)bin_filename.c_str());
    }


    assert(t.find('@') == string::npos && t.find('$') == string::npos);
    
#ifdef REV_COMP
#ifndef BWA_MEM_BUG
    t.push_back('@');
#else
    eprintln("No char placed between ref seq and reverse complement, to replicate BWA-MEM bug.");
#endif
    // appending reverse complement
    for(int64_t i=(index_t)t.size()-1-(t.back()=='@');i>=0;i--) {
#ifndef NO_DNA_ORD 
        t.push_back(dna[3-dna_ord(t[i])]);
#else 
        t.push_back(dna[3-(__lg(t[i]-'A'+2)-1)]);
#endif 
    }
#endif 
    t.push_back('$');
    assert(n == (index_t)t.size());





    vector<index_t> sa(n);
    if(numeric_limits<index_t>::max() > n && numeric_limits<index_t>::min() < 0) {
        saisxx(t.c_str(), sa.data(), (index_t)n);
    } else {
        vector<int64_t> _sa(n);
        saisxx(t.c_str(), _sa.data(), (int64_t)n);
        for(index_t i=0;i<n;i++) sa[i]=_sa[i];
    }
    

    fmi = new FMI<index_t>(t.c_str(), sa.data(), "@"+dna, bin_filename); // do not directly call load/save!
    rmi = new IPBWT_RMI<index_t, uint64_t>(t, rmi_filename, K, num_rmi_leaf_nodes, sa.data());

    // build rnk = sa^(-1)

    // build lcp
    {
        lcpp1 = new uint8_t[n+1]();
        large_lcpp1.clear();
        auto write_lcp = [&](index_t i, int64_t result) {
            if(result+1 >= LCPP1_MAX) {
                large_lcpp1.push_back({i, (index_t)(result+1)});
                lcpp1[i] = LCPP1_MAX;
            } else {
                lcpp1[i] = (uint8_t)(result+1);
            }
        };
        {
            vector<index_t> rnk(n);
            for(index_t i=0; i<n; i++) rnk[sa[i]]=i;
            for(int64_t i=0, me=0; i<n; i++) {
                if(rnk[i]==0) me=-1;
                else {
                    me = max(decltype(me)(0),me-1);
                    for(const auto mx = max((index_t)i, sa[rnk[i]-1]);
                            mx+me<n && t[i+me]==t[sa[rnk[i]-1]+me]; me++);
                }
                write_lcp(rnk[i], me);
            }
            write_lcp(n, -1);
        }
        sort(large_lcpp1.begin(), large_lcpp1.end());

        eprintln("large lcp size = %lu", large_lcpp1.size());
        eprintln("large lcp space usage = %.6fN", (double)(large_lcpp1.size()*sizeof(large_lcpp1[0])*1.0/n));
    }

    // build tree-structure of uni-lcps
    {
        // build order for traversing
        vector<index_t> dec_lcp_order(n-1); {
            int64_t max_lcp = LCPP1_MAX;
            for(auto p: large_lcpp1) max_lcp = max(max_lcp, (int64_t)p.second);
            vector<index_t> cnt(max_lcp+1, 0);
            for(index_t i=1; i<=n-1; i++) cnt[get_lcp(i)]++;
            for(index_t j=1; j<(index_t)cnt.size(); j++) {
                cnt[j] += cnt[j-1];
            }
            for(int64_t i=n-1; i>=1; i--) {
                dec_lcp_order[--cnt[get_lcp(i)]] = i;
            }
            reverse(dec_lcp_order.begin(), dec_lcp_order.end());
            assert(is_sorted(dec_lcp_order.begin(), dec_lcp_order.end(), [&](index_t i, index_t j){return make_pair(get_lcp(i), i) > make_pair(get_lcp(j),j);}));
            eprintln("dec lcp ordering ok.");
        }

        // il, ir for O(1) interval-merge algorithm
        vector<index_t> il(n);
        iota(il.begin(), il.end(), ((index_t)0));
        vector<index_t> ir=il;

        lcpi = new LcpInfo[n+1]();

        lcpi[0].s_width = lcpi[n].s_width = WID_MAX;
        b_width.push_back({0,n+1});
        b_width.push_back({n,n+1});


        for(auto itl=dec_lcp_order.begin(), itr=itl; itl!=dec_lcp_order.end(); itl=itr) {
            for(itr++; itr!=dec_lcp_order.end() &&
                       get_lcp(*itr)==get_lcp(*itl) &&
                       il[(*prev(itr))-1] == *itr; itr++);

            if(itr-itl>1) reverse(itl+1, itr);

            for(auto it=itl;it!=itr;it++) {
                const index_t i = *it;
                const index_t l = il[i-1];
                const index_t r = ir[i];
                il[r]=l;
                ir[l]=r;
                auto w = r-l+1;
                if(w>=WID_MAX) {
                    b_width.push_back({i, w});
                    lcpi[i].s_width = WID_MAX;
                } else {
                    lcpi[i].s_width = (uint16_t)w;
                }
                if(it+1==itr) {
                    uint8_t my_mask = 0;
                    for(int j=0;j<(int)dna.size();j++) {
                        auto intv = fmi->backward_extend({l, r+1}, dna[j]);
                        if(intv.low < intv.high) my_mask |= 1<<j;
                    }
                    lcpi[i].bw_ext_msk = my_mask;
                }
            }
        }
        sort(b_width.begin(), b_width.end());
        eprintln("large width size = %lu", b_width.size());
        eprintln("large width space usage = %.6fN", (double)(b_width.size()*sizeof(b_width[0])*1.0/n));
    }

    eprintln("%s build done.", (char*)bin_filename.c_str());
    save(bin_filename);
    eprintln("save done.");
}

#endif
