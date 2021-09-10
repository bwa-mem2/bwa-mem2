#include "common.h"

SMEM_out::SMEM_out(int _id, int _q_l, int _q_r, index_t _ref_l, index_t _ref_r)
{
		    id = _id;
		    q_l = _q_l;
		    q_r = _q_r;
		    ref_l = _ref_l;
		    ref_r = _ref_r;
}

vector_based_output::vector_based_output(int a)
{ 
	id  = a;

}
Output::Output(int a)
{ 
	id  = a;
}

void Info::set(int a, int b, index_t c, index_t d)
{
	l = a; r = b; intv = make_pair(c,d);
}


uint64_t Info::get_enc_str(){

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


void load(string filename, vector<char*> ptrs, vector<size_t> sizes) {
    ifstream instream(filename.c_str(), ifstream::binary);
    instream.seekg(0);
    assert(ptrs.size() == sizes.size());
    for(size_t i=0; i<ptrs.size(); i++) {
        instream.read(ptrs[i], sizes[i]);
    }
    instream.close();
}

void save(string filename, vector<char*> ptrs, vector<size_t> sizes) {
    ofstream outstream(filename.c_str(), ofstream::binary);
    outstream.seekp(0);
    assert(ptrs.size() == sizes.size());
    for(size_t i=0; i<ptrs.size(); i++) {
        outstream.write(ptrs[i], sizes[i]);
    }
    outstream.close();
}

int64_t FCLAMP(double inp, double bound) {
  if (inp < 0.0) return 0;
  return (inp > bound ? bound : (int64_t)inp);
}

