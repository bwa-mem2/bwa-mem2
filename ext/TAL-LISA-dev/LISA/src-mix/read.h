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
#include <cstdlib>
#include <zlib.h>
#include "kseq.h"
#include "utils.h"
//#include "bwa.h"
KSEQ_INIT(gzFile, gzread)


int64_t pac_seq_len(const char *fn_pac)
{
	FILE *fp;
	int64_t pac_len;
	uint8_t c;
	fp = xopen(fn_pac, "rb");
	err_fseek(fp, -1, SEEK_END);
	pac_len = err_ftell(fp);
	err_fread_noeof(&c, 1, 1, fp);
	err_fclose(fp);
	return (pac_len - 1) * 4 + (int)c;
}

void read_seq_lisa(string filename, string &reference_seq){
	uint8_t *buf2;
	int64_t i, pac_size, seq_len;
	FILE *fp;
	filename += ".pac";	
	const char *fn_pac = filename.c_str();

	// initialization
	seq_len = pac_seq_len(fn_pac);
    assert(seq_len > 0);
    assert(seq_len <= 0x7fffffffffL);
	fp = xopen(fn_pac, "rb");

	// prepare sequence
	pac_size = (seq_len>>2) + ((seq_len&3) == 0? 0 : 1);
	buf2 = (uint8_t*)calloc(pac_size, 1);
    assert(buf2 != NULL);
	err_fread_noeof(buf2, 1, pac_size, fp);
	err_fclose(fp);
	for (i = 0; i < seq_len; ++i) {
		int nt = buf2[i>>2] >> ((3 - (i&3)) << 1) & 3;
        switch(nt)
        {
            case 0:
                reference_seq += "A";
            break;
            case 1:
                reference_seq += "C";
            break;
            case 2:
                reference_seq += "G";
            break;
            case 3:
                reference_seq += "T";
            break;
            default:
                fprintf(stderr, "ERROR! Value of nt is not in 0,1,2,3!");
                exit(EXIT_FAILURE);
        }
	}
	free(buf2);
}

string read_seq(string filename) {
    static char in[123465];
    FILE* fp = fopen(filename.c_str(),"r");
    assert(fp != NULL); 
    fscanf(fp,"%123464[^\n]",in);
    assert(in[0] == '>');

    string result = "";
    while(fscanf(fp,"%123464s",in)==1) {
        result += in;
    }
    static_assert(sizeof(result.size()) * __CHAR_BIT__ == 64, "you are not using 64-bit system?!");
    eprintln("Read ref seq length = %lu", result.size());
    
    for(char &c:result) if(!strchr(dna.c_str(),c)) {
        c = dna[rand()%4];
        // error_quit("Error: ref file not in upper-case ACGT format");
    }
    fclose(fp);
    return result;
}
pair<string, int> read_query(string filename) {
    //rem: static char in[123465];
    gzFile fp = gzopen(filename.c_str(), "r");
    kseq_t *ks = kseq_init(fp);

    string ret;
    int num_queries = 0;
    int max_query_len = 0;
	int64_t m, n;
	//rem: char *seq;
	m = n = 0;
	while (kseq_read(ks) >= 0)
	{
        num_queries++;
        string cur = string(ks->seq.s);

        max_query_len = max(max_query_len, (int)cur.size());
        //for(char &c:cur) if(!strchr(dna.c_str(),c)) {
        //    c = '.';
        //    // error_quit("Error: query set file not in upper-case ACGT format");
        //}
        ret += cur;
    }
    gzclose(fp);
    eprintln("Detected num of queries = %d, max_query_len = %d", num_queries, max_query_len);

    //Fix: Memory leak
    kseq_destroy(ks);
    return {ret, max_query_len};
}



pair<string, int> read_query_separated_with_dot(string filename) {
    //rem: static char in[123465];
    gzFile fp = gzopen(filename.c_str(), "r");
    kseq_t *ks = kseq_init(fp);

    string ret;
    int num_queries = 0;
    int max_query_len = 0;
    int64_t m, n;
    //rem: char *seq;
	m = n = 0;
	while (kseq_read(ks) >= 0)
	{
        num_queries++;
        string cur = string(ks->seq.s);
        max_query_len = max(max_query_len, (int)cur.size());
        for(char &c:cur) if(!strchr(dna.c_str(),c)) {
            c = '.';
            // error_quit("Error: query set file not in upper-case ACGT format");
        }
        ret += cur + ";";
    }
    gzclose(fp);
    eprintln("Detected num of queries = %d, max_query_len = %d", num_queries, max_query_len);
    //Fix: Memory leak
    kseq_destroy(ks);
    return {ret, max_query_len};
}







