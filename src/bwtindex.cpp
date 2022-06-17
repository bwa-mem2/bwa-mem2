/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

   Modified Copyright (C) 2019 Intel Corporation, Heng Li.
   Contacts: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
   Heng Li <hli@jimmy.harvard.edu> 
*/


#include "sais.h"
#include "read.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <zlib.h>
#include "bntseq.h"
#include "bwa.h"
#include "bwt.h"
#include "utils.h"
#include "FMI_search.h"
#include "LISA_search.h"

int bwa_index(int argc, char *argv[]) // the "index" command
{
	string bin_file = argv[0];
	string mem2_home = get_abs_location(bin_file);

	argc = argc - 1;
	argv = argv + 1;
		
	int c;
	char *prefix = 0;
	int min_seed_len = 0;
	uint64_t num_rmi_leaf = 0;
	while ((c = getopt(argc, argv, "p:k:l:")) >= 0) {
		if (c == 'p') prefix = optarg;
		else if (c == 'k') min_seed_len = atoi(optarg);
		else if (c == 'l') num_rmi_leaf = atoi(optarg);
		else return 1;
	}

	if (optind + 1 > argc) {
		fprintf(stderr, "Usage: bwa-mem2 index [-p prefix] <in.fasta>\n");
		return 1;
	}
	
	assert(num_rmi_leaf >= 0 && num_rmi_leaf < UINT64_MAX);
	assert(min_seed_len >= 0 && min_seed_len < INT_MAX);
	if (min_seed_len == 0) min_seed_len = 19; // default value 
	if (prefix == 0) prefix = argv[optind];
#ifndef ENABLE_LISA 
	bwa_idx_build(argv[optind], prefix);
#else
	lisa_idx_build(argv[optind], prefix, min_seed_len, num_rmi_leaf, mem2_home);
#endif
	return 0;
}

int bwa_idx_build(const char *fa, const char *prefix)
{
	extern void bwa_pac_rev_core(const char *fn, const char *fn_rev);

	clock_t t;
	int64_t l_pac;
	{ // nucleotide indexing
		gzFile fp = xzopen(fa, "r");
		t = clock();
		fprintf(stderr, "[bwa_index] Pack FASTA... ");
		l_pac = bns_fasta2bntseq(fp, prefix, 1);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		err_gzclose(fp);
        FMI_search *fmi = new FMI_search(prefix);
        // fmi->build_index();
        fmi->build_index(0);
        delete fmi;
	}
	return 0;
}

int lisa_idx_build(const char *fa, const char *prefix, int min_seed_len, uint64_t num_rmi_leaf, string mem2_home)
{
	extern void bwa_pac_rev_core(const char *fn, const char *fn_rev);

	clock_t t;
	int64_t l_pac;
	 // nucleotide indexing
	gzFile fp = xzopen(fa, "r");
	t = clock();
	fprintf(stderr, "[lisa_index] Pack FASTA... ");
	l_pac = bns_fasta2bntseq(fp, prefix, 1);
	err_gzclose(fp);
    FMI_search *fmi = new FMI_search(prefix);
    fmi->build_index(0);
    delete fmi;

    string ref_seq_file = (string) prefix;
		
    string seq; 
	{ 
    	gzFile fp = xzopen(ref_seq_file.c_str(), "r");
    	read_seq_lisa(ref_seq_file, seq);
    	fprintf(stderr, "Read ref file done.\n");
 	}
    string path = mem2_home + "/ext/TAL";
	LISA_search<index_t> *lisa =  new LISA_search<index_t>(seq, seq.size(), ref_seq_file, min_seed_len + 1, num_rmi_leaf, path);

	delete lisa;
	fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	
	return 0;
}
