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
*/

/* Contact: Heng Li <hli@jimmy.harvard.edu> */

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

int bwa_index(int argc, char *argv[]) // the "index" command
{
	int c;
	char *prefix = 0;
	while ((c = getopt(argc, argv, "p:")) >= 0) {
		if (c == 'p') prefix = optarg;
		else return 1;
	}

	if (optind + 1 > argc) {
		fprintf(stderr, "Usage: bwa-mem2 index [-p prefix] <in.fasta>\n");
		return 1;
	}
	if (prefix == 0) prefix = argv[optind];
	bwa_idx_build(argv[optind], prefix);
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
        fmi->build_index();
        delete fmi;
	}
	return 0;
}
