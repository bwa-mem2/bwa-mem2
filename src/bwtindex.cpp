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
#include "bwtbuild.h"

#ifdef _DIVBWT
#include "divsufsort.h"
#endif

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

int bwa_index(int argc, char *argv[]) // the "index" command
{
	int c, algo_type = BWTALGO_AUTO, is_64 = 0, block_size = 10000000;
	char *prefix = 0, *str;
	while ((c = getopt(argc, argv, "6a:p:b:")) >= 0) {
		switch (c) {
		case 'a': // if -a is not set, algo_type will be determined later
			if (strcmp(optarg, "rb2") == 0) algo_type = BWTALGO_RB2;
			else if (strcmp(optarg, "bwtsw") == 0) algo_type = BWTALGO_BWTSW;
			else if (strcmp(optarg, "is") == 0) algo_type = BWTALGO_IS;
			else err_fatal(__func__, "unknown algorithm: '%s'.", optarg);
			break;
		case 'p': prefix = strdup(optarg); break;
		case '6': is_64 = 1; break;
		case 'b':
			block_size = strtol(optarg, &str, 10);
			if (*str == 'G' || *str == 'g') block_size *= 1024 * 1024 * 1024;
			else if (*str == 'M' || *str == 'm') block_size *= 1024 * 1024;
			else if (*str == 'K' || *str == 'k') block_size *= 1024;
			break;
		default: return 1;
		}
	}

	if (optind + 1 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa index [options] <in.fasta>\n\n");
		fprintf(stderr, "Options: -a STR    BWT construction algorithm: bwtsw, is or rb2 [auto]\n");
		fprintf(stderr, "         -p STR    prefix of the index [same as fasta name]\n");
		fprintf(stderr, "         -b INT    block size for the bwtsw algorithm (effective with -a bwtsw) [%d]\n", block_size);
		fprintf(stderr, "         -6        index files named as <in.fasta>.64.* instead of <in.fasta>.* \n");
		fprintf(stderr, "\n");
		fprintf(stderr,	"Warning: `-a bwtsw' does not work for short genomes, while `-a is' and\n");
		fprintf(stderr, "         `-a div' do not work not for long genomes.\n\n");
		return 1;
	}
	if (prefix == 0) {
		prefix = (char *) malloc(strlen(argv[optind]) + 4);
		strcpy(prefix, argv[optind]);
		if (is_64) strcat(prefix, ".64");
	}
	bwa_idx_build(argv[optind], prefix, algo_type, block_size);
	free(prefix);
	return 0;
}

int bwa_idx_build(const char *fa, const char *prefix, int algo_type, int block_size)
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
        build_index(prefix);
	}
	return 0;
}
