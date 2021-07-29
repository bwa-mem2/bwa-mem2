/*
This file is part of mpiBWA

NGS aligner inspired by BWA

The project was developped by Frederic Jarlier from Institut Curie and Nicolas Joly from Institut Pasteur

Copyright (C) 2016-2017  Institut Curie / Institut Pasteur

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <sys/mman.h>
#include <sys/stat.h>

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <inttypes.h> /* For PRIu64 */
#include <libgen.h>
#include <limits.h>   /* For PATH_MAX */
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#if NUMA_ENABLED
#include <numa.h>
#endif

#if AFF && (__linux__)
#include <sys/sysinfo.h>
int affy[256];
#endif

#include "bwa.h"
#include "bwamem.h"
#include "utils.h"
#include "macro.h"
#include "bandedSWA.h"
#include "profiling.h"
#include "FMI_search.h"
#include "fastmap.h"
#include "bwamem.h"
#include "parallel_aux.h"


#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "2.2.1"
#endif

//due to mpi_read_at limit buffer 1gb
#define DEFAULT_INBUF_SIZE (1024*1024*1024)
//#define DEFAULT_INBUF_SIZE 1024

/* We require 64bit offsets */
#ifndef MPI_OFFSET
#define MPI_OFFSET MPI_LONG_LONG
#endif

#ifdef TIMIMG
#define xfprintf fprintf
#else
#define xfprintf(...) /**/
#endif


#define STRIPING_FACTOR "6"
//#define STRIPING_UNIT "4194304"   // 4 MB
//#define STRIPING_UNIT "268435456"   // 256 MB
//#define STRIPING_UNIT "536870912"   // 500 MB
#define STRIPING_UNIT "1073741824"  // 1GB
//#define STRIPING_UNIT "1610612736"  // 1.5GB
//#define STRIPING_UNIT "2147483648"  // 2GB
//#define STRIPING_UNIT "2684354560"  // 2.5GB
//#define STRIPING_UNIT "3221225472"  // 3GB
//#define STRIPING_UNIT "3758096384"  // 3.5GB


#define NB_PROC  "16" //numer of threads for writing
#define CB_NODES "2" //numer of server for writing
#define CB_BLOCK_SIZE  "268435456" /* 256 MBytes - should match FS block size */
#define CB_BUFFER_SIZE  "3758096384" /* multiple of the block size by the number of proc*/
#define DATA_SIEVING_READ "enable"
#define min(a,b) (a>=b?b:a)
#define SMALL_STACK (1024*1024)
#define BIG_STACK (1024*1024*512)

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "2.2.1"
#endif


// ----------------------------------
uint64_t proc_freq, tprof[LIM_R][LIM_C], prof[LIM_R];
int nthreads;
int num_ranks = 1, myrank = 0;
int64_t reference_seq_len;
// ----------------------------------

void __cpuid(unsigned int i, unsigned int cpuid[4]) {
#ifdef _WIN32
    __cpuid((int *) cpuid, (int)i);

#else
    asm volatile
        ("cpuid" : "=a" (cpuid[0]), "=b" (cpuid[1]), "=c" (cpuid[2]), "=d" (cpuid[3])
            : "0" (i), "2" (0));
#endif
}



// --------------
// global vars
FMI_search *fmi;
extern uint64_t proc_freq, tprof[LIM_R][LIM_C];
extern unsigned char nst_nt4_table[256];
extern int num_ranks, myrank;
extern int nthreads;
uint8_t *ref_string;
int readLen, affy[256];
int64_t nreads, memSize;

int HTStatus()
{
    unsigned int cpuid[4];
    char platform_vendor[12];
    __cpuid(0, cpuid);
    ((unsigned int *)platform_vendor)[0] = cpuid[1]; // B
    ((unsigned int *)platform_vendor)[1] = cpuid[3]; // D
    ((unsigned int *)platform_vendor)[2] = cpuid[2]; // C
    std::string platform = std::string(platform_vendor, 12);

    __cpuid(1, cpuid);
    unsigned int platform_features = cpuid[3]; //D

     __cpuid(1, cpuid);
    unsigned int num_logical_cpus = (cpuid[1] >> 16) & 0xFF; // B[23:16]
     //fprintf(stderr, "#logical cpus: ", num_logical_cpus);
     unsigned int num_cores = -1;
     if (platform == "GenuineIntel") {
            __cpuid(4, cpuid);
            num_cores = ((cpuid[0] >> 26) & 0x3f) + 1; //A[31:26] + 1
            fprintf(stderr, "Platform vendor: Intel.\n");
        } else  {
            fprintf(stderr, "Platform vendor unknown.\n");
        }
        // fprintf(stderr, "#physical cpus: ", num_cores);
     int ht = platform_features & (1 << 28) && num_cores < num_logical_cpus;
     if (ht)
        fprintf(stderr, "CPUs support hyperThreading !!\n");
     return ht;
}


static void update_a(mem_opt_t *opt, const mem_opt_t *opt0)
{
	if (opt0->a) { // matching score is changed
		if (!opt0->b) opt->b *= opt->a;
		if (!opt0->T) opt->T *= opt->a;
		if (!opt0->o_del) opt->o_del *= opt->a;
		if (!opt0->e_del) opt->e_del *= opt->a;
		if (!opt0->o_ins) opt->o_ins *= opt->a;
		if (!opt0->e_ins) opt->e_ins *= opt->a;
		if (!opt0->zdrop) opt->zdrop *= opt->a;
		if (!opt0->pen_clip5) opt->pen_clip5 *= opt->a;
		if (!opt0->pen_clip3) opt->pen_clip3 *= opt->a;
		if (!opt0->pen_unpaired) opt->pen_unpaired *= opt->a;
	}
}


void memoryAlloc(ktp_aux_t *aux, worker_t &w, int32_t nreads, int32_t nthreads)
{
    mem_opt_t *opt = aux->opt;
    int32_t memSize = nreads;
    int32_t readLen = READ_LEN;

    
    w.regs = NULL; w.chain_ar = NULL; w.seedBuf = NULL;

    w.regs = (mem_alnreg_v *) calloc(memSize, sizeof(mem_alnreg_v));
    w.chain_ar = (mem_chain_v*) malloc (memSize * sizeof(mem_chain_v));
    w.seedBuf = (mem_seed_t *) calloc(sizeof(mem_seed_t),  memSize * AVG_SEEDS_PER_READ);

    assert(w.seedBuf  != NULL);
    assert(w.regs     != NULL);
    assert(w.chain_ar != NULL);

    w.seedBufSize = BATCH_SIZE * AVG_SEEDS_PER_READ;

    
    int64_t allocMem = memSize * sizeof(mem_alnreg_v) +
    memSize * sizeof(mem_chain_v) +
    sizeof(mem_seed_t) * memSize * AVG_SEEDS_PER_READ;
    fprintf(stderr, "------------------------------------------\n");
    fprintf(stderr, "1. Memory pre-allocation for Chaining: %0.4lf MB\n", allocMem/1e6);
    
     int64_t wsize = BATCH_SIZE * SEEDS_PER_READ;
    for(int l=0; l<nthreads; l++)
    {
        w.mmc.seqBufLeftRef[l*CACHE_LINE]  = (uint8_t *)
            _mm_malloc(wsize * MAX_SEQ_LEN_REF * sizeof(int8_t) + MAX_LINE_LEN, 64);
        w.mmc.seqBufLeftQer[l*CACHE_LINE]  = (uint8_t *)
            _mm_malloc(wsize * MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_LINE_LEN, 64);
        w.mmc.seqBufRightRef[l*CACHE_LINE] = (uint8_t *)
            _mm_malloc(wsize * MAX_SEQ_LEN_REF * sizeof(int8_t) + MAX_LINE_LEN, 64);
        w.mmc.seqBufRightQer[l*CACHE_LINE] = (uint8_t *)
            _mm_malloc(wsize * MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_LINE_LEN, 64);

        w.mmc.wsize_buf_ref[l*CACHE_LINE] = wsize * MAX_SEQ_LEN_REF;
        w.mmc.wsize_buf_qer[l*CACHE_LINE] = wsize * MAX_SEQ_LEN_QER;

        assert(w.mmc.seqBufLeftRef[l*CACHE_LINE]  != NULL);
        assert(w.mmc.seqBufLeftQer[l*CACHE_LINE]  != NULL);
        assert(w.mmc.seqBufRightRef[l*CACHE_LINE] != NULL);
        assert(w.mmc.seqBufRightQer[l*CACHE_LINE] != NULL);
    }

    for(int l=0; l<nthreads; l++) {
        w.mmc.seqPairArrayAux[l]      = (SeqPair *) malloc((wsize + MAX_LINE_LEN)* sizeof(SeqPair));
        w.mmc.seqPairArrayLeft128[l]  = (SeqPair *) malloc((wsize + MAX_LINE_LEN)* sizeof(SeqPair));
        w.mmc.seqPairArrayRight128[l] = (SeqPair *) malloc((wsize + MAX_LINE_LEN)* sizeof(SeqPair));
        w.mmc.wsize[l] = wsize;

        assert(w.mmc.seqPairArrayAux[l] != NULL);
        assert(w.mmc.seqPairArrayLeft128[l] != NULL);
        assert(w.mmc.seqPairArrayRight128[l] != NULL);
    }


    allocMem = (wsize * MAX_SEQ_LEN_REF * sizeof(int8_t) + MAX_LINE_LEN) * opt->n_threads * 2+
        (wsize * MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_LINE_LEN) * opt->n_threads  * 2 +
        wsize * sizeof(SeqPair) * opt->n_threads * 3;
    fprintf(stderr, "2. Memory pre-allocation for BSW: %0.4lf MB\n", allocMem/1e6);
    for (int l=0; l<nthreads; l++)
    {
        w.mmc.wsize_mem[l]     = BATCH_MUL * BATCH_SIZE *               readLen;
        w.mmc.matchArray[l]    = (SMEM *) _mm_malloc(w.mmc.wsize_mem[l] * sizeof(SMEM), 64);
        w.mmc.min_intv_ar[l]   = (int32_t *) malloc(w.mmc.wsize_mem[l] * sizeof(int32_t));
        w.mmc.query_pos_ar[l]  = (int16_t *) malloc(w.mmc.wsize_mem[l] * sizeof(int16_t));
        w.mmc.enc_qdb[l]       = (uint8_t *) malloc(w.mmc.wsize_mem[l] * sizeof(uint8_t));
        w.mmc.rid[l]           = (int32_t *) malloc(w.mmc.wsize_mem[l] * sizeof(int32_t));
        w.mmc.lim[l]           = (int32_t *) _mm_malloc((BATCH_SIZE + 32) * sizeof(int32_t), 64); // candidate not for reallocation, deferred for next round of changes.
    }

    allocMem = nthreads * BATCH_MUL * BATCH_SIZE * readLen * sizeof(SMEM) +
    nthreads * BATCH_MUL * BATCH_SIZE * readLen *sizeof(int32_t) +
    nthreads * BATCH_MUL * BATCH_SIZE * readLen *sizeof(int16_t) +
    nthreads * BATCH_MUL * BATCH_SIZE * readLen *sizeof(int32_t) +
    nthreads * (BATCH_SIZE + 32) * sizeof(int32_t);
    fprintf(stderr, "3. Memory pre-allocation for BWT: %0.4lf MB\n", allocMem/1e6);
    fprintf(stderr, "------------------------------------------\n");
}

int main(int argc, char *argv[]) {
	

	const char *mode = NULL;
	char *progname = basename(argv[0]);
	char *file_r1 = NULL, *file_r2 = NULL;
	char *buffer_r1, *buffer_r2;
	char *file_out = NULL;
	char *buffer_out;
	char *file_ref = NULL;
	char *rg_line = NULL, *hdr_line = NULL, *pg_line = NULL;
	char file_map[PATH_MAX], file_tmp[PATH_MAX];
	char *p, *q, *s, *e;
	uint8_t *a, *addr, *addr_map;
	int fd_in1, fd_in2;
	int proc_num, rank_num, rank_shr;
	int res, count;
	int files, nargs;
	int no_mt_io, i, c, copy_comment = 0;
	int ignore_alt = 0;
	double bef, aft;
	size_t localsize;
	size_t n = 0;
	off_t locoff, locsiz, *alloff, *curoff, maxsiz, totsiz, filsiz;
	struct stat stat_out, stat_r1, stat_r2;
    int dofixmate = 0;
    int  write_format = 2;
    int fixed_chunk_size = -1;
    int compression_level = 3;
	
    MPI_Aint size_shr;
	MPI_Comm comm_shr;
	MPI_File fh_r1, fh_r2, fh_out, fh_tmp;
	MPI_Offset *coff, m, size_map, size_tot;
	MPI_Status status;
	MPI_Win win_shr;

	mem_opt_t *opt, opt0;
	mem_pestat_t pes[4], *pes0 = NULL;
	bwaidx_t indix;
	bseq1_t *seqs;

	if (argc < 2) {
		fprintf(stderr, "program: %s is a MPI version of BWA MEM2\n"
			"version: %s\n"
			"\nusage : mpirun -n TOTAL_PROC %s mem -t 8 -o SAM_FILE REFERENCE_GENOME FASTQ_R1 [FASTQ_R2]\n"
            "\n\tTOTAL_PROC tells how many cores will be used by MPI to parallelize the computation.\n"
			"\nrequirements : from the reference genome index file generated with the command 'bwa index'\n"
            "\tyou need to create a reference genome map file with 'mpiBAWIdx' that comes along\n"
            "\twith this program as follows:\n"
            "\n\t\tmpiBWAIdx myReferenceGenome.fa\n\n"
			"\tIt creates a .map file that will be used in shared memory as reference genome.\n"
            "\ninput:\n"
            "\tREFERENCE_GENOME: reference genome name (e.g. myReferenceGenome.fa).\n"
            "\t\tDo not provide the '.map' extension of the file genareted with 'mpiBWAIdx'\n"
            "\n\tFASTQ_R1: fastq file for R1\n"
            "\n\tFASTQ_R2: fastq file for R2 if the data come from paired-end sequencing (optional)\n"
            "\noutput: SAM_FILE\n"
            "\noptions: 'bwa mem' options can be passed from command line (e.g. mpiBWA mem -t 8 -k 18)\n"
            "\nFor more detailed documentation visit:\n"
            "\thttps://github.com/bioinfo-pf-curie/mpiBWA\n"
            "\nCopyright (C) 2020  Institut Curie <http://www.curie.fr> \n"
            "\nThis program comes with ABSOLUTELY NO WARRANTY. \n"
            "This is free software, and you are welcome to redistribute it \n"
            "under the terms of the CeCILL License. \n"
	    	"\ncontact: Frederic Jarlier (frederic.jarlier@curie.fr) \n",
			progname, PACKAGE_VERSION, progname);
			return 1; }


	ktp_aux_t		 aux;
    worker_t         w;
	bool			 is_o	   = 0;
	int64_t			 nread_lim = 0;

	// uint64_t tim = __rdtsc();
	memset(&aux, 0, sizeof(ktp_aux_t));
	memset(pes, 0, 4 * sizeof(mem_pestat_t));
	for (i = 0; i < 4; ++i) pes[i].failed = 1;
	
	// opterr = 0;
	aux.fp = stdout;
	aux.opt = opt = mem_opt_init();
	memset(&opt0, 0, sizeof(mem_opt_t));


	if (strcmp(argv[1], "mem") == 0)
	{
		uint64_t tim = __rdtsc();
		kstring_t pg = {0,0,0};
		extern char *bwa_pg;

		fprintf(stderr, "-----------------------------\n");
#if __AVX512BW__
		fprintf(stderr, "Executing in AVX512 mode!!\n");
#endif
#if ((!__AVX512BW__) & (__AVX2__))
		fprintf(stderr, "Executing in AVX2 mode!!\n");
#endif
#if ((!__AVX512BW__) && (!__AVX2__) && (__SSE2__))
		fprintf(stderr, "Executing in SSE4.1 mode!!\n");
#endif
#if ((!__AVX512BW__) && (!__AVX2__) && (!__SSE2__))
		fprintf(stderr, "Executing in Scalar mode!!\n");
#endif
		fprintf(stderr, "-----------------------------\n");

		ksprintf(&pg, "@PG\tID:bwa\tPN:bwa\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);
		for (int i = 1; i < argc; ++i) ksprintf(&pg, " %s", argv[i]);
		ksprintf(&pg, "\n");
		bwa_pg = pg.s;
        
        //create pg_line for create_sam_header
        asprintf(&pg_line, "@PG\tID:bwa\tPN:bwa\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);
        for (int i = 1; i < argc; ++i) asprintf(&pg_line, "%s %s", pg_line, argv[i]);
        
	}	

    /*
	 * DEBUG
	 *
	 */
	fprintf(stderr, "DEBUG 0 \n");
	/* Validate provided command (first argument) */
	if (strcmp(argv[1], "mem") != 0) {
		fprintf(stderr, "%s: unsupported %s command\n", progname, argv[1]);
		return 1; }

	/*
	 * DEBUG
	 *
	 */
	fprintf(stderr, "DEBUG 1 \n");
	
	/* Parse input arguments */
    
	while ((c = getopt(argc, argv, "bg51qpaMCSPVYjk:K:c:v:s:r:t:R:A:B:O:E:U:w:L:d:T:Q:D:m:I:N:W:x:G:h:y:K:X:H:o:f:")) >= 0)
	{

        //special mpi param
        if (c == 'o') file_out = optarg;
        else if (c == 'b') write_format = 1;
        else if (c == 'g') write_format = 0;
        else if (c == 'f') dofixmate =1;
        //return to mem options
        

		else if (c == 'k') opt->min_seed_len = atoi(optarg), opt0.min_seed_len = 1;
		else if (c == '1') no_mt_io = 1;
		else if (c == 'x') mode = optarg;
		else if (c == 'w') opt->w = atoi(optarg), opt0.w = 1;
		else if (c == 'A') opt->a = atoi(optarg), opt0.a = 1, assert(opt->a >= INT_MIN && opt->a <= INT_MAX);
		else if (c == 'B') opt->b = atoi(optarg), opt0.b = 1, assert(opt->b >= INT_MIN && opt->b <= INT_MAX);
		else if (c == 'T') opt->T = atoi(optarg), opt0.T = 1, assert(opt->T >= INT_MIN && opt->T <= INT_MAX);
		else if (c == 'U')
			opt->pen_unpaired = atoi(optarg), opt0.pen_unpaired = 1, assert(opt->pen_unpaired >= INT_MIN && opt->pen_unpaired <= INT_MAX);
		else if (c == 't')
			opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1, assert(opt->n_threads >= INT_MIN && opt->n_threads <= INT_MAX);
		else if (c == 'o')
		{
			file_out = optarg;		
		}
		else if (c == 'q') {
			nread_lim = atoi(optarg);
		}
		else if (c == 'P') opt->flag |= MEM_F_NOPAIRING;
		else if (c == 'a') opt->flag |= MEM_F_ALL;
		else if (c == 'p') opt->flag |= MEM_F_PE | MEM_F_SMARTPE;
		else if (c == 'M') opt->flag |= MEM_F_NO_MULTI;
		else if (c == 'S') opt->flag |= MEM_F_NO_RESCUE;
		else if (c == 'Y') opt->flag |= MEM_F_SOFTCLIP;
		else if (c == 'V') opt->flag |= MEM_F_REF_HDR;
        else if (c == '5') opt->flag |= MEM_F_PRIMARY5 | MEM_F_KEEP_SUPP_MAPQ; // always apply MEM_F_KEEP_SUPP_MAPQ with -5
        else if (c == 'q') opt->flag |= MEM_F_KEEP_SUPP_MAPQ;
		else if (c == 'c') opt->max_occ = atoi(optarg), opt0.max_occ = 1;
		else if (c == 'd') opt->zdrop = atoi(optarg), opt0.zdrop = 1;
		else if (c == 'v') bwa_verbose = atoi(optarg);
		else if (c == 'j') ignore_alt = 1;
		else if (c == 'r')
			opt->split_factor = atof(optarg), opt0.split_factor = 1.;
		else if (c == 'D') opt->drop_ratio = atof(optarg), opt0.drop_ratio = 1.;
		else if (c == 'm') opt->max_matesw = atoi(optarg), opt0.max_matesw = 1;
		else if (c == 's') opt->split_width = atoi(optarg), opt0.split_width = 1;
		else if (c == 'G')
			opt->max_chain_gap = atoi(optarg), opt0.max_chain_gap = 1;
		else if (c == 'N')
			opt->max_chain_extend = atoi(optarg), opt0.max_chain_extend = 1;
		else if (c == 'W')
			opt->min_chain_weight = atoi(optarg), opt0.min_chain_weight = 1;
		else if (c == 'y')
			opt->max_mem_intv = atol(optarg), opt0.max_mem_intv = 1;
		else if (c == 'C') aux.copy_comment = 1;
		else if (c == 'K') fixed_chunk_size = atoi(optarg);
		else if (c == 'X') opt->mask_level = atof(optarg);
		else if (c == 'h')
		{
			opt0.max_XA_hits = opt0.max_XA_hits_alt = 1;
			opt->max_XA_hits = opt->max_XA_hits_alt = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->max_XA_hits_alt = strtol(p+1, &p, 10);
		}
		else if (c == 'Q')
		{
			opt0.mapQ_coef_len = 1;
			opt->mapQ_coef_len = atoi(optarg);
			opt->mapQ_coef_fac = opt->mapQ_coef_len > 0? log(opt->mapQ_coef_len) : 0;
		}
		else if (c == 'O')
		{
			opt0.o_del = opt0.o_ins = 1;
			opt->o_del = opt->o_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->o_ins = strtol(p+1, &p, 10);
		}
		else if (c == 'E')
		{
			opt0.e_del = opt0.e_ins = 1;
			opt->e_del = opt->e_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->e_ins = strtol(p+1, &p, 10);
		}
		else if (c == 'L')
		{
			opt0.pen_clip5 = opt0.pen_clip3 = 1;
			opt->pen_clip5 = opt->pen_clip3 = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->pen_clip3 = strtol(p+1, &p, 10);
		}
		else if (c == 'R')
		{
		if ((rg_line = bwa_set_rg(optarg)) == 0) {
			free(opt);
			if (is_o)
				fclose(aux.fp);
				return 1;
            }
		}
		else if (c == 'H')
		{
			if (optarg[0] != '@')
			{
				FILE *fp;
				if ((fp = fopen(optarg, "r")) != 0)
				{
					char *buf;
					buf = (char *) calloc(1, 0x10000);
					while (fgets(buf, 0xffff, fp))
					{
						i = strlen(buf);
						assert(buf[i-1] == '\n');
						buf[i-1] = 0;
						hdr_line = bwa_insert_header(buf, hdr_line);
					}
					free(buf);
					fclose(fp);
				}
			} else hdr_line = bwa_insert_header(optarg, hdr_line);
		}
		else if (c == 'I')
		{
			aux.pes0 = pes;
			pes[1].failed = 0;
			pes[1].avg = strtod(optarg, &p);
			pes[1].std = pes[1].avg * .1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].std = strtod(p+1, &p);
			pes[1].high = (int)(pes[1].avg + 4. * pes[1].std + .499);
			pes[1].low  = (int)(pes[1].avg - 4. * pes[1].std + .499);
			if (pes[1].low < 1) pes[1].low = 1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].high = (int)(strtod(p+1, &p) + .499);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].low  = (int)(strtod(p+1, &p) + .499);
		}
		else {
			free(opt);
			if (is_o)
				fclose(aux.fp);
			return 1;
		}
	}
	

	/*
	 * DEBUG
	 *
	 */
	fprintf(stderr, "DEBUG 3 \n");

	/* Further input parsing */
	if (mode)
	{
		fprintf(stderr, "WARNING: bwa-mem2 doesn't work well with long reads or contigs; please use minimap2 instead.\n");
		if (strcmp(mode, "intractg") == 0)
		{
			if (!opt0.o_del) opt->o_del = 16;
			if (!opt0.o_ins) opt->o_ins = 16;
			if (!opt0.b) opt->b = 9;
			if (!opt0.pen_clip5) opt->pen_clip5 = 5;
			if (!opt0.pen_clip3) opt->pen_clip3 = 5;
		}
		else if (strcmp(mode, "pacbio") == 0 || strcmp(mode, "pbref") == 0 || strcmp(mode, "ont2d") == 0)
		{
			if (!opt0.o_del) opt->o_del = 1;
			if (!opt0.e_del) opt->e_del = 1;
			if (!opt0.o_ins) opt->o_ins = 1;
			if (!opt0.e_ins) opt->e_ins = 1;
			if (!opt0.b) opt->b = 1;
			if (opt0.split_factor == 0.) opt->split_factor = 10.;
			if (strcmp(mode, "ont2d") == 0)
			{
				if (!opt0.min_chain_weight) opt->min_chain_weight = 20;
				if (!opt0.min_seed_len) opt->min_seed_len = 14;
				if (!opt0.pen_clip5) opt->pen_clip5 = 0;
				if (!opt0.pen_clip3) opt->pen_clip3 = 0;
			}
			else
			{
				if (!opt0.min_chain_weight) opt->min_chain_weight = 40;
				if (!opt0.min_seed_len) opt->min_seed_len = 17;
				if (!opt0.pen_clip5) opt->pen_clip5 = 0;
				if (!opt0.pen_clip3) opt->pen_clip3 = 0;
			}
		}
		else
		{
			fprintf(stderr, "[E::%s] unknown read type '%s'\n", __func__, mode);
			free(opt);
			if (is_o)
				fclose(aux.fp);
			return 1;
        }
	} else update_a(opt, &opt0);
	

	bwa_fill_scmat(opt->a, opt->b, opt->mat);


	files = 0;
	file_ref = argv[optind+1+0];
	fprintf(stderr, "reference = %s \n", file_ref);

    file_ref = argv[optind+1];
    if (optind + 2 < argc) {
        file_r1 = argv[optind+2]; files += 1;
    }
    if (optind + 3 < argc) {
        file_r2 = argv[optind+3]; files += 1;
        opt->flag |= MEM_F_PE;
    }

    

	/* Derived file names */
	sprintf(file_map, "%s.map", file_ref);

	/* start up MPI */
	//res = MPI_Init(&argc, &argv);
	int threads_ok;
    int provided;
    res = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    assert(res == MPI_SUCCESS);
    threads_ok = provided >= MPI_THREAD_MULTIPLE;
    res = MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    assert(res == MPI_SUCCESS);
    res = MPI_Comm_rank(MPI_COMM_WORLD, &rank_num);
    assert(res == MPI_SUCCESS);
    if (rank_num == 0) {
                fprintf(stderr, "rank %d: is multithread ok ? %d\n", rank_num, threads_ok );
                fprintf(stderr, "all ranks will use %d threads per mpi job \n", opt->n_threads);
        }

	 // some internal structures
	char *p1, *q1, *e1, *p2, *q2, *e2;
	int line_number, line_number2;
    int NUM_THREADS = opt->n_threads;	
	int64_t bases;
	double local_time_spend_mapping = 0;
	double before_local_mapping = 0;
	double after_local_mapping	= 0;
	double total_time_local_mapping  = 0;
	double total_time_mapping = 0;
	double total_time_reading_seq = 0;
	double grand_total_time_reading_seq = 0;
	double total_time_parsing = 0;
	double grand_total_time_parsing = 0;
	double total_time_writing = 0;
	double grand_total_time_writing = 0;
	size_t reads_r1, reads_r2, reads;	
	size_t offset_chunk;
	size_t offset_chunk_2;
	size_t size_chunk;
	size_t size_chunk_2;
	size_t total_local_reads_aligned= 0;
	size_t total_reads_check 		= 0;
	size_t local_num_reads = 0;
	size_t total_num_reads = 0;
	size_t grand_total_num_reads = 0;
	size_t *begin_offset_chunk 	    = NULL;
	size_t *chunk_size 		        = NULL;
	size_t *reads_in_chunk 		    = NULL;

    size_t *all_begin_offset_chunk  = NULL;
    size_t *all_chunk_size          = NULL;
    size_t *all_reads_in_chunk      = NULL;

    size_t *all_begin_offset_chunk_2  = NULL;
    size_t *all_chunk_size_2          = NULL;
    size_t *all_reads_in_chunk_2      = NULL;

    uint8_t magic[] = "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00";
    pthread_attr_t attr;
    pthread_t threads[NUM_THREADS];


	//MPI_Info finfo;
	//MPI_Info_create(&finfo);
	/*
	 * In this part you shall adjust the striping factor and unit according
	 * to the underlying filesystem.
	 * Harmless for other file system.
	 *
	 */
	//MPI_Info_set(finfo,"striping_factor", STRIPING_FACTOR);
	//MPI_Info_set(finfo,"striping_unit", STRIPING_UNIT); //2G striping
	//MPI_Info_set(finfo,"ind_rd_buffer_size", STRIPING_UNIT); //2gb buffer
	//MPI_Info_set(finfo,"romio_ds_read",DATA_SIEVING_READ);

	/*
	 * for collective reading and writing
	 * should be adapted too and tested according to the file system
	 * Harmless for other file system.
	 */
	//MPI_Info_set(finfo,"nb_proc", NB_PROC);
	//MPI_Info_set(finfo,"cb_nodes", CB_NODES);
	//MPI_Info_set(finfo,"cb_block_size", CB_BLOCK_SIZE);
	//MPI_Info_set(finfo,"cb_buffer_size", CB_BUFFER_SIZE);
	//MPI_Info_set(finfo, "collective_buffering", "true");
	//MPI_Info_set(finfo,"romio_cb_write","enable");
	//MPI_Info_set(finfo,"romio_cb_read","enable");


	/* Check that output file (-o) is not null ... */
	if (file_out == NULL) {
		fprintf(stderr, "missing mandatory output file (-o)\n");
		res = MPI_Finalize();
		assert(res == MPI_SUCCESS);
		exit(2);
	}

	/* Check R1 & R2 file sizes */
	if (file_r1 != NULL && stat(file_r1, &stat_r1) == -1) {
		fprintf(stderr, "%s: %s\n", file_r1, strerror(errno));
		res = MPI_Finalize();
		assert(res == MPI_SUCCESS);
		exit(2);
	}

	if (file_r2 != NULL && stat(file_r2, &stat_r2) == -1) {
		fprintf(stderr, "%s: %s\n", file_r2, strerror(errno));
		res = MPI_Finalize();
		assert(res == MPI_SUCCESS);
		exit(2);
	}

    if (rank_num == 0)
               fprintf(stderr, "%s: controls are done. Start analyzing fastqs it could take few minutes...\n", __func__);


    char* file_out_ext = (char*)malloc((strlen(file_out) + 40) * sizeof(char));


    if(write_format == 2){
        sprintf(file_out_ext, "%s.sam", file_out);
    }
    else if (write_format == 1){
        sprintf(file_out_ext, "%s.bam", file_out);
    }
    else {
        sprintf(file_out_ext,"%s.gz",file_out);
    }

    //test if the file output exist and remove it
    if (rank_num == 0)
        if (stat(file_out_ext, &stat_out)) remove(file_out_ext);    

    fixed_chunk_size = fixed_chunk_size > 0? fixed_chunk_size : opt->chunk_size * opt->n_threads;

	if ((file_r1 != NULL && file_r2 != NULL  && (stat_r1.st_size == stat_r2.st_size)))  {
	
	/*
	 	 * 
	     * We are in case we are paired and not trimmed
	     *
	 	 */

		/* Work around build warning in non timing case */
		aft = 0; aft++;
		bef = 0; bef++;
		/*
		 * Rank 0 estimate the size of a read
		 *
		 */
		size_t slen, blen;
		off_t tmp_sz = 1024;
		if (rank_num == 0){
			// 512 Mo are
			int fd_tmp = open(file_r1, O_RDONLY, 0666);
			char *buffer = malloc(tmp_sz  + 1);
			buffer[tmp_sz] = '\0';
			size_t read_out = read(fd_tmp, buffer, tmp_sz);
			assert(read_out);
			assert(strlen(buffer) == tmp_sz);
			assert( *buffer == '@');
			/* Estimate sequence size (bases and bytes) */
			/* With this estimation we compute approximately the chuncks size and offset*/

			s = buffer;
			e = buffer + tmp_sz;
			p = q = s;
			while (q < e && *q != '\n') q++; p = ++q;
			while (q < e && *q != '\n') q++; blen = q - p; p = ++q;
			while (q < e && *q != '\n') q++; p = ++q;
			while (q < e && *q != '\n') q++; slen = q - s + 1;

			/* Split local buffer in chunks of 100000000 bases */
			free(buffer);
		}

		//Rank O broadcast the size of a read
	 	res = MPI_Bcast(&blen, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
		assert(res == MPI_SUCCESS);

		/*
		 * Split sequence files in chunks
		 */
		MPI_File mpi_fd_in1, mpi_fd_in2;
		bef = MPI_Wtime();

		res = MPI_File_open(MPI_COMM_WORLD, file_r1,  MPI_MODE_RDONLY , MPI_INFO_NULL, &mpi_fd_in1);
		assert(res == MPI_SUCCESS);
		
		res = MPI_File_open(MPI_COMM_WORLD, file_r1,  MPI_MODE_RDONLY , MPI_INFO_NULL, &mpi_fd_in2);
		assert(res == MPI_SUCCESS);
		
		/*
		 * first we parse the buffer and see
		 * how many reads we have
		 */

		size_t read_number=0;
		assert(fd_in1 != -1);
		size_t *goff = NULL; //global offset contain the start offset in the fastq
        goff = calloc( (proc_num * NUM_THREADS + 1) , sizeof(size_t));
        size_t *goff_inter = calloc( (proc_num * NUM_THREADS + 1) , sizeof(size_t));
        bef = MPI_Wtime();
        find_process_starting_offset_mt(goff, stat_r1.st_size, file_r1, proc_num, rank_num, NUM_THREADS);
        aft = MPI_Wtime();
        fprintf(stderr, "%s: rank %d time spend in finding process start offset = (%.02f) \n", __func__, rank_num, aft - bef);

        int i12=0;
        for ( i12 = 0; i12 < proc_num * NUM_THREADS + 1; i12++)  goff_inter[i12] = goff[i12];
                
        char *current_line = NULL;
        res = MPI_Allgather(&goff_inter[rank_num*NUM_THREADS], NUM_THREADS, MPI_LONG_LONG_INT, goff , NUM_THREADS, MPI_LONG_LONG_INT, MPI_COMM_WORLD);
        assert(res == MPI_SUCCESS);
        free(goff_inter);

        size_t local_num_reads          = 0;
        size_t total_num_reads          = 0;
        size_t u1						= 0;
        size_t *local_read_offsets      = calloc(1 , sizeof(size_t));
        size_t *local_read_bytes        = calloc(1, sizeof(size_t));
        int *local_read_size            = calloc(1, sizeof(int));


        assert( local_read_bytes != NULL);
        assert( local_read_offsets != NULL);
        assert( local_read_size != NULL);
        bef = MPI_Wtime();

        pthread_attr_t attr;
        pthread_attr_init(&attr);
        pthread_attr_setstacksize(&attr, SMALL_STACK);
        pthread_attr_setdetachstate(&attr, 0);

        pthread_t threads_1[NUM_THREADS];

        struct struct_data_thread_1 *td_1 =  malloc(NUM_THREADS * sizeof(struct struct_data_thread_1));

        size_t *local_num_reads_t            = calloc(NUM_THREADS, sizeof(size_t));
        size_t *total_num_reads_t            = calloc(NUM_THREADS, sizeof(size_t));
        size_t **local_read_offsets_t        = calloc(NUM_THREADS, sizeof(size_t*));
        size_t **local_read_bytes_t          = calloc(NUM_THREADS, sizeof(size_t*));
        int **local_read_size_t              = calloc(NUM_THREADS, sizeof(int*));

        int ret_code_1 = 0;
        int goff_idx = 0;
        for ( n = 0; n < NUM_THREADS; n++){

            goff_idx = (rank_num * NUM_THREADS) + n;
            td_1[n].offset_in_file_mt         = goff[goff_idx];
            td_1[n].size2read_mt              = goff[goff_idx + 1] - goff[goff_idx];
            td_1[n].file_r1_mt                = file_r1;
            td_1[n].local_num_reads_mt        = &local_num_reads_t[n];
            td_1[n].total_num_reads_mt        = &total_num_reads_t[n];
            td_1[n].local_read_offsets_mt     = &local_read_offsets_t[n];
            td_1[n].local_read_size_mt        = &local_read_size_t[n];
            td_1[n].local_read_bytes_mt       = &local_read_bytes_t[n];
            td_1[n].proc_num_mt               = proc_num;
            td_1[n].rank_num_mt               = rank_num;
            td_1[n].thread_num_mt             = n;
            td_1[n].previous_read_num         = 0;
            ret_code_1 = pthread_create(&threads_1[n], &attr, find_reads_size_and_offsets_mt, (void *)(&td_1[n]));

        }

            int g=0;
            total_num_reads = 0;
            for (n = 0; n < NUM_THREADS; n++){
                    pthread_join(threads_1[n], (void *)(&td_1[n]));
            total_num_reads += *(td_1[n].total_num_reads_mt);
        }

         local_read_offsets  = calloc(total_num_reads, sizeof(size_t));
            local_read_size     = calloc(total_num_reads, sizeof(int));
            local_read_bytes    = calloc(total_num_reads, sizeof(size_t));

            assert(local_read_offsets);
            assert(local_read_size);
            assert(local_read_bytes);

        size_t tmp_var = 0;
            for (n = 0; n < NUM_THREADS; n++){
                    td_1[n].local_read_offsets     = local_read_offsets;
                    td_1[n].local_read_size        = local_read_size;
                    td_1[n].local_read_bytes       = local_read_bytes;
                    td_1[n].previous_read_num      = tmp_var;
                    tmp_var                        += *(td_1[n].total_num_reads_mt);

                    ret_code_1 = pthread_create(&threads_1[n], &attr, copy_local_read_info_mt, (void *)(&td_1[n]));

            }

            for (n = 0; n < NUM_THREADS; n++){
                    pthread_join(threads_1[n], (void *)(&td_1[n]));

                    free(local_read_offsets_t[n]);
                    free(local_read_bytes_t[n]);
                    free(local_read_size_t[n]);
            }


            free(local_num_reads_t);
            free(total_num_reads_t);

            free(local_read_offsets_t);
            free(local_read_bytes_t);
            free(local_read_size_t);

            pthread_attr_destroy(&attr);
            free(td_1);

         MPI_Barrier(MPI_COMM_WORLD);
        aft = MPI_Wtime();
        fprintf(stderr, "%s: rank %d num reads parsed: %zu ::: time spend reading and parsing entire fastq = (%.02f) \n", __func__, rank_num, total_num_reads, aft - bef);

		if (goff) free(goff);
	
		MPI_File_close(&mpi_fd_in2);			
		res = MPI_Reduce(&total_num_reads, &grand_total_num_reads, 1,MPI_LONG_LONG_INT, MPI_SUM, 0,MPI_COMM_WORLD);
		assert ( res == MPI_SUCCESS);	
		res = MPI_Bcast(&grand_total_num_reads, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
		assert ( res == MPI_SUCCESS);	

		size_t num_reads_by_proc[proc_num];
		res = MPI_Allgather(&total_num_reads, 1, MPI_LONG_LONG_INT, num_reads_by_proc, 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);
		//we dispatch the previous result
	
		if ( rank_num == 0 || rank_num ==1 )
			fprintf(stderr, "rank %d ::: total_num_reads = %zu \n", rank_num, grand_total_num_reads);

		/*
		 * Now each rank compute buffers offset start and end.
		 */

		local_num_reads = total_num_reads;
		//now we estimate the number of chunk per rank
		size_t chunck_num = (local_num_reads * blen) / ( fixed_chunk_size / 2);
		chunck_num += 2; //the last chunk hold the remain bases

		size_t h=0;
		for ( h = 0; h < total_num_reads; h++){
			assert( local_read_size[h] == blen );
			assert( local_read_offsets[h] >= 0 );
		}


		// we allocate vector for chunks offset
		begin_offset_chunk 	= calloc(chunck_num, sizeof(size_t));
		chunk_size      	= calloc(chunck_num, sizeof(size_t));
		reads_in_chunk 		= calloc(chunck_num, sizeof(size_t));

		assert( begin_offset_chunk != NULL );
		assert( chunk_size != NULL );
		assert( reads_in_chunk != NULL );
	
		size_t chunk_count = 0;

		maxsiz = ( opt->chunk_size * opt->n_threads) / 2; 
		MPI_Barrier(MPI_COMM_WORLD);
		fprintf(stderr,"rank %d ::: Call find_chunks_info \n", rank_num);
		// the detail of he paramters is at the function definition
		// fprintf(stderr, "rank %d ::: begin_offset_chunk = %zu \n", rank_num, begin_offset_chunk[0]);
		// fprintf(stderr, "rank %d ::: begin_offset_chunk_2 = %zu \n", rank_num, begin_offset_chunk_2[0]);
		// fprintf(stderr, "rank %d ::: chunk_size = %zu \n", rank_num, chunk_size);
		// fprintf(stderr, "rank %d ::: chunk_size_2 = %zu \n", rank_num, chunk_size_2);

		bef = MPI_Wtime();
		find_chunks_info(begin_offset_chunk,
				 chunk_size,
				 reads_in_chunk,
				 local_read_size,
				 local_read_bytes,
	    		 local_read_offsets,
				 rank_num,
				 proc_num,
				 local_num_reads,
				 grand_total_num_reads,
				 maxsiz,
				 &chunk_count);		

        MPI_Barrier(MPI_COMM_WORLD);

		aft = MPI_Wtime();
		fprintf(stderr, "%s: rank %d time spend evaluating chunks = (%.02f) \n", __func__, rank_num, aft - bef);

		free(local_read_offsets);
		free(local_read_size);

        bef = MPI_Wtime();
        

        //we get all_chunk_size, all_begin_offset_chunk, all_reads_in_chunk
        size_t total_chunks = 0;

        res = MPI_Reduce(&chunk_count, &total_chunks, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        assert ( res == MPI_SUCCESS);
        res = MPI_Bcast(&total_chunks, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        assert ( res == MPI_SUCCESS);

        //fprintf(stderr, "%s: rank %d total chunks = %zu \n", __func__, rank_num, total_chunks);
        //fprintf(stderr, "%s: rank %d local chunks = %zu \n", __func__, rank_num, chunk_count);

        size_t chunks_per_rank[proc_num];
        res = MPI_Allgather(&chunk_count, 1, MPI_LONG_LONG_INT, chunks_per_rank, 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);

        //for (i = 0; i < proc_num; i++) fprintf(stderr, "%s: rank %d chck per rank[%d] = %zu \n", __func__, rank_num, i, chunks_per_rank[i]);

        
        all_begin_offset_chunk = malloc(total_chunks * sizeof(size_t));
        all_begin_offset_chunk[0]=0;
        all_chunk_size = malloc(total_chunks * sizeof(size_t));
        all_chunk_size[0] = 0;
        all_reads_in_chunk =  malloc(total_chunks * sizeof(size_t));
        all_reads_in_chunk[0] = 0;
        
        int displ_chunk[proc_num];
        displ_chunk[0] = 0; 
       
        for (i = 1; i < proc_num; i++) displ_chunk[i] = (displ_chunk[i-1] + chunks_per_rank[i-1]);
      
        int indx=displ_chunk[rank_num];
        for (i = 0; i <  chunks_per_rank[rank_num]; i++) { 
            all_chunk_size[indx+i]=chunk_size[i];
            all_begin_offset_chunk[indx+i]=begin_offset_chunk[i];
            all_reads_in_chunk[indx+i]=reads_in_chunk[i];
        }

        if (rank_num > 0){
            res=MPI_Send(chunk_size, chunks_per_rank[rank_num], MPI_LONG_LONG_INT, 0, rank_num, MPI_COMM_WORLD);
            assert(res == MPI_SUCCESS);
            }
        else{
            for (i = 1; i < proc_num; i++){
                res=MPI_Recv(&(all_chunk_size[displ_chunk[i]]), chunks_per_rank[i], MPI_LONG_LONG_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                assert(res == MPI_SUCCESS);
            }
        }

        if (rank_num > 0){
            res=MPI_Send(begin_offset_chunk, chunks_per_rank[rank_num], MPI_LONG_LONG_INT, 0, rank_num, MPI_COMM_WORLD);
            assert(res == MPI_SUCCESS);
            }
        else{
            for (i = 1; i < proc_num; i++){
                res=MPI_Recv(&(all_begin_offset_chunk[displ_chunk[i]]), chunks_per_rank[i], MPI_LONG_LONG_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                assert(res == MPI_SUCCESS);
            }
        }


        if (rank_num > 0){
            res=MPI_Send(reads_in_chunk, chunks_per_rank[rank_num], MPI_LONG_LONG_INT, 0, rank_num, MPI_COMM_WORLD);
            assert(res == MPI_SUCCESS);
            }
        else{
            for (i = 1; i < proc_num; i++){
                res=MPI_Recv(&(all_reads_in_chunk[displ_chunk[i]]), chunks_per_rank[i], MPI_LONG_LONG_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                assert(res == MPI_SUCCESS);
            }
        }

            
        res=MPI_Bcast(all_chunk_size, total_chunks, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        res=MPI_Bcast(all_reads_in_chunk, total_chunks, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        res=MPI_Bcast(all_begin_offset_chunk, total_chunks, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

        assert(res == MPI_SUCCESS);


        for ( i = 0; i < total_chunks; i++)
                assert(all_chunk_size[i] != 0);

        free(chunk_size);
        free(begin_offset_chunk);
        free(reads_in_chunk);
      

        aft = MPI_Wtime();
        fprintf(stderr, "%s: rank %d time spend gathering chunks info = (%.02f) \n", __func__, rank_num, aft - bef);        

        #if NUMA_ENABLED
            int  deno = 1;
            int tc = numa_num_task_cpus();
            int tn = numa_num_task_nodes();
            int tcc = numa_num_configured_cpus();
            fprintf(stderr, "num_cpus: %d, num_numas: %d, configured cpus: %d\n", tc, tn, tcc);
            int ht = HTStatus();
            if (ht) deno = 2;
        
            if (nthreads < tcc/tn/deno) {
                fprintf(stderr, "Enabling single numa domain...\n\n");
                struct bitmask *mask = numa_bitmask_alloc(numa_num_possible_nodes());
                numa_bitmask_clearall(mask);
                numa_bitmask_setbit(mask, 0);
                numa_bind(mask);
                numa_bitmask_free(mask);
            }
        #endif
       
        #if AFF && (__linux__)
        {
            unsigned int cpuid[4];
            asm volatile
                ("cpuid" : "=a" (cpuid[0]), "=b" (cpuid[1]), "=c" (cpuid[2]), "=d" (cpuid[3])
                : "0" (0xB), "2" (1));
            int num_logical_cpus = cpuid[1] & 0xFFFF;

            asm volatile
                ("cpuid" : "=a" (cpuid[0]), "=b" (cpuid[1]), "=c" (cpuid[2]), "=d" (cpuid[3])
                : "0" (0xB), "2" (0));
            int num_ht = cpuid[1] & 0xFFFF;
            int num_total_logical_cpus = get_nprocs_conf();
            int num_sockets = num_total_logical_cpus / num_logical_cpus;
            fprintf(stderr, "#sockets: %d, #cores/socket: %d, #logical_cpus: %d, #ht/core: %d\n",
                num_sockets, num_logical_cpus/num_ht, num_total_logical_cpus, num_ht);
        
            for (int i=0; i<num_total_logical_cpus; i++) affy[i] = i;
            int slookup[256] = {-1};

            if (num_ht == 2 && num_sockets == 2)
            {
                for (int i=0; i<num_total_logical_cpus; i++) {
                    std::ostringstream ss;
                    ss << i;
                    std::string str = "/sys/devices/system/cpu/cpu"+ ss.str();
                    str = str +"/topology/thread_siblings_list";

                    FILE *fp = fopen(str.c_str(), "r");
                    if (fp == NULL) {
                        fprintf("Error: Cant open the file..\n");
                        break;
                    }
                    else {
                        int a, b, v;
                        char ch[10] = {'\0'};
                        fgets(ch, 10, fp);
                        v = sscanf(ch, "%u,%u",&a,&b);
                        if (v == 1) v = sscanf(ch, "%u-%u",&a,&b);
                        if (v == 1) {
                            fprintf(stderr, "Mis-match between HT and threads_sibling_list...%s\n", ch);
                            fprintf(stderr, "Continuing with default affinity settings..\n");
                            break;
                        }
                        slookup[a] = 1;
                        slookup[b] = 2;
                        fclose(fp);
                    }
                }
                int a = 0, b = num_total_logical_cpus / num_ht;
                for (int i=0; i<num_total_logical_cpus; i++) {
                    if (slookup[i] == -1) {
                    fprintf(stderr, "Unseen cpu topology..\n");
                    break;
                    }
                    if (slookup[i] == 1) affy[a++] = i;
                    else affy[b++] = i;
                }
            }
        }
        #endif
     
        uint64_t tim = __rdtsc();
        fprintf(stderr, "* Ref file: %s\n", file_ref);
        aux.fmi = new FMI_search(file_ref);
        //aux.fmi->load_index(); 
        
              
        uint8_t *shared_ref;    
        MPI_Comm comm_shr;
        int count4, rank_shr = 0;    
        size_t size_read;
 
        res = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_shr);
        assert(res == MPI_SUCCESS);

        res = MPI_Comm_rank(comm_shr, &rank_shr);
        assert(res == MPI_SUCCESS);

        char binary_seq_file[200];
        sprintf(binary_seq_file, "%s.0123", file_ref);

        char pac_file[200];
        sprintf(pac_file, "%s.pac", file_ref);

        MPI_File fh_ref_file, fh_pac_file;
        MPI_Aint size_shr_ref, size_shr_pac;
        MPI_Offset size_map_ref,size_map_pac;
        MPI_Win win_shr_ref, win_shr_pac;

        uint8_t *addr_map_ref;
        uint8_t *addr_map_pac;       


        MPI_Info win_info;
        MPI_Info_create(&win_info);
        MPI_Info_set(win_info, "alloc_shared_non_contig", "true");

        res = MPI_File_open(MPI_COMM_WORLD, binary_seq_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_ref_file);
        assert(res == MPI_SUCCESS);

        res = MPI_File_open(MPI_COMM_WORLD, pac_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_pac_file);
        assert(res == MPI_SUCCESS);

        res = MPI_File_get_size(fh_pac_file, &size_map_pac);
        assert(res == MPI_SUCCESS);
        
        res = MPI_File_get_size(fh_ref_file, &size_map_ref);
        assert(res == MPI_SUCCESS);

        size_shr_ref = (rank_shr == 0) ? size_map_ref : 0;
        res = MPI_Win_allocate_shared(size_shr_ref * sizeof(uint8_t), 8 , win_info, comm_shr, &addr_map_ref, &win_shr_ref);

        size_shr_pac = (rank_shr == 0) ? size_map_pac : 0;
        res = MPI_Win_allocate_shared(size_shr_pac * sizeof(uint8_t), 8 , win_info, comm_shr, &addr_map_pac, &win_shr_pac);

        res = MPI_Win_shared_query(win_shr_pac, MPI_PROC_NULL, &size_shr_pac, &res, &addr_map_pac);
        assert(res == MPI_SUCCESS);

        res = MPI_Win_shared_query(win_shr_ref, MPI_PROC_NULL, &size_shr_ref, &res, &addr_map_ref);
        assert(res == MPI_SUCCESS);

        uint8_t *c = addr_map_ref;
        uint8_t *d = addr_map_pac;
        
        while(rank_shr == 0) {
            res = MPI_File_read(fh_pac_file, d, size_shr_pac, MPI_UINT8_T, &status);
            assert(res == MPI_SUCCESS);
            res = MPI_Get_count(&status, MPI_UINT8_T, &count4);
            assert(res == MPI_SUCCESS);
            if (count4 == 0) break;
            d += count4;
            size_read += count4;
            size_tot += count4;
        }
        
        size_read = 0;

        while(rank_shr == 0) {
            res = MPI_File_read(fh_ref_file, c, size_shr_ref, MPI_UINT8_T, &status);
            assert(res == MPI_SUCCESS);
            res = MPI_Get_count(&status, MPI_UINT8_T, &count4);
            assert(res == MPI_SUCCESS);
            if (count4 == 0) break;
            c += count4;
            size_read += count4;
            size_tot += count4;
        }

        res = MPI_Win_fence(0, win_shr_ref);
        assert(res == MPI_SUCCESS);

        res = MPI_Win_fence(0, win_shr_pac);
        assert(res == MPI_SUCCESS);

        res = MPI_File_close(&fh_pac_file);
        assert(res == MPI_SUCCESS);

        res = MPI_File_close(&fh_ref_file);
        assert(res == MPI_SUCCESS);

        aux.ref_string = addr_map_ref;
        aux.fmi->load_index_mpi(addr_map_pac);
                
        tprof[FMI][0] += __rdtsc() - tim;

		///Create SAM header
		 /*
          * Create SAM header
          * TODO: Add line for BWA version
          */

        if (rank_num == 0) {
            if (write_format == 2)
                create_sam_header(file_out_ext, (aux.fmi->idx), &count, hdr_line, rg_line, pg_line, rank_num);
            else
                create_bam_header(file_out_ext, (aux.fmi->idx), &count, hdr_line, rg_line, pg_line, rank_num, compression_level);
        }
    
        if (hdr_line) free(hdr_line);
        if (rg_line) free(rg_line);
        if (pg_line) free(pg_line);

        res = MPI_File_open(MPI_COMM_WORLD, file_out_ext, MPI_MODE_WRONLY|MPI_MODE_APPEND, MPI_INFO_NULL, &fh_out);
        assert(res == MPI_SUCCESS);

        if (file_r1 != NULL) {
            res = MPI_File_open(MPI_COMM_WORLD, file_r1, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_r1);
            assert(res == MPI_SUCCESS);
        }
        if (file_r2 != NULL) {
            res = MPI_File_open(MPI_COMM_WORLD, file_r2, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_r2);
            assert(res == MPI_SUCCESS);
        }
                 
        
        fprintf(stderr, "* Done reading reference genome !!\n\n");

        if (ignore_alt)
        for (i = 0; i < fmi->idx->bns->n_seqs; ++i)
            fmi->idx->bns->anns[i].is_alt = 0;

		/*
		 *	In this part we initiate the worker threads  
		 *
		 */
        aux.ks2                       = 0;
        aux.ks                        = 0;
        aux.task_size = opt->chunk_size * opt->n_threads;
		fprintf(stderr, "\naux.task_size : %d\n", aux.task_size);
		worker_t	 w2;
		mem_opt_t	*opt			  = aux.opt;
		nthreads = opt->n_threads; // global variable for profiling!
        w2.nthreads = opt->n_threads;
        readLen = blen;
		//here we initialize pthreads 
		fprintf(stderr, "\nThreads used (compute): %d\n", nthreads);
        aux.actual_chunk_size = all_chunk_size[0];
        fprintf(stderr, "\naux2->actual_chunk_size: %d\n", aux.actual_chunk_size);
		nreads = aux.actual_chunk_size/ (readLen) + 10;
		fprintf(stderr, "Info: projected #read in a task: %ld\n", (long)nreads);
	
		/* All memory allocation */
		memoryAlloc(&aux, w2, nreads , nthreads);
        w2.ref_string = aux.ref_string;
        w2.fmi = aux.fmi;
        w2.nreads  = nreads;
		ktp_data_t *ret = (ktp_data_t *) calloc(1, sizeof(ktp_data_t));


       // here we loop until there's nothing to read
       // in the offset and size file
       //   before_local_mapping = MPI_Wtime();
		before_local_mapping = MPI_Wtime();       
		//we loop the chunck_count
		//for (u1 = 0; u1 < chunk_count; u1++){
        u1 = rank_num;  
        while ( u1 < total_chunks){
        
			offset_chunk = all_begin_offset_chunk[u1];
			size_chunk   = all_chunk_size[u1];
			assert(size_chunk != 0);
			/*
			 * Read sequence datas ...
			 *
			 */
			bef = MPI_Wtime();

			buffer_r1 = malloc(size_chunk+1);
			assert(buffer_r1 != NULL);
			buffer_r1[size_chunk] = 0;
		
            buffer_r2 = malloc(size_chunk+1);
            assert(buffer_r2 != NULL);
            buffer_r2[size_chunk]=0;
            
            struct struct_pread_fastq *td_pread1;
            td_pread1 = malloc (NUM_THREADS * sizeof(struct struct_pread_fastq));
            bef = MPI_Wtime();
            pthread_attr_t attr4;
            pthread_attr_init(&attr4);
            pthread_attr_setstacksize(&attr4, BIG_STACK);
            pthread_attr_setdetachstate(&attr4, 0);

            for( n = 0; n < NUM_THREADS; n++ ){
                td_pread1[n].total_thread = NUM_THREADS;
                td_pread1[n].thread_id = n;
                td_pread1[n].job_rank = rank_num;
                td_pread1[n].offset= offset_chunk;
                td_pread1[n].size = size_chunk;
                td_pread1[n].buffer = buffer_r1;
                td_pread1[n].fd  = fh_r1;
                int ret_code = pthread_create(&threads[n], &attr4, pread_fastq_chunck, (void *)(&td_pread1[n]));
                if (ret_code) {
                    fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", ret_code);
                }
           }
           for(n=0; n<NUM_THREADS; n++)
                pthread_join(threads[n], (void *)(&td_pread1[n]));

            pthread_attr_destroy(&attr4);
            free(td_pread1);

            assert(strlen(buffer_r1) == size_chunk);
            assert(*buffer_r1 == '@');

            struct struct_pread_fastq *td_pread2;
            td_pread2 = malloc (NUM_THREADS * sizeof(struct struct_pread_fastq));
            bef = MPI_Wtime();
            pthread_attr_t attr5;
            pthread_attr_init(&attr5);
            pthread_attr_setstacksize(&attr5, BIG_STACK);
            pthread_attr_setdetachstate(&attr5, 0);

            for( n = 0; n < NUM_THREADS; n++ ){
                td_pread2[n].total_thread = NUM_THREADS;
                td_pread2[n].thread_id = n;
                td_pread2[n].job_rank = rank_num;
                td_pread2[n].offset= offset_chunk;
                td_pread2[n].size = size_chunk;
                td_pread2[n].buffer = buffer_r2;
                td_pread2[n].fd  = fh_r2;

                int ret_code = pthread_create(&threads[n], &attr5, pread_fastq_chunck, (void *)(&td_pread2[n]));
                if (ret_code) {
                    fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", ret_code);
                }
           }
           for(n=0; n<NUM_THREADS; n++)
                pthread_join(threads[n], (void *)(&td_pread2[n]));

            pthread_attr_destroy(&attr5);
            free(td_pread2);

            assert(strlen(buffer_r2) == size_chunk);
            assert(*buffer_r2 == '@');

            if (u1 < total_chunks){
                MPI_File_seek(fh_r1, (MPI_Offset)all_begin_offset_chunk[u1 + proc_num], MPI_SEEK_CUR ) ;
                MPI_File_seek(fh_r2, (MPI_Offset)all_begin_offset_chunk[u1 + proc_num], MPI_SEEK_CUR ) ;
            }

			aft = MPI_Wtime();
			fprintf(stderr, "%s: read sequences (%.02f)\n", __func__, aft - bef);
			total_time_reading_seq += (aft - bef);

			bef = MPI_Wtime();
			reads_r1 = all_reads_in_chunk[u1];
			
			if (file_r2 != NULL) reads_r2 = all_reads_in_chunk[u1];

			reads = reads_r1 + reads_r2; bases = 0;
			total_local_reads_aligned += reads/2;
			assert(reads <= INT_MAX);
			fprintf(stderr, "%s: num_rank = %d :: number of paired reads = %zu \n", __func__, rank_num, reads);

			/* Parse sequences ... */
			seqs = malloc(reads * sizeof(*seqs));
			assert(seqs != NULL);

			//case one we are paired
			if (file_r1 != NULL && file_r2 !=NULL) {
				p1 = q1 = buffer_r1; e1 = buffer_r1 + size_chunk; line_number = 0;
				p2 = q2 = buffer_r2; e2 = buffer_r2 + size_chunk; 
				
				while (q1 < e1) {
					if (*q1 != '\n') { q1++; q2++; continue; }
					/* We have a full line ... process it */
					*q1 = '\0'; 
					*q2 = '\0';
					n = files * (line_number / 4);
				
					switch (line_number % 4) {
					case 0: /* Line1: Name and Comment */
						assert(*p1 == '@');
						seqs[n].name   = p1 + 1;
						seqs[n+1].name = p2 + 1;
					
						while (*p1 && !isspace((unsigned char)*p1)) {p1++;p2++;}
						if (*(p1-2) == '/' && isdigit((unsigned char)*(p1-1))) {*(p1-2) = '\0'; *(p2-2) = '\0';}
						if (*p1) {*p1++ = '\0'; *p2++ = '\0';}
						seqs[n].comment = (copy_comment != 0) ? p1 : NULL;
						seqs[n].sam = NULL;
						seqs[n+1].comment = (copy_comment != 0) ? p2 : NULL;
						seqs[n+1].sam = NULL;					
						break;
					case 1: /* Line2: Sequence */
						seqs[n].seq = p1;
						seqs[n].l_seq = q1 - p1;
						seqs[n+1].seq = p2;
						seqs[n+1].l_seq = q2 - p2;

						bases += seqs[n].l_seq;
						bases += seqs[n+1].l_seq;
						break;
					case 2: /* Line3: Ignored */
						assert(*p1 == '+');
						break;
					case 3: /* Line4: Quality */
						seqs[n].qual = p1;
						seqs[n+1].qual = p2;
						break; }
					p1 = ++q1; 
					p2 = ++q2; 
					line_number++; 
				}
				
			}
		
			aft = MPI_Wtime();
			fprintf(stderr, "%s: parsed sequences (%.02f)\n", __func__, aft - bef);
			total_time_parsing += (aft - bef);
			fprintf(stderr, "rank %d ::: [M::%s] read %zu sequences (%ld bp)...\n",rank_num , __func__, reads, (long)bases);
			/* Datas computation ... */
			bef = MPI_Wtime();
			
			ret->seqs   = seqs;
            		ret->n_seqs = reads;				

            //int64_t size = 0;
            //for (int i = 0; i < ret->n_seqs; ++i) size += ret->seqs[i].l_seq;
            //fprintf(stderr, "\t[0000][ M::%s] read %d sequences (%ld bp)...\n",
            //        __func__, ret->n_seqs, (long)size);
	
			fprintf(stderr, "[RANK %d] 1: Calling process()\n", rank_num);

			//fprintf(stderr, "Reallocating initial memory allocations!!\n");
			if (w2.regs) free(w2.regs); 
			if (w2.chain_ar) free(w2.chain_ar); 
			if (w2.seedBuf) free(w2.seedBuf);
			
			nreads = ret->n_seqs;
			w2.regs = (mem_alnreg_v *) calloc(nreads, sizeof(mem_alnreg_v));
			w2.chain_ar = (mem_chain_v*) malloc (nreads * sizeof(mem_chain_v));
			w2.seedBuf = (mem_seed_t *) calloc(sizeof(mem_seed_t),nreads * AVG_SEEDS_PER_READ);
			
            /* new call */
            w2.nreads=ret->n_seqs;
            mem_process_seqs(opt,
                             aux.n_processed,
                             ret->n_seqs,
                             ret->seqs,
                             aux.pes0,
                             w2);
            aux.n_processed += ret->n_seqs;            

			aft = MPI_Wtime();
			fprintf(stderr, "%s: computed mappings (%.02f)\n", __func__, aft - bef);

            if (dofixmate){
                bef = MPI_Wtime();
                int ret_code = 0;
                        struct thread_data *td;
                        td = malloc (NUM_THREADS * sizeof(struct thread_data));
                        bef = MPI_Wtime();
                        pthread_attr_init(&attr);
                        pthread_attr_setstacksize(&attr, SMALL_STACK);
                        pthread_attr_setdetachstate(&attr, 0);

                        for( n = 0; n < NUM_THREADS; n++ ){
                            td[n].total_thread = NUM_THREADS;
                                td[n].thread_id = n;
                                assert(seqs);
                                td[n].seqs_thr = seqs;
                                td[n].job_rank = rank_num;
                                td[n].total_reads = reads;
                                td[n].start_index_seqs = 0;
                                td[n].final_index_seqs = reads-1;
                                td[n].indix_thr = &indix;
                                td[n].total_lines = 0;
                                ret_code = pthread_create(&threads[n], &attr, call_fixmate, (void *)(&td[n]));
                                if (ret_code) {
                                    fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", ret_code);
                                }
                        }
                        for(n=0; n<NUM_THREADS; n++)
                                pthread_join(threads[n], (void *)(&td[n]));

                pthread_attr_destroy(&attr);
                        free(td);
                aft = MPI_Wtime();
                        fprintf(stderr, "rank: %d :: %s: time spend in fixmate (%.02f) \n", rank_num, __func__, aft - bef);

            }

            // write SAM format
            if ( write_format == 2){
            	bef = MPI_Wtime();
                pthread_attr_t attr;
                pthread_attr_init(&attr);
                pthread_attr_setstacksize(&attr, BIG_STACK);
                pthread_attr_setdetachstate(&attr, 0);

                struct struct_data_thread *td;
                td = calloc (opt->n_threads, sizeof(struct struct_data_thread));

                int rest = reads%NUM_THREADS;
                int quot = reads/NUM_THREADS;
                pthread_t threads[NUM_THREADS];

                int ret_code = 0;
                for ( n = 0; n < NUM_THREADS; n++){
                    td[n].seqs_thr = seqs;
                    td[n].begin_index = n * quot;
                    td[n].end_index = ((n+1) * quot);
                    if ( n == (NUM_THREADS - 1)) td[n].end_index = reads;
                    td[n].file_desc = fh_out;
                    ret_code = pthread_create(&threads[n], &attr, write_sam_mt, (void *)(&td[n]));
                }

                for (n = 0; n < NUM_THREADS; n++)
                    pthread_join(threads[n], (void *)(&td[n]));

                pthread_attr_destroy(&attr);
                free(td);

                aft = MPI_Wtime();
                total_time_writing += (aft - bef);
                free(buffer_r1);
                free(buffer_r2);

                }
                //write BGZF 
		if ( write_format == 0 ) {

                        int ret_code = 0;
                        struct thread_data_compress *tdc;
                        tdc = malloc (NUM_THREADS * sizeof(struct thread_data_compress));
                        bef = MPI_Wtime();
                        pthread_attr_init(&attr);
                        pthread_attr_setstacksize(&attr, SMALL_STACK);
                        pthread_attr_setdetachstate(&attr, 0);

                        for( n = 0; n < NUM_THREADS; n++ ){
                            tdc[n].total_thread = NUM_THREADS;
                                tdc[n].thread_id = n;
                                tdc[n].seqs_thr = seqs;
                                tdc[n].job_rank = rank_num;
                                tdc[n].total_reads = reads;
                                tdc[n].thr_comp_sz = 0;
                                tdc[n].comp_level = compression_level;
                                tdc[n].compressed_buffer_thread  = 0;
                                tdc[n].fh_out = fh_out;
                                ret_code = pthread_create(&threads[n], &attr, compress_and_write_bgzf_thread, (void *)(&tdc[n]));
                                if (ret_code) {
                                    fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", ret_code);
                                }
                        }
                        for(n=0; n<NUM_THREADS; n++)
                                pthread_join(threads[n], (void *)(&tdc[n]));
                        
                        /*
                        size_t compressed_size = 0;
                        for(n=0; n<NUM_THREADS; n++)
                                compressed_size += tdc[n].thr_comp_sz;

                        uint8_t *compressed_buff;
                        compressed_buff =  malloc( compressed_size*sizeof(uint8_t) + 1);
                        compressed_buff[compressed_size] = 0;
                        size_t total_cpy = 0;
                        for(n=0; n<NUM_THREADS; n++){
                                memcpy(compressed_buff + total_cpy, tdc[n].compressed_buffer_thread, tdc[n].thr_comp_sz);
                                total_cpy += tdc[n].thr_comp_sz;
                                free(tdc[n].compressed_buffer_thread);
                        }


                        aft = MPI_Wtime();
                        fprintf(stderr, "%s: compress time (%.02f)\n", __func__, aft - bef);


                        size_t compSize = compressed_size;
                        bef = MPI_Wtime();
                        res = MPI_File_write_shared(fh_out, compressed_buff, compSize, MPI_BYTE, &status);
                        assert(res == MPI_SUCCESS);
                        res = MPI_Get_count(&status, MPI_BYTE, &count);
                        assert(res == MPI_SUCCESS);
                        assert(count == (int)compSize);
                        */
                        pthread_attr_destroy(&attr);
                        free(tdc);
                        for (n = 0; n < reads; n++) free(seqs[n].sam);
                            free(seqs);
			aft = MPI_Wtime();
                        total_time_writing += (aft - bef);
                        free(buffer_r1);
                        free(buffer_r2);
                        
                    }

                    // write BAM
                    if ( write_format == 1 ) {

                        int ret_code = 0;
                        struct thread_data_compress *tdc;
                        tdc = malloc (NUM_THREADS * sizeof(struct thread_data_compress));
                        bef = MPI_Wtime();
                        pthread_attr_init(&attr);
                        pthread_attr_setstacksize(&attr, SMALL_STACK);
                        pthread_attr_setdetachstate(&attr, 0);

                        for( n = 0; n < NUM_THREADS; n++ ){
                            tdc[n].total_thread = NUM_THREADS;
                                tdc[n].thread_id = n;
                                tdc[n].seqs_thr = seqs;
                                tdc[n].job_rank = rank_num;
                                tdc[n].total_reads = reads;
                                tdc[n].thr_comp_sz = 0;
                                tdc[n].comp_level = compression_level;
                                tdc[n].compressed_buffer_thread  = 0;
                                tdc[n].fh_out = fh_out;
                                ret_code = pthread_create(&threads[n], &attr, compress_and_write_bam_thread, (void *)(&tdc[n]));
                                if (ret_code) {
                                    fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", ret_code);
                                }
                        }
                        for(n=0; n<NUM_THREADS; n++)
                                pthread_join(threads[n], (void *)(&tdc[n]));


                        pthread_attr_destroy(&attr);
                        free(tdc);
                        for (n = 0; n < reads; n++) free(seqs[n].sam);
                            free(seqs);
			aft = MPI_Wtime();
                        total_time_writing += (aft - bef);
                        free(buffer_r1);
                        free(buffer_r2);

                    }
                    	
                    	fprintf(stderr, "%s: wrote results (%.02f)\n", __func__, aft - bef);

            		u1 += proc_num;
		} //end for (u1 = 0; u1 < chunk_count; u1++){

        MPI_Barrier(MPI_COMM_WORLD);
        if ( (write_format == 1) && (rank_num == 0)) {

            res = MPI_File_write_shared(fh_out, magic, 28, MPI_BYTE, &status);
            assert(res == MPI_SUCCESS);
            res = MPI_Get_count(&status, MPI_BYTE, &count);
            assert(res == MPI_SUCCESS);
            assert(count == 28);
        }

        free(w2.chain_ar);
        free(w2.regs);
        free(w2.seedBuf);

        for(int l=0; l<nthreads; l++) {
            _mm_free(w2.mmc.seqBufLeftRef[l*CACHE_LINE]);
            _mm_free(w2.mmc.seqBufRightRef[l*CACHE_LINE]);
            _mm_free(w2.mmc.seqBufLeftQer[l*CACHE_LINE]);
            _mm_free(w2.mmc.seqBufRightQer[l*CACHE_LINE]);
        }

        for(int l=0; l<nthreads; l++) {
            free(w2.mmc.seqPairArrayAux[l]);
            free(w2.mmc.seqPairArrayLeft128[l]);
            free(w2.mmc.seqPairArrayRight128[l]);
        }

        for(int l=0; l<nthreads; l++) {
            _mm_free(w2.mmc.matchArray[l]);
            free(w2.mmc.min_intv_ar[l]);
            free(w2.mmc.query_pos_ar[l]);
            free(w2.mmc.enc_qdb[l]);
            free(w2.mmc.rid[l]);
            _mm_free(w2.mmc.lim[l]);
        }

        res = MPI_Win_free(&win_shr_ref);
        assert(res == MPI_SUCCESS);
        res = MPI_Win_free(&win_shr_pac);
        assert(res == MPI_SUCCESS);

        /* shall we delete?
         * because error in read_index_ele.cpp:52
        */
        //delete(aux.fmi);


	} //end if (file_r2 != NULL && stat_r1.st_size == stat_r2.st_size)
	
	
	if (file_r1 != NULL && file_r2 != NULL && stat_r1.st_size != stat_r2.st_size) {
	   
	   /*
		*	We are in the case the reads are paired and trimmed
		*/

		aft = 0; aft++;
		bef = 0; bef++;

        /*
         *  Rank 0 estimate the size of a read
         *
         */
        size_t slen, blen;
        off_t tmp_sz = 1024;
        if (rank_num == 0){
            int fd_tmp = open(file_r1, O_RDONLY, 0666);
            char *buffer = malloc(tmp_sz  + 1);
            buffer[tmp_sz] = '\0';
            size_t read_out = read(fd_tmp, buffer, tmp_sz);
            assert(read_out);
            assert(strlen(buffer) == tmp_sz);
            assert( *buffer == '@');
            s = buffer;
            e = buffer + tmp_sz;
            p = q = s;
            while (q < e && *q != '\n') q++; p = ++q;
            while (q < e && *q != '\n') q++; blen = q - p; p = ++q;
            while (q < e && *q != '\n') q++; p = ++q;
            while (q < e && *q != '\n') q++; slen = q - s + 1;
             free(buffer);
        }

        //Rank O broadcast the size of a read
        res = MPI_Bcast(&blen, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        assert(res == MPI_SUCCESS);
        
		//global offset contains the starting offset in the fastq for each process
		//TODO?: use only one vector for all the global offsets iot reduce communications
        size_t *goff    = NULL;
        size_t *goff2   = NULL;
        goff    = malloc((proc_num * NUM_THREADS + 1) * sizeof(size_t));
        goff2   = malloc((proc_num * NUM_THREADS + 1) * sizeof(size_t));

        size_t *goff_inter = calloc( (proc_num * NUM_THREADS + 1) , sizeof(size_t));
        size_t *goff_inter_2 = calloc( (proc_num * NUM_THREADS + 1) , sizeof(size_t));

		//this function is used to fill the goff vectors
        bef = MPI_Wtime();
        find_process_starting_offset_mt(goff, stat_r1.st_size, file_r1, proc_num, rank_num, NUM_THREADS);
        find_process_starting_offset_mt(goff2, stat_r2.st_size, file_r2, proc_num, rank_num, NUM_THREADS);
        aft = MPI_Wtime();
        fprintf(stderr, "%s: rank %d time spend in finding process start offset = (%.02f) \n", __func__, rank_num, aft - bef);

		//now we exchange the goff buffer between all proc
		//rank 0 gather the vector
        int i12=0;
        for ( i12 = 0; i12 < proc_num * NUM_THREADS + 1; i12++ ){
            goff_inter[i12] = goff[i12];
            goff_inter_2[i12] = goff2[i12];
        }

        res = MPI_Allgather(&goff_inter[rank_num*NUM_THREADS], NUM_THREADS, MPI_LONG_LONG_INT, goff , NUM_THREADS, MPI_LONG_LONG_INT, MPI_COMM_WORLD);
        assert(res == MPI_SUCCESS);

        res = MPI_Allgather(&goff_inter_2[rank_num*NUM_THREADS], NUM_THREADS, MPI_LONG_LONG_INT, goff2 , NUM_THREADS, MPI_LONG_LONG_INT, MPI_COMM_WORLD);
        assert(res == MPI_SUCCESS);

        free(goff_inter);
        free(goff_inter_2);

		///resources needed to find the offsets and size of each read.	
		size_t grand_total_num_reads 	= 0; 
		size_t grand_total_num_reads_2	= 0;
		size_t local_num_reads 		= 0;
		size_t local_num_reads_2	= 0;
		size_t total_num_reads 		= 0;
		size_t total_num_reads_2	= 0;

		//Here I decided to keep separated vectors because otherwise with all the reallocs they would be too big and take too much contiguous memory space
		int *local_read_size    	    = calloc(1, sizeof(int));
		int *local_read_size_2    	    = calloc(1, sizeof(int));
		size_t *local_read_bytes    	= calloc(1, sizeof(size_t));
		size_t *local_read_bytes_2    	= calloc(1, sizeof(size_t));
		size_t *local_read_offsets 	    = calloc(1, sizeof(size_t));
		size_t *local_read_offsets_2 	= calloc(1, sizeof(size_t));

		assert( local_read_offsets 	!= NULL);
		assert( local_read_size 	!= NULL);
		assert( local_read_bytes 	!= NULL);
		assert( local_read_offsets_2 	!= NULL);
		assert( local_read_size_2 	!= NULL);
		assert( local_read_bytes_2 	!= NULL);
	
		///find offsets and sizes for the first file
		bef = MPI_Wtime();
        pthread_attr_t attr;
        pthread_attr_init(&attr);
        pthread_attr_setstacksize(&attr, SMALL_STACK);
        pthread_attr_setdetachstate(&attr, 0);
        pthread_t threads_1[NUM_THREADS];

        struct struct_data_thread_1 *td_1;
        size_t *local_num_reads_t            = calloc(NUM_THREADS, sizeof(size_t));
        size_t *total_num_reads_t            = calloc(NUM_THREADS, sizeof(size_t));
        size_t **local_read_offsets_t        = calloc(NUM_THREADS, sizeof(size_t*));
        size_t **local_read_bytes_t          = calloc(NUM_THREADS, sizeof(size_t*));
        int **local_read_size_t              = calloc(NUM_THREADS, sizeof(int*));

        td_1 = calloc(NUM_THREADS, sizeof(struct struct_data_thread_1));
        int ret_code_1 = 0;
        int goff_idx = 0;
        for ( n = 0; n < NUM_THREADS; n++){
            goff_idx = (rank_num * NUM_THREADS) + n;
            td_1[n].offset_in_file_mt         = goff[goff_idx];
            td_1[n].size2read_mt              = goff[goff_idx + 1] - goff[goff_idx];
            td_1[n].file_r1_mt                = file_r1;
            td_1[n].local_num_reads_mt        = &local_num_reads_t[n];
            td_1[n].total_num_reads_mt        = &total_num_reads_t[n];
            td_1[n].local_read_offsets_mt     = &local_read_offsets_t[n];
            td_1[n].local_read_size_mt        = &local_read_size_t[n];
            td_1[n].local_read_bytes_mt       = &local_read_bytes_t[n];
            td_1[n].proc_num_mt               = proc_num;
            td_1[n].rank_num_mt               = rank_num;
            td_1[n].thread_num_mt             = n;
            td_1[n].previous_read_num         = 0;
            ret_code_1 = pthread_create(&threads_1[n], &attr, find_reads_size_and_offsets_mt, (void *)(&td_1[n]));

        }

        total_num_reads = 0;
        for (n = 0; n < NUM_THREADS; n++){
            pthread_join(threads_1[n], (void *)(&td_1[n]));
            total_num_reads += *(td_1[n].total_num_reads_mt);
        }

        local_read_offsets  = calloc(total_num_reads, sizeof(size_t));
        local_read_size     = calloc(total_num_reads, sizeof(int));
        local_read_bytes    = calloc(total_num_reads, sizeof(size_t));

        assert(local_read_offsets);
        assert(local_read_size);
        assert(local_read_bytes);

        size_t tmp_var = 0;
        for (n = 0; n < NUM_THREADS; n++){
            td_1[n].local_read_offsets     = local_read_offsets;
            td_1[n].local_read_size        = local_read_size;
            td_1[n].local_read_bytes       = local_read_bytes;
            td_1[n].previous_read_num      = tmp_var;
            tmp_var                        += *(td_1[n].total_num_reads_mt);
            ret_code_1 = pthread_create(&threads_1[n], &attr, copy_local_read_info_mt, (void *)(&td_1[n]));
        }

        for (n = 0; n < NUM_THREADS; n++){
            pthread_join(threads_1[n], (void *)(&td_1[n]));
            free(local_read_offsets_t[n]);
            free(local_read_bytes_t[n]);
            free(local_read_size_t[n]);
        }

        free(local_num_reads_t);
        free(total_num_reads_t);
        free(local_read_offsets_t);
        free(local_read_bytes_t);
        free(local_read_size_t);

        pthread_attr_destroy(&attr);
        free(td_1);

        res = MPI_Reduce(&total_num_reads, &grand_total_num_reads, 1, MPI_LONG_LONG_INT, MPI_SUM, 0,MPI_COMM_WORLD);
        assert ( res == MPI_SUCCESS);

        res = MPI_Bcast(&grand_total_num_reads, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        assert ( res == MPI_SUCCESS);

        if ( rank_num == 0 )
            fprintf(stderr, "rank %d ::: total_num_reads (1rst file)= %zu \n", rank_num, grand_total_num_reads);

        MPI_Barrier(MPI_COMM_WORLD);

        pthread_attr_t attr2;
        pthread_attr_init(&attr2);
        pthread_attr_setstacksize(&attr2, SMALL_STACK);
        pthread_attr_setdetachstate(&attr2, 0);

        pthread_t threads_2[NUM_THREADS];

        struct struct_data_thread_1 *td_2;
        size_t *local_num_reads_t2            = calloc(NUM_THREADS, sizeof(size_t));
        size_t *total_num_reads_t2            = calloc(NUM_THREADS, sizeof(size_t));
        size_t **local_read_offsets_t2        = calloc(NUM_THREADS, sizeof(size_t*));
        size_t **local_read_bytes_t2          = calloc(NUM_THREADS, sizeof(size_t*));
        int **local_read_size_t2              = calloc(NUM_THREADS, sizeof(int*));

        td_2 = calloc(NUM_THREADS, sizeof(struct struct_data_thread_1));

        int ret_code_2 = 0;
        int goff_idx2 = 0;
        for ( n = 0; n < NUM_THREADS; n++){
            goff_idx2 = (rank_num * NUM_THREADS) + n;
            td_2[n].offset_in_file_mt         = goff2[goff_idx2];
            td_2[n].size2read_mt              = goff2[goff_idx2 + 1] - goff2[goff_idx2];
            td_2[n].file_r1_mt                = file_r2;
            td_2[n].local_num_reads_mt        = &local_num_reads_t2[n];
            td_2[n].total_num_reads_mt        = &total_num_reads_t2[n];
            td_2[n].local_read_offsets_mt     = &local_read_offsets_t2[n];
            td_2[n].local_read_size_mt        = &local_read_size_t2[n];
            td_2[n].local_read_bytes_mt       = &local_read_bytes_t2[n];
            td_2[n].proc_num_mt               = proc_num;
            td_2[n].rank_num_mt               = rank_num;
            td_2[n].thread_num_mt             = n;
            td_2[n].previous_read_num         = 0;
            ret_code_2 = pthread_create(&threads_2[n], &attr2, find_reads_size_and_offsets_mt, (void *)(&td_2[n]));
        }
        total_num_reads_2 = 0;
        for (n = 0; n < NUM_THREADS; n++){
            pthread_join(threads_1[n], (void *)(&td_2[n]));
            total_num_reads_2 += *(td_2[n].total_num_reads_mt);
        }

        local_read_offsets_2  = calloc(total_num_reads_2, sizeof(size_t));
        local_read_size_2     = calloc(total_num_reads_2, sizeof(int));
        local_read_bytes_2    = calloc(total_num_reads_2, sizeof(size_t));

        assert(local_read_offsets_2);
        assert(local_read_size_2);
        assert(local_read_bytes_2);

        size_t tmp_var2 = 0;
        for (n = 0; n < NUM_THREADS; n++){
            td_2[n].local_read_offsets     = local_read_offsets_2;
            td_2[n].local_read_size        = local_read_size_2;
            td_2[n].local_read_bytes       = local_read_bytes_2;
            td_2[n].previous_read_num      = tmp_var2;
            tmp_var2                       += *(td_2[n].total_num_reads_mt);

            ret_code_2 = pthread_create(&threads_2[n], &attr2, copy_local_read_info_mt, (void *)(&td_2[n]));
        }

        for (n = 0; n < NUM_THREADS; n++){
            pthread_join(threads_2[n], (void *)(&td_2[n]));
            free(local_read_offsets_t2[n]);
            free(local_read_bytes_t2[n]);
            free(local_read_size_t2[n]);
        }

        free(local_num_reads_t2);
        free(total_num_reads_t2);

        free(local_read_offsets_t2);
        free(local_read_bytes_t2);
        free(local_read_size_t2);

        pthread_attr_destroy(&attr2);
        free(td_2);


        res = MPI_Reduce(&total_num_reads_2, &grand_total_num_reads_2, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        assert ( res == MPI_SUCCESS);

        res = MPI_Bcast(&grand_total_num_reads_2, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        assert ( res == MPI_SUCCESS);
        
        if ( rank_num == 0 )
            fprintf(stderr, "rank %d ::: total_num_reads (for 2nd file) = %zu \n", rank_num, grand_total_num_reads_2);

        assert ( grand_total_num_reads_2 == grand_total_num_reads );

		aft = MPI_Wtime();
		fprintf(stderr, "%s: rank %d num reads parsed: %zu ::: time spend reading and parsing entire fastq = (%.02f) \n", __func__, rank_num, total_num_reads, aft - bef);

		if (goff) free(goff);
        if (goff2) free (goff2);

		///Now find the required information to create the chunks, such as total number of bases
		local_num_reads 	= total_num_reads;
		local_num_reads_2 	= total_num_reads_2;
	
		//change chunck_size if necessary
		//opt->chunk_size=10000000;

		size_t length_sum	= 0;
		size_t length_sum_2	= 0;
		size_t i;
		//since all the sizes are not the same anymore, you have to loop

		size_t min_num_read = min(local_num_reads, local_num_reads_2);
		size_t bases_tmp 	= 0;
		size_t chunck_num 	= 0;
		size_t chunck_num_2	= 0;

		for(i=0; i < min_num_read; i++)   {

			bases_tmp  += (local_read_size[i] + local_read_size_2[i]);

			if ( bases_tmp > ( opt->chunk_size * opt->n_threads)){

				bases_tmp = 0;
				chunck_num++;
				chunck_num_2++;
			}
		}

		chunck_num 	+= 2;
		chunck_num_2 	+= 2;
	
		fprintf(stderr,"Rank %d :: chunk_num = %zu\n", rank_num, chunck_num);
		fprintf(stderr,"Rank %d :: chunk_num_2 = %zu\n", rank_num, chunck_num_2);

		///assert that the sizes and offsets found earlier are real.
		//changed the comparisons here for the size. they were compared to blen which has become impossible
		size_t h=0;
		/*
		for (h=0; h < total_num_reads; h++)
		{
			assert(local_read_size[h] > 0 );
			assert(local_read_offsets[h] >= 0 );
		}
		for(h=0; h < total_num_reads_2; h++)
		{
			assert(local_read_size_2[h] > 0 );
			assert(local_read_offsets_2[h] >= 0 );
		}
		*/
		///Now that we know how many chunks we need, we create the vectors associated
		///and we update their informations according to the file they are treating

		// we allocate vector for chunks offset
		size_t *begin_offset_chunk 	= calloc(chunck_num, sizeof(size_t));
		size_t *chunk_size 			= calloc(chunck_num, sizeof(size_t));
		size_t *reads_in_chunk 		= calloc(chunck_num, sizeof(size_t));

		assert( begin_offset_chunk != NULL );
		assert( chunk_size != NULL );
		assert( reads_in_chunk != NULL );
	
		size_t chunk_count = 0;

		size_t *begin_offset_chunk_2 	= calloc(chunck_num_2, sizeof(size_t));
		size_t *chunk_size_2 		    = calloc(chunck_num_2, sizeof(size_t));
		size_t *reads_in_chunk_2 	    = calloc(chunck_num_2, sizeof(size_t));

		assert( begin_offset_chunk_2 != NULL );
		assert( chunk_size_2 != NULL );
		assert( reads_in_chunk_2 != NULL );

		size_t chunk_count_2 = 0;

		maxsiz = fixed_chunk_size; 
		
		fprintf(stderr,"rank %d ::: Call find_chunks_info \n", rank_num);
		// the detail of he paramters is at the function definition
		// fprintf(stderr, "rank %d ::: begin_offset_chunk = %zu \n", rank_num, begin_offset_chunk[0]);
		// fprintf(stderr, "rank %d ::: begin_offset_chunk_2 = %zu \n", rank_num, begin_offset_chunk_2[0]);
		// fprintf(stderr, "rank %d ::: chunk_size = %zu \n", rank_num, chunk_size);
		// fprintf(stderr, "rank %d ::: chunk_size_2 = %zu \n", rank_num, chunk_size_2);

		bef = MPI_Wtime();
		find_chunks_info_trim(begin_offset_chunk,
					 begin_offset_chunk_2,
					 chunk_size,
					 chunk_size_2,
					 reads_in_chunk,
					 reads_in_chunk_2,
					 local_read_size,
					 local_read_size_2,
					 local_read_bytes,
					 local_read_bytes_2,
					 local_read_offsets,
					 local_read_offsets_2,
					 rank_num,
					 proc_num,
					 local_num_reads,
					 local_num_reads_2,
					 grand_total_num_reads,
					 grand_total_num_reads_2,
					 maxsiz,
					 &chunk_count);

		aft = MPI_Wtime();
		fprintf(stderr, "%s ::: rank %d ::: evaluating offsets chuncks and sizes (%.02f) found %zu chuncks \n", __func__, rank_num, aft - bef, chunk_count);
	
		MPI_Barrier(MPI_COMM_WORLD);

		if (local_read_offsets) 	free(local_read_offsets);
		if (local_read_size) 		free(local_read_size);
		if (local_read_offsets_2) 	free(local_read_offsets_2);
		if (local_read_size_2)		free(local_read_size_2);

		//verify the number of reads

		size_t num_reads_1 = 0;
		size_t num_reads_2 = 0;
		size_t total_num_reads_v1 = 0;
		size_t total_num_reads_v2 = 0;

		for(h=0; h < chunk_count; h++){
			 assert( reads_in_chunk[h] == reads_in_chunk_2[h]);
			 num_reads_1 += reads_in_chunk[h];
			 num_reads_2 += reads_in_chunk_2[h];
		}

		res = MPI_Reduce(&num_reads_1, &total_num_reads_v1, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		assert ( res == MPI_SUCCESS);

		res = MPI_Bcast(&total_num_reads_v1, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
		assert ( res == MPI_SUCCESS);

		if ( rank_num == 0 )
			fprintf(stderr, "rank %d ::: after finding chunks info ::: total_num_reads forward = %zu \n", rank_num, total_num_reads_v1);

		res = MPI_Reduce(&num_reads_2, &total_num_reads_v2, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		assert ( res == MPI_SUCCESS);

		res = MPI_Bcast(&total_num_reads_v2, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
		assert ( res == MPI_SUCCESS);

		if ( rank_num == 0 )
			fprintf(stderr, "rank %d ::: after finding chunks info ::: total_num_reads backward = %zu \n", rank_num, total_num_reads_v2);
		
		assert( total_num_reads_v1 == total_num_reads_v2 );
		assert( grand_total_num_reads == total_num_reads_v1);
		assert( grand_total_num_reads_2 == total_num_reads_v2);
		MPI_Barrier(MPI_COMM_WORLD);

		//display test to check that each vector was well filled with the right information
		//fprintf(stderr,"rank %d ::: begin offset %zu, %zu\n", rank_num, begin_offset_chunk[0], begin_offset_chunk_2[0]);
		//fprintf(stderr,"rank %d ::: chunk size %zu, %zu\n", rank_num, chunk_size[0], chunk_size_2[0]);
		//fprintf(stderr,"rank %d ::: reads in chunk %zu, %zu\n", rank_num, reads_in_chunk[0], reads_in_chunk_2[0]);


        size_t total_chunks = 0;

        res = MPI_Reduce(&chunk_count, &total_chunks, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        assert ( res == MPI_SUCCESS);
        res = MPI_Bcast(&total_chunks, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        assert ( res == MPI_SUCCESS);

        size_t chunks_per_rank[proc_num];
        res = MPI_Allgather(&chunk_count, 1, MPI_LONG_LONG_INT, chunks_per_rank, 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);

        all_begin_offset_chunk = malloc(total_chunks * sizeof(size_t));
        all_begin_offset_chunk[0]=0;
        all_chunk_size = malloc(total_chunks * sizeof(size_t));
        all_chunk_size[0] = 0;
        all_reads_in_chunk =  malloc(total_chunks * sizeof(size_t));
        all_reads_in_chunk[0] = 0;

        all_begin_offset_chunk_2 = malloc(total_chunks * sizeof(size_t));
        all_begin_offset_chunk_2[0]=0;
        all_chunk_size_2 = malloc(total_chunks * sizeof(size_t));
        all_chunk_size_2[0] = 0;
        all_reads_in_chunk_2 =  malloc(total_chunks * sizeof(size_t));
        all_reads_in_chunk_2[0] = 0;

        int displ_chunk[proc_num];
        displ_chunk[0] = 0;

        for (i = 1; i < proc_num; i++) displ_chunk[i] = (displ_chunk[i-1] + chunks_per_rank[i-1]);

        int indx=displ_chunk[rank_num];
        for (i = 0; i <  chunks_per_rank[rank_num]; i++) {
            all_chunk_size[indx+i]=chunk_size[i];
            all_begin_offset_chunk[indx+i]=begin_offset_chunk[i];
            all_reads_in_chunk[indx+i]=reads_in_chunk[i];

            all_chunk_size_2[indx+i]=chunk_size_2[i];
            all_begin_offset_chunk_2[indx+i]=begin_offset_chunk_2[i];
            all_reads_in_chunk_2[indx+i]=reads_in_chunk_2[i];

        }

        if (rank_num > 0){
            res=MPI_Send(chunk_size, chunks_per_rank[rank_num], MPI_LONG_LONG_INT, 0, rank_num, MPI_COMM_WORLD);
            assert(res == MPI_SUCCESS);
            }
        else{
            for (i = 1; i < proc_num; i++){
                res=MPI_Recv(&(all_chunk_size[displ_chunk[i]]), chunks_per_rank[i], MPI_LONG_LONG_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                assert(res == MPI_SUCCESS);
            }
        }

        if (rank_num > 0){
            res=MPI_Send(chunk_size_2, chunks_per_rank[rank_num], MPI_LONG_LONG_INT, 0, rank_num, MPI_COMM_WORLD);
            assert(res == MPI_SUCCESS);
            }
        else{
            for (i = 1; i < proc_num; i++){
                res=MPI_Recv(&(all_chunk_size_2[displ_chunk[i]]), chunks_per_rank[i], MPI_LONG_LONG_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                assert(res == MPI_SUCCESS);
            }
        }

        if (rank_num > 0){
            res=MPI_Send(begin_offset_chunk, chunks_per_rank[rank_num], MPI_LONG_LONG_INT, 0, rank_num, MPI_COMM_WORLD);
            assert(res == MPI_SUCCESS);
            }
        else{
            for (i = 1; i < proc_num; i++){
                res=MPI_Recv(&(all_begin_offset_chunk[displ_chunk[i]]), chunks_per_rank[i], MPI_LONG_LONG_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                assert(res == MPI_SUCCESS);
            }
        }

        if (rank_num > 0){
            res=MPI_Send(begin_offset_chunk_2, chunks_per_rank[rank_num], MPI_LONG_LONG_INT, 0, rank_num, MPI_COMM_WORLD);
            assert(res == MPI_SUCCESS);
            }
        else{
            for (i = 1; i < proc_num; i++){
                res=MPI_Recv(&(all_begin_offset_chunk_2[displ_chunk[i]]), chunks_per_rank[i], MPI_LONG_LONG_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                assert(res == MPI_SUCCESS);
            }
        }

        if (rank_num > 0){
            res=MPI_Send(reads_in_chunk, chunks_per_rank[rank_num], MPI_LONG_LONG_INT, 0, rank_num, MPI_COMM_WORLD);
            assert(res == MPI_SUCCESS);
            }
        else{
            for (i = 1; i < proc_num; i++){
                res=MPI_Recv(&(all_reads_in_chunk[displ_chunk[i]]), chunks_per_rank[i], MPI_LONG_LONG_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                assert(res == MPI_SUCCESS);
            }
        }
        
        if (rank_num > 0){
            res=MPI_Send(reads_in_chunk_2, chunks_per_rank[rank_num], MPI_LONG_LONG_INT, 0, rank_num, MPI_COMM_WORLD);
            assert(res == MPI_SUCCESS);
            }
        else{
            for (i = 1; i < proc_num; i++){
                res=MPI_Recv(&(all_reads_in_chunk_2[displ_chunk[i]]), chunks_per_rank[i], MPI_LONG_LONG_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                assert(res == MPI_SUCCESS);
            }
        }
        
        res=MPI_Bcast(all_chunk_size, total_chunks, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        res=MPI_Bcast(all_reads_in_chunk, total_chunks, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        res=MPI_Bcast(all_begin_offset_chunk, total_chunks, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

        res=MPI_Bcast(all_chunk_size_2, total_chunks, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        res=MPI_Bcast(all_reads_in_chunk_2, total_chunks, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        res=MPI_Bcast(all_begin_offset_chunk_2, total_chunks, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

        assert(res == MPI_SUCCESS);

        for ( i = 0; i < total_chunks; i++)
                assert(all_chunk_size[i] != 0);

        free(chunk_size);
        free(begin_offset_chunk);
        free(reads_in_chunk);
        free(chunk_size_2);
        free(begin_offset_chunk_2);
        free(reads_in_chunk_2);

		/// Map reference genome indexes in shared memory (by host)
		
		uint64_t tim = __rdtsc();
        fprintf(stderr, "* Ref file: %s\n", file_ref);
        aux.fmi = new FMI_search(file_ref);
        
        uint8_t *shared_ref;
        MPI_Comm comm_shr;
        int count4, rank_shr = 0;    
        size_t size_read;

        res = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_shr);
        assert(res == MPI_SUCCESS); 

        res = MPI_Comm_rank(comm_shr, &rank_shr);
        assert(res == MPI_SUCCESS);

        char binary_seq_file[200]; 
        sprintf(binary_seq_file, "%s.0123", file_ref);
    
        char pac_file[200];
        sprintf(pac_file, "%s.pac", file_ref);
    
        MPI_File fh_ref_file, fh_pac_file;
        MPI_Aint size_shr_ref, size_shr_pac;
        MPI_Offset size_map_ref,size_map_pac;
        MPI_Win win_shr_ref, win_shr_pac;

        uint8_t *addr_map_ref;
        uint8_t *addr_map_pac;


        MPI_Info win_info;
        MPI_Info_create(&win_info);
        MPI_Info_set(win_info, "alloc_shared_non_contig", "true");

        res = MPI_File_open(MPI_COMM_WORLD, binary_seq_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_ref_file);
        assert(res == MPI_SUCCESS);

        res = MPI_File_open(MPI_COMM_WORLD, pac_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_pac_file);
        assert(res == MPI_SUCCESS);

        res = MPI_File_get_size(fh_pac_file, &size_map_pac);
        assert(res == MPI_SUCCESS);

        res = MPI_File_get_size(fh_ref_file, &size_map_ref);
        assert(res == MPI_SUCCESS);

        size_shr_ref = (rank_shr == 0) ? size_map_ref : 0;
        res = MPI_Win_allocate_shared(size_shr_ref * sizeof(uint8_t), 8 , win_info, comm_shr, &addr_map_ref, &win_shr_ref);

        size_shr_pac = (rank_shr == 0) ? size_map_pac : 0;
        res = MPI_Win_allocate_shared(size_shr_pac * sizeof(uint8_t), 8 , win_info, comm_shr, &addr_map_pac, &win_shr_pac);

        res = MPI_Win_shared_query(win_shr_pac, MPI_PROC_NULL, &size_shr_pac, &res, &addr_map_pac);
        assert(res == MPI_SUCCESS);

        res = MPI_Win_shared_query(win_shr_ref, MPI_PROC_NULL, &size_shr_ref, &res, &addr_map_ref);
        assert(res == MPI_SUCCESS);

        uint8_t *c = addr_map_ref;
        uint8_t *d = addr_map_pac;

        while(rank_shr == 0) {
            res = MPI_File_read(fh_pac_file, d, size_shr_pac, MPI_UINT8_T, &status);
            assert(res == MPI_SUCCESS);
            res = MPI_Get_count(&status, MPI_UINT8_T, &count4);
            assert(res == MPI_SUCCESS);
            if (count4 == 0) break;
            d += count4;
            size_read += count4;
            size_tot += count4;
        }

        size_read = 0;

        while(rank_shr == 0) {
            res = MPI_File_read(fh_ref_file, c, size_shr_ref, MPI_UINT8_T, &status);
            assert(res == MPI_SUCCESS);
            res = MPI_Get_count(&status, MPI_UINT8_T, &count4);
            assert(res == MPI_SUCCESS);
            if (count4 == 0) break;
            c += count4;
            size_read += count4;
            size_tot += count4;
        }

        res = MPI_Win_fence(0, win_shr_ref);
        assert(res == MPI_SUCCESS);

        res = MPI_Win_fence(0, win_shr_pac);
        assert(res == MPI_SUCCESS);

        res = MPI_File_close(&fh_pac_file);
        assert(res == MPI_SUCCESS);

        res = MPI_File_close(&fh_ref_file);
        assert(res == MPI_SUCCESS);

        aux.ref_string = addr_map_ref;
        aux.fmi->load_index_mpi(addr_map_pac);


        //aux.fmi->load_index();    
        tprof[FMI][0] += __rdtsc() - tim;


		///Create SAM header
		
        if (rank_num == 0) {
            if (write_format == 2)
                create_sam_header(file_out_ext, (aux.fmi->idx), &count, hdr_line, rg_line, pg_line, rank_num);
            else
                create_bam_header(file_out_ext, (aux.fmi->idx), &count, hdr_line, rg_line, pg_line, rank_num, compression_level);
        }
        
        if (hdr_line) free(hdr_line);
        if (rg_line) free(rg_line);
        if (pg_line) free(pg_line);

		///This is a testline to stop the program wherever I want
		if(0){ MPI_Barrier(MPI_COMM_WORLD); MPI_Finalize(); return 0;}

		res = MPI_File_open(MPI_COMM_WORLD, file_out_ext, MPI_MODE_WRONLY|MPI_MODE_APPEND, MPI_INFO_NULL, &fh_out);
		assert(res == MPI_SUCCESS);

		if (file_r1 != NULL) {
			//fprintf(stderr, "rank %d ::: open file %s \n",rank_num, file_r1);
			res = MPI_File_open(MPI_COMM_WORLD, file_r1, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_r1);
			assert(res == MPI_SUCCESS);
		}
		if (file_r2 != NULL) {
			//fprintf(stderr, "rank %d ::: open file %s \n",rank_num, file_r2);
			res = MPI_File_open(MPI_COMM_WORLD, file_r2, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_r2);
			assert(res == MPI_SUCCESS);
		}

		fprintf(stderr, "* Done reading reference genome..\n");

        if (ignore_alt)
        for (i = 0; i < fmi->idx->bns->n_seqs; ++i)
            fmi->idx->bns->anns[i].is_alt = 0;

		/*
		 *	In this part we initiate the worker threads  
		 *
		 */
        aux.ks2                       = 0;
        aux.ks                        = 0;
        aux.task_size = opt->chunk_size * opt->n_threads;
		fprintf(stderr, "\naux.task_size : %d\n", aux.task_size);
		worker_t	 w2;
		mem_opt_t	*opt			  = aux.opt;
		nthreads = opt->n_threads; // global variable for profiling!
        w2.nthreads = opt->n_threads;
        
        readLen = blen;
        //here we initialize pthreads 
		fprintf(stderr, "\nThreads used (compute): %d\n", nthreads);
        aux.actual_chunk_size = all_chunk_size[0];
        fprintf(stderr, "\naux2->actual_chunk_size: %d\n", aux.actual_chunk_size);
		nreads = aux.actual_chunk_size/ (readLen) + 10;
		fprintf(stderr, "Info: projected #read in a task: %ld\n", (long)nreads);
	
		/* All memory allocation */
		memoryAlloc(&aux, w2, nreads , nthreads);
        w2.ref_string = aux.ref_string;
        w2.fmi = aux.fmi;
        w2.nreads  = nreads;
		ktp_data_t *ret = (ktp_data_t *) calloc(1, sizeof(ktp_data_t));


		buffer_r1 = buffer_r2 = NULL; seqs = NULL;
		before_local_mapping = MPI_Wtime();

		// here we loop until there's nothing to read
		//we loop the chunck_count
		//size_t u1 = 0; 
		//for (u1 = 0; u1 < chunk_count; u1++){
        size_t  u1 = rank_num;
        while ( u1 < total_chunks){


			offset_chunk 		= all_begin_offset_chunk[u1];
			size_chunk   		= all_chunk_size[u1];
			offset_chunk_2 		= all_begin_offset_chunk_2[u1];
			size_chunk_2 		= all_chunk_size_2[u1];

			assert(size_chunk 	!= 0);
			assert(size_chunk_2 != 0);
		
			bef = MPI_Wtime();
			///allocate the buffer sizes to store an entire chunk in it
			buffer_r1 = malloc(size_chunk + 1);
			assert(buffer_r1 != NULL);
			buffer_r1[size_chunk] = 0;
		
			buffer_r2 = malloc(size_chunk_2 + 1);
			assert(buffer_r2 != NULL);
			buffer_r2[size_chunk_2]=0;

			// FOR TEST //
			//MPI_Type_contiguous(size_chunk, MPI_CHAR, &arraytype_r1);
			//MPI_Type_commit(&arraytype_r1);
			//MPI_File_set_view(fh_r2, (MPI_Offset)offset_chunk, MPI_CHAR, MPI_CHAR, "native", finfo ) ; 
			//res = MPI_File_read(fh_r2, buffer_r1, 1, arraytype_r1, &status);

			// FOR TEST //
			//MPI_Type_contiguous(size_chunk, MPI_CHAR, &arraytype_r2);
			//MPI_Type_commit(&arraytype_r2);
			//MPI_File_set_view(fh_r1, (MPI_Offset)offset_chunk, MPI_CHAR, MPI_CHAR, "native", finfo ) ; 
			//res = MPI_File_read(fh_r1, buffer_r2, 1, arraytype_r2, &status);		

		
			///Read the files and fill the buffers at the right offset and for the right size
			/*res = MPI_File_read_at(fh_r1, (MPI_Offset)offset_chunk, buffer_r1, size_chunk, MPI_CHAR, &status);
			assert(res == MPI_SUCCESS);
			res = MPI_Get_count(&status, MPI_CHAR, &count);
			assert(res == MPI_SUCCESS);
			assert(count == (int)size_chunk && *buffer_r1 == '@');

			res = MPI_File_read_at(fh_r2, (MPI_Offset)offset_chunk_2, buffer_r2, size_chunk_2, MPI_CHAR, &status);
			assert(res == MPI_SUCCESS);
			res = MPI_Get_count(&status, MPI_CHAR, &count);
			assert(res == MPI_SUCCESS);
			assert(count == (int) size_chunk_2 && *buffer_r2 == '@');
            */
            
            struct struct_pread_fastq *td_pread1;
            td_pread1 = malloc (NUM_THREADS * sizeof(struct struct_pread_fastq));
            bef = MPI_Wtime();
            pthread_attr_t attr4;
            pthread_attr_init(&attr4);
            pthread_attr_setstacksize(&attr4, BIG_STACK);
            pthread_attr_setdetachstate(&attr4, 0);

            for( n = 0; n < NUM_THREADS; n++ ){
                td_pread1[n].total_thread = NUM_THREADS;
                td_pread1[n].thread_id = n;
                td_pread1[n].job_rank = rank_num;
                td_pread1[n].offset= offset_chunk;
                td_pread1[n].size = size_chunk;
                td_pread1[n].buffer = buffer_r1;
                td_pread1[n].fd  = fh_r1;
                int ret_code = pthread_create(&threads[n], &attr4, pread_fastq_chunck, (void *)(&td_pread1[n]));
                if (ret_code) {
                    fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", ret_code);
                }
           }
           for(n=0; n<NUM_THREADS; n++)
                pthread_join(threads[n], (void *)(&td_pread1[n]));

            pthread_attr_destroy(&attr4);
            free(td_pread1);

            assert(strlen(buffer_r1) == size_chunk);
            assert(*buffer_r1 == '@');

            struct struct_pread_fastq *td_pread2;
            td_pread2 = malloc (NUM_THREADS * sizeof(struct struct_pread_fastq));
            bef = MPI_Wtime();
            pthread_attr_t attr5;
            pthread_attr_init(&attr5);
            pthread_attr_setstacksize(&attr5, BIG_STACK);
            pthread_attr_setdetachstate(&attr5, 0);

            for( n = 0; n < NUM_THREADS; n++ ){
                td_pread2[n].total_thread = NUM_THREADS;
                td_pread2[n].thread_id = n;
                td_pread2[n].job_rank = rank_num;
                td_pread2[n].offset= offset_chunk_2;
                td_pread2[n].size = size_chunk_2;
                td_pread2[n].buffer = buffer_r2;
                td_pread2[n].fd  = fh_r2;

                int ret_code = pthread_create(&threads[n], &attr5, pread_fastq_chunck, (void *)(&td_pread2[n]));
                if (ret_code) {
                    fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", ret_code);
                }
           }
           for(n=0; n<NUM_THREADS; n++)
                pthread_join(threads[n], (void *)(&td_pread2[n]));

            pthread_attr_destroy(&attr5);
            free(td_pread2);

            assert(strlen(buffer_r2) == size_chunk_2);
            assert(*buffer_r2 == '@');

            if (u1 < total_chunks){
                MPI_File_seek(fh_r1, (MPI_Offset)all_begin_offset_chunk[u1 + proc_num], MPI_SEEK_CUR ) ;
                MPI_File_seek(fh_r2, (MPI_Offset)all_begin_offset_chunk[u1 + proc_num], MPI_SEEK_CUR ) ;
            }

			///some stats
			aft = MPI_Wtime();
			fprintf(stderr, "%s: num_rank = %d :: read sequences (%.02f)\n", __func__, rank_num, aft - bef);
			total_time_reading_seq += (aft - bef);

			///operations on the number of reads
			bef = MPI_Wtime();
			reads_r1 = all_reads_in_chunk[u1];
			reads_r2 = all_reads_in_chunk_2[u1];
			reads = reads_r1 + reads_r2;
			bases = 0;
			total_local_reads_aligned += reads/2;
			assert(reads <= INT_MAX);

			/* Parse sequences ... */
			seqs = malloc(reads * sizeof(*seqs));
			assert(seqs != NULL);

			if (file_r1 != NULL) {
				p1 = q1 = buffer_r1; e1 = buffer_r1 + size_chunk; line_number = 0;
							
				while (q1 < e1) {
					if (*q1 != '\n') { q1++; continue; }
					/* We have a full line ... process it */
					*q1 = '\0'; 
					//n will take 4 times in row the same value from 0 with a step of 2
					//it is because the 2 operations are using int types, so it is cast as int twice (within and then outside the brackets)
					n = files * (line_number / 4);
				
					switch (line_number % 4) {
					case 0: /* Line1: Name and Comment */
						assert(*p1 == '@');
						seqs[n].name   = p1 + 1;
						while (*p1 && !isspace((unsigned char)*p1)) {p1++;}
						if (*(p1-2) == '/' && isdigit((unsigned char)*(p1-1))) {*(p1-2) = '\0';}
						if (*p1) {*p1++ = '\0';}
						seqs[n].comment = (copy_comment != 0) ? p1 : NULL;
						seqs[n].sam = NULL;
						break;
					case 1: /* Line2: Sequence */
						seqs[n].seq = p1;
						seqs[n].l_seq = q1 - p1;
						bases += seqs[n].l_seq;
						break;
					case 2: /* Line3: Ignored */
						assert(*p1 == '+');
						break;
					case 3: /* Line4: Quality */
						seqs[n].qual = p1;
						break; }
					p1 = ++q1; 
					line_number++; 
				}
			}
			assert( (line_number/4) == reads_r1);
		
			if (file_r2 != NULL) {
				p2 = q2 = buffer_r2; e2 = buffer_r2 + size_chunk_2; line_number = 0;
							
				while (q2 < e2) {
					if (*q2 != '\n') { q2++; continue; }
					/* We have a full line ... process it */
					*q2 = '\0'; 
					//n will take 4 times in row the same value from 0 with a step of 2
					//it is because the 2 operations are using int types, so it is cast as int twice (within and then outside the brackets)
					n = files * (line_number / 4);
				
					switch (line_number % 4) {
					case 0: /* Line1: Name and Comment */
						assert(*p2 == '@');
						seqs[n+1].name   = p2 + 1;
						while (*p2 && !isspace((unsigned char)*p2)) {p2++;}
						if (*(p2-2) == '/' && isdigit((unsigned char)*(p2-1))) {*(p2-2) = '\0';}
						if (*p2) {*p2++ = '\0';}
						seqs[n+1].comment = (copy_comment != 0) ? p2 : NULL;
						seqs[n+1].sam = NULL;
						break;
					case 1: /* Line2: Sequence */
						seqs[n+1].seq = p2;
						seqs[n+1].l_seq = q2 - p2;
						bases += seqs[n+1].l_seq;
						break;
					case 2: /* Line3: Ignored */
						assert(*p2 == '+');
						break;
					case 3: /* Line4: Quality */
						seqs[n+1].qual = p2;
						break; }
					p2 = ++q2; 
					line_number++; 
				}
			}

			assert( (line_number/4) == reads_r2);
			aft = MPI_Wtime();

			fprintf(stderr, "rank %d :: %s: parsed %zu read of sequences (%ld bp) in (%.02f)\n", rank_num, __func__, reads, (long)bases, aft - bef);
			total_time_parsing += (aft - bef);
			//fprintf(stderr, "rank %d ::: [M::%s] read %zu sequences (%ld bp)...\n",rank_num , __func__, reads, (long)bases);
			/* Datas computation ... */
			bef = MPI_Wtime();
			fprintf(stderr, "rank %d ::::  Call memseqs with count sequences %zu \n", rank_num, reads);
			//mem_process_seqs(opt, indix.bwt, indix.bns, indix.pac, 0, (int)reads, seqs, pes0);
			ret->seqs   = seqs;
            ret->n_seqs     = reads;
            int64_t size = 0;
            for (int i = 0; i < ret->n_seqs; ++i) size += ret->seqs[i].l_seq;
            fprintf(stderr, "\t[0000][ M::%s] read %d sequences (%ld bp)...\n",
                    __func__, ret->n_seqs, (long)size);
            fprintf(stderr, "[RANK %d] 1: Calling process()\n", rank_num);
            if (w2.regs) free(w2.regs);
            if (w2.chain_ar) free(w2.chain_ar);
            if (w2.seedBuf) free(w2.seedBuf);

            nreads = ret->n_seqs;
            w2.regs = (mem_alnreg_v *) calloc(nreads, sizeof(mem_alnreg_v));
            w2.chain_ar = (mem_chain_v*) malloc (nreads * sizeof(mem_chain_v));
            w2.seedBuf = (mem_seed_t *) calloc(sizeof(mem_seed_t),nreads * AVG_SEEDS_PER_READ);
            w2.nreads=ret->n_seqs;
            mem_process_seqs(opt,
                             aux.n_processed,
                             ret->n_seqs,
                             ret->seqs,
                             aux.pes0,
                             w2);
            aux.n_processed += ret->n_seqs;
            aft = MPI_Wtime();
            fprintf(stderr, "rank %d :: %s: computed mappings (%.02f)\n", rank_num, __func__, aft - bef);

            total_time_mapping += (aft - bef);


            aft = MPI_Wtime();
			fprintf(stderr, "rank %d :: %s: computed mappings (%.02f)\n", rank_num, __func__, aft - bef);

			total_time_mapping += (aft - bef);

            if (dofixmate){
                bef = MPI_Wtime();
                int ret_code = 0;
                struct thread_data *td;
                td = malloc (NUM_THREADS * sizeof(struct thread_data));
                bef = MPI_Wtime();
                pthread_attr_init(&attr);
                pthread_attr_setstacksize(&attr, SMALL_STACK);
                pthread_attr_setdetachstate(&attr, 0);

                for( n = 0; n < NUM_THREADS; n++ ){
                    td[n].total_thread = NUM_THREADS;
                    td[n].thread_id = n;
                    assert(seqs);
                    td[n].seqs_thr = seqs;
                    td[n].job_rank = rank_num;
                    td[n].total_reads = reads;
                    td[n].start_index_seqs = 0;
                    td[n].final_index_seqs = reads-1;
                    td[n].indix_thr = &indix;
                    td[n].total_lines = 0;
                    ret_code = pthread_create(&threads[n], &attr, call_fixmate, (void *)(&td[n]));
                    if (ret_code) {
                        fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", ret_code);
                    }
                }
                for(n=0; n<NUM_THREADS; n++)
                    pthread_join(threads[n], (void *)(&td[n]));
                pthread_attr_destroy(&attr);
                free(td);
                aft = MPI_Wtime();
                fprintf(stderr, "rank: %d :: %s: time spend in fixmate (%.02f) \n", rank_num, __func__, aft - bef);
            }
		
            if ( write_format == 2){

                bef = MPI_Wtime();
                pthread_attr_t attr;
                pthread_attr_init(&attr);
                pthread_attr_setstacksize(&attr, BIG_STACK);
                pthread_attr_setdetachstate(&attr, 0);

                struct struct_data_thread *td;
                td = calloc (opt->n_threads, sizeof(struct_data_thread));
                
                int rest = reads%NUM_THREADS;
                int quot = reads/NUM_THREADS;
                pthread_t threads[NUM_THREADS];
                
                int ret_code = 0;
                for ( n = 0; n < NUM_THREADS; n++){
                    td[n].seqs_thr = seqs;
                    td[n].begin_index = n * quot;
                    td[n].end_index = ((n+1) * quot);
                    if ( n == (NUM_THREADS - 1)) td[n].end_index = reads;
                    td[n].file_desc = fh_out;
                    ret_code = pthread_create(&threads[n], &attr, write_sam_mt, (void *)(&td[n]));
                }

                for (n = 0; n < NUM_THREADS; n++)
                    pthread_join(threads[n], (void *)(&td[n]));

                pthread_attr_destroy(&attr);
                free(td);
                aft = MPI_Wtime();
                total_time_writing += (aft - bef);
                free(buffer_r1);
                free(buffer_r2);
 
            }

            if ( write_format == 0 ) {

                int ret_code = 0;
                struct thread_data_compress *tdc;
                tdc = malloc (NUM_THREADS * sizeof(struct thread_data_compress));
                bef = MPI_Wtime();
                pthread_attr_init(&attr);
                pthread_attr_setstacksize(&attr, SMALL_STACK);
                pthread_attr_setdetachstate(&attr, 0);

                for( n = 0; n < NUM_THREADS; n++ ){

                    tdc[n].total_thread = NUM_THREADS;
                    tdc[n].thread_id = n;
                    tdc[n].seqs_thr = seqs;
                    tdc[n].job_rank = rank_num;
                    tdc[n].total_reads = reads;
                    tdc[n].thr_comp_sz = 0;
                    tdc[n].comp_level = compression_level;
                    tdc[n].compressed_buffer_thread  = 0;
                    tdc[n].fh_out = fh_out;
                    ret_code = pthread_create(&threads[n], &attr, compress_and_write_bgzf_thread, (void *)(&tdc[n]));
                    if (ret_code) {
                        fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", ret_code);
                    }
                }
                for(n=0; n<NUM_THREADS; n++)
                   pthread_join(threads[n], (void *)(&tdc[n]));
                pthread_attr_destroy(&attr);
                free(tdc);
                for (n = 0; n < reads; n++) free(seqs[n].sam);
                free(seqs);
		aft = MPI_Wtime();
                total_time_writing += (aft - bef);
                free(buffer_r1);
                free(buffer_r2);
            }
            
            if ( write_format == 1 ) {

                int ret_code = 0;
                struct thread_data_compress *tdc;
                tdc = malloc (NUM_THREADS * sizeof(struct thread_data_compress));
                bef = MPI_Wtime();
                pthread_attr_init(&attr);
                pthread_attr_setstacksize(&attr, SMALL_STACK);
                pthread_attr_setdetachstate(&attr, 0);

                for( n = 0; n < NUM_THREADS; n++ ){

                    tdc[n].total_thread = NUM_THREADS;
                    tdc[n].thread_id = n;
                    tdc[n].seqs_thr = seqs;
                    tdc[n].job_rank = rank_num;
                    tdc[n].total_reads = reads;
                    tdc[n].thr_comp_sz = 0;
                    tdc[n].comp_level = compression_level;
                    tdc[n].compressed_buffer_thread  = 0;
                    tdc[n].fh_out = fh_out;
                    ret_code = pthread_create(&threads[n], &attr, compress_and_write_bam_thread, (void *)(&tdc[n]));
                    if (ret_code) {
                        fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", ret_code);
                    }
                }
                for(n=0; n<NUM_THREADS; n++)
                   pthread_join(threads[n], (void *)(&tdc[n]));
                pthread_attr_destroy(&attr);
                free(tdc);
                for (n = 0; n < reads; n++) free(seqs[n].sam);
                free(seqs);
		aft = MPI_Wtime();
                total_time_writing += (aft - bef);
                free(buffer_r1);
                free(buffer_r2);
            }

            
            fprintf(stderr, "%s: wrote results (%.02f)\n", __func__, aft - bef);
			fprintf(stderr, "rank: %d :: finish for chunck %zu \n", rank_num, u1);
            u1 += proc_num;

	} //end for loop on chunks


        MPI_Barrier(MPI_COMM_WORLD);
        if ( (write_format == 1) && (rank_num == 0)) {

                res = MPI_File_write_shared(fh_out, magic, 28, MPI_BYTE, &status);
                assert(res == MPI_SUCCESS);
                res = MPI_Get_count(&status, MPI_BYTE, &count);
                assert(res == MPI_SUCCESS);
                assert(count == 28);
        }

        free(w2.chain_ar);
        free(w2.regs);
        free(w2.seedBuf);

        for(int l=0; l<nthreads; l++) {
            _mm_free(w2.mmc.seqBufLeftRef[l*CACHE_LINE]);
            _mm_free(w2.mmc.seqBufRightRef[l*CACHE_LINE]);
            _mm_free(w2.mmc.seqBufLeftQer[l*CACHE_LINE]);
            _mm_free(w2.mmc.seqBufRightQer[l*CACHE_LINE]);
        }

        for(int l=0; l<nthreads; l++) {
            free(w2.mmc.seqPairArrayAux[l]);
            free(w2.mmc.seqPairArrayLeft128[l]);
            free(w2.mmc.seqPairArrayRight128[l]);
        }

        for(int l=0; l<nthreads; l++) {
            _mm_free(w2.mmc.matchArray[l]);
            free(w2.mmc.min_intv_ar[l]);
            free(w2.mmc.query_pos_ar[l]);
            free(w2.mmc.enc_qdb[l]);
            free(w2.mmc.rid[l]);
            _mm_free(w2.mmc.lim[l]);
        }

        MPI_Win_free(&win_shr_ref);
        MPI_Win_free(&win_shr_pac);

         /* shall we delete?
          * because error in read_index_ele.cpp:52
          */
        //delete(aux.fmi);
        //


	}// end else case files are trimmed

if (file_r1 != NULL && file_r2 == NULL){

		/*
		 *
		 * We are in the case reads are single
		 * We do both case trimmed or not 
		 *
		 * We use the same algo as in paired trimmed
		 */
		
        
		aft = 0; aft++;
		bef = 0; bef++;

        size_t slen;
        int32_t blen;

        MPI_File mpi_fd_in1;
        bef = MPI_Wtime();

        res = MPI_File_open(MPI_COMM_WORLD, file_r1,  MPI_MODE_RDONLY , MPI_INFO_NULL, &mpi_fd_in1);
        assert(res == MPI_SUCCESS);

		//global offset contains the starting offset in the fastq for each process
		//TODO?: use only one vector for all the global offsets iot reduce communications
		size_t *goff 	= malloc((proc_num * NUM_THREADS + 1) * sizeof(size_t));
        size_t *goff_inter = calloc( (proc_num * NUM_THREADS + 1) , sizeof(size_t));
		
		//this function is used to fill the goff vectors
		//find_process_starting_offset(goff, stat_r1.st_size, file_r1, proc_num, rank_num);
        find_process_starting_offset_mt(goff, stat_r1.st_size, file_r1, proc_num, rank_num, NUM_THREADS);		


		//now we exchange the goff buffer between all proc
		int i12=0;
        for ( i12 = 0; i12 < proc_num * NUM_THREADS + 1; i12++ )
            goff_inter[i12] = goff[i12];
				
        res = MPI_Allgather(&goff_inter[rank_num*NUM_THREADS], NUM_THREADS, MPI_LONG_LONG_INT, goff , NUM_THREADS, MPI_LONG_LONG_INT, MPI_COMM_WORLD);
        assert(res == MPI_SUCCESS);

        free(goff_inter);
						
		///resources needed to find the offsets and size of each read.	
		size_t grand_total_num_reads 	= 0; 
		size_t local_num_reads 			= 0;
		size_t total_num_reads 			= 0;
		
		//Here I decided to keep separated vectors because otherwise with all the reallocs they would be too big and take too much contiguous memory space
		int *local_read_size    		= calloc(1, sizeof(int));
		size_t *local_read_bytes    	= calloc(1, sizeof(size_t));
		size_t *local_read_offsets 		= calloc(1, sizeof(size_t));
		
		assert( local_read_offsets 		!= NULL);
		assert( local_read_size 		!= NULL);
		assert( local_read_bytes 		!= NULL);
		
		///find offsets and sizes for the first file
		bef = MPI_Wtime();
        pthread_attr_t attr;
        pthread_attr_init(&attr);
        pthread_attr_setstacksize(&attr, SMALL_STACK);
        pthread_attr_setdetachstate(&attr, 0);
        pthread_t threads_1[NUM_THREADS];

        struct struct_data_thread_1 *td_1;
        size_t *local_num_reads_t            = calloc(NUM_THREADS, sizeof(size_t));
        size_t *total_num_reads_t            = calloc(NUM_THREADS, sizeof(size_t));
        size_t **local_read_offsets_t        = calloc(NUM_THREADS, sizeof(size_t*));
        size_t **local_read_bytes_t          = calloc(NUM_THREADS, sizeof(size_t*));
        int **local_read_size_t              = calloc(NUM_THREADS, sizeof(int*));

        td_1 = calloc(NUM_THREADS, sizeof(struct struct_data_thread_1));
        int ret_code_1 = 0;
        int goff_idx = 0;
        for ( n = 0; n < NUM_THREADS; n++){
            goff_idx = (rank_num * NUM_THREADS) + n;
            td_1[n].offset_in_file_mt         = goff[goff_idx];
            td_1[n].size2read_mt              = goff[goff_idx + 1] - goff[goff_idx];
            td_1[n].file_r1_mt                = file_r1;
            td_1[n].local_num_reads_mt        = &local_num_reads_t[n];
            td_1[n].total_num_reads_mt        = &total_num_reads_t[n];
            td_1[n].local_read_offsets_mt     = &local_read_offsets_t[n];
            td_1[n].local_read_size_mt        = &local_read_size_t[n];
            td_1[n].local_read_bytes_mt       = &local_read_bytes_t[n];
            td_1[n].proc_num_mt               = proc_num;
            td_1[n].rank_num_mt               = rank_num;
            td_1[n].thread_num_mt             = n;
            td_1[n].previous_read_num         = 0;
            ret_code_1 = pthread_create(&threads_1[n], &attr, find_reads_size_and_offsets_mt, (void *)(&td_1[n]));

        }

        total_num_reads = 0;
        for (n = 0; n < NUM_THREADS; n++){
            pthread_join(threads_1[n], (void *)(&td_1[n]));
            total_num_reads += *(td_1[n].total_num_reads_mt);
        }

        local_read_offsets  = calloc(total_num_reads, sizeof(size_t));
        local_read_size     = calloc(total_num_reads, sizeof(int));
        local_read_bytes    = calloc(total_num_reads, sizeof(size_t));

        assert(local_read_offsets);
        assert(local_read_size);
        assert(local_read_bytes);

        size_t tmp_var = 0;
        for (n = 0; n < NUM_THREADS; n++){
            td_1[n].local_read_offsets     = local_read_offsets;
            td_1[n].local_read_size        = local_read_size;
            td_1[n].local_read_bytes       = local_read_bytes;
            td_1[n].previous_read_num      = tmp_var;
            tmp_var                        += *(td_1[n].total_num_reads_mt);
            ret_code_1 = pthread_create(&threads_1[n], &attr, copy_local_read_info_mt, (void *)(&td_1[n]));
        }

        for (n = 0; n < NUM_THREADS; n++){
            pthread_join(threads_1[n], (void *)(&td_1[n]));
            free(local_read_offsets_t[n]);
            free(local_read_bytes_t[n]);
            free(local_read_size_t[n]);
        }

        free(local_num_reads_t);
        free(total_num_reads_t);
        free(local_read_offsets_t);
        free(local_read_bytes_t);
        free(local_read_size_t);

        pthread_attr_destroy(&attr);
        free(td_1);

        res = MPI_Reduce(&total_num_reads, &grand_total_num_reads, 1, MPI_LONG_LONG_INT, MPI_SUM, 0,MPI_COMM_WORLD);
        assert ( res == MPI_SUCCESS);

        res = MPI_Bcast(&grand_total_num_reads, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        assert ( res == MPI_SUCCESS);

        if ( rank_num == 0 )
            fprintf(stderr, "rank %d ::: total num reads in fastq = %zu \n", rank_num, grand_total_num_reads);

        MPI_Barrier(MPI_COMM_WORLD);

		aft = MPI_Wtime();
		fprintf(stderr, "%s: rank %d num reads parsed: %zu ::: time spend reading and parsing entire fastq file = (%.02f) \n", __func__, rank_num, total_num_reads, aft - bef);

		if (goff) free(goff);

		MPI_Barrier(MPI_COMM_WORLD);
			
		///Now find the required information to create the chunks, such as total number of bases
		local_num_reads 	= total_num_reads;
        
        //now we estimate the number of chunk per rank
        size_t min_num_read = local_num_reads;
        size_t bases_tmp    = 0;
        size_t chunck_num   = 0;
        
        for(i=0; i < min_num_read; i++)   {

            bases_tmp  += (local_read_size[i]);

            if ( bases_tmp > ( opt->chunk_size * opt->n_threads)){

                bases_tmp = 0;
                chunck_num++;
            }
        }

        chunck_num  += 2;
        	
		fprintf(stderr,"Rank %d :: chunk_num = %zu\n", rank_num, chunck_num);
		
		///assert that the sizes and offsets found earlier are real.
		//changed the comparisons here for the size. they were compared to blen which has become impossible
		size_t h=0;
		
		for (h=0; h < total_num_reads; h++)
		{
			assert(local_read_size[h] > 0 );
			assert(local_read_offsets[h] >= 0 );
		}

        blen = local_read_size[0];
		
		///Now that we know how many chunks we need, we create the vectors associated
		///and we update their informations according to the file they are treating

		// we allocate vector for chunks offset
		size_t *begin_offset_chunk 	= calloc(chunck_num, sizeof(size_t));
		size_t *chunk_size 			= calloc(chunck_num, sizeof(size_t));
		size_t *reads_in_chunk 		= calloc(chunck_num, sizeof(size_t));

		assert( begin_offset_chunk != NULL );
		assert( chunk_size != NULL );
		assert( reads_in_chunk != NULL );
	
		size_t chunk_count = 0;

		maxsiz = ( opt->chunk_size * opt->n_threads); 
		MPI_Barrier(MPI_COMM_WORLD);
		fprintf(stderr,"rank %d ::: Call find_chunks_info \n", rank_num);
		// the detail of he paramters is at the function definition
		// fprintf(stderr, "rank %d ::: begin_offset_chunk = %zu \n", rank_num, begin_offset_chunk[0]);
		// fprintf(stderr, "rank %d ::: begin_offset_chunk_2 = %zu \n", rank_num, begin_offset_chunk_2[0]);
		// fprintf(stderr, "rank %d ::: chunk_size = %zu \n", rank_num, chunk_size);
		// fprintf(stderr, "rank %d ::: chunk_size_2 = %zu \n", rank_num, chunk_size_2);

		bef = MPI_Wtime();
		find_chunks_info(begin_offset_chunk,
				chunk_size,
				reads_in_chunk,
				local_read_size,
				local_read_bytes,
				local_read_offsets,
				rank_num,
				proc_num,
				local_num_reads,
				grand_total_num_reads,
				maxsiz,
				&chunk_count);

		aft = MPI_Wtime();
		fprintf(stderr, "%s ::: rank %d ::: evaluating offsets chuncks and sizes (%.02f) found %zu chuncks \n", __func__, rank_num, aft - bef, chunk_count);
	
		MPI_Barrier(MPI_COMM_WORLD);

		if (local_read_offsets) 	free(local_read_offsets);
		if (local_read_size) 		free(local_read_size);
		
		//verify the number of reads

		size_t num_reads_1 = 0;
		size_t total_num_reads_v1 = 0;
		
		for(h=0; h < chunk_count; h++){
			num_reads_1 += reads_in_chunk[h];
		}

		fprintf(stderr, "rank %d ::: after find_chunks_info ::: num_reads_1 = %zu \n", rank_num, num_reads_1);
		
		res = MPI_Reduce(&num_reads_1, &total_num_reads_v1, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		assert ( res == MPI_SUCCESS);

		res = MPI_Bcast(&total_num_reads_v1, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
		assert ( res == MPI_SUCCESS);

		if ( rank_num == 0 || rank_num ==1 )
			fprintf(stderr, "rank %d ::: after find_chunks_info ::: total_num_reads forward = %zu \n", rank_num, total_num_reads_v1);

		MPI_Barrier(MPI_COMM_WORLD);

		//display test to check that each vector was well filled with the right information
		//fprintf(stderr,"rank %d ::: begin offset %zu, %zu\n", rank_num, begin_offset_chunk[0], begin_offset_chunk_2[0]);
		//fprintf(stderr,"rank %d ::: chunk size %zu, %zu\n", rank_num, chunk_size[0], chunk_size_2[0]);
		//fprintf(stderr,"rank %d ::: reads in chunk %zu, %zu\n", rank_num, reads_in_chunk[0], reads_in_chunk_2[0]);

		/// Map reference genome indexes in shared memory (by host)
		/*bef = MPI_Wtime();
		map_indexes(file_map, &count, &indix, &ignore_alt, &win_shr);
		aft = MPI_Wtime();
		fprintf(stderr, "%s ::: rank %d ::: mapped indexes (%.02f)\n", __func__, rank_num, aft - bef);
		*/

         size_t total_chunks = 0;

        res = MPI_Reduce(&chunk_count, &total_chunks, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        assert ( res == MPI_SUCCESS);
        res = MPI_Bcast(&total_chunks, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        assert ( res == MPI_SUCCESS);

        //fprintf(stderr, "%s: rank %d total chunks = %zu \n", __func__, rank_num, total_chunks);
        //fprintf(stderr, "%s: rank %d local chunks = %zu \n", __func__, rank_num, chunk_count);

        size_t chunks_per_rank[proc_num];
        res = MPI_Allgather(&chunk_count, 1, MPI_LONG_LONG_INT, chunks_per_rank, 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);

        //for (i = 0; i < proc_num; i++) fprintf(stderr, "%s: rank %d chck per rank[%d] = %zu \n", __func__, rank_num, i, chunks_per_rank[i]);


        all_begin_offset_chunk = malloc(total_chunks * sizeof(size_t));
        all_begin_offset_chunk[0]=0;
        all_chunk_size = malloc(total_chunks * sizeof(size_t));
        all_chunk_size[0] = 0;
        all_reads_in_chunk =  malloc(total_chunks * sizeof(size_t));
        all_reads_in_chunk[0] = 0;

        int displ_chunk[proc_num];
        displ_chunk[0] = 0;

        for (i = 1; i < proc_num; i++) displ_chunk[i] = (displ_chunk[i-1] + chunks_per_rank[i-1]);

        int indx=displ_chunk[rank_num];
        for (i = 0; i <  chunks_per_rank[rank_num]; i++) {
            all_chunk_size[indx+i]=chunk_size[i];
            all_begin_offset_chunk[indx+i]=begin_offset_chunk[i];
            all_reads_in_chunk[indx+i]=reads_in_chunk[i];
        }

         if (rank_num > 0){
            res=MPI_Send(chunk_size, chunks_per_rank[rank_num], MPI_LONG_LONG_INT, 0, rank_num, MPI_COMM_WORLD);
            assert(res == MPI_SUCCESS);
            }
        else{
            for (i = 1; i < proc_num; i++){
                res=MPI_Recv(&(all_chunk_size[displ_chunk[i]]), chunks_per_rank[i], MPI_LONG_LONG_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                assert(res == MPI_SUCCESS);
            }
        }

        if (rank_num > 0){
            res=MPI_Send(begin_offset_chunk, chunks_per_rank[rank_num], MPI_LONG_LONG_INT, 0, rank_num, MPI_COMM_WORLD);
            assert(res == MPI_SUCCESS);
            }
        else{
            for (i = 1; i < proc_num; i++){
                res=MPI_Recv(&(all_begin_offset_chunk[displ_chunk[i]]), chunks_per_rank[i], MPI_LONG_LONG_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                assert(res == MPI_SUCCESS);
            }
        }


        if (rank_num > 0){
            res=MPI_Send(reads_in_chunk, chunks_per_rank[rank_num], MPI_LONG_LONG_INT, 0, rank_num, MPI_COMM_WORLD);
            assert(res == MPI_SUCCESS);
            }
        else{
            for (i = 1; i < proc_num; i++){
                res=MPI_Recv(&(all_reads_in_chunk[displ_chunk[i]]), chunks_per_rank[i], MPI_LONG_LONG_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                assert(res == MPI_SUCCESS);
            }
        }


        res=MPI_Bcast(all_chunk_size, total_chunks, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        res=MPI_Bcast(all_reads_in_chunk, total_chunks, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        res=MPI_Bcast(all_begin_offset_chunk, total_chunks, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

        assert(res == MPI_SUCCESS);

        size_t total = 0;

        for ( i = 0; i < total_chunks; i++){
                assert(all_chunk_size[i] != 0);
                total += all_reads_in_chunk[i];    
        }

        assert (total_num_reads_v1 == total);

        free(chunk_size);
        free(begin_offset_chunk);
        free(reads_in_chunk);

        aft = MPI_Wtime();
        fprintf(stderr, "%s: rank %d time spend gathering chunks info = (%.02f) \n", __func__, rank_num, aft - bef);


		uint64_t tim = __rdtsc();
        fprintf(stderr, "* Ref file: %s\n", file_ref);
        aux.fmi = new FMI_search(file_ref);

        uint8_t *shared_ref;
        MPI_Comm comm_shr;
        int count4, rank_shr = 0;
        size_t size_read;

        res = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_shr);
        assert(res == MPI_SUCCESS);

        res = MPI_Comm_rank(comm_shr, &rank_shr);
        assert(res == MPI_SUCCESS);

        char binary_seq_file[200];
        sprintf(binary_seq_file, "%s.0123", file_ref);

        char pac_file[200];
        sprintf(pac_file, "%s.pac", file_ref);

        MPI_File fh_ref_file, fh_pac_file;
        MPI_Aint size_shr_ref, size_shr_pac;
        MPI_Offset size_map_ref,size_map_pac;
        MPI_Win win_shr_ref, win_shr_pac;

        uint8_t *addr_map_ref;
        uint8_t *addr_map_pac;


        MPI_Info win_info;
        MPI_Info_create(&win_info);
        MPI_Info_set(win_info, "alloc_shared_non_contig", "true");

        res = MPI_File_open(MPI_COMM_WORLD, binary_seq_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_ref_file);
        assert(res == MPI_SUCCESS);

        res = MPI_File_open(MPI_COMM_WORLD, pac_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_pac_file);
        assert(res == MPI_SUCCESS);

        res = MPI_File_get_size(fh_pac_file, &size_map_pac);
        assert(res == MPI_SUCCESS);

        res = MPI_File_get_size(fh_ref_file, &size_map_ref);
        assert(res == MPI_SUCCESS);

        size_shr_ref = (rank_shr == 0) ? size_map_ref : 0;
        res = MPI_Win_allocate_shared(size_shr_ref * sizeof(uint8_t), 8 , win_info, comm_shr, &addr_map_ref, &win_shr_ref);

        size_shr_pac = (rank_shr == 0) ? size_map_pac : 0;
        res = MPI_Win_allocate_shared(size_shr_pac * sizeof(uint8_t), 8 , win_info, comm_shr, &addr_map_pac, &win_shr_pac);

        res = MPI_Win_shared_query(win_shr_pac, MPI_PROC_NULL, &size_shr_pac, &res, &addr_map_pac);
        assert(res == MPI_SUCCESS);

        res = MPI_Win_shared_query(win_shr_ref, MPI_PROC_NULL, &size_shr_ref, &res, &addr_map_ref);
        assert(res == MPI_SUCCESS);

        uint8_t *c = addr_map_ref;
        uint8_t *d = addr_map_pac;

        while(rank_shr == 0) {
            res = MPI_File_read(fh_pac_file, d, size_shr_pac, MPI_UINT8_T, &status);
            assert(res == MPI_SUCCESS);
            res = MPI_Get_count(&status, MPI_UINT8_T, &count4);
            assert(res == MPI_SUCCESS);
            if (count4 == 0) break;
            d += count4;
            size_read += count4;
            size_tot += count4;
        }

        size_read = 0;

        while(rank_shr == 0) {
            res = MPI_File_read(fh_ref_file, c, size_shr_ref, MPI_UINT8_T, &status);
            assert(res == MPI_SUCCESS);
            res = MPI_Get_count(&status, MPI_UINT8_T, &count4);
            assert(res == MPI_SUCCESS);
            if (count4 == 0) break;
            c += count4;
            size_read += count4;
            size_tot += count4;
        }

        res = MPI_Win_fence(0, win_shr_ref);
        assert(res == MPI_SUCCESS);

        res = MPI_Win_fence(0, win_shr_pac);
        assert(res == MPI_SUCCESS);

        res = MPI_File_close(&fh_pac_file);
        assert(res == MPI_SUCCESS);

        res = MPI_File_close(&fh_ref_file);
        assert(res == MPI_SUCCESS);

        aux.ref_string = addr_map_ref;
        aux.fmi->load_index_mpi(addr_map_pac);

        tprof[FMI][0] += __rdtsc() - tim;

		///Create SAM header
		//TODO: Add line for BWA version
		//create_sam_header(file_out, aux, &count, hdr_line, rg_line, rank_num);

        if (rank_num == 0) {
            if (write_format == 2)
                create_sam_header(file_out_ext, (aux.fmi->idx), &count, hdr_line, rg_line, pg_line, rank_num);
            else
                create_bam_header(file_out_ext, (aux.fmi->idx), &count, hdr_line, rg_line, pg_line, rank_num, compression_level);
        }


		///This is a testline to stop the program wherever I want
		if(0){ MPI_Barrier(MPI_COMM_WORLD); MPI_Finalize(); return 0;}

		res = MPI_File_open(MPI_COMM_WORLD, file_out_ext, MPI_MODE_WRONLY|MPI_MODE_APPEND, MPI_INFO_NULL, &fh_out);
		assert(res == MPI_SUCCESS);

		if (file_r1 != NULL) {
			//fprintf(stderr, "rank %d ::: open file %s \n",rank_num, file_r1);
			res = MPI_File_open(MPI_COMM_WORLD, file_r1, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_r1);
			assert(res == MPI_SUCCESS);
		}
		
		/* Matrix for SWA */
        bwa_fill_scmat(opt->a, opt->b, opt->mat);

         // reading ref string from the file
        tim = __rdtsc();
	
        fprintf(stderr, "* Done reading reference genome !!\n\n");

        if (ignore_alt)
        for (i = 0; i < fmi->idx->bns->n_seqs; ++i)
            fmi->idx->bns->anns[i].is_alt = 0;

		/*
		 *	In this part we initiate the worker threads  
		 *
		 */
        aux.ks2                       = 0;
        aux.ks                        = 0;
        aux.task_size = opt->chunk_size * opt->n_threads;
		fprintf(stderr, "\naux.task_size : %d\n", aux.task_size);
		worker_t	 w2;
		mem_opt_t	*opt			  = aux.opt;
		nthreads = opt->n_threads; // global variable for profiling!
        w2.nthreads = opt->n_threads;
        //int32_t blen=151;
        readLen = blen; 
		//here we initialize pthreads 
		fprintf(stderr, "\nThreads used (compute): %d\n", nthreads);
        aux.actual_chunk_size = all_chunk_size[0];
        fprintf(stderr, "\naux2->actual_chunk_size: %d\n", aux.actual_chunk_size);
        
		nreads = aux.actual_chunk_size/ (readLen) + 10;
		fprintf(stderr, "Info: projected #read in a task: %ld\n", (long)nreads);
	
		/* All memory allocation */
		memoryAlloc(&aux, w2, nreads , nthreads);
        w2.ref_string = aux.ref_string;
        w2.fmi = aux.fmi;
        w2.nreads  = nreads;
		ktp_data_t *ret = (ktp_data_t *) calloc(1, sizeof(ktp_data_t));

		buffer_r1 = NULL; seqs = NULL;
		before_local_mapping = MPI_Wtime();

		// here we loop until there's nothing to read
		//we loop the chunck_count
		//size_t u1 = 0; 
		//for (u1 = 0; u1 < chunk_count; u1++){
        size_t u1 = rank_num;
        while ( u1 < total_chunks){

			offset_chunk 		= all_begin_offset_chunk[u1];
			size_chunk   		= all_chunk_size[u1];
			
			assert(size_chunk 	!= 0);
					
			bef = MPI_Wtime();
			///allocate the buffer sizes to store an entire chunk in it
			buffer_r1 = malloc(size_chunk + 1);
			assert(buffer_r1 != NULL);
			buffer_r1[size_chunk] = 0;
		
			// FOR TEST //
			//MPI_Type_contiguous(size_chunk, MPI_CHAR, &arraytype_r2);
			//MPI_Type_commit(&arraytype_r2);
			//MPI_File_set_view(fh_r1, (MPI_Offset)offset_chunk, MPI_CHAR, MPI_CHAR, "native", finfo ) ; 
			//res = MPI_File_read(fh_r1, buffer_r2, 1, arraytype_r2, &status);		

		
			///Read the files and fill the buffers at the right offset and for the right size
			//res = MPI_File_read_at(fh_r1, (MPI_Offset)offset_chunk, buffer_r1, size_chunk, MPI_CHAR, &status);
			//assert(res == MPI_SUCCESS);
			//res = MPI_Get_count(&status, MPI_CHAR, &count);
			//assert(res == MPI_SUCCESS);
			//assert(count == (int)size_chunk && *buffer_r1 == '@');
			//assert(strlen( buffer_r1 ) == size_chunk);

            struct struct_pread_fastq *td_pread1;
            td_pread1 = malloc (NUM_THREADS * sizeof(struct struct_pread_fastq));
            bef = MPI_Wtime();
            pthread_attr_t attr4;
            pthread_attr_init(&attr4);
            pthread_attr_setstacksize(&attr4, BIG_STACK);
            pthread_attr_setdetachstate(&attr4, 0);

            for( n = 0; n < NUM_THREADS; n++ ){
                td_pread1[n].total_thread = NUM_THREADS;
                td_pread1[n].thread_id = n;
                td_pread1[n].job_rank = rank_num;
                td_pread1[n].offset= offset_chunk;
                td_pread1[n].size = size_chunk;
                td_pread1[n].buffer = buffer_r1;
                td_pread1[n].fd  = fh_r1;
                int ret_code = pthread_create(&threads[n], &attr4, pread_fastq_chunck, (void *)(&td_pread1[n]));
                if (ret_code) {
                    fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", ret_code);
                }
           }
           for(n=0; n<NUM_THREADS; n++)
                pthread_join(threads[n], (void *)(&td_pread1[n]));

            pthread_attr_destroy(&attr4);
            free(td_pread1);

            assert(strlen(buffer_r1) == size_chunk);
            assert(*buffer_r1 == '@');

			///some stats
			aft = MPI_Wtime();
			fprintf(stderr, "%s: num_rank = %d :: read sequences (%.02f)\n", __func__, rank_num, aft - bef);
			total_time_reading_seq += (aft - bef);

			///operations on the number of reads
			bef = MPI_Wtime();
			reads_r1 = all_reads_in_chunk[u1];
			reads = reads_r1;
            bases = 0;
			total_local_reads_aligned += reads/2;
			assert(reads <= INT_MAX);

			/* Parse sequences ... */
			seqs = malloc(reads * sizeof(*seqs));
			assert(seqs != NULL);
            files=1;
			if (file_r1 != NULL) {
				p1 = q1 = buffer_r1; e1 = buffer_r1 + size_chunk; line_number = 0;
							
				while (q1 < e1) {

					if (*q1 != '\n') { q1++; continue; }
					/* We have a full line ... process it */
					*q1 = '\0'; 
					//n will take 4 times in row the same value from 0 with a step of 2
					//it is because the 2 operations are using int types, so it is cast as int twice (within and then outside the brackets)
					n = files * (line_number / 4);
				
					switch (line_number % 4) {
					case 0: /* Line1: Name and Comment */
						assert(*p1 == '@');
						seqs[n].name   = p1 + 1;
						while (*p1 && !isspace((unsigned char)*p1)) {p1++;}
						if (*(p1-2) == '/' && isdigit((unsigned char)*(p1-1))) {*(p1-2) = '\0';}
						if (*p1) {*p1++ = '\0';}
						seqs[n].comment = (copy_comment != 0) ? p1 : NULL;
						seqs[n].sam = NULL;
						break;
					case 1: /* Line2: Sequence */
						seqs[n].seq = p1;
						seqs[n].l_seq = q1 - p1;
						bases += seqs[n].l_seq;
						break;
					case 2: /* Line3: Ignored */
						assert(*p1 == '+');
						break;
					case 3: /* Line4: Quality */
						seqs[n].qual = p1;
						break; }
					p1 = ++q1; 
					line_number++; 
				}
			}
			assert( (line_number/4) == reads_r1);
					
			aft = MPI_Wtime();

			fprintf(stderr, "rank %d :: %s: parsed %zu read of sequences (%ld bp) in (%.02f)\n", rank_num, __func__, reads, (long)bases, aft - bef);
			total_time_parsing += (aft - bef);
			//fprintf(stderr, "rank %d ::: [M::%s] read %zu sequences (%ld bp)...\n",rank_num , __func__, reads, (long)bases);
			/* Datas computation ... */
			bef = MPI_Wtime();
			fprintf(stderr, "rank %d ::::  Call memseqs with count sequences %zu \n", rank_num, reads);
			//mem_process_seqs(opt, indix.bwt, indix.bns, indix.pac, 0, (int)reads, seqs, pes0);
			ret->seqs   = seqs;
            ret->n_seqs     = reads;
            int64_t size = 0;
            for (int i = 0; i < ret->n_seqs; ++i) size += ret->seqs[i].l_seq;
            fprintf(stderr, "\t[0000][ M::%s] read %d sequences (%ld bp)...\n",
                    __func__, ret->n_seqs, (long)size);
            fprintf(stderr, "[RANK %d] 1: Calling process()\n", rank_num);
            if (w2.regs) free(w2.regs);
            if (w2.chain_ar) free(w2.chain_ar);
            if (w2.seedBuf) free(w2.seedBuf);

            nreads = ret->n_seqs;
            w2.regs = (mem_alnreg_v *) calloc(nreads, sizeof(mem_alnreg_v));
            w2.chain_ar = (mem_chain_v*) malloc (nreads * sizeof(mem_chain_v));
            w2.seedBuf = (mem_seed_t *) calloc(sizeof(mem_seed_t),nreads * AVG_SEEDS_PER_READ);
            w2.nreads=ret->n_seqs;
            mem_process_seqs(opt,
                             aux.n_processed,
                             ret->n_seqs,
                             ret->seqs,
                             aux.pes0,
                             w2);
            aux.n_processed += ret->n_seqs;
            aft = MPI_Wtime();
			fprintf(stderr, "rank %d :: %s: computed mappings (%.02f)\n", rank_num, __func__, aft - bef);

			total_time_mapping += (aft - bef);

			/* Write results ... */
            if ( write_format == 2){
                bef = MPI_Wtime();
                pthread_attr_t attr;
                pthread_attr_init(&attr);
                pthread_attr_setstacksize(&attr, BIG_STACK);
                pthread_attr_setdetachstate(&attr, 0);

                struct struct_data_thread *td;
                td = calloc (opt->n_threads, sizeof(struct_data_thread));

                int rest = reads%NUM_THREADS;
                int quot = reads/NUM_THREADS;
                pthread_t threads[NUM_THREADS];

                int ret_code = 0;
                for ( n = 0; n < NUM_THREADS; n++){
                    td[n].seqs_thr = seqs;
                    td[n].begin_index = n * quot;
                    td[n].end_index = ((n+1) * quot);
                    if ( n == (NUM_THREADS - 1)) td[n].end_index = reads;
                    td[n].file_desc = fh_out;
                    ret_code = pthread_create(&threads[n], &attr, write_sam_mt, (void *)(&td[n]));
                }

                for (n = 0; n < NUM_THREADS; n++)
                    pthread_join(threads[n], (void *)(&td[n]));

                pthread_attr_destroy(&attr);
                free(td);
                free(seqs);
                free(buffer_r1);

            }
            if ( write_format == 0 ) {

                int ret_code = 0;
                struct thread_data_compress *tdc;
                tdc = malloc (NUM_THREADS * sizeof(struct thread_data_compress));
                bef = MPI_Wtime();
                pthread_attr_init(&attr);
                pthread_attr_setstacksize(&attr, SMALL_STACK);
                pthread_attr_setdetachstate(&attr, 0);

                for( n = 0; n < NUM_THREADS; n++ ){
                    tdc[n].total_thread = NUM_THREADS;
                    tdc[n].thread_id = n;
                    tdc[n].seqs_thr = seqs;
                    tdc[n].job_rank = rank_num;
                    tdc[n].total_reads = reads;
                    tdc[n].thr_comp_sz = 0;
                    tdc[n].comp_level = compression_level;
                    tdc[n].compressed_buffer_thread  = 0;
                    tdc[n].fh_out = fh_out;
                    ret_code = pthread_create(&threads[n], &attr, compress_and_write_bgzf_thread, (void *)(&tdc[n]));
                    if (ret_code) {
                        fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", ret_code);
                    }
                 }
                 for(n=0; n<NUM_THREADS; n++)
                    pthread_join(threads[n], (void *)(&tdc[n]));
                pthread_attr_destroy(&attr);
                free(tdc);
                for (n = 0; n < reads; n++) free(seqs[n].sam);
                free(seqs);
                free(buffer_r1);
            }

            if ( write_format == 1 ) {

                int ret_code = 0;
                struct thread_data_compress *tdc;
                tdc = malloc (NUM_THREADS * sizeof(struct thread_data_compress));
                bef = MPI_Wtime();
                pthread_attr_init(&attr);
                pthread_attr_setstacksize(&attr, SMALL_STACK);
                pthread_attr_setdetachstate(&attr, 0);

                for( n = 0; n < NUM_THREADS; n++ ){
                    tdc[n].total_thread = NUM_THREADS;
                    tdc[n].thread_id = n;
                    tdc[n].seqs_thr = seqs;
                    tdc[n].job_rank = rank_num;
                    tdc[n].total_reads = reads;
                    tdc[n].thr_comp_sz = 0;
                    tdc[n].comp_level = compression_level;
                    tdc[n].compressed_buffer_thread  = 0;
                    tdc[n].fh_out = fh_out;
                    ret_code = pthread_create(&threads[n], &attr, compress_and_write_bam_thread, (void *)(&tdc[n]));
                    if (ret_code) {
                        fprintf(stderr, "ERROR; return code from pthread_create() is %d\n", ret_code);
                    }
                 }
                 for(n=0; n<NUM_THREADS; n++)
                    pthread_join(threads[n], (void *)(&tdc[n]));
                pthread_attr_destroy(&attr);
                free(tdc);
                for (n = 0; n < reads; n++) free(seqs[n].sam);
                free(seqs);
                free(buffer_r1);
            }

		aft = MPI_Wtime();
		fprintf(stderr, "rank: %d :: %s: wrote results (%.02f) \n", rank_num, __func__, aft - bef);
		total_time_writing += (aft - bef);
		fprintf(stderr, "rank: %d :: finish for chunck %zu \n", rank_num, u1);
            	u1 += proc_num;

	} //end for loop on chunks

         MPI_Barrier(MPI_COMM_WORLD);

        /* Dealloc memory allcoated in the header section */   
        if ( (write_format == 1) && (rank_num == 0)) {

            res = MPI_File_write_shared(fh_out, magic, 28, MPI_BYTE, &status);
            assert(res == MPI_SUCCESS);
            res = MPI_Get_count(&status, MPI_BYTE, &count);
            assert(res == MPI_SUCCESS);
            assert(count == 28);
        }

        free(w2.chain_ar);
        free(w2.regs);
        free(w2.seedBuf);
    
        for(int l=0; l<nthreads; l++) {
            _mm_free(w2.mmc.seqBufLeftRef[l*CACHE_LINE]);
            _mm_free(w2.mmc.seqBufRightRef[l*CACHE_LINE]);
            _mm_free(w2.mmc.seqBufLeftQer[l*CACHE_LINE]);
            _mm_free(w2.mmc.seqBufRightQer[l*CACHE_LINE]);
        }

        for(int l=0; l<nthreads; l++) {
            free(w2.mmc.seqPairArrayAux[l]);
            free(w2.mmc.seqPairArrayLeft128[l]);
            free(w2.mmc.seqPairArrayRight128[l]);
        }

        for(int l=0; l<nthreads; l++) {
            _mm_free(w2.mmc.matchArray[l]);
            free(w2.mmc.min_intv_ar[l]);
            free(w2.mmc.query_pos_ar[l]);
            free(w2.mmc.enc_qdb[l]);
            free(w2.mmc.rid[l]);
            _mm_free(w2.mmc.lim[l]);
        }

        MPI_Win_free(&win_shr_ref);
        MPI_Win_free(&win_shr_pac);

         /* shall we delete?
          * because error in read_index_ele.cpp:52
          */
        //delete(aux.fmi);
        

	}

	/*
	*
	*   Print some statistics
	*
	*/


	after_local_mapping	 = MPI_Wtime();
	total_time_local_mapping = after_local_mapping - before_local_mapping;

	res = MPI_Allreduce(&total_time_local_mapping, &total_time_mapping, 1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	assert ( res == MPI_SUCCESS);

	res = MPI_Allreduce(&total_local_reads_aligned, &total_reads_check, 1,MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
	assert ( res == MPI_SUCCESS);

	res = MPI_Allreduce(&total_time_reading_seq, &grand_total_time_reading_seq, 1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	assert ( res == MPI_SUCCESS);

	res = MPI_Allreduce(&total_time_parsing, &grand_total_time_parsing, 1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	assert ( res == MPI_SUCCESS);

	res = MPI_Allreduce(&total_time_writing, &grand_total_time_writing, 1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	assert ( res == MPI_SUCCESS);
	
	fprintf(stderr, "rank %d :::: total_time_reading_seq = %.02f seconds \n", rank_num, grand_total_time_reading_seq);
	fprintf(stderr, "rank %d :::: total_time_writing = %.02f seconds \n", rank_num, grand_total_time_writing);
	fprintf(stderr, "rank %d :::: total_time_parsing = %.02f  seconds \n", rank_num, grand_total_time_parsing);
	fprintf(stderr, "rank %d :::: total_time_local_mapping = %.02f seconds \n", rank_num, total_time_local_mapping);
	fprintf(stderr, "rank %d :::: total_time_mapping = %.02f seconds \n", rank_num, total_time_mapping);
	fprintf(stderr, "rank %d :::: total_reads_check for all ranks= %zu \n", rank_num, total_reads_check);
	fprintf(stderr, "rank %d :::: total_local_reads_aligned = %zu \n", rank_num, total_local_reads_aligned);
	

	fprintf(stderr, "rank %d :::: finish mappings for reads \n", rank_num);
	if (all_begin_offset_chunk != NULL) free(all_begin_offset_chunk);
	if (all_chunk_size != NULL) free(all_chunk_size);
	if (all_reads_in_chunk != NULL) free(all_reads_in_chunk);

    if (all_begin_offset_chunk_2 != NULL) free(all_begin_offset_chunk_2);
    if (all_chunk_size_2 != NULL) free(all_chunk_size_2);
    if (all_reads_in_chunk_2 != NULL) free(all_reads_in_chunk_2);

	if (opt != NULL) free(opt);

	if (file_r2 != NULL) {
		res = MPI_File_close(&fh_r2);
		assert(res == MPI_SUCCESS);
	}
	if (file_r1 != NULL) {
		res = MPI_File_close(&fh_r1);
		assert(res == MPI_SUCCESS);
	}

	res = MPI_File_close(&fh_out);
	assert(res == MPI_SUCCESS);

	res = MPI_Finalize();
	assert(res == MPI_SUCCESS);
	
	return 0;
}

