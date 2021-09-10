/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

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

Authors: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>.
*****************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <limits>
#include <unistd.h>
#include "bandedSWA.h"


#define DEFAULT_MATCH 1
#define DEFAULT_MISMATCH 4
#define DEFAULT_OPEN 6
#define DEFAULT_EXTEND 1
#define DEFAULT_AMBIG -1

FILE *fksw;
double freq = 2.3*1e9;
int w_match, w_mismatch, w_open, w_extend, w_ambig;
uint64_t SW_cells;
char *pairFileName;
FILE *pairFile;
int8_t h0 = 0;
double clock_freq;
uint64_t prof[10][112], data, SW_cells2;

void bwa_fill_scmat(int a, int b, int ambig, int8_t mat[25]) {
	int i, j, k;
	for (i = k = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j)
			mat[k++] = i == j? a : -b;
		mat[k++] = ambig; // ambiguous base
	}
	for (j = 0; j < 5; ++j) mat[k++] = ambig;
}


void parseCmdLine(int argc, char *argv[])
{
    int i;
    w_match = DEFAULT_MATCH;
    w_mismatch = DEFAULT_MISMATCH;
    w_open = DEFAULT_OPEN;
    w_extend = DEFAULT_EXTEND;
	w_ambig = DEFAULT_AMBIG;
	
    int pairFlag = 0;
    for(i = 1; i < argc; i+=2)
    {
        if(strcmp(argv[i], "-match") == 0)
        {
            w_match = atoi(argv[i + 1]);
        }
        if(strcmp(argv[i], "-mismatch") == 0) //penalty, +ve number
        {
            w_mismatch = atoi(argv[i + 1]);
        }
		if(strcmp(argv[i], "-ambig") == 0)
        {
            w_ambig = atoi(argv[i + 1]);
        }

        if(strcmp(argv[i], "-gapo") == 0)
        {
            w_open = atoi(argv[i + 1]);
        }
        if(strcmp(argv[i], "-gape") == 0)
        {
            w_extend = atoi(argv[i + 1]);
        }
        if(strcmp(argv[i], "-pairs") == 0)
        {
            pairFileName = argv[i + 1];
            pairFlag = 1;
        }
        if(strcmp(argv[i], "-h0") == 0)
        {
            h0 = atoi(argv[i + 1]);
        }
    }
    if(pairFlag == 0)
    {
        printf("ERROR! pairFileName not specified.\n");
        exit(0);
    }
}

int loadPairs(SeqPair *seqPairArray, uint8_t *seqBufRef, uint8_t* seqBufQer)
{
    int32_t numPairs = 0;
    while(numPairs < MAX_NUM_PAIRS_ALLOC)
    {
		int32_t h0 = 0;
		char temp[10];
		fgets(temp, 10, pairFile);
		sscanf(temp, "%d", &h0);
        //if(!fgets((char *)(seqBuf + numPairs * 2 * MAX_SEQ_LEN), MAX_SEQ_LEN, pairFile))
		if(!fgets((char *)(seqBufRef + numPairs * MAX_SEQ_LEN_REF), MAX_SEQ_LEN_REF, pairFile))
        {
            break;
        }
        //if(!fgets((char *)(seqBuf + (numPairs * 2 + 1) * MAX_SEQ_LEN), MAX_SEQ_LEN, pairFile))
		if(!fgets((char *)(seqBufQer + numPairs * MAX_SEQ_LEN_QER), MAX_SEQ_LEN_QER, pairFile))	
        {
            printf("ERROR! Odd number of sequences in %s\n", pairFileName);
            break;
        }

        SeqPair sp;
        sp.id = numPairs;
        // sp.seq1 = seqBuf + numPairs * 2 * MAX_SEQ_LEN;
        // sp.seq2 = seqBuf + (numPairs * 2 + 1) * MAX_SEQ_LEN;
        sp.len1 = strnlen((char *)(seqBufRef + numPairs * MAX_SEQ_LEN_REF), MAX_SEQ_LEN_REF) - 1;
        sp.len2 = strnlen((char *)(seqBufQer + numPairs * MAX_SEQ_LEN_QER), MAX_SEQ_LEN_QER) - 1;
		sp.h0 = h0;
        // sp.score = 0;

		//uint8_t *seq1 = seqBuf + numPairs * 2 * MAX_SEQ_LEN;
        //uint8_t *seq2 = seqBuf + (numPairs * 2 + 1) * MAX_SEQ_LEN;
		uint8_t *seq1 = seqBufRef + numPairs * MAX_SEQ_LEN_REF;
        uint8_t *seq2 = seqBufQer + numPairs * MAX_SEQ_LEN_QER;
		sp.idr =  numPairs * MAX_SEQ_LEN_REF;
		sp.idq =  numPairs * MAX_SEQ_LEN_QER;
		
		for (int l=0; l<sp.len1; l++)
			seq1[l] -= 48;
		for (int l=0; l<sp.len2; l++)
			seq2[l] -= 48;

		// sp.idr = sp.idq = 0;
		sp.seqid = sp.regid = sp.score = sp.tle = sp.gtle = sp.qle = -1;
		sp.gscore = sp.max_off = -1;
		
        seqPairArray[numPairs] = sp;
        numPairs++;
        // SW_cells += (sp.len1 * sp.len2);
    }
    // fclose(pairFile);
    return numPairs;
}
// profiling stats
uint64_t find_stats(uint64_t *val, int nt, double &min, double &max, double &avg) {
	min = 1e10;
	max = 0;
	avg = 0;
	for (int i=0; i<nt; i++) {
		avg += val[i];
		if (max < val[i]) max = val[i];
		if (min > val[i]) min = val[i];		
	}
	avg /= nt;

	return 1;
}

int main(int argc, char *argv[])
{
	fksw = fopen("fksw.txt", "w");
	if (argc < 3) {
		printf("usage: <exec> -pairs <InSeqFile>!!\n");
		exit(0);
	}
	
	clock_freq = _rdtsc();
	sleep(1);
	clock_freq = _rdtsc() - clock_freq;
	
    parseCmdLine(argc, argv);
    SeqPair *seqPairArray = (SeqPair *)_mm_malloc(MAX_NUM_PAIRS * sizeof(SeqPair), 64);
	// OutScore *outScoreArray = (OutScore *)_mm_malloc(MAX_NUM_PAIRS * sizeof(OutScore), 64);
	
	pairFile = fopen(pairFileName, "r");	
    if(pairFile == NULL)
    {
        fprintf(stderr, "Could not open file: %s\n", pairFileName);
        exit(0);
    }

    uint8_t *seqBufRef = (uint8_t *)_mm_malloc(MAX_SEQ_LEN_REF * MAX_NUM_PAIRS_ALLOC *
											sizeof(int8_t), 64);
    uint8_t *seqBufQer = (uint8_t *)_mm_malloc(MAX_SEQ_LEN_QER * MAX_NUM_PAIRS_ALLOC *
											sizeof(int8_t), 64);

	int8_t mat[25];
	bwa_fill_scmat(w_match, w_mismatch, w_ambig, mat);
	int zdrop = 100, w = 100, end_bonus = 5;
 	BandedPairWiseSW *bsw = new BandedPairWiseSW(w_open, w_extend, w_open, w_extend,
												 zdrop, end_bonus, mat,
												 w_match, w_mismatch, 1);
	
    //int64_t endTick, pTotalTicks = 0;
    int64_t startTick, totalTicks = 0, readTim = 0;
#if __AVX2__
    int64_t myTicks = 0;
	int numThreads = 1;
#endif
	// SW_cells = 0;
	int32_t numPairs = 0, totalPairs = 0;
	
	uint64_t tim = _rdtsc();
	numPairs = loadPairs(seqPairArray, seqBufRef, seqBufQer);
	readTim += _rdtsc() - tim;
	
	while(numPairs)
	{
        //int i;
		// for(i = 0; i < numPairs; i++)
		// {
		// 	SeqPair *p = seqPairArray + i;
		// 	int qle, tle, gtle, gscore, max_off;
		// 	int32_t h0 = p->h0;
	
		// 	int score = bsw.scalarBandedSWA(p->len2, p->seq2, p->len1,
		// 									p->seq1, w, h0, &qle, &tle,
		// 									&gtle, &gscore, &max_off);		
		// 	totalPairs++;
		// }
		startTick = __rdtsc();		
		
#if __AVX2__
		// printf("Executing AVX2 vector code...\n");
		uint64_t timM = _rdtsc();
		// bsw->getScores8(seqPairArray, seqBufRef, seqBufQer, numPairs, 1, w);
		bsw->getScores16(seqPairArray, seqBufRef, seqBufQer, numPairs, 1, w);
		//bsw->scalarBandedSWAWrapper(seqPairArray, seqBufRef, seqBufQer, numPairs, 1, w);
		myTicks += _rdtsc() - timM;
#else
		// printf("Executing scalar code...\n");
		bsw->scalarBandedSWAWrapper(seqPairArray, seqBufRef, seqBufQer, numPairs, 1, w);
#endif
		totalTicks += __rdtsc() - startTick;
		
		uint64_t tim = _rdtsc();
		numPairs = loadPairs(seqPairArray, seqBufRef, seqBufQer);
		readTim += _rdtsc() - tim;
		
		totalPairs += numPairs;
	}
	//totalTicks += endTick - startTick;
#if __AVX2__
		printf("Executed AVX2 vector code...\n");
#else
		printf("Executed scalar code...\n");
#endif

	tim = _rdtsc();
	sleep(1);
	freq = _rdtsc() - tim;
	
	printf("Processor freq: %0.2lf MHz\n", freq/1e6);
	
#if __AVX2__
	//int64_t myTicks = bsw->getTicks();
	printf("Read time  = %0.2lf\n", readTim/freq);
    printf("Overall SW cycles = %ld, %0.2lf\n", myTicks, myTicks*1.0/freq);
	printf("Total Pairs processed: %d\n", totalPairs);
    printf("SW cells(T)  = %ld\n", SW_cells);
	printf("SW cells(||)  = %ld\n", SW_cells2);
    printf("SW GCUPS  = %lf\n", SW_cells * freq/1e9 / myTicks);

	printf("More stats:\n");
	// double freq = 2.3*1e9;
	double min, max, avg;
	find_stats(prof[1], numThreads, min, max, avg);
	printf("Time in pre-processing: %0.2lf (%0.2lf, %0.2lf)\n",
		   avg*1.0/freq, min*1.0/freq, max*1.0/freq);
	find_stats(prof[0], numThreads, min, max, avg);
	printf("Time spent in smithWaterman(): %0.2lf (%0.2lf, %0.2lf)\n",
		   avg*1.0/freq, min*1.0/freq, max*1.0/freq);

	printf("\nDebugging info:\n");
	printf("Time taken for DP loop: %0.2lf\n", prof[DP][0]*1.0/freq);
	printf("Time taken for DP loop upper part: %0.2lf\n", prof[DP3][0]*1.0/freq);	
	printf("Time taken for DP inner loop: %0.2lf\n", prof[DP1][0]*1.0/freq);
	printf("Time taken for DP loop lower part: %0.2lf\n", prof[DP2][0]*1.0/freq);	

#else
	fprintf(stderr, "Processor is running @%0.2lf Mhz\n", clock_freq/1e6);
    //printf("SW cycles = %ld\n", pwsw->getTicks());
	fprintf(stderr, "SW cycles = %ld, time = %0.6lf Sec\n",
			totalTicks, totalTicks*1.0/clock_freq);
	fprintf(stderr, "Total pairs computed: %d\n", totalPairs);
	fprintf(stderr, "Total cells computed: %ld\n", bsw->SW_cells);
#endif
	/**** free memory *****/
	_mm_free(seqPairArray);
	_mm_free(seqBufRef);
	_mm_free(seqBufQer);
	delete bsw;
	
	fclose(pairFile);
	fclose(fksw);
	return 1;
}
