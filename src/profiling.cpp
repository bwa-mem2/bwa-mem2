/*************************************************************************************
                    GNU GENERAL PUBLIC LICENSE
           		      Version 3, 29 June 2007

BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
Copyright (C) 2019  Vasimuddin Md, Sanchit Misra, Intel Corporation, Heng Li.
    
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License at https://www.gnu.org/licenses/ for more details.


TERMS AND CONDITIONS FOR DISTRIBUTION OF THE CODE
                                             
1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 
3. Neither the name of Intel Corporation nor the names of its contributors may
   be used to endorse or promote products derived from this software without
   specific prior written permission.

Authors: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>.
*****************************************************************************************/

#include <stdio.h>
#include "macro.h"
#include <stdint.h>
#if MPI_ENABLED
#include <mpi.h>
#endif


extern uint64_t proc_freq, tprof[LIM_R][LIM_C];
extern int nthreads;
extern int myrank, num_ranks;

int find_opt(uint64_t *a, int len, uint64_t *max, uint64_t *min, double *avg)
{
	*max = 0;
	*min = 1e15;
	*avg = 0;

	int i=0;
	for (i=0; i<len; i++)
	{
		if (a[i] > *max) *max = a[i];
		if (a[i] < *min) *min = a[i];
		*avg += a[i];
	}
	*avg /= len;

	return 1;
}

int display_stats()
{
	uint64_t max, min;
	double avg;
	printf( "No. of OMP threads: %d\n", nthreads);
	printf( "Processor is runnig @%lf MHz\n", proc_freq*1.0/1e6);
	printf("Runtime profile:\n");

	printf( "\n\t Time taken for main_mem function: %0.2lf Sec\n\n",
			tprof[MEM][0]*1.0/proc_freq);

	printf("\tIO times (sec) :\n");
	find_opt(tprof[READ_IO], 1, &max, &min, &avg);
	printf( "\tReading IO time (reads) avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

	find_opt(tprof[REF_IO], 1, &max, &min, &avg);
	printf( "\tReading IO time (Reference Genome) avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

	// find_opt(tprof[CONVERT], 1, &max, &min, &avg);
	// printf( "\tTime spent in ref. seq. conversion, avg: %0.2lf, (%0.2lf, %0.2lf)\n",
	// 		avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

	find_opt(tprof[FMI], 1, &max, &min, &avg);
	printf( "\nIndex read time avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

	printf("\n\tOverall time (sec) (Excluding Index reading time):\n");
	// find_opt(tprof[PROCESS], 1, &max, &min, &avg);
	printf( "\tPROCESS() (Total compute time + (read + SAM) IO time) : %0.2lf\n",
			tprof[PROCESS][0]*1.0/proc_freq);

	find_opt(tprof[MEM_PROCESS2], 1, &max, &min, &avg);
	printf( "\tMEM_PROCESS_SEQ() (Total compute time (Kernel + SAM)), avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

	printf("\n\t SAM Processing time (sec):\n");
	find_opt(tprof[WORKER20], 1, &max, &min, &avg);
	printf( "\t--WORKER_SAM avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
	
#if HIDE
	find_opt(tprof[SAM1], 1, &max, &min, &avg);
	printf( "\t\tWORKER_SAM1 avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
	find_opt(tprof[SAM2], 1, &max, &min, &avg);
	printf( "\t\tWORKER_SAM2 avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
	
	find_opt(tprof[ALIGN1], 1, &max, &min, &avg);
	printf( "\t\t\tWORKER_ALIGN1 avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
	
	find_opt(tprof[SAM3], 1, &max, &min, &avg);
	printf( "\t\tWORKER_SAM3 avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
#endif
	
	printf("\n\tKernels' compute time (sec):\n");
	find_opt(tprof[WORKER10], 1, &max, &min, &avg);
	printf( "\tTotal kernel (smem+sal+bsw) time avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
#if HIDE
	//find_opt(tprof[MEM_ALN_M1], nthreads, &max, &min, &avg);
	//printf( "\t\tMEM_ALN_CHAIN_FLT avg: %0.2lf, (%0.2lf, %0.2lf)\n",
	//		avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
	//
	//find_opt(tprof[MEM_ALN_M2], nthreads, &max, &min, &avg);
	//printf( "\t\tMEM_ALN_CHAIN_SEED avg: %0.2lf, (%0.2lf, %0.2lf)\n",
	//		avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
#endif
	
	find_opt(tprof[MEM_COLLECT], nthreads, &max, &min, &avg);
	printf( "\t\tSMEM compute avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

#if HIDE
	//find_opt(tprof[MEM_CHAIN], nthreads, &max, &min, &avg);
	//printf( "\t\tMEM_CHAIN avg: %0.2lf, (%0.2lf, %0.2lf)\n",
	//		avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
#endif
	
	find_opt(tprof[MEM_SA_BLOCK], nthreads, &max, &min, &avg);
	printf( "\t\tSAL compute avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
#if HIDE
	find_opt(tprof[MEM_SA], nthreads, &max, &min, &avg);
	printf( "\t\t\t\tMEM_SA avg: %0.2lf, (%0.2lf, %0.2lf)\n\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
#endif
	
	// printf("\n\t BSW compute time (sec):\n");
	find_opt(tprof[MEM_ALN2], nthreads, &max, &min, &avg);
	printf( "\t\tBSW time, avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

#if HIDE
	find_opt(tprof[MEM_ALN2_UP], nthreads, &max, &min, &avg);
	printf( "\t\t\tPRE PROC avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
		
	find_opt(tprof[SORT], nthreads, &max, &min, &avg);
	printf( "\t\t\tSWA SORT avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

	find_opt(tprof[CLEFT], nthreads, &max, &min, &avg);
	printf( "\t\t\tKSW BLOCK L-EXTEND avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

	find_opt(tprof[MEM_ALN2_B], nthreads, &max, &min, &avg);
	printf( "\t\t\t\tKSW L-EXTEND 16 bit lane avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

	find_opt(tprof[MEM_ALN2_D], nthreads, &max, &min, &avg);
	printf( "\t\t\t\tKSW L-EXTEND  8 bit lane avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

	find_opt(tprof[CRIGHT], nthreads, &max, &min, &avg);
	printf( "\t\t\tKSW BLOCK R-EXTEND  avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

	find_opt(tprof[MEM_ALN2_C], nthreads, &max, &min, &avg);
	printf( "\t\t\t\tKSW R-EXTEND 16 bit lane avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

	find_opt(tprof[MEM_ALN2_E], nthreads, &max, &min, &avg);
	printf( "\t\t\t\tKSW R-EXTEND  8 bit lane avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

	find_opt(tprof[MEM_ALN2_DOWN], nthreads, &max, &min, &avg);
	printf( "\t\t\tPOST PROC avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

	find_opt(tprof[MEM_ALN2_DOWN1], nthreads, &max, &min, &avg);
	printf( "\t\t\t\tPOST PROC avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

	find_opt(tprof[POST_SWA], nthreads, &max, &min, &avg);
	printf( "\t\tAfter SWA avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
	
	
	// find_opt(tprof[PE1], nthreads, &max, &min, &avg);
	printf("\tDEBUG: Avg. numPairsLeft int8  %0.2lf\n", tprof[PE1][0]*1.0/tprof[PE2][0]);
	printf("\tDEBUG: Avg. numPairsRight int8 %0.2lf\n", tprof[PE3][0]*1.0/tprof[PE4][0]);
	printf("\tDEBUG: Avg. numPairsLeft >int8  %0.2lf\n", tprof[PE5][0]*1.0/tprof[PE6][0]);
	printf("\tDEBUG: Avg. numPairsRight >int8 %0.2lf\n\n", tprof[PE7][0]*1.0/tprof[PE8][0]);

	printf("\tL: %ld, R:%ld\n", tprof[PE1][0]+tprof[PE5][0], tprof[PE3][0]+tprof[PE7][0] );
	printf("\tL8: %ld, L16: %ld, R8: %ld, R16: %ld\n",
		   tprof[PE1][0], tprof[PE5][0], tprof[PE3][0], tprof[PE7][0] );
	printf("\tDEBUG: Avg. chain size per read: %0.2lf (%ld %ld)\n",
	       tprof[PE11][0]*1.0/tprof[PE12][0], tprof[PE11][0], tprof[PE12][0]);

	printf("\tDEBUG: PE12 (#aln): %ld, PE11 (Total chains): %ld, PE13 (max batch): %ld "
	       "PE14 (#reads): %ld, PE15 (Max chain): %ld, PE16 (#seeds): %ld, "
		   "PE17 (iters): %ld\n",
	       tprof[PE12][0], tprof[PE11][0], tprof[PE13][0], tprof[PE14][0], tprof[PE15][0],
	       tprof[PE16][0], tprof[PE17][0]);
	
	printf("\tPE18: %ld, PE19: %ld, PE20: %ld\n", tprof[PE18][0], tprof[PE19][0], tprof[PE20][0]);
	printf("\tPE21 (malloc): %ld, PE22: %ld, PE23: %ld\n", tprof[PE21][0], tprof[PE22][0], tprof[PE23][0]);
	printf("\tPE24 (malloc): %ld, PE25: %ld, PE26: %ld\n", tprof[PE24][0], tprof[PE25][0], tprof[PE26][0]);
	printf("\tFile size: %ld, Task_size: %ld, BATCH_SIZE: %d\n", tprof[MISC][0], tprof[MISC][1], BATCH_SIZE);

	printf("\nSTATSV\n");
	printf( "%0.2lf\n",	tprof[PROCESS][0]*1.0/proc_freq);
	find_opt(tprof[MEM_PROCESS2], 1, &max, &min, &avg);
	printf( "%0.2lf\n",	avg*1.0/proc_freq);
	find_opt(tprof[WORKER10], 1, &max, &min, &avg);
	printf( "%0.2lf\n",	avg*1.0/proc_freq);
	find_opt(tprof[MEM_COLLECT], nthreads, &max, &min, &avg);
	printf( "%0.2lf\n",	avg*1.0/proc_freq);
	find_opt(tprof[MEM_SA_BLOCK], nthreads, &max, &min, &avg);
	printf( "%0.2lf\n",	avg*1.0/proc_freq);
	find_opt(tprof[MEM_SA], nthreads, &max, &min, &avg);
	printf( "%0.2lf\n",	avg*1.0/proc_freq);
	double val = 0;
	find_opt(tprof[MEM_ALN2_UP], nthreads, &max, &min, &avg);
	printf( "%0.2lf\n",	avg*1.0/proc_freq);
	val += avg;
	find_opt(tprof[CLEFT], nthreads, &max, &min, &avg);
	printf( "%0.2lf\n",	avg*1.0/proc_freq);
	val += avg;
	find_opt(tprof[CRIGHT], nthreads, &max, &min, &avg);
	printf( "%0.2lf\n",	avg*1.0/proc_freq);
	val += avg;
	find_opt(tprof[MEM_ALN2_DOWN], nthreads, &max, &min, &avg);
	printf( "%0.2lf\n",	avg*1.0/proc_freq);
	val += avg;
	printf( "%0.2lf\n",	val*1.0/proc_freq);
	printf( "%0.2lf\n",	(tprof[REF_IO][0] + tprof[FMI][0])*1.0/proc_freq);
	printf( "%0.2lf\n",	tprof[READ_IO][0]*1.0/proc_freq);
#endif
	// printf("\tMemory usage (GB):\n");
	// printf("\tAvg: %0.2lf, Peak: %0.2lf\n", tprof[PE21][0]*1.0/1e9, tprof[PE22][0]*1.0/1e9);
		
	return 1;
}

#if MPI_ENABLED
void mpi_profile()
{

	tprof[MPI_TIME][2] = tprof[MEM_PROCESS2][0];
	tprof[MPI_TIME][3] = tprof[READ_IO][0];
	tprof[MPI_TIME][4] = tprof[OUTPUT][0];
	int num = 10;
	MPI_Reduce(tprof[MPI_TIME], tprof[MPI_TIME_SUM], num,
			   MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(tprof[MPI_TIME], tprof[MPI_TIME_MAX], num,
			   MPI_UINT64_T, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(tprof[MPI_TIME], tprof[MPI_TIME_MIN], num,
			   MPI_UINT64_T, MPI_MIN, 0, MPI_COMM_WORLD);
	
	if (ROOT_) {
		printf("\nNo. of MPI ranks: %d\n", num_ranks);
		int64_t avg, min, max;
		
		printf("MPI performance stats:\n");
		avg = tprof[MPI_TIME_SUM][1]/num_ranks;
		min = tprof[MPI_TIME_MIN][1];
		max = tprof[MPI_TIME_MAX][1];
		printf("\t(Threaded) Kernel execution block time: %0.2lf (%0.2lf, %0.2lf)\n",
			   avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);


		avg = tprof[MPI_TIME_SUM][3]/num_ranks;
		min = tprof[MPI_TIME_MIN][3];
		max = tprof[MPI_TIME_MAX][3];
		printf("\t\t1. Read time (not done using MPI_File_read()): %0.2lf (%0.2lf, %0.2lf)\n",
			   avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

		avg = tprof[MPI_TIME_SUM][2]/num_ranks;
		min = tprof[MPI_TIME_MIN][2];
		max = tprof[MPI_TIME_MAX][2];
		printf("\t\t2. Mem process time: %0.2lf (%0.2lf, %0.2lf)\n",
			   avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

		avg = tprof[MPI_TIME_SUM][4]/num_ranks;
		min = tprof[MPI_TIME_MIN][4];
		max = tprof[MPI_TIME_MAX][4];
		printf("\tMPI (SAM) output time: %0.2lf (%0.2lf, %0.2lf)\n",
			   avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
	}
}
#endif
