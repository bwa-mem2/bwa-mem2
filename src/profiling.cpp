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
// #if MPI_ENABLED
// #include <mpi.h>
// #endif


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
	fprintf(stderr, "No. of OMP threads: %d\n", nthreads);
	fprintf(stderr, "Processor is runnig @%lf MHz\n", proc_freq*1.0/1e6);
	fprintf(stderr, "Runtime profile:\n");

	fprintf(stderr, "\n\t Time taken for main_mem function: %0.2lf Sec\n\n",
			tprof[MEM][0]*1.0/proc_freq);

	fprintf(stderr, "\tIO times (sec) :\n");
	find_opt(tprof[READ_IO], 1, &max, &min, &avg);
	fprintf(stderr, "\tReading IO time (reads) avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

	find_opt(tprof[SAM_IO], 1, &max, &min, &avg);
	fprintf(stderr, "\tWriting IO time (SAM) avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

	find_opt(tprof[REF_IO], 1, &max, &min, &avg);
	fprintf(stderr, "\tReading IO time (Reference Genome) avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

	find_opt(tprof[FMI], 1, &max, &min, &avg);
	fprintf(stderr, "\tIndex read time avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

	fprintf(stderr, "\n\tOverall time (sec) (Excluding Index reading time):\n");
	// find_opt(tprof[PROCESS], 1, &max, &min, &avg);
	fprintf(stderr, "\tPROCESS() (Total compute time + (read + SAM) IO time) : %0.2lf\n",
			tprof[PROCESS][0]*1.0/proc_freq);

	find_opt(tprof[MEM_PROCESS2], 1, &max, &min, &avg);
	fprintf(stderr, "\tMEM_PROCESS_SEQ() (Total compute time (Kernel + SAM)), avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

	fprintf(stderr, "\n\t SAM Processing time (sec):\n");
	find_opt(tprof[WORKER20], 1, &max, &min, &avg);
	fprintf(stderr, "\t--WORKER_SAM avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
	
#if HIDE
	find_opt(tprof[SAM1], 1, &max, &min, &avg);
	fprintf(stderr, "\t\tWORKER_SAM1 avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
	find_opt(tprof[SAM2], 1, &max, &min, &avg);
	fprintf(stderr, "\t\tWORKER_SAM2 avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
	
	find_opt(tprof[ALIGN1], 1, &max, &min, &avg);
	fprintf(stderr, "\t\t\tWORKER_ALIGN1 avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
	
	find_opt(tprof[SAM3], 1, &max, &min, &avg);
	fprintf(stderr, "\t\tWORKER_SAM3 avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
#endif
	
	fprintf(stderr, "\n\tKernels' compute time (sec):\n");
	find_opt(tprof[WORKER10], 1, &max, &min, &avg);
	fprintf(stderr, "\tTotal kernel (smem+sal+bsw) time avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
	
#if HIDE
	find_opt(tprof[MEM_ALN_M1], nthreads, &max, &min, &avg);
	fprintf(stderr, "\t\tMEM_ALN_CHAIN_FLT avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
	
	find_opt(tprof[MEM_ALN_M2], nthreads, &max, &min, &avg);
	fprintf(stderr, "\t\tMEM_ALN_CHAIN_SEED avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
#endif
	
	find_opt(tprof[MEM_COLLECT], nthreads, &max, &min, &avg);
	fprintf(stderr, "\t\tSMEM compute avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

#if HIDE
	find_opt(tprof[MEM_CHAIN], nthreads, &max, &min, &avg);
	fprintf(stderr, "\t\tMEM_CHAIN avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
#endif
	
	find_opt(tprof[MEM_SA_BLOCK], nthreads, &max, &min, &avg);
	fprintf(stderr, "\t\tSAL compute avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
	
#if HIDE
	find_opt(tprof[MEM_SA], nthreads, &max, &min, &avg);
	fprintf(stderr, "\t\t\t\tMEM_SA avg: %0.2lf, (%0.2lf, %0.2lf)\n\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
#endif
	
	// printf("\n\t BSW compute time (sec):\n");
	find_opt(tprof[MEM_ALN2], nthreads, &max, &min, &avg);
	fprintf(stderr, "\t\tBSW time, avg: %0.2lf, (%0.2lf, %0.2lf)\n",
			avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

#if HIDE
	fprintf(stderr, "\nSTATSV\n");
	fprintf(stderr, "%0.2lf\n",	tprof[PROCESS][0]*1.0/proc_freq);
	find_opt(tprof[MEM_PROCESS2], 1, &max, &min, &avg);
	fprintf(stderr, "%0.2lf\n",	avg*1.0/proc_freq);
	find_opt(tprof[WORKER10], 1, &max, &min, &avg);
	fprintf(stderr, "%0.2lf\n",	avg*1.0/proc_freq);
	find_opt(tprof[MEM_COLLECT], nthreads, &max, &min, &avg);
	fprintf(stderr, "%0.2lf\n",	avg*1.0/proc_freq);
	find_opt(tprof[MEM_SA_BLOCK], nthreads, &max, &min, &avg);
	fprintf(stderr, "%0.2lf\n",	avg*1.0/proc_freq);
	find_opt(tprof[MEM_SA], nthreads, &max, &min, &avg);
	fprintf(stderr, "%0.2lf\n",	avg*1.0/proc_freq);
	double val = 0;
	find_opt(tprof[MEM_ALN2_UP], nthreads, &max, &min, &avg);
	fprintf(stderr, "%0.2lf\n",	avg*1.0/proc_freq);
	val += avg;
	find_opt(tprof[CLEFT], nthreads, &max, &min, &avg);
	fprintf(stderr, "%0.2lf\n",	avg*1.0/proc_freq);
	val += avg;
	find_opt(tprof[CRIGHT], nthreads, &max, &min, &avg);
	fprintf(stderr, "%0.2lf\n",	avg*1.0/proc_freq);
	val += avg;
	find_opt(tprof[MEM_ALN2_DOWN], nthreads, &max, &min, &avg);
	fprintf(stderr, "%0.2lf\n",	avg*1.0/proc_freq);
	val += avg;
	fprintf(stderr, "%0.2lf\n",	val*1.0/proc_freq);
	fprintf(stderr, "%0.2lf\n",	(tprof[REF_IO][0] + tprof[FMI][0])*1.0/proc_freq);
	fprintf(stderr, "%0.2lf\n",	tprof[READ_IO][0]*1.0/proc_freq);
#endif
	// printf("\tMemory usage (GB):\n");
	// printf("\tAvg: %0.2lf, Peak: %0.2lf\n", tprof[PE21][0]*1.0/1e9, tprof[PE22][0]*1.0/1e9);
		
	return 1;
}

