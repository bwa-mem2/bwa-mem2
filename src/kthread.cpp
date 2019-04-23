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

#include "kthread.h"
#include <omp.h>
#include <stdio.h>

extern uint64_t proc_freq, tprof[LIM_R][LIM_C];

void kt_for(void (*func)(void*, int, int, int), void *data, int n)
{
	/* Note: we dont need ktf_worker anymore!
	   omp threading of num threads */
#if 0
	
	int i = 0;
#pragma omp parallel for schedule(dynamic)
	for (i=0; i<n; i+= BATCH_SIZE)
	{
		// uint64_t timG = _rdtsc();
		int tid = omp_get_thread_num();
		int st = i;
		int ed = i + BATCH_SIZE < n? i + BATCH_SIZE : n;
		// set thread memory shared in data
		// printf("[%0.4d] st :%d, ed: %d\n", tid, st, ed-st);
		func(data, st, ed-st, tid);
		// tprof[MEM_SA_BLOCK][tid] += _rdtsc() - timG;
	}
	
#else

	
	int i = 0;
#pragma omp parallel
	{
		int tid = omp_get_thread_num();
#pragma omp for schedule(dynamic)
		for (i=0; i<n; i+= BATCH_SIZE)
		{
			int st = i;
			int ed = i + BATCH_SIZE < n? i + BATCH_SIZE : n;
			// set thread memory shared in data
			// printf("[%0.4d] st :%d, ed: %d\n", tid, st, ed-st);
			func(data, st, ed-st, tid);
		}
	}	
#endif
}

