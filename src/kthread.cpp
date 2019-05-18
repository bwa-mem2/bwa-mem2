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

Authors: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
         Heng Li <hli@jimmy.harvard.edu>.
*****************************************************************************************/

#include "kthread.h"
// #include <omp.h>
#include <stdio.h>

extern uint64_t proc_freq, tprof[LIM_R][LIM_C];
extern int nthreads, affy[256];
int g_itr;

static inline long steal_work(kt_for_t *t)
{
	int i, min_i = -1;
	long k, min = LONG_MAX;
	for (i = 0; i < t->n_threads; ++i)
		if (min > t->w[i].i) min = t->w[i].i, min_i = i;
	k = __sync_fetch_and_add(&t->w[min_i].i, t->n_threads);
	// return k >= t->n? -1 : k;
	return k*BATCH_SIZE >= t->n? -1 : k;
}

#if 0  // GSHARED
static void *ktf_worker(void *data)
{
	ktf_worker_t *w = (ktf_worker_t*)data;
	long i, val = 0;
	for (;;) {
		// i = __sync_fetch_and_add(&w->i, w->t->n_threads);
		i = __sync_fetch_and_add(&g_itr, 1);
		int st = i * BATCH_SIZE;
		if (st >= w->t->n) break;
		int ed = (i + 1) * BATCH_SIZE < w->t->n? (i + 1) * BATCH_SIZE : w->t->n;
		w->t->func(w->t->data, st, ed-st, w->i);
	}
	pthread_exit(0);
}

#else

static void *ktf_worker(void *data)
{
	ktf_worker_t *w = (ktf_worker_t*)data;
	long i, val = 0;
	int tid = w->i;
	// fprintf(stderr, "i: %d, CPU: %d\n", tid , sched_getcpu());
	
	for (;;) {
		i = __sync_fetch_and_add(&w->i, w->t->n_threads);
		// i = __sync_fetch_and_add(&g_itr, 1);
		int st = i * BATCH_SIZE;
		if (st >= w->t->n) break;
		int ed = (i + 1) * BATCH_SIZE < w->t->n? (i + 1) * BATCH_SIZE : w->t->n;
		w->t->func(w->t->data, st, ed-st, tid);
	}

	while ((i = steal_work(w->t)) >= 0) {
		int st = i * BATCH_SIZE;
		int ed = (i + 1) * BATCH_SIZE < w->t->n? (i + 1) * BATCH_SIZE : w->t->n;
		w->t->func(w->t->data, st, ed-st, tid);
	}
	pthread_exit(0);
}
#endif

void kt_for(void (*func)(void*, int, int, int), void *data, int n)
{
	int i;
	kt_for_t t;
	pthread_t *tid;
	t.func = func, t.data = data, t.n_threads = nthreads, t.n = n;
	t.w = (ktf_worker_t*) malloc (nthreads * sizeof(ktf_worker_t));
	tid = (pthread_t*) malloc (nthreads * sizeof(pthread_t));
	for (i = 0; i < nthreads; ++i)
		t.w[i].t = &t, t.w[i].i = i;

	pthread_attr_t attr;
    cpu_set_t cpus;
    pthread_attr_init(&attr);
	
	// printf("getcpu: %d\n", sched_getcpu());
	g_itr = 0;
	for (i = 0; i < nthreads; ++i) {
#if 1
		CPU_ZERO(&cpus);
		// CPU_SET(i, &cpus);
		CPU_SET(affy[i], &cpus);
		pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpus);	
		pthread_create(&tid[i], &attr, ktf_worker, &t.w[i]);
#else
		pthread_create(&tid[i], NULL, ktf_worker, &t.w[i]);
#endif
	}
	for (i = 0; i < nthreads; ++i) pthread_join(tid[i], 0);

    free(t.w);
	free(tid);
}
