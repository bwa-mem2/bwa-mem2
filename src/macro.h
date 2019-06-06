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

Authors: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>
*****************************************************************************************/

#ifndef _MACRO_HPP
#define _MACRO_HPP

#include <stdio.h>

#define VER 0
#define printf_(x,y...)								\
	{												\
		if(x)										\
			fprintf(stderr, y);						\
	}

/* Note: BSW-specific macros are in src/bandedSWA.h file */

#define H0_ -99
#define SEEDS_PER_READ 500           /* Avg seeds per read */
#define MAX_SEEDS_PER_READ 500       /* Max seeds per read */
#define AVG_SEEDS_PER_READ 64       /* Used for storing seeds in chains*/
#define BATCH_SIZE 512               /* Block of reads alloacted to a thread for processing*/
#define BATCH_MUL 20
#define SEEDS_PER_CHAIN 2

#define SEQ_LEN8 128
#define MAX_LINE_LEN 256

#define ROOT_ (myrank == 0)
#define LIM_R 128
#define LIM_C 128


/*** Runtime profiling macros ***/
#define INDEX 0
#define MEM 1
#define MEM2 2
#define MEM3 4
#define SAM1 5
#define SAM2 3
#define SAM3 7
#define MPI_TIME 8
#define MEM_PROCESS10 9
#define MEM_PROCESS2 10
#define READ_IO 11
#define PROCESS 12
#define REF_IO 13
#define PREPROCESS 14
#define CONVERT 15
#define MPI_TIME_SUM 16
#define OUTPUT 17
#define MPI_TIME_MIN 18
#define POST_SWA 19
#define MPI_TIME_MAX 20
#define SAM_IO 21
#define ALIGN1 22

#define KT_FOR 24
#define KTF_WORKER 26
#define WORKER20 28
#define WORKER21 30
#define WORKER10 32
#define WORKER11 34
#define MEM_ALN 36
#define MEM_CHAIN 38
#define MEM_COLLECT 40
#define BWT_REVERSE 41
#define BWA_BUILD 42
#define PACKED 43
#define SA 44
#define BWT_REVERSE_A 45
#define BWT_REVERSE_B 46
#define MEM_SA 47
#define MEM_ALN2 48
#define MEM_ALN2_A 49
#define MEM_ALN2_B 50
#define MEM_ALN2_C 51
#define EXTEND 52
#define FORWARD 53
#define MEM_CHAIN1 54
#define MEM_CHAIN2 55
#define SMEM1 56
#define SMEM2 57
#define SMEM3 58
#define BWT_FORWARD_A 59
#define STR 60
#define MISC 61
#define MEM_ALN2_UP 62
#define BWT_FORWARD_B 63
#define CLEFT 64
#define CRIGHT 65
#define MEM_ALN_M1 66
#define MEM_ALN_M2 67
#define MEM_SA_BLOCK 68
#define SEQ_FETCH 69
#define MEM_ALN2_PRE 70
#define QLEN 71
#define TLEN 72
#define CNT 73
#define WAVG 74
#define WCNT 75
#define WMAX 76
#define WMIN 77
#define KSW 78
#define PE 79
#define PESW 80
#define PESORT 81
#define INTROSORT 82
#define PE1 83
#define PE3 84
#define PE2 85
#define PE4 86
#define PE5 87
#define PE6 88
#define PE7 89
#define PE8 90
#define PE11 91
#define PE12 92
#define PE13 93
#define PE14 94
#define PE15 95
#define PE16 96
#define PE17 97
#define PE18 98
#define PE19 99
#define PE20 100
#define PE21 101
#define PE22 102
#define PE23 103
#define MEM_ALN2_DOWN 104
#define MEM_ALN2_DOWN1 105
#define SORT 106
#define FMI 107
#define MEM_ALN2_D 108
#define MEM_ALN2_E 109
#define PE24 110
#define PE25 111
#define PE26 112


#endif
