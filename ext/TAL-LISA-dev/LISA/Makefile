##/*************************************************************************************
##                           The MIT License
##
##   Copyright (C) 2020 Intel Labs.
##
##   Permission is hereby granted, free of charge, to any person obtaining a copy
##   of this software and associated documentation files (the "Software"), to deal
##   in the Software without restriction, including without limitation the rights
##   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
##   copies of the Software, and to permit persons to whom the Software is
##   furnished to do so, subject to the following conditions:
##
##   The above copyright notice and this permission notice shall be included in all
##   copies or substantial portions of the Software.
##   
##   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
##   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
##   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
##   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
##   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
##   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
##   SOFTWARE.
##
##Contacts: Sanchit Misra <sanchit.misra@intel.com>; Vasimuddin Md <vasimuddin.md@intel.com>
##*****************************************************************************************/

HP= -DNO_HUGE_PAGE
ifeq ($(huge_page), 1)
	HP= -DHUGE_PAGE
endif

LISA_CC = icpc 
LISA_CFLAGS = -DSAIS=1 -std=c++17 -xCORE-AVX512 -Ofast -fopenmp -Wall -Wshadow -Wno-char-subscripts
LISA_CPPFLAGS = -DOUTPUT ${HP} -DNO_DNA_ORD -DREV_COMP -DBWA_MEM_BUG -DNOCHUNK -DVECTORIZE  -DENABLE_PREFETCH -D_64BIT
LISA_CPPFLAGS_EXACT_SEARCH = ${HP} -DNO_DNA_ORD -DBWA_MEM_BUG -DNOCHUNK -DVECTORIZE  -DENABLE_PREFETCH -D_64BIT
LISA_BUILD_RMI_FLAGS = -DBWA_MEM_BUG -DNOCHUNK -xCORE-AVX512  -DVECTORIZE -DENABLE_PREFETCH -D_64BIT 
LISA_INCLUDE = -I src-mix/ -I ../ext/ -I ../ext/safestringlib/include/
LISA_LDLIBS = -lz -L./.. -ltal -L ../ext/safestringlib/ -lsafestring 




VTUNE = -DVTUNE_ANALYSIS -I/swtools/intel/vtune_amplifier/include/ -littnotify -L/swtools/intel/vtune_amplifier/lib64/ 

.cpp.o:                                                                                                                                                   ${LISA_CC} ${LISA_CFLAGS} ${LISA_CPPFLAGS} ${LISA_INCLUDE} -c ${LISA_LDLIBS} $< -o $@

make-obj:
	${LISA_CC} ${LISA_CFLAGS} ${LISA_CPPFLAGS} ${LISA_INCLUDE} -c src-modular/*.cpp ${LISA_LDLIBS} -DPRINT_OUTPUT -o lisalib.o

lisa: smem-lisa exact-search-lisa build-index-forward-only-lisa build-index-with-rev-complement-lisa
smem-lisa: benchmarks/bench-smem-lisa.cpp
	${LISA_CC} ${LISA_CFLAGS} ${LISA_CPPFLAGS} ${LISA_INCLUDE} src-mix/common.o src-mix/chunkEncode.o  src-mix/qbwt-rmi-batched.o benchmarks/bench-smem-lisa.cpp ${LISA_LDLIBS} -DPRINT_OUTPUT -o smem-lisa.o


exact-search-lisa: benchmarks/bench-fixed-len-e2e-match-lisa.cpp
	${LISA_CC} ${LISA_CFLAGS} ${LISA_CPPFLAGS_EXACT_SEARCH} ${LISA_INCLUDE} -DENABLE_PREFETCH_OPT benchmarks/bench-fixed-len-e2e-match-lisa.cpp ${LISA_LDLIBS} -DPRINT_OUTPUT -o exact-search-lisa.o 


build-index-forward-only-lisa: benchmarks/build-lisa-index.cpp
	${LISA_CC} ${LISA_CFLAGS} ${LISA_BUILD_RMI_FLAGS} ${LISA_INCLUDE} benchmarks/build-lisa-index.cpp ${LISA_LDLIBS} -o build-index-forward-only-lisa.o

build-index-with-rev-complement-lisa: benchmarks/build-lisa-index.cpp
	${LISA_CC} ${LISA_CFLAGS} ${LISA_BUILD_RMI_FLAGS} -DREV_COMP ${LISA_INCLUDE} benchmarks/build-lisa-index.cpp ${LISA_LDLIBS} -o build-index-with-rev-complement-lisa.o

clean:
	rm *.o 



# Experimental compilation flags
#TEMP_LISA_CFLAGS = -DSAIS=1 -std=c++17  -Ofast -fopenmp -Wall -Wshadow -Wno-char-subscripts
#TEMP_LISA_CPPFLAGS = -DREV_COMP -DBWA_MEM_BUG -DNOCHUNK   -DENABLE_PREFETCH -D_64BIT
#LISA_CPPFLAGS_SCALAR = -DOUTPUT ${HP} -DNO_DNA_ORD -DREV_COMP -DBWA_MEM_BUG -DNOCHUNK -DNo_VECTORIZE  -DENABLE_PREFETCH -D_64BIT
#LISA_CPPFLAGS_NO_OPT = -DOUTPUT -DNO_DNA_ORD -DREV_COMP -DBWA_MEM_BUG -DNOCHUNK -DNo_VECTORIZE  -DENABLE_PREFETCH -D_64BIT
#LISA_CPPFLAGS_EXACT_SEARCH_SCALAR = ${HP} -DNO_DNA_ORD -DBWA_MEM_BUG -DNOCHUNK -DNO_VECTORIZE  -DENABLE_PREFETCH -D_64BIT
#LISA_CPPFLAGS_EXACT_SEARCH_UNOPT = -DVECTORIZE  -DENABLE_PREFETCH -DNO_DNA_ORD -DBWA_MEM_BUG -DNOCHUNK -D_64BIT
#LISA_OUTPUT = qbwt-rmi-parallel-make.o
