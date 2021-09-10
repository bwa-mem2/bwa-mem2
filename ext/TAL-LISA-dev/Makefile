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
##		Saurabh Kalikar <saurabh.kalikar@intel.com>
##*****************************************************************************************/

BUILD_INDEX_FOR_ONLY= build-index-forward-only
BUILD_INDEX_WITH_RC= build-index-with-rev-complement
BENCH_SMEM=		bench-smem
BENCH_FIXED_LEN_E2E=		bench-fixed-len-e2e-match
#CXX=		icpc
ifeq ($(CXX), icpc)
    CC= icc
else ifeq ($(CXX), g++)
	CC=gcc
endif

ARCH_FLAGS=	-msse4.1
#ARCH_FLAGS=	-mavx512bw
MEM_FLAGS=	-DSAIS=1
CPPFLAGS=	-DENABLE_PREFETCH $(MEM_FLAGS) -DKSW=1
INCLUDES=   -Iext -Iext/safestringlib/include -Isrc/FMI/
LIBS=		-fopenmp -lm -lz -L. -ltal -Lext/safestringlib/ -lsafestring
OBJS=		ext/utils.o \
			ext/kstring.o ext/bntseq.o \
			src/FMI/FMI_search.o ext/bseq.o
TAL_LIB=    libtal.a
SAFE_STR_LIB=    ext/safestringlib/libsafestring.a

ifeq ($(arch),sse)
		ARCH_FLAGS=-msse4.1
else ifeq ($(arch),avx2)
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-march=core-avx2 #-xCORE-AVX2
	else	
		ARCH_FLAGS=-mavx2
	endif
else ifeq ($(arch),avx512)
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-xCORE-AVX512
	else	
		ARCH_FLAGS=-mavx512bw
	endif
else ifeq ($(arch),native)
	ARCH_FLAGS=-march=native
else ifneq ($(arch),)
## To provide a different architecture flag like -march=core-avx2.
	ARCH_FLAGS=$(arch)
endif

#CXXFLAGS=	-g -O3 -fpermissive $(ARCH_FLAGS) #-Wall ##-xSSE2
CXXFLAGS=	-g -O3 -fopenmp $(ARCH_FLAGS)

#.PHONY:all clean depend
.PHONY:all clean depend
#.PHONY:all clean
.SUFFIXES:.cpp .o

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(BUILD_INDEX_FOR_ONLY) $(BUILD_INDEX_WITH_RC) $(BENCH_SMEM) $(BENCH_FIXED_LEN_E2E)

$(BUILD_INDEX_FOR_ONLY):$(TAL_LIB) $(SAFE_STR_LIB) benchmarks/build-index-forward-only.o
	$(CXX) $(CXXFLAGS) benchmarks/build-index-forward-only.o $(LIBS) -o $@

$(BUILD_INDEX_WITH_RC):$(TAL_LIB) $(SAFE_STR_LIB) benchmarks/build-index-with-rev-complement.o
	$(CXX) $(CXXFLAGS) benchmarks/build-index-with-rev-complement.o $(LIBS) -o $@

$(BENCH_SMEM):$(TAL_LIB) $(SAFE_STR_LIB) benchmarks/bench-smem.o
	$(CXX) $(CXXFLAGS) benchmarks/bench-smem.o $(LIBS) -o $@

$(BENCH_FIXED_LEN_E2E):$(TAL_LIB) $(SAFE_STR_LIB) benchmarks/bench-fixed-len-e2e-match.o
	$(CXX) $(CXXFLAGS) benchmarks/bench-fixed-len-e2e-match.o $(LIBS) -o $@

$(TAL_LIB):$(OBJS)
	ar rcs $(TAL_LIB) $(OBJS)

$(SAFE_STR_LIB):
	cd ext/safestringlib/ && make clean && make CC=$(CC) directories libsafestring.a

clean:
	rm -fr src/*.o ext/*.o benchmarks/*.o $(TAL_LIB) $(BUILD_INDEX_FOR_ONLY) $(BUILD_INDEX_WITH_RC) $(BENCH_SMEM) $(BENCH_FIXED_LEN_E2E)
	cd ext/safestringlib && make CC=$(CC) clean

depend:
	(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CXXFLAGS) $(CPPFLAGS) -I. -- src/*.cpp)

KEY_TYPE=-DUINT64
ifeq ($(key_type),f64)
	KEY_TYPE=-DF64
else ifeq ($(key_type),uint64)
	KEY_TYPE=-DUINT64
endif

#LISA hash
lisa_hash: ./benchmarks/build-lisa-hash-index.cpp  ./src/LISA-hash/lisa_hash.h
	$(CXX) ./benchmarks/build-lisa-hash-index.cpp -o build-lisa-hash-index -Ofast -DVECTORIZE ${KEY_TYPE} -I./src/LISA-hash -march=native

#Dynamic programming based chaining benchmark
dp_chain: ./benchmarks/bench-dp-chaining.cpp ./src/dynamic-programming/parallel_chaining_32_bit.h
	$(CXX) ./benchmarks/bench-dp-chaining.cpp -I src/dynamic-programming/ -Ofast -o bench-dp-chaining -march=native -DDEBUG


# DO NOT DELETE

src/FMI/FMI_search.o: src/FMI/FMI_search.h ext/bntseq.h
src/FMI/FMI_search.o: ext/utils.h  ext/sais.h
ext/bntseq.o: ext/bntseq.h ext/utils.h ext/kseq.h
#ext/bwa.o: ext/bntseq.h ext/bwa.h ext/utils.h
#ext/bwa.o: ext/kstring.h ext/kseq.h
ext/kstring.o: ext/kstring.h
ext/utils.o: ext/utils.h ext/kseq.h
ext/bseq.o: ext/bseq.h ext/utils.h ext/kseq.h
benchmarks/bench-smem.o: src/FMI/FMI_search.h ext/bntseq.h ext/utils.h
benchmarks/bench-smem.o:  ext/sais.h
