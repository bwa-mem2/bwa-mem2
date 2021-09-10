##/*************************************************************************************
##                           The MIT License
##
##   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
##   Copyright (C) 2019  Intel Corporation, Heng Li.
##
##   Permission is hereby granted, free of charge, to any person obtaining
##   a copy of this software and associated documentation files (the
##   "Software"), to deal in the Software without restriction, including
##   without limitation the rights to use, copy, modify, merge, publish,
##   distribute, sublicense, and/or sell copies of the Software, and to
##   permit persons to whom the Software is furnished to do so, subject to
##   the following conditions:
##
##   The above copyright notice and this permission notice shall be
##   included in all copies or substantial portions of the Software.
##
##   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
##   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
##   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
##   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
##   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
##   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
##   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
##   SOFTWARE.
##
##Contacts: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
##                                Heng Li <hli@jimmy.harvard.edu> 
##*****************************************************************************************/

ifneq ($(portable),)
	STATIC_GCC=-static-libgcc -static-libstdc++
endif

EXE=		bwa-mem2
#CXX=		icpc
ifeq ($(CXX), icpc)
	CC= icc
else ifeq ($(CXX), g++)
	CC=gcc
endif		
ARCH_FLAGS=	-msse -msse2 -msse3 -mssse3 -msse4.1
MEM_FLAGS=	-DSAIS=1 -std=c++17 
CPPFLAGS+=	-DENABLE_PREFETCH -DV17=1 -DMATE_SORT=0 $(MEM_FLAGS)
LISA_SMEM_FLAGS= -DOUTPUT -DNO_DNA_ORD -DREV_COMP -DBWA_MEM_BUG -DNOCHUNK -DVECTORIZE  -DENABLE_PREFETCH -D_64BIT #-DHUGE_PAGE
CPPFLAGS+= $(LISA_SMEM_FLAGS)
EXT=         ext/TAL
EXTEXT=      ext/TAL/ext
EXTBSW=      $(EXT)/src/BSW
EXTFMI=      $(EXT)/src/FMI
EXTLISA=     $(EXT)/LISA
EXTLISASRC=  $(EXTLISA)/src-mix
INCLUDES=   -Isrc -I$(EXT)/ext/safestringlib/include -I$(EXTBSW) -I$(EXTFMI) -I$(EXTEXT) -I$(EXTLISASRC)
LIBS=		-lpthread -lm -lz -L. -lbwa -L$(EXT)/ext/safestringlib -lsafestring $(STATIC_GCC)
OBJS=		src/fastmap.o src/bwtindex.o $(EXTEXT)/utils.o src/memcpy_bwamem.o src/kthread.o \
			src/kstring.o src/ksw.o $(EXTEXT)/bntseq.o src/bwamem.o src/profiling.o $(EXTBSW)/bandedSWA.o \
			$(EXTFMI)/FMI_search.o  src/bwamem_pair.o src/kswv.o src/bwa.o $(EXTEXT)/bseq.o\
			src/bwamem_extra.o src/kopen.o 
LISA_OBJ=	$(EXTLISASRC)/chunkEncode.o $(EXTLISASRC)/common.o $(EXTLISASRC)/qbwt-rmi-batched.o
OBJS+= $(LISA_OBJ)

BWA_LIB=    libbwa.a
SAFE_STR_LIB=    $(EXT)/ext/safestringlib/libsafestring.a

ifeq ($(arch),sse41)
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-msse4.1
	else
		ARCH_FLAGS=-msse -msse2 -msse3 -mssse3 -msse4.1
	endif
else ifeq ($(arch),sse42)
	ifeq ($(CXX), icpc)	
		ARCH_FLAGS=-msse4.2
	else
		ARCH_FLAGS=-msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2
	endif
else ifeq ($(arch),avx)
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-mavx ##-xAVX
	else	
		ARCH_FLAGS=-mavx
	endif
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
# To provide a different architecture flag like -march=core-avx2.
	ARCH_FLAGS=$(arch)
else
myall:multi
endif

CXXFLAGS+=	-g -O3 -fpermissive $(ARCH_FLAGS) #-Wall ##-xSSE2

.PHONY:all clean depend multi
.SUFFIXES:.cpp .o

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(EXE)

multi:
	rm -f src/*.o $(BWA_LIB); cd $(EXT)/ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=sse41    EXE=bwa-mem2.sse41    CXX=$(CXX) all
	rm -f src/*.o $(BWA_LIB); cd $(EXT)/ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=sse42    EXE=bwa-mem2.sse42    CXX=$(CXX) all
	rm -f src/*.o $(BWA_LIB); cd $(EXT)/ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=avx    EXE=bwa-mem2.avx    CXX=$(CXX) all
	rm -f src/*.o $(BWA_LIB); cd $(EXT)/ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=avx2   EXE=bwa-mem2.avx2     CXX=$(CXX) all
	rm -f src/*.o $(BWA_LIB); cd $(EXT)/ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=avx512 EXE=bwa-mem2.avx512bw CXX=$(CXX) all
	$(CXX) -Wall -O3 src/runsimd.cpp -I$(EXT)/ext/safestringlib/include -L$(EXT)/ext/safestringlib/ -lsafestring $(STATIC_GCC) -o bwa-mem2


$(EXE):$(BWA_LIB) $(SAFE_STR_LIB) src/main.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) src/main.o $(BWA_LIB) $(LIBS) -o $@

$(BWA_LIB):$(OBJS)
	ar rcs $(BWA_LIB) $(OBJS)

$(SAFE_STR_LIB):
	cd $(EXT)/ext/safestringlib/ && $(MAKE) clean && $(MAKE) CC=$(CC) directories libsafestring.a

clean:
	rm -fr src/*.o $(BWA_LIB) $(EXE) bwa-mem2.sse41 bwa-mem2.sse42 bwa-mem2.avx bwa-mem2.avx2 bwa-mem2.avx512bw
	cd $(EXT)/ext/safestringlib/ && $(MAKE) clean
	rm -rf $(EXTBSW)/*.o
	rm -rf $(EXTFMI)/*.o
	rm -rf $(EXTEXT)/*.o

depend:
	(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CXXFLAGS) $(CPPFLAGS) -I. -- src/*.cpp)

# DO NOT DELETE

$(EXTBSW)/bandedSWA.o: $(EXTBSW)/bandedSWA.h $(EXTBSW)/bandedSWA.cpp src/macro.h
$(EXTFMI)/FMI_search.o: $(EXTFMI)/FMI_search.h $(EXTEXT)/bntseq.h 
##$(EXTFMI)/FMI_search.o: $(EXTFMI)/utils.h src/macro.h $(EXTFMI)/bwa.h src/bwt.h src/sais.h
$(EXTFMI)/FMI_search.o: $(EXTEXT)/sais.h src/bwa.h $(EXTEXT)/bseq.h
$(EXTEXT)/bntseq.o: $(EXTEXT)/bntseq.h $(EXTEXT)/utils.h src/macro.h src/kseq.h $(EXTEXT)/khash.h
src/bwa.o: $(EXTEXT)/bntseq.h src/bwa.h $(EXTEXT)/bseq.h src/bwt.h src/macro.h src/ksw.h $(EXTEXT)/utils.h
src/bwa.o: src/kstring.h $(EXTEXT)/kvec.h src/kseq.h
src/bwamem.o: src/bwamem.h src/bwt.h $(EXTEXT)/bntseq.h src/bwa.h $(EXTEXT)/bseq.h src/macro.h
src/bwamem.o: src/kthread.h $(EXTBSW)/bandedSWA.h src/kstring.h src/ksw.h
src/bwamem.o: $(EXTEXT)/kvec.h $(EXTEXT)/ksort.h $(EXTEXT)/utils.h src/profiling.h
src/bwamem.o: $(EXTFMI)/FMI_search.h  src/kbtree.h
src/bwamem_extra.o: src/bwa.h $(EXTEXT)/bseq.h $(EXTEXT)/bntseq.h src/bwt.h src/macro.h src/bwamem.h
src/bwamem_extra.o: src/kthread.h $(EXTBSW)/bandedSWA.h src/kstring.h src/ksw.h
src/bwamem_extra.o: $(EXTEXT)/kvec.h $(EXTEXT)/ksort.h $(EXTEXT)/utils.h src/profiling.h
src/bwamem_extra.o: $(EXTFMI)/FMI_search.h 
src/bwamem_pair.o: src/kstring.h src/bwamem.h src/bwt.h $(EXTEXT)/bntseq.h
src/bwamem_pair.o: src/bwa.h $(EXTEXT)/bseq.h src/macro.h src/kthread.h $(EXTBSW)/bandedSWA.h
src/bwamem_pair.o: src/ksw.h $(EXTEXT)/kvec.h $(EXTEXT)/ksort.h $(EXTEXT)/utils.h
src/bwamem_pair.o: src/profiling.h $(EXTFMI)/FMI_search.h 
src/bwamem_pair.o: src/kswv.h
src/bwtindex.o: $(EXTEXT)/bntseq.h src/bwa.h $(EXTEXT)/bseq.h src/bwt.h src/macro.h $(EXTEXT)/utils.h
src/bwtindex.o: $(EXTFMI)/FMI_search.h 
src/fastmap.o: src/fastmap.h src/bwa.h $(EXTEXT)/bseq.h $(EXTEXT)/bntseq.h src/bwt.h src/macro.h
src/fastmap.o: src/bwamem.h src/kthread.h $(EXTBSW)/bandedSWA.h src/kstring.h
src/fastmap.o: src/ksw.h $(EXTEXT)/kvec.h $(EXTEXT)/ksort.h $(EXTEXT)/utils.h src/profiling.h
src/fastmap.o: $(EXTFMI)/FMI_search.h  src/kseq.h
src/kstring.o: src/kstring.h
src/ksw.o: src/ksw.h src/macro.h
src/kswv.o: src/kswv.h src/macro.h src/ksw.h $(EXTBSW)/bandedSWA.h
src/kthread.o: src/kthread.h src/macro.h src/bwamem.h src/bwt.h $(EXTEXT)/bntseq.h
src/kthread.o: src/bwa.h $(EXTEXT)/bseq.h $(EXTBSW)/bandedSWA.h src/kstring.h src/ksw.h $(EXTEXT)/kvec.h
src/kthread.o: $(EXTEXT)/ksort.h $(EXTEXT)/utils.h src/profiling.h $(EXTFMI)/FMI_search.h
#src/kthread.o: $(EXTFMI)/read_index_ele.h
src/main.o: src/main.h src/kstring.h $(EXTEXT)/utils.h src/macro.h $(EXTBSW)/bandedSWA.h
src/main.o: src/profiling.h
src/profiling.o: src/macro.h
## $(EXTFMI)/read_index_ele.o: $(EXTFMI)/read_index_ele.h $(EXTFMI)/utils.h $(EXTFMI)/bntseq.h
##$(EXTFMI)/read_index_ele.o: $(EXTFMI)/macro.h
$(EXTEXT)/utils.o: $(EXTEXT)/utils.h $(EXTEXT)/ksort.h src/kseq.h
src/memcpy_bwamem.o: src/memcpy_bwamem.h
