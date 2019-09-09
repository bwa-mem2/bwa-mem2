##/*************************************************************************************
##                           The MIT License
##
##   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
##   Copyright (C) 2019  Vasimuddin Md, Sanchit Misra, Intel Corporation, Heng Li.
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


EXE=		fmi_test smem2_test bwt_seed_strategy_test sa2ref_test xeonbsw
CXX=		icpc
CXXFLAGS=	-std=c++11 -fopenmp -mtune=native -march=native
CPPFLAGS=	-DENABLE_PREFETCH -DBWA_OTHER_ELE=0
INCLUDES=	-I../src
LIBS=		-fopenmp

.PHONY:all clean depend
.SUFFIXES:.cpp .o

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(EXE)

FMI_search.o:../src/FMI_search.cpp
	$(CC) -c $(CXXFLAGS) $< -o $@

fmi_test:fmi_test.o FMI_search.o
	$(CC) -o $@ $^ $(LIBS)

smem2_test: smem2_test.o FMI_search.o
	$(CC) -o $@ $^ $(LIBS)

sa2ref_test: sa2ref_test.o FMI_search.o
	$(CC) -o $@ $^ $(LIBS)

bwt_seed_strategy_test:bwt_seed_strategy_test.o FMI_search.o
	$(CC) -o $@ $^ $(LIBS)

bandedSWA.o:../src/bandedSWA.cpp ../src/bandedSWA.h
	$(CXX) -c $(CXXFLAGS) $< -o $@

xeonbsw:main_banded.o bandedSWA.o
	$(CXX) -o $@ $^ $(LIBS)

clean:
	rm -fr *.o $(EXE)

depend:
	(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -- *.cpp ../src/FMI_search.cpp ../src/bandedSWA.cpp)

# DO NOT DELETE

bwt_seed_strategy_test.o: ../src/fasta_file.h ../src/FMI_search.h
bwt_seed_strategy_test.o: ../src/bntseq.h ../src/read_index_ele.h
bwt_seed_strategy_test.o: ../src/utils.h ../src/macro.h
fmi_test.o: ../src/fasta_file.h ../src/FMI_search.h ../src/bntseq.h
fmi_test.o: ../src/read_index_ele.h ../src/utils.h ../src/macro.h
main_banded.o: ../src/bandedSWA.h ../src/macro.h
sa2ref_test.o: ../src/fasta_file.h ../src/FMI_search.h ../src/bntseq.h
sa2ref_test.o: ../src/read_index_ele.h ../src/utils.h ../src/macro.h
smem2_test.o: ../src/fasta_file.h ../src/FMI_search.h ../src/bntseq.h
smem2_test.o: ../src/read_index_ele.h ../src/utils.h ../src/macro.h
../src/FMI_search.o: ../src/FMI_search.h ../src/bntseq.h
../src/FMI_search.o: ../src/read_index_ele.h ../src/utils.h ../src/macro.h
../src/bandedSWA.o: ../src/bandedSWA.h ../src/macro.h
