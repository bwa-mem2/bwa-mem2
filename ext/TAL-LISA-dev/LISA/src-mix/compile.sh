LISA_CC="icpc " 
LISA_CFLAGS="-DSAIS=1 -std=c++17 -xCORE-AVX512 -Ofast -fopenmp -Wall -Wshadow -Wno-char-subscripts"
LISA_CPPFLAGS="-DOUTPUT ${HP} -DNO_DNA_ORD -DREV_COMP -DBWA_MEM_BUG -DNOCHUNK -DVECTORIZE  -DENABLE_PREFETCH -D_64BIT"
LISA_CPPFLAGS_EXACT_SEARCH="${HP} -DNO_DNA_ORD -DBWA_MEM_BUG -DNOCHUNK -DVECTORIZE  -DENABLE_PREFETCH -D_64BIT"
LISA_BUILD_RMI_FLAGS="-DBWA_MEM_BUG -DNOCHUNK -xCORE-AVX512  -DVECTORIZE -DENABLE_PREFETCH -D_64BIT" 
LISA_INCLUDE="-I src-mix/ -I ../../ext/ -I ../../ext/safestringlib/include/ -I ../../src/FMI"
LISA_LDLIBS="-lz -L./../.. -ltal -L ../../ext/safestringlib/ -lsafestring" 


for cpp_file in ./*.cpp
do
	${LISA_CC} ${LISA_CFLAGS} ${LISA_CPPFLAGS} ${LISA_INCLUDE} ${LISA_LDLIBS} -DPRINT_OUTPUT -c ${cpp_file}
done
