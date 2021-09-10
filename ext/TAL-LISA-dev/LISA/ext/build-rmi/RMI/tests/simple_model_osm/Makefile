
result: test
	$(shell ./test > stdout) 
	echo $(.SHELLSTATUS) > result
	cat stdout >> result

rmi.cpp: ../rmi
	../rmi ../osm_cellids_200M_uint64 rmi robust_linear,linear 262144

test: main.cpp rmi.cpp
	# -lstdc++fs is required for ancient G++s
	g++ -std=c++17 -Wall -O3 -ffast-math -march=native main.cpp rmi.cpp -o test -lstdc++fs

.PHONY: clean
clean:
	rm -rf test result rmi* 
