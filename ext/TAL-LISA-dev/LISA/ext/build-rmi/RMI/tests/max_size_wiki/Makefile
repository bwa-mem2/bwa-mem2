
result: test
	$(shell ./test > stdout) 
	echo $(.SHELLSTATUS) > result
	cat stdout >> result

rmi.cpp: ../rmi
	RMI_OPTIMIZER_PROFILE=fast ../rmi ../wiki_ts_200M_uint64 rmi cubic,linear --max-size 50000000

test: main.cpp rmi.cpp
	# -lstdc++fs is required for ancient G++s
	g++ -std=c++17 -Wall -O3 -ffast-math -march=native main.cpp rmi.cpp -o test -lstdc++fs

.PHONY: clean
clean:
	rm -rf test result rmi* 
