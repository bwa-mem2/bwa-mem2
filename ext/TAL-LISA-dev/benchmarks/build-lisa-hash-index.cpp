/*************************************************************************************
MIT License

Copyright (c) 2020 Intel Labs

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Authors: Saurabh Kalikar <saurabh.kalikar@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
*****************************************************************************************/

#include<iostream>
#include <algorithm>
#include <math.h>
#include <cmath>
#include <fstream>
//#include <filesystem>
#include <iostream>
#include <stdio.h>
#include <sys/types.h>
#include<sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <x86intrin.h>
#ifdef VTUNE_ANALYSIS
#include <ittnotify.h>
#endif
#include <omp.h>
#include<vector>

#include "lisa_hash.h"
using namespace std;



int main(int argc, char* argv[]){

#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
	char *prefix;
	long num_leaf = 0;
	if(argc < 2){
		cout<<"Invalid number of arguments: 1. key-value file [2. Query file] [3. Number of rmi leaf nodes]\n";
		exit(0);
	}
	
	if(argc == 4)
		num_leaf = atoi(argv[3]);

	//lisa_hash<rmi_key_type, value_type>
	lisa_hash<uint64_t, uint64_t> lh(argv[1], prefix, num_leaf);	

	if(argc == 2) return 0;

	vector<uint64_t> minimizers;
	ifstream inputQuery(argv[2]);
	uint64_t key;
	while(inputQuery>>key){
		minimizers.push_back(key);
	}
	int numq = minimizers.size();
 	

  	int64_t startTick = __rdtsc();
	int *nn;
	uint64_t *pp = lh.get_hash_values_batched(&minimizers[0], numq, nn);

	int64_t endTick = __rdtsc();
	cout<<"Ticks per lookup = "<< (endTick- startTick)/numq<<endl;

#ifdef OUTPUT
	int count = 0;
	for(int i = 0; i < numq; i++){
		for(int j = 0; j< nn[i]; j++){
			cout<<pp[count++]<< " ";
		}
	cout<<endl;

	}
#endif

}
