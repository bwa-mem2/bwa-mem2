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

Authors: Saurabh Kalikar <saurabh.kalikar@intel.com>; Sanchit Misra <sanchit.misra@intel.com>
*****************************************************************************************/
void prepareChunkBatchVectorized_v1(Info* qPool, int qPoolSize, uint64_t* str_enc, int64_t* intv_all, int K){

	int64_t p[21][8];
	
	int64_t i = 0;
	for(i = 0; (i+8) < qPoolSize; i = i + 8){
		for(int64_t j = 0 ; j < 8; j++){
			Info &q = qPool[i+j];
			int start = q.l - K;
			//for(int k = 0; k < 21; k++){
				p[0][j] = q.p[0 + start];
				p[1][j] = q.p[1 + start];
				p[2][j] = q.p[2 + start];
				p[3][j] = q.p[3 + start];
				p[4][j] = q.p[4 + start];

				p[5][j] = q.p[5 + start];
				p[6][j] = q.p[6 + start];
				p[7][j] = q.p[7 + start];
				p[8][j] = q.p[8 + start];
				p[9][j] = q.p[9 + start];


				p[10][j] = q.p[10 + start];
				p[11][j] = q.p[11 + start];
				p[12][j] = q.p[12 + start];
				p[13][j] = q.p[13 + start];
				p[14][j] = q.p[14 + start];


				p[15][j] = q.p[15 + start];
				p[16][j] = q.p[16 + start];
				p[17][j] = q.p[17 + start];
				p[18][j] = q.p[18 + start];
				p[19][j] = q.p[19 + start];
				p[20][j] = q.p[20 + start];
		//	}
	

			intv_all[2 * (i+j)] = q.intv.first;
			intv_all[2 * (i+j) + 1] = q.intv.second;
			const char *ptr = qPool[i + j + 35].p; int offset = qPool[i + j + 35].l -  K;
        		my_prefetch((const char*)(ptr + offset) , _MM_HINT_T0);
		}
		
		__m512i vAns = _mm512_setzero_si512();	

		int64_t* ptr;// = &p[j][0];
		__m512i vP;// =  _mm512_loadu_si512(ptr);

		for(int64_t j = 0 ; j < 21; j++){
			ptr = &p[j][0];
			vAns = _mm512_slli_epi64(vAns, 2);
			vP =  _mm512_loadu_si512(ptr);
			vAns =  _mm512_or_epi64(vAns, vP);
		}
		//Store results of 8 chunks
		_mm512_store_epi64(str_enc + i, vAns);
	}

	//printf("%lld, %lld\n", i, qPoolSize);

	 for( int64_t j = i; j < qPoolSize; j++){
			Info &q = qPool[j];
			uint64_t nxt_ext = 0;
			
			for(int itr = q.l-K; itr != q.l; itr++) {
			    nxt_ext = (nxt_ext<<2) | (q.p[itr]); 
			}
		str_enc[j] = nxt_ext;		
		intv_all[2 * j] = q.intv.first;
		intv_all[2 * j + 1] = q.intv.second;
	 }


}
void prepareChunkBatchVectorized(Info* qPool, int qPoolSize, uint64_t* str_enc, int64_t* intv_all, int K){

	uint64_t offset[24] = {40,38,36,34,32,30,28,26,24,22,20,18,16,14,12,10,8,6,4,2,0,0,0,0};
	uint64_t v_intv[24] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; 
	__m512i vOff1,vOff2,vOff3;
	vOff1 =  _mm512_loadu_si512(offset);
	vOff2 =  _mm512_loadu_si512(offset+8);
	vOff3 =  _mm512_loadu_si512(offset+16);

    for(int64_t j = 0; j < qPoolSize; j++)
    {
	Info &q = qPool[j];
	intv_all[2 * j] = q.intv.first;
	intv_all[2 * j + 1] = q.intv.second;
	
	int i = q.l - K;
//	for(
	int itr = 0;// itr < 21; itr++)
	v_intv[itr++] = q.p[i++];
	v_intv[itr++] = q.p[i++];
	v_intv[itr++] = q.p[i++];
	v_intv[itr++] = q.p[i++];
	v_intv[itr++] = q.p[i++];

	v_intv[itr++] = q.p[i++];
	v_intv[itr++] = q.p[i++];
	v_intv[itr++] = q.p[i++];
	v_intv[itr++] = q.p[i++];
	v_intv[itr++] = q.p[i++];

	v_intv[itr++] = q.p[i++];
	v_intv[itr++] = q.p[i++];
	v_intv[itr++] = q.p[i++];
	v_intv[itr++] = q.p[i++];
	v_intv[itr++] = q.p[i++];

	v_intv[itr++] = q.p[i++];
	v_intv[itr++] = q.p[i++];
	v_intv[itr++] = q.p[i++];
	v_intv[itr++] = q.p[i++];
	v_intv[itr++] = q.p[i++];
	v_intv[itr++] = q.p[i++];

	__m512i vA1, vA2, vA3;
	vA1 =  _mm512_loadu_si512(v_intv);
	vA2 =  _mm512_loadu_si512(v_intv+8);
	vA3 =  _mm512_loadu_si512(v_intv+16);


	 __m512i vShift1 = _mm512_sllv_epi64(vA1, vOff1);
	 __m512i vShift2 = _mm512_sllv_epi64(vA2, vOff2);
	 __m512i vShift3 = _mm512_sllv_epi64(vA3, vOff3);

	 __m512i or1 = _mm512_or_epi64(vShift1, vShift2);
	 __m512i or2 = _mm512_or_epi64(vShift3, or1);

	uint64_t v_nxt_ext = _mm512_reduce_add_epi64(or2);
	
	str_enc[j] = v_nxt_ext;

	const char *p = qPool[j + 40].p; int off_set = qPool[j + 40].l -  K;
        my_prefetch((const char*)(p + off_set) , _MM_HINT_T0);
    }
}

void prepareChunkBatch(Info* qPool, int qPoolSize, uint64_t* str_enc, int64_t* intv_all, int K){

		    for(int64_t j = 0; j < qPoolSize; j++)
		    {
			Info &q = qPool[j];
			uint64_t nxt_ext = 0;
#ifndef NO_DNA_ORD
			
			for(int i = q.l-K; i != q.l; i++) {
			    nxt_ext = (nxt_ext<<2) | dna_ord(q.p[i]); 
			}
#else			
			for(int i = q.l-K; i != q.l; i++) {
			    nxt_ext = (nxt_ext<<2) | (q.p[i]); 
			}
#endif

//Optimization for K=21 with full loop
#if 0
#ifndef NO_DNA_ORD
			
			nxt_ext = (nxt_ext<<2) | dna_ord(q.p[i++]);
			nxt_ext = (nxt_ext<<2) | dna_ord(q.p[i++]);
			nxt_ext = (nxt_ext<<2) | dna_ord(q.p[i++]);
			nxt_ext = (nxt_ext<<2) | dna_ord(q.p[i++]);
			nxt_ext = (nxt_ext<<2) | dna_ord(q.p[i++]);
			
			nxt_ext = (nxt_ext<<2) | dna_ord(q.p[i++]);
			nxt_ext = (nxt_ext<<2) | dna_ord(q.p[i++]);
			nxt_ext = (nxt_ext<<2) | dna_ord(q.p[i++]);
			nxt_ext = (nxt_ext<<2) | dna_ord(q.p[i++]);
			nxt_ext = (nxt_ext<<2) | dna_ord(q.p[i++]);
			
			nxt_ext = (nxt_ext<<2) | dna_ord(q.p[i++]);
			nxt_ext = (nxt_ext<<2) | dna_ord(q.p[i++]);
			nxt_ext = (nxt_ext<<2) | dna_ord(q.p[i++]);
			nxt_ext = (nxt_ext<<2) | dna_ord(q.p[i++]);
			nxt_ext = (nxt_ext<<2) | dna_ord(q.p[i++]);
			
			nxt_ext = (nxt_ext<<2) | dna_ord(q.p[i++]);
			nxt_ext = (nxt_ext<<2) | dna_ord(q.p[i++]);
			nxt_ext = (nxt_ext<<2) | dna_ord(q.p[i++]);
			nxt_ext = (nxt_ext<<2) | dna_ord(q.p[i++]);
			nxt_ext = (nxt_ext<<2) | dna_ord(q.p[i++]);
			nxt_ext = (nxt_ext<<2) | dna_ord(q.p[i++]);
#else
			nxt_ext = (q.p[i++]);
			nxt_ext = (nxt_ext<<2) + (q.p[i++]);
			nxt_ext = (nxt_ext<<2) + (q.p[i++]);
			nxt_ext = (nxt_ext<<2) + (q.p[i++]);
			nxt_ext = (nxt_ext<<2) + (q.p[i++]);
			
			nxt_ext = (nxt_ext<<2) + (q.p[i++]);
			nxt_ext = (nxt_ext<<2) + (q.p[i++]);
			nxt_ext = (nxt_ext<<2) + (q.p[i++]);
			nxt_ext = (nxt_ext<<2) + (q.p[i++]);
			nxt_ext = (nxt_ext<<2) + (q.p[i++]);
			
			nxt_ext = (nxt_ext<<2) + (q.p[i++]);
			nxt_ext = (nxt_ext<<2) + (q.p[i++]);
			nxt_ext = (nxt_ext<<2) + (q.p[i++]);
			nxt_ext = (nxt_ext<<2) + (q.p[i++]);
			nxt_ext = (nxt_ext<<2) + (q.p[i++]);
			
			nxt_ext = (nxt_ext<<2) + (q.p[i++]);
			nxt_ext = (nxt_ext<<2) + (q.p[i++]);
			nxt_ext = (nxt_ext<<2) + (q.p[i++]);
			nxt_ext = (nxt_ext<<2) + (q.p[i++]);
			nxt_ext = (nxt_ext<<2) + (q.p[i++]);
			nxt_ext = (nxt_ext<<2) + (q.p[i++]);

	
#endif

#endif
			str_enc[j] = nxt_ext;

			intv_all[2 * j] = q.intv.first;
			intv_all[2 * j + 1] = q.intv.second;
			const char *p = qPool[j + 40].p; int offset = qPool[j + 40].l -  K;
                        my_prefetch((const char*)(p + offset) , _MM_HINT_T0);
		    }
}
