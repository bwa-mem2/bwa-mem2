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

Authors: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>.
*****************************************************************************************/

#include "read_index_ele.h"

indexEle::indexEle()
{
	idx = (bwaidx_fm_t*) calloc(1, sizeof(bwaidx_fm_t));
}

indexEle::~indexEle()
{
	if (idx == 0) return;
	if (idx->mem == 0)
	{
		if (idx->bns) bns_destroy(idx->bns);
		if (idx->pac) free(idx->pac);
	} else {
		free(idx->bns->anns); free(idx->bns);
		if (!idx->is_shm) free(idx->mem);
	}
	free(idx);
}

void indexEle::bwa_idx_load_ele(const char *hint, int which)
{
	char *prefix;
	//prefix = bwa_idx_infer_prefix(hint);
	//if (prefix == 0) {
	//	printf("[E::%s] fail to locate the index files\n", __func__);
	//	return;
	//}
	int l_hint = strlen(hint);
	prefix = (char *) malloc(l_hint + 3 + 4 + 1);
	strcpy(prefix, hint);

	printf("prefix: %s\n", prefix);
	
	// idx = (bwaidx_fm_t*) calloc(1, sizeof(bwaidx_fm_t));
	if (which & BWA_IDX_BNS) {
		int i, c;
		idx->bns = bns_restore(prefix);
		//if (idx->bns == 0) {
		//	printf("Error: read_index_ele:38, bns is NULL!!\n");
		//	exit(0);
		//}
		for (i = c = 0; i < idx->bns->n_seqs; ++i)
			if (idx->bns->anns[i].is_alt) ++c;
		
		printf("[M::%s] read %d ALT contigs\n", __func__, c);
		
		if (which & BWA_IDX_PAC)
		{
			idx->pac = (uint8_t*) calloc(idx->bns->l_pac/4+1, 1);
			err_fread_noeof(idx->pac, 1, idx->bns->l_pac/4+1, idx->bns->fp_pac); // concatenated 2-bit encoded sequence
			err_fclose(idx->bns->fp_pac);
			idx->bns->fp_pac = 0;
		}
	}
	free(prefix);
}

#include <sys/file.h>
char* indexEle::bwa_idx_infer_prefix(const char *hint)
{
	char *prefix;
	int l_hint;
	FILE *fp;
	l_hint = strlen(hint);
	prefix = (char *) malloc(l_hint + 3 + 4 + 1);
	strcpy(prefix, hint);
	strcpy(prefix + l_hint, ".64.bwt");
	if ((fp = fopen(prefix, "rb")) != 0)
	{
		fclose(fp);
		prefix[l_hint + 3] = 0;
		return prefix;
	} else {
		strcpy(prefix + l_hint, ".bwt");
		if ((fp = fopen(prefix, "rb")) == 0)
		{
			free(prefix);
			return 0;
		} else {
			//flock(fileno(fp), 1);
			//flock(fileno(fp), 1);  // Unlock the file
			fclose(fp);
			prefix[l_hint] = 0;
			return prefix;
		}
	}
}

#if TEST
//int main(int argc, char* argv[])
//{
//	printf("Testing read_index_ele...\n");
//	indexEle *bwaEle = new indexEle();
//	
//	bwaEle->bwa_idx_load_ele("/projects/PCL-GBB/wasim/read_and_ref_data_1/hgaa.fa",
//							BWA_IDX_ALL);
//
//	delete bwaEle;
//}
#endif
