#include <sys/types.h>
#include <sys/mman.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <stdio.h>
#include "bwa.h"
#include "fastmap.h"

int bwa_shm_stage(ktp_aux_t *aux, const char *hint)
{
    bwaidx_fm_t *idx = aux->fmi->idx;
    // reference string in aux->ref_string (unit8_t*)

	const char *name;
	uint8_t *shm, *shm_idx;
	uint16_t *cnt;
	int shmid, to_init = 0, l;
	char path[PATH_MAX + 1];

	if (hint == 0 || hint[0] == 0) return -1;
	for (name = hint + strlen(hint) - 1; name >= hint && *name != '/'; --name);
	++name;

	if ((shmid = shm_open("/bwa2ctl", O_RDWR, 0)) < 0) {
		shmid = shm_open("/bwa2ctl", O_CREAT|O_RDWR|O_EXCL, 0644);
		to_init = 1;
	}
	if (shmid < 0) return -1;
	ftruncate(shmid, BWA_CTL_SIZE);
	shm = mmap(0, BWA_CTL_SIZE, PROT_READ|PROT_WRITE, MAP_SHARED, shmid, 0);
	cnt = (uint16_t*)shm;
	if (to_init) {
		memset(shm, 0, BWA_CTL_SIZE);
		cnt[1] = 4;
	}

	// if the index object is not contiguous in memory, re-stage it
	if (aux->fmi->idx->mem == 0) {
	  aux->fmi->bwa_idx2mem(); 
	}
	
	int64_t reflen = *(aux->fmi->reference_seq_len);

	strcat(strcpy(path, "/bwaidx2-"), name);
	if ((shmid = shm_open(path, O_CREAT|O_RDWR|O_EXCL, 0644)) < 0) {
		shm_unlink(path);
		perror("shm_open()");
		return -1;
	}
	l = 16 + strlen(name) + 1;
	if (cnt[1] + l > BWA_CTL_SIZE) return -1;
	memcpy(shm + cnt[1], &idx->l_mem, 8);
	memcpy(shm + cnt[1] + 8, &reflen, 8);
	memcpy(shm + cnt[1] + 16, name, l - 16);
	cnt[1] += l; ++cnt[0];

	ftruncate(shmid, idx->l_mem + reflen);
	shm_idx = mmap(0, idx->l_mem + reflen, PROT_READ|PROT_WRITE, MAP_SHARED, shmid, 0);
	memcpy(shm_idx, idx->mem, idx->l_mem);
	memcpy(shm_idx+ (idx->l_mem), aux->ref_string, reflen);
	free(idx->mem);
	free(aux->ref_string);
	aux->fmi->bwa_mem2idx(idx->l_mem, shm_idx);
	aux->ref_string = shm_idx + idx->l_mem;
	idx->is_shm = 1;
	return 0;
	
}

void bwa_idx_load_from_shm(const char *hint, ktp_aux_t *aux)
{
    fprintf(stderr, "* Load index from shm\n");

	const char *name;
	uint8_t *shm, *shm_idx;
	uint16_t *cnt, i;
	char *p, path[PATH_MAX + 1];
	int shmid;
	int64_t l_mem, reflen;
    bwaidx_fm_t *idx = aux->fmi->idx;

	if (hint == 0 || hint[0] == 0) return 0;
	for (name = hint + strlen(hint) - 1; name >= hint && *name != '/'; --name);
	++name;
	if ((shmid = shm_open("/bwa2ctl", O_RDONLY, 0)) < 0) return 0;
	shm = mmap(0, BWA_CTL_SIZE, PROT_READ, MAP_SHARED, shmid, 0);
	cnt = (uint16_t*)shm;
	if (cnt[0] == 0) return 0;
	for (i = 0, p = (char*)(shm + 4); i < cnt[0]; ++i) {
		memcpy(&l_mem, p, 8); p += 8;
		memcpy(&reflen, p, 8); p += 8;
		if (strcmp(p, name) == 0) break;
		p += strlen(p) + 1;
	}
	if (i == cnt[0]) return 0;

    fprintf(stderr, "* Load index from shm: size %ld, reflen %ld\n",l_mem,reflen);

	strcat(strcpy(path, "/bwaidx2-"), name);
	if ((shmid = shm_open(path, O_RDONLY, 0)) < 0) return 0;
	shm_idx = mmap(0, l_mem + reflen, PROT_READ, MAP_SHARED, shmid, 0);
	aux->fmi->bwa_mem2idx(l_mem, shm_idx);
	aux->ref_string = shm_idx + l_mem;
	idx->is_shm = 1;
}

int bwa_shm_test(const char *hint)
{
	int shmid;
	uint16_t *cnt, i;
	char *p, *shm;
	const char *name;

	if (hint == 0 || hint[0] == 0) return 0;
	for (name = hint + strlen(hint) - 1; name >= hint && *name != '/'; --name);
	++name;
	if ((shmid = shm_open("/bwa2ctl", O_RDONLY, 0)) < 0) return 0;
	shm = mmap(0, BWA_CTL_SIZE, PROT_READ, MAP_SHARED, shmid, 0);
	cnt = (uint16_t*)shm;
	for (i = 0, p = shm + 4; i < cnt[0]; ++i) {
		if (strcmp(p + 16, name) == 0) return 1;
		p += strlen(p) + 17;
	}
	return 0;
}

int bwa_shm_list(void)
{
	int shmid;
	uint16_t *cnt, i;
	char *p, *shm;
	if ((shmid = shm_open("/bwa2ctl", O_RDONLY, 0)) < 0) return -1;
	shm = mmap(0, BWA_CTL_SIZE, PROT_READ, MAP_SHARED, shmid, 0);
	cnt = (uint16_t*)shm;
	for (i = 0, p = shm + 4; i < cnt[0]; ++i) {
		int64_t l_mem, reflen;
		memcpy(&l_mem, p, 8); p += 8;
		memcpy(&reflen, p, 8); p += 8;
		printf("%s\t%ld\t%ld\n", p, (long)l_mem, (long)reflen);
		p += strlen(p) + 1;
	}
	return 0;
}

int bwa_shm_destroy(void)
{
	int shmid;
	uint16_t *cnt, i;
	char *p, *shm;
	char path[PATH_MAX + 1];

	if ((shmid = shm_open("/bwa2ctl", O_RDONLY, 0)) < 0) return -1;
	shm = mmap(0, BWA_CTL_SIZE, PROT_READ, MAP_SHARED, shmid, 0);
	cnt = (uint16_t*)shm;
	for (i = 0, p = shm + 4; i < cnt[0]; ++i) {
		int64_t l_mem, reflen;
		memcpy(&l_mem, p, 8); p += 8;
		memcpy(&reflen, p, 8); p += 8;
		fprintf(stderr, "Unloading shm index: %s\n",p);
		strcat(strcpy(path, "/bwaidx2-"), p);
		shm_unlink(path);
		p += strlen(p) + 1;
	}
	munmap(shm, BWA_CTL_SIZE);
	fprintf(stderr, "Removing shm control file\n");
	shm_unlink("/bwa2ctl");
	return 0;
}


int bwa_shm_remove(const char* hint)
{
	int shmid, shmid2, kept=0;
	uint16_t *cnt, *cnt2, i, l, to_init;
	char *p, *shm, *shm2;
	const char *name;
	char path[PATH_MAX + 1];
	
	if (hint == 0 || hint[0] == 0) return 0;
	for (name = hint + strlen(hint) - 1; name >= hint && *name != '/'; --name);
	++name;
	
	// create a second control file and copy all remaining indexes there
	if ((shmid2 = shm_open("/_bwa2ctl", O_RDWR, 0)) < 0) {
		shmid2 = shm_open("/_bwa2ctl", O_CREAT|O_RDWR|O_EXCL, 0644);
		to_init = 1;
	}
	if (shmid2 < 0) return -1;
	ftruncate(shmid2, BWA_CTL_SIZE);
	shm2 = mmap(0, BWA_CTL_SIZE, PROT_READ|PROT_WRITE, MAP_SHARED, shmid2, 0);
	cnt2 = (uint16_t*)shm2;
	if (to_init) {
		memset(shm2, 0, BWA_CTL_SIZE);
		cnt2[1] = 4;
	}
	
	// open main control file, iterate over items, copy all that are not dropped
	if ((shmid = shm_open("/bwa2ctl", O_RDONLY, 0)) < 0) return -1;
	shm = mmap(0, BWA_CTL_SIZE, PROT_READ, MAP_SHARED, shmid, 0);
	cnt = (uint16_t*)shm;
	for (i = 0, p = shm + 4; i < cnt[0]; ++i) {
		int64_t l_mem, reflen;
		memcpy(&l_mem, p, 8); p += 8;
		memcpy(&reflen, p, 8); p += 8;
		if (strcmp(p, name) == 0) { 
			fprintf(stderr, "Unloading shm index: %s\n",p);
			strcat(strcpy(path, "/bwaidx2-"), p);
			shm_unlink(path);
		} else {
			fprintf(stderr, "Keeping shm index: %s\n",p);
			l = 16 + strlen(p) + 1;
			memcpy(shm2 + cnt2[1], &l_mem, 8);
			memcpy(shm2 + cnt2[1] + 8, &reflen, 8);
			memcpy(shm2 + cnt2[1] + 16, p, l - 16);
			cnt2[1] += l; ++cnt2[0];
			++kept;
		}
		p += strlen(p) + 1;	
	}

	// now copy back
	shm_unlink("/bwa2ctl");
	if (kept>0) {
		fprintf(stderr, "Replacing shm control file\n");
		if ((shmid = shm_open("/bwa2ctl", O_RDWR, 0)) < 0) {
			shmid = shm_open("/bwa2ctl", O_CREAT|O_RDWR|O_EXCL, 0644);
			to_init = 1;
		}
		if (shmid < 0) return -1;
		ftruncate(shmid, BWA_CTL_SIZE);
		shm = mmap(0, BWA_CTL_SIZE, PROT_READ|PROT_WRITE, MAP_SHARED, shmid, 0);
		cnt = (uint16_t*)shm;
		if (to_init) {
			memset(shm, 0, BWA_CTL_SIZE);
			cnt[1] = 4;
		}	
	
		for (i = 0, p = shm2 + 4; i < cnt2[0]; ++i) {
			int64_t l_mem, reflen;
			memcpy(&l_mem, p, 8); p += 8;
			memcpy(&reflen, p, 8); p += 8;
			l = 16 + strlen(p) + 1;
			memcpy(shm + cnt[1], &l_mem, 8);
			memcpy(shm + cnt[1] + 8, &reflen, 8);
			memcpy(shm + cnt[1] + 16, p, l - 16);
			cnt[1] += l; ++cnt[0];
			p += strlen(p) + 1;	
		}
	} else {
		fprintf(stderr, "No more indexes left, removing control file\n");
	}	

	munmap(shm, BWA_CTL_SIZE);
	munmap(shm2, BWA_CTL_SIZE);
	
	shm_unlink("/_bwa2ctl");
	
	return 0;
}

int main_shm(int argc, char *argv[])
{
	int c, to_list = 0, to_drop = 0, ret = 0;
	while ((c = getopt(argc, argv, "ld")) >= 0) {
		if (c == 'l') to_list = 1;
		else if (c == 'd') to_drop = 1;
	}
	++optind;

	if (optind == argc && !to_list && !to_drop) {
		fprintf(stderr, "\nUsage: bwa shm [-d|-l][idxbase]\n\n");
		fprintf(stderr, "Options: -d       destroy all indices in shared memory\n");
		fprintf(stderr, "         -l       list names of indices in shared memory\n");
		return 1;
	}
	if (optind < argc && to_list) {
		fprintf(stderr, "[E::%s] -l cannot be used when 'idxbase' is present\n", __func__);
		return 1;
	}
	if (to_list) {
		bwa_shm_list();
	} else if (!to_drop && optind < argc) {
		if (bwa_shm_test(argv[optind]) == 0) {
			ktp_aux_t aux;
			load_complete_index(argv[optind], &aux);
			if (bwa_shm_stage(&aux, argv[optind]) < 0) {
				fprintf(stderr, "[E::%s] failed to stage the index in shared memory\n", __func__);
				ret = 1;
			}
		} else fprintf(stderr, "[M::%s] index '%s' is already in shared memory\n", __func__, argv[optind]);
	} else if (to_drop) {
		if (optind<argc) { //index name of to-destroy index was supplied
			bwa_shm_remove(argv[optind]);
		} else {
			fprintf(stderr, "Unloading all shm indexes\n");
			bwa_shm_destroy();
		}
	} else {
	  fprintf(stderr, "[M::%s] not sure what you want to do\n");
	}
	return ret;
}
