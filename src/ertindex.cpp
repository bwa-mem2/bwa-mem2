#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <limits.h>
#include <pthread.h>
#include <errno.h>
#include <math.h>
#include <fstream>
#include "utils.h"
#include "ertindex.h"
#include "memcpy_bwamem.h"

#ifdef __cplusplus
extern "C" {
#endif
#include "safe_mem_lib.h"
#include "safe_str_lib.h"
#include <snprintf_s.h>
#ifdef __cplusplus
}
#endif

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#define _set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<(((l)&3)<<1))
#define _set_pac_orig(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))

inline void getNumBranchesForKmer(bwtintv_t ok[4], int* numBranches, uint8_t* uniform_bp) {
	uint8_t i;
	for (i = 0; i < 4; ++i) {
		if (ok[i].x[2] > 0) { *numBranches += 1; *uniform_bp = i; }
	}
}

static inline void kmertoquery(uint64_t x, uint8_t *a, int l) {
	int i;
	for (i = 0; i < l; ++i) {
		a[i] = (uint8_t)((x >> (i << 1)) & 0x3);
	}
}

inline uint64_t addBytesForEntry(byte_type_t type, int count, int numHits) {
	uint64_t numBits = 0;
	switch(type) {
		case CODE:
			numBits = 8;
			break;
		case LEAF_COUNT:
			numBits = (hits_count_size_in_bits);
			break;
		case LEAF_HITS:
			numBits = (ref_ptr_size_in_bits * numHits);
			break;
		case UNIFORM_COUNT: // "Uniform"
			numBits = (char_count_size_in_bits);
			break;
		case UNIFORM_BP:
			numBits = (count << 1);
			break;
		case LEAF_PTR: // "Leaf Offset Pointer"
			numBits = leaf_offset_ptr_size_in_bits;
			break;
		case OTHER_PTR:
			numBits = other_offset_ptr_size_in_bits;
			break;
		case EMPTY_NODE:
			numBits = 0;
			break;
		default :
			break;
	}
	return (numBits % 8) ? (numBits / 8 + 1) : numBits / 8;
}

void addChildNode(node_t* p, node_t* c) {
	assert(p->numChildren <= 3); 
	p->child_nodes[p->numChildren++] = c;
	c->parent_node = p;
}

void handleLeaf(const bwt_t* bwt, const bntseq_t *bns, const uint8_t *pac, bwtintv_t ik, node_t* n, int step) {
	n->type = LEAF;
	n->numHits = ik.x[2];
	n->hits = (uint64_t*) calloc(n->numHits, sizeof(uint64_t));
	assert(n->hits != NULL);
	if (step == 1) {
		uint64_t ref_pos = 0;
		int j = 0;
		for (j = 0; j < n->numHits; ++j) {
			ref_pos = bwt_sa(bwt, ik.x[0]+j);
			n->hits[j] = ref_pos;
		}
	}
}

void handleDivergence(const bwt_t* bwt, const bntseq_t *bns, const uint8_t *pac, 
		bwtintv_t ok[4], int depth, node_t* parent_node, int step, int max_depth) {
	int i;
	bwtintv_t ok_copy[4];
	bwtintv_t ik_new;
	memcpy_bwamem(ok_copy, 4*sizeof(bwtintv_t), ok, 4*sizeof(bwtintv_t), __FILE__, __LINE__);
	for (i = 3; i >= 0; --i) {
		node_t* n = (node_t*) calloc(1, sizeof(node_t));
		assert(n != NULL);
		n->numChildren = 0;
		memset_s(n->child_nodes, 4*sizeof(node_t*), 0);
		n->pos = 0;
		n->num_bp = 0;
		if (ok_copy[i].x[2] == 0) { //!< Empty node
			n->type = EMPTY;
			n->numHits = 0;
			memcpy_bwamem(n->seq, READ_LEN * sizeof(uint8_t), parent_node->seq, (depth-1)*sizeof(uint8_t), __FILE__, __LINE__);
			n->l_seq = depth;
			addChildNode(parent_node, n);
		}
		else if (ok_copy[i].x[2] > 1 && depth != max_depth) {
			ik_new = ok_copy[i]; ik_new.info = depth+1;
			memcpy_bwamem(n->seq, READ_LEN * sizeof(uint8_t), parent_node->seq, parent_node->l_seq*sizeof(uint8_t), __FILE__, __LINE__);
			assert(depth >= 0);
			n->seq[depth] = i;
			n->pos = depth;
			n->num_bp = 1;
			n->l_seq = depth+1;
			n->numHits = ok_copy[i].x[2];
			n->type = DIVERGE;
			addChildNode(parent_node, n);
			ert_build_kmertree(bwt, bns, pac, ik_new, ok, depth+1, n, step, max_depth);
		}
		else {
			memcpy_bwamem(n->seq, READ_LEN * sizeof(uint8_t), parent_node->seq, parent_node->l_seq*sizeof(uint8_t), __FILE__, __LINE__);
			n->seq[depth] = i;
			n->pos = depth;
			n->num_bp = 1;
			n->l_seq = depth+1;
			handleLeaf(bwt, bns, pac, ok_copy[i], n, step);
			addChildNode(parent_node, n);
		}
	}   
}

void ert_build_kmertree(const bwt_t* bwt, const bntseq_t *bns, const uint8_t *pac, 
		bwtintv_t ik, bwtintv_t ok[4], int curDepth, node_t* parent_node, int step, int max_depth) {

	uint8_t uniform_bp = 0;
	int numBranches = 0, depth = curDepth;
	bwt_extend(bwt, &ik, ok, 0); //!< Extend right by 1bp
	/// Check if we need to make a uniform entry
	getNumBranchesForKmer(ok, &numBranches, &uniform_bp);
	bwtintv_t ik_new;
	/// If we find a uniform entry, extend further till we diverge or hit a leaf node
	if (numBranches == 1) {
		uint8_t uniformExtend = 1;
		ik_new = ok[uniform_bp]; ik_new.info = depth+1;
		node_t* n = (node_t*) calloc(1, sizeof(node_t));
		assert(n != NULL);
		n->numChildren = 0;
		memset_s(n->child_nodes, 4*sizeof(node_t*), 0);
		memcpy_bwamem(n->seq, READ_LEN * sizeof(uint8_t), parent_node->seq, parent_node->l_seq*sizeof(uint8_t), __FILE__, __LINE__);
		assert(depth >= 0);
		n->seq[depth] = uniform_bp;
		n->numHits = ok[uniform_bp].x[2];
		n->l_seq = depth + 1;
		n->pos = depth;
		n->num_bp = 1;
		addChildNode(parent_node, n);
		if (depth < max_depth) {
			bwtintv_t ok_init;
			memcpy_bwamem(&ok_init, sizeof(bwtintv_t), &ok[uniform_bp], sizeof(bwtintv_t), __FILE__, __LINE__);
			while (uniformExtend) {
				numBranches = 0; uniform_bp = 0;
				depth += 1;
				bwt_extend(bwt, &ik_new, ok, 0); //!< Extend right by 1bp
				getNumBranchesForKmer(ok, &numBranches, &uniform_bp);
				assert(numBranches != 0);
				if (numBranches == 1) { //<! Uniform
					ik_new = ok[uniform_bp]; ik_new.info = depth+1;
					n->seq[depth] = uniform_bp;
					n->l_seq = depth + 1;
					n->num_bp++;
					/// Hit a multi-hit leaf node
					if (depth == max_depth) {
						uniformExtend = 0;
						handleLeaf(bwt, bns, pac, ok_init, n, step);
					}
				}
				else { //!< Diverge
					uniformExtend = 0;
					n->type = UNIFORM; 
					handleDivergence(bwt, bns, pac, ok, depth, n, step, max_depth);                        
				}
			}
		}
		else { //<! Uniform, depth == max_depth, multi-hit leaf node
			uniformExtend = 0; 
			handleLeaf(bwt, bns, pac, ok[uniform_bp], n, step);
		}
	} //!< End uniform entry
	else { //!< Diverge, empty or leaf, same as above
		handleDivergence(bwt, bns, pac, ok, depth, parent_node, step, max_depth);                        
	} //!< End diverge
}

void ert_build_table(const bwt_t* bwt, const bntseq_t *bns, const uint8_t *pac,
		bwtintv_t ik, bwtintv_t ok[4], uint8_t* mlt_data, uint8_t* mh_data,  
		uint64_t* size, uint64_t* mh_size, uint8_t* aq, 
		uint64_t* numHits, uint64_t* max_next_ptr, uint64_t next_ptr_width,
		int step, int max_depth) {

	uint64_t byte_idx = *size;
	uint64_t mh_byte_idx = *mh_size;
	int i,j;
	uint8_t aq1[xmerSize];
	assert(xmerSize <= 15);
	uint64_t lep1 = 0;
	uint8_t c;
	uint64_t prevHits = ik.x[2];
	bwtintv_t ik_copy = ik;
	uint64_t mlt_byte_idx = byte_idx + (numXmers << 3);
	uint64_t xmer_entry = 0;
	uint16_t xmer_data = 0;
	uint64_t mlt_offset = mlt_byte_idx;
	for (i = 0; i < numXmers; ++i) {
		kmertoquery(i, aq1, xmerSize);
		for (j = 0; j < xmerSize; ++j) {
			c = 3 - aq1[j];
			bwt_extend(bwt, &ik, ok, 0); //!< ok contains the result of BWT extension
			if (ok[c].x[2] != prevHits) { //!< hit set changes
				lep1 |= (1ULL << j);
			}
			/// Extend right till k-mer has zero hits
			if (ok[c].x[2] >= 1) { prevHits = ok[c].x[2]; ik = ok[c]; ik.info = kmerSize + j + 1; }
			else { break; }
		}
		uint64_t num_hits = ok[c].x[2];
		if (ok[c].x[2] == 0) {
			xmer_data = ((lep1 & LEP_MASK) << METADATA_BITWIDTH) | INVALID;
		}
		else if (ok[c].x[2] == 1) {
			xmer_data = ((lep1 & LEP_MASK) << METADATA_BITWIDTH) | (SINGLE_HIT_LEAF);
			if (step == 1) {
				mlt_data[mlt_byte_idx] = 0; //!< Not a multi-hit
			}
			mlt_byte_idx++;
			if (step == 1) {
				uint64_t ref_pos = 0;
				ref_pos = bwt_sa(bwt, ok[c].x[0]);
				uint64_t leaf_data = ref_pos << 1;
				memcpy_bwamem(&mlt_data[mlt_byte_idx], 5 * sizeof(uint8_t), &leaf_data, 5 * sizeof(uint8_t), __FILE__, __LINE__);
			}
			mlt_byte_idx += 5;
			*numHits += 1;
		}
		else {
			xmer_data = ((lep1 & LEP_MASK) << METADATA_BITWIDTH) | (INFREQUENT);
			node_t* n = (node_t*) calloc(1, sizeof(node_t));
			assert(n != NULL);
			n->type = DIVERGE;
			n->pos = 0;
			n->num_bp = 0;
			memcpy_bwamem(n->seq, READ_LEN * sizeof(uint8_t), aq, kmerSize, __FILE__, __LINE__);
			n->l_seq = kmerSize;
			memcpy_bwamem(&n->seq[n->l_seq], xmerSize * sizeof(uint8_t), aq1, xmerSize * sizeof(uint8_t), __FILE__, __LINE__);
			n->l_seq += xmerSize;
			n->parent_node = 0;
			n->numChildren = 0;
			memset_s(n->child_nodes, 4*sizeof(node_t*), 0);
			n->start_addr = mlt_byte_idx;
			ert_build_kmertree(bwt, bns, pac, ik, ok, kmerSize+j, n, step, max_depth);
			ert_traverse_kmertree(n, mlt_data, mh_data, &mlt_byte_idx, &mh_byte_idx, kmerSize+j, numHits, 
					max_next_ptr, next_ptr_width, step);
			ert_destroy_kmertree(n);
		}
		if (num_hits < 20) {
			xmer_entry = (mlt_offset << KMER_DATA_BITWIDTH) | (num_hits << 17) | xmer_data;
		}
		else {
			xmer_entry = (mlt_offset << KMER_DATA_BITWIDTH) | xmer_data;
		}
		uint64_t ptr_width = (next_ptr_width < 4) ? next_ptr_width : 0;
		xmer_entry |= (ptr_width << 22);
		if (step == 1) {
			memcpy_bwamem(&mlt_data[byte_idx], 8 * sizeof(uint8_t), &xmer_entry, 8 * sizeof(uint8_t), __FILE__, __LINE__);
		}
		byte_idx += 8;
		mlt_offset = mlt_byte_idx;
		ik = ik_copy;
		prevHits = ik_copy.x[2];
	}
	*size = mlt_byte_idx;
	*mh_size = mh_byte_idx;
}

void addCode(uint8_t* mlt_data, uint64_t* byte_idx, uint8_t code, int step) {
	if (step == 1) {
		memcpy_bwamem(&mlt_data[*byte_idx], sizeof(uint8_t), &code, sizeof(uint8_t), __FILE__, __LINE__);
	}
	*byte_idx += 1;
}

void addUniformNode(uint8_t* mlt_data, uint64_t* byte_idx, uint8_t count, uint8_t* uniformBases, uint64_t hitCount, int step) {
	int numBytesForBP = addBytesForEntry(UNIFORM_BP, count, 0);
	assert(numBytesForBP < 256);
	if (step == 1) {
		memcpy_bwamem(&mlt_data[*byte_idx], sizeof(uint8_t), &count, sizeof(uint8_t), __FILE__, __LINE__);
	}
	*byte_idx += 1;
	if (step == 1) {
		int j;
		uint8_t packUniformBases[numBytesForBP];
		memset_s(packUniformBases, numBytesForBP * sizeof(uint8_t), 0);
		for (j = 0; j < count; ++j) {
			_set_pac_orig(packUniformBases, j, uniformBases[j]);
		}
		memcpy_bwamem(&mlt_data[*byte_idx], numBytesForBP * sizeof(uint8_t), packUniformBases, numBytesForBP * sizeof(uint8_t), __FILE__, __LINE__); 
	}
	*byte_idx += numBytesForBP;
} 

void addLeafNode(uint8_t* mlt_data, uint64_t* byte_idx, uint64_t ref_pos, int step) {
	if (step == 1) {
		uint64_t leaf_data = (ref_pos << 1);
		memcpy_bwamem(&mlt_data[*byte_idx], 5 * sizeof(uint8_t), &leaf_data, 5 * sizeof(uint8_t), __FILE__, __LINE__);
	}
	*byte_idx += 5;
}

void addMultiHitLeafNode(uint8_t* mlt_data, uint64_t* byte_idx, uint64_t count, uint64_t* hits, int step) {
	uint16_t k = 0;
	for (k = 0; k < count; ++k) {
		if (step == 1) {
			uint64_t leaf_data = (hits[k] << 1) | 1ULL;
			memcpy_bwamem(&mlt_data[*byte_idx], 5 * sizeof(uint8_t), &leaf_data, 5 * sizeof(uint8_t), __FILE__, __LINE__);
		}
		*byte_idx += 5;
	}
} 

void addMultiHitLeafCount(uint8_t* mlt_data, uint64_t* byte_idx, uint64_t count, int step) {
	if (step == 1) {
		memcpy_bwamem(&mlt_data[*byte_idx], 2 * sizeof(uint8_t), &count, 2 * sizeof(uint8_t), __FILE__, __LINE__);
	}
	*byte_idx += 2;
}

void addMultiHitLeafPtr(uint8_t* mlt_data, uint64_t* byte_idx, uint64_t mh_byte_idx, int step) {
	if (step == 1) {
		uint64_t mh_data = (mh_byte_idx << 1) | 1ULL; 
		memcpy_bwamem(&mlt_data[*byte_idx], 5 * sizeof(uint8_t), &mh_data, 5 * sizeof(uint8_t), __FILE__, __LINE__);
	}
	*byte_idx += 5;
} 

void ert_traverse_kmertree(node_t* n, uint8_t* mlt_data, uint8_t* mh_data, uint64_t* size, uint64_t* mh_size, int depth, uint64_t* numHits, uint64_t* max_ptr, uint64_t next_ptr_width, int step) {
	int j = 0;
	int cur_depth = depth;
	uint64_t byte_idx = *size;
	uint64_t mh_byte_idx = *mh_size;
	uint8_t code = 0;
	assert(n->numChildren != 0);
	if (n->numChildren == 1) {
		node_t* child = n->child_nodes[0];
		uint8_t c = child->seq[child->pos];
		if (child->type == LEAF) {
			// 
	  // FIXME: In rare cases, when one of the occurrences of the k-mer is at the end of the reference,
	// # hits for parent node is not equal to the sum of #hits of children nodes, and we trigger the assertion below
  // This should not affect results as long as readLength > kmerSize
  // assert(child->numHits > 1);
			code |= (LEAF << (c << 1));
			addCode(mlt_data, &byte_idx, code, step);
			addMultiHitLeafPtr(mlt_data, &byte_idx, mh_byte_idx, step);
			addMultiHitLeafCount(mh_data, &mh_byte_idx, child->numHits, step);
			addMultiHitLeafNode(mh_data, &mh_byte_idx, child->numHits, child->hits, step);
			*numHits += child->numHits;
		}
		else {
			assert(child->type == UNIFORM);
			code |= (UNIFORM << (c << 1));
			addCode(mlt_data, &byte_idx, code, step);
			addUniformNode(mlt_data, &byte_idx, child->num_bp, &child->seq[child->pos], child->numHits, step);
			ert_traverse_kmertree(child, mlt_data, mh_data, &byte_idx, &mh_byte_idx, cur_depth+child->num_bp, 
					numHits, max_ptr, next_ptr_width, step);
		}
	}
	else {
		uint8_t numEmpty = 0, numLeaves = 0;
		for (j = 0; j < n->numChildren; ++j) {
			node_t* child = n->child_nodes[j];
			uint8_t c = child->seq[child->pos];
			if (child->type == EMPTY) {
				numEmpty++;
			}
			else if (child->type == LEAF) {
				numLeaves++;
				code |= (LEAF << (c << 1));
			}
			else {
				code |= (DIVERGE << (c << 1));
			}
		}
		uint8_t numPointers = ((4 - numEmpty - numLeaves) > 0) ? (4 - numEmpty - numLeaves) : 0;
		uint64_t start_byte_idx = byte_idx;
		addCode(mlt_data, &byte_idx, code, step);
		uint64_t ptr_byte_idx = byte_idx;
		uint64_t ptrToOtherNodes[numPointers + 1]; //!< These point to children. We have one more child than number of pointers
		memset_s(ptrToOtherNodes, (numPointers + 1)*sizeof(uint64_t), 0);
		uint64_t numHitsForChildren[numPointers + 1];
		memset_s(numHitsForChildren, (numPointers + 1)*sizeof(uint64_t), 0);
		uint64_t other_idx = 0;
		if (numPointers > 0) {
			byte_idx += (numPointers*next_ptr_width);
		}
		for (j = 0; j < n->numChildren; ++j) {
			node_t* child = n->child_nodes[j];
			if (child->type == LEAF) {
				if (child->numHits == 1) {
					addLeafNode(mlt_data, &byte_idx, child->hits[0], step);
				}
				else {
					addMultiHitLeafPtr(mlt_data, &byte_idx, mh_byte_idx, step);
					addMultiHitLeafCount(mh_data, &mh_byte_idx, child->numHits, step);
					addMultiHitLeafNode(mh_data, &mh_byte_idx, child->numHits, child->hits, step);
				}
			}
		}
		if (numPointers > 0) {
			ptrToOtherNodes[other_idx] = byte_idx;
		}
		for (j = 0; j < n->numChildren; ++j) {
			node_t* child = n->child_nodes[j];
			assert(child->type != UNIFORM);
			if (child->type == DIVERGE) {
				ert_traverse_kmertree(child, mlt_data, mh_data, &byte_idx, &mh_byte_idx, 
						cur_depth+1, numHits, max_ptr, next_ptr_width, step);
				numHitsForChildren[other_idx] = child->numHits; 
				other_idx++;
				ptrToOtherNodes[other_idx] = byte_idx;
			}
		}
		for (j = 0; j < numPointers; ++j) {
			uint64_t pointerToNextNode = (ptrToOtherNodes[j] - start_byte_idx);
			if (pointerToNextNode > *max_ptr) {
				*max_ptr = pointerToNextNode;
			}
			assert(pointerToNextNode < (1 << 26));
		}
		// Fill up pointers based on size of previous children
		if (step == 1) {
			for (j = 0; j < numPointers; ++j) {
				uint64_t pointerToNextNode = (ptrToOtherNodes[j] - start_byte_idx);
				assert(pointerToNextNode < (1 << 26));
				uint64_t reseed_data = 0;
				if (numHitsForChildren[j] < 20) {
					reseed_data = (pointerToNextNode << 6) | (numHitsForChildren[j]);
				}
				else {
					reseed_data = (pointerToNextNode << 6);
				}
				memcpy_bwamem(&mlt_data[ptr_byte_idx], next_ptr_width * sizeof(uint8_t), &reseed_data, next_ptr_width * sizeof(uint8_t), __FILE__, __LINE__);
				ptr_byte_idx += next_ptr_width;
			}
		}
	}
	*size = byte_idx; 
	*mh_size = mh_byte_idx; 
}

void ert_destroy_kmertree(node_t* n) {
	int j;
	if (n == NULL) {
		return;
	}
	if (n->hits) {
		free(n->hits);
	}
	for (j = 0; j < n->numChildren; ++j) {
		ert_destroy_kmertree(n->child_nodes[j]);
	}
	free(n);
}

//
// This function builds the ERT index. 
// Note on pointers to child nodes: When building the radix tree for each k-mer, 
// we try 3 values for pointers to child nodes, 2,3,4 B and choose the smallest
// one possible.
//
void* buildIndex(void *arg) {

	thread_data_t *data = (thread_data_t *)arg;
	bwtintv_t ik, ok[4];
	uint64_t idx = 0;
	uint8_t aq[kmerSize];
	int i; 
	uint8_t c;
	uint64_t lep, prevHits, numBytesPerKmer, numBytesForMh, ref_pos, total_hits = 0, ptr = 0, max_next_ptr = 0;
	uint64_t next_ptr_width = 0; 
	uint64_t nKmerSmallPtr = 0, nKmerMedPtr = 0, nKmerLargePtr = 0;
	uint16_t kmer_data = 0;

	// File to write the multi-level tree index
	char ml_tbl_file_name[PATH_MAX];
	snprintf_s_si(ml_tbl_file_name, PATH_MAX, "%s.mlt_table_%d", data->filePrefix, data->tid);

	// Log progress
	char log_file_name[PATH_MAX];
	snprintf_s_si(log_file_name, PATH_MAX, "%s.log_%d", data->filePrefix, data->tid);

	FILE *ml_tbl_fd = 0, *log_fd = 0;

	ml_tbl_fd = fopen(ml_tbl_file_name, "wb");
	if (ml_tbl_fd == NULL) {
		printf("\nCan't open file or file doesn't exist. mlt_table errno = %d\n", errno);
		pthread_exit(NULL);
	}

	if (bwa_verbose >= 4) {
		log_fd = fopen(log_file_name, "w");
		if (log_fd == NULL) {
			printf("\nCan't open file or file doesn't exist. log errno = %d\n", errno);
			pthread_exit(NULL);
		} 
		log_file(log_fd, "Start: %lu End: %lu", data->startKmer, data->endKmer);
	}

	//
	// Loop for each k-mer and compute LEP when the hit set changes
	//
	for (idx = data->startKmer; idx < data->endKmer; ++idx) {
		max_next_ptr = 0;
		next_ptr_width = 0;
		c = 0;
		lep = 0;   // k-1-bit LEP
		prevHits = 0;
		numBytesPerKmer = 0;
		numBytesForMh = 0;
		kmertoquery(idx, aq, kmerSize); // represent k-mer as uint8_t*
		assert(aq[0] >= 0 && aq[0] <= 3);
		bwt_set_intv(data->bid->bwt, aq[0], ik); // the initial interval of a single base
		ik.info = 1; 
		prevHits = ik.x[2];

		//
		// Backward search k-mer
		//
		for (i = 1; i < kmerSize; ++i) {
			c = 3 - aq[i]; 
			bwt_extend(data->bid->bwt, &ik, ok, 0); // ok contains the result of BWT extension
			if (ok[c].x[2] != prevHits) { // hit set changes
				lep |= (1ULL << (i-1));
			}
			//
			// Extend left till k-mer has zero hits
			//
			if (ok[c].x[2] >= 1) { prevHits = ok[c].x[2]; ik = ok[c]; ik.info = i + 1; }
			else { break; }
		}

		uint64_t num_hits = ok[c].x[2];
		if (ok[c].x[2] == 0) { // "Empty" - k-mer absent in the reference genome
			kmer_data = ((lep & LEP_MASK) << METADATA_BITWIDTH) | INVALID;
		}
		else if (ok[c].x[2] == 1) { // "Leaf" - k-mer has a single hit in the reference genome
			kmer_data = ((lep & LEP_MASK) << METADATA_BITWIDTH) | (SINGLE_HIT_LEAF);
			numBytesPerKmer = 6;
			uint8_t byte_idx = 0;
			uint8_t mlt_data[numBytesPerKmer];
			if (data->step == 1) { 
				mlt_data[byte_idx] = 0; // Mark that the hit is not a multi-hit
			}
			byte_idx++;
			data->numHits[idx-data->startKmer] += ok[c].x[2];
			if (data->step == 1) {
				//
				// Look up suffix array to identify the hit position
				//
				ref_pos = bwt_sa(data->bid->bwt, ok[c].x[0]);
				uint64_t leaf_data = ref_pos << 1;
				memcpy_bwamem(&mlt_data[byte_idx], 5 * sizeof(uint8_t), &leaf_data, 5 * sizeof(uint8_t), __FILE__, __LINE__);
				fwrite(mlt_data, sizeof(uint8_t), numBytesPerKmer, ml_tbl_fd);
			}
			byte_idx += 5;
		}
		//
		// If the number of hits for the k-mer does not exceed the HIT_THRESHOLD,
		// prefer a radix-tree over a multi-level table as the radix tree for the
		// k-mer is likely to be sparse.
		//
		else if (ok[c].x[2] <= HIT_THRESHOLD) {
			kmer_data = ((lep & LEP_MASK) << METADATA_BITWIDTH) | (INFREQUENT);
			node_t* n = (node_t*) calloc(1, sizeof(node_t));
			assert(n != NULL);
			n->type = DIVERGE;
			n->pos = 0;
			n->num_bp = 0;
			memcpy_bwamem(n->seq, READ_LEN * sizeof(uint8_t), aq, kmerSize * sizeof(uint8_t), __FILE__, __LINE__);
			n->l_seq = kmerSize;
			n->parent_node = 0;
			n->numChildren = 0;
			n->numHits = ok[c].x[2];
			n->child_nodes[0] = n->child_nodes[1] = n->child_nodes[2] = n->child_nodes[3] = 0;
			n->start_addr = 0;
			uint8_t* mlt_data = 0;
			next_ptr_width = 2;
			uint8_t* mh_data = 0;
			uint64_t size = 0;
			if (data->step == 1) {
				if (idx != (numKmers - 1)) {
					size = (data->byte_offsets[idx+1] >> KMER_DATA_BITWIDTH) - (data->byte_offsets[idx] >> KMER_DATA_BITWIDTH); 
					assert(size < (1 << 26));
				} 
				else { // FIXME: This is a hack. We know the size of every k-mer tree except the last-kmer 
					size = 1 << 26;
				}
				next_ptr_width = (((data->byte_offsets[idx] >> 22) & 3) == 0)? 4 : ((data->byte_offsets[idx] >> 22) & 3);  
				mlt_data = (uint8_t*) calloc(size, sizeof(uint8_t));
				assert(mlt_data != NULL);
				mh_data = (uint8_t*) calloc(size, sizeof(uint8_t));
				assert(mh_data != NULL);
			}
			ert_build_kmertree(data->bid->bwt, data->bid->bns, data->bid->pac, ik, ok, i, n, data->step, data->readLength - 1);
			//
			// Reserve space for pointer to start of multi-hit address space
			//
			numBytesPerKmer = 4; 

			// Traverse tree and place data in memory space
			ert_traverse_kmertree(n, mlt_data, mh_data, &numBytesPerKmer, &numBytesForMh, 
					i, &data->numHits[idx-data->startKmer], &max_next_ptr, next_ptr_width, data->step);

			if (data->step == 0 || data->step == 1) {
				if (max_next_ptr >= 1024 && max_next_ptr < 262144) {
					next_ptr_width = 3;
					max_next_ptr = 0;
					numBytesPerKmer = 4;
					numBytesForMh = 0;
					ert_traverse_kmertree(n, mlt_data, mh_data, &numBytesPerKmer, &numBytesForMh, i, 
							&data->numHits[idx-data->startKmer], &max_next_ptr, next_ptr_width, data->step); 
				}
				if (max_next_ptr >= 262144) {
					next_ptr_width = 4;
					max_next_ptr = 0;
					numBytesPerKmer = 4;
					numBytesForMh = 0;
					ert_traverse_kmertree(n, mlt_data, mh_data, &numBytesPerKmer, &numBytesForMh, i, 
							&data->numHits[idx-data->startKmer], &max_next_ptr, next_ptr_width, data->step); 
				}
			}
			ert_destroy_kmertree(n);
			assert(numBytesPerKmer < (1 << 26));
			// assert(numBytesForMh < (1 << 24));
			if (data->step == 1) {
				if (idx != numKmers-1) assert((numBytesPerKmer+numBytesForMh) == size);
				memcpy_bwamem(mlt_data, 4*sizeof(uint8_t), &numBytesPerKmer, 4*sizeof(uint8_t), __FILE__, __LINE__);
				fwrite(mlt_data, sizeof(uint8_t), numBytesPerKmer, ml_tbl_fd);
				free(mlt_data);
				fwrite(mh_data, sizeof(uint8_t), numBytesForMh, ml_tbl_fd);
				free(mh_data);
			}
		}
		//
		// If the number of hits for the k-mer exceeds the HIT_THRESHOLD,
		// prefer a multi-level table to encode the suffixes for the 
		// k-mer
		//
		else {
			kmer_data = ((lep & LEP_MASK) << METADATA_BITWIDTH) | (FREQUENT); 
			uint8_t* mlt_data = 0;
			uint8_t* mh_data = 0;
			next_ptr_width = 2;
			uint64_t size = 0;
			if (data->step == 1) {
				if (idx != (numKmers - 1)) {
					size = (data->byte_offsets[idx+1] >> KMER_DATA_BITWIDTH) - (data->byte_offsets[idx] >> KMER_DATA_BITWIDTH); 
					assert(size < (1 << 26));
				} 
				else { //!< FIXME: Hack. We do not store the size of the last-kmer
					size = 1 << 26;
				}
				next_ptr_width = (((data->byte_offsets[idx] >> 22) & 3) == 0)? 4 : ((data->byte_offsets[idx] >> 22) & 3);  
				mlt_data = (uint8_t*) calloc(size, sizeof(uint8_t));
				assert(mlt_data != NULL);
				mh_data = (uint8_t*) calloc(size, sizeof(uint8_t));
				assert(mh_data != NULL);
			}
			numBytesPerKmer = 4;
			ert_build_table(data->bid->bwt, data->bid->bns, data->bid->pac, ik, ok, mlt_data, mh_data, &numBytesPerKmer,
					&numBytesForMh, aq, &data->numHits[idx-data->startKmer], &max_next_ptr, next_ptr_width, data->step,
					data->readLength - 1);
			if (data->step == 0 || data->step == 1) {
				if (max_next_ptr >= 1024 && max_next_ptr < 262144) {
					next_ptr_width = 3;
					max_next_ptr = 0;
					numBytesPerKmer = 4;
					numBytesForMh = 0;
					ert_build_table(data->bid->bwt, data->bid->bns, data->bid->pac, ik, ok, mlt_data, mh_data, 
							&numBytesPerKmer, &numBytesForMh, aq, &data->numHits[idx-data->startKmer], 
							&max_next_ptr, next_ptr_width, data->step, data->readLength - 1);
				}
				if (max_next_ptr >= 262144) {
					next_ptr_width = 4;
					max_next_ptr = 0;
					numBytesPerKmer = 4;
					numBytesForMh = 0;
					ert_build_table(data->bid->bwt, data->bid->bns, data->bid->pac, ik, ok, mlt_data, mh_data, 
							&numBytesPerKmer, &numBytesForMh, aq, &data->numHits[idx-data->startKmer], 
							&max_next_ptr, next_ptr_width, data->step, data->readLength - 1);
				}
			}
			// 
			// Traverse tree and place data in memory
			//
			assert(numBytesPerKmer < (1 << 26));
			// assert(numBytesForMh < (1 << 24));
			if (data->step == 1) {
				if (idx != numKmers-1) assert((numBytesPerKmer+numBytesForMh) == size);
				memcpy_bwamem(mlt_data, 4*sizeof(uint8_t), &numBytesPerKmer, 4*sizeof(uint8_t), __FILE__, __LINE__);
				fwrite(mlt_data, sizeof(uint8_t), numBytesPerKmer, ml_tbl_fd);
				free(mlt_data);
				fwrite(mh_data, sizeof(uint8_t), numBytesForMh, ml_tbl_fd);
				free(mh_data);
			}
		}
		if (num_hits < 20) {
			data->kmer_table[idx-data->startKmer] = (ptr << KMER_DATA_BITWIDTH) | (num_hits << 17) | kmer_data;
		}
		else {
			data->kmer_table[idx-data->startKmer] = (ptr << KMER_DATA_BITWIDTH) | kmer_data;
		}

		//
		ptr += (numBytesPerKmer + numBytesForMh);

		if (next_ptr_width == 2) {
			nKmerSmallPtr++;
		}
		else if (next_ptr_width == 3) {
			nKmerMedPtr++;
		}
		else if (next_ptr_width == 4) {
			nKmerLargePtr++;
			next_ptr_width = 0;
		}

		//
		data->kmer_table[idx-data->startKmer] |= (next_ptr_width << 22);

		if (bwa_verbose >= 4) {
			if (idx == data->endKmer-1) {
				log_file(log_fd, "TotalSize:%lu\n", ptr);
			}
			if ((idx-data->startKmer) % 10000000 == 0) {
				log_file(log_fd, "%lu,%lu,%lu", idx, numBytesPerKmer, ptr);
			}
		}

		total_hits += data->numHits[idx-data->startKmer];

	}

	//
	data->end_offset = ptr;

	if (bwa_verbose >= 4) {
		log_file(log_fd, "Hits:%lu\n", total_hits);
		log_file(log_fd, "nKmersSmallPtrs:%lu", nKmerSmallPtr);
		log_file(log_fd, "nKmersMedPtrs:%lu", nKmerMedPtr);
		log_file(log_fd, "nKmersLargePtrs:%lu", nKmerLargePtr);
		fclose(log_fd);
	}
	fclose(ml_tbl_fd);
	pthread_exit(NULL);
}

void buildKmerTrees(char* kmer_tbl_file_name, bwaidx_t* bid, char* prefix, int num_threads, int readLength) {

	FILE* kmer_tbl_fd;
	pthread_t thr[num_threads];
	int i, rc;
	thread_data_t thr_data[num_threads];
	uint64_t numKmersThread = (uint64_t)ceil(((double)(numKmers))/num_threads);
	if (bwa_verbose >= 3) {
		fprintf(stderr, "[M::%s] Computing tree sizes for each k-mer\n", __func__);
	}
	// 
	// STEP 1: Create threads. Each thread builds the index for a fraction of the k-mers
	//
	for (i = 0; i < num_threads; ++i) {
		thr_data[i].tid = i;
		thr_data[i].step = 0;
		thr_data[i].readLength = readLength;
		thr_data[i].bid = bid; 
		thr_data[i].startKmer = i*numKmersThread;
		thr_data[i].endKmer = ((i + 1)*numKmersThread > numKmers) ? numKmers : (i + 1)*numKmersThread;  
		thr_data[i].end_offset = 0;
		thr_data[i].filePrefix = prefix;
		uint64_t numKmersToProcess = thr_data[i].endKmer - thr_data[i].startKmer;
		thr_data[i].kmer_table = (uint64_t*) calloc(numKmersToProcess, sizeof(uint64_t));
		assert(thr_data[i].kmer_table != NULL);
		thr_data[i].numHits = (uint64_t*) calloc(numKmersToProcess, sizeof(uint64_t));
		assert(thr_data[i].numHits != NULL);
		if ((rc = pthread_create(&thr[i], NULL, buildIndex, &thr_data[i]))) {
			fprintf(stderr, "[M::%s] error: pthread_create, rc: %d\n", __func__, rc);
			return;
		}
	}
	// 
	// block until all threads complete
	//
	for (i = 0; i < num_threads; ++i) {
		pthread_join(thr[i], NULL);
	}

	//
	// Compute absolute offsets for each kmer's tree from per-thread relative offsets
	//
	uint64_t* kmer_table = (uint64_t*) calloc(numKmers, sizeof(uint64_t));
	assert(kmer_table != NULL);
	int tidx;
	uint64_t kidx;
	uint64_t numProcessed = 0;
	uint64_t offset = 0;
	for (tidx = 0; tidx < num_threads; ++tidx) {
		uint64_t numKmersToProcess = thr_data[tidx].endKmer - thr_data[tidx].startKmer;
		for (kidx = 0; kidx < numKmersToProcess; ++kidx) {
			uint64_t rel_offset = thr_data[tidx].kmer_table[kidx] >> KMER_DATA_BITWIDTH;
			uint16_t kmer_data = thr_data[tidx].kmer_table[kidx] & KMER_DATA_MASK;
			uint64_t ptr_width = (thr_data[tidx].kmer_table[kidx] >> 22) & 3;
			uint64_t reseed_hits = (thr_data[tidx].kmer_table[kidx] >> 17) & 0x1F;
			kmer_table[numProcessed + kidx] =   ((offset + rel_offset) << KMER_DATA_BITWIDTH) 
				| (ptr_width << 22) 
				| (reseed_hits << 17) 
				| (kmer_data);        
		}
		numProcessed += numKmersToProcess;
		offset += thr_data[tidx].end_offset;
		free(thr_data[tidx].kmer_table);
		free(thr_data[tidx].numHits);
	}

	// 
	// STEP 2 : Using estimates of each k-mer's tree size from the previous step, write the index to file
	// 
	uint64_t total_size = offset + (numKmers * 8UL);
	if (bwa_verbose >= 3) {
		fprintf(stderr, "[M::%s] Total size of ERT index = %lu B (Expected). (k-mer,tree) = (%lu,%lu)\n", __func__, total_size, numKmers * 8UL, offset);
	}    

	for (i = 0; i < num_threads; ++i) {
		thr_data[i].tid = i;
		thr_data[i].step = 1;
		thr_data[i].readLength = readLength;
		thr_data[i].bid = bid; 
		thr_data[i].startKmer = i*numKmersThread;
		thr_data[i].endKmer = ((i + 1)*numKmersThread > numKmers) ? numKmers : (i + 1)*numKmersThread;  
		thr_data[i].end_offset = 0;
		thr_data[i].filePrefix = prefix; 
		uint64_t numKmersToProcess = thr_data[i].endKmer - thr_data[i].startKmer;
		thr_data[i].kmer_table = (uint64_t*) calloc(numKmersToProcess, sizeof(uint64_t));
		thr_data[i].numHits = (uint64_t*) calloc(numKmersToProcess, sizeof(uint64_t));
		thr_data[i].byte_offsets = kmer_table;
		if ((rc = pthread_create(&thr[i], NULL, buildIndex, &thr_data[i]))) {
			fprintf(stderr, "[M::%s] error: pthread_create, rc: %d\n", __func__, rc);
			return;
		}
	}

	for (i = 0; i < num_threads; ++i) {
		pthread_join(thr[i], NULL);
	}

	if (bwa_verbose >= 3) {
		fprintf(stderr, "[M::%s] Merging per-thread tables ...\n", __func__);
	}
	//
	// Compute absolute offsets for each k-mer tree's root node
	//  
	numProcessed = 0;
	offset = 0; 
	for (tidx = 0; tidx < num_threads; ++tidx) {
		uint64_t numKmersToProcess = thr_data[tidx].endKmer - thr_data[tidx].startKmer;
		for (kidx = 0; kidx < numKmersToProcess; ++kidx) {
			uint64_t rel_offset = thr_data[tidx].kmer_table[kidx] >> KMER_DATA_BITWIDTH;
			uint16_t kmer_data = thr_data[tidx].kmer_table[kidx] & KMER_DATA_MASK;
			uint64_t ptr_width = (thr_data[tidx].kmer_table[kidx] >> 22) & 3;
			uint64_t reseed_hits = (thr_data[tidx].kmer_table[kidx] >> 17) & 0x1F;
			kmer_table[numProcessed + kidx] =   ((offset + rel_offset) << KMER_DATA_BITWIDTH) 
				| (ptr_width << 22) 
				| (reseed_hits << 17) 
				| (kmer_data);        
		}
		numProcessed += numKmersToProcess;
		offset += thr_data[tidx].end_offset;
		free(thr_data[tidx].kmer_table);
		free(thr_data[tidx].numHits);
	}
	kmer_tbl_fd = fopen(kmer_tbl_file_name, "wb");
	if (kmer_tbl_fd == NULL) {
		fprintf(stderr, "[M::%s] Can't open file or file doesn't exist.\n", __func__);
		exit(1);
	}
	fwrite(kmer_table, sizeof(uint64_t), numKmers, kmer_tbl_fd);
	fclose(kmer_tbl_fd);
	free(kmer_table);

	// 
	// Merge all per-thread trees
	//
	char ml_tbl_file_name[PATH_MAX];
	strcpy_s(ml_tbl_file_name, PATH_MAX, prefix);
	strcat_s(ml_tbl_file_name, PATH_MAX, ".mlt_table");

	if (remove(ml_tbl_file_name) == 0) {
		fprintf(stderr, "[M::%s] Overwriting existing index file (tree)\n", __func__);
	}
	std::ofstream o_mlt(ml_tbl_file_name, std::ios::binary | std::ios::app);
	if (!o_mlt.is_open()) {
		fprintf(stderr, "[M::%s] Can't open output index file for writing.\n", __func__);
		exit(1);
	}
	for (uint64_t tidx = 0; tidx < num_threads; ++tidx) {
		snprintf_s_si(ml_tbl_file_name, PATH_MAX, "%s.mlt_table_%d", prefix, tidx);
		std::ifstream i_mlt(ml_tbl_file_name, std::ios::binary);
		if (!i_mlt.is_open()) {
			fprintf(stderr, "[M::%s] Can't open per-thread index file for thread %d\n", __func__, tidx);
			exit(1);
		}
		o_mlt << i_mlt.rdbuf();
		if (remove(ml_tbl_file_name) != 0) {
			fprintf(stderr, "[M::%s] Can't remove per-thread index file (tree) for thread %d\n", __func__, tidx);
			exit(1);
		}
	}

}
