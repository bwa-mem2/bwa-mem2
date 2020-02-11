#ifndef BWA_ERT_H
#define BWA_ERT_H

#include "kvec.h"
#include "macro.h"
#include "bwa.h"
#include "bwt.h"
#include "bntseq.h"

typedef struct _node_t node_t;

typedef enum {EMPTY, LEAF, UNIFORM, DIVERGE} node_type_t;

struct _node_t {
    node_type_t type;
    int pos;
    int num_bp;
    int l_seq;
    int numChildren;
    uint64_t numHits;
    uint64_t start_addr;
    uint64_t* hits;
    uint8_t seq[READ_LEN + 1];
    node_t* parent_node;
    node_t* child_nodes[4];
};

typedef kvec_t(node_t) node_v;

typedef node_t* node_ptr_t;

typedef struct {
    int tid;
    int step;
    int readLength;
    uint64_t* kmer_table;
    uint64_t startKmer;
    uint64_t endKmer;
    bwaidx_t* bid; 
    uint64_t* numHits;
    char* filePrefix;
    uint64_t* byte_offsets;
    uint64_t end_offset;
} thread_data_t;

// FIXME : Add to options later
const uint8_t char_count_size_in_bits = 8;
const uint8_t hits_count_size_in_bits = 8;
const uint8_t ref_ptr_size_in_bits = 40;
const uint8_t leaf_offset_ptr_size_in_bits = 8;
const uint8_t other_offset_ptr_size_in_bits = 32;

typedef enum { CODE, EMPTY_NODE, LEAF_COUNT, LEAF_HITS, UNIFORM_COUNT, UNIFORM_BP, LEAF_PTR, OTHER_PTR } byte_type_t;

void ert_build_kmertree(const bwt_t* bwt, const bntseq_t* bns, const uint8_t* pac, bwtintv_t ik, bwtintv_t ok[4], int curDepth, node_t* parent_node, int step, int max_depth);

void handleDivergence(const bwt_t* bwt, const bntseq_t* bns, const uint8_t* pac, bwtintv_t ok[4], int depth, node_t* parent_node, int step, int max_depth);

void handleLeaf(const bwt_t* bwt, const bntseq_t* bns, const uint8_t* pac, bwtintv_t ik, node_t* n, int step);

void ert_build_table(const bwt_t* bwt, const bntseq_t* bns, const uint8_t* pac, bwtintv_t ik, bwtintv_t ok[4], uint8_t* mlt_data, uint8_t* mh_data, uint64_t* size, uint64_t* mh_size, uint8_t* aq, uint64_t* numHits, uint64_t* max_next_ptr, uint64_t next_ptr_width, int step, int max_depth);

void ert_traverse_kmertree(node_t* n, uint8_t* mlt_data, uint8_t* mh_data, uint64_t* byte_idx, uint64_t* mh_byte_idx, int depth, uint64_t* numHits, uint64_t* max_ptr, uint64_t next_ptr_width, int step);

void ert_destroy_kmertree(node_t* n);

void buildKmerTrees(char* kmer_tbl_file_name, bwaidx_t* bid, char* prefix, int num_threads, int readLength);

#endif
