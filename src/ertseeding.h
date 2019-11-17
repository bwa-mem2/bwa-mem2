#ifndef ERTSEEDING_HPP
#define ERTSEEDING_HPP

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <vector>
#include "kstring.h"
#include "ksw.h"
#include "kvec.h"
#include "ksort.h"
#include "utils.h"
#include "bwt.h"
#include "ertindex.h"
#include "bntseq.h"
#include "bwa.h"
#include "macro.h"
#include "profiling.h"

typedef struct {
    uint8_t c;
    int i;
    int num_hits;
    uint64_t byte_idx;
} node_info_t;

typedef kvec_t(uint64_t) u64v;
typedef kvec_t(uint8_t) u8v;
typedef kvec_t(int) intv;

typedef kvec_t(node_info_t) path_v;

typedef struct {
    int c_pivot;    //!< Pivot used to generate the SMEM
    int p_pivot;    //!< Previous pivot
    int pp_pivot;   //!< Pivot before the previous pivot. Useful in reseeding
} pivot_t;

typedef struct {
    uint8_t forward; //!< RMEM or LMEM. We need this to normalize hit positions
    int start; //!< MEM start position in read
    int end; //!< MEM end position in read. [start, end)
    int rc_start;
    int rc_end;
    int type;
    int skip_ref_fetch; 
    int fetch_leaves;
    int hitbeg;
    int hitcount;
    pivot_t pt;
    // u64v hits; //!< Hits for MEM
} mem_t;

typedef kvec_t(mem_t) mem_v;

/**
 * Index-related auxiliary data structures
 */

typedef struct {
    uint64_t* kmer_offsets;     // K-mer table
    uint8_t* mlt_table;         // Multi-level ERT
    uint8_t* leaf_table;         // Multi-level ERT
    const bwt_t* bwt;           // FM-index
    const bntseq_t* bns;        // Input reads sequences
    const uint8_t* pac;         // Reference genome (2-bit encoded)
} index_aux_t;

/**
 * 'Read' auxiliary data structures
 */

typedef struct {
    int min_seed_len;               // Minimum length of seed
    int l_seq;                      // Read length
    int ptr_width;                  // Size of pointers to child nodes in ERT
    int num_hits;                   // Number of hits for each node in the ERT
    int limit;                      // Number of hits after which extension must be stopped
    uint64_t lep[5];                // FIXME: Can support up to 320bp
    uint64_t nextLEPBit;            // Index into the LEP bit-vector
    uint64_t mlt_start_addr;        // Start address of multi-level ERT
    uint64_t mh_start_addr;         // Start address of multi-hits for each k-mer
    uint64_t lft_start_addr;         // Start address of multi-hits for each k-mer
    char* read_name;                // Read name
    char* queue_buf;                // Read sequence
    char* rc_queue_buf;             // Reverse complemented read sequence
    uint8_t* unpacked_queue_buf;    // Read sequence (2-bit encoded)
    uint8_t* unpacked_rc_queue_buf; // Reverse complemented read (2-bit encoded)
    uint8_t* read_buf;              // == queue_buf (forward) and == rc_queue_buf (backward)
} read_aux_t;

/**
 * SMEM auxiliary data structures
 */

typedef struct {
    int prevMemStart;               // Start position of previous MEM in the read
    int prevMemEnd;                 // End position of previous MEM in the read
    int curr_pivot;                 // Pivot used for forward/backward search in iteration i
    int prev_pivot;                 // Pivot used in the previous iteration (i-1)
    int prev_prev_pivot;            // Pivot used in iteration i-2 (useful for reseeding)
    int stop_be;                    // Stop backward search early if no new SMEMs can be found for pivot
} smem_helper_t;

void get_seeds(index_aux_t* iaux, read_aux_t* raux, mem_t* smems, int64_t& num_smem, u64v* hits);

void reseed(index_aux_t* iaux, read_aux_t* raux, mem_t* smems, int start, int limit, pivot_t* pt, int64_t& num_smem, u64v* hits);

void last(index_aux_t* iaux, read_aux_t* raux, mem_t* smems, int limit, int64_t& num_smem, u64v* hits);

#endif
