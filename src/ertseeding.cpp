#include "ertseeding.h"

#define smem_lt(a, b) ((a).start == (b).start ? (a).end > (b).end : (a).start < (b).start)
KSORT_INIT(mem_smem_ert, mem_t, smem_lt)

#define smem_lt_2(a, b) ((a).start == (b).start ? (a).end < (b).end : (a).start < (b).start)
KSORT_INIT(mem_smem_sort_lt_ert, mem_t, smem_lt_2)

/**
 * 256 x 3 (2B,3B,4B ptrs) x 4 (ACGT) table encoding offset to next address given a code byte for a leaf node
 */
unsigned char leaf_table[256 * 3][4] = {
        {0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,2,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {2,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,0,2,0},{0,5,0,0},
        {10,5,0,0},{0,5,0,0},{0,7,2,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,0,2,0},{0,0,2,0},
        {7,0,2,0},{0,0,2,0},{0,0,4,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,2,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {2,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {2,0,0,0},{0,0,0,0},{0,0,0,0},{0,2,0,0},
        {7,2,0,0},{0,2,0,0},{0,4,0,0},{0,0,0,0},
        {2,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {4,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,0,0,2},{0,5,0,0},
        {10,5,0,0},{0,5,0,0},{0,7,0,2},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,0,0,2},{0,0,0,2},
        {7,0,0,2},{0,0,0,2},{0,0,0,4},{0,0,5,0},
        {10,0,5,0},{0,0,5,0},{0,0,7,2},{0,10,5,0},
        {15,10,5,0},{0,10,5,0},{0,12,7,2},{0,0,5,0},
        {10,0,5,0},{0,0,5,0},{0,0,7,2},{0,0,7,2},
        {12,0,7,2},{0,0,7,2},{0,0,9,4},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,0,0,2},{0,5,0,0},
        {10,5,0,0},{0,5,0,0},{0,7,0,2},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,0,0,2},{0,0,0,2},
        {7,0,0,2},{0,0,0,2},{0,0,0,4},{0,0,0,2},
        {7,0,0,2},{0,0,0,2},{0,0,0,4},{0,7,0,2},
        {12,7,0,2},{0,7,0,2},{0,9,0,4},{0,0,0,2},
        {7,0,0,2},{0,0,0,2},{0,0,0,4},{0,0,0,4},
        {9,0,0,4},{0,0,0,4},{0,0,0,6},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,2,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {2,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,0,2,0},{0,5,0,0},
        {10,5,0,0},{0,5,0,0},{0,7,2,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,0,2,0},{0,0,2,0},
        {7,0,2,0},{0,0,2,0},{0,0,4,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,2,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {2,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {2,0,0,0},{0,0,0,0},{0,0,0,0},{0,2,0,0},
        {7,2,0,0},{0,2,0,0},{0,4,0,0},{0,0,0,0},
        {2,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {4,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {2,0,0,0},{0,0,0,0},{0,0,0,0},{0,2,0,0},
        {7,2,0,0},{0,2,0,0},{0,4,0,0},{0,0,0,0},
        {2,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {4,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,2,0},
        {7,0,2,0},{0,0,2,0},{0,0,4,0},{0,7,2,0},
        {12,7,2,0},{0,7,2,0},{0,9,4,0},{0,0,2,0},
        {7,0,2,0},{0,0,2,0},{0,0,4,0},{0,0,4,0},
        {9,0,4,0},{0,0,4,0},{0,0,6,0},{0,0,0,0},
        {2,0,0,0},{0,0,0,0},{0,0,0,0},{0,2,0,0},
        {7,2,0,0},{0,2,0,0},{0,4,0,0},{0,0,0,0},
        {2,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {4,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {4,0,0,0},{0,0,0,0},{0,0,0,0},{0,4,0,0},
        {9,4,0,0},{0,4,0,0},{0,6,0,0},{0,0,0,0},
        {4,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {6,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,3,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {3,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,0,3,0},{0,5,0,0},
        {10,5,0,0},{0,5,0,0},{0,8,3,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,0,3,0},{0,0,3,0},
        {8,0,3,0},{0,0,3,0},{0,0,6,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,3,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {3,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {3,0,0,0},{0,0,0,0},{0,0,0,0},{0,3,0,0},
        {8,3,0,0},{0,3,0,0},{0,6,0,0},{0,0,0,0},
        {3,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {6,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,0,0,3},{0,5,0,0},
        {10,5,0,0},{0,5,0,0},{0,8,0,3},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,0,0,3},{0,0,0,3},
        {8,0,0,3},{0,0,0,3},{0,0,0,6},{0,0,5,0},
        {10,0,5,0},{0,0,5,0},{0,0,8,3},{0,10,5,0},
        {15,10,5,0},{0,10,5,0},{0,13,8,3},{0,0,5,0},
        {10,0,5,0},{0,0,5,0},{0,0,8,3},{0,0,8,3},
        {13,0,8,3},{0,0,8,3},{0,0,11,6},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,0,0,3},{0,5,0,0},
        {10,5,0,0},{0,5,0,0},{0,8,0,3},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,0,0,3},{0,0,0,3},
        {8,0,0,3},{0,0,0,3},{0,0,0,6},{0,0,0,3},
        {8,0,0,3},{0,0,0,3},{0,0,0,6},{0,8,0,3},
        {13,8,0,3},{0,8,0,3},{0,11,0,6},{0,0,0,3},
        {8,0,0,3},{0,0,0,3},{0,0,0,6},{0,0,0,6},
        {11,0,0,6},{0,0,0,6},{0,0,0,9},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,3,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {3,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,0,3,0},{0,5,0,0},
        {10,5,0,0},{0,5,0,0},{0,8,3,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,0,3,0},{0,0,3,0},
        {8,0,3,0},{0,0,3,0},{0,0,6,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,3,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {3,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {3,0,0,0},{0,0,0,0},{0,0,0,0},{0,3,0,0},
        {8,3,0,0},{0,3,0,0},{0,6,0,0},{0,0,0,0},
        {3,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {6,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {3,0,0,0},{0,0,0,0},{0,0,0,0},{0,3,0,0},
        {8,3,0,0},{0,3,0,0},{0,6,0,0},{0,0,0,0},
        {3,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {6,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,3,0},
        {8,0,3,0},{0,0,3,0},{0,0,6,0},{0,8,3,0},
        {13,8,3,0},{0,8,3,0},{0,11,6,0},{0,0,3,0},
        {8,0,3,0},{0,0,3,0},{0,0,6,0},{0,0,6,0},
        {11,0,6,0},{0,0,6,0},{0,0,9,0},{0,0,0,0},
        {3,0,0,0},{0,0,0,0},{0,0,0,0},{0,3,0,0},
        {8,3,0,0},{0,3,0,0},{0,6,0,0},{0,0,0,0},
        {3,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {6,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {6,0,0,0},{0,0,0,0},{0,0,0,0},{0,6,0,0},
        {11,6,0,0},{0,6,0,0},{0,9,0,0},{0,0,0,0},
        {6,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {9,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,4,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {4,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,0,4,0},{0,5,0,0},
        {10,5,0,0},{0,5,0,0},{0,9,4,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,0,4,0},{0,0,4,0},
        {9,0,4,0},{0,0,4,0},{0,0,8,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,4,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {4,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {4,0,0,0},{0,0,0,0},{0,0,0,0},{0,4,0,0},
        {9,4,0,0},{0,4,0,0},{0,8,0,0},{0,0,0,0},
        {4,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {8,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,0,0,4},{0,5,0,0},
        {10,5,0,0},{0,5,0,0},{0,9,0,4},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,0,0,4},{0,0,0,4},
        {9,0,0,4},{0,0,0,4},{0,0,0,8},{0,0,5,0},
        {10,0,5,0},{0,0,5,0},{0,0,9,4},{0,10,5,0},
        {15,10,5,0},{0,10,5,0},{0,14,9,4},{0,0,5,0},
        {10,0,5,0},{0,0,5,0},{0,0,9,4},{0,0,9,4},
        {14,0,9,4},{0,0,9,4},{0,0,13,8},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,0,0,4},{0,5,0,0},
        {10,5,0,0},{0,5,0,0},{0,9,0,4},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,0,0,4},{0,0,0,4},
        {9,0,0,4},{0,0,0,4},{0,0,0,8},{0,0,0,4},
        {9,0,0,4},{0,0,0,4},{0,0,0,8},{0,9,0,4},
        {14,9,0,4},{0,9,0,4},{0,13,0,8},{0,0,0,4},
        {9,0,0,4},{0,0,0,4},{0,0,0,8},{0,0,0,8},
        {13,0,0,8},{0,0,0,8},{0,0,0,12},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,4,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {4,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,0,4,0},{0,5,0,0},
        {10,5,0,0},{0,5,0,0},{0,9,4,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,0,4,0},{0,0,4,0},
        {9,0,4,0},{0,0,4,0},{0,0,8,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {5,0,0,0},{0,0,0,0},{0,4,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {4,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {4,0,0,0},{0,0,0,0},{0,0,0,0},{0,4,0,0},
        {9,4,0,0},{0,4,0,0},{0,8,0,0},{0,0,0,0},
        {4,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {8,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {4,0,0,0},{0,0,0,0},{0,0,0,0},{0,4,0,0},
        {9,4,0,0},{0,4,0,0},{0,8,0,0},{0,0,0,0},
        {4,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {8,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,4,0},
        {9,0,4,0},{0,0,4,0},{0,0,8,0},{0,9,4,0},
        {14,9,4,0},{0,9,4,0},{0,13,8,0},{0,0,4,0},
        {9,0,4,0},{0,0,4,0},{0,0,8,0},{0,0,8,0},
        {13,0,8,0},{0,0,8,0},{0,0,12,0},{0,0,0,0},
        {4,0,0,0},{0,0,0,0},{0,0,0,0},{0,4,0,0},
        {9,4,0,0},{0,4,0,0},{0,8,0,0},{0,0,0,0},
        {4,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {8,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {8,0,0,0},{0,0,0,0},{0,0,0,0},{0,8,0,0},
        {13,8,0,0},{0,8,0,0},{0,12,0,0},{0,0,0,0},
        {8,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {12,0,0,0},{0,0,0,0},{0,0,0,0}
};

/**
 * 256 x 3 (2B,3B,4B ptrs) x 4 (ACGT) table encoding offset to next address given a code byte for an internal 'DIVERGE' node
 */
unsigned char code_table[256 * 3][4] = {
        {0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,2,0,0},
        {0,2,0,0},{0,2,0,0},{4,2,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,2,0,0},
        {0,2,0,0},{0,2,0,0},{4,2,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,2,0,0},
        {0,2,0,0},{0,2,0,0},{4,2,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,2,0,0},
        {0,2,0,0},{0,2,0,0},{4,2,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,2,0,0},
        {0,2,0,0},{0,2,0,0},{4,2,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{2,0,0,0},{0,2,0,0},
        {0,2,0,0},{0,2,0,0},{4,2,0,0},{0,0,2,0},
        {0,0,2,0},{0,0,2,0},{4,0,2,0},{0,0,2,0},
        {0,0,2,0},{0,0,2,0},{4,0,2,0},{0,0,2,0},
        {0,0,2,0},{0,0,2,0},{4,0,2,0},{0,4,2,0},
        {0,4,2,0},{0,4,2,0},{6,4,2,0},
        {0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,3,0,0},
        {0,3,0,0},{0,3,0,0},{6,3,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,3,0,0},
        {0,3,0,0},{0,3,0,0},{6,3,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,3,0,0},
        {0,3,0,0},{0,3,0,0},{6,3,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,3,0,0},
        {0,3,0,0},{0,3,0,0},{6,3,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,3,0,0},
        {0,3,0,0},{0,3,0,0},{6,3,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{3,0,0,0},{0,3,0,0},
        {0,3,0,0},{0,3,0,0},{6,3,0,0},{0,0,3,0},
        {0,0,3,0},{0,0,3,0},{6,0,3,0},{0,0,3,0},
        {0,0,3,0},{0,0,3,0},{6,0,3,0},{0,0,3,0},
        {0,0,3,0},{0,0,3,0},{6,0,3,0},{0,6,3,0},
        {0,6,3,0},{0,6,3,0},{9,6,3,0},
        {0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,4,0,0},
        {0,4,0,0},{0,4,0,0},{8,4,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,4,0,0},
        {0,4,0,0},{0,4,0,0},{8,4,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,4,0,0},
        {0,4,0,0},{0,4,0,0},{8,4,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,4,0,0},
        {0,4,0,0},{0,4,0,0},{8,4,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,4,0,0},
        {0,4,0,0},{0,4,0,0},{8,4,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,0,0,0},
        {0,0,0,0},{0,0,0,0},{4,0,0,0},{0,4,0,0},
        {0,4,0,0},{0,4,0,0},{8,4,0,0},{0,0,4,0},
        {0,0,4,0},{0,0,4,0},{8,0,4,0},{0,0,4,0},
        {0,0,4,0},{0,0,4,0},{8,0,4,0},{0,0,4,0},
        {0,0,4,0},{0,0,4,0},{8,0,4,0},{0,8,4,0},
        {0,8,4,0},{0,8,4,0},{12,8,4,0} 
};

/**
 * Return integer key from k-mer string. 
 *
 * @param str           Query sequence
 * @param keysize       K-mer length
 * @param index         Index into query sequence
 * @param seq_len       Length of query sequence
 * @param end_flag      Flag to indicate end of query
 * @param idx_first_N   Index of first ambiguous bp in the k-mer
 *
 * @return              Integer key (hash) for k-mer
 */
inline uint32_t getHashKey (const uint8_t *str, const int keysize, int index, int seq_len, int* end_flag, int* idx_first_N) {
    int i;
    uint32_t key = 0;
    int len = keysize;
    if (index + keysize > seq_len) {
        *end_flag = 1;
        len = seq_len - index;
    }
    for (i = 0; i < len; ++i) {
        if (str[i] != 4) {
            key |= ((uint32_t)(str[i]) << (i << 1));
        }
        else {
            *idx_first_N = i;
            break;
        }
    }
    return key;
}

/**
 * Compute offset to child and return address of child node
 *
 * @param raux              read parameters
 * @param mlt_data          radix tree of k-mer
 * @param code              type of child nodes (EMPTY:00, LEAF:01, UNIFORM:10, DIVERGE:11) - 2b encoded
 * @param c                 next character in read (used to traverse tree)
 * @param byteIdx           byte index into the mlt_data radix tree
 * 
 */
void getOffsetToChildNode(read_aux_t* raux, uint8_t* mlt_data, uint8_t code, uint8_t c, uint64_t* byteIdx) {
    uint64_t nextByteIdx = *byteIdx;
    uint64_t startByteIdx = nextByteIdx - 1;
    int offset = code_table[((raux->ptr_width - 2) << 8) + code][c];
    uint32_t reseed_data = 0;
    nextByteIdx += offset;
    memcpy(&reseed_data, &mlt_data[nextByteIdx], raux->ptr_width);
    uint32_t jumpByteIdx = reseed_data >> 6;
    raux->num_hits = (reseed_data & 0x3F);
    nextByteIdx = startByteIdx + jumpByteIdx;
    *byteIdx = nextByteIdx;
}


/**
 * Compute offset to start of leaf data
 *
 * @param raux              read parameters
 * @param code              type of child nodes (EMPTY:00, LEAF:01, UNIFORM:10, DIVERGE:11) - 2b encoded
 * @param c                 next character in read (used to traverse tree)
 * 
 */
inline int getOffsetToLeafData(read_aux_t* raux, uint8_t code, uint8_t c) {
    return leaf_table[((raux->ptr_width - 2) << 8) + code][c];
}

/**
 * This routine does depth-first tree traversal starting from an internal node to obtain all hits at leaf nodes
 * 
 * @param raux              read parameters
 * @param mlt_data          radix tree of k-mer
 * @param byteIdx           byte index into the mlt_data radix tree
 * @param mem               maximal-exact-match storage 
 * @param bc                next character in read (used to traverse tree)
 * @param hits              list of hits for read
 */
void getNextByteIdx_dfs(read_aux_t* raux, uint8_t* mlt_data, uint64_t* byte_idx, mem_t* mem, uint8_t bc, u64v* hits) {

    uint64_t nextByteIdx = *byte_idx;
    uint64_t ref_pos = 0;
    uint8_t c;
    c = 3 - bc;
    mem->skip_ref_fetch = 1; // MEM cannot be extended by leaf decompression
    uint8_t code = mlt_data[nextByteIdx++];
    uint8_t code_c = (code >> (c << 1)) & 3;
    assert(code != 0);
    if (code_c == LEAF) { // Hit a leaf node
        int k;
        uint64_t leaf_data = 0;
        nextByteIdx += getOffsetToLeafData(raux, code, c);
        memcpy(&leaf_data, &mlt_data[nextByteIdx], 5);
        if (leaf_data & 1) { // Found a multi-hit leaf node
            nextByteIdx = raux->mh_start_addr + (leaf_data >> 1);
            memcpy(&raux->num_hits, &mlt_data[nextByteIdx], 2);
            nextByteIdx += 2;
            mem->hitcount += raux->num_hits;
            for (k = 0; k < raux->num_hits; ++k) {
                memcpy(&ref_pos, &mlt_data[nextByteIdx], 5);
                kv_push(uint64_t, *hits, ref_pos >> 1);
                nextByteIdx += 5;
                ref_pos = 0;
            }
        }
        else { // Single-hit leaf node
            raux->num_hits = 1;
            mem->hitcount += raux->num_hits;
            kv_push(uint64_t, *hits, leaf_data >> 1);
        }
    }
    else if (code_c == UNIFORM) { // Multi-character internal node
        int countBP = mlt_data[nextByteIdx++];
        int numBitsForBP = countBP << 1;
        int numBytesForBP = (numBitsForBP % 8) ? (numBitsForBP / 8 + 1) : (numBitsForBP / 8);
        uint8_t packedBP[numBytesForBP];
        memcpy(packedBP, &mlt_data[nextByteIdx], numBytesForBP);
        nextByteIdx += numBytesForBP;
        uint8_t k;
        uint64_t startByteIdx = nextByteIdx;
        for (k = 0; k < 4; ++k) {
            getNextByteIdx_dfs(raux, mlt_data, &startByteIdx, mem, k, hits);
            startByteIdx = nextByteIdx; 
        }
    }
    else if (code_c == DIVERGE) { // Single-character internal node
        getOffsetToChildNode(raux, mlt_data, code, c, &nextByteIdx);
        uint8_t k;
        uint64_t startByteIdx = nextByteIdx;
        for (k = 0; k < 4; ++k) {
            getNextByteIdx_dfs(raux, mlt_data, &startByteIdx, mem, k, hits);
            startByteIdx = nextByteIdx;
        }
    }
    *byte_idx = nextByteIdx;    
}

/**
 * Gather all hits for the MEM 
 * 
 * @param raux              read parameters
 * @param mlt_data          radix tree of k-mer
 * @param byteIdx           byte index into the mlt_data radix tree
 * @param mem               MEM including hits 
 * @param hits              list of hits for read
 */
void leaf_gather(read_aux_t* raux, uint8_t* mlt_data, uint64_t* byte_idx, mem_t* mem, u64v* hits) {
    uint8_t k;
    uint64_t startByteIdx = *byte_idx; 
    uint64_t tmpByteIdx = startByteIdx; 
    for (k = 0; k < 4; ++k) {
        getNextByteIdx_dfs(raux, mlt_data, &tmpByteIdx, mem, k, hits);
        tmpByteIdx = startByteIdx;
    }
}   

/**
 * Traverse tree during backward search(seeding). Similar to getNextByteIdx, except that we don't compute LEP
 * 
 * @param raux              read parameters
 * @param mlt_data          radix tree of k-mer
 * @param byteIdx           byte index into the mlt_data radix tree
 * @param i                 index into read buffer
 * @param mem               maximal-exact-match storage  
 * @param hits              list of hits for read
 */
void getNextByteIdx_backward(read_aux_t* raux, uint8_t* mlt_data, uint64_t* byte_idx, int* i, mem_t* mem, u64v* hits) {
    uint64_t nextByteIdx = *byte_idx;
    uint64_t ref_pos = 0;
    uint8_t c, code, code_c;
    if (raux->read_buf[*i] != 4) {
        c = 3 - raux->read_buf[*i];
        code = mlt_data[nextByteIdx++];
        code_c = (code >> (c << 1)) & 3;
        assert(code != 0);
    }
    else { // Terminate MEM search when we hit an 'N'
        code_c = EMPTY;
    }
    if (code_c == EMPTY) { // Gather leaves later during forward traversal
        mem->rc_end = *i;
        mem->fetch_leaves = 1; 
    }
    else if (code_c == LEAF) { // Hit a leaf node
        *i += 1;
        mem->rc_end = *i;
        int k;
        uint64_t leaf_data = 0;
        nextByteIdx += getOffsetToLeafData(raux, code, c);
        memcpy(&leaf_data, &mlt_data[nextByteIdx], 5);
        if (leaf_data & 1) { // Found a multi-hit leaf node
            nextByteIdx = raux->mh_start_addr + (leaf_data >> 1);
            memcpy(&raux->num_hits, &mlt_data[nextByteIdx], 2);
            nextByteIdx += 2;
            mem->hitcount += raux->num_hits;
            for (k = 0; k < raux->num_hits; ++k) {
                memcpy(&ref_pos, &mlt_data[nextByteIdx], 5);
                nextByteIdx += 5;
                kv_push(uint64_t, *hits, ref_pos >> 1);
                ref_pos = 0;
            }
            //
            // We found a multi-hit. But to report hits in the same order as BWA-MEM,
            // we will fetch the hits again in a forward tree traversal
            //
            mem->fetch_leaves = 1;
        } 
        else { // Single-hit leaf node
            mem->hitcount += 1;
            kv_push(uint64_t, *hits, leaf_data >> 1);
        }
    }
    else if (code_c == UNIFORM) { // Multi-character internal node
        uint32_t j;
        int countBP = mlt_data[nextByteIdx++];
        int numBitsForBP = countBP << 1;
        int numBytesForBP = (numBitsForBP % 8) ? (numBitsForBP / 8 + 1) : (numBitsForBP / 8);
        uint8_t packedBP[numBytesForBP];
        memcpy(packedBP, &mlt_data[nextByteIdx], numBytesForBP);
        nextByteIdx += numBytesForBP;
        // Unpack base pairs
        uint8_t unpackedBP[countBP];
        for (j = 0; j < countBP; ++j) {
            unpackedBP[j] = ((packedBP[j >> 2] >> ((~(j) & 3) << 1)) & 3); 
        }
        // Count number of matching base pairs with read
        for (j = 0; j < countBP; ++j) {
            if ((*i + j) >= raux->l_seq) {
                break;
            }
            if (raux->read_buf[*i+j] == 4) {
                break;
            }
            if (3 - raux->read_buf[*i+j] != unpackedBP[j]) {
                break;
            }
        }  
        *i += j;
        if (j == countBP) { // We match all bps
            if (*i < raux->l_seq) {
                getNextByteIdx_backward(raux, mlt_data, &nextByteIdx, i, mem, hits);
            }
            else {
                mem->rc_end = *i;
            }
        }
        else { // Did not match all bps. Gather leaves for MEM in a later forward traversal 
            mem->rc_end = *i;
            mem->fetch_leaves = 1;
        }
    }
    else { // Single-character internal node
        getOffsetToChildNode(raux, mlt_data, code, c, &nextByteIdx);
        *i += 1;
        if (*i < raux->l_seq) {
            getNextByteIdx_backward(raux, mlt_data, &nextByteIdx, i, mem, hits);
        }
        else {
            mem->rc_end = *i;
        }
    }
    *byte_idx = nextByteIdx;    
}

/**
 * Traverse tree during backward search (reseeding). Terminate when fewer than 'raux->limit' hits are found for 
 * an internal node  
 * 
 * @param raux              read parameters
 * @param mlt_data          radix tree of k-mer
 * @param byteIdx           byte index into the mlt_data radix tree
 * @param i                 index into read buffer
 * @param mem               maximal-exact-match storage 
 * @param hits              list of hits for each read 
 */
void getNextByteIdx_backward_wlimit(read_aux_t* raux, uint8_t* mlt_data, uint64_t* byte_idx, int* i, mem_t* mem, u64v* hits) {

    uint64_t nextByteIdx = *byte_idx;
    uint64_t ref_pos = 0;
    uint8_t c, code, code_c;
    if (raux->read_buf[*i] != 4) {
        c = 3 - raux->read_buf[*i];
        code = mlt_data[nextByteIdx++];
        code_c = (code >> (c << 1)) & 3;
        assert(code != 0);
    }
    else { // Terminate MEM search when we hit an 'N'
        code_c = EMPTY;
    }
    if (code_c == EMPTY) { // Gather leaves later during forward traversal
        mem->rc_end = *i;
        mem->fetch_leaves = 1;
    }
    else if (code_c == LEAF) { // Hit a leaf node
        int k;
        uint64_t leaf_data = 0;
        nextByteIdx += getOffsetToLeafData(raux, code, c);
        memcpy(&leaf_data, &mlt_data[nextByteIdx], 5);
        if (leaf_data & 1) { // Found a multi-hit leaf node
            nextByteIdx = raux->mh_start_addr + (leaf_data >> 1);
            memcpy(&raux->num_hits, &mlt_data[nextByteIdx], 2);
            nextByteIdx += 2;
            // Hits exceed reseeding threshold
            if (raux->num_hits >= raux->limit) {
                mem->hitcount += raux->num_hits;
                for (k = 0; k < raux->num_hits; ++k) {
                    memcpy(&ref_pos, &mlt_data[nextByteIdx], 5);
                    nextByteIdx += 5;
                    kv_push(uint64_t, *hits, ref_pos >> 1);
                    ref_pos = 0;
                }
                *i += 1;
            } 
        }
        //
        // We found a multi-hit. But to report hits in the same order as BWA-MEM,
        // we will fetch the hits again in a forward tree traversal
        //
        mem->fetch_leaves = 1;
        mem->rc_end = *i;
    }
    else if (code_c == UNIFORM) { // Multi-character internal node
        uint32_t j;
        int countBP = mlt_data[nextByteIdx++];
        int numBitsForBP = countBP << 1;
        int numBytesForBP = (numBitsForBP % 8) ? (numBitsForBP / 8 + 1) : (numBitsForBP / 8);
        uint8_t packedBP[numBytesForBP];
        memcpy(packedBP, &mlt_data[nextByteIdx], numBytesForBP);
        nextByteIdx += numBytesForBP;
        // Unpack base pairs
        uint8_t unpackedBP[countBP];
        for (j = 0; j < countBP; ++j) {
            unpackedBP[j] = ((packedBP[j >> 2] >> ((~(j) & 3) << 1)) & 3); 
        }
        // Count number of matching base pairs with read
        for (j = 0; j < countBP; ++j) {
            if ((*i + j) >= raux->l_seq) {
                break;
            }
            if (raux->read_buf[*i+j] == 4) {
                break;
            }
            if (3 - raux->read_buf[*i+j] != unpackedBP[j]) {
                break;
            }
        }  
        *i += j;
        if (j == countBP) { // We match all bps
            if (*i < raux->l_seq) {
                getNextByteIdx_backward_wlimit(raux, mlt_data, &nextByteIdx, i, mem, hits);
            }
            else {
                mem->rc_end = *i;
                mem->fetch_leaves = 1;
            }
        }
        else { // Did not match all bps. Gather leaves for MEM in a later forward traversal 
            mem->rc_end = *i;
            mem->fetch_leaves = 1;
        }
    }
    else { // Single-character internal node
        raux->num_hits = 0;
        getOffsetToChildNode(raux, mlt_data, code, c, &nextByteIdx);
        // In the internal nodes, raux->num_hits = 0 is used to reprsent # hits > 20
        if ((raux->num_hits == 0) || (raux->num_hits >= raux->limit)) {
            *i += 1;
            if (*i < raux->l_seq) {
                getNextByteIdx_backward_wlimit(raux, mlt_data, &nextByteIdx, i, mem, hits);
            }
            else {
                mem->rc_end = *i;
                mem->fetch_leaves = 1;
            }
        }
        else {
            mem->rc_end = *i;
            mem->fetch_leaves = 1;
        }
    }
    *byte_idx = nextByteIdx;    
}

/**
 * Traverse tree during forward search (seeding). 
 * 
 * @param raux              read parameters
 * @param mlt_data          radix tree of k-mer
 * @param byteIdx           byte index into the mlt_data radix tree
 * @param i                 index into read buffer
 * @param mem               maximal-exact-match storage
 * @param hits              list of hits for read 
 */
void getNextByteIdx(read_aux_t* raux, uint8_t* mlt_data, uint64_t* byte_idx, int* i, mem_t* mem, u64v* hits) {

    uint64_t nextByteIdx = *byte_idx;
    uint64_t parent_byte_idx = nextByteIdx;
    uint64_t ref_pos = 0;
    uint8_t c, code, code_c;
    if (raux->read_buf[*i] != 4) {
        c = 3 - raux->read_buf[*i];
        code = mlt_data[nextByteIdx++];
        code_c = (code >> (c << 1)) & 3;
        assert(code != 0);
    }
    else { // Terminate MEM search when we hit an 'N'
        code_c = EMPTY;
    }
    uint64_t lep_idx = 0;
    uint64_t lep_bit_idx = 0;
    if (code_c == EMPTY) {
        if (mem->start == 0) { // FIXME: Gather leaves later during forward traversal even for MEMs starting at read_pos = 0
            int mem_len = *i;
            if (mem_len >= raux->min_seed_len) { // Only gather leaves when MEM length exceeds threshold
                uint64_t startByteIdx = parent_byte_idx;
                leaf_gather(raux, mlt_data, &startByteIdx, mem, hits); 
            }
        }
        // Update LEP for backward search
        lep_idx = raux->nextLEPBit >> 6;
        lep_bit_idx = raux->nextLEPBit & (0x3FULL);
        raux->lep[lep_idx] |= (1ULL << lep_bit_idx);
        raux->nextLEPBit += 1;
    }
    else if (code_c == LEAF) { // Hit a leaf node
        int k;
        uint64_t leaf_data = 0;
        nextByteIdx += getOffsetToLeafData(raux, code, c);
        memcpy(&leaf_data, &mlt_data[nextByteIdx], 5);
        if (leaf_data & 1) { // Found a multi-hit leaf node
            nextByteIdx = raux->mh_start_addr + (leaf_data >> 1);
            memcpy(&raux->num_hits, &mlt_data[nextByteIdx], 2);
            nextByteIdx += 2;
            mem->hitcount += raux->num_hits;
            for (k = 0; k < raux->num_hits; ++k) {
                memcpy(&ref_pos, &mlt_data[nextByteIdx], 5);
                kv_push(uint64_t, *hits, ref_pos >> 1);
                nextByteIdx += 5;
                ref_pos = 0;
            }
        }
        else { // Single-hit leaf node
            raux->num_hits = 1;
            mem->hitcount += raux->num_hits;
            kv_push(uint64_t, *hits, leaf_data >> 1);
        }
        // Update LEP for backward search
        lep_idx = raux->nextLEPBit >> 6;
        lep_bit_idx = raux->nextLEPBit & (0x3FULL);
        raux->lep[lep_idx] |= (1ULL << lep_bit_idx);
        raux->nextLEPBit += 1;
        *i += 1;
    }
    else if (code_c == UNIFORM) { // Multi-character internal node
        uint32_t j;
        int countBP = mlt_data[nextByteIdx++];
        int numBitsForBP = countBP << 1;
        int numBytesForBP = (numBitsForBP % 8) ? (numBitsForBP / 8 + 1) : (numBitsForBP / 8);
        uint8_t packedBP[numBytesForBP];
        memcpy(packedBP, &mlt_data[nextByteIdx], numBytesForBP);
        nextByteIdx += numBytesForBP;
        // Unpack base pairs
        uint8_t unpackedBP[countBP];
        for (j = 0; j < countBP; ++j) {
            unpackedBP[j] = ((packedBP[j >> 2] >> ((~(j) & 3) << 1)) & 3); 
        }
        // Count number of matching base pairs with read
        for (j = 0; j < countBP; ++j) {
            if ((*i + j) >= raux->l_seq) { // Don't run past the end of the read
                break;
            }
            if (raux->read_buf[*i+j] == 4) {
                break;
            } 
            if ((3 - raux->read_buf[*i+j]) != unpackedBP[j]) {
                break;
            }
        }  
        raux->nextLEPBit += j;
        *i += j;
        if (j == countBP) { // If we match all bases of uniform entry
            if (*i == raux->l_seq) { // Check if we reached the end of the read
                if (mem->start == 0) {
                    leaf_gather(raux, mlt_data, &nextByteIdx, mem, hits); 
                }
                lep_idx = raux->nextLEPBit >> 6;
                lep_bit_idx = raux->nextLEPBit & (0x3FULL);
                raux->lep[lep_idx] |= (1ULL << lep_bit_idx);
            }
            else {
                if (*i < raux->l_seq) {
                    getNextByteIdx(raux, mlt_data, &nextByteIdx, i, mem, hits);
                }
            }
        }
        else { 
            // 
            // We did not match all bases of uniform entry
            // Fetch all hits from leaf nodes for backward extension (dfs :( )
            //
            assert(*i <= raux->l_seq);
            if (mem->start == 0) {
                int mem_len = *i;
                if (mem_len >= raux->min_seed_len) {
                    leaf_gather(raux, mlt_data, &nextByteIdx, mem, hits); 
                } 
            }
            // Update LEP to start backward search from last matching bp
            lep_idx = raux->nextLEPBit >> 6;
            lep_bit_idx = raux->nextLEPBit & (0x3FULL);
            raux->lep[lep_idx] |= (1ULL << lep_bit_idx);
        }
    }
    else { // Single-character internal node
        getOffsetToChildNode(raux, mlt_data, code, c, &nextByteIdx);
        lep_idx = raux->nextLEPBit >> 6;
        lep_bit_idx = raux->nextLEPBit & (0x3FULL);
        raux->lep[lep_idx] |= (1ULL << lep_bit_idx);
        raux->nextLEPBit += 1;
        *i += 1;
        if (*i < raux->l_seq) {
            getNextByteIdx(raux, mlt_data, &nextByteIdx, i, mem, hits);
        }
        else {
            if (mem->start == 0) {
                leaf_gather(raux, mlt_data, &nextByteIdx, mem, hits); 
            }
            raux->lep[raux->nextLEPBit >> 6] |= (1ULL << (raux->nextLEPBit & (0x3FULL)));
            raux->nextLEPBit += 1;
        }
    }
    *byte_idx = nextByteIdx;    
}

/**
 * Traverse tree during forward search (reseeding). Terminate when fewer than 'raux->limit' hits are found for an 
 * internal node    
 * 
 * @param raux              read parameters
 * @param mlt_data          radix tree of k-mer
 * @param byteIdx           byte index into the mlt_data radix tree
 * @param i                 index into read buffer
 * @param mem               maximal-exact-match storage 
 * @param visited           stack to store list of visited nodes
 * @param hits              list of hits for read
 */
void getNextByteIdx_wlimit(read_aux_t* raux, uint8_t* mlt_data, uint64_t* byte_idx, int* i, mem_t* mem, path_v* visited, u64v* hits) {

    uint64_t nextByteIdx = *byte_idx;
    uint64_t parent_byte_idx = nextByteIdx;
    uint64_t ref_pos = 0;
    uint8_t c, code, code_c;
    if (raux->read_buf[*i] != 4) {
        c = 3 - raux->read_buf[*i];
        code = mlt_data[nextByteIdx++];
        code_c = (code >> (c << 1)) & 3;
        assert(code != 0);
    }
    else { // Terminate MEM search when we hit an 'N'
        code_c = EMPTY;
    }
    uint64_t lep_idx = 0;
    uint64_t lep_bit_idx = 0;
    if (code_c == EMPTY) {
        if (mem->start == 0) {
            int mem_len = *i;
            if (mem_len >= raux->min_seed_len) {
                leaf_gather(raux, mlt_data, &parent_byte_idx, mem, hits); 
            }
        }
        // Update LEP for backward search
        lep_idx = raux->nextLEPBit >> 6;
        lep_bit_idx = raux->nextLEPBit & (0x3FULL);
        raux->lep[lep_idx] |= (1ULL << lep_bit_idx);
        raux->nextLEPBit += 1;
    }
    else if (code_c == LEAF) { // Hit a leaf node
        int k;
        uint64_t leaf_data = 0;
        nextByteIdx += getOffsetToLeafData(raux, code, c);
        memcpy(&leaf_data, &mlt_data[nextByteIdx], 5);
        if (leaf_data & 1) { // Found a multi-hit leaf node
            nextByteIdx = raux->mh_start_addr + (leaf_data >> 1);
            memcpy(&raux->num_hits, &mlt_data[nextByteIdx], 2);
            nextByteIdx += 2;
        }
        else { // Single-hit leaf node
            raux->num_hits = 1;
        }
        // Hits exceed reseeding threshold
        if (raux->num_hits >= raux->limit) {
            mem->hitcount += raux->num_hits;
            for (k = 0; k < raux->num_hits; ++k) {
                memcpy(&ref_pos, &mlt_data[nextByteIdx], 5);
                kv_push(uint64_t, *hits, ref_pos >> 1);
                nextByteIdx += 5;
                ref_pos = 0;
            }
            *i += 1;
        }
        else {
            // Do DFS traversal to gather leaves starting from parent node
            if (mem->start == 0) {
                int mem_len = *i;
                if (mem_len >= raux->min_seed_len) { 
                    node_info_t* tmp_node_info = &kv_pop(*visited);
                    leaf_gather(raux, mlt_data, &tmp_node_info->byte_idx, mem, hits);
                }
            }
        }
        // Update LEP for backward search
        lep_idx = raux->nextLEPBit >> 6;
        lep_bit_idx = raux->nextLEPBit & (0x3FULL);
        raux->lep[lep_idx] |= (1ULL << lep_bit_idx);
        raux->nextLEPBit += 1;
    }
    else if (code_c == UNIFORM) { // Multi-character internal node
        uint32_t j;
        int countBP = mlt_data[nextByteIdx++];
        int numBitsForBP = countBP << 1;
        int numBytesForBP = (numBitsForBP % 8) ? (numBitsForBP / 8 + 1) : (numBitsForBP / 8);
        uint8_t packedBP[numBytesForBP];
        memcpy(packedBP, &mlt_data[nextByteIdx], numBytesForBP);
        nextByteIdx += numBytesForBP;
        // Unpack base pairs
        uint8_t unpackedBP[countBP];
        for (j = 0; j < countBP; ++j) {
            unpackedBP[j] = ((packedBP[j >> 2] >> ((~(j) & 3) << 1)) & 3); 
        }
        // Count number of matching base pairs with read
        for (j = 0; j < countBP; ++j) {
            if ((*i + j) >= raux->l_seq) { // Don't run past the end of the read
                break;
            }
            if (raux->read_buf[*i+j] == 4) {
                break;
            } 
            if ((3 - raux->read_buf[*i+j]) != unpackedBP[j]) {
                break;
            }
        }  
        raux->nextLEPBit += j;
        *i += j;
        if (j == countBP) { // If we match all bases of uniform entry
            if (*i == raux->l_seq) { // Check if we reached the end of the read
                if (mem->start == 0) {
                    leaf_gather(raux, mlt_data, &nextByteIdx, mem, hits); 
                }
                lep_idx = raux->nextLEPBit >> 6;
                lep_bit_idx = raux->nextLEPBit & (0x3FULL);
                raux->lep[lep_idx] |= (1ULL << lep_bit_idx);
            }
            else {
                if (*i < raux->l_seq) {
                    getNextByteIdx_wlimit(raux, mlt_data, &nextByteIdx, i, mem, visited, hits);
                }
            }
        }
        else { 
            // 
            // We did not match all bases of uniform entry
            // Fetch all hits from leaf nodes for backward extension (dfs :( )
            // 
            assert(*i <= raux->l_seq);
            if (mem->start == 0) {
                int mem_len = *i;
                if (mem_len >= raux->min_seed_len) {
                    leaf_gather(raux, mlt_data, &nextByteIdx, mem, hits);
                }
            }
            lep_idx = raux->nextLEPBit >> 6;
            lep_bit_idx = raux->nextLEPBit & (0x3FULL);
            raux->lep[lep_idx] |= (1ULL << lep_bit_idx);
        }
    }
    else {
        getOffsetToChildNode(raux, mlt_data, code, c, &nextByteIdx);
        lep_idx = raux->nextLEPBit >> 6;
        lep_bit_idx = raux->nextLEPBit & (0x3FULL);
        raux->lep[lep_idx] |= (1ULL << lep_bit_idx);
        raux->nextLEPBit += 1;
        // In the internal nodes, raux->num_hits = 0 is used to reprsent # hits > 20
        if ((raux->num_hits == 0) || (raux->num_hits >= raux->limit)) {
            node_info_t nif;
            nif.byte_idx = nextByteIdx; 
            nif.num_hits = raux->num_hits;
            kv_push(node_info_t, *visited, nif);
            *i += 1;
            if (*i < raux->l_seq) {
                getNextByteIdx_wlimit(raux, mlt_data, &nextByteIdx, i, mem, visited, hits);
            }
            else {
                if (mem->start == 0) {
                    leaf_gather(raux, mlt_data, &nextByteIdx, mem, hits); 
                }
                raux->lep[raux->nextLEPBit >> 6] |= (1ULL << (raux->nextLEPBit & (0x3FULL)));
                raux->nextLEPBit += 1;
            }
        }
        else {
            // Do DFS traveral to gather leaves from parent node
            if (mem->start == 0) {
                int mem_len = *i;
                if (mem_len >= raux->min_seed_len) {
                    node_info_t* tmp_node_info = &kv_pop(*visited);
                    leaf_gather(raux, mlt_data, &tmp_node_info->byte_idx, mem, hits); 
                }
            }
        }
    }
    *byte_idx = nextByteIdx;    
}

/**
 * Traverse tree during forward search (LAST). Terminate when fewer than 'raux->limit' hits are found for an internal node
 * Minimum MEM length for LAST is opt->min_seed_len + 1    
 * 
 * @param raux              read parameters
 * @param mlt_data          radix tree of k-mer
 * @param byteIdx           byte index into the mlt_data radix tree
 * @param i                 index into read buffer
 * @param mem               maximal-exact match found 
 * @param visited           stack to store list of visited nodes
 * @param hits              list of hits for read
 */
void getNextByteIdx_last(read_aux_t* raux, uint8_t* mlt_data, uint64_t* byte_idx, int* i, mem_t* mem, u64v* hits) {

    uint64_t nextByteIdx = *byte_idx;
    uint64_t ref_pos = 0;
    uint8_t c, code, code_c;
    if (raux->read_buf[*i] != 4) {
        c = 3 - raux->read_buf[*i];
        code = mlt_data[nextByteIdx++];
        code_c = (code >> (c << 1)) & 3;
        assert(code != 0);
    }
    else { // Terminate MEM search when we hit an 'N'
        code_c = EMPTY;
    }
    if (code_c == EMPTY) {
        *i += 1;
    }
    else if (code_c == LEAF) { // Hit a leaf node
        int k;
        uint64_t leaf_data = 0;
        nextByteIdx += getOffsetToLeafData(raux, code, c);
        memcpy(&leaf_data, &mlt_data[nextByteIdx], 5);
        if (leaf_data & 1) { // multi-hit leaf
            nextByteIdx = raux->mh_start_addr + (leaf_data >> 1);
            memcpy(&raux->num_hits, &mlt_data[nextByteIdx], 2);
            nextByteIdx += 2;
            mem->hitcount += raux->num_hits;
            for (k = 0; k < raux->num_hits; ++k) {
                memcpy(&ref_pos, &mlt_data[nextByteIdx], 5);
                kv_push(uint64_t, *hits, ref_pos >> 1);
                nextByteIdx += 5;
                ref_pos = 0;
            }
        }
        else { // single-hit leaf
            raux->num_hits = 1;
            mem->hitcount += raux->num_hits;
            kv_push(uint64_t, *hits, leaf_data >> 1);
        }
        *i += 1;
    }
    else if (code_c == UNIFORM) { // Multi-character internal node
        uint32_t j;
        int countBP = mlt_data[nextByteIdx++];
        int numBitsForBP = countBP << 1;
        int numBytesForBP = (numBitsForBP % 8) ? (numBitsForBP / 8 + 1) : (numBitsForBP / 8);
        uint8_t packedBP[numBytesForBP];
        memcpy(packedBP, &mlt_data[nextByteIdx], numBytesForBP);
        nextByteIdx += numBytesForBP;
        // Unpack base pairs
        uint8_t unpackedBP[countBP];
        for (j = 0; j < countBP; ++j) {
            unpackedBP[j] = ((packedBP[j >> 2] >> ((~(j) & 3) << 1)) & 3); 
        }
        // Count number of matching base pairs with read
        for (j = 0; j < countBP; ++j) {
            if ((*i + j) >= raux->l_seq) { //!< Don't run past the end of the read
                break;
            } 
            if (raux->read_buf[*i+j] == 4) {
                break;
            }
            if ((3 - raux->read_buf[*i+j]) != unpackedBP[j]) {
                break;
            }
        }  
        *i += j;
        int len = *i - mem->start;
        // 
        // LAST stop criterion: MEM is sufficiently long and not too frequent
        //
        int stop_extension = (raux->num_hits > 0 && 
                              raux->num_hits < raux->limit && 
                              len >= (raux->min_seed_len + 1)) ? 1 : 0;
        if (stop_extension) { 
            leaf_gather(raux, mlt_data, &nextByteIdx, mem, hits);
            *i = mem->start + (raux->min_seed_len + 1);
        }
        else { 
            if (j == countBP) {
                if (*i < raux->l_seq) {
                    getNextByteIdx_last(raux, mlt_data, &nextByteIdx, i, mem, hits);
                }
            }
            else {
                *i += 1;
            }
        }
    }
    else if (code_c == DIVERGE) {
        getOffsetToChildNode(raux, mlt_data, code, c, &nextByteIdx);
        *i += 1;
        int len = *i - mem->start;
        // 
        // LAST stop criterion: MEM is sufficiently long and not too frequent
        //
        int stop_extension = (raux->num_hits > 0 && 
                              raux->num_hits < raux->limit && 
                              len >= (raux->min_seed_len + 1)) ? 1 : 0;
        if (stop_extension) {
            leaf_gather(raux, mlt_data, &nextByteIdx, mem, hits);
        }
        else {
            if (*i < raux->l_seq) {
                getNextByteIdx_last(raux, mlt_data, &nextByteIdx, i, mem, hits);
            }
        }
    }
    *byte_idx = nextByteIdx;    
}

/**
 * Main backward search function (seeding). Lookup up k-mer and/or x-mer table and identify root of ERT 
 * 
 * @param iaux              index parameters
 * @param raux              read parameters
 * @param i                 index into read buffer
 * @param mem               maximal-exact match found
 * @param hits              list of hits for read
 */
void leftExtend(index_aux_t* iaux, read_aux_t* raux, int* i, mem_t* mem, u64v* hits) {
    
    uint8_t code;
    uint8_t* mlt_data;
    uint64_t byte_idx = 0, ref_pos = 0, kmer_entry = 0, start_addr = 0;
    uint32_t hashval = 0;
    int idx_first_N = -1;
    hashval = getHashKey(&raux->read_buf[*i], kmerSize, *i, raux->l_seq, 0, &idx_first_N);
    if (idx_first_N != -1) {
        *i += (kmerSize + xmerSize);
        mem->rc_end = *i;
        return;
    }
    // index-table lookup
    kmer_entry = iaux->kmer_offsets[hashval];
    // index-table entry type
    code = kmer_entry & METADATA_MASK;
    // pointer to root of tree
    start_addr = kmer_entry >> KMER_DATA_BITWIDTH;
    // width used for internal pointers in tree
    raux->ptr_width = (((kmer_entry >> 22) & 3) == 0) ? 4 : ((kmer_entry >> 22) & 3);
    byte_idx = 0;
    if (code == INVALID) { // k-mer absent
        *i += (kmerSize + xmerSize);
        mem->rc_end = *i;
    }
    else if (code == SINGLE_HIT_LEAF) { // single-hit k-mer
        mlt_data = &iaux->mlt_table[start_addr];
        byte_idx++;
        memcpy(&ref_pos, &mlt_data[byte_idx], 5);                         
        byte_idx += 5;
        mem->hitcount += 1;
        kv_push(uint64_t, *hits, ref_pos >> 1);
        *i += kmerSize;
        mem->rc_end = *i;
    }
    else if (code == INFREQUENT) { // k-mer has fewer than 256 hits
        *i += kmerSize;
        mlt_data = &iaux->mlt_table[start_addr];
        if (*i < raux->l_seq) {
            // 
            // First 4 bytes of tree store the start of all multi-hit leaves for the k-mer
            //
            memcpy(&raux->mh_start_addr, &mlt_data[byte_idx], 4);
            byte_idx += 4;
            getNextByteIdx_backward(raux, mlt_data, &byte_idx, i, mem, hits);
        }
        else {
            mem->rc_end = *i;
        }
    }
    else { // k-mer has large, dense tree, do an additional x-mer lookup
        uint64_t xmer_entry;
        uint64_t ptr = 0;
        *i += kmerSize;
        mlt_data = &iaux->mlt_table[start_addr];
        hashval = getHashKey(&raux->read_buf[*i], xmerSize, *i, raux->l_seq, 0, &idx_first_N);
        memcpy(&raux->mh_start_addr, &mlt_data[byte_idx], 4);
        byte_idx += 4;
        memcpy(&xmer_entry, &mlt_data[byte_idx + (hashval << 3)], 8);
        code = xmer_entry & METADATA_MASK;
        ptr = xmer_entry >> KMER_DATA_BITWIDTH;
        if (idx_first_N != -1) {
            *i += xmerSize;
            mem->rc_end = *i;
            return;
        }
        if (code == INVALID) {
            *i += xmerSize;
            mem->rc_end = *i;
        }
        else if (code == SINGLE_HIT_LEAF) {
            byte_idx = ptr;
            byte_idx++;
            memcpy(&ref_pos, &mlt_data[byte_idx], 5);                         
            byte_idx += 5; 
            mem->hitcount += 1;
            kv_push(uint64_t, *hits, ref_pos >> 1);
            *i += xmerSize;
            mem->rc_end = *i;
        }
        else {
            byte_idx = ptr;
            *i += xmerSize; 
            if (*i < raux->l_seq) {
                getNextByteIdx_backward(raux, mlt_data, &byte_idx, i, mem, hits);
            }
            else {
                mem->rc_end = *i;
            }
        }
    }
}

/**
 * Main backward search function (reseeding). Lookup up k-mer and/or x-mer table and identify root of ERT.
 * Return early if root node has fewer than 'raux->limit' hits 
 * 
 * @param iaux              index parameters
 * @param raux              read parameters
 * @param i                 index into read buffer
 * @param mem               maximal-exact-match storage
 * @param hits              list of hits for read 
 */
void leftExtend_wlimit(index_aux_t* iaux, read_aux_t* raux, int* i, mem_t* mem, u64v* hits) {
    
    uint8_t code;
    uint8_t* mlt_data;
    uint64_t byte_idx = 0, kmer_entry = 0, start_addr = 0;
    uint32_t hashval = 0;
    int idx_first_N = -1;
    hashval = getHashKey(&raux->read_buf[*i], kmerSize, *i, raux->l_seq, 0, &idx_first_N);
    if (idx_first_N != -1) {
        *i += (kmerSize + xmerSize);
        mem->rc_end = *i;
        return;
    }
    // index-table lookup
    kmer_entry = iaux->kmer_offsets[hashval];
    // index-table entry type
    code = kmer_entry & METADATA_MASK;
    // pointer to root of tree
    start_addr = kmer_entry >> KMER_DATA_BITWIDTH;
    // width used for internal pointers in tree
    raux->ptr_width = (((kmer_entry >> 22) & 3) == 0) ? 4 : ((kmer_entry >> 22) & 3);
    // # hits for k-mer (0 if > 20 hits)
    raux->num_hits = (kmer_entry >> 17) & 0x1F;
    byte_idx = 0;
    if (code == INVALID) { // k-mer absent
        *i += (kmerSize + xmerSize);
        mem->rc_end = *i;
    }
    else if (code == SINGLE_HIT_LEAF) { // single-hit k-mer
        *i += (kmerSize + xmerSize);        
        mem->rc_end = *i;
    }
    else if (code == INFREQUENT) { // k-mer has fewer than 256 hits
        *i += kmerSize;
        mlt_data = &iaux->mlt_table[start_addr];
        if ((raux->num_hits == 0) || (raux->num_hits >= raux->limit)) {
            if (*i < raux->l_seq) {
                // 
                // First 4 bytes of tree store the start of all multi-hit leaves for the k-mer
                //
                memcpy(&raux->mh_start_addr, &mlt_data[byte_idx], 4);
                byte_idx += 4;
                getNextByteIdx_backward_wlimit(raux, mlt_data, &byte_idx, i, mem, hits);
            }
            else {
                // Leaf gathering to be done later
                mem->rc_end = *i;
                mem->fetch_leaves = 1;
            }
        }
        else {
            mem->rc_end = *i;
        }
    }
    else { // k-mer has large, dense tree, do an additional x-mer lookup
        uint64_t xmer_entry;
        uint64_t ptr = 0;
        *i += kmerSize;
        mlt_data = &iaux->mlt_table[start_addr];
        hashval = getHashKey(&raux->read_buf[*i], xmerSize, *i, raux->l_seq, 0, &idx_first_N);
        // 
        // First 4 bytes of tree store the start of all multi-hit leaves for the k-mer
        // This helps in decoding nodes and creates compact trees
        //
        memcpy(&raux->mh_start_addr, &mlt_data[byte_idx], 4);
        byte_idx += 4;
        memcpy(&xmer_entry, &mlt_data[byte_idx + (hashval << 3)], 8);
        code = xmer_entry & METADATA_MASK;
        ptr = xmer_entry >> KMER_DATA_BITWIDTH;
        raux->num_hits = (xmer_entry >> 17) & 0x1F;
        if (idx_first_N != -1) {
            *i += xmerSize;
            mem->rc_end = *i;
            return;
        }
        if (code == INVALID) {
            *i += xmerSize;
            mem->rc_end = *i;
        }
        else if (code == SINGLE_HIT_LEAF) {
            *i += xmerSize;            
            mem->rc_end = *i;
        }
        else {
            byte_idx = ptr;
            *i += xmerSize; 
            if ((raux->num_hits == 0) || (raux->num_hits >= raux->limit)) {
                if (*i < raux->l_seq) {
                    getNextByteIdx_backward_wlimit(raux, mlt_data, &byte_idx, i, mem, hits);
                }   
                else {
                    // Leaf gathering
                    mem->rc_end = *i;
                    mem->fetch_leaves = 1;
                }
            }
            else {
                mem->rc_end = *i;
            }
        }
    }
}

/**
 * Fetch hits for all MEMs identified after backward search (termination conditions based on reseeding)
 * 
 * @param raux              read parameters
 * @param mlt_data          radix tree of k-mer
 * @param byteIdx           byte index into the mlt_data radix tree
 * @param idx               index into read buffer
 * @param mem               maximal-exact-match storage
 * @param hits              list of hits for read 
 */
void getNextByteIdx_fetch_leaves_prefix_reseed(read_aux_t* raux, uint8_t* mlt_data, uint64_t* byte_idx, int idx, mem_t* mem, path_v* visited, u64v* hits) {

    uint64_t nextByteIdx = *byte_idx;
    uint64_t parent_byte_idx = nextByteIdx;
    uint64_t ref_pos = 0;
    uint8_t c;
    int i = idx;
    assert(raux->read_buf[i] != 4); // Should not see N in SMEMs
    c = 3 - raux->read_buf[i];
    uint8_t code = mlt_data[nextByteIdx++];
    uint8_t code_c = (code >> (c << 1)) & 3;
    assert(code != 0);
    if (code_c == EMPTY) { // Do leaf gathering
        mem->end = i;
        int mem_len = mem->end - mem->start;
        if (mem_len >= raux->min_seed_len) {
            leaf_gather(raux, mlt_data, &parent_byte_idx, mem, hits);
        }
    }
    else if (code_c == LEAF) {
        int k;
        uint64_t leaf_data = 0;
        nextByteIdx += getOffsetToLeafData(raux, code, c);
        memcpy(&leaf_data, &mlt_data[nextByteIdx], 5);
        if (leaf_data & 1) {
            nextByteIdx = raux->mh_start_addr + (leaf_data >> 1);
            memcpy(&raux->num_hits, &mlt_data[nextByteIdx], 2);
            nextByteIdx += 2;
        }
        else {
            raux->num_hits = 1;
        }
        if (raux->num_hits >= raux->limit) {
            mem->hitcount += raux->num_hits;
            for (k = 0; k < raux->num_hits; ++k) {
                memcpy(&ref_pos, &mlt_data[nextByteIdx], 5);
                kv_push(uint64_t, *hits, ref_pos >> 1);
                nextByteIdx += 5;
                ref_pos = 0;
            }
            i += 1;
            mem->end = i;
	        mem->is_multi_hit = 1; // decompress leaf node for potentially longer match
        }
        else {
            mem->end = i;
            int mem_len = mem->end - mem->start;
            if (mem_len >= raux->min_seed_len) {
                node_info_t* tmp_node_info = &kv_pop(*visited);
                leaf_gather(raux, mlt_data, &tmp_node_info->byte_idx, mem, hits);
            } 
        }
    } 
    else if (code_c == UNIFORM) {
        uint32_t j;
        int countBP = mlt_data[nextByteIdx++];
        int numBitsForBP = countBP << 1;
        int numBytesForBP = (numBitsForBP % 8) ? (numBitsForBP / 8 + 1) : (numBitsForBP / 8);
        uint8_t packedBP[numBytesForBP];
        memcpy(packedBP, &mlt_data[nextByteIdx], numBytesForBP);
        nextByteIdx += numBytesForBP;
        // Unpack base pairs
        uint8_t unpackedBP[countBP];
        for (j = 0; j < countBP; ++j) {
            unpackedBP[j] = ((packedBP[j >> 2] >> ((~(j) & 3) << 1)) & 3); 
        }
        // Count number of matching base pairs with read
        for (j = 0; j < countBP; ++j) {
            if ((i + j) >= raux->l_seq) {
                break;
            }
            if (3 - raux->read_buf[i+j] != unpackedBP[j]) {
                break;
            }
        }  
        i += j;
        if (j == countBP) {
            if (i < raux->l_seq) {
                getNextByteIdx_fetch_leaves_prefix_reseed(raux, mlt_data, &nextByteIdx, i, mem, visited, hits);
            }
            else {
		        mem->end = i;
                int mem_len = mem->end - mem->start;
                if (mem_len >= raux->min_seed_len) {
                    leaf_gather(raux, mlt_data, &nextByteIdx, mem, hits);
                }
            }
        }
        else {
            mem->end = i;
            int mem_len = mem->end - mem->start;
            if (mem_len >= raux->min_seed_len) {
                leaf_gather(raux, mlt_data, &nextByteIdx, mem, hits);
            }
        }
    }
    else if (code_c == DIVERGE) {
        raux->num_hits = 0;
        getOffsetToChildNode(raux, mlt_data, code, c, &nextByteIdx);
        if ((raux->num_hits == 0) || (raux->num_hits >= raux->limit)) {
            node_info_t nif;
            nif.byte_idx = nextByteIdx;
            nif.num_hits = raux->num_hits;
            kv_push(node_info_t, *visited, nif);
            i += 1;
            if (i < raux->l_seq) {
                getNextByteIdx_fetch_leaves_prefix_reseed(raux, mlt_data, &nextByteIdx, i, mem, visited, hits);
            }
            else {
                mem->end = i;
                int mem_len = mem->end - mem->start;
                if (mem_len >= raux->min_seed_len) {
                    leaf_gather(raux, mlt_data, &nextByteIdx, mem, hits);
                } 
            }
        }
        else {
	        mem->end = i;
            int mem_len = mem->end - mem->start;
            if (mem_len >= raux->min_seed_len) {
                node_info_t* tmp_node_info = &kv_pop(*visited);
                leaf_gather(raux, mlt_data, &tmp_node_info->byte_idx, mem, hits); 
            } 
        }
    }
    *byte_idx = nextByteIdx;    

}

/**
 * Fetch hits for all MEMs identified after backward search (extend beyond mem->end)
 * 
 * @param raux              read parameters
 * @param mlt_data          radix tree of k-mer
 * @param byteIdx           byte index into the mlt_data radix tree
 * @param idx               index into read buffer
 * @param mem               maximal-exact-match storage
 * @param hits              list of hits for read 
 */
void getNextByteIdx_fetch_leaves_prefix(read_aux_t* raux, uint8_t* mlt_data, uint64_t* byte_idx, int idx, mem_t* mem, u64v* hits) {

    uint64_t nextByteIdx = *byte_idx;
    uint64_t ref_pos = 0;
    uint8_t c;
    int i = idx;
    assert(raux->read_buf[i] != 4); // Should not see N in SMEMs
    c = 3 - raux->read_buf[i];
    uint8_t code = mlt_data[nextByteIdx++];
    uint8_t code_c = (code >> (c << 1)) & 3;
    assert(code != 0);
    if (code_c == EMPTY) { // Do leaf gathering
        mem->end = i;
        int mem_len = mem->end - mem->start;
        if (mem_len >= raux->min_seed_len) {
            nextByteIdx = *byte_idx;
            leaf_gather(raux, mlt_data, &nextByteIdx, mem, hits);
        }
    }
    else if (code_c == LEAF) {
        int k;
        uint64_t leaf_data = 0;
        nextByteIdx += getOffsetToLeafData(raux, code, c);
        memcpy(&leaf_data, &mlt_data[nextByteIdx], 5);
        if (leaf_data & 1) {
            nextByteIdx = raux->mh_start_addr + (leaf_data >> 1);
            memcpy(&raux->num_hits, &mlt_data[nextByteIdx], 2);
            nextByteIdx += 2;
            mem->hitcount += raux->num_hits;
            for (k = 0; k < raux->num_hits; ++k) {
                memcpy(&ref_pos, &mlt_data[nextByteIdx], 5);
                kv_push(uint64_t, *hits, ref_pos >> 1);
                nextByteIdx += 5;
                ref_pos = 0;
            }
        }
        else {
            raux->num_hits = 1;
            mem->hitcount += raux->num_hits;
            kv_push(uint64_t, *hits, leaf_data >> 1);
        }
        i += 1;
        mem->end = i;
    } 
    else if (code_c == UNIFORM) {
        uint32_t j;
        int countBP = mlt_data[nextByteIdx++];
        int numBitsForBP = countBP << 1;
        int numBytesForBP = (numBitsForBP % 8) ? (numBitsForBP / 8 + 1) : (numBitsForBP / 8);
        uint8_t packedBP[numBytesForBP];
        memcpy(packedBP, &mlt_data[nextByteIdx], numBytesForBP);
        nextByteIdx += numBytesForBP;
        // Unpack base pairs
        uint8_t unpackedBP[countBP];
        for (j = 0; j < countBP; ++j) {
            unpackedBP[j] = ((packedBP[j >> 2] >> ((~(j) & 3) << 1)) & 3); 
        }
        // Count number of matching base pairs with read
        for (j = 0; j < countBP; ++j) {
            if ((i + j) >= raux->l_seq) {
                break;
            }
            if (3 - raux->read_buf[i+j] != unpackedBP[j]) {
                break;
            }
        }  
        i += j;
        if (j == countBP) {
            if (i < raux->l_seq) {
                getNextByteIdx_fetch_leaves_prefix(raux, mlt_data, &nextByteIdx, i, mem, hits);
            }
            else {
		        mem->end = i;
                int mem_len = mem->end - mem->start;
                if (mem_len >= raux->min_seed_len) {
                    leaf_gather(raux, mlt_data, &nextByteIdx, mem, hits);
                }
            }
        }
        else {
            mem->end = i;
            int mem_len = mem->end - mem->start;
            if (mem_len >= raux->min_seed_len) {
                leaf_gather(raux, mlt_data, &nextByteIdx, mem, hits);
            }
        }
    }
    else if (code_c == DIVERGE) {
        raux->num_hits = 0;
        getOffsetToChildNode(raux, mlt_data, code, c, &nextByteIdx);
        i += 1;
        if (i < raux->l_seq) {
            getNextByteIdx_fetch_leaves_prefix(raux, mlt_data, &nextByteIdx, i, mem, hits);
        }
        else {
	        mem->end = i;
            int mem_len = mem->end - mem->start;
            if (mem_len >= raux->min_seed_len) {
                leaf_gather(raux, mlt_data, &nextByteIdx, mem, hits);
            }
        }
    }

}

/**
 * Forward traversal to fetch hits for all MEMs identified after backward search.
 * 
 * @param raux              read parameters
 * @param mlt_data          radix tree of k-mer
 * @param byteIdx           byte index into the mlt_data radix tree
 * @param idx               index into read buffer
 * @param mem               maximal-exact-match storage
 * @param hits              list of hits for read 
 */
void getNextByteIdx_fetch_leaves(read_aux_t* raux, uint8_t* mlt_data, uint64_t* byte_idx, int idx, mem_t* mem, u64v* hits) {

    uint64_t nextByteIdx = *byte_idx;
    uint64_t ref_pos = 0;
    uint8_t c;
    int i = idx;
    assert(raux->read_buf[i] != 4); // Should not see N in SMEMs
    c = 3 - raux->read_buf[i];
    uint8_t code = mlt_data[nextByteIdx++];
    uint8_t code_c = (code >> (c << 1)) & 3;
    assert(code != 0);
    assert(code_c != EMPTY);
    if (code_c == LEAF) {
        int k;
        uint64_t leaf_data = 0;
        nextByteIdx += getOffsetToLeafData(raux, code, c);
        memcpy(&leaf_data, &mlt_data[nextByteIdx], 5);
        if (leaf_data & 1) {
            nextByteIdx = raux->mh_start_addr + (leaf_data >> 1);
            memcpy(&raux->num_hits, &mlt_data[nextByteIdx], 2);
            nextByteIdx += 2;
            mem->hitcount += raux->num_hits;
            for (k = 0; k < raux->num_hits; ++k) {
                memcpy(&ref_pos, &mlt_data[nextByteIdx], 5);
                kv_push(uint64_t, *hits, ref_pos >> 1);
                nextByteIdx += 5;
                ref_pos = 0;
            }
        }
        else {
            raux->num_hits = 1;
            mem->hitcount += raux->num_hits;
            kv_push(uint64_t, *hits, leaf_data >> 1);
        }
        i += 1;
    } 
    else if (code_c == UNIFORM) {
        uint32_t j;
        int countBP = mlt_data[nextByteIdx++];
        int numBitsForBP = countBP << 1;
        int numBytesForBP = (numBitsForBP % 8) ? (numBitsForBP / 8 + 1) : (numBitsForBP / 8);
        uint8_t packedBP[numBytesForBP];
        memcpy(packedBP, &mlt_data[nextByteIdx], numBytesForBP);
        nextByteIdx += numBytesForBP;
        // Unpack base pairs
        uint8_t unpackedBP[countBP];
        for (j = 0; j < countBP; ++j) {
            unpackedBP[j] = ((packedBP[j >> 2] >> ((~(j) & 3) << 1)) & 3); 
        }
        // Count number of matching base pairs with read
        for (j = 0; j < countBP; ++j) {
            if ((i + j) >= raux->l_seq) {
                break;
            }
            if (3 - raux->read_buf[i+j] != unpackedBP[j]) {
                break;
            }
        }  
        i += j;
        if (j == countBP) {
            if (i < mem->end) { // only extend till end of MEM found during previous backward search
                getNextByteIdx_fetch_leaves(raux, mlt_data, &nextByteIdx, i, mem, hits);
            }
            else {
                leaf_gather(raux, mlt_data, &nextByteIdx, mem, hits);
            }
        }
        else {
            leaf_gather(raux, mlt_data, &nextByteIdx, mem, hits);
        }
    }
    else if (code_c == DIVERGE) {
        raux->num_hits = 0;
        getOffsetToChildNode(raux, mlt_data, code, c, &nextByteIdx);
        i += 1;
        if (i < mem->end) { // only extend till end of MEM found during previous backward search
            getNextByteIdx_fetch_leaves(raux, mlt_data, &nextByteIdx, i, mem, hits);
        }
        else {
            leaf_gather(raux, mlt_data, &nextByteIdx, mem, hits);
        } 
    }

}

/**
 * Extend as much as possible to the right based on reseeding criteria.
 * 
 * Return early if root node has fewer than 'raux->limit' hits 
 * 
 * @param iaux              index parameters
 * @param raux              read parameters
 * @param mem               maximal-exact-match storage
 * @param hits              list of hits for read 
 */
void rightExtend_fetch_leaves_prefix_reseed(index_aux_t* iaux, read_aux_t* raux, mem_t* mem, u64v* hits) {

    uint8_t code;
    uint8_t* mlt_data;
    uint64_t byte_idx = 0, ref_pos = 0, kmer_entry = 0, start_addr = 0;
    uint32_t hashval = 0;
    int flag = 0;
    int i = mem->start;
    int end = mem->end;
    int idx_first_N = -1;
    hashval = getHashKey(&raux->read_buf[i], kmerSize, i, raux->l_seq, &flag, &idx_first_N);
    // index-table lookup
    kmer_entry = iaux->kmer_offsets[hashval];
    // index-table entry type
    code = kmer_entry & METADATA_MASK;
    // pointer to root of tree
    start_addr = kmer_entry >> KMER_DATA_BITWIDTH;
    raux->mh_start_addr = 0;
    // width used for internal pointers in tree
    raux->ptr_width = (((kmer_entry >> 22) & 3) == 0) ? 4 : ((kmer_entry >> 22) & 3);
    // # hits for k-mer (0 if > 20 hits)
    raux->num_hits = (kmer_entry >> 17) & 0x1F;
    byte_idx = 0;
    assert(code != INVALID);
    if (code == SINGLE_HIT_LEAF) {
        mem->end = i;
    }
    else if (code == INFREQUENT) { // k-mer has fewer than 256 hits
        if ((raux->num_hits == 0) || (raux->num_hits >= raux->limit)) {
            i += kmerSize;
            mlt_data = &iaux->mlt_table[start_addr];
            // 
            // First 4 bytes of tree store the start of all multi-hit leaves for the k-mer
            //
            memcpy(&raux->mh_start_addr, &mlt_data[byte_idx], 4);
            byte_idx += 4;
            if (i < raux->l_seq) {
                path_v visited;
                kv_init(visited);
                node_info_t nif;
                nif.byte_idx = byte_idx;
                nif.num_hits = raux->num_hits;
                kv_push(node_info_t, visited, nif);
                getNextByteIdx_fetch_leaves_prefix_reseed(raux, mlt_data, &byte_idx, i, mem, &visited, hits);
                kv_destroy(visited);
            }
            else {
                mem->end = i;
                int mem_len = mem->end - mem->start;
                if (mem_len >= raux->min_seed_len) {
                    leaf_gather(raux, mlt_data, &byte_idx, mem, hits);
                }
            }
        }
        else {
            mem->end = i;
        }
    }
    else if (code == FREQUENT) { // k-mer has large, dense tree, do an additional x-mer lookup
        uint64_t xmer_entry;
        uint64_t ptr = 0;
        mlt_data = &iaux->mlt_table[start_addr];
        hashval = getHashKey(&raux->read_buf[i + kmerSize], xmerSize, i + kmerSize, raux->l_seq, 0, &idx_first_N);
        memcpy(&raux->mh_start_addr, &mlt_data[byte_idx], 4);
        byte_idx += 4;
        memcpy(&xmer_entry, &mlt_data[byte_idx + (hashval << 3)], 8);
        code = xmer_entry & METADATA_MASK;
        ptr = xmer_entry >> KMER_DATA_BITWIDTH;
        // 5 bits to encode number of hits at each node
        raux->num_hits = (xmer_entry >> 17) & 0x1F;
        if (code == INVALID) {
            mem->end = i;
        }
        else if (code == SINGLE_HIT_LEAF) {
            mem->end = i;
        }
        else {
            // When a node has greater than 20 (opt->max_mem->intv) hits, we store a 0 in the hits field 
            if ((raux->num_hits == 0) || (raux->num_hits >= raux->limit)) {
                byte_idx = ptr;
                i += (kmerSize + xmerSize); 
                if (i < raux->l_seq) {
                    path_v visited;
                    kv_init(visited);
                    node_info_t nif;
                    nif.byte_idx = byte_idx;
                    nif.num_hits = raux->num_hits;
                    kv_push(node_info_t, visited, nif);
                    getNextByteIdx_fetch_leaves_prefix_reseed(raux, mlt_data, &byte_idx, i, mem, &visited, hits);
                    kv_destroy(visited);
                }
                else {
                    mem->end = i;
                    int mem_len = mem->end - mem->start;
                    if (mem_len >= raux->min_seed_len) {
                        leaf_gather(raux, mlt_data, &byte_idx, mem, hits);
                    }
                }
            }
            else {
                mem->end = i;
            }
        }
    }
}

/**
 * Extend as much as possible to the right and gather leaves if MEM length >= opt->min_seed_len.
 * 
 * @param iaux              index parameters
 * @param raux              read parameters
 * @param mem               maximal-exact-match storage
 * @param hits              list of hits for read 
 */
void rightExtend_fetch_leaves_prefix(index_aux_t* iaux, read_aux_t* raux, mem_t* mem, u64v* hits) {

    uint8_t code;
    uint8_t* mlt_data;
    uint64_t byte_idx = 0, ref_pos = 0, kmer_entry = 0, start_addr = 0;
    uint32_t hashval = 0;
    int flag = 0;
    int i = mem->start;
    int end = mem->end;
    int idx_first_N = -1;
    hashval = getHashKey(&raux->read_buf[i], kmerSize, i, raux->l_seq, &flag, &idx_first_N);
    // index-table lookup
    kmer_entry = iaux->kmer_offsets[hashval];
    // index-table entry type
    code = kmer_entry & METADATA_MASK;
    // pointer to root of tree
    start_addr = kmer_entry >> KMER_DATA_BITWIDTH;
    raux->mh_start_addr = 0;
    // width used for internal pointers in tree, 2 bits, ptr_width = 4 is encoded as 0
    raux->ptr_width = (((kmer_entry >> 22) & 3) == 0) ? 4 : ((kmer_entry >> 22) & 3);
    byte_idx = 0;
    assert(code != INVALID);
    if (code == SINGLE_HIT_LEAF) {
        mlt_data = &iaux->mlt_table[start_addr];
        byte_idx++;
        memcpy(&ref_pos, &mlt_data[byte_idx], 5);                         
        mem->hitcount += 1;
        kv_push(uint64_t, *hits, ref_pos >> 1);
        byte_idx += 5;                                                    
        i += kmerSize;
        mem->end = i;
    }
    else if (code == INFREQUENT) { // k-mer has fewer than 256 hits
        i += kmerSize;
        mlt_data = &iaux->mlt_table[start_addr];
        // 
        // First 4 bytes of tree store the start of all multi-hit leaves for the k-mer
        //
        memcpy(&raux->mh_start_addr, &mlt_data[byte_idx], 4);
        byte_idx += 4;
        if (i < raux->l_seq) {
            getNextByteIdx_fetch_leaves_prefix(raux, mlt_data, &byte_idx, i, mem, hits);
        }
        else {
            mem->end = i;
            int mem_len = mem->end - mem->start;
            if (mem_len >= raux->min_seed_len) {
                leaf_gather(raux, mlt_data, &byte_idx, mem, hits);
            }
        }
    }
    else if (code == FREQUENT) { // k-mer has large, dense tree, do an additional x-mer lookup
        uint64_t xmer_entry;
        uint64_t ptr = 0;
        mlt_data = &iaux->mlt_table[start_addr];
        hashval = getHashKey(&raux->read_buf[i + kmerSize], xmerSize, i + kmerSize, raux->l_seq, 0, &idx_first_N);
        memcpy(&raux->mh_start_addr, &mlt_data[byte_idx], 4);
        byte_idx += 4;
        memcpy(&xmer_entry, &mlt_data[byte_idx + (hashval << 3)], 8);
        code = xmer_entry & METADATA_MASK;
        ptr = xmer_entry >> KMER_DATA_BITWIDTH;
        if (code == INVALID) {
            mem->end = i;
        }
        else if (code == SINGLE_HIT_LEAF) {
            byte_idx = ptr;
            byte_idx++;
            memcpy(&ref_pos, &mlt_data[byte_idx], 5);                         
            mem->hitcount += 1;
            kv_push(uint64_t, *hits, ref_pos >> 1);
            byte_idx += 5;                                                    
            i += (kmerSize + xmerSize);
            mem->end = i;
        }
        else {
            byte_idx = ptr;
            i += (kmerSize + xmerSize); 
            if (i < raux->l_seq) {
                getNextByteIdx_fetch_leaves_prefix(raux, mlt_data, &byte_idx, i, mem, hits);
            }
            else {
                mem->end = i;
                int mem_len = mem->end - mem->start;
                if (mem_len >= raux->min_seed_len) {
                    leaf_gather(raux, mlt_data, &byte_idx, mem, hits);
                }
            }
        }
    }
}

/**
 * Fetch hits for all MEMs identified after backward search.
 * 
 * Note that backward search functions above only perform tree traversal without gathering hits as they will be fetched in
 * a different order than required by BWA-MEM (i.e., all hits for MEM need to be sorted by right-context). We re-traverse and
 * the tree and fetch the hits in the correct order below (this has a minor performance penalty)
 * 
 * Return early if root node has fewer than 'raux->limit' hits 
 * 
 * @param iaux              index parameters
 * @param raux              read parameters
 * @param mem               maximal-exact-match storage
 * @param hits              list of hits for read 
 */
void rightExtend_fetch_leaves(index_aux_t* iaux, read_aux_t* raux, mem_t* mem, u64v* hits) {

    uint8_t code;
    uint8_t* mlt_data;
    uint64_t byte_idx = 0, kmer_entry = 0, start_addr = 0;
    uint32_t hashval = 0;
    int flag = 0;
    int i = mem->start;
    int end = mem->end;
    int idx_first_N = -1;
    hashval = getHashKey(&raux->read_buf[i], kmerSize, i, raux->l_seq, &flag, &idx_first_N);
    // index-table lookup
    kmer_entry = iaux->kmer_offsets[hashval];
    // index-table entry type
    code = kmer_entry & METADATA_MASK;
    // pointer to root of tree
    start_addr = kmer_entry >> KMER_DATA_BITWIDTH;
    raux->mh_start_addr = 0;
    // width used for internal pointers in tree
    raux->ptr_width = (((kmer_entry >> 22) & 3) == 0) ? 4 : ((kmer_entry >> 22) & 3);
    byte_idx = 0;
    assert(code != INVALID);
    assert(code != SINGLE_HIT_LEAF);
    if (code == INFREQUENT) { // k-mer has fewer than 256 hits
        i += kmerSize;
        mlt_data = &iaux->mlt_table[start_addr];
        // 
        // First 4 bytes of tree store the start of all multi-hit leaves for the k-mer
        //
        memcpy(&raux->mh_start_addr, &mlt_data[byte_idx], 4);
        byte_idx += 4;
        if (i < end) {
            getNextByteIdx_fetch_leaves(raux, mlt_data, &byte_idx, i, mem, hits);
        }
        else {
            leaf_gather(raux, mlt_data, &byte_idx, mem, hits);
        }
    }
    else if (code == FREQUENT) { // k-mer has large, dense tree, do an additional x-mer lookup
        uint64_t xmer_entry;
        uint64_t ptr = 0;
        i += kmerSize;
        mlt_data = &iaux->mlt_table[start_addr];
        hashval = getHashKey(&raux->read_buf[i], xmerSize, i, raux->l_seq, 0, &idx_first_N);
        memcpy(&raux->mh_start_addr, &mlt_data[byte_idx], 4);
        byte_idx += 4;
        memcpy(&xmer_entry, &mlt_data[byte_idx + (hashval << 3)], 8);
        code = xmer_entry & METADATA_MASK;
        ptr = xmer_entry >> KMER_DATA_BITWIDTH;
        assert(code != INVALID);
        assert(code != SINGLE_HIT_LEAF);
        byte_idx = ptr;
        i += xmerSize; 
        if (i < end) {
            getNextByteIdx_fetch_leaves(raux, mlt_data, &byte_idx, i, mem, hits);
        }
        else {
            leaf_gather(raux, mlt_data, &byte_idx, mem, hits);
        }
    }
}

/**
 * Main forward search function (seeding). Lookup up k-mer and/or x-mer table and identify root of ERT 
 * 
 * @param iaux              index parameters
 * @param raux              read parameters
 * @param i                 Index into read buffer
 * @param mem               maximal-exact-match storage
 * @param hits              list of hits for read 
 */
void rightExtend(index_aux_t* iaux, read_aux_t* raux, int* i, mem_t* mem, u64v* hits) {
    
    uint8_t code;
    uint8_t* mlt_data;
    uint64_t byte_idx = 0, ref_pos = 0, kmer_entry = 0, start_addr = 0;
    uint32_t hashval = 0;
    uint64_t lep_data = 0;
    int flag = 0;
    int idx_first_N = -1;
    hashval = getHashKey(&raux->read_buf[*i], kmerSize, *i, raux->l_seq, &flag, &idx_first_N);
    // index-table lookup
    kmer_entry = iaux->kmer_offsets[hashval];
    // index-table entry type
    code = kmer_entry & METADATA_MASK;
    lep_data = (kmer_entry >> METADATA_BITWIDTH) & LEP_MASK;
    // pointer to root of tree
    start_addr = kmer_entry >> KMER_DATA_BITWIDTH;
    uint64_t mlt_start_addr = raux->mlt_start_addr = start_addr;
    raux->mh_start_addr = 0;
    // width used for internal pointers in tree
    raux->ptr_width = (((kmer_entry >> 22) & 3) == 0) ? 4 : ((kmer_entry >> 22) & 3);
    // LEP takes up kmerSize-1 bits. Last LEP bit is at position = kmerSize-2.
    if (*i <= 64-kmerSize) {
        raux->lep[0] |= (lep_data << *i);
    }
    else if (*i < 64) { 
        raux->lep[0] |= (lep_data << *i);
        raux->lep[1] |= (lep_data >> (64-*i));
    }
    else if (*i <= 128-kmerSize) {
        raux->lep[1] |= (lep_data << (*i-64));
    }
    else if (*i < 128) {
        raux->lep[1] |= (lep_data << (*i-64));
        raux->lep[2] |= (lep_data >> (128-*i));
    }
    else if (*i <= 192-kmerSize) {
        raux->lep[2] |= (lep_data << (*i-128));
    }
    else if (*i < 192) {
        raux->lep[2] |= (lep_data << (*i-128));
        raux->lep[3] |= (lep_data >> (192-*i));
    }
    else if (*i <= 256-kmerSize) { 
        raux->lep[3] |= (lep_data << (*i-192));
    }
    else if (*i < 256) {
        raux->lep[3] |= (lep_data << (*i-192));
        raux->lep[4] |= (lep_data >> (256-*i));
    }
    else {
        raux->lep[4] |= (lep_data << (*i-256));
    }
    raux->nextLEPBit = *i + kmerSize - 1;
    byte_idx = 0;
    // We found an ambiguous base in the kmer. Stop extension at ambiguous base and record LEP
    if (idx_first_N != -1) {
        if (*i != 0) {
            raux->nextLEPBit = *i + idx_first_N - 1;
            raux->lep[raux->nextLEPBit >> 6] |= (1ULL << (raux->nextLEPBit & (0x3FULL)));
        }
        *i += idx_first_N;
        return;
    }
    if (flag) {
        raux->nextLEPBit = (raux->l_seq - 1); 
        *i = raux->l_seq;
        raux->lep[raux->nextLEPBit >> 6] |= (1ULL << (raux->nextLEPBit & (0x3FULL)));
        return;
    }
    if (code == INVALID) {
        // Do backward extension using LEP
        *i += (kmerSize + xmerSize);
    }
    else if (code == SINGLE_HIT_LEAF) {
        mlt_data = &iaux->mlt_table[start_addr];
        byte_idx++;
        memcpy(&ref_pos, &mlt_data[byte_idx], 5);                         
        mem->hitcount += 1;
        kv_push(uint64_t, *hits, ref_pos >> 1);
        byte_idx += 5;                                                    
        *i += kmerSize;
    }
    else if (code == INFREQUENT) { // k-mer has fewer than 256 hits
        *i += kmerSize;
        mlt_data = &iaux->mlt_table[start_addr];
        if (*i < raux->l_seq) {
            // 
            // First 4 bytes of tree store the start of all multi-hit leaves for the k-mer
            //
            memcpy(&raux->mh_start_addr, &mlt_data[byte_idx], 4);
            byte_idx += 4;
            getNextByteIdx(raux, mlt_data, &byte_idx, i, mem, hits);
        }
        else {
            raux->lep[raux->nextLEPBit >> 6] |= (1ULL << (raux->nextLEPBit & (0x3FULL)));
            raux->nextLEPBit += 1;
        }
    }
    else if (code == FREQUENT) { // k-mer has large, dense tree, do an additional x-mer lookup
        uint64_t xmer_entry;
        uint64_t ptr = 0;
        *i += kmerSize;
        int k;
        uint64_t lepBit = 0;
        flag = 0;
        mlt_data = &iaux->mlt_table[mlt_start_addr];
        hashval = getHashKey(&raux->read_buf[*i], xmerSize, *i, raux->l_seq, &flag, &idx_first_N);
        memcpy(&raux->mh_start_addr, &mlt_data[byte_idx], 4);
        byte_idx += 4;
        memcpy(&xmer_entry, &mlt_data[byte_idx + (hashval << 3)], 8);
        code = xmer_entry & METADATA_MASK;
        lep_data = (xmer_entry >> METADATA_BITWIDTH) & 0xF;
        ptr = xmer_entry >> KMER_DATA_BITWIDTH;
        int xmerLen = 0;
        if (raux->l_seq - *i > xmerSize) {
            xmerLen = xmerSize;
        }
        else {
            xmerLen = raux->l_seq - *i;
        }
        for (k = 0; k < xmerLen; ++k) {
            lepBit = (lep_data >> k) & 1;
            raux->lep[raux->nextLEPBit >> 6] |= (lepBit << (raux->nextLEPBit & (0x3FULL)));
            raux->nextLEPBit++;
        }
        // We found an ambiguous base in the kmer. Stop extension at ambiguous base and record LEP
        if (idx_first_N != -1) {
            raux->nextLEPBit = *i + idx_first_N - 1;
            raux->lep[raux->nextLEPBit >> 6] |= (1ULL << (raux->nextLEPBit & (0x3FULL)));
            *i += idx_first_N;
            return;
        }
        if (flag) {
            raux->nextLEPBit = (raux->l_seq - 1); 
            *i = raux->l_seq;
            raux->lep[raux->nextLEPBit >> 6] |= (1ULL << (raux->nextLEPBit & (0x3FULL)));
            return;
        }
        if (code == INVALID) {
            *i += xmerSize;
        }
        else if (code == SINGLE_HIT_LEAF) {
            byte_idx = ptr;
            byte_idx++;
            memcpy(&ref_pos, &mlt_data[byte_idx], 5);                         
            mem->hitcount += 1;
            kv_push(uint64_t, *hits, ref_pos >> 1);
            byte_idx += 5;                                                    
            *i += xmerSize;
        }
        else {
            byte_idx = ptr;
            *i += xmerSize;
            if (*i < raux->l_seq) {
                getNextByteIdx(raux, mlt_data, &byte_idx, i, mem, hits);
            }
            else {
                raux->lep[raux->nextLEPBit >> 6] |= (1ULL << (raux->nextLEPBit & (0x3FULL)));
                raux->nextLEPBit += 1;
            }
        }
    }
}

/**
 * Main forward search function (reseeding). Lookup up k-mer and/or x-mer table and identify root of ERT 
 * Return early if root node has fewer than 'raux->limit' hits 
 * 
 * @param iaux              index parameters
 * @param raux              read parameters
 * @param i                 Index into read buffer
 * @param mem               maximal-exact-match storage
 * @param hits              list of hits for read 
 */
void rightExtend_wlimit(index_aux_t* iaux, read_aux_t* raux, int* i, mem_t* mem, u64v* hits) {
    
    uint8_t code;
    uint8_t* mlt_data;
    uint64_t byte_idx = 0, kmer_entry = 0, start_addr = 0;
    uint32_t hashval = 0;
    uint64_t lep_data = 0;
    int flag = 0;
    int idx_first_N = -1;
    hashval = getHashKey(&raux->read_buf[*i], kmerSize, *i, raux->l_seq, &flag, &idx_first_N);
    // index-table lookup
    kmer_entry = iaux->kmer_offsets[hashval];
    // index-table entry type
    code = kmer_entry & METADATA_MASK;
    lep_data = (kmer_entry >> METADATA_BITWIDTH) & LEP_MASK;
    // pointer to root of tree
    start_addr = kmer_entry >> KMER_DATA_BITWIDTH;
    uint64_t mlt_start_addr = raux->mlt_start_addr = start_addr;
    raux->mh_start_addr = 0;
    // width used for internal pointers in tree, 2 bits, ptr_width = 4 is encoded as 0
    raux->ptr_width = (((kmer_entry >> 22) & 3) == 0) ? 4 : ((kmer_entry >> 22) & 3);
    raux->num_hits = (kmer_entry >> 17) & 0x1F;
    // LEP takes up kmerSize-1 bits. Last LEP bit is at position = kmerSize-2.
    if (*i <= 64-kmerSize) {
        raux->lep[0] |= (lep_data << *i);
    }
    else if (*i < 64) { 
        raux->lep[0] |= (lep_data << *i);
        raux->lep[1] |= (lep_data >> (64-*i));
    }
    else if (*i <= 128-kmerSize) {
        raux->lep[1] |= (lep_data << (*i-64));
    }
    else if (*i < 128) {
        raux->lep[1] |= (lep_data << (*i-64));
        raux->lep[2] |= (lep_data >> (128-*i));
    }
    else if (*i <= 192-kmerSize) {
        raux->lep[2] |= (lep_data << (*i-128));
    }
    else if (*i < 192) {
        raux->lep[2] |= (lep_data << (*i-128));
        raux->lep[3] |= (lep_data >> (192-*i));
    }
    else if (*i <= 256-kmerSize) { 
        raux->lep[3] |= (lep_data << (*i-192));
    }
    else if (*i < 256) {
        raux->lep[3] |= (lep_data << (*i-192));
        raux->lep[4] |= (lep_data >> (256-*i));
    }
    else {
        raux->lep[4] |= (lep_data << (*i-256));
    }
    raux->nextLEPBit = *i + kmerSize - 1;
    byte_idx = 0;
    // We found an ambiguous base in the kmer. Stop extension at ambiguous base and record LEP
    if (idx_first_N != -1) {
        if (*i != 0) {
            raux->nextLEPBit = *i + idx_first_N - 1;
            raux->lep[raux->nextLEPBit >> 6] |= (1ULL << (raux->nextLEPBit & (0x3FULL)));
        }
        *i += idx_first_N;
        return;
    }
    if (flag) {
        raux->nextLEPBit = (raux->l_seq - 1); 
        *i = raux->l_seq;
        raux->lep[raux->nextLEPBit >> 6] |= (1ULL << (raux->nextLEPBit & (0x3FULL)));
        return;
    }
    if (code == INVALID) {
        // Do backward extension using LEP
        *i += (kmerSize + xmerSize);
    }
    else if (code == SINGLE_HIT_LEAF) {
        *i += (kmerSize + xmerSize);
    }
    else if (code == INFREQUENT) { // k-mer has fewer than 256 hits
        *i += kmerSize;
        mlt_data = &iaux->mlt_table[mlt_start_addr];
        if ((raux->num_hits == 0) || (raux->num_hits >= raux->limit)) {
            if (*i < raux->l_seq) {
                path_v visited;
                kv_init(visited);
                node_info_t nif;
                nif.byte_idx = byte_idx;
                nif.num_hits = raux->num_hits;
                kv_push(node_info_t, visited, nif);
                // 
                // First 4 bytes of tree store the start of all multi-hit leaves for the k-mer
                //
                memcpy(&raux->mh_start_addr, &mlt_data[byte_idx], 4);
                byte_idx += 4;
                getNextByteIdx_wlimit(raux, mlt_data, &byte_idx, i, mem, &visited, hits);
                kv_destroy(visited);
            }
            else {
                raux->lep[raux->nextLEPBit >> 6] |= (1ULL << (raux->nextLEPBit & (0x3FULL)));
                raux->nextLEPBit += 1;
            }
        }
    }
    else if (code == FREQUENT) { // k-mer has large, dense tree, do an additional x-mer lookup
        uint64_t xmer_entry;
        uint64_t ptr = 0;
        *i += kmerSize;
        int k;
        uint64_t lepBit = 0;
        flag = 0;
        mlt_data = &iaux->mlt_table[mlt_start_addr];
        hashval = getHashKey(&raux->read_buf[*i], xmerSize, *i, raux->l_seq, &flag, &idx_first_N);
        memcpy(&raux->mh_start_addr, &mlt_data[byte_idx], 4);
        byte_idx += 4;
        memcpy(&xmer_entry, &mlt_data[byte_idx + (hashval << 3)], 8);
        code = xmer_entry & METADATA_MASK;
        lep_data = (xmer_entry >> METADATA_BITWIDTH) & 0xF;
        ptr = xmer_entry >> KMER_DATA_BITWIDTH;
        // 5 bits to encode hits for each node
        raux->num_hits = (xmer_entry >> 17) & 0x1F;
        int xmerLen = 0;
        if (raux->l_seq - *i > xmerSize) {
            xmerLen = xmerSize;
        }
        else {
            xmerLen = raux->l_seq - *i;
        }
        for (k = 0; k < xmerLen; ++k) {
            lepBit = (lep_data >> k) & 1;
            raux->lep[raux->nextLEPBit >> 6] |= (lepBit << (raux->nextLEPBit & (0x3FULL)));
            raux->nextLEPBit++;
        }
        // We found an ambiguous base in the kmer. Stop extension at ambiguous base and record LEP
        if (idx_first_N != -1) {
            raux->nextLEPBit = *i + idx_first_N - 1;
            raux->lep[raux->nextLEPBit >> 6] |= (1ULL << (raux->nextLEPBit & (0x3FULL)));
            *i += idx_first_N;
            return;
        }
        if (flag) {
            raux->nextLEPBit = (raux->l_seq - 1); 
            *i = raux->l_seq;
            raux->lep[raux->nextLEPBit >> 6] |= (1ULL << (raux->nextLEPBit & (0x3FULL)));
            return;
        }
        if (code == INVALID) {
            *i += xmerSize;
        }
        else if (code == SINGLE_HIT_LEAF) {
            *i += xmerSize;
        }
        else {
            byte_idx = ptr;
            *i += xmerSize;
            if ((raux->num_hits == 0) || (raux->num_hits >= raux->limit)) {
                if (*i < raux->l_seq) {
                    path_v visited;
                    kv_init(visited);
                    node_info_t nif;
                    nif.byte_idx = byte_idx;
                    nif.num_hits = raux->num_hits;
                    kv_push(node_info_t, visited, nif);
                    getNextByteIdx_wlimit(raux, mlt_data, &byte_idx, i, mem, &visited, hits);
                    kv_destroy(visited);
                }
                else {
                    raux->lep[raux->nextLEPBit >> 6] |= (1ULL << (raux->nextLEPBit & (0x3FULL)));
                    raux->nextLEPBit += 1;
                }
            }
        }
    }
}

/*
 * Main forward search function for LAST heuristic
 * 
 * @param iaux              index parameters
 * @param raux              read parameters
 * @param i                 Index into read buffer
 * @param mem               maximal-exact-match storage
 * @param hits              list of hits for read 
 */
void rightExtend_last(index_aux_t* iaux, read_aux_t* raux, int* i, mem_t* mem, u64v* hits) {
    
    uint8_t code;
    uint8_t* mlt_data;
    uint64_t byte_idx = 0, ref_pos = 0, kmer_entry = 0, start_addr = 0;
    int flag = 0;
    uint32_t hashval = 0;
    int idx_first_N = -1;
    hashval = getHashKey(&raux->read_buf[*i], kmerSize, *i, raux->l_seq, &flag, &idx_first_N);
    // We found an ambiguous base in the kmer. Stop extension
    if (idx_first_N != -1) {
        *i += (idx_first_N + 1);
        return;
    }
    if (flag) {
        *i = raux->l_seq;
        return;
    }
    // index-table lookup
    kmer_entry = iaux->kmer_offsets[hashval];
    code = kmer_entry & METADATA_MASK;
    start_addr = kmer_entry >> KMER_DATA_BITWIDTH;
    uint64_t mlt_start_addr = raux->mlt_start_addr = start_addr;
    raux->mh_start_addr = 0;
    raux->ptr_width = (((kmer_entry >> 22) & 3) == 0) ? 4 : ((kmer_entry >> 22) & 3);
    raux->num_hits = (kmer_entry >> 17) & 0x1F;
    byte_idx = 0;
    if (code == INVALID) {
        // Do backward extension using LEP
        *i += kmerSize;
    }
    else if (code == SINGLE_HIT_LEAF) {
        mlt_data = &iaux->mlt_table[mlt_start_addr];
        byte_idx++;
        memcpy(&ref_pos, &mlt_data[byte_idx], 5);                         
        mem->hitcount += 1;
        kv_push(uint64_t, *hits, ref_pos >> 1);
        byte_idx += 5;                                                    
        *i += kmerSize;
    }
    else if (code == INFREQUENT) {
        *i += kmerSize;
        mlt_data = &iaux->mlt_table[mlt_start_addr];
        if (*i < raux->l_seq) { // Length <= (kmer + xmer). Need not check num_hits here.
            memcpy(&raux->mh_start_addr, &mlt_data[byte_idx], 4);
            byte_idx += 4;
            getNextByteIdx_last(raux, mlt_data, &byte_idx, i, mem, hits);
        }
    }
    else if (code == FREQUENT) {
        uint64_t xmer_entry;
        uint64_t ptr = 0;
        *i += kmerSize;
        flag = 0;
        mlt_data = &iaux->mlt_table[mlt_start_addr];
        hashval = getHashKey(&raux->read_buf[*i], xmerSize, *i, raux->l_seq, &flag, &idx_first_N);
        // We found an ambiguous base in the kmer. Stop extension
        if (idx_first_N != -1) {
            *i += (idx_first_N + 1);
            return;
        }
        if (flag) {
            *i = raux->l_seq;
            return;
        }
        memcpy(&raux->mh_start_addr, &mlt_data[byte_idx], 4);
        byte_idx += 4;
        memcpy(&xmer_entry, &mlt_data[byte_idx + (hashval << 3)], 8);
        code = xmer_entry & METADATA_MASK;
        ptr = xmer_entry >> KMER_DATA_BITWIDTH;
        raux->num_hits = (xmer_entry >> 17) & 0x1F;
        if (code == INVALID) {
            *i += xmerSize;
        }
        else if (code == SINGLE_HIT_LEAF) {
            byte_idx = ptr;
            byte_idx++;
            memcpy(&ref_pos, &mlt_data[byte_idx], 5);                         
            mem->hitcount += 1;
            kv_push(uint64_t, *hits, ref_pos >> 1);
            byte_idx += 5;                                                    
            *i += xmerSize;
        }
        else {
            byte_idx = ptr;
            *i += xmerSize;
            if ((raux->num_hits == 0) || (raux->num_hits >= raux->limit) || ((*i - mem->start) < (raux->min_seed_len + 1))) {
                if (*i < raux->l_seq) {
                    getNextByteIdx_last(raux, mlt_data, &byte_idx, i, mem, hits);
                }
            }
            else { // number of hits is less than the limit and seed length >= opt->min_seed_len
                leaf_gather(raux, mlt_data, &byte_idx, mem, hits);
            }
        }
    }
}

/*
 * Initialize MEM parameters for backward extension
 *
 * @param lep               LEP bit vector
 * @param mem               maximal-exact-match
 * @param j                 Index into LEP
 * @param seq_len           Read length
 * @param min_seed_len      Minimum seed length
 *
 * @return mem_valid        valid flag to indicate if backward search needs to be performed from given position j
 */
inline int init_mem(uint64_t* lep, mem_t* mem, int j, int seq_len, int min_seed_len) {
    int lep_bit_set = (lep[j >> 6] >> (j & 63)) & 1; 
    int in_valid_range = (j >= (min_seed_len-1)) ? 1 : 0;
    int mem_valid = (lep_bit_set && in_valid_range); 
    mem->end = j + 1; // [start, end) indexing 
    mem->rc_start = seq_len - j - 1;
    mem->rc_end = mem->rc_start;
    mem->skip_ref_fetch = 0;
    mem->forward = 0;
    mem->fetch_leaves = 0;
    mem->hitbeg = mem->hitcount = 0;
    mem->end_correction = 0;
    mem->is_multi_hit = 0;
    return mem_valid;
}

/*
 * Compute final SMEMs and their hits after considering their overlaps
 *
 * Leaf expansion is also performed to get the actual length of the MEM. 
 * Note that leaf nodes store compressed suffixes by including a pointer to the reference genome
 *
 * @param iaux          index related parameters
 * @param raux          read related parameters
 * @param mem           maximal-exact-match storage
 * @param sh            helper data structure to keep track of start and end positions of previously identified MEMs
 * @param smems         List of SMEMs
 *
 * @return n            number of backward extensions to skip
 */
int check_and_add_smem_prefix_reseed(index_aux_t* iaux, read_aux_t* raux, mem_t* mem, smem_helper_t* sh, mem_v* smems, u64v* hits) {

    mem->start = raux->l_seq - mem->rc_end; // Adjust start position of LMEM
    int lmemLen = mem->end - mem->start, rmemLen = -1, next_be_point;
    if (mem->hitcount > 0 && !mem->skip_ref_fetch) {
        int64_t len;
        int64_t start_ref_pos = hits->a[mem->hitbeg] - mem->rc_start;
        int64_t end_ref_pos = hits->a[mem->hitbeg];
        // Fetch reference to check for extra matching bps on the right
        uint8_t* rseq = bns_get_seq_v2(iaux->bns->l_pac, iaux->pac, start_ref_pos, end_ref_pos, &len, 0); 
        int m, numMatchingBP = 0;
        for (m = 1; m <= len; ++m) {
            if (rseq[mem->rc_start - m] == raux->read_buf[mem->rc_start - m]) {
                numMatchingBP++;
            }
            else {
                break;
            }
        }
        mem->end += numMatchingBP;
    	mem->end_correction += numMatchingBP;
    
        start_ref_pos = hits->a[mem->hitbeg] + lmemLen;
        end_ref_pos = start_ref_pos + mem->start;
        
        // Fetch reference to check for extra matching bps on the left
        rseq = bns_get_seq_v2(iaux->bns->l_pac, iaux->pac, start_ref_pos, end_ref_pos, &len, 0); 
        numMatchingBP = 0;
        for (m = 0; m < len; ++m) {
            if (rseq[m] == raux->read_buf[mem->rc_end + m]) {
                numMatchingBP++;
            }
            else {
                break;
            }
        }
        mem->start -= numMatchingBP;
    }
    // Adjust start position of MEM by extra matching bps
    lmemLen = mem->end - mem->start;
    next_be_point = mem->end;
    if (mem->hitcount == 1) {
        if (lmemLen >= raux->min_seed_len) {
            kv_push(mem_t, *smems, *mem);
        }
        else {
            next_be_point += (raux->min_seed_len - lmemLen);
        }   
    }
    // perform forward extension only for non-empty MEMs from backward extension
    else if (mem->fetch_leaves && (mem->start <= (raux->l_seq - raux->min_seed_len))) {
        hits->n -= mem->hitcount;
        mem->hitbeg = hits->n;
        mem->hitcount = 0;
        raux->read_buf = raux->unpacked_queue_buf;
        rightExtend_fetch_leaves_prefix_reseed(iaux, raux, mem, hits);
        raux->read_buf = raux->unpacked_rc_queue_buf;
        rmemLen = mem->end - mem->start;
        next_be_point = mem->end;
        if (mem->hitcount > 0) {
            if (mem->is_multi_hit) {
                int64_t len;
		        int64_t start_ref_pos = hits->a[mem->hitbeg] + rmemLen;
		        int64_t end_ref_pos = hits->a[mem->hitbeg] + raux->l_seq - mem->start;
		        /// Fetch reference
		        uint8_t* rseq = bns_get_seq_v2(iaux->bns->l_pac, iaux->pac, start_ref_pos, end_ref_pos, &len, 0); 
		        int m;
		        int numMatchingBP = 0;
                /// Check for matching bases
                for (m = 0; m < len; ++m) {
                    if (rseq[m] == raux->unpacked_queue_buf[mem->end + m]) {
                        numMatchingBP++;
                    }
                    else {
                        break;
                    }
                }
                mem->end += numMatchingBP;
                rmemLen = mem->end - mem->start;
                next_be_point = mem->end;
            }
            if (rmemLen >= raux->min_seed_len && mem->end <= sh->mem_end_limit) {
                kv_push(mem_t, *smems, *mem);
            }
            else {
                next_be_point += (raux->min_seed_len - rmemLen);       
            }           
        }
        else { // we don't have min_seed_len match for this start position
            assert(rmemLen <= raux->min_seed_len);
            next_be_point += (raux->min_seed_len - rmemLen);       
        }
    }
    else {
	    if (lmemLen <= raux->min_seed_len) {
            next_be_point += (raux->min_seed_len - lmemLen);
        }
    }

    return next_be_point;
}

/*
 * Compute final SMEMs and their hits after considering their overlaps
 *
 * Leaf expansion is also performed to get the actual length of the MEM. 
 * Note that leaf nodes store compressed suffixes by including a pointer to the reference genome
 *
 * @param iaux          index related parameters
 * @param raux          read related parameters
 * @param mem           maximal-exact-match storage
 * @param sh            helper data structure to keep track of start and end positions of previously identified MEMs
 * @param smems         list of SMEMs
 * @param hits          list of hits for read
 *
 * @return n            number of backward extensions to skip
 */
int check_and_add_smem_prefix(index_aux_t* iaux, read_aux_t* raux, mem_t* mem, smem_helper_t* sh, mem_v* smems, u64v* hits) {

    mem->start = raux->l_seq - mem->rc_end; // Adjust start position of LMEM
    int lmemLen = mem->end - mem->start, rmemLen = -1, next_be_point;
    if (mem->hitcount > 0 && !mem->skip_ref_fetch) {
        int64_t len;
        int64_t start_ref_pos = hits->a[mem->hitbeg] - mem->rc_start;
        int64_t end_ref_pos = hits->a[mem->hitbeg];
        // Fetch reference to check for extra matching bps
        uint8_t* rseq = bns_get_seq_v2(iaux->bns->l_pac, iaux->pac, start_ref_pos, end_ref_pos, &len, 0); 
        int m, numMatchingBP = 0;
        for (m = 1; m <= len; ++m) {
            if (rseq[mem->rc_start - m] == raux->read_buf[mem->rc_start - m]) {
                numMatchingBP++;
            }
            else {
                break;
            }
        }
        mem->end += numMatchingBP;
    	mem->end_correction += numMatchingBP;
    
        start_ref_pos = hits->a[mem->hitbeg] + lmemLen;
        end_ref_pos = start_ref_pos + mem->start;
        // Fetch reference to check for extra matching bps
        rseq = bns_get_seq_v2(iaux->bns->l_pac, iaux->pac, start_ref_pos, end_ref_pos, &len, 0); 
        numMatchingBP = 0;
        for (m = 0; m < len; ++m) {
            if (rseq[m] == raux->read_buf[mem->rc_end + m]) {
                numMatchingBP++;
            }
            else {
                break;
            }
        }
        // free(rseq);
        mem->start -= numMatchingBP;
    }
    lmemLen = mem->end - mem->start;
    next_be_point = mem->end;
    // skip forward re-traversal for single-git MEMs
    if (mem->hitcount == 1) {
        if (lmemLen >= raux->min_seed_len) {
            kv_push(mem_t, *smems, *mem);
        }
        else {
            next_be_point += (raux->min_seed_len - lmemLen);
        }   
    }
    else if (mem->fetch_leaves && (mem->start <= (raux->l_seq - raux->min_seed_len))) {
        hits->n -= mem->hitcount;
        mem->hitbeg = hits->n;
        mem->hitcount = 0;
        raux->read_buf = raux->unpacked_queue_buf;
        rightExtend_fetch_leaves_prefix(iaux, raux, mem, hits);
        raux->read_buf = raux->unpacked_rc_queue_buf;
        rmemLen = mem->end - mem->start;
        next_be_point = mem->end;
        if (mem->hitcount > 0) {
            int64_t len;
            int64_t start_ref_pos = hits->a[mem->hitbeg] + rmemLen;
            int64_t end_ref_pos = hits->a[mem->hitbeg] + raux->l_seq - mem->start;
            /// Fetch reference
            uint8_t* rseq = bns_get_seq_v2(iaux->bns->l_pac, iaux->pac, start_ref_pos, end_ref_pos, &len, 0); 
            int m;
            int numMatchingBP = 0;
            /// Check for matching bases
            for (m = 0; m < len; ++m) {
                if (rseq[m] == raux->unpacked_queue_buf[mem->end + m]) {
                    numMatchingBP++;
                }
                else {
                    break;
                }
            }
            mem->end += numMatchingBP;
            rmemLen = mem->end - mem->start;
            next_be_point = mem->end;
            if (rmemLen >= raux->min_seed_len) {
                kv_push(mem_t, *smems, *mem);
            }
            else {
		        next_be_point += (raux->min_seed_len - rmemLen);    
            }           
        }
        else { // we don't have min_seed_len match for this start position
            assert(rmemLen <= raux->min_seed_len);	        
            next_be_point += (raux->min_seed_len - rmemLen);    
        }
    }
    else {
	    assert(lmemLen <= raux->min_seed_len);
        next_be_point += (raux->min_seed_len - lmemLen);       
    }

    return next_be_point;
}

/*
 * Compute final SMEMs and their hits after considering their overlaps
 *
 * Leaf expansion is also performed to get the actual length of the MEM. 
 * Note that leaf nodes store compressed suffixes by including a pointer to the reference genome
 *
 * @param iaux          index related parameters
 * @param raux          read related parameters
 * @param mem           maximal-exact-match storage
 * @param sh            helper data structure to keep track of start and end positions of previously identified MEMs
 * @param smems         list of SMEMs
 * @param hits          list of hits for read
 */
void check_and_add_smem(index_aux_t* iaux, read_aux_t* raux, mem_t* mem, smem_helper_t* sh, mem_v* smems, u64v* hits) {

    mem->start = raux->l_seq - mem->rc_end; // Adjust start position of LMEM
    int lmemLen = mem->end - mem->start;
    if (mem->hitcount > 0 && !mem->skip_ref_fetch) {
        int64_t len;
        int64_t start_ref_pos = hits->a[mem->hitbeg] + lmemLen;
        int64_t end_ref_pos = start_ref_pos + mem->start;
        // Fetch reference to check for extra matching bps
        uint8_t* rseq = bns_get_seq_v2(iaux->bns->l_pac, iaux->pac, start_ref_pos, end_ref_pos, &len, 0); 
        int m, numMatchingBP = 0;
        for (m = 0; m < len; ++m) {
            if (rseq[m] == raux->read_buf[mem->rc_end + m]) {
                numMatchingBP++;
            }
            else {
                break;
            }
        }
        // Adjust start position of MEM by extra matching bps
        mem->start -= numMatchingBP;
    }
    lmemLen = mem->end - mem->start;
    if (lmemLen >= raux->min_seed_len) {
        // Check if MEM lies completely within previously discovered MEM
        if (mem->start < sh->prevMemStart || mem->end > sh->prevMemEnd) {
            // Extra work for equivalency with BWA-MEM. 
            if (mem->fetch_leaves) {
                hits->n -= mem->hitcount;
                mem->hitbeg = hits->n;
                mem->hitcount = 0;
                raux->read_buf = raux->unpacked_queue_buf;
                rightExtend_fetch_leaves(iaux, raux, mem, hits);
                raux->read_buf = raux->unpacked_rc_queue_buf;
            }
            if (mem->hitcount > 0) {
                mem->pt.c_pivot = sh->curr_pivot;
                mem->pt.p_pivot = sh->prev_pivot;
                mem->pt.pp_pivot = sh->prev_prev_pivot;
                kv_push(mem_t, *smems, *mem);
                if (mem->start <= (sh->prev_pivot + 1)) {
                    sh->stop_be = 1;
                }
            }
            sh->prevMemStart = mem->start;
            sh->prevMemEnd = mem->end;
        }
    }
}

/*
 * This function replaces bwt_smem1() and uses ERT to generate SMEMs
 *
 * @param iaux          index related parameters
 * @param raux          read related parameters
 * @param smems         list of SMEMs
 * @param hits          list of hits for read
 */
void get_seeds_prefix(index_aux_t* iaux, read_aux_t* raux, mem_v* smems, u64v* hits) {

    smem_helper_t sh;
    memset(&sh, 0, sizeof(smem_helper_t));
    sh.prevMemStart = raux->l_seq;
    sh.prevMemEnd = 0;
    int i = 0, j = 0;
    sh.prev_pivot = -1;
    sh.prev_prev_pivot = -1;
    memset(raux->lep, 0, 5 * sizeof(uint64_t));
    while (i < raux->l_seq) { // Begin identifying RMEMs
        mem_t rm;
        memset(&rm, 0, sizeof(mem_t));
        rm.start = i; 
        rm.forward = 1;
        rm.hitbeg = hits->n;
        sh.curr_pivot = rm.start;
        raux->read_buf = raux->unpacked_queue_buf;
        rightExtend(iaux, raux, &i, &rm, hits); //!< Compute LEP.  
        // Lazy expansion of leaf nodes. 
        if (rm.hitcount > 0 && !rm.skip_ref_fetch) {
            int64_t len;
            int64_t start_ref_pos = hits->a[rm.hitbeg] + i - rm.start;
            int64_t end_ref_pos = hits->a[rm.hitbeg] + raux->l_seq - rm.start;
            // Fetch reference
            uint8_t* rseq = bns_get_seq_v2(iaux->bns->l_pac, iaux->pac, start_ref_pos, end_ref_pos, &len, 0); 
            int m;
            int numMatchingBP = 0;
            // Check for matching bases
            for (m = 0; m < len; ++m) {
                if (rseq[m] == raux->unpacked_queue_buf[i+m]) {
                    numMatchingBP++;
                }
                else {
                    raux->lep[(i+m-1) >> 6] |= (1ULL << ((i+m-1) & (0x3FULL)));
                    break;
                }
            }
            // Last base of RMEM must have LEP bit set
            if (m == len) {
                raux->lep[(i+m-1) >> 6] |= (1ULL << ((i+m-1) & (0x3FULL)));
            }
            i += numMatchingBP;
        }
        rm.end = i;
        int rmemLen = rm.end - rm.start;
        // No left-extension for position 0 in read
        // rm.start is the current pivot
        if (rm.start == 0) {
            if (rmemLen >= raux->min_seed_len) {
                if (rm.hitcount > 0) { 
                    kv_push(mem_t, *smems, rm);
                }
            }
            else {
                hits->n -= rm.hitcount;
            }
            memset(raux->lep, 0, 5 * sizeof(uint64_t));
        }
        else { // perform all backward extensions
            hits->n -= rm.hitcount;
            uint64_t* lep = raux->lep;
            int seq_len = raux->l_seq;
            int min_seed_len = raux->min_seed_len;
            sh.stop_be = 0; 
            int min_j = (rm.start > min_seed_len) ? (rm.start-1) : (min_seed_len-1);
            int max_j = rm.end - 1;
            j = min_j;
	        sh.prev_pivot = rm.start;
            while (j <= max_j) {
                mem_t m;
                int be_point;
                int mem_valid = init_mem(lep, &m, j, seq_len, min_seed_len);
                m.hitbeg = hits->n;
                int next_j = j + 1;
                if (mem_valid) {
                    be_point = j + 1;
                    if (be_point >= min_seed_len) {
                        int rc_i = seq_len - be_point; 
                        raux->read_buf = raux->unpacked_rc_queue_buf;
                        leftExtend(iaux, raux, &rc_i, &m, hits);
                        next_j = check_and_add_smem_prefix(iaux, raux, &m, &sh, smems, hits);
                    }
                }
                j = next_j;
                if (m.end > i) {
                    i = m.end;
                }
            }
        }
        raux->read_buf = raux->unpacked_queue_buf;
        // Skip all ambiguous bases
        while (i < raux->l_seq) {
            if (raux->read_buf[i] == 4) {
                ++i;
            }
            else {
                break;
            }
        }
        // Check if there other ambiguous bases within min_seed_len bases of the start of the MEM
        while ((i < raux->l_seq) && (i - rm.start) < raux->min_seed_len) {
            if (raux->read_buf[i] == 4) {
                ++i;
                break;
            }
            ++i;
        }
        sh.prev_prev_pivot = sh.prev_pivot;
        sh.prev_pivot = rm.start;
        memset(raux->lep, 0, 5 * sizeof(uint64_t));
    }
    #ifdef PRINT_SMEM
        ks_introsort(mem_smem_sort_lt_ert, smems->n, smems->a); // Sort SMEMs based on start pos in read. For DEBUG. 
        for (i = 0; i < smems->n; ++i) {
            // printf("[SMEM]:%d,%d\n", smems->a[i].start, smems->a[i].end);
            int idx;
            for (idx = 0; idx < smems->a[i].hitcount; ++idx) {
                if (smems->a[i].forward || smems->a[i].fetch_leaves) {
                    printf("[SMEM]:%d,%d,%lu\n", smems->a[i].start, smems->a[i].end, hits->a[smems->a[i].hitbeg + idx]);
                }
                else {
                    printf("[SMEM]:%d,%d,%lu\n", smems->a[i].start, smems->a[i].end, (iaux->bns->l_pac << 1) - hits->a[smems->a[i].hitbeg + idx] - (smems->a[i].end - smems->a[i].start - smems->a[i].end_correction));
                }
            }
        }
    #endif
}

/*
 * This function replaces bwt_smem1() and uses ERT to generate SMEMs
 *
 * @param iaux          index related parameters
 * @param raux          read related parameters
 * @param smems         list of SMEMs
 * @param hits          list of hits for read
 */
void get_seeds(index_aux_t* iaux, read_aux_t* raux, mem_v* smems, u64v* hits) {

    smem_helper_t sh;
    memset(&sh, 0, sizeof(smem_helper_t));
    sh.prevMemStart = raux->l_seq;
    sh.prevMemEnd = 0;
    int i = 0, j = 0;
    sh.prev_pivot = -1;
    sh.prev_prev_pivot = -1;
    memset(raux->lep, 0, 5 * sizeof(uint64_t));
    while (i < raux->l_seq) { // Begin identifying RMEMs
        mem_t rm;
        memset(&rm, 0, sizeof(mem_t));
        rm.start = i; 
        rm.forward = 1;
        rm.hitbeg = hits->n;
        sh.curr_pivot = rm.start;
        raux->read_buf = raux->unpacked_queue_buf;
        rightExtend(iaux, raux, &i, &rm, hits); // Compute LEP.  
        // Lazy expansion of leaf nodes. 
        if (rm.hitcount > 0 && !rm.skip_ref_fetch) {
            int64_t len;
            int64_t start_ref_pos = hits->a[rm.hitbeg] + i - rm.start;
            int64_t end_ref_pos = hits->a[rm.hitbeg] + raux->l_seq - rm.start;
            // Fetch reference
            uint8_t* rseq = bns_get_seq_v2(iaux->bns->l_pac, iaux->pac, start_ref_pos, end_ref_pos, &len, 0); 
            int m;
            int numMatchingBP = 0;
            // Check for matching bases
            for (m = 0; m < len; ++m) {
                if (rseq[m] == raux->unpacked_queue_buf[i+m]) {
                    numMatchingBP++;
                }
                else {
                    raux->lep[(i+m-1) >> 6] |= (1ULL << ((i+m-1) & (0x3FULL)));
                    break;
                }
            }
            // Last base of RMEM must have LEP bit set
            if (m == len) {
                raux->lep[(i+m-1) >> 6] |= (1ULL << ((i+m-1) & (0x3FULL)));
            }
            i += numMatchingBP;
        }
        rm.end = i;
        int rmemLen = rm.end - rm.start;
        // No left-extension for position 0 in read
        // rm.start is the current pivot
        if (rm.start == 0) {
            if (rmemLen >= raux->min_seed_len) {
                if (rm.hitcount > 0) { 
                    rm.pt.c_pivot = sh.curr_pivot;
                    rm.pt.p_pivot = sh.prev_pivot;
                    rm.pt.pp_pivot = sh.prev_prev_pivot; 
                    kv_push(mem_t, *smems, rm);
                }
            }
            else {
                hits->n -= rm.hitcount;
            }
            memset(raux->lep, 0, 5 * sizeof(uint64_t));
        }
        else {
            hits->n -= rm.hitcount;
            uint64_t* lep = raux->lep;
            int seq_len = raux->l_seq;
            int min_seed_len = raux->min_seed_len;
            j = rm.end-1;
            sh.stop_be = 0; 
            int min_j = (rm.start > min_seed_len) ? (rm.start-1) : (min_seed_len-1);
            while (j >= min_j) {
                mem_t m;
                int be_point;
                int mem_valid = init_mem(lep, &m, j, seq_len, min_seed_len);
                m.hitbeg = hits->n;
                if (mem_valid) {
                    be_point = j + 1;
                    if (be_point >= min_seed_len) {
                        int rc_i = seq_len - be_point; 
                        raux->read_buf = raux->unpacked_rc_queue_buf;
                        leftExtend(iaux, raux, &rc_i, &m, hits);
                        check_and_add_smem(iaux, raux, &m, &sh, smems, hits);
                        if (sh.stop_be) break;
                    }
                }
                j -= 1;
            }
        }
        raux->read_buf = raux->unpacked_queue_buf;
        // Skip all ambiguous bases
        while (i < raux->l_seq) {
            if (raux->read_buf[i] == 4) {
                ++i;
            }
            else {
                break;
            }
        }
        // Check if there other ambiguous bases within min_seed_len bases of the start of the MEM
        while ((i < raux->l_seq) && (i - rm.start) < raux->min_seed_len) {
            if (raux->read_buf[i] == 4) {
                ++i;
                break;
            }
            ++i;
        }
        sh.prev_prev_pivot = sh.prev_pivot;
        sh.prev_pivot = rm.start;
        memset(raux->lep, 0, 5 * sizeof(uint64_t));
    }
    #ifdef PRINT_SMEM
        ks_introsort(mem_smem_sort_lt_ert, smems->n, smems->a); // Sort SMEMs based on start pos in read. For DEBUG. 
        for (i = 0; i < smems->n; ++i) {
            // printf("[SMEM]:%d,%d\n", smems->a[i].start, smems->a[i].end);
            int idx;
            for (idx = 0; idx < smems->a[i].hitcount; ++idx) {
                if (smems->a[i].forward || smems->a[i].fetch_leaves) {
                    printf("[SMEM]:%d,%d,%lu\n", smems->a[i].start, smems->a[i].end, hits->a[smems->a[i].hitbeg + idx]);
                }
                else {
                    printf("[SMEM]:%d,%d,%lu\n", smems->a[i].start, smems->a[i].end, (iaux->bns->l_pac << 1) - hits->a[smems->a[i].hitbeg + idx] - (smems->a[i].end - smems->a[i].start));
                }
            }
        }
    #endif
}

/*
 * This function performs reseeding of SMEMs
 *
 * @param iaux          index related parameters
 * @param raux          read related parameters
 * @param smems         list of SMEMs
 * @param start         pivot position in read
 * @param limit         hit threshold below which tree traversal must stop
 * @param pt            track pivot information to reduce work done during reseeding
 * @param hits          list of hits for every read
 */
void reseed_prefix(index_aux_t* iaux, read_aux_t* raux, mem_v* smems, int start, int limit, pivot_t* pt, u64v* hits) {

    smem_helper_t sh;
    memset(&sh, 0, sizeof(smem_helper_t));
    sh.prevMemStart = raux->l_seq;
    sh.prevMemEnd = 0;
    int i = start, j = 0;
    int old_n = smems->n;
    memset(raux->lep, 0, 5 * sizeof(uint64_t));
    mem_t rm;
    memset(&rm, 0, sizeof(mem_t));
    rm.start = i;
    rm.forward = 1;
    rm.hitbeg = hits->n;
    sh.prev_pivot = (rm.start >= pt->c_pivot) ? pt->p_pivot : pt->pp_pivot;
    raux->read_buf = raux->unpacked_queue_buf;
    raux->limit = limit;
    rightExtend_wlimit(iaux, raux, &i, &rm, hits); // Compute LEP.  
    // Lazy expansion of leaf nodes. 
    if (rm.hitcount > 0 && !rm.skip_ref_fetch) {
        int64_t len;
        int64_t start_ref_pos = hits->a[rm.hitbeg] + i - rm.start;
        int64_t end_ref_pos = hits->a[rm.hitbeg] + raux->l_seq - rm.start;
        // Fetch reference
        uint8_t* rseq = bns_get_seq_v2(iaux->bns->l_pac, iaux->pac, start_ref_pos, end_ref_pos, &len, 0);
        int m;
        int numMatchingBP = 0;
        // Check for matching bases
        for (m = 0; m < len; ++m) {
            if (rseq[m] == raux->unpacked_queue_buf[i+m]) {
                numMatchingBP++;
            }
            else {
                raux->lep[(i+m-1) >> 6] |= (1ULL << ((i+m-1) & (0x3FULL)));
                break;
            }
        }
        // Last base of RMEM must have LEP bit set
        if (m == len) {
            raux->lep[(i+m-1) >> 6] |= (1ULL << ((i+m-1) & (0x3FULL)));
        }
        i += numMatchingBP;
    }
    rm.end = i;
    int rmemLen = rm.end - rm.start;
    if (rm.start == 0) {
        if (rmemLen >= raux->min_seed_len) {
            if (rm.hitcount > 0) { 
                kv_push(mem_t, *smems, rm);
            }
        }
        else {
            hits->n -= rm.hitcount;
        }
        memset(raux->lep, 0, 5 * sizeof(uint64_t));
    }
    // Begin left-extension, i.e., right extension on reverse complemented read
    else {
        hits->n -= rm.hitcount;
        uint64_t* lep = raux->lep;
        int seq_len = raux->l_seq;
        int min_seed_len = raux->min_seed_len;
        sh.stop_be = 0;
        int min_j = (rm.start > min_seed_len) ? (rm.start-1) : (min_seed_len-1);
        int max_j = rm.end - 1;
        j = min_j;
        sh.prev_pivot = rm.start;
        sh.mem_end_limit = rm.end;
	    while (j <= max_j) {
            mem_t m;
            int be_point;
            int mem_valid = init_mem(lep, &m, j, seq_len, min_seed_len);
            m.hitbeg = hits->n;
            int next_j = j + 1;
            if (mem_valid) {
                be_point = j + 1;
                if (be_point >= min_seed_len) {
                    int rc_i = seq_len - be_point; 
                    raux->read_buf = raux->unpacked_rc_queue_buf;
                    leftExtend_wlimit(iaux, raux, &rc_i, &m, hits);
                    next_j = check_and_add_smem_prefix_reseed(iaux, raux, &m, &sh, smems, hits);
                }
            }
            j = next_j;
        }
    }
    #ifdef PRINT_SMEM
        ks_introsort(mem_smem_sort_lt_ert, smems->n - old_n, &smems->a[old_n]); // Debug: Sort SMEMs based on start pos in read. 
        for (i = old_n; i < smems->n; ++i) {
            int idx;
            for (idx = 0; idx < smems->a[i].hitcount; ++idx) {
                if (smems->a[i].forward || smems->a[i].fetch_leaves) {
                    printf("[Reseed]:%d,%d,%lu\n", smems->a[i].start, smems->a[i].end, hits->a[smems->a[i].hitbeg + idx]);
                }
                else {
                    printf("[Reseed]:%d,%d,%lu\n", smems->a[i].start, smems->a[i].end, (iaux->bns->l_pac << 1) - hits->a[smems->a[i].hitbeg + idx] - (smems->a[i].end - smems->a[i].start - smems->a[i].end_correction));
                }
            }
        }
    #endif
}

/*
 * This function performs reseeding of SMEMs
 *
 * @param iaux          index related parameters
 * @param raux          read related parameters
 * @param smems         list of SMEMs and their hits
 * @param start         pivot position in read
 * @param limit         hit threshold below which tree traversal must stop
 * @param pt            track pivot information to reduce work done during reseeding
 * @param hits          list of hits for read
 */
void reseed(index_aux_t* iaux, read_aux_t* raux, mem_v* smems, int start, int limit, pivot_t* pt, u64v* hits) {

    smem_helper_t sh;
    memset(&sh, 0, sizeof(smem_helper_t));
    sh.prevMemStart = raux->l_seq;
    sh.prevMemEnd = 0;
    int i = start, j = 0;
    int old_n = smems->n;
    memset(raux->lep, 0, 5 * sizeof(uint64_t));
    mem_t rm;
    memset(&rm, 0, sizeof(mem_t));
    rm.start = i;
    rm.forward = 1;
    rm.hitbeg = hits->n;
    sh.prev_pivot = (rm.start >= pt->c_pivot) ? pt->p_pivot : pt->pp_pivot;
    raux->read_buf = raux->unpacked_queue_buf;
    raux->limit = limit;
    rightExtend_wlimit(iaux, raux, &i, &rm, hits); //!< Compute LEP.  
    // Lazy expansion of leaf nodes. 
    if (rm.hitcount > 0 && !rm.skip_ref_fetch) {
        int64_t len;
        int64_t start_ref_pos = hits->a[rm.hitbeg] + i - rm.start;
        int64_t end_ref_pos = hits->a[rm.hitbeg] + raux->l_seq - rm.start;
        // Fetch reference
        uint8_t* rseq = bns_get_seq_v2(iaux->bns->l_pac, iaux->pac, start_ref_pos, end_ref_pos, &len, 0);
        int m;
        int numMatchingBP = 0;
        // Check for matching bases
        for (m = 0; m < len; ++m) {
            if (rseq[m] == raux->unpacked_queue_buf[i+m]) {
                numMatchingBP++;
            }
            else {
                raux->lep[(i+m-1) >> 6] |= (1ULL << ((i+m-1) & (0x3FULL)));
                break;
            }
        }
        // Last base of RMEM must have LEP bit set
        if (m == len) {
            raux->lep[(i+m-1) >> 6] |= (1ULL << ((i+m-1) & (0x3FULL)));
        }
        i += numMatchingBP;
    }
    rm.end = i;
    int rmemLen = rm.end - rm.start;
    if (rm.start == 0) {
        if (rmemLen >= raux->min_seed_len) {
            if (rm.hitcount > 0) { 
                kv_push(mem_t, *smems, rm);
            }
        }
        else {
            hits->n -= rm.hitcount;
        }
        memset(raux->lep, 0, 5 * sizeof(uint64_t));
    }
    // Begin left-extension, i.e., right extension on reverse complemented read
    else {
        hits->n -= rm.hitcount;
        uint64_t* lep = raux->lep;
        int seq_len = raux->l_seq;
        int min_seed_len = raux->min_seed_len;
        j = rm.end-1;
        sh.stop_be = 0;
        int min_j = (rm.start > min_seed_len) ? (rm.start-1) : (min_seed_len-1);
        while (j >= min_j) {
            mem_t m;
            int be_point;
            int mem_valid = init_mem(lep, &m, j, seq_len, min_seed_len);
            m.hitbeg = hits->n;
            if (mem_valid) {
                be_point = j + 1;
                if (be_point >= min_seed_len) {
                    int rc_i = seq_len - be_point; 
                    raux->read_buf = raux->unpacked_rc_queue_buf;
                    leftExtend_wlimit(iaux, raux, &rc_i, &m, hits);
                    check_and_add_smem(iaux, raux, &m, &sh, smems, hits);
                    if (sh.stop_be) break;
                }
            }
            j -= 1;
        }
    }
    #ifdef PRINT_SMEM
        ks_introsort(mem_smem_sort_lt_ert, smems->n - old_n, &smems->a[old_n]); // Debug: Sort SMEMs based on start pos in read. 
        for (i = old_n; i < smems->n; ++i) {
            int idx;
            for (idx = 0; idx < smems->a[i].hitcount; ++idx) {
                if (smems->a[i].forward || smems->a[i].fetch_leaves) {
                    printf("[Reseed]:%d,%d,%lu\n", smems->a[i].start, smems->a[i].end, hits->a[smems->a[i].hitbeg + idx]);
                }
                else {
                    printf("[Reseed]:%d,%d,%lu\n", smems->a[i].start, smems->a[i].end, (iaux->bns->l_pac << 1) - hits->a[smems->a[i].hitbeg + idx] - (smems->a[i].end - smems->a[i].start));
                }
            }
        }
    #endif
}

/*
 * This function performs the LAST heuristic
 *
 * @param iaux          index related parameters
 * @param raux          read related parameters
 * @param smems         list of SMEMs
 * @param limit         hit threshold above which tree traversal must stop
 * @param hits          list of hits for read
 */
void last(index_aux_t* iaux, read_aux_t* raux, mem_v* smems, int limit, u64v* hits) {

    int i = 0;
    int old_n = smems->n;
    const uint8_t minSeedLen = raux->min_seed_len + 1; // LAST exits seeding only when seed length >= 20
    raux->limit = limit;
    while (i < raux->l_seq) { // Begin identifying RMEMs
        mem_t rm;
        rm.start = i; rm.end = 0;
        rm.forward = 1;
        rm.skip_ref_fetch = 0;
        rm.hitbeg = hits->n;
        rm.hitcount = 0;
        raux->read_buf = raux->unpacked_queue_buf;
        rightExtend_last(iaux, raux, &i, &rm, hits);  
        // Lazy expansion of leaf nodes. 
        if (rm.hitcount > 0 && !rm.skip_ref_fetch) {
            int64_t len;
            int64_t start_ref_pos = hits->a[rm.hitbeg] + i - rm.start;
            int64_t end_ref_pos = hits->a[rm.hitbeg] + raux->l_seq - rm.start;
            // Fetch reference
            uint8_t* rseq = bns_get_seq_v2(iaux->bns->l_pac, iaux->pac, start_ref_pos, end_ref_pos, &len, 0); 
            int m, numMatchingBP = 0;
            // Check for matching bases
            for (m = 0; m < len; ++m) {
                int seedLen = (i + m) - rm.start;
                int match_next_bp = ((seedLen < minSeedLen) || (rm.hitcount >= raux->limit)) ? 1 : 0;
                if (match_next_bp) {
                    if ((rseq[m] == raux->unpacked_queue_buf[i+m])) {
                        numMatchingBP++;
                    }
                    else {
                        ++i; // Increment i on every mismatch for LAST to match BWA-MEM
                        hits->n -= rm.hitcount;
                        rm.hitcount = 0;
                        break;
                    }
                }
                else {
                    break;
                }
            }
            i += numMatchingBP;
        }
        rm.end = i;
        int rmemLen = rm.end - rm.start;
        if (rmemLen >= minSeedLen) {
            if (rm.hitcount > 0 && rm.hitcount < raux->limit) { 
                kv_push(mem_t, *smems, rm);
            }
            else {
                hits->n -= rm.hitcount;
            }
        }
        else {
            hits->n -= rm.hitcount;
        }
        int foundN = 0;
        // Skip all ambiguous bases
        assert(i > 0);
        if (raux->read_buf[i-1] == 4) { // Found an 'N' during lookup
            foundN = 1;
        }
        if (!foundN) {
            // Check if there other ambiguous bases within min_seed_len bases of the start of the MEM
            while ((i < raux->l_seq) && (i - rm.start) < minSeedLen) {
                if (raux->read_buf[i] == 4) {
                    ++i;
                    break;
                }
                ++i;
            }
        }
    }
    #ifdef PRINT_SMEM
        ks_introsort(mem_smem_sort_lt_ert, smems->n - old_n, &smems->a[old_n]); // Debug: Sort SMEMs based on start pos in read. 
        for (i = old_n; i < smems->n; ++i) {
            int idx;
            for (idx = 0; idx < smems->a[i].hitcount; ++idx) {
                printf("[LAST]:%d,%d,%lu\n", smems->a[i].start, smems->a[i].end, hits->a[smems->a[i].hitbeg + idx]);
            } 
        }
    #endif
}

