#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "bwa.h"
#include "bwamem.h"

#define MD_MIN_QUALITY 15

/*! @abstract the read is paired in sequencing, no matter whether it is mapped in a pair */
#define BAM_FPAIRED        1
/*! @abstract the read is mapped in a proper pair */
#define BAM_FPROPER_PAIR   2
/*! @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
#define BAM_FUNMAP         4
/*! @abstract the mate is unmapped */
#define BAM_FMUNMAP        8
/*! @abstract the read is mapped to the reverse strand */
#define BAM_FREVERSE      16
/*! @abstract the mate is mapped to the reverse strand */
#define BAM_FMREVERSE     32
/*! @abstract this is read1 */
#define BAM_FREAD1        64
/*! @abstract this is read2 */
#define BAM_FREAD2       128
/*! @abstract not primary alignment */
#define BAM_FSECONDARY   256
/*! @abstract QC failure */
#define BAM_FQCFAIL      512
/*! @abstract optical or PCR duplicate */
#define BAM_FDUP        1024
/*! @abstract supplementary alignment */
#define BAM_FSUPPLEMENTARY 2048

//use to encode quality like samtools
//from htslib/hts.c
/*
const unsigned char seq_nt16_table2[256] = {
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
     1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,

    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};
*/


typedef struct readInfo {
    
    size_t pos;                    	 /**< clipped read position (coordinate in sam file) */
    size_t isize;			
    size_t mpos;			 /**< mate coordinate >**/
    size_t unclippedPos;           	 /**< unclipped read position */
    uint32_t dist2mate;			 /**< the distance to the next read>**/ 	
    
    unsigned int flag;                   /**< flag of the read */
    unsigned int pair_num;               /**< tell if the read is first(1) or second in the pair(2) */    
    unsigned int orientation;
    uint32_t mapq;			 /**<the mapping quality>**/	
    uint32_t mmapq;			 /**< MQ:Z tag (the quality of the mate) >**/     

    int tid;			         /**< chromosom id >*/
    int mtid;				 /**< chromosom mate id >*/
    int l_seq;				 /**< read length >*/
    uint32_t score;                	 /**< read phred score */
    uint32_t mscore;             	 /**< paired-end phred score */
     
    char *name;                          /**< read name, TODO:we don't need to store Qname, we use it to compte physical location for optical duplicates */
    char *cigar;                         /**< cigar string, TODO:we don't need to store cigar, we use it to compute unclipped coordinate */
    char *mcigar;                        /**< mate cigar >**/
    char *aux;				 /** < auxillary string that contains NM: XA: ZA:... > **/        
    char *seq;			         /** < nucleotids bases >**/
    uint8_t *qual;			         /*  < nucleotids bases quality>*/

} readInfo;

int fixmate( int rank, bseq1_t *seqs_1, bseq1_t *seqs_2, int *read_num_1, int *read_num_2, bwaidx_t *indix );
int readParsing (char *sam_buff, readInfo *read, bwaidx_t *indix);
int getTokenTab(char **offset, char **tokenTab);
int sync_mq_mc(readInfo *src, readInfo *dest);
int sync_mate(readInfo *a, readInfo *b);
int add_mate_score(readInfo *src, readInfo *dest);
int sam_write_unmapped_and_munmapped(readInfo *read, char **final_buffer, bwaidx_t *indix);
int sam_write_supp_and_secondary(readInfo *read, char **final_buffer, bwaidx_t *indix);
int sam_write_mate_unmapped(readInfo *read, char **final_buffer, bwaidx_t *indix);
int sam_write(readInfo *read, char **final_buffer, bwaidx_t *indix);



