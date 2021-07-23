
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <assert.h>
#include <ctype.h>
#include "fixmate.h"
#include "bwa.h"
#include "bwamem.h"

#define min(a,b) (a>=b?b:a)

#define SWAP_UINT16(x) (((x) >> 8) | ((x) << 8))
#define SWAP_UINT32(x) (((x) >> 24) | (((x) & 0x00FF0000) >> 8) | (((x) & 0x0000FF00) << 8) | ((x) << 24))


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


unsigned char seq_nt4_table2[256] = {
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};




char *flag2string(int flagValue) {
    char *valueString = calloc(5, sizeof(char));
    sprintf(valueString, "%u", flagValue);
    return valueString;
}

int getTokenTab(char **offset, char **tokenTab) {
      
    
    char *lineHead = *offset, *next;
    next = strchr(*offset, '\t');
    
    if (next) {
    	*offset = next;
        int counter = (*offset - lineHead) / sizeof(char);
        (*offset)++;
        *tokenTab = strndup(lineHead, counter);
        return counter;
    } else {
	char *lineEnd = strchr(*offset, '\n');
    	size_t remainsTextSize = lineEnd ? (lineEnd - *offset) / sizeof(char) : strlen(*offset);
    	*offset += remainsTextSize;
    	*tokenTab = remainsTextSize ? strndup(lineHead, remainsTextSize) : NULL;
    	return remainsTextSize;
    }
}

inline unsigned int readBits(unsigned int x, int k) {
    // from K&R page 49
    return ((x & ( 1 << k )) >> k);
}

inline void writeBits(unsigned int *x, int k) {
    *x |= (1 << k);
}

inline void toggleBits(unsigned int *x, int k) {
    *x ^= readBits(*x, k);
}

inline void switchBits(unsigned int *x, int k, int j){
    toggleBits(x, k);
    toggleBits(x, j);
}

inline int fastqToPhred (char ch) {
    int ich = ch;
    assert(ich > 33 || ich < 126);
    return (ich - 33);
}

size_t getReferenceLength(char *cigar) {
    int refLen = 0, len = 0;
    char *q = cigar;

    while (*q) {
        len = 0;

        /* Compute the length of bases which the operator applies */
        while (isdigit(*q)) {
            len = 10 * len + (*q - '0');
            q++;
        }

        /* Cigar operator is one of M/D/N/X/EQ */
        if (*q && (*q == 'M' || *q == 'D' || *q == 'N' || *q == 'X' || *q == '=')) {
            refLen += len;
        }

        q++;
    }

    return refLen;
}


/*
hts_pos_t bam_endpos(const bam1_t *b)
{
    hts_pos_t rlen = (b->core.flag & BAM_FUNMAP)? 0 : bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
    if (rlen == 0) rlen = 1;
    return b->core.pos + rlen;
}

hts_pos_t bam_cigar2rlen(int n_cigar, const uint32_t *cigar)
{
    int k;
    hts_pos_t l;
    for (k = l = 0; k < n_cigar; ++k)
        if (bam_cigar_type(bam_cigar_op(cigar[k]))&2)
            l += bam_cigar_oplen(cigar[k]);
    return l;
}


*/



int readParsing (char *sam_buff, readInfo *read, bwaidx_t *indix) {

    int rank, num_proc;
    char *q = sam_buff;
    char *tokenCar;
    char *u;
    unsigned int readUnmapped;
    unsigned int mateUnmapped;
    unsigned int secondaryAlignment;
    unsigned int supplementaryAlignment;
    unsigned int firstInPair;
    int i = 0;
    int j = 0;
    uint32_t score;
       
    while (1) {
        switch (i++) {

            // These fields are ignored
            
            break;
            case 0: //qname we already have it
		        q = strchr(q, '\t') + 1;	
                break;
            case 1: // Flag part
		        getTokenTab(&q, &tokenCar);
                read->flag  = atoi(tokenCar);
                assert(read->flag);
                free(tokenCar);
                break;

            case 2: // RNAME part

		        getTokenTab(&q, &tokenCar);
                for ( j = 0; j < indix->bns->n_seqs; ++j) {
                        if (strcmp(indix->bns->anns[j].name, tokenCar) == 0) {
                            read->tid = j;
                            break;
                        }
			            else read->tid = -1;
                }
                free(tokenCar);
                break;

	    case 3: // The position chromosome wise
		        getTokenTab(&q, &tokenCar);
		        if (strcmp(tokenCar, "*") == 0) read->pos = -1;
                else read->pos = atoll(tokenCar);
                free(tokenCar);
                break;

	    case 4: // The MAPQ
		
		        getTokenTab(&q, &tokenCar);
                read->mapq = atoi(tokenCar);
                free(tokenCar);
                break;

	    case 5: // the CIGAR string
	            //q = strchr(q, '\t') + 1;
		        getTokenTab(&q, &tokenCar);
                read->cigar = strdup(tokenCar);
                assert(read->cigar);
		        free(tokenCar);
		        break;

   	    case 6: // the mate chromosom

		        getTokenTab(&q, &tokenCar);

                /* mate chromosome same as read chromosome */
                if (strcmp(tokenCar, "=") == 0) {
                    read->mtid = read->tid;

                } else {
                    for ( j = 0; j < indix->bns->n_seqs; ++j) {
                        if (strcmp(indix->bns->anns[j].name, tokenCar) == 0) {
                            read->mtid = j;
                            break;
                        }
			else read->mtid = -1;
                    }
                }
                free(tokenCar);
                break;

	   case 7: // here we have the PNEXT, the position of the next read	
		        getTokenTab(&q, &tokenCar);
		        if (strcmp(tokenCar, "*") == 0) read->mpos = -1;
                else read->mpos = atol(tokenCar);
                free(tokenCar);
                break;

	   case 8: // the Mate MAPQ 		
		        getTokenTab(&q, &tokenCar);
                read->dist2mate = atol(tokenCar);
                free(tokenCar);
                break;

	   case 9: // the segment SEQ
                //q = strchr(q, '\t') + 1;
                //break; 
		        getTokenTab(&q, &tokenCar);
		        read->seq = strdup(tokenCar);
		        free(tokenCar); 
	
                break;

           case 10:

                getTokenTab(&q, &tokenCar);
		read->qual=strndup(tokenCar, strlen(read->seq));	
		
	        read->score=0;
	      	uint8_t *u = read->qual;
        	int i;
		
    		for (i = 0; i < strlen(read->seq); i++) {
        	    if ((u[i]-33) >= MD_MIN_QUALITY) read->score += (u[i]-33);
    		}

                q = strchr(q, '\t') + 1;	
                free(tokenCar);
                break;

	   case 11:
		/*  after quality string we have the 
 		 *  auxillary tags we copy it in the 
 		 *  read->aux
 		 */
                if ( *q == '\t' ) q++;
                char *start = q;
                char *end = strchr(q, '\n');
                int res = 0;
		        //replace  \n 
		        //*end='\n';		
                res = asprintf(&(read->aux), start, end - start);
                assert(res > 0);
                break;

	   default:
                break;
        }	
        if (i == 12) break;	
    }

return 0;
}

int getLine(char **offset, char **tokenLine) {

	char *lineHead = *offset, *next;
        next = strchr(*offset, '\n');
        if (next) {
	        *offset = next;
                int counter = (*offset - lineHead) / sizeof(char);
                (*offset)++;
                *tokenLine = strndup(lineHead, counter + 1);
                return counter;
         } else {
                *tokenLine = NULL;
                return 0;
         }
}

void sync_unmapped_pos_inner(readInfo *src, readInfo *dest) {
	if ((dest->flag & BAM_FUNMAP) && !(src->flag & BAM_FUNMAP)) {
        // Set unmapped read's RNAME and POS to those of its mapped mate
        // (recommended best practice, ensures if coord sort will be together)
        	dest->tid = src->tid;
                dest->pos = src->pos;
        }
}

void sync_mate_inner(readInfo *src, readInfo *dest)
{
    // sync mate pos information
    dest->mtid = src->tid; dest->mpos = src->pos;
    // sync flag info
    if (src->flag&BAM_FREVERSE)
    	dest->flag |= BAM_FMREVERSE;
    else
        dest->flag &= ~BAM_FMREVERSE;
    if (src->flag & BAM_FUNMAP) {
    	dest->flag |= BAM_FMUNMAP;
    }
}


// Returns 0 on success, -1 on failure.
int sync_mq_mc(readInfo *src, readInfo *dest)
{
	if ( (src->flag & BAM_FUNMAP) == 0 ) { // If mapped
        	// Copy Mapping Quality
        	// MQ tag in dest
	        dest->mmapq = src->mapq;
                // Copy mate cigar if either read is mapped
		if ( (src->flag & BAM_FUNMAP) == 0 || (dest->flag & BAM_FUNMAP) == 0 ) {
	        	//size_t clen = strlen( src->cigar);
			dest->mcigar = strdup(src->cigar);	
		}
      		return 0;
	}
	return 0;
}
int sync_mate(readInfo *a, readInfo *b)
{
    sync_unmapped_pos_inner(a,b);
    sync_unmapped_pos_inner(b,a);
    sync_mate_inner(a,b);
    sync_mate_inner(b,a);
    if (sync_mq_mc(a,b) < 0) return -1;
    if (sync_mq_mc(b,a) < 0) return -1;
    return 0;
}

int add_mate_score(readInfo *src, readInfo *dest)
{
    dest->mscore = src->score;
    return 0;
}

int sam_write_supp_and_secondary(readInfo *read, char **final_buffer, bwaidx_t *indix){

        char *current_line;
        int res = 0;
	    res = asprintf(&current_line, "%s\t%d\t%s\t%zu\t%d\t%s\t%s\t%lu\t%d\t%s\t%s\t%s", read->name,
                   	read->flag, indix->bns->anns[read->tid].name, read->pos, read->mapq, 
				        read->cigar, indix->bns->anns[read->mtid].name, read->mpos, read->dist2mate,
                            read->seq,  read->qual, read->aux);

        assert( res > 0);

        free(read->aux);
        free(read->cigar);
        free(read->seq);
	free(read->qual);
        size_t len1 = 0;
        if (*final_buffer) len1=strlen(*final_buffer);
        size_t len2 = strlen(current_line);
        if (len1){
        	*final_buffer=realloc(*final_buffer, (len1 + len2 + 1));
                assert(*final_buffer);
                (*final_buffer)[len1 + len2]=0;
                memcpy(*final_buffer + len1, current_line, len2);
        }
        else{ 
            free(*final_buffer);
            *final_buffer=strdup(current_line);
        }
        free(current_line);
        return 0;
  }
  
int sam_write_discordant(readInfo *read, char **final_buffer, bwaidx_t *indix){

    	char *current_line;
    	int res = 0;
	res = asprintf(&current_line, "%s\t%d\t%s\t%zu\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s", read->name,
                                read->flag, indix->bns->anns[read->tid].name, read->pos, read->mapq,
                                        read->cigar, "*", 0, 0, read->seq,  read->qual, read->aux);

    	assert(res > 0);
    	free(read->aux);
    	free(read->cigar);
    	free(read->seq);
	free(read->qual);
	size_t len1 = 0;
	if (*final_buffer) len1=strlen(*final_buffer);
        size_t len2 = strlen(current_line);
        if (len1){
                *final_buffer=realloc(*final_buffer, (len1 + len2 + 1));
                assert(*final_buffer);
                (*final_buffer)[len1 + len2]=0;
                memcpy(*final_buffer + len1, current_line, len2);
        }
	else{
        free(*final_buffer); 
        *final_buffer=strdup(current_line);
        }
        free(current_line);
        return 0;
}



int sam_write_unmapped_and_munmapped(readInfo *read, char **final_buffer, bwaidx_t *indix){

	char *current_line;
    char *mchr=NULL;
	char *chr=NULL;
    int res = 0;

	if (read->tid == -1 && read->mtid == -1) {chr = "*"; mchr = "*";}
 	else chr=indix->bns->anns[read->tid].name;

        if ( (read->tid != -1)  && (read->tid == read->mtid )) {
			chr=indix->bns->anns[read->tid].name;
			mchr="=";
			}
        if (read->mtid != -1 ) mchr = indix->bns->anns[read->mtid].name;
	
	//uint8_t mscore_t = 0;
        //memcpy(&mscore_t, (uint8_t *)&read->mscore, sizeof(uint32_t));

	//read->mscore = bswap_32(read->mscore);
	SWAP_UINT32(read->mscore);
    	assert(chr);
	assert(mchr);
        res = asprintf(&current_line, "%s\t%d\t%s\t%zu\t%d\t%s\t%s\t%zu\t%d\t%s\t%s\tms:i:%u\t%s", read->name,
                                read->flag, chr, read->pos, read->mapq, read->cigar, mchr, read->mpos, read->dist2mate,
                                                read->seq,  read->qual, read->mscore, read->aux);

        assert( res > 0);
        free(read->aux);
        free(read->cigar);
        free(read->seq);
	free(read->qual);
        size_t len1 = 0;
        if (*final_buffer) len1=strlen(*final_buffer);
        size_t len2 = strlen(current_line);

        if (len1){
                *final_buffer=realloc(*final_buffer, (len1 + len2 + 1));
                assert(*final_buffer);
		        (*final_buffer)[len1 + len2]=0;
                memcpy(*final_buffer + len1, current_line, len2);
        }
	else {
        free(*final_buffer);
        *final_buffer=strdup(current_line);
    }
	free(current_line);
	return 0;
}
int sam_write_mate_unmapped(readInfo *read, char **final_buffer, bwaidx_t *indix){
	
	char *current_line;
    int res = 0;
    char *mchr;
    if ( read->tid == read->mtid ) mchr="=";
    else mchr = indix->bns->anns[read->mtid].name;

	//uint8_t mscore_t = 0;
        //memcpy(&mscore_t, (uint8_t *)&read->mscore, sizeof(uint32_t));

	//read->mscore = bswap_32(read->mscore);
	SWAP_UINT32(read->mscore);
	if ( read->flag&BAM_FUNMAP)
        res = asprintf(&current_line, "%s\t%d\t%s\t%zu\t%d\t%s\t%s\t%zu\t%d\t%s\t%s\tMQ:i:%d\tMC:Z:%s\tms:i:%u\t%s", read->name,
                        read->flag, indix->bns->anns[read->tid].name, read->pos,
                            read->mapq, read->cigar, mchr, read->mpos, read->dist2mate,
                                read->seq,  read->qual, read->mmapq, read->mcigar, read->mscore, read->aux);

	else 
		res = asprintf(&current_line, "%s\t%d\t%s\t%zu\t%d\t%s\t%s\t%zu\t%d\t%s\t%s\tMC:Z:*\tms:i:%u\t%s", read->name,
                                read->flag, indix->bns->anns[read->tid].name, read->pos,
                                        read->mapq, read->cigar, mchr, read->mpos, read->dist2mate,
                                                read->seq,  read->qual,  read->mscore, read->aux);
 


        assert( res > 0 );
        free(read->aux);
        free(read->cigar);
        free(read->seq);
	free(read->qual);
        assert(current_line);
        size_t len1 = 0;
        if (*final_buffer) len1=strlen(*final_buffer);
        size_t len2 = strlen(current_line);

        if (len1){
                *final_buffer=realloc(*final_buffer, len1 + len2 + 1);
                assert(*final_buffer);
                (*final_buffer)[len1 + len2]=0;
                memcpy(*final_buffer + len1, current_line, len2);
        }
        else {
            free(*final_buffer);
            *final_buffer=strdup(current_line);
        }
        free(current_line);
        return 0;
}

int sam_write(readInfo *read, char **final_buffer, bwaidx_t *indix){

	//we use asprintf to write in the final buffer
	//we write in order
	//Qname, Flag
	char *current_line;
    	int res = 0;
    	char *mchr;
	if ( read->tid == read->mtid ) mchr="=";
	else mchr = indix->bns->anns[read->mtid].name;

	//uint8_t mscore_t = 0;
	//memcpy(&mscore_t, (uint8_t *)&read->mscore, sizeof(uint32_t));

	//read->mscore = bswap_32(read->mscore);
	SWAP_UINT32(read->mscore);
	res = asprintf(&current_line, "%s\t%d\t%s\t%zu\t%d\t%s\t%s\t%zu\t%d\t%s\t%s\tMQ:i:%d\tms:i:%u\t%s", read->name, 
				read->flag, indix->bns->anns[read->tid].name, read->pos, 
					read->mapq, read->cigar, mchr, read->mpos, read->dist2mate, 
						read->seq,  read->qual, read->mmapq, read->mscore, read->aux);


    	assert( res > 0 );
    	free(read->aux);
    	free(read->cigar);
    	free(read->seq);
        free(read->qual);    
	assert(current_line);
	size_t len1 = strlen(*final_buffer);
	size_t len2 = strlen(current_line);

    // fprintf(stderr, "len1 = %d\n", len1);

	if (len1){
        //free(*final_buffer);
        //*final_buffer=malloc((len1 + len2 + 1)*sizeof(char)); 
		*final_buffer=realloc(*final_buffer, len1 + len2 + 1);
		assert(*final_buffer);
		(*final_buffer)[len1 + len2 ]=0;
		memcpy(*final_buffer + len1 , current_line, len2);
	}
	else {
        free(*final_buffer); 
        *final_buffer=strdup(current_line);
        }
	free(current_line);
	return 0;
}



int fixmate( int rank, bseq1_t *seqs_1, bseq1_t *seqs_2, int *reads_1, int *reads_2, bwaidx_t *indix ){
		/*  First part of the algorithme */

		/* number of ref chr : indix.bns->n_seqs
		   name of chromosoms : indix.bns->anns[s].name
		   
		   seqs contains:   int l_seq, id;
		   		    char *name, *comment, *seq, *qual, *sam;
		*/
        int res = 0;
        char *seqs_fxmt1=NULL;
        char *seqs_fxmt2=NULL;
        seqs_fxmt1 = malloc(sizeof(char)*1);
        seqs_fxmt2 = malloc(sizeof(char)*1);
        //seqs_fxmt1[0]=(char*)malloc(sizeof(char));
        //seqs_fxmt2[0]=(char*)malloc(sizeof(char));
        seqs_fxmt1[0]=0;
        seqs_fxmt2[0]=0;


		int d = 0;
		char *m = seqs_1->sam;
		int read_num_1 = 0;
		while (d < strlen(seqs_1->sam)) {
			if (*m++ == '\n') read_num_1++;
			d++;
		}
		
		d = 0;
		m = seqs_2->sam;
        int read_num_2 = 0;
        while ( d < strlen(seqs_2->sam)  ) {
			if (*m++ == '\n') read_num_2++;
			d++;
		}
		*reads_1 = read_num_1;
		*reads_2 = read_num_2;	
		//we allocate a read_info
		//read_num_1 = read_num_2 = 10;
		readInfo **reads = calloc ( (read_num_1 + read_num_2), sizeof(readInfo *));
		int i = 0;
		//we parse the sam buffer original
		//fprintf(stderr, " read_num1 : %d ## read_num2 : %d \n", read_num_1, read_num_2);
		// fprintf(stderr, " start analysing  seqs1=%s \n seqs2=%s\n", seqs_1->sam, seqs_2->sam);
		char *q =seqs_1->sam; 
		while (*q) {
			char *tok;
        		// parse read
       	    		getLine(&q, &tok);
			//we fill first info
			reads[i] = malloc (sizeof(readInfo));
			reads[i]->name = seqs_1->name;
			//reads[i]->qual = seqs_1->qual;
			reads[i]->l_seq = seqs_1->l_seq;
			//now we get the rest
            		assert (readParsing(tok, reads[i], indix) == 0);
            		free(tok);
			i++;
		}
		
		char *q2 = seqs_2->sam;
        	while (*q2) {
            		char *tok2;
            		// parse read
            		getLine(&q2, &tok2);
			reads[i] = malloc (sizeof(readInfo));
            		//we fill first info
            		reads[i]->name = seqs_2->name;
            		//reads[i]->qual = seqs_2->qual;
			reads[i]->l_seq = seqs_2->l_seq;
            		//reads[i]->seq = seqs_2->seq;
            		assert (readParsing(tok2, reads[i], indix) == 0);
            		free(tok2);
            		i++;
       		}
              
		/*
		if ( i != ( read_num_1 + read_num_2)){ 
		    fprintf(stderr, " problem finish parsing line 1= %d :: 2 = %d :: i =%d \n", read_num_1, read_num_2, i);
		    fprintf(stderr, " seqs1= %s \n seqs2 = %s \n ", seqs_1->sam, seqs_2->sam);	
		}
		assert(i == ( read_num_1 + read_num_2));
		*/
		int curr = 0; int has_prev = 0;
		size_t pre_end, cur_end;
		
		//read from the first pair (not supplementary or secondary)
		readInfo *read1;
		readInfo *read2;
		int have_pair1 = 0; //check we have a pair of BAM_FPAIRED reads	
		int pair1_idx1 = 0;
		int pair1_idx2 = 0;	
		int have_pair2 = 0; //check we have a pair of !BAM_FPAIRED reads
		int pair2_idx1 = 0;
		int pair2_idx2 = 0;

		int both_mapped=0;
		int both_unmapped=0;
		int supplemtary=0;
		
		//now we care for the unmapped and mate unmapped		
		for  ( i=0; i< ( read_num_1 + read_num_2); i++ ){

			//we check both pair unmapped	
                        if ( (reads[i]) && (reads[i]->flag & BAM_FUNMAP) && (reads[i]->flag & BAM_FMUNMAP)){
                             
				  if  (reads[i]->flag & BAM_FREAD1) {read1 = reads[i]; have_pair1++; pair1_idx1 = i; }
                                  if  (reads[i]->flag & BAM_FREAD2) {read2 = reads[i]; have_pair1++; pair1_idx2 = i; }
				
				   if (have_pair1 == 2){
                                                read1->flag |= BAM_FPAIRED;
                                                read2->flag |= BAM_FPAIRED;
                                                assert(add_mate_score(read1, read2) == 0);
                                                assert(add_mate_score(read2, read1) == 0);
                                                assert ( sam_write_unmapped_and_munmapped(read1, &seqs_fxmt1, indix) == 0 );
                                                assert ( sam_write_unmapped_and_munmapped(read2, &seqs_fxmt2, indix) == 0 );
                                                free(reads[pair1_idx1]);
                                                free(reads[pair1_idx2]);
                                                reads[pair1_idx1] = NULL;
                                                reads[pair1_idx2] = NULL;
                                                have_pair1=0;

                                        }
                        }
		}
		have_pair1=0;	

		for  ( i=0; i< ( read_num_1 + read_num_2); i++ ){

			if ( (reads[i]) && (reads[i]->flag&BAM_FPAIRED) && !(reads[i]->flag&BAM_FSECONDARY)
                                        && !(reads[i]->flag&BAM_FSUPPLEMENTARY) && !(reads[i]->flag&BAM_FUNMAP) && !(reads[i]->flag&BAM_FMUNMAP)){

                                if ((reads[i]) && (reads[i]->flag & BAM_FPAIRED)){

                                        if  (reads[i]->flag & BAM_FREAD1) {read1 = reads[i]; have_pair1++; pair1_idx1 = i; }
                                        if  (reads[i]->flag & BAM_FREAD2) {read2 = reads[i]; have_pair1++; pair1_idx2 = i; }

                                        if (have_pair1 == 2){

                                                if ( (read1->tid != read1->mtid) && (read2->tid != read2->mtid)) {
                                                        read1->flag &= ~BAM_FPAIRED;
                                                        read2->flag &= ~BAM_FPAIRED;
                                                        read1->mtid = -1;
                                                        read2->mtid = -1;
                                                        read1->mpos = 0;
                                                        read2->mpos = 0;
                                                        read1->dist2mate = 0;
                                                        read2->dist2mate = 0;
                                                        assert ( sam_write_discordant( read1, &seqs_fxmt1, indix) == 0 );
                                                        assert ( sam_write_discordant( read2, &seqs_fxmt2, indix) == 0 );
                                                }
                                                else {
                                                        read1->flag |= BAM_FPAIRED;
                                                        read2->flag |= BAM_FPAIRED;
                                                        sync_mate( read1, read2);
                                                        assert(add_mate_score(read1, read2) == 0);
                                                        assert(add_mate_score(read2, read1) == 0);
                                                        assert ( sam_write(read1, &seqs_fxmt1, indix) == 0 );
                                                        assert ( sam_write(read2, &seqs_fxmt2, indix) == 0 );
                                                }
                                                free(reads[pair1_idx1]);
                                                free(reads[pair1_idx2]);
                                                reads[pair1_idx1] = NULL;
                                                reads[pair1_idx2] = NULL;

                                                have_pair1=0;

                                        }
                                }
                        }
		}

		for  ( i=0; i< ( read_num_1 + read_num_2); i++ ){

                        if ( (reads[i]) && ((reads[i]->flag & BAM_FSECONDARY) || (reads[i]->flag & BAM_FSUPPLEMENTARY))){
                                if  (reads[i]->flag & BAM_FREAD1)
                                        assert ( sam_write_supp_and_secondary( reads[i], &seqs_fxmt1, indix) == 0 );
                                else
                                        assert ( sam_write_supp_and_secondary( reads[i], &seqs_fxmt2, indix) == 0 );
                         
                                free(reads[i]);
                                reads[i]=NULL;
                                supplemtary++;
                        }
			if ((reads[i]) && (reads[i]->flag&BAM_FPAIRED)){
                                        if  (reads[i]->flag&BAM_FREAD1) {read1 = reads[i]; have_pair2++; pair2_idx1 = i;}
                                        if  (reads[i]->flag&BAM_FREAD2) {read2 = reads[i]; have_pair2++; pair2_idx2 = i;}
                                        if (have_pair2 == 2){
                                                read1->flag |= BAM_FPAIRED;
                                                read2->flag |= BAM_FPAIRED;
                                                sync_mate( read1, read2);
                                                assert(add_mate_score(read1, read2) == 0);
                                                assert(add_mate_score(read2, read1) == 0);
                                                assert ( sam_write_mate_unmapped(read1, &seqs_fxmt1, indix) == 0 );
                                                assert ( sam_write_mate_unmapped(read2, &seqs_fxmt2, indix) == 0 );
                                                have_pair2=0;
                                                free(reads[pair2_idx1]);
                                                free(reads[pair2_idx2]);
                                                reads[pair2_idx1] = NULL;
                                                reads[pair2_idx2] = NULL;
                                        }

                                }

                }

		for  ( i=0; i< ( read_num_1 + read_num_2); i++ ) assert(reads[i] == NULL);
		
		free(seqs_1->sam);
		res = asprintf(&(seqs_1->sam),"%s",seqs_fxmt1);
        assert ( res > 0 );    
		free(seqs_2->sam);
        res = asprintf(&(seqs_2->sam),"%s",seqs_fxmt2);
        assert(res > 0);
		free(seqs_fxmt1);
        free(seqs_fxmt2);

		if (reads) free(reads);
		return 0;
}





                                                                                                                               
