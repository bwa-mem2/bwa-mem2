


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


void init_goff(size_t *goff, MPI_File mpi_filed, size_t fsize,int numproc,int rank);
void find_process_starting_offset(size_t *goff, size_t size, char* file_to_read, int proc_num, int rank_num);
void find_process_starting_offset_mt(size_t *goff, size_t size, char* file_to_read, int proc_num, int rank_num, int nthreads);


void find_reads_size_and_offsets(size_t offset_in_file,
								size_t siz2read,
								char *file_to_read,
								size_t *p_local_num_reads,
								size_t *p_total_num_reads,
								size_t **local_read_offsets,
								int **local_read_size,
								size_t **local_read_bytes,
								int proc_num, int rank_num);
void find_chunks_info_trim(	size_t *begin_offset_chunk,
				size_t *begin_offset_chunk_2,
				size_t *chunk_size,
				size_t *chunk_size_2,
				size_t *reads_in_chunk,
				size_t *reads_in_chunk_2,
				int *local_read_size,
				int *local_read_size_2,
				size_t *local_read_bytes,
				size_t *local_read_bytes_2,
				size_t *local_read_offsets,
				size_t *local_read_offsets_2,
				int rank_num,
				int proc_num,
				size_t local_num_reads,
				size_t local_num_reads_2,
				size_t grand_total_num_reads,
				size_t grand_total_num_reads_2,
				off_t maxsiz,
				size_t *chunk_count);

void find_chunks_info(  size_t *begin_offset_chunk,
			size_t *chunk_size,
			size_t *reads_in_chunk,
			int *local_read_size,
			size_t *local_read_bytes,
			size_t *local_read_offsets,
			int rank_num,
			int proc_num,
			size_t local_num_reads,
			size_t grand_total_num_reads,
			off_t maxsiz,
			size_t *chunk_count);

void map_indexes(char *file_map, ktp_aux_t *aux, MPI_Win *win_shr_ref, MPI_Win *win_shr_pac);
void create_sam_header(char *file_out, bwaidx_fm_t *indix, int *count, char *hdr_line, char *rg_line, char *pg_line, int rank_num);

int getChr(char *str, char *chrNames[], int nbchr, char *tmp_chr);
void create_sam_header_by_chr_file(char *file_out[], 
		bwaidx_t *indix, int *count, char *hdr_line, char *rg_line, int rank_num);

void *compress_and_write_thread(void *threadarg);
void *compress_and_write_bam_thread(void *threadarg);
void *compress_and_write_bgzf_thread(void *threadarg);
void *compress_thread_by_chr(void *threadarg);
void *compress_thread_by_chr_single(void *threadarg);
void *call_fixmate(void *threadarg);
void *copy_local_read_info_mt(void *thread_arg);
void *find_reads_size_and_offsets_mt(void *thread_arg);
void *copy_buffer_write_thr(void *thread_arg);
void *compute_buffer_size_thr(void *thread_arg);
void create_bam_header(char *file_out, bwaidx_fm_t *indix, int *count, 
    char *hdr_line, char *rg_line, char *pg_line, int rank_num, int compression_level);
void create_bam_header_by_chr_file(char *file_out[], bwaidx_t *indix, int *count, 
    char *hdr_line, char *rg_line, char *pg_line, int rank_num, int compression_level, int dofixmate);
void *write_sam_mt(void *thread_arg);



void *pread_fastq_chunck(void *thread_arg );


/****************************************
 *  
 *       Thread Structures and functions
 *  
 ****************************************/

struct thread_data
{
        int  job_rank;
        int  thread_id;
        int  total_thread;
        int  total_lines;
        size_t total_reads;
        int start_index_seqs;
        int final_index_seqs;
        bseq1_t *seqs_thr;
        bwaidx_t *indix_thr;
};

struct thread_data_compress
{
        int  job_rank;
        int  thread_id;
        int  total_thread;
        size_t total_reads;
        MPI_File fh_out;
        bseq1_t *seqs_thr;
        uint8_t *compressed_buffer_thread;
        int comp_level;
        size_t thr_comp_sz;
};

struct thread_data_compress_by_chr
{
        int  job_rank;
        int  thread_id;
        int  total_thread;
	bwaidx_t *indix_thr;
        size_t total_sam_line;
        int comp_level;
        char **start_addr;
	int *add_in_disc;
        int *line_size_to_cpy;
        int *sam_buff_dest;
        int incrmnt;
        MPI_File *fh_out;
};

struct thread_data_compress_by_chr_single
{
        int  job_rank;
        int  thread_id;
        int  total_thread;
        bwaidx_t *indix_thr;
        size_t total_sam_line;
        int comp_level;
        char **start_addr;
        int *line_size_to_cpy;
        int *sam_buff_dest;
        int incrmnt;
        MPI_File *fh_out;
};


//struct used to find started offset of rank and threads in fastq 
/*struct struct_data_thread{
    bseq1_t *seqs_thr;
    int begin_index;
    int end_index;
    int size_thr;
    char *buffer_out;
    MPI_File file_desc;
    size_t offset_write;
};
*/

struct struct_data_thread{
    bseq1_t *seqs_thr;
    int begin_index;
    int end_index;
    MPI_File file_desc;
};


//struct used to find read size, offset, and bytes in multithread
struct struct_data_thread_1{
            size_t offset_in_file_mt;
            size_t size2read_mt;
            char   *file_r1_mt;
            size_t *local_num_reads_mt;
            size_t *total_num_reads_mt;
            size_t **local_read_offsets_mt;
            int    **local_read_size_mt;
            size_t **local_read_bytes_mt;
            size_t *local_read_offsets;
            int    *local_read_size;
            size_t *local_read_bytes;
            size_t previous_read_num;
            int proc_num_mt;
            int rank_num_mt;
            int thread_num_mt;
};

struct struct_pread_fastq{

            int total_thread;
            int thread_id;
            int job_rank;
            size_t offset;
            size_t size;
            char *buffer;
            MPI_File fd;
};

