/*
This file is part of mpiBWA

The project was developped by Frederic Jarlier from Institut Curie and Nicolas Joly from Institut Pasteur

NGS aligner inspired by BWA-MEM 

Copyright (C) 2016-2021  Institut Curie / Institut Pasteur

You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
   
The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND. Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data. 
   
The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.

*/




#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#include <mpi.h>
#include <fcntl.h>
#include <pthread.h>
#include <stdio.h>
#include <zlib.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <sys/stat.h>
#include "utils.h"
#include "fastmap.h"
#include "fixmate.h"
#include "tokenizer.h"
#include "parallel_aux.h"
#include "bgzf.h"
#include "bgzf.c"
#include "FMI_search.h"
#define DEFAULT_INBUF_SIZE (1024*1024*1024)
#define NB_PROC  "16" //numer of threads for writing
#define CB_NODES "2" //numer of server for writing
#define CB_BLOCK_SIZE  "268435456" /* 256 MBytes - should match FS block size */
#define CB_BUFFER_SIZE  "3758096384" /* multiple of the block size by the number of proc*/
#define DATA_SIEVING_READ "enable"
#define min(a,b) (a>=b?b:a)
#define MAX_CHAR_SIZE 2048
#define MAX_CHR_NAME_SIZE 200
#define SMALL_STACK (1024*1024)
#define BIG_STACK (1024*1024*512)

#ifdef TIMING
#define xfprintf fprintf
#else
#define xfprintf(...) /**/
#endif



void init_goff(size_t *goff, MPI_File mpi_filed, size_t fsize,int numproc,int rank){


	char * current_line = NULL;
	MPI_Status status;
	int i = 0;
	int j = 0;
	// TODO problem here when we have only one job!
	size_t lsize = fsize/numproc;
	goff[0]=0;
	for(i=1;i<numproc;i++){goff[i]=lsize*i;}
	goff[numproc]=fsize;
	for(i=1;i<numproc;i++){
		current_line =(char*)calloc(2000,sizeof(char));
		//we read 2K caracters
		MPI_File_read_at(mpi_filed, (MPI_Offset)goff[i], current_line, 2000, MPI_CHAR, &status);
		j=0;
		while (j < 2000 ){
			//we check the \n after and before
			j++;
			if (current_line[j] == '+' && current_line[j-1] == '\n') {j++;break;}}
		//then we look for return after quality
		while( j < 2000){j++;if (current_line[j] == '\n') {break;}}
		goff[i]+=(j+1); free(current_line);}

	return;
}

///Function used to find the starting reading offset in a given file according to the number of processes used.
//Parameters (ordered):	pointer on the vector that will save these offsets
//			the size of the given file
//			the studied file
//			the number or processes used
//			the rank of the actual process
//quick workflow : 	find the size to read per process
//			define the starting offset arbitrarily
//			browse a sample, from that offset, find the next beginning of a new read
//			update the starting offset in the vector
void find_process_starting_offset_mt(size_t *goff, size_t size, char* file_to_read, int proc_num, int rank_num, int nthreads)
{
	///MPI related resources
	int res; //used to assert success of MPI communications
	MPI_File mpi_fd; // file descriptor used to open and read from the right file
	MPI_Status status;
	res = MPI_File_open(MPI_COMM_SELF, file_to_read,  MPI_MODE_RDONLY , MPI_INFO_NULL, &mpi_fd); //open the wanted file
	assert(res==MPI_SUCCESS);
	
	///other resources
	off_t tmp_sz = 1024; //size of the sample from the file, big enough to contain a full read
    	char *buffer_r0 = malloc( tmp_sz + 1); //buffer used to save the sample
	buffer_r0[tmp_sz] = '\0'; 
	size_t lsize = size/(proc_num * nthreads); //proportion of the file 1 process should read
	int i; //used as an iterator
    	char *p, *q, *e; //pointers on the buffer ot find the start of a read
	
	///define the arbitrary offsets
	goff[0]=0;
	for(i = 1 ; i < (proc_num * nthreads) ; i++){ goff[i] = lsize*i;}
    	goff[proc_num * nthreads] = size;

	int k = 0;
   	size_t first_index = rank_num * nthreads;

   	while ( k < nthreads){

 	    res = MPI_File_read_at(mpi_fd,  (MPI_Offset)goff[first_index + k], buffer_r0, tmp_sz, MPI_CHAR, &status); //read the wanted part of the file nd save it into the buffer
    	assert(res == MPI_SUCCESS);
		p = buffer_r0;
    	e = buffer_r0 + tmp_sz;

		//browse the buffer to find the beginning of the next read

		while (p < e) {
            if (*p != '@') { p++; continue; }
            if (p != buffer_r0 && *(p-1) != '\n') { p++; continue; }
            q = p + 1;
            while (q < e && *q != '\n') q++; q++;
            while (q < e && *q != '\n') q++; q++;
            if (q < e && *q == '+') break;
            p++;
        }

		assert(*p == '@');
       	//we update begining offsets with the found value
       	goff[first_index + k] += p - buffer_r0;        
        memset( buffer_r0, 0, tmp_sz * sizeof(char));
       	k++;
	}
	fprintf(stderr, "finish finding offset rank %d \n", rank_num);
    ///free the resources no longer needed
    free(buffer_r0);
    res = MPI_File_close(&mpi_fd);
    assert(res == MPI_SUCCESS);

}

void *find_reads_size_and_offsets_mt(void *thread_arg){


    struct struct_data_thread_1 *my_data;
    my_data = (struct struct_data_thread_1 *) thread_arg;

    size_t offset_in_file           = my_data->offset_in_file_mt;
    size_t siz2read                 = my_data->size2read_mt;
    char   *file_to_read            = my_data->file_r1_mt;
    size_t *p_local_num_reads       = my_data->local_num_reads_mt;
    size_t *p_total_num_reads       = my_data->total_num_reads_mt;
    size_t **local_read_offsets     = my_data->local_read_offsets_mt;
    int    **local_read_size        = my_data->local_read_size_mt;
    size_t **local_read_bytes       = my_data->local_read_bytes_mt;
    int    proc_num                 = my_data->proc_num_mt;
    int    rank_num                 = my_data->rank_num_mt;
    int    thread_num               = my_data->thread_num_mt;

    MPI_File  mpi_fd;
    MPI_Status status;
    int count;
    //int fd;
    int res;
    res = MPI_File_open(MPI_COMM_WORLD, file_to_read, MPI_MODE_RDONLY, MPI_INFO_NULL, &mpi_fd);
    assert(res==MPI_SUCCESS);
    //fd = open(file_to_read, O_RDONLY);
    

    char *buffer_r;
    char *b, *r, *t, *e;
    size_t offset_end_buff;
    size_t pos_in_vect = 0;
    size_t lines = 0;
    size_t total_parsing = 0;
    size_t g=0;

    MPI_Datatype arraytype;
    MPI_Datatype arraytype0;
    MPI_Datatype arraytype_r1;
    MPI_Datatype arraytype_r2;

    *p_total_num_reads = 0;
    size_t read_buffer_sz = 0;
    if ( siz2read < DEFAULT_INBUF_SIZE ) read_buffer_sz = siz2read;
    else read_buffer_sz = DEFAULT_INBUF_SIZE;

    while (1){
        
        buffer_r = malloc(read_buffer_sz + 1);
        assert( buffer_r != NULL );
        buffer_r[read_buffer_sz] = '0';

        //pread(fd, buffer_r, read_buffer_sz, offset_in_file);
        res = MPI_File_read_at(mpi_fd, (MPI_Offset)offset_in_file, buffer_r, read_buffer_sz, MPI_CHAR, &status);
        assert(res == MPI_SUCCESS);
        res = MPI_Get_count(&status, MPI_CHAR, &count);
        assert(res == MPI_SUCCESS);
        assert(count == (int)read_buffer_sz && *buffer_r == '@');
	   

        b = buffer_r;
        r = b + read_buffer_sz;

        if ( read_buffer_sz == DEFAULT_INBUF_SIZE){
            while (r-- != b){if (*r == '\n' && *(r+1) == '+') {r--; break;}}
            while (r-- != b){if (*r == '\n') break;}
            while (r-- != b){if (*r == '@') break;}
            r--;
            offset_end_buff = (r - b);
        }
        else
            offset_end_buff = (r - b);

	    t = buffer_r ;
        e = buffer_r + offset_end_buff;
        lines = 0;
        while (t++ < e){if (*t == '\n') lines++;}


        *p_local_num_reads =  (lines/4);
        *p_total_num_reads += *p_local_num_reads;
        *local_read_size    = (int *)realloc(*local_read_size, sizeof(int) * (*p_total_num_reads));
        *local_read_offsets = (size_t *)realloc(*local_read_offsets, sizeof(size_t) * (*p_total_num_reads));
        *local_read_bytes   = (size_t *)realloc(*local_read_bytes, sizeof(size_t) * (*p_total_num_reads));

        assert( *local_read_offsets != NULL);
        assert( *local_read_size    != NULL);
        assert( *local_read_bytes   != NULL);

        t = buffer_r;

        int size=0;

        size_t lines2 = 0;
        size_t lines3 = 0;

        t = buffer_r;
        int last_size;
        int done = 0;
        size_t start_read_offset = 0;

        g = offset_in_file;

        int count = 0;

        while (t < e){

            if (lines3 == lines) break;
            assert( *t == '@');

            start_read_offset = g;
            while (*t != '\n'){ t++; g++;}
            t++;g++; lines3++;
            while (*t != '\n'){ size++; t++; g++;}
            t++;g++; lines3++;
            while (*t != '\n'){ t++; g++;}
            t++;g++; lines3++;
            while (*t != '\n'){ t++; g++;}
            lines3++;

            (*local_read_offsets)[pos_in_vect]  = start_read_offset;
            (*local_read_bytes)[pos_in_vect]    = (g - start_read_offset) + 1;
            assert((*local_read_bytes)[pos_in_vect] != 0);
            (*local_read_size)[pos_in_vect]     = size;
            assert((*local_read_size)[pos_in_vect] != 0);

            size = 0;
		    pos_in_vect++;
            t++;g++;

        }

        assert( lines == lines3 );
        total_parsing   += offset_end_buff;
        if (total_parsing == siz2read) {free(buffer_r); break;}
        if ((siz2read - total_parsing) < DEFAULT_INBUF_SIZE)
            read_buffer_sz = siz2read - total_parsing;
        else read_buffer_sz = DEFAULT_INBUF_SIZE;

        offset_in_file  += offset_end_buff + 1;
        free(buffer_r);
    }

    assert(total_parsing == siz2read);
    MPI_File_close(&mpi_fd);
    //close(fd);
}

void *copy_local_read_info_mt(void *thread_arg){
	struct struct_data_thread_1 *my_data;
    	my_data = (struct struct_data_thread_1 *) thread_arg;
        size_t *p_local_num_reads       = my_data->local_num_reads_mt;
        size_t *p_total_num_reads       = my_data->total_num_reads_mt;
        size_t **local_read_offsets     = my_data->local_read_offsets_mt;
        int    **local_read_size        = my_data->local_read_size_mt;
        size_t **local_read_bytes       = my_data->local_read_bytes_mt;
        int    proc_num                 = my_data->proc_num_mt;
        int    rank_num                 = my_data->rank_num_mt;
        int    thread_num               = my_data->thread_num_mt;
    
    
        size_t offset = my_data->previous_read_num;
    	int *p      = my_data->local_read_size + offset;
    	memmove(p, *local_read_size, *p_total_num_reads * sizeof(int));

    	size_t *p1  = my_data->local_read_offsets + offset;
    	memmove(p1, *local_read_offsets, *p_total_num_reads * sizeof(size_t));

    	size_t *p2  = my_data->local_read_bytes + offset;
    	memmove(p2, *local_read_bytes, *p_total_num_reads * sizeof(size_t));

}







///Function used to find the size and starting offset of each read in a given file
//Parameters (ordered):	the starting offset of the running process
//			the size of the file the process has to read
//			the given file
//			a pointer on the number of reads in the studied interval
//			a pointer on the number of reads in the file
//			a doubly differenciated pointer on the array used to store the offsets
//			a doubly differenciated pointer on the array used to store the sizes
//			the number of processes used
//			the rank of the running process
//quick worflow : 	read from the starting offset 1Go at a time
//			find the last read in that buffer and stop at its end
//			find the number of lines in contains and use it to realloc the size and offset vector
//			then browse the buffer to find the size and starting offset of each read
//			display stats about the time it took to do all the previous actions
void find_reads_size_and_offsets(size_t offset_in_file,
								size_t siz2read,
								char *file_to_read,
								size_t *p_local_num_reads,
								size_t *p_total_num_reads,
								size_t **local_read_offsets,
								int **local_read_size,
								size_t **local_read_bytes,
								int proc_num, int rank_num){
	
	///MPI related resources
	MPI_File  mpi_fd; //file handle for the given file
	int res; 
	res = MPI_File_open(MPI_COMM_WORLD, file_to_read, MPI_MODE_RDONLY, MPI_INFO_NULL, &mpi_fd);
	assert(res==MPI_SUCCESS);
	
	///Non MPI related resources
	char *buffer_r;
	char *b, *r, *t, *e; //used as iterators
	size_t offset_end_buff;
	size_t pos_in_vect = 0;
	size_t lines = 0;
	size_t total_parsing = 0;
	size_t g=0;
	size_t chunck_count =0;
        size_t chunck_number = 0;
        size_t i1 = 0;
        size_t i2 = 0;

	MPI_Datatype arraytype;
	MPI_Datatype arraytype0;	
	MPI_Datatype arraytype_r1;
	MPI_Datatype arraytype_r2;
	//init p_total_num_reads
	*p_total_num_reads = 0;
	//you can only read a certain size of file at a time iot reduce the processor load
	size_t read_buffer_sz = 0;
        if ( siz2read < DEFAULT_INBUF_SIZE ) {chunck_number = 1; read_buffer_sz = siz2read;}
        else {
                i1 = siz2read / (size_t)DEFAULT_INBUF_SIZE;
                i2 = siz2read % (size_t)DEFAULT_INBUF_SIZE;
                chunck_number = i1;
                if (i2 != 0) chunck_number++;
                read_buffer_sz = DEFAULT_INBUF_SIZE;
        }
	
	while (chunck_count < chunck_number){

		buffer_r = malloc(read_buffer_sz + 1);
		assert( buffer_r != NULL );
		buffer_r[read_buffer_sz] = '0';
		
		// FOR TEST //
		// MPI_Type_contiguous(read_buffer_sz, MPI_CHAR, &arraytype0);
		// MPI_Type_commit(&arraytype0);
		// MPI_File_set_view(mpi_fd_in2, (MPI_Offset)off_in_file, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL ) ; 
		// res = MPI_File_read(mpi_fd_in2, b , 1, arraytype0, &status);
		res = MPI_File_read_at(mpi_fd, (MPI_Offset)offset_in_file, buffer_r, read_buffer_sz, MPI_CHAR, MPI_STATUS_IGNORE);
		assert(res == MPI_SUCCESS);
		assert(*buffer_r == '@');	

		//we search the last read 
		b = buffer_r;
		r = b + read_buffer_sz;

		if ( read_buffer_sz == DEFAULT_INBUF_SIZE){
			r--;
			//go to last previous read
			while (r != (b+1)){if (*r == '+' && *(r-1) == '\n') break; else r--;}
            		while (r != b){if (*r == '\n') break; else r--;}
            		while (r != b){if (*r == '@') break; else r--;}
            		r--; //stop at \n
            		offset_end_buff = (r - b);
		}
		else
			offset_end_buff = (r - b);

		//we count the number of lines
		t = buffer_r ;
		e = buffer_r + offset_end_buff;
		lines = 0;		
		while (t++ < e){if (*t == '\n') lines++;}

		assert( lines%4 == 0);

		*p_local_num_reads =  (lines/4);
		*p_total_num_reads += *p_local_num_reads;
		//realloc the vectors according to the new number of reads you have
		//they are doubly differenciated because the realloc can change the vector pointer
		*local_read_size  	= (int *)realloc(*local_read_size, sizeof(int) * (*p_total_num_reads));
		*local_read_offsets	= (size_t *)realloc(*local_read_offsets, sizeof(size_t) * (*p_total_num_reads));
		*local_read_bytes	= (size_t *)realloc(*local_read_bytes, sizeof(size_t) * (*p_total_num_reads));

		assert( *local_read_offsets != NULL);
		assert( *local_read_size 	!= NULL);
		assert( *local_read_bytes 	!= NULL);
		
		//for the first first read
		t = buffer_r;
		
		int size=0;
		size_t lines2 = 0;
		size_t lines3 = 0;

		t = buffer_r;
		//we move to the next line because we test *(t-1)
		int last_size;
		int done = 0;
		size_t start_read_offset = 0;

		g = offset_in_file;
		//here we browse the buffer and each time we reach a new read we save:
		//- the size of the previous one
		//- the offset of the next one
		//the last_size resource is used in case we're treating the last read

		int count = 0;

		while (t < e){

			if (lines3 == lines) break;
			assert( *t == '@');

			//we get the read start offset
			start_read_offset = g;
			while (*t != '\n'){ t++; g++;} //the name
			t++;g++; lines3++;
			while (*t != '\n'){ size++; t++; g++;} //the read
			t++;g++; lines3++;//skip \n
			while (*t != '\n'){ t++; g++;} //the +
			t++;g++; lines3++;
			while (*t != '\n'){ t++; g++;} //the qual
			lines3++;

			(*local_read_offsets)[pos_in_vect]	= start_read_offset;
			(*local_read_bytes)[pos_in_vect]	= (g - start_read_offset) + 1;
			assert((*local_read_bytes)[pos_in_vect] != 0);
			(*local_read_size)[pos_in_vect] 	= size;
			assert((*local_read_size)[pos_in_vect] != 0);

			size = 0;
			pos_in_vect++;
			t++;g++; //we go to the next read

		}

		assert( lines == lines3 );
		chunck_count++;

                if ( (chunck_count + 1) < chunck_number) total_parsing   += offset_end_buff + 1;
                else total_parsing   += offset_end_buff;

                if ((siz2read - total_parsing) < DEFAULT_INBUF_SIZE)
                        read_buffer_sz = siz2read - total_parsing;
                else read_buffer_sz = DEFAULT_INBUF_SIZE;

                offset_in_file  += offset_end_buff + 1;
                free(buffer_r);
	}

	assert(total_parsing == siz2read);
	MPI_File_close(&mpi_fd);	

}

///Function used to evaluate the number, size, number of reads and starting offset of the chunks
//Parameters (ordered):	pointer on the vector that has the starting offset of each chunk for f1 
//			pointer on the vector that has the starting offset of each chunk for f2 
//			pointer on the vector that has the size of each chunk for f1 (in bytes)
//			pointer on the vector that has the size of each chunk for f2 (in bytes)
//			pointer on the vector that has the number of reads in each chunk for f1
//			pointer on the vector that has the number of reads in each chunk for f2
//			pointer on the vector that has the size of each read for f1
//			pointer on the vector that has the size of each read for f2
//			pointer on the vector that has the starting offset of each read for f1
//			pointer on the vector that has the starting offset of each read for f2
//			rank of the running process
//			number of processes used
//			number of reads in the size given to the running process
//			total number of reads in the file
//			maximum size of a chunk
//			number of chunks needed by the file
//quick workfllow : 	every process waits for the previous one to end him info except rank 0
//			rank 0 browses his reads and gets the chunk info
//			when he hits maxsiz, he send the current info to the next process
//			next processes does the same until every read is taken into account
//			then a last chunk is filled with the remaining bases
//			display some stats about the time it took to do all this
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
				size_t *chunk_count){
	
	//non MPI related resources
	int count 		= 0;
	int res;
	int i,bef,aft;
	int read_nb 		= 0; // have to use a buffer value instead of the iterator directly in case there are more reads handled by the process than there are in a chunk.
	int nb_reads_to_send	= 0;
	int nb_reads_to_recv	= 0;
	int reads_type_to_send	= 0; 	// use to know if we shall send forward (0) or backward (1)
	int reads_type_to_recv	= 0; 	// use to know if we shall recieve forward (0) or backward (1)
	int *sizes_to_send	= NULL;
	int *sizes_to_recv	= NULL;
	int *p_sizes  		= NULL; //reserved pointer
	int *p_sizes_2  	= NULL; //reserved pointer
	
	size_t *p_bytes  	= NULL; //reserved pointer
	size_t *p_bytes_2  	= NULL; //reserved pointer
	size_t u1 		= 0;
	size_t local_read_min;
	size_t *bytes_to_send	= NULL;
	size_t *bytes_to_recv	= NULL;
	size_t *p_offset 	= NULL; //reserved pointer
	size_t *p_offset_2 	= NULL; //reserved pointer
	size_t reads_recieved	= 0;
	size_t reads_recieved_2	= 0;
	size_t *offsets_to_send	= NULL;
	size_t *offsets_to_recv	= NULL;
	size_t counter_bases  	= 0;
	size_t bytes_in_chunk 	= 0, bytes_in_chunk_2 = 0;
	size_t index_in_chunk 	= 0;
	size_t begin_offset 		= 0, begin_offset_2	= 0;
	size_t bases_previous_chunk     = 0;
	size_t offset_previous_chunk 	= 0, offset_previous_chunk_2 = 0;
	size_t size_previous_chunk 	= 0, size_previous_chunk_2 = 0;
	size_t x 			=0,y=0; 	// iterators of the number of reads in each file
	size_t read1=0,read2=0;
	size_t size_chunk = 100;
	size_t tmp_cnt=0;

	MPI_Status status;

	if (rank_num > 0){
		// we wait the previous rank to send the final size chunk
		// and offset of the chunk
		MPI_Recv(&bases_previous_chunk,
				1,
				MPI_LONG_LONG_INT,
				rank_num - 1,
				0,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);

		MPI_Recv(&offset_previous_chunk,
				1,
				MPI_LONG_LONG_INT,
				rank_num - 1,
				1,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);

		MPI_Recv(&offset_previous_chunk_2,
				1,
				MPI_LONG_LONG_INT,
				rank_num - 1,
				2,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);

		MPI_Recv(&size_previous_chunk,
				1,
				MPI_LONG_LONG_INT,
				rank_num - 1,
				3,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);

		MPI_Recv(&size_previous_chunk_2,
				1,
				MPI_LONG_LONG_INT,
				rank_num - 1,
				4,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);

		//we send the type or reads we are going to send
		MPI_Recv(&reads_type_to_recv,
				1,
				MPI_INT,
				rank_num - 1,
				5,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);

		//we send the size of the vector for the recieving part
		MPI_Recv(&nb_reads_to_recv,
				1,
				MPI_LONG_LONG_INT,
				rank_num - 1,
				6,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);

		if ( reads_type_to_recv == 0){
			//we recieve for forward reads part
			reads_recieved 		= nb_reads_to_recv;
			reads_recieved_2 	= 0;
		}
		else{
			//we recieve for backward reads part
			reads_recieved 		= 0;
			reads_recieved_2 	= nb_reads_to_recv;
		}

		//here we malloc the size
		sizes_to_recv 		= realloc(sizes_to_recv, nb_reads_to_recv * sizeof(int));
		bytes_to_recv 		= realloc(bytes_to_recv, nb_reads_to_recv * sizeof(size_t));
		offsets_to_recv 	= realloc(offsets_to_recv, nb_reads_to_recv * sizeof(size_t));

		assert(sizes_to_recv != NULL);
		assert(bytes_to_recv != NULL);
		assert(offsets_to_recv != NULL);

		//we send the vector of offsets to the next rank
		MPI_Recv(offsets_to_recv,
				nb_reads_to_recv,
				MPI_LONG_LONG_INT,
				rank_num - 1,
				7,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);

		//we send the vector of sizes to the next rank
		MPI_Recv(sizes_to_recv,
				nb_reads_to_recv,
				MPI_INT,
				rank_num - 1,
				8,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);

		MPI_Recv(bytes_to_recv,
				nb_reads_to_recv,
				MPI_LONG_LONG_INT,
				rank_num - 1,
				9,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);

		// we send read1 and read2
		// number of reads already parsed
		MPI_Recv(&read1,
				1,
				MPI_LONG_LONG_INT,
				rank_num - 1,
				10,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);

		MPI_Recv(&read2,
				1,
				MPI_LONG_LONG_INT,
				rank_num - 1,
				11,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);

	}

	if (rank_num == 0){
		
		// here you find which file buffer of the rank 0 holds the least reads
		local_read_min = min(local_num_reads, local_num_reads_2);

		if(local_num_reads_2 == local_read_min) {
			nb_reads_to_send 	= (local_num_reads - local_read_min);
			reads_type_to_send 	= 0;
		}
		else {
			nb_reads_to_send 	= (local_num_reads_2 - local_read_min);
			reads_type_to_send	= 1;
		}
		for (u1 = 0; u1 < local_read_min; u1++){

			if (u1 >= 1){
				assert(*(local_read_offsets + x)   != *(local_read_offsets + x - 1));
				assert(*(local_read_offsets_2 + y) != *(local_read_offsets_2 + y - 1));
			}

			// we for offsets multiple of maxsize
			if (counter_bases == 0){
				begin_offset 	= *(local_read_offsets + x);
				begin_offset_2 	= *(local_read_offsets_2 + y);
			}

			assert( (*(local_read_size + x)) != 0 );
			assert( (*(local_read_bytes + x)) != 0 );

			counter_bases 	+= *(local_read_size + x);
			bytes_in_chunk  += *(local_read_bytes + x);
			x++;read1++;

			assert( (*(local_read_size_2 + y)) != 0 );
			assert( (*(local_read_bytes_2 + y)) != 0 );

			counter_bases	 += *(local_read_size_2 + y);
			bytes_in_chunk_2 += *(local_read_bytes_2 + y);
			y++;read2++;

			if ( counter_bases > maxsiz){

				assert( read1 == read2 );

				//the current chunk is full
				//then we have the number of bases wanted
				begin_offset_chunk[index_in_chunk] 		= begin_offset;
				begin_offset_chunk_2[index_in_chunk] 		= begin_offset_2;
				chunk_size[index_in_chunk]    			= bytes_in_chunk;
				chunk_size_2[index_in_chunk] 	   		= bytes_in_chunk_2;
				reads_in_chunk[index_in_chunk]	   		= read1;
				reads_in_chunk_2[index_in_chunk]   		= read2;
				*chunk_count 			   += 1;
				index_in_chunk 			   += 1;
				//we reset counter of bases
				counter_bases  		= 0;
				bytes_in_chunk		= 0;
				bytes_in_chunk_2	= 0;
				read_nb				= 0;
				read1				= 0;
				read2				= 0;
			}
		} //end for loop on local_read_min
	}//end if rank_num == 0
	else{

		//receive the info about the current chunk being filled
		counter_bases 			= bases_previous_chunk;
		bytes_in_chunk 			= size_previous_chunk;
		bytes_in_chunk_2 		= size_previous_chunk_2;
		begin_offset			= offset_previous_chunk;
		begin_offset_2			= offset_previous_chunk_2;
		///receive the nb of reads sent
		//malloc a vector to its size
		//fill it with the sent data
		//redo the calculus of the nb of reads to treat and to send
		//once you get that new number, do the rest of the treatments
		//once you're finished check that everything went smoothly and send the right nb ofreads to the next rank.
		x=0; y=0;
		size_t counter_tmp = 0;
		///here you find which file buffer of the rank 0 that holds the last reads

		//we test if for the last rank the reads are equal
		if (rank_num == (proc_num - 1))
			assert((local_num_reads + reads_recieved) == (local_num_reads_2 + reads_recieved_2));

		local_read_min = min( local_num_reads + reads_recieved, local_num_reads_2 + reads_recieved_2);

		if ( (local_num_reads_2 + reads_recieved_2) == local_read_min ) {
			nb_reads_to_send 	= (local_num_reads + reads_recieved) - local_read_min;
			reads_type_to_send 	= 0; //we send forward
		}
		else {
			nb_reads_to_send 	= (local_num_reads_2 + reads_recieved_2) - local_read_min;
			reads_type_to_send	= 1; //we send backward
		}

		// here we initialyze p_offset and p_sizes
		// according to recieve part
		if ( reads_type_to_recv == 0){
			// we recieve from forward reads part
			p_offset		= offsets_to_recv;
			p_sizes		= sizes_to_recv;
			p_bytes		= bytes_to_recv;

			p_offset_2	= local_read_offsets_2;
			p_sizes_2	= local_read_size_2;
			p_bytes_2	= local_read_bytes_2;
		}
		else{
			// we recieve from backward reads part
			p_offset	= local_read_offsets;
			p_sizes		= local_read_size;
			p_bytes		= local_read_bytes;

			p_offset_2	= offsets_to_recv;
			p_sizes_2	= sizes_to_recv;
			p_bytes_2	= bytes_to_recv;
		}

		for (u1 = 0; u1 < local_read_min; u1++){

			if (u1 >= 1){
				assert(*(p_offset + x) != *(p_offset + x - 1));
				assert(*(p_offset_2 + y) != *(p_offset_2 + y - 1));
			}
			// when we reached the number of reads recieved
			// we switch of vector and pass to the local vector
			if (counter_tmp == nb_reads_to_recv){

				if ( reads_type_to_recv == 0){
					p_offset	= local_read_offsets;
					p_sizes		= local_read_size;
					p_bytes		= local_read_bytes;
					x = 0; //only pos in forward buff are reset
				}
				else{
					p_offset_2	= local_read_offsets_2;
					p_sizes_2	= local_read_size_2;
					p_bytes_2	= local_read_bytes_2;
					y = 0; //only pos in backward buff are reset
				}
			} //end condition switch



			//we for offsets multiple of maxsize
			if (counter_bases == 0){
				begin_offset 	= *(p_offset + x);
				begin_offset_2 	= *(p_offset_2 + y);
			}

			counter_bases 	+= *(p_sizes + x);
			bytes_in_chunk  += *(p_bytes + x);
			x++; read1++;

			counter_bases	 += *(p_sizes_2 + y);
			bytes_in_chunk_2 += *(p_bytes_2 + y);
			y++; read2++;
			counter_tmp++;

			if ( counter_bases > maxsiz){

				assert( read1 == read2 );

				//the current chunk is full
				//then we have the number of bases wanted
				begin_offset_chunk[index_in_chunk] 		= begin_offset;
				begin_offset_chunk_2[index_in_chunk] 		= begin_offset_2;
				chunk_size[index_in_chunk]    			= bytes_in_chunk;
				chunk_size_2[index_in_chunk] 	   		= bytes_in_chunk_2;

				reads_in_chunk[index_in_chunk]	   		= read1;
				reads_in_chunk_2[index_in_chunk]   		= read2;
				*chunk_count		+= 1;
				index_in_chunk 		+= 1;
				//we reset counter of bases
				counter_bases  		= 0;
				bytes_in_chunk 		= 0;
				bytes_in_chunk_2	= 0;
				read_nb			= 0;
				read1			= 0;
				read2			= 0;
			}
		} 	// end for loop

		//free(offsets_to_recv);
		//free(sizes_to_recv);
	} // end else
	

	if ( rank_num == (proc_num - 1) ){

		// we complete the last chunk
		// with the last reads
		begin_offset_chunk[index_in_chunk] 	= begin_offset;
		begin_offset_chunk_2[index_in_chunk]	= begin_offset_2;
		chunk_size[index_in_chunk] 		= bytes_in_chunk;
		chunk_size_2[index_in_chunk]		= bytes_in_chunk_2;

		// the total num reads in the file is : grand_total_num_reads
		// grand total - sum from 0 to index in chunk -1 = num reads in that chunk
		int i_chunks	 	= 0;
		size_t sum_reads 	= 0;
		size_t sum_reads_2	= 0;

		for(i_chunks = 0; i_chunks < index_in_chunk; i_chunks++){
			sum_reads 	+= reads_in_chunk[i_chunks];
			sum_reads_2 += reads_in_chunk_2[i_chunks];
		}

		reads_in_chunk[index_in_chunk]		= read1;
		reads_in_chunk_2[index_in_chunk]	= read2;

		index_in_chunk 	+=1;
		*chunk_count 	+=1;

	}
	if (rank_num < (proc_num -1)){

		//we send to rank + 1
		MPI_Send(&counter_bases,
				1,
				MPI_LONG_LONG_INT,
				rank_num + 1,
				0,
				MPI_COMM_WORLD);

		MPI_Send(&begin_offset,
				1,
				MPI_LONG_LONG_INT,
				rank_num + 1,
				1,
				MPI_COMM_WORLD);

		MPI_Send(&begin_offset_2,
				1,
				MPI_LONG_LONG_INT,
				rank_num + 1,
				2,
				MPI_COMM_WORLD);

		MPI_Send(&bytes_in_chunk,
				1,
				MPI_LONG_LONG_INT,
				rank_num + 1,
				3,
				MPI_COMM_WORLD);

		MPI_Send(&bytes_in_chunk_2,
				1,
				MPI_LONG_LONG_INT,
				rank_num + 1,
				4,
				MPI_COMM_WORLD);

		//we send the type or reads we are going to send
		MPI_Send(&reads_type_to_send,
				1,
				MPI_INT,
				rank_num + 1,
				5,
				MPI_COMM_WORLD);

		// we send the size of the vector for the recieving part
		MPI_Send(&nb_reads_to_send,
				1,
				MPI_LONG_LONG_INT,
				rank_num + 1,
				6,
				MPI_COMM_WORLD);

		if ( reads_type_to_send == 0){
			offsets_to_send	= local_read_offsets + local_read_min - reads_recieved;
			sizes_to_send 	= local_read_size + local_read_min - reads_recieved;
			bytes_to_send	= local_read_bytes + local_read_min - reads_recieved;
		}
		else {
			offsets_to_send	= local_read_offsets_2 + local_read_min - reads_recieved_2;
			sizes_to_send 	= local_read_size_2 + local_read_min - reads_recieved_2;
			bytes_to_send	= local_read_bytes_2 + local_read_min - reads_recieved_2;
		}
		// we send the vector of offsets to the next rank

		MPI_Send(offsets_to_send,
				nb_reads_to_send,
				MPI_LONG_LONG_INT,
				rank_num + 1,
				7,
				MPI_COMM_WORLD);

		//we send the vector of sizes to the next rank
		MPI_Send(sizes_to_send,
				nb_reads_to_send,
				MPI_INT,
				rank_num + 1,
				8,
				MPI_COMM_WORLD);

		MPI_Send(bytes_to_send,
				nb_reads_to_send,
				MPI_LONG_LONG_INT,
				rank_num + 1,
				9,
				MPI_COMM_WORLD);

		// we send read1 and read2
		// number of reads already parsed
		MPI_Send(&read1,
				1,
				MPI_LONG_LONG_INT,
				rank_num + 1,
				10,
				MPI_COMM_WORLD);

		MPI_Send(&read2,
				1,
				MPI_LONG_LONG_INT,
				rank_num + 1,
				11,
				MPI_COMM_WORLD);

	}


	free(sizes_to_recv);
	free(bytes_to_recv);
	free(offsets_to_recv);

}


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
			size_t *chunk_count){

	//non MPI related resources	
	int i,bef,aft;
	int read_nb=0;   // have to use a buffer value instead of the iterator directly in case there are more reads handled by the process than there are in a chunk.
	int nb_reads_to_send		= 0;
	int nb_reads_to_recv		= 0;
	int reads_type_to_send		= 0; 	// use to know if we shall send forward (0) or backward (1)
	int reads_type_to_recv		= 0; 	// use to know if we shall recieve forward (0) or backward (1)
	int *sizes_to_send		= NULL;
	int *sizes_to_recv		= NULL;
	int *p_sizes  			= NULL; //reserved pointer
	int count 			= 0;
	int res 			= 0;
	size_t *p_bytes  		= NULL; //reserved pointer
	size_t u1			= 0;
	size_t local_read_min		= 0;
	size_t *bytes_to_send		= NULL;
	size_t *bytes_to_recv		= NULL;
	size_t *p_offset 		= NULL; //reserved pointer
	size_t reads_recieved		= 0;
	size_t *offsets_to_send		= NULL;
	size_t *offsets_to_recv		= NULL;
	size_t counter_bases  	 	= 0;
	size_t bytes_in_chunk 	 	= 0;
	size_t index_in_chunk 		= 0;
	size_t begin_offset 		= 0;
	size_t bases_previous_chunk 	= 0;
	size_t offset_previous_chunk 	= 0;
	size_t size_previous_chunk 	= 0;
	size_t x=0,y=0; 					// iterators of the number of reads in each file
	size_t read1=0,read2=0;
	size_t size_chunk = 100;
	size_t tmp_cnt=0;

	MPI_Status status;

	if (rank_num > 0){
		// we wait the previous rank to send the final size chunk
		// and offset of the chunk
		MPI_Recv(&bases_previous_chunk,
				1,
				MPI_LONG_LONG_INT,
				rank_num - 1,
				0,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);

		MPI_Recv(&offset_previous_chunk,
				1,
				MPI_LONG_LONG_INT,
				rank_num - 1,
				1,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);

		MPI_Recv(&size_previous_chunk,
				1,
				MPI_LONG_LONG_INT,
				rank_num - 1,
				2,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);

		//we send the type or reads we are going to send
		MPI_Recv(&reads_type_to_recv,
				1,
				MPI_INT,
				rank_num - 1,
				3,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);

		//we send the size of the vector for the recieving part
		MPI_Recv(&nb_reads_to_recv,
				1,
				MPI_LONG_LONG_INT,
				rank_num - 1,
				4,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);

		//we recieve for forward reads part
		reads_recieved 		= nb_reads_to_recv;
		
		//here we malloc the size
		sizes_to_recv 		= realloc(sizes_to_recv, nb_reads_to_recv * sizeof(int));
		bytes_to_recv 		= realloc(bytes_to_recv, nb_reads_to_recv * sizeof(size_t));
		offsets_to_recv 	= realloc(offsets_to_recv, nb_reads_to_recv * sizeof(size_t));

		assert(sizes_to_recv != NULL);
		assert(bytes_to_recv != NULL);
		assert(offsets_to_recv != NULL);

		//we send the vector of offsets to the next rank
		MPI_Recv(offsets_to_recv,
				nb_reads_to_recv,
				MPI_LONG_LONG_INT,
				rank_num - 1,
				5,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);

		//we send the vector of sizes to the next rank
		MPI_Recv(sizes_to_recv,
				nb_reads_to_recv,
				MPI_INT,
				rank_num - 1,
				6,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);

		MPI_Recv(bytes_to_recv,
				nb_reads_to_recv,
				MPI_LONG_LONG_INT,
				rank_num - 1,
				7,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);

		// we send read1 and read2
		// number of reads already parsed
		MPI_Recv(&read1,
				1,
				MPI_LONG_LONG_INT,
				rank_num - 1,
				8,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);

	}

	if (rank_num == 0){
		
		// here you find which file buffer of the rank 0 holds the least reads
		local_read_min = local_num_reads;

		// we are only in forward
		reads_type_to_send 	= 0;
		
		for (u1 = 0; u1 < local_read_min; u1++){

			if (u1 >= 1){
				assert(*(local_read_offsets + x)   != *(local_read_offsets + x - 1));
			}

			// we for offsets multiple of maxsize
			if (counter_bases == 0){
				begin_offset 	= *(local_read_offsets + x);
			}

			assert( (*(local_read_size + x)) != 0 );
			assert( (*(local_read_bytes + x)) != 0 );

			counter_bases 	+= *(local_read_size + x);
			bytes_in_chunk  += *(local_read_bytes + x);
			x++;read1++;
			
			if ( counter_bases > maxsiz){

				//the current chunk is full
				//then we have the number of bases wanted
				begin_offset_chunk[index_in_chunk] 		= begin_offset;
				chunk_size[index_in_chunk] 	   			= bytes_in_chunk;
				reads_in_chunk[index_in_chunk]	   		= read1;
				*chunk_count 			   += 1;
				index_in_chunk 			   += 1;
				//we reset counter of bases
				counter_bases  		= 0;
				bytes_in_chunk		= 0;
				read_nb				= 0;
				read1				= 0;
				read2				= 0;
			}
		} //end for loop on local_read_min
	}//end if rank_num == 0
	else{

		//receive the info about the current chunk being filled
		counter_bases 		= bases_previous_chunk;
		bytes_in_chunk 		= size_previous_chunk;
		begin_offset		= offset_previous_chunk;
		
		///receive the nb of reads sent
		//malloc a vector to its size
		//fill it with the sent data
		//redo the calculus of the nb of reads to treat and to send
		//once you get that new number, do the rest of the treatments
		//once you're finished check that everything went smoothly and send the right nb ofreads to the next rank.
		x=0; y=0;
		size_t counter_tmp = 0;
		///here you find which file buffer of the rank 0 that holds the last reads

		local_read_min = local_num_reads + reads_recieved;

		
		nb_reads_to_send 	= (local_num_reads + reads_recieved) - local_read_min;
		reads_type_to_send 	= 0; //we send forward
		
		// here we initialyze p_offset and p_sizes
		// according to recieve part
		// we recieve from forward reads part
		p_offset	= offsets_to_recv;
		p_sizes		= sizes_to_recv;
		p_bytes		= bytes_to_recv;
		
		for (u1 = 0; u1 < local_read_min; u1++){

			if (u1 >= 1){
				assert(*(p_offset + x) != *(p_offset + x - 1));
			}
			
			// when we reached the number of reads recieved
			// we switch of vector and pass to the local vector
			if (counter_tmp == nb_reads_to_recv){

					p_offset	= local_read_offsets;
					p_sizes		= local_read_size;
					p_bytes		= local_read_bytes;
					x 			= 0; //only pos in forward buff are reset
			} //end condition switch

			//we for offsets multiple of maxsize
			if (counter_bases == 0){
				begin_offset 	= *(p_offset + x);
			}

			counter_bases 	+= *(p_sizes + x);
			bytes_in_chunk  += *(p_bytes + x);
			x++; read1++;

			counter_tmp++;

			if ( counter_bases > maxsiz){

				//the current chunk is full
				//then we have the number of bases wanted
				begin_offset_chunk[index_in_chunk] 		= begin_offset;
				chunk_size[index_in_chunk] 	   			= bytes_in_chunk;
				
				reads_in_chunk[index_in_chunk]	   		= read1;
				*chunk_count		+= 1;
				index_in_chunk 		+= 1;
				//we reset counter of bases
				counter_bases  		= 0;
				bytes_in_chunk 		= 0;
				read_nb				= 0;
				read1				= 0;
				
			}
		} 	// end for loop

		//free(offsets_to_recv);
		//free(sizes_to_recv);
	} // end else
	

	if ( rank_num == (proc_num - 1) ){

		// we complete the last chunk
		// with the last reads
		begin_offset_chunk[index_in_chunk] 		= begin_offset;
		chunk_size[index_in_chunk] 				= bytes_in_chunk;
		
		// the total num reads in the file is : grand_total_num_reads
		// grand total - sum from 0 to index in chunk -1 = num reads in that chunk
		int i_chunks	 	= 0;
		size_t sum_reads 	= 0;
		
		for(i_chunks = 0; i_chunks < index_in_chunk; i_chunks++){
			sum_reads 	+= reads_in_chunk[i_chunks];
		}

		reads_in_chunk[index_in_chunk]		= read1;
		
		index_in_chunk 	+=1;
		*chunk_count 	+=1;

	}
	if (rank_num < (proc_num -1)){

		//we send to rank + 1
		MPI_Send(&counter_bases,
				1,
				MPI_LONG_LONG_INT,
				rank_num + 1,
				0,
				MPI_COMM_WORLD);

		MPI_Send(&begin_offset,
				1,
				MPI_LONG_LONG_INT,
				rank_num + 1,
				1,
				MPI_COMM_WORLD);

		MPI_Send(&bytes_in_chunk,
				1,
				MPI_LONG_LONG_INT,
				rank_num + 1,
				2,
				MPI_COMM_WORLD);

		//we send the type or reads we are going to send
		MPI_Send(&reads_type_to_send,
				1,
				MPI_INT,
				rank_num + 1,
				3,
				MPI_COMM_WORLD);

		// we send the size of the vector for the recieving part
		MPI_Send(&nb_reads_to_send,
				1,
				MPI_LONG_LONG_INT,
				rank_num + 1,
				4,
				MPI_COMM_WORLD);

		if ( reads_type_to_send == 0){
			offsets_to_send	= local_read_offsets + local_read_min - reads_recieved;
			sizes_to_send 	= local_read_size + local_read_min - reads_recieved;
			bytes_to_send	= local_read_bytes + local_read_min - reads_recieved;
		}
		// we send the vector of offsets to the next rank

		MPI_Send(offsets_to_send,
				nb_reads_to_send,
				MPI_LONG_LONG_INT,
				rank_num + 1,
				5,
				MPI_COMM_WORLD);

		//we send the vector of sizes to the next rank
		MPI_Send(sizes_to_send,
				nb_reads_to_send,
				MPI_INT,
				rank_num + 1,
				6,
				MPI_COMM_WORLD);

		MPI_Send(bytes_to_send,
				nb_reads_to_send,
				MPI_LONG_LONG_INT,
				rank_num + 1,
				7,
				MPI_COMM_WORLD);

		// we send read1 and read2
		// number of reads already parsed
		MPI_Send(&read1,
				1,
				MPI_LONG_LONG_INT,
				rank_num + 1,
				8,
				MPI_COMM_WORLD);
	}

	free(sizes_to_recv);
	free(bytes_to_recv);
	free(offsets_to_recv);
}

///Function used to map the reference genome whiwh will be used in the alignment.
//Parameters (ordered):	file containing the reference genoma
//			count of top level elements in the file
//			bwa parameter
//			i don't know what this is
//			pointer on he window of shared memory used to save the mapped indexes
void map_indexes(char *file_ref, ktp_aux_t *aux, MPI_Win *win_shr_ref, MPI_Win *win_shr_pac)
{
	   

        uint64_t tim = __rdtsc();
        fprintf(stderr, "* Ref file: %s\n", file_ref);
        aux->fmi = new FMI_search(file_ref);
        int res;
        MPI_Status status;
        uint8_t *shared_ref;    
        MPI_Comm comm_shr;
        int count4, rank_shr = 0;    
        size_t size_read, size_tot;
 
        res = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_shr);
        assert(res == MPI_SUCCESS);

        res = MPI_Comm_rank(comm_shr, &rank_shr);
        assert(res == MPI_SUCCESS);

        char binary_seq_file[200];
        sprintf(binary_seq_file, "%s.0123", file_ref);

        char pac_file[200];
        sprintf(pac_file, "%s.pac", file_ref);

        MPI_File fh_ref_file, fh_pac_file;
        MPI_Aint size_shr_ref, size_shr_pac;
        MPI_Offset size_map_ref,size_map_pac;
        
        uint8_t *addr_map_ref;
        uint8_t *addr_map_pac;       


        MPI_Info win_info;
        MPI_Info_create(&win_info);
        MPI_Info_set(win_info, "alloc_shared_non_contig", "true");

        res = MPI_File_open(MPI_COMM_WORLD, binary_seq_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_ref_file);
        assert(res == MPI_SUCCESS);

        res = MPI_File_open(MPI_COMM_WORLD, pac_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_pac_file);
        assert(res == MPI_SUCCESS);

        res = MPI_File_get_size(fh_pac_file, &size_map_pac);
        assert(res == MPI_SUCCESS);
        
        res = MPI_File_get_size(fh_ref_file, &size_map_ref);
        assert(res == MPI_SUCCESS);

        size_shr_ref = (rank_shr == 0) ? size_map_ref : 0;
        res = MPI_Win_allocate_shared(size_shr_ref * sizeof(uint8_t), 8 , win_info, comm_shr, &addr_map_ref, win_shr_ref);

        size_shr_pac = (rank_shr == 0) ? size_map_pac : 0;
        res = MPI_Win_allocate_shared(size_shr_pac * sizeof(uint8_t), 8 , win_info, comm_shr, &addr_map_pac, win_shr_pac);

        res = MPI_Win_shared_query(*win_shr_pac, MPI_PROC_NULL, &size_shr_pac, &res, &addr_map_pac);
        assert(res == MPI_SUCCESS);

        res = MPI_Win_shared_query(*win_shr_ref, MPI_PROC_NULL, &size_shr_ref, &res, &addr_map_ref);
        assert(res == MPI_SUCCESS);

        uint8_t *c = addr_map_ref;
        uint8_t *d = addr_map_pac;
        
        while(rank_shr == 0) {
            res = MPI_File_read(fh_pac_file, d, size_shr_pac, MPI_UINT8_T, &status);
            assert(res == MPI_SUCCESS);
            res = MPI_Get_count(&status, MPI_UINT8_T, &count4);
            assert(res == MPI_SUCCESS);
            if (count4 == 0) break;
            d += count4;
            size_read += count4;
            size_tot += count4;
        }
        
        size_read = 0;

        while(rank_shr == 0) {
            res = MPI_File_read(fh_ref_file, c, size_shr_ref, MPI_UINT8_T, &status);
            assert(res == MPI_SUCCESS);
            res = MPI_Get_count(&status, MPI_UINT8_T, &count4);
            assert(res == MPI_SUCCESS);
            if (count4 == 0) break;
            c += count4;
            size_read += count4;
            size_tot += count4;
        }

        res = MPI_Win_fence(0, *win_shr_ref);
        assert(res == MPI_SUCCESS);

        res = MPI_Win_fence(0, *win_shr_pac);
        assert(res == MPI_SUCCESS);

        res = MPI_File_close(&fh_pac_file);
        assert(res == MPI_SUCCESS);

        res = MPI_File_close(&fh_ref_file);
        assert(res == MPI_SUCCESS);

        aux->ref_string = addr_map_ref;
        aux->fmi->load_index_mpi(addr_map_pac);
	
}







// Function used to create the header of the .sam result file
// Parameters (ordered):	result file
//			bwa parameter
//			count of top level elements in the file
//			header string
//			rg = read group string
void create_sam_header(char *file_out, bwaidx_fm_t *indix, int *count, char *hdr_line, char *rg_line, char* pg_line, int rank_num)
{

	//MPI resources
	MPI_File fh_out;
	MPI_Status status;

	//non MPI related resources
	int bef,aft;
	int res;
	
	//the rank 0 takes care of that task
	
		int s, len;
		char *buff;
		struct stat stat_file_out;

		//We test if the output sam exists
		if ( stat(file_out, &stat_file_out)  != -1 ) {
		    res = MPI_File_delete(file_out, MPI_INFO_NULL);
		    assert(res == MPI_SUCCESS);
		}
				
		res = MPI_File_open(MPI_COMM_SELF, file_out,  MPI_MODE_CREATE|MPI_MODE_WRONLY , MPI_INFO_NULL, &fh_out);
		assert(res == MPI_SUCCESS);
		/* Add reference sequence lines */
		for (s = 0; s < (*indix).bns->n_seqs; ++s) {
		len = asprintf(&buff, "@SQ\tSN:%s\tLN:%d\n", (*indix).bns->anns[s].name, (*indix).bns->anns[s].len);
		res = MPI_File_write(fh_out, buff, len, MPI_CHAR, &status);
		assert(res == MPI_SUCCESS);
		res = MPI_Get_count(&status, MPI_CHAR, count);
		assert(res == MPI_SUCCESS);
		assert(*count == len);
		free(buff);
		}
		/* Add header lines */
		if (hdr_line != NULL) {
		len = asprintf(&buff, "%s\n", hdr_line);
		res = MPI_File_write(fh_out, buff, len, MPI_CHAR, &status);
		assert(res == MPI_SUCCESS);
		res = MPI_Get_count(&status, MPI_CHAR, count);
		assert(res == MPI_SUCCESS);
		assert(*count == len);
		free(buff);
		}
		/* Add read group line */
		if (rg_line != NULL) {
		len = asprintf(&buff, "%s\n", rg_line);
		res = MPI_File_write(fh_out, buff, len, MPI_CHAR, &status);
		assert(res == MPI_SUCCESS);
		res = MPI_Get_count(&status, MPI_CHAR, count);
		assert(res == MPI_SUCCESS);
		assert(*count == len);
		free(buff);
		}

        /* Add programm line */
        if (pg_line != NULL) {
        len = asprintf(&buff, "%s\n", pg_line);
        res = MPI_File_write(fh_out, buff, len, MPI_CHAR, &status);
        assert(res == MPI_SUCCESS);
        res = MPI_Get_count(&status, MPI_CHAR, count);
        assert(res == MPI_SUCCESS);
        assert(*count == len);
        free(buff);
        }

		res = MPI_File_close(&fh_out);
		assert(res == MPI_SUCCESS);
	
}

void create_bam_header(char *file_out, bwaidx_fm_t *indix, int *count, 
	char *hdr_line, char *rg_line, char *pg_line, int rank_num, int compression_level){

			int s = 0;
			int s1 = 0;
			int len = 0;
			struct stat stat_path;
			char *buff[(indix)->bns->n_seqs + 2];
			MPI_Status status;
			int res = 0;
			for (s = 0; s < (indix)->bns->n_seqs; s++) 
				len += asprintf(&buff[s], "@SQ\tSN:%s\tLN:%d\n", (indix)->bns->anns[s].name, (indix)->bns->anns[s].len);
			
			s1 = s;
			
			if (hdr_line != NULL) {
				len += asprintf(&buff[s++], "%s\n", hdr_line);
				s1++;	
			}
			if (rg_line != NULL) {
				len += asprintf(&buff[s++], "%s\n", rg_line);			
				s1++;
			}
             if (pg_line != NULL) {
                len += asprintf(&buff[s++], "%s\n", pg_line);
                s1++;
            }

			char *buff_header = malloc(len + 1);
            		buff_header[len] = 0;
            		int disp = 0;
			char *z = buff_header;
			char *z1;	
			
			for (s = 0; s < s1; s++){
				z1 = buff[s];
                		while (*z1) {*z++ = *z1++;}
				free(buff[s]);
			}	
            
            		BGZF *fp_header;
			fp_header = calloc(1, sizeof(BGZF));
			uint8_t *compressed_header = NULL;
			int compressed_size_header = 0;

			int block_length = MAX_BLOCK_SIZE;
			int bytes_written;
			int length = strlen(buff_header);

			fp_header->open_mode = 'w';
			fp_header->uncompressed_block_size = MAX_BLOCK_SIZE;
			fp_header->uncompressed_block = malloc(MAX_BLOCK_SIZE);
			fp_header->compressed_block_size = MAX_BLOCK_SIZE;
			fp_header->compressed_block = malloc(MAX_BLOCK_SIZE);
			fp_header->cache_size = 0;
			fp_header->block_address = 0;
			fp_header->block_offset = 0;
			fp_header->block_length = 0;
			fp_header->compress_level = compression_level < 0? Z_DEFAULT_COMPRESSION : compression_level; 
			if (fp_header->compress_level > 9) fp_header->compress_level = Z_DEFAULT_COMPRESSION;

			const bgzf_byte_t *input = (void *)buff_header;

			if (fp_header->uncompressed_block == NULL) {
				fp_header->uncompressed_block = malloc(fp_header->uncompressed_block_size);
			}

			block_length = fp_header->uncompressed_block_size;
			bytes_written = 0;
			compressed_header =  malloc(strlen(buff_header) * sizeof(uint8_t));

			while (bytes_written < length) {
				int copy_length = bgzf_min(block_length - fp_header->block_offset, length - bytes_written);
				bgzf_byte_t* buffer = fp_header->uncompressed_block;
				memcpy(buffer + fp_header->block_offset, input, copy_length);
				fp_header->block_offset += copy_length;
				input += copy_length;
				bytes_written += copy_length;
				
				while (fp_header->block_offset > 0) {
					int block_length;
					block_length = deflate_block(fp_header, fp_header->block_offset);
					memcpy(compressed_header + compressed_size_header, fp_header->compressed_block, block_length);
					compressed_size_header += block_length;
					fp_header->block_address += block_length;
				}
			}
			free(buff_header);

			if ( stat(file_out, &stat_path)  != -1 ) {
		        res = MPI_File_delete(file_out, MPI_INFO_NULL);
               	assert(res == MPI_SUCCESS);
            }

			MPI_File fh_out_h;
		    	res = MPI_File_open(MPI_COMM_SELF, file_out, MPI_MODE_CREATE|MPI_MODE_WRONLY|MPI_MODE_APPEND, MPI_INFO_NULL, &fh_out_h);
            		assert(res == MPI_SUCCESS);
 			size_t compSize = compressed_size_header;
			res = MPI_File_write(fh_out_h, compressed_header, compSize, MPI_BYTE, &status);
			assert(res == MPI_SUCCESS);
			res = MPI_Get_count(&status, MPI_BYTE, count);
			assert(res == MPI_SUCCESS);
			assert(*count == (int)compSize);
			res = MPI_File_close(&fh_out_h);
		    	assert(res == MPI_SUCCESS);

}

void create_bam_header_by_chr_file(char *file_out[], bwaidx_t *indix, int *count, char *hdr_line, char *rg_line, int rank_num, int compression_level, int dofixmate){


		
        MPI_Status status;
	int s = 0;
	int i = 0;
	int len = 0;
	int res = 0;
        char *buff;
        struct stat stat_file_out;
	int incrmnt = 0;
	if (!dofixmate) incrmnt += 1;

	if (rank_num == 0) {

		// Remove sam files if already exists      
		for (s = 0; s < (indix)->bns->n_seqs + incrmnt; ++s) {
			if ( stat(file_out[s], &stat_file_out)  != -1 ) {
			        res = MPI_File_delete(file_out[s], MPI_INFO_NULL);
                                assert(res == MPI_SUCCESS);
		        }
		}
		
		

                int s1 = 0;
                struct stat stat_path;
                char *buff[(indix)->bns->n_seqs + 2];
                MPI_Status status;
                        
                for (s = 0; s < (indix)->bns->n_seqs; s++)
                	len += asprintf(&buff[s], "@SQ\tSN:%s\tLN:%d\n", (indix)->bns->anns[s].name, (indix)->bns->anns[s].len);

                s1 = s;

                if (hdr_line != NULL) {
                	len += asprintf(&buff[s++], "%s\n", hdr_line);
                        s1++;
                }
                if (rg_line != NULL) {
                        len += asprintf(&buff[s++], "%s\n", rg_line);
                        s1++;
                }
                char *buff_header = malloc(len + 1);
                buff_header[len] = 0;
                int disp = 0;
                char *z = buff_header;
                char *z1;

                for (s = 0; s < s1; s++){
                       z1 = buff[s];
                       while (*z1) {*z++ = *z1++;}
                       free(buff[s]);
                }

		

                        BGZF *fp_header;
                        fp_header = calloc(1, sizeof(BGZF));
                        uint8_t *compressed_header = NULL;
                        int compressed_size_header = 0;

                        int block_length = MAX_BLOCK_SIZE;
                        int bytes_written;
                        int length = strlen(buff_header);

                        fp_header->open_mode = 'w';
                        fp_header->uncompressed_block_size = MAX_BLOCK_SIZE;
                        fp_header->uncompressed_block = malloc(MAX_BLOCK_SIZE);
                        fp_header->compressed_block_size = MAX_BLOCK_SIZE;
                        fp_header->compressed_block = malloc(MAX_BLOCK_SIZE);
                        fp_header->cache_size = 0;
                        fp_header->block_address = 0;
                        fp_header->block_offset = 0;
                        fp_header->block_length = 0;
                        fp_header->compress_level = compression_level < 0? Z_DEFAULT_COMPRESSION : compression_level;
			            if (fp_header->compress_level > 9) fp_header->compress_level = Z_DEFAULT_COMPRESSION;

                        const bgzf_byte_t *input = (void *)buff_header;

                        if (fp_header->uncompressed_block == NULL) {
                                fp_header->uncompressed_block = malloc(fp_header->uncompressed_block_size);
                        }

                        block_length = fp_header->uncompressed_block_size;
                        bytes_written = 0;
                        compressed_header =  malloc(strlen(buff_header) * sizeof(uint8_t));

                        while (bytes_written < length) {
                                int copy_length = bgzf_min(block_length - fp_header->block_offset, length - bytes_written);
                                bgzf_byte_t* buffer = fp_header->uncompressed_block;
                                memcpy(buffer + fp_header->block_offset, input, copy_length);
                                fp_header->block_offset += copy_length;
                                input += copy_length;
                                bytes_written += copy_length;

                                while (fp_header->block_offset > 0) {
                                        int block_length;
                                        block_length = deflate_block(fp_header, fp_header->block_offset);
                                        memcpy(compressed_header + compressed_size_header, fp_header->compressed_block, block_length);
                                        compressed_size_header += block_length;
                                        fp_header->block_address += block_length;
                                }
                        }
                        free(buff_header);

		                for (i = 0; i < (indix)->bns->n_seqs + incrmnt; ++i) {                        
                            MPI_File fh_out_h;
                            res = MPI_File_open(MPI_COMM_SELF, file_out[i], MPI_MODE_CREATE|MPI_MODE_WRONLY|MPI_MODE_APPEND, MPI_INFO_NULL, &fh_out_h);
                            assert(res == MPI_SUCCESS);
                            size_t compSize = compressed_size_header;
                            res = MPI_File_write(fh_out_h, compressed_header, compSize, MPI_BYTE, &status);
                            assert(res == MPI_SUCCESS);
                            res = MPI_Get_count(&status, MPI_BYTE, count);
                            assert(res == MPI_SUCCESS);
                            assert(*count == (int)compSize);
                            res = MPI_File_close(&fh_out_h);
                            assert(res == MPI_SUCCESS);
		                }
	}
}



/****************************************
 *  
 *       Thread Structures and functions
 *  
 ****************************************/

void *call_fixmate(void *threadarg){


        struct thread_data *my_data;
        my_data = (struct thread_data *) threadarg;
        size_t total_sam_line = 0;
        int next, i, n, m;
        char currentLine[MAX_CHAR_SIZE];
        bwaidx_t *indix_tmp = my_data->indix_thr;
        int rank_num = my_data->job_rank;
        bseq1_t *seqs =  my_data->seqs_thr;
        int current_sam_line1 = 0;
        int start = my_data->start_index_seqs;

        int current_sam_line2=0;

        n = my_data->thread_id * 2;
        m = my_data->thread_id * 2 + 1;
        int incr = my_data->total_thread * 2;
        int total_reads_parsed = 0;
	int total_reads_per_thread =  my_data->total_reads / my_data->total_thread;
        int read_num_1 = 0;
        int read_num_2 = 0;
        int tmp = 0;
        if ( my_data->thread_id == ( my_data->total_thread - 1 ))
                total_reads_per_thread = my_data->total_reads - ( total_reads_per_thread * (my_data->total_thread - 1));

        do {
                read_num_1=0;
                read_num_2=0;
		fixmate (rank_num, &(seqs[n]), &(seqs[m]), &read_num_1, &read_num_2, indix_tmp);

                tmp += read_num_1 + read_num_2;
                n = n + incr;
                m = m + incr;
                current_sam_line1=0;
                current_sam_line2=0;
                total_reads_parsed = total_reads_parsed + 2;
        } while ( ( m < my_data->total_reads) );
        my_data->total_lines = tmp;
        pthread_exit((void *)&my_data);

}

void *compress_and_write_thread(void *threadarg){

        struct thread_data_compress *my_data;
        my_data = (struct thread_data_compress *) threadarg;
        int compression_level = my_data->comp_level;
        int rank_num = my_data->job_rank;
        int res = 0;
        int count = 0;
        MPI_Status status;
        bseq1_t *seqs =  my_data->seqs_thr;
        MPI_File mpifh = my_data->fh_out;
        int start_index = (my_data->total_reads / my_data->total_thread) * my_data->thread_id;
        int end_index = 0;
        if (my_data->thread_id < my_data->total_thread)
                end_index = (my_data->total_reads / my_data->total_thread) * (my_data->thread_id + 1) - 1;
        else
                end_index = my_data->total_reads;
        int n = 0;
        size_t localsize = 0;
        for (n = start_index; n < end_index; n++) {
                seqs[n].l_seq = strlen(seqs[n].sam);
                localsize += seqs[n].l_seq;
        }
        assert(localsize <= INT_MAX);
        char *buffer_out = malloc((localsize + 1)*sizeof(char));
	    assert(buffer_out != NULL);
	    buffer_out[localsize] = 0;
        char *p = buffer_out;
        for (n = start_index; n < end_index; n++) {
                memmove(p, seqs[n].sam, seqs[n].l_seq);
                p += seqs[n].l_seq;
        }
        
        uint8_t *p2 = buffer_out;
        size_t compressed_size =0;
        BGZF *fp;
        fp = calloc(1, sizeof(BGZF));
        int block_length = MAX_BLOCK_SIZE;
        int bytes_written;
        int length = strlen(buffer_out);
        fp->open_mode = 'w';
        fp->uncompressed_block_size = MAX_BLOCK_SIZE;
        fp->uncompressed_block = malloc(MAX_BLOCK_SIZE);
        fp->compressed_block_size = MAX_BLOCK_SIZE;
        fp->compressed_block = malloc(MAX_BLOCK_SIZE);
        fp->cache_size = 0;
        fp->cache = kh_init(cache);
        fp->block_address = 0;
        fp->block_offset = 0;
        fp->block_length = 0;
        fp->compress_level = compression_level < 0 ? Z_DEFAULT_COMPRESSION : compression_level;
        if (fp->compress_level > 9) {
                fp->compress_level = Z_DEFAULT_COMPRESSION;
        }

        const bgzf_byte_t *input = (void *)buffer_out;

        if (fp->uncompressed_block == NULL) {
                fp->uncompressed_block = malloc(fp->uncompressed_block_size);
        }
        input = (void *)buffer_out;
        block_length = fp->uncompressed_block_size;
        bytes_written = 0;

        while (bytes_written < length) {
                int copy_length = bgzf_min(block_length - fp->block_offset, length - bytes_written);
                bgzf_byte_t *buffer = fp->uncompressed_block;
                memcpy(buffer + fp->block_offset, input, copy_length);
                fp->block_offset += copy_length;
                input += copy_length;
                bytes_written += copy_length;

                while (fp->block_offset > 0) {
                        int block_length;
                        block_length = deflate_block(fp, fp->block_offset);

                        memcpy(p2 + compressed_size, fp->compressed_block, block_length);
                        compressed_size += block_length;
                        fp->block_address += block_length;
                }
        }

        size_t compSize = compressed_size;
        res = MPI_File_write_shared(mpifh, buffer_out, compSize, MPI_BYTE, &status);
        assert(res == MPI_SUCCESS);
        res = MPI_Get_count(&status, MPI_BYTE, &count);
        assert(res == MPI_SUCCESS);
        assert(count == (int)compSize);

        free(fp->uncompressed_block);
        free(fp->compressed_block);
        my_data->thr_comp_sz = compressed_size;
        kh_destroy(cache, fp->cache);
        free(buffer_out);
}

void *compress_thread_by_chr(void *threadarg){

        struct thread_data_compress_by_chr *my_data;
        my_data = (struct thread_data_compress_by_chr *) threadarg;
        int compression_level = my_data->comp_level;
        int rank_num = my_data->job_rank;
	    bwaidx_t *indix = my_data->indix_thr;
	    int total_thread = my_data->total_thread;
	    size_t total_sam_line = my_data->total_sam_line;
	    int comp_level = my_data->comp_level;
        char **start_addr = my_data->start_addr;
        int *line_size_to_cpy = my_data->line_size_to_cpy;
        int *sam_buff_dest = my_data->sam_buff_dest;
        int incrmnt = my_data->incrmnt;
	    int *add_in_disc = my_data->add_in_disc;
	    MPI_File *fh_out = my_data->fh_out;
	    MPI_Status status;
	    int count;
	    int res;
	    double bef, aft;
	    int nbchr = indix->bns->n_seqs + incrmnt;
	    int start_index = (total_sam_line / total_thread) * my_data->thread_id;
        int end_index = 0;
        if (my_data->thread_id < my_data->total_thread)
                end_index = (my_data->total_sam_line / my_data->total_thread) * (my_data->thread_id + 1) - 1;
        else
                end_index = my_data->total_sam_line;
        int n = 0;
	    size_t *chr_buff_size = calloc( indix->bns->n_seqs + incrmnt, sizeof(size_t));
	    for (n = start_index; n < end_index; n++){
		    chr_buff_size[sam_buff_dest[n]] += line_size_to_cpy[n];  
		    if ( add_in_disc[n] ) chr_buff_size[nbchr -2] += line_size_to_cpy[n];
	    }
	    char *buffer_out_vec[indix->bns->n_seqs + incrmnt];
	    for (n = 0; n < (indix->bns->n_seqs + incrmnt); n++) {
            assert(chr_buff_size[n] <= INT_MAX);
            if (chr_buff_size[n]){
                buffer_out_vec[n] = calloc( chr_buff_size[n] + 1, sizeof(char));
                assert(buffer_out_vec[n] != NULL);
                buffer_out_vec[n][chr_buff_size[n]] = '\0';
            }
        }
        size_t *actual_size = calloc(nbchr, sizeof(size_t));
        char *p_temp2;
        int u = 0;
        for (n = start_index; n < end_index; n++) {

  		p_temp2 = buffer_out_vec[sam_buff_dest[n]] + actual_size[sam_buff_dest[n]];
        memmove(p_temp2, start_addr[n], line_size_to_cpy[n]);
        actual_size[sam_buff_dest[n]] += line_size_to_cpy[n];
        if ( add_in_disc[n] ){
		    p_temp2 = buffer_out_vec[nbchr - 2 ] + actual_size[nbchr - 2 ];
            memmove(p_temp2, start_addr[n], line_size_to_cpy[n]);
            actual_size[nbchr -2] += line_size_to_cpy[n];
            }
        }
    	free(actual_size);
	
	    for (n = 0; n < (indix->bns->n_seqs + incrmnt); n++) {

		    if (chr_buff_size[n]){

        		char *buffer_out = buffer_out_vec[n];
        		uint8_t *compressed_buff = malloc(chr_buff_size[n] * sizeof(uint8_t));
        		uint8_t *p2 = compressed_buff;
        		size_t compressed_size =0;
        		BGZF *fp;
        		fp = calloc(1, sizeof(BGZF));
        		int block_length = MAX_BLOCK_SIZE;
        		int bytes_written;
        		int length = strlen(buffer_out);
        		fp->open_mode = 'w';
        		fp->uncompressed_block_size = MAX_BLOCK_SIZE;
        		fp->uncompressed_block = malloc(MAX_BLOCK_SIZE);
        		fp->compressed_block_size = MAX_BLOCK_SIZE;
        		fp->compressed_block = malloc(MAX_BLOCK_SIZE);
        		fp->cache_size = 0;
        		fp->cache = kh_init(cache);
        		fp->block_address = 0;
        		fp->block_offset = 0;
        		fp->block_length = 0;
        		fp->compress_level = compression_level < 0 ? Z_DEFAULT_COMPRESSION : compression_level;
        		if (fp->compress_level > 9) {
                		fp->compress_level = Z_DEFAULT_COMPRESSION;
        		}

	        	const bgzf_byte_t *input = (void *)buffer_out;

        		if (fp->uncompressed_block == NULL) {
                		fp->uncompressed_block = malloc(fp->uncompressed_block_size);
        		}
        		input = (void *)buffer_out;
        		block_length = fp->uncompressed_block_size;
        		bytes_written = 0;

        		while (bytes_written < length) {
                		int copy_length = bgzf_min(block_length - fp->block_offset, length - bytes_written);
                		bgzf_byte_t *buffer = fp->uncompressed_block;
                		memcpy(buffer + fp->block_offset, input, copy_length);
                		fp->block_offset += copy_length;
                		input += copy_length;
                		bytes_written += copy_length;

                		while (fp->block_offset > 0) {
                        		int block_length;
                        		block_length = deflate_block(fp, fp->block_offset);

                        		memcpy(p2 + compressed_size, fp->compressed_block, block_length);
                        		compressed_size += block_length;
                        		fp->block_address += block_length;
                		}
        		}
		
			size_t compSize = compressed_size;
                	bef = MPI_Wtime();
                	res = MPI_File_write_shared(fh_out[n], compressed_buff, compSize, MPI_BYTE, &status);
                	assert(res == MPI_SUCCESS);
                	res = MPI_Get_count(&status, MPI_BYTE, &count);
                	assert(res == MPI_SUCCESS);
                	assert(count == (int)compSize);

 	       		free(fp->uncompressed_block);
        		free(fp->compressed_block);
			free(buffer_out_vec[n]);        	
        		kh_destroy(cache, fp->cache);
        		free(compressed_buff);
		}
	}
	free(chr_buff_size);
	//free(buffer_out_vec);
}

void *compress_thread_by_chr_single(void *threadarg){

        struct thread_data_compress_by_chr_single *my_data;
        my_data = (struct thread_data_compress_by_chr_single *) threadarg;
        int compression_level = my_data->comp_level;
        int rank_num = my_data->job_rank;
        bwaidx_t *indix = my_data->indix_thr;
        int total_thread = my_data->total_thread;
        size_t total_sam_line = my_data->total_sam_line;
        int comp_level = my_data->comp_level;
        char **start_addr = my_data->start_addr;
        int *line_size_to_cpy = my_data->line_size_to_cpy;
        int *sam_buff_dest = my_data->sam_buff_dest;
        int incrmnt = my_data->incrmnt;
        MPI_File *fh_out = my_data->fh_out;
        MPI_Status status;
        int count;
        int res;
        double bef, aft;
        int nbchr = indix->bns->n_seqs + incrmnt;
        int start_index = (total_sam_line / total_thread) * my_data->thread_id;
        int end_index = 0;
        if (my_data->thread_id < my_data->total_thread)
                end_index = (my_data->total_sam_line / my_data->total_thread) * (my_data->thread_id + 1) - 1;
        else
                end_index = my_data->total_sam_line;
        int n = 0;
        size_t *chr_buff_size = calloc( indix->bns->n_seqs + incrmnt, sizeof(size_t));
        for (n = start_index; n < end_index; n++){
                chr_buff_size[sam_buff_dest[n]] += line_size_to_cpy[n];
        }
        char *buffer_out_vec[indix->bns->n_seqs + incrmnt];
        for (n = 0; n < (indix->bns->n_seqs + incrmnt); n++) {
                 assert(chr_buff_size[n] <= INT_MAX);
                 if (chr_buff_size[n]){
                           buffer_out_vec[n] = calloc( chr_buff_size[n] + 1, sizeof(char));
                           assert(buffer_out_vec[n] != NULL);
                           buffer_out_vec[n][chr_buff_size[n]] = '\0';
                 }
        }
        size_t *actual_size = calloc(nbchr, sizeof(size_t));
        char *p_temp2;
        int u = 0;
        for (n = start_index; n < end_index; n++) {

                p_temp2 = buffer_out_vec[sam_buff_dest[n]] + actual_size[sam_buff_dest[n]];
                memmove(p_temp2, start_addr[n], line_size_to_cpy[n]);
                actual_size[sam_buff_dest[n]] += line_size_to_cpy[n];
        }
        free(actual_size);

        for (n = 0; n < (indix->bns->n_seqs + incrmnt); n++) {

                if (chr_buff_size[n]){

                        char *buffer_out = buffer_out_vec[n];
                        uint8_t *compressed_buff = malloc(chr_buff_size[n] * sizeof(uint8_t));
                        uint8_t *p2 = compressed_buff;
                        size_t compressed_size =0;
                        BGZF *fp;
                        fp = calloc(1, sizeof(BGZF));
                        int block_length = MAX_BLOCK_SIZE;
                        int bytes_written;
                        int length = strlen(buffer_out);
                        fp->open_mode = 'w';
                        fp->uncompressed_block_size = MAX_BLOCK_SIZE;
                        fp->uncompressed_block = malloc(MAX_BLOCK_SIZE);
                        fp->compressed_block_size = MAX_BLOCK_SIZE;
                        fp->compressed_block = malloc(MAX_BLOCK_SIZE);
                        fp->cache_size = 0;
                        fp->cache = kh_init(cache);
                        fp->block_address = 0;
                        fp->block_offset = 0;
                        fp->block_length = 0;
                        fp->compress_level = compression_level < 0 ? Z_DEFAULT_COMPRESSION : compression_level;
                        if (fp->compress_level > 9) {
                                fp->compress_level = Z_DEFAULT_COMPRESSION;
                        }

                        const bgzf_byte_t *input = (void *)buffer_out;

                        if (fp->uncompressed_block == NULL) {
                                fp->uncompressed_block = malloc(fp->uncompressed_block_size);
                        }
                        input = (void *)buffer_out;
                        block_length = fp->uncompressed_block_size;
                        bytes_written = 0;

                        while (bytes_written < length) {
                                int copy_length = bgzf_min(block_length - fp->block_offset, length - bytes_written);
                                bgzf_byte_t *buffer = fp->uncompressed_block;
                                memcpy(buffer + fp->block_offset, input, copy_length);
                                fp->block_offset += copy_length;
                                input += copy_length;
                                bytes_written += copy_length;

                                while (fp->block_offset > 0) {
					int block_length;
                                        block_length = deflate_block(fp, fp->block_offset);

                                        memcpy(p2 + compressed_size, fp->compressed_block, block_length);
                                        compressed_size += block_length;
                                        fp->block_address += block_length;
                                }
                        }

                        size_t compSize = compressed_size;
                        bef = MPI_Wtime();
                        res = MPI_File_write_shared(fh_out[n], compressed_buff, compSize, MPI_BYTE, &status);
                        assert(res == MPI_SUCCESS);
                        res = MPI_Get_count(&status, MPI_BYTE, &count);
                        assert(res == MPI_SUCCESS);
                        assert(count == (int)compSize);

                        free(fp->uncompressed_block);
                        free(fp->compressed_block);
                        free(buffer_out_vec[n]);
                        kh_destroy(cache, fp->cache);
                        free(compressed_buff);
                }
        }
        free(chr_buff_size);
        //free(buffer_out_vec);
}
 

int getChr(char *str, char *chrNames[], int nbchr, char *tmp_chr) {
    int i = 0, found = 0, size;
    char *str1 = str, *str2;

    str2 = str1 + 1;

    for (; *str2 != '\t'; str2++);

    size = strlen(str1) - strlen(str2);
    assert(size != 0);

    for (i = 0; i < size; i++) {
        tmp_chr[i] = str1[i + 1];
    }

    tmp_chr[size - 1] = '\0';
    assert(strlen(tmp_chr) != 0);
  
    for (i = 0, found = 0; i < nbchr && !found; i++) {
        found = !strcmp(tmp_chr, chrNames[i]);
    }

    return i - 1;
}

void create_sam_header_by_chr_file(char *file_out[], 
		bwaidx_t *indix, int *count, char *hdr_line, char *rg_line, int rank_num)
{

        //MPI resources
        MPI_File fh_out[(*indix).bns->n_seqs];
        MPI_Status status;
        //non MPI related resources
        int bef, i, aft;
        int res;
        //the rank 0 takes care of that task
        if (rank_num == 0) {
             int s, len;
             char *buff;
	     struct stat stat_file_out;

	     // Remove sam files if already exists 	
             for (s = 0; s < (*indix).bns->n_seqs; ++s) {		
             	
		if ( stat(file_out[s], &stat_file_out)  != -1 ) {
        		res = MPI_File_delete(file_out[s], MPI_INFO_NULL);
                        assert(res == MPI_SUCCESS);
                }
		res = MPI_File_open(MPI_COMM_SELF, file_out[s], MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_out[s]);
             	assert(res == MPI_SUCCESS);
	     }

	      for (i = 0; i < (*indix).bns->n_seqs; ++i) {	
             	// Add reference sequence lines 
             	for (s = 0; s < (*indix).bns->n_seqs; ++s) {
             		len = asprintf(&buff, "@SQ\tSN:%s\tLN:%d\n", (*indix).bns->anns[s].name, (*indix).bns->anns[s].len);
                	res = MPI_File_write(fh_out[i], buff, len, MPI_CHAR, &status);
                	assert(res == MPI_SUCCESS);
                	res = MPI_Get_count(&status, MPI_CHAR, count);
                	assert(res == MPI_SUCCESS);
                	assert(*count == len);
                	free(buff);
              	}
	      	// Add header lines 
 	      	if (hdr_line != NULL) {
              		len = asprintf(&buff, "%s\n", hdr_line);
                	res = MPI_File_write(fh_out[i], buff, len, MPI_CHAR, &status);
                	assert(res == MPI_SUCCESS);
                	res = MPI_Get_count(&status, MPI_CHAR, count);
                	assert(res == MPI_SUCCESS);
                	assert(*count == len);
                	free(buff);
               	}
       	       	if (rg_line != NULL) {
               		len = asprintf(&buff, "%s\n", rg_line);
                        res = MPI_File_write(fh_out[i], buff, len, MPI_CHAR, &status);
                        assert(res == MPI_SUCCESS);
                        res = MPI_Get_count(&status, MPI_CHAR, count);
                        assert(res == MPI_SUCCESS);
                        assert(*count == len);
                        free(buff);
			res = MPI_File_close(&fh_out[i]);
	                assert(res == MPI_SUCCESS);
                	}
		}
        }

        bef = MPI_Wtime();
        res = MPI_Barrier(MPI_COMM_WORLD);
        assert(res == MPI_SUCCESS);
        aft = MPI_Wtime();
        xfprintf(stderr, "%s: synched processes (%.02f)\n", __func__, aft - bef);
}
/*
void compute_buffer_size_thr(void *thread_arg){

    struct struct_data_thread *my_data;
    my_data = (struct struct_data_thread *) thread_arg;
    int i = 0;
    for (i = my_data->begin_index; i < my_data->end_index; i++){
        my_data->seqs_thr[i].l_seq = strlen(my_data->seqs_thr[i].sam);
        my_data->size_thr += my_data->seqs_thr[i].l_seq;
    }
}

void copy_buffer_write_thr(void *thread_arg){
    struct struct_data_thread *my_data;
    my_data = (struct struct_data_thread *) thread_arg;

    my_data->buffer_out = malloc( my_data->size_thr + 1);
    assert(my_data->buffer_out);
    my_data->buffer_out[my_data->size_thr] = '0';
    char *p = my_data->buffer_out;
    int i = 0;
    int count;
    MPI_Status status;
    for (i = my_data->begin_index; i < my_data->end_index; i++){
        assert( strlen(my_data->seqs_thr[i].sam)> 0 );
        memmove(p, my_data->seqs_thr[i].sam, strlen(my_data->seqs_thr[i].sam) * sizeof(char));
        p += my_data->seqs_thr[i].l_seq;
        free(my_data->seqs_thr[i].sam);
    }
    //fprintf(stderr, "in copy_buffer_write_thr \n");
    //size_t res = pwrite(my_data->file_desc, my_data->buffer_out, my_data->size_thr, 0);
    //assert(res == my_data->size_thr);
      
    int res = MPI_File_write_shared(my_data->file_desc, my_data->buffer_out, my_data->size_thr, MPI_CHAR, &status);
    assert(res == MPI_SUCCESS);
    res = MPI_Get_count(&status, MPI_CHAR, &count);
    assert(res == MPI_SUCCESS);
    assert(count == (int)my_data->size_thr);
    free(my_data->buffer_out);
}
*/
void write_sam_mt(void *thread_arg){

    struct struct_data_thread *my_data;
    my_data = (struct struct_data_thread *) thread_arg;
    int i = 0;
    int j = 0;
    int count;
    MPI_Status status;
    size_t total_to_wrt = 0;
    size_t *l_seq = malloc((my_data->end_index - my_data->begin_index + 1)*sizeof(size_t));

    for (i = my_data->begin_index; i < my_data->end_index; i++)  {
        l_seq[j]= strlen(my_data->seqs_thr[i].sam);
        total_to_wrt += l_seq[j];
        j++;
    }    
    char *buffer_out = malloc( total_to_wrt + 1);
    buffer_out[total_to_wrt] = '0';
    char *p = buffer_out;
    j = 0;
    for (i = my_data->begin_index; i < my_data->end_index; i++){
   
        memmove(p, my_data->seqs_thr[i].sam, l_seq[j] * sizeof(char));
        p += l_seq[j];
        j++;
        free(my_data->seqs_thr[i].sam);
    }
 
    int res = MPI_File_write_shared(my_data->file_desc, buffer_out, total_to_wrt, MPI_CHAR, &status);
    assert(res == MPI_SUCCESS);
    res = MPI_Get_count(&status, MPI_CHAR, &count);
    assert(res == MPI_SUCCESS);
    assert(count == (int)total_to_wrt);

    free(buffer_out);   
    free(l_seq);

}



void pread_fastq_chunck(void *thread_arg ){

    struct struct_pread_fastq *my_data;
    my_data = (struct struct_pread_fastq *) thread_arg;


    int num_thread      = my_data->total_thread;
    int thread_id       = my_data->thread_id;
    size_t begin_offset = my_data->offset;
    char* buffer        = my_data->buffer;
    size_t total_size   = my_data->size;
    MPI_File fd        = my_data->fd;
    MPI_Status status;
    int res, count;

    size_t size2read = total_size / (size_t)num_thread;
    size_t rest      = total_size % (size_t)size2read;

    size_t offset_th_read   = begin_offset + (thread_id)*size2read;
    size_t offset_in_buffer = (thread_id)*size2read;

    size_t offset_thr = begin_offset + (thread_id)*size2read;
    
    if (thread_id == (num_thread - 1)) size2read +=rest;
   
    //pread(fd, buffer + offset_in_buffer, size2read, offset_thr);
    res = MPI_File_read_at(fd, offset_thr, buffer + offset_in_buffer, size2read, MPI_CHAR, &status);      
    assert(res == MPI_SUCCESS);
    res = MPI_Get_count(&status, MPI_CHAR, &count);
    assert(res == MPI_SUCCESS);
    assert(count == (int)size2read);

}





/*
void find_reads_size_and_offsets_mt(void *thread_arg){


    struct struct_data_thread_1 *my_data;
    my_data = (struct struct_data_thread_1 *) thread_arg;

    size_t offset_in_file           = my_data->offset_in_file_mt;
    size_t siz2read                 = my_data->size2read_mt;
    char   *file_to_read            = my_data->file_r1_mt;
    size_t *p_local_num_reads       = my_data->local_num_reads_mt;
    size_t *p_total_num_reads       = my_data->total_num_reads_mt;
    size_t **local_read_offsets     = my_data->local_read_offsets_mt;
    int    **local_read_size        = my_data->local_read_size_mt;
    size_t **local_read_bytes       = my_data->local_read_bytes_mt;
    int    proc_num                 = my_data->proc_num_mt;
    int    rank_num                 = my_data->rank_num_mt;
    int    thread_num               = my_data->thread_num_mt;


    //fprintf(stderr, "in find_reads_size_and_offsets_mt thread %d file =%s  \n", thread_num, file_to_read);
    //fprintf(stderr, "in find_reads_size_and_offsets_mt thread %d size2read =%zu  \n", thread_num, siz2read);
    //fprintf(stderr, "in find_reads_size_and_offsets_mt thread %d offset =%zu  \n", thread_num, offset_in_file);
    //fprintf(stderr, "in find_reads_size_and_offsets_mt rank %d thread %d step 0 \n", rank_num, thread_num);

    MPI_File  mpi_fd;
    int fd;
    int res;
    //res = MPI_File_open(MPI_COMM_WORLD, file_to_read, MPI_MODE_RDONLY, MPI_INFO_NULL, &mpi_fd);
    //assert(res==MPI_SUCCESS);
    fd = open(file_to_read, O_RDONLY);
    //fprintf(stderr, "in find_reads_size_and_offsets_mt rank %d thread %d step 1 \n", rank_num, thread_num);


    char *buffer_r;
    char *b, *r, *t, *e;
    size_t offset_end_buff;
    size_t pos_in_vect = 0;
    size_t lines = 0;
    size_t total_parsing = 0;
    size_t g=0;

    MPI_Datatype arraytype;
    MPI_Datatype arraytype0;
    MPI_Datatype arraytype_r1;
    MPI_Datatype arraytype_r2;

    *p_total_num_reads = 0;
    size_t read_buffer_sz = 0;
    if ( siz2read < DEFAULT_INBUF_SIZE ) read_buffer_sz = siz2read;
    else read_buffer_sz = DEFAULT_INBUF_SIZE;

     while (1){

        buffer_r = malloc(read_buffer_sz + 1);
        assert( buffer_r != NULL );
        buffer_r[read_buffer_sz] = '0';

        ssize_t res = pread(fd, buffer_r, read_buffer_sz, offset_in_file);
        assert( read_buffer_sz == res );
        //res = MPI_File_read_at(mpi_fd, (MPI_Offset)offset_in_file, buffer_r, read_buffer_sz, MPI_CHAR, MPI_STATUS_IGNORE);
        //assert(res == MPI_SUCCESS);
        assert(*buffer_r == '@');

        b = buffer_r;
        r = b + read_buffer_sz;

        if ( read_buffer_sz == DEFAULT_INBUF_SIZE){
            while (r-- != b){if (*r == '\n' && *(r+1) == '+') {r--; break;}}
            while (r-- != b){if (*r == '\n') break;}
            while (r-- != b){if (*r == '@') break;}
            r--;
            offset_end_buff = (r - b);
        }
        else
            offset_end_buff = (r - b);


        //fprintf(stderr, "in find_reads_size_and_offsets_mt rank %d thread %d step 2 \n", rank_num, thread_num);


        t = buffer_r ;
        e = buffer_r + offset_end_buff;
        lines = 0;
        while (t++ < e){if (*t == '\n') lines++;}

        *p_local_num_reads =  (lines/4);
        *p_total_num_reads += *p_local_num_reads;
        *local_read_size    = (int *)realloc(*local_read_size, sizeof(int) * (*p_total_num_reads));
        *local_read_offsets = (size_t *)realloc(*local_read_offsets, sizeof(size_t) * (*p_total_num_reads));
        *local_read_bytes   = (size_t *)realloc(*local_read_bytes, sizeof(size_t) * (*p_total_num_reads));

        assert( *local_read_offsets != NULL);
        assert( *local_read_size    != NULL);
        assert( *local_read_bytes   != NULL);

        t = buffer_r;

        int size=0;

        size_t lines2 = 0;
        size_t lines3 = 0;

        t = buffer_r;
        int last_size;
        int done = 0;
        size_t start_read_offset = 0;

        g = offset_in_file;

        int count = 0;

        while (t < e){

            if (lines3 == lines) break;
            assert( *t == '@');

            start_read_offset = g;
            while (*t != '\n'){ t++; g++;}
            t++;g++; lines3++;
            while (*t != '\n'){ size++; t++; g++;}
            t++;g++; lines3++;
            while (*t != '\n'){ t++; g++;}
            t++;g++; lines3++;
            while (*t != '\n'){ t++; g++;}
            lines3++;

 			(*local_read_offsets)[pos_in_vect]  = start_read_offset;
            (*local_read_bytes)[pos_in_vect]    = (g - start_read_offset) + 1;
            assert((*local_read_bytes)[pos_in_vect] != 0);
            (*local_read_size)[pos_in_vect]     = size;
            assert((*local_read_size)[pos_in_vect] != 0);

            size = 0;
            pos_in_vect++;
            t++;g++;

        }

        assert( lines == lines3 );
        total_parsing   += offset_end_buff;
        if (total_parsing == siz2read) {free(buffer_r); break;}
        if ((siz2read - total_parsing) < DEFAULT_INBUF_SIZE)
            read_buffer_sz = siz2read - total_parsing;
        else read_buffer_sz = DEFAULT_INBUF_SIZE;

        offset_in_file  += offset_end_buff + 1;
        free(buffer_r);
    }

    assert(total_parsing == siz2read);
    //MPI_File_close(&mpi_fd);
    close(fd);
}

void copy_local_read_info_mt(void *thread_arg){

    struct struct_data_thread_1 *my_data;
    my_data = (struct struct_data_thread_1 *) thread_arg;

    size_t *p_local_num_reads       = my_data->local_num_reads_mt;
    size_t *p_total_num_reads       = my_data->total_num_reads_mt;
    size_t **local_read_offsets     = my_data->local_read_offsets_mt;
    int    **local_read_size        = my_data->local_read_size_mt;
    size_t **local_read_bytes       = my_data->local_read_bytes_mt;
    int    proc_num                 = my_data->proc_num_mt;
    int    rank_num                 = my_data->rank_num_mt;
    int    thread_num               = my_data->thread_num_mt;



    size_t offset = my_data->previous_read_num;

    //fprintf(stderr, " rank %d thread %d copy %zu reads from %zu \n", rank_num, thread_num, *p_total_num_reads, offset);

    int *p      = my_data->local_read_size + offset;
    memmove(p, *local_read_size, *p_total_num_reads * sizeof(int));

    size_t *p1  = my_data->local_read_offsets + offset;
    memmove(p1, *local_read_offsets, *p_total_num_reads * sizeof(size_t));

    size_t *p2  = my_data->local_read_bytes + offset;
    memmove(p2, *local_read_bytes, *p_total_num_reads * sizeof(size_t));
}
*/

