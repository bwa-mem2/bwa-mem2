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

Authors:  Sanchit Misra <sanchit.misra@intel.com>; Vasimuddin Md <vasimuddin.md@intel.com>;
*****************************************************************************************/

#ifndef FASTA_FILE_H
#define FASTA_FILE_H

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<assert.h>

long long int read_multi_fasta_file_batch(FILE *input, int64_t first, int64_t last, int64_t cur, char* str_db)
{
        char * line = NULL;
        size_t len = 0;
        ssize_t read;

        if(input==NULL){
                return -1;
        }
        int64_t i = cur;

        char name[1000]="", content[1000]="";
		
        name[0]='\0',content[0]='\0';
        while( (read=getline(&line,&len,input))!=-1 ) {
            int line_len = strnlen(line, 1000); 
            assert(line_len < 1000);
			if( line_len==0 || line[0] == '>' ) { /* Identifier marker*/
				if( strlen(name)!=0 ){ 

					if(i >= first)
					{
						memcpy(str_db,content,(strlen(content)));
						str_db+=(strlen(content)-1);
					}
					content[0]='\0';
					i++;
					if(i == last)
						break;
					name[0]='\0';
				}
				if(( line_len !=0) ) {
					memcpy(name,&line[1], line_len);
				}
				content[0]='\0';
			} else if( strlen(name)!=0 ) {
				if( strstr(line," ") != NULL ){ /* Invalid sequence--no spaces allowed*/
					name[0]='\0';
					content[0]='\0';
				} else {
                    assert(strlen(content) + line_len < 1000);
					strcat(content,line);
				}
			}
        }
        if(strlen(content) > 0)
        {
            memcpy(str_db,content,strlen(content));
            i++;
        }
        return i - first;
}
long long int read_multi_fasta_file(char *file_path, char* str_db)
{
        FILE *input;
        char * line = NULL;
        size_t len = 0;
        ssize_t read;

        input=fopen(file_path,"r");

        if(input==NULL){
                return -1;
        }
        int i=0;
        char name[1000]="", content[1000]="";

        name[0]='\0',content[0]='\0';

        while( (read=getline(&line,&len,input))!=-1 ){
            int line_len = strnlen(line, 1000); 
            assert(line_len < 1000);			
			if(line_len==0 || line[0] == '>' ){ /* Identifier marker*/
				if( strlen(name)!=0 ){ 

					memcpy(str_db,content,(strlen(content)));
					str_db+=(strlen(content)-1);
					i++;
					name[0]='\0';
				}
				if(( line_len!=0) ){
					if (strlen(line) >= 1000) {
						printf("Error: fasta_file.h:92, name buffer of size 1000 overflowed!!\n");
						exit(0);
					}
					memcpy(name,&line[1], line_len);
				}
				content[0]='\0';
			} else if( strlen(name)!=0 ){
				if( strstr(line," ") != NULL ){ /* Invalid sequence--no spaces allowed*/
					name[0]='\0';
					content[0]='\0';
				} else {
					if (strlen(content) + line_len >= 1000) {
						printf("Error: fasta_file.h:108, name buffer of size 1000 overflowed!!\n");
						exit(0);
					}					
					strcat(content,line);
				}
			}
        }
        memcpy(str_db,content,strlen(content));
        i++;
		fclose(input);
        return i;
}
long long int count_fasta_file(char *file_path, char *str_db)
{
        /*assume file has just one sequence*/
        FILE *input;
        char * line = NULL;
        size_t len = 0;
        ssize_t read;

        input=fopen(file_path,"r");

        if(input==NULL){
                return -1;
        }
        long long int i=0;
        char name[1000], content[1000];
        name[0]='\0',content[0]='\0';

        while( (read=getline(&line,&len,input))!=-1 ){
			if( strlen(line)==0 || line[0] == '>' ){ /* Identifier marker*/
				if( strnlen(name, 1000)!=0 ){

					i+=strlen(content)-1;
					name[0]='\0';
				}
				if(( strlen(line)!=0) ) {
					if (strlen(line) > 1000) {
						printf("Error: fasta_file.h:137, name buffer of size 1000 overflowed!!\n");
						exit(0);
					}
					memcpy(name,&line[1],strlen(line)-1);
				}
				content[0]='\0';
			} else if( strnlen(name, 1000)!=0 ){
				if( strstr(line," ") != NULL ){ /* Invalid sequence--no spaces allowed*/
					name[0]='\0';
					content[0]='\0';
				} else {
					i+=strlen(line)-1;
				}
			}
        }
		fclose(input);
        return i;
}

#endif
