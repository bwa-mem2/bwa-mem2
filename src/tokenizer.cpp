/*
   This file is part of mpiSORT
   
   Copyright Institut Curie 2020
   
   This software is a computer program whose purpose is to sort SAM file.
   
   You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
   
   The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND. Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data. 
   
   The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.
*/

/*
   Module:
     tokenizer.c

   Authors:
    Frederic Jarlier, 	Institut Curie
	Nicolas Joly, 		Institut Pasteur
	Nicolas Fedy,		Institut Curie
	Leonor Sirotti,	 	Institut Curie
	Thomas Magalhaes,	Institut Curie
	Paul Paganiban,		Institut Curie
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

#include "tokenizer.h"

// does NOT allocate memory : had to be done elsewhere
int tokenizer(char* str, const char delim, char* token){

	static char* last;
	int found = 0, length = 0;
	char* check;
	int i;

	if(str)
		last = str;

	else if (!last || *last == 0)
		return 0;

	check = last;

	while (check && !found && last[length]){

		if (*check == delim)
			found = 1;
		else{
			check++;
			length++;
		}
	}

	if (!found)
		return 0;

	for(i = 0; i < length; i++){
		token[i] = *last++;
	}

	token[length] = 0;
	last++;

	return 1;
}
