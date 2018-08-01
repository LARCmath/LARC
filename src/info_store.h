//                            info_store.h
//                  Matrix Information Store for LARC
/******************************************************************
 *                                                                *
 * Copyright 2014, Institute for Defense Analyses                 *
 * 4850 Mark Center Drive, Alexandria, VA; 703-845-2500           *
 * This material may be reproduced by or for the US Government    *
 * pursuant to the copyright license under the clauses at DFARS   *
 * 252.227-7013 and 252.227-7014.                                 *
 *                                                                *
 * LARC : Linear Algebra via Recursive Compression                *
 * Authors:                                                       *
 *   - Steve Cuccaro (IDA-CCS)                                    *
 *   - John Daly (LPS)                                            *
 *   - John Gilbert (UCSB, IDA adjunct)                           *
 *   - Jenny Zito (IDA-CCS)                                       *
 *                                                                *
 * Additional contributors are listed in "LARCcontributors".      *
 *                                                                *
 * POC: Jennifer Zito <jszito@super.org>                          *
 * Please contact the POC before disseminating this code.         *
 *                                                                *
 *****************************************************************/


#ifndef LARC_INFO_STORE_H
#define LARC_INFO_STORE_H

#include <stdint.h>
#include "larc.h"


/* lists out all the info_names and description of them 
   to help users interface with info chain functions  */
void list_info_names();

// This function returns the name associated with enum info_type value
char *return_info_name(uint64_t info_type);

// return an array of the info_names to the user
// char** info_names_array();

/* Prints a summary of info_store usage */
void info_store_report(char *outfilepath);

/* Returns 1 on successfully creating the info_store */
int create_info_store(size_t exponent);

struct larc_info_t *info_find(info_type_t info_type, int64_t matID, int create_flag);

/* Returns 0 on success, or error code if failed  */
int info_set(info_type_t info_type, int64_t my_matID, const char *info_data);

/* Returns the matID previously stored, or MATRIX_ID_INVALID if not found */
char *info_get(info_type_t info_type, int64_t my_matID);

/* translates a string into the appropriate enum value */
info_type_t get_info_type_from_string_name(const char *info_name);

/* returns the hashID value associated with the info record */
int64_t info_hashID_by_matrixIDs(int64_t my_matID, const char *info_name);

// prints a hash chain from the infostore to a file.
int info_hash_chain_info_to_file(uint64_t hash, char *outfilepath, char *comment);

// prints a hash chain from the infostore to screen.
int info_hash_chain_info_to_screen(uint64_t hash, char *comment);

// Cleans a hash chain in the info store for a given hash
int clean_info_hash_chain(uint64_t hashID);


#endif
