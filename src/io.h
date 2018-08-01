//                          io.h
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


#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "json.h"
#include "larc.h"
#include "info_store.h"
#include "matrix_store.h"


/* Reads a file starting with the dimensions, followed by all entries 
   The entries are listed in the order you would get by reading each row in turn */
mat_add_t read_row_major_matrix_from_file ( char * file_path); 

/* Python interface version of read_row_major_matrix_from_file */
int64_t read_row_major_matrix_from_file_matrixID(char * file_path);



// Take a row major list descibing a matrix and put it in the store.
mat_add_t row_major_list_to_store(ScalarType * dense_mat, 
				  mat_level_t    row_level, 
				  mat_level_t    col_level, 
				  int64_t        dim_whole 
				  );

// Python interface version of  row_major_list_to_store
int64_t row_major_list_to_store_matrixID(ScalarType  * dense_mat, 
				       mat_level_t        current_row_level, 
				       mat_level_t        current_col_level, 
				       int64_t        orig_num_cols 
				       );


/* Prints every entry of a SMALL matrix to screen by address*/
void print_matrix_naive(mat_add_t A_ptr);

/* Python Interface: Prints every entry of a SMALL matrix to screen by matrixID */
void print_matrix_naive_by_matrixID(int64_t A_mID);



/* Stores a json formatted compressed matrix to a file */
int matrix_write_json_file(mat_add_t id, char *path);

/* Python Interface: Stores a json formatted compressed matrix to a file */
int matrix_write_json_file_matrixID(int64_t m_mID, char *path);

/* Prints every entry of a SMALL matrix to a file */
void print_matrix_to_file_naive(mat_add_t, char *);

/* Python interface: Prints every entry of a SMALL matrix to a file */
void print_matrix_to_file_naive_by_matrixID(int64_t, char *);



/* Prints nonzero entries of a SMALL or VERY SPARSE matrix to a file */
void print_matrix_nonzeros_to_file(mat_add_t, char *);  

/* Python interface: Prints nonzero entries of a SMALL or VERY SPARSE matrix to a file */
void print_matrix_nonzeros_to_file_by_matrixID(int64_t, char *);  



/* Reads a json formatted compressed matrix from a file */
mat_add_t matrix_read_json_file(char *path);

/* Python interface: Reads a json formatted compressed matrix from a file, returns matrixID */
int64_t matrix_read_json_file_matrixID(char *path);

/* Python interface: Reads a two json formatted compressed matrices
   from files, returns the matrixID if they are equal 
   matrices, and 0 otherwise */
int64_t equal_matrices_in_json_files(char *path1, char *path2);


mat_add_t matrix_read_json_file_zeroize_small_scalars(char *path, ScalarType threshold);

int64_t matrix_read_json_file_zeroize_small_scalars_matrixID(char *path, ScalarType threshold);




mat_add_t matrix_read_json_file_flatten_small_scalars(char *path, ScalarType threshold, ScalarType flatValue);

int64_t matrix_read_json_file_flatten_small_scalars_matrixID(char *path, ScalarType threshold,  ScalarType flatValue);


#endif
