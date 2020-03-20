//                          io.h
/******************************************************************
 *                                                                *
 * Copyright (C) 2014, Institute for Defense Analyses             *
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
 * Questions: larc@super.org                                      *
 *                                                                *
 * All rights reserved.                                           *
 *                                                                *
 * Redistribution and use in source and binary forms, with or     *
 * without modification, are permitted provided that the          *
 * following conditions are met:                                  *
 *   - Redistribution of source code must retain the above        *
 *     copyright notice, this list of conditions and the          *
 *     following disclaimer.                                      *
 *   - Redistribution in binary form must reproduce the above     *
 *     copyright notice, this list of conditions and the          *
 *     following disclaimer in the documentation and/or other     *
 *     materials provided with the distribution.                  *
 *   - Neither the name of the copyright holder nor the names of  *
 *     its contributors may be used to endorse or promote         *
 *     products derived from this software without specific prior *
 *     written permission.                                        *
 *                                                                *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND         *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,    *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF       *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE       *
 * DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER NOR        *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,   *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT   *
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;   *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)       *
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN      *
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR   *
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, *
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.             *
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


/*!
 * \ingroup larc
 *
 * \brief Reads a file starting with the dimensions, followed by all entries.
 *
 * The entries are listed in the order you would get by reading each row in
 * turn
 *
 * \param file_path The path to the file to be read.
 * \return A pointer to the stored matrix.
 */
mat_ptr_t read_row_major_matrix_from_file ( char * file_path); 

/*!
 * \ingroup larc
 *
 * \brief Take a row major list descibing a matrix and put it in the store.
 *
 * This is a recursive routine, and requires the original number of columns
 * (the paramater dim_whole) at lower levels of the recursion.
 *
 * \param dense_mat A pointer to scalarType; the matrix to be stored
 * \param row_level The base2 log of the number of rows in the current matrix
 * \param col_level The base2 log of the number of columns in the current matrix
 * \param dim_whole The number of columns of the matrix originally passed to the recursion
 * \return A pointer to the stored matrix.
 */
mat_ptr_t row_major_list_to_store(scalarType * dense_mat, 
				  mat_level_t    row_level, 
				  mat_level_t    col_level, 
				  int64_t        dim_whole 
				  );

/*!
 * \ingroup larc
 * \brief Stores a formatted compressed LARCmatrix to a json file and applies
 * func to scalar entries before they are written.
 *
 * The function func() is applied to each scalar as it is written out to the
 * json file. The value of the scalar stored in the matrix store is unchanged.
 *
 * \param m_ptr The matrix pointer for the matrix to be written.
 * \param path The file to which the matrix will be written.
 * \param func A function that will be applied to each scalar before writing
 * \return A value indicating success or the type of error.
 */
int write_and_alter_vals_larcMatrix_file(mat_ptr_t m_ptr, char *path, void (*func)(scalarType*, const scalarType));

/*!
 * \ingroup larc
 * \brief Reads a LARCMatrix formatted compressed matrix from a json file and
 * applies func to scalar entries as they are read in.
 *
 * The function func() is applied to each scalar in the file before passing
 * it to the routines which put it into the matrix store (or find that the
 * modified value is already there...)
 *
 * \param path The filename of the LARCMatrix file to be read into LARC.
 * \param func A function that will be applied to each scalar before storing it
 * \return The pointer to the stored matrix.
 */
mat_ptr_t read_and_alter_vals_larcMatrix_file_return_matPTR(char *path, void (*func)(scalarType*, const scalarType));

#ifndef SWIG
/*!
 * \ingroup larc
 * \brief Probably-obsolete function which reads in a matrix from a LARCMatrix json file, turning small scalars to zero before storing the matrix
 * \param path The location of the json file to be read
 * \param threshold The smallest scalarType value to not be zeroed
 * \return The pointer to the stored matrix
 */
mat_ptr_t matrix_read_larcMatrix_file_zeroize_small_scalars(char *path, scalarType threshold);

/*!
 * \ingroup larc
 * \brief Probably-obsolete function which reads in a matrix from a LARCMatrix json file and performs a flattening function to the data before storage
 * \param path The location of the json file to be read
 * \param threshold 
 * \param flatValue 
 * \return The pointer to the stored matrix
 */
mat_ptr_t matrix_read_larcMatrix_file_flatten_small_scalars(char *path, scalarType threshold, scalarType flatValue);
#endif // ifndef SWIG

/*!
 * \ingroup larc
 * \brief Prints every entry of a SMALL matrix to screen by address
 * 
 * \param A_ptr The pointer to a matrix to be printed to the screen.
 */
void print_naive_by_matPTR(mat_ptr_t A_ptr);

/*!
 * \ingroup larc
 * \brief Stores a LARCMatrix formatted compressed matrix to a json file
 * 
 * \param id The pointer to a matrix to be stored.
 * \param path The filename the matrix will be written to.
 * \return A value indicating success or the type of error.
 */
int write_larcMatrix_file_by_matPTR(mat_ptr_t id, char *path);

/*!
 * \ingroup larc
 * \brief Stores a SMALL uncompressed matrix to a file
 * 
 * \param id The pointer to a matrix to be stored.
 * \param path The filename the matrix will be written to.
 */
void write_naive_by_matPTR(mat_ptr_t id, char *path);

/*!
 * \ingroup larc
 * \brief Stores the nonzero elements of a matrix and their positions to a file.
 * 
 * This routine should only be used for SMALL or VERY SPARSE matrices
 *
 * \param id The pointer to a matrix to be stored.
 * \param path The filename the matrix will be written to.
 */
void write_matrix_nonzeros_by_matPTR(mat_ptr_t id, char *path);  

/*!
 * \ingroup larc
 * \brief Reads a LARCMatrix formatted compressed matrix from a json file 
 * 
 * \param path The json file from which the matrix will be read.
 * \return The pointer to the stored matrix.
 */
mat_ptr_t read_larcMatrix_file_return_matPTR(char *path);

/*!
 * \ingroup larc
 * \brief Reads a LARCMatrix formatted compressed matrix from a json file 
 * 
 * This routine assumes the LARCMatrix file to be in an older format, in which
 * the scalar values are not stored as strings. It is in LARC to maintain
 * backwards compatibility.
 *
 * \param path The json file from which the matrix will be read.
 * \return The pointer to the stored matrix.
 */
mat_ptr_t read_larcMatrix_file_legacy_return_matPTR(char *path);

/*!
 * \ingroup larc
 * \brief Python interface version of read_row_major_matrix_from_file
 * 
 * \param file_path The file from which the matrix will be read.
 * \return The matrixID of the stored matrix.
 */
int64_t read_row_major_matrix_from_file_matrixID(char * file_path);

/*!
 * \ingroup larc
 * \brief Python interface version of row_major_list_to_store
 * 
 * The parameters for this routine mirror those of row_major_list_to_store
 * (somewhat unnecessarily, since the _matrixID routine is not recursive).
 *
 * \param dense_mat The matrix to be stored, with character strings as elements
 * \param current_row_level The row level of the matrix
 * \param current_col_level The column level of the matrix
 * \param orig_num_cols The number of columns of the matrix
 *
 * \return The matrixID of the stored matrix.
 */
int64_t row_major_list_to_store_matrixID(char **dense_mat, 
				       mat_level_t        current_row_level, 
				       mat_level_t        current_col_level, 
				       int64_t        orig_num_cols 
				       );

/*!
 * \ingroup larc
 * \brief Python Interface: Prints every entry of a SMALL matrix to screen
 *
 * \param A_mID The matrix ID of the matrix to be printed to the screen.
 */
void print_naive_by_matID(int64_t A_mID);

/*!
 * \ingroup larc
 * \brief Python Interface: Stores a LARCMatrix formatted compressed matrix to a json file
 *
 * \param m_mID The matrixID of the matrix to be stored in compressed format.
 * \param path The path to the file that will be created.
 * \return A value indicating success or the type of error.
 */
int write_larcMatrix_file_by_matID(int64_t m_mID, char *path);

/*!
 * \ingroup larc
 * \brief Python interface: Prints every entry of a SMALL matrix to a file
 * 
 * \param id The matrixID of the matrix to be stored (uncompressed!)
 * \param path The path to the file that will be created.
 */
void write_naive_by_matID(int64_t id, char *path);

/*!
 * \ingroup larc
 * \brief Python interface: Prints nonzero entries of a SMALL or VERY SPARSE matrix to a file
 *
 * \param id The matrixID of the matrix to be stored.
 * \param path The path to the file that will be created.
 */
void write_matrix_nonzeros_by_matID(int64_t id, char *path);  

/*!
 * \ingroup larc
 * \brief Python interface: Reads a larcMatrix formatted compressed matrix from a json file
 *
 * \param path The path to the json file to be read into the matrix store
 * \return The matrixID of the stored file
 */
int64_t read_larcMatrix_file_return_matID(char *path);

/*!
 * \ingroup larc
 * \brief Python interface: Reads a LARCMatrix formatted compressed matrix from a json file
 *
 * This routine assumes the LARCMatrix file to be in an older format, in which
 * the scalar values are not stored as strings. It is in LARC to maintain
 * backwards compatibility.
 *
 * \param path The path to the json file to be read into the matrix store
 * \return The matrixID of the stored file
 */
int64_t read_larcMatrix_file_legacy_return_matID(char *path);

/*!
 * \ingroup larc
 * \brief Python interface: checks two LARCMatrix json files to see if they contain the same matrix
 *
 * The compressed LARCMatrix file format refers to matrices by matrixID, which even
 * for identical matrices may vary depending on the order in which matrices
 * were stored for a particular program run. This routine reads both matrices
 * into the same matrix store to see if the returned matrixIDs are the same.
 *
 * \param path1 The path to the first json file
 * \param path2 The path to the second json file
 * \return The matrixID for the matrices if they are the same, 0 otherwise
 */
int64_t equal_matrices_in_larcMatrix_files(char *path1, char *path2);



// THE FOLLOWING FUNCTIONS NEED DOXYGEN HEADERSL

int64_t matrix_read_larcMatrix_file_zeroize_small_scalars_matrixID(char *path, char *threshold);

int64_t matrix_read_larcMatrix_file_flatten_small_scalars_matrixID(char *path, char *threshold, char *flatValue);

mat_ptr_t read_matrixMarketExchange_file_return_matPTR(char * file_path);

size_t write_larcMatrix_file_return_larcSize(mat_ptr_t m_ptr, char *path,
							 void (*func)(scalarType*,
							 const scalarType));


size_t write_larcMatrix_file_return_larcSize_by_matID(int64_t m_mID, char *path);

size_t write_larcMatrix_file_return_larcSize_by_matPTR(mat_ptr_t m_ptr, char *path);


#endif
