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
 *   - Mark Pleszkoch (IDA-CCS)                                   *
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
#include "scalars.h"
#include "info_store.h"
#include "matrix_store.h"

/*!
 * \ingroup larc
 * \brief Reads a matrix represented in row major format from a file and adds it to the MatrixStore
 * 
 * \param file_path The file from which the matrix will be read.
 * \return The packedID of the stored matrix.
 */
int64_t read_row_major_matrix_from_file(char * file_path);

/*!
 * \ingroup larc
 * \brief Adds an array of values (in string form) held in memory to the MatrixStore
 * 
 * The parameters for this routine mirror those of
 * recursive_row_major_list_to_store except that the matrix is made up of
 * character strings rather than scalarType values
 *
 * \param dense_mat The matrix to be stored, with character strings as elements
 * \param current_row_level The row level of the matrix
 * \param current_col_level The column level of the matrix
 * \param orig_num_cols The number of columns of the matrix
 *
 * \return The packedID of the stored matrix.
 */
int64_t row_major_list_to_store(char **dense_mat, 
				       mat_level_t  current_row_level, 
				       mat_level_t  current_col_level, 
				       int64_t  orig_num_cols );

/*!
 * \ingroup larc
 * \brief Prints every entry of a SMALL matrix to screen
 *
 * \param A_pID The packedID of the matrix to be printed to the screen.
 */
void print_naive(int64_t A_pID);

/*!
 * \ingroup larc
 * \brief Prints every entry of a SMALL matrix to a file
 * 
 * \param A_pID The packedID of the matrix to be stored (uncompressed!)
 * \param path The path to the file that will be created.
 */
void fprint_naive(int64_t A_pID, char *path);

/*!
 * \ingroup larc
 * \brief Prints nonzero entries of a SMALL or VERY SPARSE matrix to a file
 *
 * \param A_pID The packedID of the matrix to be stored.
 * \param path The path to the file that will be created.
 */
void fprint_matrix_nonzeros(int64_t A_pID, char *path);  

/*!
 * \ingroup larc
 * \brief Reads a larcMatrix formatted compressed matrix from a json file and adds it to the MatrixStore
 *
 * This routine assumes that the LARCMatrix file is line-by-line,
 * character-by-character equivalent
 * to what the fprint_larcMatrixFile() function produces.
 *
 * \param path The path to the json file to be read into the MatrixStore
 * \return The packedID of the stored file
 */
int64_t read_larcMatrixFile(char *path);

/*!
 * \ingroup larc
 * \brief Reads a larcMatrix formatted compressed matrix from a json file and adds it to the MatrixStore
 *
 * This routine can handle a JSON file that is *not* line-by-line,
 * character-by-character equivalent
 * to what the fprint_larcMatrixFile() function produces. Because of its
 * generality, it requires substantially more resources than the default
 * read_larcMatrixFile() routine.
 *
 * \param path The path to the json file to be read into the MatrixStore
 * \return The packedID of the stored file
 */
int64_t read_anyFormat_larcMatrixFile(char *path);

/*!
 * \ingroup larc
 * \brief Reads a LARCMatrix formatted compressed matrix from a json file and adds it to the MatrixStore
 *
 * This routine assumes the LARCMatrix file to be in an older format, in which
 * the scalar values are not stored as strings. It is in LARC to maintain
 * backwards compatibility. Like the read_anyFormat_larcMatrixFile() routine,
 * it uses a general JSON reader and may require substantially more resources
 * than one would expect for the read_larcMatrixFile() routine.
 *
 * \param path The path to the json file to be read into the MatrixStore
 * \return The packedID of the stored file
 */
int64_t read_legacy_larcMatrixFile(char *path);

/*!
 * \ingroup larc
 * \brief checks two LARCMatrix json files to see if they contain the same matrix
 *
 * The legacy compressed LARCMatrix file format refers to matrices by matrixID,
 * which even for identical matrices may vary depending on the order in
 * which matrices were stored for a particular program run. This routine reads
 * both matrices into the same MatrixStore to see if the returned packedIDs
 * are the same.
 *
 * If the files were written in the newer sequentially-numbered format, then
 * they will have identical matrix records (there may be some differences in
 * InfoStore fields) and this program is not needed.
 *
 * \param path1 The path to the first json file
 * \param path2 The path to the second json file
 * \return The packedID for the matrices if they are the same, 0 otherwise
 */
int64_t equal_matrices_in_larcMatrix_files(char *path1, char *path2);

/*!
 * \ingroup larc
 * \brief Reads a matrix represented in Matrix Market Exchange format from a file and adds it to the MatrixStore
 * 
 * \param file_path The file from which the matrix will be read.
 * \return The packedID of the stored matrix.
 */
int64_t read_matrixMarketExchange_file(char * file_path);

/*!
 * \ingroup larc
 * \brief Writes a LARCmatrix to a file and returns the LARCsize
 *
 * \param m_pID The packedID of the matrix to be written in LARC format
 * \param path The location of the new larcMatrix json file
 * \return The larcSize of the matrix with packedID m_pID
 */
size_t fprint_larcMatrixFile(int64_t m_pID, char *path);


/*!
 * \ingroup larc
 * \brief Writes a file with unique scalars in a matrix and a file with infostore metadata.
 *
 * \param m_pID The packedID whose scalars will be listed in one file, and metadata in other.
 * \param path File for list of scalars is path. File for meta data has ".info" suffix to path.
 * \return The number of unique scalars in the matrix.
 */
size_t fprint_uniqScalar_file(int64_t m_pID, char *path);



#endif
