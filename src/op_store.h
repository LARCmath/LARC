//                              op_store.h
//                Matrix Operations Stores for LARC
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


#ifndef LARC_OP_STORE_H
#define LARC_OP_STORE_H

#include <stdint.h>
#include "larc.h"

/*!
 * \ingroup larc
 * \brief Produces a hash value from the three components of an operation
 * \param in1_pID The packedID of the first matrix/scalar
 * \param in2_pID The packedID of the second matrix/scalar
 * \param op_type The enum type for the operation to be stored
 * \return The desired hash value
 */
uint64_t
hash_from_op(uint64_t in1_pID, uint64_t in2_pID, op_type_t op_type);

/*!
 * \ingroup larc
 * \brief Returns enum value corresponding to operation, for hashing
 * \param op_name String value for operation, e.g. "SUM".
 * \return The integer enum value corresponding to the operation
 */
op_type_t get_op_type_from_string_name(char *op_name);

/*!
 * \ingroup larc
 * \brief Adds an operation to the operations store
 *
 * \param op The type of operation to be added
 * \param in1_pID packedID of the first matrix of the operation
 * \param in2_pID packedID of the second matrix of the operation
 * \param out_pID packedID of the matrix created by the operation
 * \param hash Identifies the hash chain to which this operation will be added
 *
 * \return 0 on success, or error code if failed 
 */
int op_set(op_type_t op, uint64_t in1_pID, uint64_t in2_pID, uint64_t out_pID,
	uint64_t hash);

/*!
 * \ingroup larc
 * \brief Retrieves the packedID of the result of an operation stored in the operations store
 *
 * \param op The type of the operation
 * \param in1_pID packedID of the first matrix of the operation
 * \param in2_pID packedID of the second matrix of the operation
 * \param hash Identifies the hash chain to be searched for the operation
 * \return The packedID for the stored scalar or matrix, or MATRIX_ID_INVALID if not found
 */
int64_t op_get(op_type_t op, uint64_t in1_pID, uint64_t in2_pID,
             uint64_t hash);

/*!
 * \ingroup larc
 * \brief Lists out all the op_names and description of them to help users interface with op chain functions
 */
void list_op_names();

#ifndef SWIG
/*!
 * \ingroup larc
 * \brief Creates the operations store
 *
 * \param exponent The base2 log of the number of hash chains in the operations store
 * \return 1 on successfully creating the operations store
 */
int create_op_store(size_t exponent);

/*!
 * \ingroup larc
 * \brief Frees memory allocated to the operations store
 *
 * This routine is called by shutdown_larc. It should not be called from
 * any other routine.
 */
void free_op_store(void);
#endif // #ifndef SWIG

/*!
 * \ingroup larc
 * \brief Retrieve the exponent for the op store
 * \return The exponent used to create the op store hash table
 */
size_t get_op_store_exp(void);

/*!
 * \ingroup larc
 * \brief Prints a summary of operations store usage
 * \param outfilepath The location for the report (can be "stdout")
 */
void op_store_report(char *outfilepath);

//int64_t op_hashID_by_matrixIDs(int64_t in1_pID, int64_t in2_pID, char *op_name);

/*!
 * \ingroup larc
 * \brief prints a hash chain from the opstore to a file
 *
 * \param hash The index into the operations store
 * \param outfilepath The file to which the report will be written (can be "stdout")
 * \param comment A user-specified string printed to the file
 * \return 1 on success, otherwise error code
 */
int op_hash_chain_info_to_file(uint64_t hash, char *outfilepath, char *comment);

/*!
 * \ingroup larc
 * \brief prints a hash chain from the opstore to screen
 * 
 * \param hash The index into the operations store
 * \param comment A user-specified string printed to the file
 * \return 1 on success, otherwise error code
 */
int op_hash_chain_info_to_screen(uint64_t hash, char *comment);

/*!
 * \ingroup larc
 * \brief Cleans a hash chain in the op store for a given hash
 *
 * This removes any links in the hash chain that refer to matrices that have
 * been removed from the matrix store.
 *
 * \param hashID The hash chain which will have
 * \return 1 on success, otherwise error code
 */
int clean_op_hash_chain(uint64_t hashID);

/*!
 * \ingroup larc
 * \brief Cleans all hash chains in the op store
 *
 * This removes any hash chain links in the operations store that refer to
 * matrices that have been removed from the matrix store.
 *
 * \return 1 on success, otherwise error code
 */
int clean_op_store();

/*!
 * \ingroup larc
 * \brief Empties all hash chains in the op store
 *
 * This removes all records from the operations store
 * \return 1 on success, otherwise error code
 */
int empty_op_store();

#ifdef HASHSTATS
/*!
 * \ingroup larc
 * \brief Prints hash table statistics for the operations store
 *
 * prints hash table statistics for op store, creating files with
 * accesses and nodes for each bucket (hash_val).  It also produces
 * as standard hash_report to stdout or file. This works only if
 * HASHSTATS is defined.
 *
 * \param accesses_file
 * \param nodes_file
 * \param report_file
 */
void op_hashstats(char *accesses_file,  char *nodes_file, char *report_file);
#endif // HASHSTATS

#endif
