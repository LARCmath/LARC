//                       hash.h
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

#ifndef LARC_HASH_H
#define LARC_HASH_H

/******************************************************************************
 *     Hash functions for Global Recursive Compression Scheme                 *
 *     for rapid find and insert functions for both the matrix and            *
 *     operations stores                                                      *
 *     last updated:   Fall 2015                                              *
 ******************************************************************************/

#include <stdint.h> // uint64_t
#include "matrix_store.h"
#include "mar.h"

#ifndef SWIG
/*!
 * \ingroup larc
 * \brief Performs the initial setup of a hash table for LARC
 * \param exponent The log base 2 of the size of the hash table
 * \return A pointer to the hash table
 */
hash_table_t *alloc_hash(int exponent);

/*!
 * \ingroup larc
 * \brief Reverses the allocations of alloc_hash
 * \param table A pointer to the hash table to be freed
 */
void dealloc_hash(hash_table_t *table);

/*!
 * \ingroup larc
 * \brief Creates and initializes a hash node 
 * \return A pointer to the hash node
 */
hash_node_t *alloc_hash_node();

/*
 * lookup/insert/delete an entry
 */
/*!
 * \ingroup larc
 * \brief This function currently does nothing
 * \param table A pointer to a hash table
 * \param record_ptr A record supposedly in this hash table
 * \param packedID A packedID supposedly in this hash table
 *
 * If the hash chain is in the matrixStore or ScalarStore, the data passed must
 * be a packedID (with record_ptr == RECORD_PTR_INVALID); if the hash chain is
 * in the infoStore or operationsStore, the data passed must
 * be a record_ptr (with packedID == MATRIX_ID_INVALID).
 *
 * \param hash_val Indicates the hash chain to be searched
 * \return A NULL pointer
 */
hash_node_t *hash_lookup(hash_table_t *table, record_ptr_t record_ptr,
        int64_t packedID, uint64_t hash_val);

/*!
 * \ingroup larc
 * \brief Inserts a hash node at the tail of a specified hash chain
 * \param table A pointer to a hash table
 * \param record_ptr A record to be added to a hash chain
 * \param packedID A packedID to be added to a hash chain
 *
 * If the hash chain is in the matrixStore or ScalarStore, the data passed must
 * be a packedID (with record_ptr == RECORD_PTR_INVALID); if the hash chain is
 * in the infoStore or operationsStore, the data passed must
 * be a record_ptr (with packedID == MATRIX_ID_INVALID).
 *
 * \param hash_val Indicates the hash chain to be modified
 * \return A pointer to the new hash node
 */
hash_node_t *hash_insert_at_tail(hash_table_t *table, record_ptr_t record_ptr,
        int64_t packedID, uint64_t hash_val);

/*!
 * \ingroup larc
 * \brief Inserts a hash node at the head of a specified hash chain
 * \param table A pointer to a hash table
 * \param record_ptr A record to be added to a hash chain
 * \param packedID A packedID to be added to a hash chain
 *
 * If the hash chain is in the matrixStore or ScalarStore, the data passed must
 * be a packedID (with record_ptr == RECORD_PTR_INVALID); if the hash chain is
 * in the infoStore or operationsStore, the data passed must
 * be a record_ptr (with packedID == MATRIX_ID_INVALID).
 *
 * \param hash_val Indicates the hash chain to be modified
 * \return A pointer to the new hash node
 */
hash_node_t *hash_insert_at_head(hash_table_t *table, record_ptr_t record_ptr,
        int64_t packedID, uint64_t hash_val);

/*!
 * \ingroup larc
 * \brief Removes a hash node from a specified hash chain, assuming that the data passed to the routine matches what is in the hash node
 *
 * If the hash chain is in the matrixStore or ScalarStore, the data passed must
 * be a packedID (with record_ptr == RECORD_PTR_INVALID); if the hash chain is
 * in the infoStore or operationsStore, the data passed must
 * be a record_ptr (with packedID == MATRIX_ID_INVALID).
 *
 * \param table A pointer to a hash table
 * \param record_ptr A record to be removed from a hash chain
 * \param packedID A packedID to be removed from a hash chain
 * \param hash_val Indicates the hash chain to be modified
 * \return 1 if the hash node was successfully removed, 0 otherwise
 */
int hash_node_remove(hash_table_t *table, record_ptr_t record_ptr,
        int64_t packedID, uint64_t hash_val);

/*!
 * \ingroup larc
 * \brief Removes a specified hash node from the hash table
 * \param table A pointer to a hash table
 * \param n A record to be removed from a hash chain
 * \param hash_val Indicates the hash chain to be modified
 * \return the node that followed n in the hash chain
 */
hash_node_t * hash_node_remove_node(hash_table_t *table, hash_node_t *n, uint64_t hash_val);

/*
 * report out
 */
/*!
 * \ingroup larc
 * \brief Prints a report of the current contents of a hash table
 * \param table A pointer to a hash table
 * \param fp A file pointer indicating where the output of this routine will be
 * \param verbose Currently not used
 * \return The size of the hash table in bytes
 */
uint64_t hash_report(hash_table_t *table, FILE *fp, int verbose);

/*!
 * \ingroup larc
 * \brief Returns the head node of a hash chain in a hash table
 * \param table A hash table
 * \param hash_val The hash value for the desired hash chain
 * \return The head node of the specified hash chain
 */
inline hash_node_t * hash_get_chain(hash_table_t *table, uint64_t hash_val)
{
  return table->heads[hash_val&table->mask];
}


#ifndef MAR   // SPRmode
/*!
 * \ingroup larc
 * \brief Obtains a hash value obtained from a single scalar value 
 * \param scalar The value to be hashed
 * \param exponent The log base 2 of the size of the hash table
 * \return The hash value
 */
uint64_t
region_hash_from_scalarType(scalarType scalar, uint64_t exponent);
#endif  // SPRmode

/*!
 * \ingroup larc
 * \brief Obtains a hash value from the four packedIDs which make up a matrix
 * \param sub An array of four packedIDs
 * \param mat_type One of (MATRIX, ROW_VECTOR, COL_VECTOR)
 * \param exponent The log base 2 of the size of the hash table
 * \return The hash value
 */
uint64_t
hash_from_matrix_subMatList(int64_t sub[4], matrix_type_t mat_type, uint64_t exponent);
#endif

/*!
 * \ingroup larc
 * \brief Hashes a single integer to the range [0,2**exponent)
 * \param key A 64-bit unsigned integer to be hashed
 * \param exponent Determines the size of the hash table (2**exponent)
 * \return The hash value
 */
uint64_t mult_golden_hash(uint64_t key, uint64_t exponent);

/*!
 * \ingroup larc
 * \brief Returns a hash value between 0 and 2^exponent-1
 * \param int1 One 64-bit unsigned integer to be hashed
 * \param int2 Another 64-bit unsigned integer to be hashed
 * \param exponent The log base 2 of the size of the hash table
 * \return The hash value
 */
uint64_t 
recursive_hash_from_two_integers
   (uint64_t int1, uint64_t int2, uint64_t exponent);

/*!
 * \ingroup larc
 * \brief Returns a hash value between 0 and 2^exponent-1
 * \param int1 One 64-bit integer to be hashed
 * \param int2 Another 64-bit integer to be hashed
 * \param int3 A third 64-bit integer to be hashed
 * \param exponent The log base 2 of the size of the hash table
 * \return The hash value
 */
uint64_t 
recursive_hash_from_three_integers
   (uint64_t int1, uint64_t int2, uint64_t int3, uint64_t exponent);

/*!
 * \ingroup larc
 * \brief Returns a hash value between 0 and 2^exponent-1
 * \param int1 One 64-bit integer to be hashed
 * \param int2 Another 64-bit integer to be hashed
 * \param int3 A third 64-bit integer to be hashed
 * \param int4 A fourth 64-bit integer to be hashed
 * \param exponent The log base 2 of the size of the hash table
 * \return The hash value
 */
uint64_t 
recursive_hash_from_four_integers 
   (uint64_t int1, uint64_t int2, uint64_t int3, uint64_t int4, uint64_t exponent);

/*!
 * \ingroup larc
 * \brief Returns a hash value between 0 and 2^exponent-1
 * \param tlist A list of 64-bit integers to be recursively hashed
 * \param tlength The number of values in that list
 * \param exponent The log base 2 of the size of the hash table
 * \return The hash value
 */
uint64_t
recursive_hash_from_int_list (uint64_t *tlist, uint64_t tlength, uint64_t exponent);

/*!
 * \ingroup larc
 * \brief Hashes one long double to an integer in range [0,2**exponent)
 * \param ldoub The long-double-precision number
 * \param exponent Determines the size of the hash table (2**exponent)
 * \return The hash value
 */
uint64_t
hash_from_one_longdouble (long double ldoub, uint64_t exponent);

/*!
 * \ingroup larc
 * \brief Hashes two long doubles to an integer in range [0,2**exponent)
 * \param d1l The first long-double-precision number
 * \param d2l The first long-double-precision number
 * \param exponent Determines the size of the hash table (2**exponent)
 * \return The hash value
 */
uint64_t
hash_from_two_longdoubles (long double d1l, long double d2l, uint64_t exponent);


#ifndef SWIG

#ifdef HASHSTATS
// hash table statistics produced when HASHSTATS is defined 
/*!
 * \ingroup larc
 * \brief Outputs statistics about a hash table to specified files
 * \param hash_table A pointer to a hash table
 * \param accesses_file The file reporting on hash chain accesses
 * \param nodes_file The file reporting on the number of nodes in the hash chains
 * \param report_file The file containing the hash report
 */
void hashstats_to_files( hash_table_t *hash_table, 
                         char *accesses_file, 
                         char *nodes_file, 
                         char *report_file
			);
#endif   // HASHSTATS


#ifdef MAR   // MARmode
/*!
 * \ingroup larc
 * \brief Returns a hash table index into the matrix store calculated from the MAR tile index
 * \param target_tile_index_PTR  The structure
 */
uint64_t  hash_tile_index(MAR_tile_index_t *target_tile_index_PTR);

/*!
 * \ingroup larc
 * \brief The hash filter is a short hash function that is independent
 *        of the standard hash into MatrixStore.
 * \param target_tile_index_PTR  The structure
 * \param hashfilter_bits  Number of bits for hashfilter
 */
uint64_t  hashfilter_of_tileindex(const MAR_tile_index_t *target_tile_index_PTR, size_t hashfilter_bits);

#endif   // MARmode


#endif  // SWIG


#endif  // 
