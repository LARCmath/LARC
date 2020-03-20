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

/*
 * allocate a new hash table
 */
/*!
 * \ingroup larc
 * \brief Performs the initial setup of a hash table for LARC
 * \param exponent The log base 2 of the size of the hash table
 * \return A pointer to the hash table
 */
hash_table_t *alloc_hash(int exponent);

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
 * \param hash_val Indicates the hash chain to be searched
 * \return A NULL pointer
 */
hash_node_t *hash_lookup(hash_table_t *table, record_ptr_t record_ptr, uint64_t hash_val);

/*!
 * \ingroup larc
 * \brief Inserts a hash node at the tail of a specified hash chain
 * \param table A pointer to a hash table
 * \param record_ptr A record to be added to a hash chain
 * \param hash_val Indicates the hash chain to be modified
 * \return A pointer to the new hash node
 */
hash_node_t *hash_insert_at_tail(hash_table_t *table, record_ptr_t record_ptr, uint64_t hash_val);

/*!
 * \ingroup larc
 * \brief Inserts a hash node at the head of a specified hash chain
 * \param table A pointer to a hash table
 * \param record_ptr A record to be added to a hash chain
 * \param hash_val Indicates the hash chain to be modified
 * \return A pointer to the new hash node
 */

hash_node_t *hash_insert_at_head(hash_table_t *table, record_ptr_t record_ptr, uint64_t hash_val);

/*!
 * \ingroup larc
 * \brief Removes a hash node from a specified hash chain, assuming that the data passed to the routine matches what is in the hash node
 * \param table A pointer to a hash table
 * \param record_ptr A record to be removed from a hash chain
 * \param hash_val Indicates the hash chain to be modified
 * \return 1 if the hash node was successfully removed, 0 otherwise
 */
int hash_node_remove(hash_table_t *table, record_ptr_t record_ptr, uint64_t hash_val);

/*!
 * \ingroup larc
 * \brief Removes a specified hash node from the hash table
 * \param table A pointer to a hash table
 * \param n A record to be removed from a hash chain
 * \param hash_val Indicates the hash chain to be modified
 */
void hash_node_remove_node(hash_table_t *table, hash_node_t *n, uint64_t hash_val);

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

/*
 * Utility
 */
/*!
 * \ingroup larc
 * \brief Returns the size of a hash table
 * \param table The hash table to be sized
 * \return The size of the hash table
 */
static inline uint64_t
hash_max(hash_table_t *table)
{
  return 1<<table->exponent;
}

// creates a mask that grabs the low exponent number of bits
// of the hash_val and then references that entry in the heads array
/*!
 * \ingroup larc
 * \brief Returns the head node of a hash chain in a hash table
 * \param table A hash table
 * \param hash_val The hash value for the desired hash chain
 * \return The head node of the specified hash chain
 */
static inline hash_node_t *
hash_get_chain(hash_table_t *table, uint64_t hash_val)
{
  return table->heads[hash_val&((1<<table->exponent)-1)];
}


#ifndef SWIG
/* Returns a 64 bit hash value, where the input is of scalarType and
 * exponent is the hash table exponent */
/*!
 * \ingroup larc
 * \brief Obtains a hash value obtained from a single scalar value 
 * \param scalar The value to be hashed
 * \param mat_type Should be SCALAR
 * \param exponent The log base 2 of the size of the hash table
 * \return The hash value
 */
uint64_t
hash_from_matrix_scalar(scalarType scalar, matrix_type_t mat_type, uint64_t exponent);
#endif

/* Returns a 64 bit hash value, where the input is an array of 4 mat_ptr_t 
 * and exponent is the hash table exponent */
/*!
 * \ingroup larc
 * \brief Obtains a hash value from the four pointers which make up a matrix panel
 * \param sub An array of four matrix pointers
 * \param mat_type One of (MATRIX, ROW_VECTOR, COL_VECTOR)
 * \param exponent The log base 2 of the size of the hash table
 */
uint64_t
hash_from_matrix_panel(mat_ptr_t sub[4], matrix_type_t mat_type, uint64_t exponent);

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
 * \brief Hashes one double to an integer in range [0,2**exponent)
 * \param doub The double-precision number
 * \param exponent Determines the size of the hash table (2**exponent)
 * \return The hash value
 */
uint64_t
hash_from_one_double (double doub, uint64_t exponent);

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
 * \brief Hashes two doubles to an integer in range [0,2**exponent)
 * \param d1 The first double-precision number
 * \param d2 The first double-precision number
 * \param exponent Determines the size of the hash table (2**exponent)
 * \return The hash value
 */
uint64_t
hash_from_two_doubles (double d1, double d2, uint64_t exponent);

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

#endif


#endif
