//                            info_store.h
//                  Matrix Information Store for LARC
/******************************************************************
 *                                                                *
 * Copyright (C) 2014-2024, Institute for Defense Analyses        *
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


#ifndef LARC_INFO_STORE_H
#define LARC_INFO_STORE_H

#include <stdint.h>
#include "larc.h"

/*!
 * \ingroup larc
 * \brief Generates a hash value based on an input packedID and an info_type_t value
 * \param pID A packedID
 * \param info_type An enum value indicating the type of metadata stored
 * \return The hash value
 */
uint64_t
hash_from_info(int64_t pID, info_type_t info_type);

/*!
 * \ingroup larc
 * \brief Print out the names understood by the info store
 *
 * lists out all the info_names and description of them 
 * to help users interface with info chain functions 
 */
void list_info_names();

/*!
 * \ingroup larc
 * \brief Outputs the name associated to an enum info_type value
 * \param info_type An enum value 
 * \return A string containing the name of the enum value
 *
 * This function returns the name associated with enum info_type value
 */

char *return_info_name(uint64_t info_type);

// return an array of the info_names to the user
// char** info_names_array();

/*!
 * \ingroup larc
 * \brief Prints out information about the info store
 * \param outfilepath The location to which the function prints
 */
void info_store_report(char *outfilepath);

#ifndef SWIG
/*!
 * \ingroup larc
 * \brief Creates the LARC infoStore
 * \param exponent The log base 2 of the size of the infoStore hash table
 * \return 1 if infoStore created successfully
 */
int create_info_store(size_t exponent);

/*!
 * \ingroup larc
 * \brief Frees memory allocated to LARC's infoStore
 *
 * This routine is called by shutdown_larc. It should not be called by
 * any other routine.
 */
void free_info_store(void);
#endif // #ifndef SWIG

/*!
 * \ingroup larc
 * \brief Creates a new larc_info_t record and fills it with the variables passed to this function
 * \param info_type The type of metadata to be stored
 * \param my_pID The packedID to be stored
 * \param info_data Metadata about the matrix
 * \return 0 on success, or error code if failed 
 */
int info_set(info_type_t info_type, int64_t my_pID, const char *info_data);

/*!
 * \ingroup larc
 * \brief Finds any stored metadata of a particular type about a matrix
 * \param info_type The type of metadata
 * \param my_pID The matrix of interest
 * \return A string containing the metadata
 */
char *info_get(info_type_t info_type, int64_t my_pID);

/*!
 * \ingroup larc
 * \brief Given a string that corresponds to the name of an enum value for the infoStore, returns the enum value
 * \param info_name The input string
 * \return The enum value if it exists, otherwise returns INVALID_INFO
 */
info_type_t get_info_type_from_string_name(const char *info_name);

/*!
 * \ingroup larc
 * \brief Returns the hash value associated to a packedID and the type of metadata to be stored
 * \param my_pID The packedID
 * \param info_name The type of metadata to be stored
 * \return The hashID value associated with the info record 
 */
int64_t get_info_hash(int64_t my_pID, const char *info_name);

/*!
 * \ingroup larc
 * \brief Prints information about the infoStore hash chain corresponding to a specified hash value
 * \param hash The hash value for an infoStore hash chain
 * \param outfilepath The file to which the information is printed
 * \param comment A string that the user can have printed to that file
 * \return 1 on success
 */
int fprint_info_hash_chain_info(uint64_t hash, char *outfilepath, char *comment);

/*!
 * \ingroup larc
 * \brief Same as fprint_info_hash_chain_info except that the print is to stdout
 * \param hash The hash value for an infoStore hash chain
 * \param comment A string that the user can have printed to stdout
 * \return 1 on success
 */
int print_info_hash_chain_info(uint64_t hash, char *comment);

/*!
 * \ingroup larc
 * \brief Cleans a hash chain in the info store
 * \param hashID The hash specifing the hash chain to be cleaned
 * \return 1 on success
 */
int clean_info_store_hash_chain(uint64_t hashID);

#endif
