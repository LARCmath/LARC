//                       organize.h
/******************************************************************
 *                                                                *
 * Copyright (C) 2014, Institute for Defense Analyses             *
 * 4850 Mark Center Drive, Alexandria, VA; 703-845-2500           *
 * This material may be reproduced by or for the US Government    *
 * pursuant to the copyright license under the clauses at DFARS   *
 * 252.227-7013 and 252.227-7014.                                 *
 *                                                                *
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

#ifndef ORGANIZE_H
#define ORGANIZE_H

// Standard Libaries
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <pthread.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \ingroup larc
 * \brief Spawns a reporting thread
 * \param period The number of seconds between reports
 */
void create_report_thread(uint64_t period); 

/*!
 * \ingroup larc
 * \brief Manages the reporting thread
 * \param arg The time in seconds between reports
 */
void *report(void *arg);

/*!
 * \ingroup larc
 * \brief Writes various LARC statistics to a designated location
 * \param tabs determines whether some report lines begin with tabs
 * \param outfilepath The location for the report output (can be "stdout")
 */
void rusage_report(int tabs, char *outfilepath);

/*!
 * \ingroup larc
 * \brief Makes an immediate statistics report to a designated location 
 * \param outfilepath The path to the reporing location
 */
void report_now(char *outfilepath);

/*!
 * \ingroup larc
 * \brief Tests available memory to determine whether the LARC process should terminate itself
 * \param arg The number of seconds between memory tests
 */
void *seppuku(void *arg);

/*!
 * \ingroup larc
 * \brief Creates the matrix and ops stores, sets parameters for neighborhoods
 *
 * \param matrix_store_exp The base2 log of the desired size of the matrix store
 * \param op_store_exp The base2 log of the desired size of the operations store
 * \param max_level The base2 log of the largest dimension of any matrix to be created
 * \param sighash The number of bits of precision kept in neighborhood
 * \param zerobitthresh The negative base2 log of the smallest number not in the neighborhood of zero
 * \param verbose Determines which informational messages are printed
 */
void initialize_larc(size_t matrix_store_exp, size_t op_store_exp, mat_level_t max_level, int sighash, int zerobitthresh, int verbose);

/*!
 * \ingroup larc
 * \brief A debugging tool; prints out the precision of double and long double on the current processor
 */
void precision_testing();

/*!
 * \ingroup larc
 * \brief Prints paragraph explaining matrix and operations store (see also explain_locality_sensitive_hashing).
 */
void explain_matrix_and_operation_storage();
  
/*!
 * \ingroup larc
 * \brief Prints paragraph explaining locality sensitive hashing and identification of scalars in same region.
 */
void explain_locality_sensitive_hashing();
  
/*!
 * \ingroup larc
 * \brief Prints explanation of level.
 */
void explain_level();
  
/*!
 * \ingroup larc
 * \brief Prints explanation of scalarType and lists the options.
 */
void explain_scalarType();  

/*!
 * \ingroup larc
 * \brief Prints the allowed values for verbose and their meanings.
 */
void explain_verbosity();

/*!
 * \ingroup larc
 * \brief Looks in /proc/meminfo and returns the availableMemory in kB, if unable to open /proc/meminfo returns 0.
 */
int64_t memory_available_GiB();
  

#ifdef __cplusplus
}
#endif

#endif
