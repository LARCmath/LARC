//                       organize.h
/******************************************************************
 *                                                                *
 * Copyright (C) 2014-2024, Institute for Defense Analyses        *
 * 4850 Mark Center Drive, Alexandria, VA; 703-845-2500           *
 * This material may be reproduced by or for the US Government    *
 * pursuant to the copyright license under the clauses at DFARS   *
 * 252.227-7013 and 252.227-7014.                                 *
 *                                                                *
 * LARC (Linear Algebra via Recursive Compression)                *
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

#ifndef ORGANIZE_H
#define ORGANIZE_H

// Standard Libaries
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
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
 *
 * The thread calls report() every "period" seconds
 * to print current LARC statistics to the screen.
 *
 * \param period The number of seconds between reports
 */
void create_report_thread(uint64_t period); 

#ifndef SWIG
/*!
 * \ingroup larc
 * \brief Stops a thread
 *
 * This routine is used to stop the reporting thread started by
 * create_report_thread(), or to stop the stopHogging thread started
 * by initialize_larc(). It calls pthread_cancel() and pthread_join()
 * to do so.
 *
 * \param thread The handle of the thread to be stopped
 */
void stop_thread(pthread_t *thread);
#endif // #ifndef SWIG

/*!
 * \ingroup larc
 * \brief Writes various LARC statistics to a designated location
 *
 * Memory statistics are reported from the standard location /proc/meminfo,
 * and time statistics from the rusage structure.
 * In addition to creating this report, it tests to see if the program is
 * using too much memory (defined as 95% of memory allocated, along with
 * the program itself using 25% of this memory) and if so exits gracefully.
 *
 * \param tabs controls initial tabs for ome report lines (0: no tabs)
 * \param outfilepath The location for the report output (can be "stdout")
 */
void memory_and_time_report(int tabs, char *outfilepath);

/*!
 * \ingroup larc
 * \brief Makes an immediate statistics report to a designated location 
 *
 * This program may be called directly by the user to get an immediate
 * report of the matrixStore, operationsStore, and time/memory usage
 * statistics. More commonly, it is called periodically by the reporting
 * thread via the report() function.
 *
 * \param outfilepath The path to the reporing location
 */
void report_now_on_all_stores(char *outfilepath);

/*!
 * \ingroup larc
 * \brief Prints the allowed values for verbose and their meanings.
 */
void explain_verbosity();


/*!
 * \ingroup larc
 * \brief Prints out the names of various explanatory files included with the LARC and MyPyLARC packages
 */
void list_explanatory_resources();

/*!
 * \ingroup larc
 * \brief Prints overview of LARC and MyPyLARC packages.
 */
void introduce_LARC_and_MyPyLARC();


/*!
 * \ingroup larc
 * \brief Prints paragraph explaining matrix and operations store (see also explain_collapsingScalars).
 */
void explain_matrix_and_operation_storage();
  
/*!
 * \ingroup larc
 * \brief Prints explanation of level.
 */
void explain_level();

/*!
 * \ingroup larc
 * \brief Prints explanation of why LARC collapses nearly identical scalars together.
 *
 *  LARC only saves the first appearing scalar from a set of nearly identical scalars in order to address finite precision issues and allow pseudo symbolic computation. 
 */
void explain_scalar_techniques();


/*!
 * \ingroup larc
 * \brief Prints overview of ways hash functions are used in LARC and lists explanatory functions.
 *
 *  LARC uses hash functions to implement many of its capabilities: 
 *  unique storage of each matrix and submatrix,
 *  quick retrieval of a matrix record without knowing its MatrixID, 
 *  memoizing  operations so they only need be carried out once, 
 *  using a hash filter for  fast traversal of hash chains, and snapping a scalar to
 *  to previously-stored sufficiently-close scalar.
 */
void explain_hashing_uses_in_LARC();


/*!
 * \ingroup larc
 * \brief Prints explanation of multiplicative Fibonacci hash and LARC implementation.
 *
 *  LARC has a short and elegant implementation of a multiplicative Fibonacci hash for 64-bit integers.  LARC uses this fundamental hash to build all other hash functions.
 */
void explain_multiplicative_Fibonacci_hash();

    
/*!
 * \ingroup larc
 * \brief Prints explanation of how hashes of lists are used to access MatrixStore and OperationsStore.
 *
 *  This function describes how LARC hashes a list.
 *  List hashing is used in the MatrixStore to retrieve a MatrixRecord
 *  without knowing the MatrixID.  A clever choice of inputs to the hash
 *  function allows this, while ensuring that only a single MatrixRecord
 *  is made for each unique matrix.
 *  List hashing is used in the OperationsStore to retrieve the results of
 *  operations that have already been carried out before.
 */
void explain_hashing_lists();


/*!
 * \ingroup larc
 * \brief Prints explanation of various scalarTypes are hashed in LARC.
 *
 * This function describes how LARC hashes
 *  the various scalarTypes in a way that preserves the entropy
 *  of the underlying scalars.
 */
void explain_hashing_scalarTypes();

/*!
 * \ingroup larc
 * \brief Prints explanation of how LARC snaps scalars to sufficiently-close previously-stored scalars.
 *
 *  This function describes how LARC implements
 *  snapping to a previously stored scalar that is 'close enough' to a
 *  target scalar.   In SPR mode, the locality sensitive hash will, with
 *  high probablility, snap to a good previously stored scalar.
 *  In MAR mode, the locality sensitive hash is guaranteed to find a
 *  scalar to snap to (if any that are sufficiently nearby are stored),
 *  however, it might not snap to the one which is closest to the target
 * (if there is more than one stored scalar nearby).
 */
void explain_hashing_for_snapping_scalars();


/*!
 * \ingroup larc
 * \brief Prints description of how hash filters are used to speed searches along hash chains.
 *
 *  This function describes how LARC uses hash filters
 *  LARC hash filters are small integers saved in the hash chains of
 *  the ScalarStore.  When traversing the hash chain searching for a
 *  previously stored scalar in the right region, an initial check is made
 *  to see whether the value of the filter is correct, before a more 
 *  expensive part of the search is carried out.
 */
void explain_hash_filters();


/*!
 * \ingroup larc
 * \brief Prints description of how to see hash chain statistics and what should be expected.
 *
 *  This function  describes how to get statistics on hash 
 *  chain length and what expected longest chain length should be.
 */
void explain_hash_statistics();
    
    
  
/*!
 * \ingroup larc
 * \brief Prints explanation of scalarType and lists the options.
 */
void explain_scalarType();

/*!
 * \ingroup larc
 * \brief Returns information about the current scalarType.
 */
const char* get_string_scalarType();

/*!
 * \ingroup larc
 * \brief Returns information about the LARC version.
 * \return A string with the value of the current LARC software releasae.
 */
const char* get_string_larc_version();

/*!
 * \ingroup larc
 * \brief Prints information about the LARC version.
 */
void print_larc_version();

/*!
 * \ingroup larc
 * \brief Prints information about the LARC settings.
 */
void print_larc_parameters_and_compile_options();

/*!
 * \ingroup larc
 * \brief Looks in /proc/meminfo and returns the availableMemory in GB, if unable to open /proc/meminfo returns 0.
 * \return the number of Gigabytes of memory used by main() and its subroutines
 */
int64_t memory_available_GiB();

/*!
 * \ingroup larc
 * \brief Creates the matrix and ops stores, sets parameters for locality hash regions.
 * 
 * LARC has two types of locality hashes SPR (Single-tile Probabilitistic Retrieval hash) and MAR
 * (Multi-tile Assured Retrieval hash). LARC only saves at most one scalar in each
 * locality region.  The width of the locality hash regions is specified to be
 * 1/(2^(regionbitparam)) for those regions not "near zero" (meaning that for
 * for real scalarTypes the region contains the point zero, and for complex 
 * scalarTypes contains a point on the real or imaginary axes). The width of the
 * regions "near zero" is 1/(2^(zeroregionbitparam)).
 * In MAR, zeroregionbitparam is constrained to be regionbitparam-1, which creates
 * a region that is the optimal width of a MAR super-region.
 * (see MAR explanation in papers for more detailed explanation).
 *
 * \param matrix_store_exp The base2 log of the desired size of the matrix store
 * \param op_store_exp The base2 log of the desired size of the operations store
 * \param max_level The base2 log of the largest dimension of any matrix to be created
 * \param regionbitparam The width of locality hash regions not "near zero" will be 1/(2^(regionbitparam))
 * \param zeroregionbitparam The width of the locality hash regions "near zero" is 1/2^(zeroregionbitparam).
 * \param verbose Determines which informational messages are printed
 */
void initialize_larc(size_t matrix_store_exp, size_t op_store_exp, mat_level_t max_level, int regionbitparam, int zeroregionbitparam, int verbose);

/*!
 * \brief This routine deallocates all memory used by larc and also terminates all threads started by initialize_larc.
 *
 * The intent for this routine is to allow a user to get a fresh copy of
 * LARC. One use case would be to have both initialize_larc and shutdown_larc
 * in a loop for testing RSH parameters for a particular problem.
 */
void shutdown_larc(void);

/*!
 * \ingroup larc
 * \brief A debugging tool; prints out the precision of double and long double on the current processor
 */
void precision_testing();


#ifdef __cplusplus
}
#endif

#endif
