//                        organize.c 
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

/*!
 * \file organize.c
 * \brief This file contains LARC initialization routines, functions which
 * print explanatory text, and routines which print LARC statistics.
 *
 * These functions print explanatory text:
 *    introduce_LARC_and_MyPyLARC() 
 *       a quick overview of the ideas for LARC and MyPyLARC
 *   explain_verbosity()
 *      what verbosity levels are available and what they do
 *   explain_matrix_and_operation_storage()
 *      How LARC stores matrix records and operations records using hashes.
 *   explain_level()
 *      matrix size is 2^row_level * 2^col_level,max_level
 *   explain_scalar_techniques()
 *      LARC saves only the first appearing of a set of nearly identical scalars 
 *      it uses regional representation and locality hashing to implement this.
 *      LARC uses this technique to  mimic symbolic computation
 *      and reduce numerical precision issues.
 *      Compares SPR and MAR methods.
 *   explain_scalarType()
 *      User-selected scalarTypes are selected at compile time and include:
 *      Real, Integer, Complex, and Multiprecision options
 *   explain_multiplicative_Fibonacci_hash()
 *       This function explains the fundamental hash mult_golden_hash() 
 *       from which all other hash functions in LARC are built.
 *       mult_golden_hash() is a multiplicative Fibonacci hash which acts 
 *       on 64-bit integers.
 *
 * These functions print information that helps with reproduciblity of experiments:
 *   print_larc_version()
 *     Larc version (major/minor/patch), GIT date, GCC version, and scalarType,
 *     plus info on multiprecision packages where relevant
 *   create_output_strings_to_print_settings(,,,,)
 *     Matrix/Operations store hash exponents/Region hash parameters/max_level
 *   print_larc_parameters_and_compile_options()
 *     Get above 5 parameters and pass to "print_larc_parameters_and_compile_options()"
 *
 *  This function initializes the storage using the 5 parameters plus verbosity
 *      initialize_larc(,,,,,)
 *       Create+preload storage/ Check parameter+abort conditions/ verbosity print
 *
 * These routines provide for threaded reporting of running LARC statistics
 *   create_report_thread(uint64_t period)
 *     Calls pthread_create to create a report_thread with a timeout parameter
 *   stop_thread(pthread_t *thread)
 *     Calls pthread_cancel and pthread_join for thread control
 *   memory_and_time_report(int tabs, char *outfilepath)
 *     Report stats. Mem from /proc/meminfo.Time from rusage struct.(Self-kill)
 *   report(void *arg)
 *     Handles sleep time and periodic report.  Calls "report_now_on_all_stores" below.
 *   report_now_on_all_stores(char *outfilepath)
 *     Calls memory_and_time_report and report routines from matrix_store.c, op_store.c
 *
 *Additional routines for Safety checks
 *  stopHogging(void *arg)
 *     Background checking for Memory/Timing Self-kill conditions
 *  precision_testing()
 *     Simply prints out the mantissa information for Doubles and Long Doubles
 *  memory_available_GiB()
 *     Gets memory statistics from standard location /proc/meminfo in GigBytes
 *  shutdown_larc(void)
 *     Clears/frees all LARC memory if initializing LARC without exiting.
 */


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
#include <inttypes.h> //for printing 
#include <curses.h>

// Our header files structures and functions
#include "larc.h"
#include "op_store.h"
#include "global.h"
#include "organize.h"
#include "version.h"
#include "scalars.h"
#include "show_scalars.h"
#include "info_store.h"
#include "matmath.h"
#include "matrix_store.h"

// Global verbosity parameter.  See larc.h for explanation.
verbose_type_t VERBOSE;

/*
 * "local global" variables needed for the two thread processes spawned by
 * LARC. The stopHogging process kills the program that initialized LARC if
 * it consumes too large a percentage of total machine resources. The 
 * report_thread process periodically prints out statistics on the resources
 * that LARC is using.
 */
pthread_t stopHogging_thread, report_thread;
static char report_thread_exists = 0;

/*!
 * \ingroup larc
 * \brief Tests available memory to determine whether the LARC process should terminate itself
 *
 * This function is called periodically by the stopHogging thread to determine
 * whether the current program should be killed due to excessive consumption
 * of system resources (at least 95% of total memory in use, of which at least
 * 25% is in use by this program). A full usage report is printed before exit.
 *
 * \param arg The number of seconds between memory tests
 */
static void *stopHogging(void *arg)
{
  unsigned int timeout = *(unsigned int*)arg;
  struct rusage rusage;
  int tabs = 0;
  
  while (1) {
    sleep(timeout);
    printf("Running stopHogging check\n");

    // get memory statistics from standard location /proc/meminfo
    FILE *fp = fopen("/proc/meminfo", "r");
    uint64_t mem_total=0, mem_free=0, cached=0;
    if (fp)
      {
	char buf[100];
        // read first line of file /proc/meminfo
	if (fgets(buf, 100, fp))
	  {
	    sscanf(buf, "MemTotal: %" SCNu64 " kB", &mem_total);
	  }
        // read second line of file /proc/meminfo
	if (fgets(buf, 100, fp))
	  {
	    sscanf(buf, "MemFree: %" SCNu64 " kB", &mem_free);
	  }
	(void)! fgets(buf, 100, fp);
        // read fourth line of file /proc/meminfo
	if (fgets(buf, 100, fp))
	  {
	    sscanf(buf, "Cached: %" SCNu64 " kB", &cached);
	  }
	fclose(fp);
      }
    else if (VERBOSE>SILENT)
      {
        printf("WARNING in stopHogging check: could not open /proc/meminfo\n");
        printf("No checking will be performed!\n");
      }

    // determine whether necessary to kill process
    if (
	(mem_free || cached) && 
	mem_total &&
	(((1.0*mem_free)+(1.0*cached))/mem_total < 0.05) && 
	(1.0*rusage.ru_maxrss/mem_total > 0.25)
	)
      {
	if (!getrusage(RUSAGE_SELF, &rusage))
	  {
	    double cputime = (rusage.ru_utime.tv_sec + rusage.ru_stime.tv_sec
			      +(rusage.ru_utime.tv_usec/1000000.0)
			      + (rusage.ru_stime.tv_usec/1000000.0));
            int cpumins = ((int) cputime + 30) % 60;
            int cpuhours = ((int) cputime + 1800) % 3600;
	    printf("%sCPU Time: %3.2f sec (roughly %d mins) (roughly %d hours)\n", 
                   tabs?"\t":"", cputime, cpumins,cpuhours);
	    printf("%sMax RSS size: %gMB\n", tabs?"\t":"", rusage.ru_maxrss/1000.0);
	  }
        report_now_on_all_stores("stdout");
        printf("%s Mem Free=%" PRIu64 ", Cached=%" PRIu64 ", Mem Total=%" PRIu64 "\n", tabs?"\t":"", mem_free,cached, mem_total);
        printf("\n%sTotal Free/Memory: %g/%gMB\n", tabs?"\t":"", (mem_free+cached)/1000.0, mem_total/1000.0);
        printf(">>>Too little \"real\" memory left, exiting<<<\n");
        exit(1);
      }
  }
}

/*!
 * \ingroup larc
 * \brief Manages the reporting thread
 *
 * This function is called periodically by the LARC reporting thread. It
 * calls report_now_on_all_stores() to output statistics on the matrixStore, operationsStore
 * and on total time/memory usage to stdout. report_now_on_all_stores() also may kill the
 * process that called initialize_larc() in the same way and under the same
 * conditions as the stopHogging thread.
 *
 * \param arg The time in seconds between reports
 */
static void *report(void *arg)
{
    uint64_t timeout = *(uint64_t *)arg;
    printf("\n\n=== Starting Periodic Reporting with Interval of %" PRIu64 " Seconds ===\n",timeout);

    while (1) {
	sleep(timeout);
	printf("\n\n==================== Periodic Report (%" PRIu64 " sec) ====================\n", timeout);
        printf("NOTE: the main thread continues to run while the periodic report is\n");
        printf("generated, so number of stored matrices may increase!");
	printf("\n========================================================================\n\n");
	report_now_on_all_stores("stdout");
	// matrix_store_report();
	// op_store_report();
	// memory_and_time_report(0,"stdout");
	printf("\n========================================================================\n\n");
    }
    exit(1);
}

// These routines provide for threaded reporting of LARC statistics during
// operation.

void create_report_thread(uint64_t period) 
{
  // this routine spawns a thread which provides periodic reporting of LARC
  // matrixStore, operationStore, and time/memory statistics.
  static uint64_t timeout;
  timeout = period;
  int s;
  s = pthread_create(&report_thread, NULL, report, &timeout);
  if (s != 0) { errno = s; perror("pthread_create"); exit(EXIT_FAILURE); }
  report_thread_exists = 1;
}

void stop_thread(pthread_t *thread)
{
  // This routine safely terminates the pthread given as its argument.
  int s;
  void *res;
  s = pthread_cancel(*thread);
  if (s != 0) { errno = s; perror("pthread_create"); exit(EXIT_FAILURE); }
  // the join function does not return until the thread has terminated
  s = pthread_join(*thread, &res);
  if (s != 0) { errno = s; perror("pthread_join"); exit(EXIT_FAILURE); }
  if (res != PTHREAD_CANCELED)
  {
    fprintf(stderr,"ERROR in stop_thread: result of pthread_join() was\n");
    fprintf(stderr,"not PTHREAD_CANCELED!\n");
    exit(EXIT_FAILURE);
  }
}

// This routine provides for threaded reporting of time and memory statistics.
// It also causes the program to exit if it detects excessive memory is being
// consumed.
void memory_and_time_report(int tabs, char *outfilepath)
{
  FILE *f; 

  if (strcmp(outfilepath,"stdout")) {
     printf("Printing rusage report to file %s\n", outfilepath);
     f = fopen(outfilepath, "a"); 

  }
  else {
    printf("Printing rusage report to screen\n");
    f = stdout;
  }
	
  // STILL CONTAINS stopHogging function
  struct rusage rusage;
  FILE *fp = fopen("/proc/meminfo", "r");
  uint64_t mem_total=0, mem_free=0, cached=0;
  
  if (!getrusage(RUSAGE_SELF, &rusage))
    {
      double cputime = (rusage.ru_utime.tv_sec + rusage.ru_stime.tv_sec
			+(rusage.ru_utime.tv_usec/1000000.0)
			+ (rusage.ru_stime.tv_usec/1000000.0));
      fprintf(f,"%sCPU Time: %3.2f sec\n", tabs?"\t":"", cputime);
      fprintf(f,"%sMax RSS size: %gMB\n", tabs?"\t":"", rusage.ru_maxrss/1000.0);
#ifdef HASHSTATS
      if (VERBOSE>SILENT)
        fprintf(f,"\nWARNING: HASHSTATS on, using more time and memory!\n");
      if (VERBOSE>BASIC) {
        fprintf(f," When HASHSTATS is defined LARC will generate statistics\n");
        fprintf(f," on the numbers of visits to each hash bucket and\n");
        fprintf(f," the lengths of the hash chain from each bucket.\n");
	fprintf(f," To change: comment out define HASHSTATS in src/larc.h.\n");
        }
#else
      if (VERBOSE>BASIC) {
        fprintf(f,"\nHASHSTATS is not defined in larc.h, uncomment the\n");
        fprintf(f,"  define statement, to generate statistics on\n");
        fprintf(f,"  the numbers of visits to each hash bucket and the\n");
        fprintf(f,"  lengths of the hash chain from each bucket.\n");
        }
#endif      
    }
  if (fp)
    {
      char buf[100];
      
      // read first line of file /proc/meminfo
      if (fgets(buf, 100, fp))
	{
	  sscanf(buf, "MemTotal: %" SCNu64 " kB", &mem_total);
	}
      // read second line of file /proc/meminfo
      if (fgets(buf, 100, fp))
	{
	  sscanf(buf, "MemFree: %" SCNu64 " kB", &mem_free);
	}
      (void)! fgets(buf, 100, fp);
      // read fourth line of file /proc/meminfo
      if (fgets(buf, 100, fp))
	{
	  sscanf(buf, "Cached: %" SCNu64 " kB", &cached);
	}
      fclose(fp);
    }
  fprintf(f,"%sMem Free=%" PRIu64 ", Cached=%" PRIu64 ", Mem Total=%" PRIu64 "\n", tabs?"\t":"", mem_free,cached, mem_total);
  fprintf(f,"%sTotal Free/Memory: %g/%gMB\n", tabs?"\t":"", (mem_free+cached)/1000.0, mem_total/1000.0);
  fflush(f);
  
  // static int kill_me_too_big = 1;
  //  if (kill_me_too_big && mem_free && mem_total &&

  if (
      (mem_free || cached) && 
       mem_total &&
      (((1.0*mem_free)+(1.0*cached))/mem_total < 0.05) && 
       (1.0*rusage.ru_maxrss/mem_total > 0.25)
     )
    {
      fprintf(f,">>>Too little \"real\" memory left, exiting<<<\n");
      fflush(f);
      if  (strcmp(outfilepath,"stdout")) fclose(f);   // if not stdout then close
      exit(1);
    }
}

// This routine is called by the report() function, in turn called
// periodically by the reporting thread. It may also be called directly
// by the user.
void report_now_on_all_stores(char *outfilepath)
{
		matrix_store_report(outfilepath);
		op_store_report(outfilepath);
		memory_and_time_report(0,outfilepath);
}

// The following functions print information to stdout. Some are merely
// explanatory information, while others provide information about the
// current LARC version or the parameters chosen at compile time or during
// initialization.

void explain_verbosity() {
    // This function prints an explanation of the different levels of verbosity
    fprintf(stdout,"\nThe different levels of verbosity are:\n");
    fprintf(stdout,"\tSILENT: verbose = 0   only prints errors;\n");
    fprintf(stdout,"\tBASIC:  verbose = 1   prints minimal information\n");
    fprintf(stdout,"\t                      including errors and warnings;\n");
    fprintf(stdout,"\tCHATTY: verbose = 2   adds informative comments;\n");
    fprintf(stdout,"\tDEBUG:  verbose = 3   adds debugging comments;\n");
    fprintf(stdout,"\tALL:    verbose = 4   simply won't shut up.\n");
}

void introduce_LARC_and_MyPyLARC() {
    fprintf(stdout,"OVERVIEW:                    \n");
    fprintf(stdout,"LARC (Linear Algebra via Recursive Compression)                    \n");
    fprintf(stdout,"is a software package that stores matrices in a recursively compressed \n"); 
    fprintf(stdout,"format and performs operations on them while they are compressed.   \n");
    fprintf(stdout,"LARC assigns each scalar, vector, and matrix a unique integer MatrixID. \n");
    fprintf(stdout,"A 2^m by 2^n matrix is then recursively defined by a list of the four \n");
    fprintf(stdout,"MatrixIDs of its quadrant submatrices.						    \n");
    fprintf(stdout,"Functions depending on these MatrixIDs are used to produce hash	   \n");
    fprintf(stdout,"values accessing three hash tables which store scalars, matrices and   \n");
    fprintf(stdout,"previously computed operations.							    \n");
    fprintf(stdout,"These tables enable LARC to quickly retrieve previously stored matrices  \n");
    fprintf(stdout,"and operations, and thus to never store the same matrix twice or repeat \n");
    fprintf(stdout,"an operation.						    \n");
    fprintf(stdout,"											    \n");
    fprintf(stdout,"LARC contains a method called Regional Representative Retrieval            \n");
    fprintf(stdout,"that uses a locality sensitive hash and other techniques to snap	       \n");
    fprintf(stdout,"nearby scalars together.  Scalars that are important to identities	       \n");
    fprintf(stdout,"can be preloaded with precision, enabling LARC to mimic the        \n");
    fprintf(stdout,"behavior of symbolic computations over a ring.		       \n");
    fprintf(stdout,"									       \n");
    fprintf(stdout,"The LARC matrix math package and a companion Python package       \n");
    fprintf(stdout,"MyPyLARC containing tutorials and sample applications are available	  \n");
    fprintf(stdout,"on GitHub/LARCmath. 					       \n");
    fprintf(stdout,"An explanatory companion paper to the associated GitHub     \n");
    fprintf(stdout,"code repository supplying motivations and explanations for LARC	       \n");
    fprintf(stdout,"methodology as well as examples of applications for LARC and      \n");
    fprintf(stdout,"MyPyLARC is in MyPyLARC/doc/LARCandMyPyLARC_ver2.1_2021.pdf.\n");
}

void explain_matrix_and_operation_storage() {
    // This function prints an explanation of how LARC stores matrix
    // records and operations records .
    // Also see explain_scalar_techniques() 
    fprintf(stdout,"\nStoring Matrix Records and Operation Records in LARC:\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"The LARC (linear algebra via recursive compression) package\n");
    fprintf(stdout,"represents 2^r by 2^c matrices in a recursively compressed format\n");
    fprintf(stdout,"based on a decomposition of the matrix into its quadrant submatrices.\n");
    fprintf(stdout,"LARC assigns an integer MatrixID to each unique matrix and submatrix.\n");
    fprintf(stdout,"Matrices are stored recursively by recording their MatrixID and a list\n");
    fprintf(stdout,"of the MatrixIDs of their four half-dimension quadrant submatrices.\n");
    fprintf(stdout,"For example, if a matrix M had the four quadrant submatrices A, B, C, D,\n");
    fprintf(stdout,"given in the order upper-left, upper-right, lower-left, lower-right,\n");
    fprintf(stdout,"then the MatrixRecord for M would contain the MatrixIDs:\n");
    fprintf(stdout,"                         M_ID;       [A_ID, B_ID, C_ID, D_ID].\n");
    fprintf(stdout,"LARC was designed to work on applications with some quadrant submatrix\n");
    fprintf(stdout,"reuse.  Thus only the first copy of each matrix is stored.  In the example, if\n");
    fprintf(stdout,"matrices A and D were equal, then the MatrixRecord for M would contain:\n"); 
    fprintf(stdout,"                         M_ID;       [A_ID, B_ID, C_ID, A_ID].\n");
    fprintf(stdout,"The list of the four MatrixIDs of the quadrant submatrices is called the\n");
    fprintf(stdout,"SubMatList.  The MatrixRecord for matrix M is stored in the MatrixStore\n"); 
    fprintf(stdout,"hash table at the index location given by a hash of SubMatList(M).\n");
    fprintf(stdout,"Vectors of length 2^m are stored in a similar way by subdivision in half.\n");
    fprintf(stdout,"When the dimension of a matrix is only 1 by 1 (a single scalar value) it no\n");
    fprintf(stdout,"longer has a recursive definition.  The scalar is given a ScalarRecord in the\n");
    fprintf(stdout,"ScalarStore hash table containing its MatrixID and the value of the scalar.\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"OPERATIONS:\n");
    fprintf(stdout,"Matrix operations are performed recusively and take MatrixIDs as their\n");
    fprintf(stdout,"input arguments, and return a MatrixID. LARC performs operations on\n");
    fprintf(stdout,"matrices recursively without leaving the compressed format. During a recursive\n");
    fprintf(stdout,"operation that returns a matrix, the SubMatList of the result of the operation\n");
    fprintf(stdout,"is produced.  A check is made to see whether this SubMatList corresponds\n");
    fprintf(stdout,"to a matrix that has already been stored in the MatrixStore. This check is\n");
    fprintf(stdout,"fast, because hashing the SubMatList gives the location in the MatrixStore\n");
    fprintf(stdout,"where such a MatrixRecord would have been stored.   If the matrix was\n");
    fprintf(stdout,"already stored its MatrixID is returned, otherwise a new MatrixRecord is\n");
    fprintf(stdout,"created and its MatrixID is returned.\n");
    fprintf(stdout,"Matrix operations are memoized in the OperationsStore hash table and \n");
    fprintf(stdout,"indexed by a hash of the operation name and the list of MatrixIDs of the\n");
    fprintf(stdout,"inputs to the operation.  In this way operations never need to be repeated.\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"The storage and retrieval of scalars from the ScalarStore in LARC uses\n");
    fprintf(stdout,"some special techniques which allow LARC to snap nearby scalars together\n");
    fprintf(stdout,"and mimic symbolic computation for selected identities.   An explanation\n");
    fprintf(stdout,"of these techniques is given in explain_scalar_techniques().\n");
    //    fprintf(stdout,"A list of supported scalarTypes is given in explain_scalarTypes().\n");
    //    fprintf(stdout,"Parameters for LARC initialization are described in explain_initialization().\n");
    // fprintf(stdout,"Supported formats for input and output of matrices, and methods of matrix\n");
    //  fprintf(stdout,"creation are described in explain_matrix_input_output_and_creation.\n");
}

void explain_level() {
    // this function just prints a message to explain what a 'level' means
    // in LARC
    fprintf(stdout,"\nMatrix Dimensions:\n");
    fprintf(stdout,"Because LARC defines matrices recursively in terms of quadrant submatrices,\n");
    fprintf(stdout,"each dimension of a matrix in LARC must be a power of two. We\n");
    fprintf(stdout,"use 'level' to indicate the exponent of two describing the dimension\n");
    fprintf(stdout,"of the matrix; the size of a matrix is 2^(row_level) by 2^(col_level).\n");
    fprintf(stdout,"When LARC reads in a matrix from an external format such as sparsely\n");
    fprintf(stdout,"formated Matrix Market Exchange Format, it will pad to the nearest\n");
    fprintf(stdout,"power of 2 dimension.\n");
    fprintf(stdout,"In the current version of LARC, the max_level of matrices is selected at\n");
    fprintf(stdout,"initialization, restricting max size to 2^(max_level) by 2^(max_level).\n");
}

void explain_scalar_techniques() {
   //  LARC saves only the first appearing of a set of nearly identical scalars
   //
   // NOTE:  This explanation is very similar to explain_hashing_for_snapping_scalars()
   //             This other one has more on hash function specifics, while this one has
   //             a section on mimicking symbolic computation.
   // 
   fprintf(stdout,"\n");
   fprintf(stdout,"Scalar Techniques: Finite Precision Issues, Regional Representative\n");
   fprintf(stdout,"Retrieval (snapping scalars together) and Pseudo Symbolic Computation:\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"The earliest version of LARC, when used on large applications, required\n");
   fprintf(stdout,"large amounts of memory to save sets of scalars that were\n");
   fprintf(stdout,"almost identical.  If the math operations had been carried out\n");
   fprintf(stdout,"symbolically the scalars would have been identical, but the finite\n");
   fprintf(stdout,"precision of floating point operations created tiny differences.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"The LARC team has developed a special retrieval method for scalars in LARC\n");
   fprintf(stdout,"called *regional representative retrieval* that addresses these finite\n");
   fprintf(stdout,"precision issues and allows scalars differing only by a tiny amount to\n");
   fprintf(stdout,"be collapsed to a single value. A newly calculated scalar will be replaced\n");
   fprintf(stdout,"with (snapped to) an almost identical previously stored value.\n");
   fprintf(stdout,"LARC saves the first scalar in a region, then uses a fast hash-table based\n");
   fprintf(stdout,"search method (described below) to allow replacement of a newly\n");
   fprintf(stdout,"calculated scalar with the previously-stored (almost identical) scalar.\n");
   fprintf(stdout,"The distance for snapping is tuned by an initialization parameter.\n");
   fprintf(stdout,"This technique allows LARC to handle finite precision issues and also to\n");
   fprintf(stdout,"mimic symbolic computation.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"In order to perform a fast hash-table-based search for nearly identical scalars,\n");
   fprintf(stdout,"LARC divides the space of scalars into small tiles.  For complex scalars,\n");
   fprintf(stdout,"the tiles are two dimensional and for real scalars they are one dimensional.\n");
   fprintf(stdout,"The LARC scalar-storage hash table uses a *locality sensitive hash* that\n");
   fprintf(stdout,"combines two functions, a function T which calculates a tile label for\n");
   fprintf(stdout,"a scalar S, and a function H which hashes tile labels. Then H(T(S)) is a\n");
   fprintf(stdout,"locality hash that hashes all scalars within a single tile to the same value.\n"); 
   fprintf(stdout,"LARC uses this locality hash to identify cases where a newly calculated\n");
   fprintf(stdout,"scalar need not be stored, because it can be replaced by a previously\n");
   fprintf(stdout,"stored scalar.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"LARC has two different schemes for storing and snapping scalars:\n");
   fprintf(stdout,"   - SPR (short for Single-tile Probabilistic Retrieval), and\n");
   fprintf(stdout,"   - MAR (short for Multi-tile Assured Retrieval).\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"SPR replaces a newly calculated scalar with a previously stored scalar\n");
   fprintf(stdout,"in the same tile if one exists, or otherwise stores the new scalar.\n");
   fprintf(stdout,"SPR succeeds with high probability in identifying a newly calculated\n");
   fprintf(stdout,"scalar with an almost-identical previously stored scalar, only failing\n");
   fprintf(stdout,"to collapse two close scalars that are on opposite sides of a tile\n");
   fprintf(stdout,"boundary.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"MAR uses a more complicated scheme in which adjacent tiles can be\n");
   fprintf(stdout,"grouped together to form a region.  This region will contain a single\n"); 
   fprintf(stdout,"stored scalar which lies someplace in the central portion of the region.\n");
   fprintf(stdout,"This trick is carried out by pasting together more than one tile.\n");
   fprintf(stdout,"MAR has a guarantee of snapping a newly calculated scalar to some\n");
   fprintf(stdout,"sufficiently-nearby previously-stored scalar (if it exists), but it may\n");
   fprintf(stdout,"not be the closest previously-stored scalar.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"For either SPR or MAR we call the unique stored scalar in a region\n");
   fprintf(stdout,"the *representative* of the region (for SPR a region consists of a\n");
   fprintf(stdout,"single tile).\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"More details on these two schemes can be found in the explanatory\n");
   fprintf(stdout,"paper in MyPyLARC/doc/LARCandMyPyLARC_ver2.1_2021.pdf.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"Mimicking of Symbolic Computation:\n");
   fprintf(stdout,"LARC uses its locality-sensitive hash to mimic symbolic computation\n");
   fprintf(stdout,"for selected mathematical identities. This is achieved by preloading\n");
   fprintf(stdout,"(in full precision) a selection of scalars involved in these identities.\n");
   fprintf(stdout,"By preloading we mean to place scalars in the store immediately\n");
   fprintf(stdout,"after initialization, so that they become the representatives of their\n");
   fprintf(stdout,"regions.   For example, in our implementation of a discrete Fourier\n");
   fprintf(stdout,"transform we preload the nth-roots of unity 1, w, w^2, ....w^(n-1)\n");
   fprintf(stdout,"where w = e^(2 pi i / n).\n");
   fprintf(stdout,"This assures that when LARC calls:\n");
   fprintf(stdout,"       product((matrixID(w^j ),matrixID(w^k )) \n");
   fprintf(stdout,"the function will return matrixID(w^m ) where m = j+k (mod n).\n");
   fprintf(stdout,"The scalars 0 and 1 are preloaded into LARC during initialization\n");
   fprintf(stdout,"because of their importance as identity elements in mathematical\n");
   fprintf(stdout,"operations.   Similarly, we also preload all identity and zero\n");
   fprintf(stdout,"matrices and set flags in their matrix records so that we can quickly\n");
   fprintf(stdout,"identify and short cut any matrix operation identities.\n");
}

// HASH EXPLANATION #1
void explain_hashing_uses_in_LARC() {
   //  This function explains the different kinds of hashes that are needed
   //  to give LARC its capabilities: unique storage of each matrix and submatrix,
   //  quick retrieval of a matrix record without knowing its MatrixID, memoizing
   //  operations so they only need be carried out once, using a hash filter for
   //  fast traversal of hash chains, and snapping new scalars to
   //  to previously stored scalars that are "close-enough". 
   fprintf(stdout,"\n");
   fprintf(stdout,"                       How LARC uses hashing\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"Hash functions are central to LARC's capabilities.  They are used for:\n");
   fprintf(stdout,"ensuring unique storage of each matrix and submatrix,\n");
   fprintf(stdout,"quick retrieval of a MatrixRecord without knowing its MatrixID,\n");
   fprintf(stdout,"memoizing operations so they only need be carried out once, \n");
   fprintf(stdout,"fast search traversal of hash chains by using a hash filter, and\n");
   fprintf(stdout,"snapping new scalars to sufficiently close previously stored scalars.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"LARC is designed to work on large problems and can create huge\n");
   fprintf(stdout,"hash tables for storing scalars, matrices, and operations; so it is\n");
   fprintf(stdout,"very important that LARC has very well-behaved hash functions.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"* explain_multiplicative_Fibonacci_hash - describes the fundamental\n");
   fprintf(stdout,"     hash that is used to build all other hash functions in LARC.\n");
   fprintf(stdout,"     This is a multiplicative Fibonacci hash has optimal distribution\n");
   fprintf(stdout,"     and is elegantly implemented as an entirely integer calculation.\n");
   fprintf(stdout,"     The hash is often used to produce an E-bit value used as an index\n");
   fprintf(stdout,"     into a hash table of size 2^E.    LARC also uses this hash function\n");
   fprintf(stdout,"     to output a 64-bit integer in order to build other hash functions.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"* explain_hashing_lists - describes how LARC hashes a list.  \n");
   fprintf(stdout,"     List hashing is used in the MatrixStore to retrieve a MatrixRecord\n");
   fprintf(stdout,"     without knowing the MatrixID.  A clever choice of inputs to the hash\n");
   fprintf(stdout,"     function allows this, while ensuring that only a single MatrixRecord\n");
   fprintf(stdout,"     is made for each unique matrix.\n");
   fprintf(stdout,"     List hashing is used in the OperationsStore to retrieve the results of\n");
   fprintf(stdout,"     operations that have already been carried out before.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"* explain_hashing_scalarTypes - describes how LARC hashes\n");
   fprintf(stdout,"     the various scalarTypes in a way that preserves the entropy\n");
   fprintf(stdout,"     of the underlying scalars.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"* explain_hashing_for_snapping_scalars - describes how LARC implements\n");
   fprintf(stdout,"     snapping to a previously stored scalar that is 'close enough' to a\n");
   fprintf(stdout,"     target scalar.   In SPR mode, the locality sensitive hash will, with\n");
   fprintf(stdout,"     high probablility, snap to a good previously stored scalar. \n");
   fprintf(stdout,"     In MAR mode, the locality sensitive hash is guaranteed to find a \n");
   fprintf(stdout,"     scalar to snap to (if any that are sufficiently nearby are stored), \n");
   fprintf(stdout,"     however, it might not snap to the one which is closest to the target\n");
   fprintf(stdout,"     (if there is more than one stored scalar nearby).\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"* explain_hash_filters - describes how LARC uses hash filters. \n");
   fprintf(stdout,"     LARC hash filters are small integers saved in the hash chains of\n");
   fprintf(stdout,"     the ScalarStore.  When traversing the hash chain searching for a\n");
   fprintf(stdout,"     previously stored scalar in the right region, an initial check is made\n");
   fprintf(stdout,"     to see whether the value of the filter is correct, before a more \n");
   fprintf(stdout,"     expensive part of the search is carried out.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"* explain_hash_statistics - describes how to get statistics on hash \n");
   fprintf(stdout,"     chain length and what expected longest chain length should be.\n");
   fprintf(stdout,"\n");
}

// HASH EXPLANATION #2
void explain_multiplicative_Fibonacci_hash() {
   //  This function explains the fundamental hash function that all other hash
   //  functions in LARC are built from.
   //  TODO: When writing this in tex get the accent on Sos and her 2 references from 
   //      Knuth as well as the 2 or 3 references on botanists and the golden mean.
   //      Also make some sort of image to show the distribution of the hash like in Knuth.
   fprintf(stdout,"\n");
   fprintf(stdout,"       LARC's multiplicative Fibonacci hash - its fundamental hash\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"LARC uses hash functions for a variety of reasons including the\n");
   fprintf(stdout,"quick retrieval of previously stored information in potentially huge\n");
   fprintf(stdout,"hash tables such as the ScalarStore and the MatrixStore.  So it is\n");
   fprintf(stdout,"very important that LARC has a very well-behaved hash function.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"LARC has a simple optimally-distributed hash function from which all\n");
   fprintf(stdout,"its other hash functions are built.  This hash function is a multiplicative\n");
   fprintf(stdout,"Fibonacci hash, mult_golden_hash(), which acts on 64-bit integers\n");
   fprintf(stdout,"and is implemented in LARC using entirely integer mathematics.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"A multiplicative hash H(key) uses a constant M (called the multiplier) \n");
   fprintf(stdout,"which is between 0 and 1 and outputs an integer hash value \n");
   fprintf(stdout,"between 0 and S-1 (in our case  S = hash table size) as follows:\n");
   fprintf(stdout,"        H(key) = integer part of [S * (fractional part of (M * key))]. \n");
   fprintf(stdout,"See more in Corman, Leiserson, Rivest and Steine's \n");
   fprintf(stdout,"Introduction to Algorithms (3rd Edition, pp. 263-264).\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"Fibonacci hashes are particularly well-behaved multiplicative hashes\n");
   fprintf(stdout,"whose multiplicative constant is close to the golden ratio\n");
   fprintf(stdout,"        M = (sqrt(5) - 1)/2 = .6180339886... .  \n");
   fprintf(stdout,"See Knuth's Sorting and Searching (2nd Edition, pp.516-518). \n");
   fprintf(stdout,"As Knuth mentions, it was shown by Vera Sos that the golden ratio\n");
   fprintf(stdout,"is one of a general class of numbers which give multiplicative hashes\n");
   fprintf(stdout,"with the best possible separations when given sequential values as\n");
   fprintf(stdout,"input.  This property of the golden ratio was known to botanists as\n");
   fprintf(stdout,"early as 1837 and seen in the rotational spacing between plant parts.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"We use a few tricks to implement the multiplicative Fibonacci hash\n"); 
   fprintf(stdout,"in LARC mult_golden_hash().  It will act on 64-bit integers and output\n"); 
   fprintf(stdout,"E-bit integers where E is the exponent of the hash table size S = 2^E.\n"); 
   fprintf(stdout,"We expect E to be no larger than 64, which means that when\n");
   fprintf(stdout,"implementing a multiplicative hash we only need to know at most\n");
   fprintf(stdout,"the first 64-bits to the right of the decimal place in our computations.\n");
   fprintf(stdout,"In particular, we only need to know the first 64-bits to the right of the\n");
   fprintf(stdout,"decimal of the golden ratio M = (sqrt(5) - 1)/2.  We precompute these\n");
   fprintf(stdout,"64 bits with multiprecision and store them as an integer\n");
   fprintf(stdout,"              11400714819323198485ULL .    \n");
   fprintf(stdout,"This is equal to the integer part of  [2^{64} * (sqrt(5) - 1)/2].\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"We were able to use this prestored value to implement LARC's\n");
   fprintf(stdout,"multiplicative Fibonacci hash entirely in terms of integers:\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"uint64_t mult_golden_hash(uint64_t key, uint64_t E  \n");
   fprintf(stdout,"{\n");
   fprintf(stdout,"       return (11400714819323198485ULL * key) >> (64 - E);\n");
   fprintf(stdout,"}\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"It takes a minute to convince oneself that this is equivalent to the\n");
   fprintf(stdout,"multiplicative Fibonancci hash H(key) which was described above as\n");
   fprintf(stdout,"   H(key) =  the integer part of [(2^E) * (fractional part of (M * key))] \n");
   fprintf(stdout,"                    where M =  (sqrt(5) - 1)/2].\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"First, remember that: \n");
   fprintf(stdout,"   11400714819323198485 =  the integer part of [2^{64} * M], \n");
   fprintf(stdout,"        where the golden mean M = (sqrt(5) - 1)/2] = .61803... , thus\n");
   fprintf(stdout,"   11400714819323198485 is less than 2^{64}, and we can fit it in\n");
   fprintf(stdout,"an unsigned 64-bit integer.\n");
   fprintf(stdout,"This corresponds to the first 64-bits of the fractional part of M.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"Because all the arguments are unsigned 64-bit integers, when we\n");
   fprintf(stdout,"compute the product P = (11400714819323198485ULL * key)\n");
   fprintf(stdout,"the result is only the low 64-bits of the integer result.\n");
   fprintf(stdout,"Thus it follows that this product as a uint64_t corresponds to\n");
   fprintf(stdout,"the first 64-bits of the fractional part of (M * key).\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"In LARC we assume  the hash exponent E can not be larger than 64,\n");
   fprintf(stdout,"so the number (64-E) is between 0 and 64.\n");
   fprintf(stdout,"The right shift operation >> (64-E) extracts the top E bits of P\n");
   fprintf(stdout,"which are the first E bits of the fractional part of (M * key).\n");
   fprintf(stdout,"Therefore,  the routine mult_golden_hash returns the intended value,\n");
   fprintf(stdout,"   =   the integer part of [(2^E) * (fractional part of (M * key))] \n");
   fprintf(stdout,"         where M is the golden mean (sqrt(5) - 1)/2].\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"Without having to use any floating point math, we have calculated\n");
   fprintf(stdout,"the multiplicative Fibonacci hash.  \n");
   fprintf(stdout,"\n");
   fprintf(stdout,"LARC can use mult_golden_hash(key, E) as an index into an \n");
   fprintf(stdout,"E-bit hash table.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"LARC often uses mult_golden_hash(key, 64), a hash function which\n");
   fprintf(stdout,"outputs a unsigned 64-bit integer, to build other hash functions.\n");
   fprintf(stdout,"\n");
}

// HASH EXPLANATION #3
void explain_hashing_lists() {
   /* This function describes how LARC hashes a list. */
   /* List hashing is used in the MatrixStore to retrieve a MatrixRecord */
   /* without knowing the MatrixID.  A clever choice of inputs to the hash */
   /* function allows this, while ensuring that only a single MatrixRecord */
   /* is made for each unique matrix. */
   /* List hashing is used in the OperationsStore to retrieve the results of */
   /* operations that have already been carried out before. */
   fprintf(stdout,"\n");
   fprintf(stdout,"       Hashing Lists with applications to MatrixStore and OperationsStore\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"LARC has an elegant and well-distributed fundamental hash called \n");
   fprintf(stdout,"mult_golden_hash(key, E) which takes as input a 64-bit integer\n");
   fprintf(stdout,"integer key, and using only integer math outputs an E-bit integer.\n");
   fprintf(stdout,"The details are given in explain_multiplicative_Fibonacci_hash().\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"LARC can use mult_golden_hash(key, E) as an index into an \n");
   fprintf(stdout,"E-bit hash table.  \n");
   fprintf(stdout,"\n");
   fprintf(stdout,"LARC uses mult_golden_hash(key, 64) to build a hash acting on\n");
   fprintf(stdout,"lists of 64-bit integers by alternating the hash with a bit-wise mod two\n");
   fprintf(stdout,"addition of the next integer in the list with the last hash result.\n");
   fprintf(stdout,"The final step is to use mult_golden_hash(last_sum,E) to produce\n");
   fprintf(stdout,"an E-bit hash index.  Here is an example:\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"uint64_t recursive_hash_from_three_integers\n");
   fprintf(stdout,"        (uint64_t int1, uint64_t int2, uint64_t int3, uint64_t E)\n");
   fprintf(stdout,"{\n");
   fprintf(stdout,"     uint64_t ping,pong;  \n");
   fprintf(stdout,"     ping = int1; 		  \n");
   fprintf(stdout,"     pong = mult_golden_hash(ping,64);   \n");
   fprintf(stdout,"     ping = pong ^ int2;	  \n");
   fprintf(stdout,"     pong = mult_golden_hash(ping,64);  \n");
   fprintf(stdout,"     ping = pong ^ int3;  \n");
   fprintf(stdout,"     pong = mult_golden_hash(ping,E); \n");			   
   fprintf(stdout,"     return pong;                       \n");
   fprintf(stdout,"}\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"A straightforward application of this hash of three integers is used\n");
   fprintf(stdout,"to memoize the results of operations in the OperationsStore.\n"); 
   fprintf(stdout,"\n");
   fprintf(stdout,"For example, say we want to multiply matrix A by matrix B.\n");
   fprintf(stdout,"We first check to see if the answer to this operation is already.\n");
   fprintf(stdout,"known by looking in the OperationsStore in the hash chain given by:\n");
   fprintf(stdout,"index = recursive_hash_from_three_integers(A_ID, B_ID, Num('mult')\n");
   fprintf(stdout,"       where A_ID is the MatrixID of matrix A,\n");
   fprintf(stdout,"                  B_ID is the MatrixID of matrix B, and\n");
   fprintf(stdout,"                  Num('mult') is an integer designating the multiply operation\n");
   fprintf(stdout,"If the appropriate record is found in the hash chain, then we return\n");
   fprintf(stdout,"the MatrixID of the result of this operation.   If no record is found,\n");
   fprintf(stdout,"then the computation is performed and a record is created.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"LARC also uses a hash of a list to index the MatrixStore. \n");
   fprintf(stdout,"This is a very powerful technique for LARC and allows us to\n");
   fprintf(stdout,"to ensure that we only save one copy of a matrix, and allow us\n");
   fprintf(stdout,"to identify quickly whether the result of a calculation is some\n");
   fprintf(stdout,"previously stored matrix without yet knowing the MatrixID of the result.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"The trick is based on the fact that LARC carries out its operations\n");
   fprintf(stdout,"recursively and matrices (larger than 1 by 1) are stored recursively\n");
   fprintf(stdout,"by simply listing the four MatrixIDs of the quadrant submatrices.\n");
   fprintf(stdout,"For example, a matrix M whose four quadrant submatrices are\n");
   fprintf(stdout,"[A, B; C, D] would have in its matrix record five different MatrixIDs:\n");
   fprintf(stdout,"     Record for M:    M_ID     [A_ID, B_ID, C_ID, D_ID]  \n");
   fprintf(stdout,"\n");
   fprintf(stdout,"LARC uses recursive_hash_from_four_integers(A_ID, B_ID, C_ID, D_ID)\n");
   fprintf(stdout,"as the index into the MatrixStore where the record for M will be stored.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"When an operation is carried out recursively and it is almost done\n");
   fprintf(stdout,"it will have identified that the solution S to the operation has four\n");
   fprintf(stdout,"quadrant submatrices with MatrixIDs  S1_ID, S2_ID, S3_ID, S4_ID\n");
   fprintf(stdout,"so it hashes this list and looks in the MatrixStore to see if such a\n");
   fprintf(stdout,"matrix already has a record.   If it does, the MatrixID of S is returned.\n");
   fprintf(stdout,"If there is no record, a new record is made and its MatrixID returned.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"Storing and retrieving matrices consisting of a single scalar value\n");
   fprintf(stdout,"require a whole different set of LARC tricks which are discussed in \n");
   fprintf(stdout,"explain_hashing_scalarTypes() and explain_hashing_for_snapping_scalars.\n");
   fprintf(stdout,"\n");
}

// HASH EXPLANATION #4
void explain_hashing_scalarTypes() {
   // This function describes how LARC hashes
   //  the various scalarTypes in a way that preserves the entropy
   //  of the underlying scalars.
   fprintf(stdout,"\n");
   fprintf(stdout,"              How LARC hashes various scalarTypes\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"We want well-behaved hashes for each scalarType that preserves\n");
   fprintf(stdout,"all the entropy of its scalars.  By this, we mean that modifying any\n");
   fprintf(stdout,"part of a scalar should essentially randomize the new hash value.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"LARC supports many different types of scalars including standard\n");
   fprintf(stdout,"types like integer, complex, rational, and multiprecision versions;\n");
   fprintf(stdout,"and some unusual types discussed in explain_scalarTypes() and \n");
   fprintf(stdout,"the paper in MyPyLARC/doc/LARCandMyPyLARC_ver2.1_2021.pdf.\n");
   fprintf(stdout,"The particular scalarType is selected when compiling, for example:\n");
   fprintf(stdout,"\tmake TYPE=COMPLEX\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"To illustrate the basic strategies used to hash the various scalarTypes\n");
   fprintf(stdout,"we will explain how to hash three scalarTypes:\n");
   fprintf(stdout,"INTEGER (via mult_golden_hash), REAL (via entropy extraction),\n");
   fprintf(stdout,"and multiprecision integer, MPINTEGER, (via list hashing).\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"Hashing INTEGER:\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"As was discussed in explain_multiplicative_Fibonacci_hash(),\n");
   fprintf(stdout,"LARC has an elegant and well-distributed fundamental hash called \n");
   fprintf(stdout,"mult_golden_hash(key, E) which takes as input a 64-bit integer\n");
   fprintf(stdout,"integer key, and using only integer math outputs an E-bit integer.\n");
   fprintf(stdout,"Usually this output is used as the index into an 2^E size hash table.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"Hashing REAL scalarType (via Extraction of Full Entropy):\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"LARC hashes a real longdouble by converting all the bits from \n");
   fprintf(stdout,"the mantissa, sign and exponent, into a single 64-bit integer, and then\n");
   fprintf(stdout,"applying mult_golden_hash(key, E) to return an E-bit integer.\n");
   fprintf(stdout,"This procedure implies that as long as our underlying hash is strong\n");
   fprintf(stdout,"that modifying any part of the real number input should essentially\n");
   fprintf(stdout,"randomize the output of the hash.  \n");
   fprintf(stdout,"The hash_from_one_longdouble() routine is in larc/src/hash.c.\n");
   fprintf(stdout,"It also handles the strange fact that zero in long double can have\n");
   fprintf(stdout,"either a 'positive' or 'negative' sign, hence different sign bits, and \n");
   fprintf(stdout,"we want both versions of zero to be treated the same way by the hash.\n");
   // fprintf(stdout,"Here is the routine:\n");
   // fprintf(stdout,"uint64_t hash_from_longdouble (long double ldoub, uint64_t E)\n");
   // fprintf(stdout,"{				   \n");
   // fprintf(stdout,"  uint64_t int_of_doub;   \n");
   // fprintf(stdout,"			    \n");
   // fprintf(stdout,"  // ensure that we aren't dealing with a negative zero    \n");
   // fprintf(stdout,"  if (fpclassify(ldoub) == FP_ZERO) ldoub = 0.0;	    \n");
   // fprintf(stdout,"										  \n");
   // fprintf(stdout,"  // interpret the bits of the double to be an integer	  \n");
   // fprintf(stdout,"  // int_of_doub = *((uint64_t*)&doub);	    \n");
   // fprintf(stdout,"  memcpy((void *) &int_of_doub, (void *) &ldoub, sizeof(int_of_doub));\n");
   // fprintf(stdout,"								    \n");
   // fprintf(stdout,"  return mult_golden_hash(int_of_doub, E);	    \n");
   //fprintf(stdout,"}								   \n");
   fprintf(stdout,"\n");
   fprintf(stdout,"Hashing multiprecision integers, MPINTEGER (via List Hashing): \n");
   fprintf(stdout,"\n");
   fprintf(stdout,"The scalarType MPINTEGER is a GMP multiprecision integer.\n");
   fprintf(stdout,"GMP integers are stored in (32 or) 64 bit words called `limbs'.\n");
   fprintf(stdout,"LARC converts the information from the MPINTEGER limbs into a\n");
   fprintf(stdout,"list of 64-bit integers.  It also uses one additional list item to handle\n");
   fprintf(stdout,"the sign of the MPINTEGER.\n");
   fprintf(stdout,"As we discussed in explain_hashing_lists(), LARC hashes lists by rounds\n");
   fprintf(stdout,"of applying mult_golden_hash(key, 64) to a 64-bit integer, alternated\n");
   fprintf(stdout,"with bit-wise mod 2 addition with the next 64-bit integer. On the final\n");
   fprintf(stdout,"step use mult_golden_hash(key, E) so we return an E-bit integer.\n");
   fprintf(stdout,"We apply list hashing to the list of integers that hold the information\n");
   fprintf(stdout,"from the MPINTEGER producing the desired hash function.\n");
   fprintf(stdout,"The larc_mpz_hash() routine is in larc/src/scalars.c.\n");
   fprintf(stdout,"This routine does enforce that different representations of zero\n");
   fprintf(stdout,"return the same value.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"Hashing Other ScalarTypes:\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"All other scalarTypes are hashed with combining three strategies:\n");
   fprintf(stdout,"entropy preservation, application of previously created hash functions,\n");
   fprintf(stdout,"and list hashing.  Here are a few examples.\n");
   fprintf(stdout,"A COMPLEX scalarType consists of two long doubles for the real and\n");
   fprintf(stdout,"imaginary parts, which we convert to two 64-bit integers preserving,\n");
   fprintf(stdout,"all bits of the long doubles, then we apply list hashing to these two\n");
   fprintf(stdout,"integers to create the desired COMPLEX scalarType hash.\n");
   fprintf(stdout,"For MPRATIONAL we use the larc_mpz_hash on the numerator and\n");
   fprintf(stdout,"denominator to produce two 64-bit integers, then apply list hashing.\n");
   fprintf(stdout,"For MPCOMPLEX we use larc_mpz_hash on the real and imaginary\n");
   fprintf(stdout,"parts and then apply list hashing.\n");
   fprintf(stdout,"The special LARC scalarTypes Clifford, Upper, and Lower\n");
   fprintf(stdout,"depend on structures, but the strategy is still the same.\n");
   fprintf(stdout,"We either translate information directly into 64-bit integers\n");
   fprintf(stdout,"or apply existing hash functions on some part of the information\n");
   fprintf(stdout,"to produce 64-bit integers, then take the full list of 64-bit integers\n");
   fprintf(stdout,"and apply list hashing.\n");
   fprintf(stdout,"\n");
/*  * ** INTEGER: do int_hash -> mult_golden_hash */
/*  * ** REAL: do hash_from_one_longdouble */
/*  * ** COMPLEX: do hash_from_two_longdoubles */
/*  * ** MPINTEGER: do larc_mpz_hash -> calls mult_golden_hash via list on a sign bit, */
/*  *    loops overXORing each limb in then hashing */
/*  * ** MPREAL: do larc_mpfr_hash -> tests for zero with mpfr_regular_p() */
/*  *    if zero, returns mult_golden_hash of zero, */
/*  *    otherwise, behaves as MPINTEGER */
/*  * ** MPCOMPLEX: do larc_mpfr_hash on real and imaginary parts, pass results */
/*  *    to recursive_hash_from_two_integers() */
/*  * ** MPRATIONAL: do larc_mpz_hash on numerator, denominator, pass results to */
/*  *    recursive_hash_from_two_integers() */
/*  * ** MPRATCOMPLEX: do mprational hash on real and imaginary parts, pass */
/*  *    results to recursive_hash_from_two_integers() */
/*    fprintf(stdout,"\tClifford:  a LARC structure with array of mpq_t basis coefficients\n"); */
/*    fprintf(stdout,"\tUpper and Lower:  a LARC structure used for approximate\n"); */
/*    fprintf(stdout,"\t                              calculations on probabilities.\n"); */
}


// HASH EXPLANATION #5
void explain_hashing_for_snapping_scalars() {
   // This function describes how LARC implements
   // snapping to a previously stored scalar that is 'close enough' to a
   // target scalar.   In SPR mode, the locality sensitive hash will, with
   // high probablility, snap to a good previously stored scalar.
   // In MAR mode, the locality sensitive hash is guaranteed to find a
   // scalar to snap to (if any that are sufficiently nearby are stored),
   // however, it might not snap to the one which is closest to the target
   // (if there is more than one stored scalar nearby).
   // 
   // NOTE:  This explanation is very similar to explain_scalar_techniques().
   //             This one has more on hash function specifics, the other has
   //             a section on mimicking symbolic computation.
   // 
   fprintf(stdout,"\n");
   fprintf(stdout,"             Locality Sensitive Hashing and Snapping Scalars\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"LARC reduces the spawning of nearly identical scalars (mostly caused\n");
   fprintf(stdout,"by finite precision issues) by 'snapping'. Snapping replaces a newly\n");
   fprintf(stdout,"calculated scalar with a sufficiently-nearby previously-stored scalar.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"The earliest version of LARC, when used on large applications,\n");
   fprintf(stdout,"required large amounts of memory in order to save sets of scalars that\n");
   fprintf(stdout,"were almost identical.  If the math operations had been carried out\n");
   fprintf(stdout,"symbolically the scalars would have been identical, but the finite\n");
   fprintf(stdout,"precision of floating point operations created tiny differences.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"We realized that if we could find a previously stored scalar that was\n");
   fprintf(stdout,"almost identical to a value that we had just calculated, that we could\n");
   fprintf(stdout,"substitute the previously stored scalar as a replacement for the new\n");
   fprintf(stdout,"value and prevent the spawning of nearly identical scalars (and the\n");
   fprintf(stdout,"corresponding spawning of nearly identical larger matrices).\n");
   fprintf(stdout,"Our first attempts to solve this problem required too much search time\n");
   fprintf(stdout,"since we might be searching a substantial portion of a large ScalarStore.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"We finally came up with a scheme that, with high probability, hashed \n");
   fprintf(stdout,"two nearly identical scalars to the same hash chain.  This meant\n");
   fprintf(stdout,"that we could greatly reduce spawning, simple by searching a single\n");
   fprintf(stdout,"hash chain in the ScalarStore.\n");
   fprintf(stdout,"We had rediscovered a hashing idea called 'locality sensitive hashing',\n");
   fprintf(stdout,"and then applied it to LARC as a way to 'snap' a newly calculated scalar\n");
   fprintf(stdout,"to a previously almost identical scalar (if one existed).\n");
   fprintf(stdout,"The basic idea is that the space of scalars is divided into tiny regions\n");
   fprintf(stdout,"and we only save the first scalar we find in a region in the scalar store.\n");
   fprintf(stdout,"This scalar becomes the *representative* of that region.  \n");
   fprintf(stdout,"Subsequently calculated scalars that would fall in a represented region\n");
   fprintf(stdout,"are not stored, and instead 'snapped  to' the value of the representative.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"We call this method Regional Representative Retrieval (RRR).\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"When LARC is initialized a parameter determines how close scalars\n");
   fprintf(stdout,"must be to each other in order to possibly be snapped together.\n");
   fprintf(stdout,"The parameter determines the size of tiny tiles. We have one method \n");
   fprintf(stdout,"that uses a single tile per region and another method that uses more. \n");
   fprintf(stdout,"\n");
   fprintf(stdout,"The LARC ScalarStore hash table uses a locality sensitive hash that\n");
   fprintf(stdout,"combines two functions, a function T which calculates a tile label for\n");
   fprintf(stdout,"a scalar S, and a function H which hashes tile labels. Then H(T(S)) is a\n");
   fprintf(stdout,"locality hash that hashes all scalars within a single tile to the same index.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"For example, if the tile function T truncates the value of S after 2 digits\n");
   fprintf(stdout,"to the right of the decimal place then the locality hash function H(T(S))\n");
   fprintf(stdout,"would store S= 3.1415 in hash chain C = H(3.14).  Then if a calculation\n");
   fprintf(stdout,"yielded S'=  3.1468, LARC would search hash chain C=H(3.14), find the\n");
   fprintf(stdout,"record containing S, and snap S' to this previously stored value.\n");
   fprintf(stdout,"The ScalarID of the stored record (if found) or new record (if created)\n");
   fprintf(stdout,"is returned as the 'answer' of the computation.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"LARC has two different schemes for storing and snapping scalars:\n");
   fprintf(stdout,"   - SPR (short for Single-tile Probabilistic Retrieval), and\n");
   fprintf(stdout,"   - MAR (short for Multi-tile Assured Retrieval).\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"SPR replaces a newly calculated scalar with a previously stored scalar\n");
   fprintf(stdout,"in the same tile if one exists, or otherwise stores the new scalar.\n");
   fprintf(stdout,"SPR succeeds with high probability in identifying a newly calculated\n");
   fprintf(stdout,"scalar with an almost-identical previously stored scalar, only failing\n");
   fprintf(stdout,"to collapse two close scalars that are on opposite sides of a tile\n");
   fprintf(stdout,"boundary.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"MAR uses a more complicated scheme in which adjacent tiles can be\n");
   fprintf(stdout,"grouped together to form a region.  This region will contain a single\n"); 
   fprintf(stdout,"stored scalar which lies someplace in the central portion of the region.\n");
   fprintf(stdout,"This trick is carried out by pasting together more than one tile.\n");
   fprintf(stdout,"MAR has a guarantee of snapping a newly calculated scalar to some\n");
   fprintf(stdout,"sufficiently-nearby previously-stored scalar (if it exists), but it may\n");
   fprintf(stdout,"not be the closest previous-stored scalar.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"We mentioned above that the locality hash function H(T()) is a \n");
   fprintf(stdout,"composition of a hash function H, and a tile label function T. \n");
   fprintf(stdout,"In SPR mode, the tile label will be the center point of the tile and \n");
   fprintf(stdout,"the hash function is the appropriate type of hash for the given\n");
   fprintf(stdout,"scalarType that LARC is compiled in. These the hash functions for\n");
   fprintf(stdout,"different scalarTypes are described in explain_hashing_scalarTypes().\n");
   fprintf(stdout,"In MAR mode, the tile label will an integer (or pair of integers for\n"); 
   fprintf(stdout,"complex types); you can think of these as indexing a grid.\n");
   fprintf(stdout,"Thus the hash function for MAR mode is the hash for MPINTEGER\n");
   fprintf(stdout,"or a hash built up from list hashing a pair of MPINTEGER hashes.\n");
   fprintf(stdout,"More details regional representive retrieval and thes two modes are in\n");
   fprintf(stdout,"the paper in MyPyLARC/doc/LARCandMyPyLARC_ver2.1_2021.pdf.\n");
   fprintf(stdout,"\n");
// left out - mimic symbolic computation
}

// HASH EXPLANATION #6
void explain_hash_filters() {
    //  This function describes how LARC uses hash filters
    //   LARC hash filters are small integers saved in the hash chains of
    //   the ScalarStore.  When traversing the hash chain searching for a
    //   previously stored scalar in the right region, an initial check is made
    //   to see whether the value of the filter is correct, before a more 
    //   expensive part of the search is carried out.
   fprintf(stdout,"\n");
   fprintf(stdout,"               Hash Filters to speed up ScalarStore searches\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"LARC uses a technique called 'hash filters' to speed searches\n"); 
   fprintf(stdout,"of the ScalarStore hash table for previously stored scalars.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"As described in explain_hashing_for_snapping_scalars(), LARC has a\n");
   fprintf(stdout,"scheme that divides up scalar space into tiny regions where it will \n");
   fprintf(stdout,"save at most one scalar called the representative of that region.\n");
   fprintf(stdout,"Any new scalar that would fall in the same region will be snapped to\n");
   fprintf(stdout,"the value of the representative scalar instead of being placed in a\n");
   fprintf(stdout,"new ScalarRecord. \n");
   fprintf(stdout,"\n");
   fprintf(stdout,"These regions are made up of tiles.  In SPR mode each region has a\n");
   fprintf(stdout,"single tile and if there is a representative of that tile, then the tile has\n");
   fprintf(stdout,"an associated ScalarRecord saved in the ScalarStore in the hash chain\n");
   fprintf(stdout,"indexed by the hash H(A) of the tile label A.  In MAR mode the regions \n");
   fprintf(stdout,"can have multiple tiles and each tile will have a ScalarRecord\n");
   fprintf(stdout,"stored in the hash chain determined by the hash of its tile label.\n");
   fprintf(stdout,"The one tile of the MAR region that contains the representative scalar\n");
   fprintf(stdout,"is called the primary region and its ScalarRecord has the saved scalar\n");
   fprintf(stdout,"as well as information that determines the other tiles in the region.\n");
   fprintf(stdout,"The tiles that are not primary tiles have ScalarRecords that contain\n");
   fprintf(stdout,"pointers back to the primary tile's record.  Thus whether we are in\n");
   fprintf(stdout,"MAR or SPR mode, when a new scalar S is calculated, we determine\n");
   fprintf(stdout,"the label L of the tile containing S, then we search the hash chain \n");
   fprintf(stdout,"labelled by the hash H(L) to see if there is a ScalarRecord for that\n");
   fprintf(stdout,"tile already saved in the ScalarStore. If there is such a ScalarRecord\n");
   fprintf(stdout,"then LARC will snap the new scalar to the stored scalar value that is\n");
   fprintf(stdout,"associated with the region containing this tile (by returning the ScalarID\n");
   fprintf(stdout,"of the stored scalar).  If there is not a ScalarRecord associated with L\n");
   fprintf(stdout,"in the store, then a new region is claimed with S as its representative\n");
   fprintf(stdout,"(and the newly assigned ScalarID is returned).    \n");
   fprintf(stdout,"    \n");
   fprintf(stdout,"When in MAR mode we will use the additional technique of hash filters.\n");
// is it only in MAR mode that the tile labels are multiprecision
// grid labels for the tiles and hence bad to store?							       
   fprintf(stdout,"    \n");
   fprintf(stdout,"Although we only need to look in the ScalarRecords along the \n");
   fprintf(stdout,"single hash chain with index H(L), this might involve a lot of work.\n");
   fprintf(stdout,"When the ScalarStore gets very full there may be many\n");
   fprintf(stdout,"tile labels A, B, C, etc.  that have ScalarRecords in the same\n");
   fprintf(stdout,"hash chain H(A) = H(B) = H(C).	\n");
   fprintf(stdout,"Furthermore, the check to see whether a ScalarRecord\n");
   fprintf(stdout,"has the same tile label as L actually takes a bit of computation\n");
   fprintf(stdout,"because the labels are not saved and must be recalculated.\n");
   fprintf(stdout,"This is where the hash filter check saves us time.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"As LARC traverses the hash chain before making the check to see \n");
   fprintf(stdout,"whether the current ScalarRecord has the same tile label as the \n");
   fprintf(stdout,"tile label L we are looking for, LARC first does an initial check to \n");
   fprintf(stdout,"see whether the ScalarRecord has the right hash filter value.\n");
   fprintf(stdout,"    \n");
   fprintf(stdout,"The hash filter value F(A) associated with a tile labeled A is a short\n");
   fprintf(stdout,"16-bit hash function that was calculated and saved when we first\n");
   fprintf(stdout,"made the ScalarRecord for the tile labelled A.  We check to see if\n");
   fprintf(stdout,"F(A) = F(L) and if it does not, then this is not the ScalarRecord we\n");
   fprintf(stdout,"are seeking and we can move to the next ScalarRecord in the chain.\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"Note that the hash function H() which determines the index of\n");
   fprintf(stdout,"the hash chain in the ScalarStore and the hash function F() which\n");
   fprintf(stdout,"determines the hash filter must be independent hash functions.\n");
   fprintf(stdout,"Every tile label A, B, C, ... with a ScalarRecord in the hash chain\n");
   fprintf(stdout,"(as we are searching for a tile with label L) has the property\n");
   fprintf(stdout,"that H(L) = H(A) = H(B) = H(C).  We would like it to be the case,\n");
   fprintf(stdout,"that the hash filter has the property that if A and B are different tiles\n");
   fprintf(stdout,"with H(A) = H(B), that there is only a one in 2^{16}-th chance\n");
   fprintf(stdout,"that A and B have the same hash filter values F(A) = F(B).\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"Thus LARC ScalarStore searches only have to search a single hash chain\n");
   fprintf(stdout,"and can skip the secondary work of checking whether tile labels \n");
   fprintf(stdout,"match (A = L) except for the rare cases in which the hash filter F(L)\n");
   fprintf(stdout,"for our new scalar in tile L, matches the prestored hash filter value\n");
   fprintf(stdout,"F(A) for the a ScalarRecord.\n");
   fprintf(stdout,"\n");
}


// HASH EXPLANATION #7
void explain_hash_statistics() {
    // explain_hash_statistics - describes how to get statistics on hash 
    //  chain length and what expected longest chain length should be.
   fprintf(stdout,"\n");
   fprintf(stdout,"   Getting Statistics on Hashing, What should you expect?\n");
   fprintf(stdout,"\n");
   fprintf(stdout,"LARC has a routine that will track statistics on the hash tables.\n");
   fprintf(stdout,"When the user edits larc/src/larc.h so that HASHSTATS is defined,\n");
   fprintf(stdout,"a program hashstats_to_files() will create files with statistics on\n");
   fprintf(stdout,"the LARC hash tables.  It tracks various statistics:\n");
   fprintf(stdout,"   (1) accesses[hash_val] = An array of the number of accesses \n");
   fprintf(stdout,"                 each hash bucket (hash_val) has so far. This  \n");
   fprintf(stdout,"                 is incremented by both hits and misses. \n");
   fprintf(stdout,"   (2) nodes[hash_val] = An array of the number of nodes (records)  \n");
   fprintf(stdout,"                currently in each hash chain.  These values are\n");
   fprintf(stdout,"                incremented by inserts, and decremented\n");
   fprintf(stdout,"                by removes.                     \n");
   fprintf(stdout,"and a standard hash report is sent to stdout or a file.    \n");
   fprintf(stdout,"\n");
   fprintf(stdout,"For users who would like to know more about the theoretical expected\n");
   fprintf(stdout,"length of hash chains in a well-behaved hash (which is what one \n");
   fprintf(stdout,"would expect for a truely random function), we highly recommend\n");
   fprintf(stdout,"   'Balls into Bins - A Simple and Tight Analysis'\n");
   fprintf(stdout,"    by Martin Raab and Angelika Steger,\n");
   fprintf(stdout,"which is in\n");
   fprintf(stdout,"   Springer Lecture Notes in Computer Science 1518\n");
   fprintf(stdout,"   'Randomization and Approximation Techniques in Computer Science'\n");
   fprintf(stdout,"    edited by Luby, Rolim, and Serna.\n");
   fprintf(stdout,"   \n");
/*
The statistic we were most
interested in was the maximum count. This should correlate to the statistic called maximum
occupancy in the problem in which m balls are tossed with uniform probability into n bins.
For this experiment m and n were equal. We verified that the maximum occupancy was
approximately log log log n n , which is the value predicted for random data for the case m = n [?].
We also used a random number generator to implement an experiment tossing m balls into
n bins (matching the m, n values from our hash experiment). We visually compared the two
distributions and saw a very close match.
*/
}



void explain_scalarType() {
    // this function explains the different scalarTypes available in LARC.
    fprintf(stdout,"\nSupported scalarTypes for LARC:\n");
    fprintf(stdout,"The user selects what type of scalars will be used inside of LARC\n");
    fprintf(stdout,"matrices by compiling with a command such as:\n");
    fprintf(stdout,"\tmake TYPE=INTEGER\n");
    fprintf(stdout,"This scalarType also determines the appropriate scalar arithmetic.\n");
    fprintf(stdout,"The scalarTypes supported by LARC are:\n");
    fprintf(stdout,"\tReal (default): C long double\n");
    fprintf(stdout,"\tInteger:        C int64_t \n");
    fprintf(stdout,"\tBoolean:        C int64_t \n");
    fprintf(stdout,"\tComplex:        C long double complex \n");
    fprintf(stdout,"\tMPInteger:      GMP multiprecisionf mpz_t \n");
    fprintf(stdout,"\tMPRational:     GMP multiprecison  mpq_t \n");
    fprintf(stdout,"\tMPReal:         GMP multiprecision mpfr_t \n");
    fprintf(stdout,"\tMPcomplex:      GMP multiprecison mpc_t \n");
    fprintf(stdout,"\tMPRatComplex:   a LARC structure with real and imag in mpq_t\n");
    fprintf(stdout,"\tClifford:  a LARC structure with array of mpq_t basis coefficients\n");
    fprintf(stdout,"\tUpper and Lower:  a LARC structure used for approximate\n");
    fprintf(stdout,"\t                              calculations on probabilities.\n");
}

void list_explanatory_resources() {
       fprintf(stdout,"**********************************************************************\n");
       fprintf(stdout,"*  RESOURCES:                                                        *\n");
       fprintf(stdout,"*   A detailed explanatory paper on the LARC (Linear Algebra via     *\n");
       fprintf(stdout,"*   Recursive Compression) package and the MyPyLARC Tutorial         *\n");
       fprintf(stdout,"*   and Sample Applications package is available at:                 *\n");
       fprintf(stdout,"*   GitHub.com/LARCmath/MyPyLARC/doc/LARCandMyPyLARC_version2_June2021.pdf *\n");
       fprintf(stdout,"*                                                                    *\n");
       fprintf(stdout,"*   Additional information can be found at the following locations:  *\n");
       fprintf(stdout,"*     MyPyLARC/Tutorial/README:  short tutorial routine descriptions *\n");
       fprintf(stdout,"*     MyPyLARC/Tutorial/newuser_instructions:                        *\n");
       fprintf(stdout,"*                                suggested order to do tutorials     *\n");
       fprintf(stdout,"*     MyPyLARC/README.md:        overview of MyPyLARC and LARC,      *\n");
       fprintf(stdout,"*                                compiling, matrix operations list   *\n");
       fprintf(stdout,"*     MyPyLARC/doc:              explanatory paper and poster        *\n");
       fprintf(stdout,"*     MyPyLARC/larc/doc:         intro slides, contributors list     *\n");
       fprintf(stdout,"*     MyPyLARC/html/index.html:  doxygen documentation to view       *\n");
       fprintf(stdout,"*                                in browser                          *\n");
       fprintf(stdout,"**********************************************************************\n");
}

const char* get_string_scalarType() {
    static char buf[1024];

    // write the scalar type to the info_store
#if defined(USE_INTEGER)
    snprintf(buf, sizeof(buf), "INTEGER");
#elif defined(USE_BOOLEAN)
    snprintf(buf, sizeof(buf), "BOOLEAN");
#elif defined(USE_COMPLEX)
    snprintf(buf, sizeof(buf), "COMPLEX");
#elif defined(USE_REAL)
    snprintf(buf, sizeof(buf), "REAL");
#elif defined(USE_MPINTEGER)
    snprintf(buf, sizeof(buf), "MPINTEGER");
#elif defined(USE_MPRATIONAL)
    snprintf(buf, sizeof(buf), "MPRATIONAL");
#elif defined(USE_MPRATCOMPLEX)
    snprintf(buf, sizeof(buf), "MPRATCOMPLEX");
#elif defined(USE_MPREAL)
    snprintf(buf, sizeof(buf), "MPREAL");
#elif defined(USE_MPCOMPLEX)
    snprintf(buf, sizeof(buf), "MPCOMPLEX");
#elif defined(USE_CLIFFORD)
#ifdef IS_COMPLEX
    snprintf(buf, sizeof(buf), "COMPLEX CLIFFORD: %s", clifford_description);
#else
    snprintf(buf, sizeof(buf), "REAL CLIFFORD: %s", clifford_description);
#endif  // #ifdef IS_COMPLEX
#elif defined(USE_LOWER)
    snprintf(buf, sizeof(buf), "LOWER %d", NUM_EXPONENTS);
#elif defined(USE_UPPER)
    snprintf(buf, sizeof(buf), "UPPER %d", NUM_EXPONENTS);
#else
    snprintf(buf, sizeof(buf), "Unknown Scalar Type");
#endif  // #if defined(USE_INTEGER)

    return buf;
}

const char* get_string_larc_version() {
    static char ver_buf[20000];
    ver_buf[0] = '\0';

    char line_buf[1024];
    snprintf(line_buf, sizeof(line_buf), "  LARC v%d.%d.%d\n", MAJOR_VERSION, MINOR_VERSION, PATCH_VERSION);
    strcat(ver_buf, line_buf);
#ifdef GIT_COMMIT_DATE
    snprintf(line_buf, sizeof(line_buf), "  Commit Date: %s\n", GIT_COMMIT_DATE);
    strcat(ver_buf, line_buf);
#endif // #ifdef GIT_COMMIT_DATE
#ifdef __GNUC__
    snprintf(line_buf, sizeof(line_buf), "  GCC version used: %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
    strcat(ver_buf, line_buf);
#endif // #ifdef __GNUC__
#ifdef __GLIBC__
    snprintf(line_buf, sizeof(line_buf), "  GLIBC version used: %d.%d\n", __GLIBC__, __GLIBC_MINOR__);
    strcat(ver_buf, line_buf);
#endif // #ifdef __GLIBC__
    snprintf(line_buf, sizeof(line_buf), "  NCURSES version used: %s\n", NCURSES_VERSION);
    strcat(ver_buf, line_buf);

    snprintf(line_buf, sizeof(line_buf), "  The scalarType is currently %s.\n", get_string_scalarType());
    strcat(ver_buf, line_buf);

    strcat(ver_buf, get_string_gmp_versions());

#ifdef MAR
    snprintf(line_buf, sizeof(line_buf), "  The locality hash used is: MAR (Multi-tile Assured Retrieval).\n");
    strcat(ver_buf, line_buf);
#else
    snprintf(line_buf, sizeof(line_buf), "  The locality hash used is: SPR (Single-tile Probabilistic Retrieval).\n");
    strcat(ver_buf, line_buf);
#endif  // #ifdef MAR

#ifdef HASHSTATS
    strcat(ver_buf, "  The HASHSTATS flag is enabled; detailed hash statistics are being gathered.\n");
#else
    strcat(ver_buf, "  The HASHSTATS flag is disabled; detailed hash statistics are NOT being gathered.\n");
#endif  // #ifdef HASHSTATS

    return ver_buf;
}

void print_larc_version() {
    // This function prints to screen information about the version of LARC 
    // being used and the current scalarType (set at compile time)
    printf("%s", get_string_larc_version());
}

/*!
 * \ingroup larc
 * \brief A utility to print LARC initialization parameters 
 *
 * This routine outputs the LARC initialization parameters to screen, along
 * with a small amount of information on what each parameter influences.
 *
 * \param matrix_store_exp The base2 log of the size of the matrix store
 * \param op_store_exp The base2 log of the size of the operations store
 * \param max_level The base2 log of the largest dimension of any matrix to be created
 * \param regionbitparam The number of bits of precision kept in SPR regions
 * \param zeroregionbitparam The negative base2 log of the smallest number not in the SPR region of zero. MAR ignores this and sets it automatically.
 * 
 */
static void create_output_strings_to_print_settings(size_t matrix_store_exp, size_t op_store_exp, mat_level_t max_level, int regionbitparam, int zeroregionbitparam) {
    printf("  ================================================================\n");
    printf("  Initializing LARC stores with:\n");
    printf("    Matrix store hash exponent %zd\n",matrix_store_exp);
    printf("    Operations store hash exp  %zd\n",op_store_exp);
    printf("    Maximum matrix level %d\n",max_level);
    printf("    locality-sensitive hashing parameters are:\n");
    printf("       regionbitparam = %d; ",regionbitparam);
    printf("       zeroregionbitparam = %d\n",zeroregionbitparam);
    printf("    (see documentation for what these do)\n");
}

void print_larc_parameters_and_compile_options() {
  // This function calls five access functions to retrieve the parameters
  // passed to LARC at initialization, then passes this information to 
  // create_output_strings_to_print_settings to print them to the screen. See the
  // latter function for descriptions of the parameters.
  create_output_strings_to_print_settings(
      get_nonscalar_store_exp(),
      get_op_store_exp(),
      max_level_allowed_matrixStore(),
      get_regionbitparam(),
      get_zeroregionbitparam());
}

int64_t memory_available_GiB()
{
  // get memory statistics from standard location /proc/meminfo
  FILE *fp = fopen("/proc/meminfo", "r");
  int64_t mem_total=0, mem_free=0, mem_avail=0, cached=0, mem_avail_GiB=0;
  
  if (fp)
    {
      char buf[100];
      // read first line of file /proc/meminfo
      if (fgets(buf, 100, fp))
	{
	  if (sscanf(buf, "MemTotal: %" SCNd64 " kB", &mem_total) < 1) {
	    fprintf(stderr,"WARNING in %s not reading MemTotal\n",__func__);
	  }
	}
      // read second line of file /proc/meminfo
      if (fgets(buf, 100, fp))
	{
	  if (sscanf(buf, "MemFree: %" SCNd64 " kB", &mem_free) < 1) {
	    fprintf(stderr,"WARNING in %s not reading MemFree\n",__func__);
	  }
	}
       // read third line of file /proc/meminfo
      if (fgets(buf, 100, fp))
	{
	  if (sscanf(buf, "MemAvailable: %" SCNd64 " kB", &mem_avail) < 1) {
	    fprintf(stderr,"WARNING in %s not reading MemAvailable\n",__func__);
	  }
	}
      // skip fourth line 
      (void)! fgets(buf, 100, fp);
      // read fifth line of file /proc/meminfo
      if (fgets(buf, 100, fp))
	{
	  if (sscanf(buf, "Cached: %" SCNd64 " kB", &cached) < 1) {
	    fprintf(stderr,"WARNING in %s not reading Cached\n",__func__);
	  }
	}
      fclose(fp);
    }
  else if (VERBOSE>SILENT)
    {
      printf("WARNING in %s: could not open /proc/meminfo\n",__func__);
      printf("Could not retrieve information on available memory.\n");
    }

  // printf("The memory available is %ld\n", mem_avail);

  // if (VERBOSE==ALL) {
  // printf("Mem Total = %ld Mem Free=%ld, Cached=%ld, Mem Avail=%ld\n",mem_total,mem_free,cached, mem_avail);
  // }

  // convert from kB to Gigibytes by 2^20 factor, that is 20 bits.
  mem_avail_GiB = mem_avail >> 20;
  
  return(mem_avail_GiB);
}


// The main initialization function for LARC; it creates storage and preloads
// select scalars and matrices into the matrixStore. See Doxygen documentation
// in organize.h for descriptions of the parameters.
void initialize_larc(size_t matrix_store_exp, size_t op_store_exp, mat_level_t max_level, int regionbitparam, int zeroregionbitparam, int verbose)
{
  if (verbose > SILENT) {
#ifdef MAR
    fprintf(stdout,"LARC will be running in MAR mode.\n");
#else 
    fprintf(stdout,"LARC will be running in SPR mode.\n");
#endif
  }
  
  // Error  checking
  // region parameters must be in set {-1 (default), 0, 1, 2, ...}
  if (regionbitparam < -1) {
    fprintf(stderr,"Invalid entry for regionbitparam. Default of %d used.\n",
	    regionbitparam = DEFAULT_REGIONBITPARAM);
  }
  if (zeroregionbitparam < -1) {
    fprintf(stderr,"Invalid entry for zeroregionbitparam. Default used.\n");
  }
  
  // Region parameter defaults:
  // regionbitparam default 
  if (regionbitparam == -1) {
    regionbitparam = DEFAULT_REGIONBITPARAM;
  }
  // zeroregionbitparam for SPR
#ifndef MAR
  if (zeroregionbitparam == -1) {
    zeroregionbitparam = regionbitparam;
  }
#endif
// zeroregionbitparam for MAR
#ifdef MAR
  zeroregionbitparam = regionbitparam;
#endif
  
  // Error checking for SPR
  // regions near real or imag zero must be as large or larger than other regions
  // hence zeroregionbitparam must be <= regionbitparam. 
#ifndef MAR
  if (zeroregionbitparam > regionbitparam) {
    zeroregionbitparam = regionbitparam;
  }
#endif
 
 // warnings about machine precison   
#if defined(USE_COMPLEX) || defined(USE_REAL)  
 if (regionbitparam > LDBL_MANT_DIG)
   {
     fprintf(stderr,"WARNING: regionbitparam value is %d which is greater\n",regionbitparam);
     fprintf(stderr,"than machine precision of %d, so every scalar will\n",LDBL_MANT_DIG);
     fprintf(stderr,"be in its own locality hash region for hashing. Machine\n");
     fprintf(stderr,"precision may cause many nearly equal scalars\n");
     fprintf(stderr,"to be stored, and LARC not to work well.\n");
   }
#endif
   
   // if verbose parameter outside allowed range, set VERBOSE to default level
   // (BASIC); otherwise set VERBOSE accordingly
   verbose_type_t vbot = SILENT; // should be min
   verbose_type_t vtop = ALL;  // should be max
   if (verbose>=(int)vbot && verbose <=(int)vtop )
     { VERBOSE = (verbose_type_t) verbose; }
   else {VERBOSE = BASIC; }
   
   if (VERBOSE==CHATTY || VERBOSE ==ALL) {
     introduce_LARC_and_MyPyLARC();
     explain_matrix_and_operation_storage();
     explain_level();
     explain_scalar_techniques();
     explain_scalarType();
     explain_verbosity();
     list_explanatory_resources();
  }  

  if (VERBOSE>SILENT) {
    print_larc_version();
    create_output_strings_to_print_settings(matrix_store_exp, op_store_exp, max_level,
                                    regionbitparam, zeroregionbitparam);
    printf("\n");
  }
  
  // Initialize scalar operations.
  init_scalarOps(verbose);

  // Initialize common use variables in globals.c
  // This must come after scalar ops are set. 
  scratchVars_init();

  if (verbose == DEBUG)  printf("about to create matrix store\n");
  // Create the matrix store
  if (!create_matrix_store(matrix_store_exp, max_level, regionbitparam, zeroregionbitparam)) {
    fprintf(stderr,"Unable to create matrix store. Try a smaller exponent than %zd.\n", matrix_store_exp);
    exit(1);
  }

  // Create the operations store
  if (verbose == DEBUG) printf("about to create ops store\n");
  if (!create_op_store(op_store_exp)) {
    fprintf(stderr,"Unable to create operations store. Try a smaller exponent than %zd.\n", op_store_exp);
    exit(1);
  }

  // Create the info store for meta data (use sparingly)
  size_t info_store_exp = 10;   // default for now
  if (verbose == DEBUG)  printf("about to create info store\n");
  if (!create_info_store(info_store_exp)) {
    fprintf(stderr,"Unable to create info store. Try a smaller exponent than %zd.\n", info_store_exp);
    exit(1);
  }

  // Preload basic scalars, matrices for matrix math+linear algebra+circuits
  // preload zero and identity matrices (requires matrix store)
  // preload integer Hadamard matrices (also requires op store)
  if (verbose == DEBUG) printf("about to preload scalars\n");
  if (!preload_matrix_and_scalar_stores()) {
    fprintf(stderr,"Unable to preload matrix store.\n");
    exit(1);
  }
  if (verbose == DEBUG)  printf("finished preload of scalars\n");

  if (VERBOSE>BASIC) {
    printf("    Created matrix store and operations store,\n");
    printf("       and preloaded basic scalars and matrices.\n");
    printf("  ================================================================\n\n");
  }

  // Play nicely with other computer users with create_stopHogging_thread
  // which will kill off our program if we are hogging almost all of the memory
  // for some reason timeout must be static for this to work properly
  static unsigned int timeout = 600; // in seconds: 
  int s;
  s = pthread_create(&stopHogging_thread, NULL, stopHogging, &timeout);
  if (s != 0) { errno = s; perror("pthread_create"); exit(EXIT_FAILURE); }
  if (verbose==DEBUG)  printf("end of initialize_larc()\n");  
}

void shutdown_larc(void)
{
        // This routine clears all data and frees all memory that was
        // allocated for LARC. It is only needed for programs which
        // want to initialize a fresh copy of LARC without exiting first.
        stop_thread(&stopHogging_thread);
	if (report_thread_exists)
	{
	        stop_thread(&report_thread);
		report_thread_exists = 0;
	}
        scratchVars_clear();
        free_info_store();
        free_op_store();
        free_matrix_store();
}
      
/* This test was written when we were having an issue with setting
 * the -z parameter for FFT 2^3. We needed -z 53
 * or smaller, however this value zeros out very small scalars, e.g.
 * 1/2^60.  */
void precision_testing() {
  printf("\n**************************************\n");
  printf("The value of DBL_MANT_DIG is %d\n",DBL_MANT_DIG);
  printf("The value of LDBL_MANT_DIG is %d\n",LDBL_MANT_DIG);
  printf("**************************************\n\n");
  // scalarType smallScalarType = pow(2,-60);
  // long double smallLongDouble = powl(2,-61);
  // int64_t f_ID = create_FFTMat(3);
  // print_naive(f_ID);
  return;
}

