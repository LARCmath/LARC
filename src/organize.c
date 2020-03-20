//                        organize.c 
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

// Our header files structures and functions
#include "larc.h"
#include "matrix_store.h"
#include "op_store.h"
#include "global.h"
#include "organize.h"
#include "version.h"
#include "scalars.h"
#include "info_store.h"
#include "matmath.h"


void explain_verbosity() {
  // this function just prints a message to explain the different levels of verbosity
    fprintf(stdout,"\nThe different levels of verbosity are:\n");
    fprintf(stdout,"\tSILENT:  verbose = 0    only prints errors;\n");
    fprintf(stdout,"\tBASIC:   verbose = 1    prints errors, warnings, and minimal information;\n");
    fprintf(stdout,"\tCHATTY:  verbose = 2    adds informative comments;\n");
    fprintf(stdout,"\tDEBUG:   verbose = 3    adds debugging comments;\n");
    fprintf(stdout,"\tALL:     verbose = 4    simply won't shut up.\n");
}


void explain_matrix_and_operation_storage() {
  // this function just prints a message to explain how LARC
  // stores matrix records and operations records
  // also see explain_locality_sensitive_hashing()
    fprintf(stdout,"\nStoring Matrix Records and Operation Records:\n");
    fprintf(stdout,"In the LARC (linear algebra via recursive compression) package, each\n");
    fprintf(stdout,"scalar, vector, and matrix is stored exactly once as a matrix-record\n");
    fprintf(stdout,"in the LARC matrix-record hash table and is assigned a unique\n");
    fprintf(stdout,"matrixID.  In matrix records, a matrix is defined recursively by a\n");
    fprintf(stdout,"list of the four matrixIDs corresponding to its four quadrant\n");
    fprintf(stdout,"submatrices.  The matrix-store hash table has length 2^(mat_exp) and\n");
    fprintf(stdout,"is indexed by a hash of the four quadrant submatrix matrixIDs (from\n");
    fprintf(stdout,"its recursive definition), except when the matrix is a scalar, in\n");
    fprintf(stdout,"which case a special locality-sensitive hash is used.  The\n");
    fprintf(stdout,"operations-memoization hash table has length 2^(op_exp) and is indexed\n");
    fprintf(stdout,"by a hash containing all the matrixIDs of inputs to the operation and\n");
    fprintf(stdout,"also the name of the operation.\n");
}

  
void explain_locality_sensitive_hashing() {
  // this function just prints a message to explain locality sensitive
  // hashing and identification of scalars within regions.
    fprintf(stdout,"\nLocality Sensitive Hashing and Identifying Scalars within each Region:\n");
    fprintf(stdout,"In LARC the space of scalars is divided into regions of size\n");
    fprintf(stdout,"1/2^(region_bits), and the first scalar stored in LARC that lies in a\n");
    fprintf(stdout,"region becomes the unique representative in that region.  That is, any\n");
    fprintf(stdout,"new value that would lie in an already occupied region is treated as\n");
    fprintf(stdout,"equivalent to the originally stored representative and is not stored.\n");
    fprintf(stdout,"Since zero is particularly important to mathematical identities it is\n");
    fprintf(stdout,"pre-stored as the representative of its region and any scalar within a\n");
    fprintf(stdout,"distance of 1/2^(max{bits_to_zero,region_bits}) of zero will be\n");
    fprintf(stdout,"treated as zero. The locality-sensitive hash insures that scalars that\n");
    fprintf(stdout,"would be in the same region hash to the same value, making it fast to\n");
    fprintf(stdout,"find a previously stored representative.\n");
}


void explain_level() {
  // this function just prints a message to explain levels
    fprintf(stdout,"\nMatrix Dimensions:\n");
    fprintf(stdout,"Because LARC defines matrices recursively in terms of quadrant\n");
    fprintf(stdout,"submatrices, each dimension of a matrix in LARC must be a power of\n");
    fprintf(stdout,"two. We use level to be the exponent of two describing the dimension\n");
    fprintf(stdout,"of the matrix; the size of a matrix is 2^(row_level) by 2^(col_level).\n");
    fprintf(stdout,"For efficiency, the max_level of matrices is selected at initialization,\n");
    fprintf(stdout,"so that matrices can be no bigger than 2^(max_level) by 2^(max_level).\n");
}


void explain_scalarType() {
  // this function just prints a message to explain locality sensitive
  // hashing and identification of scalars within regions.
    fprintf(stdout,"\nLARC scalarType:\n");
    fprintf(stdout,"The user selects what type of scalars will be used inside of LARC\n");
    fprintf(stdout,"matrices by compiling with a command such as:\n");
    fprintf(stdout,"\tmake TYPE=INTEGER\n");
    fprintf(stdout,"This scalarType also determines the appropriate scalar arithmetic.\n");
    fprintf(stdout,"The scalarTypes supported by LARC are:\n");
    fprintf(stdout,"\tReal (default): C long double\n");
    fprintf(stdout,"\tInteger:        C int64_t \n");
    fprintf(stdout,"\tComplex:        C long double complex \n");
    fprintf(stdout,"\tMPInteger:      multiprecisionf mpz_t \n");
    fprintf(stdout,"\tMPRational:     multiprecison  mpq_t \n");
    fprintf(stdout,"\tMPReal:         multiprecision mpfr_t \n");
    fprintf(stdout,"\tMPcomplex:      multiprecison mpc_t \n");
    fprintf(stdout,"\tMPRatComplex:   a LARC structure with real and imag in mpq_t\n");
}



void initialize_larc(size_t matrix_store_exp, size_t op_store_exp, mat_level_t max_level, int sighash, int zerobitthresh, int verbose)
{
  if (sighash == -1) {sighash = SIGHASH_DEFAULT;}
#if defined(USE_COMPLEX) || defined(USE_REAL)  
  if (sighash > LDBL_MANT_DIG)
  {
    fprintf(stderr,"WARNING: sighash value is %d which is greater\n",sighash);
    fprintf(stderr,"than machine precision of %d, so no approximation\n",LDBL_MANT_DIG);
    fprintf(stderr,"occurs via neighborhood hashing and machine\n");
    fprintf(stderr,"precision may cause many nearly equal scalars\n");
    fprintf(stderr,"to be stored, and LARC not to work well.\n");
  }
#endif  
  if (zerobitthresh == -1) {zerobitthresh = ZEROBITTHRESH_DEFAULT;}
   
  // if verbose parameter outside allowed range, set VERBOSE to default level
  // (BASIC); otherwise set VERBOSE accordingly
  verbose_type_t vbot = SILENT; // should be min
  verbose_type_t vtop = DEBUG;  // should be max
  if (verbose>=(int)vbot && verbose <=(int)vtop )
    { VERBOSE = (verbose_type_t) verbose; }
  else {VERBOSE = BASIC; }

  if (VERBOSE==CHATTY || VERBOSE ==ALL) {
    explain_matrix_and_operation_storage();
    explain_locality_sensitive_hashing();
    explain_level();
    explain_scalarType();
    explain_verbosity();
  }  

  if (VERBOSE>SILENT) {
    printf("  LARC v%d.%d.%d\n", MAJOR_VERSION, MINOR_VERSION, PATCH_VERSION);
    printf("  ================================================================\n");
    printf("  Initializing LARC stores with:\n");
    printf("    Matrix store hash exponent %zd\n",matrix_store_exp);
    printf("    Operations store hash exp  %zd\n",op_store_exp);
    printf("    Maximum matrix level %d\n",max_level);
    printf("    locality-approximation constants are:\n");
    printf("       round to %d significant bits\n",sighash);
    printf("       truncate to zero if within %d bits of zero\n",zerobitthresh);
#ifdef USE_INTEGER 
    printf("  The scalarType is currently INTEGER.\n\n");
#endif
#ifdef USE_COMPLEX 
    printf("  The scalarType is currently COMPLEX.\n\n");
#endif
#ifdef USE_REAL 
   printf("  The scalarType is currently REAL.\n\n");
#endif
#ifdef USE_MPINTEGER
   printf("  The scalarType is currently MPINTEGER.\n\n");
#endif
#ifdef USE_MPRATIONAL
   printf("  The scalarType is currently MPRATIONAL.\n\n");
#endif
#ifdef USE_MPRATCOMPLEX
   printf("  The scalarType is currently MPRATCOMPLEX.\n\n");
#endif
#ifdef USE_MPREAL
   printf("  The scalarType is currently MPREAL.\n\n");
#endif
#ifdef USE_MPCOMPLEX
   printf("  The scalarType is currently MPCOMPLEX.\n\n");
#endif
  }
  
  // Check to see if scalar operations have already been set by user. If yes,
  // continue. If partially, throw error and exit. Otherwise, use a default
  // set of operations. 
  int ops_set = check_scalarOps();
  if (ops_set == 0){
      if (VERBOSE>BASIC)
          printf("    No scalar ops defined. Defaulting.\n    ");
      init_arithmetic_scalarOps(verbose);
  }
  else if (ops_set == 1){
      if (VERBOSE>BASIC)
          printf("    Scalar ops predefined by user.\n");
  }
  else {
      fprintf(stderr,"\nERROR: some, but not all, scalar ops defined.\n");
      exit(1);
  }

  // Initialize common use variables in globals.c
  // This must come after scalar ops are set. 
  scratchVars_init();

  // Create the matrix store
  if (!create_matrix_store(matrix_store_exp, max_level, sighash, zerobitthresh)) {
    fprintf(stderr,"Unable to create matrix store. Try a smaller exponent than %zd.\n", matrix_store_exp);
    exit(1);
  }

  // Create the operations store
  if (!create_op_store(op_store_exp)) {
    fprintf(stderr,"Unable to create operations store. Try a smaller exponent than %zd.\n", op_store_exp);
    exit(1);
  }

  // Create the info store for meta data (use sparingly)
  size_t info_store_exp = 10;   // default for now
  if (!create_info_store(info_store_exp)) {
    fprintf(stderr,"Unable to create info store. Try a smaller exponent than %zd.\n", info_store_exp);
    exit(1);
  }

  // Preload basic scalars, matrices for matrix math+linear algebra+circuits
  // preload zero and identity matrices (requires matrix store)
  // preload integer Hadamard matrices (also requires op store)
  if (!preload_matrix_store()) {
    fprintf(stderr,"Unable to preload matrix store.\n");
    exit(1);
  }

  if (VERBOSE>BASIC) {
    printf("    Created matrix store and operations store,\n");
    printf("       and preloaded basic scalars and matrices.\n");
    printf("  ================================================================\n\n");
  }

  // Play nicely with other computer users with create_seppuku_thread
  // which will kill off our program if we are hogging almost all of the memory
  // for some reason timeout must be static for this to work properly
  static unsigned int timeout = 600; // in seconds: 
  static pthread_t thd;
  pthread_create(&thd, NULL, seppuku, &timeout);
}


void create_report_thread(uint64_t period) 
{
  static pthread_t thd;
  static uint64_t timeout;
  timeout = period;
  pthread_create(&thd, NULL, report, &timeout);
}


void rusage_report(int tabs, char *outfilepath)
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
	
  // STILL CONTAINS Seppuku function (which is also in seppuku)
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
#else
      if (VERBOSE>BASIC)
        fprintf(f,"\nHASHSTATS is off, turn on for lengths of hash chains.\n");
#endif      
    }
  if (fp)
    {
      char buf[100];
      
      // read first line of file /proc/meminfo
      if (fgets(buf, 100, fp))
	{
	  sscanf(buf, "MemTotal: %ld kB", &mem_total);
	}
      // read second line of file /proc/meminfo
      if (fgets(buf, 100, fp))
	{
	  sscanf(buf, "MemFree: %ld kB", &mem_free);
	}
      fgets(buf, 100, fp);
      // read fourth line of file /proc/meminfo
      if (fgets(buf, 100, fp))
	{
	  sscanf(buf, "Cached: %ld kB", &cached);
	}
      fclose(fp);
    }
  fprintf(f,"%sMem Free=%ld, Cached=%ld, Mem Total=%ld\n", tabs?"\t":"", mem_free,cached, mem_total);
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

void *report(void *arg)
{
	uint64_t timeout = *(uint64_t *)arg;
        printf("\n\n=== Starting Periodic Reporting with Interval of %" PRIu64 " Seconds ===\n",timeout);

	while (1) {
		sleep(timeout);
		printf("\n\n==================== Periodic Report (%" PRIu64 " sec) ====================\n", timeout);
                printf("NOTE: the main thread continues to run while the periodic report is\n");
                printf("generated, so number of stored matrices may increase!");
		printf("\n========================================================================\n\n");
		report_now("stdout");
		// matrix_store_report();
		// op_store_report();
		// rusage_report(0,"stdout");
		printf("\n========================================================================\n\n");
	}
}

void report_now(char *outfilepath)
{
		matrix_store_report(outfilepath);
		op_store_report(outfilepath);
		rusage_report(0,outfilepath);
}

void *seppuku(void *arg)
{
  unsigned int timeout = *(unsigned int*)arg;
  struct rusage rusage;
  int tabs = 0;
  
  while (1) {
    sleep(timeout);
    printf("Running seppuku check\n");

    // get memory statistics from standard location /proc/meminfo
    FILE *fp = fopen("/proc/meminfo", "r");
    uint64_t mem_total=0, mem_free=0, cached=0;
    if (fp)
      {
	char buf[100];
        // read first line of file /proc/meminfo
	if (fgets(buf, 100, fp))
	  {
	    sscanf(buf, "MemTotal: %ld kB", &mem_total);
	  }
        // read second line of file /proc/meminfo
	if (fgets(buf, 100, fp))
	  {
	    sscanf(buf, "MemFree: %ld kB", &mem_free);
	  }
	fgets(buf, 100, fp);
        // read fourth line of file /proc/meminfo
	if (fgets(buf, 100, fp))
	  {
	    sscanf(buf, "Cached: %ld kB", &cached);
	  }
	fclose(fp);
      }
    else if (VERBOSE>SILENT)
      {
        printf("WARNING in seppuku check: could not open /proc/meminfo\n");
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
        report_now("stdout");
        printf("%s Mem Free=%ld, Cached=%ld, Mem Total=%ld\n", tabs?"\t":"", mem_free,cached, mem_total);
        printf("\n%sTotal Free/Memory: %g/%gMB\n", tabs?"\t":"", (mem_free+cached)/1000.0, mem_total/1000.0);
        printf(">>>Too little \"real\" memory left, exiting<<<\n");
        exit(1);
      }
  }
}


/* We are having an issue with setting the -z parameter
 * for FFT 2^3 with current precision, we need -z 53
 * or smaller, however this value zeros out very small scalars, e.g.
 * 1/2^60.  Thus we would like to develop a higher 
 * precision capability that would allow both situations
 * to correctly compute with the same -z parameter. */
void precision_testing() {
  printf("\n**************************************\n");
  printf("The value of DBL_MANT_DIG is %d\n",DBL_MANT_DIG);
  printf("The value of LDBL_MANT_DIG is %d\n",LDBL_MANT_DIG);
  printf("**************************************\n\n");
  // scalarType smallScalarType = pow(2,-60);
  // long double smallLongDouble = powl(2,-61);
  // mat_ptr_t f_PTR = create_fft_matrix(3);
  // print_naive_by_matPTR(f_PTR);
  return;
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
	  if (sscanf(buf, "MemTotal: %ld kB", &mem_total) < 1) {
	    fprintf(stderr,"WARNING in %s not reading MemTotal\n",__func__);
	  }
	}
      // read second line of file /proc/meminfo
      if (fgets(buf, 100, fp))
	{
	  if (sscanf(buf, "MemFree: %ld kB", &mem_free) < 1) {
	    fprintf(stderr,"WARNING in %s not reading MemFree\n",__func__);
	  }
	}
       // read third line of file /proc/meminfo
      if (fgets(buf, 100, fp))
	{
	  if (sscanf(buf, "MemAvailable: %ld kB", &mem_avail) < 1) {
	    fprintf(stderr,"WARNING in %s not reading MemAvailable\n",__func__);
	  }
	}
      // skip fourth line 
      fgets(buf, 100, fp);
      // read fifth line of file /proc/meminfo
      if (fgets(buf, 100, fp))
	{
	  if (sscanf(buf, "Cached: %ld kB", &cached) < 1) {
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


