//                        organize.c 
/******************************************************************
 *                                                                *
 * Copyright 2014, Institute for Defense Analyses                 *
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
 * POC: Jennifer Zito <jszito@super.org>                          *
 * Please contact the POC before disseminating this code.         *
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
#include "info_store.h"


void initialize_larc(size_t matrix_store_exp, size_t op_store_exp, mat_level_t max_level, 
                     int sighash, int zerobitthresh)
{
  if (sighash == -1) {sighash = SIGHASH_DEFAULT;}
  if (sighash > DBL_MANT_DIG)
  { printf("WARNING: sighash value of %d will do nothing\n",sighash); }
  if (zerobitthresh == -1) {zerobitthresh = ZEROBITTHRESH_DEFAULT;}
  int verbose = 1;
  if (verbose) {
    printf("  LARC v%d.%d.%d\n", MAJOR_VERSION, MINOR_VERSION, PATCH_VERSION);
    printf("  ===============================================\n");
    printf("  Initializing LARC stores with:\n");
    printf("    Matrix store hash exponent %zd\n",matrix_store_exp);
    printf("    Operations store hash exp  %zd\n",op_store_exp);
    printf("    Maximum matrix level %d\n",max_level);
    printf("    locality-approximation constants are:\n");
    printf("       round to %d significant bits\n",sighash);
    printf("       truncate to zero if within %d bits of zero\n",zerobitthresh);
#ifdef USE_INTEGER 
    printf("  The scalarType is currently INTEGER.\n");
#endif
#ifdef USE_COMPLEX 
    printf("  The scalarType is currently COMPLEX.\n");
#endif
#ifdef USE_REAL 
   printf("  The scalarType is currently REAL.\n");
#endif
  }
  
  // Create the matrix store
  if (!create_matrix_store(matrix_store_exp, max_level, sighash, zerobitthresh)) {
    printf("Unable to create matrix store. Try a smaller exponent than %zd.\n", matrix_store_exp);
    exit(1);
  }

  // Create the operations store
  if (!create_op_store(op_store_exp)) {
    printf("Unable to create operations store. Try a smaller exponent than %zd.\n", op_store_exp);
    exit(1);
  }

  // Create the info store for meta data (use sparingly)
  size_t info_store_exp = 10;   // default for now
  if (!create_info_store(info_store_exp)) {
    printf("Unable to create info store. Try a smaller exponent than %zd.\n", info_store_exp);
    exit(1);
  }

  // Preload basic scalars, matrices for matrix math, linear algebra and circuits
  // preload zero and identity matrices (requires matrix store)
  // preload integer Hadamard matrices (also requires op store)
  if (!preload_matrix_store()) {
    printf("Unable to preload matrix store.\n");
    exit(1);
  }


  if (verbose) {
    printf("    Created matrix store and operations store,\n");
    printf("       and preloaded basic scalars and matrices.\n");
    printf("  ===============================================\n\n");
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
      fprintf(f,"\nWARNING: HASHSTATS on, using more time and memory!\n");
#else
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
  // mat_add_t f_PTR = create_fft_matrix(3);
  // print_matrix_naive(f_PTR);
  return;
}
