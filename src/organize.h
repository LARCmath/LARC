//                       organize.h
/******************************************************************
 *                                                                *
 * Copyright 2014, Institute for Defense Analyses                 *
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
 * POC: Jennifer Zito <jszito@super.org>                          *
 * Please contact the POC before disseminating this code.         *
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

void create_report_thread(uint64_t period); 
void *report(void *arg);
void rusage_report(int tabs, char *outfilepath);
void report_now(char *outfilepath);

void *seppuku(void *arg);

void initialize_larc(size_t matrix_store_exp, size_t op_store_exp, mat_level_t max_level, 
                     int sighash, int zerobitthresh);


void precision_testing();

#ifdef __cplusplus
}
#endif

#endif
