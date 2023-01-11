//                         logging.c
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


// Standard Libaries
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/types.h>

// Our header files structures and functions
#include "larc.h"
#include "matrix_store.h"
#include "op_store.h"
#include "global.h"

/*!
 * \file logging.c
 * \brief Routines for user logging of projects - still preliminary
 */

/*!
 * \ingroup larc
 * \brief Returns the current working directory.
 *
 * Code copied from: https://www.gnu.org/software/libc/manual/html_node/Working-Directory.html
 */
char * gnu_getcwd_copy ()
{
  size_t size = 100;

  while (1)
    {
      char *buffer = (char *) malloc (size);
      if (getcwd (buffer, size) == buffer)
        return buffer;
      free (buffer);
      if (errno != ERANGE)
        return 0;
      size *= 2;
    }
}

/* create_log_dir
   Authors: Jenny and Steve
   This function creates a directory named <log_dir>/YYYYMMDD/YYYYMMDD.hhmmss
   to be used to store a copy of the input and all the output for this run.
   It allocates space for a directory name string, then returns a pointer to 
   the string.  The top level directory <log_dir> is placed inside the location
   the main program is called from and is by default "log" if this
   function is called with empty string instead of log_name being 
   passed as an argument.
 */

/*!
 * \ingroup larc
 * \brief Creates in the current directory a directory named [log_name]/YYYYMMDD/YYYYMMDD.hhmmss
 * 
 * The directory is created within the directory that the main program is 
 * called from. If one or both of the directories above it ([log_name] and
 * [log_name]/YYYYMMDD) do not exist, they are created first.
 *
 * \param log_name A string with the name to be used for the top-level directory (if null, the default name is 'log')
 * \return A string containing the path to the bottom-level log directory
 */
char * create_log_dir(char * log_name)
{

  char *log_dir;
  if ((strcmp(log_name, "") == 0) || (log_name == NULL)) log_dir = "log";
  else     log_dir = log_name;


        // Grab the current directory path
        // char *current_path =  malloc(1024);
        char *current_path = gnu_getcwd_copy();
        // printf("LOG: The current path is %s, and has length %d\n",
        //        current_path,strlen(current_path));

        // make the path to the top log directory
        // char *log_path =  malloc(strlen(current_path) + strlen(log_dir) + 1);
        char *log_path =  malloc(strlen(current_path) + strlen(log_dir) + 2);
        if (log_path == NULL) { ALLOCFAIL(); }
	strcpy(log_path,current_path);
        strcat(log_path,"/");
        strcat(log_path,log_dir);
        // printf("LOG: The log directory path is %s, and has length %d\n",
        //        log_path,strlen(log_path));
        mkdir(log_path,0777);

        // Create a time stamp for YYMMDD
        struct tm time_struct;
        time_t epoch_time;
        time(&epoch_time);
        gmtime_r(&epoch_time,&time_struct);
        unsigned yyyy,mo,hh,dd,mm,ss;
        yyyy = 1900 + time_struct.tm_year;
        mo = 1 + time_struct.tm_mon;
        dd = time_struct.tm_mday;
        hh = time_struct.tm_hour;
        mm = time_struct.tm_min;
        ss = time_struct.tm_sec;

        // create name of day-stamp subdirectory inside main log directory
        int size_time_buffers = 256;
        char *time_buffer_day = malloc(256);
        if (time_buffer_day == NULL) { ALLOCFAIL(); }
        char *time_format_day="%4u%02u%02u";
        snprintf(time_buffer_day,size_time_buffers,time_format_day,yyyy,mo,dd);
        // printf("LOG: The time string is %s, and has length %d\n",
        //        time_buffer_day,strlen(time_buffer_day));


        // make a sub directory in the log, named with the day stamp        
        char *day_log_path =  malloc(strlen(log_path) + strlen(time_buffer_day) + 2);
        if (day_log_path == NULL) { ALLOCFAIL(); }
	day_log_path = strcpy(day_log_path,log_path);
        strcat(day_log_path,"/");
        strcat(day_log_path,time_buffer_day);
        // printf("LOG: The log sub-directory path is %s, and has length %d\n",
        //       day_log_path,strlen(day_log_path));
        mkdir(day_log_path,0777);

  
        // Create a bottom level logging subdirectory with complete time stamp
        char *time_buffer_full = malloc(256);
        if (time_buffer_full == NULL) { ALLOCFAIL(); }
        char *time_format_full="%4u%02u%02u.%02u%02u%02u";
        snprintf(time_buffer_full,size_time_buffers,time_format_full,yyyy,mo,dd,hh,mm,ss);
        // printf("LOG: The time string is %s, and has length %d\n",
        //       time_buffer_full,strlen(time_buffer_full));

        // make a sub sub directory in the log, named with the time stamp        
        char *time_log_path =  malloc(strlen(day_log_path) + strlen(time_buffer_full) + 2);
        if (time_log_path == NULL) { ALLOCFAIL(); }
	time_log_path = strcpy(time_log_path,day_log_path);
        strcat(time_log_path,"/");
        strcat(time_log_path,time_buffer_full);
        // printf("LOG: The log sub-sub-directory path is %s, and has length %d\n",
        //       time_log_path,strlen(time_log_path));
        mkdir(time_log_path,0777);


        // Another way to set permissions
	//        int dir_status = mkdir(log_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
 

        // Free everything
        free(current_path); 
        free(log_path); 
        free(time_buffer_day);
        free(day_log_path); 
        free(time_buffer_full);
        
        return(time_log_path);

}
