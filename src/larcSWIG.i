%module larcSWIG

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

/*
  The following functions all return char * (string) type from malloc'd
  memory and are declared in headers in this .i file. This is how you
  tell SWIG it's ok for the wrapper functions to free the memory.
*/
%newobject count_entries_larcMatrixFile;              /* experimental.h */
%newobject info_get;                                  /* info_store.h */
%newobject create_log_dir;                            /* larc.h */
%newobject sca_get_str;                               /* larc.h */
%newobject matrix_maxnorm_matrixID;                   /* matmath.h */
%newobject matrix_l2norm_matrixID;                    /* matmath.h */
%newobject matrix_tracenorm_matrixID;                 /* matmath.h */
%newobject matrix_trace_from_matrixID;                /* matmath.h */
%newobject scalar_string_from_matrixID;               /* matmath.h */
%newobject matrix_count_entries_matrixID;             /* matmath.h */
%newobject matrix_list_scalars_larcMatrix;            /* matmath.h */
%newobject get_valString_from_matID_and_coords;       /* matrix_store.h */
%newobject matrix_trace_matrixID;                     /* matrix_store.h */
/* 
  The following functions return char * (string) type, but not from
  malloc'd memory, so SWIG should not try to free it!

info_store.h: return_info_name;

  The following functions return char * (string) type from malloc'd
  memory, but their headers are not visible to SWIG and they should
  not be included in the %newobject declarations.

json.h: j_make_rel_filename
logging.c: get_current_dir_name (this is a unistd.h function)
*/

%{
#include "type.h"
#include "matmath.h"
#include "larc.h"
#include "io.h"
#include "global.h" 
#include "matrix_store.h"
#include "op_store.h"
#include "organize.h"
#include "info_store.h"
#include "fft.h"
#include "hash.h"
#include "scalars.h"
#include "experimental.h"
#include <complex.h>
#include <gmp.h>
#include <pthread.h>
%}

#define SWIGWORDSIZE64
// #enddef

%include "stdint.i"
%include "complex.i"
%include "cstring.i"
%include "carrays.i"
%include "cpointer.i"


/* This tells SWIG to output int64_t ** as list of length two integer lists. */
/*
%typemap(out) int64_t ** {
    int len, i;
    len = 0;
    while ($1[len]) len++;
    $result = PyList_New(len);
    for (i = 0; i < len; i++){
        PyObject *l = PyList_New(2);
        PyList_SetItem(l, 0, PyLong_FromLong($1[i][0]));
        PyList_SetItem(l, 1, PyLong_FromLong($1[i][1]));
        PyList_SetItem($result, i, l);
    }
}
*/

/* This tells SWIG to output int64_t ** as list of integer lists. 
   In this setup, -1's mark the end of each list and a NULL marks the end
   of the final list */
/* As of now, locate_entries_larcMatrixFile is the only routine that makes use of this.*/
%typemap(out) int64_t ** {
    int len;
    int i, j;
    len = 0;
    while (NULL != $1[len]) len++;
    $result = PyList_New(len);
    for (i = 0; i < len; i++){
        j = 0;
        PyObject *l = PyList_New(0);
        while (0 <= $1[i][j]){
            PyList_Append(l, PyLong_FromLong($1[i][j]));
            j ++;
        }
        PyList_SetItem($result, i, l);
    }
}

/* This tells SWIG to input list of strings as char **. */
%typemap(in) char ** {
    /* Check if is a list */
    if (PyList_Check($input)){
        int size = PyList_Size($input);
        int i = 0;
        $1 = (char **) malloc((size+1)*sizeof(char *));
        for (i = 0; i < size; i++){
            PyObject *o = PyList_GetItem($input, i);
            if (PyUnicode_Check(o))
                $1[i] = PyBytes_AsString(PyUnicode_AsEncodedString(o, "utf-8", "strict"));
            else {
                PyErr_SetString(PyExc_TypeError, "list must contain strings");
                free($1);
                return NULL;
            }
        }
        $1[i] = 0;
    }
    else {
        PyErr_SetString(PyExc_TypeError, "argument is not a list");
        return NULL;
    }
}

/* This cleans up the char ** array we malloc'd before the function call. */
%typemap(freearg) char ** {
    free((char *) $1);
}


// # NOTE: global.h uses the USE_INTEGER/REAL/COMPLEX that is #defined in
// # type.h, so type.h must be before it in the list

%include "type.h"
%include "matmath.h"
%include "io.h"
%include "larc.h"
%include "global.h"
%include "matrix_store.h"
%include "op_store.h"
%include "organize.h"
%include "fft.h"
%include "hash.h"
%include "scalars.h"
%include "info_store.h"
%include "experimental.h"

%array_class(int, intArray);
%array_class(long int, int64Array); // # works because SWIGWORDSIZE64 defined
%array_class(double, realArray);
%array_class(complex, complexArray);

/*
# Why buildArray doesn't work anymore and what we did:
# we think we were having problems with scalarType not being known
# by SWIG generated python, and so we had problems with the routine
# scalarTypeArray_setitem(newarr, i, val) since val was scalarType
# so we starting working around doing buildArray, by using map and
# other functions that python knows and we removed the following line:
# %array_functions(scalarType, scalarTypeArray);
# This caused all the jupyter notebook files to stop working, but
# since only a few people use them it was not noticed until 5-2019.
*/
