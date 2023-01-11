#!/usr/bin/env python3

 #*################################################################
 #                                                                #
 # Copyright (C) 2014, Institute for Defense Analyses             #
 # 4850 Mark Center Drive, Alexandria, VA; 703-845-2500           #
 # This material may be reproduced by or for the US Government    #
 # pursuant to the copyright license under the clauses at DFARS   #
 # 252.227-7013 and 252.227-7014.                                 #
 #                                                                #
 # LARC : Linear Algebra via Recursive Compression                #
 # Authors:                                                       #
 #   - Steve Cuccaro (IDA-CCS)                                    #
 #   - John Daly (LPS)                                            #
 #   - John Gilbert (UCSB, IDA adjunct)                           #
 #   - Mark Pleszkoch (IDA-CCS)                                   #
 #   - Jenny Zito (IDA-CCS)                                       #
 #                                                                #
 # Additional contributors are listed in "LARCcontributors".      #
 #                                                                #
 # Questions: larc@super.org                                      #
 #                                                                #
 # All rights reserved.                                           #
 #                                                                #
 # Redistribution and use in source and binary forms, with or     #
 # without modification, are permitted provided that the          #
 # following conditions are met:                                  #
 #   - Redistribution of source code must retain the above        #
 #     copyright notice, this list of conditions and the          #
 #     following disclaimer.                                      #
 #   - Redistribution in binary form must reproduce the above     #
 #     copyright notice, this list of conditions and the          #
 #     following disclaimer in the documentation and/or other     #
 #     materials provided with the distribution.                  #
 #   - Neither the name of the copyright holder nor the names of  #
 #     its contributors may be used to endorse or promote         #
 #     products derived from this software without specific prior #
 #     written permission.                                        #
 #                                                                #
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND         #
 # CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,    #
 # INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF       #
 # MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE       #
 # DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER NOR        #
 # CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,   #
 # SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT   #
 # NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;   #
 # LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)       #
 # HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN      #
 # CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR   #
 # OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, #
 # EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.             #
 #                                                                #
 #*################################################################

from __future__ import print_function

import os
import sys
import numpy as np
import random
sys.path.append(os.path.join(os.path.dirname(__file__),"../../src"))
import pylarc
from ctypes import *

# Define string for using in formating filenames
scalarTypeStr = pylarc.cvar.scalarTypeStr

def rand_real(range_start, range_stop):
    """
    Generates a random real number by using random.random to generate a random x in [0,1)
    and transforming it to be y between range_start and range_stop: 
        y = x * (range_stop - range_start) + range_start
    """
    return random.random()*(range_stop - range_start) + range_start

def rand_complex(range_start, range_stop):
    """
    Generates a random complex number of the form a+bj where a,b are real numbers
    in [range_start, range_stop). 
    """
    return np.complex(rand_real(range_start, range_stop), rand_real(range_start, range_stop))

def matrix_random_matrixID(scalar_type, row_level, col_level, range_start, range_stop, sparsity=None):
    """
    Generate a matrix of size m x n, where m = 2^row_level and n = 2^col_level, with entries 
    in [range_start, range_stop). If a sparsity value in [0,1] is given, rounds 
    sparsity * m * n to the closest integer value and generates the matrix such that there 
    are exactly that many (randomly placed) zeros in the matrix. 
    """
    #size of matrix
    n = (2**row_level)*(2**col_level)

    # COMPLEX CASE
    if scalar_type.lower() in ["c", "complex"]:
        rand_f = rand_complex

    # INTEGER CASE
    elif scalar_type.lower() in ["i", "integer"]:
        rand_f = random.randrange
        #check that start and stop are integers
        if type(range_start) != int or type(range_stop) != int:
            raise ValueError("range_start and range_stop must be integers when scalarType = INTEGER!")

    # REAL CASE
    elif scalar_type.lower() in ["r", "real"]:
        rand_f = rand_real

    else:
        raise Exception("Do not know how to build matrix for type %s." %scalar_type)

    #apply sparsity requirements if present
    if sparsity is not None:
        #calculate required number of zeros from sparsity
        numZerosReq = int(round(sparsity * n))
        #randomly choose zero locations
        zeroPos = random.sample(range(n), numZerosReq)
        #generate random values in desired range, placing zeros
        randVals = []
        for i in range(n):
            if i in zeroPos:
                randVals.append(0)
            else:
                val = 0
                while val == 0:
                    val = rand_f(range_start, range_stop)
                randVals.append(val)
    #otherwise, just generate random values in desired range
    else: 
        randVals = [rand_f(range_start, range_stop) for i in range(n)]
    #convert list to a matrix and add to store
    try:
        randMat_ID = row_major_list_to_store(pylarc.map_to_str(randVals, scalarTypeStr), row_level, col_level, 2**col_level)
    except TypeError:
        raise TypeError("argument 1 of type 'ScalarType'")

    #return matrix id
    return randMat_ID

# If you want a zero matrix, there is a routine in matrix_store.c for that
# int64_t get_zero_pID(mat_level_t row_level, mat_level_t col_level);
# Don't use the old code from Andy for matrix_zero_matrixID.

# def matrix_zero_matrixID(row_level, col_level):
#     """
#     Builds the zero matrix of the requested size and returns the matrixID.
#     """
#     n = (2**row_level)*(2**col_level)
#     vals = [0 for i in range(n)]

#     try: 
#         zeroMat_ID = row_major_list_to_store(pylarc.map_to_str(vals, scalarTypeStr), row_level, col_level, 2**col_level)
#     except TypeError:
#         raise TypeError("argument 1 of type `ScalarType'")

#     return zeroMat_ID


# If you want an identity  matrix, there is a routine in matrix_store.c for that
# int64_t get_identity_pID(mat_level_t level);
# Don't use the old code from Andy for matrix_identity_matrixID.

# def matrix_identity_matrixID(level):
#     """
#     Builds a square identity matrix of requested level and returns the matrixID.
#     """
#     n = (2**level)*(2**level)
#     vals = [1*(0==(i%(2**level+1))) for i in range(n)]

#     try:
#         idMat_ID = row_major_list_to_store(pylarc.map_to_str(vals, scalarTypeStr), level, level, 2**level)
#     except TypeError:
#         raise TypeError("argument 1 of type 'ScalarType'")

#     return idMat_ID

