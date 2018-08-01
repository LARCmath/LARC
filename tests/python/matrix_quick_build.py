#################################################################
#                                                               #
# Copyright 2014, Institute for Defense Analyses                #
# 4850 Mark Center Drive, Alexandria, VA; 703-845-2500          #
# This material may be reproduced by or for the US Government   #
# pursuant to the copyright license under the clauses at DFARS  #
# 252.227-7013 and 252.227-7014.                                #
#                                                               #
# LARC : Linear Algebra via Recursive Compression               #
# Authors:                                                      #
#   - Steve Cuccaro (IDA-CCS)                                   #
#   - John Daly (LPS)                                           #
#   - John Gilbert (UCSB, IDA adjunct)                          #
#   - Jenny Zito (IDA-CCS)                                      #
#                                                               #
# Additional contributors are listed in "LARCcontributors".     #
#                                                               #
# POC: Jennifer Zito <jszito@super.org>                         #
# Please contact the POC before disseminating this code.        #
#                                                               #
#################################################################

import os
import sys
import numpy as np
import random
sys.path.append(os.path.join(os.path.dirname(__file__),"../../src"))
from pylarc import *
from ctypes import *

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
        randMat_ID = row_major_list_to_store_matrixID(buildArray(randVals), row_level, col_level, 2**col_level)
    except TypeError:
        raise TypeError("argument 1 of type 'ScalarType'")

    #return matrix id
    return randMat_ID

def matrix_zero_matrixID(row_level, col_level):
    """
    Builds the zero matrix of the requested size.
    """
    n = (2**row_level)*(2**col_level)
    vals = [0 for i in range(n)]

    try: 
        zeroMat_ID = row_major_list_to_store_matrixID(buildArray(vals), row_level, col_level, 2**col_level)
    except TypeError:
        raise TypeError("argument 1 of type `ScalarType'")

    return zeroMat_ID

def matrix_identity_matrixID(level):
    """
    Builds a square identity matrix of requested level.
    """
    n = (2**level)*(2**level)
    vals = [1*(0==(i%(2**level+1))) for i in range(n)]

    try:
        idMat_ID = row_major_list_to_store_matrixID(buildArray(vals), level, level, 2**level)
    except TypeError:
        raise TypeError("argument 1 of type 'ScalarType'")

    return idMat_ID

