#              larc_utilities.py
#*################################################################
#                                                                #
# Copyright (C) 2014-2024, Institute for Defense Analyses        #
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

from __future__ import print_function, division

## \file larc_utilities.py
#  \brief This module contains Python utility functions for LARC.
#
#  We also use this file to dynamically import the module of SWIG-generated
#  Python wrappers for the LARC C code. It's necessary to import the wrappers
#  here as some of the utilities need to call LARC functions. We do the import
#  dynamically because each package built on top of LARC will have its own set
#  of SWIG-generated Python wrappers.
#
#  We have found that including more than one module of these wrappers in a
#  project build can cause a bug where the LARC MatrixStore (and other static
#  memory) gets duplicated, causing hard-to-diagnose errors when one
#  MatrixStore structure is initialized and the other is accessed. The dynamic
#  import of a single complete set of wrappers from each individual project
#  prevents this bug.
#

import os
import sys
import random
import argparse
current_directory = os.path.dirname(__file__)
# gv.py is in larc/src, one directory above this one
sys.path.append(os.path.join(current_directory,'..'))
import gv
import importlib
from ctypes import *
import numpy

# the following code emulates 'from module import *' for a dynamically loaded
# module (with the module name expressed as a string), assuming that the SWIG
# module being imported has no '__all__' attribute and no submodules. If either
# of these assumptions becomes invalid, the more complicated form below may be
# necessary.
sys.path.append(gv.PATH_TO_SWIG)
mdl = importlib.import_module(gv.SWIG_TO_USE)
for attr in dir(mdl):
    if not attr.startswith('_'):
        globals()[attr] = getattr(mdl, attr)

#mdl = __import__(gv.SWIG_TO_USE, fromlist=['*'])
# is there an __all__?
#if "__all__" in mdl.__dict__:
#    # if so respect it
#    names = mdl.__dict__["__all__"]
#else:
#    # otherwise we import all names that don't begin with _
#    names = [x for x in mdl.__dict__ if not x.startswith("_")]
#now drag names into global namespace
#globals().update({k: getattr(mdl, k) for k in names})

# utility functions begin here

## \brief Construct a string containing the exact binary represention
#         (mantissa and exp) of a double value.
#  \param dbl The double value to express as a string.
#  \result The string representing the double value 'dbl'.
def double_to_hex_string(dbl):
    """
    See Doxygen documentation.
    """
    assert type(dbl) == type(1.0)
    return dbl.hex()

## \brief Construct a string containing the exact rational representation
#         (numerator and denominator) of a double value.
#         Note that the denominator will always be a power of two.
#  \param dbl The double value to express as a string.
#  \result The string representing the double value 'dbl'.
def double_to_frac_string(dbl):
    """
    See Doxygen documentation.
    """
    assert type(dbl) == type(1.0)
    num, den = dbl.as_integer_ratio()
    return "{0}/{1}".format(num, den)

## \brief This routine converts Python values to strings suitable for parsing
#         by the LARC C code.
#
#  This is most needed for complex types, as Python complex numbers have the
#  format (1+2j) and C complex numbers the format (1+I*2). All floating point
#  numbers are also converted to a binary format to ensure that machine
#  precision is preserved. The type() command is used to determine how to do
#  this conversion. For exact types (e.g INTEGER, MPRATIONAL) binary is not
#  necessary.
#
#  \param val A Python variable
#  \param scalar_type A string that determines the format of the output
#  \result A string (suitable to pass to C routines through SWIG) which 
#         is an exact representation of the Python variable value
def value_to_string(val, scalar_type):
    """
    Returns a string representation of a Python value in a format
    appropriate for the given scalar_type; converts floating point
    numbers to a binary format to preserve machine precision, and changes
    python complex format (1+2j) to something more appropriate for C (1+I*2).
    """
    # COMPLEX CASE
    if scalar_type == "Complex":
        if type(val) == type(1j):
            return "{0}+I*{1}".format(double_to_hex_string(val.real), double_to_hex_string(val.imag))
        elif type(val) == type(1.0):
            return "{0}".format(double_to_hex_string(val))
        elif type(val) == type(1):
            return "{0}".format(val)
        else:
            assert False, "Invalid value {0} of type {1} for COMPLEX".format(val, type(val))

    # REAL CASE
    elif scalar_type == "Real":
        if type(val) == type(1.0):
            return "{0}".format(double_to_hex_string(val))
        elif type(val) == type(1):
            return "{0}".format(val)
        else:
            assert False, "Invalid value {0} of type {1} for REAL".format(val, type(val))

    # INTEGER CASE
    elif scalar_type == "Integer":
        if type(val) == type(1):
            return "{0}".format(val)
        else:
            assert False, "Invalid value {0} of type {1} for INTEGER".format(val, type(val))

    # BOOLEAN CASE
    elif scalar_type == "Boolean":
        if type(val) == type(True):
            return "{0}".format(1 if val else 0)
        elif type(val) == type(1):
            return "{0}".format(val)
        else:
            assert False, "Invalid value {0} of type {1} for BOOLEAN".format(val, type(val))

    # MPINTEGER CASE
    elif scalar_type == "MPInteger":
        if type(val) == type(1):
            return "{0}".format(val)
        else:
            assert False, "Invalid value {0} of type {1} for MPINTEGER".format(val, type(val))

    # MPRATIONAL CASE
    elif scalar_type == "MPRational":
        if type(val) == type(1.0):
            return "{0}".format(double_to_frac_string(val))
        elif type(val) == type(1):
            return "{0}".format(val)
        else:
            assert False, "Invalid value {0} of type {1} for MPRATIONAL".format(val, type(val))

    # MPRATCOMPLEX CASE
    elif scalar_type == "MPRatComplex":
        if type(val) == type(1j):
            return "{0}+I*{1}".format(double_to_frac_string(val.real), double_to_frac_string(val.imag))
        elif type(val) == type(1.0):
            return "{0}".format(double_to_frac_string(val))
        elif type(val) == type(1):
            return "{0}".format(val)
        else:
            assert False, "Invalid value {0} of type {1} for MPRATCOMPLEX".format(val, type(val))

    # MPREAL CASE
    elif scalar_type == "MPReal":
        if type(val) == type(1.0):
            return "{0}".format(double_to_hex_string(val))
        elif type(val) == type(1):
            return "{0}".format(val)
        else:
            assert False, "Invalid value {0} of type {1} for MPREAL".format(val, type(val))

    # MPCOMPLEX CASE
    elif scalar_type == "MPComplex":
        if type(val) == type(1j):
            return "{0}+I*{1}".format(double_to_hex_string(val.real), double_to_hex_string(val.imag))
        elif type(val) == type(1.0):
            return "{0}".format(double_to_hex_string(val))
        elif type(val) == type(1):
            return "{0}".format(val)
        else:
            assert False, "Invalid value {0} of type {1} for MPCOMPLEX".format(val, type(val))

    # CLIFFORD CASE
    elif scalar_type == "Clifford":
        if type(val) == type(1j):
            return "{0}+I*{1}".format(double_to_frac_string(val.real), double_to_frac_string(val.imag))
        elif type(val) == type(1.0):
            return "{0}".format(double_to_frac_string(val))
        elif type(val) == type(1):
            return "{0}".format(val)
        else:
            assert False, "Invalid value {0} of type {1} for CLIFFORD".format(val, type(val))

    # UPPER or LOWER CASE
    elif scalar_type in ("Upper", "Lower"):
        if type(val) == type(1.0):
            return "{0}".format(double_to_frac_string(val))
        elif type(val) == type(1):
            return "{0}".format(val)
        else:
            assert False, "Invalid value {0} of type {1} for BOUNDING type".format(val, type(val))

    # OTHER
    else:
        assert False, "Unhandled scalar type: {0}".format(scalar_type)

## \brief Apply value_to_string function to each element of a list.
# When in a floating point format (e.g. REAL, COMPLEX),
# this routine will use binary format in order to preserve
# full machine precision.
# For exact types (e.g INTEGER, MPRATIONAL) binary is not necessary.
#  \param vlist A list of python variables, all of the same type (integer, floating point, or complex)
#  \param scalar_type A string that determines the format of the output
#  \return A list of strings (see value_to_string function for details)
def map_to_str(vlist, scalar_type):
    """
    Apply value_to_string to each element of a list.
    """
    return [value_to_string(val, scalar_type) for val in vlist]

## \brief Generate a random real number in [range_start,range_stop)
#  \param range_start The smaller end of the real number range
#  \param range_stop The larger end of the real number range
#  \return A random real number in the desired range
def rand_real(range_start, range_stop):
    """
    Generates a random real number by using random.random to generate a random
    x in [0,1) and transforming it to be y between range_start and range_stop: 
        y = x * (range_stop - range_start) + range_start
    """
    return random.random()*(range_stop - range_start) + range_start

## \brief Generate a random rational number in [range_start,range_stop)
#  \param range_start The smaller end of the rational number range
#  \param range_stop The larger end of the rational number range
#  \return A random rational number in the desired range
def rand_mprational(range_start, range_stop):
    """
    Generates a random rational number by using rand_real and then making
    sure the denominator is 2^10.
    """
    return int((2**10)*rand_real(range_start, range_stop)) / (2**10)

## \brief Generate a random complex number with each of real and imaginary parts in [range_start,range_stop)
#  \param range_start The smaller end of the real number range
#  \param range_stop The larger end of the real number range
#  \return A random complex number in the desired range
def rand_complex(range_start, range_stop):
    """
    Generates a random complex number of the form a+bj where a,b are
    real numbers in [range_start, range_stop). 
    """
    return rand_real(range_start, range_stop) + 1j * rand_real(range_start, range_stop)

## \brief Generate a random complex rational number with each of real and imaginary parts in [range_start,range_stop)
#  \param range_start The smaller end of the rational number range
#  \param range_stop The larger end of the rational number range
#  \return A random rational number in the desired range
def rand_mpratcomplex(range_start, range_stop):
    """
    Generates a random mpratcomplex number of the form a+bj where a,b are
    real numbers in [range_start, range_stop). 
    """
    return rand_mprational(range_start, range_stop) + 1j * rand_mprational(range_start, range_stop)

## \brief Generate a single random number suitable for the matrix_random_entries function
#  \param scalar_type One of the scalarTypes supported by LARC
#  \param range_start The smallest possible return value
#  \param range_stop The top of the range of return values
#  \return A random number
def matrix_random_entry(scalar_type, range_start, range_stop):
    """
    Generates a single matrix entry for matrix_random_entries 
    that is determined by the scalarType.
    """
    # COMPLEX CASE
    if scalar_type == "Complex":
        rand_val = rand_complex(range_start, range_stop)

    # REAL CASE
    elif scalar_type == "Real":
        rand_val = rand_real(range_start, range_stop)

    # INTEGER CASE
    elif scalar_type == "Integer":
        if type(range_start) != int or type(range_stop) != int:
            raise ValueError("range_start and range_stop must be integers when scalarType is INTEGER!")
        rand_val = random.randrange(range_start, range_stop)

    # BOOLEAN CASE
    elif scalar_type == "Boolean":
        rand_val = random.randrange(0, 1)

    # MPINTEGER CASE
    elif scalar_type == "MPInteger":
        if type(range_start) != int or type(range_stop) != int:
            raise ValueError("range_start and range_stop must be integers when scalarType is MPINTEGER!")
        rand_val = random.randrange(range_start, range_stop)

    # MPRATIONAL CASE
    elif scalar_type == "MPRational":
        rand_val = rand_mprational(range_start, range_stop)

    # MPRATCOMPLEX CASE
    elif scalar_type == "MPRatComplex":
        rand_val = rand_mpratcomplex(range_start, range_stop)

    # MPREAL CASE
    elif scalar_type == "MPReal":
        rand_val = rand_real(range_start, range_stop)

    # MPCOMPLEX CASE
    elif scalar_type == "MPComplex":
        rand_val = rand_complex(range_start, range_stop)

    # CLIFFORD CASE
    elif scalar_type == "Clifford":
        # Just do Real Clifford for now.
        # TBD: Figure out how to detect Complex Clifford from Python side.
        rand_val = rand_mprational(range_start, range_stop)

    # UPPER or LOWER CASE
    elif scalar_type in ("Upper", "Lower"):
        rand_val = rand_mprational(0.0, 1.0)

    # OTHER
    else:
        raise Exception("Do not know how to build matrix entry for type %s." % scalar_type)

    return rand_val

## \brief Does the entry selection for matrix_random_matrixID.
#  \param scalar_type A string (from LARC) determining how to convert python variables into understandable strings
#  \param numEntries The size of the returned list
#  \param range_start The smallest random number returnable
#  \param range_stop The largest random number returnable
#  \param sparsity Determines the number of zeros
#  \return A list of random values to pass to matrix_random_matrixID
def matrix_random_entries(scalar_type, numEntries, range_start, range_stop, sparsity=None):
    """
    Does the entry selection for matrix_random_matrixID.
    """
    #apply sparsity requirements if present
    if sparsity is not None:
        if sparsity > 1 or sparsity < 0:
            raise Exception("Sparsity {} is not a value in [0,1]".format(sparsity))
        #calculate required number of zeros from sparsity
        numZerosReq = int(round(sparsity * numEntries))
        #randomly choose zero locations
        zeroPos = random.sample(range(numEntries), numZerosReq)
        #generate random values in desired range, placing zeros
        randVals = []
        for i in range(numEntries):
            if i in zeroPos:
                randVals.append(0)
            else:
                val = 0
                count = 0
                while val == 0:
                    val = matrix_random_entry(scalar_type, range_start, range_stop)
                    count += 1
                    if count > 100:
                        raise Exception("Invalid range for desired sparsity - must allow nonzero values for sparsity != 1.")
                randVals.append(val)
    # otherwise, just generate random values in desired range
    else: 
        randVals = [matrix_random_entry(scalar_type, range_start, range_stop) for i in range(numEntries)]
    return randVals

## \brief creates a random matrix and stores it in the LARC MatrixStore
#  \param scalar_type A string (from LARC) determining how to convert python variables into understandable strings
#  \param row_level The log-based-2 of the number of rows of the matrix
#  \param col_level The log-based-2 of the number of columns of the matrix
#  \param range_start The smallest random number returnable
#  \param range_stop The largest random number returnable
#  \param sparsity Determines the number of zeros
#  \return The LARC MatrixID for the random matrix
def matrix_random_matrixID(scalar_type, row_level, col_level, range_start, range_stop, sparsity=None):
    """
    Generate a matrix of size m x n, where m = 2^row_level and n = 2^col_level, 
    with entries in [range_start, range_stop). If a sparsity value in [0,1] is
    given, rounds sparsity * m * n to the closest integer value and generates
    the matrix such that there are exactly that many (randomly placed) zeros
    in the matrix. 
    """
    #size of matrix
    n = (2**row_level)*(2**col_level)
    randVals = matrix_random_entries(scalar_type, n, range_start, range_stop, sparsity)
    # convert list to a matrix and add to store
    # python lists of strings are automatically converted by SWIG to char ** 
    randMat_ID = row_major_list_to_store(map_to_str(randVals, scalar_type), row_level, col_level, 2**col_level)

    #return matrix id
    return randMat_ID

## \brief creates a LARC matrix with every element equal
#  \param value_str A scalar value in string format
#  \param row_level The log-based-2 of the number of rows of the matrix
#  \param col_level The log-based-2 of the number of columns of the matrix
#  \return The LARC MatrixID for the constant matrix
def matrix_constant_entry_matrixID(value_str, row_level, col_level):
    """
    Generate a matrix of size m x n, where m = 2^row_level and n = 2^col_level,
    where each entry has the value specified by value_str.
    """
    MATRIX_ID_INVALID = -1
    if (row_level == 0) and (col_level == 0):
        # Just make and return a 1 x 1 matrix.
        return get_valID_from_valString(value_str)

    elif row_level == 0:
       # For row vectors, recurse on columns.
       sub_answer = matrix_constant_entry_matrixID(value_str, row_level, col_level - 1)
       panel = [sub_answer, sub_answer, MATRIX_ID_INVALID, MATRIX_ID_INVALID]
       return get_pID_from_four_sub_pIDs(*panel, row_level, col_level)

    elif col_level == 0:
       # For column vectors, recurse on rows.
       sub_answer = matrix_constant_entry_matrixID(value_str, row_level - 1, col_level)
       panel = [sub_answer, MATRIX_ID_INVALID, sub_answer, MATRIX_ID_INVALID]
       return get_pID_from_four_sub_pIDs(*panel, row_level, col_level)

    else:
       # Otherwise, recurse on both rows and columns.
       sub_answer = matrix_constant_entry_matrixID(value_str, row_level - 1, col_level - 1)
       panel = [sub_answer, sub_answer, sub_answer, sub_answer]
       return get_pID_from_four_sub_pIDs(*panel, row_level, col_level)


##############################################################
##   Redirect stdout from C at python level
##############################################################

#from contextlib import contextmanager

# Taken from stackoverflow question titled 'Redirect stdout to a file in Python?' and answer by 'jfs' on Mar 16 '14. 

#def fileno(file_or_fd):
#    fd = getattr(file_or_fd, 'fileno', lambda: file_or_fd)()
#    if not isinstance(fd, int):
#        raise ValueError("Expected a file ('.fileno()') or a file descriptor")
#    return fd
#
#@contextmanager
#def stdout_redirected(to=os.devnull, stdout=None):
#    if stdout is None:
#        stdout = sys.stdout
#
#    stdout_fd = fileno(stdout)
#    #
#    #
#    with os.fdopen(os.dup(stdout_fd), 'wb') as copied:
#        stdout.flush()
#        try:
#            os.dup2(fileno(to), stdout_fd)
#        except ValueError:
#            with open(to, 'wb') as to_file:
#                os.dup2(to_file.fileno(), stdout_fd)
#        try:
#            yield stdout
#        finally:
#            stdout.flush()
#            os.dup2(copied.fileno(), stdout_fd)

#EXAMPLE that captures stdout and writes to file. 
#can call without arguments if you're not interested
#in writing to file. 

#stdout_fd = sys.stdout.fileno()
#with open('output.txt', 'w') as f, stdout_redirected(f):
#    print('redirected to a file')
#    os.write(stdout_fd, b'it is redirected now\n')
#    os.system('echo this is also redirected')
#print('this goes back to stdout')

## \brief Returns the row and column levels from a matrix stored on disk in a LARCmatrix json file
#  \param matrix_path The path to the file
#  \return The maximum row level and the maximum column level of the matrix
def get_levels_from_larcMatrix_file(matrix_path):
    """
    Obtains the max row_level and col_level of a matrix by reading until the
    end of the table section and observing the final entry in the table.
    input:  matrix_path path to LARCMatrix json file
    output: min_max_level
    """

    end_key   = "end"       #json table entries are numbers and 'end'
    table_key = "table"     #json table is keyed with 'table'

    # open as json file and read in LARCMatrix data 
    with open(matrix_path, 'r') as f:
        lines = ['', '']
        passed_table = False
        while (end_key not in lines[1]) or (not passed_table):
            lines.append(f.readline())
            lines.pop(0)
            # set marker that we've entered the 'table' section
            # since other sections have 'end' markers.
            if table_key in lines[0]:
                passed_table = True

    # ex:
    # '    "221":[4, 2, 169, 191, 207, 220],'
    # next line obtains [4,2]
    levels = lines[0].split(':[')[1].split(',')[:2]
    return map(int, levels)


## \brief Returns the most likely scalarType of a matrix stored on disk in a LARCmatrix json file
#  \param matrix_path The path to the file
#  \return A predicted scalarType for the data in the file
def get_scalarTypeStr_from_larcMatrix_file(matrix_path):
    """
    Predict the scalarType of a matrix by counting the number of entries in 
    lines with scalars in the LARCMatrix json file. 
    """

    end_key   = "end"       #json table entries are numbers and 'end'
    table_key = "table"     #json table is keyed with 'table'

    # open as json file and read in LARCMatrix data 
    with open(matrix_path, 'r') as f:
        line = ""
        while table_key not in line:
            line = f.readline()
        # first entry in table is always a scalar

        line = f.readline()
        while end_key not in line:
            while len(line.rstrip(",\n").split(",")) == 6:
                line = f.readline()
            if end_key in line:
                break
            scalar = line.rstrip("],\n").split(",")[2]
            scalar = scalar.strip(' \"')
            if '+' in scalar or 'i' in scalar: 
                return "Complex"
            if '/' in scalar:
                return "MPRational"
            try:
                scalar = int(scalar)
            except ValueError as e:
                return "Real"
            if scalar >= 2**64: 
                return "MPInteger"
            line = f.readline()

    return "Integer"

##############################################################################

##
# \ingroup larc
# \class LarcSettings
# \brief The LarcSettings class provides one way to initialize LARC
#
class LarcSettings:
    
    def __init__(self):
        self.mat_hash_exponent = None
        self.ops_hash_exponent = None 
        self.max_mat_level = None
        self.round_bits = -1
        self.trunc_bits = -1
        self.screen_output_level =  1

    def addArgsToParser(self, parser):
        parser.add_argument("-l", "--level", type=int, help="set 'Maximum matrix level' for initializing LARC. otherwise, automatically detected. ")
        parser.add_argument("--mat_exp", type=int, help="set 'Matrix store hash exponent' for initializing LARC. Defaults to 22 (for levels <= 10) or 30 (for levels > 10).")
        parser.add_argument("--ops_exp", type=int, help="set 'Operations store hash exponent' for initializing LARC. defaults to 19 (for levels <= 10) or 31 (for levels > 10).")
        parser.add_argument("--round_bits", type=int, default=-1, help="set number of significant bits for rounding for initializing LARC. defaults to LARC defaults.")
        parser.add_argument("--trunc_bits", type=int, default=-1, help="set number of significant bits for trucation for initializing LARC. defaults to LARC defaults.")
        parser.add_argument("-i", "--reportInterval", type=int, default=0, help="set stdout matrix/op store report time interval in minutes")
        parser.add_argument("-s", "--screen_output_level", type=int, default=1, help="set whether initialization output is verbose or quiet")

    def parseArgs(self, args):
        if args.level is not None:
            self.max_mat_level = args.level

        self.round_bits = args.round_bits
        self.trunc_bits = args.trunc_bits
        
        if args.mat_exp is not None:
            self.mat_hash_exponent = args.mat_exp
        if args.ops_exp is not None:
            self.ops_hash_exponent = args.ops_exp
        if args.screen_output_level is not None:
            self.screen_output_level = args.screen_output_level
        

    def getLevelFromFiles(self, *file_paths):
        if self.max_mat_level is None:
            self.max_mat_level = max([max(get_levels_from_larcMatrix_file(json_file)) for json_file in file_paths])
            print("\nAutomatically detected minimum 'Maximum matrix level' for initialzing LARC to be {}.".format(self.max_mat_level))
        else:
            print("\nUser provided 'Maximum matrix level' for initializing LARC to be {}.".format(self.max_mat_level))

    def initLarcArgs(self):
        if self.max_mat_level is None:
            raise Exception("Maximum matrix level needed to initialize LARC.") 
        if self.mat_hash_exponent is None:
            if self.max_mat_level <= 10:
                self.mat_hash_exponent = 22
            else:
                self.mat_hash_exponent = 30
        if self.ops_hash_exponent is None:
            if self.max_mat_level <= 10:
                self.ops_hash_exponent = 19
            else:
                self.ops_hash_exponent = 31
        if self.screen_output_level is None:
            self.screen_output_level = 1

        return [self.mat_hash_exponent, self.ops_hash_exponent,
                self.max_mat_level, self.round_bits, self.trunc_bits,
                self.screen_output_level]

    def initLarc(self):
        # Note: if you made a build on top of LARC, then the LARC that this
        # routine initializes will probably be SEPARATE from the LARC build you
        # want to be using. In this case, it's still safe to use initLarcArgs. 
        args = self.initLarcArgs()
        # THIS USED TO BE larcSWIG.
        initialize_larc(*args)

    def check_scalarType_matches_files(self, *file_paths):
        # THIS USED TO BE larcSWIG.
        scalarTypeStr = cvar.scalarTypeDef

        for json_file in file_paths:
            jScalarType = get_scalarTypeStr_from_larcMatrix_file(json_file).lower()[0]
            if scalarTypeStr != jScalarType:
                print("WARNING: current larc make scalarType ({}) might not match LARCMatrix file scalarType ({}). file path is {}".format(
                          scalarTypeStr, jScalarType, json_file) )

##
# \ingroup inlarc
# \brief Adds the value of a python variable to the matrix store and puts a hold on it
# \param value A python variable
# \param scalar_type A string which determines the conversion of the input into a LARC-understood string
# \result The MatrixID for the stored scalar
#
# A function to allow a python user to store and hold scalar values. The
# intent is to allow the user to preload the matrix store with their own
# selected scalar values, but the function may be used at any time during
# a calculation.
#
def load_and_hold_scalar_to_matrixStore(value,scalar_type):
    """
    Adds a scalar value to the LARC matrix store (or finds its matrixID if
    it has already been stored), and places a hold on it so that LARC will
    not remove it until the hold is released. The LARC function which removes
    the hold (release_hold_matrix) can be called directly from
    python.
    """
    value_string = value_to_string(value,scalar_type)
    s_ID = get_valID_from_valString(value_string)
    set_hold_matrix(s_ID)
    return s_ID

##
# \ingroup larc
# \brief Takes a numpy.matrix and converts it to a list of strings, then adds it to the matrixStore
# \param pmat The numpy matrix to be put into the MatrixStore
# \result The matrixID of the new LARCmatrix
#
def add_numpy_matrix_to_matrix_store(pmat):
    """
    This utility is intended to make it easy for a python user
    to add a *small* dense matrix to the matrix store. For
    larger matrices, memory may be an issue...
    """

    # find dimensions of input numpy matrix
    # shape is an attribute of numpy.matrix, so an error
    # will be thrown if pmat is not a numpy matrix
    dimensions = list(pmat.shape)

    #print("pmat = ",pmat)
    #print("type(pmat) = ", type(pmat))
    #print("type(pmat[0]) = ", type(pmat[0]))
    #print("type(pmat[0][0]) = ", type(pmat[0][0]))

    # when matrix dimensions are not 2**j by 2**k, pad matrix by
    # adding the smallest number of zero rows to bottom and/or
    # zero columns to right to make them so
    [row_lev, col_lev] = map(int,numpy.ceil(numpy.log2(dimensions)))
    row_size = 2**row_lev
    col_size = 2**col_lev
    row_pad = row_size - dimensions[0]
    col_pad = col_size - dimensions[1]
    if ((row_pad,col_pad)==(0,0)):
       padded = pmat
       # turn the matrix into a 1-D list array
       # note that padded is type matrix, so flatten makes a list 
       #   containing one list of lots of numbers, so to get rid of the
       #   outer list we index the entry we want
       plist = padded.flatten().tolist()[0]
    else:
       print("in add_numpy_matrix_to_matrix_store:")
       print("input matrix padded from size ", dimensions, "to size ",
             [row_size,col_size])
       padded = numpy.pad(pmat, ((0,row_pad), (0,col_pad)), mode='constant')
       # turn the matrix into a 1-D list array
       # note that padded is type ndarray, so flatten makes a list
       #   containing lots of numbers
       plist = padded.flatten().tolist()

    # print("padded = ",padded)
    # print("type(padded) = ", type(padded))
    # print("type(padded[0]) = ", type(padded[0]))
    # print("type(padded[0][0]) = ", type(padded[0][0]))
    # print("plist = ", plist)
    # print("type(plist) = ", type(plist))
    # print("type(plist[0]) = ", type(plist[0]))

    # convert to string format appropriate to scalarType
    plist_strs = map_to_str(plist,cvar.scalarTypeStr)
    del plist

    # add matrix to matrixStore
    matID = row_major_list_to_store(plist_strs, row_lev,
                   col_lev, col_size)
    return matID
