#!/usr/bin/env python

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
sys.path.append(os.path.join(os.path.dirname(__file__),"../../src"))
import pylarc
import numpy as np
from ctypes import *

# Python Interface version of test_math.py
# Uses matrixIDs instead of addresses

if __name__ == '__main__':

    print "This code tests some basic matrix building and reading routines\n"

    # initialize larc
    mat_store_exp = 26
    op_store_exp = 24
    max_level = 10
    rnd_sig_bits = -1   # default value
    trunc_to_zero_bits = -1  # default value
    pylarc.create_report_thread(1800)
    pylarc.initialize_larc(mat_store_exp,op_store_exp,max_level,rnd_sig_bits,trunc_to_zero_bits)

    # Define string for using in formating filenames
    if pylarc.cvar.scalarTypeDef == 'i':
        scalarType = "Integer"
    elif pylarc.cvar.scalarTypeDef == 'c':
        scalarType = "Complex"
    elif pylarc.cvar.scalarTypeDef == 'r':
        scalarType = "Real"
    else:
        raise Exception('scalarTypeDef %s was not handled.'%(pylarc.cvar.scalarTypeDef,))

    #  PLAYING WITH PREWRITTEN NONSQUARE MATRIX
    filename = "../dat/in/sample.1.2.%s.json" %scalarType
    print "About to test read %s\n" %filename
    samp_matrixID = pylarc.matrix_read_json_file_matrixID(os.path.join(os.path.dirname(__file__),filename))
    print "We read in the json file\n"
    pylarc.print_matrix_naive_by_matrixID(samp_matrixID)
    
    print "does scalarM1_val print?"
    scalarM1_val = -1
    scalarM1_matrixID = pylarc.matrix_get_matrixID_from_scalar(scalarM1_val)
    pylarc.print_matrix_naive_by_matrixID(scalarM1_matrixID)
    
    print "testing scalar_mult:"
    samp2_matrixID = pylarc.scalar_mult_matrixID(scalarM1_matrixID,samp_matrixID)
    pylarc.print_matrix_naive_by_matrixID(samp2_matrixID)
    
    print "testing addition:"
    samp3_matrixID = pylarc.matrix_add_matrixID(samp_matrixID,samp2_matrixID)
    pylarc.print_matrix_naive_by_matrixID(samp3_matrixID)
    
    print "testing adjoint:"
    samp3_matrixID = pylarc.matrix_adjoint_matrixID(samp_matrixID)
    pylarc.print_matrix_naive_by_matrixID(samp3_matrixID)
    
    print "testing non-square matrix mult:"
    samp4_matrixID = pylarc.matrix_mult_matrixID(samp3_matrixID,samp_matrixID)
    pylarc.print_matrix_naive_by_matrixID(samp4_matrixID)
    print ""
    samp4_matrixID = pylarc.matrix_mult_matrixID(samp_matrixID,samp3_matrixID)
    pylarc.print_matrix_naive_by_matrixID(samp4_matrixID)
    print "testing kron product:"
    samp4_matrixID = pylarc.kronecker_product_matrixID(samp_matrixID,samp_matrixID)
    pylarc.print_matrix_naive_by_matrixID(samp4_matrixID)
    print "testing join:"
    samp4_matrixID = pylarc.join_matrixID(samp_matrixID,samp_matrixID)
    pylarc.print_matrix_naive_by_matrixID(samp4_matrixID)
    print "testing stack:"
    samp4_matrixID = pylarc.stack_matrixID(samp_matrixID,samp_matrixID)
    pylarc.print_matrix_naive_by_matrixID(samp4_matrixID)
