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
import math

if __name__ == '__main__':

    print "This code tests matrixID interface\n"

    # initialize larc
    mat_store_exp = 26
    op_store_exp = 24
    max_level = 3
    rnd_sig_bits = 15
    trunc_to_zero_bits = 20
    pylarc.initialize_larc(mat_store_exp,op_store_exp,max_level,rnd_sig_bits,trunc_to_zero_bits)
    
    pylarc.create_report_thread(1800)

    # Define string for using in formating filenames
    if pylarc.cvar.scalarTypeDef == 'i':
        raise Exception('The globbing test does not work for scalarType = Integer')
    elif pylarc.cvar.scalarTypeDef == 'c':
        scalarType = "Complex"
    elif pylarc.cvar.scalarTypeDef == 'r':
        scalarType = "Real"
    else:
        raise Exception('scalarTypeDef %s was not handled.'%(pylarc.cvar.scalarTypeDef,))

    # MAKE a matrix to test Globbing
    # two by two tests
    level = 1
    dim_whole = 2**level

    arr_a = pylarc.buildArray([.4, 0, 0, .3])
    a_mID = pylarc.row_major_list_to_store_matrixID(arr_a,level,level,dim_whole)
    arr_b = pylarc.buildArray([.8, 0, 0, .6])
    b_mID = pylarc.row_major_list_to_store_matrixID(arr_b,level,level,dim_whole)
    arr_c = pylarc.buildArray([-.4, 0, 0, -.3])
    c_mID = pylarc.row_major_list_to_store_matrixID(arr_c,level,level,dim_whole)
    print "The matrixIDs of a, b, and c are %d %d %d\n" %(a_mID,b_mID,c_mID)
    
    d_mID = pylarc.matrix_add_matrixID(a_mID,a_mID)
    
    print "matrix a:"
    pylarc.print_matrix_naive_by_matrixID(a_mID)
     
    print "\nMatrix id of a + a: %d, should be that of b: %d \n" %(d_mID,b_mID)
    if (d_mID == b_mID): 
		print "a + a PASSED: \n" 
		pylarc.print_matrix_naive_by_matrixID(d_mID)
    else: print "FAILED: [.4,0,0,.3] + [.4,0,0,.3] = [.8,0,0,.6]\n"

    arr_prod1 = pylarc.buildArray([.16, 0, 0, .09])
    prod1_mID = pylarc.row_major_list_to_store_matrixID(arr_prod1,level,level,dim_whole)
    arr_e = pylarc.buildArray([1, -1, -1, 1])
    e_mID = pylarc.row_major_list_to_store_matrixID(arr_e,level,level,dim_whole)
    print "\nmatrix e:"
    pylarc.print_matrix_naive_by_matrixID(e_mID)
    
    arr_prod2 = pylarc.buildArray([.4, -.4, -.3, .3])
    prod2_mID = pylarc.row_major_list_to_store_matrixID(arr_prod2,level,level,dim_whole)
    print "\nThe matrixIDs of prod1, e, and prod2 are %d %d %d\n" %(prod1_mID,e_mID,prod2_mID)

    m_mID = pylarc.matrix_mult_matrixID(a_mID,a_mID)
    n_mID = pylarc.matrix_mult_matrixID(a_mID,e_mID)
    print "The matrixIDs of m and n are %d %d\n" %(m_mID,n_mID)
    print "Matrix id of a * a: %d, should be that of prod1: %d \n" %(m_mID,prod1_mID)
    if (m_mID == prod1_mID): 
		print "a * a PASSED:"
		pylarc.print_matrix_naive_by_matrixID(m_mID)
    else: print "FAILED: [.4,0,0,.3] * [.4,0,0,.3] = [.16, 0, 0, .09]\n"

    print "\nMatrix id of a * e: %d, should be that of prod2: %d \n" %(n_mID,prod2_mID)
    if (n_mID == prod2_mID): 
		print "a * e PASSED:"
		pylarc.print_matrix_naive_by_matrixID(n_mID)
    else: print "FAILED: [.4,0,0,.3] * [1, -1, -1, 1] = [.4, -.4, -.3, .3]\n"

    # scalar tests
    level = 0
    dim_whole = 2**level

    arr_f = pylarc.buildArray([.4])
    f_mID = pylarc.row_major_list_to_store_matrixID(arr_f,level,level,dim_whole)
    print "We input the value into f (%d) of .4, the matrix is:" %f_mID
    pylarc.print_matrix_naive_by_matrixID(f_mID)

    arr_g = pylarc.buildArray([.4+1e-5])
    g_mID = pylarc.row_major_list_to_store_matrixID(arr_g,level,level,dim_whole)
    print "We input the value into g (%d) of .4+1e-5, the matrix is:" %g_mID
    pylarc.print_matrix_naive_by_matrixID(g_mID)

    arr_h = pylarc.buildArray([.4-1e-10])
    h_mID = pylarc.row_major_list_to_store_matrixID(arr_h,level,level,dim_whole)
    print "We input the value into h (%d) of .4-1e-10, the matrix is:" %h_mID
    pylarc.print_matrix_naive_by_matrixID(h_mID)

# remaining tests in test_globbing.py rely on index_* values that are pointers
# to scalars (defined in global.c). To get these running, we need to create new
# variables for matrixIDs and make them global.
