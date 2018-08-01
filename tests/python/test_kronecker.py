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

if __name__ == '__main__':
	# This version references matrices by matrixID instead of pointers
	
    print "This code tests the kronecker product routine"

    # initialize larc
    mat_store_exp = 26
    op_store_exp = 24
    max_level = 10
    rnd_sig_bits = -1   # default value
    trunc_to_zero_bits = -1  # default value
    pylarc.create_report_thread(1800)
    pylarc.initialize_larc(mat_store_exp,op_store_exp,max_level,rnd_sig_bits,trunc_to_zero_bits)


    # In the Makefile you can compile with "-DTYPE=REAL, -DTYPE=COMPLEX, -DTYPE=INTEGER"
    # Define string for using in formating filenames
    if pylarc.cvar.scalarTypeDef == 'i':
        scalarType = "Integer"
    elif pylarc.cvar.scalarTypeDef == 'c':
        scalarType = "Complex"
    elif pylarc.cvar.scalarTypeDef == 'r':
        scalarType = "Real"
    else:
        raise Exception('scalarTypeDef %s was not handled.'%(pylarc.cvar.scalarTypeDef,))

    # build array in C from Python list of scalars
    print "Using row_major_list_to_store on data entered from python\n"

    # create a matrix in python
    if pylarc.cvar.scalarTypeDef == 'i':
        a = np.matrix([[1, 2, 3, 4],[5, 6, 7, 8]])
        b = np.matrix([[9, 10, 11, 12]])
#        a = np.matrix([[1, 3, 5, 6],
#                       [8, 6, 3, 1],
#                       [-9, 11, 13, 15],
#                       [16, 13, 12, 10]])
    elif pylarc.cvar.scalarTypeDef == 'c':
        a = np.matrix([[1+2j, 2+3j, 3+4j, 4+5j],[5+1j, 6+2j, 7+3j, 8+4j]])
        b = np.matrix([[9-1j, 10-2j, 11-3j, 12-4j]])
#        a = np.matrix([[1+2j, 3+4j, 5+6j, 7+8j],
#                       [8+7j, 6+5j, 3+4j, 1+2j],
#                       [9+10j, 11+12j, 13+14j, 15+16j],
#                       [16+15j, 14+13j, 12+11j, 10+9j]])
    elif pylarc.cvar.scalarTypeDef == 'r':
        a = np.matrix([[1, 2, 3, 4],[5, 6, 7, 8]])
        b = np.matrix([[9, 10, 11, 12]])
#        a = np.matrix([[1, 3, .5, 6],
#                       [8, 6, 3, .1],
#                       [-9, 11, 13, 1.5],
#                       [16, 13, 12, 10]])
    else:
        raise Exception('Do not know how to build matrix for type %s.'%(pylarc.cvar.scalarTypeDef,))

    print "COLUMN_COLUMN CASE"
    # turn the matrix into an array by reading off each row in turn (row major format)
    alist = a.reshape(-1).tolist()[0]
    aarr = pylarc.buildArray(alist)
    print 'array A:', pylarc.str_scalarTypeArray(aarr, len(alist))

    # creating or finding the matrix associated with the array
    serial_a = pylarc.row_major_list_to_store_matrixID(aarr, 3, 0, 1)
    pylarc.print_matrix_naive_by_matrixID(serial_a)
    print "\n"

    # turn the matrix into an array by reading off each row in turn (row major format)
    blist = b.reshape(-1).tolist()[0]
    barr = pylarc.buildArray(blist)
    print 'array B:', pylarc.str_scalarTypeArray(barr, len(blist))

    # creating or finding the matrix associated with the array
    serial_b = pylarc.row_major_list_to_store_matrixID(barr, 2, 0, 1)
    pylarc.print_matrix_naive_by_matrixID(serial_b)
    print "\n"

    print "kronecker product A \otimes B:"
    serial_c = pylarc.kronecker_product_matrixID(serial_a,serial_b)
    pylarc.print_matrix_naive_by_matrixID(serial_c)
    print "\n"

    print "kronecker product B \otimes A:"
    serial_c = pylarc.kronecker_product_matrixID(serial_b,serial_a)
    pylarc.print_matrix_naive_by_matrixID(serial_c)
    print "\n"

    print "ROW_COLUMN CASE"
    # turn the matrix into an array by reading off each row in turn (row major format)
    alist = a.reshape(-1).tolist()[0]
    aarr = pylarc.buildArray(alist)
    print 'array A:', pylarc.str_scalarTypeArray(aarr, len(alist))

    # creating or finding the matrix associated with the array
    serial_a = pylarc.row_major_list_to_store_matrixID(aarr, 0, 3, 8)
    pylarc.print_matrix_naive_by_matrixID(serial_a)
    print "\n"

    # turn the matrix into an array by reading off each row in turn (row major format)
    blist = b.reshape(-1).tolist()[0]
    barr = pylarc.buildArray(blist)
    print 'array B:', pylarc.str_scalarTypeArray(barr, len(blist))

    # creating or finding the matrix associated with the array
    serial_b = pylarc.row_major_list_to_store_matrixID(barr, 2, 0, 1)
    pylarc.print_matrix_naive_by_matrixID(serial_b)
    print "\n"

    print "kronecker product A \otimes B:"
    serial_c = pylarc.kronecker_product_matrixID(serial_a,serial_b)
    pylarc.print_matrix_naive_by_matrixID(serial_c)
    print "\n"

    print "COLUMN_ROW CASE"
    # turn the matrix into an array by reading off each row in turn (row major format)
    alist = a.reshape(-1).tolist()[0]
    aarr = pylarc.buildArray(alist)
    print 'array A:', pylarc.str_scalarTypeArray(aarr, len(alist))

    # creating or finding the matrix associated with the array
    serial_a = pylarc.row_major_list_to_store_matrixID(aarr, 3, 0, 1)
    pylarc.print_matrix_naive_by_matrixID(serial_a)
    print "\n"

    # turn the matrix into an array by reading off each row in turn (row major format)
    blist = b.reshape(-1).tolist()[0]
    barr = pylarc.buildArray(blist)
    print 'array B:', pylarc.str_scalarTypeArray(barr, len(blist))

    # creating or finding the matrix associated with the array
    serial_b = pylarc.row_major_list_to_store_matrixID(barr, 0, 2, 4)
    pylarc.print_matrix_naive_by_matrixID(serial_b)
    print "\n"

    print "kronecker product A \otimes B:"
    serial_c = pylarc.kronecker_product_matrixID(serial_a,serial_b)
    pylarc.print_matrix_naive_by_matrixID(serial_c)
    print "\n"

    print "ROW_ROW CASE"
    # turn the matrix into an array by reading off each row in turn (row major format)
    alist = a.reshape(-1).tolist()[0]
    aarr = pylarc.buildArray(alist)
    print 'array A:', pylarc.str_scalarTypeArray(aarr, len(alist))

    # creating or finding the matrix associated with the array
    serial_a = pylarc.row_major_list_to_store_matrixID(aarr, 0, 3, 8)
    pylarc.print_matrix_naive_by_matrixID(serial_a)
    print "\n"

    # turn the matrix into an array by reading off each row in turn (row major format)
    blist = b.reshape(-1).tolist()[0]
    barr = pylarc.buildArray(blist)
    print 'array B:', pylarc.str_scalarTypeArray(barr, len(blist))

    # creating or finding the matrix associated with the array
    serial_b = pylarc.row_major_list_to_store_matrixID(barr, 0, 2, 4)
    pylarc.print_matrix_naive_by_matrixID(serial_b)
    print "\n"

    print "kronecker product A \otimes B:"
    serial_c = pylarc.kronecker_product_matrixID(serial_a,serial_b)
    pylarc.print_matrix_naive_by_matrixID(serial_c)
    print "\n"

    print "kronecker product B \otimes A:"
    serial_c = pylarc.kronecker_product_matrixID(serial_b,serial_a)
    pylarc.print_matrix_naive_by_matrixID(serial_c)
    print "\n"

    print "MATRIX_ROW CASE"
    # turn the matrix into an array by reading off each row in turn (row major format)
    alist = a.reshape(-1).tolist()[0]
    aarr = pylarc.buildArray(alist)
    print 'array A:', pylarc.str_scalarTypeArray(aarr, len(alist))

    # creating or finding the matrix associated with the array
    serial_a = pylarc.row_major_list_to_store_matrixID(aarr, 1, 2, 4)
    pylarc.print_matrix_naive_by_matrixID(serial_a)
    print "\n"

    # turn the matrix into an array by reading off each row in turn (row major format)
    blist = b.reshape(-1).tolist()[0]
    barr = pylarc.buildArray(blist)
    print 'array B:', pylarc.str_scalarTypeArray(barr, len(blist))

    # creating or finding the matrix associated with the array
    serial_b = pylarc.row_major_list_to_store_matrixID(barr, 0, 2, 4)
    pylarc.print_matrix_naive_by_matrixID(serial_b)
    print "\n"

    print "kronecker product A \otimes B:"
    serial_c = pylarc.kronecker_product_matrixID(serial_a,serial_b)
    pylarc.print_matrix_naive_by_matrixID(serial_c)
    print "\n"

    print "MATRIX_COLUMN CASE"
    # turn the matrix into an array by reading off each row in turn (row major format)
    alist = a.reshape(-1).tolist()[0]
    aarr = pylarc.buildArray(alist)
    print 'array A:', pylarc.str_scalarTypeArray(aarr, len(alist))

    # creating or finding the matrix associated with the array
    serial_a = pylarc.row_major_list_to_store_matrixID(aarr, 1, 2, 4)
    pylarc.print_matrix_naive_by_matrixID(serial_a)
    print "\n"

    # turn the matrix into an array by reading off each row in turn (row major format)
    blist = b.reshape(-1).tolist()[0]
    barr = pylarc.buildArray(blist)
    print 'array B:', pylarc.str_scalarTypeArray(barr, len(blist))

    # creating or finding the matrix associated with the array
    serial_b = pylarc.row_major_list_to_store_matrixID(barr, 2, 0, 1)
    pylarc.print_matrix_naive_by_matrixID(serial_b)
    print "\n"

    print "kronecker product A \otimes B:"
    serial_c = pylarc.kronecker_product_matrixID(serial_a,serial_b)
    pylarc.print_matrix_naive_by_matrixID(serial_c)
    print "\n"

    print "ROW_MATRIX CASE"
    # turn the matrix into an array by reading off each row in turn (row major format)
    alist = a.reshape(-1).tolist()[0]
    aarr = pylarc.buildArray(alist)
    print 'array A:', pylarc.str_scalarTypeArray(aarr, len(alist))

    # creating or finding the matrix associated with the array
    serial_a = pylarc.row_major_list_to_store_matrixID(aarr, 1, 2, 4)
    pylarc.print_matrix_naive_by_matrixID(serial_a)
    print "\n"

    # turn the matrix into an array by reading off each row in turn (row major format)
    blist = b.reshape(-1).tolist()[0]
    barr = pylarc.buildArray(blist)
    print 'array B:', pylarc.str_scalarTypeArray(barr, len(blist))

    # creating or finding the matrix associated with the array
    serial_b = pylarc.row_major_list_to_store_matrixID(barr, 0, 2, 4)
    pylarc.print_matrix_naive_by_matrixID(serial_b)
    print "\n"

    print "kronecker product B \otimes A:"
    serial_c = pylarc.kronecker_product_matrixID(serial_b,serial_a)
    pylarc.print_matrix_naive_by_matrixID(serial_c)
    print "\n"

    print "COLUMN_MATRIX CASE"
    # turn the matrix into an array by reading off each row in turn (row major format)
    alist = a.reshape(-1).tolist()[0]
    aarr = pylarc.buildArray(alist)
    print 'array A:', pylarc.str_scalarTypeArray(aarr, len(alist))

    # creating or finding the matrix associated with the array
    serial_a = pylarc.row_major_list_to_store_matrixID(aarr, 1, 2, 4)
    pylarc.print_matrix_naive_by_matrixID(serial_a)
    print "\n"

    # turn the matrix into an array by reading off each row in turn (row major format)
    blist = b.reshape(-1).tolist()[0]
    barr = pylarc.buildArray(blist)
    print 'array B:', pylarc.str_scalarTypeArray(barr, len(blist))

    # creating or finding the matrix associated with the array
    serial_b = pylarc.row_major_list_to_store_matrixID(barr, 2, 0, 1)
    pylarc.print_matrix_naive_by_matrixID(serial_b)
    print "\n"

    print "kronecker product B \otimes A:"
    serial_c = pylarc.kronecker_product_matrixID(serial_b,serial_a)
    pylarc.print_matrix_naive_by_matrixID(serial_c)
    print "\n"
