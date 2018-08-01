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

if 2 == len(sys.argv):
    sys.path.append(sys.argv[1])
    sys.argv = [sys.argv[0]]

sys.path.append(os.path.join(os.path.dirname(__file__),"../../src"))
from pylarc import *
import numpy as np
from ctypes import *
import unittest

# THIS IS AN EXAMPLE OF USING THE PYTHON UNIT TEST LIBRARY.
# DOCUMENTATION LIVES ON THE OUTSIDE WEB AT
# https://docs.python.org/2/library/unittest.html .

class TestLARCmath(unittest.TestCase):
    def testAdd(self):
        ### SCALAR ADDITION ###
        ### THIS FAILS IF WE ARE NOT IN REAL  ##
        # 0 + 0 = 0
        self.assertEqual(matrix_add_matrixID(cvar.matID_scalar0,cvar.matID_scalar0), cvar.matID_scalar0) 
        
        # 1 + 0 = 1
        self.assertEqual(matrix_add_matrixID(cvar.matID_scalar1, cvar.matID_scalar0), cvar.matID_scalar1)
        
        # 0 + 1 = 1
        self.assertEqual(matrix_add_matrixID(cvar.matID_scalar0, cvar.matID_scalar1), cvar.matID_scalar1)
        
        # 1 + 0 = 0 + 1
        self.assertEqual(matrix_add_matrixID(cvar.matID_scalar1, cvar.matID_scalar0), matrix_add_matrixID(cvar.matID_scalar0, cvar.matID_scalar1))
        
        # 1 + 1 = 2
        #constructing scalar 2 and adding it to the matrix store
# NOTE: the following has been replaced, partly because it's overly complicated,
# mostly because it assumes Complex type 
        # a = np.matrix([[2+0j]])
        # alist = a.reshape(-1).tolist()[0]
        # arr = buildComplexArray(alist)
        # level = 0
        # dim = 2**level
        # aser = row_major_list_to_store_matrixID(arr, level, level, dim)
        aser = matrix_get_matrixID_from_scalar(2)
        
        self.assertEqual(matrix_add_matrixID(cvar.matID_scalar1, cvar.matID_scalar1), aser)
        
        ### SQUARE MATRIX ADDITION ###
        # check that error is thrown when dimensions do not match;
        # (would need to implement Python functionality to throw exception)
        # e.g. self.assertRaises(exception, callable, args, keywords)
        
        # 4x4 identity matrix + 4x4 zero matrix = 4x4 identity matrix
        self.assertEqual(matrix_add_matrixID(I2_matrixID, Z2_matrixID), I2_matrixID)
        
        
        ### NON-SQUARE MATRIX ADDITION ###
        # check that error is thrown when dimensions do not match


    def testScalarMult(self):
		# test scalar multiplication
		# 0 * 0
		self.assertEqual(scalar_mult_matrixID(cvar.matID_scalar0,cvar.matID_scalar0), cvar.matID_scalar0) 
		
		# 1 * 0 = 0
		self.assertEqual(scalar_mult_matrixID(cvar.matID_scalar1, cvar.matID_scalar0), cvar.matID_scalar0)
		
		# 1 * 1 = 1
		self.assertEqual(scalar_mult_matrixID(cvar.matID_scalar1, cvar.matID_scalar1), cvar.matID_scalar1)
        
		# 0 * 4x4 identity matrix = 4x4 0 matrix
		self.assertEqual(scalar_mult_matrixID(cvar.matID_scalar0,I2_matrixID), Z2_matrixID)
		
		# 1 * 4x4 identity matrix = 4x4 identity matrix
		self.assertEqual(scalar_mult_matrixID(cvar.matID_scalar1,I2_matrixID), I2_matrixID)
		
    def testMult(self):
		# test matrix multiplication
		# 4x4 identity matrix * 4x4 identity matrix = 4x4 identity matrix
		self.assertEqual(matrix_mult_matrixID(I2_matrixID,I2_matrixID), I2_matrixID) 	
		
		
#    def testKronecker(self)
#    def testAdjoint(self)
    
        

if __name__ == '__main__':

    print "Running LARC math unit tests..."
    initialize_larc(26,24,10,-1,-1)
    create_report_thread(1800)
    print "matrixID for 0 ", cvar.matID_scalar0
    print "matrixID for 1 ", cvar.matID_scalar1
    
        # Define string for using in formating filenames
    if cvar.scalarTypeDef == 'i':
        scalarType = "Integer"
    elif cvar.scalarTypeDef == 'c':
        scalarType = "Complex"
    elif cvar.scalarTypeDef == 'r':
        scalarType = "Real"
    else:
        raise Exception('scalarTypeDef %s was not handled.'%(pylarc.cvar.scalarTypeDef,))
        
    # build Zero matrices
    level = 2
    dim_whole = 2**level
    Z2_arr = buildArray([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    Z2_matrixID = row_major_list_to_store_matrixID(Z2_arr, level, level, dim_whole)
    print_matrix_naive_by_matrixID(Z2_matrixID)
    print "\nmatrixID for  the level 2 zero matrix\n", Z2_matrixID

    # build Identity matrices
    I2_arr = buildArray([1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1])
    I2_matrixID = row_major_list_to_store_matrixID(I2_arr, level, level, dim_whole)
    print_matrix_naive_by_matrixID(I2_matrixID)
    print "\nmatrixID for  the level 2 identity matrix\n", I2_matrixID
    
    unittest.main()

