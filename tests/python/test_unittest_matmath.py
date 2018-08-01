
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

# NOTE: to run on command line
#   python test_unittest_matmath.py -v

import os 
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../../src"))
from pylarc import *
import numpy as np
from ctypes import *

import matrix_quick_build
import random
import unittest

#class larcSetup():
#    """
#    In case we don't want to run initialize_larc for every test, 
#    we build a class to save initialized values one time.
#    """
#    def __init__(self):
#        initialize_larc(26,24,10,-1,-1)
## If we decide that we don't even want to set up the matrices every time
## a test is run, we can put all that in this larcSetup function. 
## then put the next two lines in the setUp() function for the test. 
## then to access a value in a test, you can do something like
##      self.init.zeroMatID
# to access persistent things to test from init. This is like unittest.setUpClass
# but I don't think we have that in this python version? 
#self.__class__.init = larcSetup()

class TestMatrixMaxnorm(unittest.TestCase):
    # Define a class variable that determines if larc initialization has been run. 
    ClassIsSetup = False

    def setUp(self):
        #NOTE: at this point, initialize_larc is run ones per TestCase subclass. 
        # If larc initialization has not been run yet, do it.
        if not self.ClassIsSetup:
            initialize_larc(26,24,10,-1,-1)
            self.__class__.ClassIsSetup = True
        # Define string for using in formating filenames
        # (cvar is in pylarc)
        if cvar.scalarTypeDef == 'i':
            self.scalarType = "Integer"
        elif cvar.scalarTypeDef == 'c':
            self.scalarType = "Complex"
        elif cvar.scalarTypeDef == 'r':
            self.scalarType = "Real"
        else:
            raise Exception('scalarTypeDef %s was not handled.'%(cvar.scalarTypeDef,))
        # set testing levels
        self.row_level = 3
        self.col_level = 3
        self.val_range= [-100, 100]
        self.sparsity = 0.3
        # build Zero matrix
        self.zeroMatID = matrix_quick_build.matrix_zero_matrixID(self.row_level, self.col_level)
        # generate random matrix
        self.randMatAID = matrix_quick_build.matrix_random_matrixID(self.scalarType, self.row_level, self.col_level, self.val_range[0], self.val_range[1], self.sparsity)
        # build -1
        self.neg1 = matrix_get_matrixID_from_scalar(-1)

    def test_maxnorm_zero(self):
        """
        maxnorm(0) = 0
        """
        self.assertEqual(matrix_maxnorm_matrixID(self.zeroMatID), 0)

    def test_maxnorm_id(self):
        """
        maxnorm(id) = id
        """
        idMatID = matrix_quick_build.matrix_identity_matrixID(self.row_level)
        self.assertEqual(matrix_maxnorm_matrixID(idMatID), 1)

    def test_maxnorm_negated(self):
        """
        maxnorm(-A) = maxnorm(A)
        """
        negAMatID = scalar_mult_matrixID(self.neg1, self.randMatAID)
        self.assertEqual(matrix_maxnorm_matrixID(self.randMatAID), matrix_maxnorm_matrixID(negAMatID))

    def test_maxnorm_scalar(self):
        """
        maxnorm(scalar) = abs(scalar)
        """
        scalMatID = matrix_quick_build.matrix_random_matrixID(self.scalarType, 0, 0, self.val_range[0], self.val_range[1])
        norm = abs(matrix_trace_matrixID(scalMatID))
        self.assertEqual(matrix_maxnorm_matrixID(scalMatID), norm)


class TestMatrixSparsity(unittest.TestCase):
    # Define a class variable that determines if larc initialization has been run. 
    ClassIsSetup = False

    def setUp(self):
        #NOTE: at this point, initialize_larc is run ones per TestCase subclass. 
        # If larc initialization has not been run yet, do it.
        if not self.ClassIsSetup:
            print "\n"
            initialize_larc(26,24,10,-1,-1)
            self.__class__.ClassIsSetup = True
        # Define string for using in formating filenames
        # (cvar is in pylarc)
        if cvar.scalarTypeDef == 'i':
            self.scalarType = "Integer"
        elif cvar.scalarTypeDef == 'c':
            self.scalarType = "Complex"
        elif cvar.scalarTypeDef == 'r':
            self.scalarType = "Real"
        else:
            raise Exception('scalarTypeDef %s was not handled.'%(cvar.scalarTypeDef,))
        # set testing levels
        self.row_level = 3
        self.col_level = 3
        self.val_range= [-100, 100]

    def test_sparse_val(self):
        """
        Verify sparsity(A) = sparsity value for generating random matrix A. 
        """
        sparsity = 0
        n = (2**self.row_level)*(2**self.col_level)
        while sparsity == 0:
            sparsity = random.random()
        randMatID = matrix_quick_build.matrix_random_matrixID(self.scalarType, self.row_level, self.col_level, self.val_range[0], self.val_range[1], sparsity)
        # convert to number of zeros
        # requested sparsity may be inaccurate by less than one zero 
        self.assertEqual(int(round(sparsity * n)), int(round(matrix_sparsity_matrixID(randMatID)*n)))

    def test_sparse_1(self):
        """
        Matrix generated with sparcity value 0 has sparsity 0.
        """
        onesMatID = matrix_quick_build.matrix_random_matrixID(self.scalarType, self.row_level, self.col_level, 1, 2, 0)
        self.assertEqual(matrix_sparsity_matrixID(onesMatID), 0)

    def test_sparse_id(self):
        """
        Identity matrix of order n has sparsity (n-1)/n.
        """
        idMatID = matrix_quick_build.matrix_identity_matrixID(self.row_level)
        self.assertEqual(matrix_sparsity_matrixID(idMatID), (2**self.row_level-1.0)/(2**self.row_level))

    def test_sparse_0(self):
        """
        Zero matrix has sparsity 0.
        """
        zeroMatID = matrix_quick_build.matrix_zero_matrixID(self.row_level, self.col_level)
        self.assertEqual(matrix_sparsity_matrixID(zeroMatID), 1)

class TestMatrixSum(unittest.TestCase):
    # Define a class variable that determines if larc initialization has been run. 
    ClassIsSetup = False

    def setUp(self):
        #NOTE: at this point, initialize_larc is run ones per TestCase subclass. 
        # If larc initialization has not been run yet, do it.
        if not self.ClassIsSetup:
            print "\n"
            initialize_larc(26,24,10,-1,-1)
            self.__class__.ClassIsSetup = True
        # Define string for using in formating filenames
        # (cvar is in pylarc)
        if cvar.scalarTypeDef == 'i':
            self.scalarType = "Integer"
        elif cvar.scalarTypeDef == 'c':
            self.scalarType = "Complex"
        elif cvar.scalarTypeDef == 'r':
            self.scalarType = "Real"
        else:
            raise Exception('scalarTypeDef %s was not handled.'%(cvar.scalarTypeDef,))
        # set testing levels
        self.row_level = 3
        self.col_level = 3
        self.val_range= [-100, 100]
        self.sparsity = .3
        # build Zero matrix
        self.zeroMatID = matrix_quick_build.matrix_zero_matrixID(self.row_level, self.col_level)
        # generate two random matrices of same size
        self.randMatAID = matrix_quick_build.matrix_random_matrixID(self.scalarType, self.row_level, self.col_level, self.val_range[0], self.val_range[1], self.sparsity)
        self.randMatBID = matrix_quick_build.matrix_random_matrixID(self.scalarType, self.row_level, self.col_level, self.val_range[0], self.val_range[1], self.sparsity)
        # build -1
        self.neg1 = matrix_get_matrixID_from_scalar(-1)

    def test_add_additive_inverse(self):
        """
        A + (-A) = 0
        """
        a = scalar_mult_matrixID(self.neg1, self.randMatAID)
        self.assertEqual(matrix_add_matrixID(self.randMatAID, a), self.zeroMatID)

    def test_add_add_zero(self):
        """
        A + 0 = A
        """
        self.assertEqual(matrix_add_matrixID(self.randMatAID, self.zeroMatID), self.randMatAID)

    def test_add_double(self):
        """
        A + A = 2A
        """
        scalMatID = matrix_get_matrixID_from_scalar(2)
        prodMatID = scalar_mult_matrixID(scalMatID, self.randMatAID)
        self.assertEqual(matrix_add_matrixID(self.randMatAID, self.randMatAID), prodMatID)

    def test_add_commutativity(self):
        """
        A + B = B + A
        """
        a = matrix_add_matrixID(self.randMatAID, self.randMatBID)
        b = matrix_add_matrixID(self.randMatBID, self.randMatAID)
        self.assertEqual(a,b)

    def test_add_commutativity_op_shortcut(self):
        """
        If A+B in op store, calculating B+A shouldn't require additional calculations.
        """
        #NOTE: this requires a reset of the store to test properly
        # if we decide to run setUp everywhere with tearDowns, remove next line. 
        self.__class__.ClassIsSetup = False
        self.setUp()
        a = matrix_add_matrixID(self.randMatAID, self.randMatBID)
        num_matrices_made = num_matrices_created()
        b = matrix_add_matrixID(self.randMatBID, self.randMatAID)
        num_matrices_made2 = num_matrices_created()
        self.assertEqual(num_matrices_made, num_matrices_made2)




class TestMatrixDifference(unittest.TestCase):
    # Define a class variable that determines if larc initialization has been run. 
    ClassIsSetup = False

    def setUp(self):
        #NOTE: at this point, initialize_larc is run ones per TestCase subclass. 
        # If larc initialization has not been run yet, do it.
        if not self.ClassIsSetup:
            print "\n"
            initialize_larc(26,24,10,-1,-1)
            self.__class__.ClassIsSetup = True
        # Define string for using in formating filenames
        # (cvar is in pylarc)
        if cvar.scalarTypeDef == 'i':
            self.scalarType = "Integer"
        elif cvar.scalarTypeDef == 'c':
            self.scalarType = "Complex"
        elif cvar.scalarTypeDef == 'r':
            self.scalarType = "Real"
        else:
            raise Exception('scalarTypeDef %s was not handled.'%(cvar.scalarTypeDef,))
        # set testing levels
        self.row_level = 3
        self.col_level = 3
        self.val_range= [-100, 100]
        self.sparsity = .3
        # build Zero matrix
        self.zeroMatID = matrix_quick_build.matrix_zero_matrixID(self.row_level, self.col_level)
        # generate two random matrices of same size
        self.randMatAID = matrix_quick_build.matrix_random_matrixID(self.scalarType, self.row_level, self.col_level, self.val_range[0], self.val_range[1], self.sparsity)
        self.randMatBID = matrix_quick_build.matrix_random_matrixID(self.scalarType, self.row_level, self.col_level, self.val_range[0], self.val_range[1], self.sparsity)
        # build -1
        self.neg1 = matrix_get_matrixID_from_scalar(-1)

    def test_diff_additive_inverse(self):
        """
        A - A = 0
        """
        self.assertEqual(self.zeroMatID, matrix_diff_matrixID(self.randMatAID, self.randMatAID))

    def test_diff_subtract_zero(self):
        """
        A - 0 = A
        """
        self.assertEqual(self.randMatAID, matrix_diff_matrixID(self.randMatAID, self.zeroMatID))

    def test_diff_uni_invert(self):
        """
        0 - A = -A
        """
        self.assertEqual(matrix_diff_matrixID(self.zeroMatID, self.randMatAID), scalar_mult_matrixID(self.neg1, self.randMatAID))

    def test_diff_double(self):
        """
        A - (-A) = A + A
        """
        a = matrix_diff_matrixID(self.randMatAID, scalar_mult_matrixID(self.neg1, self.randMatAID))
        b = matrix_add_matrixID(self.randMatAID, self.randMatAID)
        self.assertEqual(a,b)

    def test_diff_bin_invert(self):
        """
        A - B = -(B - A)
        """
        a = matrix_diff_matrixID(self.randMatAID, self.randMatBID)
        b = matrix_diff_matrixID(self.randMatBID, self.randMatAID)
        self.assertEqual(a, scalar_mult_matrixID(self.neg1, b))

    ##TODO we want to catch a failure here, but instead we get a system exit. 
    #we could intentionally skip the test but not in python 2.6. 
    #@unittest.skip("not supported") # *.skip is python 2.7 and above
    #def test_diff_size_mismatch(self):
    #    """
    #    Fail if size mismatch.
    #    """
    #    # generate matrix with one less row
    #    randMatCID = matrix_quick_build.matrix_random_matrixID(self.scalarType, self.row_level-1, self.col_level, self.val_range[0], self.val_range[1], self.sparsity)
    #    # generate matrix with one less col
    #    randMatCID = matrix_quick_build.matrix_random_matrixID(self.scalarType, self.row_level-1, self.col_level, self.val_range[0], self.val_range[1], self.sparsity)
    #    #with self.assertRaises(
    #    self.assertEqual(self.randMatAID, matrix_diff_matrixID(self.randMatAID, randMatCID))


#needed for python 2.6 (without unittest2)
if __name__ == "__main__":
    unittest.main()


