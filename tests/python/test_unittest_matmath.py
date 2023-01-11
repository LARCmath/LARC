#!/usr/bin/env python3

# NOTE: to run on command line
#   python3 -m unittest -v test_unittest_matmath
# To run individual test classes from the module, do (for example) 
#   python3 -m unittest -v test_unittest_matmath.TestMatrixMaxnorm

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

 # This test set looks at routines in matmath.c and compares
 # their output with an output calculated in python.

from __future__ import print_function

import os 
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../../src"))
from pylarc import *
from ctypes import *
import random
import unittest
import math


# Python routine to compute the complex conjugate of a number:
def conj(val):
    if type(val) == type(1j):
        return val.conjugate()
    elif type(val) == type(1.0):
        return val
    elif type(val) == type(1):
        return val
    else:
        assert False, "Unknown value type {0} of {1}.".format(type(val), val)


# Python routine to take matrix entry list and square each entry
# then multiply each entry by scale_factor
def squareEntries(entries, scale_factor):
    return [scale_factor*(conj(entry)*entry) for entry in entries]


# tests the matmath.c routine  matrix_entrySquared
class TestMatrixEntrySquared(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # in fear of memory leaks, we'll run this one per class. 
        # with stdout_redirected():
        #    initialize_larc(26,24,10,-1,-1,1)
        initialize_larc(26,24,10,-1,-1,0)  # run quiet

    def setUp(self):
        # Define string for use in formating filenames
        # (cvar is in pylarc)
        self.scalarType = cvar.scalarTypeStr
        # set testing levels
        self.row_level = 2
        self.col_level = 2
        self.val_range= [-100, 100]
        self.sparsity = 0.3

    def tearDown(self):
        # clean the matrix store after every test
        clean_matrix_storage()
        #clean_hash_store()
        #clean_op_store()

    def test_entrysquared_random_zeroscale(self):
        """
        0 * entrysquared(A) = 0 should work for random A.
        """
        scale_factor = 0
        n = (2**self.row_level)*(2**self.col_level)

        vals = matrix_random_entries(self.scalarType, n, self.val_range[0], self.val_range[1], self.sparsity)
        squared_vals = [0 for val in vals]

        vals_matID = row_major_list_to_store(map_to_str(vals, self.scalarType), self.row_level, self.col_level, 2**self.col_level)
        squared_vals_matID = row_major_list_to_store(map_to_str(squared_vals, self.scalarType), self.row_level, self.col_level, 2**self.col_level)
        self.assertEqual(matrix_entrySquared(vals_matID, str(scale_factor)), squared_vals_matID)

    @unittest.skipIf(cvar.scalarTypeStr in ("Real", "Complex", "MPReal", "MPComplex"), "ScalarType must not be Real, Complex, MPReal or MPComplex")
    def test_entrysquared_random_noscale(self):
        """
        1 * entrysquared(A) should work for random A.
        """
        scale_factor = 1
        n = (2**self.row_level)*(2**self.col_level)

        vals = matrix_random_entries(self.scalarType, n, self.val_range[0], self.val_range[1], self.sparsity)
        squared_vals = squareEntries(vals, 1)

        vals_matID = row_major_list_to_store(map_to_str(vals, self.scalarType), self.row_level, self.col_level, 2**self.col_level)
        squared_vals_matID = row_major_list_to_store(map_to_str(squared_vals, self.scalarType), self.row_level, self.col_level, 2**self.col_level)
        squared_matID = matrix_entrySquared(vals_matID, "1") 
        failure_msg = "matrices displayed above"
        if squared_matID != squared_vals_matID:
            print("Failure in test_entrysquared_random_noscale - printing relevant matrices.")
            print("Original matrix:")
            print_naive(vals_matID)
            print("Squared (by python) matrix:")
            print_naive(squared_vals_matID)
            print("Squared (by LARC) matrix:")
            print_naive(squared_matID)
        self.assertEqual(squared_matID, squared_vals_matID, failure_msg)

    @unittest.skipIf(cvar.scalarTypeStr in ("Real", "Complex", "MPReal", "MPComplex"), "ScalarType must not be Real, Complex, MPReal or MPComplex")
    def test_entrysquared_random_scale(self):
        """
        s * entrysquared(A) should work for random A.
        """
        scale_factor = matrix_random_entry(self.scalarType, -10, 10)
        n = (2**self.row_level)*(2**self.col_level)

        vals = matrix_random_entries(self.scalarType, n, self.val_range[0], self.val_range[1], self.sparsity)
        squared_vals = squareEntries(vals, scale_factor)

        vals_matID = row_major_list_to_store(map_to_str(vals, self.scalarType), self.row_level, self.col_level, 2**self.col_level)
        squared_vals_matID = row_major_list_to_store(map_to_str(squared_vals, self.scalarType), self.row_level, self.col_level, 2**self.col_level)
        squared_matID = matrix_entrySquared(vals_matID, value_to_string(scale_factor, self.scalarType))
        failure_msg = "matrices displayed above"
        if squared_matID != squared_vals_matID:
            print("Failure in test_entrysquared_random_scale - printing relevant matrices.")
            print("Scale factor: {0} ({1})".format(scale_factor, scale_factor.hex()))
            print("Original matrix:")
            print_naive(vals_matID)
            print("Squared (by python) matrix:")
            print_naive(squared_vals_matID)
            print("Squared (by LARC) matrix:")
            print_naive(squared_matID)
        self.assertEqual(squared_matID, squared_vals_matID, failure_msg)


class TestMatrixMaxnorm(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # in fear of memory leaks, we'll run this one per class. 
        # with stdout_redirected():
        #    initialize_larc(26,24,10,-1,-1,1)
        initialize_larc(26,24,10,-1,-1,0)  # run quiet

    def setUp(self):
        # Define string for use in formating filenames
        # (cvar is in pylarc)
        self.scalarType = cvar.scalarTypeStr
        # set testing levels
        self.row_level = 3
        self.col_level = 3
        self.val_range= [-100, 100]
        self.sparsity = 0.3

    def tearDown(self):
        # clean the matrix store after every test
        clean_matrix_storage()
        #clean_hash_store()
        #clean_op_store()

    def test_maxnorm_zero(self):
        """
        maxnorm(0) = 0
        """
        #* grab Zero matrix
        zeroMatID = get_zero_pID(self.row_level, self.col_level)
        maxnormID = normID(zeroMatID,0)
        maxnorm = traceID(maxnormID)
        if self.scalarType in ("Complex", "MPRatComplex", "MPComplex",
                                "Clifford"):
#            maxnorm = sum(map(float, maxnorm.split("+I*")))
            maxnormList = list(map(float, maxnorm.split("+I*")))
            if (len(maxnormList)==2):
                self.assertEqual(maxnormList[1],0)
            maxnorm = maxnormList[0]
        else:
            maxnorm = float(maxnorm)
        self.assertEqual(maxnorm, 0)

    def test_maxnorm_id(self):
        """
        maxnorm(id) = 1
        """
        idMatID = get_identity_pID(self.row_level)
        maxnormID = normID(idMatID,0)
        maxnorm = traceID(maxnormID)
        if self.scalarType in ("Complex", "MPRatComplex", "MPComplex",
                               "Clifford"):
#            maxnorm = sum(map(float, maxnorm.split("+I*")))
            maxnormList = list(map(float, maxnorm.split("+I*")))
            if (len(maxnormList)==2):
                self.assertEqual(maxnormList[1],0.0)
            maxnorm = maxnormList[0]
        else:
            maxnorm = float(maxnorm)
        self.assertEqual(maxnorm, 1)

    def test_maxnorm_negated(self):
        """
        maxnorm(-A) = maxnorm(A)
        """
        # build -1
        neg1 = get_valID_from_valString("-1")
        # generate random matrix
        randMatAID = matrix_random_matrixID(self.scalarType, self.row_level, self.col_level, self.val_range[0], self.val_range[1], self.sparsity)
        negAMatID = scalar_mult(neg1, randMatAID)
        self.assertEqual(normID(randMatAID,0), normID(negAMatID,0))

    def test_maxnorm_random(self):
        """
        placing != entry 5 in zero matrix to make A
        maxnorm(A) = 5
        """
        # get a zero matrix and place a 5 some where in the matrix
        zeroMatID = get_zero_pID(self.row_level, self.col_level)
        # TO BE VERY CONFUSING we need an (i,j) index in the matrix we use:
        # for i: we know row_level <= 2^row_level so the index i=self.row_level is a valid row index
        # for j: we know col_level <= 2^col_level so the index j=self.col_level is a valid col index
        matAID = replace_scalar_in_matrix_by_string_and_coords(zeroMatID, self.row_level, self.col_level, "5");
        maxnormID = normID(matAID,0)
        maxnorm = traceID(maxnormID)
        if self.scalarType in ("Complex", "MPRatComplex", "MPComplex",
                               "Clifford"):
#            maxnorm = sum(map(float, maxnorm.split("+I*")))
            maxnormList = list(map(float, maxnorm.split("+I*")))
            if (len(maxnormList)==2):
                self.assertEqual(maxnormList[1],0.0)
            maxnorm = maxnormList[0]
        else:
            maxnorm = float(maxnorm)
        self.assertEqual(maxnorm, 5)


class TestMatrixL2norm(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # in fear of memory leaks, we'll run this one per class. 
        # with stdout_redirected():
        #    initialize_larc(26,24,10,-1,-1,1)
        initialize_larc(26,24,10,-1,-1,0)  # run quiet

    def setUp(self):
        # Define string for use in formating filenames
        # (cvar is in pylarc)
        self.scalarType = cvar.scalarTypeStr
        # set testing levels
        self.row_level = 3
        self.col_level = 3
        self.val_range= [-100, 100]
        self.sparsity = 0.3

    def tearDown(self):
        # clean the matrix store after every test
        clean_matrix_storage()
        #clean_hash_store()
        #clean_op_store()

    def test_L2norm(self):
        """
        l2Norm(0) = 0
        """
        #* grab Zero matrix
        zeroMatID = get_zero_pID(self.row_level, self.col_level)
#        maxnormID = normID(zeroMatID,0)
#        maxnorm = traceID(maxnormID)
        l2NormID = normID(zeroMatID,2)
        l2Norm = traceID(l2NormID)
        if self.scalarType in ("Complex", "MPRatComplex", "MPComplex",
                               "Clifford"):
#            l2Norm = sum(map(float, l2Norm.split("+I*")))
            l2NormList = list(map(float, l2Norm.split("+I*")))
            if (len(l2NormList)==2):
                self.assertEqual(l2NormList[1],0)
            l2Norm = l2NormList[0]
        else:
            l2Norm = float(l2Norm)
        self.assertEqual(l2Norm, 0)

    @unittest.skipIf(cvar.scalarTypeStr in ("Integer", "MPInteger"), "ScalarType must not be Integer or MPInteger")
    def test_l2Norm_id(self):
        """
        l2Norm(id) = sqrt(8.0)
        """
        idMatID = get_identity_pID(self.row_level)
        l2NormID = normID(idMatID,2)
        l2Norm = traceID(l2NormID)
        if self.scalarType in ("Complex", "MPComplex"):
            # l2Norm = sum(map(float, l2Norm.split("+I*")))
            l2NormList = list(map(float, l2Norm.split("+I*")))
            # norm should be real
            if (len(l2NormList) == 2):
                self.assertEqual(l2NormList[1],0.0)
            l2Norm = l2NormList[0]
        elif self.scalarType in ("MPRational", "MPRatComplex","Clifford"):
            if self.scalarType in ("MPRatComplex", "Clifford"):
                #    need to split complex into real and imaginary
                #    and confirm imaginary part is zero
                l2NormList = l2Norm.split("+I*")
                if (len(l2NormList) == 2):
                    l2ImagList = l2NormList[1].split("/") 
                    if (len(l2ImagList) == 2):
                        l2Imag = float(l2ImagList[0])/float(l2ImagList[1])
                    else:
                        l2Imag = float(l2ImagList[0])
                    self.assertEqual(l2Imag,0.0)
                l2Norm = l2NormList[0]
            #    need to split real rational into numerator and denominator
            l2NormList = l2Norm.split("/")
            l2Norm = float(l2NormList[0])/float(l2NormList[1])
        else:
            l2Norm = float(l2Norm)
        self.assertEqual(l2Norm, math.sqrt(8.0))

    def test_l2Norm_negated(self):
        """
        maxnorm(-A) = maxnorm(A)
        """
        # build -1
        neg1 = get_valID_from_valString("-1")
        # generate random matrix
        randMatAID = matrix_random_matrixID(self.scalarType, self.row_level, self.col_level, self.val_range[0], self.val_range[1], self.sparsity)
        negAMatID = scalar_mult(neg1, randMatAID)
        self.assertEqual(normID(randMatAID,2), normID(negAMatID,2))

    def test_l2Norm_random(self):
        """
        placing != entry 5 in zero matrix to make A
        maxnorm(A) = 5
        """
        # get a zero matrix and place a 5 some where in the matrix
        zeroMatID = get_zero_pID(self.row_level, self.col_level)
        # TO BE VERY CONFUSING we need an (i,j) index in the matrix we use:
        # for i: we know row_level <= 2^row_level so the index i=self.row_level is a valid row index
        # for j: we know col_level <= 2^col_level so the index j=self.col_level is a valid col index
        matAID = replace_scalar_in_matrix_by_string_and_coords(zeroMatID, self.row_level, self.col_level, "5");
        l2NormID = normID(matAID,2)
        l2Norm = traceID(l2NormID)
        if self.scalarType in ("Complex", "MPRatComplex", "MPComplex",
                               "Clifford"):
#            l2Norm = sum(map(float, l2Norm.split("+I*")))
            l2NormList = list(map(float, l2Norm.split("+I*")))
            if (len(l2NormList)==2):
                self.assertEqual(l2NormList[1],0)
            l2Norm = l2NormList[0]
        else:
            l2Norm = float(l2Norm)
        self.assertEqual(l2Norm, 5)


class TestMatrixSparsity(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # in fear of memory leaks, we'll run this one per class. 
        # with stdout_redirected():
        #    initialize_larc(26,24,10,-1,-1,1)
        initialize_larc(26,24,10,-1,-1,0)  # run quiet

    def setUp(self):
        # Define string for using in formating filenames
        # (cvar is in pylarc)
        self.scalarType = cvar.scalarTypeStr
        # set testing levels
        self.row_level = 3
        self.col_level = 3
        self.val_range= [-100, 100]

    def tearDown(self):
        # clean the matrix store after every test
        clean_matrix_storage()
        #clean_hash_store()
        #clean_op_store()

    def test_sparse_val(self):
        """
        Verify sparsity(A) = sparsity value for generating random matrix A. 
        """
        sparsity = 0
        n = (2**self.row_level)*(2**self.col_level)
        # get a nonzero random sparsity
        while sparsity == 0:
            sparsity = random.random()
        randMatID = matrix_random_matrixID(self.scalarType, self.row_level, self.col_level, self.val_range[0], self.val_range[1], sparsity)
        # convert to number of zeros
        # requested sparsity may be inaccurate by less than one zero 
        self.assertEqual(int(round(sparsity * n)), int(matrix_count_entries(randMatID, "0")))

    def test_sparse_0(self):
        """
        Matrix generated with sparcity value 0 has sparsity 0.
        """
        onesMatID = matrix_random_matrixID(self.scalarType, self.row_level, self.col_level, 1, 2, 0)
        self.assertEqual(int(matrix_count_entries(onesMatID, "0")), 0)

    def test_sparse_id(self):
        """
        Identity matrix of order n has sparsity (n-1)/n.
        """
        idMatID = get_identity_pID(self.row_level)
        self.assertEqual(int(matrix_count_entries(idMatID, "0")), (2**self.row_level-1)*(2**self.row_level))

    def test_sparse_0_mat(self):
        """
        Zero matrix has sparsity 0.
        """
        zeroMatID = get_zero_pID(self.row_level, self.col_level)
        self.assertEqual(int(matrix_count_entries(zeroMatID, "0")), (2**self.row_level) * (2**self.col_level))

class TestMatrixSum(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # in fear of memory leaks, we'll run this one per class. 
        # with stdout_redirected():
        #    initialize_larc(26,24,10,-1,-1,1)
        initialize_larc(26,24,10,-1,-1,0)  # run quiet

    def setUp(self):
        # Define string for using in formating filenames
        # (cvar is in pylarc)
        self.scalarType = cvar.scalarTypeStr
        # set testing levels
        self.row_level = 3
        self.col_level = 3
        self.val_range= [-100, 100]
        self.sparsity = .3
        # get Zero matrix
        self.zeroMatID = get_zero_pID(self.row_level, self.col_level)
        # generate two random matrices of same size
        self.randMatAID = matrix_random_matrixID(self.scalarType, self.row_level, self.col_level, self.val_range[0], self.val_range[1], self.sparsity)
        self.randMatBID = matrix_random_matrixID(self.scalarType, self.row_level, self.col_level, self.val_range[0], self.val_range[1], self.sparsity)
        # build -1
        self.neg1 = get_valID_from_valString("-1")

    def tearDown(self):
        # clean the matrix store after every test
        clean_matrix_storage()
        #clean_hash_store()
        #clean_op_store()

    def test_add_additive_inverse(self):
        """
        A + (-A) = 0
        """
        a = scalar_mult(self.neg1, self.randMatAID)
        self.assertEqual(matrix_add(self.randMatAID, a), self.zeroMatID)

    def test_add_add_zero(self):
        """
        A + 0 = A
        """
        self.assertEqual(matrix_add(self.randMatAID, self.zeroMatID), self.randMatAID)

    def test_add_double(self):
        """
        A + A = 2A
        """
        scalMatID = get_valID_from_valString("2")
        prodMatID = scalar_mult(scalMatID, self.randMatAID)
        self.assertEqual(matrix_add(self.randMatAID, self.randMatAID), prodMatID)

    def test_add_commutativity(self):
        """
        A + B = B + A
        """
        a = matrix_add(self.randMatAID, self.randMatBID)
        b = matrix_add(self.randMatBID, self.randMatAID)
        self.assertEqual(a,b)

    def test_add_commutativity_op_shortcut(self):
        """
        If A+B in op store, calculating B+A shouldn't require additional calculations.
        """
        a = matrix_add(self.randMatAID, self.randMatBID)
        num_matrices_made = num_matrices_created()
        b = matrix_add(self.randMatBID, self.randMatAID)
        num_matrices_made2 = num_matrices_created()
        self.assertEqual(num_matrices_made, num_matrices_made2)


class TestMatrixDifference(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # in fear of memory leaks, we'll run this one per class. 
        # with stdout_redirected():
        #    initialize_larc(26,24,10,-1,-1,1)
        initialize_larc(26,24,10,-1,-1,0)  # run quiet

    def setUp(self):
        # Define string for using in formating filenames
        # (cvar is in pylarc)
        self.scalarType = cvar.scalarTypeStr
        # set testing levels
        self.row_level = 3
        self.col_level = 3
        self.val_range= [-100, 100]
        self.sparsity = .3
        # get a Zero matrix
        self.zeroMatID = get_zero_pID(self.row_level, self.col_level)
        # generate two random matrices of same size
        self.randMatAID = matrix_random_matrixID(self.scalarType, self.row_level, self.col_level, self.val_range[0], self.val_range[1], self.sparsity)
        self.randMatBID = matrix_random_matrixID(self.scalarType, self.row_level, self.col_level, self.val_range[0], self.val_range[1], self.sparsity)
        # build -1
        self.neg1 = get_valID_from_valString("-1")

    def tearDown(self):
        # clean the matrix store after every test
        clean_matrix_storage()
        #clean_hash_store()
        #clean_op_store()

    def test_diff_additive_inverse(self):
        """
        A - A = 0
        """
        self.assertEqual(self.zeroMatID, matrix_diff(self.randMatAID, self.randMatAID))

    def test_diff_subtract_zero(self):
        """
        A - 0 = A
        """
        self.assertEqual(self.randMatAID, matrix_diff(self.randMatAID, self.zeroMatID))

    def test_diff_uni_invert(self):
        """
        0 - A = -A
        """
        self.assertEqual(matrix_diff(self.zeroMatID, self.randMatAID), scalar_mult(self.neg1, self.randMatAID))

    def test_diff_double(self):
        """
        A - (-A) = A + A
        """
        a = matrix_diff(self.randMatAID, scalar_mult(self.neg1, self.randMatAID))
        b = matrix_add(self.randMatAID, self.randMatAID)
        self.assertEqual(a,b)

    def test_diff_bin_invert(self):
        """
        A - B = -(B - A)
        """
        a = matrix_diff(self.randMatAID, self.randMatBID)
        b = matrix_diff(self.randMatBID, self.randMatAID)
        self.assertEqual(a, scalar_mult(self.neg1, b))

    @unittest.skip("not supported") # *.skip is python 2.7 and above
    def test_diff_size_mismatch(self):
        """
        Fail if size mismatch.
        """
        # generate matrix with one less row
        randMatCID = matrix_random_matrixID(self.scalarType, self.row_level-1, self.col_level, self.val_range[0], self.val_range[1], self.sparsity)
        # generate matrix with one less col
        randMatCID = matrix_random_matrixID(self.scalarType, self.row_level-1, self.col_level, self.val_range[0], self.val_range[1], self.sparsity)
        # could try using self.assertRaises() to show an exception is raised
        self.assertEqual(self.randMatAID, matrix_diff(self.randMatAID, randMatCID))


class TestMatrixQuotient(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # in fear of memory leaks, we'll run this one per class. 
        # with stdout_redirected():
        #    initialize_larc(26,24,10,-1,-1,1)
        initialize_larc(26,24,10,-1,-1,0)  # run quiet

    def setUp(self):
        # Define string for using in formating filenames
        # (cvar is in pylarc)
        self.scalarType = cvar.scalarTypeStr

        # Choose non-zero value
        if self.scalarType in ("Complex", "MPComplex",):
            self.nz_value_str = "1.7384043-I*3.2171"
        elif self.scalarType in ("MPRatComplex",):
            self.nz_value_str = "348734/4457457-I*48484/55341"
        elif self.scalarType in ("Real", "MPReal",):
            self.nz_value_str = "1.7384043"
        elif self.scalarType in ("MPRational",):
            self.nz_value_str = "348734/4457457"
        elif self.scalarType in ("Integer",):
            self.nz_value_str = "348734"
        elif self.scalarType in ("MPInteger",):
            self.nz_value_str = "-34847385924385794875999277366366478588861875091875757281766734"
        elif self.scalarType in ("Clifford",):
            # here is a default string that works for all cases
            self.nz_value_str = "(32/7)*{1}"
#            # the following string works for Q[S2,S3]
#            self.nz_value_str = "(32/7)*{1}-(15/16)*{S2}+(7/22)*{S3}"
#            # the following string works for Q[i,S2,S3]
#            self.nz_value_str = "(32/7)*{1}+(1/17)*I-(15/16)*{S2}+2*{S2}*I+(7/22)*{S3}+(9/7)*{S6}*I"
#            # the following string works for Q[i,C2]
#            self.nz_value_str = "(32/7)*{1}+(1/17)*I-(15/16)*{C2}+2*{C4}*I"
        else:
            self.nz_value_str = "0"

        self.nzID = get_valID_from_valString(self.nz_value_str)

    def tearDown(self):
        # clean the matrix store after every test
        clean_matrix_storage()
        #clean_hash_store()
        #clean_op_store()

    def test_one(self):
        mat34 = matrix_constant_entry_matrixID(self.nz_value_str, 3, 4)
        one34 = matrix_constant_entry_matrixID("1", 3, 4)
        ans34 = scalar_divide(mat34, self.nzID)
        self.assertEqual(ans34, one34)

    def test_two(self):
        mat00 = matrix_constant_entry_matrixID(self.nz_value_str, 0, 0)
        one00 = matrix_constant_entry_matrixID("1", 0, 0)
        ans00 = scalar_divide(mat00, self.nzID)
        self.assertEqual(ans00, one00)

    def test_three(self):
        mat22 = matrix_constant_entry_matrixID(self.nz_value_str, 2, 2)
        one22 = matrix_constant_entry_matrixID("1", 2, 2)
        ans22 = scalar_divide(mat22, self.nzID)
        self.assertEqual(ans22, one22)


class TestMatrixValueCount(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # in fear of memory leaks, we'll run this one per class. 
        # with stdout_redirected():
        #    initialize_larc(26,24,10,-1,-1,1)
        initialize_larc(26,24,10,-1,-1,0)  # run quiet

    def setUp(self):
        # Define string for using in formating filenames
        # (cvar is in pylarc)
        self.scalarType = cvar.scalarTypeStr
        # set testing levels
        self.row_level = 3
        self.col_level = 3
        self.val_range= [-100, 100]
        self.sparsity = .3

    def tearDown(self):
        # clean the matrix store after every test
        clean_matrix_storage()
        #clean_hash_store()
        #clean_op_store()

    def test_count0_sparsity_rand(self):
        """
        Verify count0(A) = sparsity*size value for generating random matrix A. 
        """
        sparsity = 0
        n = (2**self.row_level)*(2**self.col_level)
        # get a nonzero random sparsity
        while sparsity == 0:
            sparsity = random.random()
        randMatID = matrix_random_matrixID(self.scalarType, self.row_level, self.col_level, self.val_range[0], self.val_range[1], sparsity)
        # convert to number of zeros
        # requested sparsity may be inaccurate by less than one zero 
        self.assertEqual(int(round(sparsity * n)), int(matrix_count_entries(randMatID, "0")))

    def test_count0_sparse_0(self):
        """
        Matrix generated with sparcity value 0 has 0 zero entries.
        """
        onesMatID = matrix_random_matrixID(self.scalarType, self.row_level, self.col_level, 1, 2, 0)
        self.assertEqual(0, int(matrix_count_entries(onesMatID, "0")))

    def test_count0_id_mat(self):
        """
        nxn identiy matrix has n^2 - n zeros.
        """
        idMatID = get_identity_pID(self.row_level)
        self.assertEqual(int(matrix_count_entries(idMatID, "0")), ((2**self.row_level) - 1)*(2**self.row_level))

    def test_count1_id_mat(self):
        """
        nxn identiy matrix n ones.
        """
        idMatID = get_identity_pID(self.row_level)
        self.assertEqual(int(matrix_count_entries(idMatID, "1")), 2**self.row_level)

    def test_count0_0_mat(self):
        """
        nxn zero matrix has n^2 zeros.
        """
        n = (2**self.row_level)*(2**self.col_level);
        zeroMatID = get_zero_pID(self.row_level, self.col_level)
        self.assertEqual(n, int(matrix_count_entries(zeroMatID, "0")))


#needed for python 2.6 (without unittest2)
if __name__ == "__main__":
    unittest.main()


