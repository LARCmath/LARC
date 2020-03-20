#!/usr/bin/env python3

# NOTE: to run on command line
#   python3 -m unittest -v test_unittest_matrix_store

 ##################################################################
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
 ##################################################################

from __future__ import print_function

import os 
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../../src"))
from pylarc import *
from ctypes import *
import unittest

class TestMacroFunctions(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # in fear of memory leaks, we'll run this once per class. 
        # with pylarc.stdout_redirected():
        #    initialize_larc(26,24,10,-1,-1,1)
        initialize_larc(26,24,10,-1,-1,0)  # run quiet
        
    def setUp(self):
        # Define string for use in formating filenames
        # (cvar is in pylarc)
        self.scalarType = cvar.scalarTypeStr

    def test_row_level(self):
        """
        Row level of matrix with row level 3 should be 3.
        """
        row_level = 3
        col_level = 4
        randMatID = pylarc.matrix_random_matrixID(self.scalarType, row_level, col_level, 0, 100)
        self.assertEqual(matrix_row_level_matrixID(randMatID), row_level)

    def test_col_level(self):
        """
        Col level of matrix with col level 4 should be 4.
        """
        row_level = 3
        col_level = 4
        randMatID = pylarc.matrix_random_matrixID(self.scalarType, row_level, col_level, 0, 100)
        self.assertEqual(matrix_row_level_matrixID(randMatID), row_level)

    def test_sub_matrix(self):
        """
        One of the submatrices should have one less level. 
        """
        row_level = 3
        col_level = 3
        randMatID = pylarc.matrix_random_matrixID(self.scalarType, row_level, col_level, 0, 100, .5)
        subMatID = matrix_sub_matrixID(randMatID, 0)
        #print_naive_by_matID(randMatID)
        #print_naive_by_matID(subMatID)
        self.assertEqual(matrix_row_level_matrixID(subMatID), matrix_row_level_matrixID(randMatID)-1)




class TestMatrixGetMatrixIDFromPanel(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # in fear of memory leaks, we'll run this once per class. 
        # with pylarc.stdout_redirected():
        #    initialize_larc(26,24,10,-1,-1,1)
        initialize_larc(26,24,10,-1,-1,0)  # run quiet

    def setUp(self):
        # Define string for use in formating filenames
        # (cvar is in pylarc)
        self.scalarType = cvar.scalarTypeStr

    def test_zero_already_present(self):
        """
        Building a zero matrix of level 4 from level 3 zero matrices
        should return the zero matrix of level 4.
        """
        level = 3
        zeroMatID = pylarc.get_zero_matrixID(level, level)
        bigZeroMatID = get_matID_from_four_subMatIDs(zeroMatID, zeroMatID, zeroMatID, zeroMatID, level+1, level+1)
        self.assertEqual(bigZeroMatID, get_zero_matrixID(level+1, level+1))

    def test_rebuild_rand_matrix(self):
        """
        Rebuilding a matrix from it's submatrices should return the original 
        matrix.
        """
        level = 4
        randMatID = pylarc.matrix_random_matrixID(self.scalarType, level, level, 0, 100, .5)
        subMatIDs = [matrix_sub_matrixID(randMatID, i) for i in range(4)]
        rebuildRandMatID = get_matID_from_four_subMatIDs(*(subMatIDs + [level, level]))
        self.assertEqual(randMatID, rebuildRandMatID)

    def test_row_vector(self):
        """
        check to see if level 3 zero vectors together are equiv to level 4 zero vector
        """
        level = 4
        MATRIX_ID_INVALID = -1
        zeroHalfVecID = pylarc.get_zero_matrixID(0, level-1)
        zeroVecID = pylarc.get_zero_matrixID(0, level)
        targetVecID = get_matID_from_four_subMatIDs(zeroHalfVecID, zeroHalfVecID, MATRIX_ID_INVALID, MATRIX_ID_INVALID, 0, level)
        self.assertEqual(zeroVecID, targetVecID)

    def test_col_vector(self):
        """
        Make a zero row vector. 
        """
        level = 4
        MATRIX_ID_INVALID = -1
        zeroHalfVecID = pylarc.get_zero_matrixID(level-1, 0)
        zeroVecID = pylarc.get_zero_matrixID(level, 0)
        targetVecID = get_matID_from_four_subMatIDs(zeroHalfVecID, MATRIX_ID_INVALID, zeroHalfVecID, MATRIX_ID_INVALID, level, 0)
        self.assertEqual(zeroVecID, targetVecID)


