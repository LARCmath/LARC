#!/usr/bin/env python3

# NOTE: to run on command line
#   python3 -m unittest -v test_unittest_matrix_store

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
sys.path.append(os.path.join(os.path.dirname(__file__),"../../src"))
from pylarc import *
from ctypes import *
import unittest

class TestOpStoreCleaning(unittest.TestCase):

    def setUp(self):
        # in fear of memory leaks, we'll run this once per class. 
        # Define string for use in formating filenames
        # (cvar is in pylarc)
        self.scalarType = cvar.scalarTypeStr
        self.level = 3
        self.val_range = [-100, 100]
        self.sparsity = 0.5
        # with stdout_redirected():
        #    initialize_larc(26,24,10,-1,-1,1)
        initialize_larc(26,24,10,-1,-1,0)    # this should have run quiet effect
        self.randMatAID = matrix_random_matrixID(self.scalarType, self.level, self.level, self.val_range[0], self.val_range[1], self.sparsity)
        self.randMatBID = matrix_random_matrixID(self.scalarType, self.level, self.level, self.val_range[0], self.val_range[1], self.sparsity)
        op_store_report("stdout");
        self.prodMatABID = matrix_mult(self.randMatAID, self.randMatBID)
        self.sumMatABID  = matrix_add(self.randMatAID, self.randMatBID)
        self.sumMatAAID  = matrix_add(self.randMatAID, self.randMatAID)

    @unittest.skip("hide") 
    def test_clean_ops_no_need(self):
        """
        If we just run clean_op_store, there shouldn't be any operations with
        invalid matrices to remove. 
        """
        op_store_report("stdout");
        clean_op_store();
        op_store_report("stdout");
        empty_op_store();
        op_store_report("stdout");
        self.assertEqual(1, 1)

    @unittest.skip("hide") 
    def test_clean_ops_mat_removal(self):
        """
        If we delete matrix B from the op store, there should be additions left
        but no multiplications. 
        """
        op_store_report("stdout");
        remove_matrix_from_store(self.randMatBID)
        clean_op_store();
        op_store_report("stdout");
        empty_op_store();
        op_store_report("stdout");
        self.assertEqual(1, 1)

    @unittest.skip("hide") 
    def test_clean_ops_all_removal(self):
        """
        If we delete matrix B from the op store, there should be additions left
        but no multiplications. 
        """
        op_store_report("stdout");
        clean_matrix_storage();
        clean_op_store();
        op_store_report("stdout");
        empty_op_store();
        op_store_report("stdout");
        self.assertEqual(1, 1)

    @unittest.skip("hide") 
    def test_resume_post_clean(self):
        """
        If we clean the matrix store, we should see operations come back if we
        do more operations afterward.
        """
        op_store_report("stdout");
        empty_op_store();
        op_store_report("stdout");
        # now we need to recalculate these
        self.prodMatABID = matrix_mult(self.randMatAID, self.randMatBID)
        self.sumMatABID  = matrix_add(self.randMatAID, self.randMatBID)
        self.sumMatAAID  = matrix_add(self.randMatAID, self.randMatAID)
        # note that there are now (KRO +-19) twice as many records created 
        # as current records since we recreated them all. 
        op_store_report("stdout");
        self.assertEqual(1, 1)
