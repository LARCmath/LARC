#!/usr/bin/env python3

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
import pylarc
import numpy as np
from ctypes import *

if __name__ == '__main__':
	# This version references matrices by matrixID instead of pointers
	
    print("This code tests some basic matrix building and reading routines\n")

    # initialize larc
    mat_store_exp = 26
    op_store_exp = 24
    max_level = 10
    rnd_sig_bits = -1   # default value
    trunc_to_zero_bits = -1  # default value
    pylarc.create_report_thread(1800)
    verbose = 1
    pylarc.initialize_larc(mat_store_exp,op_store_exp,max_level,rnd_sig_bits,trunc_to_zero_bits,verbose)


    ##############################################################
    ##  In the Makefile you can compile with:                   ##
    ##     TYPE=INTEGER, TYPE=REAL, TYPE=COMPLEX,               ##
    ##  or with multiprecision types:                           ##
    ##     TYPE=MPINTEGER, TYPE=MPRATIONAL, or TYPE=MPRatComplex   ##
    ##  Define string for using in formating filenames
    scalarTypeStr = pylarc.cvar.scalarTypeStr

    # Calculate number of matrices created, then print part of matrix store
    num_matrices_made = pylarc.num_matrices_created()
    print("\n%d matrices have been created" %num_matrices_made)
    end = num_matrices_made - 1
    filename = "../dat/out/preload.%s.store" %scalarTypeStr
    pylarc.matrix_store_info_to_file(0,end,os.path.join(os.path.dirname(__file__),filename),"After preload with parameters: 26, 24, 10.")

    # build array in C from Python list of scalars
    print("Using row_major_list_to_store on data entered from python\n")

    # create a matrix in python
    if scalarTypeStr in ("Integer", "MPInteger"):
        a = np.matrix([[1, 3, 0, 0],
                       [8, 6, 0, 0],
                       [0, 0, 13, 15],
                       [0, 0, 12, 10]])
    elif scalarTypeStr in ("Complex", "MPComplex", "MPRatComplex"):
        a = np.matrix([[1+2j, 3+4j, 5+6j, 7+8j],
                       [8+7j, 6+5j, 3+4j, 1+2j],
                       [9+10j, 11+12j, 13+14j, 15+16j],
                       [16+15j, 14+13j, 12+11j, 10+9j]])
    elif scalarTypeStr in ("Real", "MPReal", "MPRational"):
        a = np.matrix([[1, 3, .5, 6],
                       [8, 6, 3, .1],
                       [-9, 11, 13, 1.5],
                       [16, 13, 12, 10]])
    else:
        raise Exception('Do not know how to build matrix for type %s.' %scalarTypeStr)

    # test our new function matrix_is_zero_matrixID(z_matID))
    z = np.matrix([[0, 0],[0, 0]])
    zlist = z.reshape(-1).tolist()[0]
    zlist_strs = pylarc.map_to_str(zlist, scalarTypeStr)
    z_matID = pylarc.row_major_list_to_store_matrixID(zlist_strs, 1, 1, 2)
    pylarc.print_naive_by_matID(z_matID)
    print("\nThe fucntion matrix_is_zero_matrixID returns:")
    print(pylarc.matrix_is_zero_matrixID(z_matID))
    print("\n")
    
    # turn the matrix into an array by reading off each row in turn (row major format)
    alist = a.reshape(-1).tolist()[0]
    alist_strs = pylarc.map_to_str(alist, scalarTypeStr)
    # arr = pylarc.buildArray(alist)
    # print 'arr:', pylarc.str_scalarTypeArray(arr, len(alist))
    print("alist_str:", alist_strs)

    # parameters for entering the python array into the store
    level = 2
    dim_whole = 2**level

    # creating or finding the matrix associated with the array
    serial = pylarc.row_major_list_to_store_matrixID(alist_strs, level, level, dim_whole)
    pylarc.print_naive_by_matID(serial)
    print("\n")

    # Make a parent matrix from four copies of the serial matrix
    print("Creating matrix from get_matID_from_four_subMatIDs on panel input and writing LARCMatrix file\n")

    zero2 = pylarc.get_zero_matrixID(2,2)
    
    # panel = [serial,zero2,zero2,serial]
    top_matID = pylarc.get_matID_from_four_subMatIDs(serial,zero2,zero2,serial,3,3)
    pylarc.print_naive_by_matID(top_matID)
    filename = "../dat/out/blocktest.%s.json" %scalarTypeStr
    pylarc.write_larcMatrix_file_by_matID(top_matID,os.path.join(os.path.dirname(__file__), filename))


    verbose = 1    # 0 run quiet, 1 warn if not block diagonal, 2 be chatty

    # run our test code for finding the matrix IDs of block diagonals
    blockLevel = 4
    print("\nTEST with blocklevel %d" %blockLevel)
    block_array = pylarc.list_block_diagonal_matrixIDs(top_matID,blockLevel,verbose)
    print(block_array)

    # run our test code for finding the matrix IDs of block diagonals
    blockLevel = 3
    print("\nTEST with blocklevel %d" %blockLevel)
    block_array = pylarc.list_block_diagonal_matrixIDs(top_matID,blockLevel,verbose)
    print(block_array)

    # run our test code for finding the matrix IDs of block diagonals
    blockLevel = 2
    print("\nTEST with blocklevel %d" %blockLevel)
    block_array = pylarc.list_block_diagonal_matrixIDs(top_matID,blockLevel,verbose)
    print(block_array)
    
    # run our test code for finding the matrix IDs of block diagonals
    blockLevel = 1
    print("\nTEST with blocklevel %d" %blockLevel)
    block_array = pylarc.list_block_diagonal_matrixIDs(top_matID,blockLevel,verbose)
    print(block_array)

    # run our test code for finding the matrix IDs of block diagonals
    blockLevel = 0
    print("\nTEST with blocklevel %d" %blockLevel)
    block_array = pylarc.list_block_diagonal_matrixIDs(top_matID,blockLevel,verbose)
    print(block_array)


