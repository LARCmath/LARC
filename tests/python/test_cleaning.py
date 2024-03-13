#!/usr/bin/env python3

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

from __future__ import print_function

import os
import sys 
sys.path.append(os.path.join(os.path.dirname(__file__),"../../src"))
import pylarc
import numpy as np
from ctypes import *

if __name__ == '__main__':
	# This version references matrices by matrixIDs instead of pointers
	
    print("This code tests some basic matrix building and reading routines\n")

    # initialize larc
    mat_store_exp = 15
    op_store_exp = 5
    max_level = 10
    regionbitparam = -1   # default value
    zeroregionbitparam = -1  # default value
    pylarc.create_report_thread(1800)
    verbose = 1
    pylarc.initialize_larc(mat_store_exp,op_store_exp,max_level,regionbitparam,zeroregionbitparam,verbose)


    # Define string for using in formating filenames
    scalarTypeStr = pylarc.cvar.scalarTypeStr

    # Calculate number of matrices created, then print part of matrix store
    num_matrices_made = pylarc.num_matrices_created()
    print("\n%d matrices have been created" %num_matrices_made)
    end = num_matrices_made - 1
    filename = "../dat/out/preload.%s.store" %scalarTypeStr
    pylarc.fprint_store_info_for_matrixID_range(0,end,os.path.join(os.path.dirname(__file__),filename),"After preload with parameters: 26, 24, 10.")

    #  PLAYING WITH PREWRITTEN NONSQUARE MATRIX
    filename = "../dat/in/sample.1.2.%s.json" %scalarTypeStr
    print("About to test read %s\n" %filename)
    samp_mID = pylarc.read_larcMatrixFile(os.path.join(os.path.dirname(__file__),filename))

    print("We read in the LARCMatrix file\n")
    pylarc.print_naive(samp_mID)
    print("\n")

    # build array in C from Python list of scalars
    print("Using row_major_list_to_store on data entered from python\n")

    # create a matrix in python
    if scalarTypeStr in ("Integer", "MPInteger"): 
        a = np.matrix([[1, 3, 5, 6],
                       [8, 6, 3, 1],
                       [-9, 11, 13, 15],
                       [16, 13, 12, 10]])
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

    mID = pylarc.add_numpy_matrix_to_matrix_store(a)
    
    # parameters for entering the python array into the store
    level = 2
    dim_whole = 2**level

    # creating or finding the matrix associated with the array
    #mID = pylarc.row_major_list_to_store(arr, level, level, dim_whole)
    pylarc.print_naive(mID)
    print("\n")

    # make a parent matrix from four copies of the matrixID matrix
    print("Creating matrix from get_pID_from_four_sub_pIDs on panel input and writing LARCMatrix file\n")
    panel = [mID]*4   # alternatively panel=[mID,mID,mID,mID]
    mID_parent = pylarc.get_pID_from_four_sub_pIDs(mID,mID,mID,mID,3,3)
    pylarc.print_naive(mID_parent)
    filename = "../dat/out/testfile.%s.json" %scalarTypeStr
    pylarc.fprint_larcMatrixFile(mID_parent,os.path.join(os.path.dirname(__file__), filename))
    
    #  PLAYING WITH PREWRITTEN REVERSIBLE NAND LARCMatrix FILE
    print("About to test read LARCMatrix file\n")
    filename = "../dat/in/nand.%s.json" %scalarTypeStr
    nand_mID = pylarc.read_larcMatrixFile(os.path.join(os.path.dirname(__file__),filename))
    print("We read in the LARCMatrix nand file\n")
    pylarc.print_naive(nand_mID)
    print("\n")

    # TESTING READING AND WRITING OF MATRICES
    print("Testing reading row major matrix format and writing files in LARCMatrix and naive format.\n")
    filename_rmm = "../dat/in/sample.1.1.%s.rmm" %scalarTypeStr
    filename_naive = "../dat/out/sample.1.1.%s.naive" %scalarTypeStr
    filename_json = "../dat/out/sample.1.1.%s.json" %scalarTypeStr
    
    sample_mID = pylarc.read_row_major_matrix_from_file(os.path.join(os.path.dirname(__file__),filename_rmm)) 
    pylarc.fprint_naive(sample_mID,os.path.join(os.path.dirname(__file__),filename_naive))
    pylarc.fprint_larcMatrixFile(sample_mID,os.path.join(os.path.dirname(__file__),filename_json))
    print("Printing out the rrm sample matrix in naive format to screen\n")
    pylarc.print_naive(sample_mID)

    print("Testing reading row major nonsquare matrices and writing files in LARCMatrix and naive format.\n")
    filename_rmm = "../dat/in/sample.1.2.%s.rmm" %scalarTypeStr
    filename_naive = "../dat/out/sample.1.2.%s.naive" %scalarTypeStr
    filename_json = "../dat/out/sample.1.2.%s.json" %scalarTypeStr
    
    sample_mID = pylarc.read_row_major_matrix_from_file(os.path.join(os.path.dirname(__file__),filename_rmm)) 
    pylarc.fprint_naive(sample_mID,os.path.join(os.path.dirname(__file__),filename_naive))
    pylarc.fprint_larcMatrixFile(sample_mID,os.path.join(os.path.dirname(__file__),filename_json))
    print("Printing out the nonsquare rrm sample matrix\n")
    pylarc.print_naive(sample_mID)

    print("Testing printing nonzeros to file.\n")
    filename_rmm = "../dat/in/sample.1.3.%s.rmm" %scalarTypeStr
    filename_nonzeros = "../dat/out/sample.1.3.%s.nonzeros" %scalarTypeStr
    filename_json = "../dat/out/sample.1.3.%s.json" %scalarTypeStr
    
    sample_mID = pylarc.read_row_major_matrix_from_file(os.path.join(os.path.dirname(__file__), filename_rmm))
    pylarc.fprint_matrix_nonzeros(sample_mID,os.path.join(os.path.dirname(__file__),filename_nonzeros))
    pylarc.fprint_larcMatrixFile(sample_mID,os.path.join(os.path.dirname(__file__),filename_json))
    print("Here is the matrix we are testing for printing out nonzero values")
    pylarc.print_naive(sample_mID)
    
    print("\n")

    # make CNOT
    print("\nHere is the CNOT (reversible XOR) matrix\n")
    CNOT_arr = list(map(str,[1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0]))
    CNOT_mID = pylarc.row_major_list_to_store(CNOT_arr,level,level,dim_whole)
    pylarc.print_naive(CNOT_mID)

    # Calculate number of matrices created, then print part of matrix store
    num_matrices_made = pylarc.num_matrices_created()
    print("\n%d matrices have been created" %num_matrices_made)
    start = end + 1
    end = num_matrices_made - 1
    filename = "../dat/out/cnot.%s.store" %scalarTypeStr
    pylarc.fprint_store_info_for_matrixID_range(start,end,os.path.join(os.path.dirname(__file__), filename),"Loaded CNOT")

    # build Zero matrices
    print("\nHere is the level 2 zero matrix\n")
    Z2_arr = list(map(str,[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]))
    Z2_mID = pylarc.row_major_list_to_store(Z2_arr,level,level,dim_whole)
    pylarc.print_naive(Z2_mID)

    # build Identity matrices
    print("\nHere is the level 2 identity matrix\n")
    I2_arr = list(map(str,[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]))
    I2_mID = pylarc.row_major_list_to_store(I2_arr,level,level,dim_whole)
    pylarc.print_naive(I2_mID)

    # build a doubly-controlled NOT (base unit for reversible computing)
    print("\nThe 3 bit reversible AND (Toffoli) matrix with target 3rd.\n")
    # get_pID_from_four_sub_pIDs is under construction
    TOFFOLI_mID= pylarc.get_pID_from_four_sub_pIDs(I2_mID,Z2_mID,Z2_mID,CNOT_mID,3,3)
    pylarc.print_naive(TOFFOLI_mID)
    filename = "../dat/out/toffoli.%s.naive" %scalarTypeStr
    pylarc.fprint_naive(TOFFOLI_mID,os.path.join(os.path.dirname(__file__),filename))

    # use CCNOT to build an NAND
    print("\nHere is the 3 bit reversible NAND matrix with target 3rd.\n")
    NOT_matrixID = pylarc.cvar.packedID_NOT;
    I1_matrixID = pylarc.get_identity_pID(1);
    not3_matrixID = pylarc.kronecker_product(I1_matrixID,pylarc.kronecker_product(I1_matrixID, NOT_matrixID));
    nand_from_Toff_matrixID = pylarc.matrix_mult(not3_matrixID,TOFFOLI_mID);
    pylarc.print_naive(nand_from_Toff_matrixID)
    filename = "../dat/out/nandfromtoff.%s.naive" %scalarTypeStr
    pylarc.fprint_naive(nand_from_Toff_matrixID,os.path.join(os.path.dirname(__file__),filename))

    #  PLAYING WITH PREWRITTEN NONSQUARE MATRIX
    filename = "../dat/in/sample.1.2.%s.json" %scalarTypeStr
    print("About to test read %s\n" %filename)
    samp_mID = pylarc.read_larcMatrixFile(os.path.join(os.path.dirname(__file__),filename))
    print("We read in the LARCMatrix file\n")
    pylarc.print_naive(samp_mID)
    
    print("does scalarM1_val print?")
    scalarM1_val = '-1'
    scalarM1_mID = pylarc.get_valID_from_valString(scalarM1_val)
    pylarc.print_naive(scalarM1_mID)
    
    print("testing scalar_mult:")
    samp2_mID = pylarc.scalar_mult(scalarM1_mID,samp_mID)
    pylarc.print_naive(samp2_mID)
    
    print("testing addition:")
    samp3_mID = pylarc.matrix_add(samp_mID,samp2_mID)
    pylarc.print_naive(samp3_mID)
    
    # save input matrixIDs for testing op store hash chains later
    in1_test_sum_mID = samp_mID
    in2_test_sum_mID = samp2_mID
    
    print("testing adjoint:")
    samp3_mID = pylarc.adjoint(samp_mID)
    pylarc.print_naive(samp3_mID)
    adj_mID = samp3_mID
    
    print("testing non-square matrix mult:")
    samp4_mID = pylarc.matrix_mult(samp_mID,samp3_mID)
    pylarc.print_naive(samp4_mID)
    # print("")
    # samp4_mID = pylarc.matrix_mult(samp3_mID,samp_mID)
    # pylarc.print_naive(samp4_mID)

    print("testing kron product:")
    samp4_mID = pylarc.kronecker_product(samp_mID,samp_mID)
    pylarc.print_naive(samp4_mID)
    

    print("testing join:")
    samp4_mID = pylarc.join(samp_mID,samp_mID)
    pylarc.print_naive(samp4_mID)
    print("testing stack:")
    samp4_mID = pylarc.stack(samp_mID,samp_mID)
    pylarc.print_naive(samp4_mID)


    #*  TESTING DELETION
    print("\nPreparing to delete a matrix from the store.\n")
    filename = "../dat/out/temp.%s.json" %scalarTypeStr
    pylarc.fprint_larcMatrixFile(nand_mID, os.path.join(os.path.dirname(__file__),filename))
    pylarc.read_larcMatrixFile(os.path.join(os.path.dirname(__file__),filename))

    # Calculate number of matrices created, then print part of matrix store
    num_matrices_made = pylarc.num_matrices_created()
    print("\n%d matrices have been created" %num_matrices_made)
    print("Previous range printed ended with matrixID %d\n" %end)
    if (end == num_matrices_made-1) :
        print("Nothing new since last matrix store print\n")
    else :
        start = end + 1
        end = num_matrices_made - 1
        filename = "../dat/out/nand.%s.store" %scalarTypeStr
        pylarc.fprint_store_info_for_matrixID_range(start,end,os.path.join(os.path.dirname(__file__),filename),"Loaded NAND")

    # get the hashID and print the hash chain corresponding to a matrix we are about to delete
    hashID = pylarc.hash_pID(nand_mID)
    comment = "hash chain before removal"
    filename = "../dat/out/hashChain.beforeMatrixRemove"
    out_path =  os.path.join(os.path.dirname(__file__),filename)
    pylarc.fprint_nonscalar_hash_chain_info(hashID, out_path, comment)
    
    # Test deletion of a matrix
    print("Testing removal of matrix from the matrix store\n")
    num_matrices_made =  pylarc.num_matrices_created()
    end = num_matrices_made - 1
    filename = "../dat/out/nandYES.%s.store" %scalarTypeStr
    pylarc.fprint_store_info_for_matrixID_range(0,end,os.path.join(os.path.dirname(__file__),filename),"Before Removed NAND")

    pylarc.remove_matrix_from_store(nand_mID)
	
    filename = "../dat/out/nandNO.%s.store" %scalarTypeStr
    filename_json = "../dat/out/temp.%s.json" %scalarTypeStr
    print("\nDeleting the NAND matrix with matrixID", nand_mID,"from store, which had been read from %s\n"  %filename_json)
    pylarc.fprint_store_info_for_matrixID_range(0,end,os.path.join(os.path.dirname(__file__),filename),"Removed NAND")

    comment = "hash chain after removal"
    filename = "../dat/out/hashChain.afterMatrixRemove"
    out_path =  os.path.join(os.path.dirname(__file__),filename)
    pylarc.fprint_nonscalar_hash_chain_info(hashID, out_path, comment)


    pylarc.list_op_names()
    
    # Test op store hash chains after deletion
	
    # print ops store report before deletion of a matrix
    pylarc.op_store_report("stdout")
    
    # print op store hash chain for "SUM"
    op_name = "SUM"
    op_type = pylarc.get_op_type_from_string_name(op_name)
    
    # get hash for a operation record and print the hash chain
    sum_hashID = pylarc.hash_from_op(in1_test_sum_mID,in2_test_sum_mID,op_type)
    if (sum_hashID == -1):
        print("invalid matrixID requested for op hash chain")
    else:
        pylarc.op_hash_chain_info_to_screen(sum_hashID, "op hash chain before deletion")
        
    
    # delete the first input matrix
    pylarc.remove_matrix_from_store(in1_test_sum_mID)
    print("deleting matrix with matrixID", in1_test_sum_mID)
    
    # traverse the op_hash_chain by trying some new sums
    # test1_mID = pylarc.matrix_add(samp4_mID,samp4_mID)
    # test2_mID = pylarc.matrix_add(test1_mID,samp4_mID)
    # test3_mID = pylarc.matrix_add(test2_mID,test1_mID)
    
    # set a hold on a matrix by matrixID to see if it is immune to cleaning
    print("The matrixID of matrix to be held is", adj_mID)
    pylarc.set_hold_matrix(adj_mID)
    
    # clean the matrix store and print it again
    pylarc.clean_matrix_storage()
    filename = "../dat/out/nandNOcleaned.%s.store" %scalarTypeStr
    pylarc.fprint_store_info_for_matrixID_range(0,end,os.path.join(os.path.dirname(__file__),filename),"Removed NAND and cleaned matrix store")

    # clean the op store 
    for hash in range(1<<op_store_exp):
        pylarc.clean_op_hash_chain(hash) 

    # print ops store report after deletion of a matrix
    pylarc.op_store_report("stdout")

    # print same op hash chain again after deleting a matrix, holding a matrix and cleaning
    if (sum_hashID != -1):
        pylarc.op_hash_chain_info_to_screen(sum_hashID, "op hash chain after deletion and cleaning")
    else:
        print("invalid matrixID requested for op hash chain")
    
