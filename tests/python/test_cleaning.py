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
	# This version references matrices by matrixIDs instead of pointers
	
    print "This code tests some basic matrix building and reading routines\n"

    # initialize larc
    mat_store_exp = 15
    op_store_exp = 5
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

    # Calculate number of matrices created, then print part of matrix store
    num_matrices_made = pylarc.num_matrices_created()
    print "\n%d matrices have been created" %num_matrices_made
    end = num_matrices_made - 1
    filename = "../dat/out/preload.%s.store" %scalarType
    pylarc.matrix_store_info_to_file(0,end,os.path.join(os.path.dirname(__file__),filename),"After preload with parameters: 26, 24, 10.")

    #  PLAYING WITH PREWRITTEN NONSQUARE MATRIX
    filename = "../dat/in/sample.1.2.%s.json" %scalarType
    print "About to test read %s\n" %filename
    samp_mID = pylarc.matrix_read_json_file_matrixID(os.path.join(os.path.dirname(__file__),filename))

    print "We read in the json file\n"
    pylarc.print_matrix_naive_by_matrixID(samp_mID)
    print "\n"

    # build array in C from Python list of scalars
    print "Using row_major_list_to_store on data entered from python\n"

    # create a matrix in python
    if pylarc.cvar.scalarTypeDef == 'i':
        a = np.matrix([[1, 3, 5, 6],
                       [8, 6, 3, 1],
                       [-9, 11, 13, 15],
                       [16, 13, 12, 10]])
    elif pylarc.cvar.scalarTypeDef == 'c':
        a = np.matrix([[1+2j, 3+4j, 5+6j, 7+8j],
                       [8+7j, 6+5j, 3+4j, 1+2j],
                       [9+10j, 11+12j, 13+14j, 15+16j],
                       [16+15j, 14+13j, 12+11j, 10+9j]])
    elif pylarc.cvar.scalarTypeDef == 'r':
        a = np.matrix([[1, 3, .5, 6],
                       [8, 6, 3, .1],
                       [-9, 11, 13, 1.5],
                       [16, 13, 12, 10]])
    else:
        raise Exception('Do not know how to build matrix for type %s.'%(pylarc.cvar.scalarTypeDef,))

    # turn the matrix into an array by reading off each row in turn (row major format)
    alist = a.reshape(-1).tolist()[0]
    arr = pylarc.buildArray(alist)
    print 'arr:', pylarc.str_scalarTypeArray(arr, len(alist))

    # parameters for entering the python array into the store
    level = 2
    dim_whole = 2**level

    # creating or finding the matrix associated with the array
    mID = pylarc.row_major_list_to_store_matrixID(arr, level, level, dim_whole)
    pylarc.print_matrix_naive_by_matrixID(mID)
    print "\n"

    # make a parent matrix from four copies of the matrixID matrix
    print "Creating matrix from matrix_get_matrixID_from_panel on panel input and writing json file\n"
    #  **matrix_get_matrixID_from_scalar and matrix_get_matrixID_from_panel are under construction**
    panel = [mID]*4   # alternatively panel=[mID,mID,mID,mID]
    mID_parent = pylarc.matrix_get_matrixID_from_panel(mID,mID,mID,mID,3,3)
    pylarc.print_matrix_naive_by_matrixID(mID_parent)
    filename = "../dat/out/testfile.%s.json" %scalarType
    pylarc.matrix_write_json_file_matrixID(mID_parent,os.path.join(os.path.dirname(__file__), filename))
    
    #  PLAYING WITH PREWRITTEN REVERSIBLE NAND JSON FILE
    print "About to test read json file\n"
    filename = "../dat/in/nand.%s.json" %scalarType
    nand_mID = pylarc.matrix_read_json_file_matrixID(os.path.join(os.path.dirname(__file__),filename))
    print "We read in the json nand file\n"
    pylarc.print_matrix_naive_by_matrixID(nand_mID)
    print "\n"

    # TESTING READING AND WRITING OF MATRICES
    print "Testing reading row major matrix format and writing files in json and naive format.\n"
    filename_rmm = "../dat/in/sample.1.1.%s.rmm" %scalarType
    filename_naive = "../dat/out/sample.1.1.%s.naive" %scalarType
    filename_json = "../dat/out/sample.1.1.%s.json" %scalarType
    
    sample_mID = pylarc.read_row_major_matrix_from_file_matrixID(os.path.join(os.path.dirname(__file__),filename_rmm)) 
    pylarc.print_matrix_to_file_naive_by_matrixID(sample_mID,os.path.join(os.path.dirname(__file__),filename_naive))
    pylarc.matrix_write_json_file_matrixID(sample_mID,os.path.join(os.path.dirname(__file__),filename_json))
    print "Printing out the rrm sample matrix in naive format to screen\n"
    pylarc.print_matrix_naive_by_matrixID(sample_mID)

    print "Testing reading row major nonsquare matrices and writing files in json and naive format.\n"
    filename_rmm = "../dat/in/sample.1.2.%s.rmm" %scalarType
    filename_naive = "../dat/out/sample.1.2.%s.naive" %scalarType
    filename_json = "../dat/out/sample.1.2.%s.json" %scalarType
    
    sample_mID = pylarc.read_row_major_matrix_from_file_matrixID(os.path.join(os.path.dirname(__file__),filename_rmm)) 
    pylarc.print_matrix_to_file_naive_by_matrixID(sample_mID,os.path.join(os.path.dirname(__file__),filename_naive))
    pylarc.matrix_write_json_file_matrixID(sample_mID,os.path.join(os.path.dirname(__file__),filename_json))
    print "Printing out the nonsquare rrm sample matrix\n"
    pylarc.print_matrix_naive_by_matrixID(sample_mID)

    print "Testing printing nonzeros to file.\n"
    filename_rmm = "../dat/in/sample.1.3.%s.rmm" %scalarType
    filename_nonzeros = "../dat/out/sample.1.3.%s.nonzeros" %scalarType
    filename_json = "../dat/out/sample.1.3.%s.json" %scalarType
    
    sample_mID = pylarc.read_row_major_matrix_from_file_matrixID(os.path.join(os.path.dirname(__file__), filename_rmm))
    pylarc.print_matrix_nonzeros_to_file_by_matrixID(sample_mID,os.path.join(os.path.dirname(__file__),filename_nonzeros))
    pylarc.matrix_write_json_file_matrixID(sample_mID,os.path.join(os.path.dirname(__file__),filename_json))
    print "Here is the matrix we are testing for printing out nonzero values"
    pylarc.print_matrix_naive_by_matrixID(sample_mID)
    
    print "\n"

    # make CNOT
    print("\nHere is the CNOT (reversible XOR) matrix\n")
    CNOT_arr = pylarc.buildArray([1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0])
    CNOT_mID = pylarc.row_major_list_to_store_matrixID(CNOT_arr,level,level,dim_whole)
    pylarc.print_matrix_naive_by_matrixID(CNOT_mID)

    # Calculate number of matrices created, then print part of matrix store
    num_matrices_made = pylarc.num_matrices_created()
    print "\n%d matrices have been created" %num_matrices_made
    start = end + 1
    end = num_matrices_made - 1
    filename = "../dat/out/cnot.%s.store" %scalarType
    pylarc.matrix_store_info_to_file(start,end,os.path.join(os.path.dirname(__file__), filename),"Loaded CNOT")

    # build Zero matrices
    print("\nHere is the level 2 zero matrix\n")
    Z2_arr = pylarc.buildArray([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    Z2_mID = pylarc.row_major_list_to_store_matrixID(Z2_arr,level,level,dim_whole)
    pylarc.print_matrix_naive_by_matrixID(Z2_mID)

    # build Identity matrices
    print("\nHere is the level 2 identity matrix\n")
    I2_arr = pylarc.buildArray([1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1])
    I2_mID = pylarc.row_major_list_to_store_matrixID(I2_arr,level,level,dim_whole)
    pylarc.print_matrix_naive_by_matrixID(I2_mID)

    # build a Toffoli (base unit for reversible computing)
    print("\nThe 3 bit reversible AND (Toffoli) matrix with target 3rd.\n")
    # matrix_get_matrixID_from_panel is under construction
    TOFFOLI_mID= pylarc.matrix_get_matrixID_from_panel(I2_mID,Z2_mID,Z2_mID,CNOT_mID,3,3)
    pylarc.print_matrix_naive_by_matrixID(TOFFOLI_mID)
    filename = "../dat/out/toffoli.%s.naive" %scalarType
    pylarc.print_matrix_to_file_naive_by_matrixID(TOFFOLI_mID,os.path.join(os.path.dirname(__file__),filename))

    # use Toffoli to build an NAND
    print("\nHere is the 3 bit reversible NAND matrix with target 3rd.\n")
    NOT_matrixID = pylarc.cvar.matID_NOT;
    I1_matrixID = pylarc.get_identity_matrixID(1);
    not3_matrixID = pylarc.kronecker_product_matrixID(I1_matrixID,pylarc.kronecker_product_matrixID(I1_matrixID, NOT_matrixID));
    nand_from_Toff_matrixID = pylarc.matrix_mult_matrixID(not3_matrixID,TOFFOLI_mID);
    pylarc.print_matrix_naive_by_matrixID(nand_from_Toff_matrixID)
    filename = "../dat/out/nandfromtoff.%s.naive" %scalarType
    pylarc.print_matrix_to_file_naive_by_matrixID(nand_from_Toff_matrixID,os.path.join(os.path.dirname(__file__),filename))

    #  PLAYING WITH PREWRITTEN NONSQUARE MATRIX
    filename = "../dat/in/sample.1.2.%s.json" %scalarType
    print "About to test read %s\n" %filename
    samp_mID = pylarc.matrix_read_json_file_matrixID(os.path.join(os.path.dirname(__file__),filename))
    print "We read in the json file\n"
    pylarc.print_matrix_naive_by_matrixID(samp_mID)
    
    print "does scalarM1_val print?"
    scalarM1_val = -1
    scalarM1_mID = pylarc.matrix_get_matrixID_from_scalar(scalarM1_val)
    pylarc.print_matrix_naive_by_matrixID(scalarM1_mID)
    
    print "testing scalar_mult:"
    samp2_mID = pylarc.scalar_mult_matrixID(scalarM1_mID,samp_mID)
    pylarc.print_matrix_naive_by_matrixID(samp2_mID)
    
    print "testing addition:"
    samp3_mID = pylarc.matrix_add_matrixID(samp_mID,samp2_mID)
    pylarc.print_matrix_naive_by_matrixID(samp3_mID)
    
    # save input matrixIDs for testing op store hash chains later
    in1_test_sum_mID = samp_mID
    in2_test_sum_mID = samp2_mID
    
    print "testing adjoint:"
    samp3_mID = pylarc.matrix_adjoint_matrixID(samp_mID)
    pylarc.print_matrix_naive_by_matrixID(samp3_mID)
    adj_mID = samp3_mID
    
    print "testing non-square matrix mult:"
    samp4_mID = pylarc.matrix_mult_matrixID(samp3_mID,samp_mID)
    pylarc.print_matrix_naive_by_matrixID(samp4_mID)
    print ""
    samp4_mID = pylarc.matrix_mult_matrixID(samp_mID,samp3_mID)
    pylarc.print_matrix_naive_by_matrixID(samp4_mID)
    print "testing kron product:"
    samp4_mID = pylarc.kronecker_product_matrixID(samp_mID,samp_mID)
    pylarc.print_matrix_naive_by_matrixID(samp4_mID)
    

    print "testing join:"
    samp4_mID = pylarc.join_matrixID(samp_mID,samp_mID)
    pylarc.print_matrix_naive_by_matrixID(samp4_mID)
    print "testing stack:"
    samp4_mID = pylarc.stack_matrixID(samp_mID,samp_mID)
    pylarc.print_matrix_naive_by_matrixID(samp4_mID)


    ##  TESTING DELETION
    print("\nPreparing to delete a matrix from the store.\n")
    filename = "../dat/out/temp.%s.json" %scalarType
    pylarc.matrix_write_json_file_matrixID(nand_mID, os.path.join(os.path.dirname(__file__),filename))
    pylarc.matrix_read_json_file_matrixID(os.path.join(os.path.dirname(__file__),filename))

    # Calculate number of matrices created, then print part of matrix store
    num_matrices_made = pylarc.num_matrices_created()
    print "\n%d matrices have been created" %num_matrices_made
    print "Previous range printed ended with matrixID %d\n" %end
    if (end == num_matrices_made-1) :
        print "Nothing new since last matrix store print\n"
    else :
        start = end + 1
        end = num_matrices_made - 1
        filename = "../dat/out/nand.%s.store" %scalarType
        pylarc.matrix_store_info_to_file(start,end,os.path.join(os.path.dirname(__file__),filename),"Loaded NAND")

    # get the hashID and print the hash chain corresponding to a matrix we are about to delete
    hashID = pylarc.matrix_hashID_from_matrixID(nand_mID)
    comment = "hash chain before removal"
    filename = "../dat/out/hashChain.beforeMatrixRemove"
    out_path =  os.path.join(os.path.dirname(__file__),filename)
    # os.path.dirname(__file__) =  /.ccs/u01/jszito/LARC/tests/python
    # os.path.join("/.ccs/u01/jszito/LARC/tests/python","../dat/out/hashChain.beforeMatrixRemove")
    # out_path = "/.ccs/u01/jszito/LARC/tests/dat/out/hashChain.beforeMatrixRemove"
    pylarc.matrix_hash_chain_info_to_file(hashID, out_path, comment)
    
    # Test deletion of a matrix
    print "Testing removal of matrix from the matrix store\n"
    num_matrices_made =  pylarc.num_matrices_created()
    end = num_matrices_made - 1
    filename = "../dat/out/nandYES.%s.store" %scalarType
    pylarc.matrix_store_info_to_file(0,end,os.path.join(os.path.dirname(__file__),filename),"Before Removed NAND")

    pylarc.remove_matrix_from_mat_store_by_matrixID(nand_mID)
	
    filename = "../dat/out/nandNO.%s.store" %scalarType
    filename_json = "../dat/out/temp.%s.json" %scalarType
    print "\nDeleting the NAND matrix with matrixID", nand_mID,"from store, which had been read from %s\n"  %filename_json
    pylarc.matrix_store_info_to_file(0,end,os.path.join(os.path.dirname(__file__),filename),"Removed NAND")

    comment = "hash chain after removal"
    filename = "../dat/out/hashChain.afterMatrixRemove"
    out_path =  os.path.join(os.path.dirname(__file__),filename)
    pylarc.matrix_hash_chain_info_to_file(hashID, out_path, comment)


    pylarc.list_op_names()
    
    # Test op store hash chains after deletion
	
    # print ops store report before deletion of a matrix
    pylarc.op_store_report("stdout")
    
    # print op store hash chain for "SUM"
    op_name = "SUM"
    
    # get hash for a operation record and print the hash chain
    sum_hashID = pylarc.op_hashID_by_matrixIDs(in1_test_sum_mID,in2_test_sum_mID,op_name)
    if (sum_hashID == -1):
	print "invalid matrixID requested for op hash chain"
    else:
	pylarc.op_hash_chain_info_to_screen(sum_hashID, "op hash chain before deletion")
        
    
    # delete the first input matrix
    pylarc.remove_matrix_from_mat_store_by_matrixID(in1_test_sum_mID)
    print "deleting matrix with matrixID", in1_test_sum_mID
    
    # traverse the op_hash_chain by trying some new sums
    # test1_mID = pylarc.matrix_add_matrixID(samp4_mID,samp4_mID)
    # test2_mID = pylarc.matrix_add_matrixID(test1_mID,samp4_mID)
    # test3_mID = pylarc.matrix_add_matrixID(test2_mID,test1_mID)
    
    # set a hold on a matrix by matrixID to see if it is immune to cleaning
    print "The matrixID of matrix to be held is", adj_mID
    pylarc.set_hold_matrix_from_matrixID(adj_mID)
    
    # clean the matrix store and print it again
    pylarc.clean_matrix_store()
    filename = "../dat/out/nandNOcleaned.%s.store" %scalarType
    pylarc.matrix_store_info_to_file(0,end,os.path.join(os.path.dirname(__file__),filename),"Removed NAND and cleaned matrix store")

    # clean the op store 
    for hash in range(1<<op_store_exp):
		pylarc.clean_op_hash_chain(hash) 

    # print ops store report after deletion of a matrix
    pylarc.op_store_report("stdout")

    # print same op hash chain again after deleting a matrix, holding a matrix and cleaning
    if (sum_hashID != -1):
        pylarc.op_hash_chain_info_to_screen(sum_hashID, "op hash chain after deletion and cleaning")
    else:
        print "invalid matrixID requested for op hash chain"
    
