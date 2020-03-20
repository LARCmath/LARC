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
    verbose = 1
    pylarc.create_report_thread(1800)
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
    panel = [serial]*4   # alternatively panel=[serial,serial,serial,serial]
    serial_parent = pylarc.get_matID_from_four_subMatIDs(serial,serial,serial,serial,3,3)
    pylarc.print_naive_by_matID(serial_parent)
    filename = "../dat/out/testfile.%s.json" %scalarTypeStr
    pylarc.write_larcMatrix_file_by_matID(serial_parent,os.path.join(os.path.dirname(__file__), filename))
    
    #  PLAYING WITH PREWRITTEN REVERSIBLE NAND LARCMatrix FILE
    print("About to test read LARCMatrix file\n")
    filename = "../dat/in/nand.%s.json" %scalarTypeStr
    nand_matrixID = pylarc.read_larcMatrix_file_return_matID(os.path.join(os.path.dirname(__file__),filename))
    print("We read in the LARCMatrix nand file\n")
    pylarc.print_naive_by_matID(nand_matrixID)
    print("\n")

    # TESTING READING AND WRITING OF MATRICES
    print("Testing reading row major matrix format and writing files in LARCMatrix and naive format.\n")
    filename_rmm = "../dat/in/sample.1.1.%s.rmm" %scalarTypeStr
    filename_naive = "../dat/out/sample.1.1.%s.naive" %scalarTypeStr
    filename_json = "../dat/out/sample.1.1.%s.json" %scalarTypeStr

    print(os.path.join(os.path.dirname(__file__),filename_rmm))
    
    sample_matrixID = pylarc.read_row_major_matrix_from_file_matrixID(os.path.join(os.path.dirname(__file__),filename_rmm)) 
    pylarc.write_naive_by_matID(sample_matrixID,os.path.join(os.path.dirname(__file__),filename_naive))
    pylarc.write_larcMatrix_file_by_matID(sample_matrixID,os.path.join(os.path.dirname(__file__),filename_json))
    print("Printing out the rrm sample matrix in naive format to screen\n")
    pylarc.print_naive_by_matID(sample_matrixID)

    print("Testing reading row major nonsquare matrices and writing files in LARCMatrix and naive format.\n")
    filename_rmm = "../dat/in/sample.1.2.%s.rmm" %scalarTypeStr
    filename_naive = "../dat/out/sample.1.2.%s.naive" %scalarTypeStr
    filename_json = "../dat/out/sample.1.2.%s.json" %scalarTypeStr
    
    sample_matrixID = pylarc.read_row_major_matrix_from_file_matrixID(os.path.join(os.path.dirname(__file__),filename_rmm)) 
    pylarc.write_naive_by_matID(sample_matrixID,os.path.join(os.path.dirname(__file__),filename_naive))
    pylarc.write_larcMatrix_file_by_matID(sample_matrixID,os.path.join(os.path.dirname(__file__),filename_json))
    print("Printing out the nonsquare rrm sample matrix\n")
    pylarc.print_naive_by_matID(sample_matrixID)

    #  PLAYING WITH PREWRITTEN NONSQUARE MATRIX
    filename = "../dat/in/sample.1.2.%s.json" %scalarTypeStr
    print("About to test read %s\n" %filename)
    samp_matrixID = pylarc.read_larcMatrix_file_return_matID(os.path.join(os.path.dirname(__file__),filename))

    print("We read in the LARCMatrix file\n")
    pylarc.print_naive_by_matID(samp_matrixID)
    print("\n")

    print("Testing printing nonzeros to file.\n")
    filename_rmm = "../dat/in/sample.1.3.%s.rmm" %scalarTypeStr
    filename_nonzeros = "../dat/out/sample.1.3.%s.nonzeros" %scalarTypeStr
    filename_json = "../dat/out/sample.1.3.%s.json" %scalarTypeStr
    
    sample_matrixID = pylarc.read_row_major_matrix_from_file_matrixID(os.path.join(os.path.dirname(__file__), filename_rmm))
    pylarc.write_matrix_nonzeros_by_matID(sample_matrixID,os.path.join(os.path.dirname(__file__),filename_nonzeros))
    pylarc.write_larcMatrix_file_by_matID(sample_matrixID,os.path.join(os.path.dirname(__file__),filename_json))
    print("Here is the matrix we are testing for printing out nonzero values")
    pylarc.print_naive_by_matID(sample_matrixID)
    
    print("\n")

    # make CNOT
    print("\nHere is the CNOT (reversible XOR) matrix\n")

    #  CNOT_arr = pylarc.buildArray([1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0])
    CNOT_arr = list(map(str,[1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0]))
    CNOT_matrixID = pylarc.row_major_list_to_store_matrixID(CNOT_arr,level,level,dim_whole)
    pylarc.print_naive_by_matID(CNOT_matrixID)

    # Calculate number of matrices created, then print part of matrix store
    num_matrices_made = pylarc.num_matrices_created()
    print("\n%d matrices have been created" %num_matrices_made)
    start = end + 1
    end = num_matrices_made - 1
    filename = "../dat/out/cnot.%s.store" %scalarTypeStr
    pylarc.matrix_store_info_to_file(start,end,os.path.join(os.path.dirname(__file__), filename),"Loaded CNOT")

    # build Zero matrices
    print("\nHere is the level 2 zero matrix\n")
    # Z2_arr = pylarc.buildArray([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    Z2_arr = list(map(str,[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]))
    Z2_matrixID = pylarc.row_major_list_to_store_matrixID(Z2_arr,level,level,dim_whole)
    pylarc.print_naive_by_matID(Z2_matrixID)

    # build Identity matrices
    print("\nHere is the level 2 identity matrix\n")
    # I2_arr = pylarc.buildArray([1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1])
    I2_arr = list(map(str,[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]))
    I2_matrixID = pylarc.row_major_list_to_store_matrixID(I2_arr,level,level,dim_whole)
    pylarc.print_naive_by_matID(I2_matrixID)

    # build a Toffoli (base unit for reversible computing)
    print("\nThe 3 bit reversible AND (Toffoli) matrix with target 3rd.\n")
    # get_matID_from_four_subMatIDs is under construction
    TOFFOLI_matrixID= pylarc.get_matID_from_four_subMatIDs(I2_matrixID,Z2_matrixID,Z2_matrixID,CNOT_matrixID,3,3)
    pylarc.print_naive_by_matID(TOFFOLI_matrixID)
    filename = "../dat/out/toffoli.%s.naive" %scalarTypeStr
    pylarc.write_naive_by_matID(TOFFOLI_matrixID,os.path.join(os.path.dirname(__file__),filename))

    # use Toffoli to build a NAND
    print("\nHere is the 3 bit reversible NAND matrix with target 3rd.\n")
    # build 8x8 matrix I tensor I tensor X
    NOT_matrixID = pylarc.cvar.matID_NOT;
    I1_matrixID = pylarc.get_identity_matrixID(1);
    not3_matrixID = pylarc.kronecker_product_matrixID(I1_matrixID,pylarc.kronecker_product_matrixID(I1_matrixID, NOT_matrixID));
    nand_from_Toff_matrixID = pylarc.matrix_mult_matrixID(not3_matrixID,TOFFOLI_matrixID);
    pylarc.print_naive_by_matID(nand_from_Toff_matrixID)
    filename = "../dat/out/nandfromtoff.%s.naive" %scalarTypeStr
    pylarc.write_naive_by_matID(nand_from_Toff_matrixID,os.path.join(os.path.dirname(__file__),filename))
    # build from panels
    print("\nHere is the 3 bit reversible NAND matrix with target 3rd, constructed a different way.\n")
    Z1_matrixID = pylarc.get_zero_matrixID(1,1);
    topleft_matrixID = pylarc.get_matID_from_four_subMatIDs(NOT_matrixID,Z1_matrixID,Z1_matrixID,NOT_matrixID,2,2);
    botright_matrixID = pylarc.get_matID_from_four_subMatIDs(NOT_matrixID,Z1_matrixID,Z1_matrixID,I1_matrixID,2,2);
    nand_from_panels_matrixID = pylarc.get_matID_from_four_subMatIDs(topleft_matrixID,Z2_matrixID,Z2_matrixID,botright_matrixID,3,3);
    pylarc.print_naive_by_matID(nand_from_panels_matrixID)
    filename = "../dat/out/nandfrompanels.%s.naive" %scalarTypeStr
    pylarc.write_naive_by_matID(nand_from_panels_matrixID,os.path.join(os.path.dirname(__file__),filename))
    # compare with other nands
    if (nand_from_Toff_matrixID != nand_matrixID):
        print("problem with building NAND from Toffoli")
    if (nand_from_panels_matrixID != nand_matrixID):
        print("problem with building NAND from panels")


    # read the (1,2) entry from the nand matrix        
    scalarID = pylarc.get_valID_from_matID_and_coords(nand_matrixID, 1, 2);
    print("matrix ID of scalar in position (1,2) of NAND is %d" %scalarID)
    print("The value of scalar in position (1,2) of NAND is")
    pylarc.print_naive_by_matID(scalarID)

    # read the (0,1) entry from the nand matrix        
    scalarID = pylarc.get_valID_from_matID_and_coords(nand_matrixID, 0, 1);
    print("matrix ID of scalar in position (0,1) of NAND is %d" %scalarID)
    print("The value of scalar in position (0,1) of NAND is")
    pylarc.print_naive_by_matID(scalarID)
    print("The same value accessed directly using indices and get_valString_from_matID_and_coords")
    myStr = pylarc.get_valString_from_matID_and_coords(nand_matrixID,0,1);
    print("%s" %myStr);
