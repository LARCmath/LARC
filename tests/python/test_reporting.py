#!/usr/bin/env python3

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
import glob
import sys
import random # needed for toeplitz
import numpy as np # needed for toeplitz
sys.path.append(os.path.join(os.path.dirname(__file__),"../../src"))
import pylarc
from ctypes import *

if __name__ == '__main__':

    #*###############################
    #*   SET THESE PARAMETERS      #*
    #*###############################
    max_level = 8          #*  problem_size is always power of two!

    #*#####################################
    #*    Print baseline usage report    #*
    #*#####################################
    pylarc.memory_and_time_report(0, "stdout")


    #*##################################################################
    #*    LARC  Initialization of Matrix Store and Operation Stores   #*
    #*##################################################################
    #* The routine initialize_larc() does the following:              #*
    #* * creates the matrix and op stores                             #*
    #* * preloads matrix store with: standard scalars and gates,      #*
    #*   and with all zero, identity, and (integer) Hadamard matrices #*
    #*   left to max matrix size                                      #*
    #*##################################################################

    #*##################################################################
    #*    Testing to see what LSH function parameters we should use.  #*
    #*    The zeroregionbitparam values that work for F_3 are:        #*
    #*      -z 53 and smaller                                         #*
    #*    The zeroregionbitparam values that fail for F_3 are:        #*
    #*      -z 54 and larger                                          #*
    #*                                                                #*
    #*    The regionbitparam value has little effect, it              #*
    #*    works in ranges -s 10 to -s 1000                            #*
    #*                                                                #*
    #*    The version of LARC post-April 2020 will ignore -z values   #*
    #*    unless they are less than -s values.                        #*
    #*##################################################################



    #* SMALL STORES for working on desktop
    if max_level <= 8:    
        matrix_exponent = 22
        op_exponent = 19   


        zeroregionbitparam = 54
        regionbitparam = 1000

        #*####################################################
        #*  Sample failure values for LARC LSH parameters   #*
        #*  are to set         regionbitparam = 60            #*
        #*  and                zeroregionbitparam = 60      #*
        #*####################################################
        #*  Default values for LARC LSH parameters are      #*
        #*  both equal to LDBL_MANT_DIG -2 when scalarType is Real #*
        #*  regionbitparam = -1    #                        #*
        #*  zeroregionbitparam = -1 # default is regionbitparam #*
        #*  NOTE:  testing shows -z 47 will work            #*
        #*####################################################

        verbose = 1
        pylarc.initialize_larc(matrix_exponent,op_exponent,max_level,regionbitparam,zeroregionbitparam,verbose)
        pylarc.create_report_thread(180)
        print_naive = 0
        print_nonzeros = 0
        print("Problem size is small enough to run on desktop")
        if print_naive:
           print("  will print files of naive matrices")
        else: 
           print("  not printing files of naive matrices")
        if print_nonzeros:
           print("  will print files of nonzero matrices\n")
        else: 
           print("  not printing files of nonzero matrices\n")
    #* LARGE STORES
    else:      
        #* matrix_exponent = 26
        #* op_exponent = 24
        matrix_exponent = 30
        op_exponent = 31   
        regionbitparam = -1 # default value
        zeroregionbitparam = -1 # default value
 
        verbose = 1
        pylarc.initialize_larc(matrix_exponent,op_exponent,max_level,regionbitparam,zeroregionbitparam,verbose)
        pylarc.create_report_thread(3600)   # once per hour 3600
        print_naive = 0      
        print_nonzeros = 0
        print("Problem size is NOT small enough to run on desktop")
        if print_naive:
           print("  WARNING: will try to print files of naive matrices!!!")
        else: 
           print("  not printing files of naive matrices")
        if print_nonzeros:
           print("  WARNING: will print files of nonzero matrices!!!\n")
        else: 
           print("  not printing files of nonzero matrices\n")

    print("Finished creating LARC matrix and op stores and loading basic matrices.\n")
    print("stopHogging check to see if program is too large, to occur once every 10 minutes.\n")


    
    # Define string for using in formating filenames
    scalarTypeStr = pylarc.cvar.scalarTypeStr


    #*  START OF CODE STOLEN FROM toeplitz.py

    # build array in C from Python list of scalars
    # print("Using row_major_list_to_store on data entered from python\n")

    # parameters for entering the python array into the store
    level = 8
    dim_whole = 2**level

    if scalarTypeStr in ("Real", "MPReal", "MPRational"):
        randVals = [ random.random() for i in range(2*dim_whole-1)]
    elif scalarTypeStr in ("Complex", "MPComplex", "MPRatComplex"):
        randVals = [ np.complex(random.random(),random.random()) for i in range(2*dim_whole-1)]
    elif scalarTypeStr in ("Integer", "MPInteger"):
        randVals = [ random.randrange(0,10001) for i in range(2*dim_whole-1)]
    else:
        raise Exception('Do not know how to build matrix for type %s.' %scalarTypeStr)

    # create toeplitz matrix from the 2*dim_whole-1 random numbers
    a = []
    b = randVals
    # print("b: ", b)
    for i in range(dim_whole):
        a.append(list(b[dim_whole-i-1:2*dim_whole-i-1]))
    # print("a: ", a)
    amat = np.matrix(a)
    top_ID = pylarc.add_numpy_matrix_to_matrix_store(amat)

    # filename_json = "../dat/out/toeplitz.lev%d.%s.json" %(level,scalarTypeStr)
    # pylarc.print_naive(serial)
    # print("\n")
    # pylarc.fprint_larcMatrixFile(serial,os.path.join(os.path.dirname(__file__),filename_json))

    #*  END OF CODE STOLEN FROM toeplitz.py

    # some operations to put in the store
    A_ID = pylarc.matrix_add(top_ID,top_ID)
    B_ID = pylarc.matrix_mult(A_ID,top_ID)
    



    # #*###########################
    # # test reporting
    # #*###########################
    # pylarc.matrix_store_report("../dat/out/temp.mat_report")
    # pylarc.op_store_report("../dat/out/temp.op_report")
    # pylarc.memory_and_time_report(0,"../dat/out/temp.rusage0_report")
    # pylarc.memory_and_time_report(1,"../dat/out/temp.rusage1_report")

    print("\n")
    pylarc.matrix_store_report("stdout")
    print("\n")
    pylarc.op_store_report("stdout")
    print("\n")
    pylarc.memory_and_time_report(0,"stdout")
    print("\n")
    pylarc.memory_and_time_report(1,"stdout")
