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



#*########################################################################
#* test_infoStore2 reads a small file (infoTest1.json) produced by
#* test_infoStore1.py and saved in the current directory, then outputs
#* the metadata found in that LARCMatrix file
#*########################################################################

from __future__ import print_function

import os
import glob
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../../src"))
import pylarc
import numpy as np
from ctypes import *


if __name__ == '__main__':

   #*#################################
   #*     SET THESE PARAMETERS     #*
   #*#################################
    log_psize = 3         #*  problem_size is always power of two!

   #*###############################
   #*    Calculated Parameters   #*
   #*###############################
    problem_size = 2**(log_psize) #* level of the matrices (2^level by 2^level matrices)
    max_level = problem_size      #* the maximum size of matrices in the matrix store

   #*#####################################
   #*    Print baseline usage report   #*
   #*#####################################
    pylarc.memory_and_time_report(0, "stdout")

   #*##################################################################
   #*    LARC  Initialization of Matrix Store and Operation Stores  #*
   #*##################################################################
   #* The routine initialize_larc() does the following:             #*
   #* * creates the matrix and op stores                            #*
   #* * preloads matrix store with: standard scalars and gates,     #*
   #*   and with all zero, identity, and (integer) Hadamard matrices#*
   #*   left to max matrix size                                     #*
   #*##################################################################
   #* SMALL STORES for working on desktop 
    if log_psize <= 3:    
        matrix_exponent = 22
        op_exponent = 19            
        regionbitparam = -1 # default value
        zeroregionbitparam = -1 # default value is 20
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

   #*#################################
   #*  Read test file
   #*#################################
    local_name = "../dat/out/infoTest1.json"
    in_name = os.path.join(os.path.dirname(__file__),local_name)
    inFile_ID = pylarc.read_larcMatrixFile(in_name)


   #*#################################
   #*  Look at info store contents
   #*#################################
    pylarc.list_info_names()

    # info_name = "DATE"
    # info_type = pylarc.get_info_type_from_string_name(info_name)
    # info_data = pylarc.info_get(info_type, inFile_ID)
    # print("info_data = *%s* for info_name = %s\n"  %(info_data,info_name))

    # info_name = "COMPUTER"
    # info_type = pylarc.get_info_type_from_string_name(info_name)
    # info_data = pylarc.info_get(info_type, inFile_ID)
    # print("info_data = *%s* for info_name = %s\n"  %(info_data,info_name))

    max_info_type = pylarc.get_info_type_from_string_name("INVALID_INFO")
    for info_type in range(max_info_type): 
      info_name = pylarc.return_info_name(info_type)
      info_data = pylarc.info_get(info_type,inFile_ID)
      if (info_data != ""):
        print("info_data = *%s* for info_name = %s\n"  %(info_data,info_name))

   #*#################################
   #*  Write copy of test file
   #*#################################
    local_name = "../dat/out/infoTest2.json"
    out_name =  os.path.join(os.path.dirname(__file__),local_name)
    pylarc.fprint_larcMatrixFile(inFile_ID,out_name)
    print("Printing the file %s\n" %out_name)

