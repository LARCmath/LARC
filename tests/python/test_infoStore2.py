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


##########################################################################
## test_infoStore2 reads a small file (infoTest1.json) produced by
## test_infoStore1.py and saved in the current directory, then outputs
## the metadata found in that JSON file
##########################################################################

import os
import glob
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../../src"))
import pylarc
import numpy as np
from ctypes import *


if __name__ == '__main__':

    ###################################
    ##     SET THESE PARAMETERS      ##
    ###################################
    log_psize = 3          ##  problem_size is always power of two!

    #################################
    ##    Calculated Parameters    ##
    #################################
    problem_size = 2**(log_psize)  ## level of the matrices (2^level by 2^level matrices)
    max_level = problem_size       ## the maximum size of matrices in the matrix store

    #######################################
    ##    Print baseline usage report    ##
    #######################################
    pylarc.rusage_report(0, "stdout")

    ####################################################################
    ##    LARCt Initialization of Matrix Store and Operation Stores   ##
    ####################################################################
    ## The routine initialize_larc() does the following:              ##
    ## * creates the matrix and op stores                             ##
    ## * preloads matrix store with: standard scalars and gates,      ##
    ##   and with all zero, identity, and (integer) Hadamard matrices ##
    ##   left to max matrix size                                      ##
    ####################################################################
    ## SMALL STORES for working on desktop 
    if log_psize <= 3:    
        matrix_exponent = 22
        op_exponent = 19            
        rnd_sig_bits = -1 # default value
        trunc_to_zero_bits = -1 # default value is 20
        pylarc.initialize_larc(matrix_exponent,op_exponent,max_level,rnd_sig_bits,trunc_to_zero_bits)
        pylarc.create_report_thread(180)
        print_naive = 0
        print_nonzeros = 0
        print "Problem size is small enough to run on desktop"
        if print_naive:
           print "  will print files of naive matrices"
        else: 
           print "  not printing files of naive matrices"
        if print_nonzeros:
           print "  will print files of nonzero matrices\n"
        else: 
           print "  not printing files of nonzero matrices\n"
    ## LARGE STORES for cs1l,cs4l,cs9l
    else:      
        ## matrix_exponent = 26
        ## op_exponent = 24
        matrix_exponent = 30
        op_exponent = 31   
        rnd_sig_bits = -1 # default value
        trunc_to_zero_bits = -1 # default value
        ## trunc_to_zero_bits = 20 # truncate to zero if value is less than 2**(-threshold)
        ## trunc_to_zero_bits = 16 # truncate to zero if value is less than 2**(-threshold)
        pylarc.initialize_larc(matrix_exponent,op_exponent,max_level,rnd_sig_bits,trunc_to_zero_bits)
        pylarc.create_report_thread(3600)   # once per hour 3600
        print_naive = 0      
        print_nonzeros = 0
        print "Problem size is NOT small enough to run on desktop"
        if print_naive:
           print "  WARNING: will try to print files of naive matrices!!!"
        else: 
           print "  not printing files of naive matrices"
        if print_nonzeros:
           print "  WARNING: will print files of nonzero matrices!!!\n"
        else: 
           print "  not printing files of nonzero matrices\n"

    print "Finished creating LARC matrix and op stores and loading basic matrices.\n"
    print "Seppuku check to see if program is to large to occur once every 10 minutes.\n"

    ###################################
    ##  Read test file
    ###################################
    local_name = "../dat/out/infoTest1.json"
    in_name = os.path.join(os.path.dirname(__file__),local_name)
    inFile_ID = pylarc.matrix_read_json_file_matrixID(in_name)


    ###################################
    ##  Look at info store contents
    ###################################
    pylarc.list_info_names()

    # info_name = "DATE"
    # info_type = pylarc.get_info_type_from_string_name(info_name)
    # info_data = pylarc.info_get(info_type, inFile_ID)
    # print "info_data = *%s* for info_name = %s\n"  %(info_data,info_name)

    # info_name = "COMPUTER"
    # info_type = pylarc.get_info_type_from_string_name(info_name)
    # info_data = pylarc.info_get(info_type, inFile_ID)
    # print "info_data = *%s* for info_name = %s\n"  %(info_data,info_name)

    max_info_type = pylarc.get_info_type_from_string_name("INVALID_INFO")
    for info_type in range(max_info_type): 
      info_name = pylarc.return_info_name(info_type)
      info_data = pylarc.info_get(info_type,inFile_ID)
      if (info_data != ""):
        print "info_data = *%s* for info_name = %s\n"  %(info_data,info_name)

    ###################################
    ##  Write copy of test file
    ###################################
    local_name = "../dat/out/infoTest2.json"
    out_name =  os.path.join(os.path.dirname(__file__),local_name)
    pylarc.matrix_write_json_file_matrixID(inFile_ID,out_name)
    print "Printing the file %s\n" %out_name

