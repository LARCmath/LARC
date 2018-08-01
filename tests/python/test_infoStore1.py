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
## test_infoStore1 writes a small JSON file with metadata to the current
## directory
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

    # print "The globbing constants are: SIGHASH %d, ZEROBITTHRESH %d" %(cvar.py_SIGHASH, cvar.py_ZEROBITTHRESH)

   
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

    #############################
    # create a matrix in python #
    #############################
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

    #########################################################
    # turn the matrix into an array by reading off each row #
    # in turn (row major format)                            #
    #########################################################
    alist = a.reshape(-1).tolist()[0]
    arr = pylarc.buildArray(alist)
    print 'arr:', pylarc.str_scalarTypeArray(arr, len(alist))

    # parameters for entering the python array into the store
    level = 2
    dim_whole = 2**level

    # creating or finding the matrix associated with the array
    serial = pylarc.row_major_list_to_store_matrixID(arr, level, level, dim_whole)
    #########################
    # add metadata to store #
    #########################
    info_name = "DATE"
    info_type = pylarc.get_info_type_from_string_name(info_name)
    info_data = "20170707"
    print "Loading %s of %s to info_store" %(info_name,info_data)
    fail = pylarc.info_set(info_type, serial, info_data)
    if fail:
       print "Unable to write data to info store"
    else:
       print "Wrote data to info_store"
    info_name = "COMPUTER"
    info_type = pylarc.get_info_type_from_string_name(info_name)
    info_data = "a226"
    print "Loading %s of %s to info_store" %(info_name,info_data)
    fail = pylarc.info_set(info_type, serial, info_data)
    if fail:
       print "Unable to write data to info store"
    else:
       print "Wrote data to info_store"

    # write out JSON file
    local_name = "../dat/out/infoTest1.json"
    out_name = os.path.join(os.path.dirname(__file__),local_name)
    print "Writing json file with info attached to %s"  %out_name
    pylarc.matrix_write_json_file_matrixID(serial,out_name)

