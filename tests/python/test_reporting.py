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
import glob
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../../src"))
import pylarc
import numpy as np
from ctypes import *


if __name__ == '__main__':

    #################################
    ##   SET THESE PARAMETERS      ##
    #################################
    max_level = 8          ##  problem_size is always power of two!


    ####################################################
    ##   Find out if machine is desktop workstation   ##
    ##   or a CPU-cycle servers (cs1-cs6)             ##
    ####################################################
    machine = os.uname()[1]
    cs = 0        # on desktop workstation, with smaller memory
    if (machine.find('cs') >= 0):
        cs = 1    # on CPU-cycle server cs1-cs6, with larger memory
        print "This machine is a CPU-cycle server"
    else:
        print "This machine is a desktop work station"


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

    ####################################################################
    ##    Testing to see what approximation functions.
    ##    The zerobit thresh parameters that work for F_3 are:
    ##      -z 53 and smaller
    ##    The zerobit thresh parameters that fail for F_3 are:
    ##      -z 54 and larger  
    ## 
    ##    The rounding function does not effect whether it
    ##    works in ranges -s 10 to -s 1000
    ## 
    ## 
    ####################################################################


    ## SMALL STORES for working on desktop
    if max_level <= 8:    
        matrix_exponent = 22
        op_exponent = 19   


        trunc_to_zero_bits = 54
        rnd_sig_bits = 1000

        ######################################################
        ##  Sample failure values for LARC approximation    ##
        ##  are to set         rnd_sig_bits = 60            ##
        ##  and                trunc_to_zero_bits = 60      ##
        ######################################################
        ##  Default values for LARC approximation are       ##
        ##  both equal to DBL_MANT_DIG -2                   ##
        ##  rnd_sig_bits = -1    # default is 53 bits       ##
        ##  trunc_to_zero_bits = -1 # OLD default is 1074 bits  ##
        ##  trunc_to_zero_bits = -1 # OLD default is 1074 bits  ##
        ##  NOTE:  testing shows -z 47 will work            ##
        ######################################################
        # trunc_to_zero_bits = 52


        ######################################################

        ##  TODO: find out the space in which this test fails!!!

        ##  DBL_MANT_DIG is the number of digits in FLT_MANT  ##
        ##  why aren't we using DBL_MANT_BITS  ??????? the number of bits
        ##  used in the mantissa


        ######################################################

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



    ################################
    # inverse permutation matrices #
    ################################
    PI_3 = pylarc.create_perm_inv_matrixID(3)
    ###############################
    # create D matrices in python #
    ###############################
    D_3 = pylarc.create_fft_D_matrixID(3)
    ###############################
    # create C matrices in python #
    ###############################
    C_3 = pylarc.create_fft_C_matrixID(3)
    #################################
    # create FFT matrices in python #
    #################################
    print "\nF_3 matrix is:"
    F_3 = pylarc.create_fft_matrix_matrixID(3)
    pylarc.print_matrix_naive_by_matrixID(F_3)
    ##########################
    # take  fft of a vector  #
    ##########################
    A_arr = pylarc.buildArray([1,0,0,0,0,1,0,0])
    rowLevel = 3
    colLevel = 0
    dimWhole = 1 << colLevel
    A_mID = pylarc.row_major_list_to_store_matrixID(A_arr,rowLevel,colLevel,dimWhole)
    B_mID = pylarc.matrix_mult_matrixID(F_3,A_mID)

    #########################
    # print roots of unity  #
    #########################
    print "\n"
    pylarc.print_pow2_roots_unity(3)

    #########################
    # precision testing     #
    #########################
    pylarc.precision_testing()


    # #############################
    # # test reporting
    # #############################
    # pylarc.matrix_store_report("../dat/out/temp.mat_report")
    # pylarc.op_store_report("../dat/out/temp.op_report")
    # pylarc.rusage_report(0,"../dat/out/temp.rusage0_report")
    # pylarc.rusage_report(1,"../dat/out/temp.rusage1_report")

    print "\n"
    pylarc.matrix_store_report("stdout")
    print "\n"
    pylarc.op_store_report("stdout")
    print "\n"
    pylarc.rusage_report(0,"stdout")
    print "\n"
    pylarc.rusage_report(1,"stdout")
