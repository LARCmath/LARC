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
## FFT implements a sparse block recursive FFT (see Van Loan)
   #   F_k = C_k * (I_2 @ F_(k-1) ) * PI_k
   # where 
   #   PI_k is the 2^k by 2^k inverse shuffle matrix 
   #       (its transverse is its inverse and would shuffle a column vector)
   #       PI_2 = (1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1) is CNOT
   #       PI_3 = (10000000; 00100000; 00001000; 00000010; 
   #              01000000; 00010000; 00000100; 00000001) 
   #   C_k is the matrix constructed by having 
   #       UL = I_(k-1), UR = D_(k-1); 
   #       LL = I_(k-1), LR = -D_(k-1) 
   #   D_k is the diagonal matrix which had 1, w, w^2, ....w^(2^k - 1) 
   #       on the diagonal, where w is the root of unity cooresponding 
   #       to the next largest size fourier transform, that is 
   #       the 2^(k+1)th root unity.   e.g. D_1 = (1, 0: 0 , i);
   # We end up getting a formula for the FFT in terms of a product of 
   # block matrices:

   #   F_k =    (product_for j= 0   to k-1)  I_j @ C_(k-j)  
   #                                       * I_k @ F_0
   #          * (product_for j= k-1 to 0  )  I_j @ PI_(k-j) 

   # Since F_0 = I_0 and PI_1 = I_1  two of these terms are identity 
   # matrices (also interesting to note:  C_1 = H), so we have:

   #   F_k =    (product_for j= 0   to k-1)  I_j @ C_(k-j)  
   #          * (product_for j= k-2 to 0  )  I_j @ PI_(k-j) 

   # where @ is the tensor product, and * is matrix multiply.

##########################################################################

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

        ## THESE VALUES WORK
        # trunc_to_zero_bits = 54
        # rnd_sig_bits = 1000

        ## THESE VALUES WORK
        # trunc_to_zero_bits = 50
        # rnd_sig_bits = 50

        ## THESE VALUES FAIL
        trunc_to_zero_bits = 60
        rnd_sig_bits = 60

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
    print "\nPI_0 matrix is:"
    PI_0 = pylarc.create_perm_inv_matrixID(0)
    pylarc.print_matrix_naive_by_matrixID(PI_0)

    print "\nPI_1 matrix is:"
    PI_1 = pylarc.create_perm_inv_matrixID(1)
    pylarc.print_matrix_naive_by_matrixID(PI_1)

    print "\nPI_2 matrix is:"
    PI_2 = pylarc.create_perm_inv_matrixID(2)
    pylarc.print_matrix_naive_by_matrixID(PI_2)

    print "\nPI_3 matrix is:"
    PI_3 = pylarc.create_perm_inv_matrixID(3)
    pylarc.print_matrix_naive_by_matrixID(PI_3)

    # print "\nPI_4 matrix is:"
    # PI_4 = pylarc.create_perm_inv_matrixID(4)
    # pylarc.print_matrix_naive_by_matrixID(PI_4)


    #########################
    # print roots of unity  #
    #########################
    print "\n"
    pylarc.print_pow2_roots_unity(1)
    pylarc.print_pow2_roots_unity(2)
    pylarc.print_pow2_roots_unity(3)


    ###############################
    # create D matrices in python #
    ###############################
    print "\nD_0 matrix is:"
    D_0 = pylarc.create_fft_D_matrixID(0)
    pylarc.print_matrix_naive_by_matrixID(D_0)

    print "\nD_1 matrix is:"
    D_1 = pylarc.create_fft_D_matrixID(1)
    pylarc.print_matrix_naive_by_matrixID(D_1)

    print "\nD_2 matrix is:"
    D_2 = pylarc.create_fft_D_matrixID(2)
    pylarc.print_matrix_naive_by_matrixID(D_2)

    print "\nD_3 matrix is:"
    D_3 = pylarc.create_fft_D_matrixID(3)
    pylarc.print_matrix_naive_by_matrixID(D_3)


    ###############################
    # create C matrices in python #
    ###############################
    print "\nC_1 matrix is:"
    C_1 = pylarc.create_fft_C_matrixID(1)
    pylarc.print_matrix_naive_by_matrixID(C_1)

    print "\nC_2 matrix is:"
    C_2 = pylarc.create_fft_C_matrixID(2)
    pylarc.print_matrix_naive_by_matrixID(C_2)

    print "\nC_3 matrix is:"
    C_3 = pylarc.create_fft_C_matrixID(3)
    pylarc.print_matrix_naive_by_matrixID(C_3)


    #################################
    # create FFT matrices in python #
    #################################
    print "\nF_1 matrix is:"
    F_1 = pylarc.create_fft_matrix_matrixID(1)
    pylarc.print_matrix_naive_by_matrixID(F_1)

    print "\nF_2 matrix is:"
    F_2 = pylarc.create_fft_matrix_matrixID(2)
    pylarc.print_matrix_naive_by_matrixID(F_2)

    print "\nF_3 matrix is:"
    F_3 = pylarc.create_fft_matrix_matrixID(3)
    pylarc.print_matrix_naive_by_matrixID(F_3)


    #################################
    # create FFT matrices in python #
    #################################
    print("\nCreate a vector\n")
    A_arr = pylarc.buildArray([1,0,0,0,0,1,0,0])
    rowLevel = 3
    colLevel = 0
    dimWhole = 1 << colLevel
    A_mID = pylarc.row_major_list_to_store_matrixID(A_arr,rowLevel,colLevel,dimWhole)
    pylarc.print_matrix_naive_by_matrixID(A_mID)

    print("Now multiply FFT matrix by the vector to get the result\n")
    B_mID = pylarc.matrix_mult_matrixID(F_3,A_mID)
    pylarc.print_matrix_naive_by_matrixID(B_mID)




    #################################
    # precision test                #
    #################################

    k = 3
    print "\nCalculating the principal 2^%d-th root of unity:" %k
    root_matID = pylarc.principal_pow2_root_unity(k)

    print "\nPrinting 2^%d-th root of unity:" %k
    pylarc.print_matrix_naive_by_matrixID(root_matID)

    exp2Root_mID = pylarc.matrix_mult_matrixID(root_matID,root_matID)

    print "\nPrinting 2^%d-th root of unity squared:" %k
    pylarc.print_matrix_naive_by_matrixID(exp2Root_mID)

    print "\nIf the value just printed is not i, then the precision test FAILED"
    print "because the approximation method is currently using the values:"
    print "     trunc_to_zero_bits = %d" %trunc_to_zero_bits
    print "     rnd_sig_bits = %d" %rnd_sig_bits
