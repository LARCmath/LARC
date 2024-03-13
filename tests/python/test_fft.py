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
import glob
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../../src"))
import pylarc
from ctypes import *


if __name__ == '__main__':

    verbose=1

    #*##################################################
    #*   Print description sparse block algorithm     #*
    #*   from Cooley-Tukey Radix-2 Factorization      #*
    #*   see Van Loan, Computational Frameworks for   #*
    #*   the Fast Fourier Transform, p.21             #*
    #*##################################################
    print("\nWe will create the matrices used in the Cooley-Tukey")
    print("radix-2 sparse block recursive FFT that is from Van Loan,")
    print("Computational Frameworks for the Fast Fourier Transform. p.21")
    if verbose:
        print("The level k, 2**k by 2**k Fourier matrix F_k can be")
        print("generated recursively by the equation")
        print("(subscripts represent levels, not size=2^level):")
        print("   F_k = C_k * (I_1 @ F_(k-1) ) * PI_k")
        print("where ")
        print("   PI_k is the 2^k by 2^k inverse shuffle matrix ")
        print("     note: PI_k's transverse is its inverse and ")
        print("           would shuffle a column vector")
        print("     examples:")
        print("       PI_2 = (1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1)")  # swap 
        print("       PI_3 = (10000000; 00100000; 00001000; 00000010; ")
        print("               01000000; 00010000; 00000100; 00000001)")
        print("   C_k is the matrix constructed by having quadrant ")
        print("       submatrices [UL, UR; LL, LR] = ")
        print("       [I_(k-1), D_(k-1); I_(k-1), -D_(k-1)] ")   
        print("   D_k is the diagonal matrix which had ")
        print("       1, w, w^2, ....w^(2^k - 1) on the diagonal,") 
        print("       where w is the root of unity cooresponding")
        print("       to the next largest size fourier transform, ")
        print("       that is the 2^(k+1)th root unity.")
        print("       e.g. D_1 = (1, 0: 0 , i);")
        print("We end up getting a formula for the FFT in terms of")
        print("a product of block matrices:")
        print("F_k = (product_for j= 0 to k-1) I_j @ C_(k-j) * I_k @ F_0")
        print("       * (product_for j= k-1 to 0) I_j @ PI_(k-j)")
        print("Since F_0 = I_0 and PI_1 = I_1, two of these terms are identity")
        print("matrices (also interesting to note:  C_1 = H), so we have:")
        print("F_k = (product_for j= 0 to k-1) I_j @ C_(k-j) * ")
        print("      (product_for j= k-2 to 0) I_j @ PI_(k-j) ")
        print("where @ is the tensor product, and * is matrix multiply.")


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
        #*  Sample failure values for LARC LSH parameters #*
        #*  are to set         regionbitparam = 60            #*
        #*  and                zeroregionbitparam = 60      #*
        #*####################################################
        #*  Default values for LARC LSH parameters are       #*
        #*  both equal to LDBL_MANT_DIG -2 when scalarType is Real #*
        #*  regionbitparam = -1    # 	                    #*
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


    #*##############################
    # inverse permutation matrices #
    #*##############################
    print("\nPI_0 matrix is:")
    PI_0 = pylarc.create_invShufMat(0)
    pylarc.print_naive(PI_0)

    print("\nPI_1 matrix is:")
    PI_1 = pylarc.create_invShufMat(1)
    pylarc.print_naive(PI_1)

    print("\nPI_2 matrix is:")
    PI_2 = pylarc.create_invShufMat(2)
    pylarc.print_naive(PI_2)

    print("\nPI_3 matrix is:")
    PI_3 = pylarc.create_invShufMat(3)
    pylarc.print_naive(PI_3)

    # print("\nPI_4 matrix is:")
    # PI_4 = pylarc.create_invShufMat(4)
    # pylarc.print_naive(PI_4)


    #*#######################
    # print roots of unity  #
    #*#######################
    print("\n")
    pylarc.print_pow2_roots_unity(1)
    pylarc.print_pow2_roots_unity(2)
    pylarc.print_pow2_roots_unity(3)


    #*#############################
    # create D matrices in python #
    #*#############################
    print("\nD_0 matrix is:")
    D_0 = pylarc.create_FFT_DMat(0)
    pylarc.print_naive(D_0)

    print("\nD_1 matrix is:")
    D_1 = pylarc.create_FFT_DMat(1)
    pylarc.print_naive(D_1)

    print("\nD_2 matrix is:")
    D_2 = pylarc.create_FFT_DMat(2)
    pylarc.print_naive(D_2)

    print("\nD_3 matrix is:")
    D_3 = pylarc.create_FFT_DMat(3)
    pylarc.print_naive(D_3)


    #*#############################
    # create C matrices in python #
    #*#############################
    print("\nC_1 matrix is:")
    C_1 = pylarc.create_FFT_CMat(1)
    pylarc.print_naive(C_1)

    print("\nC_2 matrix is:")
    C_2 = pylarc.create_FFT_CMat(2)
    pylarc.print_naive(C_2)

    print("\nC_3 matrix is:")
    C_3 = pylarc.create_FFT_CMat(3)
    pylarc.print_naive(C_3)


    #*###############################
    # create FFT matrices in python #
    #*###############################
    print("\nF_1 matrix is:")
    F_1 = pylarc.create_FFTMat(1)
    pylarc.print_naive(F_1)

    print("\nF_2 matrix is:")
    F_2 = pylarc.create_FFTMat(2)
    pylarc.print_naive(F_2)

    print("\nF_3 matrix is:")
    F_3 = pylarc.create_FFTMat(3)
    pylarc.print_naive(F_3)


    #*###############################
    # create FFT matrices in python #
    #*###############################
    print("\nCreate a vector\n")
    A_arr = list(map(str,[1,0,0,0,0,1,0,0]))
    rowLevel = 3
    colLevel = 0
    dimWhole = 1 << colLevel
    A_mID = pylarc.row_major_list_to_store(A_arr,rowLevel,colLevel,dimWhole)
    pylarc.print_naive(A_mID)

    print("Now multiply FFT matrix by the vector to get the result\n")
    B_mID = pylarc.matrix_mult(F_3,A_mID)
    pylarc.print_naive(B_mID)


