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
import glob
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../../src"))
import pylarc
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
        print("This machine is a CPU-cycle server")
    else:
        print("This machine is a desktop work station")


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

        trunc_to_zero_bits = 50
    ##  trunc_to_zero_bits = 52
    ##  trunc_to_zero_bits = 54
    ##  rnd_sig_bits = 1000
        rnd_sig_bits = 30

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
        verbose = 1
        pylarc.initialize_larc(matrix_exponent,op_exponent,max_level,rnd_sig_bits,trunc_to_zero_bits,verbose)
        pylarc.create_report_thread(180)
        print_naive = 1
        print_nonzeros = 1
        print("Problem size is small enough to run on desktop")
        if print_naive:
           print("  will print files of naive matrices")
        else: 
           print("  not printing files of naive matrices")
        if print_nonzeros:
           print("  will print files of nonzero matrices\n")
        else: 
           print("  not printing files of nonzero matrices\n")
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
        verbose = 0
        pylarc.initialize_larc(matrix_exponent,op_exponent,max_level,rnd_sig_bits,trunc_to_zero_bits,verbose)
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
    print("Seppuku check to see if program is to large to occur once every 10 minutes.\n")


    #########################
    # print roots of unity  #
    #########################
    # verbose = 0
    # n = 4
    # print("\nRunning code to produce the %d-th roots of unity." %n)
    # pylarc.print_n_th_roots_of_unity(4,verbose)
    # # pylarc.print_n_th_roots_of_unity(24,verbose)

    ####################################################
    # load principal root of unity and return matrixID #
    ####################################################
    # verbose = 0
    # print("\nRunning code to produce the matrixID for the n-th principal root of unity.")
    # p1_mID = pylarc.principal_n_th_root_of_unity_matID(1, verbose)
    # print("principal nth root of unity for n=1 has matrixID %d" %p1_mID)
    # p2_mID = pylarc.principal_n_th_root_of_unity_matID(2, verbose)
    # print("principal nth root of unity for n=2 has matrixID %d" %p2_mID)
    # p3_mID = pylarc.principal_n_th_root_of_unity_matID(3, verbose)
    # print("principal nth root of unity for n=3 has matrixID %d" %p3_mID)
    # p4_mID = pylarc.principal_n_th_root_of_unity_matID(4, verbose)
    # print("principal nth root of unity for n=4 has matrixID %d" %p4_mID)
    # # p24_mID = pylarc.principal_n_th_root_of_unity_matID(24, verbose)
    # # print("principal nth root of unity for n=24 has matrixID %d" %p24_mID)


    #################################################################
    # fill an array with the matrixIDs of the nth roots of unity   ##
    # verbose = 0 is run quiet, verbose = 1 subroutines chatty, and  ##
    # verbose = 2 both local and subroutine calls chatty as well.   ##
    #################################################################
    ## verbose = 2
    verbose = 1
    if (verbose>1):
        call_verbose = 1
    else:    
        call_verbose = 0
    n = 24
    print("\nRunning code to produce a list of matrixIDs for all the %d-th roots of unity." %n)
    n_array = [0]*(n)
    if (verbose>1):
        print("The length of the array is %d" %len(n_array))
    for k in range(n):
        n_array[k]  = pylarc.k_th_power_of_n_th_root_of_unity_matID(k,n,call_verbose)
    print("\nHere is the array of the %dth roots of unity:" %n)   
    print(n_array)
    if (verbose>1):
        print("\nThe stored values of these roots are:")
        print("")
        for k in range(n):
           pylarc.print_naive_by_matID(n_array[k])
           print("")
    print("\nNow we can look to see if multiplication is closed, by looking")
    print("at the matrixIDs of products of pairs of these roots.")
    success = 1
    for k in range(n):
        for j in range(k+1):
            my_matrixID = pylarc.matrix_mult_matrixID(n_array[j],n_array[k])
            flag = 0  # initialize as though there was a closure failure
            for i in range(n):
                if (my_matrixID == n_array[i]):
                    flag = 1   # we found the matrixID of product in the preloaded roots
            if (flag == 0):
                success = 0
                print("\nMatrix multiplication of the %d-th roots of unity is not closed" %n)
                print("since the product of the %d-th power and %d-th power" %(j,k,))
                m = (j+k) % n
                print("which should have been the %d-th power which had value" %m )
                pylarc.print_naive_by_matID(n_array[m])

                print("The computed result has matrixID %g instead of %g" %(my_matrixID,
                                                                            n_array[m]))
                print("Which differs from the expected value by")
                # complex_dif = 0.1 + I*0.0
                # print(complex_dif)
    print("\n")
    if success:
        print("All products of pairs of roots of unity produce existing roots.")
    else:
        print("At least one product of roots of unity produced a value")
        print("that was not a preloaded root of unity, so nbhd hash is not optimal.")
        

    # sys.exit(0)

    ## OLDER VERSION OF CODE
    # verbose = 1
    # if (verbose>1):
    #     call_verbose = 1
    # else:    
    #     call_verbose = 0
    # n = 5
    # print("\nRunning code to produce a list of matrixIDs for all the %d-th roots of unity." %n)
    # n5_array = [0]*(n)
    # if (verbose>1):
    #     print("The length of the array is %d" %len(n5_array))
    # for k in range(n):
    #     n5_array[k]  = pylarc.k_th_power_of_n_th_root_of_unity_matID(k,n,call_verbose)
    # print("\nHere is the array of matrixIDs of the %dth roots of unity:" %n)
    # print(n5_array)
    # if (verbose>1):
    #     print("\nThe stored values of these roots are:")
    #     print("")
    #     for k in range(n):
    #        print("The matrixID is %d" %n5_array[k])
    #        print("The %d-th power of the %d-th root of unity is" %(k,n))
    #        pylarc.print_naive_by_matID(n5_array[k])
    #        print("")
    # print("\nNow we can look to see if multiplication is closed, by looking")
    # print("at the matrixIDs of products of pairs of these roots.")
    # for k in range(n):
    #     for j in range(k+1):
    #         my_matrixID = pylarc.matrix_mult_matrixID(n5_array[j],n5_array[k])
    #         print("The (%d,%d)th product has matrix ID" %(j,k))
    #         print(my_matrixID)
    # print("\n")

    # print("\nTODO: Fix these non desirable:")
    # print("\nFor n = 5 we see that the product of the 1st and 4th roots,")
    # my_matrixID = pylarc.matrix_mult_matrixID(n5_array[1],n5_array[4])
    # print("the (%d,%d)th product has matrix ID" %(1,4))
    # print(my_matrixID)
    # print("which has stored value")
    # pylarc.print_naive_by_matID(my_matrixID)
    # print("")
    # print("Whereas if the nbhd hash was perfect we expected to see")
    # print(n5_array[0])
    # print("which has stored value")
    # pylarc.print_naive_by_matID(n5_array[0])


    # mID = pylarc.matrix_diff_matrixID(int64_t A_mID, int64_t B_mID)
    # get numerical value mID
    # if val
    # i = 1
    # while value > 1/(2^i)for i in 
    #     see if the value is < 
    #         my_matrixID = pylarc.matrix_mult_matrixID(n5_array[j],n5_array[k])
    #         print("The (%d,%d)th product has matrix ID" %(j,k))
    #         print(my_matrixID)
    # print("\n")
