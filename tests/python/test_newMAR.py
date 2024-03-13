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
from pylarc import *
from ctypes import *
import math


if __name__ == '__main__':

    verbose = 0

    #*##################################################
    #*   For MARmode testing, we will set up some initial regions
    #*   that we grab by entering a few points into the store
    #*   then try some points near by and see if they collapse
    #*   into the regions.
    #*##################################################
    if (cvar.MARmode):
        regionbitparam = 3
        zeroregionbitparam = regionbitparam
    else:
        regionbitparam = 2
        zeroregionbitparam = 2 

    # run with warnings
    # initialize_larc(26, 24, 10, regionbitparam, zeroregionbitparam,0)
    
    # give warnings
    initialize_larc(26, 24, 10, regionbitparam, zeroregionbitparam,1)

    print("During initialization the following scalars are loaded\n")
    print("   0, 1, -1, 2, .5, -.5, -.707=-1/sqrt(2)\n")

    if (cvar.MARmode):
        print("In MARmode, with scalarType =%s" %cvar.scalarTypeStr)
        print("  With regionbitparam = %d and" %regionbitparam)
        print("  MARregions center at multiples of 1/2^(%d)" %regionbitparam)
        print("  and are usually 1/2^(%d) wide\n" %(regionbitparam-1))
    else:
        print("In SPRmode, with scalarType =%s" %cvar.scalarTypeStr)
        print("  With regionbitparam = %d and" %regionbitparam)
        print("  zeroregionbitparam = %d and" %zeroregionbitparam)
        print("  SPRregions center at multiples of 1/2^(%d)" %regionbitparam)
        print("  and are usually 1/2^(%d) wide\n" %regionbitparam)


    
    # Set up some regions by entering some scalars into the store
    print("Entering some initial scalars into the matrix store")
    # 0 is already in the store
    zeroID = get_zero_pID(0, 0)
    print("   0 is in the matrix store, with matID %d" %zeroID)
    # enter 1+i   (or 1)
    if (cvar.scalarTypeStr in ("Complex","MPComplex","MPRatComplex")):
        # at present require C-style complex number in string
        notNearZeroID = get_valID_from_valString("1+I*1")
        p = 1+1j
        print("   %d is in the matrix store, with matID %d" %(p,notNearZeroID))
    else:
        notNearZeroID = get_valID_from_valString("1")
        p = 1    # small = 1/64
        print("   %d is in the matrix store, with matID %d" %(p,notNearZeroID))

    sqrt2ID = get_valID_from_valString(value_to_string(math.sqrt(2),
                                                       cvar.scalarTypeStr))
    print("   sqrt(2) is in the matrix store, with matID %d" %sqrt2ID)
    twoID = get_valID_from_valString("2")
    print("   2 is in the matrix store, iwth matID %d" %twoID)
    print(" ")
        

    # full region width for MAR if grabbed the nbhrs it wanted
    fullW = 1.0/(2**2)   
    # 2 full region width
    full2W = 2.0/(2**2)
    # half region width for MAR if it only got the original regionTile
    halfW = 1.0/(2**3)
    # quarter region width
    quarW = 1.0/(2**4)
    print(fullW,full2W,halfW,quarW)
    print(" ")

    # negative full region width
    mfullW = -1.0/(2**2)
    # 2 full region width
    mfull2W = -2.0/(2**2)
    # half region width
    mhalfW = -1.0/(2**3)
    # quarter region width
    mquarW = -1.0/(2**4)
    print(mfullW,mfull2W,mhalfW,mquarW)
    print(" ")

    
    # full region width
    print("OUT: A pt a full regionWidth(1/4) from 0 should be in other region")
    EfID = get_valID_from_valString(value_to_string(fullW, cvar.scalarTypeStr))
    print("   Moving 1/4 east from zero we have matID %d" %EfID)
    print("   and has scalar value in the store of %s"
          %get_readableString_scalar_from_pID_and_coords(EfID,0,0))

    WfID = get_valID_from_valString(value_to_string((mfullW), cvar.scalarTypeStr))
    print("   Moving 1/4 west from zero we have matID %d" %WfID)
    print(" ")
    print("   and has scalar value in the store of %s"
          %get_readableString_scalar_from_pID_and_coords(WfID,0,0))        

    # two full region width
    print("NEXTOUT: Pts 2 full regionWidth(1/4) from 0")
    E2fID = get_valID_from_valString(value_to_string(full2W,
                                                    cvar.scalarTypeStr))
    print("   Moving 1/2 east from zero we have matID %d" %E2fID)
    print("   and has scalar value in the store of %s"
          %get_readableString_scalar_from_pID_and_coords(E2fID,0,0))

    W2fID = get_valID_from_valString(value_to_string((mfull2W), cvar.scalarTypeStr))
    print("   Moving 1/2 west from zero we have matID %d" %W2fID)
    print("   and has scalar value in the store of %s"
          %get_readableString_scalar_from_pID_and_coords(W2fID,0,0))
    print(" ")
        

    
    # quarter region width
    print("IN: A pt a quarter of a region from 0 should be in the region")
    EqID = get_valID_from_valString(value_to_string(quarW, cvar.scalarTypeStr))
    print("   Moving 1/16 east from zero we have matID %d" %EqID)
    print("   and has scalar value in the store of %s"
          %get_readableString_scalar_from_pID_and_coords(EqID,0,0))

    WqID = get_valID_from_valString(value_to_string((mquarW), cvar.scalarTypeStr))
    print("   Moving 1/16 west from zero we have matID %d" %WqID)
    print("   and has scalar value in the store of %s"
           %get_readableString_scalar_from_pID_and_coords(WqID,0,0))

    print(" ")
        

    # half region width
    print("BORDER: A pt a half a region from 0 should be on the border")
    EhID = get_valID_from_valString(value_to_string(halfW, cvar.scalarTypeStr))
    print("   Moving 1/8 east from zero we have matID %d" %EhID)

    WhID = get_valID_from_valString(value_to_string((mhalfW), cvar.scalarTypeStr))
    print("   Moving 1/8 west from zero we have matID %d" %WhID)
    print(" ")



          

    
    print("Stepping 1/32 to the east from zero ")
    i=0
    delta = 0
    while (i<=20):
        delta = delta + 1.0/32
        eID = get_valID_from_valString(value_to_string(delta,
                                                       cvar.scalarTypeStr))
        print("%d \t%s \t%g"
              %(eID,get_readableString_scalar_from_pID_and_coords(eID,0,0),delta))
        i=i+1
    print(" ")

    print("Stepping -1/32 to the east from zero ")
    i=0
    delta = 0
    while (i<=64):
        delta = delta - 1.0/32
        eID = get_valID_from_valString(value_to_string(delta,
                                                       cvar.scalarTypeStr))
        print("%d \t%s \t%g"
              %(eID,get_readableString_scalar_from_pID_and_coords(eID,0,0),delta))
        i=i+1
    print(" ")
        


        
sys.exit(0)

                                   
    # # m = p + 1/16
    # m = p + 1.0/(2**4)
    
    # mID = get_valID_from_valString(value_to_string(m, self.scalarType))
    
    # # Try putting a point right on the inclusive border
    # # m = p - 1/32
    # m = p - 1.0/(2**5)
    # if (cvar.MARmode):
    #     m = p - 1.0/(2**7)   # from a nbhr region ...
        
    # mID = get_valID_from_valString(value_to_string(m, self.scalarType))
    # print("The scalar values stored in notNearZero and m are:")
    # print_naive(notNearZeroID)
    # print_naive(mID)
    
    
    
    # self.assertEqual(matrix_mult(sqrt2ID, sqrt2ID), twoID)
    







    
    # #*##############################
    # # inverse permutation matrices #
    # #*##############################
    # print("\nPI_0 matrix is:")
    # PI_0 = pylarc.create_invShufMat(0)
    # pylarc.print_naive(PI_0)

    # print("\nPI_1 matrix is:")
    # PI_1 = pylarc.create_invShufMat(1)
    # pylarc.print_naive(PI_1)

    # print("\nPI_2 matrix is:")
    # PI_2 = pylarc.create_invShufMat(2)
    # pylarc.print_naive(PI_2)

    # print("\nPI_3 matrix is:")
    # PI_3 = pylarc.create_invShufMat(3)
    # pylarc.print_naive(PI_3)

    # # print("\nPI_4 matrix is:")
    # # PI_4 = pylarc.create_invShufMat(4)
    # # pylarc.print_naive(PI_4)


    # #*#######################
    # # print roots of unity  #
    # #*#######################
    # print("\n")
    # pylarc.print_pow2_roots_unity(1)
    # pylarc.print_pow2_roots_unity(2)
    # pylarc.print_pow2_roots_unity(3)


    # #*#############################
    # # create D matrices in python #
    # #*#############################
    # print("\nD_0 matrix is:")
    # D_0 = pylarc.create_FFT_DMat(0)
    # pylarc.print_naive(D_0)

    # print("\nD_1 matrix is:")
    # D_1 = pylarc.create_FFT_DMat(1)
    # pylarc.print_naive(D_1)

    # print("\nD_2 matrix is:")
    # D_2 = pylarc.create_FFT_DMat(2)
    # pylarc.print_naive(D_2)

    # print("\nD_3 matrix is:")
    # D_3 = pylarc.create_FFT_DMat(3)
    # pylarc.print_naive(D_3)


    # #*#############################
    # # create C matrices in python #
    # #*#############################
    # print("\nC_1 matrix is:")
    # C_1 = pylarc.create_FFT_CMat(1)
    # pylarc.print_naive(C_1)

    # print("\nC_2 matrix is:")
    # C_2 = pylarc.create_FFT_CMat(2)
    # pylarc.print_naive(C_2)

    # print("\nC_3 matrix is:")
    # C_3 = pylarc.create_FFT_CMat(3)
    # pylarc.print_naive(C_3)


    # #*###############################
    # # create FFT matrices in python #
    # #*###############################
    # print("\nF_1 matrix is:")
    # F_1 = pylarc.create_FFTMat(1)
    # pylarc.print_naive(F_1)

    # print("\nF_2 matrix is:")
    # F_2 = pylarc.create_FFTMat(2)
    # pylarc.print_naive(F_2)

    # print("\nF_3 matrix is:")
    # F_3 = pylarc.create_FFTMat(3)
    # pylarc.print_naive(F_3)


    # #*###############################
    # # create FFT matrices in python #
    # #*###############################
    # print("\nCreate a vector\n")
    # A_arr = list(map(str,[1,0,0,0,0,1,0,0]))
    # rowLevel = 3
    # colLevel = 0
    # dimWhole = 1 << colLevel
    # A_mID = pylarc.row_major_list_to_store(A_arr,rowLevel,colLevel,dimWhole)
    # pylarc.print_naive(A_mID)

    # print("Now multiply FFT matrix by the vector to get the result\n")
    # B_mID = pylarc.matrix_mult(F_3,A_mID)
    # pylarc.print_naive(B_mID)




    # #*###############################
    # # precision test                #
    # #*###############################

    # k = 3
    # print("\nCalculating the principal 2^%d-th root of unity:" %k)
    # #root_matID = pylarc.principal_pow2_root_unity_pID(k)
    # root_matID = pylarc.k_th_power_of_n_th_root_of_unity_pID(1, 1<<k, 0)

    # print("\nPrinting 2^%d-th root of unity:" %k)
    # pylarc.print_naive(root_matID)

    # exp2Root_mID = pylarc.matrix_mult(root_matID,root_matID)

    # print("\nPrinting 2^%d-th root of unity squared:" %k)
    # pylarc.print_naive(exp2Root_mID)

    # print("\nIf the value just printed is not =I, then the precision test FAILED")
    # print("The locality-sensitive hashing is currently using the values:")
    # print("     zeroregionbitparam = %d" %zeroregionbitparam)
    # print("     regionbitparam = %d" %regionbitparam)
    # if ((zeroregionbitparam!=50) or (regionbitparam!=50)):
    #     print("You might want to try 50 and 50 for 2^3 = 8th roots of unity")
