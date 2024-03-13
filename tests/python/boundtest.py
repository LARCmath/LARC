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
import sys 
sys.path.append(os.path.join(os.path.dirname(__file__),"../../src"))
import pylarc
import numpy as np
import random
from ctypes import *

if __name__ == '__main__':
	# This version references matrices by matrixID instead of pointers
	
    print("This code tests some basic matrix building and reading routines\n", flush=True)

    # initialize larc
    mat_store_exp = 26
    op_store_exp = 24
    max_level = 10
    regionbitparam = -1   # default value
    zeroregionbitparam = -1  # default value
    pylarc.create_report_thread(1800)
    verbose = 1
    pylarc.initialize_larc(mat_store_exp,op_store_exp,max_level,regionbitparam,zeroregionbitparam,verbose)

    pylarc.print_naive(pylarc.row_major_list_to_store(["0.65625"], 0, 0, 1))
    print("\n", flush=True)

    # In the Makefile you can compile with different scalarType values
    # Define string for using in formating filenames
    scalarTypeStr = pylarc.cvar.scalarTypeStr
    print("Scalar Type is: ", scalarTypeStr, flush=True)
    if scalarTypeStr not in ("Upper", "Lower"):
        print("This is only useful for bounding types.", flush=True)
        sys.exit(1)

    # Matrix creation - a
    avals = ["0.75", "0.5", "0.25", "0.5"]
    a = pylarc.row_major_list_to_store(avals, 1, 1, 2)
    print("\nA Matrix:", flush=True)
    pylarc.print_naive(a)
    npa = np.array([[0.75, 0.5], [0.25, 0.5]])
    print(npa, flush=True)

    # Matrix creation - b
    bvals = ["0.6875", "0.625", "0.3125", "0.375"]
    b = pylarc.row_major_list_to_store(bvals, 1, 1, 2)
    print("\nTrue B Matrix:", flush=True)
    pylarc.print_naive(b)
    npb = npa @ npa
    print(npb, flush=True)

    # Matrix creation - c
    cvals = ["0.671875", "0.65625", "0.328125", "0.34375"]
    c = pylarc.row_major_list_to_store(cvals, 1, 1, 2)
    print("\nTrue C Matrix:", flush=True)
    pylarc.print_naive(c)
    npc = npa @ npb
    print(npc, flush=True)

    # Matrix multiply
    bx = pylarc.matrix_mult(a,a)
    print("\nCalc B Matrix:", flush=True)
    pylarc.print_naive(bx)

    # Matrix multiply
    cx1 = pylarc.matrix_mult(bx,a)
    print("\nCalc C1=B@A Matrix:", flush=True)
    pylarc.print_naive(cx1)

    # Matrix multiply
    cx2 = pylarc.matrix_mult(a,bx)
    print("\nCalc C2=A@B Matrix:", flush=True)
    pylarc.print_naive(cx2)

    # New problem
    print("\nNew Problem\n")
    np.set_printoptions(linewidth=200)
    nptm = np.array([[0, 0.25, 0.25,    0, 0.25,    0, 0.25,    0],
                     [0, 0.25, 0.25,    0, 0.25,    0, 0.25,    0],
                     [0,    0, 0.25, 0.25,    0, 0.25,    0, 0.25],
                     [0, 0.25,    0, 0.25, 0.25,    0, 0.25,    0],
                     [0,    0, 0.25,    0, 0.25, 0.25,    0, 0.25],
                     [0, 0.25,    0, 0.25,    0, 0.25, 0.25,    0],
                     [0,    0, 0.25,    0, 0.25,    0, 0.25, 0.25],
                     [0, 0.25,    0, 0.25,    0, 0.25,    0, 0.25] ])
    print(nptm, flush=True)
    print(flush=True)
    npv = np.array([[1, 0, 0, 0, 0, 0, 0, 0]])
    print(npv[0], flush=True)
    for _1 in range(20):
        npv = npv @ nptm
        print(npv[0], flush=True)

    print(flush=True)
    oneseventh = pylarc.row_major_list_to_store(["1/7"], 0, 0, 1)
    pylarc.print_naive(oneseventh)

    print("\n", flush=True)
    vals = [ "0", "0.25", "0.25",    "0", "0.25",    "0", "0.25",    "0",
             "0", "0.25", "0.25",    "0", "0.25",    "0", "0.25",    "0",
             "0",    "0", "0.25", "0.25",    "0", "0.25",    "0", "0.25",
             "0", "0.25",    "0", "0.25", "0.25",    "0", "0.25",    "0",
             "0",    "0", "0.25",    "0", "0.25", "0.25",    "0", "0.25",
             "0", "0.25",    "0", "0.25",    "0", "0.25", "0.25",    "0",
             "0",    "0", "0.25",    "0", "0.25",    "0", "0.25", "0.25",
             "0", "0.25",    "0", "0.25",    "0", "0.25",    "0", "0.25" ]
    tm = pylarc.row_major_list_to_store(vals, 3, 3, 8)
    pylarc.print_naive(tm)
    print(flush=True)
    v = pylarc.row_major_list_to_store(["1"] + 7*["0"], 0, 3, 8)
    pylarc.print_naive(v)
    for _1 in range(20):
        v = pylarc.matrix_mult(v,tm)
        pylarc.print_naive(v)


