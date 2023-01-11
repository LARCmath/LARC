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
import sys 
sys.path.append(os.path.join(os.path.dirname(__file__),"../../src"))
import pylarc
import numpy as np
import random
from ctypes import *

if __name__ == '__main__':
	# This version references matrices by matrixID instead of pointers
	
    print("This code tests some basic matrix building and reading routines\n")

    # initialize larc
    mat_store_exp = 26
    op_store_exp = 24
    max_level = 10
    # regionbitparam = -1   # default value
    # zeroregionbitparam = -1  # default value
    regionbitparam = 3   # default value
    zeroregionbitparam = 3  # default value
    verbose = 1
    pylarc.create_report_thread(1800)
    pylarc.initialize_larc(mat_store_exp,op_store_exp,max_level,regionbitparam,zeroregionbitparam,verbose)

    # In the Makefile you can compile with different scalarType values
    scalarTypeStr = pylarc.cvar.scalarTypeStr

    # Calculate number of matrices created, then print part of matrix store
    num_matrices_made = pylarc.num_matrices_created()
    print("\n%d matrices have been created" %num_matrices_made)
    end = num_matrices_made - 1

    # build array in C from Python list of scalars
    print("Using row_major_list_to_store on data entered from python\n")

    # parameters for entering the python array into the store
    level = 2 
    dim_whole = 2**level
    randSize = 3*dim_whole-2
    print("randSize = ", randSize)

    # Create list of randSize random numbers
    if scalarTypeStr in ("Real", "MPRational", "MPReal"):
        randVals = [ random.random() for i in range(randSize)]
    elif scalarTypeStr in ("Complex", "MPRatComplex", "MPComplex"):
        randVals = [ np.complex(random.random(),random.random()) for i in range(randSize)]
    elif scalarTypeStr in ("Integer", "MPInteger"):
        randVals = [ random.randrange(0,10001) for i in range(randSize)]
    else:
        raise Exception('Do not know how to build matrix for type %s.'%scalarTypeStr)

    # create symmetric or Hermitian matrix from the
    # dim_whole*(dim_whole+1)/2 random numbers
    a = []
    b = randVals
    print(b)

    c = []
    c.extend(list(b[0:2]))
    c.extend([0]*(dim_whole-2))
    a.append(c)
    counter = 2
    for i in range(dim_whole-2):
        c = [0]*i
        c.extend(list(b[counter:counter+3]))
        c.extend([0]*(dim_whole-i-3))
        a.append(c)
        counter = counter + 3
    c =[0]*(dim_whole-2)
    c.extend(list(b[counter:counter+2]))
    a.append(c)
    
    #print("\na: ", a)
    #sys.exit(0)

    amat = np.matrix(a)
    serial = pylarc.add_numpy_matrix_to_matrix_store(amat)

    filename_json = "../dat/out/tridiag.lev%d.%s.json" %(level,scalarTypeStr)
    # pylarc.print_naive(serial)
    # print("\n")
    pylarc.fprint_larcMatrixFile(serial,os.path.join(os.path.dirname(__file__),filename_json))