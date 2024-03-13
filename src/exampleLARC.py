#        exampleLARC.py
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
from ctypes import *

# initialize larc and make sure it's set up for integers
pylarc.initialize_larc(26, 24, 15, -1, -1, 1)
#if pylarc.cvar.scalarTypeDef != 'z':
#    print("ERROR: remake with type MPINTEGER")
#    sys.exit()

# manipulate some numbers in python
typeChar = pylarc.cvar.scalarTypeDef
if typeChar == 'i':
        x = 2**5
        more = [x, 2**10, 2**15, 2**20]
elif typeChar == 'z':
        x = 2**25
        more = [x, 2**100, 2**50, 2**75]
elif typeChar == 'r':
        x = 3.77*(2**5)
        more = [x, x*2**10, x*2**15, x*2**20]
elif typeChar == 'c': # reformat to C style strings
        # note that python represents complex #s as a+b*1j
        # and we've set up our C code to expect them as a+I*b
        # also note that we can have no spaces between characters!
        x = "3.77+I*2.34"
        more = [x, "27.1+I*0", "0+I*27.1", "3.77+I*-2.34"]
elif typeChar == 'q': # need a better way to do this
        x = "32/3" 
        more = [x, "57/5", "22/7", "16/33"]
else: # in case we define another type
        x = 2**25
        more = [x, 2**100, 2**50, 2**75]

myFirstMatID = pylarc.row_major_list_to_store(list(map(str, more)), 1, 1, 2)
print("my first matrix is")
pylarc.print_naive(myFirstMatID)

print("\nentry squared and then multiplied by 10 is")
squareID = pylarc.matrix_entrySquared(myFirstMatID, "10")
pylarc.print_naive(squareID)

myfile = "../tests/dat/out/temp_"+typeChar+"_matrix.json"
print(myfile)

pylarc.fprint_larcMatrixFile(squareID, myfile)

print("\n entry squared done a different way (no factor of 10):")
sqID2 = pylarc.square_matrix_elements(myFirstMatID,pylarc.FUNC_A)
pylarc.print_naive(sqID2)
print("")


# lets look at the result of building an identity and zero matrix 
# from scratch using the row major reader and print
vals = [0 for i in range(64)]
vals_string = pylarc.map_to_str(vals, "Integer")
zeroID_3_3 = pylarc.row_major_list_to_store(vals_string,3,3,8)
print("The 8 by 8 Zero matrix is")
pylarc.print_naive(zeroID_3_3)

vals = [1*(0==(i%9)) for i in range(64)]
vals_string = pylarc.map_to_str(vals, "Integer")
identityID_3 = pylarc.row_major_list_to_store(vals_string,3,3,8)
print("The 8 by 8 identity matrix is")
pylarc.print_naive(identityID_3)

# Usually one would grab the matrixID's for the zero and
# identity matrices from the preloaded matrix store, using
# the internal faster C code as follows.
normalway_zeroID_3_3 = pylarc.get_zero_pID(3,3)
normalway_identityID_3 = pylarc.get_identity_pID(3)

# We verify these two methods produced the same matrices
if (normalway_zeroID_3_3 == zeroID_3_3):
        print("The zero matrix built using a row major read method and the")
        print("zero matrix retrieved using standard C routine are the same.")
else:
        print("Something is horribly wrong, zero matrix built with row major read")
        print("is not equal to zero matrix retrieved by standard C routine.")
        

# We verify these two methods produced the same matrices
if (normalway_identityID_3 == identityID_3):
        print("The identity matrix built using a row major read method and the")
        print("identity matrix retrieved using standard C routine are the same.")
else:
        print("Something is horribly wrong, identity matrix built with row major read")
        print("is not equal to identity matrix retrieved by standard C routine.")

pylarc.clean_matrix_storage()

v2ID = pylarc.read_larcMatrixFile(myfile)
pylarc.print_naive(v2ID)







