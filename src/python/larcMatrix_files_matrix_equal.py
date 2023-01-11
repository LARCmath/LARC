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
import glob
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../"))
import pylarc
import argparse
from timeit import default_timer as timer

## \file larcMatrix_files_matrix_equal.py
#  \brief Determines if two LARCMatries stored in json format are equal.
#
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Determines if two LARCmatrices stored as json files are equal.")
    parser.add_argument("path1", help="the path to the first json file to compare")
    parser.add_argument("path2", help="the path to the second json file to compare")
    parser.add_argument("-v", "--verbose", action="store_true", help="set to verbose printing")

    larc_starter = pylarc.LarcSettings()
    larc_starter.addArgsToParser(parser)

    args = parser.parse_args()
    larc_starter.parseArgs(args)

    larc_starter.getLevelFromFiles(args.path1, args.path2)
    larc_starter.initLarc()
#    if args.verbose:
#        larc_starter.getLevelFromFiles(args.path1, args.path2)
#        larc_starter.initLarc()
#    else:
#        with pylarc.stdout_redirected():
#            larc_starter.getLevelFromFiles(args.path1, args.path2)
#            larc_starter.initLarc()

    larc_starter.check_scalarType_matches_files(args.path1, args.path2)

    # read matrix into larc from file
    if args.verbose:
        print("Reading in matrices...")
        mat1ID = pylarc.read_larcMatrixFile(args.path1)
        mat2ID = pylarc.read_larcMatrixFile(args.path2)
        print()
    else:
#        with pylarc.stdout_redirected():
        mat1ID = pylarc.read_larcMatrixFile(args.path1)
        mat2ID = pylarc.read_larcMatrixFile(args.path2)

    MATRIX_ID_INVALID = -1
    if mat1ID == MATRIX_ID_INVALID:
        print("Matrix at {} is invalid.".format(args.path1))
    if mat2ID == MATRIX_ID_INVALID:
        print("Matrix at {} is invalid.".format(args.path2))
    # compare matrices
    print("The two matrices are {}equal.".format((mat1ID != mat2ID)*"not "))
    if args.verbose:
        print("mat1ID = {}\nmat2ID = {}".format(mat1ID, mat2ID))
        print()

    sys.exit(0)

