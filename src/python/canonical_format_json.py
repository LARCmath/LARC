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
sys.path.append(os.path.join(os.path.dirname(__file__),"../"))
import pylarc
import argparse
from timeit import default_timer as timer



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="This routine reads a LARCMatrix json file into the matrix store and then \n outputs a new json file containing the same matrix.  If the input json was in the legacy format \n (where scalar values were not strings) then you must set the flag --legacy to correctly \n read the file and convert to the new format.\n The output json comes from a nearly empty matrix store, hence will \n have matrixIDs that are usually smaller than those in the input file.")
    parser.add_argument("matrix_in", help="the path to the input json file containing the LARC matrix")
    parser.add_argument("matrix_out", help="the path to place the \"canonical\" output json file for the same LARC matrix")
    parser.add_argument("--legacy", action="store_true", help="set if reading in json files with legacy format.")
    parser.add_argument("-v", "--verbose", action="store_true", help="set to verbose printing.")

    larc_starter = pylarc.LarcSettings()
    larc_starter.addArgsToParser(parser)

    args = parser.parse_args()
    larc_starter.parseArgs(args)

    larc_starter.getLevelFromFiles(args.matrix_in)
    larc_starter.initLarc()
#    if args.verbose:
#        larc_starter.getLevelFromFiles(args.matrix_in)
#        larc_starter.initLarc()
#    else:
#        with pylarc.stdout_redirected():
#            larc_starter.getLevelFromFiles(args.matrix_in)
#            larc_starter.initLarc()

    larc_starter.check_scalarType_matches_files(args.matrix_in)

    # read matrix into larc from file
    if args.legacy:
        read_routine = pylarc.read_larcMatrix_file_legacy_return_matID
    else:
        read_routine = pylarc.read_larcMatrix_file_return_matID
    if args.verbose:
        print("Reading in {}matrix...".format(args.legacy*"legacy "))
        matID = read_routine(args.matrix_in)
    else:
#        with pylarc.stdout_redirected():
        matID = read_routine(args.matrix_in)

    if args.verbose:
        print("Writing matrix...")
        pylarc.write_larcMatrix_file_by_matID(matID, args.matrix_out)
    else:
#        with pylarc.stdout_redirected():
        pylarc.write_larcMatrix_file_by_matID(matID, args.matrix_out)


    sys.exit(0)

