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
sys.path.append(os.path.join(os.path.dirname(__file__),"../"))
import pylarc
import argparse
from timeit import default_timer as timer


## \file larcMatrix_file_count_specified_scalar.py
#  \brief Count the number of entries in a LARCMatrix json file with a given
#  scalar value
#
if __name__ == '__main__':

    # Returns the count of a specified scalar in a json matrix file 
    parser = argparse.ArgumentParser(description="Count the number of entries in a LARCmatrix json file with a given scalar value.")
    parser.add_argument("matrix_path", help="the path to the json file containing the LARC matrix")
    parser.add_argument("-x", "--scalar", default="0", help="the scalar value to count (default 0)")
    parser.add_argument("-v", "--verbose", action="store_true", help="set to verbose printing")

    larc_starter = pylarc.LarcSettings()
    larc_starter.addArgsToParser(parser)

    args = parser.parse_args()
    larc_starter.parseArgs(args)

    # Exit gracefully if matrix has more than 2^32 entries
    [mat_row_level, mat_col_level] = pylarc.get_levels_from_larcMatrix_file(args.matrix_path)
    if (mat_row_level + mat_col_level > 32):
        print("Sorry, matrix is too big to process.  Levels = [{0}, {1}]".format(mat_row_level, mat_col_level))
        sys.exit(1)

    larc_starter.getLevelFromFiles(args.matrix_path)
    larc_starter.initLarc()
#    if args.verbose:
#        larc_starter.getLevelFromFiles(args.matrix_path)
#        larc_starter.initLarc()
#    else:
#        with pylarc.stdout_redirected():
#            larc_starter.getLevelFromFiles(args.matrix_path)
#            larc_starter.initLarc()

    # Count the number of entries in a LARC json matrix file with a given scalar value.
    if args.verbose:
        print("\nCounting {} in matrix...\n".format(args.scalar))

    start = timer()
    pID = pylarc.read_larcMatrixFile(args.matrix_path)
    num = int(pylarc.matrix_count_entries(pID, args.scalar))
    end = timer()

    total_entries = 2**(mat_row_level + mat_col_level)
    if args.verbose:
        print("Computation time: {}\n".format(end - start))
    l = len(args.scalar)
    print("Matrix has \n  value {0} entries:     {1}\n  value not {0} entries: {2}\n  total entries:{5}      {3}\n  value {0} sparsity:    {4}\n".format(
              args.scalar, num, total_entries - num, total_entries, 1.0*num/total_entries, " "*l) )

    sys.exit(0)

