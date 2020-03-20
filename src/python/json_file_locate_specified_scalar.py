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
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../"))
import pylarc
import argparse
from timeit import default_timer as timer


if __name__ == '__main__':

    # Returns the number and locations of specific scalar, given a json matrix file
    parser = argparse.ArgumentParser(description="Count the number of entries in a LARCMatrix json file\n with a given scalar value and reports their locations.")
    parser.add_argument("matrix_path", help="the path to the json file containing the LARC matrix")
    parser.add_argument("-x", "--scalar", default="0", help="the scalar value to count (default 0)")
    parser.add_argument("-v", "--verbose", action="store_true", help="set to verbose printing")
    parser.add_argument("-n", "--max_to_report", default="0", help="the maximum number of locations to report (default is 100) ")

    larc_starter = pylarc.LarcSettings()
    larc_starter.addArgsToParser(parser)

    args = parser.parse_args()
    larc_starter.parseArgs(args)

    larc_starter.getLevelFromFiles(args.matrix_path)
    larc_starter.initLarc()
#    if args.verbose:
#        larc_starter.getLevelFromFiles(args.matrix_path)
#        larc_starter.initLarc()
#    else:
#        with pylarc.stdout_redirected():
#            larc_starter.getLevelFromFiles(args.matrix_path)
#            larc_starter.initLarc()

    if (args.max_to_report=="0"):
        if args.verbose:
            print("WARNING: user did not set the maximum number of locations to report: defaulting to 100", file=sys.stderr)
        args.max_to_report = "100"
    if args.verbose:
        print("max_to_report is %s" %args.max_to_report)

    # count number of zeros in matrix
    if args.verbose:
        print("\nLocating {} in matrix...\n".format(args.scalar))
    start = timer()
    locations = pylarc.locate_entries_larcMatrixFile(args.matrix_path, args.scalar, args.max_to_report)
    end = timer()
    num = len(locations)
    [mat_row_level, mat_col_level] = pylarc.get_levels_from_larcMatrix_file(args.matrix_path)
    total_entries = 2**(mat_row_level + mat_col_level)
    if args.verbose:
        print("Computation time: {}\n".format(end - start))
    print("Matrix has \n  value {0} entries:     {1}\n  value not {0} entries: {2}\n  total entries:       {3}\n  value {0} sparsity:    {4}\n".format(
              args.scalar, num, total_entries - num, total_entries, 1.0*num/total_entries) )
    if (num):
       print("Locations:")
       rw = max([len(str(l[0])) for l in locations])
       cw = max([len(str(l[1])) for l in locations])
       for i in range(len(locations)):
           l = locations[i]
           print("{:{}}: ({:{}}, {:{}})".format(i, len(str(len(locations))), l[0], rw, l[1], cw))
    sys.exit(0)

