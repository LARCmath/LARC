#################################################################
#                                                               #
# Copyright 2014, Institute for Defense Analyses                #
# 4850 Mark Center Drive, Alexandria, VA; 703-845-2500          #
# This material may be reproduced by or for the US Government   #
# pursuant to the copyright license under the clauses at DFARS  #
# 252.227-7013 and 252.227-7014.                                #
#                                                               #
# LARC : Linear Algebra via Recursive Compression               #
# Authors:                                                      #
#   - Steve Cuccaro (IDA-CCS)                                   #
#   - John Daly (LPS)                                           #
#   - John Gilbert (UCSB, IDA adjunct)                          #
#   - Jenny Zito (IDA-CCS)                                      #
#                                                               #
# Additional contributors are listed in "LARCcontributors".     #
#                                                               #
# POC: Jennifer Zito jszito@super.org or jennyzito@gmail.com    #
# Please contact the POC before disseminating this code.        #
#                                                               #
#################################################################

This is the LARC library implementing
LARC: Linear Algebra via Recursive Compression
which is written primarily in C, and 
has a Swig-generated Python interface.  

To get started, compile a fresh copy of the code in the 
top-level project directory.
  make clean
  make

To see examples of how the code works in C and to test
functionality after modifications you can type
  bin/exampleLARC

To see examples of how the code works in Python and to test
whether code is functioning well you can type (from the
top-level project directory):
  ./tests/python/test_matrixBuild.py
  ./tests/python/test_math.py
  ./tests/python/test_localHash.py
  etc.

It is currently the state that LARC has the ScalarType 
determined by the following Makefile commands:  
  make TYPE=COMPLEX
  make TYPE=REAL
  make TYPE=INTEGER
If no command is given the default is REAL (which is double).



