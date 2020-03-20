##                   README.md
###################################################################
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

Welcome to the Version 1.0 LARC library implementing 
   LARC = Linear Algebra via Recursive Compression.
The core code of LARC is written in C and there is a
a SWIG-generated Python3 interface for ease of use.
The python module containing both SWIG-generated python
wrappers for the C code and additional routines that
are written directly in python is called pylarc.

Python testing routines are in the directory tests/python.
There are also commandline python routines which are
not in pylarc and are mostly in the directory src/python.

There is a set of slides in the doc file "about_LARC_2017.pdf"
which describe the ideas behind LARC.
For a summary of changes since the paper see the files:
    doc/LARCimprovements.Updated20171006
    doc/LARCimprovements.Updated20190523 

We have just begun creating Doxygen documentation for LARC in HTML
and LaTeX formats. To create this documentation, say
    make docs
This will create two subdirectories, html and latex. To view the
HTML version, open in a browser window the file "html/index.html".
This command will work in Linux:
    xdg-open html/index.html

To get started, compile a fresh copy of the code in the 
top-level project directory.
  make clean
  make

To see examples of how the code works in C and to test
functionality after modifications:
  cd bin
  ./exampleLARC

To see examples of how the code works in Python and to test
whether code is functioning well you can type (from the
top-level project directory):
  ./tests/python/test_matrixBuild.py
  ./tests/python/test_math.py
  ./tests/python/test_localHash.py
  etc.
To run all the python tests say:
  ./tests/python/script

It is currently the state that LARC has the ScalarType 
determined by the following Makefile commands:  
  make TYPE=COMPLEX
  make TYPE=REAL
  make TYPE=INTEGER
as well as five multiprecision Scalartypes (mpz_t, mpfr_t, mpc_t, mpq_t, larc_mpratcomplex_t): 
  make TYPE=MPINTEGER
  make TYPE=MPREAL
  make TYPE=MPCOMPLEX
  make TYPE=MPRATIONAL
  make TYPE=MPRATCOMPLEX
the last type is constructed from two copies of MPRATIONAL
for the real and imaginary parts of a complex rational number.

If no command is given the default is REAL (which is long double).

Before commiting and pushing a new version of code to the
git repository, all the unit tests should be run by saying:
  make clean
  make unittests

The different ways to initialize LARC are described in the Doxygen
documentations and illustrated various places such as
in tests/python code.

The analog of Hello World for LARC is to ...
