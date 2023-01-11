#              pylarc.py
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

from __future__ import print_function, division

## \file pylarc.py
#  \brief This module provides access to all the LARC C and Python code for
#  those who run LARC directly (i.e., without a project package).
#
#  The module is set up so that all the LARC functions are in the pylarc.py
#  namespace. It is mostly useful for LARC developers who want to run the
#  unittests before making changes to LARC. Developers who create their own
#  project packages will use an equivalent to pylarc.py in their project
#  directories, which will give them single-namespace access to the C and
#  Python code for their project as well as that for LARC.

import os
import sys
current_directory = os.path.dirname(__file__)
# gv is in the current directory
import gv
# set the PATH_TO_SWIG variable to the current directory
gv.PATH_TO_SWIG = current_directory
# set SWIG_TO_USE to the name of the SWIG-generated python
gv.SWIG_TO_USE = 'larcSWIG'
# larc_utilities is in larc/src/python
sys.path.append(os.path.join(current_directory,'python'))
from larc_utilities import *
