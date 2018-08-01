#!/usr/bin/env python

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
# POC: Jennifer Zito <jszito@super.org>                         #
# Please contact the POC before disseminating this code.        #
#                                                               #
#################################################################

from larc_py import *

# standard types using %array_class macros
def buildIntegerArray(arr):
    newarr = intArray(len(arr))

    for i,val in enumerate(arr):
        newarr[i] = val

    return newarr

def buildInt64Array(arr):
    newarr = int64Array(len(arr))

    for i,val in enumerate(arr):
        newarr[i] = val

    return newarr

def buildRealArray(arr):
    newarr = realArray(len(arr))

    for i,val in enumerate(arr):
        newarr[i] = val

    return newarr

def buildComplexArray(arr):
    newarr = complexArray(len(arr))

    for i,val in enumerate(arr):
        newarr[i] = val

    return newarr

# non-standard type using %array_functions macros
def buildArray(arr):
    newarr = new_scalarTypeArray(len(arr))

    for i,val in enumerate(arr):
        scalarTypeArray_setitem(newarr, i, val)

    return newarr
    
def str_scalarTypeArray(arr, length):
    return '['+', '.join([str(scalarTypeArray_getitem(arr, i)) for i in range(length)])+']'
