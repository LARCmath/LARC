#!/usr/bin/env python3

# NOTE: to run on command line
#   python3 -m unittest -v test_unittest_scalars
# To run individual test classes from the module, do (for example) 
#   python3 -m unittest -v test_unittest_scalars.TestLocalityRegions

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
from pylarc import *
from ctypes import *
import math
import unittest


class TestLocalityRegions(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pass

    def setUp(self):
        # Define string for use in formating filenames
        # (cvar is in pylarc)
        self.scalarType = cvar.scalarTypeStr

    def test_twobytwo_doubling(self):
        level = 1
        dim_whole = 2**level

        # with stdout_redirected():
        #    initialize_larc(26,24,10,-1,-1,1)
        initialize_larc(26,24,10,-1,-1,0)   # run quiet

        if self.scalarType in ("Integer", "MPInteger"):
            arr_a = [4, 0, 0, 3]
            arr_b = [8, 0, 0, 6]
        else:
            arr_a = [.75, 0, 0, .125]
            arr_b = [1.5, 0, 0, .25]
        a_ID = row_major_list_to_store(map_to_str(arr_a, self.scalarType), level, level, dim_whole)
        b_ID = row_major_list_to_store(map_to_str(arr_b, self.scalarType), level, level, dim_whole)
        self.assertEqual(matrix_add(a_ID, a_ID), b_ID)

    def test_twobytwo_subtraction(self):
        level = 1
        dim_whole = 2**level

        # with stdout_redirected():
        #    initialize_larc(26,24,10,-1,-1,1)
        initialize_larc(26,24,10,-1,-1,0)    # run quiet

        if self.scalarType in ("Integer", "MPInteger"):
            arr_a = [4, 0, 0, 3]
            arr_b = [8, 0, 0, 6]
        else:
            arr_a = [.75, 0, 0, .125]
            arr_b = [1.5, 0, 0, .25]
        a_ID = row_major_list_to_store(map_to_str(arr_a, self.scalarType), level, level, dim_whole)
        b_ID = row_major_list_to_store(map_to_str(arr_b, self.scalarType), level, level, dim_whole)
        self.assertEqual(matrix_diff(b_ID, a_ID), a_ID)

    def test_twobytwo_negation(self):
        level = 1
        dim_whole = 2**level

        # with stdout_redirected():
        #     initialize_larc(26,24,10,-1,-1,1)
        initialize_larc(26,24,10,-1,-1,0)    # run quiet

        if self.scalarType in ("Integer", "MPInteger"):
            arr_a = [4, 0, 0, 3]
            arr_b = [8, 0, 0, 6]
        else:
            arr_a = [.75, 0, 0, .125]
            arr_b = [1.5, 0, 0, .25]
        a_ID = row_major_list_to_store(map_to_str(arr_a, self.scalarType), level, level, dim_whole)
        zero_ID = get_zero_pID(level, level)
        self.assertEqual(matrix_diff(a_ID, a_ID), zero_ID)

        
    # #################### #
    #    SCALAR TESTS START HERE       #
    # #################### #

    def test_subtract_self(self):
        # with stdout_redirected():
        #    initialize_larc(26,24,10,-1,-1,1)
        initialize_larc(26,24,10,-1,-1,0)   # run quiet

        if self.scalarType in ("Integer", "MPInteger"):
            mID = get_valID_from_valString("40")
        else:
            mID = get_valID_from_valString("0.4")
        scalarM1 = get_valID_from_valString("-1")
        neg_mID = matrix_mult(scalarM1, mID)
        zero = get_zero_pID(0, 0)
        self.assertEqual(matrix_add(neg_mID, mID), zero)

    @unittest.skipIf(cvar.scalarTypeStr in ("Boolean", "Integer", "MPInteger"), "ScalarType can not be Integer or Boolean")
    def test_zeroregionbitparam_within_region(self):
        """
        We choose our parameters to make the zero region have size (-1/32,1/32)
        for SPR mode and [-1/32,1/32) in MAR mode (for complex numbers, this
        size applies to both real and imaginary directions).
        """
        if (cvar.MARmode):
            # In MARmode, zero is at the left (and bottom) boundaries of its
            # region, and will claim the neighbor regions on those sides.
            regionbitparam = 5
            zeroregionbitparam = 5 
        else:
            # In SPR mode, zero is at the center of its region.
            regionbitparam = 4
            zeroregionbitparam = 4   # was 2

        # run quiet
        initialize_larc(26, 24, 10, regionbitparam, zeroregionbitparam,0)
        # give warnings
        # initialize_larc(26, 24, 10, regionbitparam, zeroregionbitparam,1)

        # small = 1/64
        small = 1.0/(2**6) 
        smallID = get_valID_from_valString(value_to_string(small,
                  self.scalarType))
        zeroID = get_zero_pID(0, 0)
        # print("The small value is", small)
        # print("The threshold is", get_zerorealthresh())
        self.assertEqual(smallID, zeroID)

        negsmall = -small
        negsmallID = get_valID_from_valString(value_to_string(negsmall,
                  self.scalarType))
        self.assertEqual(negsmallID, zeroID)

        if (cvar.scalarTypeStr in ("Complex","MPComplex","MPRatComplex")):
            # small = 1/64
            small = 1j/(2**6) 
            smallID = get_valID_from_valString(value_to_string(small,
                      self.scalarType))
            zeroID = get_zero_pID(0, 0)
            # print("The small value is", small)
            # print("The threshold is", get_zerorealthresh())
            self.assertEqual(smallID, zeroID)

            negsmall = -small
            negsmallID = get_valID_from_valString(value_to_string(negsmall,
                      self.scalarType))
            self.assertEqual(negsmallID, zeroID)

    @unittest.skipIf(cvar.scalarTypeStr in ("Boolean", "Integer", "MPInteger"), "ScalarType can not be Integer or Boolean")
    def test_zeroregionbitparam_region_boundaries(self):
        """
        We choose our parameters to make the zero region have size (-1/32,1/32)
        for SPR mode and [-1/32,1/32) in MAR mode (for complex numbers, this
        size and these boundaries apply to both real and imaginary directions
        for any x+iy with either x or y equal to zero).
        """
        if (cvar.MARmode):
            # In MARmode, zero is at the left (and bottom) boundaries of its
            # region, and will claim the neighbor regions on those sides.
            regionbitparam = 5
            zeroregionbitparam = 5 
        else:
            # In SPR mode, zero is at the center of its region.
            regionbitparam = 4
            zeroregionbitparam = 4   # was 2

        # run quiet
        initialize_larc(26, 24, 10, regionbitparam, zeroregionbitparam,0)
        # give warnings
        # initialize_larc(26, 24, 10, regionbitparam, zeroregionbitparam,1)

	# notSmall = 1/32 should be on the region boundary
        notSmall = 1.0/(2**5)
        notSmallID = get_valID_from_valString(value_to_string(notSmall, self.scalarType))
        zero = get_zero_pID(0, 0)
        self.assertNotEqual(notSmallID, zero)

        negnotSmall = -notSmall
        negnotSmallID = get_valID_from_valString(value_to_string(negnotSmall,
             self.scalarType))

        if (cvar.MARmode):
            # negative boundary is inclusive
            self.assertEqual(negnotSmallID, zero)
        else:
            # negative boundary is exclusive
            self.assertNotEqual(negnotSmallID, zero)

        if (cvar.scalarTypeStr in ("Complex","MPComplex","MPRatComplex")):
            notSmall = 1j/(2**5)
            notSmallID = get_valID_from_valString(
                value_to_string(notSmall, self.scalarType))
            self.assertNotEqual(notSmallID, zero)

            negnotSmall = -notSmall
            negnotSmallID = get_valID_from_valString(
                value_to_string(negnotSmall, self.scalarType))

            if (cvar.MARmode):
                # negative boundary is inclusive
                self.assertEqual(negnotSmallID, zero)
            else:
                # negative boundary is exclusive
                self.assertNotEqual(negnotSmallID, zero)

    @unittest.skipIf(cvar.scalarTypeStr in ("Boolean", "Integer", "MPInteger"), "ScalarType can not be Integer or Boolean")
    def test_zeroregionbitparam_outside_boundaries(self):
        """
        when zeroregionbitparam = t anything y with |y|>2^t will be nonzero.
        """
        if (cvar.MARmode):
            # In MARmode, zero is at the left (and bottom) boundaries of its
            # region, and will claim the neighbor regions on those sides.
            regionbitparam = 5
            zeroregionbitparam = 5 
        else:
            # In SPR mode, zero is at the center of its region.
            regionbitparam = 4
            zeroregionbitparam = 4   # was 2

        # run quiet
        initialize_larc(26, 24, 10, regionbitparam, zeroregionbitparam,0)
        # give warnings
        # initialize_larc(26, 24, 10, regionbitparam, zeroregionbitparam,1)

	# notSmall = 1/16 should be outside of boundary
        notSmall = 1.0/(2**4)
        notSmallID = get_valID_from_valString(value_to_string(notSmall, self.scalarType))
        zero = get_zero_pID(0, 0)
        self.assertNotEqual(notSmallID, zero)

        negnotSmall = -notSmall
        negnotSmallID = get_valID_from_valString(value_to_string(negnotSmall, self.scalarType))
        self.assertNotEqual(negnotSmallID, zero)

        if (cvar.scalarTypeStr in ("Complex","MPComplex","MPRatComplex")):
            notSmall = 1j/(2**4)
            notSmallID = get_valID_from_valString(
                value_to_string(notSmall, self.scalarType))
            self.assertNotEqual(notSmallID, zero)

            negnotSmall = -notSmall
            negnotSmallID = get_valID_from_valString(
                value_to_string(negnotSmall, self.scalarType))
            self.assertNotEqual(negnotSmallID, zero)

    @unittest.skipIf(cvar.scalarTypeStr in ("Boolean", "Integer", "MPInteger"), "ScalarType can not be Integer or Boolean")
    def test_regionbitparam_within_region(self):
        """
        Set t = regionbitparam. For non-complex types, the integer scalar k>0
        has a region with boundaries [k-1/2^(t+1),k+1/2^(t+1)) when LARC is
        operating in SPR mode, and a super-region with boundaries
        [k-1/2^t,k+1/2^t) when operating in MAR mode.
        In contrast, the region for the scalar -k has boundaries
        (-k-1/2^(t+1),-k+1/2^(t+1)] in SPR mode and a super-region with
        boundaries [-k-1/2^t,-k+1/2^t) when in MAR mode (note that the
        inclusive boundary shifts for SPR to be near the axis, but is always
        to the more negative side for MAR).
        For complex types, we consider scalars k + I*\ell. Since the real and
        imaginary parts are handled separately, the MAR mode boundaries for
        any integers k,\ell are [k-1/2^t,k+1/2^t) in the real direction and 
        [\ell-1/2^t,\ell+1/2^t) in the imaginary direction. The regions for the
        SPR mode are similar to those described above, noting that the signs of
        k and \ell determine which boundaries are inclusive.
	The largest power of two that we can add to k without being on or
        beyond a region boundary is 1/2^{t+2} in SPR mode and 1/2^{t+1} in MAR
        mode.
	The reason for using 1+I with complex types is to ensure that we are
	away from the near-zero regions.
        """
        # Set region parameters so that the region around integer k>0 is
        # [k-1/32,k+1/32) and is 1/16 wide (in both real and imag directions)
        if (cvar.MARmode):
            regionbitparam = 5
        else:
            regionbitparam = 4

        # run quiet
        initialize_larc(26, 24, 10, regionbitparam,-1,0)  

        # create a scalar value that is not near zero for either the
        # real or imag direction. Since we are inside the boundary, the
        # test will look the same for either MAR or SPR.
        if (cvar.scalarTypeStr in ("Complex","MPComplex","MPRatComplex")):
            # at present require C-style complex number in string
            notNearZeroID = get_valID_from_valString("1+I*1")
            p = 1+1j
        else:
            notNearZeroID = get_valID_from_valString("1")
            p = 1

        # m = 1 + I \pm 1/64   for complex (no I for real)           
        m = p  + 1.0/(2**6) 
        mID = get_valID_from_valString(value_to_string(m, self.scalarType))
        self.assertEqual(mID, notNearZeroID)
        m = p  - 1.0/(2**6) 
        mID = get_valID_from_valString(value_to_string(m, self.scalarType))
        self.assertEqual(mID, notNearZeroID)
        if (cvar.scalarTypeStr in ("Complex","MPComplex","MPRatComplex")):
           m = p  + 1j/(2**6)
           mID = get_valID_from_valString(value_to_string(m, self.scalarType))
           self.assertEqual(mID, notNearZeroID)
           m = p  - 1j/(2**6) 
           mID = get_valID_from_valString(value_to_string(m, self.scalarType))
           self.assertEqual(mID, notNearZeroID)


    @unittest.skipIf(cvar.scalarTypeStr in ("Boolean", "Integer", "MPInteger"), "ScalarType can not be Integer or Boolean")
    def test_regionbitparam_outside_region(self):
        """
        See comments for test_regionbitparam_inside_region. In this case,
        we avoid the complications with the boundaries by adding or subtracting
        a number big enough to guarantee we are outside the region.
        """
        # Set region parameters so that the region around zero is
        # [-1/32,1/32) and is 1/16 wide (in both real and imag directions)
        if (cvar.MARmode):
            regionbitparam = 5
        else:
            regionbitparam = 4

        # with stdout_redirected():
        #    initialize_larc(26, 24, 10, regionbitparam,-1,1)
        initialize_larc(26, 24, 10, regionbitparam,-1,0)  # run quiet

        if (cvar.scalarTypeStr in ("Complex","MPComplex","MPRatComplex")):
            # at present require C-style complex number in string
            notNearZeroID = get_valID_from_valString("1+I*1")
            p = 1+1j
        else:
            notNearZeroID = get_valID_from_valString("1")
            p = 1

        # m = p + 1/16
        m = p + 1.0/(2**4)
        mID = get_valID_from_valString(value_to_string(m, self.scalarType))
        self.assertNotEqual(mID, notNearZeroID)
        m = p - 1.0/(2**4)
        mID = get_valID_from_valString(value_to_string(m, self.scalarType))
        self.assertNotEqual(mID, notNearZeroID)

        if (cvar.scalarTypeStr in ("Complex","MPComplex","MPRatComplex")):
            m = p + 1j/(2**4)
            mID = get_valID_from_valString(value_to_string(m, self.scalarType))
            self.assertNotEqual(mID, notNearZeroID)
            m = p - 1j/(2**4)
            mID = get_valID_from_valString(value_to_string(m, self.scalarType))
            self.assertNotEqual(mID, notNearZeroID)

    @unittest.skipIf(cvar.scalarTypeStr in ("Boolean", "Integer", "MPInteger"), "ScalarType can not be Integer or Boolean")
    def test_regionbitparam_region_boundaries(self):
        """
        See comments for test_regionbitparam_inside_region. In the real number
        case, SPR regions and MAR superregions of the same size containing the
        integer k>0 have the same exclusive and inclusive boundaries; however,
        when k<0 the inclusiveness of the boundaries is switched for SPR. For
        complex numbers, we must deal with four possibilities. In MAR, each of 
        \pm|k| + I*\pm|\ell| has inclusive boundaries on the side closer to
        negative infinity, but for SPR the inclusiveness of each boundary is
        determined by the sign of k or \ell. (k==0 and \ell==0 are tested
        elsewhere; for these cases, SPR has no inclusive boundary perpendicular
        to the axis contained within the SPR region.)
        """
        # Set region parameters so that the region around zero is
        # [-1/32,1/32) and is 1/16 wide (in both real and imag directions)
        if (cvar.MARmode):
            regionbitparam = 5
        else:
            regionbitparam = 4

        # with stdout_redirected():
        #    initialize_larc(26, 24, 10, regionbitparam,-1,1)
        initialize_larc(26, 24, 10, regionbitparam,-1,0)  # run quiet

        if (cvar.scalarTypeStr in ("Complex","MPComplex","MPRatComplex")):
            # test for ++ quadrant
            # at present require C-style complex number in string
            notNearZeroID = get_valID_from_valString("1+I*1")
            p = 1+1j
            m = p - 1/32 # on inclusive border
            mID = get_valID_from_valString(value_to_string(m, self.scalarType))
            self.assertEqual(mID, notNearZeroID)
            m = p + 1/32 # on exclusive border
            mID = get_valID_from_valString(value_to_string(m, self.scalarType))
            self.assertNotEqual(mID, notNearZeroID)
            m = p - 1j/32 # on inclusive border
            mID = get_valID_from_valString(value_to_string(m, self.scalarType))
            self.assertEqual(mID, notNearZeroID)
            m = p + 1j/32 # on exclusive border
            mID = get_valID_from_valString(value_to_string(m, self.scalarType))
            self.assertNotEqual(mID, notNearZeroID)
            # test for +- quadrant
            # at present require C-style complex number in string
            notNearZeroID = get_valID_from_valString("1-I*1")
            p = 1-1j
            m = p - 1/32 # on inclusive border
            mID = get_valID_from_valString(value_to_string(m, self.scalarType))
            self.assertEqual(mID, notNearZeroID)
            m = p + 1/32 # on exclusive border
            mID = get_valID_from_valString(value_to_string(m, self.scalarType))
            self.assertNotEqual(mID, notNearZeroID)
            m = p - 1j/32 # on inclusive border for MAR, exclusive for SPR
            mID = get_valID_from_valString(value_to_string(m, self.scalarType))
            if (cvar.MARmode):
                self.assertEqual(mID, notNearZeroID)
            else:
                self.assertNotEqual(mID, notNearZeroID)
            m = p + 1j/32 # on exclusive border for MAR, inclusive for SPR
            mID = get_valID_from_valString(value_to_string(m, self.scalarType))
            if (cvar.MARmode):
                self.assertNotEqual(mID, notNearZeroID)
            else:
                self.assertEqual(mID, notNearZeroID)
            # test for -+ quadrant
            # at present require C-style complex number in string
            notNearZeroID = get_valID_from_valString("-1+I*1")
            p = -1+1j
            m = p - 1/32 # on inclusive border for MAR, exclusive for SPR
            mID = get_valID_from_valString(value_to_string(m, self.scalarType))
            if (cvar.MARmode):
                self.assertEqual(mID, notNearZeroID)
            else:
                self.assertNotEqual(mID, notNearZeroID)
            m = p + 1/32 # on exclusive border for MAR, inclusive for SPR
            mID = get_valID_from_valString(value_to_string(m, self.scalarType))
            if (cvar.MARmode):
                self.assertNotEqual(mID, notNearZeroID)
            else:
                self.assertEqual(mID, notNearZeroID)
            m = p - 1j/32 # on inclusive border
            mID = get_valID_from_valString(value_to_string(m, self.scalarType))
            self.assertEqual(mID, notNearZeroID)
            m = p + 1j/32 # on exclusive border
            mID = get_valID_from_valString(value_to_string(m, self.scalarType))
            self.assertNotEqual(mID, notNearZeroID)

            # test for -+ quadrant
            # at present require C-style complex number in string
            notNearZeroID = get_valID_from_valString("-1-I*1")
            p = -1-1j
            m = p - 1/32 # on inclusive border for MAR, exclusive for SPR
            mID = get_valID_from_valString(value_to_string(m, self.scalarType))
            if (cvar.MARmode):
                self.assertEqual(mID, notNearZeroID)
            else:
                self.assertNotEqual(mID, notNearZeroID)
            m = p + 1/32 # on exclusive border for MAR, inclusive for SPR
            mID = get_valID_from_valString(value_to_string(m, self.scalarType))
            if (cvar.MARmode):
                self.assertNotEqual(mID, notNearZeroID)
            else:
                self.assertEqual(mID, notNearZeroID)
            m = p - 1j/32 # on inclusive border for MAR, exclusive for SPR
            mID = get_valID_from_valString(value_to_string(m, self.scalarType))
            if (cvar.MARmode):
                self.assertEqual(mID, notNearZeroID)
            else:
                self.assertNotEqual(mID, notNearZeroID)
            m = p + 1j/32 # on exclusive border for MAR, inclusive for SPR
            mID = get_valID_from_valString(value_to_string(m, self.scalarType))
            if (cvar.MARmode):
                self.assertNotEqual(mID, notNearZeroID)
            else:
                self.assertEqual(mID, notNearZeroID)
        else: # not complex number
            # test for positive numbers
            notNearZeroID = get_valID_from_valString("1")
            p = 1
            m = p - 1/32 # on inclusive border
            mID = get_valID_from_valString(value_to_string(m, self.scalarType))
            self.assertEqual(mID, notNearZeroID)
            m = p + 1/32 # on exclusive border
            mID = get_valID_from_valString(value_to_string(m, self.scalarType))
            self.assertNotEqual(mID, notNearZeroID)
            # test for negative numbers
            notNearZeroID = get_valID_from_valString("-1")
            p = -1
            m = p - 1/32 # on inclusive border for MAR, exclusive for SPR
            mID = get_valID_from_valString(value_to_string(m, self.scalarType))
            if (cvar.MARmode):
                self.assertEqual(mID, notNearZeroID)
            else:
                self.assertNotEqual(mID, notNearZeroID)
            m = p + 1/32 # on exclusive border for MAR, inclusive for SPR
            mID = get_valID_from_valString(value_to_string(m, self.scalarType))
            if (cvar.MARmode):
                self.assertNotEqual(mID, notNearZeroID)
            else:
                self.assertEqual(mID, notNearZeroID)


    @unittest.skipIf(cvar.scalarTypeStr in ("Boolean", "Integer", "MPInteger"), "ScalarType can not be Integer or Boolean")  
    # @unittest.skipIf(True, "Not sure why sqrt2*sqrt2 = 2.0 in floating point")
    def test_sqrt2_squared_is_2(self):
        """
        sqrt(2) squared is the same as 2.
        """

        # Set region parameters so that the region around zero is
        # (-1/64,1/64) and is 1/32 wide (in both real and imag directions)
        if (cvar.MARmode):
            regionbitparam = 6
        else:
            regionbitparam = 5

        # with stdout_redirected():
        #    initialize_larc(26, 24, 10, regionbitparam,zeroregionbitparam,1)
        initialize_larc(26, 24, 10, regionbitparam,-1,0)  # run quiet

        sqrt2ID = get_valID_from_valString(value_to_string(math.sqrt(2), self.scalarType))
        twoID = get_valID_from_valString("2")

        self.assertEqual(matrix_mult(sqrt2ID, sqrt2ID), twoID)
        

    @unittest.skipIf(cvar.scalarTypeDef in ("i", "z"), "ScalarType must not be Integer or MPInteger")
    def test_sqrt2_formula_works(self):
        """
        (sqrt2+sqrt2)/2 is the same as sqrt2.
        """

        # with stdout_redirected():

        #    initialize_larc(26, 24, 10, -1,-1,1)
        initialize_larc(26, 24, 10, -1,-1,0)  # run quiet

        halfID = get_valID_from_valString(value_to_string(0.5, self.scalarType))
        sqrt2ID = get_valID_from_valString(value_to_string(math.sqrt(2), self.scalarType))
        doubleID = matrix_add(sqrt2ID, sqrt2ID)
        resultID = matrix_mult(halfID, doubleID)
        failure_msg = "matrices displayed above"
        if resultID != sqrt2ID:
            print("Failure in test_sqrt2_formula_works - printing relevant matrices.")
            print("Half matrix:")
            print_naive(halfID)
            print("Sqrt2 matrix:")
            print_naive(sqrt2ID)
            print("Double Sqrt2 matrix:")
            print_naive(doubleID)
            print("Result matrix:")
            print_naive(resultID)
        self.assertEqual(resultID, sqrt2ID, failure_msg)

#    @unittest.skipIf(cvar.scalarTypeStr in ("Boolean", "Integer", "MPInteger", "Clifford"), "ScalarType can not be Integer or Boolean or Clifford(for now)")
    @unittest.skipIf(cvar.scalarTypeStr in ("Boolean", "Integer", "MPInteger"), "ScalarType can not be Integer or Boolean")
    def test_inv_sqrt2_formula_works(self):
        """
        ((1/sqrt2)^2) is the same as 1/2.
        """

        # with stdout_redirected():
        #    initialize_larc(26, 24, 10, -1,-1,1)
        initialize_larc(26, 24, 10, -1, -1, 0)  # run quiet

        halfID = get_valID_from_valString(value_to_string(0.5, self.scalarType))
        rootID = get_pID_for_enum_const(SCALAR_ENUM_INV_SQRT2)
        resultID = matrix_mult(rootID, rootID);
        self.assertEqual(resultID, halfID)


