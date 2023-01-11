//                       typeCLIFFORD.h
/******************************************************************
 *                                                                *
 * Copyright (C) 2014, Institute for Defense Analyses             *
 * 4850 Mark Center Drive, Alexandria, VA; 703-845-2500           *
 * This material may be reproduced by or for the US Government    *
 * pursuant to the copyright license under the clauses at DFARS   *
 * 252.227-7013 and 252.227-7014.                                 *
 *                                                                *
 * Authors:                                                       *
 *   - Steve Cuccaro (IDA-CCS)                                    *
 *   - John Daly (LPS)                                            *
 *   - John Gilbert (UCSB, IDA adjunct)                           *
 *   - Mark Pleszkoch (IDA-CCS)                                   *
 *   - Jenny Zito (IDA-CCS)                                       *
 *                                                                *
 * Additional contributors are listed in "LARCcontributors".      *
 *                                                                *
 * Questions: larc@super.org                                      *
 *                                                                *
 * All rights reserved.                                           *
 *                                                                *
 * Redistribution and use in source and binary forms, with or     *
 * without modification, are permitted provided that the          *
 * following conditions are met:                                  *
 *   - Redistribution of source code must retain the above        *
 *     copyright notice, this list of conditions and the          *
 *     following disclaimer.                                      *
 *   - Redistribution in binary form must reproduce the above     *
 *     copyright notice, this list of conditions and the          *
 *     following disclaimer in the documentation and/or other     *
 *     materials provided with the distribution.                  *
 *   - Neither the name of the copyright holder nor the names of  *
 *     its contributors may be used to endorse or promote         *
 *     products derived from this software without specific prior *
 *     written permission.                                        *
 *                                                                *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND         *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,    *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF       *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE       *
 * DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER NOR        *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,   *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT   *
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;   *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)       *
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN      *
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR   *
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, *
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.             *
 *                                                                *
 *****************************************************************/

#ifndef LARC_TYPE_H
#define LARC_TYPE_H


#define USE_CLIFFORD

/*
 * The following Clifford algebras are currently implemented:
 *
 *  - MPR_CLIFFORD_S2: Q(sqrt(2))
 *  - MPC_CLIFFORD_S2: Q(i, sqrt(2))
 *  - MPR_CLIFFORD_S2_S3: Q(sqrt(2), sqrt(3))
 *  - MPC_CLIFFORD_S2_S3: Q(i, sqrt(2), sqrt(3))
 *  - MPR_CLIFFORD_C2: Q(2^(1/3))
 *  - MPC_CLIFFORD_C2: Q(i, 2^(1/3))
 */

// #define MPR_CLIFFORD_S2
#define MPC_CLIFFORD_S2
// #define MPR_CLIFFORD_S2_S3
// #define MPC_CLIFFORD_S2_S3
// #define MPR_CLIFFORD_C2
// #define MPC_CLIFFORD_C2

#define IS_RATIONAL
#define IS_MP

#if defined(MPC_CLIFFORD_S2) || defined(MPC_CLIFFORD_S2_S3) || defined(MPC_CLIFFORD_C2)
#define IS_COMPLEX
#endif


 /*!
  * \brief The CLIFFORD_DIMENSION is the number of new constant values
  * that are adjoined to the rational numbers (including the
  * constant value for 1).  This does not count the imaginary
  * constant values, but just the real ones.
  */
#if defined(MPR_CLIFFORD_S2) || defined(MPC_CLIFFORD_S2)
     #define CLIFFORD_DIMENSION 2

#elif defined(MPR_CLIFFORD_S2_S3) || defined(MPC_CLIFFORD_S2_S3)
     #define CLIFFORD_DIMENSION 4

#elif defined(MPR_CLIFFORD_C2) || defined(MPC_CLIFFORD_C2)
     #define CLIFFORD_DIMENSION 3

#endif  // defined(MPR_CLIFFORD_S2) || defined(MPC_CLIFFORD_S2)

#endif   // #define LARC_TYPE_H
