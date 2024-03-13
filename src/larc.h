//                         larc.h
/******************************************************************
 *                                                                *
 * Copyright (C) 2014-2024, Institute for Defense Analyses        *
 * 4850 Mark Center Drive, Alexandria, VA; 703-845-2500           *
 * This material may be reproduced by or for the US Government    *
 * pursuant to the copyright license under the clauses at DFARS   *
 * 252.227-7013 and 252.227-7014.                                 *
 *                                                                *
 * LARC : Linear Algebra via Recursive Compression                *
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

#ifndef LARC_LARC_H
#define LARC_LARC_H

#include <stdio.h> // FILE
#include <stddef.h> // size_t
#include <stdint.h> // int8_t, int32_t, uint32_t, uint64_t, uintptr_t
#include <float.h> // LDBL_MANT_DIG
#include <complex.h> // complex numbers

#include <gmp.h>
#include <mpfr.h>
#include <mpc.h>
#include "type.h"


/*!
 * \include high-level-doc.d
 */

/*!
 * \file larc.h
 * \brief Defines data structures and macros used in LARC
 */

// The Multi-tile Assured Retrieval Hash (MAR) will try to steal
// neighboring regions to also use in its expanded region
// whenever those regions are within epsilon of the scalar
// and are not already occupied.
#define MAR

// To print out alerts on when MAR or SPR snaps occur,
// set the value of ALERT_ON_SNAP to 1.  To disable alerts,
// set the value of ALERT_ON_SNAP to 0.
#define ALERT_ON_SNAP 0

// To save the tile indices for scalars inside each scalar record,
// which will also take more memory, but not require recomputing the
// tile index every time, enable the following definition:
// #define STORE_TILE_INDEX

// The global verbosity parameter is used to set the minimum number of
// messages printed, and will eventually be used everywhere to determine
// what is printed. The default value is BASIC. It may be overridden
// by local variable in individual routines or when LARC is initialized.
// 
// Most of the levels of verbosity are self-explanatory, but the DEBUG level
// deserves further detail. In routines which are actively being debugged,
// the debugging messages should be protected by a statement
//      if (VERBOSE >= DEBUG)
// and once debugging is completed these should be changed to 
//      if (VERBOSE > DEBUG).
// Thus, routines will contain print statements which have been useful for
// debugging in the past, but which are not printed unless VERBOSE == ALL.
typedef enum verbosity {
  SILENT,   // supresses everything but error messages
  BASIC,  // allows warnings
  CHATTY,   // informational messages are printed
  DEBUG,    // situational use (see above)
  ALL       // all print statements will produce output
} verbose_type_t;

extern verbose_type_t VERBOSE;

typedef enum norm {
  L_infty,
  L_1,
  L_2
} norm_type_t;

enum sca_constant_spec {
    SCALAR_ENUM_SQRT2,                 // sqrt(2)
    SCALAR_ENUM_INV_SQRT2,             // 1/sqrt(2)
    SCALAR_ENUM_SQRT3,                 // sqrt(3)
    SCALAR_ENUM_INV_SQRT3,             // 1/sqrt(3)
    SCALAR_ENUM_SQRT6,                 // sqrt(6)
    SCALAR_ENUM_INV_SQRT6,             // 1/sqrt(6)
    SCALAR_ENUM_CUBERT2,               // 2^(1/3)
    SCALAR_ENUM_INV_CUBERT2,           // 2^(-1/3)
    SCALAR_ENUM_CUBERT4,               // 2^(2/3)
    SCALAR_ENUM_INV_CUBERT4,           // 2^(-2/3)
};

// The "ifndef SWIG" keep SWIG from creating access functions for the elements
// of the scalarType structure. SWIG attempts to assign the array pointer which
// represents the multiprecision type to a scalar value, resulting in errors.
// This is not a problem because we should never directly access the scalarType
// from Python; we instead have functions that convert the value to a string
// before it is returned.

#ifndef SWIG
#ifdef USE_INTEGER
  typedef int64_t scalarType;
#endif
#ifdef USE_BOOLEAN
  typedef int64_t scalarType;
#endif
#ifdef USE_COMPLEX
  typedef long double complex scalarType;
#endif
#ifdef USE_REAL
  typedef long double scalarType;
#endif
#ifdef USE_MPINTEGER
  typedef mpz_t scalarType;
#endif
#ifdef USE_MPRATIONAL
  typedef mpq_t scalarType;
#endif

#ifdef USE_MPRATCOMPLEX
  // The mpratcomplex type used to be defined directly as the larc_mpratcomplex_s
  // structure instead of as a one element array of that structure as currently
  // defined below.  The current defintion mirrors the definitions used by the
  // GMP, MPFR, and MPC libraries.  This change, along with a few other changes,
  // was made as part of an attempt to get rid of some segfaults we were having.
  typedef struct larc_mpratcomplex_s
  {
      mpq_t real;
      mpq_t imag;
  } larc_mpratcomplex_t[1];
  typedef larc_mpratcomplex_t scalarType;
#endif

#ifdef USE_MPREAL
  typedef mpfr_t scalarType;
#endif
#ifdef USE_MPCOMPLEX
  typedef mpc_t scalarType;
#endif
#ifdef USE_CLIFFORD
#include "clifford.h"
  typedef clifford_t scalarType;
#endif
#ifdef IS_BOUNDING
#define NUM_EXPONENTS 3
typedef uint32_t bexp_type;
  typedef struct larc_exponent_scalar_s
  {
      bexp_type explist[NUM_EXPONENTS];
      unsigned int is_zero   :1;
      unsigned int is_one    :1;
      unsigned int is_nan    :1;
      unsigned int is_exact  :1;
      unsigned int xflags    :12;
  } larc_exponent_scalar_t[1];
  typedef larc_exponent_scalar_t scalarType;
#endif

// if adding additional types, also see scalarTypeDef in global.h, global.c


#endif // ifndef SWIG

/*  These types are new names for uintptr_t.
    In this way, the type is actually a pointer to the record
    disguised as an integer.  The record_ptr_t is used in 
    hash chain function definitions, but when hash functions
    are called it will be either a matns_ptr_t or mats_ptr_t if we are dealing
    with matrix records, or a op_ptr_t if we are dealing with
    op_store records.   
    Note that structures for these records are local to matrix_store.c
    or op_store.c respectively, so they can not be cast outside
    of these files.
*/
typedef uintptr_t* record_ptr_t; // is used in definition of hash functions
                                // to stand for either matns_ptr_t, mats_ptr_t
                                //  or op_ptr_t

/* We currently use -1 to indicate that nothing is in the matrix store 
   (see variable mat_level_t largest), so
   it can't be changed to a uint without rethinking this. */
typedef int32_t mat_level_t;

// LARC matrices include hash_nodes, so this definition must come first.
typedef struct hash_node
{
  union 
  {
     // both types are 64-bit
     // The matrixID for this packedID could be in either the scalarStore or 
     // matrixStore. As part of a node entry, it shouldn't be -1 (though the
     // associated record pointer could later be nulled)
     int64_t packedID;
     // This is a generic pointer that actually holds either an op_ptr_t or
     // info_ptr_t.
     // (Formerly, it could also have held matns_ptr_t or mats_ptr_t, but
     // nodes in the matrixStore or scalarStore should now only use packedIDs.)
     record_ptr_t record_ptr;
  };
  uint32_t record_hits;       // # of times that this node was the path to a record retrieval
  uint16_t hashfilter;            //  short hash checked as primary test for correct tile
  uint16_t tile_offset_flag  :4;  // 0000 for primary tile, xxyy for real yy and imaginary xx, 01=plus 1, 10=minus 1
  uint16_t hits_maxxed       :1;  // 0 if hits is in range, 1 if maxed out
  uint16_t xflags           :11;  // extra/unused flag positions
  struct hash_node *next;        // pointer to next hash node in chain or to "head"
  struct hash_node *prev;        // pointer to last hash node in chain or to "tail"
} hash_node_t;


#ifndef SWIG
// In LARC MAR space is divided into tiles, if the type is complex
// this tile is indexed by two MPintegers, if the type is real
// then the tile is indexed by one MPinteger.
typedef struct MAR_tile_index
{
#ifdef IS_COMPLEX
  mpz_t real_index;
  mpz_t imag_index;
#else
  mpz_t index;
#endif  
} MAR_tile_index_t;
#endif // SWIG


// The structure of a LARC matrix must be protected from SWIG. The reason is
// that SWIG tries to create access functions for 'scalarType trace_element'
// or when the matrix is 1 by 1 by 'scalarType scalar_value'
// and, when scalarType is multiprecision, attempts to assign the array pointer
// which represents the multiprecision type to a scalar value.
// This is not a problem because we should never directly access the scalarType
// from Python; we instead have functions that convert the value to a string
// before it is returned.
#ifndef SWIG
/*****************************************************************************
 *                                                                           *
 *                   Struct for a compressed matrix                          *
 *                                                                           *
 *  A LARC ((Linear Algebra via Recursive Compression)) matrixVal                  *
 *  either contains four submatrix pointers                                  *
 *  or a single scalar number.                                               *
 *  At level = 0 larc_matrix_t contains a single scalar number and           *
 *  at level > 0 larc_matrix_t contains pointers to the four submatrices.    *
 *  Space in units of int64_t:                                               *
 *         row, level, flags+ref_ct, packedID,                               *
 *         trace or scalar value (2 values if complex, 1 if real)            *
 *         subs (4)                                                          * 
 *                                                                           *
 *****************************************************************************/
// matric record that is for a matrix that is larger than 1 by 1
struct larc_nsmatrix
{
  // UNIQUE IDENTIFICATION since matns_ptr_t and mats_ptr_t can be released
  // and reassigned
  int64_t    packedID; // a unique integer, containing a sequentially-numbered
                       // matrixID as well as certain flag bits

  // LEVELS   matrix is 2^row_level by 2^col_level
  mat_level_t row_level;         
  mat_level_t col_level;

  // FLAGS  - are bits inside a half word
  uint16_t      info;        // 0 if nothing in info store, flip bit for each enum type in info_store used
  unsigned int  hold    :13; // incremented/decremented by user, can remove associated matrix only if 0
  //  unsigned int  MARnhbr :1;  // if set, submatrix[0] is the pID of MARrep
  unsigned int  lock    :1;  // if set, do not EVER remove this matrix from matrix store
  unsigned int  iszero  :1;  // is this a zero matrix  (was unint32_t)
  unsigned int  isid    :1;  // is this an identity matrix (was unint32_t)

  // NUMBER OF APPEARS_AS_SUB MATRICES  - uses other half of the word with flags
  uint32_t    appears_as_sub_count; 

  //   trace_element is the trace for square matrices with level > 0
  scalarType   trace_element; // for level >0 this is the matrix trace

  // for non-scalar matrices (either row or column level>0) in both MAR and
  // SPR modes, we need to keep track of the four submatrices in order
  // [0 1; 2 3]. The row level of these submatrices is one less than the
  // row level of the matrix containing it (to a minimum of 0), and the same
  // is true for their column levels.
  int64_t subMatList[4];

}; 


// scalar records
struct larc_smatrix
{
  // UNIQUE IDENTIFICATION since matns_ptr_t and mats_ptr_t can be released
  // and reassigned
  int64_t    packedID; // a unique integer, containing a sequentially-numbered
                       // matrixID as well as certain flag bits

#ifdef MAR
#ifdef STORE_TILE_INDEX
  MAR_tile_index_t *tile;
#endif // #ifdef STORE_TILE_INDEX
#endif // #ifdef MAR
  
  // FLAGS  - are bits inside a half word
  uint16_t      info;        // 0 if nothing in info store, flip bit for each enum type in info_store used
  unsigned int  hold    :13; // incremented/decremented by user, can remove associated matrix only if 0
  //  unsigned int  MARnhbr :1;  // if set, submatrix[0] is the pID of MARrep
  unsigned int  lock    :1;  // if set, do not EVER remove this matrix from matrix store
  unsigned int  iszero  :1;  // is this a zero matrix  (was unint32_t)
  unsigned int  isid    :1;  // is this an identity matrix (was unint32_t)

  // NUMBER OF APPEARS_AS_SUB MATRICES  - uses other half of the word with flags
  uint32_t    appears_as_sub_count; 

  scalarType   scalar_value;  // the scalar value being stored

  // for scalar matrices (level==0), we keep track of the hash nodes that are
  // relevant to this record. In SPR mode, tile_node[0] contains the pointer
  // of the hash node that points to this scalar record.
  // In MAR mode, tile_node[0] contains the pointer to the hash node for the
  // primary tile, and the other three are either NULL pointers or
  // pointers to the hash nodes for the neighbor tiles.
  hash_node_t *tile_node[4];

} ; 

typedef struct larc_nsmatrix larc_nsmatrix_t;
typedef struct larc_smatrix larc_smatrix_t;

typedef larc_nsmatrix_t *matns_ptr_t;
typedef larc_smatrix_t *mats_ptr_t;

#endif // ifndef SWIG


/*****************************************************************************
 *                 The different types of matrices                           *
 *****************************************************************************/
typedef enum matrix_type {
  SCALAR,               // both row_level and col_level are 0
  ROW_VECTOR,           // 1 row, multiple columns (row_level =0)
  COL_VECTOR,           // multiple rows, 1 column (col_level =0)
  MATRIX                // both row_level and col_level are greater than 0
} matrix_type_t;

/* Creates and returns a log directory for input and output from this run */
char * create_log_dir(char * log_name);


// This will print out a message whenever we have a memory allocation problem.
#define ALLOCFAIL() \
	fprintf(stderr,"*** MEMORY ALLOCATION FAIL --> %s : %s %d\n", __FILE__, __func__, __LINE__)

// The default exponent for hashing matrices and operation results
// This can be overridden by the -x command line option for the matrix store and
// the -X command line option for the operations store.
#define DEFAULT_MATRIX_STORE_EXPONENT 26
#define DEFAULT_OP_STORE_EXPONENT 24
#define DEFAULT_MAX_LEVEL 32

// The size of the locality sensitive hashing regions 
// are determined by a parameter regionbitparam that
// is passed in during initialization and then saved
// as a matrix_store structure parameter.
// See more details in the matrix_store.c comments, or
// various documentation and tutorial functions.
//
// This parameter was set by running some tests for
// determining when the regions would be too small, then
// subtracting a few bits from the parameter to make the
// regions somewhat larger than this minimum successful size.
// -> long double, long double complex: the parameter
// should be smaller than or equal to 59 (defaulted below to 56)
// -> 256-bit multiprecision floating point (and 
// indefinite precision rational rounding for e.g.
// square root): the parameter should be smaller than
// or equal to 249 (defaulted below to 245)
// 
#ifdef IS_MP
#define DEFAULT_REGIONBITPARAM 245
#else
#define DEFAULT_REGIONBITPARAM 56
#endif

// After testing to verify, I think we want the code to always run in
// HASH_CHAIN_GROWS_AT_TAIL modes, because zero and identity matrices
// are stored first.
#define HASH_CHAIN_GROWS_AT_TAIL
// #define HASH_CHAIN_GROWS_AT_HEAD


//////////////////////////////////////////////////////
// To allow the size of all the hash chains to be printed to a file
// and to allow the size of the largest hash chain to be calculated
// define HASHSTATS
// but be warned that this will create four arrays of size  2^hash_exponent
//////////////////////////////////////////////////////
// #define HASHSTATS

typedef struct hash_table  // holds the general info on the hash table and chain pointers
{
  uint64_t exponent;         // create uses exponent: nentries = 2^exponent
  uint64_t nentries;         // nentries = number of possible hash chains (num buckets)
  uint64_t mask;	     // will be nentries-1; used to mask indices into table
  uint64_t active_chains;    // number of chains with something in them ??
  uint64_t active_entries;   // number of hash nodes in the entire hash table
  uint64_t inserts;          // number of hash nodes that were ever created
  uint64_t hits;             // # of times that this hash table successfully retrieved a record
  uint64_t misses;           // number of times we looked for nodes without success
  uint64_t record_size;      // the size of an op_store record or matrix record  
  hash_node_t **heads;       // array of pointers to node at head of hash chain
  hash_node_t **tails;       // array of pointers to node at tail of hash chain
#ifdef HASHSTATS
  uint64_t *num_accesses;    // array indexed by hash_val of number of times chain was accessed
  uint64_t *num_nodes;       // array indexed by hash_val of number of entries in chain
  // uint64_t *chain_length;    // the length of each hash chain
#endif
} hash_table_t;


/* The types of operations to store. */
/* We use the X macro trick - the three arguments will become the constant
 * used in the enum, the string it translates to, and the verbose description
 * of the operation including that string. The macro X is defined and undefined
 * locally to return the value we need at that time. */
/* Make sure INVALID_OP remains the last operation if more are added! */

/*!
 * \ingroup larc
 * \brief uses X macro to keep together the enum OPERATIONS and their short and long descriptions
 */
#define OPERATIONS \
        X(SUM, "SUM", "SUM (matrix addition)") \
        X(DIFF, "DIFFERENCE", "DIFFERENCE (matrix subtraction)") \
        X(PRODUCT, "PRODUCT", "PRODUCT (matrix multiplication)") \
        X(QUOTIENT_SCALAR, "QUOTIENT_SCALAR", "QUOTIENT by SCALAR (matrix-scalar elementwise division)") \
        X(MATRIX_OR, "OR", "OR (non-exclusive matrix or)") \
	X(BINARY_PROD, "BINARY_PRODUCT", "BINARY_PRODUCT (binary matrix multiplication)") \
        X(KRONECKER, "KRONECKER", "KRONECKER (tensor product)") \
        X(ENTRYSQUARE, "ENTRYSQUARE", "ENTRYSQUARE (square each value in matrix)") \
        X(JOIN, "JOIN", "JOIN (adjoin matrices side by side)") \
        X(STACK, "STACK", "STACK (adjoin matrices vertically)") \
        X(ADJOINT, "ADJOINT", "ADJOINT (conjugate transpose of a matrix)") \
        X(ZEROCOUNT, "ZEROCOUNT", "ZEROCOUNT (count of zero entries of a matrix)") \
        X(BASIS_CHANGE, "BASIS_CHANGE", "BASIS_CHANGE: calculate B*A*Adjoint(B). (When B is unitary and A is an endomorphism, this is a basis change for A)") \
	X(MAX_ELEMENT, "MAX_ELEMENT", "MAX_ELEMENT: the 1x1 submatrix within a matrix with maximum norm value") \
	X(NORM_VAL, "NORM_VAL", "NORM_VAL: L1, (L2)^2, or L\\infty norm of matrix (depending on arguments)") \
        X(FUNC_A, "FUNCTION_A", "FUNCTION_A (generic function for user)") \
        X(FUNC_B, "FUNCTION_B", "FUNCTION_B (generic function for user)") \
        X(FUNC_C, "FUNCTION_C", "FUNCTION_C (generic function for user)") \
        X(INVALID_OP, "INVALID_OP", "INVALID_OP (undefined)") 
#define X(a,b,c) a,
enum operation_types { OPERATIONS };
#undef X

typedef enum operation_types op_type_t;

struct larc_op_t {
        op_type_t op_type;
        int64_t in1_pID;
        int64_t in2_pID;
        int64_t out_pID;
};

typedef struct larc_op_t *op_ptr_t;



/* The types of information to store. */
/* We use the X macro trick - the three arguments will become the constant
 * used in the enum, the string it translates to, and the verbose description
 *   for now we call the string it translates to, the "info_name"
 *   and the list of these is called the "info_names".
 * of the information type including that string. The macro X is defined and undefined
 * locally to return the value we need at that time. */
/* Make sure INVALID_INFO remains the last operation if more are added! */
/* Note there is an info flag in each larc matrix, with a bit flipped for the enum info_types used */
// in info_store.c it is expected that the info_name is at most 64 characters long

/*!
 * \ingroup larc
 * \brief uses X macro to keep together the enum INFO and their short and long descriptions
 */
#define INFO \
        X(SCALELOG, "SCALELOG", "SCALELOG (scalars in saved matrix are 2^SCALELOG * true values)") \
        X(NUMITERATIONS, "NUMITERATIONS", "NUMITERATIONS (the number of iterations completed)") \
        X(DATE, "DATE", "DATE (the date YYYYMMDD that computation occurred)") \
        X(COMPUTER, "COMPUTER", "COMPUTER (which computational platform was used)") \
        X(INPUTVEC, "INPUTVEC", "INPUTVEC (path to LARCMatrix file for input vector)") \
        X(OTHERINFO, "OTHERINFO", "OTHERINFO (string describing what is in OTHERMATRIX)") \
        X(OTHERMATRIX, "OTHERMATRIX", "OTHERMATRIX (path to LARCMatrix file for something else)") \
        X(COMMENT, "COMMENT", "COMMENT (string of your choosing)") \
        X(LARC_SIZE, "LARC_SIZE", "LARC_SIZE (number of unique quadrant submatrices in matrix quadtree)") \
        X(UNIQ_SCALARS, "UNIQ_SCALARS", "UNIQ_SCALARS (number of unique scalars in matrix)") \
        X(SCALARTYPE, "SCALARTYPE", "SCALARTYPE (the scalar type that LARC was compiled with)") \
        X(REGIONTYPE, "REGIONTYPE", "REGIONTYPE (whether MAR or SPR mode was selected at compile time)") \
        X(REGIONBITPARAM, "REGIONBITPARAM", "REGIONBITPARAM (region bit parameter)") \
        X(ZEROREGIONBITPARAM, "ZEROREGIONBITPARAM", "ZEROREGIONBITPARAM (zero region bit parameter)") \
        X(INVALID_INFO, "INVALID_INFO", "INVALID_INFO (undefined)") 

#define X(a,b,c) a,
enum info_types { INFO };
#undef X

typedef enum info_types info_type_t;

struct larc_info_t {
  int64_t my_pID;  // The matrix, which has associated meta data
  info_type_t info_type;  // the type of meta data that is recorded (see enum above)
  char info_data[256];  // the string containing the meta information that we want to keep
                          // e.g. a path to a LARCMatrix file, a string with an integer, etc.
};

typedef struct larc_info_t *info_ptr_t;



#endif
