//                         larc.h
/******************************************************************
 *                                                                *
 * Copyright (C) 2014, Institute for Defense Analyses             *
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
 * \mainpage High-Level LARC Documentation
 * \section data_sec Data Structure
 *
 * The data structure used by LARC is recursive.
 *
 * \section python_sec Python Interface
 * There is a C base and a Python interface.
 *
 */


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

verbose_type_t VERBOSE;

#ifdef USE_INTEGER
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
// The "ifndef SWIG" keep SWIG from creating access functions for the elements
// of the scalarType structure. SWIG attempts to assign the array pointer which
// represents the multiprecision type to a scalar value, resulting in errors.
// This is not a problem because we should never directly access the scalarType
// from Python; we instead have functions that convert the value to a string
// before it is returned.
#ifndef SWIG
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
#endif // ifndef SWIG
#ifdef USE_MPREAL
  typedef mpfr_t scalarType;
#endif
#ifdef USE_MPCOMPLEX
  typedef mpc_t scalarType;
#endif

// if adding additional types, also see scalarTypeDef in global.h, global.c

/*******************************************************************
 *  Prototypes for scalar operations that are parameter based.     *
 *   must be initialized before preloading matrix store            *
 *   If adding additional scalar operations, update scalars.c      *
 *   especially check_scalarOps().                                 *
 ******************************************************************/

/****************************************************************************
 * These functions ensure that multiprecision types (which usually allocate *
 * memory for structures) are properly created, initialized and freed. When *
 * using standard C scalar types, these functions are usually just wrappers *
 * for C operations.                                                       */

/*!
 * \ingroup larc
 * \brief Allocates memory for a new variable of type scalarType (if needed for that scalarType)
 * \param s_ptr A pointer to scalarType (to hold the new variable)
 */
void     (*sca_init)    (scalarType* s_ptr);
/*!
 * \ingroup larc
 * \brief Deallocates memory for a variable of type scalarType (if needed for that scalarType)
 * \param s_ptr A pointer to scalarType to be freed
 */
void     (*sca_clear)   (scalarType* s_ptr);
/*!
 * \ingroup larc
 * \brief Copies the value of one scalarType to another
 * \param d_ptr A pointer to a scalarType value
 * \param s_ptr A pointer to a scalarType value
 */
void     (*sca_set)     (scalarType* d_ptr, const scalarType s_ptr);
/*!
 * \ingroup larc
 * \brief Converts a string into the appropriate scalarType and stores that value
 * \param s_ptr The scalarType value to be written to the store
 * \param input_string The string-form representation to the scalarType value
 */
void     (*sca_set_str) (scalarType* s_ptr, const char *input_str);
/*!
 * \ingroup larc
 * \brief Converts two double precision floats into the appropriate scalarType and stores that value
 * \param real_val The real part of the value
 * \param imag_val The imaginary part of the value (if the scalarType is not a complex type, this must be zero)
 */
void     (*sca_set_2ldoubles) (scalarType*, long double real_val, long double imag_val);
/*!
 * \ingroup larc
 * \brief Converts a scalarType into string format and returns that string; this routine uses malloc()
 * \param s A scalarType variable
 * \return The strning form of the value of s
 */
char    *(*sca_get_str) (const scalarType s);
/*!
 * \ingroup larc
 * \brief Hashes a scalarType to a 64-bit integer according to a predefined function
 * \param s A scalarType variable
 * \param exp The log base 2 of the size of the hash table
 * \return The hash value of the input scalar
 */
uint64_t (*sca_hash)    (const scalarType s, uint64_t exp); 
/*!
 * \ingroup larc
 * \brief Adds two scalarType variables and returns the result
 * \param s A pointer to the sum of the two variables
 * \param s1 The first of two scalarType variables to be added
 * \param s2 The second of two scalarType variables to be added
 */
void     (*sca_add)   (scalarType* s, const scalarType s1, const scalarType s2);
/*!
 * \ingroup larc
 * \brief Multiplies two scalarType variables and returns the result
 * \param s A pointer to the product of the two variables
 * \param s1 The first of two scalarType variables to be multiplied
 * \param s2 The second of two scalarType variables to be multiplied
 */
void     (*sca_mult)  (scalarType* p, const scalarType s1, const scalarType s2);
/*!
 * \ingroup larc
 * \brief Divides one scalarType variable by another and returns the result
 * \param s A pointer to the result of division 
 * \param n The scalarType variable to be the numerator
 * \param d The scalarType variables to be the denominator
 */
void     (*sca_divide)  (scalarType* s, const scalarType n, const scalarType d);
/*!
 * \ingroup larc
 * \brief Returns the neighborhood approximation to a scalarType value
 * \param output The neighborhood representative of the input scalarType
 * \param input The input scalarType
 */
void     (*sca_nbhd_approx)(scalarType* output, const scalarType input);
/**********************************************************************
 * sca_cmp: The function’s return value is equal to zero, if the two  *
 * scalars are equal; and nonzero otherwise. Negative vs positive     *
 * nonzero values aren't required for LARC yet.                       *
 * If they were, return a positive value if first is greater than     *
 * second, and a negative value if first is less than second.         *
***********************************************************************/
/*!
 * \ingroup larc
 * \brief Compares two scalarType values
 * \param a The first scalarType variable
 * \param b The second scalarType variable
 * \return -1/0/1 if the first is less than/equal to/greater than the second
 */
int      (*sca_cmp)     (const scalarType a, const scalarType b); //
/*!
 * \ingroup larc
 * \brief Checks for exact equality of two scalarType variables
 * \param a The first scalarType variable
 * \param b The second scalarType variable
 * \return 1 if the the two scalarType values are equal, 0 otherwise
 * 
 */
int      (*sca_eq)   (const scalarType a, const scalarType b); //
/*!
 * \ingroup larc
 * \brief Checks for equality of the neighborhood approximations for two scalarType variables
 * \param a The first scalarType variable
 * \param b The second scalarType variable
 * \return 1 if the the two scalarType values have the same neighborhood approximation, 0 otherwise
 */
int      (*sca_eq_approx)   (const scalarType a, const scalarType b); //
/*!
 * \ingroup larc
 * \brief Returns the sqrt of the scalarType variable
 * \param a A pointer to the sqrt of the input scalarType variable
 * \param b The input scalarType variable
 */
void     (*sca_sqrt)    (scalarType* a, const scalarType b);
/*!
 * \ingroup larc
 * \brief Returns the norm of the scalarType variable
 * \param a A pointer to the norm of the input scalarType variable
 * \param b The input scalarType variable
 */
void     (*sca_norm)    (scalarType* a, const scalarType b);
/*!
 * \ingroup larc
 * \brief Returns the conjugate of the scalarType variable
 * \param a A pointer to the conjugate of the input scalarType variable
 * \param b The input scalarType variable
 */
void     (*sca_conj)    (scalarType* a, const scalarType b);
/*!
 * \ingroup larc
 * \brief Tests whether the (complex) input has an imaginary part
 * \param a The input scalarType variable
 * \return 1 if the scalarType variable has no imaginary part, 0 otherwise
 */
int      (*sca_is_real)  (const scalarType a);

/*  These types are new names for uintptr_t.
    In this way, the type is actually a pointer to the record
    disguised as an integer.  The record_ptr_t is used in 
    hash chain function definitions, but when hash functions
    are called it will be either a mat_ptr_t if we are dealing
    with matrix records, or a op_ptr_t if we are dealing with
    op_store records.   
    Note that structures for these records are local to matrix_store.c
    or op_store.c respectively, so they can not be cast outside
    of these files.
*/
typedef uintptr_t* record_ptr_t; // is used in definition of hash functions
                                // to stand for either mat_ptr_t or op_ptr_t

/* We currently use -1 to indicate that nothing is in the matrix store 
   (see variable mat_level_t largest), so
   it can't be changed to a uint without rethinking this. */
typedef int32_t mat_level_t;

// The structure of a LARC matrix must be protected from SWIG. The reason is
// that SWIG tries to create access functions for 'scalarType trace_element'
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
 *         row, level, flags+ref_ct, matrixID, trace (2 if complex, 1 if real) *
 *         subs (4)                                                          * 
 *                                                                           *
 *****************************************************************************/
struct larc_matrix
{
  // LEVELS   matrix is 2^row_level by 2^col_level
  mat_level_t row_level;         
  mat_level_t col_level;         

  // FLAGS  - are bits inside a half word
  uint32_t    info    :16;  // 0 if nothing in info store, flip bit for each enum type in info_store used
  uint32_t    hold    :13;  // incremented/decremented by user, can remove associated matrix only if 0
  uint32_t    lock    :1;  // if set, do not EVER remove this matrix from matrix store
  uint32_t    iszero  :1;  // is this a zero matrix  (was unint32_t)
  uint32_t    isid    :1;  // is this an identity matrix (was unint32_t)

  // NUMBER OF APPEARS_AS_SUB MATRICES  - uses other half of the word with flags
  uint32_t    appears_as_sub_count; 

  // UNIQUE IDENTIFICATION since mat_ptr_t can be released and reassigned
  int64_t    matrixID; // a unique integer

  // SCALAR / TRACE
  //   trace_element is scalar element if levels are 0 
  //   trace_element is the trace for square matrices with level > 0
  scalarType   trace_element; // scalar value if row_level == col_level == 0

  // RECURSIVELY DEFINE USING CHILDREN QUADRANT SUBMATRICES
  struct larc_matrix *submatrix[4];  // for level > 0, the one level down
                             // submatrices in order [0 1 ; 2 3]
                             // are specified as indices in storage_table
} ; 

typedef struct larc_matrix larc_matrix_t;

typedef larc_matrix_t *mat_ptr_t;

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

// After testing to verify, I think we want the code to always run in
// HASH_CHAIN_GROWS_AT_TAIL modes, because zero and identity matrices
// are stored first.
#define HASH_CHAIN_GROWS_AT_TAIL
// #define HASH_CHAIN_GROWS_AT_HEAD

// SIGHASH determines the number of significant bits used
// when checking to see if two scalar values are identical.
// this must be used in both the store retrieval and in the
// comparision function.
// ZEROBITTHRESH determines a threshold below which small numbers
// are set to zero, e.g. if < 2^(-ZEROBITTHRESH) then set to zero.
// These values are in bits 
// CHANGE AS OF 22Aug2017: default is for both these parameters to be set so
// that there is no locality approximation for the hashing
// *there are LDBL_MANT_DIG bits in the significand, so this keeps them all
// *the smallest positive nonzero number representable in double precision is
//  2^{-1074} which is 0x1 as represented as a long unsigned int, so any
//  nonnegative number smaller than this is by definition zero
//  In python code we often use rnd_sig_bits to specify the SIGHASH
//      and use trunc_to_zero_bits to specify the ZEROBITTHRESH
#define SIGHASH_DEFAULT (LDBL_MANT_DIG-2) // formerly 15
#define ZEROBITTHRESH_DEFAULT (LDBL_MANT_DIG-2) // formerly 20

//////////////////////////////////////////////////////
// To allow the size of all the hash chains to be printed to a file
// and to allow the size of the largest hash chain to be calculated
// define HASHSTATS
// but be warned that this will create four arrays of size  2^hash_exponent
//////////////////////////////////////////////////////
// #define HASHSTATS

typedef struct hash_node
{
  record_ptr_t record_ptr;  // either a mat_ptr_t or op_ptr_t
  uint64_t hits;            // number of times node's record has been retrieved
  struct hash_node *next;   // pointer to next hash node in chain or to "head"
  struct hash_node *prev;   // pointer to last hash node in chain or to "tail"
} hash_node_t;

typedef struct hash_table  // holds the general info on the hash table and chain pointers
{
  uint64_t exponent;         // create uses exponent: nentries = 2^exponent
  uint64_t nentries;         // nentries = number of possible hash chains (num buckets)
  uint64_t active_chains;    // number of chains with something in them ??
  uint64_t active_entries;   // number of hash nodes in the entire hash table
  uint64_t inserts;          // number of hash nodes that were ever created
  uint64_t hits;             // number of times a hash node was successfully retrieved
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
        int64_t in1_matID;
        int64_t in2_matID;
        int64_t out_matID;
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
#define INFO \
        X(SCALELOG, "SCALELOG", "SCALELOG (scalars in saved matrix are 2^SCALELOG * true values)") \
        X(NUMROUNDS, "NUMROUNDS", "NUMROUNDS (the number of rounds completed)") \
        X(DATE, "DATE", "DATE (the date YYYYMMDD that computation occurred)") \
        X(COMPUTER, "COMPUTER", "COMPUTER (which computational platform was used)") \
        X(INPUTVEC, "INPUTVEC", "INPUTVEC (path to LARCMatrix file for input vector)") \
        X(OTHERINFO, "OTHERINFO", "OTHERINFO (string describing what is in OTHERMATRIX)") \
        X(OTHERMATRIX, "OTHERMATRIX", "OTHERMATRIX (path to LARCMatrix file for something else)") \
        X(COMMENT, "COMMENT", "COMMENT (string of your choosing)") \
        X(INVALID_INFO, "INVALID_INFO", "INVALID_INFO (undefined)") 

#define X(a,b,c) a,
enum info_types { INFO };
#undef X

typedef enum info_types info_type_t;

struct larc_info_t {
  int64_t my_matID;  // The matrix, which has associated meta data
  info_type_t info_type;  // the type of meta data that is recorded (see enum above)
  char info_data[256];  // the string containing the meta information that we want to keep
                          // e.g. a path to a LARCMatrix file, a string with an integer, etc.
};

typedef struct larc_info_t *info_ptr_t;



#endif
