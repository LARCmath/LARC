//                         larc.h
/******************************************************************
 *                                                                *
 * Copyright 2014, Institute for Defense Analyses                 *
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
 * POC: Jennifer Zito <jszito@super.org>                          *
 * Please contact the POC before disseminating this code.         *
 *                                                                *
 *****************************************************************/

#ifndef LARC_LARC_H
#define LARC_LARC_H

#include <stddef.h> // size_t
#include <stdint.h> // int8_t, int32_t, uint32_t, uint64_t, uintptr_t
#include <float.h> // DBL_MANT_DIG
#ifdef __cplusplus
#include <complex> // complex numbers
#else
#include <complex.h> // complex numbers
#endif

#include "type.h"

#ifdef USE_INTEGER
  typedef int64_t ScalarType;
#endif
#ifdef USE_COMPLEX
#ifdef __cplusplus
  typedef std::complex<double> ScalarType;
#else
  typedef complex ScalarType;
#endif
#endif
#ifdef USE_REAL
  typedef double ScalarType;
#endif
/*!
 * \cond
 */

// if adding additional types, also see scalarTypeDef in global.h, global.c


/*  These types are new names for uintptr_t.
    In this way, the type is actually a pointer to the record
    disguised as an integer.  The record_ptr_t is used in 
    hash chain function definitions, but when hash functions
    are called it will be either a mat_add_t if we are dealing
    with matrix records, or a op_ptr_t if we are dealing with
    op_store records.   
    Note that structures for these records are local to matrix_store.c
    or op_store.c respectively, so they can not be cast outside
    of these files.
*/
typedef uintptr_t* record_ptr_t; // is used in definition of hash functions
                                // to stand for either mat_add_t or op_add_t

/* We currently use -1 to indicate that nothing is in the matrix store 
   (see variable mat_level_t largest), so
   it can't be changed to a uint without rethinking this. */
typedef int32_t mat_level_t;

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
 *         subs (4), adjoint                                                 *
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

  // UNIQUE IDENTIFICATION since mat_add_t can be released and reassigned
  int64_t    matrixID; // a unique integer

  // SCALAR / TRACE
  //   trace_element is scalar element if levels are 0 
  //   trace_element is the trace for square matrices with level > 0
  ScalarType   trace_element; // scalar value if row_level == col_level == 0

  // RECURSIVELY DEFINE USING CHILDREN QUADRANT SUBMATRICES
  struct larc_matrix *submatrix[4];  // for level > 0, the one level down
                             // submatrices in order [0 1 ; 2 3]
                             // are specified as indices in storage_table
  // POINTER TO ADJOINT OR TRANSPOSE, WHEN KNOWN
  struct larc_matrix *adjoint;       // adjoint of this matrix

} ; 

typedef struct larc_matrix larc_matrix_t;

typedef larc_matrix_t *mat_add_t;



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
	printf("*** MEMORY ALLOCATION FAIL --> %s : %s %d\n", __FILE__, __func__, __LINE__)

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
// *there are DBL_MANT_DIG bits in the significand, so this keeps them all
// *the smallest positive nonzero number representable in double precision is
//  2^{-1074} which is 0x1 as represented as a long unsigned int, so any
//  nonnegative number smaller than this is by definition zero
#define SIGHASH_DEFAULT DBL_MANT_DIG-2 // formerly 15
#define ZEROBITTHRESH_DEFAULT DBL_MANT_DIG-2 // formerly 20

//////////////////////////////////////////////////////
// To allow the size of all the hash chains to be printed to a file
// and to allow the size of the largest hash chain to be calculated
// define HASHSTATS
// but be warned that this will create four arrays of size  2^hash_exponent
//////////////////////////////////////////////////////
// #define HASHSTATS

typedef struct hash_node
{
  record_ptr_t record_ptr;  // either a mat_add_t or op_add_t
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
        X(KRONECKER, "KRONECKER", "KRONECKER (tensor product)") \
        X(ENTRYSQUARE, "ENTRYSQUARE", "ENTRYSQUARE (square each value in matrix)") \
        X(JOIN, "JOIN", "JOIN (adjoin matrices side by side)") \
        X(STACK, "STACK", "STACK (adjoin matrices vertically)") \
        X(ADJOINT, "ADJOINT", "ADJOINT (conjugate transpose of a matrix)") \
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

typedef struct larc_op_t *op_add_t;



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
        X(INPUTVEC, "INPUTVEC", "INPUTVEC (path to json file for input vector)") \
        X(OTHERINFO, "OTHERINFO", "OTHERINFO (string describing what is in OTHERMATRIX)") \
        X(OTHERMATRIX, "OTHERMATRIX", "OTHERMATRIX (path to json file for something else)") \
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
                          // e.g. a path to a json file, a string with an integer, etc.
};

typedef struct larc_info_t *info_add_t;

/*!
 * \endcond
 */


#endif
