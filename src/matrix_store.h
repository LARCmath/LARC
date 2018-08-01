//                   matrix_store.h
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


#ifndef LARC_MATRIX_STORE_H
#define LARC_MATRIX_STORE_H

#ifdef __cplusplus
extern "C" {
#endif // #ifdef __cplusplus

#include <stddef.h>      // size_t
#include "larc.h"

/************************************************************************
 *               MACRO FUNCTIONS acting on matrix id (mat_add_t)       *
 ***********************************************************************/
// This value for matrix ptr is also for the fake panels in row and col vectors
#define MATRIX_PTR_INVALID ((mat_add_t)0)
// This value for matrix id 
#define MATRIX_ID_INVALID ((int64_t)-1)
// Boolean functions
#define matrix_is_invalid(m_ptr) (m_ptr == MATRIX_PTR_INVALID)
#define matrix_is_id(m_ptr) (m_ptr && m_ptr->isid)
#define matrix_is_zero(m_ptr) (m_ptr && m_ptr->iszero)
#define matrix_has_lock(m_ptr) (m_ptr && m_ptr->lock)
#define matrix_has_hold(m_ptr) (m_ptr && m_ptr->hold)
// Functions retrieving values inside a larc_matrix structure

/*!
 * \ingroup larc
 *
 * \brief Return the trace of a matrix.
 *
 * \param m_ptr A handle to the matrix for which we want the trace.
 *
 * \return The trace of the matrix.
 */
#define matrix_trace(m_ptr) (m_ptr->trace_element)
#define matrix_row_level(m_ptr) (m_ptr->row_level)
#define matrix_col_level(m_ptr) (m_ptr->col_level)
#define matrix_sub(m_ptr,s) (m_ptr->submatrix[s])
#define matrix_appears_as_sub_count(m_ptr) (m_ptr->appears_as_sub_count)
// the adjoint_field can contain MATRIX_PTR_INVALID
#define matrix_adjoint_field(m_ptr) (m_ptr->adjoint)

/******************************************************************
 *         ACCESSORS acting on matrix id (mat_add_t)              *
 *****************************************************************/

/*!
 * \ingroup larc
 *
 * \brief Get a value within a matrix.
 *
 * \param m A handle to the matrix from which to retrieve a value.
 * \param row The row within the matrix.
 * \param col The column within the matrix.
 *
 * \return The value within the matrix m at the location specified by row and col.
 */
ScalarType matrix_get_value ( mat_add_t m , mat_level_t row , mat_level_t col );

/*!
 * \ingroup larc
 *
 * \brief Set a value within a matrix.
 *
 * \param m   A handle to the matrix for which we set a value.
 * \param row The row within the matrix.
 * \param col The column within the matrix.
 * \param v   The value to set.
 *
 * \return A new matrix that is the same as the input matrix, except
 * with the value at the specified row and column set to the value v.
 */
mat_add_t  matrix_set_value ( mat_add_t m , mat_level_t row , mat_level_t col , ScalarType v );

/*************************************************************************
 *         INLINE FUNCTIONS acting on matrix id (mat_add_t)              *
 ************************************************************************/
// functions return enum matrix_type_t 
inline static matrix_type_t matrix_type_from_levels(mat_level_t row_level, mat_level_t col_level)
{
  if ((row_level == 0) && (col_level == 0)) { return SCALAR;}
  if (row_level == 0) { return ROW_VECTOR;}
  if (col_level == 0) { return COL_VECTOR;}
  return MATRIX;
}

inline static matrix_type_t matrix_type(mat_add_t m_ptr) 
{
  mat_level_t row_level = matrix_row_level(m_ptr);
  mat_level_t col_level = matrix_col_level(m_ptr);
  return matrix_type_from_levels(row_level,col_level);
}

inline static int64_t get_matrixID_from_ptr(mat_add_t m_ptr) {
  return ((m_ptr==MATRIX_PTR_INVALID) ? MATRIX_ID_INVALID : m_ptr->matrixID);
}

/* retrieve the maximum level of row and columns allowed in the store */
mat_level_t maximum_level(void);

/* retrieve the locality-approximation parameters, number of significant bits to round */
int get_sighash(void);

/* retreive the locality-approximation parameters, distance from zero to collapse to zero */
double get_zerorealthresh(void);

/* Prints a summary of the matrix store usage */
void matrix_store_report(char *outfilepath);

/* Returns 1 on successfully creating the matrix store, 0 on failure */
int create_matrix_store(size_t exponent, mat_level_t max_level, int sighash, int zerobitthresh);

int lock_matrix(mat_add_t m_ptr);
int set_hold_matrix(mat_add_t m_ptr);
int release_hold_matrix(mat_add_t m_ptr);

/* python interface versions using matrixID */
int set_hold_matrix_from_matrixID(int64_t m_mID);
int release_hold_matrix_from_matrixID(int64_t m_mID);

/* Preload all matrices Jenny and Steve have decided to preload */
int preload_matrix_store(void);

/* access to commonly used preloaded matrices */
int64_t get_zero_matrixID(mat_level_t row_level, mat_level_t col_level);
mat_add_t get_zero_matrix_ptr(mat_level_t row_level, mat_level_t col_level);
int64_t get_identity_matrixID(mat_level_t level);
mat_add_t get_identity_matrix_ptr(mat_level_t level);
int64_t get_iHadamard_matrixID(mat_level_t level);
mat_add_t get_iHadamard_matrix_ptr(mat_level_t level);


int matrix_appears_as_sub_count_increment(mat_add_t id);
int matrix_appears_as_sub_count_decrement(mat_add_t id);

// python interface versions - needed to avoid passing void * type
int64_t matrix_get_matrixID_from_panel(int64_t A_mID, int64_t B_mID,
        int64_t C_mID, int64_t D_mID, mat_level_t row_level, 
        mat_level_t col_level);
int64_t matrix_get_matrixID_from_scalar(ScalarType val);

mat_add_t matrix_get_ptr_scalar(ScalarType scalar);
mat_add_t matrix_get_ptr_panel(mat_add_t panel[4], mat_level_t row_level, mat_level_t col_level);

/* Returns the number of scalars currently in memory */
size_t matrix_store_scalarCount(void);

/* Returns the number of matrices (excluding scalars) in memory */
size_t matrix_store_matrixCount(void);

int64_t matrix_adjoint_matrixID(int64_t m_mID);
/*!
 * \ingroup larc
 *
 * \param id A handle to the matrix for which we would like the adjoint.
 *
 * \return A handle to the adjoint matrix.
 */
mat_add_t matrix_adjoint(mat_add_t id);

/* Useful for some linear algebra applications */
// int preload_iHadamard(void);

/* The total number of matrices created and stored in  */
/* the matrix store (the current number of stored     */
/* may be less if there was matrix store cleaning).     */
uint64_t  num_matrices_created(void);

uint64_t num_matrices_in_store (void);

mat_add_t mat_ptr_from_matrixID(int64_t matrixID, const char* arg_no,
        const char* calling_routine, int print_flag);

int matrix_store_info_to_file(uint64_t start, uint64_t end, 
			      char *outfilepath, char *comment);

int matrix_hash_chain_info_to_file(uint64_t hash, char *outfilepath, char *comment);

int matrix_hash_chain_info_to_screen(uint64_t hash, char *comment);

int remove_matrix_from_mat_store_by_matrixID (int64_t  m_mID);

int remove_matrix_from_mat_store (mat_add_t  m_ptr);

int table_mat_ptr_by_matrixID_Remove_entry(mat_add_t m_ptr);

ScalarType locality_approx(ScalarType scalar);

int64_t matrix_hashID_from_matrixID(int64_t m_mID);

int clean_matrix_store();

int clean_matrix_hash_chain(uint64_t hash);

int64_t get_diag_hist(uint64_t i);

int64_t matrix_hash_chain_length(uint64_t i);

inline static ScalarType matrix_trace_matrixID(int64_t m_ID)
{
  mat_add_t m_ptr = mat_ptr_from_matrixID(m_ID,"only",__func__,0);
  return matrix_trace(m_ptr);
}

#ifdef HASHSTATS
/* prints hash table statistics for matrix store, creating files with 
   accesses and nodes for each bucket (hash_val).  It also produces
   as standard hash_report to stdout or file.
   This works only if HASHSTATS is defined.     */
void matrix_hashstats(char *accesses_file,  char *nodes_file, char *report_file);
#endif

#ifdef __cplusplus
}
#endif // #ifdef __cplusplus

#endif    
