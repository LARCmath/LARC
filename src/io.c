//                            io.c 
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


#include <inttypes.h> //for printing uint64_t and int64_t
#include <curses.h>
#include <ctype.h>
#include "larc.h"
#include "global.h"
#include "organize.h"
#include "io.h"
#include <errno.h>

/*!
 * \file io.c
 * \brief This file contains the routines which transfer data between disk
 * and LARC.
 *
 * Uncompressed data can be read in either row-major-matrix or Matrix Market
 * Exchange formats. Compressed data in LARCMatrix (json) format can be both
 * read and written; sparse data can also be saved to disk as a list of 
 * nonzero locations. There are also functions which change LARCMatrix data as
 * it is being read or written, which may be of use on occasion. Small
 * matrices can also be written to file or screen in uncompressed format.
 */

// This structure groups together all the information related to the 
// quadtree-like structure used to change (row,col,val) information into
// a properly compressed LARCMatrix.
struct dataForQLS
{
        int64_t **quad;
        mat_level_t *internal_node_level;
        int64_t num_internal_nodes;
        mat_level_t *leaf_level;
};


// This structure groups together information about a sparse matrix which
// was read in as (row,col,val) triplets. In MatrixMarketExchange format,
// this type of data is labeled as 'coordinate'.
struct dataByCoordinate
{
        unsigned *row, *col; 
        scalarType *val;
};

// This structure groups together all the information needed for
// writing a LARCmatrix json output file, or a unique list of scalars.
typedef struct uniq_submatrix_array_struct
{
  char *array;
  size_t   sizeArray;
  size_t   larcSizeCounter;
  size_t   numUniqScalars;
} usas_t;


/*!
 * \ingroup larc
 * \brief allocates memory for a uniq_submatrix_array structure, used when writing a LARCmatrix json output file or a unique list of scalars
 *
 * \param packedID the matrix (either to be written or to have its unique scalars extracted and written)
 * \result A pointer to the structure allocated
 */
static usas_t *create_usasPTR(int64_t packedID) {
    if (packedID == MATRIX_ID_INVALID)
    {
        fprintf(stderr,"in %s, invalid matrixID passed\n",__func__);
        return NULL;
    }
    usas_t *usasPTR = calloc(1, sizeof(usas_t));
    //    size_t   larcSizeCounter= 0;
    //    size_t   numUniqScalars = 0;
    int64_t matrixID = MID_FROM_PID(packedID);
    usasPTR->sizeArray = matrixID+1;
    usasPTR->array = calloc(matrixID+1, sizeof(char));
    if (usasPTR->array == NULL) {
        ALLOCFAIL();
        return 0;
    }
    return(usasPTR);
}
  
/*!
 * \ingroup larc
 * \brief frees memory for a uniq_submatrix_array structure
 *
 * \param usasPTR The pointer to the structure to be freed
 */
static void free_usasPTR(usas_t *usasPTR) {
    free(usasPTR->array);
    free(usasPTR);
}

/*!
 * \ingroup larc
 * \brief This routine recursively fills the uniq_submatrix_array structure for the given packedID.
 *
 * \param usasPTR the (previously allocated) structure which will hold the data
 * \param packedID the matrix (either to be written or to have its unique scalars extracted and written)
 */
static void fill_usasPTR(usas_t *usasPTR, int64_t packedID)
{
    int local_verbose = 0;

    // In row and column vectors skip over nonfilled quadrants
    if (packedID == MATRIX_ID_INVALID) return;
  
    int64_t matrixID = MID_FROM_PID(packedID);

    if (local_verbose)
    {
	printf("In function %s the matrixID is %ld\n",__func__,matrixID);
    }
    
    // simple error checking
    if (matrixID >= usasPTR->sizeArray) {
        printf("This can never happen, in function %s!\n", __func__);
        exit(0);
    }

    if (local_verbose)
    {
	printf("Checking to see if array[matrixID] is nonzero \n");
    }
    
    // if we have seen this matrix before
    if (usasPTR->array[matrixID] != 0) return;
    
    if (local_verbose)
    {
	printf("The value in array[matrixID] is %d\n",usasPTR->array[matrixID]);
    }
    
    // have not seen this submatrix yet
    // (usas_PTR.array[packedID] == 0)
    // scalar matrix
    if (IS_SCALAR(packedID))
    {
        usasPTR->array[matrixID] = 1;
        ++usasPTR->numUniqScalars;
    }
    // nonscalar matrix	
    else
    {
        usasPTR->array[matrixID] = 2;

        // recursively check the children of m_ID
        matns_ptr_t m_PTR = (matns_ptr_t)get_recordPTR_from_pID(packedID,
                       "", __func__, 0);
        for (int i=0; i<4; ++i)
            fill_usasPTR(usasPTR,m_PTR->subMatList[i]);
    }
    // larcSizeCounter will be the number of unique quadrant submatrices
    ++usasPTR->larcSizeCounter;

    if (local_verbose)
    {
        printf(" \n");
    }

   return;
}

    
/****************************************************************
 *              recursive_row_major_list_to_store
 *    jszito note: John Gilbert taught me this neat recursive method;
 *                 he says comes from the old Fortran BLAS
 *                 (basic linear algebra subprograms).
 *    It takes a dense m by n matrix which is given as a 1-dim list of 
 *    all its entries in row major order (reading each row in order):
 *                   a00 a01 ... a0(n-1)
 *                   a10 a11 ... a1(n-1)
 *                    ...  
 *                   a(m-1)0 a(m-1)1 ... a(m-1)(n-1)
 *         where m = 2^row_level, n= 2^col_level, and 
 *         orig_num_cols = n
 *         aij are scalars of whatever scalarType we are using.
 *    We specify a submatrix by giving its top left corner aij
 *    and its dimension (and the orginal number of columns).
 *
 *    This routine works recursively on the 4 submatrixIDs.
 *    It loads the matrix into the matrix store c using
 *    get_pID_from_array_of_four_sub_pIDs
 *    which returns the matrixID of the matrix in the store.
 *
 *    Any quadrant submatrix in the original matrix can be specified by 
 *    passing: 
 *       the size of the original matrix orig_num_cols 
 *             which gives the number of elements in each row
 *       the current_row_level of the submatrix which is
 *           one less than the last matrix, unless we 
 *           are in the ROW_VECTOR case.
 *       the current_col_level of the submatrix which is
 *           one less than the last matrix, unless we 
 *           are in the COL_VECTOR case.
 *       the position in dense_mat at which the submatrix starts
 *             e.g. in a 4 by 4 matrix, submatrix[0] starts at 0
 *                  submatrix[1] starts at 2, submatrix[2] starts at 8,
 *                  and submatrix[3] starts at 10, and
 *       the subarray of dense_mat that starts at that index.
 *
 *    The recursive nature of the call takes care of everything else.
 *    When get_pID_from_array_of_four_sub_pIDs has the MatrixID of
 *    each of the four quadrant submatrices
 *    (in the vector cases two of these are MATRIX_ID_INVALID),
 *    then it can be used to find the mat_ptr of a matrix.
 *    
 **********************************************************************/

/*!
 * \ingroup larc
 *
 * \brief Recursive routine used to add a matrix to the MatrixStore
 *
 * This recursive routine requires the original number of columns
 * (the paramater dim_whole) at lower levels of the recursion.
 *
 * \param dense_mat A pointer to scalarType; the matrix to be stored
 * \param current_row_level The base2 log of the number of rows in the current matrix
 * \param current_col_level The base2 log of the number of columns in the current matrix
 * \param orig_num_cols The number of columns of the matrix originally passed to the recursion
 * \return The packedID of the stored matrix.
 */
static int64_t recursive_row_major_list_to_store(scalarType  * dense_mat, 
                    mat_level_t current_row_level, 
                    mat_level_t current_col_level, 
                    int64_t     orig_num_cols)
{
#ifdef DEBUG_IO_C
  printf("--> %s : %s: %d\n", __FILE__, __func__, __LINE__);
#endif // #ifdef DEBUG_IO_C

  int verbose = 0;

  // SCALAR CASE
  if ((current_row_level == 0) && (current_col_level == 0)){
    return (get_scalarPTR_for_scalarVal(*dense_mat))->packedID;
  }

  // Need this for remaining cases
  int64_t subMatList[4];  // indices of submatrices


  // ROW VECTOR
  if (current_row_level == 0) 
  {
    uint64_t new_num_cols = (uint64_t)1 << (current_col_level-1);   // divide current number by two

    // top corner
    subMatList[0] = recursive_row_major_list_to_store(dense_mat, 0,
                    current_col_level-1,orig_num_cols);

    // step right dim_new cols
    subMatList[1] = recursive_row_major_list_to_store(dense_mat + new_num_cols,
                    0, current_col_level-1,orig_num_cols);

    // step down dim_new rows
    subMatList[2] = MATRIX_ID_INVALID;

    // step down and over dim_new rows and cols
    subMatList[3] = MATRIX_ID_INVALID;

#ifdef DEBUG_IO_C
  printf("<--%s : %s: %d\n", __FILE__, __func__,  __LINE__);
#endif // #ifdef DEBUG_IO_C 

  }

  // COLUMN VECTOR
  else if (current_col_level == 0)
  {
    uint64_t new_num_rows = (uint64_t)1 << (current_row_level-1);   // divide current number by two

    // top corner
    subMatList[0] = recursive_row_major_list_to_store(dense_mat,
                    current_row_level-1, 0, orig_num_cols);

    // step right dim_new cols
    subMatList[1] = MATRIX_ID_INVALID;

    // step down dim_new rows
    subMatList[2] = recursive_row_major_list_to_store(dense_mat +
           new_num_rows*orig_num_cols, current_row_level-1, 0, orig_num_cols); 

    // step down and over dim_new rows and cols
    subMatList[3] = MATRIX_ID_INVALID;

#ifdef DEBUG_IO_C
  printf("<--%s : %s: %d\n", __FILE__, __func__,  __LINE__);
#endif // #ifdef DEBUG_IO_C

  }

  // MATRIX CASE
  else 
  {
    uint64_t new_num_cols = (uint64_t)1 << (current_col_level-1);   // divide current number by two
    uint64_t new_num_rows = (uint64_t)1 << (current_row_level-1);   // divide current number by two

    if (verbose) {
      printf("\nIn recursive call with orig_num_cols = %"PRId64"\n",orig_num_cols);
      printf("  The four subsubMatLists will be at row level %"PRId32", col level %"PRId32", ",
                current_row_level-1, current_col_level-1);
      printf("and will have offsets: \n");
      printf("     subMatList[0]: 0\n");
      printf("     subMatList[1]: new_num_cols = %" PRIu64 "\n",new_num_cols);
      printf("     subMatList[2]: new_num_rows*orig_num_cols = %" PRIu64 "\n",
                new_num_rows*orig_num_cols);
      printf("     subMatList[3]: new_num_rows*orig_num_cols+new_num_cols = %" PRIu64 "\n",
                new_num_rows*orig_num_cols+new_num_cols);
    }  // end verbose

    // top corner
    subMatList[0] = recursive_row_major_list_to_store(dense_mat,
           current_row_level-1, current_col_level-1, orig_num_cols);

    // step right dim_new cols
    subMatList[1] = recursive_row_major_list_to_store(dense_mat +
           new_num_cols, current_row_level-1, current_col_level-1,
           orig_num_cols); 

    // step down dim_new rows
    subMatList[2] = recursive_row_major_list_to_store(dense_mat +
           new_num_rows*orig_num_cols, current_row_level-1, 
           current_col_level-1, orig_num_cols); 

    // step down and over dim_new rows and cols
    subMatList[3] = recursive_row_major_list_to_store(dense_mat +
           new_num_rows*orig_num_cols+new_num_cols,
           current_row_level-1, current_col_level-1, orig_num_cols); 
#ifdef DEBUG_IO_C
  printf("<--%s : %s: %d\n", __FILE__, __func__, __LINE__);
#endif // #ifdef DEBUG_IO_C 

  }
  return get_pID_from_array_of_four_sub_pIDs(subMatList,
         current_row_level, current_col_level);
}

// TWO ROUTINES FOR READING IN FROM ROW_MAJOR_LIST FORMAT:
//
// This routine works with multiprecision data.
// int64_t row_major_list_to_store(char **dense_mat, 
//                    mat_level_t        current_row_level, 
//                    mat_level_t        current_col_level, 
//                    int64_t        orig_num_cols 
//                    )
//
// This routine does NOT work with multiprecision data.
// int64_t read_row_major_matrix_from_file(char * file_path)
//
// row_major_list_to_store():
// This wrapper routine takes as its first argument an array of strings,
// which will be converted to the proper scalarType value before calling 
// the core function recursive_row_major_list_to_store().
int64_t row_major_list_to_store(char **dense_mat, 
                    mat_level_t        current_row_level, 
                    mat_level_t        current_col_level, 
                    int64_t        orig_num_cols 
                    )
{

  // allocate sufficient space to store the matrix "dense" in memory
  uint64_t N = ((uint64_t)1 << current_row_level) 
               * ((uint64_t)1 << current_col_level);
  scalarType *dense = malloc(N*sizeof(scalarType));
  if (dense == NULL) {
    fprintf(stderr,"ERROR in %s, line %d: Out of memory.\n", __func__, __LINE__);
    exit(1);
  }

  // convert each string value in dense_mat to the proper scalar type
  // and put it into our "dense" matrix
  for (int64_t i = 0; i < N; i++){
    sca_init(&(dense[i]));
    sca_set_str(&(dense[i]), dense_mat[i]);
  }

  // call recursive_row_major_list_to_store to put "dense" into the MatrixStore
  // (or possibly ScalarStore)
  int64_t C_pID = recursive_row_major_list_to_store(dense,
                  current_row_level, current_col_level, orig_num_cols); 

  // clear the memory used to store "dense"
  for (int64_t i = 0; i < N; i++)
    sca_clear(&(dense[i]));
  free(dense);

  return C_pID;
}

/****************************************************************
 *  python interface: read_row_major_matrix_from_file
 ****************************************************************/
int64_t read_row_major_matrix_from_file(char * file_path)
{

#ifdef DEBUG_IO_C
  printf("--> %s : %s: %d\n", __FILE__, __func__, __LINE__);
#endif // #ifdef DEBUG_IO_C 
  int verbose = 0;
  int error_code = 0;

  if (verbose) printf("In %s\n", __func__);

  FILE *fp;
  if((fp = fopen(file_path, "r")) == NULL) {
    fprintf(stderr,"No such file: %s\n",file_path);
    exit(1);    
  }
  if (verbose) printf("Opened the data file %s for reading:\n",file_path);

  mat_level_t row_level;
  mat_level_t col_level;
  // Read in the level
  int ret = fscanf(fp, "%u %u", &row_level, &col_level);
  if (ret == 2) {
    if (verbose) printf("\nReading the row and column levels: %u %u\n", (mat_level_t)row_level, (mat_level_t)col_level);
    fflush(stdout);
  } 
  else if(errno != 0) {
    perror("fscanf:");
    error_code = 1;
    fprintf(stderr,"ERROR scanning levels in first line of file %s\n", file_path);
  }
  else if(ret == EOF) {
    error_code = 2;
    fprintf(stderr,"ERROR: file %s ended early in %s\n", file_path,  __func__);
  }
  else {
    fprintf(stderr,"ERROR while reading levels in %s (routine %s).\n", file_path, __func__);
    error_code = 2;
  }
  if (error_code) {
#ifdef DEBUG_IO_C
    printf("<-- %s : fail: %d\n", __FILE__, __LINE__);
#endif // #ifdef DEBUG_IO_C 
    return MATRIX_ID_INVALID; 
  } 

  // if the file is open but nothing is left to read print EOF
  if(feof(fp)){
    fputs("EOF",stderr);
    fprintf(stderr,"ERROR: File stopped earlier than expected.\n");
    fclose(fp);
#ifdef DEBUG_IO_C
    printf("<-- %s : fail: %d\n", __FILE__, __LINE__);
#endif // #ifdef DEBUG_IO_C 
    return MATRIX_ID_INVALID; 
  }

  // successfully read the level
  if (verbose)   printf("  levels are %d %d\n",(int)row_level, (int)col_level);

  // calculate size of the array
  uint64_t array_size = (uint64_t)1 << (row_level+col_level);
  if (verbose)  printf("  array length is %" PRIu64 "\n", array_size);

  // Now we will handle the entries of the matrix
  scalarType * matrix_array;

  // Allocate space for the array of scalarType
  if (NULL==(matrix_array = calloc(array_size,sizeof(scalarType)) )) {
    error_code = 1;  // unable to allocate space for array
    fprintf(stderr,"ERROR: Unable to allocate space for array.\n");
#ifdef DEBUG_IO_C
    printf("<-- %s : fail: %d\n", __FILE__, __LINE__);
#endif // #ifdef DEBUG_IO_C 
    return MATRIX_ID_INVALID; 
  } 

  // Read the array of numbers 
  int max_size = 1000;
  char a[max_size];
  int64_t i;
  for (i=0; i<array_size; ++i) {
     // format for numbers must not include whitespace!
     int ret = fscanf(fp,"%s\n",a);
     if (ret == EOF) {
        fprintf(stderr,"ERROR: fscanf() failed to read a line (i=%" PRId64 ")",i);
        break;
     }
     else if (errno != 0) {
        perror("fscanf");
        fprintf(stderr,"(i=%" PRId64 ")\n",i);
        break;
     }
     else if (ret == max_size) {
        fprintf(stderr,"ERROR: line read by fscanf() may exceed ");
        fprintf(stderr,"%d characters (i=%" PRId64 ")\n",max_size,i);
     }
     sca_init(&(matrix_array[i]));
     sca_set_str(&(matrix_array[i]), a);
  }

  // if the file is open but nothing is left to read print EOF
  if(feof(fp)){
    fputs("EOF row major reader.\n",stderr);
  }

  fclose(fp);

  int64_t m_pID;
  if (array_size == i) {
#ifdef DEBUG_IO_C
    printf("<-- %s : %s success: %d\n", __FILE__, __func__,__LINE__);
#endif // #ifdef DEBUG_IO_C 
    m_pID = recursive_row_major_list_to_store(matrix_array,
                row_level,col_level,(uint64_t)1<<col_level);
  }
  else {
    if (verbose) printf("Failed to read row major matrix from file.\n");
#ifdef DEBUG_IO_C
    printf("<-- %s : %s fail: %d\n", __FILE__, __func__, __LINE__);
#endif // #ifdef DEBUG_IO_C 
    m_pID = MATRIX_ID_INVALID;
  }

  // clean up any space allocated by sca_init
  for (int64_t j=0; j<i; ++j) sca_clear(&(matrix_array[j]));
  // clean up scalarType array
  free(matrix_array);

  return m_pID;
}


// TWO NAIVE MATRIX PRINTING ROUTINES
// * these routines print out the full matrix with no compression, and must
// * be used carefully...
// void print_naive(int64_t A_pID)
// void fprint_naive(int64_t A_pID, char * filename)

/*******************************************************************
*                       (naive) print_naive       *
* accepts a matrixID for matrix to print                           *
* Converts to a pointer before calling                                                                *
********************************************************************/
void print_naive(int64_t A_pID)
{
  // The following checks to ensure that m_pID
  // has valid record pointer (eg, not RECORD_PTR_INVALID)
  record_ptr_t m_Rptr;
  check_validity_one_input(A_pID, __func__, &m_Rptr);

  char *entry_string;

  if (IS_SCALAR(A_pID))
  {
    mats_ptr_t A_ptr = (mats_ptr_t)m_Rptr;
    entry_string = sca_get_readable_approx_str(A_ptr->scalar_value);
    printf("%s\n", entry_string);
    free(entry_string);
    return;
  }

  // get the matrix pointer from the matrixID, and see if still in store
  matns_ptr_t A_ptr = (matns_ptr_t)m_Rptr;

  int i,j;
  int row_level =  A_ptr->row_level;
  int col_level =  A_ptr->col_level;

  // set the maximum level for printing a matrix
  // these numbers are intentionally set larger than what we think
  // would be reasonable for a matrix printed to a terminal window
#ifdef IS_MP
  int max_print_level = 6;
#elif defined(IS_COMPLEX)
  int max_print_level = 7;
#else
  int max_print_level = 10;
#endif

  if (row_level>max_print_level || col_level>max_print_level)
  {
        fprintf(stderr,"in %s, matrix has row level %d",__func__,row_level);
        fprintf(stderr," and column level %d,\n",col_level);
        fprintf(stderr,"which is too large to naive print to screen. ");
        fprintf(stderr,"If you must do this anyway,\nuse write_naive routine");
        fprintf(stderr," with 'stdout' as output filename\n");
        return;
  }

  // before continuing, calculate number of characters to be printed
  // and check that the total number that will be printed is less than
  // MAX_TO_PRINT. (this could double the time required for printing,
  // we decided that was OK)
#define MAX_TO_PRINT 100000
  uint64_t c_counter, max_columns, max_chars;

  unsigned row_dim = 1 << row_level;
  unsigned col_dim = 1 << col_level;
  max_chars = max_columns = 0;
  for (i = 0;i<row_dim;i++) {
    c_counter = 0;
    for (j = 0;j<col_dim; j++) {
      mats_ptr_t e_ptr = get_scalarPTR_from_pID_and_coords(i, j, A_pID);
      entry_string = sca_get_readable_approx_str(e_ptr->scalar_value);
      c_counter += strlen(entry_string)+1;
      free(entry_string);
    }
    c_counter += 2;
    max_chars += c_counter;
    max_columns = MAX(max_columns,c_counter);
    if (max_chars>MAX_TO_PRINT) break;
  }

  // limit naive printing to reasonable size
  if (max_chars>MAX_TO_PRINT)
  {
        fprintf(stderr,"in %s, matrix has dimensions",__func__);
        fprintf(stderr," %u x %u\n",row_dim,col_dim);
        fprintf(stderr,"and requires more than %" PRIu64 " printed characters,\n",
                max_chars);
        fprintf(stderr,"which is too large to naive print to screen. ");
        fprintf(stderr,"If you must do this anyway,\nuse write_naive routine");
        fprintf(stderr," with 'stdout' as output filename\n");
        return;
  }

  if (VERBOSE>BASIC) {
    printf("%s given matrix which, if printed, would require\n", __func__);
    printf("%" PRIu64 " characters total, with a maximum of \n", max_chars);
    printf("%" PRIu64 " characters per row of the matrix\n", max_columns);
  }

  // determine size of terminal using curses library
  int maxrow, maxcol;
  initscr();
  getmaxyx(stdscr,maxrow,maxcol);
  endwin();
  // printf("number of rows in terminal window is %d\n",maxrow);
  // printf("number of columns in terminal window is %d\n",maxcol);

  c_counter = maxrow; // avoids compiler warning that maxrow set, unused
  for(i = 0; i<row_dim; i++) {
    // find the value to be printed in string form
    mats_ptr_t e_ptr = get_scalarPTR_from_pID_and_coords(i, 0, A_pID);
    entry_string = sca_get_readable_approx_str(e_ptr->scalar_value);
    // always print first element in row without checking line length
    printf("%s ", entry_string);
    c_counter = strlen(entry_string)+1;
    free(entry_string);

    for(j = 1; j<col_dim; j++) {
      // find the value to be printed in string form
      mats_ptr_t e_ptr = get_scalarPTR_from_pID_and_coords(i, j, A_pID);
      entry_string = sca_get_readable_approx_str(e_ptr->scalar_value);
      // Determine if printing this string will cause the number of characters
      // printed for this line to be greater than the terminal width, and if it
      // would, start a new indented line.
      // For a matrix whose entries produce long strings, this can result in
      // one matrix element per line.
      c_counter += strlen(entry_string)+1; // includes space
      if (c_counter>maxcol) { 
        printf("\n      ");
        c_counter = 6+strlen(entry_string)+1;
        }
      printf("%s ", entry_string);
      free(entry_string);
    }
    // always put newline at end of matrix row
    printf("\n");
  }
  fflush(stdout);
}


/********************************************************************
* Python interface version of (naive) print_matrix_to_file          *
*   this could be implemented less naively ...                      *
********************************************************************/
void fprint_naive(int64_t A_pID, char * filename)
{
  // The following checks to ensure that m_pID
  // has valid record pointer (eg, not RECORD_PTR_INVALID)
  record_ptr_t m_Rptr;
  check_validity_one_input(A_pID, __func__, &m_Rptr);

  FILE *f;
  // remember strcmp returns 0 on equality
  if (strcmp(filename,"stdout")) f = fopen(filename, "w");
  else f = stdout;

  char *entry_string;

  if (IS_SCALAR(A_pID))
  {
     mats_ptr_t A_ptr = (mats_ptr_t)m_Rptr;
     entry_string = sca_get_readable_approx_str(A_ptr->scalar_value);
     fprintf(f,"%s\n", entry_string);
     free(entry_string);
     if (strcmp(filename,"stdout")) fclose(f);
     return;
  }

  // get the matrix pointer from the packedID, and see if still in store 
  matns_ptr_t A_ptr = (matns_ptr_t)m_Rptr;

  int i,j;
  int row_level =  A_ptr->row_level;
  int col_level =  A_ptr->col_level;
  unsigned row_dim = 1 << row_level;
  unsigned col_dim = 1 << col_level;

  for(i = 0;i<row_dim;i++){
    for(j = 0;j<col_dim; j++){
      mats_ptr_t e_ptr = get_scalarPTR_from_pID_and_coords(i, j, A_pID);
      char *entry_string = sca_get_readable_approx_str(e_ptr->scalar_value);
      fprintf(f,"%s ", entry_string);
      free(entry_string);
    }
    fprintf(f," \n");
  }
  if (strcmp(filename,"stdout")) fclose(f);
}

// SOMEWHAT LESS NAIVE PRINTING ROUTINE
// void fprint_matrix_nonzeros(int64_t A_pID, char* filename)
// * there is no limit to the number of nonzero values printed...

/********************************************************************
* (Python interface) fprint_matrix_nonzeros            *
********************************************************************/
void fprint_matrix_nonzeros(int64_t A_pID, char * filename)
{
  // The following checks to ensure that A_pID
  // has valid record pointer (eg, not RECORD_PTR_INVALID)
  record_ptr_t m_Rptr;
  check_validity_one_input(A_pID, __func__, &m_Rptr);

  if (IS_SCALAR(A_pID))
  {
     // get the matrix pointer from the matrixID, and see if still in store 
     mats_ptr_t A_ptr = (mats_ptr_t)m_Rptr;

     FILE *f = fopen(filename, "w");
     if (!(A_ptr->iszero))
     {
         char *entry_string = sca_get_readable_approx_str(A_ptr->scalar_value);
         fprintf(f, "%s\t(0,0)\n", entry_string);
         free(entry_string);
     }
     fclose(f);
     return;
  }

  // get the matrix pointer from the matrixID, and see if still in store 
  matns_ptr_t A_ptr = (matns_ptr_t)m_Rptr;

  FILE *f = fopen(filename, "w");
  mat_level_t row_level =  A_ptr->row_level;
  mat_level_t col_level =  A_ptr->col_level;
  uint64_t row_dim = (uint64_t)1 << row_level;
  uint64_t col_dim = (uint64_t)1 << col_level;

  for(uint64_t i = 0; i < row_dim; i++){
    for(uint64_t j = 0; j < col_dim; j++){
      mats_ptr_t e_ptr = get_scalarPTR_from_pID_and_coords(i, j, A_pID);
      if (!(e_ptr->iszero))
      {
        char *entry_string = sca_get_readable_approx_str(e_ptr->scalar_value);
        fprintf(f, "%s\t(%" PRIu64 ",%" PRIu64 ")\n", entry_string, i, j);
        free(entry_string);
      }
    }
  }
  fclose(f);
}


// A COMPARISON BETWEEN THE MATRICES STORED IN TWO COMPRESSED JSON FILES
// This (apparently) is in io.c since data is read from files...
// The routine works with multiprecision types and can be called from the
// python interface. It puts both matrices
// into the matrix store, and returns 0 if the matrices are not equal or the
// matrixID if they are equal.
int64_t equal_matrices_in_larcMatrix_files(char *path1, char *path2)
{
  if (!strcmp(path1,path2)) /* paths are equal */
    {
      fprintf(stderr,"Warning %s : %s fail: line %d, path1 and path2 are the same.\n",
              __FILE__, __func__, __LINE__);
    }
  int64_t m1_pID = read_larcMatrixFile(path1);
  int64_t m2_pID = read_larcMatrixFile(path2);
  if (( m1_pID==MATRIX_ID_INVALID) || (m2_pID==MATRIX_ID_INVALID)) {
      fprintf(stderr,"Warning %s : %s fail: line %d, invalid matrix\n",
              __FILE__, __func__, __LINE__);
      return(0);
  }
  if (m1_pID == m2_pID) {
    // If packedID were zero, could not tell the values were equal. However,
    // at initialization the first matrixID is assigned to 0, a scalar, and
    // with the is_scalar bit set the assigned packedID will be nonzero,
    // so packedID==0 cannot happen.
    return m1_pID;
  }
  else return 0;
}

// ROUTINES WHICH EXECUTE READING AND WRITING FROM/TO JSON COMPRESSED FORMAT
// 
// THE WRITING ROUTINES:
// int fprint_larcMatrixFile(int64_t m_pID, char *path)
// static void update_infoStore_automatic_entries(int64_t m_pID)
// static int write_infoStore_to_larcMatrix_file(int64_t m_pID, FILE *f)

/*!
 * \ingroup larc
 * \brief Puts automatic entries into info store.
 *
 * The automatic entries are SCALARTYPE, REGIONTYPE, REGIONBITPARAM, and ZEROREGIONBITPARAM.
 * \param m_pID The packedID of the matrix which will have its infoStore entries added or updated
 */
static void update_infoStore_automatic_entries(int64_t m_pID)
{
    char buf[256];
    int err_flag;

    // write the scalar type to the InfoStore
    snprintf(buf, sizeof(buf), get_string_scalarType());

    err_flag =  info_set(SCALARTYPE, m_pID, buf);
    if (err_flag) {
        fprintf(stderr,"Error writing scalar type to info store in %s\n", __func__);
        exit(1);
    }

    // write the region type to the InfoStore
#ifdef MAR
    snprintf(buf, sizeof(buf), "MAR");
#else
    snprintf(buf, sizeof(buf), "SPR");
#endif  // #ifdef MAR

    err_flag =  info_set(REGIONTYPE, m_pID, buf);
    if (err_flag) {
        fprintf(stderr,"Error writing region type to info store in %s\n", __func__);
        exit(1);
    }

    // write the region bit param to the InfoStore
    snprintf(buf, sizeof(buf), "%d", get_regionbitparam());

    err_flag =  info_set(REGIONBITPARAM, m_pID, buf);
    if (err_flag) {
        fprintf(stderr,"Error writing region bit param to info store in %s\n", __func__);
        exit(1);
    }

    // write the zero region bit param to the InfoStore
    snprintf(buf, sizeof(buf), "%d", get_zeroregionbitparam());

    err_flag =  info_set(ZEROREGIONBITPARAM, m_pID, buf);
    if (err_flag) {
        fprintf(stderr,"Error writing zero region bit param to info store in %s\n", __func__);
        exit(1);
    }
} // end update_infoStore_automatic_entries

/* This function should retrieve all information out of the InfoStore
   for the given matrix associated with m_pID and write it to the
   json file associated with the output file pointer f.  */

/*!
 * \ingroup larc
 * \brief Writes all InfoStore information about a given matrix to a json file
 * \param m_pID PackedID of the matrix to be stored in compressed json format
 * \param f File pointer for the json file where the metadata will be written
 * \return 1 on success
 */
static int write_infoStore_to_larcMatrix_file(int64_t m_pID, FILE *f)
{
  int verbose = 0;

  if (!f)
    return -1;
  if (m_pID==MATRIX_ID_INVALID)
    return -1;

  enum info_types i;
  int was_info = 0;
  char* info_data;
  char* info_name;
  for (i=0;i<INVALID_INFO;++i) {
    info_data = info_get(i,m_pID);
    if (verbose) printf("info_data is %s, info_enum is %d\n",info_data,i);
    // if info_data has an entry then the following test is true
    if (strcmp(info_data,"")) {
      if (was_info == 0)  {
        fprintf(f, "\n  \"info\":{\n"); 
        was_info = 1;
      }
      info_name = return_info_name(i);
      fprintf(f, "   \"%s\":\"%s\",\n",info_name,info_data);
    }
    free(info_data);
  }
  if (was_info) {
    // fprintf(f, "      \"end\":0 },");
    fprintf(f, "      \"end\":\"\" },");
  }
  return(0);
}


// THE READING ROUTINES:
// int64_t read_larcMatrixFile(char *path)
// int64_t read_anyFormat_larcMatrixFile(char *path)
// int64_t read_legacy_larcMatrixFile(char *path)
//
// Similar to but less complicated than the writing case. These functions read
// string-formatted data from a compressed json matrix file, converts strings
// to the appropriate scalarType, then puts the values into the matrixStore.
//
// The default version of the read function, read_larcMatrixFile, assumes that
// the format of the input file is line-by-line
// equivalent to the fprint_larcMatrixFile() output format, whereas the
// _anyFormat_ version only assumes that the input file is valid JSON and
// contains the correct fields.
//
// The _legacy version of the read function will correctly interpret
// old-style compressed json in which the data fields are C integers,
// doubles or complex numbers rather than string representations of these
// numbers. (The switch to strings was necessary to enable GMP multiprecision.)
// For obvious reasons, there are no legacy writing routines.

// Simple JSON reader, expects the file read to be in the format produced
// by fprint_larcMatrixFile().
int64_t read_larcMatrixFile(char *path)
{
    int verbose = 0;
    FILE *f = fopen(path, "r");
    if (verbose) {
        printf("VERBOSE: %s: path = %s\n",__func__, path);
    }
    if (NULL == f) {
        fprintf(stderr,"%s: no file found at\n\t%s:\n\texiting\n", __func__, path);
        exit(1);
    }

    // These are the fields we are looking for.
    // const char *str_matmax = "\"matrixID_max\":";
    // const char *str_matID = "\"matid\":";
    // const char *str_info = "\"info\":";
    // const char *str_table = "\"table\":";
    // const char *str_end = "\"end\":";

#if defined(IS_RATIONAL) || defined(USE_MPINTEGER)
    char *linebuf, *val_str;
#else
    const int buffer_length = 4096;
    char linebuf[buffer_length];
    char val_str[1024];
#endif
    uint64_t max_matrixID;
    uint64_t matid = 0;
    int64_t *mapID = NULL;
    long long int info_seek = -1;
    int info_count = 0;
    char string2compare[6];

    // Here we loop over the input file line by line.
    int state = 0;
    long long int linenum = 0;
    while (1) {
#if defined(IS_RATIONAL) || defined(USE_MPINTEGER)
        if (-1 == fscanf(f, " %m[^\n] ", &linebuf))
#else
        if (NULL == fgets(linebuf, buffer_length, f))
#endif
        {
              fprintf(stderr,"%s: finished with EOF\n",__func__);
              break;
        }
        ++linenum;
        if (verbose) printf("linebuf = %s\n",linebuf);
        if (verbose) printf("linenum = %lld\n",linenum);
        if (1 == linenum) {
            // The first line of the file just contains an open brace.
        } else if (2 == linenum) {
            // The second line of the file contains the max matrixID.
            if (1==sscanf(linebuf, " \" matrixID_max \" : %" PRIu64 ,
			&max_matrixID) ) {
                if (verbose) {
                    printf("VERBOSE: %s: max matrixID = %" PRIu64 " on line %lld\n", __func__, max_matrixID, linenum);
                }
                // check that we can allocate an array this large
                // (i.e., total number of bytes does not exceed SIZE_MAX)
                uint64_t size_check = SIZE_MAX/(max_matrixID + 1)/sizeof(int64_t);
                if (size_check==0) {
                    fprintf(stderr,"Error is %s: matrixID_max too large - try\n", __func__);
                    fprintf(stderr,"renumbering the matrix_IDs in %s\n", path);
                    fprintf(stderr,"to reduce this value to something more reasonable.\n");
                    exit(1);
                }
                // the mapID array is indexed by the matrixIDs found in the input file
                // but contains the packedIDs assigned when added to the MatrixStore
                mapID = (int64_t *) calloc(max_matrixID+1, sizeof(int64_t));
                if (mapID == NULL){
                    fprintf(stderr,"Error in %s: failed to allocate array.\n", __func__);
                    exit(1);
                }
            } else {
                fprintf(stderr,"%s(%d): Error on line %lld of %s.\nexiting\n", __func__, __LINE__, linenum, path);
                exit(1);
            }
        } else if (3 == linenum) {
            // The third line of the file contains the matID.
            if (1==sscanf(linebuf, " \" matid \" : %" PRIu64 , &matid) ) {
                if (verbose) {
                    printf("VERBOSE: %s: matid = %" PRIu64 " on line %lld\n", __func__, matid, linenum);
                 }
            } else {
                fprintf(stderr,"%s(%d): Error on line %lld of %s.\nexiting\n", __func__, __LINE__, linenum, path);
                exit(1);
            }
        } else if (4 == linenum) {
            // The fourth line of the file is either "info" or "table".
            if (1!=sscanf(linebuf, " \" %5[a-z] \" : \{", string2compare)) {
                fprintf(stderr,"%s(%d): Error on line %lld of %s.\nexiting\n", __func__, __LINE__, linenum, path);
                exit(1);
            }
	    if (!strcmp(string2compare,"info")) {
                state = 1;  // for "info"
                info_seek = ftell(f);
                if (verbose) {
                    printf("VERBOSE: %s: Found \"info\" on line %lld\n", __func__, linenum);
		    printf("\tinfo_seek set to %lld, state to 1\n", info_seek);
                }
            } else if (!strcmp(string2compare,"table")) {
                state = 3;  // for "table"
                if (verbose) {
                    printf("VERBOSE: %s: Found \"table\" on line %lld\n", __func__, linenum);
                }
            }
        } else if (1 == state) {
            // This line is an entry under "info".
            sscanf(linebuf, " \" %5[a-z] \"",string2compare);
            if (!strcmp(string2compare,"end")) {
                state = 2;  // next line must be "table"
                if (verbose) {
                    printf("VERBOSE: reached end of info, state = 2\n");
                }
            } else {
                // Skip over "info" entries for now.  We will process them later.
                if (verbose) printf("skipping %s\n",linebuf);
                ++info_count;
            }
        } else if (2 == state) {
            // line is just after the "info" section, so it must be "table"
            sscanf(linebuf, " \" %5[a-z] \" ",string2compare);
            if (!strcmp(string2compare,"table")) {
                state = 3;
                if (verbose) {
                    printf("VERBOSE: %s: Found \"table\" on line %lld\n", __func__, linenum);
		    printf("state now 3\n");
                }
            } else {
                fprintf(stderr,"%s(%d): Error on line %lld of %s.\nexiting\n", __func__, __LINE__, linenum, path);
                exit(1);
            }
        } else if (3 == state) {
            // This line is an entry under "table".
            sscanf(linebuf, " \" %5[a-z] \" ",string2compare);
            if (!strcmp(string2compare,"end")) {
                state = 4; // all subsequent lines ignored
            } else {
                // Scan the matrix index and the row and column levels.
                uint64_t m_index;
                int row_level, col_level;
                int rv = sscanf(linebuf, " \" %" SCNu64 " \" : [ %d, %d,",
                                &m_index, &row_level, &col_level);
                if (3 != rv) {
                    fprintf(stderr,"%s(%d): Error on line %lld of %s.\nexiting\n",
                                   __func__, __LINE__, linenum, path);
                    exit(1);
                }

                if ((0 == row_level) && (0 == col_level)) {
                    // Scalar case.  Only the scalar string should be remaining.
                    int rv = sscanf(linebuf,
#if defined(IS_RATIONAL) || defined(USE_MPINTEGER)
		" \" %*u \" : [ %*d, %*d, \" %m[^]\"] \" ]", &val_str
#else
		" \" %*u \" : [ %*d, %*d, \" %1023[^]\"] \" ]", val_str
#endif
                );
		    if (rv != 1) {
                        fprintf(stderr,"%s(%d): Error on line %lld of %s.\nexiting\n",
                                       __func__, __LINE__, linenum, path);
                        exit(1);
                    }
                    mapID[m_index] = get_valID_from_valString(val_str);
                    if (verbose) {
                        printf("VERBOSE: %s: Scalar m_index = %" PRIu64 " for \"%s\""
                               " maps to %" PRId64 " on line %lld\n",
                               __func__, m_index, val_str, mapID[m_index], linenum);
                    } // end verbose
#if defined(IS_RATIONAL) || defined(USE_MPINTEGER)
                    free(val_str);
#endif
                } else {
                    // Non-scalar case.  The 4 panel indices should be remaining.
                    int64_t panel_index[4];
                    int rv = sscanf(linebuf,
" \" %*u \" : [ %*d, %*d, %" SCNd64 ", %" SCNd64 ", %" SCNd64 ", %" SCNd64 " ]",
		panel_index, panel_index+1, panel_index+2, panel_index+3);
                    if (4 != rv) {
                        fprintf(stderr,"%s(%d): Error on line %lld of %s.\nexiting\n",
                                       __func__, __LINE__, linenum, path);
                        exit(1);
                    }
                    int64_t content_pID[4];
                    for (int k = 0; k < 4; ++k) {
                       // The matrix ID -1 is used whenever we have no matrix,
                       // e.g. in row and col place keepers
                       if (-1 == panel_index[k]) {
                           content_pID[k] = MATRIX_ID_INVALID;
                       }
                       else {
                         content_pID[k] = mapID[panel_index[k]];
                       }
                    }
                    mapID[m_index] = get_pID_from_array_of_four_sub_pIDs(
                         content_pID, row_level, col_level);
                    if (verbose) {
                        printf("VERBOSE: %s: Non-scalar m_index = %" PRIu64
                               " for row_level = %d, col_level = %d,"
                               " panels [%" PRId64 ", %" PRId64 ", %" PRId64 ", %" PRId64 "]"
                               " maps to %" PRId64 " on line %lld\n",
                               __func__, m_index, row_level, col_level,
                               panel_index[0], panel_index[1], panel_index[2], panel_index[3],
                               mapID[m_index], linenum);
                    }
                } // end of non-scalar case
            } // end of "table" entry processing
        } // end of state 3 section 
#if defined(IS_RATIONAL) || defined(USE_MPINTEGER)
        free(linebuf);
#endif
    }  // end of while loop over lines in input file

    int64_t packedID =  mapID[matid];
    matns_ptr_t ret = (matns_ptr_t) get_recordPTR_from_pID(packedID,"",__func__,0);
    printf("\nRead %d,%d-level matrix from json file %s\n",
           ret->row_level, ret->col_level, path);
    printf ("  stored matrix now has address %p", ret);
    printf ("  and packedID %" PRId64 "\n", packedID);
    free(mapID);

    const char *scalar_type_string = get_string_scalarType();

    if (-1 != info_seek) {
        if (verbose) {
            printf("VERBOSE: %s: Starting \"info\" entry processing for %d entries.\n",
                   __func__, info_count);
        }
        int rv = fseek(f, info_seek, SEEK_SET);
        if (rv < 0) {
            fprintf(stderr,"%s(%d): Error from file seek of %s.\nexiting\n",
                           __func__, __LINE__, path);
            exit(1);
        }
        for (int j = 0; j < info_count; ++j) {
#if defined(IS_RATIONAL) || defined(USE_MPINTEGER)
            if (1 != fscanf(f, " %m[^\n] ", &linebuf))
#else
            if (NULL == fgets(linebuf, buffer_length, f))
#endif
            {
                fprintf(stderr,"%s(%d): Error re-reading \"info\" entry #%d of %s.\nexiting\n",
                               __func__, __LINE__, j+1, path);
                exit(1);
            }
            if (verbose) {
                printf("VERBOSE: %s: Processing \"info\" entry #%d: %s",
                       __func__, j+1, linebuf);
            }

            char info_name[256];
            char info_data[4096];
	    sscanf(linebuf, " \" %255[^\"] \" : \" %4095[^\"] \",",
		info_name, info_data);
            info_type_t info_type = get_info_type_from_string_name(info_name);
            if (SCALARTYPE == info_type) {
                if (0 != strcmp(info_data, scalar_type_string)) {
                    fprintf(stderr,"Warning %s : %s : line %d, "
                            "Scalar type is different than when larcMatrix file was written.\n",
                            __FILE__, __func__, __LINE__);
                    fprintf(stderr,"Scalar type of larcMatrix file was: %s.\n", info_data);
                    fprintf(stderr,"Current scalar type is: %s.\n", scalar_type_string);
                }
            }
            int fail = info_set(info_type, packedID, info_data);
            if (fail) {
                fprintf(stderr,"panic in %s\n",__func__);
                exit(1);
            }

            if (verbose) {
                printf("VERBOSE: %s: For packedID %" PRId64 ", added key"
                       " %s to InfoStore with value \"%s\"\n",
                       __func__, packedID, info_name, info_data);
            }
#if defined(IS_RATIONAL) || defined(USE_MPINTEGER)
            free(linebuf);
#endif
        }  // end of loop through "info" entries
    }  // end of test for "info" section

    fclose(f);
    return packedID;
}

int64_t read_anyFormat_larcMatrixFile(char *path)
{
// This function is used to read a compressed matrix stored in our json format
// into the matrix store, optionally changing it to a different matrix by 
// applying a function to all scalar values in the matrix. When the function
// passed is sca_set(), the data is left unchanged.

  int verbose = 0;
  FILE *f = fopen(path, "r");
  // printf("  path = %s\n",path);
  if (f == NULL)
  {
    fprintf(stderr,"%s: no file found at\n\t%s:\n\texiting\n", __func__, path);
    exit(1);
  }
  json_t *j = j_parse_file(f);
  // printf("  parsed json\n");

  uint64_t max_matrixID = j_lookup_num64(j, "matrixID_max");
  // printf("  max matrixID is %ld\n", max_matrixID);

  // loop index for reading info and table
  uint64_t i;

  json_t *t = j_key_lookup(j, "table");

  // check that we can allocate an array this large
  // (i.e., total number of bytes does not exceed SIZE_MAX)
  i = SIZE_MAX/(max_matrixID + 1)/sizeof(int64_t);
  if (i==0) {
      fprintf(stderr,"Error is %s: matrixID_max too large - try\n", __func__);
      fprintf(stderr,"renumbering the matrix_IDs in %s\n", path);
      fprintf(stderr,"to reduce this value to something more reasonable.\n");
      exit(1);
  }

  // the mapID array is indexed by the matrixIDs found in the input file
  // but contains the packedIDs assigned when added to the MatrixStore
  int64_t *mapID = calloc(max_matrixID+1, sizeof(int64_t));
  if (mapID == NULL){
    fprintf(stderr,"Error in %s: failed to allocate array.\n", __func__);
    exit(1);
  }
  
  int64_t len = j_key_count(t);

  if (verbose) {printf("In %s, before the table line reader\n", __func__);}
  
  for (i=0; i<len; i++)
    {
      json_t *p = j_key_index(t, i);
      // index == current matrixID
      int64_t index = (int64_t) atoll(p->name); // maybe atol suffices? 
//      if (index>=max_matrixID)
//      {
//         fprintf(stderr,"index > max_matrixID! exiting\n"); exit(1);
//      }
      if (verbose) {printf("For i=%"PRId64", and index=%"PRId64"\n",i,index);}
      if (j_is_array(p))
	{
	  int len1 = j_array_get_length(p);
          if (verbose) {printf("The length of this line is %d\n",len1);}
	  switch(len1)
	    {
	    case 3:			/* expect a scalar link */
	      {
                if (verbose) {printf("In scalar case \n");}
		int row_level = j_get_num64(j_array_index(p, 0));
		int col_level = j_get_num64(j_array_index(p, 1));
                if (!j_is_string(j_array_index(p,2))){
                 fprintf(stderr,"ERROR in %s:\nExpected scalar values to be given as strings.\nIf LARCMatrix file is in 'legacy' format, use utils/canonical_format_json.py\nto convert to LARCMatrix format with scalars saved as strings.\n", __func__);
                    exit(1);
                }
                if (row_level || col_level) {
		  fprintf(stderr,"error in %s:\n\texpected ", __func__ );
                  fprintf(stderr,"scalar, but levels greater than zero!\n");
                }
                const char *val_str = j_get_string(j_array_index(p,2));
                mapID[index] = get_valID_from_valString((char *)val_str); 
                if (verbose) {
		  printf("%"PRId64" %" PRId64 " %d %d %s\n", index,
			 mapID[index], row_level, col_level, val_str);
		}
	      }
              break;
	    case 6:			/* expect a 4-tuple */
	      {
                if (verbose) {printf("In non scalar case \n");}
		int row_level = j_get_num64(j_array_index(p, 0));
		int col_level = j_get_num64(j_array_index(p, 1));
                if (verbose) {printf("  levels are %d %d \n",row_level,col_level);}
                int64_t content_pID[4];
                int64_t temp;
                for (int k=0; k< 4; ++k) {
		   temp = j_get_num64(j_array_index(p, k+2));
                   if (verbose) {printf("temp for k = %d is %" PRId64 " \n",k,temp);}
                   // The matrix ID -1 is used whenever we have no matrix, 
                   // e.g. in row and col place keepers
                   if (temp == -1) {
                       content_pID[k] = MATRIX_ID_INVALID;
                   }
                   else {
                     content_pID[k] = mapID[temp];
		   }
		}
                mapID[index] = get_pID_from_array_of_four_sub_pIDs(
                     content_pID, row_level, col_level);
                if (verbose){
		  printf("ind %" PRId64 " mID %" PRId64 " (%d,%d)",
                        index, mapID[index], row_level, col_level);
		  printf(" [%ld %ld %ld %ld]\n", content_pID[0],
                        content_pID[1], content_pID[2], content_pID[3]);
		}
	      }
	      break; 
	      
	    default:
              fprintf(stderr,"Error in %s(%s): unexpected number of entries per line\n",
                __func__, path);
              fprintf(stderr,"Expected num entries to be 3 or 6, but had %d entries\n",len1);
              exit(1);
	      break;
	    }
	}
    }

  int64_t packedID =  mapID[j_lookup_num64(j, "matid")];
  matns_ptr_t ret = (matns_ptr_t)get_recordPTR_from_pID(packedID,"",__func__,0);
  printf("\nRead %d,%d-level matrix from json file %s\n", 
	 ret->row_level, ret->col_level, path);
  printf ("  stored matrix now has address %p", ret);
  printf ("  and packedID %" PRId64 "\n", packedID);
  free(mapID);

  const char *scalar_type_string = get_string_scalarType();

  if (j_key_exists(j, "info")) {
    if (verbose) printf("Found info key in %s\n", __func__);
    json_t *t_info = j_key_lookup(j, "info");
    int64_t len_info = j_key_count(t_info);
    
    if (verbose) {printf("In %s, before the info line reader\n", __func__);}
    
    for (i=0; i<len_info-1; i++)
      {
	json_t *p = j_key_index(t_info, i);
	const char *info_name = j_get_name(p);
	const char *info_data = j_get_string(p);
	if (verbose) {printf("For i=%"PRId64", and key=%s\n", i, info_name);}
	info_type_t info_type = get_info_type_from_string_name((char *) info_name);
        if (SCALARTYPE == info_type) {
          if (strcmp(info_data, scalar_type_string) != 0) {
            fprintf(stderr,"Warning %s : %s : line %d, Scalar type is different than when larcMatrix file was written.\n",
                    __FILE__, __func__, __LINE__);
            fprintf(stderr,"Scalar type of larcMatrix file was: %s.\n", info_data);
            fprintf(stderr,"Current scalar type is: %s.\n", scalar_type_string);
          }
        }
	int fail =  info_set(info_type, packedID, info_data);
	if (fail) {fprintf(stderr,"panic in %s\n",__func__); exit(1);}
	if (verbose) {
	  printf("Added info field %s to InfoStore\n", info_name);
	  printf("  with value %s\n", info_data);
	}
      }
  }  // end of if j_key_exists for "info"

  fclose(f);
  j_set_null(j); free(j);

  return packedID;
}

// This routine is called by canonical_format_json.py when the python routine
// is used to convert old-style json compressed matrices (with standard C data
// types) into new-style json compressed matrices (with scalars expressed as
// strings)
int64_t read_legacy_larcMatrixFile(char *path)
{
  
#if defined(IS_MP) || defined(IS_BOUNDING)
  // exit program if LARC is compiled for any type which is not a standard C
  // type (multiprecision types, bounding types)
  fprintf(stderr,"ERROR in %s: MP types have no legacy json format.\n", __func__);
  exit(1);
#endif

  int verbose = 0;
  FILE *f = fopen(path, "r");
  // printf("  path = %s\n",path);
  if (f == NULL)
  {
    fprintf(stderr,"%s: no file found at\n\t%s:\n\texiting\n", __func__, path);
    exit(1);
  }
  json_t *j = j_parse_file(f);
  // printf("  parsed json\n");

  uint64_t max_matrixID = j_lookup_num64(j, "matrixID_max");
  // printf("  max matrixID is %ld\n", max_matrixID);

  // loop index for reading info and table
  uint64_t i;

  json_t *t = j_key_lookup(j, "table");

  // check that we can allocate an array this large
  // (i.e., total number of bytes does not exceed SIZE_MAX)
  i = SIZE_MAX/(max_matrixID + 1)/sizeof(int64_t);
  if (i==0) {
      fprintf(stderr,"Error is %s: matrixID_max too large - try\n", __func__);
      fprintf(stderr,"renumbering the matrix_IDs in %s\n", path);
      fprintf(stderr,"to reduce this value to something more reasonable.\n");
      exit(1);
  }

  // the mapID array is indexed by the matrixIDs found in the input file
  // but contains the packedIDs assigned when added to the MatrixStore
  int64_t *mapID = calloc(max_matrixID + 1, sizeof(int64_t));
  if (mapID == NULL){
    fprintf(stderr,"Error in %s: failed to allocate array.\n", __func__);
    exit(1);
  }
  
  int64_t len = j_key_count(t);

  if (verbose) {printf("In %s, before the table line reader\n", __func__);}
  
  for (i=0; i<len; i++)
    {
      json_t *p = j_key_index(t, i);
      // index == current matrixID
      int64_t index = (int64_t) atoll(p->name); // maybe atol suffices? 
      if (verbose) {printf("For i=%"PRId64", and index=%"PRId64"\n",i,index);}
      if (j_is_array(p))
        {
          int len1 = j_array_get_length(p);
          if (verbose) {printf("The length of this line is %d\n",len1);}
          switch(len1)
            {
#ifdef USE_COMPLEX
            case 4:                        /* expect a scalar link */
              {
                if (verbose) {printf("In complex scalar case \n");}
                int row_level = j_get_num64(j_array_index(p, 0));
                int col_level = j_get_num64(j_array_index(p, 1));
                // double complex_val[2];
                // complex_val[0] = j_get_double(j_array_index(p,2));
                // complex_val[1] = j_get_double(j_array_index(p,3));
                // scalarType val = complex_val[0]+I*complex_val[1];
                scalarType val;
                sca_init(&val);
                val =  j_get_double(j_array_index(p,2)) +
                        I*j_get_double(j_array_index(p,3));
                if (row_level || col_level) {
                        fprintf(stderr,"error in %s:\n\texpected complex ", __func__ );
                        fprintf(stderr,"scalar, but levels greater than zero!\n");
                }
                mats_ptr_t s_ptr = get_scalarPTR_for_scalarVal(val);
                mapID[index] = s_ptr->packedID;
                if (verbose) {
                  printf("%"PRId64" %p %" PRId64 " 0 0 %Lg+%Lgi\n", index,
                         s_ptr, mapID[index], creall(val), cimagl(val));
                }
                sca_clear(&val);
              }
              break;
#endif
#ifdef USE_REAL
            case 3:                        /* expect a scalar link */
              {
                if (verbose) {printf("In real scalar case \n");}
                int row_level = j_get_num64(j_array_index(p, 0));
                int col_level = j_get_num64(j_array_index(p, 1));
                scalarType mat_val;
                sca_init(&mat_val);
                mat_val = j_get_double(j_array_index(p,2));
                if (row_level || col_level) {
                        fprintf(stderr,"error in %s:\n\texpected real ", __func__ );
                        fprintf(stderr,"scalar, but levels greater than zero!\n");
                }
                mats_ptr_t s_ptr = get_scalarPTR_for_scalarVal(mat_val);
                mapID[index] = s_ptr->packedID;
                if (verbose) {
                  printf("%"PRId64" %p %" PRId64 " 0 0 %Lg\n", index, 
                         s_ptr, mapID[index], mat_val);
                }
                sca_clear(&mat_val);
              }
              break;
#endif
#ifdef USE_INTEGER
            case 3:                        /* expect a scalar link */
              {
                if (verbose) {printf("In integer scalar case \n");}
                int row_level = j_get_num64(j_array_index(p, 0));
                int col_level = j_get_num64(j_array_index(p, 1));
                scalarType mat_val;
                sca_init(&mat_val);
                mat_val = j_get_num64(j_array_index(p,2));
                if (row_level || col_level) {
                        fprintf(stderr,"error in %s:\n\texpected integer ", __func__ );
                        fprintf(stderr,"scalar, but levels greater than zero!\n");
                }
                mats_ptr_t s_ptr = get_scalarPTR_for_scalarVal(mat_val);
                mapID[index] = s_ptr->packedID;
                if (verbose) {
                  printf("%"PRId64" %p %" PRId64 " 0 0 %ld\n", index, 
                         s_ptr, mapID[index], mat_val);
                }
                sca_clear(&mat_val);
              }
              break;
#endif
#ifdef USE_BOOLEAN
            case 3:                        /* expect a scalar link */
              {
                if (verbose) {printf("In boolean scalar case \n");}
                int row_level = j_get_num64(j_array_index(p, 0));
                int col_level = j_get_num64(j_array_index(p, 1));
                scalarType mat_val;
                sca_init(&mat_val);
                mat_val = j_get_num64(j_array_index(p,2));
                if (row_level || col_level) {
                        fprintf(stderr,"error in %s:\n\texpected boolean ", __func__ );
                        fprintf(stderr,"scalar, but levels greater than zero!\n");
                }
                mats_ptr_t s_ptr = get_scalarPTR_for_scalarVal(mat_val);
                mapID[index] = s_ptr->packedID;
                if (verbose) {
                  printf("%" PRId64 " %p %" PRId64 " 0 0 %ld\n", index, 
                         s_ptr, mapID[index], mat_val);
                }
                sca_clear(&mat_val);
              }
              break;
#endif
              
            case 6:                        /* expect a 4-tuple */
              {
                if (verbose) {printf("In non scalar case \n");}
                int row_level = j_get_num64(j_array_index(p, 0));
                int col_level = j_get_num64(j_array_index(p, 1));
                if (verbose) {printf("  levels are %d %d \n",row_level,col_level);}
                int64_t content_pID[4];
                int64_t temp;
                for (int k=0; k< 4; ++k) {
                   temp = j_get_num64(j_array_index(p, k+2));
                   if (verbose) {printf("temp for k = %d is %" PRId64 " \n",k,temp);}
                   // The matrix ID -1 is used whenever we have no matrix, 
                   // e.g. in row and col place keepers
                   if (temp == -1) { content_pID[k] = MATRIX_ID_INVALID;}
                   else {
                     content_pID[k] = mapID[temp];
                   }
                }
                mapID[index] = get_pID_from_array_of_four_sub_pIDs(
                    content_pID, row_level,col_level);
                if (verbose){
                  printf("ind %" PRId64 " mID %" PRId64 " (%d,%d)",
                        index, mapID[index], row_level, col_level);
                  printf(" [%" PRId64 " %" PRId64 " %" PRId64 " %" PRId64 "]\n",
                        content_pID[0], content_pID[1], content_pID[2],
                        content_pID[3]);
                }
              }
              break; 
              
            default:
              fprintf(stderr,"Error in %s(%s): unexpected number of entries per line\n",
                __func__, path);
#ifdef USE_COMPLEX
              fprintf(stderr,"In complex mode expected num entries to be 4 or 6, but had %d entries\n",len1);
#else
              fprintf(stderr,"In real or integer mode expected num entries to be 3 or 6, but had %d entries\n",len1);
#endif
              exit(1);
              break;
            }
        }
    }

  int64_t packedID = mapID[j_lookup_num64(j, "matid")];
  free(mapID);

  if (packedID == MATRIX_ID_INVALID)
  {
    fprintf(stderr,"in %s, input packedID is invalid\n",__func__);
    exit(0);
  }

  if (IS_SCALAR(packedID))
  {
      printf("\nRead scalar matrix from json file %s\n", path);
      mats_ptr_t scaptr = (mats_ptr_t)get_recordPTR_from_pID(packedID,"",__func__,0);
      printf ("  stored matrix now has address %p", scaptr);
  }
  else
  {
     matns_ptr_t matptr = (matns_ptr_t)get_recordPTR_from_pID(packedID,"",__func__,0);
     printf("\nRead %d,%d-level matrix from json file %s\n",
         matptr->row_level, matptr->col_level, path);
     printf ("  stored matrix now has address %p", matptr);
  }
  printf ("  and packedID %" PRId64 "\n", packedID);

  if (j_key_exists(j, "info")) {
    if (verbose) printf("Found info key in %s\n", __func__);
    json_t *t_info = j_key_lookup(j, "info");
    int64_t len_info = j_key_count(t_info);
    
    if (verbose) {printf("In %s, before the info line reader\n", __func__);}
    
    for (i=0; i<len_info-1; i++)
      {
        json_t *p = j_key_index(t_info, i);
        const char *info_name = j_get_name(p);
        const char *info_data = j_get_string(p);
        if (verbose) {printf("For i=%"PRId64", and key=%s\n", i, info_name);}
	info_type_t info_type = get_info_type_from_string_name((char *) info_name);
	int fail =  info_set(info_type, packedID, info_data);
        if (fail) {fprintf(stderr,"panic in %s\n",__func__); exit(1);}
        if (verbose) {
          printf("Added info field %s to InfoStore\n", info_name);
          printf("  with value %s\n", info_data);
        }
      }
  }  // end of if j_key_exists for "info"

  fclose(f);
  j_set_null(j); free(j);
  return packedID;
}

/************************************************************************
 * We will be creating a structure which looks like a quadtree and 
 * contains the information needed to create the LARC matrix quadtree.
 * We will use the indices 0 to num_nonzeros-1 for matrices containing
 * a single nonzero entry.  We will call these "leaf nodes" since they
 * will be leafs in our quadtree-like structure.
 * We will use indices from num_nonzeros to num_nonzeros + internal_nodes-1
 * for nodes that will represent matrices with more than one nonzero entry.  
 * We call the nodes with more than one nonzero entry "internal nodes"
 * and the internal nodes each have four children, their quads,  which will
 * either contain a -1 if they correspond to a zero quadrant submatrix
 * or contain a leaf node index or contain another internal node index.
 * Note that this quadtree-like structure is not a full quadtree: many
 * branches terminate early with a leaf labeled -1 to indicate the entire
 * branch below this is zero; also, many branches terminate early with a
 * leaf labeled with a small index (these correspond to a matrix with a
 * single nonzero entry).
 ************************************************************************/
/*!
 * \ingroup larc
 * \brief Constructs a quadtree-like structure from input coordinate/value data, to be used to generate the compressed LARCMatrix for the matrix described by this data
 * \param num_nonzeros The number of (i,j,val_{ij}) pairs in the input data structure
 * \param max_level The log of the size of the LARCMatrix to be generated
 * \param data A structure containing the (i,j,val_{ij}) pairs
 * \param qls A structure containing the information which will be used to construct the compressed LARCMatrix.
 */
static void create_QLS_from_position_value_pairs(unsigned num_nonzeros,
mat_level_t max_level, struct dataByCoordinate *data, struct dataForQLS *qls)
{
  int verbose = 1;

  /************************************************************************
   * For the worst case, the number of internal nodes is smaller than
   *     (num_nonzeros) * (the base_2 log of the dimension of the matrix).
   * We will be lazy and just create arrays big enough to hold the worst case
   * stuff. And we add an extra num_nonzeros to our array so that we can be
   * even more lazy and use distinct numbers for all the leaf nodes and
   * internal nodes.
   * Each internal node has 4 children, quad[0], quad[1], quad[2], and quad[3]
   ************************************************************************/
  const int64_t max_internal_nodes = num_nonzeros * (max_level+1);

  // first index that is beyond the leaf_node indices will be the internal
  // node with maximal level
  const int64_t top_internal_node_index = num_nonzeros;

  if (verbose==DEBUG) {
    fprintf(stderr,"we set max_internal_nodes to %" PRId64 "\n",max_internal_nodes);
    fprintf(stderr,"the first internalnode index is set to %" PRId64 "\n",
        top_internal_node_index);
  }

  /************************************************************************
   * We initialize the value of each quad to be -1 
   * indicating they correspond to zero matrices.
   ************************************************************************/
  for (int64_t node_index = top_internal_node_index ;
       node_index < max_internal_nodes ; node_index++)
  {
      for (int quad_index = 0 ; quad_index <4 ;quad_index++)
        qls->quad[node_index][quad_index] = -1;
  }
  if (verbose==DEBUG) fprintf(stderr,"all quads are initialized to -1\n");


  /************************************************************************
   * We create the first internal node in our quadtree like structure and
   * which will be the top node (the root). 
   * It corresponds to a 2^max_level by 2^max_level matrix.
   * Each child node in the four quads of a node will have level -1 of 
   * it's parent, and correspond to the four quadrant submatrics of the 
   * parent matrix.
   * This top node will have the first MatrixMarket entry with leaf_index =0 
   * as one of its quads initially. 
   * Whenever a new leaf node wants to join our quadtree like structure in
   * the same position as an old leaf, we will create a new internal node
   * and try to insert the colliding leaves in as children of the internal
   * node.  As we go down the tree, the level of the matrices corresponding
   * to nodes decreases.
   ************************************************************************/

  // set the root node in the QLS to max_level
  qls->internal_node_level[top_internal_node_index] = max_level;
  // the root  is our first internal node
  qls->num_internal_nodes = 1;

  if (verbose==DEBUG)
  {
    fprintf(stderr,
        "internal_node_level[top_internal_node_index] set to %d\n",max_level);
  }

  /************************************************************************
  * Loop over the MM nonzero entries and insert each into our own quadtree like
  * structure, creating internal nodes as needed, until each nonzero entry
  * is on a "leaf node" of some level. The leaf node will correspond to a 
  * submatrix of the given leaf_level, whose scalars are all zero except for 
  * a single nonzero entry given by val[leaf_node_index] whose
  * position is calculated using row[leaf_node_index] and col[leaf_node_index]
  * and the leaf_level of the submatrix.
  ************************************************************************/


  /************************************************************************
  * We loop over all the other nonzero entries in the MatrixMarket making
  * a leaf node for each nonzero entry and pushing the leaf node into our
  * structure by starting at the first internal node and branching downward
  * until we either can substitute into a position with a -1 (corresponding
  * to an all zero submatrix) or we have a collision with another leaf node,
  * in which case we create a new internal node at that location and push the
  * existing object under the new internal node, then try to add the new node
  * under the new internal node (which is easy if in a quad labeled -1,
  * otherwise we must continue to recurse.
  *
  * The first leaf node is inserted into the top internal node and given the
  * same level as that internal node (it will be lowered later when other
  * nodes are added). We loop over the other leaf nodes, modifying the
  * internal node tree structure as needed. Algorithmically, we have three
  * possible cases when inserting a new leaf L: 
  *
  * 1. current location holds -1,
  *    ACTION: replace -1 with our leaf node L
  *    break out of inner loop to deal with next leaf node
  *
  * 2. current location holds an internal node
  *                 (number >= top_internal_node_index)
  *    ACTION descend to this internal node and re-evaluate situation
  *    this is handled by a nested while() loop
  *
  * 3. current location holds a leaf node K
  *                 (number < top_internal_node_index)
  *    ACTION:  a. create new internal node B and replace K with B
  *                            this new node will have level one less than
  *                            the internal node contaning it
  *             b. K becomes a quad in B (since B is empty this is simple)
  *                     the leaf level of K is reduced by 1
  *             c. reduce the leaf level of L by 1
  *             d. decend to B and re-evaluate situation
  *                   this is handled by a nested while() loop
  ************************************************************************/

  // initial values for the loop
  mat_level_t current_level = max_level;

  // initialize the quadtree-like structure with the first element of the
  // matrix market data (index 0). This node's leaf_level is set to one less
  // than the max_level.
  int quad_index = get_quad_in_submatrix(data->row[0],
        data->col[0],current_level);
  qls->quad[top_internal_node_index][quad_index] = 0;  
  qls->leaf_level[0] = max_level-1;

  if (verbose==DEBUG)
  {
    fprintf(stderr,"For leaf_node_index 0, we set\n");
    fprintf(stderr,"quad_index = %d, ",quad_index);
    fprintf(stderr,"top_internal_node_index = %" PRId64 "\n",
        top_internal_node_index);
    fprintf(stderr,"quad[%" PRId64 "][%d] set to 0 and has level %d\n",
        top_internal_node_index, quad_index,
        qls->internal_node_level[top_internal_node_index]);
    fprintf(stderr,"setting leaf_level[0] to %d\n",qls->leaf_level[0]);
  } 

  // loop to place the leaf_nodes into our quadtree like structure
  // and create any internal nodes necessary to tack the leafs together
  for (int64_t leaf_node_index=1; leaf_node_index<top_internal_node_index;
        ++leaf_node_index) 
  {
    // the algorithm first attempts to insert a leaf node into the
    // internal node at top node of tree
    int64_t node_index = top_internal_node_index;
    current_level = max_level; // the level of the internal node at top of tree
    if (verbose==DEBUG)
    {
      fprintf(stderr,"\ntop of loop: leaf_node_index is now %" PRId64 "\n",
                leaf_node_index);
      fprintf(stderr,"node_index set to %" PRId64 "\n",node_index);
    }

    // while loop exits when leaf node has been successfully inserted
    while(1) {
      if (verbose==DEBUG)
      {
        fprintf(stderr,"top of while(1) loop:\n");
        fprintf(stderr,"calling get_quad_in_submatrix\n");
      }
      // find the quad index for the current leaf. This depends on the
      // current matrix level, which will decrease as loop repeats
      quad_index = get_quad_in_submatrix(data->row[leaf_node_index],
                data->col[leaf_node_index], current_level);
      if (verbose==DEBUG)
        fprintf(stderr,"returned quad_index value is %d\n",quad_index);

      // case 1: The quad we would like to put the leaf node in corresponds
      // to a zero matrix (which is indicated when the value in the quad is
      // -1). We can replace the zero matrix with the leaf node.
      if (qls->quad[node_index][quad_index] == -1) {
        // zero matrix is currently in desired location, so we can insert
        // leaf node and continue the leaf_node_index loop
        qls->quad[node_index][quad_index] = leaf_node_index;
        // this leaf's level is set to one level below internal node level
        qls->leaf_level[leaf_node_index] = current_level-1;   
        if (verbose==DEBUG)
        {
          fprintf(stderr,"case 1 executed:\n");
          fprintf(stderr,"quad[%" PRId64 "][%d] set to %" PRId64 " with leaf level %d\n",
                node_index,quad_index,leaf_node_index,current_level-1);
          fprintf(stderr,"this node is level %d\n",
                qls->internal_node_level[node_index]);
        }
        break; // out of while(1), go to next leaf node
      }

      // case 2: The quad we would like to put the leaf node in corresponds
      // to an internal node (indicated when the value in the quad is a node
      // number that is greater than or equal to top_internal_node_index).
      // Since an internal node was in desired position to insert the leaf
      // node, we need to descend further down into the tree from the internal
      // node and stay in the while loop to attempt insertion again.
      if (qls->quad[node_index][quad_index] >= top_internal_node_index) {
        node_index = qls->quad[node_index][quad_index];
        --current_level;
        if (verbose==DEBUG)
        {
          fprintf(stderr,"case 2 executed: node index changed to ");
          fprintf(stderr,"%" PRId64 ",\n current_level changed to %d\n",
                node_index, current_level);
          fprintf(stderr,"this node is level %d\n",
                qls->internal_node_level[node_index]);
        }
        continue; // while(1)
      }

      // case 3: The quad value corresponds to a previously stored leaf node
      // (indicated by the value in the quad being between 0 and
      // leaf_node_index). In this case we have a collision with a previously
      // inserted leaf node, so we will create a new internal node and push
      // the previous leaf into one of the quads of this new internal node,
      // then continue in the while loop trying to insert our leaf node
      // starting from the new internal node.

      if (verbose==DEBUG) fprintf(stderr,"case 3 executing:\n");
      int64_t old_leaf_node_index = qls->quad[node_index][quad_index];

      // case 3': Make sure that the new and old leaf nodes do not have the
      // same (i,j) values. This will not happen in correct MM data, but we
      // need to protect against it.
      if (data->row[old_leaf_node_index]==data->row[leaf_node_index] &&
          data->col[old_leaf_node_index]==data->col[leaf_node_index])
      {
        // if data values are different, we can't go on; if data is repeated,
        // though technically invalid, we could continue, but right now we
        // don't have the code written to fix up our data structure after this
        // happens.
        // So, either way we exit the program
        //
        // The error message adds 1 to the stored (i,j) values to compensate
        // for our previously subtracting 1 to become zero-based
        fprintf(stderr,"ERROR: input matrix market data has two (or more) ");
        fprintf(stderr,"entries at position (%d, %d)\n",
          1+data->row[leaf_node_index],1+data->col[leaf_node_index]);
        // add additional info about whether data is same or different
        if (sca_eq(data->val[old_leaf_node_index],data->val[leaf_node_index]))
          fprintf(stderr,"with the same value!\n");
        else
          fprintf(stderr,"with different values!\n");
        fprintf(stderr,"exiting due to bad input data...\n");
        exit(0);
      }

      // a) create new internal node, level one less than current internal node
      current_level--;
      int64_t new_internal_node_index =
                top_internal_node_index + qls->num_internal_nodes;
      qls->internal_node_level[new_internal_node_index] = current_level;
      qls->num_internal_nodes++;
      // b) replace old leaf with new internal node
      qls->quad[node_index][quad_index] = new_internal_node_index;
      if (verbose==DEBUG)
      {
        fprintf(stderr,"create new internal node with index %" PRId64 " and level %d\n",
                      new_internal_node_index, current_level);
        fprintf(stderr,"quad[%" PRId64 "][%d] changed to %" PRId64 "\n",
                      node_index, quad_index, new_internal_node_index);
        fprintf(stderr,"node %" PRId64 " is level %d\n",
                node_index, qls->internal_node_level[node_index]);
      }
      // c) place old leaf in a quad of the new internal node, and also
      //    change the node level of this leaf
      if (verbose==DEBUG)
      {
        fprintf(stderr,"calling get_quad_in_submatrix, with ");
        fprintf(stderr,"row value %d, column value %d, and ",
                      data->row[old_leaf_node_index], data->col[old_leaf_node_index]);
        fprintf(stderr,"level %d\n",current_level);
      }
      quad_index = get_quad_in_submatrix(data->row[old_leaf_node_index],
                      data->col[old_leaf_node_index],current_level);
      qls->quad[new_internal_node_index][quad_index] = old_leaf_node_index;
      qls->leaf_level[old_leaf_node_index] = current_level-1;
      if (verbose==DEBUG)
      {
        fprintf(stderr,"resulting quad_index is %d\n",quad_index);
        fprintf(stderr,"quad[%" PRId64 "][%d] set to %" PRId64 "\n",
                new_internal_node_index,quad_index,old_leaf_node_index);
        fprintf(stderr,"node %" PRId64 " is level %d\n", new_internal_node_index,
                qls->internal_node_level[new_internal_node_index]);
        fprintf(stderr,"leaf level of leaf %" PRId64 " is now %d\n",
                old_leaf_node_index, qls->leaf_level[old_leaf_node_index]);
      }
      // d) set the node_index to the new internal node then continue in the
      // while loop
      node_index = new_internal_node_index;
      continue;
    } // while (1)
  } // for (leaf_node_index)

  if (verbose==DEBUG)
        fprintf(stderr,"completed loops over leaf node index, while(1)\n");
} // end create_QLS_from_position_value_pairs

/*!
 * \ingroup larc
 * \brief Uses a quadtree-like structure and the value field from input coordinate/value data to generate the compressed LARCMatrix for the matrix described by this data
 * \param num_nonzeros The number of (i,j,val_{ij}) pairs in the input data structure
 * \param data A structure containing the (i,j,val_{ij}) pairs
 * \param qls A structure containing the information which will be used to construct the compressed LARCMatrix.
 * \return The packedID of the LARCMatrix generated
 */
static int64_t build_larcMatrix_using_QLS(int64_t num_nonzeros,
        struct dataByCoordinate *data, struct dataForQLS *qls)
{
  const int verbose = 1;

  // Now we will create the LARC matrix from our quadtree like reprentation
  // of the MatrixMarket matrix and the input data values
  
  // First we go through all the leaf nodes of our structure and load a
  // matrix into LARC that has the right level and has a single nonzero entry
  // in the appropriate position.

  // Create an array of matrix pointers for LARC matrices
  int64_t *node_ID_ptr;
  node_ID_ptr = (int64_t *)malloc(
        (num_nonzeros + qls->num_internal_nodes)*sizeof(int64_t));
  if (node_ID_ptr == NULL) {
    fprintf(stderr,"ERROR in %s, line %d: Out of memory.\n",__func__,__LINE__);
    exit(1);
  }

  // Matrix Pointers are created first for all the maximally sized
  // quadrant submatrices in the LARC matrix with only a single nonzero
  // value.  These are precisely those corresponding to the leaf nodes
  // of our quadtreelike structure, which corresponded to the nonzero
  // values we read in from the MatrixMarket file.

  // Loop through the leaf_index for all the nonzero scalars
  // passing their value, level, and row, col coordinates
  // and construct the appropriate matrix in LARC and return the ptr.
  if (verbose==DEBUG) fprintf(stderr,"start of loop over leaves\n");

  for (int leaf_index=0;leaf_index<num_nonzeros;++leaf_index) {
    if (verbose==DEBUG)
    {
      fprintf(stderr,"leaf_index = %d\n",leaf_index);
      fprintf(stderr," leaf_level[%d] = %d\n",
                    leaf_index,qls->leaf_level[leaf_index]);
      fprintf(stderr," row[%d] = %d\n",leaf_index,data->row[leaf_index]);
      fprintf(stderr," col[%d] = %d\n",leaf_index,data->col[leaf_index]);
      fprintf(stderr," val[%d] = %s\n",leaf_index,
	sca_get_readable_approx_str(data->val[leaf_index]));
    }

    int64_t valID=get_valID_from_valString(sca_get_readable_approx_str(data->val[leaf_index]));
    node_ID_ptr[leaf_index] = get_matrix_with_single_nonzero_at_coords(
           qls->leaf_level[leaf_index], data->row[leaf_index],
           data->col[leaf_index], valID);

    // The scalar value is now stored in the matrixStore, so we no longer
    // need it as a scalarType. Clearing it will free up memory when
    // scalarType is multiprecision (and otherwise is a no-op)
    sca_clear(data->val+leaf_index);

  } // loop over leaf_index
  if (verbose==DEBUG) fprintf(stderr,"end of loop over leaves\n");

  // at this point, we can free the memory from several arrays, namely
  // data->row, data->col, data->val, and qls->leaf_level.
  free(data->val); data->val = (scalarType *)NULL;
  free(data->row); data->row = (unsigned *)NULL;
  free(data->col); data->col = (unsigned *)NULL;
  free(qls->leaf_level); qls->leaf_level = (mat_level_t *)NULL;

  if (verbose) 
  {
     fprintf(stderr,"matrix market reader:");
     fprintf(stderr,"finished making single-nonzero submatrices\n");
     fprintf(stderr,"now building LARC matrix \n");
  }

  // PackedIDs are created for all the internal nodes by going through
  // the indices in reverse order to their insertion. We will use the set of
  // quad values for each internal_node to create a list of packedIDs
  // to build the LARC matrix for that internal node.
  int64_t subMatList[4];

  // #nodes = #leaves + #internal nodes
  const int64_t num_nodes = num_nonzeros + qls->num_internal_nodes;

  if (verbose==DEBUG) fprintf(stderr,"start of loop over internal nodes\n");
  // since the node with maximum level is the first internal node, we loop
  // backward over the internal nodes to fill up the LARCMatrix structure;
  // the lower-level submatrices we create first will be used in building up
  // the higher-level submatrices.

  // The top internal node index == num_nonzeros;
  for (int64_t internal_node_index = num_nodes - 1;
       internal_node_index>=num_nonzeros; --internal_node_index)
  {
    // internal_node_index
    const mat_level_t level = qls->internal_node_level[internal_node_index];
    if (verbose==DEBUG)
    {
      fprintf(stderr,"\tinternal node index is %" PRId64 "\n",internal_node_index);
      fprintf(stderr,"\tinternal node level is %d\n",level);
    }

    // since we assuming a sparse matrix we retrieve the packedID
    // for the zero matrix of the correct level for quads
    int64_t zero_packedID = get_zero_pID(level-1,level-1);

    // Build a subMatList from packedIDs associated with each quad value,
    if (verbose==DEBUG)
        fprintf(stderr,"start of loop over internal node indices\n");

    for (int index=0;index<4;++index) {
      int64_t quad_value = qls->quad[internal_node_index][index];
      if (verbose==DEBUG)
      {
        fprintf(stderr,"index = %d\n",index);
        fprintf(stderr,"quad_value = %" PRId64 "\n",quad_value);
      }

      // get the packedID for the zero or nonzero submatrix
      if (quad_value == -1) subMatList[index] = zero_packedID;  
      else subMatList[index] = node_ID_ptr[quad_value];
      if (verbose==DEBUG) fprintf(stderr,"going to next index value\n");
    }

    if (verbose==DEBUG)
        fprintf(stderr,"end of loop over internal node indices\n");

    // build a matrix from the four subMatIDs
    node_ID_ptr[internal_node_index] = get_pID_from_array_of_four_sub_pIDs(
                subMatList, level, level);
  } // end loop over internal_node_index
  if (verbose==DEBUG)
  {
     fprintf(stderr,"end of loop over internal nodes\n");
     fprintf(stderr,"num_nonzeros = %" PRId64 "\n",num_nonzeros);
  }

  // The pointer to the LARC matrix is the pointer associated with the top node
  // of our node_ptr list, which is at num_nonzeros (just after the last leaf
  // node pointer)
//  matns_ptr_t ret_ptr = node_ptr[num_nonzeros];
  int64_t ret_pID = node_ID_ptr[num_nonzeros];
  if (verbose==DEBUG)
     fprintf(stderr,"about to free node_ptr array\n");

  free(node_ID_ptr);

  if (verbose==DEBUG) 
     fprintf(stderr,"freed node_ptr array\n");

  return ret_pID;
} // end build_larcMatrix_using_QLS


/********************************************************
 * routine: read_matrixMarketExchange_file()
 * Algorithm for converting from Matrix Market coordinate to LARC format 
 * created by Jenny Zito and coded by Jenny and Steve Cuccaro.
 * to read a Matrix Market Exchange "coordinate" formated
 * file which looks like the following.
 *********************************************************
 * %% MatrixMarket matrix coordinate integer general
 * 
 * % Comments: The matrix described in this file is 
 * %           coordinate format, with integer values,
 * %           and the matrix has no special symmetry type.
 * %           It is 8 by 8 and has 5 nonzero entries; this
 * %           information about the matrix is the first
 * %           uncommented line of the file, followed by 
 * %           a line for each nonzero entry, specifying:
 * %           row column value
 * 
 * 8 8 5
 * 2 2 10
 * 6 2 250
 * 3 5 -1
 * 1 6 33
 * 7 4 -15
 *****************************************************/
int64_t read_matrixMarketExchange_file(char * file_path)
{

#ifdef DEBUG_IO_C
  printf("--> %s : %s: %d\n", __FILE__, __func__, __LINE__);
#endif // #ifdef DEBUG_IO_C 
  const int verbose = 1;
  int error_code = 0;

  if (verbose) printf("In %s\n", __func__);

  
  /************************************************************************
  * Open the Matrix Market matrix file.
  ************************************************************************/
  FILE *fp;
  if((fp = fopen(file_path, "r")) == NULL) {
    fprintf(stderr,"No such file: %s\n",file_path);
    exit(1);   
  }
  if (verbose) printf("Opened the data file %s for reading:\n",file_path);
  
  // read in header line
#define MM_MAX_LINE_LENGTH 1025
#define MM_MAX_TOKEN_LENGTH 64
  char headerMM[MM_MAX_LINE_LENGTH];
  if (fgets(headerMM,MM_MAX_LINE_LENGTH, fp) == NULL)
  {
     fprintf(stderr,"matrix market reader: input file has no header!\n");
     return MATRIX_ID_INVALID;
  }

  // parse header line, check component values
  char bannerMM[MM_MAX_TOKEN_LENGTH];
  char mtxMM[MM_MAX_TOKEN_LENGTH];
  char crdMM[MM_MAX_TOKEN_LENGTH];
  char data_typeMM[MM_MAX_TOKEN_LENGTH];
  char storage_schemeMM[MM_MAX_TOKEN_LENGTH];

  int ival;
  if ((ival=sscanf(headerMM, "%s %s %s %s %s", bannerMM, mtxMM,
        crdMM, data_typeMM, storage_schemeMM)) != 5)
  {
     fprintf(stderr,"matrix market reader: input file has incorrect header!\n");
     fprintf(stderr,"only read in %d terms of 5.\n",ival);
     return MATRIX_ID_INVALID;
  }

  // check banner (first field) for validity
  if (strcmp(bannerMM,"%%MatrixMarket"))
  {
     fprintf(stderr,"matrix market reader: 'banner' ('%s') is incorrect!\n",
        bannerMM);
     return MATRIX_ID_INVALID;
  }

  // check second field ('matrix')
  char *p;
  for (p=mtxMM; *p!='\0'; *p=tolower(*p),p++); // convert to lower case
  if (strcmp(mtxMM,"matrix"))
  {
     fprintf(stderr,"matrix market reader: 'mtx' (%s)is incorrect!\n", mtxMM);
     return MATRIX_ID_INVALID;
  }

  // check third field ('coordinate' or 'array')
  for (p=crdMM; *p!='\0'; *p=tolower(*p),p++); // convert to lower case
  if (0==strcmp(crdMM,"coordinate")); // sparse format
  else // LARC only processes sparse format MM files
  {
     fprintf(stderr,"matrix market reader: 'crd' is '%s';\n",crdMM);
     fprintf(stderr,"LARC does not currently support this format.\n");
     return MATRIX_ID_INVALID;
  }

  // check fourth field (datatype)
  for (p=data_typeMM; *p!='\0'; *p=tolower(*p),p++); // convert to lower case
  // valid types are integer, real, double, complex, pattern
  if (strcmp(data_typeMM,"integer") && strcmp(data_typeMM,"real") &&
        strcmp(data_typeMM,"double") && strcmp(data_typeMM,"complex")  &&
        strcmp(data_typeMM,"pattern") )
  {
     fprintf(stderr,"matrix market reader: 'data_type' is '%s';\n",data_typeMM);
     fprintf(stderr,"LARC does not currently support this type.\n");
     return MATRIX_ID_INVALID;
  }

  // check fifth field
  for (p=storage_schemeMM; *p!='\0'; *p=tolower(*p),p++); //convert to lowercase
  // valid schemes are general, symmetric, skew-symmetric, hermitian
  if (   strcmp(storage_schemeMM,"general")
      && strcmp(storage_schemeMM,"symmetric")
      && strcmp(storage_schemeMM,"skew-symmetric")
      && strcmp(storage_schemeMM,"hermitian") )
  {
     fprintf(stderr,"matrix market reader: 'storage_scheme' is '%s';\n",
        storage_schemeMM);
     fprintf(stderr,"LARC does not currently support this scheme.\n");
     return MATRIX_ID_INVALID;
  }

  // check what type LARC is compiled for, and give warning if there might
  // be a problem 
  int typeflag = 0;
  if (0==strcmp(data_typeMM,"complex"))
  {
     typeflag = strcmp(scalarTypeStr,"Integer");
     typeflag |= strcmp(scalarTypeStr,"MPInteger");
     typeflag |= strcmp(scalarTypeStr,"Real");
     typeflag |= strcmp(scalarTypeStr,"MPReal");
     typeflag |= strcmp(scalarTypeStr,"MPRational");
  }
  else if ( (0==strcmp(data_typeMM,"real")) ||
        (0==strcmp(data_typeMM,"double")) )
  {
     typeflag = strcmp(scalarTypeStr,"Integer");
     typeflag |= strcmp(scalarTypeStr,"MPInteger");
  }
  if (typeflag)
  {
     fprintf(stderr,"matrix market reader: data type of file is\n");
     fprintf(stderr,"%s, but LARC is compiled with scalarType %s\n",
               data_typeMM, scalarTypeStr);
     fprintf(stderr,"This could be a problem...\n");
  }

  /************************************************************************
  * Reading the Matrix Market coordinate file: 
  *   COMMENTS AND BLANK LINES:
  *       * ignore comment lines starting with "% "
  *       * ignore blank lines
  ************************************************************************/
  char line[MM_MAX_LINE_LENGTH];
  do
  {
      if (fgets(line,MM_MAX_LINE_LENGTH,fp)==NULL)
      {
         fprintf(stderr,"matrix market reader: no data in file %s\n",file_path);
         return MATRIX_ID_INVALID;
      }
  } while (line[0] == '%');
  
  /************************************************************************
  * Reading the Matrix Market coordinate file:
  *   MATRIX STATISTICS:
  *       * the first line of the MMfile that does not start with a "%" contains
  *         information about the size of the MMmatrix and number of nonzero entries.
  *         MM_num_rows  MM_num_cols num_nonzeros
  ************************************************************************/
  unsigned MM_num_rows;
  unsigned MM_num_cols;
  unsigned MM_num_nonzeros;
  
  // read in the first line from the MatrixMarket file not starting with "%"
  int ret = sscanf(line, "%u %u %u", &MM_num_rows, &MM_num_cols,
         &MM_num_nonzeros);
  if (ret == 3) {
    if (verbose)
       printf("\nMM_num_rows, MM_num_cols, MM_num_nonzeros: %u %u %u\n",
                        MM_num_rows, MM_num_cols, MM_num_nonzeros);
    fflush(stdout);
  } 
  else if(errno != 0) {
    perror("fscanf:");
    error_code = 1;
    fprintf(stderr,"ERROR scanning levels in first line of file %s\n", file_path);
  }
  else if(ret == EOF) {
    error_code = 2;
    fprintf(stderr,"ERROR: file %s ended early in %s\n", file_path,  __func__);
  }
  else {
    fprintf(stderr,"return value is %d\n",ret);
    fprintf(stderr,"ERROR while reading levels in %s (routine %s).\n", file_path, __func__);
    error_code = 2;
  }
  if (error_code) {
#ifdef DEBUG_IO_C
    printf("<-- %s : fail: %d\n", __FILE__, __LINE__);
#endif // #ifdef DEBUG_IO_C 
    return MATRIX_ID_INVALID; 
  } 

  // if the file is open but nothing is left to read print EOF
  if(feof(fp)){
    fputs("EOF",stderr);
    fprintf(stderr,"ERROR: File stopped earlier than expected.\n");
    fclose(fp);
#ifdef DEBUG_IO_C
    printf("<-- %s : fail: %d\n", __FILE__, __LINE__);
#endif // #ifdef DEBUG_IO_C 
    return MATRIX_ID_INVALID; 
  }

  // if matrix is not 'general', it must be square
  // exit if input data dimensions are inconsistent with storage scheme
  if ( strcmp(storage_schemeMM,"general") && (MM_num_rows != MM_num_cols) )
  {
    fprintf(stderr,"ERROR in matrix market input data: non-square matrix\n");
    fprintf(stderr,"has storage scheme '%s': exiting...\n",storage_schemeMM);
    exit(1);
  }

  /************************************************************************
  * To store the MM_num_rows by MM_num_columns Matrix Market matrix, 
  * LARC will use a 2^row_level by 2^col_level matrix
  ************************************************************************/
  mat_level_t row_level;
  mat_level_t col_level;

  // calculating the row and column levels
  row_level = 0;
  int64_t dim = 1;
  while (dim < MM_num_rows) {
    dim <<= 1;
    row_level++;
  }

  col_level = 0;
  dim = 1;
  while (dim < MM_num_cols) {
    dim <<= 1;
    col_level++;
  }

  if (verbose) printf("  levels are %d %d\n",(int)row_level, (int)col_level);


  // EVENTUALLY HANDLE NONSQUARE MATRICES
  // Temporary hack: make the LARC matrix square
  mat_level_t max_level = MAX(row_level,col_level);
  if (row_level != col_level) {
    row_level = col_level = max_level;
  }

  /************************************************************************
  * Reading the Matrix Market coordinate file:
  *    NONZERO MATRIX ENTRIES:
  *       * Each nonzero matrix entry will have a line in the MMfile with
  *         row and col coordinates followed by the value of scalar in that
  *         position:
  *         I J A(I,J)
  *         If the data_typeMM is "complex" then A(I,J) has two entries REAL
  *         and IMAG.
  *    SYMMETRY: If the matrix is
  *        * 'general': there are MM_num_nonzeros nonzero entries.
  *        * 'skew-symmetric': there are 2*MM_num_nonzeros nonzero
  *           entries, since diagonal elements must be zero.
  *        * 'symmetric' or 'hermitian': there are at most 2*MM_num_nonzeros
  *          entries, because diagonal elements are not duplicated
  ************************************************************************/
  unsigned max_num_nonzeros = MM_num_nonzeros;
  if ( strcmp(storage_schemeMM,"general") ) max_num_nonzeros *= 2;

  struct dataByCoordinate data;
  data.row = (unsigned *)malloc(max_num_nonzeros*sizeof(unsigned));
  if (data.row == NULL) {
    fprintf(stderr,"ERROR in %s, line %d: Out of memory.\n", __func__, __LINE__);
    exit(1);
  }
  data.col = (unsigned *)malloc(max_num_nonzeros*sizeof(unsigned));
  if (data.col == NULL) {
    fprintf(stderr,"ERROR in %s, line %d: Out of memory.\n", __func__, __LINE__);
    exit(1);
  }
  data.val = (scalarType *)malloc(max_num_nonzeros*sizeof(scalarType));
  if (data.val == NULL) {
    fprintf(stderr,"ERROR in %s, line %d: Out of memory.\n", __func__, __LINE__);
    exit(1);
  }

  if (verbose==DEBUG) {
     fprintf(stderr,"have allocated data arrays, using ");
     size_t mem_alloc = max_num_nonzeros*(2*sizeof(unsigned) + sizeof(scalarType));
     fprintf(stderr,"%zd bytes of memory\n",mem_alloc);
  }

  // Now we will handle the entries of the matrix
  // The following code is all that is needed for scheme = 'general';
  int i;
  if (!strcmp(data_typeMM,"complex")) {
    long double val_real, val_imag;
    for (i=0;i<MM_num_nonzeros;++i) {
      ret = fscanf(fp, "%u %u %Lg %Lg", data.row+i, data.col+i,
                                      &val_real, &val_imag);
      sca_init(data.val+i);
      sca_set_2ldoubles(data.val+i,val_real,val_imag);
    }
  }
  else if ( (!strcmp(data_typeMM,"real")) || (!strcmp(data_typeMM,"double")) )
  {
    long double val_real;
    for (i=0;i<MM_num_nonzeros;++i) {
      ret = fscanf(fp, "%u %u %Lg", data.row+i, data.col+i, &val_real);
      sca_init(data.val+i);
      sca_set_2ldoubles(data.val+i,val_real,(long double)0.0);
    }
  }
  else if (!strcmp(data_typeMM,"integer")) {
    int64_t val_integer;
    for (i=0;i<MM_num_nonzeros;++i) {
      ret = fscanf(fp, "%u %u %" SCNd64, data.row+i,data.col+i,&val_integer);
      sca_init(data.val+i);
      sca_set_2ldoubles(data.val+i, (long double)val_integer,(long double)0.0);
    }
  }
  else if (!strcmp(data_typeMM,"pattern")) {
    ret = fscanf(fp, "%u %u",data.row,data.col);
    sca_init(data.val);
    sca_set_2ldoubles(data.val, (long double)1.0,(long double)0.0);
    for (i=1;i<MM_num_nonzeros;++i) {
      ret = fscanf(fp, "%u %u",data.row+i,data.col+i);
      sca_init(data.val+i);
      sca_set(data.val+i,data.val[0]);
    }
  }
  else {
    fprintf(stderr,"FAIL in %s reading arrays: Unknown data_typeMM\n",__func__);
    exit(0);
  }
  fclose(fp);

  unsigned num_nonzeros = max_num_nonzeros;
  // deal with non-general cases: any of the other schemes taking advantange
  // of symmetry [symmetric, skew-symmetric, hermitian] would put the value
  // read into position (i,j) and the [same, negated, adjoint] value into
  // position (j,i), and would only be valid for square matrices
  if ( strcmp(storage_schemeMM,"general") )
  {
    for (i=MM_num_nonzeros;i<max_num_nonzeros;++i)
    {
      int i0 = i - MM_num_nonzeros;
      if (data.row[i0] == data.col[i0]) { --num_nonzeros; }
      else
      {
        data.row[i] = data.col[i0];
        data.col[i] = data.row[i0];
        sca_init(data.val+i);
        if (!strcmp(storage_schemeMM,"symmetric"))
        {
          sca_set(data.val+i,data.val[i0]);
        }
        else if (!strcmp(storage_schemeMM,"skew-symmetric"))
        {
          sca_mult(data.val+i,data.val[i0],scalarM1);
        }
        else if (!strcmp(storage_schemeMM,"hermitian"))
        {
          sca_conj(data.val+i,data.val[i0]);
        }
        // other options (skew-hermitian?) would go here
        else
        {
           fprintf(stderr,"ERROR in matrix market reader: non-implemented\n");
           fprintf(stderr,"storage_scheme found in processing\n");
           exit(2);
        }
      }
    }
  }

  // change all indices to be zero based instead
  for (i=0;i<num_nonzeros;++i) {
    --data.row[i];
    --data.col[i];
    if (verbose==DEBUG) 
    {
        fprintf(stderr,"\trow[%d] decremented, is now %d\n",i,data.row[i]);
        fprintf(stderr,"\tcol[%d] decremented, is now %d\n",i,data.col[i]);
    }
  }

  // data has been read in
  if (verbose) 
  {
     fprintf(stderr,"matrix market reader: data read in\n");
     fprintf(stderr,"now adding data to Quadtree-Like Structure\n");
  }

  // TRIVIAL CASE: If there is only one nonzero element in the MatrixMarket
  // matrix then we return the packedID of the LARC matrix with a single
  // entry given by this nonzero element.
  if (verbose==DEBUG) fprintf(stderr,"num_nonzeros is %d\n",num_nonzeros);
  if (num_nonzeros == 1) {
    int64_t sca_pID = get_scalarPTR_for_scalarVal(data.val[0])->packedID;
    if ((row_level==0) && (col_level==0)) // SCALAR matrix
       return sca_pID;
    int64_t mat_pID = get_matrix_with_single_nonzero_at_coords(
           max_level, data.row[0], data.col[0], sca_pID);
    return mat_pID;
  }

  // NON-TRIVIAL CASE: we first construct a quadtree-like structure from
  // the row and column coordinates, then use this to determine the largest
  // submatrices containing only one nonzero value, and build the LARCMatrix
  // from these submatrices. We have to allocate a fair amount of memory for
  // the structures generated, and some of this memory is freed within the
  // routine build_larcMatrix_using_QLS() once it is no longer needed.

  // allocate arrays needed for quadtree-like structure
  struct dataForQLS qls;
  qls.leaf_level = (mat_level_t *)malloc(num_nonzeros*sizeof(mat_level_t));
  if (qls.leaf_level == NULL) {
    fprintf(stderr,"ERROR in %s, line %d: Out of memory.\n", __func__, __LINE__);
    exit(1);
  }
  int64_t max_internal_nodes = num_nonzeros * (max_level+1);
  qls.internal_node_level = (mat_level_t *)
        malloc(max_internal_nodes*sizeof(mat_level_t));
  if (qls.internal_node_level == NULL) {
    fprintf(stderr,"ERROR in %s, line %d: Out of memory.\n", __func__, __LINE__);
    exit(1);
  }
  // allocate memory for the quad doubly-dimensioned array
  qls.quad = (int64_t **)malloc(max_internal_nodes*sizeof(int64_t *));
  if (qls.quad == NULL) {
    fprintf(stderr,"ERROR in %s, line %d: Out of memory.\n", __func__, __LINE__);
    exit(1);
  }
  for (i=0;i<max_internal_nodes;++i) {
      qls.quad[i] = (int64_t *)malloc(4*sizeof(int64_t));
      if (qls.quad[i] == NULL) {
        fprintf(stderr,"ERROR in %s, line %d: Out of memory.\n", __func__, __LINE__);
        exit(1);
      }
  }

  // create a quadtree-like structure which provides information on the
  // largest submatrix for which each non-zero entry is the only non-zero
  // This also determines the number of internal nodes needed for the qls
  create_QLS_from_position_value_pairs(num_nonzeros, max_level, &data, &qls);

  if (verbose==DEBUG)
  {
     size_t mem_alloc = 0;
     // memory needed for qls structure
     mem_alloc += num_nonzeros*sizeof(mat_level_t);
     mem_alloc += max_internal_nodes*sizeof(mat_level_t);
     mem_alloc += max_internal_nodes*sizeof(int64_t *);
     mem_alloc += max_internal_nodes*4*sizeof(int64_t);
     // memory allocated internally in build_larcMatrix_using_QLS
     mem_alloc += (num_nonzeros + qls.num_internal_nodes)*sizeof(matns_ptr_t);
     fprintf(stderr,"total memory needed for QLS and scratch space, not\n");
     fprintf(stderr,"counting memory allocated for multiprecision\n");
     fprintf(stderr,"variables: %zd bytes\n",mem_alloc);
  }
  if (verbose) 
  {
     fprintf(stderr,"matrix market reader: QLS created\n");
     fprintf(stderr,"now making LARC submatrices with single nonzeros\n");
  }

  // We build the submatrices of the LARC matrix which 
  // 1) each have a single nonzero and
  // 2) are as large as possible for this matrix
  // By design, the remaining submatrices are all-zero

  int64_t m_pID = build_larcMatrix_using_QLS(num_nonzeros, &data, &qls);

  // Release any memory allocated to the structures initialized
  // by this routine that have not already been freed by the previous
  // subroutine.
  free(qls.internal_node_level);
  for (int64_t i=0; i< max_internal_nodes; ++i) free(qls.quad[i]);
  free(qls.quad);

  if (verbose) fprintf(stderr,"matrix market reader: finished\n");
  return m_pID;

}  // end read_matrixMarketExchange_file


/*!
 * \ingroup larc
 * \brief This is an internal subroutine to write a LARCmatrix entry to a file.
 *   If needed, sub-component entries will be written out first.
 *
 * \param f The open file pointer for the output file
 * \param s The matrixID of the entry to be written out
 * \param usasPTR The unique submatrix array structure for the matrix
 * \param matID_map The array of renumbered matrix IDs
 * \param map_counter A pointer to the next number to use for renumbering
 * \param do_renumber If non-zero, renumber the matrix IDs on output
 */
static void recursive_larcMatrixFile_element_output(
                FILE *f, int64_t s, usas_t *usasPTR,
                int64_t *matID_map, int64_t *map_counter, int do_renumber)
{
    // if the array value is zero then we don't do anything for this s
    if (usasPTR->array[s] == 0) return;

    if (do_renumber) {
        // map the value s to the smallest non-neg integer that we have not used yet
        // then increment the map_counter
        matID_map[s] = (*map_counter)++;
    } else {
        // Keep the matrix ID unchanged.
        matID_map[s] = s;
    }

    if (usasPTR->array[s] == 1)
    {
        mats_ptr_t s_ptr =(mats_ptr_t)get_recordPTR_from_pID(
            PID_FROM_SCALAR_MID(s),"",__func__,0);
        // in COMPLEX case, write routine to output "a+I*b" as "a, b"
        char *val_string =sca_get_readable_approx_str(s_ptr->scalar_value);
        fprintf(f, "    \"%" PRId64 "\":[0, 0, \"%s\"],\n",
            matID_map[s], val_string);
        free(val_string);

        // Set to zero so as not to print this element out twice.
        usasPTR->array[s] = 0;
    } // end scalar branch
    else   // NONSCALAR
    {
        matns_ptr_t s_ptr = (matns_ptr_t)get_recordPTR_from_pID(
            PID_FROM_NONSCALAR_MID(s),"",__func__,0);

        // recursively print out sub-components first,
        // if they haven't already been printed out yet.
        for (int i = 0; i < 4; ++i) {
            if (s_ptr->subMatList[i]!=MATRIX_ID_INVALID) {
                int64_t submatrixID = MID_FROM_PID(s_ptr->subMatList[i]);
                recursive_larcMatrixFile_element_output(f, submatrixID, usasPTR, matID_map, map_counter, do_renumber);
            }
        }

        // print the panel, which contains matrixIDs for these matrices
        // For ROW_VECTORS, COL_VECTORS and MATRICES there is a panel[0]
        fprintf(f, "    \"%" PRId64 "\":[%d, %d, %" PRId64 ", ",
            matID_map[s], s_ptr->row_level, s_ptr->col_level,
            matID_map[MID_FROM_PID(s_ptr->subMatList[0])]);

        // For ROW_VECTORS, MATRICES there is a panel[1]
        // For COL_VECTORS print -1
        if (s_ptr->subMatList[1]!=MATRIX_ID_INVALID)
            fprintf(f, "%" PRId64 ", ",
                matID_map[MID_FROM_PID(s_ptr->subMatList[1])]);
        else { fprintf(f, "-1, "); }

        // For COL_VECTORS, MATRICES there is a panel[2]
        // For ROW_VECTORS print -1
        if (s_ptr->subMatList[2]!=MATRIX_ID_INVALID)
            fprintf(f, "%" PRId64 ", ",
                matID_map[MID_FROM_PID(s_ptr->subMatList[2])]);
        else { fprintf(f, "-1, "); }

        // For MATRICES there is a panel[3]
        // For ROW_VECTORS or COL_VECTORS print -1
        if (s_ptr->subMatList[3]!=MATRIX_ID_INVALID)
            fprintf(f, "%" PRId64 "],\n",
                matID_map[MID_FROM_PID(s_ptr->subMatList[3])]);
        else { fprintf(f, "-1],\n"); }

        // Set to zero so as not to print this entry out twice.
        usasPTR->array[s] = 0;
    } // end nonscalar branch
}


/*!
 * \ingroup larc
 * \brief This is an internal subroutine to write a LARCmatrix to a file
 *   with options to control the output order and whether or not the matrix IDs
 *   are renumbered, and returns the LARCsize.
 *
 * \param m_pID The packedID of the matrix to be written in LARC format
 * \param path The location of the new larcMatrix json file
 * \param use_topLeft If non-zero, use top left order
 * \param do_renumber If non-zero, renumber the matrix IDs on output
 * \return The larcSize of the matrix with packedID m_pID
 */
static size_t fprint_larcMatrixFile_common(int64_t m_pID, char *path, int use_topLeft, int do_renumber)
{
    if (m_pID == MATRIX_ID_INVALID)
    {
        fprintf(stderr,"in %s, invalid matrixID passed\n",__func__);
        return 0;
    }

    int local_debug = 0;
    if (local_debug) 
        fprintf(stdout, "Inside %s, about to call create_usasPTR\n",__func__);

    // create a uniq submatrix array structure (usas) for that matrix
    usas_t *usasPTR = create_usasPTR(m_pID);
    if (local_debug) 
        fprintf(stdout, "Inside %s, about to call fill_usasPTR\n",__func__);
    fill_usasPTR(usasPTR,m_pID);

    int64_t matrixID = MID_FROM_PID(m_pID);

    // if (VERBOSE >= DEBUG)
    if (local_debug) 
        fprintf(stdout, "About to load info store\n");

    // load the info store with counts from the usas, Note: the format
    // for matrixIDs limits the number of uniq submatrices to 2^50 < 10^20
    char larcSize_str[20];
    sprintf(larcSize_str,"%ld", usasPTR->larcSizeCounter);
    char uniqScalarCount_str[20];
    sprintf(uniqScalarCount_str,"%ld", usasPTR->numUniqScalars);
    info_set(LARC_SIZE, m_pID, larcSize_str);
    info_set(UNIQ_SCALARS, m_pID, uniqScalarCount_str);

    // add meta data to the info store on compile time scalarType , width info
    update_infoStore_automatic_entries(m_pID);

    // the matID_map array is indexed by the matrixIDs found in the input file
    // but contains the packedIDs assigned when added to the MatrixStore
    int64_t *matID_map = (int64_t *) calloc(matrixID + 1, sizeof(int64_t));
     
    // open a file to write into
    FILE *f = fopen(path, "w");
  
    if (!f)
    {
        fprintf(stderr,"ERROR in %s, could not open %s for writing\n",
            __func__,path);
        return 0;
    }

    if (VERBOSE > BASIC)
    {
        fprintf(stdout,"\nWriting matrix with matrixID %" PRId64 " to %s\n",
            matrixID, path);
        fprintf(stdout,"The larcSize is %zd and numUniqScalars is %zd.\n",
                 usasPTR->larcSizeCounter, usasPTR->numUniqScalars);
    }
  
    // json file contains: matrixID_max, matid, info struct, table struct

    // print the matrixID_max, and matid to the file
    if (use_topLeft) {
        fprintf(f, "{\n  \"matrixID_max\":%" PRIu64 ",\n  \"matid\":%" PRIu64 ",",
              matrixID+1, matrixID);
    } else {
        fprintf(f, "{\n  \"matrixID_max\":%" PRIu64 ",\n  \"matid\":%" PRIu64 ",",
              usasPTR->larcSizeCounter, usasPTR->larcSizeCounter-1);
    }

    // print the InfoStore to the file
    int info_ret = write_infoStore_to_larcMatrix_file(m_pID, f);
    if (info_ret != 0) {
        fprintf(stderr,"ERROR in %s, failure in function\n",__func__);
        fprintf(stderr,"\t write_infoStore_to_larcMatrix_file\n");
    }

    // print header for the table of submatrices contained in m
    fprintf(f, "\n  \"table\":{\n");

    // map_counter is the next matID that can be used in our downsized mapping
    int64_t map_counter = 0;

    // print every submatrix s of matrix m (including itself)
    // this corresponds to everything where array value is 1 (for scalar)
    // or 2 (for matrix).
    // Note that the array values are set back to 0 upon output, to avoid
    // writing out the same entry twice.
    if (use_topLeft) {
        recursive_larcMatrixFile_element_output(f, matrixID, usasPTR, matID_map, &map_counter, do_renumber);
    } else {
        // s is a matrixID not a packedID, and we will need to keep that in mind
        for (int64_t s=0; s <= matrixID; ++s) {
            recursive_larcMatrixFile_element_output(f, s, usasPTR, matID_map, &map_counter, do_renumber);
        }
    }

    // end the table structure and the json file
    fprintf(f, "      \"end\":0 }\n}\n");
    fclose(f);

    if (VERBOSE>SILENT) printf("wrote file %s\n",path);

    // grab the larcSize then free the usas
    size_t larcSize = usasPTR->larcSizeCounter;
    free_usasPTR(usasPTR);
    free(matID_map);

    return(larcSize);
}

// See io.h for doxygen comment.
size_t fprint_larcMatrixFile(int64_t m_pID, char *path)
{
    size_t larcSize;

    larcSize = fprint_larcMatrixFile_common(m_pID, path, 0, 1);
    return(larcSize);
}


// See io.h for doxygen comment.
size_t fprint_larcMatrixFile_topLeftOrder(int64_t m_pID, char *path)
{
    size_t larcSize;

    larcSize = fprint_larcMatrixFile_common(m_pID, path, 1, 0);
    return(larcSize);
}


size_t fprint_uniqScalar_file(int64_t m_pID, char *path)
{

#ifndef IS_COMPLEX
    fprintf(stderr,"The routine %s is only usable for complex types\n",__func__);
    return(0);
#endif    

    if (m_pID == MATRIX_ID_INVALID)
    {
        fprintf(stderr,"in %s, invalid matrixID passed\n",__func__);
        return 0;
    }

    int64_t matrixID = MID_FROM_PID(m_pID);
     
    int local_debug = 0;
    if (local_debug) 
          fprintf(stdout, "Inside %s, about to call create_usasPTR\n",__func__);

    // create a uniq submatrix array structure (usas) for that matrix
    usas_t *usasPTR = create_usasPTR(m_pID);
    fill_usasPTR(usasPTR,m_pID);

     // open a file to write into
     FILE *f = fopen(path, "w");
  
     if (!f)
    {
        fprintf(stderr,"ERROR in %s, could not open %s for writing\n",__func__,path);
        return 0;
    }

     if (VERBOSE > BASIC) {
         fprintf(stdout,"\nWriting matrix with matrixID %" PRId64 " to %s\n",
               matrixID, path);
         fprintf(stdout,"The larcSize is %zd and numUniqScalars is %zd.\n",
                 usasPTR->larcSizeCounter, usasPTR->numUniqScalars);
    }

    for (int64_t s=0; s <= matrixID;s++) {
         // if the array value is zero then this matID is not in top matrix
         // if the array value is two then this matID is not for a scalar submatrix
         if (usasPTR->array[s] != 1)  continue;

	 // now we are only in the case where s corresponds to a scalar

	    // TODO: eventually write scalars.c code to do 
	    //  char *val_string = sca_get_column_format_str(s_ptr->scalar_value);
            //  fprintf(f, "%s\n", val_string);
	    //  free(val_string);
	    // CAN COMPARE TO non col format
            // char *val_string = sca_get_readable_approx_str(s_ptr->scalar_value);
	    //  fprintf(f, "%s\n", val_string);
	    //  free(val_string);
           #ifdef USE_COMPLEX	    
                 mats_ptr_t s_ptr = (mats_ptr_t)get_recordPTR_from_pID(
                      PID_FROM_SCALAR_MID(s),"",__func__,0);
                 long double real_part =  creall(s_ptr->scalar_value);
		 long double imag_part =  cimagl(s_ptr->scalar_value);
		 fprintf(f, "%Lg %Lg\n",real_part,imag_part);
           #else	   //  not a handled type
	          fprintf(stderr,"This code only works for type COMPLEX\n");
	          exit(0);
           #endif   // finished cases

  }  // end usasPTR->array loop
  fclose(f);

     // create a new path with .info at the end.
     char info_path[300];
     sprintf(info_path,"%s.info",path);
  
     // open a file to write into
     FILE *f_info = fopen(info_path, "w");
  
     if (!f_info)
    {
        fprintf(stderr,"ERROR in %s, could not open %s for writing\n",__func__,info_path);
        return 0;
    }

     if (VERBOSE > BASIC) {
         fprintf(stdout,"\nWriting associated infostore metadata to %s\n", info_path);
    }

     // if (VERBOSE >= DEBUG)
     if (local_debug)    fprintf(stdout, "About to load info store\n");
     // load the info store with counts from the usas, Note: the format
     // for matrixIDs limits the number of uniq submatrices to 2^50 < 10^20
     char larcSize_str[20];
     sprintf(larcSize_str,"%ld",usasPTR->larcSizeCounter);
     char uniqScalarCount_str[20];
     sprintf(uniqScalarCount_str,"%ld",usasPTR->numUniqScalars);
     info_set(LARC_SIZE, m_pID,larcSize_str);
     info_set(UNIQ_SCALARS, m_pID,uniqScalarCount_str);

     // add meta data to the info store on compile time scalarType , width info
     update_infoStore_automatic_entries(m_pID);

    // 
    fprintf(f_info,"The scalar file can be sorted using unix: 'sort -k1g -k2g'\n");
    fprintf(f_info,"Here is the associated metadata from the InfoStore:\n");

     
     
    // print the InfoStore to the file
    int info_ret = write_infoStore_to_larcMatrix_file(m_pID, f_info);
    if (info_ret != 0) {
        fprintf(stderr,"ERROR in %s, failure in function\n",__func__);
        fprintf(stderr,"\t write_infoStore_to_larcMatrix_file\n");
    }
    fclose(f_info);

     // Save off the number of unique scalars and then free the usas
     size_t nus = usasPTR->numUniqScalars;
     free_usasPTR(usasPTR);

     return(nus);
}
