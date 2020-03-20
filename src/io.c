//                            io.c 
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


#include <inttypes.h> //for printing uint64_t and int64_t
#include <curses.h>
#include "global.h"
#include "larc.h"
#include "io.h"


// FOUR ROUTINES FOR READING IN FROM ROW_MAJOR_LIST FORMAT:
//
// These two routines work with multiprecision data.
// int64_t row_major_list_to_store_matrixID(char **dense_mat, 
//                    mat_level_t        current_row_level, 
//                    mat_level_t        current_col_level, 
//                    int64_t        orig_num_cols 
//                    )
// mat_ptr_t row_major_list_to_store(scalarType  * dense_mat, 
//                    mat_level_t        current_row_level, 
//                    mat_level_t        current_col_level, 
//                    int64_t        orig_num_cols 
//                    )
//
// These two routines do NOT work with multiprecision data.
// int64_t read_row_major_matrix_from_file_matrixID(char * file_path)
// mat_ptr_t read_row_major_matrix_from_file(char * file_path)
//
// The *_matrixID routines are wrappers for the ones that return type
// mat_ptr_t (pointers to matrices). They are used in our python interface
// so we can avoid passing structure pointers, including those which hold
// GMP multiprecision data.

// row_major_list_to_store_matrixID():
// This wrapper routine takes as its first argument an array of strings,
// which will be converted to the proper scalarType value before calling 
// the core function row_major_list_to_store().
int64_t row_major_list_to_store_matrixID(char **dense_mat, 
                    mat_level_t        current_row_level, 
                    mat_level_t        current_col_level, 
                    int64_t        orig_num_cols 
                    )
{

  // allocate sufficient space to store the matrix "dense" in memory
  int64_t N = (1L << current_row_level) * (1L << current_col_level);
  scalarType *dense = malloc(N*sizeof(scalarType));

  // convert each string value in dense_mat to the proper scalar type
  // and put it into our "dense" matrix
  for (int64_t i = 0; i < N; i++){
    sca_init(&(dense[i]));
    sca_set_str(&(dense[i]), dense_mat[i]);
  }

  // call row_major_list_to_store to put "dense" into the matrixStore
  mat_ptr_t C_ptr = row_major_list_to_store(dense, current_row_level, current_col_level, orig_num_cols); 
  int64_t C_mID = get_matID_from_matPTR(C_ptr);

  // clear the memory used to store "dense"
  for (int64_t i = 0; i < N; i++)
    sca_clear(&(dense[i]));
  free(dense);

  return C_mID;
}


/****************************************************************
 *                     row_major_list_to_store
 *    jszito note: John Gilbert taught me this neat recursive method;
 *                 he says comes from the old Fortran blas routines.
 *    It takes a dense m by n matrix which is given as a list of all 
 *    its entries in row major order: 
 *                   a00 a01 ... a0(n-1) a10 a11 ... a1(n-1) ...  
 *                        a(m-1)0 a(m-1)1 ... a(m-1)(n-1)
 *    And specifies a submatrix by giving its top left corner aij
 *    and its dimension.
 *    This routine works recursively on the 4 subpanels
 *    It loads the matrix into the matrix store c using get_matPTR_from_array_of_four_subMatPTRs
 *    which returns the matrix index of the matrix in the store.
 *
 *    This function gets called from the top level of an m x n 
 *    matrix and acts on dense_mat which is an m*n long array 
 *    of complex values
 *    and is the row major ordering of a complex matrix
 *                   a00 a01 ... a0(n-1) a10 a11 ... a1(n-1) ...  
 *                        a(m-1)0 a(m-1)1 ... a(m-1)(n-1)
 *    where m = 2^row_level, n= 2^col_level, and 
 *    orig_num_cols = n
 *
 *    Any sub matrix in the original matrix can be specified by 
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
 *                  submatrix[1] starts at 8, submatrix[2] starts at 2,
 *                  and submatrix[3] starts at 10, and
 *       the subarray of dense_mat that starts at that index.
 *
 *    The recursive nature of the call takes care of everything else.
 *    When get_matPTR_from_array_of_four_subMatPTRs has the mat_ptr of each of the four submatrices
 *    (in the vector cases two of these are MATRIX_PTR_INVALID)
 *    then it can be used to find the mat_ptr of a matrix.
 *    
 **********************************************************************/
mat_ptr_t row_major_list_to_store(scalarType  * dense_mat, 
                    mat_level_t        current_row_level, 
                    mat_level_t        current_col_level, 
                    int64_t        orig_num_cols 
                    )
{
#ifdef DEBUG
  printf("--> %s : %s: %d\n", __FILE__, __func__, __LINE__);
#endif

  // do we want to add a check to see if dense_mat list has
  // 2^col_level*2^row_level entries?
  // SAC: no, this is a recursive routine and the check should be outside

  
  // This printing for verbose has been updated for nonsquare matrices
  // but not for scalarType handling. Note that this is inside a recursive
  // routine, so the printing is likely to be messy.
  int verbose = 0;
#if 0
  if (verbose) { 
    printf("\nIn subroutine row_major_list_to_store, with row level ");
    printf("%d and column level %d\n",current_row_level,current_col_level);
    printf("and orig_num_cols %ld\n",orig_num_cols);
    int i,j;
    for (i=0; i < (1L << current_row_level); ++i) {
      for (j=0; j < (1L << current_col_level); ++j) {
#ifdef USE_COMPLEX
        printf("%4.1lf + %4.1lf i, \t ",creall(dense_mat[i*orig_num_cols+j]),
                               cimagl(dense_mat[i*orig_num_cols+j]));
#else
        printf("%4.1lf, \t ",dense_mat[i*orig_num_cols+j]);
#endif
      }
      printf("\n");
    }
  }  // end verbose
#endif
    
  // SCALAR CASE
  if ((current_row_level == 0) && (current_col_level == 0)){
    return get_valMatPTR_from_val(*dense_mat);
  }

  // Need this for remaining cases
  mat_ptr_t panel[4];  // indices of submatrices


  // ROW VECTOR
  if ((current_row_level == 0) && (current_col_level > 0)){
    int64_t new_num_cols = 1L << (current_col_level-1);   // divide current number by two

    // top corner
    panel[0] = row_major_list_to_store(dense_mat,0,current_col_level-1,orig_num_cols);

    // step right dim_new cols
    panel[1] = row_major_list_to_store(dense_mat+new_num_cols,0, current_col_level-1,orig_num_cols); 

    // step down dim_new rows
    panel[2] = MATRIX_PTR_INVALID;

    // step down and over dim_new rows and cols
    panel[3] = MATRIX_PTR_INVALID;

#ifdef DEBUG
  printf("<--%s : %s: %d\n", __FILE__, __func__,  __LINE__);
#endif  

  }


  // COLUMN VECTOR
  if ((current_row_level > 0) && (current_col_level == 0)){
    int64_t new_num_rows = 1L << (current_row_level-1);   // divide current number by two

    // top corner
    panel[0] = row_major_list_to_store(dense_mat,current_row_level-1,0,orig_num_cols);

    // step right dim_new cols
    panel[1] = MATRIX_PTR_INVALID;

    // step down dim_new rows
    panel[2] = row_major_list_to_store(dense_mat+new_num_rows*orig_num_cols,current_row_level-1,0,orig_num_cols); 

    // step down and over dim_new rows and cols
    panel[3] = MATRIX_PTR_INVALID;

#ifdef DEBUG
  printf("<--%s : %s: %d\n", __FILE__, __func__,  __LINE__);
#endif  

  }

  // MATRIX CASE
  if ((current_row_level > 0) && (current_col_level > 0)){
    int64_t new_num_cols = 1L << (current_col_level-1);   // divide current number by two
    int64_t new_num_rows = 1L << (current_row_level-1);   // divide current number by two

    if (verbose) {
      printf("\nIn recursive call with orig_num_cols = %"PRId64"\n",orig_num_cols);
      printf("  The four subpanels will be at row level %"PRId32", col level %"PRId32", ",
                current_row_level-1, current_col_level-1);
      printf("and will have offsets: \n");
      printf("     panel[0]: 0\n");
      printf("     panel[1]: new_num_cols = %ld\n",new_num_cols);
      printf("     panel[2]: new_num_rows*orig_num_cols = %ld\n",
                new_num_rows*orig_num_cols);
      printf("     panel[3]: new_num_rows*orig_num_cols+new_num_cols = %ld\n",
                new_num_rows*orig_num_cols+new_num_cols);
    }  // end verbose

    // top corner
    panel[0] = row_major_list_to_store(dense_mat, current_row_level-1,
                                       current_col_level-1, orig_num_cols);

    // step right dim_new cols
    panel[1] = row_major_list_to_store(dense_mat+new_num_cols, current_row_level-1, 
                                       current_col_level-1, orig_num_cols); 

    // step down dim_new rows
    panel[2] = row_major_list_to_store(dense_mat+new_num_rows*orig_num_cols, current_row_level-1, 
                                       current_col_level-1, orig_num_cols); 

    // step down and over dim_new rows and cols
    panel[3] = row_major_list_to_store(dense_mat+new_num_rows*orig_num_cols+new_num_cols,
                                       current_row_level-1, current_col_level-1, orig_num_cols); 
#ifdef DEBUG
  printf("<--%s : %s: %d\n", __FILE__, __func__, __LINE__);
#endif  

  }

  return get_matPTR_from_array_of_four_subMatPTRs(panel, current_row_level, current_col_level);
}

/****************************************************************
 *  python interface: read_row_major_matrix_from_file_matrixID
 ****************************************************************/
int64_t read_row_major_matrix_from_file_matrixID(char * file_path)
{
  mat_ptr_t m_ptr = read_row_major_matrix_from_file(file_path);
  int64_t m_mID = get_matID_from_matPTR(m_ptr);
  return m_mID;
}


/****************************************************************
 *  read_row_major_matrix_from_file
 *  This routine should now handle multiprecision data!
 ****************************************************************/
mat_ptr_t read_row_major_matrix_from_file(char * file_path)
{
#ifdef DEBUG
  printf("--> %s : %s: %d\n", __FILE__, __func__, __LINE__);
#endif  
  int verbose = 0;
  int error_code = 0;
  int errno = 0; 

  if (verbose) printf("In %s\n", __func__);

  FILE *fp;
  if((fp = fopen(file_path, "r+")) == NULL) {
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
#ifdef DEBUG
    printf("<-- %s : fail: %d\n", __FILE__, __LINE__);
#endif  
    return (MATRIX_PTR_INVALID); 
  } 

  // if the file is open but nothing is left to read print EOF
  if(feof(fp)){
    fputs("EOF",stderr);
    fprintf(stderr,"ERROR: File stopped earlier than expected.\n");
    fclose(fp);
#ifdef DEBUG
    printf("<-- %s : fail: %d\n", __FILE__, __LINE__);
#endif  
    return (MATRIX_PTR_INVALID); 
  }

  // successfully read the level
  if (verbose)   printf("  levels are %d %d\n",(int)row_level, (int)col_level);

  // calculate size of the array
  int64_t array_size = (1L << (row_level+col_level));
  if (verbose)  printf("  array length is %ld\n", array_size);

  // Now we will handle the entries of the matrix
  scalarType * matrix_array;

  // Allocate space for the array of scalarType
  if (!(matrix_array = calloc(array_size,sizeof(scalarType)) )) {
    error_code = 1;  // unable to allocate space for array
    fprintf(stderr,"ERROR: Unable to allocate space for array.\n");
#ifdef DEBUG
    printf("<-- %s : fail: %d\n", __FILE__, __LINE__);
#endif  
    return (MATRIX_PTR_INVALID); 
  } 

  // Read the array of numbers 
  int max_size = 1000;
  char a[max_size];
  int64_t i;
  for (i=0; i<array_size; ++i) {
     // format for numbers must not include whitespace!
     int ret = fscanf(fp,"%s\n",a);
     if (ret == EOF) {
        fprintf(stderr,"ERROR: fscanf() failed to read a line (i=%ld)",i);
        break;
     }
     else if (errno != 0) {
        perror("fscanf");
        fprintf(stderr,"(i=%ld)\n",i);
        break;
     }
     else if (ret == max_size) {
        fprintf(stderr,"ERROR: line read by fscanf() may exceed ");
        fprintf(stderr,"%d characters (i=%ld)\n",max_size,i);
     }
     sca_init(&(matrix_array[i]));
     sca_set_str(&(matrix_array[i]), a);
  }

  // if the file is open but nothing is left to read print EOF
  if(feof(fp)){
    fputs("EOF",stderr);
  }

  fclose(fp);

  mat_ptr_t ret_ptr;
  if (array_size == i) {
#ifdef DEBUG
    printf("<-- %s : %s success: %d\n", __FILE__, __func__,__LINE__);
#endif  
    ret_ptr = row_major_list_to_store(matrix_array,row_level,col_level,
						1L<<col_level);
  }
  else {
    if (verbose) printf("Failed to read row major matrix from file.\n");
#ifdef DEBUG
    printf("<-- %s : %s fail: %d\n", __FILE__, __func__, __LINE__);
#endif  
    ret_ptr = MATRIX_PTR_INVALID;
  }

  // clean up any space allocated by sca_init
  for (int64_t j=0; j<i; ++j) sca_clear(&(matrix_array[j]));
  // clean up scalarType array
  free(matrix_array);

  return ret_ptr;
}  // end read_row_major_matrix_from_file


// FOUR NAIVE MATRIX PRINTING ROUTINES
// * these routines print out the full matrix with no compression, and must
// * be used carefully...
// void print_naive_by_matID(int64_t A_mID)
// void print_naive_by_matPTR(mat_ptr_t ptr_A)
// void write_naive_by_matID(int64_t A_mID, char * filename)
// void write_naive_by_matPTR(mat_ptr_t ptr_A, char * filename)

/*******************************************************************
*                       (naive) print_naive_by_matID       *
* Python Interface for print_naive_by_matPTR      			*
* accepts a matrixID for matrix to print                           *
* Converts to a pointer before calling								*
********************************************************************/
void print_naive_by_matID(int64_t A_mID)
{  
  // get the matrix pointer from the matrixID, and see if still in store 
  mat_ptr_t A_ptr = get_matPTR_from_matID(A_mID, "", __func__,0);

  if (A_ptr == MATRIX_PTR_INVALID) { 
      fprintf(stderr,"MatrixID %ld is INVALID.\n", A_mID);
      exit(1); 
  }

  // should have a valid matrix pointer, to call internal routine for printing
  print_naive_by_matPTR(A_ptr);
}


/********************************************************************
*                       (naive) print_matrix  (to screen)           *
********************************************************************/
void print_naive_by_matPTR(mat_ptr_t ptr_A)
{
  int i,j;
  int row_level =  matrix_row_level(ptr_A);
  int col_level =  matrix_col_level(ptr_A);
  unsigned row_dim = 1 << row_level;
  unsigned col_dim = 1 << col_level;

  // limit naive printing to reasonable sizes
  if (!row_dim || !col_dim)
  {
        fprintf(stderr,"in print_naive_by_matPTR, dimension of matrix too big\n");
        fprintf(stderr,"to fit into a 32-bit integer - we're not printing that!\n");
        return;
  }

  // before continuing, calculate number of characters to be printed
  // (this could double the time required for printing, we decided that was OK)
  scalarType *entry;
  char *entry_string;
  uint64_t c_counter, max_columns, max_chars;

  max_chars = max_columns = 0;
  for(i = 0;i<row_dim;i++) {
    c_counter = 0;
    for(j = 0;j<col_dim; j++) {
      entry = &matrix_trace(get_valMatPTR_from_matPTR_and_coords(i, j, ptr_A));
      entry_string = sca_get_str(*entry);
      c_counter += strlen(entry_string)+1;
      free(entry_string);
    }
    c_counter += 2;
    max_chars += c_counter;
    max_columns = MAX(max_columns,c_counter);
  }

  if (VERBOSE>BASIC) {
    printf("%s given matrix which, if printed, would require\n", __func__);
    printf("%lu characters total, with a maximum of \n", max_chars);
    printf("%lu characters per row of the matrix\n", max_columns);
  }

  // set the maximum level for printing a matrix
  // these numbers are intentionally set larger than what we think
  // would be reasonable for a matrix printed to a terminal window
# if defined(USE_INTEGER) || defined(USE_REAL)
  int max_print_level = 10;
# elif defined(USE_COMPLEX)
  int max_print_level = 7;
# else // multiprecision type
  int max_print_level = 6;
#endif

  // limit naive printing to reasonable size
  if (row_level>max_print_level || col_level>max_print_level ||
                max_chars>100000)
  {
        fprintf(stderr,"in %s, matrix has dimensions",__func__);
        fprintf(stderr," %u x %u\n",row_dim,col_dim);
        fprintf(stderr,"and would require %lu printed characters,\n",max_chars);
        fprintf(stderr,"which is too large to naive print to screen. ");
        fprintf(stderr,"If you must do this anyway,\nuse write_naive routine");
        fprintf(stderr," with 'stdout' as output filename\n");
        return;
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
    entry = &matrix_trace(get_valMatPTR_from_matPTR_and_coords(i, 0, ptr_A));
    entry_string = sca_get_str(*entry);
    // always print first element in row without checking line length
    printf("%s ", entry_string);
    c_counter = strlen(entry_string)+1;
    free(entry_string);

    for(j = 1; j<col_dim; j++) {
      // find the value to be printed in string form
      entry = &matrix_trace(get_valMatPTR_from_matPTR_and_coords(i, j, ptr_A));
      entry_string = sca_get_str(*entry);
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
}



/********************************************************************
* Python interface version of (naive) print_matrix_to_file          *
*   this could be implemented less naively ...                      *
********************************************************************/
void write_naive_by_matID(int64_t A_mID, char * filename)
{
  // get the matrix pointer from the matrixID, and see if still in store 
  mat_ptr_t A_ptr = get_matPTR_from_matID(A_mID, "", __func__,0);
  if (A_ptr == MATRIX_PTR_INVALID) { exit(1); }

  // should have a valid matrix pointer, to call internal routine for printing matrix to a file
  write_naive_by_matPTR(A_ptr, filename);

}



/********************************************************************
*                       (naive) print_matrix_to_file                *
*   this could be implemented less naively ...                      *
********************************************************************/
void write_naive_by_matPTR(mat_ptr_t ptr_A, char * filename)
{
  FILE *f;
  // remember strcmp returns 0 on equality
  if (strcmp(filename,"stdout")) f = fopen(filename, "w");
  else f = stdout;
  int i,j;
  int row_level =  matrix_row_level(ptr_A);
  int col_level =  matrix_col_level(ptr_A);
  int row_dim = 1 << row_level;
  int col_dim = 1 << col_level;
  scalarType *entry;

  for(i = 0;i<row_dim;i++){
    for(j = 0;j<col_dim; j++){
      entry = &matrix_trace(get_valMatPTR_from_matPTR_and_coords(i, j, ptr_A));
      char *entry_string = sca_get_str(*entry);
      fprintf(f,"%s ", entry_string);
      free(entry_string);
    }
    fprintf(f," \n");
  }
  if (strcmp(filename,"stdout")) fclose(f);
}

// SOMEWHAT LESS NAIVE PRINTING ROUTINES
// void write_matrix_nonzeros_by_matID(int64_t A_mID, char* filename)
// void write_matrix_nonzeros_by_matPTR(mat_ptr_t ptr_A, char * filename)
// * there is no limit to the number of nonzero values printed...

/********************************************************************
* (Python interface) write_matrix_nonzeros_by_matID            *
********************************************************************/
void write_matrix_nonzeros_by_matID(int64_t A_mID, char * filename)
{
  // get the matrix pointer from the matrixID, and see if still in store 
  mat_ptr_t A_ptr = get_matPTR_from_matID(A_mID, "", __func__,0);
  if (A_ptr == MATRIX_PTR_INVALID) { exit(1); }

  write_matrix_nonzeros_by_matPTR(A_ptr, filename);

}


/********************************************************************
*                write_matrix_nonzeros_by_matPTR               *
********************************************************************/
void write_matrix_nonzeros_by_matPTR(mat_ptr_t ptr_A, char * filename)
{
  FILE *f = fopen(filename, "w");
  mat_level_t row_level =  matrix_row_level(ptr_A);
  mat_level_t col_level =  matrix_col_level(ptr_A);
  int64_t row_dim = 1L << row_level;
  int64_t col_dim = 1L << col_level;
  scalarType *entry;

  for(int64_t i = 0; i < row_dim; i++){
    for(int64_t j = 0; j < col_dim; j++){
      entry = &matrix_trace(get_valMatPTR_from_matPTR_and_coords(i, j, ptr_A));
      if (0 == sca_eq(*entry, scalar0)){
        char *entry_string = sca_get_str(*entry);
        fprintf(f, "%s\t(%ld,%ld)\n", entry_string, i, j);
        free(entry_string);
      }
    }
  }
  fclose(f);
}

// A COMPARISON BETWEEN THE MATRICES STORED IN TWO COMPRESSED JSON FILES
// This (apparently) is in io.c since data is read from files...
// Since the routine calls read_larcMatrix_file_return_matPTR, it works with multiprecision
// types and can be called from the python interface. It puts both matrices
// into the matrix store, and returns 0 if the matrices are not equal or the
// matrixID if they are equal.
int64_t
equal_matrices_in_larcMatrix_files(char *path1, char *path2)
{
  if (!strcmp(path1,path2)) /* paths are equal */
    {
      fprintf(stderr,"Warning %s : %s fail: line %d, path1 and path2 are the same.\n",
              __FILE__, __func__, __LINE__);
    }
  mat_ptr_t m1_ptr = read_larcMatrix_file_return_matPTR(path1);
  mat_ptr_t m2_ptr = read_larcMatrix_file_return_matPTR(path2);
  if ( matrix_is_invalid(m1_ptr) || matrix_is_invalid(m2_ptr) ) {
      fprintf(stderr,"Warning %s : %s fail: line %d, invalid matrix\n",
              __FILE__, __func__, __LINE__);
      return(0);
  }
  int64_t m1_mID =  get_matID_from_matPTR(m1_ptr);
  int64_t m2_mID =  get_matID_from_matPTR(m2_ptr);
  if (m1_mID == m2_mID) {
    return m1_mID;
  }
  else return 0;
}

// ROUTINES WHICH EXECUTE READING AND WRITING FROM/TO JSON COMPRESSED FORMAT
// 
// THE WRITING ROUTINES:
// int write_larcMatrix_file_by_matID(int64_t m_mID, char *path)
// int write_larcMatrix_file_by_matPTR(mat_ptr_t m_ptr, char *path)
// static int recursive_write_larcMatrix_file_by_matPTR(
//    mat_ptr_t m_ptr, FILE *f, char *output_flags,
//    void (*func)(scalarType*, const scalarType))
// static int write_infoStore_to_larcMatrix_file(mat_ptr_t m_ptr, FILE *f)
// int write_and_alter_vals_larcMatrix_file(mat_ptr_t m_ptr, char *path,
//      void (*func)(scalarType*, const scalarType))

// The core routine is write_and_alter_vals_larcMatrix_file(),  which
// does some testing, writes the json file header, calls
// write_infoStore_to_larcMatrix_file() to write out any info_store information, then
// calls the recursive routine recursive_write_larcMatrix_file_by_matPTR() to output the
// data in the matrix. The matrixID values are output unchanged. Each scalar
// value seen has func() applied to it, and the resulting scalarType value is
// converted to a string before being output. Once that is completed, the file
// footer is written and the file handle closed, and finally the function
// returns to its calling routine.
//
// We have not written a matrixID version for this routine, though it would be
// trivial, because it should always be called by some other C routine which
// defines a value for the function pointer
//    void *(func)(scalarType *a, const scalarType b)
// and the calling routine would pass the pointer to the matrix rather than a
// matrixID.
//
// func() is applied to every scalar (b) that is in the matrix before the
// modified matrix is written. To put the unchanged scalar b in the location a,
// the function we would pass is sca_set(), which does *a=b. This is hardcoded
// into the "default" writing routines write_larcMatrix_file_by_matPTR() and its
// *_matrixID version.

int write_larcMatrix_file_by_matID(int64_t m_mID, char *path)
{
  // get the matrix pointer from the matrixID, and see if still in store 
  mat_ptr_t m_ptr = get_matPTR_from_matID(m_mID, "", __func__,0);
  if (m_ptr == MATRIX_PTR_INVALID) { exit(1); }

  // calculate matrix pointer version of function
  
  return write_larcMatrix_file_by_matPTR(m_ptr, path);
}

int write_larcMatrix_file_by_matPTR(mat_ptr_t m_ptr, char *path)
{
    // the function sca_set as the last argument ensures the written data
    // is the same as that in the matrixStore
    return write_and_alter_vals_larcMatrix_file(m_ptr, path, sca_set);
}

/*!
 * \ingroup larc
 * \brief A worker routine for writing a json formatted compressed matrix recursively
 * \param m_ptr The pointer to the matrix to be written
 * \param f A file pointer
 * \param output_flags An array of flags used to track already-written matrixIDs
 * \param func A function that will be applied to each scalar before writing
 * \return 1 on success
 */
static int
recursive_write_larcMatrix_file_by_matPTR(mat_ptr_t m_ptr, FILE *f, char *output_flags, void (*func)(scalarType*, const scalarType))
{
  if (matrix_is_invalid(m_ptr)) {
    fprintf(stderr,"%s: invalid matrix pointer passed\n", __func__);
    return -1;
  }

  /* output_flags keeps track of which matrixIDs have already appeared */
  if (output_flags[get_matID_from_matPTR(m_ptr)] > 0)
    return 1;

  output_flags[get_matID_from_matPTR(m_ptr)] = 1;
  matrix_type_t mtype = matrix_type(m_ptr);

  if (mtype == SCALAR)
    {
      // in COMPLEX case, write routine to output "a+I*b" as "a, b"
      scalarType *val = &scratchVars.submit_to_store;
      // the function func() is applied to the scalar value in m_ptr, and
      // the result put into val; if the function is sca_set, the stored value
      // is merely copied
      func(val, matrix_trace(m_ptr));
      char *val_string = sca_get_str(*val);
      fprintf(f, "    \"%ld\":[%d, %d, \"%s\"],\n", get_matID_from_matPTR(m_ptr), matrix_row_level(m_ptr), 
              matrix_col_level(m_ptr), val_string);
      free(val_string);
    }
  else   // NONSCALAR
    {
      int ret;
     
      mat_ptr_t panel[4];
      for (int i = 0; i < 4; ++i) {
        panel[i] = matrix_sub(m_ptr,i);
      }

      // obtain the pointers for the matrices in the panel
      ret = recursive_write_larcMatrix_file_by_matPTR(panel[0], f, output_flags, func);
      if (ret < 0) return ret;

      if (mtype!=COL_VECTOR) {
        ret = recursive_write_larcMatrix_file_by_matPTR(panel[1], f, output_flags, func);
        if (ret < 0) return ret;
      }

      if (mtype!=ROW_VECTOR) {
        ret = recursive_write_larcMatrix_file_by_matPTR(panel[2], f, output_flags, func);
        if (ret < 0) return ret;
      }

      if (mtype==MATRIX) {
        ret = recursive_write_larcMatrix_file_by_matPTR(panel[3], f, output_flags, func);
        if (ret < 0) return ret;
      }

      // print the panel, which contains matrixIDs for these matrices 
      fprintf(f, "    \"%ld\":[%d, %d, %ld, ", get_matID_from_matPTR(m_ptr), 
                matrix_row_level(m_ptr), matrix_col_level(m_ptr),
	        get_matID_from_matPTR(panel[0]));

      if (mtype!=COL_VECTOR) {
                fprintf(f, "%ld, ", get_matID_from_matPTR(panel[1]));
      }
      else { fprintf(f, "-1, "); }

      if (mtype!=ROW_VECTOR) {
                fprintf(f, "%ld, ", get_matID_from_matPTR(panel[2]));
      }
      else { fprintf(f, "-1, "); }

      if (mtype==MATRIX) {
                fprintf(f, "%ld],\n", get_matID_from_matPTR(panel[3]));
      }
      else { fprintf(f, "-1],\n"); }

      //fprintf(f, "    \"%ld\":[%d, %d, %ld, %ld, %ld, %ld],\n", get_matID_from_matPTR(m_ptr), 
      //      matrix_row_level(m_ptr),matrix_col_level(m_ptr),
      //      get_matID_from_matPTR(panel[0]),get_matID_from_matPTR(panel[1]),
      //      get_matID_from_matPTR(panel[2]),get_matID_from_matPTR(panel[3]));
    }
  return 1;
}


/* This function should retrieve all information out of the info_store
   for the given matrix associated with m_ptr and write it to the
   json file associated with the output file pointer f.
*/
/*!
 * \ingroup larc
 * \brief Writes all info_store information about a given matrix to a json file
 * \param m_ptr Pointer to the matrix to be stored in compressed json format
 * \param f File pointer for the json file where the metadata will be written
 * \return 1 on success
 */
static int write_infoStore_to_larcMatrix_file(mat_ptr_t m_ptr, FILE *f) {
  int verbose = 0;

  if (!f)
    return -1;
  if (matrix_is_invalid(m_ptr)) 
    return -1;

  int64_t m_ID = get_matID_from_matPTR(m_ptr);

  enum info_types i;
  int was_info = 0;
  char* info_data;
  char* info_name;
  for (i=0;i<INVALID_INFO;++i) {
    info_data = info_get(i,m_ID);
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

int write_and_alter_vals_larcMatrix_file(mat_ptr_t m_ptr, char *path,
        void (*func)(scalarType*, const scalarType))
{
  /* output_flags keeps track of which matrixIDs have already appeared */

  FILE *f = fopen(path, "w");
  int ret;

  if (!f)
    {
        fprintf(stderr,"in %s, could not open %s for writing\n",__func__,path);
        return -1;
    }
  if (matrix_is_invalid(m_ptr)) 
    {
        fprintf(stderr,"in %s, passed matrix pointer is invalid\n",__func__);
        return -1;
    }
  char *output_flags = calloc(num_matrices_created(), sizeof(char));
  if (output_flags == NULL) {
    ALLOCFAIL();
  }

  if (!output_flags) 
    return -1;

  printf("\nWriting %d,%d-level matrix with matrixID %" PRId64 " to %s\n", 
	 matrix_row_level(m_ptr), matrix_col_level(m_ptr),
         get_matID_from_matPTR(m_ptr), path);
  // printf("      the matrix id is %lu, with level %d %d\n",id, matrix_row_level(m_ptr), matrix_col_level(m_ptr));

  // fprintf(f, "{\n  \"matrixID_max\":%" PRIu64 ",\n  \"matid\":%" PRIu64 ",\n  \"table\":{\n", num_matrices_created(), get_matID_from_matPTR(m_ptr));


  // json file contains: matriID_max, matid, optional info struct, table struct
  fprintf(f, "{\n  \"matrixID_max\":%" PRIu64 ",\n  \"matid\":%" PRIu64 ",",
	  num_matrices_created(), get_matID_from_matPTR(m_ptr));

  // the info structure is printed when it contains information
  ret = write_infoStore_to_larcMatrix_file(m_ptr, f);
  if (ret != 0) {
    fprintf(stderr,"WARNING in %s, function write_infoStore_to_larcMatrix_file failed\n",__func__);
  }

  // always print the matrix table recursively
  fprintf(f, "\n  \"table\":{\n"); 
  ret = recursive_write_larcMatrix_file_by_matPTR(m_ptr, f, output_flags, func);
  free (output_flags);

  // end the table structure and the json file
  fprintf(f, "      \"end\":0 }\n}\n");
  fclose(f);

  return ret;
}

// THE READING ROUTINES:
// int64_t read_larcMatrix_file_return_matID(char *path)
// mat_ptr_t read_larcMatrix_file_return_matPTR(char *path)
// mat_ptr_t read_and_alter_vals_larcMatrix_file_return_matPTR(
//      char *path, void (*func)(scalarType*, const scalarType))
// int64_t read_larcMatrix_file_legacy_return_matID(char *path)
// mat_ptr_t read_larcMatrix_file_legacy_return_matPTR(char *path)
//
// Similar to but less complicated than the writing case. The core function is
// read_and_alter_vals_larcMatrix_file_return_matPTR(). This function reads in the
// string-formatted data from a compressed json matrix file, converts the string
// to the appropriate scalarType, then applies func() to the data before the
// value is put into the matrixStore. The core function is called by the
// "default" reading functions read_larcMatrix_file_return_matPTR and its *_matrixID version
// with the function func() hardcoded to be sca_set(), which does not change
// the data.
//
// There are also _legacy versions of the read function which will correctly
// interpret old-style compressed json in which the data fields are C integers,
// doubles or complex numbers rather than string representations of these
// numbers. (The switch to strings was necessary to enable GMP multiprecision.)
// For obvious reasons, there are no legacy writing routines.

/* Python interface to return a matrixID after reading a json file */
int64_t read_larcMatrix_file_return_matID(char *path)
{
  
  mat_ptr_t m_ptr = read_larcMatrix_file_return_matPTR(path);
  int64_t m_mID =  get_matID_from_matPTR(m_ptr);
  return m_mID;
  
}

mat_ptr_t read_larcMatrix_file_return_matPTR(char *path)
{
// This function is the default version for reading in json formatted
// matrix data, converting the strings in the json to the correct scalarType.
// It makes use of the more general function which allows the user to
// specify a function which is applied to each scalar as it is read in. The
// function sca_set(scalarType *a, const scalarType b) performs *a=b and thus
// leaves the input data unchanged.
  return read_and_alter_vals_larcMatrix_file_return_matPTR(path, sca_set);
}

mat_ptr_t read_and_alter_vals_larcMatrix_file_return_matPTR(char *path, void (*func)(scalarType*, const scalarType))
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
  i = SIZE_MAX/(max_matrixID + 1)/sizeof(mat_ptr_t);
  if (i==0) {
      fprintf(stderr,"Error is %s: matrixID_max too large - try\n", __func__);
      fprintf(stderr,"renumbering the matrix_IDs in %s\n", path);
      fprintf(stderr,"to reduce this value to something more reasonable.\n");
      exit(1);
  }

  mat_ptr_t *map = calloc(max_matrixID+1, sizeof(mat_ptr_t));
  if (!map){
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
	    case 3:			/* expect a scalar link */
	      {
                if (verbose) {printf("In scalar case \n");}
		int row_level = j_get_num64(j_array_index(p, 0));
		int col_level = j_get_num64(j_array_index(p, 1));
                if (!j_is_string(j_array_index(p,2))){
                 fprintf(stderr,"ERROR in %s:\nExpected scalar values to be given as strings.\nIf LARCMatrix file is in 'legacy' format, use utils/canonical_format_json.py\nto convert to LARCMatrix format with scalars saved as strings.\n", __func__);
                    exit(1);
                }
                scalarType *mat_val = &scratchVars.submit_to_store;
                const char *val_str = j_get_string(j_array_index(p,2));
                sca_set_str(mat_val, val_str);
// This is the sole place in this routine where func() is called. If the line
// is omitted or func is set to sca_set, then the scalar that was read in is
// left unchanged.
                func(mat_val, *mat_val);
                if (row_level || col_level) {
		  fprintf(stderr,"error in %s:\n\texpected ", __func__ );
                  fprintf(stderr,"scalar, but levels greater than zero!\n");
                }
		map[index] = get_valMatPTR_from_val(*mat_val);
                if (verbose) {
		  // printf("%"PRId64" %zd %" PRId64 " %d %d %g\n", index, map[index], 
                  char *mat_val_string = sca_get_str(*mat_val);
		  printf("%"PRId64" %p %" PRId64 " %d %d %s\n", index, map[index], 
			 get_matID_from_matPTR((map[index])),
			 row_level, col_level, mat_val_string);
                  free(mat_val_string);
		}
	      }
              break;
	    case 6:			/* expect a 4-tuple */
	      {
                if (verbose) {printf("In non scalar case \n");}
		int row_level = j_get_num64(j_array_index(p, 0));
		int col_level = j_get_num64(j_array_index(p, 1));
                if (verbose) {printf("  levels are %d %d \n",row_level,col_level);}
		mat_ptr_t content[4];
                int64_t temp;
                for (int k=0; k< 4; ++k) {
		   temp = j_get_num64(j_array_index(p, k+2));
                   if (verbose) {printf("temp for k = %d is %ld \n",k,temp);}
                   // The matrix ID -1 is used whenever we have no matrix, 
                   // e.g. in row and col place keepers
                   if (temp == -1) { content[k] = MATRIX_PTR_INVALID;}
                   else {
                     content[k] = map[temp];
		   }
		}
		/* content[0] = map[j_get_num64(j_array_index(p, 2))]; */
		/* content[1] = map[j_get_num64(j_array_index(p, 3))]; */
		/* content[2] = map[j_get_num64(j_array_index(p, 4))]; */
		/* content[3] = map[j_get_num64(j_array_index(p, 5))]; */
		map[index] = get_matPTR_from_array_of_four_subMatPTRs(content,row_level,col_level);
                if (verbose){
		  printf("ind %" PRId64 " map %p mID %" PRId64 " (%d,%d)",
                        index, map[index], get_matID_from_matPTR(map[index]), 
			 row_level, col_level);
		  printf(" [%p %p %p %p]\n", 
                        content[0], content[1], content[2], content[3]);
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

  mat_ptr_t ret = map[j_lookup_num64(j, "matid")];
  int64_t matid = get_matID_from_matPTR(ret);
  printf("\nRead %d,%d-level matrix from json file %s\n", 
	 matrix_row_level(ret), matrix_col_level(ret), path);
  printf ("  stored matrix now has address %p", ret);
  printf ("  and matrixID %" PRId64 "\n", matid);
  free(map);

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
	
	int fail =  info_set(info_type, matid, info_data);
	if (fail) {fprintf(stderr,"panic in %s\n",__func__); exit(1);}
	if (verbose) {
	  printf("Added info field %s to info_store\n", info_name);
	  printf("  with value %s\n", info_data);
	}
      }
  }  // end of if j_key_exists for "info"

  fclose(f);
  j_set_null(j); free(j);

  return ret;
}

// This routine is called by canonical_format_json.py when the python routine
// is used to convert old-style json compressed matrices (with standard C data
// types) into new-style json compressed matrices (with scalars expressed as
// strings)
int64_t read_larcMatrix_file_legacy_return_matID(char *path)
{
  
  mat_ptr_t m_ptr = read_larcMatrix_file_legacy_return_matPTR(path);
  int64_t m_mID =  get_matID_from_matPTR(m_ptr);
  return m_mID;
  
}

// At this point, this routine should only be called by
// read_larcMatrix_file_legacy_return_matID(), since it assumes the old style
// json compressed matrices and cannot handle multiprecision scalarTypes. 
mat_ptr_t read_larcMatrix_file_legacy_return_matPTR(char *path)
{

#if defined(USE_MPINTEGER) || defined(USE_MPRATIONAL) || defined(USE_MPRATCOMPLEX) || defined(USE_MPREAL) || defined(USE_MPCOMPLEX)
  // exit program if LARC is compiled for multiprecision types
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
  i = SIZE_MAX/(max_matrixID + 1)/sizeof(mat_ptr_t);
  if (i==0) {
      fprintf(stderr,"Error is %s: matrixID_max too large - try\n", __func__);
      fprintf(stderr,"renumbering the matrix_IDs in %s\n", path);
      fprintf(stderr,"to reduce this value to something more reasonable.\n");
      exit(1);
  }

  mat_ptr_t *map = calloc(max_matrixID+1, sizeof(mat_ptr_t));
  if (!map){
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
	    case 4:			/* expect a scalar link */
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
		map[index] = get_valMatPTR_from_val(val);
                if (verbose) {
		  printf("%"PRId64" %p %" PRId64 " %d %d %Lg+%Lgi\n", index, map[index], 
			 get_matID_from_matPTR((map[index])),
			 row_level, col_level, creall(val), cimagl(val));
		}
                sca_clear(&val);
	      }
	      break;
#endif
#ifdef USE_REAL
	    case 3:			/* expect a scalar link */
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
		map[index] = get_valMatPTR_from_val(mat_val);
                if (verbose) {
		  // printf("%"PRId64" %zd %" PRId64 " %d %d %Lg\n", index, map[index], 
		  printf("%"PRId64" %p %" PRId64 " %d %d %Lg\n", index, map[index], 
			 get_matID_from_matPTR((map[index])),
			 row_level, col_level, mat_val);
		}
                sca_clear(&mat_val);
	      }
              break;
#endif
#ifdef USE_INTEGER
	    case 3:			/* expect a scalar link */
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
		map[index] = get_valMatPTR_from_val(mat_val);
                if (verbose) {
		  printf("%"PRId64" %p %" PRId64 " %d %d %ld\n", index, map[index], 
			 get_matID_from_matPTR(map[index]),
			 row_level, col_level, mat_val);
		}
                sca_clear(&mat_val);
	      }
              break;
#endif
	      
	    case 6:			/* expect a 4-tuple */
	      {
                if (verbose) {printf("In non scalar case \n");}
		int row_level = j_get_num64(j_array_index(p, 0));
		int col_level = j_get_num64(j_array_index(p, 1));
                if (verbose) {printf("  levels are %d %d \n",row_level,col_level);}
		mat_ptr_t content[4];
                int64_t temp;
                for (int k=0; k< 4; ++k) {
		   temp = j_get_num64(j_array_index(p, k+2));
                   if (verbose) {printf("temp for k = %d is %ld \n",k,temp);}
                   // The matrix ID -1 is used whenever we have no matrix, 
                   // e.g. in row and col place keepers
                   if (temp == -1) { content[k] = MATRIX_PTR_INVALID;}
                   else {
                     content[k] = map[temp];
		   }
		}
		/* content[0] = map[j_get_num64(j_array_index(p, 2))]; */
		/* content[1] = map[j_get_num64(j_array_index(p, 3))]; */
		/* content[2] = map[j_get_num64(j_array_index(p, 4))]; */
		/* content[3] = map[j_get_num64(j_array_index(p, 5))]; */
		map[index] = get_matPTR_from_array_of_four_subMatPTRs(content,row_level,col_level);
                if (verbose){
		  printf("ind %" PRId64 " map %p mID %" PRId64 " (%d,%d)",
                        index, map[index], get_matID_from_matPTR(map[index]), 
			 row_level, col_level);
		  printf(" [%p %p %p %p]\n", 
                        content[0], content[1], content[2], content[3]);
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

  mat_ptr_t ret = map[j_lookup_num64(j, "matid")];
  int64_t matid = get_matID_from_matPTR(ret);
  printf("\nRead %d,%d-level matrix from json file %s\n", 
	 matrix_row_level(ret), matrix_col_level(ret), path);
  printf ("  stored matrix now has address %p", ret);
  printf ("  and matrixID %" PRId64 "\n", matid);
  free(map);

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
	
	int fail =  info_set(info_type, matid, info_data);
	if (fail) {fprintf(stderr,"panic in %s\n",__func__); exit(1);}
	if (verbose) {
	  printf("Added info field %s to info_store\n", info_name);
	  printf("  with value %s\n", info_data);
	}
      }
  }  // end of if j_key_exists for "info"

  fclose(f);
  j_set_null(j); free(j);

  return ret;
}

// SPECIFIC EXAMPLES OF ALTERING A MATRIX BEFORE STORING IT
// These functions are left over from before the read_and_alter_scalars
// functionality was implemented. 
//
// They are unittested, so simply removing them is a bad idea.

/*
This function reads in the matrix stored in a json file format
at location path, and modifies the small scalars in the matrices
to zeros, then stores the new matrix into the matrix store
and returns the matrixID of the resulting collapsed matrix.
The small scalars are any scalars that are less than or equal
to the value threshold.
TODO:
Keep in the info file any information about which matrix this
came from and another record with how much it was zeroized.
*/
mat_ptr_t matrix_read_larcMatrix_file_zeroize_small_scalars(char *path, scalarType threshold)
{
    mat_ptr_t mat_ptr = matrix_read_larcMatrix_file_flatten_small_scalars(path, threshold, scalar0);
    
    return mat_ptr;
}

/* Python interface to return a matrixID after reading a json file */
int64_t matrix_read_larcMatrix_file_zeroize_small_scalars_matrixID(char *path, char *threshold)
{
    scalarType thresh;
    sca_init(&thresh);
    sca_set_str(&thresh, threshold);

    mat_ptr_t m_ptr = matrix_read_larcMatrix_file_zeroize_small_scalars(path, thresh);
    sca_clear(&thresh);

    return get_matID_from_matPTR(m_ptr);
}

typedef struct flattening_vars_s {
    int initialized;
    scalarType threshold;
    scalarType flatValue;
} flattening_vars_t;

flattening_vars_t flatVars = {0};

// not for general use
void flatten_scalar(scalarType *flat, const scalarType flatten_this)
{
    sca_norm(flat, flatten_this);
    // set flat to flatValue if 
    //      |flatten_this| < threshold
    // and flatten_this != 0. 
    if ((sca_cmp(*flat, flatVars.threshold) < 0)
           && (sca_eq(flatten_this, scalar0) == 0))
        sca_set(flat, flatVars.flatValue);
    else
        sca_set(flat, flatten_this);
}

/*
This function reads in the matrix stored in a json file format
at location path, and modifies the small scalars in the matrices
to zeros, then stores the new matrix into the matrix store
and returns the matrixID of the resulting collapsed matrix.
The small scalars are any scalars that are less than or equal
to the value threshold.
TODO:
Keep in the info file any information about which matrix this
came from and another record with how much it was flattend.
*/
mat_ptr_t matrix_read_larcMatrix_file_flatten_small_scalars(char *path, scalarType threshold, scalarType flatValue)
{
    if (!flatVars.initialized){
        sca_init(&flatVars.threshold);
        sca_init(&flatVars.flatValue);
        flatVars.initialized = 1;
    }
    sca_set(&flatVars.threshold, threshold);
    sca_set(&flatVars.flatValue, flatValue);
    // flatten_scalar routine checks that values are nonzero before flattening.
    return read_and_alter_vals_larcMatrix_file_return_matPTR(path, flatten_scalar);
}


/* Python interface to return a matrixID after reading a json file */
int64_t matrix_read_larcMatrix_file_flatten_small_scalars_matrixID(char *path, char *threshold, char *flatValue)
{
  scalarType thresh, flatVal;
  sca_init(&thresh);
  sca_init(&flatVal);
  sca_set_str(&thresh, threshold);
  sca_set_str(&flatVal, flatValue);

  mat_ptr_t m_ptr = matrix_read_larcMatrix_file_flatten_small_scalars(path, thresh, flatVal);
  sca_clear(&thresh);
  sca_clear(&flatVal);

  return get_matID_from_matPTR(m_ptr);
}




/********************************************************
 * routine: read_matrixMarketExchange_file_return_matPTR()
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
mat_ptr_t read_matrixMarketExchange_file_return_matPTR(char * file_path)
{

#ifdef USE_INTEGER
  
#ifdef DEBUG
  printf("--> %s : %s: %d\n", __FILE__, __func__, __LINE__);
#endif  
  int verbose = 1;
  int error_code = 0;
  int errno = 0; 

  if (verbose) printf("In %s\n", __func__);

  
  /************************************************************************
  * Open the Matrix Market matrix file.
  ************************************************************************/
  FILE *fp;
  if((fp = fopen(file_path, "r+")) == NULL) {
    fprintf(stderr,"No such file: %s\n",file_path);
    exit(1);   
  }
  if (verbose) printf("Opened the data file %s for reading:\n",file_path);

  
  /************************************************************************
  * Reading the Matrix Market coordinate file: 
  *   HEADER:
  *      For example: %% MatrixMarket matrix coordinate integer general
  *      * Read first line starting with "%%" and grab the space-separated strings.
  *      * Verify that contents of this line have first three strings:
  *        "MatrixMarket", "matrix", and  "coordinate".
  *      * The fourth string should be:
  *        "integer" "real" or "complex"
  *        Assign this string to MM_scalarType.
  *      * The fifth string should be "general" "symmetric", "skew-symmetric" or "Hermetian"
  *        Assign this string to MM_symmetry
  ************************************************************************/
  // ADD CODE HERE and
  // REPLACE THIS HACK:
  char* MM_scalarType = "integer";  // TODO: handle real and complex eventually
//  char* MM_symmetry = "general";    // TODO: handle other symmetry types for which MM
//                                    //       is only giving us the lower triangular matrix.

  /************************************************************************
  * Reading the Matrix Market coordinate file: 
  *   COMMENTS AND BLANK LINES:
  *       * ignore comment lines starting with "% "
  *       * ignore blank lines
  ************************************************************************/
  // ADD CODE HERE

  
  /************************************************************************
  * Reading the Matrix Market coordinate file:
  *   MATRIX STATISTICS:
  *       * the first line of the MMfile that does not start with a "%" contains
  *         information about the size of the MMmatrix and number of nonzero entries.
  *         MM_num_rows  MM_num_cols num_nonzeros
  ************************************************************************/
  unsigned MM_num_rows;
  unsigned MM_num_cols;
  unsigned num_nonzeros;
  
  // read in the first line from the MatrixMarket file which doesn't start with "%"
  int ret = fscanf(fp, "%u %u %u", &MM_num_rows, &MM_num_cols, &num_nonzeros);
  if (ret == 3) {
    if (verbose) printf("\nRead in MM_num_rows, MM_num_cols, num_nonzeros: %u %u %u\n",
			MM_num_rows, MM_num_cols, num_nonzeros);
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
#ifdef DEBUG
    printf("<-- %s : fail: %d\n", __FILE__, __LINE__);
#endif  
    return (MATRIX_PTR_INVALID); 
  } 

  // if the file is open but nothing is left to read print EOF
  if(feof(fp)){
    fputs("EOF",stderr);
    fprintf(stderr,"ERROR: File stopped earlier than expected.\n");
    fclose(fp);
#ifdef DEBUG
    printf("<-- %s : fail: %d\n", __FILE__, __LINE__);
#endif  
    return (MATRIX_PTR_INVALID); 
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

  if (verbose)   printf("  levels are %d %d\n",(int)row_level, (int)col_level);

  // EVENTUALLY HANDLE NONSQUARE MATRICES
  // Temporay hack: is make the LARC matrix square
  mat_level_t max_level = MAX(row_level,col_level);
  if (row_level != col_level) {
    row_level = col_level = max_level;
  }

  /************************************************************************
  * Reading the Matrix Market coordinate file:
  *    NONZERO MATRIX ENTRIES:
  *       * Each nonzero matrix entry will have a line in the MMfile with
  *         row and col coordinates followed by the value of scalar in that position:
  *         I J A(I,J)
  *         If the MM_scalarType is "complex" then A(I,J) has two entries REAL and IMAG.
  ************************************************************************/
  unsigned row[num_nonzeros], col[num_nonzeros];
  int64_t val_integer[num_nonzeros];
  long double val_real[num_nonzeros];
  long double val_imag[num_nonzeros];
  scalarType val[num_nonzeros];
  mat_level_t leaf_level[num_nonzeros];
  
  // make three arrays of length num_nonzeros
  // Now we will handle the entries of the matrix
  int i;
  if (!strcmp(MM_scalarType,"complex")) {
    for (i=0;i<num_nonzeros;++i) {
      ret = fscanf(fp, "%u %u %Lg %Lg", &row[i], &col[i],
		                      &val_real[i], &val_imag[i]);
      sca_init(&val[i]);
      sca_set_2ldoubles(&val[i],val_real[i],val_imag[i]);
    }
  }
  else if (!strcmp(MM_scalarType,"real")) {
    for (i=0;i<num_nonzeros;++i) {
      ret = fscanf(fp, "%u %u %Lg", &row[i], &col[i], &val_real[i]);
      sca_init(&val[i]);
      sca_set(&val[i],val_real[i]);
    }
  }
  else if (!strcmp(MM_scalarType,"integer")) {
    for (i=0;i<num_nonzeros;++i) {
      ret = fscanf(fp, "%u %u %ld", &row[i], &col[i], &val_integer[i]);
      sca_init(&val[i]);
      sca_set(&val[i],val_integer[i]);
    }
  }
  else {
    fprintf(stderr,"FAIL in %s reading arrays: Unknown MM_scalarType\n",__func__);
    exit(0);
  }

  // change all indices to be zero based instead
  for (i=0;i<num_nonzeros;++i) {
    --row[i];
    --col[i];
  }

  // ADD CODE below to actually handle the case when MM_scalarType is complex or real.
  // Our hack is to only handle integer for now.


  

  // TRIVIAL CASE: If there is only one nonzero element in the MatrixMarket matrix
  //               then we return the matrix ptr of the LARC matrix with a single
  //               entry given by this nonzero element.
  leaf_level[0] = max_level;
  if (num_nonzeros == 1) {
    return get_matPTR_single_nonzero_using_valPTR_at_coords(leaf_level[0],row[0],col[0],&val[0]);  
  }


  /************************************************************************
   * We will be creating a structure which looks like a quadtree and 
   * contains the information needed to create the LARC matrix quadtree.
   * We will use the indices 0 to num_nonzeros-1 for matrices containing
   * a single nonzero entry.  We will call these "leaf nodes" since they
   * will be leafs in our quadtree structure.
   * We will use indices from num_nonzeros to num_nonzeros + internal_nodes-1
   * for nodes that will represent matrices with more than one nonzero entry.  
   * We call the nodes with more than one nonzero entry "internal nodes"
   * and the internal nodes each have four children, their quads,  which will
   * either contain a -1 if they correspond to a zero quadrant submatrix
   * or contain a leaf node index or contain another internal node index.
   ************************************************************************/


  /************************************************************************
   * For the worst case, the number of internal nodes is smaller than: 
   *     (num_nonzeros) * (the base_2 log of the dimension of the matrix)
   * We will be lazy and just create arrays big enough to hold the worst case stuff.
   * And we add an extra an num_nonzeros to our array so that we can be even more
   * lazy and use distinct numbers for all the leaf nodes and internal nodes.
   * Each internal node has 4 children, quad[0], quad[1], quad[2], and quad[3]
   ************************************************************************/
  int64_t max_internal_nodes = num_nonzeros * (max_level+1);
  int64_t first_internal_node_index = num_nonzeros;   // first index that is beyond
                                                      // the leaf_node indices
  int64_t internal_node_index;
  int64_t quad[max_internal_nodes][4];
  int quad_index;

  
  /************************************************************************
   * We initialize the value of each quad to be -1 
   * indicating they correspond to zero matrices.
   ************************************************************************/
  for (internal_node_index = first_internal_node_index ;
       internal_node_index < max_internal_nodes ;
       internal_node_index++)
    {
      for (quad_index = 0 ; quad_index <4 ;quad_index++) {
	quad[internal_node_index][quad_index] = -1;
      }
    }


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
  int64_t num_internal_nodes = 1;   // our first internal node is the root
  int64_t top_internal_node_index = first_internal_node_index;

  mat_level_t internal_node_level[max_internal_nodes];
  internal_node_level[top_internal_node_index] = max_level;

  
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
  * a leaf node for each nonzero entry and pushing the leaf node into our structure
  * by starting at the top_internal node and branching downward until
  * we either can substitute into a position with a -1 (corresponding to an
  * all zero submatrix) or we have a collision another leaf node, in which
  * case we create
  * a new internal node at that location and push the existing object
  * under the new internal node, then try to push or current
  * we have three possible situations when inserting a new leaf L
  * 1. current location holds -1,
  *    ACTION replace -1 with our leaf node, and go to next node to insert
  * 2. current location holds an internal node (number >= first_internal_node_index)
  *    ACTION descend to proper quad and re-evaluate situation
  * 3. current location holds an leaf node K (number < first_internal_node_index)
  *    ACTION:  a. create new internal node B and replace K with B
  *             b. K becomes a quad in B (since B is empty this is simple)
  *             c. starting from B, try to place L
  ************************************************************************/

  // initial values for the loop
  mat_level_t current_level = max_level;

  // Find the quad in which to place the leaf_node 0 
  int64_t leaf_node_index = 0;
  quad_index = get_quad_in_submatrix(row[leaf_node_index],col[leaf_node_index],current_level);
  quad[internal_node_index][quad_index] = leaf_node_index;  
  
  // loop to place the leaf_nodes into our quadtree like structure
  // and create any internal nodes necessary to tack the leafs together
  for (leaf_node_index=1;leaf_node_index<first_internal_node_index;++leaf_node_index) {

    // first attempt at insertion leaf node is in relation to internal node at top node of tree
    int64_t node_index = top_internal_node_index;
    current_level = max_level;  // the level of the internal node at top of tree

    // while loop exits  when leaf node has been successfully inserted
    while(1) {
      // find the quad index for the current leaf. This depends on the
      // current matrix level, which will decrease as loop repeats
      quad_index = get_quad_in_submatrix(row[leaf_node_index],col[leaf_node_index],current_level);

      // case 1: The quad we would like to put the leaf node in corresponds to a zero matrix
      //         (which is indicated when the value in the quad is -1).
      //         We can replace the zero matrix with the leaf node.
      if (quad[node_index][quad_index] == -1) {
        // zero matrix is currently in desired location, so we can insert leaf node
	// and continue the leaf_node_index loop
        quad[node_index][quad_index] = leaf_node_index;
	// this leaf's level is one level below internal node level
	leaf_level[leaf_node_index] = current_level-1;   
        break;
      }

      // case 2: The quad we would like to put the leaf node in corresponds to an internal node
      //         (indicated when the value in the quad is a node number that 
      //          greater than or equal to the first_internal_node_index).
      //         Since an internal node was in desired position to insert 
      //         the leaf node, we need to descend further down into the tree 
      //         from the internal node and stay in the while loop to
      //         attempt insertion again.
      if (quad[node_index][quad_index] >= first_internal_node_index) {
        node_index = quad[node_index][quad_index];
        --current_level;
        continue;
      }

      // case 3: The quad value corresponds to a previously stored leaf node
      //         (indicated by the value in the quad being between for the
      //          old_leaf_node_index being between 0 and first_internal_node_index - 1).
      //         In this case we have a collision with a previously inserted leaf node
      //         so we will create a new internal node and push the previous leaf into one
      //         of this new internal nodes quads, then continue in the while loop 
      //         trying to insert our leaf node starting from the new internal node.
      int64_t old_leaf_node_index = quad[node_index][quad_index];
      // a) create new internal node, level one less than current internal node
      current_level--;
      int64_t new_internal_node_index = first_internal_node_index + num_internal_nodes;
      internal_node_level[new_internal_node_index] = current_level;
      num_internal_nodes++;
      // b) replace old leaf with new internal node
      quad[node_index][quad_index] = new_internal_node_index;
      // c) place old leaf in a quad of the new internal node
      quad_index = get_quad_in_submatrix(row[old_leaf_node_index],
                                      col[old_leaf_node_index],current_level);
      quad[new_internal_node_index][quad_index] = old_leaf_node_index;
      // c) set the nod_index to the new internal node then continue in the while node
      //    trying to insert the leaf node.
      node_index = new_internal_node_index;
      continue;
    } // while (1)
  } // for (leaf_node_index)
  // This completes making our internal quadtree like representation
  // of the matrix.  

  // make a matrix pointer to point to the LARC matrix we will build
  mat_ptr_t ret_ptr;

  
  // Now we will create the LARC matrix from our quadtree like reprentation
  // of the MatrixMarket matrix.
  
  // First we go through all the leaf nodes of our structure and load a
  // matrix into LARC that has the right level and has a single nonzero entry
  // in the appropriate position:

  // Create an array of matrix pointers for LARC matrices
  mat_ptr_t node_ptr[num_nonzeros + num_internal_nodes];

  // Matrix Pointers are created first for all the maximally sized
  // quadrant submatrices in the LARC matrix with only a single nonzero
  // value.  These are precisely those cooresponding to the leaf nodes
  // of our quadtreelike structure, which cooresponded to the nonzero
  // values we read in from the MatrixMarket file.
  // Loop through the leaf_index for all the nonzero scalars
  // passing their value, level, and row, col coordinates
  // and construct the appropriate matrix in LARC and return the ptr.
  for (int leaf_index=0;leaf_index<num_nonzeros;++leaf_index) {
    node_ptr[leaf_index] = get_matPTR_single_nonzero_using_valPTR_at_coords(
	   leaf_level[leaf_index],
	   row[leaf_index],
	   col[leaf_index],
	   &val[leaf_index]);
  }

  // Matrix pointers are created for all the internal nodes by 
  // going through the indices in reverse order to their insertion.
  // We will use the set of quad values for each internal_node
  // to create a panel of matrix pointers to build the LARC matrix
  // for that internal node.
  mat_ptr_t panel[4];
  mat_level_t level;
  int quad_value;
  int64_t num_nodes = num_nonzeros + num_internal_nodes;  // leafs + internal nodes

  for (internal_node_index = num_nodes -1;
       internal_node_index>=first_internal_node_index;
       --internal_node_index) {

    level = internal_node_level[max_internal_nodes];

    // since we assuming a sparse matrix we retrieve the pointer
    // for the zero matrix of the correct level for quads
    mat_ptr_t zero_matrix = get_zero_matrix_ptr(level-1,level-1);

    // Build a panel from matrix pointer associate with each quad value,
    for (int index=0;index<4;++index) {
      quad_value = quad[internal_node_index][index];
      if (quad_value == -1) panel[index] = zero_matrix;  
      else panel[index] = node_ptr[quad_value];  // get the pointer to nonzero submatrix
    }

    // build a matrix from the four panels
    node_ptr[internal_node_index] = get_matPTR_from_array_of_four_subMatPTRs(panel, level, level);
  }

  // release the space from the array of scalarType
  for (leaf_node_index=0; leaf_node_index<num_nonzeros; ++leaf_node_index) {
    sca_clear(&(val[leaf_node_index]));
  }

  // the pointer to the LARC matrix is the the pointer associated with the top node
  ret_ptr = node_ptr[top_internal_node_index];

  return ret_ptr;

#endif    // USE_INTEGER

  return MATRIX_PTR_INVALID;
  
  
}  // end read_row_major_matrix_from_file


///////////////////////////////////////////////////////////////////////


static int   
recursive_write_larcMatrix_file_return_larcSize_by_matPTR(mat_ptr_t m_ptr, FILE *f, char *output_flags, size_t *larcSize_PTR, void (*func)(scalarType*, const scalarType))
{
  if (matrix_is_invalid(m_ptr)) {
    fprintf(stderr,"%s: invalid matrix pointer passed\n", __func__);
    return -1;
  }

  /* output_flags keeps track of which matrixIDs have already appeared */
  if (output_flags[get_matID_from_matPTR(m_ptr)] > 0)
    return 1;

  output_flags[get_matID_from_matPTR(m_ptr)] = 1;
  *larcSize_PTR +=1;
  matrix_type_t mtype = matrix_type(m_ptr);

  if (mtype == SCALAR)
    {
      // in COMPLEX case, write routine to output "a+I*b" as "a, b"
      scalarType *val = &scratchVars.submit_to_store;
      // the function func() is applied to the scalar value in m_ptr, and
      // the result put into val; if the function is sca_set, the stored value
      // is merely copied
      func(val, matrix_trace(m_ptr));
      char *val_string = sca_get_str(*val);
      fprintf(f, "    \"%ld\":[%d, %d, \"%s\"],\n", get_matID_from_matPTR(m_ptr), matrix_row_level(m_ptr), 
              matrix_col_level(m_ptr), val_string);
      free(val_string);
    }
  else   // NONSCALAR
    {
      int ret;
     
      mat_ptr_t panel[4];
      for (int i = 0; i < 4; ++i) {
        panel[i] = matrix_sub(m_ptr,i);
      }

      // obtain the pointers for the matrices in the panel
      // For ROW_VECTORS, COL_VECTORS and MATRICES there is a panel[0]
      ret = recursive_write_larcMatrix_file_return_larcSize_by_matPTR(panel[0], f, output_flags, larcSize_PTR, func);
      if (ret < 0) return ret;

      // For ROW_VECTORS, MATRICES there is a panel[1]
      if (mtype!=COL_VECTOR) {
        ret = recursive_write_larcMatrix_file_return_larcSize_by_matPTR(panel[1], f, output_flags, larcSize_PTR, func);
        if (ret < 0) return ret;
      }

      // For COL_VECTORS, MATRICES there is a panel[2]
      if (mtype!=ROW_VECTOR) {
        ret = recursive_write_larcMatrix_file_return_larcSize_by_matPTR(panel[2], f, output_flags, larcSize_PTR, func);
        if (ret < 0) return ret;
      }

      // For MATRICES there is a panel[3]
      if (mtype==MATRIX) {
        ret = recursive_write_larcMatrix_file_return_larcSize_by_matPTR(panel[3], f, output_flags, larcSize_PTR, func);   
        if (ret < 0) return ret;
      }

      // print the panel, which contains matrixIDs for these matrices 
      // For ROW_VECTORS, COL_VECTORS and MATRICES there is a panel[0]
      fprintf(f, "    \"%ld\":[%d, %d, %ld, ", get_matID_from_matPTR(m_ptr), 
                matrix_row_level(m_ptr), matrix_col_level(m_ptr),
	        get_matID_from_matPTR(panel[0]));

      // For ROW_VECTORS, MATRICES there is a panel[1]
      // For COL_VECTORS print -1
      if (mtype!=COL_VECTOR) {
                fprintf(f, "%ld, ", get_matID_from_matPTR(panel[1]));
      }
      else { fprintf(f, "-1, "); }

      // For COL_VECTORS, MATRICES there is a panel[2]
      // For ROW_VECTORS print -1
      if (mtype!=ROW_VECTOR) {
                fprintf(f, "%ld, ", get_matID_from_matPTR(panel[2]));
      }
      else { fprintf(f, "-1, "); }

      // For MATRICES there is a panel[3]
      // For ROW_VECTORS or COL_VECTORS print -1
      if (mtype==MATRIX) {
                fprintf(f, "%ld],\n", get_matID_from_matPTR(panel[3]));
      }
      else { fprintf(f, "-1],\n"); }

      //fprintf(f, "    \"%ld\":[%d, %d, %ld, %ld, %ld, %ld],\n", get_matID_from_matPTR(m_ptr), 
      //      matrix_row_level(m_ptr),matrix_col_level(m_ptr),
      //      get_matID_from_matPTR(panel[0]),get_matID_from_matPTR(panel[1]),
      //      get_matID_from_matPTR(panel[2]),get_matID_from_matPTR(panel[3]));
    }
  return 1;
}



size_t write_larcMatrix_file_return_larcSize(mat_ptr_t m_ptr, char *path,
        void (*func)(scalarType*, const scalarType))
{
// return zero if no file is written for whatever reason
/* output_flags keeps track of which matrixIDs have already appeared */

  int verbose = 0;
  
  FILE *f = fopen(path, "w");

  size_t larcSize = 0;

  if (!f)
    {
        fprintf(stderr,"ERROR in %s, could not open %s for writing\n",__func__,path);
        return 0;
    }
  if (matrix_is_invalid(m_ptr)) 
    {
        fprintf(stderr,"ERROR in %s, passed matrix pointer is invalid\n",__func__);
        return 0;
    }
  char *output_flags = calloc(num_matrices_created(), sizeof(char));
  if (output_flags == NULL) {
    ALLOCFAIL();
  }

  if (!output_flags)  {
        fprintf(stderr,"ERROR in %s, unable to allocate space.\n",__func__);
        return 0;
  }

  if (verbose) {
      printf("\nWriting %d,%d-level matrix with matrixID %" PRId64 " to %s\n", 
	     matrix_row_level(m_ptr), matrix_col_level(m_ptr),
	     get_matID_from_matPTR(m_ptr), path);
    }
  
  // printf("      the matrix id is %lu, with level %d %d\n",id, matrix_row_level(m_ptr), matrix_col_level(m_ptr));

  // fprintf(f, "{\n  \"matrixID_max\":%" PRIu64 ",\n  \"matid\":%" PRIu64 ",\n  \"table\":{\n", num_matrices_created(), get_matID_from_matPTR(m_ptr));


  // json file contains: matriID_max, matid, optional info struct, table struct
  fprintf(f, "{\n  \"matrixID_max\":%" PRIu64 ",\n  \"matid\":%" PRIu64 ",",
	  num_matrices_created(), get_matID_from_matPTR(m_ptr));

  // the info structure is printed when it contains information

  // TODO
  //    * retrieve the larcSize from infoStore if it exists
  //    * after writing the the recursive file, we have calculated a larcSize
  //    * check to see if these are the same
  //    * enter the new larcSize into the matrix store
  //    * insert the larcSize either at the end of the larcFile or at the beginning.
  int info_ret = write_infoStore_to_larcMatrix_file(m_ptr, f);
  if (info_ret != 0) {
    fprintf(stderr,"ERROR in %s, failure in function\n",__func__);
    fprintf(stderr,"\t write_infoStore_to_larcMatrix_file\n");
  }

  int recur_ret;
  // always print the matrix table recursively
  fprintf(f, "\n  \"table\":{\n"); 
  recur_ret = recursive_write_larcMatrix_file_return_larcSize_by_matPTR(
	                           m_ptr, f, output_flags, &larcSize, func);
  if (recur_ret == -1) {
    fprintf(stderr,"ERROR in %s, failure at some level of recursive function\n",__func__);
    fprintf(stderr,"\t recursive_write_larcMatrix_file_return_larcSize_by_matPTR\n");
  }
    
  free (output_flags);

  // print out the larcSize
  if (verbose) {
    fprintf(stdout,"The larcSize is %zd for the larcMatrix in file %s.\n",larcSize,path);
  }

  // end the table structure and the json file
  fprintf(f, "      \"end\":0 }\n}");
  fclose(f);

  return larcSize;
}


size_t write_larcMatrix_file_return_larcSize_by_matID(int64_t m_mID, char *path)
{
  // get the matrix pointer from the matrixID, and see if still in store 
  mat_ptr_t m_ptr = get_matPTR_from_matID(m_mID, "", __func__,0);
  if (m_ptr == MATRIX_PTR_INVALID) { exit(1); }

  // calculate matrix pointer version of function
  
  return write_larcMatrix_file_return_larcSize_by_matPTR(m_ptr, path);
}

size_t write_larcMatrix_file_return_larcSize_by_matPTR(mat_ptr_t m_ptr, char *path)
{
    // the function sca_set as the last argument ensures the written data
    // is the same as that in the matrixStore
    return write_larcMatrix_file_return_larcSize(m_ptr, path, sca_set);
}


