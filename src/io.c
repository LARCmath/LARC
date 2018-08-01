//                            io.c 
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


#include <inttypes.h> //for printing uint64_t and int64_t
#include "io.h"


int64_t row_major_list_to_store_matrixID(ScalarType  * dense_mat, 
                    mat_level_t        current_row_level, 
                    mat_level_t        current_col_level, 
                    int64_t        orig_num_cols 
                    )
{
  mat_add_t C_ptr = row_major_list_to_store(dense_mat,current_row_level,current_col_level,orig_num_cols); 
  int64_t C_mID = get_matrixID_from_ptr(C_ptr);
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
 *    It loads the matrix into the matrix store c using matrix_get_ptr_panel
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
 *                  and submatrix[3] starts at 10,
               and
 *       the subarray of dense_mat that starts at that index.
 *
 *    The recursive nature of the call takes care of everything else.
 *    When matrix_get_ptr_panel has the mat_ptr of each of the four submatrices
 *    (in the vector cases two of these are MATRIX_PTR_INVALID)
 *    then it can be used to find the mat_ptr of a matrix.
 *    
 **********************************************************************/

mat_add_t row_major_list_to_store(ScalarType  * dense_mat, 
                    mat_level_t        current_row_level, 
                    mat_level_t        current_col_level, 
                    int64_t        orig_num_cols 
                    )
{
#ifdef DEBUG
  printf("--> %s : %s: %d\n", __FILE__, __func__, __LINE__);
#endif  

  // Note that no check to verify if dense_mat is sufficiently long for size. 

  // This has not been updated for nonsquare matrices
  //int verbose = 0;
/*   if (verbose) { */
/*     printf("\nIn subroutine row_major_list_to_store with level %ld ",level); */
/*     printf("and orig_num_cols %ld\n",orig_num_cols); */
/*     int i,j;  */
/*     for (i =0; i < (1L << level); ++i) { */
/*       for (j =0; j < (1L << level); ++j) { */
/* #ifdef USE_COMPLEX */
/*         printf("%4.1lf + %4.1lf i, \t ",creal(dense_mat[i*orig_num_cols+j]), */
/*                                cimag(dense_mat[i*orig_num_cols+j])); */
/* #else */
/*         printf("%4.1lf, \t ",dense_mat[i*orig_num_cols+j]); */
/* #endif */
/*       } */
/*       printf("\n"); */
/*     } */
/*   }  // end verbose */
    
  // SCALAR CASE
  if ((current_row_level == 0) && (current_col_level == 0)){
    return matrix_get_ptr_scalar(*dense_mat);
  }

  // Need this for remaining cases
  mat_add_t panel[4];  // indices of submatrices


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

#if 0
    if (verbose) {
      printf("\nIn recursive call portion with dim_new %ld\n",dim_new);
      printf("  The four subpanels will be at level %ld and orig_num_cols %ld\n",level-1,orig_num_cols);
      printf("  and will have offsets: \n");
      printf("     panel[0]: 0\n");
      printf("     panel[1]: dim_new = %ld\n",dim_new);
      printf("     panel[2]: dim_new*orig_num_cols = %ld\n",dim_new*orig_num_cols);
      printf("     panel[3]: dim_new*orig_num_cols+dim_new = %ld\n",dim_new*orig_num_cols+dim_new);
    }  // end verbose
#endif

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

  return matrix_get_ptr_panel(panel, current_row_level, current_col_level);
}


/****************************************************************
 *  python interface: read_row_major_matrix_from_file_matrixID
 ****************************************************************/
int64_t read_row_major_matrix_from_file_matrixID(char * file_path)
{
  mat_add_t m_ptr = read_row_major_matrix_from_file(file_path);
  int64_t m_mID = get_matrixID_from_ptr(m_ptr);
  return m_mID;
}



/****************************************************************
 *  read_row_major_matrix_from_file
 ****************************************************************/
mat_add_t read_row_major_matrix_from_file(char * file_path)
{
#ifdef DEBUG
  printf("--> %s : %s: %d\n", __FILE__, __func__, __LINE__);
#endif  
  int verbose = 0;
  int error_code = 0;
  int errno = 0; // TEMPORARY FIX

  if (verbose) printf("In %s\n", __func__);

  FILE *fp;
  if((fp = fopen(file_path, "r+")) == NULL) {
    printf("No such file: %s\n",file_path);
    exit(1);    // LOOK UP MAGIC WORD
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
    perror("scanf:");
    error_code = 1;
    printf("ERROR scanning levels in first line of file %s\n", file_path);
  }
  else if(ret == EOF) {
    error_code = 2;
    printf("ERROR: file %s ended early in %s\n", file_path,  __func__);
  }
  else {
    printf("ERROR while reading levels in %s (routine %s).\n", file_path, __func__);
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
    puts("EOF");
    printf("ERROR: File stopped earlier than expected.\n");
    fclose(fp);
#ifdef DEBUG
    printf("<-- %s : fail: %d\n", __FILE__, __LINE__);
#endif  
    return (MATRIX_PTR_INVALID); 
  }

  // successfully read the level
  if (verbose)   printf("  levels are %d %d\n",(int) row_level, (int) col_level);

  // calculate size of the array
  int64_t array_size = (1 << row_level) * (1<<col_level);
  if (verbose)  printf("  array length is %d\n",(int) array_size);

  // Now we will handle the entries of the matrix
  ScalarType * matrix_array;

  // Allocate space for the array of ScalarType
  if (!(matrix_array = calloc(array_size,sizeof(ScalarType)) )) {
    error_code = 1;  // unable to allocate space for array
    printf("ERROR: Unable to allocate space for array.\n");
#ifdef DEBUG
    printf("<-- %s : fail: %d\n", __FILE__, __LINE__);
#endif  
    return (MATRIX_PTR_INVALID); 
  } 

  // Read the array of numbers 
#ifdef USE_INTEGER 
  int64_t input_value;
#endif
#ifdef USE_COMPLEX
  double real, imag;
#endif
#ifdef USE_REAL
  double real;
#endif
  int i=0;

  for (i=0;i<array_size;++i) {
#ifdef USE_COMPLEX
    int ret = fscanf(fp, "%lf+%lfi", &real, &imag);
    if (ret == 2) {
      if (verbose) printf("\ni is %d\n",i);
      if (verbose) printf("  read in: %lf \t %lf\n", real, imag);
      fflush(stdout);
      matrix_array[i] = real + imag*I;
    } 
#endif
#ifdef USE_REAL
    int ret = fscanf(fp, "%lf", &real);
    if (ret == 1) {
      if (verbose) printf("\ni is %d\n",i);
      if (verbose) printf("  read in: %lf \n", real);
      fflush(stdout);
      matrix_array[i] = real;
    } 
#endif
#ifdef USE_INTEGER
    int ret = fscanf(fp, "%ld", &input_value);
    if (ret == 1) {
      if (verbose) printf("\ni is %d\n",i);
      if (verbose) printf("  read in: %ld \n", input_value);
      fflush(stdout);
      matrix_array[i] = input_value;
    } 
#endif
    else if(errno != 0) {
      printf("ERROR: Error while reading array entry.\n");
      perror("scanf:");
      error_code =2;
      if (verbose) printf("  Error was %d\n",error_code);
    }
    else if(ret == EOF) {
      printf("ERROR: EOF Error while reading array entry.\n");
      error_code = 3;
      if (verbose) printf("  Error was %d\n",error_code);
    }
    else {
      printf("ERROR: Data not in expected format\n");
      error_code = 4;
      if (verbose) printf("  Error was %d\n",error_code);
    }
  } // end for loop
  if (verbose) printf("\n");

  if (error_code) {
#ifdef DEBUG
  printf("<-- %s : %s fail: %d\n", __FILE__, __func__, __LINE__);
#endif  
  return (MATRIX_PTR_INVALID); 
  }

  // if the file is open but nothing is left to read print EOF
  if(feof(fp)){
    puts("EOF");
  }

  fclose(fp);

  if (array_size == i) {
#ifdef DEBUG
    printf("<-- %s : %s success: %d\n", __FILE__, __func__,__LINE__);
#endif  
    mat_add_t ret_ptr = row_major_list_to_store(matrix_array, 
                                          row_level,col_level,
						1<<col_level);
    return ret_ptr;

  }
  else {
    if (verbose) printf("Failed to read column major matrix from file.\n");
#ifdef DEBUG
    printf("<-- %s : %s fail: %d\n", __FILE__, __func__, __LINE__);
#endif  
    return MATRIX_PTR_INVALID;
  }

}  // end read_row_major_matrix_from_file



/*********************************************************************
 *                  get_i_j_th_entry_of_matrix                       *
 *                                                                   *
 *   This algorithm finds the ith jth entry of a matrix, by          *
 *   looking at the binary representation of i and j and             *
 *   recursively calling the function on the appropriate sub matrix  *
 *   depending on the (high bit of i, high bit of j):                *
 *     (0,0) -> sub_matrix[0],  (0,1) -> sub_matrix[1],              *
 *     (1,0) -> sub_matrix[1],  (1,1)-> sub_matrix[3]                *
 *                                                                   *
 *   NOTE:  This function returns a mat_add_t return_ptr             *
 *          To get the actual value of the i_jth entry  use          *
 *          matrix_trace(return_ptr)                                 *
 *********************************************************************/
static
mat_add_t get_ijth_matrix_entry(long int row_i, long int col_j, mat_add_t m_ptr) 
{
  matrix_type_t mat_type = matrix_type(m_ptr);

  // Case: SCALAR  return the matrix_ptr of that scalar
  if (mat_type == SCALAR) {return m_ptr;
  } 

  // In all other cases the calculation will proceed recursively
  // using a subpanel of the matrix and new_i and new_j
  long int half_i_dim;  // half size i dimension
  int i_adjust;    // 0 or half_i_dim depending on size of i
  int i_top_bit = 0;
  long int new_i = 0;  // mod off half dim
  long int half_j_dim;  // half size j dimension
  int j_adjust; // (0 or half_j_dim) depending on size of j
  int j_top_bit = 0;
  long int new_j = 0; // mod off half dim

  // Case: mat_type  MATRIX / COL_VECTOR, calculate new_i, and i_top_bit
  int row_level = matrix_row_level(m_ptr);
  if (row_level > 0) {
    half_i_dim = 1 << (row_level-1);  // half size i dimension
    i_adjust = half_i_dim & row_i;   // 0 or half_i_dim depending on size of i
    i_top_bit = i_adjust >> (row_level-1);
    new_i = row_i - i_adjust; // mod off half dim
  } 

  // Case: mat_type  MATRIX / ROW_VECTOR, calculate new_j, and j_top_bit
  int col_level = matrix_col_level(m_ptr);
  if  (col_level > 0) {
    half_j_dim = 1 << (col_level-1);  // half size j dimension
    j_adjust = half_j_dim & col_j;   // (0 or half_j_dim) depending on size of j
    j_top_bit = j_adjust >> (col_level-1);
    new_j = col_j - j_adjust; // mod off half dim
  }
  
  // Calculate the correct submatrix to use in recursive call
  int panel_index = i_top_bit * 2 + j_top_bit;
  return get_ijth_matrix_entry(new_i,new_j,matrix_sub(m_ptr,panel_index));

}
 
/*******************************************************************
*                       (naive) print_matrix_naive_by_matrixID       *
* Python Interface for print_matrix_naive      			*
* accepts a matrixID for matrix to print                           *
* Converts to a pointer before calling								*
********************************************************************/
void print_matrix_naive_by_matrixID(int64_t A_mID)
{  
  // get the matrix pointer from the matrixID, and see if still in store 
  mat_add_t A_ptr = mat_ptr_from_matrixID(A_mID, "first", __func__,0);
  if (A_ptr == MATRIX_PTR_INVALID) { exit(1); }

  // should have a valid matrix pointer, to call internal routine for printing
  print_matrix_naive(A_ptr);
}


/********************************************************************
*                       (naive) print_matrix  (to screen)           *
********************************************************************/
void print_matrix_naive(mat_add_t ptr_A)
{
  int i,j;
  int row_level =  matrix_row_level(ptr_A);
  int col_level =  matrix_col_level(ptr_A);
  int row_dim = 1 << row_level;
  int col_dim = 1 << col_level;
  ScalarType entry;

  // debugging
  for(i = 0;i<row_dim;i++){
    for(j = 0;j<col_dim; j++){
      entry = matrix_trace(get_ijth_matrix_entry(i, j, ptr_A));

#ifdef USE_COMPLEX
      printf("%lg+%lgi ",creal(entry),cimag(entry));
#endif
#ifdef USE_REAL
      printf("%lg ",entry);
#endif
#ifdef USE_INTEGER
      printf("%ld ",entry);
#endif
      }
   printf(" \n");
   }
}



/********************************************************************
* Python interface version of (naive) print_matrix_to_file          *
*   this could be implemented less naively ...                      *
********************************************************************/
void print_matrix_to_file_naive_by_matrixID(int64_t A_mID, char * filename)
{
  // get the matrix pointer from the matrixID, and see if still in store 
  mat_add_t A_ptr = mat_ptr_from_matrixID(A_mID, "first", __func__,0);
  if (A_ptr == MATRIX_PTR_INVALID) { exit(1); }

  // should have a valid matrix pointer, to call internal routine for printing matrix to a file
  print_matrix_to_file_naive(A_ptr, filename);

}



/********************************************************************
*                       (naive) print_matrix_to_file                *
*   this could be implemented less naively ...                      *
********************************************************************/
void print_matrix_to_file_naive(mat_add_t ptr_A, char * filename)
{
  FILE *f = fopen(filename, "w");
  int i,j;
  int row_level =  matrix_row_level(ptr_A);
  int col_level =  matrix_col_level(ptr_A);
  int row_dim = 1 << row_level;
  int col_dim = 1 << col_level;
  ScalarType entry;

  for(i = 0;i<row_dim;i++){
    for(j = 0;j<col_dim; j++){
      entry = matrix_trace(get_ijth_matrix_entry(i, j, ptr_A));
#ifdef USE_COMPLEX
      if (cimag(entry) == 0)
	fprintf(f,"%lg ",creal(entry));
      else
	fprintf(f,"%lg+%lgi ",creal(entry),cimag(entry));
#endif
#ifdef USE_REAL
      fprintf(f,"%lg ",entry);
#endif
#ifdef USE_INTEGER
      fprintf(f,"%ld ",entry);
#endif
    }
    fprintf(f," \n");
  }
  fclose(f);
}



/********************************************************************
* (Python interface) print_matrix_nonzeros_to_file_by_matrixID            *
********************************************************************/
void print_matrix_nonzeros_to_file_by_matrixID(int64_t A_mID, char * filename)
{
  // get the matrix pointer from the matrixID, and see if still in store 
  mat_add_t A_ptr = mat_ptr_from_matrixID(A_mID, "first", __func__,0);
  if (A_ptr == MATRIX_PTR_INVALID) { exit(1); }

  print_matrix_nonzeros_to_file(A_ptr, filename);

}

 


/********************************************************************
*                print_matrix_nonzeros_to_file               *
********************************************************************/
void print_matrix_nonzeros_to_file(mat_add_t ptr_A, char * filename)
{
  FILE *f = fopen(filename, "w");
  int i,j;
  int row_level =  matrix_row_level(ptr_A);
  int col_level =  matrix_col_level(ptr_A);
  int row_dim = 1 << row_level;
  int col_dim = 1 << col_level;
  ScalarType entry;

  for(i = 0;i<row_dim;i++){
    for(j = 0;j<col_dim; j++){
      entry = matrix_trace(get_ijth_matrix_entry(i, j, ptr_A));
      if (entry != 0) {
#ifdef USE_COMPLEX
        if (cimag(entry) == 0)
     	  fprintf(f,"%lg\t(%d,%d)\n",creal(entry),i,j);
        else
	  fprintf(f,"%lg+%lgi\t(%d,%d)\n",creal(entry),cimag(entry),i,j);
#endif
#ifdef USE_REAL
     	fprintf(f,"%lg\t(%d,%d)\n",entry,i,j);
#endif
#ifdef USE_INTEGER
     	fprintf(f,"%ld\t(%d,%d)\n",entry,i,j);
#endif
      }
    }
  }
  fclose(f);
}


static int
matrix_write_json_file_1(mat_add_t m_ptr, FILE *f, char *output_flags)
{
  if (matrix_is_invalid(m_ptr)) {
    printf("%s: invalid matrix pointer passed\n", __func__);
    return -1;
  }

  /* output_flags keeps track of which matrixIDs have already appeared */
  if (output_flags[get_matrixID_from_ptr(m_ptr)] > 0)
    return 1;

  output_flags[get_matrixID_from_ptr(m_ptr)] = 1;
  matrix_type_t mtype = matrix_type(m_ptr);

  if (mtype == SCALAR)
    {
#ifdef USE_COMPLEX
      fprintf(f, "    \"%ld\":[%d, %d, %.25g, %.25g],\n", get_matrixID_from_ptr(m_ptr), matrix_row_level(m_ptr),
	      matrix_col_level(m_ptr),creal(matrix_trace(m_ptr)), cimag(matrix_trace(m_ptr)));
#endif
#ifdef USE_REAL
      fprintf(f, "    \"%ld\":[%d, %d, %.25g],\n", get_matrixID_from_ptr(m_ptr), matrix_row_level(m_ptr), 
              matrix_col_level(m_ptr),matrix_trace(m_ptr));
#endif
#ifdef USE_INTEGER
      fprintf(f, "    \"%ld\":[%d, %d, %ld],\n", get_matrixID_from_ptr(m_ptr), matrix_row_level(m_ptr), 
              matrix_col_level(m_ptr),matrix_trace(m_ptr));
#endif
    }
  else   // NONSCALAR
    {
      int ret;
     
      mat_add_t panel[4];
      for (int i = 0; i < 4; ++i) {
        panel[i] = matrix_sub(m_ptr,i);
      }
      ret = matrix_write_json_file_1(panel[0], f, output_flags);
      if (ret < 0) return ret;
      if (mtype!=COL_VECTOR) {
        ret = matrix_write_json_file_1(panel[1], f, output_flags);
        if (ret < 0) return ret;
      }
      if (mtype!=ROW_VECTOR) {
        ret = matrix_write_json_file_1(panel[2], f, output_flags);
        if (ret < 0) return ret;
      }
      if (mtype==MATRIX) {
        ret = matrix_write_json_file_1(panel[3], f, output_flags);
        if (ret < 0) return ret;
      }
      fprintf(f, "    \"%ld\":[%d, %d, %ld, ", get_matrixID_from_ptr(m_ptr), 
                matrix_row_level(m_ptr), matrix_col_level(m_ptr),
	        get_matrixID_from_ptr(panel[0]));
      if (mtype!=COL_VECTOR) {
                fprintf(f, "%ld, ", get_matrixID_from_ptr(panel[1]));
      }
      else { fprintf(f, "-1, "); }
      if (mtype!=ROW_VECTOR) {
                fprintf(f, "%ld, ", get_matrixID_from_ptr(panel[2]));
      }
      else { fprintf(f, "-1, "); }
      if (mtype==MATRIX) {
                fprintf(f, "%ld],\n", get_matrixID_from_ptr(panel[3]));
      }
      else { fprintf(f, "-1],\n"); }
      //fprintf(f, "    \"%ld\":[%d, %d, %ld, %ld, %ld, %ld],\n", get_matrixID_from_ptr(m_ptr), 
      //      matrix_row_level(m_ptr),matrix_col_level(m_ptr),
      //      get_matrixID_from_ptr(panel[0]),get_matrixID_from_ptr(panel[1]),
      //      get_matrixID_from_ptr(panel[2]),get_matrixID_from_ptr(panel[3]));
    }
  return 1;
}


int
matrix_write_json_file_matrixID(int64_t m_mID, char *path)
{
  // get the matrix pointer from the matrixID, and see if still in store 
  mat_add_t m_ptr = mat_ptr_from_matrixID(m_mID, "first", __func__,0);
  if (m_ptr == MATRIX_PTR_INVALID) { exit(1); }

  // calculate matrix pointer version of function
  
  return matrix_write_json_file(m_ptr, path);
}


/* This function should retrieve all information out of the info_store
   for the given matrix associated with m_ptr and write it to the
   json file associated with the output file pointer f.
*/
static int matrix_write_json_file_info(mat_add_t m_ptr, FILE *f) {
  int verbose = 0;

  if (!f)
    return -1;
  if (matrix_is_invalid(m_ptr)) 
    return -1;

  int64_t m_ID = get_matrixID_from_ptr(m_ptr);

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
      free(info_data);
    }
  }
  if (was_info) {
    // fprintf(f, "      \"end\":0 },");
    fprintf(f, "      \"end\":\"\" },");
  }
  return(0);
}


int
matrix_write_json_file(mat_add_t m_ptr, char *path)
{
  /* output_flags keeps track of which matrixIDs have already appeared */
  char *output_flags = calloc(num_matrices_created(), sizeof(char));
  if (output_flags == NULL) {
    ALLOCFAIL();
  }
  FILE *f = fopen(path, "w");
  int ret;

  if (!output_flags)
    return -1;
  if (!f)
    return -1;
  if (matrix_is_invalid(m_ptr)) 
    return -1;

  printf("\nWriting %d,%d-level matrix with matrixID %" PRId64 " to %s\n", 
	 matrix_row_level(m_ptr), matrix_col_level(m_ptr),
         get_matrixID_from_ptr(m_ptr), path);
  // printf("      the matrix id is %lu, with level %d %d\n",id, matrix_row_level(m_ptr), matrix_col_level(m_ptr));

  // fprintf(f, "{\n  \"matrixID_max\":%" PRIu64 ",\n  \"matid\":%" PRIu64 ",\n  \"table\":{\n", num_matrices_created(), get_matrixID_from_ptr(m_ptr));


  // json file contains: matriID_max, matid, optional info struct, table struct
  fprintf(f, "{\n  \"matrixID_max\":%" PRIu64 ",\n  \"matid\":%" PRIu64 ",",
	  num_matrices_created(), get_matrixID_from_ptr(m_ptr));


  /* // TODO: take this fake out, write the routine matrix_write_json_file_info *\/ */
  /* //       and call it instead *\/ */
  /* fprintf(f, "\n  \"info\":{\n");   */
  /* fprintf(f, "    \"NUMROUNDS\":\"8\",\n");  */
  /* fprintf(f, "      \"end\":\"\" },");  */

  // the info structure is printed when it contains information
  ret = matrix_write_json_file_info(m_ptr, f);
  if (ret != 0) {
    printf("WARNING in %s, function matrix_write_json_file_info  failed\n",__func__);
  }

  // always print the matrix table recursively
  fprintf(f, "\n  \"table\":{\n"); 
  ret = matrix_write_json_file_1(m_ptr, f, output_flags);
  free (output_flags);

  // end the table structure and the json file
  fprintf(f, "      \"end\":0 }\n}");
  fclose(f);

  return ret;
}


/* Python interface to return a matrixID after reading a json file */
int64_t
matrix_read_json_file_matrixID(char *path)
{
  
  mat_add_t m_ptr = matrix_read_json_file(path);
  int64_t m_mID =  get_matrixID_from_ptr(m_ptr);
  return m_mID;
  
}


/* Python interface: Reads a two json formatted compressed matrices
   from files, returns their matrixID if they are equal 
   matrices, and 0 otherwise */
int64_t
equal_matrices_in_json_files(char *path1, char *path2)
{
  if (!strcmp(path1,path2)) /* paths are equal */
    {
      printf("Warning %s : %s fail: line %d, path1 and path2 are the same.\n",
              __FILE__, __func__, __LINE__);
    }
  mat_add_t m1_ptr = matrix_read_json_file(path1);
  mat_add_t m2_ptr = matrix_read_json_file(path2);
  if ( matrix_is_invalid(m1_ptr) || matrix_is_invalid(m2_ptr) ) {
      printf("Warning %s : %s fail: line %d, invalid matrix\n",
              __FILE__, __func__, __LINE__);
      return(0);
  }
  int64_t m1_mID =  get_matrixID_from_ptr(m1_ptr);
  int64_t m2_mID =  get_matrixID_from_ptr(m2_ptr);
  if (m1_mID == m2_mID) {
    return m1_mID;
  }
  else return 0;
}





mat_add_t
matrix_read_json_file(char *path)
{
  int verbose = 0;
  FILE *f = fopen(path, "r");
  // printf("  path = %s\n",path);
  if (f == NULL)
  {
    printf("%s: no file found at\n\t%s:\n\texiting\n", __func__, path);
    exit(1);
  }
  json_t *j = j_parse_file(f);
  // printf("  parsed json\n");

  int64_t max_matrixID = j_lookup_num64(j, "matrixID_max");
  // printf("  max matrixID is %ld\n", max_matrixID);

  // loop index for reading info and table
  int64_t i;

  json_t *t = j_key_lookup(j, "table");

  mat_add_t *map = calloc(max_matrixID+1, sizeof(mat_add_t));
  int64_t len = j_key_count(t);

  if (verbose) {printf("In %s, before the table line reader\n", __func__);}
  
  for (i=0; i<len; i++)
    {
      json_t *p = j_key_index(t, i);
      int64_t index = atoi(p->name);
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
                // ScalarType val = complex_val[0]+I*complex_val[1];
                ScalarType val =  j_get_double(j_array_index(p,2)) +
		        I*j_get_double(j_array_index(p,3));
                if (row_level || col_level) {
                        printf("error in %s:\n\texpected complex ", __func__ );
                        printf("scalar, but levels greater than zero!\n");
                }
		map[index] = matrix_get_ptr_scalar(val);
                if (verbose) {
		  printf("%"PRId64" %p %" PRId64 " %d %d %g+%gi\n", index, map[index], 
			 get_matrixID_from_ptr((map[index])),
			 row_level, col_level, creal(val), cimag(val));
		}
	      }
	      break;
#endif
#ifdef USE_REAL
	    case 3:			/* expect a scalar link */
	      {
                if (verbose) {printf("In real scalar case \n");}
		int row_level = j_get_num64(j_array_index(p, 0));
		int col_level = j_get_num64(j_array_index(p, 1));
		ScalarType mat_val = j_get_double(j_array_index(p,2));
                if (row_level || col_level) {
                        printf("error in %s:\n\texpected real ", __func__ );
                        printf("scalar, but levels greater than zero!\n");
                }
		map[index] = matrix_get_ptr_scalar(mat_val);
                if (verbose) {
		  // printf("%"PRId64" %zd %" PRId64 " %d %d %g\n", index, map[index], 
		  printf("%"PRId64" %p %" PRId64 " %d %d %g\n", index, map[index], 
			 get_matrixID_from_ptr((map[index])),
			 row_level, col_level, mat_val);
		}
	      }
              break;
#endif
#ifdef USE_INTEGER
	    case 3:			/* expect a scalar link */
	      {
                if (verbose) {printf("In integer scalar case \n");}
		int row_level = j_get_num64(j_array_index(p, 0));
		int col_level = j_get_num64(j_array_index(p, 1));
		ScalarType mat_val = j_get_num64(j_array_index(p,2));
                if (row_level || col_level) {
                        printf("error in %s:\n\texpected integer ", __func__ );
                        printf("scalar, but levels greater than zero!\n");
                }
		map[index] = matrix_get_ptr_scalar(mat_val);
                if (verbose) {
		  printf("%"PRId64" %p %" PRId64 " %d %d %ld\n", index, map[index], 
			 get_matrixID_from_ptr(map[index]),
			 row_level, col_level, mat_val);
		}
	      }
              break;
#endif
	      
	    case 6:			/* expect a 4-tuple */
	      {
                if (verbose) {printf("In non scalar case \n");}
		int row_level = j_get_num64(j_array_index(p, 0));
		int col_level = j_get_num64(j_array_index(p, 1));
                if (verbose) {printf("  levels are %d %d \n",row_level,col_level);}
		mat_add_t content[4];
                int temp;
                for (int k=0; k< 4; ++k) {
		   temp = j_get_num64(j_array_index(p, k+2));
                   if (verbose) {printf("temp for k = %d is %d \n",k,temp);}
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
		map[index] = matrix_get_ptr_panel(content,row_level,col_level);
                if (verbose){
		  printf("ind %" PRId64 " map %p mID %" PRId64 " (%d,%d)",
                        index, map[index], get_matrixID_from_ptr(map[index]), 
			 row_level, col_level);
		  printf(" [%p %p %p %p]\n", 
                        content[0], content[1], content[2], content[3]);
		}
	      }
	      break; 
	      
	    default:
              printf("Error in %s(%s): unexpected number of entries per line\n",
                __func__, path);
#ifdef USE_COMPLEX
              printf("In complex mode expected num entries to be 4 or 6, but had %d entries\n",len1);
#else
              printf("In real or integer mode expected num entries to be 3 or 6, but had %d entries\n",len1);
#endif
              exit(1);
	      break;
	    }
	}
    }

  mat_add_t ret = map[j_lookup_num64(j, "matid")];
  int64_t matid = get_matrixID_from_ptr(ret);
  printf("\nRead %d,%d-level matrix from json file %s\n", 
	 matrix_row_level(ret), matrix_col_level(ret),path);
  printf ("  stored matrix now has address %p", ret);
  printf ("  and matrixID %" PRId64 "\n", matid);
  free(map);

  // YOU ARE HERE
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
	if (verbose) {printf("For i=%"PRId64", and key=%s\n",i,info_name);}
	info_type_t info_type = get_info_type_from_string_name((char *)info_name);
	
	int fail =  info_set(info_type, matid, info_data);
	if (fail) {printf("panic in %s\n",__func__); exit(1);}
	if (verbose) {
	  printf("Added info field %s to info_store\n", info_name);
	  printf("  with value %s\n", info_data);
	}
      }
  }  // end of if j_key_exists for "info"

  fclose(f);
  return ret;
}




/* Python interface to return a matrixID after reading a json file */
int64_t
matrix_read_json_file_zeroize_small_scalars_matrixID(char *path, ScalarType threshold)
{
  
  mat_add_t m_ptr = matrix_read_json_file_zeroize_small_scalars(path,threshold);
  int64_t m_mID =  get_matrixID_from_ptr(m_ptr);
  return m_mID;
  
}


/*
This function reads in the matrix stored in a json file format
at location path, and modifies the small scalars in the matrices
to zeros, then stores the new matrix into the matrix store
and returns the matrixID of the resulting collapsed matrix.
The small scalars are any scalars that are less than or equal
to the value threshold.
TODO:
This function eventually should be fixed to run on COMPLEX, 
check to see if threshold has no real component, then zeroize
any scalars whose absolute value is less than threshold
TODO:
Keep in the info file any information about which matrix this
came from and another record with how much it was zeroized.
TODO: 
Change it to return a matID
*/

mat_add_t
matrix_read_json_file_zeroize_small_scalars(char *path, ScalarType threshold)
{
  int verbose = 0;
  FILE *f = fopen(path, "r");
  // printf("  path = %s\n",path);
  if (f == NULL)
  {
    printf("%s: no file found at\n\t%s:\n\texiting\n", __func__, path);
    exit(1);
  }
  json_t *j = j_parse_file(f);
  // printf("  parsed json\n");

  int64_t max_matrixID = j_lookup_num64(j, "matrixID_max");
  // printf("  max matrixID is %ld\n", max_matrixID);

  json_t *t = j_key_lookup(j, "table");

  mat_add_t *map = calloc(max_matrixID+1, sizeof(mat_add_t));
  int64_t len = j_key_count(t);
  int64_t i;

  if (verbose) {printf("In %s, before the line reader\n", __func__);}
  
  for (i=0; i<len; i++)
    {
      json_t *p = j_key_index(t, i);
      int64_t index = atoi(p->name);
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
                printf("The function %s, needs to be modified to handle COMPLEX case\n",__func__);
                exit(1);
                if (verbose) {printf("In complex scalar case \n");}
		int row_level = j_get_num64(j_array_index(p, 0));
		int col_level = j_get_num64(j_array_index(p, 1));
		// double complex_val[2];
		// complex_val[0] = j_get_double(j_array_index(p,2));
		// complex_val[1] = j_get_double(j_array_index(p,3));
                // ScalarType val = complex_val[0]+I*complex_val[1];
                ScalarType val =  j_get_double(j_array_index(p,2)) +
		        I*j_get_double(j_array_index(p,3));
                if (row_level || col_level) {
                        printf("error in %s:\n\texpected complex ", __func__ );
                        printf("scalar, but levels greater than zero!\n");
                }
		map[index] = matrix_get_ptr_scalar(val);
                if (verbose) {
		  printf("%"PRId64" %p %" PRId64 " %d %d %g+%gi\n", index, map[index], 
			 get_matrixID_from_ptr((map[index])),
			 row_level, col_level, creal(val), cimag(val));
		}
	      }
	      break;
#endif
#ifdef USE_REAL
	    case 3:			/* expect a scalar link */
	      {
                if (verbose) {printf("In real scalar case \n");}
		int row_level = j_get_num64(j_array_index(p, 0));
		int col_level = j_get_num64(j_array_index(p, 1));
		ScalarType mat_val = j_get_double(j_array_index(p,2));
                if (row_level || col_level) {
                        printf("error in %s:\n\texpected real ", __func__ );
                        printf("scalar, but levels greater than zero!\n");
                }

                // zeroize small scalars
                if (mat_val <= threshold) {mat_val = (ScalarType) 0;}

		map[index] = matrix_get_ptr_scalar(mat_val);
                if (verbose) {
		  // printf("%"PRId64" %zd %" PRId64 " %d %d %g\n", index, map[index], 
		  printf("%"PRId64" %p %" PRId64 " %d %d %g\n", index, map[index], 
			 get_matrixID_from_ptr((map[index])),
			 row_level, col_level, mat_val);
		}
	      }
              break;
#endif
#ifdef USE_INTEGER
	    case 3:			/* expect a scalar link */
	      {
                if (verbose) {printf("In integer scalar case \n");}
		int row_level = j_get_num64(j_array_index(p, 0));
		int col_level = j_get_num64(j_array_index(p, 1));
		ScalarType mat_val = j_get_num64(j_array_index(p,2));
                if (row_level || col_level) {
                        printf("error in %s:\n\texpected integer ", __func__ );
                        printf("scalar, but levels greater than zero!\n");
                }

                // zeroize small scalars
                if (mat_val <= threshold) {mat_val = (ScalarType) 0;}

		map[index] = matrix_get_ptr_scalar(mat_val);
                if (verbose) {
		  printf("%"PRId64" %p %" PRId64 " %d %d %ld\n", index, map[index], 
			 get_matrixID_from_ptr(map[index]),
			 row_level, col_level, mat_val);
		}
	      }
              break;
#endif
	      
	    case 6:			/* expect a 4-tuple */
	      {
                if (verbose) {printf("In non scalar case \n");}
		int row_level = j_get_num64(j_array_index(p, 0));
		int col_level = j_get_num64(j_array_index(p, 1));
                if (verbose) {printf("  levels are %d %d \n",row_level,col_level);}
		mat_add_t content[4];
                int temp;
                for (int k=0; k< 4; ++k) {
		   temp = j_get_num64(j_array_index(p, k+2));
                   if (verbose) {printf("temp for k = %d is %d \n",k,temp);}
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
		map[index] = matrix_get_ptr_panel(content,row_level,col_level);
                if (verbose){
		  printf("ind %" PRId64 " map %p mID %" PRId64 " (%d,%d)",
                        index, map[index], get_matrixID_from_ptr(map[index]), 
			 row_level, col_level);
		  printf(" [%p %p %p %p]\n", 
                        content[0], content[1], content[2], content[3]);
		}
	      }
	      break; 
	      
	    default:
              printf("Error in %s(%s): unexpected number of entries per line\n",
                __func__, path);
#ifdef USE_COMPLEX
              printf("In complex mode expected num entries to be 4 or 6, but had %d entries\n",len1);
#else
              printf("In real or integer mode expected num entries to be 3 or 6, but had %d entries\n",len1);
#endif
              exit(1);
	      break;
	    }
	}
    }

  mat_add_t ret = map[j_lookup_num64(j, "matid")];
  printf("\nRead %d,%d-level matrix from json file %s\n", 
	 matrix_row_level(ret), matrix_col_level(ret),path);
  printf ("  stored matrix now has address %p", ret);
  printf ("  and matrixID %" PRId64 "\n", get_matrixID_from_ptr(ret));
  free(map);
  fclose(f);
  return ret;
}










/* Python interface to return a matrixID after reading a json file */
int64_t
matrix_read_json_file_flatten_small_scalars_matrixID(char *path, ScalarType threshold, ScalarType flatValue)
{
  
  mat_add_t m_ptr = matrix_read_json_file_flatten_small_scalars(path,threshold,flatValue);
  int64_t m_mID =  get_matrixID_from_ptr(m_ptr);
  return m_mID;
  
}


/*
This function reads in the matrix stored in a json file format
at location path, and modifies the small scalars in the matrices
to zeros, then stores the new matrix into the matrix store
and returns the matrixID of the resulting collapsed matrix.
The small scalars are any scalars that are less than or equal
to the value threshold.
TODO:
This function eventually should be fixed to run on COMPLEX, 
check to see if threshold has no real component, then flatten
any scalars whose absolute value is less than threshold
TODO:
Keep in the info file any information about which matrix this
came from and another record with how much it was flattend.
TODO: 
Change it to return a matID
*/

mat_add_t
matrix_read_json_file_flatten_small_scalars(char *path, ScalarType threshold, ScalarType flatValue)
{
  int verbose = 0;
  FILE *f = fopen(path, "r");
  // printf("  path = %s\n",path);
  if (f == NULL)
  {
    printf("%s: no file found at\n\t%s:\n\texiting\n", __func__, path);
    exit(1);
  }
  json_t *j = j_parse_file(f);
  // printf("  parsed json\n");

  int64_t max_matrixID = j_lookup_num64(j, "matrixID_max");
  // printf("  max matrixID is %ld\n", max_matrixID);

  json_t *t = j_key_lookup(j, "table");

  mat_add_t *map = calloc(max_matrixID+1, sizeof(mat_add_t));
  int64_t len = j_key_count(t);
  int64_t i;

  if (verbose) {printf("In %s, before the line reader\n", __func__);}
  
  for (i=0; i<len; i++)
    {
      json_t *p = j_key_index(t, i);
      int64_t index = atoi(p->name);
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
                printf("The function %s, needs to be modified to handle COMPLEX case\n",__func__);
                exit(1);
                if (verbose) {printf("In complex scalar case \n");}
		int row_level = j_get_num64(j_array_index(p, 0));
		int col_level = j_get_num64(j_array_index(p, 1));
		// double complex_val[2];
		// complex_val[0] = j_get_double(j_array_index(p,2));
		// complex_val[1] = j_get_double(j_array_index(p,3));
                // ScalarType val = complex_val[0]+I*complex_val[1];
                ScalarType val =  j_get_double(j_array_index(p,2)) +
		        I*j_get_double(j_array_index(p,3));
                if (row_level || col_level) {
                        printf("error in %s:\n\texpected complex ", __func__ );
                        printf("scalar, but levels greater than zero!\n");
                }
		map[index] = matrix_get_ptr_scalar(val);
                if (verbose) {
		  printf("%"PRId64" %p %" PRId64 " %d %d %g+%gi\n", index, map[index], 
			 get_matrixID_from_ptr((map[index])),
			 row_level, col_level, creal(val), cimag(val));
		}
	      }
	      break;
#endif
#ifdef USE_REAL
	    case 3:			/* expect a scalar link */
	      {
                if (verbose) {printf("In real scalar case \n");}
		int row_level = j_get_num64(j_array_index(p, 0));
		int col_level = j_get_num64(j_array_index(p, 1));
		ScalarType mat_val = j_get_double(j_array_index(p,2));
                if (row_level || col_level) {
                        printf("error in %s:\n\texpected real ", __func__ );
                        printf("scalar, but levels greater than zero!\n");
                }

                // flatten small scalars
                if ((mat_val <= threshold && mat_val != (ScalarType) 0)) {mat_val = flatValue;}

		map[index] = matrix_get_ptr_scalar(mat_val);
                if (verbose) {
		  // printf("%"PRId64" %zd %" PRId64 " %d %d %g\n", index, map[index], 
		  printf("%"PRId64" %p %" PRId64 " %d %d %g\n", index, map[index], 
			 get_matrixID_from_ptr((map[index])),
			 row_level, col_level, mat_val);
		}
	      }
              break;
#endif
#ifdef USE_INTEGER
	    case 3:			/* expect a scalar link */
	      {
                if (verbose) {printf("In integer scalar case \n");}
		int row_level = j_get_num64(j_array_index(p, 0));
		int col_level = j_get_num64(j_array_index(p, 1));
		ScalarType mat_val = j_get_num64(j_array_index(p,2));
                if (row_level || col_level) {
                        printf("error in %s:\n\texpected integer ", __func__ );
                        printf("scalar, but levels greater than zero!\n");
                }

                // flatten small scalars
                if (mat_val <= threshold) {mat_val = flatValue;}

		map[index] = matrix_get_ptr_scalar(mat_val);
                if (verbose) {
		  printf("%"PRId64" %p %" PRId64 " %d %d %ld\n", index, map[index], 
			 get_matrixID_from_ptr(map[index]),
			 row_level, col_level, mat_val);
		}
	      }
              break;
#endif
	      
	    case 6:			/* expect a 4-tuple */
	      {
                if (verbose) {printf("In non scalar case \n");}
		int row_level = j_get_num64(j_array_index(p, 0));
		int col_level = j_get_num64(j_array_index(p, 1));
                if (verbose) {printf("  levels are %d %d \n",row_level,col_level);}
		mat_add_t content[4];
                int temp;
                for (int k=0; k< 4; ++k) {
		   temp = j_get_num64(j_array_index(p, k+2));
                   if (verbose) {printf("temp for k = %d is %d \n",k,temp);}
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
		map[index] = matrix_get_ptr_panel(content,row_level,col_level);
                if (verbose){
		  printf("ind %" PRId64 " map %p mID %" PRId64 " (%d,%d)",
                        index, map[index], get_matrixID_from_ptr(map[index]), 
			 row_level, col_level);
		  printf(" [%p %p %p %p]\n", 
                        content[0], content[1], content[2], content[3]);
		}
	      }
	      break; 
	      
	    default:
              printf("Error in %s(%s): unexpected number of entries per line\n",
                __func__, path);
#ifdef USE_COMPLEX
              printf("In complex mode expected num entries to be 4 or 6, but had %d entries\n",len1);
#else
              printf("In real or integer mode expected num entries to be 3 or 6, but had %d entries\n",len1);
#endif
              exit(1);
	      break;
	    }
	}
    }

  mat_add_t ret = map[j_lookup_num64(j, "matid")];
  printf("\nRead %d,%d-level matrix from json file %s\n", 
	 matrix_row_level(ret), matrix_col_level(ret),path);
  printf ("  stored matrix now has address %p", ret);
  printf ("  and matrixID %" PRId64 "\n", get_matrixID_from_ptr(ret));
  free(map);
  fclose(f);
  return ret;
}





