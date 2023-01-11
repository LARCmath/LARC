//                       experimental.c 
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

// Standard Libaries
#include <gmp.h>

// Our header files structures and functions
#include "larc.h"
#include "matrix_store.h"
#include "op_store.h"
#include "global.h"
#include "organize.h"
#include "version.h"
#include "info_store.h"
#include "experimental.h"
#include "scalars.h"
#include "matmath.h"

/*!
 * \file experimental.c
 * \brief This file contains several utility routines which are 
 * either not fully tested or scale poorly. Use with caution.
 */

// This routine is experimental and the routine 
//  matrix_count_entries
// in matmath.c is preferred. 
void mpz_count_entries_larcMatrixFile(mpz_t sca_count, char *larcMatrixFile_path, scalarType scalar)
{
    FILE *file = fopen(larcMatrixFile_path, "r");
    int64_t i;

    if (file == NULL)
    {
        fprintf(stderr,"%s: no file found at\n\t%s:\n\texiting\n", __func__, larcMatrixFile_path);
        exit(1);
    }
    json_t *jData = j_parse_file(file);
    //printf("  parsed json\n");

    int64_t matrixID = j_lookup_num64(jData, "matid");
    //printf("  matrixID is %ld\n", matrixID);
    json_t *jTable = j_key_lookup(jData, "table");

    // check that we can allocate an array this large
    // (i.e., total number of bytes does not exceed SIZE_MAX)
    if (SIZE_MAX/(matrixID + 1)/sizeof(mpz_t) == 0) {
        fprintf(stderr,"Error is %s: matrixID_max too large - try renumbering.\n", __func__);
        exit(1);
    }
    // originally written for counting zeros - numZeros will actually count
    // number of entries of scalar. 
    mpz_t *numZeros = (mpz_t *) malloc((matrixID + 1) * sizeof(mpz_t));
    if (numZeros == NULL){
        fprintf(stderr,"Error in %s: failed to allocate array.\n", __func__);
        exit(1);
    }
    //printf("declared an array of mpz_t\n");
    // init and set the mpz_t's to 0's. 
    for (i = 0; i < matrixID + 1; i++){
        mpz_init_set_ui(numZeros[i], 0);
    }
    //printf("  array set to zeros\n");

    // for each key/matrixID in table, count number of zeros
    int64_t jTableLen = j_key_count(jTable);
    scalarType temp;
    sca_init(&temp);
    int found = 0; //to record whether we have found table entry for scalar yet. 
//    float next = 5;
    for (i = 0; i < jTableLen; i++){
//        if (100*i / jTableLen >= 5 && 100*i/jTableLen >= next){
//            printf("now %02f through table.\n", next);
//            next += 5;
//        }
        json_t *entry = j_key_index(jTable, i);
        int64_t matID = (int64_t) atoll(entry->name);

        if (j_is_array(entry)){
            // If we have found the entry in the matrix holding scalar, then 
            // we don't need to consider SCALAR entries anymore. 
            if (!found){
                // SCALAR entries have length 3
                if (3 == j_array_get_length(entry)){
                    // get string value of table entry
                    if (!j_is_string(j_array_index(entry,2))){
                        fprintf(stderr,"ERROR in %s(%s): entry 2 of line not string (%" PRId64 " of table)\n",
                                __func__, larcMatrixFile_path, i);
                        exit(1);
                    }
                    // set temp scalarType with entry string
                    sca_set_str(&temp, j_get_string(j_array_index(entry, 2)));
                    //printf("  entry %ld: scalar %s\n", i, sca_get_readable_approx_str(temp));
                    // if temp == scalar, we found the entry so set count to one
                    // for this matID. 
                    if (0 != sca_eq(temp, scalar)){
                        mpz_set_ui(numZeros[matID], 1);
                        // record that we don't need to check scalars anymore. 
                        found = 1;
                    }
                }
            }
            // If we have not found the entry in the matrix holding scalar, 
            // then the counts of any matrices from table entries must be zero
            // so not worth looking at them. 
            else {
                // SUBMATRIX entries have length 6
                if (6 == j_array_get_length(entry)){
                    for (int j = 0; j < 4; j++){
                        // get the number of zeros for each submatrix and 
                        // add to the total count for current matID. 
                        int64_t subMatID = j_get_num64(j_array_index(entry, j+2));
                        if (subMatID != MATRIX_ID_INVALID)
                            mpz_add(numZeros[matID], numZeros[matID], numZeros[subMatID]);
                    }
                }
            }
        } // end (if j_is_array(entry))
    } // end loop over len

//    //Debug
//    for (i = 0; i < matrixID + 1; i ++){
//        gmp_printf("MatrixID %ld has %Zd entries\n", i, numZeros[i]);
//    }
    
    //printf("Preparing output\n");
    mpz_set(sca_count, numZeros[matrixID]);

    // cleaning
    sca_clear(&temp);
    for(i = 0; i < matrixID + 1; i ++)
        mpz_clear(numZeros[i]);
    free(numZeros);
}

// right now this only works for scalarType MPINTEGER
char *count_entries_larcMatrixFile(char *larcMatrixFile_path, char *scalar_str)
{
    scalarType scalar;
    mpz_t ret;

    // interpret scalar_str as scalarType
    sca_init(&scalar);
    sca_set_str(&scalar, scalar_str);
    
    // run count entries routine
    mpz_init(ret);
    mpz_count_entries_larcMatrixFile(ret, larcMatrixFile_path, scalar);

    // prepare string output from mpz_t ret
    size_t out_str_size = mpz_sizeinbase(ret, 10) + 2;
    char *out_str = calloc(out_str_size, sizeof(char));
    if (NULL == out_str){
        fprintf(stderr,"ERROR: allocation failure in %s.\n", __func__);
        exit(1);
    }
    gmp_snprintf(out_str, out_str_size, "%Zd", ret);

    // clean scalarType and mpz_t variables
    sca_clear(&scalar);
    mpz_clear(ret);

    return out_str;
}

struct count_info {
    mpz_t count;
    int64_t subMatIDs[4];
};

struct sca_info {
    int64_t matID;
    int64_t row;
    int64_t col;
};

// This routine is in development and does not have a mat_ptr_t equivalent
// in matmath.c
// It works for small matrices but according to the TODO in the code, doesn't 
// always work. 
void mpz_locate_entries_larcMatrixFile(int64_t ***locations_ptr, mpz_t count, char *path, scalarType scalar, unsigned maxno)
{
    FILE *file = fopen(path, "r");
    int64_t i; 

    if (file == NULL)
    {
        fprintf(stderr,"%s: no file found at\n\t%s:\n\texiting\n", __func__, path);
        exit(1);
    }
    json_t *jData = j_parse_file(file);
    //printf("  parsed json\n");

    int64_t matrixID = j_lookup_num64(jData, "matid");
    //printf("  matrixID is %ld\n", matrixID);
    json_t *jTable = j_key_lookup(jData, "table");

    // check that we can allocate an array this large
    // (i.e., total number of bytes does not exceed SIZE_MAX)
    if (SIZE_MAX/(matrixID + 1)/sizeof(mpz_t) == 0) {
        fprintf(stderr,"Error is %s: matrixID_max too large - try renumbering.\n", __func__);
        exit(1);
    }

    struct count_info *records = malloc((matrixID + 1) * sizeof(struct count_info));
    if (records == NULL) { ALLOCFAIL(); }
    // matrixID + 1 because we'll be indexing by [0, matrixID]. 
    for (i = 0; i < matrixID + 1; i ++){
        mpz_init(records[i].count);
        mpz_set_ui(records[i].count, 0);
    }

    // for each key/matrixID in table, count number of zeros
    int64_t jTableLen = j_key_count(jTable);
    int64_t scaID = -1; //to record whether we have found table entry for scalar yet. 
    scalarType temp;
    sca_init(&temp);
    for (i = 0; i < jTableLen; i++){
        json_t *entry = j_key_index(jTable, i);
        int64_t matID = (int64_t) atoll(entry->name);

        if (j_is_array(entry)){
            // If we have found the entry in the matrix holding scalar, then 
            // we don't need to consider SCALAR entries anymore. 
            if (scaID == -1){
                // SCALAR entries have length 3
                // if scalar value == scalar, record 1 for entry matID
                if (3 == j_array_get_length(entry)){
                    // get string value of table entry
                    if (!j_is_string(j_array_index(entry,2))){
                        fprintf(stderr,"ERROR in %s(%s): entry 2 of line not string (%" PRId64 " of table)\n",
                                __func__, path, i);
                        exit(1);
                    }
                    // set temp scalarType with entry string
                    sca_set_str(&temp, j_get_string(j_array_index(entry, 2)));
                    //printf("  entry %ld: scalar %s\n", i, sca_get_readable_approx_str(temp));
                    // if temp == scalar, we found the entry so set count to one
                    // for this matID. 
                    if (0 != sca_eq(temp, scalar)){
                        // record that we don't need to check scalars anymore. 
                        // by recording matrixID of scalar
                        scaID = matID;
                        // record count of one
                        mpz_set_ui(records[matID].count, 1);
                    }
                }
            }
            // If we have not found the entry in the matrix holding scalar, 
            // then the counts of any matrices from table entries must be zero
            // so not worth looking at them. 
            else {
                // SUBMATRIX entries have length 6
                if (6 == j_array_get_length(entry)){
                    for (int j = 0; j < 4; j++){
                        // get the number of zeros for each submatrix and 
                        // add to the total count for current matID. 
                        int64_t subMatID = j_get_num64(j_array_index(entry, j+2));
                        if (subMatID != MATRIX_ID_INVALID)
                            mpz_add(records[matID].count, records[matID].count, records[subMatID].count);
                        // record the submatrices that contributed 
                        records[matID].subMatIDs[j] = subMatID;
                    }
                }
            }
        } // end (if j_is_array(entry))
    } // end loop over len
    sca_clear(&temp);
    j_set_null(jData); free(jData);


//    //Debug
//    for (i = 0; i < matrixID + 1; i ++){
//        gmp_printf("MatrixID %ld has %Zd entries\n", i, numZeros[i]);
//    }
    
    
    // if there were no entries with scalar, return NULL for empty list 
    if (mpz_cmp_ui(records[matrixID].count, 0) == 0){
        mpz_set_ui(count, 0);
        *locations_ptr = NULL;
        free(records);
        return;
    }
    
    // (~arbitrarily) we'll say that if there are more than 2*32 appearances
    // of a number, we'll refuse to locate them all. In particular, we'll 
    // protect against overflowing <<total>>. 
    if (mpz_cmp_ui(records[matrixID].count, maxno) > 0){
        fprintf(stderr,"error in %s: too many appearances of scalar to return.\n", __func__);
        exit(0);
    }

    // allocate an array of info for each scalar occurence 
    int64_t total = mpz_get_si(records[matrixID].count);
    //TODO: in my example, this total value is coming out WAY high. for rd11.json
    // error in mpz_locate_entries_larcMatrixFile: failed to allocate loc array with total -172825635700342774.
    // create and initialize array of length numZeros[matrixID] and width 3
    // [current matrixID, partial row bin expansion, partial col bin expansion]
    struct sca_info *locs = calloc(total, sizeof(struct sca_info));
    if (locs == NULL){
        fprintf(stderr,"error in %s: failed to allocate location info array with size %" PRId64 ".\n", __func__, total);
        exit(0);
    }
    for (i = 0; i < total; i++){
        // set current matrixID to final matrix matrixID for each scalar
        locs[i].matID = matrixID;
    }

    // until we hit the scalar level, make passes through the list of scalar
    // locations
    //int p = 0; // temp var for debug statement below
    while (locs[0].matID != scaID){
        for (i = 0; i < total; ){ // incrementing happens inside loop
            // for a given scalar location, consider the matrixID to do next
            int64_t curID = locs[i].matID;
            int64_t j,k;
            // go through each of the 4 quadrants 
            for (j = 0; j < 4; j++){ 
                // make edits to / replace matrixID-to-do-next with quadrant
                int64_t subMatID = records[curID].subMatIDs[j];
                // and do this to exactly the number of scalars that the  
                // quadrant contributes. Note that since total fits in an
                // int64_t, records[subMatID].count must be able to as well. 
                if (subMatID != MATRIX_ID_INVALID){
                    for (k = 0; mpz_cmp_ui(records[subMatID].count, k) > 0; k++){
                        // reset next matrixID
                        locs[i+k].matID = subMatID;
                        // include row/col location info (add and shift)
                        locs[i+k].row = (locs[i+k].row << 1) + (j>1);
                        locs[i+k].col = (locs[i+k].col << 1) + (j%2);
                        // debug: 
                        //printf("touched ind %ld on pass %d: subMatID = %ld, curID = %ld: j=%ld -> add %d to r, add %ld to c\n", i+k, p, subMatID, curID, j, j>1, j%2);
                    }
                    i += k; //increment i by the progress we've made, namely k=records[subMatID]
                }
            }
        }
        //p += 1; //part of debug print statement
    }

//    // debug
//    for (i = 0; i < total; i++){
//        printf("scalar occ %ld: row %ld col %ld\n", i, locs[i].row, locs[i].col);
//    }

    //printf("Preparing output\n");
    mpz_set(count, records[matrixID].count);
    *locations_ptr = malloc(total * sizeof(int64_t*));
    if (locations_ptr == NULL) { ALLOCFAIL(); }
    for (i = 0; i < total; i++){
        (*locations_ptr)[i] = malloc(2 * sizeof(int64_t));
        (*locations_ptr)[i][0] = locs[i].row;
        (*locations_ptr)[i][1] = locs[i].col;
    }

    for(i = 0; i < matrixID + 1; i ++)
        mpz_clear(records[i].count);
    free(records);
    free(locs);
}

// Not to be used outside the SWIG interface (i.e. with Python). 
// If programming in C, use mpz_locate_entries_larcMatrixFile. 
int64_t **locate_entries_larcMatrixFile(char *path, char *scalar_str, char *maxno_str)
{
    scalarType scalar;
    mpz_t ret;
    int64_t **locations;
    int64_t **locations2;

    // interpret scalar_str as scalarType
    sca_init(&scalar);
    sca_set_str(&scalar, scalar_str);

    // interpret maxno as unsigned int
    char *ptr; // needed for call, not used otherwise
    unsigned maxno = (unsigned)strtoul(maxno_str,&ptr,10);

    // run count entries routine
    mpz_init(ret);
    mpz_locate_entries_larcMatrixFile(&locations, ret, path, scalar, maxno);

    // prepare for swig interface
    // as of now, the ends of the r/c lists are marked with -1
    // and the end of the array is marked with NULL
    int64_t i;
    locations2 = malloc((mpz_get_ui(ret) + 1) * sizeof(int64_t *));
    if (locations2 == NULL) { ALLOCFAIL(); }
    for (i = 0; mpz_cmp_ui(ret, i) > 0; i++){
        locations2[i] = malloc(3*sizeof(int64_t));
        if (locations2[i] == NULL) { ALLOCFAIL(); }
        locations2[i][0] = locations[i][0];
        locations2[i][1] = locations[i][1];
        locations2[i][2] = -1;
    }
    locations2[i] = NULL;

    // clean scalarType and mpz_t variables
    sca_clear(&scalar);
    mpz_clear(ret);
    free(locations); 

    return locations2;
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
  if (map==NULL){
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
                  char *mat_val_string = sca_get_readable_approx_str(*mat_val);
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
                   if (verbose) {printf("temp for k = %d is %" PRId64 " \n",k,temp);}
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
	 reti->row_level, ret->col_level, path);
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
mat_ptr_t matrix_read_larcMatrix_file_zeroize_small_scalars_PTR(char *path, scalarType threshold)
{
    mat_ptr_t mat_ptr = matrix_read_larcMatrix_file_flatten_small_scalars_PTR(path, threshold, scalar0);

    return mat_ptr;
}

/* Python interface to return a matrixID after reading a json file */
int64_t matrix_read_larcMatrix_file_zeroize_small_scalars_matrixID(char *path, char *threshold)
{
    scalarType thresh;
    sca_init(&thresh);
    sca_set_str(&thresh, threshold);

    mat_ptr_t m_ptr = matrix_read_larcMatrix_file_zeroize_small_scalars_PTR(path, thresh);
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
static void flatten_scalar(scalarType *flat, const scalarType flatten_this)
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
mat_ptr_t matrix_read_larcMatrix_file_flatten_small_scalars_PTR(char *path, scalarType threshold, scalarType flatValue)
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

  mat_ptr_t m_ptr = matrix_read_larcMatrix_file_flatten_small_scalars_PTR(path, thresh, flatVal);
  sca_clear(&thresh);
  sca_clear(&flatVal);

  return get_matID_from_matPTR(m_ptr);
}

/*!
 * \ingroup larc
 * \brief Writes all info_store information about a given matrix to a json file
 *
 * This code is an exact copy of a static function in io.c as it existed on
 * 3 September 2020. It is needed here for
 * writePTR_and_alter_vals_larcMatrix_file(), which was formerly in io.c but
 * was moved here as a more suitable location.
 *
 * \param m_ptr Pointer to the matrix to be stored in compressed json format
 * \param f File pointer for the json file where the metadata will be written
 * \return 1 on success
 */
static int write_infoStore_to_larcMatrix_file(mat_ptr_t m_ptr, FILE *f) {
  int verbose = 0;

  if (!f)
    return -1;
  if (m_ptr==MATRIX_PTR_INVALID) 
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

/*!
 * \ingroup larc
 * \brief A worker routine for writePTR_and_alter_vals_larcMatrix_file() 
 * \param m_ptr The pointer to the matrix to be written
 * \param f A file pointer
 * \param output_flags An array of flags used to track already-written matrixIDs
 * \param func A function that will be applied to each scalar before writing
 * \return 1 on success
 */
static int
recursive_write_and_alter_larcMatrix_file_by_matPTR(mat_ptr_t m_ptr, FILE *f,
   char *output_flags, void (*func)(scalarType*, const scalarType))
{
  if (m_ptr==MATRIX_PTR_INVALID) {
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
      func(val, m_ptr->scalar_value);
      char *val_string = sca_get_readable_approx_str(*val);
      fprintf(f, "    \"%" PRId64 "\":[%d, %d, \"%s\"],\n",
         get_matID_from_matPTR(m_ptr), 0, 0, val_string);
      free(val_string);
    }
  else   // NONSCALAR
    {
      int ret;
     
      // would need to be fixed:
      // mat_ptr_t -> either matns_ptr_t or mats_ptr_t, depending on the
      // level of the matrix m_ptr; code would need to work with either type
      mat_ptr_t panel[4];
      for (int i = 0; i < 4; ++i) {
        record_ptr_t rec_ptr;
        check_validity_one_input(m_ptr->subMatList[i],__func__,&rec_ptr);
        panel[i] = (mat_ptr_t)rec_ptr;
      }

      // obtain the pointers for the matrices in the panel
      ret = recursive_write_and_alter_larcMatrix_file_by_matPTR(panel[0], f,
		output_flags, func);
      if (ret < 0) return ret;

      if (mtype!=COL_VECTOR) {
        ret = recursive_write_and_alter_larcMatrix_file_by_matPTR(panel[1], f,
		output_flags, func);
        if (ret < 0) return ret;
      }

      if (mtype!=ROW_VECTOR) {
        ret = recursive_write_and_alter_larcMatrix_file_by_matPTR(panel[2], f,
		output_flags, func);
        if (ret < 0) return ret;
      }

      if (mtype==MATRIX) {
        ret = recursive_write_and_alter_larcMatrix_file_by_matPTR(panel[3], f,
		output_flags, func);
        if (ret < 0) return ret;
      }

      // print the panel, which contains matrixIDs for these matrices 
      fprintf(f, "    \"%" PRId64 "\":[%d, %d, %" PRId64 ", ", get_matID_from_matPTR(m_ptr), 
                m_ptr->row_level, m_ptr->col_level,
	        get_matID_from_matPTR(panel[0]));

      if (mtype!=COL_VECTOR) {
                fprintf(f, "%" PRId64 ", ", get_matID_from_matPTR(panel[1]));
      }
      else { fprintf(f, "-1, "); }

      if (mtype!=ROW_VECTOR) {
                fprintf(f, "%" PRId64 ", ", get_matID_from_matPTR(panel[2]));
      }
      else { fprintf(f, "-1, "); }

      if (mtype==MATRIX) {
                fprintf(f, "%" PRId64 "],\n", get_matID_from_matPTR(panel[3]));
      }
      else { fprintf(f, "-1],\n"); }

      //fprintf(f, "    \"%ld\":[%d, %d, %ld, %ld, %ld, %ld],\n", get_matID_from_matPTR(m_ptr), 
      //      m_ptr->row_level,m_ptr->col_level,
      //      get_matID_from_matPTR(panel[0]),get_matID_from_matPTR(panel[1]),
      //      get_matID_from_matPTR(panel[2]),get_matID_from_matPTR(panel[3]));
    }
  return 1;
}

int writePTR_and_alter_vals_larcMatrix_file(mat_ptr_t m_ptr, char *path,
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
  if (m_ptr==MATRIX_PTR_INVALID) 
    {
        fprintf(stderr,"in %s, passed matrix pointer is invalid\n",__func__);
        return -1;
    }
  char *output_flags = calloc(num_matrices_created(), sizeof(char));
  if (output_flags == NULL) { ALLOCFAIL(); }

  if (!output_flags) 
    return -1;

  printf("\nWriting %d,%d-level matrix with matrixID %" PRId64 " to %s\n", 
	 m_ptr->row_level, m_ptr->col_level,
         get_matID_from_matPTR(m_ptr), path);
  // printf("      the matrix id is %lu, with level %d %d\n",id, m_ptr->row_level, m_ptr->col_level);

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
  ret = recursive_write_and_alter_larcMatrix_file_by_matPTR(m_ptr, f,
	output_flags, func);
  free (output_flags);

  // end the table structure and the json file
  fprintf(f, "      \"end\":0 }\n}\n");
  fclose(f);

  return ret;
}



