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
    if (!numZeros){
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
                        fprintf(stderr,"ERROR in %s(%s): entry 2 of line not string (%ld of table)\n", __func__, larcMatrixFile_path, i);
                        exit(1);
                    }
                    // set temp scalarType with entry string
                    sca_set_str(&temp, j_get_string(j_array_index(entry, 2)));
                    //printf("  entry %ld: scalar %s\n", i, sca_get_str(temp));
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
                        fprintf(stderr,"ERROR in %s(%s): entry 2 of line not string (%ld of table)\n", __func__, path, i);
                        exit(1);
                    }
                    // set temp scalarType with entry string
                    sca_set_str(&temp, j_get_string(j_array_index(entry, 2)));
                    //printf("  entry %ld: scalar %s\n", i, sca_get_str(temp));
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
    if (!locs){
        fprintf(stderr,"error in %s: failed to allocate location info array with size %ld.\n", __func__, total);
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
    for (i = 0; mpz_cmp_ui(ret, i) > 0; i++){
        locations2[i] = malloc(3*sizeof(int64_t));
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


