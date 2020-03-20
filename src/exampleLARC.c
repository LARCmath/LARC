//                    exampleLARC.c
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
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <pthread.h>
#include <math.h>
#include <inttypes.h> //for printing int64_t and uint64_t

// Our header files structures and functions
#include "global.h"
#include "larc.h"
#include "hash.h"
#include "io.h"
#include "json.h"
#include "matmath.h"
#include "matrix_store.h"
#include "op_store.h"
#include "organize.h"
#include "scalars.h"
#include "experimental.h"


// TODO: does -1 have a global name
// how to handle examples




/* A routine to negate scalar of any scalarType.
 * Result is passed by <negated> parameter.*/ 
void negate(scalarType *negated, const scalarType val){
    scalarType neg1;
    sca_init(&neg1);
    sca_set_str(&neg1, "-1");

    sca_mult(negated, neg1, val);
    sca_clear(&neg1);
}

/* Rounds input to nearest integer and then returns string version of 
 * integer. */
char *real_to_str(const long double val){
    char *out = calloc(21, sizeof(char));
    int64_t val_rounded = rint(val);
    snprintf(out, 20, "%ld", val_rounded);
    return out;
}


// this is a goofy function, that illustrates what you can do
char *zero_to_one_else_two(const scalarType val){
    // make a scalarType with value 0.
    scalarType zero;
    sca_init(&zero);
    sca_set_str(&zero, "0");

    // compare the incoming val. if it is 0, make it 1.
    int val_is_zero = sca_eq(val, zero);
    sca_clear(&zero);

    int size_of_output = 2;
    char *out = calloc(size_of_output + 1, sizeof(char));
    if (val_is_zero)
        snprintf(out, size_of_output, "1");
    else
        snprintf(out, size_of_output, "2");
    return out;
}



int main (int argc, char *argv[])
{
  initialize_larc(25, 26, 10, -1, -1, 1);

    /*
     * Example 0:
     * [type=ANY]
     * -----------
     *
     * Create some matices, do some operations on them, print some store
     * reports. 
     *
     * Step 1: we create some matrices from scratch. 
     *
     * Step 2: we add the matrices together, square the entries of the result, 
     * and multiply the entries of the result by 7. 
     *
     * Step 3: we check the size of the matrix store, print out a matrix store
     * report, remove all but the original two matrices from the store, 
     * clean the store, and print out a report. 
     *
     */

    printf("Running Example 1.\n");
    // Step 1: build some matrices.
    // one way to create a matrix is to build the scalarType values 
    // from scratch using the scalar functions.
    scalarType scalar0;
    sca_init(&scalar0);
    sca_set_str(&scalar0, "7"); // 3 exists in all TYPES

    // produce a matrix pointer for scalar 7
    mat_ptr_t seven_mat = get_valMatPTR_from_val(scalar0);
    printf("Here is a matrix with one entry:\n");
    print_naive_by_matPTR(seven_mat);
    
    // alternatively, if you don't care to keep your code `TYPE agnostic'
    // you can declare your scalar as you normally would (using the actual
    // type or scalarType)
    // Eg. REAL, INTEGER, and COMPLEX don't need to be initialized
    // so you could skip that step. 
#ifdef USE_REAL
    long double scalar1 = 2.333;
    scalarType scalar2 = -3;
#elif defined(USE_INTEGER)
    int64_t scalar1 = 2;
    scalarType scalar2 = -3;
#elif defined(USE_COMPLEX)
    long double complex scalar1 = 2+I*1;
    scalarType scalar2 = -3;
#elif defined(USE_MPINTEGER)
    mpz_t scalar1;
    mpz_init(scalar1);
    mpz_ui_pow_ui(scalar1, 2, 65);
    scalarType scalar2;  
    //notice that mpz_init is called on mpz_t but sca_init is called on pointer 
    //to mpz_t
    sca_init(&scalar2);  
    mpz_set_si(scalar2, -3);
#elif defined(USE_MPRATIONAL)
    mpq_t scalar1;
    mpq_init(scalar1);
    mpq_set_ui(scalar1, 3, 5);
    scalarType scalar2;
    sca_init(&scalar2);
    mpq_set_si(scalar2, -3, 1);
#elif defined(USE_MPRATCOMPLEX)
    larc_mpratcomplex_t scalar1;
    sca_init(&scalar1);
    mpq_set_ui(scalar1->real, 3, 5);
    mpq_set_ui(scalar1->imag, 2, 5);
    scalarType scalar2;
    sca_init(&scalar2);
    mpq_set_si(scalar2->real, -3, 1);
    mpq_set_si(scalar2->imag, -1, 1);
#elif defined(USE_MPREAL)
    mpfr_t scalar1;
    sca_init(&scalar1);
    sca_set_2ldoubles(&scalar1, 3.0, 0.0);
    scalarType scalar2;
    sca_init(&scalar2);
    sca_set_2ldoubles(&scalar2, 4.7, 0.0);
#elif defined(USE_MPCOMPLEX)
    mpc_t scalar1;
    sca_init(&scalar1);
    sca_set_2ldoubles(&scalar1, 3.0, 1.0);
    scalarType scalar2;
    sca_init(&scalar2);
    sca_set_2ldoubles(&scalar2, 4.7, -1.0);
#endif

    // lastly, we could grab an already initialized scalar from a short list
    // used in preloading matrices
    mat_ptr_t mat_ptr_scalarM1 = get_matPTR_from_matID(matID_scalarM1, "", __func__, 0);
    scalarType scalar3;
    sca_init(&scalar3);
    sca_set(&scalar3, matrix_trace(mat_ptr_scalarM1));   // so scalar4 = -1.
    // I think you could also do this:
    // scalarType *scalar3 = &matrix_trace(mat_ptr_scalarM1);
    // but then you (probably) don't want to change it's value or clear it
    // because that would effect the one in the store (which intuitively should
    // not change... -1 should stay -1). 

    // the matrix can be built from a list in row major form
    scalarType *vals = malloc(4 * sizeof(scalarType));
    // apparently you can't do something like:
    //vals[0] = scalar0; 
    // with mpz_t and mpq_t so we have to reproduce these scalar values
    for (int i = 0; i < 4; i++)
        sca_init(&(vals[i]));
    sca_set(&(vals[0]), scalar0); 
    sca_set(&(vals[1]), scalar1); 
    sca_set(&(vals[2]), scalar2); 
    sca_set(&(vals[3]), scalar3);
    mat_ptr_t mat_ptr = row_major_list_to_store(vals, 1, 1, 2);

    // now that the scalars are incorporated into matrices, we can clean them
    // up.
    sca_clear(&scalar0);
    sca_clear(&scalar1);
    sca_clear(&scalar2);
    sca_clear(&scalar3);
    for (int i = 0; i < 4; i++)
        sca_clear(&(vals[i]));
    free(vals);

    // then we can print the matrix to screen
    printf("Here is our first matrix:\n");
    print_naive_by_matPTR(mat_ptr);

    // it might be easier to build from a list of character strings
    // (from the python interface)
    char **vals2 = calloc(4, sizeof(char*));
    vals2[0] = vals2[3] = "1";
    vals2[1] = "2";
    vals2[2] = "0"; 
    int64_t matID = row_major_list_to_store_matrixID(vals2, 1, 1, 2);
    free(vals2);
    mat_ptr_t mat_ptr2 = get_matPTR_from_matID(matID, "", __func__, 0);
    printf("Here is another matrix:\n");
    print_naive_by_matPTR(mat_ptr2);

    // Step 2: 
    // we add the matrices together, 
    mat_ptr_t sum = matrix_add(mat_ptr, mat_ptr2);
    printf("When you add them together you get:\n");
    print_naive_by_matPTR(sum);
    // square the entries of the result, 
    // and multiply the entries of the result by 7. 
    mat_ptr_t squared_sum = matrix_entrySquared(sum, seven_mat);
    printf("Squaring entries and multiplying by seven gives:\n");
    print_naive_by_matPTR(squared_sum);

    // Step 3:
    // check the size of the matrix store now
    printf("The size of the matrix store is %ld.\n", matrix_store_matrixCount() + matrix_store_scalarCount());

    // print a report on the matrix store to screen
    matrix_store_report("stdout");

    // hold the two matrices we started with
    set_hold_matrix(mat_ptr);
    set_hold_matrix(mat_ptr2);
    // clear the matrix store of everything else and repair the op store
    // to account for the lost matrices. 
    clean_matrix_store();
    empty_op_store();
    printf("Now the size of the matrix store is %ld.\n", matrix_store_matrixCount() + matrix_store_scalarCount());
    // for good measure, let's delete the first matrix. 
    release_hold_matrix(mat_ptr);
    remove_matrix_from_mat_store_by_matrix(mat_ptr);
    printf("Now the size of the matrix store is %ld.\n", matrix_store_matrixCount() + matrix_store_scalarCount());

    /*
     * Example 1: 
     * [type=REAL]
     * ------------
     * #json #reading #writing #convertingtypes
     *
     * Reading and writing matrices from/to json files: 
     *
     * Step 1: We create a REAL valued matrix and write the matrix
     * to a json file. 
     *
     * Step 2: we use a custom function that negates real values and read in 
     * the json file such that we get a new matrix where all the entries are 
     * negated.
     *
     * Step 3: we use a custom function that rounds real values to integers
     * and write out the matrix so that it can be read in with type INTEGER. 
     */

#ifdef USE_REAL 
    printf("Running Example 1.\n");
    // Step 1: 
    // I'll just pick some values to use as entries to our sample matrix:
    long double valsE1[4] = {0.5L, 0.3333L, -3.0L, 8.2223L};
    mat_ptr_t mat_ptrE1 = row_major_list_to_store(valsE1, 1, 1, 2);

    char *filepath = "../tests/dat/out/example1_original.json";
    write_larcMatrix_file_by_matPTR(mat_ptrE1, filepath);
    printf("The original matrix is:\n");
    print_naive_by_matPTR(mat_ptrE1);
    printf("\n\n");

    // Step 2: 
    // function <negate> defined above multiplies a scalar by -1 (for any
    // scalar type- for REAL in this case). 
    mat_ptr_t neg_mat_ptr = read_and_alter_vals_larcMatrix_file_return_matPTR(filepath, negate);
    // verify that the following file is as expected. 
    write_larcMatrix_file_by_matPTR(neg_mat_ptr, "../tests/dat/out/example1_negated.json");
    printf("The negated matrix is:\n");
    print_naive_by_matPTR(neg_mat_ptr);
    printf("\n\n");

    // Step 3: 
    // At any given time, we can set custom functions in for the scalar 
    // operations package.
    // In this case, we change the sca_get_str routine to be the real_to_str
    // routine written above that first rounds real numbers before converting
    // them to integer strings. 
    define_sca_get_str(real_to_str);
    write_larcMatrix_file_by_matPTR(neg_mat_ptr, "../tests/dat/out/example1_negAndRounded.json");
    printf("The integer negated matrix is:\n");
    print_naive_by_matPTR(neg_mat_ptr);
    // reset the scalar operations to normal
    init_arithmetic_scalarOps(1);
#else
    printf("Example 1 only works with TYPE=REAL.\n");
#endif


    /*
     * Example 2: [any type]
     * ------------
     * #json #random #sparsity #count #locate #scalars
     *
     * Creating random (binary) matrices. 
     *
     * Step 1: We create a random matrix by cooking up the entries and building
     * the matrix from scratch. 
     *
     * Step 2: We also use a specialized routine that produces a random binary
     * matrix with a specified sparsity. (Particularly useful when used to
     * compare against structured binary matrices with matching sparsity -
     * difference in size should be due to 'structure'.)
     * 
     * Step 3: We count 1's in the matrix and locate the 1's. 
     *
     * Step 4: We further manipulate a binary random matrix. 
     */

    printf("Running Example 2.\n");

    // Step 1: Create a random matrix from scratch. 
    // To avoid checking the type, I'll choose some 'random' values common
    // to all types. Trust me, they're random. 
    char **vals3 = calloc(4, sizeof(char*));
    vals3[0] = "1";
    vals3[1] = "9";
    vals3[2] = "11"; 
    vals3[3] = "2";
    matID = row_major_list_to_store_matrixID(vals3, 1, 1, 2);
    free(vals3);
    mat_ptr = get_matPTR_from_matID(matID, "", __func__, 0);

    printf("A small 'random' matrix:\n");
    print_naive_by_matPTR(mat_ptr);
    printf("\n\n");

    // Step 2: This routine produces a matrix for us: we input how many
    // of a single nonzero scalar we want. In the matrixID version, the
    // scalar has to be "1", but we will use "2" to be fancy. 
    scalarType scalar;
    sca_init(&scalar);
    sca_set_str(&scalar, "2");
    mat_ptr_t scalar_ptr = get_valMatPTR_from_val(scalar);
    // now we have a pointer to the scalar 2
    // we'll ask that the matrix be level 3 square, 
    mat_level_t row_level = 3;
    mat_level_t col_level = 3;
    // we'll ask that the matrix have ten 2's appear in it. 
    mpz_t scalarNum;
    mpz_init(scalarNum);
    mpz_set_ui(scalarNum, 10);

    // create and write random matrix to file. 
    mat_ptr_t randmat = random_bool_matrix_from_count(row_level, col_level, scalarNum, scalar_ptr);
    char *filepathE2 = "../tests/dat/out/example2_randmat.json";
    write_larcMatrix_file_by_matPTR(randmat, filepathE2);

    // Step 3: we count 0's in the random matrix and make sure it makes
    // sense. 
    // This mpz_t type will hold the count. 
    mpz_t zeroCount;
    mpz_init(zeroCount);
    // We need a scalarType version of 0, the scalar we'll be counting. 
    scalarType zero;
    sca_init(&zero);
    sca_set_str(&zero, "0");
    
    matrix_count_entries(zeroCount, randmat, zero);
    char *zero_string = sca_get_str(zero);
    gmp_printf("That random matrix we made has %Zd entries with value %s.\n", zeroCount, zero_string);
    free(zero_string);
    sca_clear(&zero);

    // Incidentally, if we add the number of zeros to the number of the 
    // sparse scalar, we should get the matrix size. 
    mpz_t zsum;
    mpz_init(zsum);
    mpz_add(zsum, zeroCount, scalarNum);
    if (0 == mpz_cmp_ui(zsum, (1L << row_level)*(1L << col_level)))
        printf("Test passed!\n");
    else
        printf("Test failed!\n");

    // We could locate where those sparse scalars were placed:
    // but we have to use the json file, not the matrix pointer. 
    int64_t **locations;
    mpz_locate_entries_larcMatrixFile(&locations, scalarNum, filepathE2, scalar, 100);
    char *scalar_string = sca_get_str(scalar);
    printf("The scalar %s can be found at the following coordinates:\n", scalar_string);
    free(scalar_string);
    for (int64_t i = 0; mpz_cmp_ui(scalarNum, i) > 0; i++) {
        printf("  [%ld, %ld]\n", locations[i][0], locations[i][1]);
        free(locations[i]);
    }
    free(locations);

    // Step 4: just for fun, we could hook the sca_get_str routine again
    // (see Example 1) and make it so when we write out the random matrix
    // again, the 0's are converted to 1's giving us a matrix of 1's and
    // 2's instead of 0's and 2's. 
    define_sca_get_str(zero_to_one_else_two);
    write_larcMatrix_file_by_matPTR(randmat, "../tests/dat/out/example2_newrandmat.json");
    // reset the scalar operations to normal
    

    // clean up and return operations to normal
    mpz_clear(zsum);
    mpz_clear(zeroCount);
    mpz_clear(scalarNum);
    sca_clear(&scalar);
    init_arithmetic_scalarOps(1);
}

