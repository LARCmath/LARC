//                    exampleLARC.c
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

// This version uses matrixIDs instead of pointers

static inline void print_usage (char *program)
{
      printf("Usage: %s [-h] [-r] [-v] [-t <int>] [-m <int>] [-o <dir>] [-s <int>] [-x <int>] [-X <int>] [-z <int>]\n", program);
      printf("Options:\n");
      printf("  -h           Display this usage information\n");
      printf("  -r           Print a report of matrix and op stores when program completes\n");
      printf("  -v           Print verbose information\n");
      printf("  -t <int>     Display reporting information every <int> seconds>\n");
      printf("  -m <int>     Specify maximum size of matrices is 2^<int> by 2^<int> (default %d)\n", DEFAULT_MAX_LEVEL);
      printf("  -o <dir>     Save output data in <dir>\n");
      printf("  -x <int>     Use <int> bits for hashing matrix_store (default %zd)\n", (size_t)DEFAULT_MATRIX_STORE_EXPONENT);
      printf("  -X <int>     Use <int> bits for hashing op_store (default %zd)\n", (size_t)DEFAULT_OP_STORE_EXPONENT);
      printf("  -s <int>     Locality-approximation parameter rounds to <int> significant bits (default %d)\n", SIGHASH_DEFAULT);
      printf("  -z <int>     Locality-approximation parameter: value within <int> bits of zero equates to zero  (default %d)\n", ZEROBITTHRESH_DEFAULT);
      printf("\n");
      exit(1);
}

static inline double wallsec (void)
{
      struct timeval current_time;
      gettimeofday(&current_time, NULL);
      return current_time.tv_sec + 1e-6*current_time.tv_usec;
}


int main (int argc, char *argv[])
{
      // Default hash exponent size
      size_t hash_exponent_matrix = DEFAULT_MATRIX_STORE_EXPONENT;
      size_t hash_exponent_op = DEFAULT_OP_STORE_EXPONENT;
      int zerobitthresh = ZEROBITTHRESH_DEFAULT;
      int sighash = SIGHASH_DEFAULT;
      mat_level_t max_level = DEFAULT_MAX_LEVEL;

      // Command line parameters
      char *output_dir = NULL;
      //      char *approxtype = NULL, *value = NULL;
      int opt = 0;
      int final_report = 0;   // report on table sizes and usage at end of run
      int verbose = 1;
      unsigned int report_period = 0;     // periodic reporting

      // Completion statistics
      double start_time, stop_time;


#define COMMAND_LINE
#ifdef COMMAND_LINE
/*       // Parse command line options */
       while ((opt = getopt(argc, argv, "hrvt:m:o:s:X:x:z:")) != -1) { 
             switch (opt) { 
             case 'r':
                   final_report = 1;
                   break;
             case 'v':
                   verbose++;
                   break;
             case 't':
                   report_period = atoi(optarg);
                   break;
             case 'm':
                   max_level = atoi(optarg);
                   break;
	     case 'o':
                   if (output_dir != NULL) {
                         printf("Error: You can only specify one output directory\n");
                         print_usage(argv[0]);
                   }
                   output_dir = optarg;
                   break;
             case 's':
                   sighash = atoi(optarg);
                   break;
             case 'X':
                   hash_exponent_op = atoi(optarg);
                   break;
             case 'x':
                   hash_exponent_matrix = atoi(optarg);
                   break;
             case 'z':
                   zerobitthresh = atoi(optarg);
                   break;
             default:
                   print_usage(argv[0]);
                   break;
             }
       }

       if (verbose) {
           printf ("\n==================================================================\n");
             printf ("|  The hash table exponent for matrix store is %zd and the hash table size is %zd\n", 
                hash_exponent_matrix, (size_t)1<<hash_exponent_matrix);
             printf ("|  The hash table exponent for each op store is %zd and the hash table size is %zd\n", 
           hash_exponent_op, (size_t)1<<hash_exponent_op);
	     printf( " The value of sighash is %d\n",sighash);
	     printf( " The value of zerobitthresh is %d\n",zerobitthresh);
	     printf( " The value of max_level is %d\n",max_level);
             printf ("==================================================================\n\n");
       }
#endif

      // if -t option sets report_period then a report on tables and 
      // memory given every report_period seconds
      if (report_period) create_report_thread(report_period);

      struct stat st = {0};
      if (stat("./out", &st) == -1) {
	mkdir("./out", 0700);
      }

      // Initialize the matrix store and op stores
      // Preloads basic matrices: Zeros, Identities, integer Hadamards, and other common 2x2 and 4x4 matrices
      // Start the seppuku thread to kill program if it is using too much memory
      initialize_larc(hash_exponent_matrix,hash_exponent_op,max_level,sighash,zerobitthresh);

      start_time = wallsec();
       
      printf("The matrix store and operation stores have been created.\n");
      printf("Now we will preload the matrix store with all-zero, identity, and other useful matrices.\n");

      // Calculate number of matrices created, then print part of matrix store
      int num_matrices_made = num_matrices_created();
      printf("%d matrices have been created\n",num_matrices_made);
      int last_matrix_index = num_matrices_made - 1;
      
      char comment[1024];
      sprintf(comment,"Matrix store after preload; parameters=%zd,%zd,%d,%d,%d.\n",
		    (size_t)hash_exponent_matrix,(size_t)hash_exponent_op,max_level,sighash,zerobitthresh);
      char file_name[1024] = "./out/matrix_store_after_preload";
      matrix_store_info_to_file(0,last_matrix_index,file_name,comment);
      

printf("Now we give some examples of:\n");
      printf("  * how to use input/output,\n");
      printf("  * and use math functions.\n");
#ifndef USE_INTEGER
      printf("  * We also illustrate the effect of the locality-approximation parameters,\n");
      printf("    which effect numerical precision issues.\n\n");
#endif 
      printf("First we show how to build a matrix, read it into the larc store\n");
      printf("and store it in a file or print it to the screen in various formats.\n");


      /*****************************
       * EXAMPLES OF I/O FUNCTIONS *
       *****************************/


      //We support integer, real, and complex numbers.
#ifndef USE_COMPLEX
      ScalarType scalar7 = 7;
#else 
      //If complex is used, be sure to add 0*I for values without an imaginary component.
      ScalarType scalar7 = 7 + 0*I;
#endif
   
      //storing scalar 7 
      //the last two arguments are the level of the matrix: the matrix is of size 2^level x 2^level 
      int64_t mID_scalar7 = matrix_get_matrixID_from_scalar(scalar7);
      
      //If you anticipate using a value many times in the course of a computation,
      //you can "hold" it in the matrix store and release it later
      set_hold_matrix_from_matrixID(mID_scalar7);
      //If you would like to re-use a particular matrix, you can store the name using extern. 
      //See global.c and global.h for an example of how to do this.

      printf("We have stored scalar 7 and are printing it to the screen.\n");
      print_matrix_naive_by_matrixID(mID_scalar7);

       //Now we construct a square 2x2 matrix of all 7's.
      //The 0th entry in the submatrix corresponds to the upper left quadrant of the matrix,
      //the 1st to the upper right, the 2nd to the lower left, and the 3rd to the lower right.
      int64_t panel_square7[4];
      panel_square7[0] = mID_scalar7;
      panel_square7[1] = mID_scalar7;
      panel_square7[2] = mID_scalar7;
      panel_square7[3] = mID_scalar7;
      //the matrix is of size 2^1 x 2^1, so both the row and column levels are 1
      int64_t mID_square7 = matrix_get_matrixID_from_panel(mID_scalar7,mID_scalar7,mID_scalar7,mID_scalar7, 1, 1);
 
      printf("We now write square 7 to screen and various file formats in the out directory.\n");
      print_matrix_naive_by_matrixID(mID_square7);
      //writing to a JSON file
      matrix_write_json_file_matrixID(mID_square7, "./out/square7.json");
      //writing to a naive format
      print_matrix_to_file_naive_by_matrixID(mID_square7,"./out/square7.naive");
     

      //NB: THIS IS COMMENTED OUT TO RESOLVE AN UNUSED VAR WARNING.
      //writing a 2x2 matrix of 7's and 0's in row-major format to the store
      //ScalarType rmm[4] = {7,0,0,7};
      //here, we specify the row and the column levels, as well as the length of the row major matrix
      //mat_add_t mID_rmm = row_major_list_to_store_matrixID(rmm, 2, 2, 4);
      
      //Now we construct a 2x1 matrix of all 7's.
      int64_t mID_rect7 = stack_matrixID(mID_scalar7, mID_scalar7);
      printf("We are printing 2x1 matrix of all 7's to the screen.\n");
      print_matrix_naive_by_matrixID(mID_rect7);

      //Now we print the part of the matrix store that has been created since the preload.
      num_matrices_made = num_matrices_created();
      printf("\n%d matrices have been created.\n", num_matrices_made);
      int start = last_matrix_index + 1;
      int end = num_matrices_made - 1;
      matrix_store_info_to_file(start,end,"./out/ALLsevens.store","Matrices added after preload");


      // get hash value for a matrix we are about to remove
      int64_t hashID = matrix_hashID_from_matrixID(mID_rect7);

      // print the hash chain
      matrix_hash_chain_info_to_file(hashID, "./out/hashChain1", "hash chain before removal");
      matrix_hash_chain_info_to_screen(hashID, "hash chain before removal");

      //We can also delete matrices (that are not held or locked) from the store.
      printf("\nDeleting the column of 7s from the store\n");
      remove_matrix_from_mat_store_by_matrixID(mID_rect7);
      matrix_store_info_to_file(start,end,"./out/SOMEsevens.store","Removed column of sevens");

      // print the hash chain
      matrix_hash_chain_info_to_file(hashID, "./out/hashChain2", "hash chain after removal");
	  matrix_hash_chain_info_to_screen(hashID, "hash chain after removal");

      /******************************
       * EXAMPLES OF MATH FUNCTIONS *
       ******************************/
      //syntax is fairly self-explanatory: operation(mID_matrix_A, mID_matrix_B)
      //NB: adjoint only has one argument, since it is not a binary operation
      
      ScalarType scalarM1 = -1;
#ifndef USE_INTEGER
      ScalarType scalar0 = 0;
//      int64_t mID_scalar0 = matrix_get_matrixID_from_scalar(scalar0);
#endif
#ifdef USE_COMPLEX
      ScalarType scalar0i1 = 1*I;
      int64_t mID_scalar0i1 = matrix_get_matrixID_from_scalar(scalar0i1);
#endif
      int64_t mID_scalarM1 = matrix_get_matrixID_from_scalar(scalarM1);

      printf("testing scalar_mult:\n");
#ifdef USE_COMPLEX
      int64_t samp_mID = scalar_mult_matrixID(mID_scalar0i1,mID_square7);
#else
      int64_t samp_mID = scalar_mult_matrixID(mID_scalarM1,mID_square7);
#endif
      print_matrix_naive_by_matrixID(samp_mID);
      printf("testing addition:\n");
      int64_t samp1_mID = matrix_add_matrixID(samp_mID,samp_mID);
      print_matrix_naive_by_matrixID(samp1_mID);
      printf("testing adjoint:\n");
      int64_t samp2_mID = matrix_adjoint_matrixID(samp_mID);
      print_matrix_naive_by_matrixID(samp2_mID);
      printf("testing matrix mult:\n");
      int64_t samp3_mID = matrix_mult_matrixID(samp2_mID,samp_mID);
      print_matrix_naive_by_matrixID(samp3_mID);
      
      // test printing op store hash chain
      char *op_name = "PRODUCT";
      hashID = op_hashID_by_matrixIDs(samp2_mID, samp_mID, op_name);
      printf("About to test ability to print op hash chains to file\n");
      op_hash_chain_info_to_screen(hashID, "print hash chain including matrix multiplication");
      
      printf("testing kron product:\n");
      int64_t samp4_mID = kronecker_product_matrixID(samp_mID,samp_mID);
      print_matrix_naive_by_matrixID(samp4_mID);
      printf("testing join:\n");
      int64_t samp5_mID = join_matrixID(samp_mID,samp_mID);
      print_matrix_naive_by_matrixID(samp5_mID);
      printf("testing stack:\n");
      int64_t samp6_mID = stack_matrixID(samp_mID,samp_mID);
      print_matrix_naive_by_matrixID(samp6_mID);

#ifndef USE_INTEGER
      /*******************************
       * EXAMPLES OF LOCALITY-APPROXIMATION EFFECT *
       *******************************/

      //As is, all of the tests will pass in this section of the code.
      //However, if you modify the parameters zerobitthresh and 
      //sighash, you may start to see performance issues or numerical
      //precision issues if you set these values to high or low values
      //respectively.

      // Make matrices to test locality-approximation
      // two by two tests
      int level = 1;
      int dim_whole = pow(2.0,(double)level);

      ScalarType arr_a[4] = {.4, 0, 0, .3};
      int64_t a_mID = row_major_list_to_store_matrixID(arr_a,level,level,dim_whole);
      ScalarType arr_b[4] = {.8, 0, 0, .6};
      int64_t b_mID = row_major_list_to_store_matrixID(arr_b,level,level,dim_whole);
      ScalarType arr_c[4] = {-.4, 0, 0, -.3};
      int64_t c_mID = row_major_list_to_store_matrixID(arr_c,level,level,dim_whole); 
      // printf("The matrixIDs of a, b, and c are %zd %zd %zd\n", 
      printf("The matrixIDs of a, b, and c are %" PRIu64 " %" PRIu64 " %" PRIu64 " \n",
		(uint64_t)a_mID, (uint64_t)b_mID, (uint64_t)c_mID);
     
      int64_t d_mID = matrix_add_matrixID(a_mID,a_mID);
      printf("MatrixID of a + a %" PRIu64 ", should be that of b %" PRIu64 " \n", (uint64_t)d_mID, (uint64_t)b_mID);
      if (d_mID == b_mID){
	printf("PASSED: [.4,0,0,.3] + [.4,0,0,3.] = [.8,0,0,.6]\n");
      }
      else { 
	printf("FAILED: [.4,0,0,.3] + [.4,0,0,3.] = [.8,0,0,.6]\n");
      }
      
      int64_t e_mID = matrix_add_matrixID(b_mID,c_mID);
      printf("MatrixID of b + c %" PRIu64 ", should be that of a %" PRIu64 " \n", (uint64_t)e_mID, (uint64_t)a_mID);
      if (e_mID == a_mID) { 
	printf("PASSED: [.8,0,0,.6] + [-.4,0,0,-.3] = [.4,0,0,3.]\n");
      }
      else { 
	printf("FAILED: [.8,0,0,.6] + [-.4,0,0,-.3] = [.4,0,0,3.]\n");
      }
    
    
      // scalar tests
      level = 0;
      dim_whole = pow(2.0,(double)level);
      
      ScalarType arr_f[1] = {.4};
      int64_t f_mID = row_major_list_to_store_matrixID(arr_f,level,level,dim_whole);
      //printf("We input the value into f (%zd) of .4, the matrix is:\n",(uint64_t)f_mID);
      // print_matrix_naive_by_matrixID(f_mID);

      ScalarType arr_g[1] = {.4+1e-5};
      int64_t g_mID = row_major_list_to_store_matrixID(arr_g,level,level,dim_whole);
      // printf("We input the value into g (%zd) of .4+1e-5, the matrix is:\n" ,(uint64_t)g_mID);
      // print_matrix_naive_by_matrixID(g_mID);

      ScalarType arr_h[1] = {.4-1e-10};
      int64_t h_mID = row_major_list_to_store_matrixID(arr_h,level,level,dim_whole);
      //printf("We input the value into h (%zd) of .4-1e-10, the matrix is:\n", (uint64_t)h_mID);
      //print_matrix_naive_by_matrixID(h_mID);

      int64_t i_mID = matrix_mult_matrixID(h_mID, mID_scalarM1);
      printf("The calculated negative of h (%" PRIu64 ") is i (%" PRIu64 ") with matrix\n" ,(uint64_t)h_mID,(uint64_t)i_mID);
      print_matrix_naive_by_matrixID(i_mID);

      int64_t z_mID = matrix_get_matrixID_from_scalar(scalar0);
      printf("The prestored value for zero is z (%" PRIu64 ") with matrix\n" , (uint64_t)z_mID);
      print_matrix_naive_by_matrixID(z_mID);

      int64_t j_mID = matrix_add_matrixID(f_mID,i_mID);
      printf("We calculate f + i = j (%" PRIu64 "), with matrix\n" ,(uint64_t)j_mID);
      print_matrix_naive_by_matrixID(j_mID);

      if (j_mID == z_mID) {
	printf("PASSED: .4 + (-1)*(.4) = 0\n");
      }
      else {
	printf("FAILED: .4 + (-1)*(.4) = 0\n");
      }

      int64_t k_mID = matrix_mult_matrixID(g_mID,mID_scalarM1);
      printf("The calculated negative of g (%" PRIu64 ") is k (%" PRIu64 ") with matrix\n" ,(uint64_t)g_mID,(uint64_t)k_mID);
      print_matrix_naive_by_matrixID(k_mID);

      int64_t l_mID = matrix_add_matrixID(f_mID,k_mID);
      printf("We calculate f + k = l (%" PRIu64 "), with matrix\n" ,(uint64_t)l_mID);
      print_matrix_naive_by_matrixID(l_mID);

      if (g_mID == l_mID){
	printf("ZERO THRESHOLD at least 1/10^5: SO (-1)*(.4 + 1e-5) + (.4) = 0\n");
      }
      else { 
	printf("ZERO THRESHOLD greater than 1/10^5: SO (-1)*(.4 + 1e-5) + (.4) != 0\n");
      }

      ScalarType arr_w[1] = {1e-20};
      int64_t w_mID = row_major_list_to_store_matrixID(arr_w,level,level,dim_whole);
      if (w_mID == z_mID) {
	printf("APPROXIMATION: 1e-20 = 0\n");
      }
      else { 
	printf("APPROXIMATION: 1e-20 != 0\n");
      }

      ScalarType arr_y[1] = {1e-10};
      int64_t y_mID = row_major_list_to_store_matrixID(arr_y,level,level,dim_whole);
      if (y_mID == z_mID){ 
	printf("APPROXIMATION: 1e-10 = 0\n");
      }
      else {
	printf("APPROXIMATION: 1e-10 != 0\n");
      }

      ScalarType arr_v[1] = {7};
      int64_t v_mID = row_major_list_to_store_matrixID(arr_v,level,level,dim_whole);
      ScalarType arr_x[1] = {7+1e-20};
      int64_t x_mID = row_major_list_to_store_matrixID(arr_x,level,level,dim_whole);
      int64_t cc_mID = matrix_add_matrixID(w_mID,v_mID);
      if (x_mID == cc_mID) { 
	printf("C vs LARC local-approx addition: larc(7)+larc(1e-20) = larc(7+1e-20)\n");
      }
      else {
	printf("C vs LARC local-approx addition: larc(7)+larc(1e-20) != larc(7+1e-20)\n");
      }

      ScalarType arr_aa[1] = {7+1e-10};
      int64_t aa_mID = row_major_list_to_store_matrixID(arr_aa,level,level,dim_whole);
      int64_t bb_mID = matrix_add_matrixID(y_mID,v_mID);
      if (bb_mID == aa_mID) { 
	printf("C vs LARC local-approx addition: larc(7)+larc(1e-10) = larc(7+1e-10)\n");
      }
      else { 
	printf("C vs LARC local-approx addition: larc(7)+larc(1e-10) != larc(7+1e-10)\n");
      }
      
      level = 0;
      dim_whole = pow(2.0,(double)level);
      ScalarType arr_m[1] = {sqrt(2)};
      int64_t m_mID = row_major_list_to_store_matrixID(arr_m,level,level,dim_whole);
      printf("We input the value into m (%" PRIu64 ") of sqrt(2), the matrix is:\n" ,(uint64_t)m_mID);
      print_matrix_naive_by_matrixID(m_mID);

      ScalarType arr_n[1] = {2};
      int64_t n_mID = row_major_list_to_store_matrixID(arr_n,level,level,dim_whole);
      printf("We input the value into n (%" PRIu64 ") of 2, the matrix is:\n" ,(uint64_t)n_mID);
      print_matrix_naive_by_matrixID(n_mID);
	  printf("and the matrixID is% " PRId64 ":\n" ,n_mID);

      int64_t p_mID = matrix_mult_matrixID(m_mID,m_mID);
      if (p_mID == n_mID){ 
	printf("PASSED: sqrt(2) squared is the same as 2\n");
      }
      else { 
	printf("FAILED: sqrt(2) squared is the same as 2\n");
      }

      int64_t q_mID = matrix_add_matrixID(m_mID,m_mID);
      ScalarType arr_r[1] = {0.5};
      int64_t r_mID = row_major_list_to_store_matrixID(arr_r,level,level,dim_whole);
      int64_t s_mID = matrix_mult_matrixID(q_mID,r_mID);
      if (m_mID == s_mID) { 
	printf("PASSED: (sqrt2+sqrt2)/2 is the same as sqrt2\n");
      }
      else {
	printf("FAILED: (sqrt2+sqrt2)/2 is the same as sqrt2\n");
      }

      ScalarType arr_t[1] = {1.0/sqrt(2)};
      int64_t t_mID = row_major_list_to_store_matrixID(arr_t,level,level,dim_whole);
      int64_t u_mID = matrix_mult_matrixID(t_mID,t_mID);
      if (u_mID == r_mID) { 
	printf("PASSED Hadamard test1: ((1/sqrt2)^2 is the same as .5\n");
      }
      else { 
	printf("FAILED Hadamard test1: ((1/sqrt2)^2 is the same as .5\n");
      }

      matrix_store_info_to_file(0,num_matrices_created()-1,"./out/local_approx.store","finished locality approximation testing");
#endif
	

      stop_time = wallsec();
      printf("\nElapsed time = %.4g secs\n", stop_time - start_time);

      if (final_report || verbose) {
            matrix_store_report("stdout");
            op_store_report("stdout");
            rusage_report(0,"stdout");
      }

      return 0;
}
