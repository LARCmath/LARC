#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

#include "larc.h"
#include "scalars.h"
#include "show_layout.h"


/*!
 * \file show_layout.c
 * \brief This file contains debugging routines that allow the user to print out
 *        the internal layout of the main LARC structures.
 *
 */

static void show_byte_array_header (size_t num_bytes) {
    int i;

    printf("  ");
    for (i = 0; i < num_bytes; ++i) {
        printf(" %2d", i);
    }
    printf("\n");

    printf("  ");
    for (i = 0; i < num_bytes; ++i) {
        printf(" --");
    }
    printf("\n");
}

static void reinit_layout_data (size_t num_bytes, unsigned char *all_zeros,
                         unsigned char *all_ones) {
    int i;

    for (i = 0; i < num_bytes; ++i) {
        all_zeros[i] = '\0';
        all_ones[i] = '\xff';
    }
}

static void calc_layout_answer (size_t num_bytes, unsigned char *all_zeros,
                         unsigned char *all_ones, unsigned char *answer) {
    int i;

    for (i = 0; i < num_bytes; ++i) {
        answer[i] = ((unsigned char) '\xff') ^ all_zeros[i] ^ all_ones[i];
    }
}

static void print_layout_answer (size_t num_bytes, unsigned char *answer,
                          const char *name) {
    int i;

    printf("  ");
    for (i = 0; i < num_bytes; ++i) {
        printf(" %02X", answer[i]);
    }
    printf("   %s\n", name);
}

#define HANDLE_FIELD(typ, fld) \
    reinit_layout_data(num_bytes, all_zeros, all_ones); \
    ((typ *) all_zeros)->fld = 0; \
    ((typ *) all_ones)->fld = 0; \
    calc_layout_answer(num_bytes, all_zeros, all_ones, answer); \
    print_layout_answer(num_bytes, answer, #fld); \

#define HANDLE_SCALAR_FIELD(typ, fld) \
    reinit_layout_data(num_bytes, all_zeros, all_ones); \
    sca_init(&(((typ *) all_zeros)->fld)); \
    sca_init(&(((typ *) all_ones)->fld)); \
    calc_layout_answer(num_bytes, all_zeros, all_ones, answer); \
    print_layout_answer(num_bytes, answer, #fld); \

void show_layout_hash_node () {
    size_t num_bytes = sizeof(hash_node_t);
    unsigned char *all_zeros = (unsigned char *) malloc(num_bytes);
    unsigned char *all_ones = (unsigned char *) malloc(num_bytes);
    unsigned char *answer = (unsigned char *) malloc(num_bytes);

    printf("\nSize of 'hash_node_t' = %zd bytes.\n\n", num_bytes);
    show_byte_array_header(num_bytes);

    HANDLE_FIELD(hash_node_t, record_ptr);
    HANDLE_FIELD(hash_node_t, record_hits);
    HANDLE_FIELD(hash_node_t, hashfilter);
    HANDLE_FIELD(hash_node_t, tile_offset_flag);
    HANDLE_FIELD(hash_node_t, hits_maxxed);
    HANDLE_FIELD(hash_node_t, xflags);
    HANDLE_FIELD(hash_node_t, next);
    HANDLE_FIELD(hash_node_t, prev);
}

void show_layout_larc_matrix () {
    size_t num_bytes = sizeof(larc_nsmatrix_t);
    unsigned char *all_zeros = (unsigned char *) malloc(num_bytes);
    unsigned char *all_ones = (unsigned char *) malloc(num_bytes);
    unsigned char *answer = (unsigned char *) malloc(num_bytes);

    printf("\nSize of 'larc_nsmatrix_t' = %zd bytes.\n\n", num_bytes);
    show_byte_array_header(num_bytes);

    HANDLE_FIELD(larc_nsmatrix_t, row_level);
    HANDLE_FIELD(larc_nsmatrix_t, col_level);
#ifdef MAR
#ifdef STORE_TILE_INDEX
    HANDLE_FIELD(larc_nsmatrix_t, tile);
#endif // #ifdef STORE_TILE_INDEX
#endif // #ifdef MAR
    HANDLE_FIELD(larc_nsmatrix_t, info);
    HANDLE_FIELD(larc_nsmatrix_t, hold);
    HANDLE_FIELD(larc_nsmatrix_t, lock);
    HANDLE_FIELD(larc_nsmatrix_t, iszero);
    HANDLE_FIELD(larc_nsmatrix_t, isid);
    HANDLE_FIELD(larc_nsmatrix_t, appears_as_sub_count);
    HANDLE_FIELD(larc_nsmatrix_t, packedID);
    HANDLE_SCALAR_FIELD(larc_nsmatrix_t, trace_element);
    HANDLE_SCALAR_FIELD(larc_smatrix_t, scalar_value);
    HANDLE_FIELD(larc_nsmatrix_t, subMatList[0])
    HANDLE_FIELD(larc_nsmatrix_t, subMatList[1]);
    HANDLE_FIELD(larc_nsmatrix_t, subMatList[2]);
    HANDLE_FIELD(larc_nsmatrix_t, subMatList[3]);
    HANDLE_FIELD(larc_smatrix_t, tile_node[0]);
    HANDLE_FIELD(larc_smatrix_t, tile_node[1]);
    HANDLE_FIELD(larc_smatrix_t, tile_node[2]);
    HANDLE_FIELD(larc_smatrix_t, tile_node[3]);
}


