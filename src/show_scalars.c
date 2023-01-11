#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpc.h>

#include "larc.h"
#include "show_scalars.h"

/*!
 * \file show_scalars.c
 * \brief This file contains debugging routines that allow the user to print out the internal details of the C, GMP, and MPC types used by LARC.
 *
 */

const char* get_string_gmp_versions () {
    static char gmp_buf[5000];
    gmp_buf[0] = '\0';

    char line_buf[1000];
    snprintf(line_buf, 1000, "    GMP version is: %d.%d.%d\n", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);
    strcat(gmp_buf, line_buf);
    snprintf(line_buf, 1000, "    MPFR version is: %s\n", MPFR_VERSION_STRING);
    strcat(gmp_buf, line_buf);
    snprintf(line_buf, 1000, "    MPC version is: %s\n", MPC_VERSION_STRING);
    strcat(gmp_buf, line_buf);

    return gmp_buf;
}

void print_gmp_versions () {
    printf("%s", get_string_gmp_versions());
}

void debug_out_mpz (const mpz_t v, const char * title) {
    char *indent = "";
    if (strncmp(title, "RatComplex -", 12) == 0) {
        indent = "   ";
    } else if (strncmp(title, "mpq -", 5) == 0) {
        indent = "   ";
    } else {
        printf("\n");
    }
    printf("%sTitle: %s, MPZ Addr: %p\n", indent, title, v);
    if (v == NULL) {
        printf("=> ERROR: Pointer to mpz structure is NULL.\n");
        return;
    }

    char *s10 = mpz_get_str(NULL, 10, v);
    char *s16 = mpz_get_str(NULL, 16, v);
    printf("%sBase10: %s\n", indent, s10);
    printf("%sBase16: %s\n", indent, s16);
    free(s10);
    free(s16);

    int num_of_limbs = abs(v->_mp_size);
    if (num_of_limbs > 0) {
        printf("%sLimbs[%d..%d]:  0x", indent, num_of_limbs-1, 0);
        printf("%016lx", v->_mp_d[num_of_limbs - 1]);
        int i;
        for (i = num_of_limbs - 2; i >= 0; --i) {
            printf("-%016lx", v->_mp_d[i]);
        }
        printf("\n");
    } else {
        printf("%sLimbs[-1..0]:  NO LIMBS USED (ZERO).\n", indent);
    }

    int total_alloc_limbs = v->_mp_alloc;
    if (num_of_limbs < total_alloc_limbs) {
        printf("%sExtra[%d..%d]:  0x", indent, total_alloc_limbs-1, num_of_limbs);
        printf("%016lx", v->_mp_d[total_alloc_limbs - 1]);
        int i;
        for (i = total_alloc_limbs - 2; i >= num_of_limbs; --i) {
            printf("-%016lx", v->_mp_d[i]);
        }
        printf("\n");
    }

    printf("%sMPZ Header: limbptr=%p, size=%ld, alloc=%ld\n",
        indent, v->_mp_d, (long) v->_mp_size, (long) v->_mp_alloc);
}

void debug_out_mpq (const mpq_t v, const char * title) {
    printf("\nTitle: %s, MPQ Addr: %p\n", title, v);
    if (v == NULL) {
        printf("=> ERROR: Pointer to mpq structure is NULL.\n");
        return;
    }
    debug_out_mpz(mpq_numref(v), "mpq - Numerator");
    printf("   --\n");
    debug_out_mpz(mpq_denref(v), "mpq - Denominator");
}

#ifdef USE_MPRATCOMPLEX
void debug_out_RatComplex (const larc_mpratcomplex_t v, const char * title) {
    printf("\nTitle: %s, MPRATCOMP Addr: %p\n", title, v);
    if (v == NULL) {
        printf("=> ERROR: Pointer to larc_mpratcomplex_s structure is NULL.\n");
        return;
    }
    debug_out_mpz(mpq_numref(v->real), "RatComplex - Numerator of Real Part");
    printf("   --\n");
    debug_out_mpz(mpq_denref(v->real), "RatComplex - Denominator of Real Part");
    printf("   --\n");
    debug_out_mpz(mpq_numref(v->imag), "RatComplex - Numerator of Imaginary Part");
    printf("   --\n");
    debug_out_mpz(mpq_denref(v->imag), "RatComplex - Denominator of Imaginary Part");
}
#endif // #ifdef USE_MPRATCOMPLEX

#ifdef USE_CLIFFORD
void debug_out_Clifford (const clifford_t v, const char * title) {
    printf("\nTitle: %s, CLIFFORD Addr: %p\n", title, v);
    if (v == NULL) {
        printf("=> ERROR: Pointer to clifford_s structure is NULL.\n");
        return;
    }
    int i;
    char sbuf[1024];
    for (i = 0; i < CLIFFORD_DIMENSION; ++i) {
        if (i != 0) {
            printf("   --\n");
        }
        (void) snprintf(sbuf, sizeof(sbuf), "Clifford - Numerator of Real Coeff # %d - %s",
            i, clifford_consts[i]);
        debug_out_mpz(mpq_numref(v->real_coeffs[i]), sbuf);
        printf("   --\n");
        (void) snprintf(sbuf, sizeof(sbuf), "Clifford - Denominator of Real Coeff # %d - %s",
            i, clifford_consts[i]);
        debug_out_mpz(mpq_denref(v->real_coeffs[i]), sbuf);
    }
#ifdef IS_COMPLEX
    for (i = 0; i < CLIFFORD_DIMENSION; ++i) {
        printf("   --\n");
        (void) snprintf(sbuf, sizeof(sbuf), "Clifford - Numerator of Imaginary Coeff # %d - %s",
            i, clifford_consts[i]);
        debug_out_mpz(mpq_numref(v->imag_coeffs[i]), sbuf);
        printf("   --\n");
        (void) snprintf(sbuf, sizeof(sbuf), "Clifford - Denominator of Imaginary Coeff # %d - %s",
            i, clifford_consts[i]);
        debug_out_mpz(mpq_denref(v->imag_coeffs[i]), sbuf);
    }
#endif // #ifdef IS_COMPLEX
}
#endif // #ifdef USE_CLIFFORD

void debug_out_mpfr (const mpfr_t v, const char * title) {
    char *indent = "";
    if (strncmp(title, "mpc -", 5) == 0) {
        indent = "   ";
    } else {
        printf("\n");
    }
    printf("%sTitle: %s, MPFR Addr: %p\n", indent, title, v);
    if (v == NULL) {
        printf("=> ERROR: Pointer to mpfr structure is NULL.\n");
        return;
    }

    mpfr_exp_t v10exp, v16exp;
    char *s10 = mpfr_get_str(NULL, &v10exp, 10, 0, v, MPFR_RNDN);
    char *s16 = mpfr_get_str(NULL, &v16exp, 16, 0, v, MPFR_RNDN);
    if (s10[0] == '-') {
        printf("%sBase10: -0.%s E%ld\n", indent, &s10[1], v10exp);
    } else {
        printf("%sBase10: 0.%s E%ld\n", indent, s10, v10exp);
    }
    if (s16[0] == '-') {
        printf("%sBase16: -0.%s E%ld\n", indent, &s16[1], v16exp);
    } else {
        printf("%sBase16: 0.%s E%ld\n", indent, s16, v16exp);
    }
    mpfr_free_str(s10);
    mpfr_free_str(s16);

    int num_of_limbs = ceil(((double) v->_mpfr_prec) / ((double) GMP_NUMB_BITS));
    printf("%sLimbs[%d..%d]:  0x", indent, num_of_limbs-1, 0);
    printf("%016lx", v->_mpfr_d[num_of_limbs - 1]);
    int i;
    for (i = num_of_limbs - 2; i >= 0; --i) {
        printf("-%016lx", v->_mpfr_d[i]);
    }
    printf("\n");

    printf("%sMPFR Header: limbptr=%p, prec=%ld, sign=%ld, exp=%ld",
        indent, v->_mpfr_d, (long) v->_mpfr_prec,
        (long) v->_mpfr_sign, (long) v->_mpfr_exp );
    if (v->_mpfr_exp == __MPFR_EXP_ZERO) {
        printf(" (ZERO)\n");
    } else if (v->_mpfr_exp == __MPFR_EXP_NAN) {
        printf(" (NaN)\n");
    } else if (v->_mpfr_exp == __MPFR_EXP_INF) {
        printf(" (INF)\n");
    } else {
        printf("\n");
    }
}

void debug_out_mpc (const mpc_t v, const char * title) {
    printf("\nTitle: %s, MPC Addr: %p\n", title, v);
    if (v == NULL) {
        printf("=> ERROR: Pointer to mpc structure is NULL.\n");
        return;
    }
    debug_out_mpfr(mpc_realref(v), "mpc - Real Part");
    printf("   --\n");
    debug_out_mpfr(mpc_imagref(v), "mpc - Imaginary Part");
}

void debug_out_integer (const int64_t v, const char * title) {
    int num_bytes = sizeof(v);
    printf("\nTitle: %s, Integer (size = %d bytes)\n", title, num_bytes);
    printf("Base10: %" PRId64 "\n", v);

    printf("Bytes[%d..%d]:  0x", num_bytes-1, 0);
    printf("%02hhx", ((char *) &v)[num_bytes - 1]);
    int i;
    for (i = num_bytes - 2; i >= 0; --i) {
        printf("-%02hhx", ((char *) &v)[i]);
    }
    printf("\n");
}

void debug_out_real (const long double v, const char * title) {
    int num_bytes = sizeof(v);
    printf("\nTitle: %s, Real (size = %d bytes)\n", title, num_bytes);
    printf("Value: %.20Lg\n", v);
    printf("Hex Value: %La\n", v);

    printf("Bytes[%d..%d]:  0x", num_bytes-1, 0);
    printf("%02hhx", ((char *) &v)[num_bytes - 1]);
    int i;
    for (i = num_bytes - 2; i >= 0; --i) {
        printf("-%02hhx", ((char *) &v)[i]);
    }
    printf("\n");
}

void debug_out_complex (const long double complex v, const char * title) {
    int num_bytes = sizeof(v);
    printf("\nTitle: %s, Complex (size = %d bytes)\n", title, num_bytes);
    printf("Value: %.20Lg +I* %.20Lg\n", creall(v), cimagl(v));
    printf("Hex Value: %La +I* %La\n", creall(v), cimagl(v));

    printf("Bytes[%d..%d]:  0x", (num_bytes/2)-1, 0);
    printf("%02hhx", ((char *) &v)[(num_bytes/2) - 1]);
    int i;
    for (i = (num_bytes/2) - 2; i >= 0; --i) {
        printf("-%02hhx", ((char *) &v)[i]);
    }
    printf("\n");

    printf("Bytes[%d..%d]:  0x", num_bytes-1, num_bytes/2);
    printf("%02hhx", ((char *) &v)[num_bytes - 1]);
    for (i = num_bytes - 2; i >= num_bytes/2; --i) {
        printf("-%02hhx", ((char *) &v)[i]);
    }
    printf("\n");
}

void debug_out_scalar (const scalarType v, const char * title) {

#ifdef USE_MPINTEGER
    debug_out_mpz(v, title);
#elif defined(USE_MPRATIONAL)
    debug_out_mpq(v, title);
#elif defined(USE_MPRATCOMPLEX)
    debug_out_RatComplex(v, title);
#elif defined(USE_CLIFFORD)
    debug_out_Clifford(v, title);
#elif defined(USE_MPREAL)
    debug_out_mpfr(v, title);
#elif defined(USE_MPCOMPLEX)
    debug_out_mpc(v, title);
#elif defined(USE_INTEGER) || defined(USE_BOOLEAN)
    debug_out_integer(v, title);
#elif defined(USE_REAL)
    debug_out_real(v, title);
#elif defined(USE_COMPLEX)
    debug_out_complex(v, title);
#else
    printf("\nTitle: %s, UNKNOWN SCALAR TYPE!!!\n", title);
#endif // #ifdef USE_MPINTEGER

}
