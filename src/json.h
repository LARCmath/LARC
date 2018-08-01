//                         json.h
/*
 *
 * json.h
 * POC: Michael Schneider <maschn2@super.org>
 *
 * Copyright 2012, Institute for Defense Analyses
 * 4850 Mark Center Drive, Alexandria, VA; 703-845-2500
 * This material may be reproduced by or for the US Government
 * pursuant to the copyright license under the clauses at DFARS
 * 252.227-7013 and 252.227-7014.
 *
 */

#ifndef JSON_H
#define JSON_H

#include <gmp.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <stdarg.h>
#include <stdlib.h>
#include <setjmp.h>

/*!
 * \cond
 */

////////////////////////////////////////////////////////////////////////////////
// structures
////////////////////////////////////////////////////////////////////////////////

typedef struct json {
    enum {J_STRING, J_NUM64, J_BIGNUM, J_BIGNUM_DOUBLE, J_DOUBLE, J_OBJECT, J_ARRAY, J_TRUE, J_FALSE, J_NULL} type;
    struct json *parent;
    struct json *list, *list_end, *list_alloc_end;
    char *name; // only for struct in list
    char *string, *string_end, *string_alloc_end;
    int64_t num64; // for objects and arrays, the number of structs allocated in "list"
    double num_double;
    int bignum_initialized;
    mpz_t bignum;

    char *orig_filename;
    int longjmp_set;
    jmp_buf longjmp_env;
    char *longjmp_string;
} json_t;

////////////////////////////////////////////////////////////////////////////////
// prototypes for static inline functions
////////////////////////////////////////////////////////////////////////////////

// type inquiry

static inline int j_is_string(json_t *j);
static inline int j_is_num64(json_t *j);
static inline int j_is_bignum(json_t *j);
static inline int j_is_double(json_t *j);
static inline int j_is_object(json_t *j);
static inline int j_is_array(json_t *j);
static inline int j_is_boolean(json_t *j);
static inline int j_is_true(json_t *j);
static inline int j_is_false(json_t *j);
static inline int j_is_null(json_t *j);

// gets

static inline const char *j_get_string(json_t *j);
static inline int64_t j_get_num64(json_t *j);
static inline mpz_t *j_get_bignum(json_t *j);
static inline double j_get_double(json_t *j);
static inline int j_get_bool(json_t *j);
static inline int j_get_true(json_t *j);
static inline int j_get_false(json_t *j);
static inline void j_get_null(json_t *j);
static inline const char *j_get_name(json_t *j);

// sets

static inline void j_set_string(json_t *j, const char *format, ...) __attribute__ ((format (printf, 2, 3)));
static inline void j_set_num64(json_t *j, int64_t num);
static inline void j_set_bignum(json_t *j, const mpz_t num);
static inline void j_set_double(json_t *j, double num);
static inline json_t *j_set_object(json_t *j);
static inline json_t *j_set_array(json_t *j);
static inline void j_set_boolean(json_t *j, int bool);
static inline void j_set_true(json_t *j);
static inline void j_set_false(json_t *j);
static inline void j_set_null(json_t *j);

// objects

static inline json_t *j_key_lookup(json_t *j, const char *name);
static inline json_t *j_key_add(json_t *j, const char *name);
static inline int     j_key_exists(json_t *j, const char *name);
static inline int64_t j_key_count(json_t *j);
static inline json_t *j_key_index(json_t *j, int64_t index);

// Convenience to add a key and a value.
static inline void j_add_string(json_t *j, const char *name, const char *format, ...) __attribute__ ((format (printf, 3, 4)));
static inline void j_add_num64(json_t *j, const char *name, int64_t num);
static inline void j_add_bignum(json_t *j, const char *name, const mpz_t num);
static inline void j_add_double(json_t *j, const char *name, double num);
static inline json_t *j_add_object(json_t *j, const char *name);
static inline json_t *j_add_array(json_t *j, const char *name);
static inline void j_add_boolean(json_t *j, const char *name, int bool);

// arrays

static inline json_t *j_array_index(json_t *j, int64_t index);
static inline json_t *j_array_append(json_t *j);
static inline void    j_array_set_length(json_t *j, int64_t length);
static inline int64_t j_array_get_length(json_t *j);

////////////////////////////////////////////////////////////////////////////////
// standard prototypes
////////////////////////////////////////////////////////////////////////////////

json_t *j_new(json_t *parent);

json_t *j_parse_string(const char *c);
json_t *j_parse_file(FILE *fp);
json_t *j_parse_filename(const char *filename);
json_t *j_parse_filename_into_null_object(const char *filename, json_t *root);

char *j_make_rel_filename(const json_t *j, const char *filename); // must free() non-NULL result

int j_snprint(char *dst, size_t n, json_t *j);

void j_throw(json_t *j, const char *format, ...) __attribute__ ((format (printf, 2, 3)));

////////////////////////////////////////////////////////////////////////////////
// error functions
////////////////////////////////////////////////////////////////////////////////

#define j_error(...) { \
      printf("j_error() at %s:%d in %s() ", __FILE__, __LINE__, __func__); \
      printf(__VA_ARGS__);                                                 \
      printf("\n");                                                        \
      abort(); }

#define j_check(error_if_false, ...) if (!(error_if_false)) { \
      printf("j_check(%s) failed at %s:%d in %s()", #error_if_false, __FILE__, __LINE__, __func__); \
      printf(__VA_ARGS__);                                                                          \
      printf("\n");                                                                                 \
      abort(); }

#define SPRINTF(out, arg) { \
        va_list a1, a2; \
        va_start(a1, arg); \
        /*va_copy(a2, a1);*/ \
        va_start(a2, arg); \
        int s1 = vsnprintf(NULL, 0, arg, a1); \
        va_end(a1); \
        j_check(s1 >= 0, "vsnprintf output error"); \
        out = malloc(s1+1); \
        j_check(out, "malloc(%u) failed",s1+1); \
        int s2 = vsnprintf(out, s1+1, arg, a2); \
        j_check(s2 == s1,"string length mismatch"); \
        va_end(a2); }

////////////////////////////////////////////////////////////////////////////////
// internal prototypes (only for use in this file)
////////////////////////////////////////////////////////////////////////////////

json_t *_j_parse_string(const char *c, json_t *root);
json_t *_j_parse_file(FILE *fp, json_t *root);
json_t *_j_parse_filename(const char *filename, json_t *root);

////////////////////////////////////////////////////////////////////////////////
// implementation for static inline functions
////////////////////////////////////////////////////////////////////////////////

// type inquiry

static inline int j_is_string(json_t *j) { return j && j->type == J_STRING; }
static inline int j_is_num64(json_t *j)  { return j && j->type == J_NUM64; }
static inline int j_is_bignum(json_t *j) { return j && ((j->type == J_NUM64) || (j->type == J_BIGNUM) || (j->type == J_BIGNUM_DOUBLE)); }
static inline int j_is_double(json_t *j) { return j && ((j->type == J_NUM64) || (j->type == J_DOUBLE) || (j->type == J_BIGNUM_DOUBLE)); }
static inline int j_is_object(json_t *j) { return j && j->type == J_OBJECT; }
static inline int j_is_array(json_t *j)  { return j && j->type == J_ARRAY; }
static inline int j_is_boolean(json_t *j){ return j && ((j->type == J_TRUE) || (j->type == J_FALSE)); }
static inline int j_is_true(json_t *j)   { return j && j->type == J_TRUE; }
static inline int j_is_false(json_t *j)  { return j && j->type == J_FALSE; }
static inline int j_is_null(json_t *j)   { return !j || j->type == J_NULL; }

// gets

static inline const char *j_get_string(json_t *j) {
    if (!j_is_string(j)) {
        j_throw(j, "expected string");
        return NULL;
    }
    return j->string;
}

static inline int64_t j_get_num64(json_t *j) {
    if (!j_is_num64(j)) {
        j_throw(j, "expected num64");
        return 0;
    }
    return j->num64;
}

static inline mpz_t *j_get_bignum(json_t *j) {
    if ((j->type == J_NUM64) && !j->bignum_initialized) {
        mpz_init(j->bignum);
        mpz_set_ui(j->bignum, j->num64);
        j->bignum_initialized = 1;
    }
    if (!j_is_bignum(j)) {
        j_throw(j, "expected bignum");
        return NULL;
    }
    return &j->bignum;
}

static inline double j_get_double(json_t *j) {
    if (!j_is_double(j)) {
        j_throw(j, "expected double");
        return 0.0;
    }
    return j->num_double;
}

static inline int j_get_bool(json_t *j) {
    if (j_is_true(j)) return 1;
    if (j_is_false(j)) return 0;
    j_throw(j, "expected true or false");
    return 0;
}

static inline int j_get_true(json_t *j) {
    if (j_is_true(j)) return 1;
    j_throw(j, "expected true");
    return 0;
}

static inline int j_get_false(json_t *j) {
    if (j_is_false(j)) return 1;
    j_throw(j, "expected false");
    return 0;
}

static inline void j_get_null(json_t *j) {
    if (j_is_null(j)) return;
    j_throw(j, "expected null");
}

static inline const char *j_get_name(json_t *j) {
    if (!j) return NULL;
    return j->name;
}

// sets

static inline void j_set_string(json_t *j, const char *format, ...) {
    j_set_null(j);
    j->type = J_STRING;

    SPRINTF(j->string, format);
}

static inline void j_set_num64(json_t *j, int64_t num) {
    j_set_null(j);
    j->type = J_NUM64;
    j->num64 = num;
}

static inline void j_set_bignum(json_t *j, const mpz_t num) {
    j_set_null(j);
    j->type = J_BIGNUM;
    mpz_init(j->bignum);
    j->bignum_initialized = 1;
    mpz_set(j->bignum, num);
}

static inline void j_set_double(json_t *j, double num) {
    j_set_null(j);
    j->type = J_DOUBLE;
    j->num_double = num;
}

static inline void j_set_boolean(json_t *j, int boolean) {
    if (boolean) j_set_true(j); else j_set_false(j);
}

static inline json_t *j_set_object(json_t *j) {
    j_set_null(j);
    j->type = J_OBJECT;
    return j;
}

static inline json_t *j_set_array(json_t *j) {
    j_set_null(j);
    j->type = J_ARRAY;
    return j;
}

static inline void j_set_true(json_t *j) {
    j_set_null(j);
    j->type = J_TRUE;
}

static inline void j_set_false(json_t *j) {
    j_set_null(j);
    j->type = J_FALSE;
}

static inline void j_set_null(json_t *j) {
    struct json *i;
    if (j->list) {
      for (i=j->list; i<j->list_end; i++) {
	j_set_null(i);
	if (i->name) { free(i->name); i->name = NULL;}
      }
      free(j->list);
      j->list = j->list_end = j->list_alloc_end = NULL;
    }

    if (j->string) { free(j->string); j->string = j->string_end = j->string_alloc_end = NULL; }
    if (j->bignum_initialized) { mpz_clear(j->bignum); j->bignum_initialized = 0; }
    if (j->orig_filename) { free(j->orig_filename); j->orig_filename = NULL; }
    j->type = J_NULL;
}

// objects

static inline json_t *_j_lookup(json_t *j, const char *name, int create) {
    if (j_is_object(j)) {
	json_t *p;
        for (p = j->list; p < j->list_end; p++) {
            if (!strcmp(name, p->name)) return p;
        }
        if (create) {p = j_new(j); p->name = strdup(name); return p;}
        j_throw(j, "field '%s' not found", name);
    }
    j_throw(j, "expected object");
    return NULL;
}

static inline json_t *j_key_lookup(json_t *j, const char *name) {
    return _j_lookup(j, name, 0);
}

static inline json_t *j_key_add(json_t *j, const char *name) {
    return _j_lookup(j, name, 1);
}

static inline int j_key_exists(json_t *j, const char *name) {
    if (!j_is_object(j)) {
        j_throw(j, "expected object");
        return 0;
    }
    json_t *p;
    for (p = j->list; p < j->list_end; p++) {
        if (!strcmp(name, p->name)) return 1;
    }
    return 0;
}

static inline int64_t j_key_count(json_t *j) {
    if (!j_is_object(j)) {
        j_throw(j, "expected object");
        return 0;
    }
    return j->list_end - j->list;
}

static inline json_t *j_key_index(json_t *j, int64_t index) {
    if (!j_is_object(j)) {
        j_throw(j, "expected object");
        return NULL;
    }
    if (index >= j->list_end - j->list) {
        j_throw(j, "element %ld out of bounds [0,%ld)",index,j->num64);
        return NULL;
    }
    return &j->list[index];
}

// Convenience to add a key and a value.
static inline void j_add_string(json_t *j, const char *name, const char *format, ...) {
  json_t *j2 = j_key_add(j, name);

  j2->type = J_STRING;
  SPRINTF(j2->string, format);
}

static inline void j_add_num64(json_t *j, const char *name, int64_t num) {
  j_set_num64(j_key_add(j, name), num);
}

static inline void j_add_bignum(json_t *j, const char *name, const mpz_t num) {
  j_set_bignum(j_key_add(j, name), num);
}

static inline void j_add_double(json_t *j, const char *name, double num) {
  j_set_double(j_key_add(j, name), num);
}

static inline json_t *j_add_object(json_t *j, const char *name) {
  json_t *j_sub = j_key_add(j, name);
  return j_set_object(j_sub);
}

static inline json_t *j_add_array(json_t *j, const char *name) {
  json_t *j_sub = j_key_add(j, name);
  return j_set_array(j_sub);
}

static inline void j_add_boolean(json_t *j, const char *name, int bool) {
  j_set_boolean(j_key_add(j, name), bool);
}


// arrays

static inline json_t *j_array_index(json_t *j, int64_t index) {
    if (!j_is_array(j)) {
        j_throw(j, "expected array");
        return NULL;
    }
    if (index >= j->list_end - j->list) {
        j_throw(j, "element %ld out of array bounds [0,%ld)",index,j->num64);
        return NULL;
    }
    return &j->list[index];
}

static inline void j_array_set_length(json_t *j, int64_t length) {
    if (!j_is_array(j)) {
        j_throw(j, "expected array");
    }
    while (length >= j->list_end - j->list) {
        j_new(j);
    }
    while (length < j->list_end - j->list) {
        j_set_null(j->list_end--);
    }
}

static inline json_t *j_array_append(json_t *j) {
    if (!j_is_array(j)) {
        j_throw(j, "expected array");
        return NULL;
    }
    return j_new(j);
}

static inline int64_t j_array_get_length(json_t *j) {
    if (!j_is_array(j)) {
        j_throw(j, "expected array");
        return 0;
    }
    return j->list_end - j->list;
}

////////////////////////////////////////////////////////////////////////////////
// magic parse
////////////////////////////////////////////////////////////////////////////////

#define j_magic_parse_filename(f, e) ( { \
        json_t *j=j_new(NULL); \
        j_check(j, "j_new() failed"); \
        j->longjmp_set=1; \
        char **ee = e; \
        if (setjmp(j->longjmp_env)) { \
            if (ee!=NULL) { \
                *(ee) = j->longjmp_string; \
            } else { \
                printf("EXCEPTION: %s\n",j->longjmp_string); \
                free(j->longjmp_string); \
            } \
            j_set_null(j); \
            free(j); \
            j=NULL; \
        } else { \
            j=_j_parse_filename(f,j); \
        } \
    j; })

////////////////////////////////////////////////////////////////////////////////
// convenience combinations for lookup + get
////////////////////////////////////////////////////////////////////////////////

static inline const char* j_lookup_string(json_t *j, const char *name) { return j_get_string(j_key_lookup(j,name)); }
static inline int64_t j_lookup_num64(json_t *j, const char *name) { return j_get_num64(j_key_lookup(j,name)); }
static inline double j_lookup_double(json_t *j, const char *name) { return j_get_double(j_key_lookup(j,name)); }
static inline mpz_t *j_lookup_bignum(json_t *j, const char *name) { return j_get_bignum(j_key_lookup(j,name)); }
static inline json_t *j_lookup_index(json_t *j, const char *name, int index) { return j_array_index(j_key_lookup(j, name),index); }
static inline int j_lookup_bool(json_t *j, const char *name) { return j_get_bool(j_key_lookup(j,name)); }

/*!
 * \endcond
 */

#endif

/*
 */
