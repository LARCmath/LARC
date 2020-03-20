//                          json.c
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

// Author: Michael Schneider

#include "json.h"
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <limits.h>
#include <strings.h>

#define PRETTY_WRAPPING

#ifdef PRETTY_WRAPPING
#include <sys/ioctl.h>
static int wrap_width() {
    static struct winsize w;
    static int did_ioctl = 0;
    if (!did_ioctl) {
        memset(&w, 0, sizeof(struct winsize));
        ioctl(0, TIOCGWINSZ, &w);
        if (w.ws_col == 0) w.ws_col = SHRT_MAX;
        did_ioctl = 1;
    }
    return w.ws_col;
}
#define WRAP_WIDTH wrap_width()
#else
#define WRAP_WIDTH 0
#endif

static inline size_t _j_snprint_path(char *str, size_t size, json_t *j) {
    size_t used = 0;
    int r = -1;
    char *str2 = (size>used) ? str+used : NULL;
    size_t size2 = (size>used) ? size-used : 0;
    if (!j) {
        r = snprintf(str2, size2, "<null>");
    } else if (j->parent) {
        used += _j_snprint_path(str, size, j->parent);
        str2 = (size>used) ? str+used : NULL; // since "used" was updated
        size2 = (size>used) ? size-used : 0; // since "used" was updated
        if (j->parent->type == J_OBJECT) {
            r = snprintf(str2, size2, "->%s", j->name); // parent is object
        } else if (j->parent->type == J_ARRAY) {
            int64_t index = j - j->parent->list;
            r = snprintf(str2, size2, "[%ld]", index); // parent is array
        } else {
            j_error("unexpected parent type");
        }
    } else if (j->type == J_OBJECT) {
        if (j->orig_filename) {
            r = snprintf(str2, size2, "\"%s\"", j->orig_filename);
        } else {
            r = snprintf(str2, size2, "<json_object>"); // top level object
        }
    } else if (j->type == J_ARRAY) {
        if (j->orig_filename) {
            r = snprintf(str2, size2, "\"%s\"", j->orig_filename);
        } else {
            r = snprintf(str2, size2, "<json_array>"); // top level object
        }
    } else {
        if (j->orig_filename) {
            r = snprintf(str2, size2, "\"%s\"", j->orig_filename);
        } else {
            r = snprintf(str2, size2, "<json_value>"); // top level value
        }
    }
    j_check(r >= 0, "snprintf output error");
    used += r;
    return used;
}

static inline char *_j_get_name(json_t *j) {
    size_t s1 = _j_snprint_path(NULL, 0, j);
    char *rv = malloc(s1+1);
    j_check(rv, "malloc(%lu) failed",s1+1);
    size_t s2 = _j_snprint_path(rv, s1+1, j);
    j_check(s2 == s1,"string length mismatch");
    return rv;
}

void j_throw(json_t *j, const char *format, ...) {
    const char *throw_message = "%s in %s";
    char *msg;
    SPRINTF(msg, format);
    char *name = _j_get_name(j);
    int s1 = snprintf(NULL,0,throw_message,msg,name);
    j_check(s1 >= 0, "snprintf output error"); \
    char *msg2 = malloc(s1+1);
    j_check(msg2, "malloc(%u) failed",s1+1); \
    int s2 = snprintf(msg2,s1+1,throw_message,msg,name);
    j_check(s2 == s1,"string length mismatch"); \
    free(name);
    free(msg);

    if (0) printf("(debug) %s\n",msg2); // debug

    json_t *root = j;
    if (root) while (root->parent) root = root->parent;
    if (root && root->longjmp_set) {
        if (root->longjmp_string) {
            fprintf(stderr,"WARNING! lost %s\n",msg2);
            free(msg2);
        } else {
            root->longjmp_string = msg2;
        }
        longjmp(root->longjmp_env, 1);
        j_error("longjmp returned?!");
    } else {
        fprintf(stderr,"(no exception handler) %s\n",msg2);
        free(msg2);
        abort();
    }
}

static json_t *_parse_error(json_t *current, int linecount, const char *c, const char *message) {

    j_throw(current, "\"%s\" parsing line %d",message,linecount);

    /*
    char *name = _j_get_name(current);
    printf("parse error '%s' on line %d while parsing %s\n",message,linecount,name);
    free(name);
    */
    while (current->parent) current = current->parent;
    j_set_null(current);
    /*
       free(current);
    return NULL;
    */
    return current;
}

static int _is_double(const char *p)
{
  int i;
  if (*p == '-') p++; // skip any leading '-'
  i = strspn(p, "0123456789");
  if (p[i] == '.' || p[i] == 'e' || p[i] == 'E')
    return 1;
  return 0;
}

static char *_parse_number(json_t *j, const char *p) {
    char *end;
    if (_is_double(p)) {
        j->num_double = strtod(p, &end);
        j->type = J_DOUBLE;
    } else {
        j->num64 = strtoll(p, &end, 0);
        //if (j->num64 == LLONG_MAX || j->num64 == LLONG_MIN) {
        // alternate test that doesn't use LLONG_* constants
        if (((j->num64 > 0) && (j->num64+1 < 0)) || ((j->num64 < 0) && (j->num64-1 > 0))) {
            if (!j->bignum_initialized) {
                mpz_init(j->bignum);
                j->bignum_initialized = 1;
            }

            int neg = (*p == '-') ? 1 : 0;
            int i = strspn(p+neg, "0123456789") + neg;
            char p1[i+1];
            strncpy(p1, p, i);
            p1[i] = 0;
            end = (char *)p+i;
            mpz_set_str(j->bignum, p1, 10);
            j->type = J_BIGNUM;
            j->num_double = mpz_get_d(j->bignum);
        } else {
            j->type = J_NUM64;
            j->num_double = j->num64;
        }
    }
    return end;
}

static void string_append(json_t *j, char c) {
    if (!j->string_end) { // initial allocation
        size_t alloc = 64;
        j->string = malloc(alloc);
        j->string_end = j->string + 0;
        j->string_alloc_end = j->string + alloc;
    } else if (j->string_end+1 >= j->string_alloc_end) { // extend allocation
        size_t size = j->string_end - j->string;
        size_t alloc = j->string_alloc_end - j->string;
        alloc *= 2;
        j->string = realloc(j->string, alloc);
        j->string_end = j->string + size;
        j->string_alloc_end = j->string + alloc;
    }
    if (c) {
        *j->string_end++ = c;
        *j->string_end = 0;
    }
}

static json_t *j_new_root() {
    json_t *j = malloc(sizeof(json_t));
    memset(j, 0, sizeof(json_t));
    j->type = J_NULL;
    return j;
}

static json_t *new_json(json_t *parent) {
    j_check(parent, "null parent?");
    j_check((parent->type == J_OBJECT) || (parent->type == J_ARRAY), "bad parent type");
    if (!parent->list_end) { // initial allocation
        size_t alloc = 16;
        parent->list = malloc(alloc * sizeof(json_t));
        parent->list_end = parent->list + 0;
        parent->list_alloc_end = parent->list + alloc;
    } else if (parent->list_end >= parent->list_alloc_end) { // extend allocation
        size_t size = parent->list_end - parent->list;
        size_t alloc = parent->list_alloc_end - parent->list;
        alloc *= 2;
        parent->list = realloc(parent->list, alloc * sizeof(json_t));
        parent->list_end = parent->list + size;
        parent->list_alloc_end = parent->list + alloc;
    }
    memset(parent->list_end, 0, sizeof(json_t));
    parent->list_end->parent = parent;
    return parent->list_end++;
}

json_t *j_new(json_t *parent) {
    json_t *p = parent ? new_json(parent) : j_new_root();
    p->type = J_NULL;
    return p;
}

json_t *_j_parse_string(const char *c, json_t *root) {
    json_t *j = NULL;
    j_check(root, "null root?");
    j_check(root->type == J_NULL, "root reused?");
    int linecount = 1;
    enum {NORMAL, COMMENT, STRING} lex_state = NORMAL;
    enum {INITIAL, OBJ_NAME, COLON, VALUE, COMMA} parse_state = INITIAL;
    for (; *c; c++) {
        //printf("char is %c, lex_state is %d, parse_state is %d\n",*c,lex_state,parse_state);
        if (*c == '\n') linecount++;
        if (lex_state == COMMENT) {
            if (*c == '\n') lex_state = NORMAL;
        } else if (lex_state == STRING) {
            if (*c == '\\') {
                c++;
                switch (*c) {
                    case 0: return _parse_error(j,linecount,c,"unterminated string");
                    case '"':
                    case '\\':
                    case '/': string_append(j, *c); break;
                    case 'b': string_append(j, '\b'); break;
                    case 'f': string_append(j, '\f'); break;
                    case 'n': string_append(j, '\n'); break;
                    case 'r': string_append(j, '\r'); break;
                    case 't': string_append(j, '\t'); break;
                    case 'u': return _parse_error(j,linecount,c,"unicode not supported");
                    default: return _parse_error(j,linecount,c,"bad escape sequence in string");
                }
            } else if (*c == '"') {
                string_append(j, 0); // ensure string exists even if zero length
                lex_state = NORMAL;
                if (parse_state == OBJ_NAME) {
                    j->name = j->string;
                    j->string = j->string_end = j->string_alloc_end = NULL;
                    parse_state = COLON;
                } else if (parse_state == VALUE) {
                    j_check(j->type == J_STRING, "parser bug");
                    j = (j == root) ? NULL : j->parent; parse_state = j ? COMMA : INITIAL;
                } else j_error("parser bug");
            } else {
                string_append(j, *c);
            }
        } else if (isspace(*c)) {
            // do nothing for whitespace outside a comment or string
        } else if (*c == '/') {
            if (c[1] == '/') { c++; lex_state = COMMENT; }
            else return _parse_error(j,linecount,c,"unexpected /");
        } else if (parse_state == INITIAL) {
            // if j, it's really FINAL
            if (root->type != J_NULL) {
                return _parse_error(j,linecount,c,"extra text at end");
            } else {
                j = root;
                if (*c == '{') {
                    parse_state = OBJ_NAME;
                    root->type = J_OBJECT;
                    j = new_json(root);
                } else if (*c == '[') {
                    parse_state = VALUE;
                    root->type = J_ARRAY;
                    j = new_json(root);
                } else return _parse_error(j,linecount,c,"expected [ or {");
            }
        } else if (parse_state == OBJ_NAME) {
            if (*c == '"') { lex_state = STRING; string_append(j, 0); }
            else if (*c == '}') {
                j = j->parent;
                if ((j->type == J_OBJECT) && (j->list+1 == j->list_end)) {
                    // allow empty object
                    j->list_end--;
                    j = (j == root) ? NULL : j->parent; parse_state = j ? COMMA : INITIAL;
                } else return _parse_error(j,linecount,c,"unexpected }");
            } else return _parse_error(j,linecount,c,"expected object name");
        } else if (parse_state == COLON) {
            if (*c == ':') parse_state = VALUE;
            else return _parse_error(j,linecount,c,"expected colon");
        } else if (parse_state == VALUE) {
            switch (*c) {
                case '"': lex_state = STRING; j->type = J_STRING; string_append(j, 0); break;

                case '0': case '1': case '2': case '3': case '4':
                case '5': case '6': case '7': case '8': case '9':
                case '-':
                          c = _parse_number(j, c)-1;
                          j = (j == root) ? NULL : j->parent; parse_state = j ? COMMA : INITIAL;
                          break;

                case '{': parse_state = OBJ_NAME; j->type = J_OBJECT; j = new_json(j); break;
                case '[': parse_state = VALUE; j->type = J_ARRAY; j = new_json(j); break;
                case '}': return _parse_error(j,linecount,c,"unexpected }");
                case ']': j = j->parent;
                          j_check(j, "null parent?");
                          if ((j->type == J_ARRAY) && (j->list+1 == j->list_end)) {
                              // allow empty array
                              j->list_end--;
                              j = (j == root) ? NULL : j->parent; parse_state = j ? COMMA : INITIAL;
                          } else return _parse_error(j,linecount,c,"unexpected ]");
                          break;
                case 't': if (!strncmp(c,"true",4)) {
                              j->type = J_TRUE;
                              c += 3;
                              j = (j == root) ? NULL : j->parent; parse_state = j ? COMMA : INITIAL;
                          } else return _parse_error(j,linecount,c,"bad value, hoping for true");
                          break;
                case 'f': if (!strncmp(c,"false",5)) {
                              j->type = J_FALSE;
                              c += 4;
                              j = (j == root) ? NULL : j->parent; parse_state = j ? COMMA : INITIAL;
                          } else return _parse_error(j,linecount,c,"bad value, hoping for false");
                          break;
                case 'n': if (!strncmp(c,"null",4)) {
                              j->type = J_NULL;
                              c += 3;
                              j = (j == root) ? NULL : j->parent; parse_state = j ? COMMA : INITIAL;
                          } else return _parse_error(j,linecount,c,"bad value, hoping for null");
                          break;
                default: return _parse_error(j,linecount,c,"unrecognized unquoted string");

            }
        } else if (parse_state == COMMA) {  // actually  "," or "}" or "]"
            if (*c == ',') {
                if (j->type == J_ARRAY) parse_state = VALUE;
                else if (j->type == J_OBJECT) parse_state = OBJ_NAME;
                else j_error("parser bug");
                j = new_json(j);
            } else if (*c == '}') {
                if (j->type == J_OBJECT) {
                    j = (j == root) ? NULL : j->parent; parse_state = j ? COMMA : INITIAL;
                } else if (j->type == J_ARRAY) {
                    return _parse_error(j,linecount,c,"expected comma or ]");
                } else {
                    j_error("parser bug");
                }
            } else if (*c == ']') {
                if (j->type == J_ARRAY) {
                    j = (j == root) ? NULL : j->parent; parse_state = j ? COMMA : INITIAL;
                } else if (j->type == J_OBJECT) {
                    return _parse_error(j,linecount,c,"expected comma or }");
                } else {
                    j_error("parser bug");
                }
            } else return _parse_error(j,linecount,c,"expected comma");
        }
    }
    if (j) return _parse_error(j,linecount,c,"unterminated object or array");
    return root;
}

json_t *_j_parse_file(FILE *fp, json_t *root) {
    if (!fp) return NULL;
    int rv = fseek(fp, 0, SEEK_END);
    if (rv < 0) { perror("j_parse_file: fseek: "); return NULL; }
    long end = ftell(fp);
    if (end < 0) { perror("j_parse_file: ftell: "); return NULL; }
    rv = fseek(fp, 0, SEEK_SET);
    if (rv < 0) { perror("j_parse_file: fseek: "); return NULL; }
    // printf("  End is %ld\n",end);
    // char buf[end+1];
    char *buf = malloc(end+1);
    size_t size = fread(buf,1,end+1,fp);
    j_check(size == end, "incomplete fread");
    //printf("j_parse_file read %ld bytes\n", size);
    buf[size] = 0;
    //printf("Buffer was:\n%s", buf);
    //return _j_parse_string(buf, root);
    json_t *rvs = _j_parse_string(buf, root);
    free(buf);
    return rvs;
}

json_t *_j_parse_filename(const char *filename, json_t *root) {
    FILE *fp = fopen(filename, "r");
    if (!fp) { perror("_j_parse_filename: Could not open file.\n"); return root; }
    json_t *rv = _j_parse_file(fp, root);
    fclose(fp);
    rv->orig_filename = strdup(filename);
    return rv;
}

json_t *j_parse_string(const char *c) {
    json_t *root = j_new_root();
    return _j_parse_string(c, root);
}

json_t *j_parse_file(FILE *fp) {
    json_t *root = j_new_root();
    return _j_parse_file(fp, root);
}

json_t *j_parse_filename(const char *filename) {
    json_t *root = j_new_root();
    return _j_parse_filename(filename, root);
}

json_t *j_parse_filename_into_null_object(const char *filename, json_t *root) {
  return _j_parse_filename(filename, root);
}


int snprint_spaces(char *dst, size_t n, int count) {
    int i;
    for (i=0; i<count; i++) {
        if (i < n) dst[i] = ' ';
    }
    if (i < n) dst[i] = 0;
    else if (n > 0) dst[n-1] = 0;
    return count;
}

#define P0(dst,pos,n,format)     { int rv=snprintf(dst+pos,(n>pos)?n-pos:0,format);         if (rv<0) return rv; else pos+=rv; }
#define P1(dst,pos,n,format,arg) { int rv=snprintf(dst+pos,(n>pos)?n-pos:0,format,arg);     if (rv<0) return rv; else pos+=rv; }
#define PG(dst,pos,n,format,arg) { int rv=gmp_snprintf(dst+pos,(n>pos)?n-pos:0,format,arg); if (rv<0) return rv; else pos+=rv; }
#define PS(dst,pos,n,count)      { int rv=snprint_spaces(dst+pos,(n>pos)?n-pos:0,count);    if (rv<0) return rv; else pos+=rv; }

int _j_snprint(char *dst, int pos, size_t n, json_t *j, int indent) {
    int rv;
    if (!j) {
        P0(dst,pos,n,"null");
    } else {
        if (j->name) {
            P1(dst,pos,n,"\"%s\":",j->name);
        }
        switch (j->type) {
            case J_STRING: P1(dst,pos,n,"\"%s\"",j->string); break;
            case J_OBJECT: {
                               json_t *p;
                               if (j->name) {
                                   P0(dst,pos,n,"\n");
                                   PS(dst,pos,n,indent);
                               }
                               P0(dst,pos,n,"{ ");
                               int initial_pos = pos;
                               int bottom_level = 1;
                               for (p=j->list; p<j->list_end; p++) {
                                   if (j_is_object(p) || j_is_array(p)) {
                                       bottom_level = 0;
                                       break;
                                   }
                               }
                               for (p=j->list; p<j->list_end; p++) {
                                   rv = _j_snprint(dst,pos,n,p,indent+2);
                                   if (rv < 0) return rv; else pos = rv;
                                   if (p+1 < j->list_end) {
                                       if (bottom_level && (pos - initial_pos + 20 < WRAP_WIDTH)) {
                                           P0(dst,pos,n,", ");
                                       } else {
                                           P0(dst,pos,n,",\n");
                                           PS(dst,pos,n,indent+2);
                                           initial_pos = pos;
                                       }
                                      /*
                                       P0(dst,pos,n,",\n");
                                       PS(dst,pos,n,indent+2);
                                       */
                                   }
                               }
                               P0(dst,pos,n," }");
                               break;
                           }
            case J_ARRAY: {
                              json_t *p;
                              if (j->name) {
                                  P0(dst,pos,n,"\n");
                                  PS(dst,pos,n,indent);
                              }
                              P0(dst,pos,n,"[ ");
                              int initial_pos = pos;
                              int bottom_level = 1;
                              for (p=j->list; p<j->list_end; p++) {
                                  if (j_is_object(p) || j_is_array(p)) {
                                      bottom_level = 0;
                                      break;
                                  }
                              }
                              for (p=j->list; p<j->list_end; p++) {
                                  rv = _j_snprint(dst,pos,n,p,indent+2);
                                  if (rv < 0) return rv; else pos = rv;
                                  if (p+1 < j->list_end) {
                                      if (bottom_level && (pos - initial_pos + 20 < WRAP_WIDTH)) {
                                          P0(dst,pos,n,", ");
                                      } else {
                                          P0(dst,pos,n,",\n");
                                          PS(dst,pos,n,indent+2);
                                          initial_pos = pos;
                                      }
                                  }
                              }
                              P0(dst,pos,n," ]");
                              break;
                          }
            case J_TRUE:   P0(dst,pos,n,"true"); break;
            case J_FALSE:  P0(dst,pos,n,"false"); break;
            case J_NULL:   P0(dst,pos,n,"null"); break;
            case J_NUM64:  P1(dst,pos,n,"%ld", j->num64); break;
            case J_DOUBLE: P1(dst,pos,n,"%lg", j->num_double); break;
            case J_BIGNUM: PG(dst,pos,n,"%Zd", j->bignum); break;
            default:       P0(dst,pos,n,"<something_else>"); break;
        }
    }
    return pos;
}

int j_snprint(char *dst, size_t n, json_t *j) {
    int pos = _j_snprint(dst, 0, n, j, 0);
    if (pos >= 0) P0(dst,pos,n,"\n");
    return pos;
}

char *j_make_rel_filename(const json_t *j, const char *filename) {
    if (!filename) return NULL;
    if (filename[0] == '/') return strdup(filename);
    while (j && !j->orig_filename && j->parent) j = j->parent;
    const char *slash = (j && j->orig_filename) ? rindex(j->orig_filename, '/') : (char *)NULL;
    int len1 = slash ? (slash - j->orig_filename + 1) : 0;
    int len2 = strlen(filename);
    char *rv = malloc(len1+len2+1);
    strncpy(rv, j->orig_filename, len1);
    strncpy(rv+len1, filename, len2);
    rv[len1+len2] = 0;
    return rv;
}

#ifdef NEED_JSON_MAIN
int main() {
    char buf[1024*1024];
    size_t size = fread(buf,1,1024*1024,stdin);
    buf[size] = 0;
    json_t *j = j_parse_string(buf);
    int s = j_snprint(buf, 0, j);
    j_check(s >= 0, "j_snprint returned error");
    char buf2[s+1];
    int s2 = j_snprint(buf2, s+1, j);
    j_check(s == s2, "string length mismatch");
    printf("%s",buf2);
    printf("\n");
    printf("s2 was %d\n",s2);
    return 0;
}
#endif

/*
 */
