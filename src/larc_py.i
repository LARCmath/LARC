%module larc_py

%{
#include "type.h"
#include "matmath.h"
#include "larc.h"
#include "io.h"
#include "global.h" 
#include "matrix_store.h"
#include "op_store.h"
#include "organize.h"
#include "info_store.h"
#include "fft.h"
#include "hash.h"
#include <complex.h>
#include <pthread.h>
%}

#define SWIGWORDSIZE64
#enddef

%include "stdint.i"
%include "complex.i"
%include "cstring.i"
%include "carrays.i"
%include "cpointer.i"

# NOTE: global.h uses the USE_INTEGER/REAL/COMPLEX that is #defined in
# type.h, so type.h must be before it in the list

%include "type.h"
%include "matmath.h"
%include "larc.h"
%include "io.h"
%include "global.h"
%include "matrix_store.h"
%include "op_store.h"
%include "organize.h"
%include "info_store.h"
%include "hash.h"
%include "fft.h"

%array_class(int, intArray);
%array_class(long int, int64Array); # works because SWIGWORDSIZE64 defined
%array_class(double, realArray);
%array_class(complex, complexArray);
%array_functions(ScalarType, scalarTypeArray);
