##                   README.md

                                                                 
  Copyright (C) 2014, Institute for Defense Analyses             
  4850 Mark Center Drive, Alexandria, VA; 703-845-2500           
  This material may be reproduced by or for the US Government    
  pursuant to the copyright license under the clauses at DFARS   
  252.227-7013 and 252.227-7014.                                 
                                                                 
  LARC : Linear Algebra via Recursive Compression                
  Authors:                                                       
    - Steve Cuccaro (IDA-CCS)                                    
    - John Daly (LPS)                                            
    - John Gilbert (UCSB, IDA adjunct)                           
    - Mark Pleszkoch (IDA-CCS)                                     
    - Jenny Zito (IDA-CCS)                                       
                                                                 
  Additional contributors are listed in "LARCcontributors".      
                                                                 
  Questions: larc@super.org                                      
                                                                 
  All rights reserved.                                           
                                                                 
  Redistribution and use in source and binary forms, with or     
  without modification, are permitted provided that the          
  following conditions are met:                                  
    - Redistribution of source code must retain the above        
      copyright notice, this list of conditions and the          
      following disclaimer.                                      
    - Redistribution in binary form must reproduce the above     
      copyright notice, this list of conditions and the          
      following disclaimer in the documentation and/or other     
      materials provided with the distribution.                  
    - Neither the name of the copyright holder nor the names of  
      its contributors may be used to endorse or promote         
      products derived from this software without specific prior 
      written permission.                                        
                                                                 
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND         
  CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,    
  INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF       
  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE       
  DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER NOR        
  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,   
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT   
  NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;   
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)       
  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN      
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR   
  OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.             
                                                                 
 

Welcome to the Version 2.1 LARC library implementing

    LARC = Linear Algebra via Recursive Compression.

LARC is useful for carrying out linear algebra on very
large power of 2 dimensioned matrices and vectors which
have some internal repetition or structure.  LARC recursively
compresses matrices and carries out matrix operations
while in the compressed format.

The core code of LARC is written in C and there is a
a SWIG-generated Python3 interface for ease of use.
The python module pylarc.py imports SWIG-generated python
wrappers for the C code and contains additional
routines that are written directly in python.

LARC contains some Python testing routines that are in
the directory tests/python.
There are also commandline python routines which are
not in pylarc and are mostly in the directory src/python.

There is a set of slides in the doc file "about\_LARC.pdf"
which describe the ideas behind LARC.

**We strongly recommend that those using LARC for the first time download from GitHub/LARCmath the package MyPyLARC.
MyPyLARC is a tutorial and example package which uses
LARC as its math package, and can serve as a template
for writing one's own package using LARC. It will
automatically load a copy of LARC into a larc subdirectory.
MyPyLARC is written primarily in Python with user
friendly example routines, including a few in
Jupyter Notebook.  There is an explanatory paper on
LARC and MyPyLARC in MyPyLARC/doc/LARCandMyPyLARC.pdf,
and also a conference poster at
MyPyLARC/doc/LARCposterSIAMAnnualMeeting2020.pdf.**

LARC version 2.1 has not been optimized for any particular
application and MyPyLARC version 2.1 is still evolving.

We have Doxygen documentation for LARC in HTML
and LaTeX formats. To view this documentation:
open the file "html/index.html" in a browser
window.  This command will work in Linux:

    gio open html/index.html

To get started, compile a fresh copy of the code
in the top-level LARC directory by typing "make".
On the first attempt to compile, you will get a message
instructing you to create a file in LARC/local
called Makefile.conf. Several examples are in
LARC/local which can be copied to your Makefile.conf
and modified so that the paths reflect the paths
in your environment to various software packages.

As of August 2022, LARC uses:
\*   Python 3.6 or later
\*   GNU GMP 5.0 or later
\*   GNU MPFR 4.0.0 or later
\*   GNU MPC 1.1.0 or later
\*   SWIG 3.0.0 or later (to create Python interface for C routines)
\*   Doxygen-1.8.13 exact version (for documentation)
\*	* 1.8.x works but with warnings that can be ignored

We also recommend
\*   GCC 10 or later

Python3 and the multiprecision libraries can be obtained from the
Anaconda package, version 2020.02 or later, or may be downloaded and
installed individually. These should be available on line for upload
at several places including: https://anaconda.com, https://www.swig.org,
https://www.doxygen.nl, and https://gcc.gnu.org.

If you are not downloading MyPyLARC you can see some
examples of how to run LARC in these places:
\* examples of how the LARC code works from C
  (this is also a good test before doing git commits):
\*    cd bin
\*    ./exampleLARC
\* examples of how the to use LARC from Python
  (or checking the code before git commit):
    type (from the top-level project directory):
\*  ./tests/python/test\_matrixBuild.py
\*  ./tests/python/test\_math.py
\*  ./tests/python/test\_localHash.py
  etc.
  To run all the python tests say:
\*  ./tests/python/script

The following Makefile commands tell LARC to define its
scalarType to be one of the pre-defined C datatypes:
\*    make TYPE=COMPLEX (long double complex)
\*    make TYPE=REAL (long double)
\*    make TYPE=INTEGER (int64\_t)

as well as five multiprecision Scalartypes
\*    make TYPE=MPINTEGER    (mpz\_t)
\*    make TYPE=MPREAL       (256-bit mpfr\_t)
\*    make TYPE=MPCOMPLEX    (256-256-bit mpc\_t)
\*    make TYPE=MPRATIONAL   (mpq\_t)
\*    make TYPE=MPRATCOMPLEX (larc\_mpratcomplex\_t)

or one of the specialty types:
\*    make TYPE=BOOLEAN      (int64\_t)
\*    make TYPE=CLIFFORD     (clifford\_t)
\*    make TYPE=UPPER        (larc\_exponent\_scalar\_t)
\*    make TYPE=LOWER        (larc\_exponent\_scalar\_t)

the MPRATCOMPLEX type is constructed from two copies of mpq\_t
for the real and imaginary parts of a complex rational number.
If no TYPE option is specified, the default is REAL (C long double).

Before commiting a new version of code to your local
git repository, all the unit tests should be run
by typing:
\*  make clean
\*  make unittests

The different ways to set initialization parameters
for LARC are described in the Doxygen
documentations and illustrated various places such as
in tests/python code.

As we said above, those using LARC for the first time may find it helpful
to download from GitHub/LARCmath the package MyPyLARC and use its 
tutorials and examples. The MyPyLARC README.md gives guidance on
getting the package working, and the file Tutorial/newuser\_instructions
suggests which order to do the tutorials.

#### List of math operations from src/matmath.h
          (For details, also see the doxygen page on larc/matmath.h)

void get\_array\_of\_scalars\_in\_larcMatrixFile(scalarType \*\*scalars\_ptr,
                                   int64\_t \*numScalars, char \*path);
				   
int64\_t matrix\_add(int64\_t A\_pID, int64\_t B\_pID);

int64\_t matrix\_diff(int64\_t A\_pID, int64\_t B\_pID);

int64\_t matrix\_mult(int64\_t A\_pID, int64\_t B\_pID);

int64\_t matrix\_mult\_clean(int64\_t A\_pID, int64\_t B\_pID,
                                   mat\_level\_t cleanThresh);
				   
int64\_t scalar\_mult(int64\_t A\_pID, int64\_t B\_pID);

int64\_t scalar\_divide(int64\_t A\_pID, int64\_t B\_pID);

int64\_t kronecker\_product(int64\_t A\_pID, int64\_t B\_pID);

int64\_t join(int64\_t A\_pID, int64\_t B\_pID);

int64\_t stack(int64\_t A\_pID, int64\_t B\_pID);

int64\_t matrix\_entrySquared(int64\_t m\_pID, char \*scale\_factor);

int64\_t iHadamard\_times\_matrix(int64\_t A\_pID);

int64\_t matrix\_basischange\_A\_by\_B(int64\_t B\_pID, int64\_t A\_pID);

int64\_t matrix\_saxpy(int64\_t A\_pID, int64\_t B\_pID, 
                      int64\_t scalar\_a\_pID, int64\_t scalar\_b\_pID);
			      
int64\_t vector\_dot\_product(int64\_t A\_pID, int64\_t B\_pID, int verbose);

int64\_t matrix\_times\_iHadamard(int64\_t A\_pID);

char \*tracenorm(int64\_t m\_pID, char \*scale\_factor);

char \*traceID(int64\_t m\_pID);

char \*get\_scalar\_value\_string(int64\_t m\_pID);

char \*matrix\_count\_entries(int64\_t mat\_pID, char \*scalar\_str);

char \*get\_list\_\of\_scalars\_in\_larcMatrixFile(char \*path);

int64\_t random\_bool\_matrix\_from\_count(mat\_level\_t row\_level,
	mat\_level\_t col\_level, char\* numOnes);
					
int64\_t matrix\_element\_with\_maxNorm(int64\_t mat\_pID);

int64\_t normID(int64\_t mat\_pID, int whichNorm);

int64\_t create\_const\_matrix(int64\_t constID, mat\_level\_t row\_level,
	mat\_level\_t col\_level);

int64\_t apply\_function\_to\_matrix\_values(int64\_t m\_pID,
        void (*func)(scalarType*, const scalarType), op\_type\_t op\_memz);

