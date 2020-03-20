##                   Makefile
###################################################################
 #                                                                #
 # Copyright (C) 2014, Institute for Defense Analyses             #
 # 4850 Mark Center Drive, Alexandria, VA; 703-845-2500           #
 # This material may be reproduced by or for the US Government    #
 # pursuant to the copyright license under the clauses at DFARS   #
 # 252.227-7013 and 252.227-7014.                                 #
 #                                                                #
 # LARC : Linear Algebra via Recursive Compression                #
 # Authors:                                                       #
 #   - Steve Cuccaro (IDA-CCS)                                    #
 #   - John Daly (LPS)                                            #
 #   - John Gilbert (UCSB, IDA adjunct)                           #
 #   - Jenny Zito (IDA-CCS)                                       #
 #                                                                #
 # Additional contributors are listed in "LARCcontributors".      #
 #                                                                #
 # Questions: larc@super.org                                      #
 #                                                                #
 # All rights reserved.                                           #
 #                                                                #
 # Redistribution and use in source and binary forms, with or     #
 # without modification, are permitted provided that the          #
 # following conditions are met:                                  #
 #   - Redistribution of source code must retain the above        #
 #     copyright notice, this list of conditions and the          #
 #     following disclaimer.                                      #
 #   - Redistribution in binary form must reproduce the above     #
 #     copyright notice, this list of conditions and the          #
 #     following disclaimer in the documentation and/or other     #
 #     materials provided with the distribution.                  #
 #   - Neither the name of the copyright holder nor the names of  #
 #     its contributors may be used to endorse or promote         #
 #     products derived from this software without specific prior #
 #     written permission.                                        #
 #                                                                #
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND         #
 # CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,    #
 # INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF       #
 # MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE       #
 # DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER NOR        #
 # CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,   #
 # SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT   #
 # NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;   #
 # LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)       #
 # HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN      #
 # CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR   #
 # OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, #
 # EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.             #
 #                                                                #
 ##################################################################
BINDIR = bin
OBJDIR = obj
SRCDIR = src

# Get the local paths 
include local/Makefile.conf
#    GMPIDIR (gmp-6 include),
#    GMPLDIR (gmp-6 lib),
#    MPIDIR (Anaconda3 include),
#    MPLDIR (Anaconda3 lib)

LIBDIR = lib
TESTDIR = tests
DOCDIR_HTML = html

# the following determines which scalarType LARC is compiled with:
# e.g. USE_REAL, USE_INTEGER, USE_COMPLEX, USE_MPINTEGER, ...
ifeq ($(TYPE),REAL)
   $(warning specified REAL type...)
   SOURCE_FILE = $(SRCDIR)/typeREAL.h
else ifeq ($(TYPE),INTEGER)
   $(warning specified INTEGER type...)
   SOURCE_FILE = $(SRCDIR)/typeINTEGER.h
else ifeq ($(TYPE),COMPLEX)
   $(warning specified COMPLEX type...)
   SOURCE_FILE = $(SRCDIR)/typeCOMPLEX.h
else ifeq ($(TYPE),MPREAL)
   $(warning specified MPREAL type...)
   SOURCE_FILE = $(SRCDIR)/typeMPREAL.h
else ifeq ($(TYPE),MPCOMPLEX)
   $(warning specified MPCOMPLEX type...)
   SOURCE_FILE = $(SRCDIR)/typeMPCOMPLEX.h
else ifeq ($(TYPE),MPINTEGER)
   $(warning specified MPINTEGER type...)
   SOURCE_FILE = $(SRCDIR)/typeMPINTEGER.h
else ifeq ($(TYPE),MPRATIONAL)
   $(warning specified MPRATIONAL type...)
   SOURCE_FILE = $(SRCDIR)/typeMPRATIONAL.h
else ifeq ($(TYPE),MPRATCOMPLEX)
   $(warning specified MPRATCOMPLEX type...)
   SOURCE_FILE = $(SRCDIR)/typeMPRATCOMPLEX.h
else
$(warning Using default real type...)
   SOURCE_FILE = $(SRCDIR)/typeREAL.h
endif


OPTS = -O2 -g -pg -Wall -fPIC
CFLAGS = $(OPTS) -I$(SRCDIR) -I$(MPIDIR) -I$(GMPIDIR) -std=gnu99
LIBS = -L$(MPLDIR) -L$(GMPLDIR) -lcurses -lm -lrt -lmpc -lmpfr -lgmp -lpthread -ltinfo
CC = gcc

# The file larcSWIG_wrap.c is created by SWIG, and therefore may or may not
# be present in the SRCDIR. We use SRCSNOWRAP to make sure we're not expecting
# the *wrap.o file to be in the OBJDIR or in OBJS.
SRCS = $(wildcard $(SRCDIR)/*.c)
SRCSNOWRAP = $(filter-out %_wrap.c, $(SRCS))
OBJS = $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(SRCSNOWRAP))
INCS = $(wildcard $(SRCDIR)/*.h)

PYINC = -I$(PYINCDIR)
PYINCDIR := $(shell python3 -c 'import sys; print(sys.prefix+"/include/python"+str(sys.version_info[0])+"."+str(sys.version_info[1])+sys.abiflags)')

TEST_SRCS = $(wildcard $(TESTDIR)/*_tests.c)
TESTS = $(patsubst %.c, %, $(TEST_SRC))

BIN_LIST = exampleLARC 
BINS = $(patsubst %,$(BINDIR)/%,$(BIN_LIST))

TARGET = lib/liblarc.a
SO_TARGET = $(patsubst %.a, %.so, $(TARGET))

DOCS = $(DOCDIR_HTML)/index.html

.PHONY: all build clean unittests # tests install

all: NEWTYPE $(SRCS) $(BINS) $(SO_TARGET) $(SRCDIR)/_larcSWIG.so $(TARGET) $(DOCS) # tests

# This is intended to check whether (re)compliation is necessary
NEWTYPE: FORCE
	@cmp -s $(SOURCE_FILE) $(SRCDIR)/type.h || cp -f $(SOURCE_FILE) $(SRCDIR)/type.h

FORCE:

# since there is no file called build, this was redoing the library every time
#$(TARGET): build $(OBJS)
# the .so and .a libraries need not include the python wrapper object
# (they also don't need the BIN_LIST objects, but this does no harm)
$(TARGET): $(OBJS)
	ar rcs $@ $(OBJS)
	ranlib $@

$(SO_TARGET): $(OBJS) # $(TARGET)
	$(CC) -shared -o $@ $(OBJS)

build: 
	@mkdir -p obj
	@mkdir -p bin
	@mkdir -p lib

# tests: CFLAGS += $(TARGET)
# tests: $(TESTS)
#	sh ./tests/runtests.sh

clean: 
	@# again, there is no file called build, so no need to remove it
	@# rm -rf build $(OBJS) $(TESTS) $(BINS) 
	rm -rf $(OBJS) $(TESTS) $(BINS) 
	rm -f $(TESTDIR)/python/*.pyc
	rm -f $(SO_TARGET) $(TARGET)
	rm -f $(SRCDIR)/*.pyc $(SRCDIR)/_pylarc.so
	rm -f $(SRCDIR)/larcSWIG.py $(SRCDIR)/_larcSWIG.so
	rm -f $(SRCDIR)/*_wrap.o $(SRCDIR)/*_wrap.c
	rm -rf $(DOCDIR_HTML)
	rm -f tests/tests.log
	find . -name "*.gc*" -exec rm {} \;
	rm -rf `find . -name "*.dSYM" -print`
	rm -rf src/__pycache__
	rm -rf src/python/__pycache__
	rm -rf tests/python/__pycache__
	rm -f test_output

#install: all
#	install -d $(DESTDIR)/$(PREFIX)/lib/
#	install $(TARGET) $(DESTDIR)/$(PREFIX)/lib/  

# we keep all the SWIG dependent stuff in SRCDIR to avoid problems later
$(SRCDIR)/_larcSWIG.so: $(SRCDIR)/larcSWIG_wrap.o $(OBJS)
	#$(CC) $(CFLAGS) -shared $^ -o $@ $(LIBS)
	$(CC) $(CFLAGS) -shared $^ -o $@ $(LIBS) -Wl,-rpath=$(MPLDIR),-rpath=$(GMPLDIR)

$(SRCDIR)/larcSWIG_wrap.o: $(SRCDIR)/larcSWIG_wrap.c
	$(CC) $(CFLAGS) $(PYINC) -c $< -o $@

# using SRCSNOWRAP avoids a dependency loop (larcSWIG_wrap.c shouldn't depend
# on itself)
$(SRCDIR)/larcSWIG_wrap.c: $(SRCDIR)/larcSWIG.i $(SRCSNOWRAP) $(INCS)
	swig -python $<

# the binary executables do not depend on the python wrapper
$(BINS): $(OBJS) $(BINDIR)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS) -Wl,-rpath=$(MPLDIR),-rpath=$(GMPLDIR)

$(BINDIR): 
	@test -d $(BINDIR)/ || mkdir -p $(BINDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.c $(INCS)
	@test -d $(OBJDIR)/ || mkdir -p $(OBJDIR)
	$(CC) -c -o $@ $< $(CFLAGS)

$(DOCS): $(SRCS) $(SRCDIR)/Doxyfile.in
	doxygen $(SRCDIR)/Doxyfile.in



unittests:
	echo "Running unittests for all types:" > test_output
	#python2.7 -m unittest discover -p 'test_unittest*.py' -s tests/python 2>> test_output

	make TYPE=INTEGER
	echo "TYPE=INTEGER" >> test_output
	python3 -m unittest discover -p 'test_unittest*.py' -s tests/python 2>> test_output
	#tests/python/run_unittests.py 2>> test_output

	make TYPE=REAL
	echo "TYPE=REAL" >> test_output
	python3 -m unittest discover -p 'test_unittest*.py' -s tests/python 2>> test_output
	#tests/python/run_unittests.py 2>> test_output

	make TYPE=COMPLEX
	echo "TYPE=COMPLEX" >> test_output
	python3 -m unittest discover -p 'test_unittest*.py' -s tests/python 2>> test_output
	#tests/python/run_unittests.py 2>> test_output

	make TYPE=MPINTEGER
	echo "TYPE=MPINTEGER" >> test_output
	python3 -m unittest discover -p 'test_unittest*.py' -s tests/python 2>> test_output
	#tests/python/run_unittests.py 2>> test_output

	make TYPE=MPREAL
	echo "TYPE=MPREAL" >> test_output
	python3 -m unittest discover -p 'test_unittest*.py' -s tests/python 2>> test_output
	#tests/python/run_unittests.py 2>> test_output

	make TYPE=MPCOMPLEX
	echo "TYPE=MPCOMPLEX" >> test_output
	python3 -m unittest discover -p 'test_unittest*.py' -s tests/python 2>> test_output
	#tests/python/run_unittests.py 2>> test_output

	make TYPE=MPRATIONAL
	echo "TYPE=MPRATIONAL" >> test_output
	python3 -m unittest discover -p 'test_unittest*.py' -s tests/python 2>> test_output
	#tests/python/run_unittests.py 2>> test_output

	make TYPE=MPRATCOMPLEX
	echo "TYPE=MPRATCOMPLEX" >> test_output
	python3 -m unittest discover -p 'test_unittest*.py' -s tests/python 2>> test_output
	#tests/python/run_unittests.py 2>> test_output

	echo "ALL TESTS COMPLETED" >> test_output
	cat test_output
#test:
#	make clean
#	make all
#	cd $(TESTDIR) && tcsh -f larc_test

