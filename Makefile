#################################################################
#                                                               #
# Copyright 2014, Institute for Defense Analyses                #
# 4850 Mark Center Drive, Alexandria, VA; 703-845-2500          #
# This material may be reproduced by or for the US Government   #
# pursuant to the copyright license under the clauses at DFARS  #
# 252.227-7013 and 252.227-7014.                                #
#                                                               #
# LARC : Linear Algebra via Recursive Compression               #
# Authors:                                                      #
#   - Steve Cuccaro (IDA-CCS)                                   #
#   - John Daly (LPS)                                           #
#   - John Gilbert (UCSB, IDA adjunct)                          #
#   - Jenny Zito (IDA-CCS)                                      #
#                                                               #
# Additional contributors are listed in "LARCcontributors".     #
#                                                               #
# POC: Jennifer Zito <jszito@super.org>                         #
# Please contact the POC before disseminating this code.        #
#                                                               #
#################################################################
BINDIR = bin
OBJDIR = obj
SRCDIR = src
LIBDIR = lib
TESTDIR = tests

# the following determines whether LARC is compiled with USE_REAL, USE_INTEGER
# or USE_COMPLEX.
ifeq ($(TYPE),REAL)
   SOURCE_FILE = $(SRCDIR)/.typeREAL.h
else ifeq ($(TYPE),INTEGER)
   SOURCE_FILE = $(SRCDIR)/.typeINTEGER.h
else ifeq ($(TYPE),COMPLEX)
   SOURCE_FILE = $(SRCDIR)/.typeCOMPLEX.h
else
$(warning Using default real type...)
   SOURCE_FILE = $(SRCDIR)/.typeREAL.h
endif


OPTS = -O2 -Wall -fPIC
CFLAGS = $(OPTS) -I$(SRCDIR) -std=gnu99
LIBS = -lm -lrt -lgmp -lpthread 
CC = gcc

# The file larc_py_wrap.c is created by SWIG, and therefore may or may not
# be present in the SRCDIR. We use SRCSNOWRAP to make sure we're not expecting
# the *wrap.o file to be in the OBJDIR or in OBJS.
SRCS = $(wildcard $(SRCDIR)/*.c)
SRCSNOWRAP = $(filter-out %_wrap.c, $(SRCS))
OBJS = $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(SRCSNOWRAP))
INCS = $(wildcard $(SRCDIR)/*.h)

PYINC = -I$(PYINCDIR)
PYINCDIR := $(shell python2 -c 'import sys; print sys.prefix+"/include/python"+str(sys.version_info[0])+"."+str(sys.version_info[1])')

TEST_SRCS = $(wildcard $(TESTDIR)/*_tests.c)
TESTS = $(patsubst %.c, %, $(TEST_SRC))

BIN_LIST = exampleLARC 
BINS = $(patsubst %,$(BINDIR)/%,$(BIN_LIST))

TARGET = lib/liblarc.a
SO_TARGET = $(patsubst %.a, %.so, $(TARGET))

.PHONY: all build clean #tests install

all: NEWTYPE $(SRCS) $(BINS) $(SO_TARGET) $(SRCDIR)/_larc_py.so $(TARGET) # tests

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
	rm -f $(SO_TARGET) $(TARGET)
	rm -f $(SRCDIR)/larc_py.pyc $(SRCDIR)/_larc_py.so $(SRCDIR)/larc_py_wrap.o $(SRCDIR)/larc_py_wrap.c $(SRCDIR)/larc_py.py
	rm -f tests/tests.log
	find . -name "*.gc*" -exec rm {} \;
	rm -rf `find . -name "*.dSYM" -print`

#install: all
#	install -d $(DESTDIR)/$(PREFIX)/lib/
#	install $(TARGET) $(DESTDIR)/$(PREFIX)/lib/  

# we keep all the SWIG dependent stuff in SRCDIR to avoid problems later
$(SRCDIR)/_larc_py.so: $(SRCDIR)/larc_py_wrap.o $(OBJS)
	$(CC) $(CFLAGS) -shared $^ -o $@ $(LIBS)

$(SRCDIR)/larc_py_wrap.o: $(SRCDIR)/larc_py_wrap.c
	$(CC) $(CFLAGS) $(PYINC) -c $< -o $@

# using SRCSNOWRAP avoids a dependency loop (larc_py_wrap.c shouldn't depend
# on itself)
$(SRCDIR)/larc_py_wrap.c: $(SRCDIR)/larc_py.i $(SRCSNOWRAP) $(INCS)
	swig -python $<

# the binary executables do not depend on the python wrapper
$(BINS): $(OBJS) $(BINDIR)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS)

$(BINDIR): 
	@test -d $(BINDIR)/ || mkdir -p $(BINDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.c $(INCS)
	@test -d $(OBJDIR)/ || mkdir -p $(OBJDIR)
	$(CC) -c -o $@ $< $(CFLAGS)


#test:
#	make clean
#	make all
#	cd $(TESTDIR) && tcsh -f larc_test

