##                   Makefile
###################################################################
 #                                                                #
 # Copyright (C) 2014-2024, Institute for Defense Analyses        #
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
 #   - Mark Pleszkoch (IDA-CCS)                                   #
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
LARCDIR = ../..
LARCSRCDIR = $(LARCDIR)/src
LARCLIBDIR = $(LARCDIR)/lib

# Get the local paths 
include $(LARCDIR)/local/Makefile.conf
#    GMPIDIR (gmp-6 include),
#    GMPLDIR (gmp-6 lib),
#    MPIDIR (Anaconda3 include),
#    MPLDIR (Anaconda3 lib)

OPTS = -O2 -g -pg -Wall -fPIC
CFLAGS = $(OPTS) -I$(LARCSRCDIR) -I$(MPIDIR) -I$(GMPIDIR) -std=gnu99
LIBS = -L$(LARCLIBDIR) -L$(MPLDIR) -L$(GMPLDIR) -lcurses -lm -llarc -lrt -lmpc -lmpfr -lgmp -lpthread -ltinfo
CC = gcc

.PHONY: clean

test1: test1.o
	$(CC) $(CFLAGS) test1.o -o $@ $(LIBS) -Wl,-rpath=$(MPLDIR):$(GMPLDIR):$(LARCLIBDIR)

test1.o: test1.c
	$(CC) -c $< $(CFLAGS)

test2: test2.o
	$(CC) $(CFLAGS) test2.o -o $@ $(LIBS) -Wl,-rpath=$(MPLDIR):$(GMPLDIR):$(LARCLIBDIR)

test2.o: test2.c
	$(CC) -c $< $(CFLAGS)

clean:
	rm -f *.o test1 test2
