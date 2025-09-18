
# File: Makefile
# Date: 
# Author: T. Quinn Smith
# Principal Investigator: Dr. Zachary A. Szpiech
# Purpose: Build privateD.

CC?=gcc
CFLAGS = -c -Wall -g -I lib
LFLAGS = -g -o

bin/privateD: lib src
	mkdir -p bin
	$(CC) $(LFLAGS) bin/privateD src/*.o lib/*.o lib/gsl/*.o -lz -lm

.PHONY: src
src: src/Main.o 

src/Main.o: src/Interface.o src/HaplotypeEncoder.o src/PrivateD.o
	$(CC) $(CFLAGS) src/Main.c -o src/Main.o

src/Interface.o:
	$(CC) $(CFLAGS) src/Interface.c -o src/Interface.o

src/PrivateD.o: src/BlockList.o
	$(CC) $(CFLAGS) -DHAVE_INLINE src/PrivateD.c -o src/PrivateD.o

src/BlockList.o:
	$(CC) $(CFLAGS) src/BlockList.c -o src/BlockList.o

src/HaplotypeEncoder.o: src/VCFLocusParser.o
	$(CC) $(CFLAGS) src/HaplotypeEncoder.c -o src/HaplotypeEncoder.o

src/VCFLocusParser.o:
	$(CC) $(CFLAGS) src/VCFLocusParser.c -o src/VCFLocusParser.o

.PHONY: lib 
lib: lib/kstring.o lib/gsl

lib/kstring.o:
	$(CC) $(CFLAGS) lib/kstring.c -o lib/kstring.o

.PHONY: lib/gsl
lib/gsl:
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/error.c -o lib/gsl/error.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/message.c -o lib/gsl/message.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/stream.c -o lib/gsl/stream.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/default.c -o lib/gsl/default.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/rng.c -o lib/gsl/rng.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/mt.c -o lib/gsl/mt.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/types.c -o lib/gsl/types.o

.PHONY: clean
clean:
	rm src/*.o bin/* lib/*.o lib/gsl/*.o