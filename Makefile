
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
	$(CC) $(LFLAGS) bin/privateD src/*.o lib/*.o -lz -lm

.PHONY: src
src: src/Main.o 

src/Main.o: src/Interface.o src/HaplotypeEncoder.o src/PrivateD.o
	$(CC) $(CFLAGS) src/Main.c -o src/Main.o

src/Interface.o:
	$(CC) $(CFLAGS) src/Interface.c -o src/Interface.o

src/PrivateD.o: src/BlockList.o
	$(CC) $(CFLAGS) src/PrivateD.c -o src/PrivateD.o

src/BlockList.o:
	$(CC) $(CFLAGS) src/BlockList.c -o src/BlockList.o

src/HaplotypeEncoder.o: src/VCFLocusParser.o
	$(CC) $(CFLAGS) src/HaplotypeEncoder.c -o src/HaplotypeEncoder.o

src/VCFLocusParser.o:
	$(CC) $(CFLAGS) src/VCFLocusParser.c -o src/VCFLocusParser.o

.PHONY: lib 
lib: lib/kstring.o

lib/kstring.o:
	$(CC) $(CFLAGS) lib/kstring.c -o lib/kstring.o

.PHONY: clean
clean:
	rm src/*.o bin/* lib/*.o