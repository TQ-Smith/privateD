
# File: Makefile
# Date: 8 July 2025
# Author: T. Quinn Smith
# Principal Investigator: Dr. Zachary A. Szpiech
# Purpose: Build dplus.

CC?=gcc
CFLAGS = -c -Wall -g -I lib
LFLAGS = -g -o

bin/dplus: lib src
	mkdir -p bin
	$(CC) $(LFLAGS) bin/dplus src/*.o lib/*.o -lz -lm

.PHONY: src
src: src/Main.o 

src/Main.o: src/Interface.o src/VCFLocusParser.o src/Dplus.o
	$(CC) $(CFLAGS) src/Main.c -o src/Main.o

src/Interface.o:
	$(CC) $(CFLAGS) src/Interface.c -o src/Interface.o

src/Dplus.o: src/BlockList.o
	$(CC) $(CFLAGS) src/Dplus.c -o src/Dplus.o

src/BlockList.o:
	$(CC) $(CFLAGS) src/BlockList.c -o src/BlockList.o

src/VCFLocusParser.o:
	$(CC) $(CFLAGS) src/VCFLocusParser.c -o src/VCFLocusParser.o

.PHONY: lib 
lib: lib/kstring.o

lib/kstring.o:
	$(CC) $(CFLAGS) lib/kstring.c -o lib/kstring.o

.PHONY: clean
clean:
	rm src/*.o bin/* lib/*.o