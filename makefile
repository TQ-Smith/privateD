
# File: Makefile
# Date: 23 December 2024
# Author: T. Quinn Smith
# Principal Investigator: Dr. Zachary A. Szpiech
# Purpose: Build privateD.

CC?=gcc
CFLAGS = -c -Wall -g
LFLAGS = -g -o

bin/privateD: src/Main.o
	mkdir -p bin
	$(CC) $(LFLAGS) bin/privateD src/*.o -lz -lm

src/Main.o: src/HaplotypeEncoder.o src/PrivateD.o
	$(CC) $(CFLAGS) src/Main.c -o src/Main.o

src/PrivateD.o: src/BlockList.o
	$(CC) $(CFLAGS) src/PrivateD.c -o src/PrivateD.o

src/BlockList.o:
	$(CC) $(CFLAGS) src/BlockList.c -o src/BlockList.o

src/HaplotypeEncoder.o: src/VCFLocusParser.o
	$(CC) $(CFLAGS) src/HaplotypeEncoder.c -o src/HaplotypeEncoder.o

src/VCFLocusParser.o:
	$(CC) $(CFLAGS) src/VCFLocusParser.c -o src/VCFLocusParser.o

.PHONY: clean
clean:
	rm src/*.o bin/* lib/lapack/*.o