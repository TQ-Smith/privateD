
// File: main.c
// Date: 23 December 2024
// Version 1: 
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Three population introgression test using private allelic richness.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "VCFLocusParser.h"
#include "HaplotypeEncoder.h"
#include "../lib/ketopt.h"
#include "../lib/klist.h"
#include "../lib/kstring.h"

// Create a string-to-string hash table.
KHASH_SET_INIT_STR(str)

// Create an integer array associating sample name to population name.
int* labelSamples() {
    gzFile sampleToPopFile = gzopen(argv[2], "r");
    int errnum;
    gzerror(sampleToPopFile, &errnum);
    if (errnum != Z_OK) {
        printf("<sampleToPop.csv> does not exist or not compressed using gzip. Exiting!\n");
        return 1;
    }
    // Create the file buffer to read in <sampleToPop.csv.gz>.
    kstream_t* stream = ks_init(sampleToPopFile);
    kstring_t* buffer = init_kstring(NULL);
    
    // Get the three population names.
    int n;
    int* fields = ksplit(s, ',', &n);
	for (int i = 0; i < n; i++) {
		
    }

    gzclose(sampleToPopFile);
    ks_destroy(stream);
    destroy_kstring(buffer);
}

// Print the help menu for privateD.
// Accepts: void.
// Returns: void.
void print_help() {
    printf("\n");
    printf("privateD v1.0 January 2025\n");
    printf("----------------------\n\n");
    printf("Written by T. Quinn Smith\n");
    printf("Principal Investigator: Zachary A. Szpiech\n");
    printf("The Pennsylvania State University\n\n");
    printf("Usage: privateD [options] <inFile.vcf.gz> <sampleToPop.csv.gz> <popList>\n\n");
    printf("<inFile.vcf.gz>             The input VCF file.\n");
    printf("<sampleToPop.csv.gz>        Comma seperate file associating each sample with a population.\n");
    printf("<popList>                   Names of the three populations to test seperated by commas.\n\n");
    printf("Options:\n");
    printf("    -g                      The maximum standardized sample size used for rarefaction. Default 2.\n");
    printf("    -b                      Block size for jackknife. Default 10,000 KB.\n");
    printf("    -h                      Haplotype size in number of loci. Default 1.\n");
    printf("\n");
}

// Long options.
static ko_longopt_t long_options[] = { {NULL, 0, 0} };

int main (int argc, char *argv[]) {
    
    // If no arguments are given, print the help menu.
    if (argc == 1) {
        print_help();
        return 0;
    }

    // Default values for options.
    int g = 2;
    int blockSize = 10000;
    int h = 1;

    // Parse options
    ketopt_t options = KETOPT_INIT;
    int c; 
    while ((c = ketopt(&options, argc, argv, 1, "g:b:h:", long_options)) >= 0) {
        if (c == 'g') { g = (int) strtol(options.arg, (char**) NULL, 10); }
        if (c == 'b') { blockSize = (int) strtol(options.arg, (char**) NULL, 10); }
        if (c == 'h') { h = (int) strtol(options.arg, (char**) NULL, 10); }
        if (c == '?') { printf("Error! \"%s\" is unknown! Exiting ...\n", argv[options.i - 1]); return 1; }
	}

    // Check options.
    if (g < 2) {
        printf("g must be >= 2. Exiting!\n");
        return 1;
    }
    if (blockSize < 1) {
        printf("b must be >= 1. Exiting!\n");
        return 1;
    }
    if (h < 1) {
        printf("h must be >= 1. Exiting!\n");
        return 1;
    }

    // Open VCF file for reading. If valid, create the haplotype encoder.
    VCFLocusParser_t* vcfFile = init_vcf_locus_parser(argv[1], 0, 1, false);
    if (vcfFile == NULL) {
        printf("Supplied VCF does not exist. Exiting!\n");
        return 1;
    }
    HaplotypeEncoder_t* encoder = init_haplotype_encoder(vcfFile -> numSamples);

    // Free all used memory.
    destroy_vcf_locus_parser(vcfFile);
    destroy_haplotype_encoder(encoder);

    return 0;
}