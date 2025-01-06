
// File: Main.c
// Date: 23 December 2024
// Version 1: 6 January 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Three population introgression test using private allelic richness.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "VCFLocusParser.h"
#include "HaplotypeEncoder.h"
#include "PrivateD.h"
#include "../lib/ketopt.h"
#include "../lib/kstring.h"

// Create a string-to-string hash table.
KHASH_MAP_INIT_STR(str, int)

// Create an integer array associating sample index to population label.
// Accepts:
//  kstring_t** sampleNames -> The sample names from the VCF file.
//  int numSamples -> The number of samples in the VCF file.
//  char* sampleToPopFileName -> The file associating samples to populations.
//  char* popList -> List of the three populations we are testing.
// Returns: int*, An array associating the sample index to the population they belong from popList.
//                  Labels are 1,2,3 for the first, second, third population in popList. -1 for not in popList.
int* labelSamples(kstring_t** sampleNames, int numSamples, char* sampleToPopFileName, char* popList) {

    // Check that three population names were supplied.
    int n = 0;
    for (int i = 0; i < strlen(popList); i++)
        if (popList[i] == ',')
            n++;
    if (n != 2) {
        printf("<popList> must contain three populations. Exiting!\n");
        return NULL;
    }

    // Get the three population names.
    char* token = strtok(popList, ",");
    kstring_t* first = init_kstring(token);
    token = strtok(NULL, ",");
    kstring_t* second = init_kstring(token);
    token = strtok(NULL, ",");
    kstring_t* third = init_kstring(token);

    // Create the file buffer to read in <sampleToPop.csv.gz>.
    gzFile sampleToPopFile = gzopen(sampleToPopFileName, "r");
    int errnum;
    gzerror(sampleToPopFile, &errnum);
    if (errnum != Z_OK) {
        destroy_kstring(first); destroy_kstring(second); destroy_kstring(third);
        gzclose(sampleToPopFile);
        printf("<sampleToPop.csv> does not exist or not compressed using gzip. Exiting!\n");
        return NULL;
    }
    kstream_t* stream = ks_init(sampleToPopFile);
    kstring_t* buffer = init_kstring(NULL);
    kstring_t* sampleName = init_kstring(NULL);
    kstring_t* popName = init_kstring(NULL);

    // Setup hash table.
    khash_t(str) *h;
    khint_t k;
    h = kh_init(str);

    // Parse <sampleToPop.csv>.
    int dret, absent, commaPosition;
    while ( ks_getuntil(stream, '\n', buffer, &dret) && dret != 0 ) {

        // Find the comma on the line.
        commaPosition = -1;
        for (int i = 0; i < ks_len(buffer); i++)
            if(ks_str(buffer)[i] == ',')
                commaPosition = i; 
            
        // If the comma does not exist, throw error and exit.
        if (commaPosition == -1) {
            destroy_kstring(first); destroy_kstring(second); destroy_kstring(third); destroy_kstring(sampleName); destroy_kstring(popName);
            gzclose(sampleToPopFile);
            ks_destroy(stream);
            destroy_kstring(buffer);
            kh_destroy(str, h);
            printf("<sampleToPop.csv> not formatted properly. Exiting!\n");
            return NULL;
        }

        // Insert (sampleName, popName) into the hash table.
        ks_overwriten(ks_str(buffer), commaPosition, sampleName);
        ks_overwrite(ks_str(buffer) + commaPosition + 1, popName);

        // If the name is not in <popList>, label it as -1.
        int popLabel = -1;
        if (strcmp(ks_str(popName), ks_str(first)) == 0)  {  popLabel = 1;  }
        if (strcmp(ks_str(popName), ks_str(second)) == 0) {  popLabel = 2;  }
        if (strcmp(ks_str(popName), ks_str(third)) == 0)  {  popLabel = 3;  }

        k = kh_put(str, h, ks_str(sampleName), &absent);
        if (absent) kh_key(h, k) = strdup(ks_str(sampleName));
        kh_val(h, k) = popLabel;

    }

    // Create association array.
    int* samplesToLabel = calloc(numSamples, sizeof(int));
    for (int i = 0; i < numSamples; i++) {
        k = kh_get(str, h, ks_str(sampleNames[i]));
        // Assign -1 if sample name does not have a pop label.
        if (k == kh_end(h))
            samplesToLabel[i] = -1;
        else
            samplesToLabel[i] = kh_val(h, k);
    }

    // Free used memory.
    destroy_kstring(first); destroy_kstring(second); destroy_kstring(third); destroy_kstring(sampleName); destroy_kstring(popName);
    gzclose(sampleToPopFile);
    ks_destroy(stream);
    destroy_kstring(buffer);
    for (k = 0; k < kh_end(h); k++)
        if (kh_exist(h, k))
            free((char*) kh_key(h, k));
    kh_destroy(str, h);

    return samplesToLabel;
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
    printf("    -g                      The maximum standardized sample size used for rarefaction. Default 1.\n");
    printf("    -b                      Block size for weighted jackknife. Default 2 MB.\n");
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
    int g = 1;
    int blockSize = 2000000;
    int h = 1;

    // Parse options
    ketopt_t options = KETOPT_INIT;
    int c; 
    while ((c = ketopt(&options, argc, argv, 1, "g:b:h:", long_options)) >= 0) {
        if (c == 'g') { g = (int) strtol(options.arg, (char**) NULL, 10); }
        else if (c == 'b') { blockSize = (int) strtol(options.arg, (char**) NULL, 10); }
        else if (c == 'h') { h = (int) strtol(options.arg, (char**) NULL, 10); }
        else { printf("Error! \"%s\" is unknown! Exiting ...\n", argv[options.i - 1]); return 1; }
	}

    // Check options.
    if (g < 1) {
        printf("g must be >= 1. Exiting!\n");
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
    // QUESTION: Can we drop monomorphic sites/haplotypes? We can for sites.
    VCFLocusParser_t* vcfFile = init_vcf_locus_parser(argv[argc - 3], 0, 1, true);
    if (vcfFile == NULL) {
        printf("Supplied VCF does not exist. Exiting!\n");
        return 1;
    }
    HaplotypeEncoder_t* encoder = init_haplotype_encoder(vcfFile -> numSamples);

    // Get array associating each sample with a population label.
    char* popList = strdup(argv[argc - 1]);
    int* samplesToLabel = labelSamples(vcfFile -> sampleNames, vcfFile -> numSamples, argv[argc - 2], popList);
    if (samplesToLabel == NULL) {
        free(popList);
        destroy_vcf_locus_parser(vcfFile);
        destroy_haplotype_encoder(encoder);
        return 1;
    }

    // The number of haplotypes belonging to samples in the three populations.
    int maxNumOfHaps = 0;
    for (int i = 0; i < vcfFile -> numSamples; i++)
        if (samplesToLabel[i] != -1)
            maxNumOfHaps += 2;
    
    // Make sure we have a sufficient number of chromosomes for g.
    if (g > maxNumOfHaps) {
        printf("g is less than the number of chromosomes belonging to the three populations. Exiting!\n");
        free(popList);
        destroy_vcf_locus_parser(vcfFile);
        destroy_haplotype_encoder(encoder);
        return 1;
    }

    // Allocate memory for resulting values.
    double* D = calloc(g, sizeof(double));
    double* stdError = calloc(g, sizeof(double));
    double* pvals = calloc(g, sizeof(double));

    // Calculate our statistics.
    privateD(vcfFile, encoder, samplesToLabel, maxNumOfHaps, g, blockSize, h, D, stdError, pvals);

    // Echo command.
    for (int i = 0; i < argc; i++) {
        printf("%s ", argv[i]);
    }
    // Print results.
    printf("\n%5s\t%10s\t%10s\t%10s\n", "g", "D^g", "stderr", "p-vals");
    for (int i = 0; i < g; i++) {
        printf("%5d\t%10.5f\t%10.5f\t%10.5f\n", i + 1, D[i], stdError[i], pvals[i]);
    }

    // Free all used memory.
    destroy_vcf_locus_parser(vcfFile);
    destroy_haplotype_encoder(encoder);
    free(popList);
    free(samplesToLabel);
    free(D);
    free(stdError);
    free(pvals);

    return 0;
}