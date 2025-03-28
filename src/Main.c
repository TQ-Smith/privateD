
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
        printf("<pop1>,<pop2>,<pop3> must contain three populations. Exiting!\n");
        return NULL;
    }

    // The three population names.
    char* tok;
    char* first;
    char* second;
    char* third;

    // Get the three population names.
    tok = strtok(popList, ",");
    first = strdup(tok);
    tok = strtok(NULL, ",");
    second = strdup(tok);
    tok = strtok(NULL, ",");
    third = strdup(tok);

    // Create the file buffer to read in <sampleToPop.tsv>.
    FILE* sampleToPopFile = fopen(sampleToPopFileName, "r");
    if (sampleToPopFile == NULL) {
        printf("<sampleToPop.tsv> does not exist. Exiting!\n");
        return NULL;
    }
    
    // Used for reading in sampleToPopFile.
    char sampleName[512];
    char popName[512];

    // Setup hash table.
    khash_t(str) *h;
    khint_t k;
    h = kh_init(str);

    int numFields = 0;

    // Parse <sampleToPop.tsv>.
    int absent;
    while ( !feof(sampleToPopFile) ) {

        numFields = fscanf(sampleToPopFile, "%s\t%s\n", sampleName, popName);

        // Invalid file. Free memory and return.
        if (numFields != 2) {
            for (k = 0; k < kh_end(h); k++)
                if (kh_exist(h, k))
                    free((char*) kh_key(h, k));
            kh_destroy(str, h);
            fclose(sampleToPopFile);
            printf("<sampleToPop.tsv> improperly formatted. Exiting!\n");
            return NULL;
        }

        // If the name is not in <popList>, label it as -1.
        int popLabel = -1;
        if (strcmp(popName, first) == 0)  {  popLabel = 1;  }
        if (strcmp(popName, second) == 0) {  popLabel = 2;  }
        if (strcmp(popName, third) == 0)  {  popLabel = 3;  }

        // Insert into hash table.
        k = kh_put(str, h, sampleName, &absent);
        if (absent) kh_key(h, k) = strdup(sampleName);
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
    free(first); free(second); free(third);
    fclose(sampleToPopFile);
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
    printf("Usage: privateD [options] <inFile.vcf.gz> <sampleToPop.tsv> <pop1>,<pop2>,<pop3>\n\n");
    printf("<inFile.vcf.gz>             The input VCF file.\n");
    printf("<sampleToPop.tsv>           Tab seperate file associating each sample with a population.\n");
    printf("<pop1>,<pop2>,<pop3>        Names of the three populations in <sampleToPop.csv.gz>.\n\n");
    printf("Options:\n");
    printf("    -g                      The maximum standardized sample size used for rarefaction. Default 1.\n");
    printf("    -b                      Block size for weighted jackknife. Default 2 MB.\n");
    printf("    -h                      Haplotype size in number of loci. Default 1.\n");
    printf("    -o                      The outputbase name.\n");
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
    char* outputBasename = NULL;

    // Parse options
    ketopt_t options = KETOPT_INIT;
    int c; 
    while ((c = ketopt(&options, argc, argv, 1, "g:b:h:o:", long_options)) >= 0) {
        if (c == 'g') { g = (int) strtol(options.arg, (char**) NULL, 10); }
        else if (c == 'b') { blockSize = (int) strtol(options.arg, (char**) NULL, 10); }
        else if (c == 'h') { h = (int) strtol(options.arg, (char**) NULL, 10); }
        else if (c == 'o') {outputBasename = options.arg;}
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
    if (outputBasename == NULL) {
        printf("-o must be given an argument. Exiting!\n");
        return -1;
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
        printf("g is more than the number of chromsomes belonging to the three populations. Exiting!\n");
        free(popList);
        destroy_vcf_locus_parser(vcfFile);
        destroy_haplotype_encoder(encoder);
        return 1;
    }

    printf("Partitioning genome into blocks ...\n");
    // Block our genome.
    BlockList_t* blockList = privateD(vcfFile, encoder, samplesToLabel, maxNumOfHaps, g, blockSize, h);

    printf("Executing weight block jackknife ...\n");
    // Execute jackknife.
    weighted_block_jackknife(blockList);

    // Print results to the files.
    printf("Printing results to files ...\n");

    // Open our three files.
    kstring_t* output = init_kstring(outputBasename);
    kputs("_block_pvals.tsv", output);
    FILE* blocksPvals = fopen(ks_str(output), "w");
    ks_overwrite(outputBasename, output);
    kputs("_block_privateD.tsv", output);
    FILE* blocksD = fopen(ks_str(output), "w");
    ks_overwrite(outputBasename, output);
    kputs("_global.tsv", output);
    FILE* global = fopen(ks_str(output), "w");
    destroy_kstring(output);

    // Print header information.
    fprintf(blocksPvals, "#Command: ");
    fprintf(blocksD, "#Command: ");
    fprintf(global, "#Command: ");
    for (int i = 0; i < argc; i++) {
        fprintf(blocksPvals, "%s ", argv[i]);
        fprintf(blocksD, "%s ", argv[i]);
        fprintf(global, "%s ", argv[i]);
    }
    fprintf(blocksPvals, "\n");
    fprintf(blocksD, "\n");
    fprintf(global, "\n");
    fprintf(blocksPvals, "BlockNum\tBlockNumOnChrom\tChrom\tStart\tEnd\tNumHaps");
    fprintf(blocksD, "BlockNum\tBlockNumOnChrom\tChrom\tStart\tEnd\tNumHaps");
    fprintf(global, "Value\tNumBlocks");
    for (int i = 1; i <= g; i++) {
        fprintf(blocksPvals, "\tg%d", i);
        fprintf(blocksD, "\tg%d", i);
        fprintf(global, "\tg%d", i);
    }
    fprintf(blocksPvals, "\n");
    fprintf(blocksD, "\n");
    fprintf(blocksD, "\n");

    // Print block information.
    Block_t* curBlock = blockList -> head;
    for (int i = 0; i < blockList -> numBlocks; i++) {
        fprintf(blocksPvals, "%d\t%d\t%s\t%d\t%d\t%d", curBlock -> blockNum, curBlock -> blockNumOnChrom, curBlock -> chrom, curBlock -> startCoordinate, curBlock -> endCoordinate, curBlock -> numHaps);
        fprintf(blocksD, "%d\t%d\t%s\t%d\t%d\t%d", curBlock -> blockNum, curBlock -> blockNumOnChrom, curBlock -> chrom, curBlock -> startCoordinate, curBlock -> endCoordinate, curBlock -> numHaps);
        for (int j = 0; j < g; j++) {
            fprintf(blocksPvals, "\t%lf", curBlock -> rarefactCounts[j].p);
            fprintf(blocksPvals, "\t%lf", curBlock -> rarefactCounts[j].num / (double) curBlock -> rarefactCounts[j].denom);
        }
        fprintf(blocksPvals, "\n");
        fprintf(blocksD, "\n");
        curBlock = curBlock -> next;
    }
    // Print global information.
    fprintf(global, "privateD\t%d", blockList -> numBlocks);
    for (int i = 0; i < g; i++) {
        fprintf(global, "\t%lf", blockList -> rarefactCounts[i].num / (double) blockList -> rarefactCounts[i].denom);
    }
    fprintf(global, "\n");
    fprintf(global, "pvals\t%d", blockList -> numBlocks);
    for (int i = 0; i < g; i++) {
        fprintf(global, "\t%lf", blockList -> rarefactCounts[i].p);
    }
    fprintf(global, "\n");
    fprintf(global, "stderrs\t%d", blockList -> numBlocks);
    for (int i = 0; i < g; i++) {
        fprintf(global, "\t%lf", blockList -> stderrs[i]);
    }
    fprintf(global, "\n");

    printf("Done!\n");
    fclose(blocksPvals);
    fclose(blocksD);
    fclose(global);
    // Free all used memory.
    destroy_vcf_locus_parser(vcfFile);
    destroy_haplotype_encoder(encoder);
    free(popList);
    free(samplesToLabel);
    destroy_block_list(blockList);

    return 0;
}