
// File: Main.c
// Date: 23 December 2024
// Version 1: 19 June 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Three population introgression test using private allelic richness.

#include "Interface.h"
#include "VCFLocusParser.h"
#include "HaplotypeEncoder.h"
#include "PrivateD.h"
#include "kstring.h"

// Create a string-to-string hash table.
KHASH_MAP_INIT_STR(str, int)

int* labelSamples(char** sampleNames, int numSamples, char* samplesToPopFileName, char* threePopList) {
    
    // Check that three population names were supplied.
    int n = 0;
    for (int i = 0; i < strlen(threePopList); i++)
        if (threePopList[i] == ',')
            n++;
    if (n != 2)
        return NULL;

    // The three population names.
    char* tok = NULL;
    char* first = NULL;
    char* second = NULL;
    char* third = NULL;

    // Get the three population names.
    tok = strtok(threePopList, ",");
    first = strdup(tok);
    tok = strtok(NULL, ",");
    second = strdup(tok);
    tok = strtok(NULL, ",");
    third = strdup(tok);

    FILE* samplesToPopFile = fopen(samplesToPopFileName, "r");

    // Used for reading in sampleToPopFile.
    char* line = calloc(1024, sizeof(char));
    char* sampleName = calloc(512, sizeof(char));
    char* popName = calloc(512, sizeof(char));

    // Setup hash table.
    khash_t(str) *h;
    khint_t k;
    h = kh_init(str);

    // Parse <sampleToPop.tsv>.
    int absent;
    size_t length = 0;
    while (getline(&line, &length, samplesToPopFile) != -1) {
        
        int numFields = sscanf(line, "%s\t%s\n", sampleName, popName);

        // If invalid formatting, free all memory and exit.
        if (numFields != 2) {
            free(line); free(sampleName); free(popName);
            free(first); free(second); free(third);
            fclose(samplesToPopFile);
            for (k = 0; k < kh_end(h); k++)
                if (kh_exist(h, k))
                    free((char*) kh_key(h, k));
            kh_destroy(str, h);
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
    if (line)
        free(line); 
    free(sampleName); free(popName);

    // Create association array.
    int* samplesToLabel = calloc(numSamples, sizeof(int));
    for (int i = 0; i < numSamples; i++) {
        k = kh_get(str, h, sampleNames[i]);
        // Assign -1 if sample name does not have a pop label.
        if (k == kh_end(h))
            samplesToLabel[i] = -1;
        else
            samplesToLabel[i] = kh_val(h, k);
    }

    // Free used memory.
    free(first); free(second); free(third);
    fclose(samplesToPopFile);
    for (k = 0; k < kh_end(h); k++)
        if (kh_exist(h, k))
            free((char*) kh_key(h, k));
    kh_destroy(str, h);

    return samplesToLabel;
}

int main (int argc, char *argv[]) {

    // Print help if no arguments were given.
    if (argc == 1) {
        print_help();
        return 0;
    }

    // Get configuration from user and exit if error.
    PrivateDConfig_t* config = init_privated_config(argc, argv);
    if (config == NULL)
        return -1;

    // Create the VCF parser and haplotype encoer.
    VCFLocusParser_t* vcfFile = init_vcf_locus_parser(config -> inputFileName, config -> MAF, config -> missingAF, true);
    HaplotypeEncoder_t* encoder = init_haplotype_encoder(vcfFile -> numSamples);

    // Associate the index of a sample in the VCF file with a population number 1,2,3. 
    //  Samples not used are marked with -1.
    int* samplesToLabel = labelSamples(vcfFile -> sampleNames, vcfFile -> numSamples, config -> samplesToPopFileName, config -> threePopList);
    if (samplesToLabel == NULL) {
        fprintf(stderr, "Improper formatting to assign samples to populations. Exiting!\n");
        destroy_privated_config(config);
        destroy_vcf_locus_parser(vcfFile);
        destroy_haplotype_encoder(encoder);
        return -1;
    }

    // Count the number of samples in the three populations.
    int numSamples = 0;
    for (int i = 0; i < vcfFile -> numSamples; i++) {
        if (samplesToLabel[i] != -1) {
            numSamples++;
        }
    }

    // Check to make sure we hace enough samples.
    if (config -> sampleSize > 2 * numSamples) {
        fprintf(stderr, "--sampleSize exceeds the number of lineages in the VCF. Exiting!\n");
        free(samplesToLabel);
        destroy_privated_config(config);
        destroy_vcf_locus_parser(vcfFile);
        destroy_haplotype_encoder(encoder);
        return -1;
    }

    // Compute privateD in each block and genome-wide.
    BlockList_t* blocks = privateD(vcfFile, encoder, samplesToLabel, numSamples, config -> sampleSize, config -> blockSize, config -> haplotypeSize);

    // Execute our weighted jackknife.
    weighted_block_jackknife(blocks);

    // Open our three output files.
    kstring_t* out = calloc(1, sizeof(kstring_t));
    ksprintf(out, "%s%s", config -> outBaseName, "_global.tsv");
    FILE* global = fopen(out -> s, "w");
    free(out -> s); free(out);
    out = calloc(1, sizeof(kstring_t));
    ksprintf(out, "%s%s", config -> outBaseName, "_pvals.tsv");
    FILE* blockPVals = fopen(out -> s, "w");
    free(out -> s); free(out);
    out = calloc(1, sizeof(kstring_t));
    ksprintf(out, "%s%s", config -> outBaseName, "_dvals.tsv");
    FILE* blockDVals = fopen(out -> s, "w");
    free(out -> s); free(out);

    // Print command used in header for convinence.
    fprintf(global, "#%s\n", config -> cmd);
    fprintf(blockPVals, "#%s\n", config -> cmd);
    fprintf(blockDVals, "#%s\n", config -> cmd);

    // Global contains a line for each of the number of haplotypes, privateD, variance, and pvalues genome-wide.
    fprintf(global, "Obs");
    for (int i = 0; i < blocks -> sampleSize; i++)
        fprintf(global, "\tg=%d", i + 1);
    fprintf(global, "\n");
    fprintf(global, "Num Haps");
    for (int i = 0; i < blocks -> sampleSize; i++)
        fprintf(global, "\t%d", blocks -> rarefactCounts[i].denom);
    fprintf(global, "\n");
    fprintf(global, "privateD");
    for (int i = 0; i < blocks -> sampleSize; i++)
        fprintf(global, "\t%lf", blocks -> rarefactCounts[i].num / (double) blocks -> rarefactCounts[i].denom);
    fprintf(global, "\n");
    fprintf(global, "Variance");
    for (int i = 0; i < blocks -> sampleSize; i++)
        fprintf(global, "\t%lf", blocks -> stderrs[i] * blocks -> stderrs[i]);
    fprintf(global, "\n");
    fprintf(global, "P-Vals");
    for (int i = 0; i < blocks -> sampleSize; i++)
        fprintf(global, "\t%lf", blocks -> rarefactCounts[i].p);
    fprintf(global, "\n");

    // Print the pvalues for each block.
    fprintf(blockPVals, "Block Num\tBlock Num on Chr\tChromosome\tStart Position\tEnd Position\tNum Haps");
    for (int i = 0; i < blocks -> sampleSize; i++)
        fprintf(blockPVals, "\tg=%d", i + 1);
    fprintf(blockPVals, "\n");
    for (Block_t* temp = blocks -> head; temp != NULL; temp = temp -> next) {
        fprintf(blockPVals, "%d\t%d\t%s\t%d\t%d\t%d", temp -> blockNum, temp -> blockNumOnChrom, temp -> chrom, temp -> startCoordinate, temp -> endCoordinate, temp -> numHaps);
        for (int i = 0; i < blocks -> sampleSize; i++)
            fprintf(blockPVals, "\t%lf", temp -> rarefactCounts[i].p);
        fprintf(blockPVals, "\n");
    }

    // Print the privateD values for each block.
    fprintf(blockDVals, "Block Num\tBlock Num on Chr\tChromosome\tStart Position\tEnd Position\tNum Haps");
    for (int i = 0; i < blocks -> sampleSize; i++)
        fprintf(blockDVals, "\tg=%d", i + 1);
    fprintf(blockDVals, "\n");
    for (Block_t* temp = blocks -> head; temp != NULL; temp = temp -> next) {
        fprintf(blockDVals, "%d\t%d\t%s\t%d\t%d\t%d", temp -> blockNum, temp -> blockNumOnChrom, temp -> chrom, temp -> startCoordinate, temp -> endCoordinate, temp -> numHaps);
        for (int i = 0; i < blocks -> sampleSize; i++)
            fprintf(blockDVals, "\t%lf", temp -> rarefactCounts[i].num / (double) temp -> rarefactCounts[i].denom);
        fprintf(blockDVals, "\n");
    }

    // Free all used memory.
    fclose(global);
    fclose(blockPVals);
    fclose(blockDVals);
    free(samplesToLabel);
    destroy_privated_config(config);
    destroy_vcf_locus_parser(vcfFile);
    destroy_haplotype_encoder(encoder);
    destroy_block_list(blocks);

    return 0;
}