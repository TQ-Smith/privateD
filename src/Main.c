
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
#include <math.h>

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
    for (int i = 0; i < vcfFile -> numSamples; i++)
        if (samplesToLabel[i] != -1) 
            numSamples++;
    
    // fprintf(stderr, "\nBlocking Genome ...\n");
    // Compute privateD in each block and genome-wide.
    BlockList_t* blocks = privateD(vcfFile, encoder, samplesToLabel, numSamples, config -> sampleSize, config -> blockSize, config -> haplotypeSize);
    // fprintf(stderr, "Finished Blocking Genome ...\n");

    // Calculate p-values.
    if (config -> replicates > 0) {
        // fprintf(stderr, "\nCalculating Bootstrap ...\n");
        bootstrap(blocks, config -> replicates);
        // fprintf(stderr, "Finished Bootstrap ...\n");
    } else {
        // fprintf(stderr, "\nCalculating Weighted Jackknife ...\n");
        weighted_block_jackknife(blocks);
        // fprintf(stderr, "Finished Weighted Jackknife ...\n");
    }

    // Default is to write to stdout.
    FILE* output = stdout;
    if (config -> outBaseName != NULL) {
        kstring_t* out = calloc(1, sizeof(kstring_t));
        ksprintf(out, "%s%s", config -> outBaseName, "_privateD.tsv");
        output = fopen(out -> s, "w");
        free(out -> s); free(out);
    }

    // Output values.
    // fprintf(stderr, "\nPrinting Output ...\n\n");
    fprintf(output, "#%s\n", config -> cmd);
    fprintf(output, "#Block_Num\tBlock_Num_on_Chr\tChromosome\tStart_Position\tEnd_Position\tNum_Haps\tNum_Alleles\tprivateD\tpvalue\n");
    for (Block_t* temp = blocks -> head; temp != NULL; temp = temp -> next)
        fprintf(output, "%d\t%d\t%s\t%d\t%d\t%d\t%lf\t%lf\t%lf\n", temp -> blockNum, temp -> blockNumOnChrom, temp -> chrom, temp -> startCoordinate, temp -> endCoordinate, temp -> numHaps, temp -> denominatorPrivateD, temp -> numeratorPrivateD / temp -> denominatorPrivateD, temp -> p);
    fprintf(output, "%d\t%d\t%s\t%d\t%d\t%d\t%lf\t%lf\t%lf\n", 0, 0, "Global", 0, 0, blocks -> numHaps, blocks -> denominatorPrivateD, blocks -> numeratorPrivateD / blocks -> denominatorPrivateD, blocks -> p);
    // fprintf(stderr, "Finished Printing Output ...\n\n");
    // fprintf(stderr, "Done!\n");

    // Free all used memory.
    if (config -> outBaseName != NULL)
        fclose(output);
    free(samplesToLabel);
    destroy_privated_config(config);
    destroy_vcf_locus_parser(vcfFile);
    destroy_haplotype_encoder(encoder);
    destroy_block_list(blocks);

    return 0;
}