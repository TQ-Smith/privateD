
// File: PrivateD.c
// Date: 26 December 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Compute privateD with weighted jackknife.

#include <math.h>
#include "PrivateD.h"
#include "Matrix.h"
#include "../lib/kvec.h"
#include "../lib/khash.h"

MATRIX_INIT(int, int)

// Define a block used in the jackknife.
typedef struct {
    // The numerator and denominator for D^g within the given block.
    int num;
    int denom;
} Block_t;

// Calculate D^g within the next block.
// Accepts:
//  VCFLocusParser_t* vcfFile -> The VCF file we are reading.
//  HaplotypeEncoder_t* encoder -> Encode the haplotypes from the VCF file.
//  double** hapCounts -> A table used to count haplotypes per population.
//  int g -> The max standardized sample size.
//  int blockSize -> The block size in base pairs.
//  int h -> The size of the haplotypes.
//  int* startOfNextBlock -> Sets the start basepair of next block.
// Returns: Block_t*, the next block or NULL if end of VCF.
Block_t* get_next_block(VCFLocusParser_t* vcfFile, HaplotypeEncoder_t* encoder, int** hapCounts, int g, int blockSize, int h, int* startOfNextBlock) {

    // No block to return if end of file reached.
    if (vcfFile -> isEOF)
        return NULL;

    Block_t* temp = calloc(1, sizeof(Block_t));

    if (*startOfNextBlock < vcfFile -> nextCoord)
        *startOfNextBlock = (int) ((vcfFile -> nextCoord - 1) / blockSize) * blockSize + 1;

    int endOfBlock = *startOfNextBlock + blockSize - 1;

    bool isOnSameChrom = true;

    printf("Block %d to ", vcfFile -> nextCoord);
    while ((isOnSameChrom = get_next_haplotype(vcfFile, encoder, h)) && vcfFile -> nextCoord <= endOfBlock) {

    }
    printf("%d\n", encoder -> endCoord);

    if (!isOnSameChrom)
        *startOfNextBlock = 1;
    else 
        *startOfNextBlock = endOfBlock + 1;

    return temp;

}

void privateD(VCFLocusParser_t* vcfFile, HaplotypeEncoder_t* encoder, int* samplesToLabel, int g, int blockSize, int h, double* D, double* stdError, double* pvals) {

    // The max number of unique haplotypes is 2 * the number of samples in pops 1,2,3.
    int maxNumberOfHaps = 0;
    for (int i = 0; i < vcfFile -> numSamples; i++)
        if (samplesToLabel[i] != -1)
            maxNumberOfHaps++;
    maxNumberOfHaps *= 2;

    // Our haplotype count table for rarefaction.
    int** hapCounts = create_int_matrix(maxNumberOfHaps, 3);

    // Our genome-wide block.
    Block_t* global = calloc(g, sizeof(Block_t));

    // Our list of blocks.
    kvec_t(Block_t*) blocks;
	kv_init(blocks);

    // Accumulate our blocks.
    int startOfNextBlock = 1;
    Block_t* temp = NULL;
    while ((temp = get_next_block(vcfFile, encoder, hapCounts, g, blockSize, h, &startOfNextBlock)) != NULL) {
        kv_push(Block_t*, blocks, temp);
        for (int i = 0; i < g; i++) {
            global[i].num += temp -> num;
            global[i].denom += temp -> denom;
        }
    }

    // Free all used memory.
    destroy_int_matrix(hapCounts, maxNumberOfHaps);
    free(global);
    for (int i = 0; i < kv_size(blocks); i++)
        free(kv_A(blocks, i));
    kv_destroy(blocks);
}