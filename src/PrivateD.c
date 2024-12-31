
// File: PrivateD.c
// Date: 26 December 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Compute privateD with weighted jackknife.

#include <math.h>
#include "PrivateD.h"
#include "Matrix.h"
#include "../lib/kvec.h"

// Our epsilong for floating point comparison.
#define EPS 1e-10

MATRIX_INIT(int, int)

// Define a block used in the jackknife.
typedef struct {
    // The numerator and denominator for D^g within the given block.
    int num;
    int denom;
} Block_t;

// Calculate Equation (1) from Szpiech et al. 2008.
// Accepts:
//  int Nj -> The size of pop j.
//  int Nji -> The number of copies of hap i in pop j.
//  int g -> The standardized sample size.
// Returns: double, the prob of not observing i in pop j from a sample of size g.
static inline double Q_ijg(int Nj, int Nji, int g) {
    if (Nj - Nji < g)
        return 0;
    double result = 0;
    for (int u = 0; u <= g - 1; u++) {
        result += log((Nj - Nji - u) / (double) (Nj - u));
    }
    return exp(result);
}

// Calculate contribution of locus to privateD.
// Accepts:
//  Block_t* block -> The block the locus is in.
//  int** hapCounts -> The haplotype counts used for rarefaction.
//  int numUniqueHaps -> The number of unique haplotypes in hapCounts.
//  int G -> The standardized sample size.
// Returns: void.
void locus_privateD(Block_t* block, int** hapCounts, int numUniqueHaps, int G) {

    // All of our variables used in the rarefaction computation.
    double pi13, pi23, temp;

    // We do this for all our standardized sample sizes.
    for (int g = 1; g <= G; g++) {
        pi13 = 0;
        pi23 = 0;
        // For each haplotype at the locus.
        for (int i = 1; i <= numUniqueHaps; i++) {
            temp = log(1 - Q_ijg(hapCounts[0][1], hapCounts[i][1], g)) + log(1 - Q_ijg(hapCounts[0][2], hapCounts[i][2], g)) + log(Q_ijg(hapCounts[0][0], hapCounts[i][0], g));
            pi13 += exp(temp);
            temp = log(1 - Q_ijg(hapCounts[0][0], hapCounts[i][0], g)) + log(1 - Q_ijg(hapCounts[0][2], hapCounts[i][2], g)) + log(Q_ijg(hapCounts[0][1], hapCounts[i][1], g));
            pi23 += exp(temp);
        }
        // If pi13 and pi23 are NOT equal, the locus contributes to D^g.
        if (fabs(pi13 - pi23) > EPS) {
            if (pi23 > pi13)
                block[g - 1].num += 1;
            else 
                block[g - 1].num += 1;
            block[g - 1].denom += 1;
        }
    }
 
}

// Count number of haplotypes to each population.
// Accepts:
//  HaplotypeEncoder_t* encoder -> Encodes the haplotypes.
//  int* sampleIndices -> The indices of the samples in the three populations.
//  int* samplesToLabel -> The population labels of the samples.
//  int** hapCounts -> The haplotype counts used for rarefaction.
//  int numSamplesInPops -> The number of samples in the three populations.
// Returns: int, the number of unique haplotypes.
int accumulate_counts(HaplotypeEncoder_t* encoder, int* sampleIndices, int* samplesToLabel, int** hapCounts, int numSamplesInPops) {
    
    // We use the encoder's hash table to do the relabeling.
    
    // The number of unique haplotypes.
    int numUniqueHaps = 0;

    khiter_t k;
    int ret;

    // Reset previous labels in hash table.
    for (int i = 0; i < numSamplesInPops; i++) {
        k = kh_get(haplotype, encoder -> labelMap, encoder -> genotypes[sampleIndices[i]].left);
        if (k == kh_end(encoder -> labelMap)) {
            k = kh_put(haplotype, encoder -> labelMap, encoder -> genotypes[sampleIndices[i]].left, &ret);
        }
        kh_value(encoder -> labelMap, k) = MISSING;
        k = kh_get(haplotype, encoder -> labelMap, encoder -> genotypes[sampleIndices[i]].right);
        if (k == kh_end(encoder -> labelMap)) {
            k = kh_put(haplotype, encoder -> labelMap, encoder -> genotypes[sampleIndices[i]].right, &ret);
        }
        kh_value(encoder -> labelMap, k) = MISSING;
    }

    // Zero totals for the three pops.
    hapCounts[0][0] = hapCounts[0][1] = hapCounts[0][2] = 0;

    // We relabel and populate hapCounts.
    for (int i = 0; i < numSamplesInPops; i++) {
        // If the genotype is missing, then we skip relabeling.
        if (encoder -> genotypes[sampleIndices[i]].left != MISSING) {
            k = kh_get(haplotype, encoder -> labelMap, encoder -> genotypes[sampleIndices[i]].left);
            // If the hap label needs to be relabeled, then zero counts and assign label to hap.
            if (kh_value(encoder -> labelMap, k) == MISSING) {
                kh_value(encoder -> labelMap, k) = numUniqueHaps;
                hapCounts[numUniqueHaps + 1][0] = hapCounts[numUniqueHaps + 1][1] = hapCounts[numUniqueHaps + 1][2] = 0;
                // Generate next label.
                numUniqueHaps++;
            }
            // Increment counter.
            hapCounts[kh_value(encoder -> labelMap, k) + 1][samplesToLabel[sampleIndices[i]] - 1]++;
            // Increment number in that population.
            hapCounts[0][samplesToLabel[sampleIndices[i]] - 1]++;
        }
        if (encoder -> genotypes[sampleIndices[i]].right != MISSING) {
            k = kh_get(haplotype, encoder -> labelMap, encoder -> genotypes[sampleIndices[i]].left);
            if (kh_value(encoder -> labelMap, k) == MISSING) {
                kh_value(encoder -> labelMap, k) = numUniqueHaps;
                hapCounts[numUniqueHaps + 1][0] = hapCounts[numUniqueHaps + 1][1] = hapCounts[numUniqueHaps + 1][2] = 0;
                numUniqueHaps++;
            }
            hapCounts[kh_value(encoder -> labelMap, k) + 1][samplesToLabel[sampleIndices[i]] - 1]++;
            hapCounts[0][samplesToLabel[sampleIndices[i]] - 1]++;
        }
    }

    return numUniqueHaps;

}

// Calculate D^g within the next block.
// Accepts:
//  VCFLocusParser_t* vcfFile -> The VCF file we are reading.
//  HaplotypeEncoder_t* encoder -> Encode the haplotypes from the VCF file.
//  int* sampleIndices -> The indices of the samples in the three populations.
//  int* samplesToLabel -> Associates sample index with population label.
//  int** hapCounts -> A table used to count haplotypes per population.
//  int numSamplesInPops -> The number of samples in the three pops.
//  int g -> The max standardized sample size.
//  int blockSize -> The block size in base pairs.
//  int h -> The size of the haplotypes.
//  int* startOfNextBlock -> Sets the start base pair of next block.
// Returns: Block_t*, the next block or NULL if end of VCF.
Block_t* get_next_block(VCFLocusParser_t* vcfFile, HaplotypeEncoder_t* encoder, int* sampleIndices, int* samplesToLabel, int** hapCounts, int numSamplesInPops, int g, int blockSize, int h, int* startOfNextBlock) {

    // No block to return if end of file reached.
    if (vcfFile -> isEOF)
        return NULL;
    
    // Create our block.
    Block_t* temp = calloc(g, sizeof(Block_t));

    // Find start of the block the next record is in.
    if (*startOfNextBlock < vcfFile -> nextCoord)
        *startOfNextBlock = (int) ((vcfFile -> nextCoord - 1) / blockSize) * blockSize + 1;
    
    // End of the current block.
    int endOfBlock = *startOfNextBlock + blockSize - 1;

    // While we read in records within the block.
    bool isOnSameChrom = true;
    while (isOnSameChrom && vcfFile -> nextCoord <= endOfBlock) {
        isOnSameChrom = get_next_haplotype(vcfFile, encoder, h);
        // QUESTION: There is a possibility that a D^g is 0 within the block.
        //  Should we exclude these blocks? Ask Zach at some point.

        // Create table of haplotype counts.
        int numUniqueHaps = accumulate_counts(encoder, sampleIndices, samplesToLabel, hapCounts, numSamplesInPops);

        // Calculate privateD for locus.
        locus_privateD(temp, hapCounts, numUniqueHaps, g);
    }

    printf("D =");
    for (int i = 0; i < g; i++)
        printf("\t%lf", temp[i].num / (double) temp[i].denom);
    printf("\n");

    // Set the start position of the next block.
    if (!isOnSameChrom)
        *startOfNextBlock = 1;
    else 
        *startOfNextBlock = endOfBlock + 1;

    // Return the block.
    return temp;

}

void privateD(VCFLocusParser_t* vcfFile, HaplotypeEncoder_t* encoder, int* samplesToLabel, int maxNumOfHaps, int g, int blockSize, int h, double* D, double* stdError, double* pvals) {

    // Indices of our samples in the three populations.
    int* sampleIndices = calloc(maxNumOfHaps / 2, sizeof(int));
    int counter = 0;
    for (int i = 0; i < vcfFile -> numSamples; i++) {
        if (samplesToLabel[i] != -1) {
            sampleIndices[counter] = i;
            counter++;
        }
    }

    // Our haplotype count table for rarefaction.
    //  NOTE: The first row contains the totals for each population.
    //        This is meant to accomodate missing data.
    int** hapCounts = create_int_matrix(maxNumOfHaps + 1, 3);

    // Our genome-wide block.
    Block_t* global = calloc(g, sizeof(Block_t));

    // Our list of blocks.
    kvec_t(Block_t*) blocks;
	kv_init(blocks);

    // Accumulate our blocks.
    int startOfNextBlock = 1;
    Block_t* temp = NULL;
    while ((temp = get_next_block(vcfFile, encoder, sampleIndices, samplesToLabel, hapCounts, maxNumOfHaps / 2, g, blockSize, h, &startOfNextBlock)) != NULL) {
        kv_push(Block_t*, blocks, temp);
        for (int i = 0; i < g; i++) {
            global[i].num += temp -> num;
            global[i].denom += temp -> denom;
        }
    }

    // Free all used memory.
    free(sampleIndices);
    destroy_int_matrix(hapCounts, maxNumOfHaps + 1);
    free(global);
    for (int i = 0; i < kv_size(blocks); i++)
        free(kv_A(blocks, i));
    kv_destroy(blocks);
}