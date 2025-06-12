
// File: PrivateD.c
// Date: 26 December 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Compute privateD with weighted block jackknife.

#include "PrivateD.h"
#include <math.h>
#include <stdio.h>

/*
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
*/

void locus_privateD(Block_t* block, int** hapCounts, int numUniqueHaps, int sampleSize) {
 
}

int accumulate_counts(HaplotypeEncoder_t* encoder, int haplotypeSize, int* samplesToLabel, int** hapCounts, int numSamples) {

    int numUniqueHaps = 0;

    khiter_t k;
    int ret;

    for (int i = 0; i < encoder -> numSamples; i++) {
        if (samplesToLabel[i] != -1) {
            k = kh_get(haplotype, encoder -> labelMap, encoder -> genotypes[i].left);
            if (k == kh_end(encoder -> labelMap))
                k = kh_put(haplotype, encoder -> labelMap, encoder -> genotypes[i].left, &ret);
            kh_value(encoder -> labelMap, k) = MISSING;
            k = kh_get(haplotype, encoder -> labelMap, encoder -> genotypes[i].right);
            if (k == kh_end(encoder -> labelMap)) 
                k = kh_put(haplotype, encoder -> labelMap, encoder -> genotypes[i].right, &ret);
            kh_value(encoder -> labelMap, k) = MISSING;
        }
    }

    hapCounts[0][0] = hapCounts[0][1] = hapCounts[0][2] = 0;

    for (int i = 0; i < encoder -> numSamples; i++) {
        if (samplesToLabel[i] != -1) {
            if (encoder -> genotypes[i].left != MISSING) { 
                k = kh_get(haplotype, encoder -> labelMap, encoder -> genotypes[i].left);
                if (kh_value(encoder -> labelMap, k) == MISSING) {
                    kh_value(encoder -> labelMap, k) = numUniqueHaps;
                    hapCounts[0][numUniqueHaps + 1] = hapCounts[1][numUniqueHaps + 1] = hapCounts[2][numUniqueHaps + 1] = 0;
                    numUniqueHaps++;
                }
                hapCounts[samplesToLabel[i] - 1][kh_value(encoder -> labelMap, k) + 1]++;
                hapCounts[samplesToLabel[i] - 1][0]++;
            }
            if (encoder -> genotypes[i].right != MISSING) { 
                k = kh_get(haplotype, encoder -> labelMap, encoder -> genotypes[i].right);
                if (kh_value(encoder -> labelMap, k) == MISSING) {
                    kh_value(encoder -> labelMap, k) = numUniqueHaps;
                    hapCounts[0][numUniqueHaps + 1] = hapCounts[1][numUniqueHaps + 1] = hapCounts[2][numUniqueHaps + 1] = 0;
                    numUniqueHaps++;
                }
                hapCounts[samplesToLabel[i] - 1][kh_value(encoder -> labelMap, k) + 1]++;
                hapCounts[samplesToLabel[i] - 1][0]++;
            }
        }
    }

    return numUniqueHaps;
}

Block_t* get_next_block(
    VCFLocusParser_t* vcfFile, 
    HaplotypeEncoder_t* encoder, 
    int* samplesToLabel, 
    int numSamples, 
    int sampleSize, 
    int blockSize, 
    int haplotypeSize,
    int endOfBlock,
    int** hapCounts
) {

    if (isEOF(vcfFile))
        return NULL;

    Block_t* block = init_block(vcfFile -> nextChrom, vcfFile -> nextCoord, sampleSize);

    bool isOnSameChrom = true;
    while (isOnSameChrom && vcfFile -> nextCoord <= endOfBlock) {
        isOnSameChrom = get_next_haplotype(vcfFile, encoder, haplotypeSize);
        int numUniqueHaps = accumulate_counts(encoder, haplotypeSize, samplesToLabel, hapCounts, numSamples);
        locus_privateD(block, hapCounts, numUniqueHaps, sampleSize);
        block -> numHaps++;
    }
    block -> endCoordinate = encoder -> endCoord;

    return block;

}

BlockList_t* privateD(VCFLocusParser_t* vcfFile, HaplotypeEncoder_t* encoder, int* samplesToLabel, int numSamples, int sampleSize, int blockSize, int haplotypeSize) {

    BlockList_t* globalList = init_block_list(sampleSize);

    int** hapCounts = (int**) calloc(3, sizeof(int*));
    hapCounts[0] = (int*) calloc(2 * numSamples + 1, sizeof(int));
    hapCounts[1] = (int*) calloc(2 * numSamples + 1, sizeof(int));
    hapCounts[2] = (int*) calloc(2 * numSamples + 1, sizeof(int));

    while (true) {
        int endOfBlock = ((int) ((vcfFile -> nextCoord - 1) / blockSize) + 1) * blockSize - 1;
        Block_t* temp = get_next_block(vcfFile, encoder, samplesToLabel, numSamples, sampleSize, blockSize, haplotypeSize, endOfBlock, hapCounts);
        if (temp == NULL)
            break;

        append_block(globalList, temp);
        for (int i = 0; i < sampleSize; i++) {
            globalList -> rarefactCounts[i].num += temp -> rarefactCounts[i].num;
            globalList -> rarefactCounts[i].denom += temp -> rarefactCounts[i].denom;
        }
    }

    free(hapCounts[0]); free(hapCounts[1]); free(hapCounts[2]);
    free(hapCounts);

    return globalList;

}