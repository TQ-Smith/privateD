
// File: PrivateD.c
// Date: 26 December 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Compute privateD with weighted block jackknife.

#include <math.h>
#include "PrivateD.h"
#include "Matrix.h"
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

Block_t* get_next_block(VCFLocusParser_t* vcfFile, HaplotypeEncoder_t* encoder, int* sampleIndices, int* samplesToLabel, int** hapCounts, int numSamplesInPops, int g, int blockSize, int h, int* startOfNextBlock) {

    // No block to return if end of file reached.
    if (vcfFile -> isEOF)
        return NULL;

    // Create our block.
    Block_t* block = init_block(ks_str(vcfFile -> nextChrom), vcfFile -> nextCoord, g);
    if (*startOfNextBlock == 1)
        block -> blockNumOnChrom = 1;
    else 
        block -> blockNumOnChrom++;

    // Find start of the block the next record is in.
    if (*startOfNextBlock < vcfFile -> nextCoord)
        *startOfNextBlock = (int) ((vcfFile -> nextCoord - 1) / blockSize) * blockSize + 1;
    
    // End of the current block.
    int endOfBlock = *startOfNextBlock + blockSize - 1;

    // While we read in records within the block.
    bool isOnSameChrom = true;
    while (isOnSameChrom && vcfFile -> nextCoord <= endOfBlock) {
        isOnSameChrom = get_next_haplotype(vcfFile, encoder, h);

        // Create table of haplotype counts.
        int numUniqueHaps = accumulate_counts(encoder, sampleIndices, samplesToLabel, hapCounts, numSamplesInPops);

        // Calculate privateD for locus.
        locus_privateD(block, hapCounts, numUniqueHaps, g);

        block -> numHaps++;
    }

    block -> endCoordinate = encoder -> endCoord;

    // Set the start position of the next block.
    if (!isOnSameChrom)
        *startOfNextBlock = 1;
    else
        *startOfNextBlock = endOfBlock + 1;

    // Return the block.
    return block;

}
*/

Block_t* get_next_block(
    VCFLocusParser_t* vcfFile, 
    HaplotypeEncoder_t* encoder, 
    int* samplesToLabel, 
    int numSamples, 
    int sampleSize, 
    int blockSize, 
    int haplotypeSize,
    int endOfBlock
) {

    if (isEOF(vcfFile))
        return NULL;

    Block_t* block = init_block(vcfFile -> nextChrom, vcfFile -> nextCoord, sampleSize);

    bool isOnSameChrom = true;
    while (isOnSameChrom && vcfFile -> nextCoord <= endOfBlock) {
        isOnSameChrom = get_next_haplotype(vcfFile, encoder, haplotypeSize);
        block -> numHaps++;
    }
    block -> endCoordinate = encoder -> endCoord;

    return block;

}

BlockList_t* privateD(VCFLocusParser_t* vcfFile, HaplotypeEncoder_t* encoder, int* samplesToLabel, int numSamples, int sampleSize, int blockSize, int haplotypeSize) {

    BlockList_t* globalList = init_block_list(sampleSize);

    while (true) {

        int endOfBlock = ((int) ((vcfFile -> nextCoord - 1) / blockSize) + 1) * blockSize - 1;
        Block_t* temp = get_next_block(vcfFile, encoder, samplesToLabel, numSamples, sampleSize, blockSize, haplotypeSize, endOfBlock);
        if (temp == NULL)
            break;

        append_block(globalList, temp);
        for (int i = 0; i < sampleSize; i++) {
            globalList -> rarefactCounts[i].num += temp -> rarefactCounts[i].num;
            globalList -> rarefactCounts[i].denom += temp -> rarefactCounts[i].denom;
        }
    }

    return globalList;

}