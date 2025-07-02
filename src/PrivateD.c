
// File: PrivateD.c
// Date: 26 December 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Compute privateD with weighted block jackknife.

#include "PrivateD.h"
#include <math.h>

#define EPS 1e-8

// P-values from jackknife with mean 0 and given std. dev. 
double get_p_val(double d, double std) {
    double Z = fabs(d / std);
    double pnorm = 0.5 * (1 + erf(Z / sqrt(2)));
    return 2 * (1 - pnorm);
}

void weighted_block_jackknife(BlockList_t* blocks) {
    
    // Calculate the jackknife estimator.
    double sum = 0;
    for (Block_t* temp = blocks -> head; temp != NULL; temp = temp -> next) {
        int n = blocks -> denom;
        double dropped = (blocks -> num - temp -> num) / (double) (blocks -> denom - temp -> denom);
        sum += (n - temp -> denom) * dropped / (double) n;
    }
    double est = blocks -> num / (double) blocks -> denom;

    // Our jackknife estimator.
    double jack = blocks -> numBlocks * est - sum;

    // Calculate the standard error.
    sum = 0;
    for (Block_t* temp = blocks -> head; temp != NULL; temp = temp -> next) {
        double h = blocks -> denom / (double) temp -> denom;
        double dropped = (blocks -> num - temp -> num) / (double) (blocks -> denom - temp -> denom);
        double est = blocks -> num / (double) blocks -> denom;
        double pseudo = h * est - (h - 1) * dropped;
        sum += (pseudo - jack) * (pseudo - jack) / (h - 1);
    }
    blocks -> stderr = sqrt(sum / blocks -> numBlocks);

    // Calculate our pvalues genome-wide.
    est = blocks -> num / (double) blocks -> denom;
    blocks -> p = get_p_val(est, blocks -> stderr);

}

// From Szpiech et al 2008.
double Q_gji(int N_j, int N_ji, int g) {
    if (N_j - N_ji < g)
        return 0.0;
    else if (g == 1)
        return (N_j - N_ji) / (double) N_j;
    double Q = 0;
    for (int i = 0; i < g; i++)
        Q += log((N_j - N_ji - i) / (double) (N_j - i));
    return exp(Q);
}

// Calculate privateD at the locus.
void locus_privateD(Block_t* block, int** hapCounts, int numUniqueHaps, int sampleSize) {
    double pi13 = 0, pi23 = 0;
    // Calculate the private allelic richness for the two combinations.
    for (int i = 0; i < numUniqueHaps; i++) {
        pi13 += exp(log(1 - Q_gji(hapCounts[0][0], hapCounts[0][i + 1], sampleSize)) + log(1 - Q_gji(hapCounts[2][0], hapCounts[2][i + 1], sampleSize)) + log(Q_gji(hapCounts[1][0], hapCounts[1][i + 1], sampleSize)));
        pi23 += exp(log(1 - Q_gji(hapCounts[1][0], hapCounts[1][i + 1], sampleSize)) + log(1 - Q_gji(hapCounts[2][0], hapCounts[2][i + 1], sampleSize)) + log(Q_gji(hapCounts[0][0], hapCounts[0][i + 1], sampleSize)));
    }
    // If they are not equal
    if (fabs(pi13 - pi23) > EPS) {
        // Increment our counts.
        if (pi23 > pi13)
            block -> num += 1;
        else 
            block -> num -= 1;
        block -> denom += 1;
    }
}

// Fill hapCounts.
//  Returns the number of unique haplotyes.
int accumulate_counts(HaplotypeEncoder_t* encoder, int haplotypeSize, int* samplesToLabel, int** hapCounts, int numSamples) {

    int numUniqueHaps = 0;

    khiter_t k;
    int ret;

    // We use the hash table in the encoder to accumulate haplotype counts.
    //  Reset all values in the hash table.
    for (int i = 0; i < encoder -> numSamples; i++) {
        if (samplesToLabel[i] != -1) {
            if (encoder -> genotypes[i].left != MISSING) {
                k = kh_get(haplotype, encoder -> labelMap, encoder -> genotypes[i].left);
                if (k == kh_end(encoder -> labelMap))
                    k = kh_put(haplotype, encoder -> labelMap, encoder -> genotypes[i].left, &ret);
                kh_value(encoder -> labelMap, k) = MISSING;
            }
            if (encoder -> genotypes[i].right != MISSING) {
                k = kh_get(haplotype, encoder -> labelMap, encoder -> genotypes[i].right);
                if (k == kh_end(encoder -> labelMap)) 
                    k = kh_put(haplotype, encoder -> labelMap, encoder -> genotypes[i].right, &ret);
                kh_value(encoder -> labelMap, k) = MISSING;
            }
        }
    }

    // Set global counts to 0.
    hapCounts[0][0] = hapCounts[1][0] = hapCounts[2][0] = 0;

    // Accumulate counts.
    for (int i = 0; i < encoder -> numSamples; i++) {
        if (samplesToLabel[i] != -1) {
            // If the left genotype is not missing.
            if (encoder -> genotypes[i].left != MISSING) { 
                k = kh_get(haplotype, encoder -> labelMap, encoder -> genotypes[i].left);
                // If the haplotype is not in the hash table
                if (kh_value(encoder -> labelMap, k) == MISSING) {
                    // Insert it into the hashtable
                    kh_value(encoder -> labelMap, k) = numUniqueHaps;
                    // Set corresponding row of counts to 0.
                    hapCounts[0][numUniqueHaps + 1] = hapCounts[1][numUniqueHaps + 1] = hapCounts[2][numUniqueHaps + 1] = 0;
                    numUniqueHaps++;
                }
                // Increment count of observed haplotype and population count.
                hapCounts[samplesToLabel[i] - 1][kh_value(encoder -> labelMap, k) + 1]++;
                hapCounts[samplesToLabel[i] - 1][0]++;
            }
            // If the right genotype is not missing.
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

    Block_t* block = init_block(vcfFile -> nextChrom, vcfFile -> nextCoord);

    bool isOnSameChrom = true;
    while (isOnSameChrom && vcfFile -> nextCoord <= endOfBlock) {
        isOnSameChrom = get_next_haplotype(vcfFile, encoder, haplotypeSize);
        
        // Fill hapCounts.
        int numUniqueHaps = accumulate_counts(encoder, haplotypeSize, samplesToLabel, hapCounts, numSamples);
        
        // If sample size was not set by the user, then we take the min of the three populations.
        int minNumLineages = (int) fmin(hapCounts[0][0], fmin(hapCounts[1][0], hapCounts[2][0]));
        if (sampleSize == -1)
            sampleSize = minNumLineages;
        // Otherwise, if it was set and we do not have the appropriate number of lineages, we skip the haplotype.
        else if (sampleSize > minNumLineages)
            continue;

        // Use hapCounts for rarefaction calculations.
        locus_privateD(block, hapCounts, numUniqueHaps, sampleSize);
        block -> numHaps++;
    }
    block -> endCoordinate = encoder -> endCoord;
    return block;

}

BlockList_t* privateD(VCFLocusParser_t* vcfFile, HaplotypeEncoder_t* encoder, int* samplesToLabel, int numSamples, int sampleSize, int blockSize, int haplotypeSize) {
    
    BlockList_t* globalList = init_block_list(sampleSize);


    // Count the number of haplotyes in each population.
    //  The first row holds the global counts.
    int** hapCounts = (int**) calloc(3, sizeof(int*));
    hapCounts[0] = (int*) calloc(2 * numSamples + 1, sizeof(int));
    hapCounts[1] = (int*) calloc(2 * numSamples + 1, sizeof(int));
    hapCounts[2] = (int*) calloc(2 * numSamples + 1, sizeof(int));

    while (true) {
        // Get the end position of the block for the next record.
        int endOfBlock = ((int) ((vcfFile -> nextCoord - 1) / (double) blockSize) + 1) * blockSize;
        
        // Get the next block.
        Block_t* temp = get_next_block(vcfFile, encoder, samplesToLabel, numSamples, sampleSize, blockSize, haplotypeSize, endOfBlock, hapCounts);
        if (temp == NULL)
            break;

        // Append to global list of blocks.
        if (globalList -> numBlocks > 0 && strcmp(temp -> chrom, globalList -> tail -> chrom) == 0)
            temp -> blockNumOnChrom = globalList -> tail -> blockNumOnChrom + 1;
        else
            temp -> blockNumOnChrom = 1;
        append_block(globalList, temp);

        // Accumulate genome-wide privateD.
        globalList -> num += temp -> num;
        globalList -> denom += temp -> denom;

        globalList -> numHaps += temp -> numHaps;
    }

    free(hapCounts[0]); free(hapCounts[1]); free(hapCounts[2]);
    free(hapCounts);

    return globalList;

}