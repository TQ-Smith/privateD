
// File: PrivateD.c
// Date: 26 December 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Compute privateD with weighted block jackknife.

#include <math.h>
#include "PrivateD.h"
#include "Matrix.h"

// Our epsilon for floating point comparison.
#define EPS 1e-10

MATRIX_INIT(int, int)

// Calculate the pvalue from the estimator and the standard error
//  given by a jackknife.
static inline double normCDF(double x) {
    return 0.05 * (1 + erf(x / sqrt(2)));
}
static inline double pval(double val, double stderr) {
    double Z = fabs(val) / stderr;
    return 2 * (1 - normCDF(Z));
}

// Calculated weighted block jackknife standard errors and pvalues for D^g.
//  References are Nick Patterson's notes.
// Accepts:
//  BlockList_t* globalList -> Our list of blocks along the genome.
// Returns: void.
void weighted_block_jackknife(BlockList_t* globalList) {
    int g = globalList -> g;

    // Compute the jackknifed estimator for each g.
    double leaveOutBlock, est;
    double* jackEst = calloc(g, sizeof(double));
    Block_t* curBlock = globalList -> head;
    for (int i = 0; i < globalList -> numBlocks; i++) {
        for (int j = 0; i < g; i++) {
            leaveOutBlock = (globalList -> rarefactCounts[j].num - curBlock -> rarefactCounts[j].num) / (double) (globalList -> rarefactCounts[j].denom - curBlock -> rarefactCounts[j].denom);
            est = globalList -> rarefactCounts[j].num / (double) globalList -> rarefactCounts[j].denom;
            jackEst[i] += (est - leaveOutBlock) + ((curBlock -> rarefactCounts[j].denom * leaveOutBlock) / (double) globalList -> rarefactCounts[j].denom);
        }
        curBlock = curBlock -> next;
    }

    // Compute the variance of the pseudovalues.
    double h, pseudo;
    curBlock = globalList -> head;
    for (int i = 0; i < globalList -> numBlocks; i++) {
        for (int j = 0; i < g; i++) {
            leaveOutBlock = (globalList -> rarefactCounts[j].num - curBlock -> rarefactCounts[j].num) / (double) (globalList -> rarefactCounts[j].denom - curBlock -> rarefactCounts[j].denom);
            h = globalList -> rarefactCounts[j].denom / (double) curBlock -> rarefactCounts[j].denom;
            est = globalList -> rarefactCounts[j].num / (double) globalList -> rarefactCounts[j].denom;
            pseudo = h * est - (h - 1) * leaveOutBlock;
            globalList -> stderrs[j] += (pseudo - est) * (pseudo - est) / (double) (h - 1);
        }
        curBlock = curBlock -> next;
    }

    // Calculate standard errors.
    for (int i = 0; i < g; i++) {
        globalList -> stderrs[i] = sqrt(globalList -> stderrs[i] / globalList -> numBlocks);
    }

    // Finally, we calculate the p-values.
    for (int i = 0; i < g; i++) {
        curBlock = globalList -> head;
        for (int j = 0; j < globalList -> numBlocks; j++) {
            curBlock -> rarefactCounts[i].p = pval(curBlock -> rarefactCounts[i].num / (double) curBlock -> rarefactCounts[i].denom, globalList -> stderrs[i]);
            curBlock = curBlock -> next;
        }
        globalList -> rarefactCounts[i].p = pval(globalList -> rarefactCounts[i].num / (double) globalList -> rarefactCounts[i].denom, globalList -> stderrs[i]);
    }

    free(jackEst);

}

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
            temp = log(1 - Q_ijg(hapCounts[0][0], hapCounts[i][0], g)) + log(1 - Q_ijg(hapCounts[0][2], hapCounts[i][2], g)) + log(Q_ijg(hapCounts[0][1], hapCounts[i][1], g));
            pi13 += exp(temp);
            temp = log(1 - Q_ijg(hapCounts[0][1], hapCounts[i][1], g)) + log(1 - Q_ijg(hapCounts[0][2], hapCounts[i][2], g)) + log(Q_ijg(hapCounts[0][0], hapCounts[i][0], g));
            pi23 += exp(temp);
        }
        // If pi13 and pi23 are NOT equal, the locus contributes to D^g.
        if (fabs(pi13 - pi23) > EPS) {
            if (pi23 > pi13)
                block -> rarefactCounts[g - 1].num += 1;
            else 
                block -> rarefactCounts[g - 1].num -= 1;
            block -> rarefactCounts[g - 1].denom += 1;
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

BlockList_t* privateD(VCFLocusParser_t* vcfFile, HaplotypeEncoder_t* encoder, int* samplesToLabel, int maxNumOfHaps, int g, int blockSize, int h) {
    
    // Indices of our samples in the three populations in the VCF file.
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

    // Our list of blocks.
    BlockList_t* globalList = init_block_list(g);

    // We use this to track start coordinate of the next block.
    int startOfNextBlock = 1;
    // Accumulate our blocks.
    Block_t* temp = NULL;
    while ((temp = get_next_block(vcfFile, encoder, sampleIndices, samplesToLabel, hapCounts, maxNumOfHaps / 2, g, blockSize, h, &startOfNextBlock)) != NULL) {
        appendBlock(globalList, temp);
        // Add counts from block to global.
        for (int i = 0; i < g; i++) {
            globalList -> rarefactCounts[i].num += temp -> rarefactCounts[i].num;
            globalList -> rarefactCounts[i].denom += temp -> rarefactCounts[i].denom;
        }
    }

    // Free all used memory.
    free(sampleIndices);
    destroy_int_matrix(hapCounts, maxNumOfHaps + 1);

    return globalList;
}