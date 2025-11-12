
// File: DSTAR.c
// Date: 26 December 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Compute privateD with weighted block jackknife.

#include "DSTAR.h"
#include <math.h>
#include <time.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

#define EPS 1e-8

void bootstrap(BlockList_t* blockList, int replicates) {

    // Allocate the distribution.
    double* dis = calloc(replicates, sizeof(double));

    // Allocate our RNG.
    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));

    // Array for easy block look up.
    Block_t** blocks = calloc(blockList -> numBlocks, sizeof(Block_t*));
    int k = 0;
    for (Block_t* temp = blockList -> head; temp != NULL; temp = temp -> next)
        blocks[k++] = temp;

    // Create our distribution.
    for (int i = 0; i < replicates; i++) {
        double num = 0;
        double denom = 0;
        for (int j = 0; j < blockList -> numBlocks; j++) {
            int randomBlock = (int) (blockList -> numBlocks * (double) gsl_rng_uniform(r));
            num += blocks[randomBlock] -> numeratorDSTAR;
            denom += blocks[randomBlock] -> denominatorDSTAR;
        }
        dis[i] = num / denom;
    }

    // Calculate pvalue for each block and global.
    int numGreater;
    for (Block_t* temp = blockList -> head; temp != NULL; temp = temp -> next) {
        numGreater = 0;
        for (int i = 0; i < replicates; i++) {
            if (temp -> numeratorDSTAR / temp -> denominatorDSTAR > dis[i])
                numGreater++;
        }
        temp -> p = 1 - 2 * fabs(0.5 - numGreater / (double) replicates);
    }
    numGreater = 0;
    for (int i = 0; i < replicates; i++) {
        if (blockList -> numeratorDSTAR / blockList -> denominatorDSTAR > dis[i])
            numGreater++;
    }
    blockList -> p = 1 - 2 * fabs(0.5 - numGreater / (double) replicates);

    free(dis);
    free(blocks);
    gsl_rng_free(r);
}

/*
// P-values from jackknife with mean 0 and given std. dev. 
double get_p_val(double d, double std) {
    double Z = fabs(d / std);
    double pnorm = 0.5 * (1 + erf(Z / sqrt(2)));
    return 2 * (1 - pnorm);
}

void weighted_block_jackknife(BlockList_t* blocks) {
    
    // Calculate the jackknife estimator.
    double sum = 0;
    int n = blocks -> numHaps;
    for (Block_t* temp = blocks -> head; temp != NULL; temp = temp -> next) {
        double dropped = (blocks -> numeratorPrivateD - temp -> numeratorPrivateD) /  (blocks -> denominatorPrivateD - temp -> denominatorPrivateD);
        sum += (n - temp -> numHaps) * dropped / (double) n;
    }
    double est = blocks -> numeratorPrivateD / (double) blocks -> denominatorPrivateD;

    // Our jackknife estimator.
    double jack = blocks -> numBlocks * est - sum;

    // Calculate the standard error.
    sum = 0;
    for (Block_t* temp = blocks -> head; temp != NULL; temp = temp -> next) {
        double h = n / (double) temp -> numHaps;
        double dropped = (blocks -> numeratorPrivateD - temp -> numeratorPrivateD) / (double) (blocks -> denominatorPrivateD - temp -> denominatorPrivateD);
        double pseudo = h * est - (h - 1) * dropped;
        sum += (pseudo - jack) * (pseudo - jack) / (h - 1);
    }
    blocks -> stder = sqrt(sum / (double) blocks -> numBlocks);

    // Calculate our pvalue for every block and global.
    blocks -> p = get_p_val(est, blocks -> stder);
    for (Block_t* temp = blocks -> head; temp != NULL; temp = temp -> next) {
        est = temp -> numeratorPrivateD / (double) temp -> denominatorPrivateD;
        temp -> p = get_p_val(est, blocks -> stder);
    }
}
*/

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

// Calculate DSTAR at the locus.
void locus_dstar(Block_t* block, int** alleleCounts, int numAlleles, int sampleSize) {
    double pi13 = 0, pi23 = 0, alpha = 0;
    // Total number of lineages.
    int totalLineages = alleleCounts[0][0] + alleleCounts[1][0] + alleleCounts[2][0];
    // Calculate the private allelic richness for the two combinations.
    for (int i = 0; i < numAlleles; i++) {
        pi13 += exp(log(1 - Q_gji(alleleCounts[0][0], alleleCounts[0][i + 1], sampleSize)) + log(1 - Q_gji(alleleCounts[2][0], alleleCounts[2][i + 1], sampleSize)) + log(Q_gji(alleleCounts[1][0], alleleCounts[1][i + 1], sampleSize)));
        pi23 += exp(log(1 - Q_gji(alleleCounts[1][0], alleleCounts[1][i + 1], sampleSize)) + log(1 - Q_gji(alleleCounts[2][0], alleleCounts[2][i + 1], sampleSize)) + log(Q_gji(alleleCounts[0][0], alleleCounts[0][i + 1], sampleSize)));
        alpha += (1 - Q_gji(totalLineages, alleleCounts[0][i + 1] + alleleCounts[1][i + 1] + alleleCounts[2][i + 1], sampleSize));
    }

    if (fabs(pi13) > EPS || fabs(pi23) > EPS) {
        block -> pi13 += pi13;
        block -> pi23 += pi23;
        block -> numeratorDSTAR += (pi23 - pi13);
        block -> denominatorDSTAR += alpha;
        block -> numLoci++;
    }
}

Block_t* get_next_block(
    VCFLocusParser_t* vcfFile,
    int* samplesToLabel, 
    int numSamples, 
    int sampleSize, 
    int blockSize,
    int endOfBlock,
    int** alleleCounts
) {

    if (isEOF(vcfFile))
        return NULL;

    Block_t* block = init_block(vcfFile -> nextChrom, vcfFile -> nextCoord);

    char* chrom = NULL;
    int coord;
    Locus* loci = calloc(numSamples, sizeof(Locus));
    int numAlleles;

    bool isOnSameChrom = true;
    while (isOnSameChrom && vcfFile -> nextCoord <= endOfBlock) {
        
        // Get the locus.
        get_next_locus(vcfFile, &chrom, &coord, &numAlleles, &loci);
        isOnSameChrom = !isEOF(vcfFile) && strcmp(chrom, vcfFile -> nextChrom) == 0;
        
        // Fill alleleCounts.
        for (int i = 0; i < 3; i++)
            for (int j = 0; j <= numAlleles; j++)
                alleleCounts[i][j] = 0;

        for (int i = 0; i < numSamples; i++) {
            if (samplesToLabel[i] != -1) {
                if (((int) LEFT_ALLELE(loci[i])) != numAlleles) {
                    alleleCounts[samplesToLabel[i] - 1][((int) (LEFT_ALLELE(loci[i]))) + 1]++;
                    alleleCounts[samplesToLabel[i] - 1][0]++;
                }
                if (((int) RIGHT_ALLELE(loci[i])) != numAlleles) {
                    alleleCounts[samplesToLabel[i] - 1][((int) (RIGHT_ALLELE(loci[i]))) + 1]++;
                    alleleCounts[samplesToLabel[i] - 1][0]++;
                }
            }
        }

        // If sample size was not set by the user, then we take the min of the three populations.
        int minNumLineages = (int) fmin(alleleCounts[0][0], fmin(alleleCounts[1][0], alleleCounts[2][0]));
        if (sampleSize == -1)
            sampleSize = minNumLineages;
        // Otherwise, if it was set and we do not have the appropriate number of lineages, we skip the block.
        else if (sampleSize > minNumLineages)
            continue;

        // Use alleleCounts for rarefaction calculations.
        locus_dstar(block, alleleCounts, numAlleles, sampleSize);
    }
    block -> endCoordinate = coord;

    free(chrom);
    free(loci);

    return block;

}

BlockList_t* dstar(VCFLocusParser_t* vcfFile, int* samplesToLabel, int numSamples, int sampleSize, int blockSize) {
    
    BlockList_t* globalList = init_block_list(sampleSize);


    // Count the number of haplotyes in each population.
    //  The first row holds the global counts.
    int** alleleCounts = (int**) calloc(3, sizeof(int*));
    alleleCounts[0] = (int*) calloc(2 * numSamples + 1, sizeof(int));
    alleleCounts[1] = (int*) calloc(2 * numSamples + 1, sizeof(int));
    alleleCounts[2] = (int*) calloc(2 * numSamples + 1, sizeof(int));

    while (true) {
        // Get the end position of the block for the next record.
        int endOfBlock = ((int) ((vcfFile -> nextCoord - 1) / (double) blockSize) + 1) * blockSize;
        
        // Get the next block.
        Block_t* temp = get_next_block(vcfFile, samplesToLabel, numSamples, sampleSize, blockSize, endOfBlock, alleleCounts);
        if (temp == NULL)
            break;

        // Append to global list of blocks.
        if (globalList -> numBlocks > 0 && strcmp(temp -> chrom, globalList -> tail -> chrom) == 0)
            temp -> blockNumOnChrom = globalList -> tail -> blockNumOnChrom + 1;
        else
            temp -> blockNumOnChrom = 1;
        append_block(globalList, temp);

        // Accumulate genome-wide DSTAR.
        globalList -> pi13 += temp -> pi13;
        globalList -> pi23 += temp -> pi23;
        globalList -> numeratorDSTAR += temp -> numeratorDSTAR;
        globalList -> denominatorDSTAR += temp -> denominatorDSTAR;

        globalList -> numLoci += temp -> numLoci;
    }

    free(alleleCounts[0]); free(alleleCounts[1]); free(alleleCounts[2]);
    free(alleleCounts);

    return globalList;

}