
// File: Dplus.c
// Date: 8 July 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Compute D and D+ statistics.

#include "Dplus.h"
#include <math.h>
#include <time.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

double pnorm(double x) {
    return 0.5 * (1.0 + erf(x / sqrt(2.0)));
}

void bootstrap(BlockList_t* blockList, int replicates, bool standard) {

    // Allocate the distribution.
    double* dDis = calloc(replicates, sizeof(double));
    double* dplusDis = calloc(replicates, sizeof(double));

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
        double dNum = 0, dplusNum = 0;
        double dDenom = 0, dplusDenom = 0;
        for (int j = 0; j < blockList -> numBlocks; j++) {
            int randomBlock = (int) (blockList -> numBlocks * (double) gsl_rng_uniform(r));
            dNum += blocks[randomBlock] -> dNum;
            dDenom += blocks[randomBlock] -> dDenom;
            dplusNum += blocks[randomBlock] -> dplusNum;
            dplusDenom += blocks[randomBlock] -> dplusDenom;
        }
        dDis[i] = dNum / dDenom;
        dplusDis[i] = dplusNum / dplusDenom;
    }

    // Calculate pvalue for each block and global.
    if (standard) {
        double meanD = 0, meanDplus = 0;
        for (int i = 0; i < replicates; i++) {
            meanD += dDis[i] / replicates;
            meanDplus += dplusDis[i] / replicates;
        }
        double stddevD = 0, stddevDplus = 0;
        for (int i = 0; i < replicates; i++) {
            stddevD += (dDis[i] - meanD) * (dDis[i] - meanD);
            stddevDplus += (dplusDis[i] - meanDplus) * (dplusDis[i] - meanDplus);
        }
        stddevD /= (replicates - 1);
        stddevD = sqrt(stddevD);
        stddevDplus /= (replicates - 1);
        stddevDplus = sqrt(stddevDplus);

        for (Block_t* temp = blockList -> head; temp != NULL; temp = temp -> next) {
            temp -> dP = 1 - pnorm((temp -> dNum / temp -> dDenom) / stddevD);
            temp -> dplusP = 1 - pnorm((temp -> dplusNum / temp -> dplusDenom) / stddevDplus);
        }
        blockList -> dP = 1- pnorm((blockList -> dNum / blockList -> dDenom) / stddevD);
        blockList -> dplusP = 1- pnorm((blockList -> dplusNum / blockList -> dplusDenom) / stddevDplus);
    } else {
        int dNumGreater, dplusNumGreater;
        for (Block_t* temp = blockList -> head; temp != NULL; temp = temp -> next) {
            dNumGreater = 0;
            dplusNumGreater = 0;
            for (int i = 0; i < replicates; i++) {
                if (temp -> dNum / temp -> dDenom > dDis[i])
                    dNumGreater++;
                if (temp -> dplusNum / temp -> dplusDenom > dplusDis[i])
                    dplusNumGreater++;
            }
            temp -> dP = 1 - 2 * fabs(0.5 - dNumGreater / (double) replicates);
            temp -> dplusP = 1 - 2 * fabs(0.5 - dplusNumGreater / (double) replicates);
        }
        dNumGreater = 0;
        dplusNumGreater = 0;
        for (int i = 0; i < replicates; i++) {
            if (blockList -> dNum / blockList -> dDenom > dDis[i])
                dNumGreater++;
            if (blockList -> dplusNum / blockList -> dplusDenom > dplusDis[i])
                dplusNumGreater++;
        }
        blockList -> dP = 1 - 2 * fabs(0.5 - dNumGreater / (double) replicates);
        blockList -> dplusP = 1 - 2 * fabs(0.5 - dplusNumGreater / (double) replicates);
    }

    free(dDis);
    free(dplusDis);
    free(blocks);
    gsl_rng_free(r);
}

Block_t* get_next_block(
    VCFLocusParser_t* vcfFile, 
    int* samplesToLabel,
    int blockSize, 
    int endOfBlock,
    int** alleleCounts
) {

    if (isEOF(vcfFile))
        return NULL;

    Block_t* block = init_block(vcfFile -> nextChrom, vcfFile -> nextCoord);

    char* chrom = NULL;
    int coord;
    Locus* loci = calloc(vcfFile -> numSamples, sizeof(Locus));
    int numAlleles;

    bool isOnSameChrom = true;
    while (isOnSameChrom && vcfFile -> nextCoord <= endOfBlock) {

        get_next_locus(vcfFile, &chrom, &coord, &numAlleles, &loci);
        isOnSameChrom = !isEOF(vcfFile) && strcmp(chrom, vcfFile -> nextChrom) == 0;
        
        if (numAlleles != 2)
            continue;

        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 3; j++)
                alleleCounts[i][j] = 0;

        for (int i = 0; i < vcfFile -> numSamples; i++) {
            if (samplesToLabel[i] != -1) {
                if (LEFT_ALLELE(loci[i]) != numAlleles) {
                    alleleCounts[samplesToLabel[i] - 1][(int) (LEFT_ALLELE(loci[i])) + 1]++;
                    alleleCounts[samplesToLabel[i] - 1][0]++;
                }
                if (RIGHT_ALLELE(loci[i]) != numAlleles) {
                    alleleCounts[samplesToLabel[i] - 1][RIGHT_ALLELE(loci[i]) + 1]++;
                    alleleCounts[samplesToLabel[i] - 1][0]++;
                }
            }
        }

        double ABBA = 0, BABA = 0, BAAA = 0, ABAA = 0, ADDA = 0, BDBD = 0, BBAA = 0;
        double p1 = 0, p2 = 0, p3 = 0, p4 = 0;
        if (alleleCounts[0][0] == 1 && alleleCounts[1][0] == 1 && alleleCounts[2][0] == 1 && alleleCounts[3][0] == 1) {
            if (alleleCounts[0][1] == 1 && alleleCounts[1][2] == 1 && alleleCounts[2][2] == 1 && alleleCounts[3][1] == 1)
                ABBA += 1;
            if (alleleCounts[0][2] == 1 && alleleCounts[1][1] == 1 && alleleCounts[2][2] == 1 && alleleCounts[3][1] == 1)
                BABA += 1;
            if (alleleCounts[0][2] == 1 && alleleCounts[1][1] == 1 && alleleCounts[2][1] == 1 && alleleCounts[3][1] == 1)
                BAAA += 1;
            if (alleleCounts[0][1] == 1 && alleleCounts[1][2] == 1 && alleleCounts[2][1] == 1 && alleleCounts[3][1] == 1)
                ABAA += 1;
            if (alleleCounts[0][2] == 1 && alleleCounts[1][2] == 1 && alleleCounts[2][1] == 1 && alleleCounts[3][1] == 1)
                BBAA += 1;
        } else {
            p1 = alleleCounts[0][2] / (double) alleleCounts[0][0];
            p2 = alleleCounts[1][2] / (double) alleleCounts[1][0];
            p3 = alleleCounts[2][2] / (double) alleleCounts[2][0];
            p4 = alleleCounts[3][2] / (double) alleleCounts[3][0];
            ABBA = (1 - p1) * p2 * p3 * (1 - p4);
            BABA = p1 * (1 - p2) * p3 * (1 - p4);
            BAAA = p1 * (1 - p2) * (1 - p3) * (1 - p4);
            ABAA = (1 - p1) * p2 * (1 - p3) * (1 - p4);

            double pd = fmax(p2, p3);
            ADDA = (1 - p1) * pd * pd * (1 - p4);
            BDBD = p1 * (1 - pd) * pd * (1 - p4);
            BBAA = p1 * p2 * (1 - p3) * (1 - p4);
        }

        block -> dNum += (ABBA - BABA);
        block -> dDenom += (ABBA + BABA);
        block -> dplusNum += (ABBA - BABA + BAAA - ABAA);
        block -> dplusDenom += (ABBA + BABA + BAAA + ABAA);

        block -> fdNum += (ABBA - BABA);
        block -> fdDenom += (ADDA - BDBD);
        block -> dfNum += (ABBA - BABA);
        block -> dfDenom += (ABBA + BBAA + BABA + BBAA);

        // Calculate nucleotide diversity.
        int total = alleleCounts[0][0] + alleleCounts[1][0] + alleleCounts[2][0];
        block -> nucleotideDiversity += 2 * (alleleCounts[0][1] + alleleCounts[1][1] + alleleCounts[2][1]) * (alleleCounts[0][2] + alleleCounts[1][2] + alleleCounts[2][2]) / (total * (total -1));

        block -> numHaps++;
    }
    block -> endCoordinate = coord;

    free(chrom);
    free(loci);

    return block;

}

BlockList_t* dplus(VCFLocusParser_t* vcfFile, int* samplesToLabel, int blockSize) {
    
    BlockList_t* globalList = init_block_list();


    // Count the number of haplotyes in each population.
    //  The first row holds the global counts.
    int** alleleCounts = (int**) calloc(4, sizeof(int*));
    alleleCounts[0] = (int*) calloc(3, sizeof(int));
    alleleCounts[1] = (int*) calloc(3, sizeof(int));
    alleleCounts[2] = (int*) calloc(3, sizeof(int));
    alleleCounts[3] = (int*) calloc(3, sizeof(int));

    while (true) {
        // Get the end position of the block for the next record.
        int endOfBlock = ((int) ((vcfFile -> nextCoord - 1) / (double) blockSize) + 1) * blockSize;
        
        // Get the next block.
        Block_t* temp = get_next_block(vcfFile, samplesToLabel, blockSize, endOfBlock, alleleCounts);
        if (temp == NULL)
            break;

        // Append to global list of blocks.
        if (globalList -> numBlocks > 0 && strcmp(temp -> chrom, globalList -> tail -> chrom) == 0)
            temp -> blockNumOnChrom = globalList -> tail -> blockNumOnChrom + 1;
        else
            temp -> blockNumOnChrom = 1;
        append_block(globalList, temp);

        globalList -> dNum += temp -> dNum;
        globalList -> dDenom += temp -> dDenom;
        globalList -> dplusNum += temp -> dplusNum;
        globalList -> dplusDenom += temp -> dplusDenom;
        globalList -> fdNum += temp -> fdNum;
        globalList -> fdDenom += temp -> fdDenom;
        globalList -> dfNum += temp -> dfNum;
        globalList -> dfDenom += temp -> dfDenom;
        globalList -> nucleotideDiversity += globalList -> nucleotideDiversity;

        globalList -> numHaps += temp -> numHaps;

    }

    free(alleleCounts[0]); free(alleleCounts[1]); free(alleleCounts[2]); free(alleleCounts[3]);
    free(alleleCounts);

    return globalList;

}