
// File: Dplus.c
// Date: 8 July 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Compute D and D+ statistics.

#include "Dplus.h"
#include <math.h>

Block_t* get_next_block(
    VCFLocusParser_t* vcfFile, 
    int* samplesToLabel, 
    int numSamples, 
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

        get_next_locus(vcfFile, &chrom, &coord, &numAlleles, &loci);
        isOnSameChrom = !isEOF(vcfFile) && strcmp(chrom, vcfFile -> nextChrom) == 0;
        
        if (numAlleles != 2)
            continue;

        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 3; j++)
                alleleCounts[i][j] = 0;

        for (int i = 0; i < numSamples; i++) {
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
                ABBA += 1;
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

        block -> numHaps++;
    }
    block -> endCoordinate = coord;

    free(chrom);
    free(loci);

    return block;

}

BlockList_t* dplus(VCFLocusParser_t* vcfFile, int* samplesToLabel, int numSamples, int blockSize) {
    
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
        Block_t* temp = get_next_block(vcfFile, samplesToLabel, numSamples, blockSize, endOfBlock, alleleCounts);
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

        globalList -> numHaps += temp -> numHaps;

    }

    free(alleleCounts[0]); free(alleleCounts[1]); free(alleleCounts[2]); free(alleleCounts[3]);
    free(alleleCounts);

    return globalList;

}