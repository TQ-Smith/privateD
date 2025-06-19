
// File: BlockList.c
// Date: 28 March 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Define the list of blocks along the genome.

#include "BlockList.h"
#include <string.h>
#include <stdlib.h>

Block_t* init_block(char* chrom, int startCoordinate, int sampleSize) {
    Block_t* block = calloc(1, sizeof(Block_t));
    block -> chrom = strdup(chrom);
    block -> startCoordinate = startCoordinate;
    block -> rarefactCounts = calloc(sampleSize, sizeof(Counts_t));
    block -> numHaps = 0;
    block -> next = NULL;
    return block;
}

BlockList_t* init_block_list(int sampleSize) {
    BlockList_t* blockList = calloc(1, sizeof(BlockList_t));
    blockList -> sampleSize = sampleSize;
    blockList -> rarefactCounts = calloc(sampleSize, sizeof(Counts_t));
    blockList -> numBlocks = 0;
    blockList -> head = NULL;
    blockList -> tail = NULL;
    blockList -> stderrs = (double*) calloc(sampleSize, sizeof(double));
    return blockList;
}

void append_block(BlockList_t* blockList, Block_t* block) {
    if (blockList -> numBlocks == 0) {
        blockList -> head = block;
        blockList -> tail = block;
        blockList -> numBlocks = 1;
        block -> blockNum = 1;
        return;
    }
    blockList -> tail -> next = block;
    blockList -> tail = block;
    blockList -> numBlocks++;
    block -> blockNum = blockList -> numBlocks;
}

void destroy_block(Block_t* block) {
    free(block -> chrom);
    free(block -> rarefactCounts);
    free(block);
}

void destroy_block_list(BlockList_t* blockList) {
    free(blockList -> rarefactCounts);
    Block_t* temp = NULL;
    for (int i = 0; i < blockList -> numBlocks; i++) {
        temp = blockList -> head;
        blockList -> head = blockList -> head -> next;
        destroy_block(temp);
    }
    free(blockList -> stderrs);
    free(blockList);
}

