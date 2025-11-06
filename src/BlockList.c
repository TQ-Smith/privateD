
// File: BlockList.c
// Date: 28 March 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Define the list of blocks along the genome.

#include "BlockList.h"
#include <string.h>
#include <stdlib.h>

Block_t* init_block(char* chrom, int startCoordinate) {
    Block_t* block = calloc(1, sizeof(Block_t));
    block -> chrom = strdup(chrom);
    block -> startCoordinate = startCoordinate;
    block -> numLoci = 0;
    block -> p = -1;
    block -> next = NULL;
    return block;
}

BlockList_t* init_block_list(int sampleSize) {
    BlockList_t* blockList = calloc(1, sizeof(BlockList_t));
    blockList -> sampleSize = sampleSize;
    blockList -> numBlocks = 0;
    blockList -> head = NULL;
    blockList -> tail = NULL;
    blockList -> numLoci = 0;
    blockList -> stder = 0;
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
    free(block);
}

void destroy_block_list(BlockList_t* blockList) {
    Block_t* temp = NULL;
    for (int i = 0; i < blockList -> numBlocks; i++) {
        temp = blockList -> head;
        blockList -> head = blockList -> head -> next;
        destroy_block(temp);
    }
    free(blockList);
}

