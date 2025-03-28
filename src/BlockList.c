
// File: BlockList.c
// Date: 28 March 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Define the list of blocks along the genome.

#include "BlockList.h"
#include <string.h>
#include <stdlib.h>

Block_t* init_block(char* chrom, int startCoordinate, int g) {
    Block_t* block = calloc(1, sizeof(Block_t));
    block -> chrom = strdup(chrom);
    block -> startCoordinate = startCoordinate;
    block -> rarefactCounts = calloc(g, sizeof(Counts_t));
    block -> g = g;
    block -> numHaps = 0;
    block -> next = NULL;
    return block;
}

BlockList_t* init_block_list(int g) {
    BlockList_t* blockList = calloc(1, sizeof(BlockList_t));
    blockList -> g = g;
    blockList -> rarefactCounts = calloc(g, sizeof(Counts_t));
    blockList -> stderrs = calloc(g, sizeof(double));
    blockList -> numBlocks = 0;
    blockList -> head = NULL;
    blockList -> tail = NULL;
    return blockList;
}

void appendBlock(BlockList_t* blockList, Block_t* block) {
    if (blockList -> numBlocks == 0) {
        blockList -> head = block;
        blockList -> tail = block;
        block -> blockNum = 1;
        return;
    }
    blockList -> tail -> next = block;
    blockList -> tail = block;
    // Increment the number of blocks in the list and 
    //  assign the appened block its number in the list.
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
    free(blockList -> stderrs);
    Block_t* temp = NULL;
    for (int i = 0; i < blockList -> numBlocks; i++) {
        temp = blockList -> head;
        blockList -> head = blockList -> head -> next;
        destroy_block(temp);
    }
    free(blockList);
}

