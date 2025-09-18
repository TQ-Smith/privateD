
// File: BlockList.h
// Date: 28 March 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Define the list of blocks along the genome.

#ifndef _BLOCK_LIST_H_
#define _BLOCK_LIST_H_

// A node in the BlockList.
typedef struct Block {
    int blockNum;
    int blockNumOnChrom;
    char* chrom;
    int startCoordinate;
    int endCoordinate;
    int numHaps;

    // Block counts.
    double dNum;
    double dDenom;
    double dplusNum;
    double dplusDenom;
    double fdNum;
    double fdDenom;
    double dfNum;
    double dfDenom;

    struct Block* next;
} Block_t;

// A list of blocks.
typedef struct BlockList {
    int numSamples;
    int sampleSize;
    
    // Global counts.
    double dNum;
    double dDenom;
    double dplusNum;
    double dplusDenom;
    double fdNum;
    double fdDenom;
    double dfNum;
    double dfDenom;

    int numHaps;

    int numBlocks;
    Block_t* head;
    Block_t* tail;

    double stderr;
} BlockList_t;

// Creates a list of blocks.
// Accepts: void.
// Returns: An empty block list.
BlockList_t* init_block_list();

// Creates a block.
// Acccepts:
//  char* chrom -> The chromosome the block is on.
//  int startCoordinate -> The coordinate of the firs record in the block.
Block_t* init_block(char* chrom, int startCoordinate);

// Adds a block to the block list.
// Accepts:
//  BlockList_t* blockList -> The list to add the block to.
//  Block_t* block -> Adds the block to the end of the list.
// Returns: void.
void append_block(BlockList_t* blockList, Block_t* block);

// Frees the memory occupied by the block list.
// Accepts:
//  BlockList_t* blockList -> The block list to destroy.
// Returns: void.
void destroy_block_list(BlockList_t* blockList);

#endif