
// File: BlockList.h
// Date: 28 March 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Define the list of blocks along the genome.

#ifndef _BLOCK_LIST_H_
#define _BLOCK_LIST_H_

// Track the counts for a given sample size in a block.
//  Also, store p value for convience.
typedef struct Counts {
    int num;
    int denom;
    double p;
} Counts_t;

// A node in the BlockList.
typedef struct Block {
    int blockNum;
    int blockNumOnChrom;
    char* chrom;
    int startCoordinate;
    int endCoordinate;
    int numHaps;
    int g;
    // An array of size g to hold counts.
    Counts_t* rarefactCounts;
    struct Block* next;
} Block_t;

// A list of blocks.
typedef struct BlockList {
    int numSamples;
    int g;
    // Array of size g. Holds global counts.
    Counts_t* rarefactCounts;
    // Array to hold standard errors from jackknife.
    double* stderrs;
    int numBlocks;
    // Head and tail pointer for our linked list.
    Block_t* head;
    Block_t* tail;
} BlockList_t;

// Creates a list of blocks.
// Accepts:
//  int g -> The standardized samples size. Redundant but good to track.
// Returns: An empty block list.
BlockList_t* init_block_list(int g);

// Creates a block.
// Acccepts:
//  char* chrom -> The chromosome the block is on.
//  int startCoordinate -> The coordinate of the firs record in the block.
//  int g -> The standardized samples size.
Block_t* init_block(char* chrom, int startCoordinate, int g);

// Adds a block to the block list.
// Accepts:
//  BlockList_t* blockList -> The list to add the block to.
//  Block_t* block -> Adds the block to the end of the list.
// Returns: void.
void appendBlock(BlockList_t* blockList, Block_t* block);

// Frees the memory occupied by the block list.
// Accepts:
//  BlockList_t* blockList -> The block list to destroy.
// Returns: void.
void destroy_block_list(BlockList_t* blockList);

#endif