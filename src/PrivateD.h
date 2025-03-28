
// File: PrivateD.h
// Date: 26 December 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Compute privateD with weighted jackknife.

#ifndef _PRIVATE_D_H_
#define _PRIVATE_D_H_

#include "VCFLocusParser.h"
#include "HaplotypeEncoder.h"
#include "BlockList.h"

// Accumulate privateD counts for each block along the genome.
// Accepts:
//  VCFLocusParser_t* vcfFile -> The VCF to read.
//  HaplotypeEncoder_t* encoder -> The encoder used for haplotypes.
//  int* samplesToLabel -> Associates the sample with a population.
//  int maxNumOfHaps -> The maximum number of possible haplotypes.
//  int g -> Ranges from 1 to g for rarefaction sample size.
//  int blockSize -> The block size used in the jackknife.
//  int h -> The haplotype size.
// Returns: BlockList_t*, a list of the counts for all the blocks along the genome.
BlockList_t* privateD(VCFLocusParser_t* vcfFile, HaplotypeEncoder_t* encoder, int* samplesToLabel, int maxNumOfHaps, int g, int blockSize, int h);

// Calculates jackknife for all blocks. Results are stored in list.
// Accetps:
//  BlockList_t* blockList -> The list of blocks.
// Returns: void.
void weighted_block_jackknife(BlockList_t* blockList);

#endif