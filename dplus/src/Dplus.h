
// File: Dplus.h
// Date: 26 December 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Compute privateD with weighted jackknife.

#ifndef _D_PLUS_H_
#define _D_PLUS_H_

#include "VCFLocusParser.h"
#include "BlockList.h"

// Partition the genome into blocks. Calculate privateD in each block.
// Accepts:
//  VCFLocusParsder_t* vcfFile -> The VCF file to parse.
//  HaplotypeEncoder_t* encoder -> Used to encode haplotypes from the VCF files.
//  int* samplesToLabels -> Associate the index of a sample with its population.
//  int numSamples -> The total number of samples in the four populations.
//  int blockSize -> The block size in base-pairs of our genome.
// Returns: BlockList_t*, the list of genome blocks.
BlockList_t* dplus(VCFLocusParser_t* vcfFile, int* samplesToLabel, int numSamples, int blockSize);

void bootstrap(BlockList_t* blockList, int replicates);

#endif