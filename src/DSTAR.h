
// File: DSTAR.h
// Date: 26 December 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Compute DSTAR with weighted jackknife.

#ifndef _DSTAR_H_
#define _DSTAR_H_

#include "VCFLocusParser.h"
#include "BlockList.h"

// Partition the genome into blocks. Calculate DSTAR in each block.
// Accepts:
//  VCFLocusParsder_t* vcfFile -> The VCF file to parse.
//  int* samplesToLabels -> Associate the index of a sample with its population.
//  int numSamples -> The total number of samples in the three populations.
//  int samplesSize -> The maximum sample size.
//  int blockSize -> The block size in base-pairs of our genome.
// Returns: BlockList_t*, the list of genome blocks.
BlockList_t* dstar(VCFLocusParser_t* vcfFile, int* samplesToLabel, int numSamples, int sampleSize, int blockSize);

// Compute the bootstrapped pvalues for global and all blocks.
void bootstrap(BlockList_t* blockList, int replicates);

#endif