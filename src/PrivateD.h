
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

// Partition the genome into blocks. Calculate privateD in each block.
// Accepts:
//  VCFLocusParsder_t* vcfFile -> The VCF file to parse.
//  HaplotypeEncoder_t* encoder -> Used to encode haplotypes from the VCF files.
//  int* samplesToLabels -> Associate the index of a sample with its population.
//  int numSamples -> The total number of samples in the three populations.
//  int samplesSize -> The maximum sample size.
//  int blockSize -> The block size in base-pairs of our genome.
//  int haplotypeSize -> The size of a haplotype in loci.
// Returns: BlockList_t*, the list of genome blocks.
BlockList_t* privateD(VCFLocusParser_t* vcfFile, HaplotypeEncoder_t* encoder, int* samplesToLabel, int numSamples, int sampleSize, int blockSize, int haplotypeSize);

// Compute the weighted jackknife pvalues for global and all blocks.
void weighted_block_jackknife(BlockList_t* blockList);

#endif