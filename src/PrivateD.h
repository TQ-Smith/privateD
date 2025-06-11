
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

BlockList_t* privateD(VCFLocusParser_t* vcfFile, HaplotypeEncoder_t* encoder, int* samplesToLabel, int numSamples, int sampleSize, int blockSize, int haplotypeSize);

void weighted_block_jackknife(BlockList_t* blockList);

#endif