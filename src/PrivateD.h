
// File: Block.h
// Date: 26 December 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Compute privateD with weighted jackknife.

#ifndef _BLOCK_H_
#define _BLOCK_H_

#include "VCFLocusParser.h"
#include "HaplotypeEncoder.h"

// Our main method that computes D^g for each g along with the weighted jackknife.
// Accepts:
//  VCFLocusParser_t* vcfFile -> The VCF to read.
//  HaplotypeEncoder_t* encoder -> The encoder used for haplotypes.
//  int g -> Ranges from 2 to g for rarefaction sample size.
//  int blockSize -> The block size used in the jackknife.
//  int h -> The haplotype size.
//  double* D -> Sets D^g for each g.
//  double* lowerCI -> Sets the lower 95% confidence bound for D^g.
//  double* upperCI -> Sets the upper 95% confidence bound for D^g.
//  double* pvals -> Sets the p-value for each g.
// Returns: void.
void privateD(VCFLocusParser_t* vcfFile, HaplotypeEncoder_t* encoder, int g, int blockSize, int h, double* D, double* lowerCI, double* upperCI, double* pvals);

#endif