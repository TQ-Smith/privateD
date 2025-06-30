
# File: dstats.py
# Date: 30 June 2025
# Author: T. Quinn Smith
# Principal Investigator: Zachary A. Szpiech 
# Purpose: Calculate D and Dplus in windows.

# NOTE: There is no error checking in this script or significance testing.

import numpy as np
import gzip
import sys

# Use a VCF reader to make my life a little better.
# Import PyVCF
import vcf

# Open up the tsv file associating sample name to population number.
def create_samples_to_pops(inputFile, samplesToPopsFile, pop1, pop2, pop3, pop4):
    samplesToPops = {}
    sampleIndexes = {}
    for sample, i in zip(inputFile.samples, range(len(inputFile.samples))):
        sampleIndexes[sample] = i
    for line in open(samplesToPops, 'r'):
        name, pop = line.split()
        if pop == pop1:
            samplesToPops[sampleIndexes[name]] = 1
        elif pop == pop2:
            samplesToPops[sampleIndexes[name]] = 2
        elif pop == pop3:
            samplesToPops[sampleIndexes[name]] = 3
        elif pop == pop4:
            samplesToPops[sampleIndexes[name]] = 4
        else:
            samplesToPops[sampleIndexes[name]] = -1
    return samplesToPops

# Our structure to hold the contents of a block.
class Block_t:
    self.chrom = None
    self.start = -1
    self.end = -1
    self.d_num = 0
    self.d_denom = 0
    self.dplus_num = 0
    self.dplus_denom = 0

def calc_block_n1(records, samplesToPops):
    block = Block_t()
    block.chrom = records[0][0]
    block.start = records[0][1]
    block.end = records[len(records) - 1][1]

    
def calc_block_freq(records, samplesToPops):
    block = Block_t()
    block.chrom = records[0][0]
    block.start = records[0][1]
    block.end = records[len(records) - 1][1]

def dstat(inputFile, samplesToPops, blockSize):
    # Decide if we are calculating the statistics off of counts or frequencies.
    n1 = False
    if samplesToPops.values().count(1) == 1 and samplesToPops.values().count(2) == 1 and samplesToPops.values().count(3) == 1 and samplesToPops.values().count(4) == 1:
        n1 = True

    # Our genome-wide value.
    globalBlock = Block_t()

    # Hold our local blocks.
    blocks = []

    # Used to track local blocks.
    prevChrom = ''
    prevEndPosition = 0

    # We append each record within a block to a list. This is really lazy.
    records = []

    # For each record in the VCF
    for record in inputFile:
        # Skip non-biallelic sites.
        if record.ALT.count(",") >= 1:
            continue
        endPos = (int((record.POS - 1) / blockSize) + 1) * blockSize
        # If record is within the current block, append to list.
        if prevChrom != '' and prevEndPosition != 0 and record.CHROM == prevChrom and record.POS < prevEndPosition:
            records.append((record.CHROM, record.POS, record.samples))
        else:
            if n1:
                block = calc_block_n1(records, samplesToPops)
            else:
                block = calc_block_freq(records, samplesToPops)
            # Add values to global.
            globalBlock.d_num += block.d_num
            globalBlock.d_denom += block.d_denom
            globalBlock.dplus_num += block.dplus_num
            globalBlock.dplus_denom += block.dplus_denom
            # Reset for next block.
            prevChrom = record.CHROM 
            prevEndPosition = endPos
            records = []
    return [globalBlock] + blocks

if __name__ == "__main__":

    if len(sys.argv) == 1:
        print("usage: python3 dstats in.vcf.gz samplesToPops.tsv blockSize pop1,pop2,pop3,pop4")
        exit()

    # Parse our command line options.
    inputFile = vcf.Reader(gzip.open(sys.argv[1], 'rt'))
    samplesToPopsFile = sys.argv[2]
    blockSize = int(sys.argv[3])
    pop1, pop2, pop3, pop4 = sys.argv[4].split(',')

    # Associate each sample with a population number.
    samplesToPops = create_samples_to_pops(inputFile, samplesToPopsFile, pop1, pop2, pop3, pop4)

    print(samplesToPops)
    del samplesToPops
    # First block is the global block. The rest are the local blocks.
    #   This is a list of tuples.
    # blocks = dstat(inputFile, samplesToPops, pop1, pop2, pop3, pop4)
    