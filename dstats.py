
# File: dstats.py
# Date: 30 June 2025
# Author: T. Quinn Smith
# Principal Investigator: Zachary A. Szpiech 
# Purpose: Calculate D and Dplus in windows.

# NOTE: There is no error checking in this script or significance testing.

import numpy as np
import gzip
import sys
import io
import pandas as pd

# Taken from https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744
def read_vcf(path):
    with gzip.open(path, 'rt') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

# Open up the tsv file associating sample name to population number.
def create_samples_to_pops(inputFile, samplesToPopsFile, pop1, pop2, pop3, pop4):
    samplesToPops = {}
    sampleIndexes = {}
    for i in range(9, len(inputFile.columns)):
        sampleIndexes[inputFile.columns[i]] = i - 9
    for line in open(samplesToPopsFile, 'r'):
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

# We spoof a struct with a dictionary.
def new_block():
    return {'blockNum': 0, 'blockNumOnChrom': 0, 'chr' : '', 'start' : 0, 'end' : 0, 'd_num': 0, 'd_denom': 0, 'dplus_num': 0, 'dplus_denom': 0}

# Parse a genotype entry in a VCF file.
def parse_geno(geno):
    if '|' in geno:
        left, right = geno.split('|')
        if left == '.':
            return (None, int(right))
        elif right == '.':
            return (int(left), None)
        else:
            return (None, None)
    elif '/' in geno:
        left, right = geno.split('/')
        if left == '.':
            return (None, int(right))
        elif right == '.':
            return (int(left), None)
        else:
            return (None, None)
    else:
        return (None, None)

def calc_n1(startIndex, numRecords, inputFile, samplesToPops):
    block = new_block()
    for i in range(startIndex, numRecords + 1):
        # We skip multiallelic sites.
        if inputFile.iloc[i]['ALT'].count(',') > 2:
            continue
        ancestral = None
        states = {}
        # For each sample, we get the non-null genotype.
        for j in range(len(samplesToPops)):
            if samplesToPops[j] == -1:
                continue
            left, right = parse_geno(inputFile.iloc[i, j + 9])
            if left == None:
                states[samplesToPops[j]] = right
            elif right == None:
                states[samplesToPops[j]] = left
            else:
                states[samplesToPops[j]] = None
        
        # If one sample is missing, we skip the record.
        if None in states.values():
            continue

        # Ancestral allele is chosen by the outgroup.
        ancestral = states[4]

        # Count indicator variables.
        ABBA = 0
        BABA = 0
        denom = 0
        if states[1] == ancestral and states[2] != ancestral and states[3] != ancestral:
            ABBA = 1
            denom = 1
        if states[1] != ancestral and states[2] == ancestral and states[3] != ancestral:
            BABA = 1
            denom = 1
        # Accumulate counts.
        block['d_num'] += ABBA - BABA 
        block['d_denom'] += denom

        # Count indicator variables.
        BAAA = 0
        ABAA = 0
        if states[1] != ancestral and states[2] == ancestral and states[3] == ancestral:
            BAAA = 1
            denom = 1
        if states[1] == ancestral and states[2] != ancestral and states[3] == ancestral:
            ABAA = 1
            denom = 1
        # Accumulate counts.
        block['dplus_num'] += ABBA - BABA + BAAA - ABAA
        block['dplus_denom'] += denom 
    
    return block
        
        
# Calculate our statistics based off of allele frequencies.
def calc_freq(startIndex, numRecords, inputFile, samplesToPops):
    block = new_block()
    for i in range(startIndex, numRecords + 1):
        # We only consider biallelic sites.
        if inputFile.iloc[i]['ALT'].count(',') > 2:
            continue
        ancestral = None
        states = {}
        numSamples = {}
        # For each sample, accumulate alternative allele and sample size totals.
        for j in range(len(samplesToPops)):
            if samplesToPops[j] == -1:
                continue
            left, right = parse_geno(inputFile.iloc[i, j + 9])
            if right == 1:
                states[samplesToPops[j]] += 1
                numSamples[samplesToPops[j]] += 1
            if left == 1:
                states[samplesToPops[j]] += 1
                numSamples[samplesToPops[j]] += 1

        # Calculate alternative allele frequencies.
        states[1] /= numSamples[1]
        states[2] /= numSamples[2]
        states[3] /= numSamples[3]
        states[4] /= numSamples[4]

        # Calculate our values.
        ABBA = (1 - states[1]) * states[2] * states[3] * (1 - states[4])
        BABA = states[1] * (1 - states[2]) * (1 - states[3]) * states[4]
        BAAA = states[1] * (1 - states[2]) * (1 - states[3]) * (1 - states[4])
        ABAA = (1 - states[1]) * states[2] * (1 - states[3]) * (1 - states[4])

        # Accumulate our values.
        block['d_num'] += ABBA - BABA 
        block['d_denom'] += ABBA + BABA 
        block['dplus_num'] += (ABBA - BABA) + (BAAA - ABAA)
        block['dplus_denom'] + (ABBA + BABA) + (BAAA + ABAA)

    return block

# Given a coordinate on a chromosome and a block size, return the end coordinate of the current block.
def get_end_coord(coord, blockSize):
    return int((coord - 1) / blockSize + 1) * blockSize

# Our method that calculate the D statistics in blocks along the genome.
def dstats(isN1, inputFile, samplesToPops, blockSize):
    # Our global values.
    globalBlock = new_block()
    # Our list of local blocks.
    blocks = []

    # Used to track the current window
    startIndex = 0
    numRecords = 0
    numBlocks = 1
    numBlocksOnChrom = 1
    while numRecords <= inputFile.shape[0]:
        # If we are on a new chromosome or if we are in a new block
        if numRecords != 0 and (numRecords == inputFile.shape[0] or inputFile.iloc[numRecords]['CHROM'] != inputFile.iloc[numRecords - 1]['CHROM'] or get_end_coord(inputFile.iloc[numRecords]['POS'], blockSize) != get_end_coord(inputFile.iloc[startIndex]['POS'], blockSize)):
            # Calculate the appropriate version of the statistic.
            if isN1:
                block = calc_n1(startIndex, numRecords - 1, inputFile, samplesToPops)
            else:
                block = calc_freq(startIndex, numRecords - 1, inputFile, samplesToPops)
            
            # Set block information.
            block['blockNum'] = numBlocks
            numBlocks += 1
            block['blockNumOnChrom'] = numBlocksOnChrom
            numBlocksOnChrom += 1
            if numRecords != inputFile.shape[0] and inputFile.iloc[numRecords]['CHROM'] != inputFile.iloc[numRecords - 1]['CHROM']:
                numBlocksOnChrom = 1
            block['chr'] = inputFile.iloc[numRecords - 1]['CHROM']
            block['start'] = inputFile.iloc[startIndex]['POS']
            block['end'] = inputFile.iloc[numRecords - 1]['POS']

            # Accumulate global values
            globalBlock['d_num'] += block['d_num']
            globalBlock['d_denom'] += block['d_denom']
            globalBlock['dplus_num'] += block['dplus_num']
            globalBlock['dplus_denom'] += block['dplus_denom']
            blocks.append(block)
            startIndex = numRecords
        numRecords += 1
    # Global counts are in the last element.
    return blocks + [globalBlock]

if __name__ == "__main__":

    if len(sys.argv) == 1:
        print("usage: python3 dstats.py in.vcf.gz samplesToPops.tsv blockSize pop1,pop2,pop3,pop4")
        exit()

    # Parse our command line options.
    inputFile = read_vcf(sys.argv[1])
    samplesToPopsFile = sys.argv[2]
    blockSize = int(sys.argv[3])
    pop1, pop2, pop3, pop4 = sys.argv[4].split(',')

    # Associate sample indices with population label.
    samplesToPops = create_samples_to_pops(inputFile, samplesToPopsFile, pop1, pop2, pop3, pop4)

    # Decide if we are calculating n1 stat or not.
    #   This is lazy. If there are 4 samples, then we calculate n1. This is not always the case. We 
    #   are assuming that the four samples only each have one lineage. 
    isN1 = True
    n = 0
    for sample in samplesToPops:
        if samplesToPops[sample] != -1:
            n += 1
    if n > 4:
        isN1 = False
    
    # Compute statistics
    blocks = dstats(isN1, inputFile, samplesToPops, blockSize)

    # Output values.
    print("#python3 dstats.py", sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    print('BLOCK\tBLOCK_ON_CHROM\tCHROM\tSTART\tEND\tNUM_LOCI\tD\tD_PLUS')
    for i in range(len(blocks) - 1):
        print(blocks[i]['blockNum'], blocks[i]['blockNumOnChrom'], blocks[i]['chr'], blocks[i]['start'], blocks[i]['end'], blocks[i]['dplus_denom'], blocks[i]['d_num'] / blocks[i]['d_denom'], blocks[i]['dplus_num'] / blocks[i]['dplus_denom'], sep='\t')
    i = len(blocks) - 1
    print('0', '0', 'Global', blocks[i]['start'], blocks[i]['end'], blocks[i]['dplus_denom'], blocks[i]['d_num'] / blocks[i]['d_denom'], blocks[i]['dplus_num'] / blocks[i]['dplus_denom'], sep='\t')
