// File: Interface.h
// Date: 19 June 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Parser command line options.

#ifndef _INTERFACE_H_
#define _INTERFACE_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// All the possible values that define a run of DSTAR.
typedef struct {
    int sampleSize;
    int blockSize;
    double MAF;
    double missingAF;
    int replicates;
    char* inputFileName;
    char* samplesToPopFileName;
    char* threePopList;
    char* cmd;
    char* outBaseName;
} DSTARConfig_t;

// Parse commnad line arguments.
DSTARConfig_t* init_dstar_config(int argc, char* argv[]);

void destroy_dstar_config(DSTARConfig_t* config);

void print_help();

#endif