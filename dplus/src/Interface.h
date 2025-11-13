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
#include <stdbool.h>

// All the possible values that define a run of privateD.
typedef struct {
    int blockSize;
    double MAF;
    double missingAF;
    int replicates;
    bool standard;
    char* inputFileName;
    char* samplesToPopFileName;
    char* fourPopList;
    char* cmd;
    char* outBaseName;
} PrivateDConfig_t;

// Parse commnad line arguments.
PrivateDConfig_t* init_privated_config(int argc, char* argv[]);

void destroy_privated_config(PrivateDConfig_t* config);

void print_help();

#endif