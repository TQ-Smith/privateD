
#ifndef _INTERFACE_H_
#define _INTERFACE_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef struct {
    int sampleSize;
    int haplotypeSize;
    int blockSize;
    double MAF;
    double missingAF;
    char* inputFileName;
    char* samplesToPopFileName;
    char* threePopList;
    char* cmd;
    char* outBaseName;
} PrivateDConfig_t;

PrivateDConfig_t* init_privated_config(int argc, char* argv[]);

void destroy_privated_config(PrivateDConfig_t* config);

void print_help();

#endif