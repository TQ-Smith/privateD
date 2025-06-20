// File: Interface.c
// Date: 19 June 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Parser command line options.

#include "Interface.h"
#include "kstring.h"
#include "ketopt.h"
#include <unistd.h>

void print_help() {
    fprintf(stderr, "\n");
    fprintf(stderr, "privateD v1.0\n");
    fprintf(stderr, "----------------------\n\n");
    fprintf(stderr, "Written by T. Quinn Smith\n");
    fprintf(stderr, "Principal Investigator: Zachary A. Szpiech\n");
    fprintf(stderr, "The Pennsylvania State University\n\n");
    fprintf(stderr, "Usage: privateD [options] <inFile.vcf.gz> <sampleToPop.tsv> <pop1>,<pop2>,<pop3>\n\n");
    fprintf(stderr, "<inFile.vcf.gz>                    The input VCF file.\n");
    fprintf(stderr, "<sampleToPop.tsv>                  Tab seperate file associating each sample with a population.\n");
    fprintf(stderr, "<pop1>,<pop2>,<pop3>               Names of the three populations to test.\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -g,--sampleSize        INT         The maximum standardized sample size.\n");
    fprintf(stderr, "                                           Default 1.\n");
    fprintf(stderr, "    -b,--blockSize         INT         Block size for weighted jackknife.\n");
    fprintf(stderr, "                                           Default 2 MB.\n");
    fprintf(stderr, "    -h,--haplotypeSize     INT         Haplotype size in number of loci.\n");
    fprintf(stderr, "                                           Default 1.\n");
    fprintf(stderr, "    -m,--MAF               DOUBLE      Biallelic sites with MAF <= DOUBLE are dropped.\n");
    fprintf(stderr, "                                           Default 0; monomorphic sites are dropped.\n");
    fprintf(stderr, "    -n,--missingAF         DOUBLE      Sites with proportion of missing genotype >= DOUBLE are dropped.\n");
    fprintf(stderr, "                                           Default 1.\n");
    fprintf(stderr, "    -o,--out               STR         The output basename of the files.\n");
    fprintf(stderr, "\n");
}

// Our options.
static ko_longopt_t long_options[] = {
    {"sampleSize",      ko_required_argument,         'g'},
    {"blockSize",       ko_required_argument,         'b'},
    {"haplotypeSize",   ko_required_argument,         'h'},
    {"MAF",             ko_required_argument,         'm'},
    {"missingAF",       ko_required_argument,         'n'},
    {"out",             ko_required_argument,         'o'},
    {0, 0, 0}
};

// Check that user supplied values are valid.
int check_configuration(PrivateDConfig_t* config) {
    if (config -> sampleSize < 1) {
        fprintf(stderr, "--sampleSize must be given an integer >= 1. Exiting!\n");
        return -1;
    }
    if (config -> blockSize < 1) {
        fprintf(stderr, "--blockSize must be given an integer >= 1. Exiting!\n");
        return -1;
    }
    if (config -> haplotypeSize < 1) {
        fprintf(stderr, "-haplotypeSize must be given an integer >= 1. Exiting!\n");
        return -1;
    }
    if (config -> MAF < 0 || config -> MAF >= 1) {
        fprintf(stderr, "-MAF must be given a real number in [0, 1). Exiting!\n");
        return -1;
    }
    if (config -> missingAF <= 0 || config -> missingAF > 1) {
        fprintf(stderr, "-missingAF must be given a real number in (0, 1]. Exiting!\n");
        return -1;
    }
    if (config -> inputFileName != NULL && access(config -> inputFileName, F_OK) != 0) {
        fprintf(stderr, "%s does not exist. Exiting!\n", config -> inputFileName);
        return -1;
    }
    if (config -> samplesToPopFileName != NULL && access(config -> samplesToPopFileName, F_OK) != 0) {
        fprintf(stderr, "%s does not exist. Exiting!\n", config -> samplesToPopFileName);
        return -1;
    }
    return 0;
}

PrivateDConfig_t* init_privated_config(int argc, char* argv[]) {

    const char *opt_str = "g:b:h:m:n:o:";
    ketopt_t options = KETOPT_INIT;
    int c;

    while ((c = ketopt(&options, argc, argv, 1, opt_str, long_options)) >= 0) {
        switch (c) {
            case ':': fprintf(stderr, "Error! Option %s is missing an argument! Exiting ...\n", argv[options.i - 1]); return NULL;
            case '?': fprintf(stderr, "Error! \"%s\" is unknown! Exiting ...\n", argv[options.i - 1]); return NULL;
        }
	}

    // Set defaults.
    PrivateDConfig_t* config = calloc(1, sizeof(PrivateDConfig_t));
    config -> sampleSize = 1;
    config -> haplotypeSize = 1;
    config -> blockSize = 2000000;
    config -> MAF = 0;
    config -> missingAF = 1;
    config -> inputFileName = NULL;
    config -> samplesToPopFileName = NULL;
    config -> threePopList = NULL;
    config -> cmd = NULL;
    config -> outBaseName = NULL;
    
    // Get user arguments.
    options = KETOPT_INIT;
    while ((c = ketopt(&options, argc, argv, 1, opt_str, long_options)) >= 0) {
        switch (c) {
            case 'g': config -> sampleSize = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 'b': config -> blockSize = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 'h': config -> haplotypeSize = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 'm': config -> MAF = (double) strtod(options.arg, (char**) NULL); break;
            case 'n': config -> missingAF = (double) strtod(options.arg, (char**) NULL); break;
            case 'o': config -> outBaseName = strdup(options.arg); break;
        }
	}

    // If the last three files were not given, then exit.
    if (argc - options.ind != 3) {
        fprintf(stderr, "Incomplete arguments. Exiting!\n");
        destroy_privated_config(config);
        return NULL;
    }

    config -> inputFileName = strdup(argv[options.ind]);
    config -> samplesToPopFileName = strdup(argv[options.ind + 1]);
    config -> threePopList = strdup(argv[options.ind + 2]);

    // Check configuration.
    if (check_configuration(config) != 0) {
        destroy_privated_config(config);
        return NULL;
    }

    // We save the long form of the command to output in the header of files.
    kstring_t* cmd = (kstring_t*) calloc(1, sizeof(kstring_t));
    ksprintf(cmd ,"privated ");
    ksprintf(cmd, "--sampleSize %d ", config -> sampleSize);
    ksprintf(cmd, "--blockSize %d ", config -> blockSize);
    ksprintf(cmd, "--haplotypeSize %d ", config -> haplotypeSize);
    ksprintf(cmd, "--MAF %lf ", config -> MAF);
    ksprintf(cmd, "--missingAF %lf ", config -> missingAF);
    ksprintf(cmd, "%s ", config -> inputFileName);
    ksprintf(cmd, "%s ", config -> samplesToPopFileName);
    ksprintf(cmd, "%s", config -> threePopList);
    config -> cmd = cmd -> s;
    free(cmd);

    if (config -> outBaseName == NULL) {
        int endPos = strcspn(config -> inputFileName, ".vcf");
        config -> outBaseName = strndup(config -> inputFileName, endPos - 1);
    }
    
    return config;
}

void destroy_privated_config(PrivateDConfig_t* config) {
    if (config == NULL)
        return;
    if (config -> inputFileName != NULL)
        free(config -> inputFileName);
    if (config -> samplesToPopFileName != NULL)
        free(config -> samplesToPopFileName);
    if (config -> threePopList != NULL)
        free(config -> threePopList);
    if (config -> cmd != NULL)
        free(config -> cmd);
    if (config -> outBaseName != NULL)
        free(config -> outBaseName);
    free(config);
}