# privateD v1.0

Three population introgression test using private allelic richness.

## Building privateD

The resulting executable will be in the **privateD/bin** directory.

```
git clone https://github.com/TQ-Smith/privateD.git 
cd privateD
make
```

## Help
```
Usage: privateD [options] <inFile.vcf.gz> <sampleToPop.tsv> <popList>

<inFile.vcf.gz>             The input VCF file.
<sampleToPop.tsv>           Comma seperate file associating each sample with a population.
<popList>                   Names of the three populations to test seperated by commas.

Options:
    -g                      The maximum standardized sample size used for rarefaction. Default 2.
    -b                      Block size for jackknife. Default 2 MB.
    -h                      Haplotype size in number of loci. Default 1.
```