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
Usage: privateD [options] <inFile.vcf.gz> <sampleToPop.tsv> <pop1>,<pop2>,<pop3>

<inFile.vcf.gz>             The input VCF file.
<sampleToPop.tsv>           Comma seperate file associating each sample with a population.
<pop1>,<pop2>,<pop3>        Names of the three populations in <sampleToPop.csv.gz>.

Options:
    -g                      The maximum standardized sample size used for rarefaction. Default 2.
    -b                      Block size for jackknife. Default 2 MB.
    -h                      Haplotype size in number of loci. Default 1.
```