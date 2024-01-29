# barcodeIBD

Run `barcodeIBD.R` to calculate pairwise comparisons of genetic SNP-barcodes using [hmmIBD](https://github.com/glipsnort/hmmIBD).

The script contains the global variables that you need to modify to fit your own data.

# Install

The pipeline was run with R version 4.1.3 and requires the packages:

- reshape2
- stringr
- ape
- readr
- dplyr
- lubridate
- igraph

# Input

Two files are required to compute the IBD between genetic SNP-barcodes.

## *barcodefile*

**The location of the SNP-barcodes file.**

This tabulated file contains the unique identifier of each SNP-barcode in each sample followed by the base called at each loci separated in individual columns.

Multi-allelic SNPs can be either 'A', 'C', 'G', 'T', 'N' or 'X'. 'N' represents a mixed call of several bases at a locus and 'X' represents an unknown base at a locus.

The first line of the file should be the column names, with no restriction for the first column but with names under the form 'chromosome.position' for the loci columns.

Example:

|ID|chr1.123456|chr2.456789|chrX.789123|
|---|---|---|---|
|SAMPLE1|A|T|G|
|SAMPLE2|X|N|T|

## *metadatafile*

**The location of the SNP-barcodes metadata file.**

This tabulated file contains the identifier of each barcode followed by various information about sample collection. At least four columns are expected in the following order:

1. Unique identifier for each SNP-barcode
2. Unique individual identifier
3. Date of collection (YYYY-MM-DD)
4. Location of collection

Additional columns can also be added.

The file should contain a header with no restriction on the names of any column.

Example:

|ID|individual|date|location|
|---|---|---|---|
|SAMPLE1|IND1|2024-01-28| Montpellier |
|SAMPLE2|IND1|2024-01-29|Montpellier|

# Options

## *chrOrder*

**The order of chromosomes when transforming them into integers.** 

If you omit a chromosome in this list, all of its loci will be discarded when building the relatedness network. You can set this variable to 'NA' if you do not require the chromosome to have a specific order.

## *outDir*

**The path to the directory where the output of this pipeline will be saved.**

## *hmmIBD*

**The location of the compiled hmmIBD binary.**

You need to install [hmmIBD](https://github.com/glipsnort/hmmIBD) prior to setting this variable.

## *runInBackground*

**Whether hmmIBD should run in a background process.**

As hmmIBD can take a substantial amount of time to run, especially with many samples or SNPs, you may want this process to run in background so that it will still run if the R session crashes. If that happened, you would need to manually zip hmmIBD output files so that the rest of pipeline can read them.

```bash
zip -j "barcodes-IBD.hmm_fract.txt.zip" "barcodes-IBD.hmm_fract.txt" && rm "barcodes-IBD.hmm_fract.txt"
zip -j "barcodes-IBD.hmm.txt.zip" "barcodes-IBD.hmm.txt" && rm "barcodes-IBD.hmm.txt"
```

## *minsites*

**The minimal number of sites to consider the IBD value of a pair.**

## *clusterIBDmin*

**The value of IBD above which SNP-barcodes can be clustered.**

For example, if you set this variable to 0.5, clusters will contain SNP-barcodes related with an IBD of 0.5 or more only. Clusters are not complete subgraphs (cliques) but rather graph components of SNP-barcodes.
