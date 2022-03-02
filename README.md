# CHEERS
 
CHEERS (Chromatin Element Enrichment Ranking by Specificity) is a method to determine enrichment of annotations in GWAS significant loci. In addition to SNP-peak overlap, CHEERS takes into account peak properties as reflected by quantitative changes in read counts within peaks. 
 
## Requirements
 
CHEERS is written for **Python 3.8**. It uses the following modules:

1. `os`
2. `sys`
3. `subprocess`
4. `time`
5. `numpy`
6. `pandas`
7. `argparse`
8. `functools`
9. `scipy`
10. `math`
 
## How to use CHEERS

### create_LD_blocks.py

Description:

For a list of SNPs, it finds all the SNPs in LD (r2 > 0.8).

*Usage*:

```
usage: create_LD_blocks.py [-h] snp_map outdir tabix_dir

Tagging SNP identification for CHEERS disease enrichment

positional arguments:
  snp_map     File containing SNPs to query for LD information
  outdir      Directory for result output
  tabix_dir   Directory containing LD tabix files
```

We uploaded the GRCh38 LD files to `ftp://ngs.sanger.ac.uk/production`. They are in the `trynka/LD_GRCh38/` directory.

*Outputs*:

Parent directory is the name of the trait; within that directory there are subdirectories for all the chromosomes (chr1, chr2...) and within each of these there are `.txt` files named after the lead SNP (for example `results_ld_rs9989735.txt`). Each file contains all the SNPs in the LD with the lead (python indexing: snp name is at position 3, chr at position 0 and base pair information is at position 4).

*Example*:

Create a directory named after the trait. Then move into the directory and run the command.

```
mkdir trait_name
cd trait_name
python create_LD_blocks.py SNP_LIST OUTPUT_DIR LD_DIR
```

-----

### CHEERS_normalize.py
 
*Description*:

This script is used to normalize read counts within peaks. It:

1. loads output of featureCounts (each txt file is a sample that contains 4 tab delimited columns without header: chr, start, end, count)
2. scales read counts to the largest library size
3. removes the bottom 10th percentile of peaks with the lowest read counts
4. quantile normalizes the library size-corrected peak counts
5. performs Euclidean normalization to obtain a cell type specificity score

*Usage*:

```
usage: CHEERS_normalize.py [-h] prefix outdir [input [input ...]]

Data processing for CHEERS disease enrichment

positional arguments:
  prefix      File prefix
  outdir      Directory for result output
  input       Pattern of files with read counts per peak

optional arguments:
  -h, --help  show this help message and exit
```

*Outputs*:

1. `prefix_counts_normToMax.txt`
2. `prefix_counts_normToMax_quantileNorm.txt`
3. `prefix_counts_normToMax_quantileNorm_euclideanNorm.txt`
 
*Example*:

```
python CHEERS_normalize.py prefix ~/output/directory ~/peak/counts/per/sample/*.txt
```

-----

### CHEERS_computeEnrichment.py
 
*Description*:  
 
This script computes the disease enrichment of a provided set of SNPs and writes the output to tab-delimited files. It takes 4 inputs to calculate the enrichment p-value:

1. Output from CHEERS_normalize.py. This is the text file containing peak coordinates and specificity scores for each of the analyzed samples(`prefix_count_normToMax_quantileNorm_euclideanNorm.txt`).
2. Output from `create_LD_blocks.py` - directory with LD information for each SNP. Alternatively provide txt file with fine-mapped SNP set.
3. trait name that will be used for output prefix
4. Directory where to output results

*Usage*:

```
usage: CHEERS_computeEnrichment.py [-h] (--ld LD | --snp_list SNP_LIST) trait outdir input

CHEERS computes the disease enrichment within functional annotations by taking into account quantitative changes in read counts within peaks

positional arguments:
  trait                Name of the analyzed trait
  outdir               Directory for result output
  input                Text file containing peak coordinates and specificity scores for each of the analyzed samples

optional arguments:
  -h, --help           show this help message and exit
  --ld LD              Directory with LD information for each SNP if CHEERS is used on all SNPs in LD
  --snp_list SNP_LIST  List of SNPs if CHEERS is used on fine-mapped set
```

*Outputs*:

1. `trait_uniquePeaks.txt` - list of unique peaks and their ranks per sample
2. `trait_SNPsOverlappingPeaks.txt` - list of all overlapping SNPs and peak ranks per sample
3. `trait_disease_enrichment_pValues.txt` - enrichment p-values per sample
4. `trait _disease_enrichment_observedMeanRank.txt` - observed mean ranks per sample
5. `trait.log` - log file containing run information
 
*Usage 1*:

```
python CHEERS_computeEnrichment.py trait_name ~/output/directory/ normToMax_quantileNorm_euclideanNorm.txt --ld ~/LD/trait/
```

*Usage 2*:

```
python CHEERS_computeEnrichment.py trait_name ~/output/directory/ normToMax_quantileNorm_euclideanNorm.txt --snp_list snp_list.txt
```
