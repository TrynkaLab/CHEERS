# CHEERS
 
CHEERS (Chromatin Element Enrichment Ranking by Specificity) is a method to determine enrichment of annotations in GWAS significant loci. In addition to SNP-peak overlap, CHEERS takes into account peak properties as reflected by quantitative changes in read counts within peaks.

Main code by Blagoje Soskic.
LD calculation by Lara Bossini-Castillo, Python 3.8 version and test dataset courtesy of [Nikhil Milind](https://github.com/NMilind). 
 
## Requirements
 
This is a rewrite of CHEERS for **Python 3.8**. Please see the `main` branch for a **Python 2.7** version. It uses the following base modules:

1. `argparse`
2. `functools`
3. `gzip`
4. `io`
5. `math`
6. `os`
7. `re`
8. `sys`
9. `time`

In addition, CHEERS requires the following modules:

1. `numpy >= 1.21.3`
2. `pandas >= 1.3.4`
3. `requests >= 2.26.0`
4. `scipy >= 1.8.0`

You can install these using the following command:

```
python3 -m pip install numpy==1.21.3 pandas==1.3.4 requests==2.26.0 scipy==1.8.0
```
 
## How to use CHEERS

### generate_mock_dataset.py

#### Description

Generates a mock ATAC-Seq dataset to run CHEERS on. The output data is in the exact format that CHEERS expects as input. The script pulls data from a publicly-available [ATAC-Seq experiment](https://doi.org/10.1038/s41588-019-0505-9) that has been deposited to the NCBI's SRA.

**Note**: Generating the mock dataset requires an internet connection.

**Note**: The mock data uses pseudorandom sampling of the original data. Set `--seed` to generate reproducible mock data.

#### Usage

The number of peaks (`n_peaks`) is restricted between 1 and 100,000. The number of SNPs (`n_snps`) is restricted between 1 and `n_peaks`. The proportion of SNPs in peaks (`prop_snps_in_peaks`) is a decimal value between 0.0 and 1.0.

```
usage: generate_mock_dataset.py [-h] [--seed SEED] n_peaks n_snps prop_snps_in_peaks out_dir

Generate mock dataset for CHEERS analysis

positional arguments:
  n_peaks             Number of peaks to sample from a real dataset
  n_snps              Number of SNPs to generate within the dataset
  prop_snps_in_peaks  Proportion of SNPs that are in peaks
  out_dir             Output directory

optional arguments:
  -h, --help          show this help message and exit
  --seed SEED         Set a seed for reproducible results
```

#### Outputs

1. `*_ReadsInPeaks.txt` - Read count files that can be generated from an ATAC-Seq experiment over a set of consensus peaks. Each file represents a unique condition.
2. `snp_list.txt` - A list of SNPs, some of which overlap peaks in the mock data.

#### Example

```
mkdir test_input/
python3 generate_mock_dataset.py 10000 1000 0.25 test_input/ --seed 42
```

-----

### CHEERS_normalize.py
 
#### Description

This script is used to normalize read counts within peaks. It:

1. Loads an input file containing read counts per peak. The format is provided below and can be generated using `generate_mock_dataset.py`.
2. Scales read counts to the largest library size.
3. Removes the bottom 10th percentile of peaks with the lowest read counts.
4. Quantile normalizes the library size-corrected peak counts.
5. Performs Euclidean normalization to obtain a cell type specificity score.

#### Usage

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

#### Outputs

1. `prefix_counts_normToMax.txt` - Counts normalized to the largest library size.
2. `prefix_counts_normToMax_quantileNorm.txt` - Quantile-normalized counts on the top 90th percentile of peaks.
3. `prefix_counts_normToMax_quantileNorm_euclideanNorm.txt` - Euclidean normalized counts that represent cell type specificity scores.
 
#### Example

```
mkdir test_input/
python3 generate_mock_dataset.py 10000 1000 0.25 test_input/ --seed 42
python3 CHEERS_normalize.py trait test_input/ test_input/*_ReadsInPeaks.txt
```

#### Input Data Format

This is the expected format of the file containing read counts per peak for a specific condition. The columns are chromosome, start, end, and read counts per peak.

```
chr1  100 200 20
chr2  400 450 15
chr3  500 600 10
```

-----

### CHEERS_computeEnrichment.py

#### Description

This script computes the disease enrichment of a provided set of SNPs and writes the output to tab-delimited files. It takes 4 inputs to calculate the enrichment p-value:

1. Trait name that will be used as a prefix for the output.
2. Directory for result output.
3. Cell type specificity scores that are generated after running `CHEERS_normalize.py`.
4. A fine-mapped SNP set to test for enrichment. The format is provided below and can be generated using `generate_mock_dataset.py`.

In the original CHEERS implementation, peaks were ranked from 0 to N-1 and the standard deviation of the uniform distribution was calculated using the floor of a ratio. This implementation ranks peaks from 1 to N and calculates the ratio for the standard deviation as a floating point value. If you want to generate the same p-values that were generated from the original CHEERS implementation, use the `--old_p_values` flag.

#### Usage

```
usage: CHEERS_computeEnrichment.py [-h] [--old_p_values] trait outdir input snp_list

CHEERS computes the disease enrichment within functional annotations by taking into account quantitative changes in read counts within peaks

positional arguments:
  trait           Name of the analyzed trait
  outdir          Directory for result output
  input           Text file containing peak coordinates and specificity scores for each of the analyzed samples
  snp_list        List of SNPs

optional arguments:
  -h, --help      show this help message and exit
  --old_p_values  Calculate p-values as was done in the original implementation of CHEERS
```

#### Outputs

1. `trait_uniquePeaks.txt` - List of unique peaks and their ranks per sample.
2. `trait_SNPsOverlappingPeaks.txt` - List of all overlapping SNPs and peak ranks per sample.
3. `trait_disease_enrichment_pValues.txt` - Enrichment p-values per sample.
4. `trait_disease_enrichment_observedMeanRank.txt` - Observed mean ranks per sample.
5. `trait.log` - Log file containing run information.
 
#### Example

```
mkdir test_input/
python3 generate_mock_dataset.py 10000 1000 0.25 test_input/ --seed 42
python3 CHEERS_normalize.py trait test_input/ test_input/*_ReadsInPeaks.txt

mkdir test_output/
python3 CHEERS_computeEnrichment.py trait test_output/ test_input/trait_counts_normToMax_quantileNorm_euclideanNorm.txt test_input/snp_list.txt
```

#### Input Data Format

This is the expected format of the file containing the fine-mapped set of SNPs to test for enrichment. The columns are SNP name, chromosome, and position.

```
rs001 chr1  20
rs002 chr2  50
rs003 chr3  100
```
