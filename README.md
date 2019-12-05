# CHEERS
 
CHEERS (Chromatin Element Enrichment Ranking by Specificity) is a method to determine enrichment of annotations in GWAS significant loci. In addition to SNP-peak overlap, CHEERS takes into account peak properties as reflected by quantitative changes in read counts within peaks. 
 
## **Requirements**
 
CHEERS is written for **Python 2.7**. It uses the following modules:
os  
sys  
glob  
re  
subprocess  
time  
numpy  
pandas  
argparse  
scipy  
math  
 
## **How to use CHEERS**

### **create_LD_blocks.py**

Description:  
For a list of SNPs, it finds all the SNPs in LD (r2>0.8).  

*Usage*:  
Create a directory named after the trait. Then move into the directory and run the command.    

mkdir trait_name  
cd trait_name  
python create_LD_blocks.py SNP_LIST OUTPUT_DIR LD_DIR  

We uploaded the GRCh38 LD files to ftp://ngs.sanger.ac.uk/production. They are in the trynka/LD_GRCh38/ directory.

*Output*:  
Parent directory is the name of the trait; within that directory there are subdirectories for all the chromosomes (chr1, chr2...) and within each of these there are .txt files named after the lead SNP (for example results_ld_rs9989735.txt). Each file contains all the SNPs in the LD with the lead (python indexing: snp name is at position 3, chr at position 0 and base pair information is at position 4)

 
### **CHEERS_normalize.py**
 
*Description*:  
This script is used to normalize read counts within peaks. It :  
1)    loads output of featureCounts (*_ReadsInPeaks.txt - each txt file is a sample that contains 4 tab delimited columns without header: chr, start, end, count )
2)    scales read counts to the largest library size
3)    removes the bottom 10th percentile of peaks with the lowest read counts
4)    quantile normalizes the library size-corrected peak counts
5)    performs Euclidean normalization to obtain a cell type specificity score
 
 
*Arguments*:  
  -h, --help   	show this help message and exit  
  --input INPUT	path to files with read counts per peak (default: None)  
  --prefix PREFIX  file prefix (default: None)  
  --outdir OUTDIR  directory where to output results (default: None)  
 
*Outputs*:  
prefix_normToMax.txt  
prefix_normToMax_quantileNorm.txt  
prefix_normToMax_quantileNorm_euclideanNorm.txt  
 
*Usage*:  
python CHEERS_normalize.py --input ~/peak/counts/per/sample/ --prefix prefix --outdir ~/output/directory/
 
 
### **CHEERS_computeEnrichment.py**
 
*Description*:  
 
This script computes the disease enrichment of a provided set of SNPs and writes the output to tab-delimited files. It takes 4 inputs to calculate the enrichment p-value:  
1)    Output from CHEERS_normalize.py. This is text file containing peak coordinates and specificity scores for each of the analyzed samples (output_normToMax_quantileNorm_euclideanNorm.txt). 
2)    Output from create_LD_blocks.py - directory with LD information for each SNP. Alternatively provide txt file with fine-mapped SNP set 
3)    trait name that will be used for output prefix  
4)    Directory where to output results  
 
*Arguments*:  
  -h, --help   	show this help message and exit  
  --input INPUT	Text file containing peak coordinates and specificity
               	scores for each of the analyzed samples (default: None)    
  --ld LD      	Directory with LD information for each SNP (default: None)  
  --snp_list SNP_LIST  list of SNPs if CHEERS is used on finemapped set
                       (default: None)  
  --trait TRAIT	Name of the analyzed trait (default: None)  
  --outdir OUTDIR  Directory where to output results (default: None)  
 
*Outputs*:  
trait_uniquePeaks.txt - list of unique peaks and their ranks per sample  
trait _SNPsOverlappingPeaks.txt - list of all overlapping SNPs and peak ranks per sample  
trait _disease_enrichment_pValues.txt - enrichment p-values per sample  
trait _disease_enrichment_observedMeanRank.txt - observed mean ranks per sample  
trait.log - log file containing run information  
 
*Usage*:  
python CHEERS_computeEnrichment.py --input normToMax_quantileNorm_euclideanNorm.txt --ld ~/LD/trait/ --trait trait_name --outdir ~/output/directory/  
or   
python CHEERS_computeEnrichment.py --input data.txt --snp_list snp_list.txt --trait trait_name --outdir ~/output/directory/  


