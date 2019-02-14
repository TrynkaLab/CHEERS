'''
:File: CHEERS_normalize.py
:Author: Blagoje Soskic, Wellcome Sanger Institute, <bs11@sanger.ac.uk>
:Last updated: 20/1/2019

This script is used to normalize counts within peaks in the dataset:
1) load output of featureCounts *_ReadsInPeaks.txt (each txt file is an analyzed sample that contains 4 columns without header: chr   start   end count )
2) correct counts per peak for library size by scaling it to the sample with the largest sum of counts
3) remove the bottom 10th percentile of peaks with the lowest read counts
4) quantile normalize the library size-corrected peak counts
5) Euclidean normalization to obtain a cell type specificity score

outputs :
outputPrefix_normToMax.txt
outputPrefix_normToMax_quantileNorm.txt
outputPrefix_normToMax_quantileNorm_euclideanNorm.txt

Usage:
    python CHEERS_normalize.py --input ~/peak/counts/per/sample/ --prefix prefix --outdir ~/output/directory/


Copyright (C) 2019  Blagoje Soskic

This file is part of CHEERS code.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

BY USING THE SOFTWARE YOU ACKNOWLEDGE THAT YOU HAVE READ AND UNDERSTAND THE
TERMS OF USE BELOW. 

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

THIS SOFTWARE IS TO BE USED AS A RESEARCH TOOL ONLY. THE SOFTWARE TOOL SHALL
NOT BE USED AS A DIAGNOSTIC DECISION MAKING SYSTEM AND MUST NOT BE USED TO
MAKE A CLINICAL DIAGNOSIS OR REPLACE OR OVERRULE A LICENSED HEALTH CARE
PROFESSIONAL'S JUDGMENT OR CLINICAL DIAGNOSIS. ANY ACTION THE RESEARCHER TAKES
IN RESPONSE TO THE INFORMATION CONTAINED WITHIN IS AT THE RESEARCHER'S
DISCRETION ONLY.

'''

#import libraries
import os
import glob
import re
import numpy as np
import pandas as pd
import argparse


#parse arguments
parser = argparse.ArgumentParser(description = "data processing for CHEERS disease enrichment", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--input", help = "path to files with read counts per peak")
parser.add_argument("--prefix", help = "file prefix")
parser.add_argument("--outdir", help = "directory where to output results")
args = parser.parse_args()

#define functions
def quantile_norm(normData):
    regions = normData.loc[:, :'end']
    normDataClean = normData.drop(['chr', 'start', 'end'], axis=1)
    rank_mean = normDataClean.stack().groupby(normDataClean.rank(method='first').stack().astype(int)).mean()
    normDataCleanQuantile = normDataClean.rank(method='min').stack().astype(int).map(rank_mean).unstack()
    normDataCleanQuantileFinal = pd.concat([regions, normDataCleanQuantile], axis=1)
    return normDataCleanQuantileFinal

def euclidean_norm(normDataCleanQuantileFinal):
    normDataQuantile = normDataCleanQuantileFinal.drop(['chr', 'start', 'end'], axis=1)
    normDataSquare = np.square(normDataQuantile)
    normDataQuantile['eucl_norm'] = np.sqrt(normDataSquare.sum(axis=1))
    normDataQuantileEucl = normDataQuantile.iloc[:, 0:].div(normDataQuantile.eucl_norm, axis=0)
    dim = normDataQuantileEucl.shape[1] - 1
    final = normDataQuantileEucl.iloc[:, :dim]
    return final

#set variables
filesNames = []
filesSums = []
filesData = {}
peakSums = []

#chrs to keep for the analysis
chrList = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']

# open and load all counts per files
for filename in glob.glob(os.path.join(args.input, '*_ReadsInPeaks.txt')):
    file = open(filename, "r")
    fileText = file.read()
    file.close()
    fileTextLines = fileText.split("\n") #split lines in file
    countSum = 0 # define count sum
    for line in fileTextLines:
        if (line): #check if line exists
            lineParts = line.split("\t")
            if (lineParts[0] in chrList):
                key = lineParts[0] + '-' + lineParts[1] + '-' + lineParts[2] #create key based on chr, start and end
                count = float(lineParts[3])
                if not(key in filesData):
                    filesData[key] = [] #if key doesn't exist create it
                countSum += count
                filesData[key].append(count) #create dictionary with key being location and value being read counts
    filesNames.append(os.path.basename(filename)) #create list of sums
    filesSums.append(countSum) #create list of names

'''
Normalize counts to account for sequencing depth.
'''
#create norm factor
maxSum = max(filesSums) #find max sum
filesSums[:] = [maxSum / x for x in filesSums] #divide every max sum with sum of reads per condition

#normalize each count
for key in filesData:
    counts = filesData[key]
    normalizedCounts = []
    for i in range(0, len(counts)):
        normalizedCounts.append(counts[i] * filesSums[i])
    filesData[key] = normalizedCounts

'''
remove bottom 10th percentile of the peaks
'''
for key in filesData:
    peakCounts = filesData[key]
    sumPeakCounts = sum(filesData[key])
    peakSums.append(sumPeakCounts)
peakSumsArray = np.array(peakSums)
bottom = np.percentile(peakSumsArray, 10)

for (key, value) in filesData.items():
    if sum(value) < bottom:
        del filesData[key]

#create headers for the outputs
normalizedCountsOutput = 'chr' + '\t' + 'start' + '\t' + 'end' + '\t'
for name in filesNames:
    name = re.sub('_ReadsInPeaks\.txt$', '', name)
    normalizedCountsOutput += str(name) + '\t'
normalizedCountsOutput += '\n'

#fill table with the normalized counts
for key in filesData:
    keyParts = key.split('-')
    normalizedCountsOutput += keyParts[0] + '\t' + keyParts[1] + '\t' + keyParts[2] + '\t'
    normalizedCounts = filesData[key]
    for count in normalizedCounts:
        if count == normalizedCounts[:-1]:
            normalizedCountsOutput += str(count)
        else:
            normalizedCountsOutput += str(count) + '\t'
    normalizedCountsOutput += '\n'

#create output for counts normalised for library size and bottom 10th percentile peaks removed
fileName1 = str(args.outdir) + str(args.prefix) + '_counts_normToMax.txt'

text_file = open(fileName1, "w")
text_file.write(normalizedCountsOutput)
text_file.close()

'''
quantile normalize the library size-corrected peak counts to make comparisons between different peaks possible
'''
#quantile normalisation and output
normData = pd.read_csv(fileName1, sep='\t')
normDataQuantile = quantile_norm(normData)
fileName2 = str(args.outdir) + str(args.prefix) + '_counts_normToMax_quantileNorm.txt'
normDataQuantile.to_csv(fileName2, header=1, index=None, sep='\t')

'''
To assess specificity score perform euclidean normalization
'''
#euclidean normalization and output
euclideanNorm = euclidean_norm(normDataQuantile)
regions = normDataQuantile.loc[:, :'end']
euclideanNorm = euclideanNorm.round(4)
euclideanNormFinal = pd.concat([regions, euclideanNorm], axis=1)

fileName3 = str(args.outdir) + str(args.prefix) + '_counts_normToMax_quantileNorm_euclideanNorm.txt'

euclideanNormFinal.to_csv(fileName3, header=1, index=None, sep='\t')
