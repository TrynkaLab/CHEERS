'''
:File: CHEERS_computeEnrichment.py
:Author: Blagoje Soskic, Wellcome Sanger Institute, <bs11@sanger.ac.uk>
:Last updated: 20/1/2019

This code computes the disease enrichment of provided set of SNPs and writes the output to tab-delimited files.
It takes 4 inputs to calculate the enrichment p-value
inputs:
1) Text file containing peak coordinates and specificity scores for each of the analyzed samples
2) Directory with LD information for each SNP
3) Name of the analyzed trait
4) Directory where to output results

usage:
python CHEERS_computeEnrichment.py --input data.txt --ld ~/LD/trait/ --trait trait_name --outdir ~/output/directory/
'''

import os
from scipy.stats import norm
import numpy
import math
import argparse
import time
import sys

#parse arguments
parser = argparse.ArgumentParser(description = "CHEERS computes the disease enrichment within functional annotations by taking into account quantitative changes in read counts within peaks", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--input", help = "Text file containing peak coordinates and specificity scores for each of the analyzed samples")
parser.add_argument("--ld", help = "Directory with LD information for each SNP if CHEERS is used on all SNPs in LD")
parser.add_argument("--snp_list", help = "list of SNPs if CHEERS is used on finemapped set")
parser.add_argument("--trait", help = "Name of the analyzed trait")
parser.add_argument("--outdir", help = "Directory where to output results")
args = parser.parse_args()

#set the timer
start_time1 = time.time()

#stop the execution if all the parameters are not there
if args.input is None:
    print ("Please provide the file that is output of CHEERS_normalize.py.")
    sys.exit()

if (args.ld is None and args.snp_list is None):
    print ("Please provide path to LD or SNP-list.")
    sys.exit(1)

if (args.ld is not None and args.snp_list is not None):
    print ("Please provide only path to LD or SNP-list.")
    sys.exit(1)

if args.trait is None:
    print("Please provide the trait name for the output files.")
    sys.exit(1)

if args.outdir is None:
    print("Please provide the path to directory to output files.")
    sys.exit(1)

#set variables
normValues = {}
normValuesSorted = {}
samplesList = []
snpsList = []
peaks = {}
overlappedPeaks ={}
overlappedSNPs = {}
uniquePeaks = {}
pValue = {}

'''
open file containing normalized counts per region. 
output of the code: CHEERS_normalize.py
'''

file = open(args.input, "r")
counts = file.read()
file.close()
countsLines = counts.split("\n")

#load data
'''
load data into dictionary that has sample as the key and list of peaks (each being the separate dictionary) as a value
example: 
{'Sample 1': [{'start': S1, 'chr1': 'chr1', 'end': E1, 'value': V1}, {'start': S2, 'chr': 'chr1', 'end': E2, 'value':V2} ... }
'''
headerNotRead = True #set whether to read the first row
for countsLine in countsLines:
    if (countsLine): #only load if string is not empty
        countsLineParts = countsLine.split("\t")
        if (not (headerNotRead)):
            chr = countsLineParts[0]
            start = countsLineParts[1]
            end = countsLineParts[2]
            for i in range(3, len(countsLineParts)):
                normValues[samplesList[i-3]].append({
                    'chr': chr,
                    'start': float(start),
                    'end': float(end),
                    'value': float(countsLineParts[i].replace('\r', ''))
                })
        else:
            headerNotRead = False  # now first raw is already read
            for i in range(3, len(countsLineParts)):
                stripped_line = countsLineParts[i].replace('\r', '')
                normValues[stripped_line] = []
                samplesList.append(stripped_line)

''' 
Ranking step
1. go over the list of samples
2. sort euclidean norm values from the lowest to the highest within each sample
3. rank sorted values (take values from the iterator in the for cycle)
4. if the value already exist rank them all the same with the lowest rank i.e (1, 2, 3, 3, 3, 4, 5 ) -> (1, 2, 3, 3, 3, 6, 7)
 '''
for sample in samplesList:
    normValues[sample] = sorted(normValues[sample], key=lambda parameter: parameter['value'])
    lastValue = None
    lastRank = None
    for rank in range(0, len(normValues[sample])):
        if (not (lastValue)):
            lastValue = normValues[sample][rank]["value"]
            lastRank = float(rank)
            normValues[sample][rank]['rank'] = lastRank
        else:
            if (lastValue == normValues[sample][rank]['value']):
                normValues[sample][rank]['rank'] = lastRank
            else:
                lastValue = normValues[sample][rank]["value"]
                lastRank = float(rank)
                normValues[sample][rank]['rank'] = lastRank


#read all snp files containing the LD information

'''
- if using LD 
- directory structure: parent directory is a name of the trait; within that directory there are directories for all the chrs and within
each chr directory there are .txt files named after lead SNP (for example results_ld_rs9989735.txt)
- each file contains all the SNPs in the LD (snp name is at position 3, chr at position 0 and bp info is at position 4)
- load it into the list of dictionaries with the structure: [{'chr': 'chr1', 'name': 'SNP1', 'pos': P1}, {'chr': 'chr1', 'name': 'SNP2', 'pos': P2}, ... ]
'''

if args.ld is not None:
    ldPerChr = os.listdir(args.ld)
    for chrDir in ldPerChr:
        leadSnps = os.listdir(os.path.join(args.ld, chrDir))
        for leadSnp in leadSnps:
            leadSnpFile = open(os.path.join(args.ld, chrDir, leadSnp), "r")
            leadSnp = leadSnpFile.read()
            leadSnpFile.close()
            leadSnpLines = leadSnp.split("\n")
            chr = None
            for snpInLD in leadSnpLines:
                if (snpInLD):
                    snpInLDParts = snpInLD.split("\t")
                    if not(chr):
                        chr = snpInLDParts[0]
                    position = float(snpInLDParts[4])
                    # create dictionary containing name, position and chr
                    snpsList.append({
                        'name': snpInLDParts[3],
                        'pos': position,
                        'chr': snpInLDParts[0]
                    })

if args.snp_list is not None:
    ldFile = open(args.snp_list, "r")
    ld = ldFile.read()
    ldFile.close()
    ldParts = ld.split("\n")
    chr = None
    for snpLd in ldParts:
        if (snpLd):
            ldParts = snpLd.split("\t")
            if not (chr):
                chr = ldParts[1]
            position = float(ldParts[2])
            snpsList.append({
                'name': ldParts[0],
                'pos': position,
                'chr': ldParts[1]
            })

# restructure the data
'''
to speed up the overlaps convert list of normalised values into dictionary with key being chr and values being a list of peaks
current structure: {'Sample 1': [{'start': S1, 'chr1': 'chr1', 'end': E1, 'value': V1}, {'start': S2, 'chr': 'chr1', 'end': E2, 'value':V2} ... }
new structure: {'chr5': [{'start': S1, 'chr': 'chr5', 'end': E1, 'rank': R1, 'sampleName': 'S1'}, {'start': S2, 'chr': 'chr5', 'end': E2, 'rank': R2, 'sampleName': 'S2'}, ...]...}
this will allow to only search peaks that are on the same chr as the SNP
'''
for sampleName in samplesList:
    for peak in normValues[sampleName]:
        chr = peak['chr']
        if not (chr in peaks):  # if the key doesn't exist create it
            peaks[chr] = []
        peaks[chr].append({
            'chr': chr,
            'start': peak["start"],
            'end': peak["end"],
            'rank': peak["rank"],
            'sampleName': sampleName
        })

#  overlap snps with peaks
'''
first match chr (key of dictionary peaks) with the chr for SNP, and only iterate through the peaks on the same chr
output is the dictionary with key being the SNP and value being the list of peak values for the sample
output example: {'SNP1': [{'start': S1, 'chr': 'chr4', 'end': E1, 'rank': R1, 'sampleName': 'S1'}, {'start': S1, 'chr': 'chr4', 'end': E1, 'rank': R1, 'sampleName': 'S2'}...]}
'''

for snp in snpsList:
    chr = snp['chr']
    pos = snp['pos']
    name = snp['name']
    if (chr in peaks):
        peaksPerChr = peaks[chr] # gives all the peaks that are on the same chr as SNP
        for peak in peaksPerChr:
            if ((pos >= peak["start"]) and (pos <= peak["end"])):
                if not (name in overlappedPeaks):  # if the key doesn't exist create it
                    overlappedPeaks[name] = []
                    overlappedSNPs[name] = snp
                overlappedPeaks[name].append(peak)

# create headers for the outputs
'''
SNPsOverlappingPeaks: snp pos	chr	start	end Sample1	Sample2	...
peaksOverlappingSNP: chr	start	end Sample1	Sample2	...
'''
SNPsOverlappingPeaks = 'snp' + '\t' + 'pos' + '\t' + 'chr' + '\t' + 'start' + '\t' + 'end' + '\t' # create output containing all SNPs overlapping peaks
peaksOverlappingSNP = 'chr' + '\t' + 'start' + '\t' + 'end' + '\t' # create the output for all the peaks overlapping SNP (multiple SNPs can hit the same peak so each peak is otput only once)

for name in samplesList:
    SNPsOverlappingPeaks += name + '\t'
    peaksOverlappingSNP += name + '\t'
SNPsOverlappingPeaks += '\n'
peaksOverlappingSNP += '\n'

'''
fill in the positions and normalised ranks for SNPsOverlappingPeaks
'''

for snp in overlappedPeaks: #take each SNP (key in the overlappedPeaks) that overlaps peaks
    overlappedPeak = overlappedPeaks[snp] #take peak that overlapps SNP
    SNPsOverlappingPeaks += snp + '\t' #add first SNP name
    SNPsOverlappingPeaks += str(int(overlappedSNPs[snp]['pos'])) + '\t'
    SNPsOverlappingPeaks += overlappedPeak[0]['chr'] + '\t' + str(int(overlappedPeak[0]['start'])) + '\t' + str(int(overlappedPeak[0]['end'])) + '\t' #since coordinates are the same for all the samples only add them once
    for peak in overlappedPeak:
        SNPsOverlappingPeaks += str(int(peak['rank'])) + '\t' #add ranks
    SNPsOverlappingPeaks += '\n'


'''
create the output of all the SNPs and the overlapping peaks
'''

name = str(args.outdir) + str(args.trait) + '_SNPsOverlappingPeaks.txt'
textFile = open(name, "w")
textFile.write(SNPsOverlappingPeaks)
textFile.close()

# fill in the positions and normalised ranks for peaksOverlappingSNP
'''
multiple SNPs in LD can hit the same peak so each peak would contribute to the score multiple times.
to remove "duplicated peaks" create key chr_start_end and restructure the data so the key is coordinate and value is the list of all ranks
for example: {'chr1_S1_E1': [Sample1_rank, Sample2_rank, Sample3_rank...], 'chr1_S2_E2': [Sample1_rank, Sample2_rank, Sample3_rank...]... }
'''
for snp in overlappedPeaks:
    overlappedPeak = overlappedPeaks[snp]
    start = overlappedPeak[0]['start']
    end = overlappedPeak[0]['end']
    chr = overlappedPeak[0]['chr']
    key = chr + '_' + str(start) + "_" + str(end)
    if not (key in uniquePeaks):
        uniquePeaks[key] = []
        for data in overlappedPeak:
            uniquePeaks[key].append(data['rank'])

'''
go over unique peaks
'''
for key in uniquePeaks:
    overlappedPeak = uniquePeaks[key]
    keyList = key.split("_")
    chr = keyList[0]
    start = float(keyList[1])
    end = float(keyList[2])
    peaksOverlappingSNP += chr + '\t' + str(int(start)) + '\t' + str(int(end)) + '\t'
    for rank in overlappedPeak:
        peaksOverlappingSNP += str(int(rank)) + '\t'
    peaksOverlappingSNP += '\n'

name2 = str(args.outdir) + str(args.trait) + '_uniquePeaks.txt'

'''
create the output of all the peaks that are overlapping SNP
'''
textFile2 = open(name2, "w")
textFile2.write(peaksOverlappingSNP)
textFile2.close()

'''
get the mean of ranks and convert it to dictionary
'''
observedMean = []
for value in (zip(*list(uniquePeaks.values()))):
    observedMean.append(numpy.mean(value))
observedMean = dict(zip(samplesList, observedMean))

'''
define parameters for descrete uniform distribution and calculate p-value
'''
N = len(normValues[sample])
n = len(uniquePeaks)
mean_sd = math.sqrt((N**2-1)/(12*n))
mean_mean = (1+N)/2

for sample in samplesList:
    observed = observedMean[sample]
    if not (sample in pValue):
        pValue[sample] = ()
    pValue[sample] = 1-norm.cdf(observed, loc=mean_mean, scale=mean_sd)

'''
create the txt file with the p-values
'''
pValueName = str(args.outdir) + str(args.trait) + '_disease_enrichment_pValues.txt'
f = open(pValueName, 'a')
for key, value in pValue.iteritems():
    f.write(key + '\t' + str(value) + '\n')
f.close()

'''
create the txt file with the mean ranks
'''
meanName = str(args.outdir) + str(args.trait) + '_disease_enrichment_observedMeanRank.txt'
f = open(meanName, 'w')
for key,value in observedMean.iteritems():
    f.write(key + '\t' + str(value) + '\n')
f.close()

'''
output in the log file
'''

end_time1 = time.time()
running_time = (end_time1 - start_time1)

logfileName = str(args.outdir) + str(args.trait) + ".log"
logfile = open(logfileName, "w")
print >> logfile, 'Total number of peaks\t%s' % (str(N))
print >> logfile, 'Number of overlapping peaks\t%s' % (str(n))
print >> logfile, 'Number of SNPs overlapping peaks\t%s' % (str(len(overlappedPeaks)))
print >> logfile, 'Distribution mean\t%s' %  (str(mean_mean))
print >> logfile, 'Distribution sd\t%s' %  (str(mean_sd))
print >> logfile, 'Running time in seconds\t%s' % (running_time)