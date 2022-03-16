#!/usr/bin/env python

'''
:File: CHEERS_computeEnrichment.py
:Author: Blagoje Soskic, Wellcome Sanger Institute, <bs11@sanger.ac.uk>
:Updated By: Nikhil Milind, Wellcome Sanger Institute, <nm18@sanger.ac.uk>
:Last Updated: 16 March 2022

This script computes the disease enrichment of a provided set of SNPs and
writes the output to tab-delimited files. It takes 4 inputs to calculate the
enrichment p-value:

1. Trait name that will be used as a prefix for the output.
2. Directory for result output.
3. Cell type specificity scores that are generated after running 
    CHEERS_normalize.py.
4. A fine-mapped SNP set to test for enrichment. The format is provided below
    and can be generated using generate_mock_dataset.py.

In the original CHEERS implementation, peaks were ranked from 0 to N-1 and the
standard deviation of the uniform distribution was calculated using the floor
of a ratio. This implementation ranks peaks from 1 to N and calculates the
ratio for the standard deviation as a floating point value. If you want to
generate the same p-values that were generated from the original CHEERS
implementation, use the --old_p_values flag.

Usage:
    python3 CHEERS_computeEnrichment.py trait test_output/ test_input/trait_counts_normToMax_quantileNorm_euclideanNorm.txt test_input/snp_list.txt


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

import argparse
import math
import os
import sys
import time

import numpy as np
import pandas as pd
import scipy.stats


def parse_arguments():

    '''
    Parse and validate command line arguments.

    :return: A namespace containing the command line argument values.
    '''

    parser = argparse.ArgumentParser(description='CHEERS computes the disease enrichment within functional annotations by taking into account quantitative changes in read counts within peaks')
    
    parser.add_argument('trait', help='Name of the analyzed trait')
    parser.add_argument('outdir', help='Directory for result output')
    parser.add_argument('input', help='Text file containing peak coordinates and specificity scores for each of the analyzed samples')
    parser.add_argument('snp_list', help='List of SNPs')
    parser.add_argument('--old_p_values', action='store_true', help='Calculate p-values as was done in the original implementation of CHEERS')

    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f'{args.input} is not a valid path to a file', file=sys.stderr)
        exit(2)
    
    if not os.path.isfile(args.snp_list):
        print(f'{args.snp_list} is not a valid path to a file', file=sys.stderr)
        exit(2)

    if not os.path.isdir(args.outdir):
        print(f'{args.outdir} is not a valid path to a directory', file=sys.stderr)
        exit(2)

    return args


def load_data(file_path):

    '''
    Reads the output from the CHEERS normalization process.

    :param file_path: File path to the Euclidean-normalized data.
    :return: A data frame of the peak specificity scores.
    '''

    norm_matrix = pd.read_csv(file_path, sep='\t')
    norm_matrix.sort_values(['chr', 'start', 'end'], ignore_index=True, inplace=True)

    return norm_matrix


def rank_norm_data(norm_data):

    '''
    For each condition, rank peaks from least specific to most specific.

    :param norm_data: The matrix of peak specificity scores.
    :return: A matrix of peak specificity ranks.
    '''

    peak_info = norm_data[['chr', 'start', 'end']]

    rank_matrix = norm_data.drop(['chr', 'start', 'end'], axis=1).rank(method='min').astype(int)

    return pd.merge(peak_info, rank_matrix, left_index=True, right_index=True)


def load_snp_list(file_path):

    '''
    Load SNPs from a file.

    :param file_path: A path to the file containing the list of SNPs.
    :return: A data frame of SNPs.
    '''

    snp_list = pd.read_csv(file_path, sep='\t', header=None)
    snp_list.columns = ['name', 'chr', 'pos']
    snp_list.sort_values(['chr', 'pos'], inplace=True)

    return snp_list


def overlap_peaks_with_snps(peak_info, snp_list):

    '''
    Determine which SNPs overlap with peaks.

    :param peak_info: Peak rank matrix.
    :param snp_list: SNP list to overlap with peak matrix.
    :return: A tuple of peak-SNP overlaps and unique peaks.
    '''

    peak_index = 0
    snp_index = 0

    peak_overlaps = list()
    snp_overlaps = list()

    peak_iterator = peak_info.iterrows()
    snp_iterator = snp_list.iterrows()

    peak_index, peak = next(peak_iterator, (None, None))
    snp_index, snp = next(snp_iterator, (None, None))

    # Iterate until either the peaks or SNPs are all considered for overlap
    while peak_index is not None and snp_index is not None:

        # If the peak occurs before the SNP by genomic coordinates, increment the peak
        # If the SNP occurs before the peak by genomic coordinates, increment the SNP
        if peak.chr == snp.chr:
            if peak.start <= snp.pos <= peak.end:
                peak_overlaps.append(peak_index)
                snp_overlaps.append(snp_index)
                snp_index, snp = next(snp_iterator, (None, None))
            elif snp.pos < peak.start:
                snp_index, snp = next(snp_iterator, (None, None))
            else:
                peak_index, peak = next(peak_iterator, (None, None))
        elif peak.chr < snp.chr:
            peak_index, peak = next(peak_iterator, (None, None))
        else:
            snp_index, snp = next(snp_iterator, (None, None))
    
    snp_list_overlaps = snp_list.drop('chr', axis=1).iloc[snp_overlaps, :]
    snp_list_overlaps.reset_index(drop=True, inplace=True)

    peak_info_overlaps = peak_info.iloc[peak_overlaps, :]
    peak_info_overlaps.reset_index(drop=True, inplace=True)

    overlap_matrix = pd.merge(snp_list_overlaps, peak_info_overlaps, left_index=True, right_index=True)
    overlap_matrix.rename({'name': 'snp'}, axis=1, inplace=True)

    unique_peak_indices = np.unique(peak_overlaps)
    unique_peaks = peak_info.iloc[unique_peak_indices, :]

    return overlap_matrix, unique_peaks


def calc_enrichment(unique_peaks, num_peaks, old_p_values):

    '''
    Use the statistical model to calculate enrichment in each condition.

    :param unique_peaks: Unique peaks with overlapping SNPs.
    :param num_peaks: The number of peaks used to rank the peaks.
    :param old_p_values: Calculate p-values to be consistent with the original CHEERS implementation.
    :return: A tuple of observed rank means, p-values, N, n, grand mean, and grand s.d.
    '''

    # Paper mentions ranks from 1...N but implementation uses 0...N-1
    # I'm subtracting 1 from the mean to be consistent with the original CHEERS
    if old_p_values:
        observed_rank_means = np.mean(unique_peaks.drop(['chr', 'start', 'end'], axis=1) - 1, axis=0)
    else:
        observed_rank_means = np.mean(unique_peaks.drop(['chr', 'start', 'end'], axis=1), axis=0)

    n = len(unique_peaks)

    # Due to the use of Python 2, integer division returns an integer value
    # I'm casting the result to an integer to be consistent with the original CHEERS
    if old_p_values:
        mean_sd = math.sqrt(int((num_peaks**2 - 1) / (12 * n)))
    else:
        mean_sd = math.sqrt((num_peaks**2 - 1) / (12 * n))
    
    mean_mean = (1 + num_peaks) / 2

    p_values = 1 - scipy.stats.norm.cdf(observed_rank_means.values, loc=mean_mean, scale=mean_sd)
    p_values = pd.Series(p_values)
    p_values.index = observed_rank_means.index

    return observed_rank_means, p_values, num_peaks, n, mean_mean, mean_sd


def write_log_file(file_path, start_time, end_time, n_snps, N, n, mean_mean, mean_sd):

    '''
    Writes a log file with relevant information.

    :param file_path: The path of the log file.
    :param start_time: The time the program was started.
    :param end_time: The time the program ended.
    :param n_snps: The number of SNPs overlapping peaks.
    :param N: The number of peaks.
    :param n: The number of peaks with overlapping SNPs.
    :param mean_mean: The grand mean of the ranks.
    :param mean_sd: The grand s.d. of the ranks.
    '''

    running_time = end_time - start_time

    with open(file_path, 'w') as f_out:

        print(f'Total number of peaks\t{N}', file=f_out)
        print(f'Number of overlapping peaks\t{n}', file=f_out)
        print(f'Number of SNPs overlapping peaks\t{n_snps}', file=f_out)
        print(f'Distribution mean\t{mean_mean}', file=f_out)
        print(f'Distribution sd\t{mean_sd}', file=f_out)
        print(f'Running time in seconds\t{running_time}', file=f_out)


def main():

    # Read command line arguments
    args = parse_arguments()

    # Record the start time of the program
    start_time = time.time()

    # Load normalized data generated in a previous step
    norm_data = load_data(args.input)

    # Rank peaks within each condition
    rank_data = rank_norm_data(norm_data)

    # Load SNP list
    snp_list = load_snp_list(args.snp_list)
    
    # Identify overlapping peak-SNP pairs and unique peaks
    peak_overlaps, unique_peaks = overlap_peaks_with_snps(rank_data, snp_list)
    
    # Write overlap information to disk
    peak_overlaps.to_csv(os.path.join(args.outdir, f'{args.trait}_SNPsOverlappingPeaks.txt'), sep='\t', index=None)
    unique_peaks.to_csv(os.path.join(args.outdir, f'{args.trait}_uniquePeaks.txt'), sep='\t', index=None)

    # Calculate p-values based on the statistical model
    observed_rank_means, p_values, N, n, mean_mean, mean_sd = calc_enrichment(unique_peaks, len(norm_data), args.old_p_values)

    # Write observed rank means and p-values to disk
    p_values.to_csv(os.path.join(args.outdir, f'{args.trait}_disease_enrichment_pValues.txt'), sep='\t', header=False)
    observed_rank_means.to_csv(os.path.join(args.outdir, f'{args.trait}_disease_enrichment_observedMeanRank.txt'), sep='\t', header=False)

    # Store the end time of the program
    end_time = time.time()

    # Write relevant information to disk
    write_log_file(os.path.join(args.outdir, f'{args.trait}.log'), start_time, end_time, len(peak_overlaps), N, n, mean_mean, mean_sd)


if __name__ == '__main__':

    main()
