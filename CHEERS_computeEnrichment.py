#!/usr/bin/env python

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
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--ld', help='Directory with LD information for each SNP if CHEERS is used on all SNPs in LD')
    group.add_argument('--snp_list', help='List of SNPs if CHEERS is used on fine-mapped set')
    
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f'{args.input} is not a valid path to a file', file=sys.stderr)
        exit(2)
    
    if args.ld and not os.path.isdir(args.ld):
        print(f'{args.ld} is not a valid path to a directory', file=sys.stderr)
        exit(2)
    
    if args.snp_list and not os.path.isfile(args.snp_list):
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

    return pd.read_csv(file_path, sep='\t').sort_values(['chr', 'start', 'end'], ignore_index=True)


def rank_norm_data(norm_data):

    '''
    For each condition, rank peaks from least specific to most specific.

    :param norm_data: The matrix of peak specificity scores.
    :return: A matrix of peak specificity ranks.
    '''

    peak_info = norm_data[['chr', 'start', 'end']]

    rank_matrix = norm_data.drop(['chr', 'start', 'end'], axis=1).rank(method='min').astype(int)

    return pd.merge(peak_info, rank_matrix, left_index=True, right_index=True)


def load_ld_snps(ld_dir):

    '''
    Load SNPs that are in LD with the tag SNPs.

    :param ld_dir: A path to the directory containing the LD information.
    :return: A data frame of SNPs.
    '''

    ld_dfs = list()

    chr_dirs = os.listdir(ld_dir)

    for chr_dir in chr_dirs:

        lead_snp_files = os.listdir(os.path.join(ld_dir, chr_dir))

        for lead_snp_file in lead_snp_files:

            lead_snp_path = os.path.join(ld_dir, chr_dir, lead_snp_file)
            ld_df = pd.read_csv(lead_snp_path, sep='\t')
            ld_dfs.append(ld_df.iloc[:, [3, 0, 4]])

    snp_list = pd.concat(ld_dfs, axis=1)
    snp_list.columns = ['name', 'chr', 'pos']

    return snp_list


def load_snp_list(file_path):

    '''
    Load SNPs from a file.

    :param file_path: A path to the file containing the list of SNPs.
    :return: A data frame of SNPs.
    '''

    snp_list = pd.read_csv(file_path, sep='\t', header=None)
    snp_list = snp_list.iloc[:, [0, 1, 2]]
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

    while peak_index is not None and snp_index is not None:
    
        if peak.chr == snp.chr:
            if snp.pos >= peak.start and snp.pos <= peak.end:
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


def calc_enrichment(unique_peaks, num_peaks):

    '''
    Use the statistical model to calculate enrichment in each condition.

    :param unique_peaks: Unique peaks with overlapping SNPs.
    :param num_peaks: The number of peaks used to rank the peaks.
    :return: A tuple of observed rank means, p-values, N, n, grand mean, and grand s.d.
    '''

    # TODO: Paper mentions ranks from 1...N but implementation uses 0...N-1
    # I'm subtracting 1 from the mean to be consistent with the original CHEERS
    # Check if this is valid with the statistical model
    observed_rank_means = np.mean(unique_peaks.drop(['chr', 'start', 'end'], axis=1) - 1, axis=0)

    n = len(unique_peaks)

    # TODO: Due to the use of Python 2, integer division returns an integer value
    # I'm casting the result to an integer to be consistent with the original CHEERS
    # Check if this is valid with the statistical model
    mean_sd = math.sqrt(int((num_peaks**2 - 1) / (12 * n)))
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

    # Load SNP list from either LD information or a flat file
    if args.ld:
        snp_list = load_ld_snps(args.ld)
    else:
        snp_list = load_snp_list(args.snp_list)
    
    # Identify overlapping peak-SNP pairs and unique peaks
    peak_overlaps, unique_peaks = overlap_peaks_with_snps(rank_data, snp_list)
    
    # Write overlap information to disk
    peak_overlaps.to_csv(os.path.join(args.outdir, f'{args.trait}_SNPsOverlappingPeaks.txt'), sep='\t', index=None)
    unique_peaks.to_csv(os.path.join(args.outdir, f'{args.trait}_uniquePeaks.txt'), sep='\t', index=None)

    # Calculate p-values based on the statistical model
    observed_rank_means, p_values, N, n, mean_mean, mean_sd = calc_enrichment(unique_peaks, len(norm_data))

    # Write observed rank means and p-values to disk
    p_values.to_csv(os.path.join(args.outdir, f'{args.trait}_disease_enrichment_pValues.txt'), sep='\t', header=False)
    observed_rank_means.to_csv(os.path.join(args.outdir, f'{args.trait}_disease_enrichment_observedMeanRank.txt'), sep='\t', header=False)

    # Store the end time of the program
    end_time = time.time()

    # Write relevant information to disk
    write_log_file(os.path.join(args.outdir, f'{args.trait}.log'), start_time, end_time, len(peak_overlaps), N, n, mean_mean, mean_sd)


if __name__ == '__main__':

    main()
