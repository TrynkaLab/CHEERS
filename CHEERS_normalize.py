#!/usr/bin/env python

import argparse
import functools
import glob
import os
import sys

import numpy as np
import pandas as pd


def parse_arguments():

    '''
    Parse and validate command line arguments.

    :return: A namespace containing the command line argument values.
    '''

    parser = argparse.ArgumentParser(description='Data processing for CHEERS disease enrichment')
    parser.add_argument('prefix', help='File prefix')
    parser.add_argument('outdir', help='Directory for result output')
    parser.add_argument('input', nargs='*', help='Pattern of files with read counts per peak')
    args = parser.parse_args()

    if not args.input:
        print(f'Input did not match any valid paths to files', file=sys.stderr)
        exit(2)
    
    if not os.path.isdir(args.outdir):
        print(f'{args.outdir} is not a valid path to a directory', file=sys.stderr)
        exit(2)

    return args


def read_files(files):

    '''
    Read input count files and generate a data frame.

    :param files: A list of file paths.
    :return: A count matrix, where each column is for one condition and each row represents a peak.
    '''

    data_in_dfs = list()

    for file_name in files:

        if not os.path.isfile(file_name):
            print(f'Could not open file {file_name}', file=sys.stderr)
            exit(1)
        
        condition_name = os.path.basename(file_name).split('.')[0]

        data_in = pd.read_csv(file_name, sep='\t', header=None)
        data_in.columns = ['chr', 'start', 'end', condition_name]

        data_in_dfs.append(data_in)

    count_data = functools.reduce(
        lambda left, right: pd.merge(left, right, on=['chr', 'start', 'end'], how='inner'),
        data_in_dfs
    )

    count_data.set_index(count_data.chr + '-' + count_data.start.map(str) + '-' + count_data.end.map(str), inplace=True)
   
    count_data.drop(['chr', 'start', 'end'], axis=1, inplace=True)

    return count_data


def normalize_counts(count_data):

    '''
    Normalize counts to account for sequencing depth.

    :param count_data: A data frame from Peak ID to list of counts across conditions.
    :return: Normalized count data.
    '''

    count_sums = count_data.sum(axis=0)

    norm_factors = count_sums.max() / count_sums

    return count_data * norm_factors


def drop_peaks(norm_data):

    '''
    Drops the bottom 10th percentile of peaks by total count.

    :param norm_data: A normalized peak count matrix.
    :return: A normalized peak count matrix with the top 90% of peaks based on total count.
    '''

    peak_sums = norm_data.sum(axis=1)

    tenth_percentile = np.percentile(peak_sums, 10)

    return norm_data[peak_sums >= tenth_percentile]


def quantile_norm(norm_data):

    '''
    Quantile normalize each condition's peak count distributions.

    :param norm_data: A normalized peak count matrix.
    :return: A quantile-normalized peak count matrix.
    '''

    rank_mean = norm_data.stack().groupby(norm_data.rank(method='first').stack().astype(int)).mean()

    return norm_data.rank(method='min').stack().astype(int).map(rank_mean).unstack()


def euclidean_norm(norm_data):

    '''
    Normalize each peak's count by the euclidean norm across all conditions.

    :param norm_data: A normalized peak count matrix.
    :return: A Euclidean-normalized peak count matrix.
    '''

    norm_data_squared = np.square(norm_data)

    euclidean_norm = np.sqrt(norm_data_squared.sum(axis=1))

    return norm_data.div(euclidean_norm, axis=0)


def write_matrix(matrix, file_path):

    '''
    Writes any matrix generated in this file. The matrix should have an index of peak IDs.

    :param matrix: A matrix to write out.
    :param file_path: The file path for the matrix output.
    '''

    out_mtx = matrix.copy()

    peak_info = pd.Series(out_mtx.index).str.split('-', expand=True)
    peak_info.columns = ['chr', 'start', 'end']
    peak_info.start = peak_info.start.astype(int)
    peak_info.end = peak_info.end.astype(int)
    peak_info.index = out_mtx.index

    out_mtx = pd.merge(peak_info, out_mtx, left_index=True, right_index=True)

    out_mtx.sort_values(['chr', 'start', 'end'], inplace=True)

    out_mtx.to_csv(file_path, index=None, sep='\t')


def main():

    # Read command line arguments
    args = parse_arguments()

    # Generate one unified count matrix from input files
    count_data = read_files(args.input)

    # Normalize data based on library size
    norm_data = normalize_counts(count_data)

    # Drop bottom 10th percentile of peaks
    norm_data = drop_peaks(norm_data)

    # Write library-normalized data to disk
    norm_data_path = os.path.join(args.outdir, f'{args.prefix}_counts_normToMax.txt')
    write_matrix(norm_data, norm_data_path)

    # Perform quantile normalization
    quantile_norm_data = quantile_norm(norm_data)

    # Write quantile-normalized data to disk
    quantile_norm_data_path = os.path.join(args.outdir, f'{args.prefix}_counts_normToMax_quantileNorm.txt')
    write_matrix(quantile_norm_data, quantile_norm_data_path)

    # Perform Euclidean normalization to generate specificity scores for peaks
    euclidean_norm_data = euclidean_norm(quantile_norm_data)

    # Write Euclidean-normalized data to disk
    euclidean_norm_data_path = os.path.join(args.outdir, f'{args.prefix}_counts_normToMax_quantileNorm_euclideanNorm.txt')
    write_matrix(euclidean_norm_data.round(4), euclidean_norm_data_path)


if __name__ == '__main__':

    main()
