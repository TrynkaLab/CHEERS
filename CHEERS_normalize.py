#!/usr/bin/env python

'''
:File: CHEERS_normalize.py
:Author: Blagoje Soskic, Wellcome Sanger Institute, <bs11@sanger.ac.uk>
:Updated By: Nikhil Milind, Wellcome Sanger Institute, <nm18@sanger.ac.uk>
:Last Updated: 16 March 2022

This script is used to normalize read counts within peaks. It:

1. Loads an input file containing read counts per peak. The format is provided
    below and can be generated using generate_mock_dataset.py.
2. Scales read counts to the largest library size.
3. Removes the bottom 10th percentile of peaks with the lowest read counts.
4. Quantile normalizes the library size-corrected peak counts.
5. Performs Euclidean normalization to obtain a cell type specificity score.

Outputs:

1. prefix_counts_normToMax.txt - Counts normalized to the largest library size.
2. prefix_counts_normToMax_quantileNorm.txt - Quantile-normalized counts on the
    top 90th percentile of peaks.
3. prefix_counts_normToMax_quantileNorm_euclideanNorm.txt - Euclidean 
    normalized counts that represent cell type specificity scores.

Usage:
    python3 CHEERS_normalize.py trait test_input/ test_input/*_ReadsInPeaks.txt


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
import functools
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
    count_data = count_data.apply(pd.to_numeric)

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
