#!/usr/bin/env python

'''
:File: generate_mock_dataset.py
:Author: Nikhil Milind, Wellcome Sanger Institute, <nm18@sanger.ac.uk>
:Last Updated: 16 March 2022

This script generates a mock dataset to run through CHEERS. This includes read
count data for peaks called in different cell conditions and a SNP list to test
for enrichment.

Usage:
    python3 generate_mock_dataset.py 10000 1000 0.25 test_input/ --seed 42


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
import gzip
import io
import os
import re
import requests

import numpy as np
import pandas as pd


def parse_arguments():

    '''
    Parse and validate command line arguments.

    :return: A namespace containing the command line argument values.
    '''

    parser = argparse.ArgumentParser(description='Generate mock dataset for CHEERS analysis')

    parser.add_argument('n_peaks', help='Number of peaks to sample from a real dataset')
    parser.add_argument('n_snps', help='Number of SNPs to generate within the dataset')
    parser.add_argument('prop_snps_in_peaks', help='Proportion of SNPs that are in peaks')
    parser.add_argument('out_dir', help='Output directory')
    parser.add_argument('--seed', help='Set a seed for reproducible results', default=42)

    args = parser.parse_args()

    if not args.n_peaks.isdigit():
        print(f'{args.n_peaks} is not a positive integer', file=sys.stderr)
        exit(2)
    
    args.n_peaks = int(args.n_peaks)

    if args.n_peaks > 100000 or args.n_peaks < 1:
        print('Please set n_peaks between 1 and 100000', file=sys.stderr)
        exit(2)
    
    if not args.n_snps.isdigit():
        print(f'{args.n_snps} is not a positive integer', file=sys.stderr)
        exit(2)
    
    args.n_snps = int(args.n_snps)

    if args.n_snps > args.n_peaks or args.n_snps < 1:
        print('Please set n_snps between 1 and the number of peaks', file=sys.stderr)
        exit(2)
    
    try:
        args.prop_snps_in_peaks = float(args.prop_snps_in_peaks)
    except ValueError:
        print(f'{args.prop_snps_in_peaks} is not a float', file=sys.stderr)
        exit(2)
    
    if args.prop_snps_in_peaks > 1.0 or args.prop_snps_in_peaks < 0.0:
        print('Please set prop_snps_in_peaks between 0.0 and 1.0', file=sys.stderr)
        exit(2)

    if not os.path.isdir(args.out_dir):
        print(f'{args.out_dir} is not a valid path to a directory', file=sys.stderr)
        exit(2)
    
    try:
        args.seed = int(args.seed)
    except ValueError:
        print(f'{args.seed} is not an integer', file=sys.stderr)
        exit(2)

    return args


def main():

    # Read command line arguments
    args = parse_arguments()

    # Set seed for reproducibility
    np.random.seed(args.seed)

    # Calculate number of SNPs inside and outside peaks
    n_snps_in_peaks = int(np.floor(args.n_snps * args.prop_snps_in_peaks))
    n_snps_outside_peaks = int(args.n_snps - n_snps_in_peaks)

    # ATAC-Seq data that is publicly available
    #   Reference DOI: 10.1038/s41588-019-0505-9
    #   GEO Accession: GSE118189
    url = 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE118189&format=file&file=GSE118189%5FATAC%5Fcounts%2Etxt%2Egz'

    # Load data directly into memory
    file_content = requests.get(url).content
    file_io = io.BytesIO(file_content)

    # Read peak count matrix content
    with gzip.GzipFile(fileobj=file_io, mode='rb') as data_in:
        peak_counts = pd.read_csv(io.StringIO(data_in.read().decode('UTF-8')), sep='\t')

    # Sample a subset of the peaks
    peak_counts = peak_counts.sample(args.n_peaks)

    # Sum together counts from the same conditions
    peak_counts = peak_counts.transpose()
    peak_counts['Condition'] = [re.sub('^[0-9]*-', '', x) for x in peak_counts.index.values]
    condition_counts = peak_counts.groupby('Condition').sum().transpose()

    # Generate peak information
    peak_info = pd.Series(condition_counts.index).str.split('_', expand=True)
    peak_info.columns = ['Chr', 'Start', 'End']
    peak_info = peak_info.set_index(condition_counts.index)

    # Write out each condition's counts
    for i in range(condition_counts.shape[1]):
        out_counts = pd.concat((peak_info, condition_counts.iloc[:, i]), axis=1)
        out_counts = out_counts.sort_values(['Chr', 'Start', 'End'])
        out_path = os.path.join(args.out_dir, f'{condition_counts.columns[i]}_ReadsInPeaks.txt')
        out_counts.to_csv(out_path, sep='\t', index=False, header=False)

    # Create SNP list
    snp_list = dict()
    snp_list['Chr'] = list()
    snp_list['Position'] = list()

    # Sample peaks to generate SNPs within peaks
    peaks_with_snps = peak_info.sample(n_snps_in_peaks)
    for i, row in peaks_with_snps.iterrows():
        snp_list['Chr'].append(str(row['Chr']))
        start = int(row['Start'])
        end = int(row['End'])
        snp_list['Position'].append(int(((end - start) / 2) + start))

    # Generate SNPs greater than 400Mb to guarantee no overlap with peaks
    for i in range(n_snps_outside_peaks):
        snp_list['Chr'].append(peak_info['Chr'].sample().values[0])
        snp_list['Position'].append(np.random.randint(400e6, 450e6))

    # Generate SNP names
    snp_list['SNP'] = [f'rs{i}' for i in range(args.n_snps)]

    # Write out the SNP list
    snp_list_df = pd.DataFrame(snp_list)
    snp_list_df = snp_list_df.sort_values(['Chr', 'Position'])
    snp_list_df = snp_list_df[['SNP', 'Chr', 'Position']]
    out_path = os.path.join(args.out_dir, 'snp_list.txt')
    snp_list_df.to_csv(out_path, sep='\t', index=False, header=False)


if __name__ == '__main__':

    main()
