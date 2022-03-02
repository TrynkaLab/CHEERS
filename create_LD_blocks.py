'''
:File: create_LD_blocks.py
:Author: Lara Bossini-Castillo, Wellcome Sanger Institute, <lbc@sanger.ac.uk>
:Last updated: 20/1/2019

usage: python create_LD_blocks.py SNP_LIST OUTPUT_DIR LD_DIR


Copyright (C) 2019  Lara Bossini-Castillo

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
import os
import sys
from subprocess import Popen, PIPE


def parse_arguments():

    '''
    Parse and validate command line arguments.

    :return: A namespace containing the command line argument values.
    '''

    parser = argparse.ArgumentParser(description='Tagging SNP identification for CHEERS disease enrichment')
    parser.add_argument('snp_map', help='File containing SNPs to query for LD information')
    parser.add_argument('outdir', help='Directory for result output')
    parser.add_argument('tabix_dir', help='Directory containing LD tabix files')
    args = parser.parse_args()

    if not os.path.isfile(args.snp_map):
        print(f'{args.snp_map} is not a valid path to a file', file=sys.stderr)
        exit(2)
    
    if not os.path.isdir(args.outdir):
        print(f'{args.outdir} is not a valid path to a directory', file=sys.stderr)
        exit(2)
    
    if not os.path.isdir(args.tabix_dir):
        print(f'{args.tabix_dir} is not a valid path to a directory', file=sys.stderr)
        exit(2)

    return args


def snp_ld_tabix(snp, chrom, bp, tabix_dir, window, r2min):

    '''
    Retrieve LD information from the tabix file for a single SNP.

    :param snp: The SNP to query in the file.
    :param chrom: The chromosome the SNP is on.
    :param bp: The base-pair position of the SNP.
    :param tabix_dir: The directory containing the tabix files.
    :param window: The window around the SNP to query for tagging SNPs.
    :param r2min: The minimum R^2 value to include as a tagging SNP.
    :return: Returns a process that has run the tabix command to retrieve LD information.
    '''

    tabix_file = os.path.join(tabix_dir, f'{chrom}_GRCh38.EUR.ld.bgz')

    start = bp - window
    if start < 0:
        start = 0
    end = bp + window

    query = f"tabix {tabix_file} {chrom}:{start}-{end} | awk '$6 >= {r2min} {{print $0}}' | grep -w {snp}"
    process = Popen(query, shell=True, stdout=PIPE)

    return process


def read_snp_map_by_chr(file):

    '''
    Reads in file with SNP mappings.

    :param file: The path to the file.
    :return: A dictionary mapping chromosomes to SNPs and SNP positions.
    '''

    snp_map = dict()
    i = 0

    with open(file, 'r') as f_in:

        snp_number = 0
        header = next(f_in).rstrip().split('\t')
        chr_i = getCol('Chrom', header, file)
        snp_i = getCol('SNP', header, file)
        bp_i = getCol('BP', header, file)

        for line in f_in:

            i += 1
            line = line.rstrip().split('\t')
            chrom = line[chr_i]
            snp = line[snp_i]

            if snp == 'excluded': continue

            bp = int(line[bp_i])
            snp_map.setdefault(chrom, dict())
            snp_map[chrom].setdefault(snp, dict())
            snp_map[chrom][snp]['bp'] = bp
            snp_number += 1
            
    print(f'Read {snp_number} SNPs from {file}')

    return snp_map


def get_col(name, line, file):

    '''
    Retrieves the index of the column that contains the queried term.

    :param name: The term to query for.
    :param line: The first line of the table.
    :param file: The path to the file.
    :return: The index of the column.
    '''

    try:

        col = line.index(name)

    except ValueError:

        print(f'File: {file} missing {name} value in the header line. STOPPING', file=sys.stderr)
        exit(1)
    
    return col


def write_ld(snp_map, tabix_dir, window, r2min, dout):

    '''
    Takes SNP map and tabix files. Identifies SNPs in LD and writes into single
    files with LD info.

    :param snp_map: A map from chromosomes to SNPs and SNP positions.
    :param tabix_dir: The directory containing LD tabix files.
    :param window: The window around the SNP position to check for tagging SNPs.
    :param r2min: The minimum R^2 value to include as a tagging SNP.
    :param dout: The output directory.
    '''

    for chrom in snp_map:

        dout_snp = os.path.join(dout, chrom)

        if not os.path.exists(dout_snp):

            print(f'{dout_snp} does not exist, creating')
            os.mkdir(dout_snp)

        for snp in snp_map[chrom]:

            bp = snp_map[chrom][snp]['bp']

            # Get the LD information from tabix file
            proc = snp_ld_tabix(snp, chrom, bp, tabix_dir, window, r2min)

            # Name of the file to write results to
            f_name = f'results_ld_{snp}.txt'
            f_out = os.path.join(dout_snp, f_name)

            # Initiate for LD output
            ld_info = dict()

            snp_in_file = False
            for line in proc.stdout:

                fields = line.rstrip().split()
                snp1 = fields[2]
                bp1 = int(fields[1])
                snp2 = fields[4]
                bp2 = int(fields[3])
                r2 = float(fields[5])
                dprim = float(fields[6])

                if snp1 == snp:

                    ld_info.setdefault(snp2, dict())
                    ld_info[snp2]['bp2'] = bp2
                    ld_info[snp2]['r2'] = r2
                    ld_info[snp2]['d'] = dprim
                    snp_in_file = True

                elif snp2 == snp:

                    ld_info.setdefault(snp1, dict())
                    ld_info[snp1]['bp2'] = bp1
                    ld_info[snp1]['r2'] = r2
                    ld_info[snp1]['d'] = dprim
                    snp_in_file = True

            if not in_file:

                print(f'{snp} on {chrom} was not found in the tabix file.')

            else:

                print(f'Writing results to {f_name}.')

                with open(f_out, 'w') as f:
                    
                    line_format = '\t'.join(['{}'] * 7) + '\n'

                    for ld_snp, ld_snp_info in ld_info.items():

                        ld_bp = ld_snp_info['bp2']
                        ld_r2 = ld_snp_info['r2']
                        ld_d = ld_snp_info['d']

                        f.write(f'{chrom}\t{snp}\t{bp}\t{ld_snp}\t{ld_bp}\t{ld_r2}\t{ld_d}\n')
    
    print('Analysis finished.')


def main():

    # Parse command line arguments
    args = parse_arguments()

    # Window and dprime parameters
    window = 5e5
    r2min = 0.8

    snp_map = read_snp_map_by_chr(args.snp_map)
    write_ld(snp_map, args.tabix_dir, window, r2min, args.outdir)


if __name__ == "__main__":
    
    main()
