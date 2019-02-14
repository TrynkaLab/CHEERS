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

import os
import sys
from subprocess import Popen, PIPE


def snpLdTabix(snp, chrom, bp, tabixDir, window, r2min):
    """
    Retrive LD info from the tabix file for a single SNP.
    """
    file = os.path.join(tabixDir, '{}_GRCh38.EUR.ld.bgz')
    tabixFile = file.format(chrom)
    print 'tabixFile'
    st = bp - window
    if st < 0:
        st = 0
    en = bp + window

    query = "tabix {tabixFile} {chrom}:{st}-{en} | awk '$6 >= {r2min} {{print $0}}' | grep -w {snp}".format(**locals())
    proc = Popen(query, shell=True, stdout=PIPE)

    return proc


def readSnpMapByChr(file):
    """
    Reads in file with SNP mappings.
    Returns dict, key:chrom, value:snp, bp.
    """
    snpMap = {}
    i = 0
    with open(file, 'r') as f:
        snpNu = 0
        header = next(f).rstrip().split("\t")
        chr_i = getCol('Chrom', header, file)
        snp_i = getCol('SNP', header, file)
        bp_i = getCol('BP', header, file)

        for line in f:
            i += 1
            line = line.rstrip().split("\t")
            # checkMappingsFormat(file,line,i)
            chr = line[chr_i]
            snp = line[snp_i]
            if snp == 'excluded': continue
            bp = int(line[bp_i])
            snpMap.setdefault(chr, {})
            snpMap[chr].setdefault(snp, {})
            snpMap[chr][snp]['bp'] = bp
            snpNu += 1
    print 'Read {snpNu} SNPs from {file}'.format(**locals())
    return snpMap


def getCol(name, line, file):
    """
    Retrives the index of column that contains queried term
    """
    try:
        col = line.index(name)
    except ValueError:
        sys.exit("File: {file} missing {name} value in the header line. STOPPING".format(**locals()))
    return col


def write_ld(snp_map, tabixDir, window, r2min, dout):
    """ Takes SNP map and tabix files, identifies SNPs in LD
        and writes into single files with LD info """

    line = '\t'.join(['{}'] * 7) + '\n'
    for chrom in snp_map:
        dout_snp = os.path.join(dout, chrom)
        if not os.path.exists(dout_snp):
            print "{} does not exist, creating".format(dout_snp)
            os.mkdir(dout_snp)

        for snp in snp_map[chrom]:
            bp = snp_map[chrom][snp]['bp']

            # get the LD info from tabix file
            proc = snpLdTabix(snp, chrom, bp, tabixDir, window, r2min)

            # name of the file to write results to
            fname = "results_ld_" + snp + ".txt"
            fout = os.path.join(dout_snp, fname)

            # initiate for LD output
            ldInfo = {}

            infile = 0
            for line in proc.stdout:
                fields = line.rstrip().split()
                snp1 = fields[2]
                bp1 = int(fields[1])
                snp2 = fields[4]
                bp2 = int(fields[3])
                r2 = float(fields[5])
                dprim = float(fields[6])

                if snp1 == snp:
                    ldInfo.setdefault(snp2, {})
                    ldInfo[snp2]['bp2'] = bp2
                    ldInfo[snp2]['r2'] = r2
                    ldInfo[snp2]['d'] = dprim
                    infile = 1
                elif snp2 == snp:
                    ldInfo.setdefault(snp1, {})
                    ldInfo[snp1]['bp2'] = bp1
                    ldInfo[snp1]['r2'] = r2
                    ldInfo[snp1]['d'] = dprim
                    infile = 1

            if not infile:
                print "{snp} on {chrom} was not found in the tabix file." \
                    .format(**locals())
            else:
                print "writing results to {}".format(fname)
                with open(fout, 'w') as f:
                    # format string
                    line_format = '\t'.join(['{}'] * 7) + '\n'

                    for ldsnp in ldInfo:
                        ldbp = ldInfo[ldsnp]['bp2']
                        ldr2 = ldInfo[ldsnp]['r2']
                        ldd = ldInfo[ldsnp]['d']
                        f.write(line_format.format(chrom, snp, bp, ldsnp, ldbp, ldr2, ldd))
    end_msg = "Analysis finished"
    print end_msg


if __name__ == "__main__":
    # file with snp map
    f_snp = sys.argv[1]
    # output directory
    d_out = sys.argv[2]
    # tabix directory
    d_tabix = sys.argv[3]

    # window and dprime parameters
    window = 5e5
    r2min = 0.8

    snp_map = readSnpMapByChr(f_snp)
    write_ld(snp_map, d_tabix, window, r2min, d_out)
