#!/usr/bin/env python3
"""
  This script takes the average CGI methylation bed file, plots the distribution of average CGI methylation levels,
  and save the plot as <averageCGI methylation bed basename>._distribution.png

usage: analyze_CGI_methylation.py -b /home/chase/code/cmatkhan/genomics/assignment5/WGBS_CGI_methylation.bed -o .
"""
import argparse
import pandas as pd
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns

def main(argv):
    # main method
    # Args: cmd line input. Use -h at cmd line to
    # display
    # Output:

    # read cmd line input arguments
    args = parseArgs(argv)
    bed_file = args.bed_path
    bed_file_basename = os.path.splitext(os.path.basename(bed_file))[0]
    output_dir = args.output_dir

    # check the bed file and output directory to make sure both exist
    if not os.path.isfile(bed_file):
        print('The .bed you entered cannot be found. Please check the path and re-launch this script.')
        sys.exit(1)
    if not os.path.isdir(output_dir):
        print('The output directory you entered does not exist. Please check the path and re-launch the script')
        sys.exit(1)

    avg_cgi_methyl = pd.read_table(bed_file, sep='\t', header = None)
    plotHist(avg_cgi_methyl,4, args.output_dir, bed_file_basename)
    print(avg_cgi_methyl)

def parseArgs(argv):
    # cmd line input
    parser = argparse.ArgumentParser(description="This script generates sbatch script and submits sbatch job.")
    parser.add_argument("-b", "--bed_path", required=True,
                        help="[Required] Path to WGS_CGI_methylation.bed")
    parser.add_argument('-o', "--output_dir", required=True,
                        help="[Required] Path to output directory")
    # remove script call from list (this is passed as list of cmd line input with script at position 0)
    args = parser.parse_args(argv[1:])
    return args

def plotHist(df,col, output_dir, bed_file_basename):
    #
    df[col].hist(grid=True, bins=20, color='#607c8e')
    plt.title('Avg CGI methylation')
    plt.xlabel('Avg CGI methylation')
    plt.ylabel('Count')
    plt.grid(axis='y', alpha=0.75)
    plt.savefig(os.path.join(output_dir, bed_file_basename + '_methylation_distribution.png'))
    plt.close()

# call main method
if __name__ == "__main__":
    main(sys.argv)