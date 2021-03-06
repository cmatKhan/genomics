#!/usr/bin/env python3
"""
that takes a WGBS bed file and:1.Calculates the methylation level for ​each CpG​ using the formula:

CpG_methylation_level = C_base_calls / C_base_calls + T_base_calls

2.Writes a bed-like file of CpG methylation. The columns will be
                1) chromosome
                2) start
                3) stop
                4) methylation level.
 Do not output CpGs that have 0X coverage. Save thefile as ​<WGBS bed basename>​_CpG_methylation.bed​.

 3.Plots the distribution, e.g., with a histogram, of CpG methylation levels as ​<WGBS bedbasename>​_methylation_distribution.png​.
 4.Plots the distribution of read coverage for ​all​ CpGs​ for coverages between 0X and 100Xas ​<WGBS bed basename>​_CpG_coverage_distribution.png​.
 5.Calculates and print the fraction of CpGs that have 0X coverage

 usage: analyze_WGBS_methylation.py -b <path_to>/BGM_WGBS.bed -o <path_to_output_directory>
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
    # Args: cmd line input. Use -h at cmd line to display
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

    # this function both returns the bed + methyl level as a df and writes the bed + methyl to output_dir
    print('\n...Calculating methylation level and coverage...')
    methyl_coverage_df = createBGSbed(bed_file, output_dir)
    print("Methylation calculations complete. See {}_CpG_methylation.bed​ in output directory for methylation level.".format(bed_file_basename))

    # create csv to plot methylation in R
    #writeCSV(methyl_coverage_df, output_dir, 'methylation_coverage.csv​')

    # plot coverage and methylation

    plotHist(methyl_coverage_df, output_dir, bed_file_basename)

    # plot with R scripts
    print("\n...Plotting distributions of methylation levels and coverage. See the output directory...")

    # print % 0x coverage
    print("\nThe fraction of 0x CpG sites is: {}".format(len(methyl_coverage_df[methyl_coverage_df['coverage'] == 0]) / len(methyl_coverage_df)))


def parseArgs(argv):
    # cmd line input
    parser = argparse.ArgumentParser(description="This script generates sbatch script and submits sbatch job.")
    parser.add_argument("-b", "--bed_path", required=True,
                        help="[Required] Path to BGM_WGS.bed")
    parser.add_argument('-o', "--output_dir", required=True,
                        help="[Required] Path to output directory")
    # remove script call from list (this is passed as list of cmd line input with script at position 0)
    args = parser.parse_args(argv[1:])
    return args

def createBGSbed(bed_file, output_dir):
    # remove records with zero coverage
    # Args: a bed file

    bed_file_basename = os.path.splitext(os.path.basename(bed_file))[0]

    # read in the .bed as a dataframe and filter out rows with 0s
    input_bed_df = pd.read_table(bed_file, sep='\t', header = None)
    # create new column with methyl level
    input_bed_df['methyl_level'] = input_bed_df.apply(lambda row: cpgMethylLevel(row[3], row[4]), axis=1)
    input_bed_df['coverage'] = input_bed_df.apply(lambda row: row[3] + row[4], axis = 1)
    # filter for 0 coverage and write chrm name, start, stop and methylation le of methyl_bed
    methyl_bed = input_bed_df.where(input_bed_df['coverage'] != 0).dropna()
    writeBed(methyl_bed[[0, 1, 2, 'methyl_level']], output_dir, bed_file_basename, '_CpG_methylation.bed​')

    # return df for graphing purposes
    return input_bed_df[['methyl_level', 'coverage']]

def cpgMethylLevel(num_c_base_calls, num_t_base_calls):
    # calculate CpG methylation
    # Args: number of C base calls and number of T base calls
    # Returns: CpG methylation level calculated by num_c_base_calls / (num_c_base_calls + num_t_base_calls)

    # to avoid division by zero
    if num_t_base_calls + num_t_base_calls == 0:
        return 0
    else:
        return num_c_base_calls / float(num_c_base_calls + num_t_base_calls)

def writeBed(df, output_dir, bed_base_name, filename):
    # write .bed style table in output directory
    # Args: a dataframe with the following columns, an output directory, and a filename
    # Returns: none
    # Output: prints a file called _CpG_methylation.bed as bed (tsv) to specified directory

    # cast position columns to int and sort
    df = df.astype({1: 'int64', 2:'int64'})
    df =df.sort_values(by=[1])

    output_path = os.path.join(output_dir, bed_base_name + filename)

    df.to_csv(output_path, index = False, header = False, sep='\t')

def writeCSV(df, output_dir, filename):
    # write a csv of methylation_level and coverage for plotting in R
    # Args: dataframe with methylation level and coverage, output directory and filename
    # Output: prints methylation_coverage.csv to output_dir

    output_path = os.path.join(output_dir, filename)
    df.to_csv(output_path, index=False)

def plotHist(df, output_dir, bed_file_basename):
    df['methyl_level'].hist(grid=True, bins=20, color='#607c8e')
    plt.title('Methylation Level')
    plt.xlabel('CpG level')
    plt.ylabel('Count')
    plt.grid(axis='y', alpha=0.75)
    plt.savefig(os.path.join(output_dir, bed_file_basename + '_methylation_distribution.png'))
    plt.close()

    sns.distplot(df['coverage'], hist=True, kde=True,
                 bins=10000, color='#607c8e',
                 hist_kws={'edgecolor': 'black'})
    # Add labels
    plt.xlim(0, 100)
    plt.title('Coverage')
    plt.xlabel('Coverage')
    plt.ylabel('Frequency')
    plt.grid(axis='y', alpha=0.75)
    plt.savefig(os.path.join(output_dir, bed_file_basename + '_CpG_coverage_distribution.png'))
    #df['coverage'].hist(grid=True, bins=10000, color='#607c8e', )
    #plt.title('Coverage')
    #plt.xlabel('Coverage')
    #plt.ylabel('Frequency')
    #plt.grid(axis='y', alpha=0.75)
    #plt.savefig(os.path.join(output_dir, bed_file_basename + '_CpG_coverage_distribution.png'))

# call main method
if __name__ == "__main__":
    main(sys.argv)
