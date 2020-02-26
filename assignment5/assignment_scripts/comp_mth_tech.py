#!/usr/bin/env python3
"""
This script Creates the following scatter plots

    a.MeDIP-seq RPKM vs. MRE-seq RPKM.
    Save the plot as ​<MeDIP-seq RPKM bedbasename>​_vs_​<MRE-seq RPKM bed basename>​.png​

    b.MeDIP-seq RPKM vs. WGBS average DNA methylation level.
    Save the plot as<MeDIP-seq RPKM bed basename>​_vs_​<WGBS methylation level bedbasename>​.png​

    c.MRE-seq RPKM vs. WGBS average DNA methylation level.
    Save the plot as<MRE-seq RPKM bed basename>​_vs_​<WGBS methylation level bedbasename>​.png​

2.Calculate the correlation for each comparison. Print the correlation to stdout.

usage: compare_methylome_technologies.py -md MeDIP_RPKM_methyl_levels.bed -mr MRE_RPKM_methyl_levels.bed -bs BGM_RPKM_methyl_levels.bed -o .
"""
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
import os
import sys

def main(argv):
    # main method
    # Ags: cmd line
    # Output: comparisons of correlation between methylation techniques of measuring methylation

    # parse cmd line arguments
    args = parseArgs(argv)
    md = args.medip
    mre = args.mre
    bs = args.wgbs
    output_dir = args.output_dir
    if not os.path.isfile(md) or not os.path.isfile(mre) or not os.path.isfile(bs):
        print('One of the inputted beds cannot be found. Please check your cmd line input and try again.')
        sys.exit(1)
    if not os.path.isdir(output_dir):
        print('ckeck the path to the output directory and try again.')
        sys.exit(1)

    # create dataframes
    medip_df = pd.read_table(md, sep = '\t', header = None)
    mre_df = pd.read_table(mre, sep = '\t', header = None)
    bs_df = pd.read_table(bs, sep = '\t', header = None)

    plotCorrelations(medip_df[[4]], mre_df[[4]], 'MeDIP_RPKM', 'MRE_RPKM', output_dir, 'MeDIP_Seq_RPKM_vs_MRE_Seq_RPKM')
    plotCorrelations(medip_df[[4]], bs_df[[4]], 'MeDIP_RPKM', 'WGBS_RPKM', output_dir, 'MeDIP_Seq_RPKM_vs_WGBS_Seq_RPKM')
    plotCorrelations(mre_df[[4]], bs_df[[4]], 'MRE_RPKM', 'WGBS_RPKM', output_dir, 'MRE_Seq_RPKM_vs_WGBS_Seq_RPKM')

    # remove outlier
    medip_mod_df = medip_df[(medip_df[1] != 9825442) & (medip_df[2] != 9826296)]
    mre_mod_df = mre_df[(mre_df[1] != 9825442) & (mre_df[2] != 9826296)]
    bs_mod_df = bs_df[(bs_df[1] != 9825442) & (bs_df[2] != 9826296)]

    # replot
    plotCorrelations(medip_mod_df[[4]], mre_mod_df[[4]], 'MeDIP_RPKM', 'MRE_RPKM', output_dir, 'No_Outlier_MeDIP_Seq_RPKM_vs_MRE_Seq_RPKM_')
    plotCorrelations(medip_mod_df[[4]], bs_mod_df[[4]], 'MeDIP_RPKM', 'WGBS_RPKM', output_dir, 'No_Outlier_MeDIP_Seq_RPKM_vs_WGBS_Seq_RPKM')
    plotCorrelations(mre_mod_df[[4]], bs_mod_df[[4]], 'MRE_RPKM', 'WGBS_RPKM', output_dir, 'No_Outlier_MRE_Seq_RPKM_vs_WGBS_Seq_RPKM')

def parseArgs(argv):
    # cmd line input
    parser = argparse.ArgumentParser(description="This script generates sbatch script and submits sbatch job.")
    parser.add_argument("-md", "--medip", required=True,
                        help="[Required] Path to MeDIP_RPKM_methyl_levels.bed")
    parser.add_argument("-mr", "--mre", required=True,
                        help="[Required] Path to MRE_RPKM_methyl_levels.bed")
    parser.add_argument("-bs", "--wgbs", required=True,
                        help="[Required] Path to WGBS_RPKM_methyl_levels.bed")
    parser.add_argument('-o', "--output_dir", required=True,
                        help="[Required] Path to output directory")

    # remove script call from list (this is passed as list of cmd line input with script at position 0)
    args = parser.parse_args(argv[1:])
    return args

def plotCorrelations(col1_rpkm, col2_rpkm, x_lab, y_lab, output_dir, main):
    # comparing sequencing technologies
    # Args: the RPKM cols of related dataframes as named matricies, output directory and main heading in filename format (used as both title and filename)
    # Output: .png with same name as main, print out of correlation to stdout

    # create output path
    output_path = os.path.join(output_dir, main + '.png')

    # create two column df and plot as a scatter plot
    scatter_df = pd.concat([col1_rpkm, col2_rpkm], axis=1)
    scatter_df.columns = [x_lab, y_lab]
    scatter_df.plot.scatter(x=x_lab, y=y_lab)
    plt.savefig(output_path)
    plt.close()

    print('\nThe pearson R coefficient for {}\nis {}\n'. format(main, scipy.stats.pearsonr(scatter_df[x_lab], scatter_df[y_lab])))


if __name__ == "__main__":
    main(sys.argv)
