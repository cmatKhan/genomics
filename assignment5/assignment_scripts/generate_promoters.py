#!/usr/bin/env python3
"""
Lastly, we want to explore the CpG methylation profiles in “promoter-CGIs” versus“non-promoter-CGIs.”
Write a script called ​generate_promoters.py​ that takes a bed file of gene coordinates and creates a bed file of their promoters.
The columns will be
    1) chromosome   2) start    3) stop     4)gene name     5) strand.

usage: generate_promoters.py​​ -b refGene.bed -o . -bn refGene_promoter.bed
"""
import argparse
import pandas as pd
import os
import sys


def main(argv):
    # main method
    # Args: cmd line input
    # Output: reGene_promoter.bed in output_dir

    # parse cmd line arguments
    args = parseArgs(argv)
    ref_gene_bed = args.bed_path
    output_dir = args.output_dir
    output_bed_name = args.output_bed_name

    # verify cmd line input
    if not os.path.isfile(ref_gene_bed):
        print('refGene.bed cannot be found. Please check the path and try again')
    if not os.path.isdir(output_dir):
        print('The output directory does not exist. Please check the path and try again')

    # read in refGene bed
    ref_gene_df = pd.read_csv(ref_gene_bed, sep='\t', header=None)

    # create a dataframe with chromosome, promoter start, promoter stop, gene name, strand
    promoter_bed = createPromoterBed(ref_gene_df)

    # write promoter_bed to file in .bed format
    writeBed(promoter_bed, output_dir, output_bed_name)


def parseArgs(argv):
    # cmd line input
    parser = argparse.ArgumentParser(description="This script generates sbatch script and submits sbatch job.")
    parser.add_argument("-b", "--bed_path", required=True,
                        help="[Required] Path to refGene.bed")
    parser.add_argument('-o', "--output_dir", required=True,
                        help="[Required] Path to output directory")
    parser.add_argument('-bn', "--output_bed_name", required=True,
                        help="[Required] The name of the .bed")

    # remove script call from list (this is passed as list of cmd line input with script at position 0)
    args = parser.parse_args(argv[1:])
    return args


def createPromoterBed(cds_df):
    # Creates a .bed with cols 1) chromosome   2) start    3) stop     4)gene name     5) strand from a bed_df with coding sequences in col2 and col3
    # Args: A df derived from a .bed with TSS in col2 and the rest of the info above
    # Returns: A .bed with the col structure in description above

    # create promoter_df initially with chromosome, gene name and strand columns
    promoter_df = cds_df[[0, 3, 5]]

    # I found this trick to calculate both columns at once here: https://stackoverflow.com/a/53666187/9708266
    def f_plus(gene_start):
        # returns a 1x2 series. Used to create promoter start and stop columns for plus strand
        return pd.Series([gene_start - 1000, gene_start - 30])

    def f_minus(gene_start):
        # simple function that returns a 1x2 series. Used to create promoter start and stop columns for minus strand
        return pd.Series([gene_start + 1000, gene_start + 30])

    # masks used to apply above functions to correct rows in cds_df in calculations below
    mask_plus = promoter_df[5] == '+'
    mask_minus = promoter_df[5] == '-'

    # calculate promoter start/stop
    pos_strand, neg_strand = pd.DataFrame(), pd.DataFrame() # create two dataframes
    # calc promoter start/stop on + strand
    pos_strand[['promoter_start', 'promoter_stop']] = cds_df[1].loc[mask_plus].apply(f_plus)
    # calc promoter start/stop on - strand
    neg_strand[['promoter_start', 'promoter_stop']] = cds_df[2].loc[mask_minus].apply(f_minus)
    # combine plus and minus into single df
    promoter_start_stop = pos_strand.combine_first(neg_strand)
    # concat onto df with cols chr, gene_name_strand
    promoter_df = pd.concat([promoter_df, promoter_start_stop], axis=1)


    # re-order columns
    cols = [0, 'promoter_start', 'promoter_stop', 3, 5]
    promoter_df = promoter_df[cols]

    return promoter_df


def writeBed(df, output_dir, filename):
    # write .bed style table in output directory
    # Args: a dataframe with the following columns and an output directory
    # Returns: none
    # Output: prints a file called _CpG_methylation.bed as bed (tsv) to specified directory

    # cast position columns to int and sort
    df = df.astype({0: 'str', 'promoter_start': 'int64', 'promoter_stop': 'int64', 3: 'str', 5:'str'})
    df = df.sort_values(by=['promoter_start'])

    output_path = os.path.join(output_dir, filename)

    df.to_csv(output_path, index=False, header=False, sep='\t')


# call main method
if __name__ == "__main__":
    main(sys.argv)
