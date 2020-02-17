#!/usr/bin/env python3
"""
Lastly, we want to explore the CpG methylation profiles in “promoter-CGIs” versus“non-promoter-CGIs.”
Write a script called ​generate_promoters.py​ that takes a bed file of gene coordinates and creates a bed file of their promoters.
The columns will be
    1) chromosome   2) start    3) stop     4)gene name     5) strand.

usage: generate_promoters.py​​<bed of gene coordinates>

Run ​generate_promoters.py​ on ​refGene.bed​.
Save the output file as refGene_promoters.bed​.
In your README, justify how you defined promoter and paste your command for creating this file.
Generate a bed file of promoter-CGIs called ​promoter_CGI.bed​
and
a bed file of non-promoter-CGIs called ​non_promoter_CGI.bed​.
Promoter-CGIs are defined as CGIs that overlap with a promoter.
In your README, justify your overlapping criteria and paste your commands for creating these files.

Hint: use ​bedtools intersect​. Calculate the average CpG methylation level for each promoter-CGI and non-promoter-CGI.
Save these files as ​average_promoter_CGI_methylation.bed​ and average_non_promoter_CGI_methylation.bed​, respectively.
Paste your commands for creating these files in your README.Hint: see your commands for part 1.1.
"""
import argparse
import pandas as pd
import os
import sys


def main(argv):
    # main method

    # parse cmd line arguments
    args = parseArgs(argv)
    try:
        ref_gene_bed = args.bed_path
        output_dir = args.output_dir
    except AttributeError:
        print('One or both of the required inputs is missing. Please see script usage by entering generate_promoters.py -h')

    # read in refGene bed
    ref_gene_df = pd.read_csv(ref_gene_bed, sep='\t', header=None)

    # create a dataframe with chromosome, promoter start, promoter stop, gene name, strand
    promoter_bed = createPromoterBed(ref_gene_df)

    # write promoter_bed to file in .bed format
    writeBed(promoter_bed, output_dir, 'promoter_CGI.bed​')


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


def createPromoterBed(cds_df):
    # Creates a .bed with cols 1) chromosome   2) start    3) stop     4)gene name     5) strand from a bed_df with coding sequences in col2 and col3
    # Args: A df derived from a .bed with TSS in col2 and the rest of the info above
    # Returns: A .bed with the col structure in description above

    # create promoter_df initially with chromosome, gene name and strand columns
    promoter_df = cds_df[[0, 3, 5]]

    # I found this trick to calculate both columns at once here: https://stackoverflow.com/a/53666187/9708266
    def f(tss):
        # simple function that returns a 1x2 series. Used to create promoter start and promoter stop columns
        # by subtracting 1000 and 40 from the cds start value
        return pd.Series([tss - 1000, tss - 30])

    # add promoter start and stop site columns using function f above
    promoter_df[['promoter_start', 'promoter_stop']] = cds_df[1].apply(f)

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
    df = df.astype({'promoter_start': 'int64', 'promoter_stop': 'int64'})
    df = df.sort_values(by=['promoter_start'])

    output_path = os.path.join(output_dir, filename)

    df.to_csv(output_path, index=False, header=False, sep='\t')


# call main method
if __name__ == "__main__":
    main(sys.argv)
