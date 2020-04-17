#! /usr/bin/env python3

"""
    analyze_MPRA.py takes in both pDNA and cDNA count files together with filtered_variant_to_barcode.txt​.
    This script will:
    1.Calculate a normalized expression(NE) level for each barcode. For each oligo-GFP construct represented by a barcode,
    the normalized expression level is calculated by log transformed ratio of cDNA to pDNA barcode count:
            NE = log2( #ofoccurrenceincDNA / #ofoccurrenceinpDNA2)
    2. Assign NE calculated by barcode back to variant by using variant-barcode assignment file.
    3.Calculate log transformed fold change(FC) in normalized expression level for each variant pair by using the formula:
            log2FC = (avg alt allele) - (avg ref allele)
    4.Determine significance of the fold change by using Mann-Whitney U test
    5.Output to file ​variant_fold_change.txt​

    usage: analyze_mpra.py -c cDNA_count.tsv -p pDNA_count.tsv -f filtered_variant_to_barcode.tsv
"""
import sys
import argparse
import logging
import pandas as pd
from scipy.stats import mannwhitneyu

def main(argv):
    """ main method
    :param argv: cmd line arguments
    """
    # create logger
    global logger
    logger = createLogger('/home/chase/code/cmatkhan/genomics/assignment11/log/assignment11.log', __name__)

    # parse cmd line arguments
    args = parseArgs(argv)
    c_dna_count_path = args.cDNA_counts
    p_dna_count_path = args.pDNA_counts
    filtered_barcode_variant_path = args.filtered_barcode_to_variant
    eqtl_df = pd.read_table(args.eQTL, sep='\t')

    # create variant_barcode_dict with structure
    ref_variant_barcode_df, alt_variant_barcode_df = createVariantBarcodeDataframe(filtered_barcode_variant_path)

    # merge barcode count files into df with columns barcode, c_dna_count, p_dna_count
    ref_variant_barcode_count_df, alt_variant_barcode_count_df = createCombinedCountDataframe(ref_variant_barcode_df,
                                                                                              alt_variant_barcode_df,
                                                                                              c_dna_count_path,
                                                                                              p_dna_count_path)
    # quantify/summarize
    variant_fc_df = analyzeCounts(ref_variant_barcode_count_df, alt_variant_barcode_count_df)
    # write output
    writeOut(variant_fc_df, 'variant_fold_change.tsv')

    # read in variant_to_eQTL_results and print out info for question 4
    df = variant_fc_df.merge(eqtl_df, how = 'inner', on='Variant_ID')
    df = df.groupby('Significant_eQTL_gene').agg('count')[['Variant_ID']].rename(columns={'Variant_ID': 'Count_by_eQTL'})
    print(df)

def parseArgs(argv):
    parser = argparse.ArgumentParser(
        description="description here")
    parser.add_argument("-c", "--cDNA_counts", required=True,
                        help='path to cDNA_count.tsv')
    parser.add_argument("-p", "--pDNA_counts", required=True,
                        help='path to pDNA_count.tsv')
    parser.add_argument("-f", "--filtered_barcode_to_variant", required=True,
                        help='path to filtered_barcode_to_variant.tsv')
    parser.add_argument("-e", "--eQTL", required=True,
                        help='path to variant_eQTL_results.txt')
    args = parser.parse_args(argv[1:])
    return args

def createBarcodeDict(variant_to_barcode_path):
    """
    read in filtered variant_to_barcode, return dict
    :param variant_to_barcode_path: path to variant_to_barcode.txt
    :return: dict of structure {variant_id_1: [[ref_barcodes...], [alt_barcodes]], variant_id_2: ...}
    """

    barcode_dict = {}

    with open(variant_to_barcode_path, 'r') as variant_to_barcode_file:
        line = variant_to_barcode_file.readline()  # skip the header

        if not line == 'Variant_ID\tREF_barcode\tALT_barcode\n':
            raise Exception('the variant_to_barcode dict does not have the proper header of:\n'
                            'Variant_ID\tREF_barcode\tALT_barcode\n'
                            'check the formatting of the file and try again')

        line = variant_to_barcode_file.readline().strip().split('\t')
        while line and not line == ['']:
            # extract header
            variant_id = line[0]
            # put list of ref barcodes and alt barcodes in list structure [[list_ref_barcodes], [list_alt_barcodes]]
            barcode_list = [line[1].split(':'), line[2].split(':')]
            # add entry to dict
            barcode_dict.setdefault(variant_id, barcode_list)
            # next line
            line = variant_to_barcode_file.readline().strip().split('\t')

    return barcode_dict

def createVariantBarcodeDataframe(variant_barcode_path):
    """
    create dataframe of variant_id and barcode split by ref and alt
    :param variant_barcode_path:
    :return: two dataframes, ref_df and alt_df with index column variant_id and one column Barcodes
    """
    barcode_dict = createBarcodeDict(variant_barcode_path)
    # data frame in with columns variant_id, ref_barcodes (a list), alt_barcodes (a list)
    barcode_df = pd.DataFrame.from_dict(barcode_dict, orient='index')
    # take columns variant_id and ref_barcodes, and then expand the list into individual rows
    ref_df = barcode_df[[0]].explode(0).reset_index().rename(columns={'index': 'Variant_ID', 0: 'Barcodes'})
    # (continued from prev. comment) st the variant_id (repeated) is associated with each ref_barcode in the list
    alt_df = barcode_df[[1]].explode(1).reset_index().rename(columns={'index': 'Variant_ID', 1: 'Barcodes'})

    return ref_df, alt_df

def createCombinedCountDataframe(ref_variant_barcode_df, alt_variant_barcode_df, c_dna_count_path, p_dna_count_path):
    """
    create dataframes with variant, barcode and count data for ref and alt
    :param ref_variant_barcode_df: see createVariantBarcodeDataframe()
    :param alt_variant_barcode_df: see createVariantBarcodeDataframe()
    :param c_dna_count_path: path to .tsv in two columns with headers Barcodes\tCount
    :param p_dna_count_path: path to .tsv in two columns with headers Barcodes\tCount
    :return: two dataframes with columns Variant_ID, Barcodes, cDNA_Count, pDNA_Count
    """
    # read in c_dna and p_dna counts as dataframes
    c_dna_count_df = pd.read_csv(c_dna_count_path, sep='\t')
    p_dna_count_df = pd.read_csv(p_dna_count_path, sep='\t')

    # merge c_dna and p_dna counts with ref_variant_barcode_df
    ref_variant_barcode_count_df = ref_variant_barcode_df.merge(c_dna_count_df, how = 'left', on='Barcodes')
    ref_variant_barcode_count_df.rename(columns={'Count': 'cDNA_Count'}, inplace=True)
    ref_variant_barcode_count_df = ref_variant_barcode_count_df.merge(p_dna_count_df, how = 'left', on='Barcodes')
    ref_variant_barcode_count_df.rename(columns={'Count': 'pDNA_Count'}, inplace=True)

    # do the same with the alt_variant_barcode_df
    alt_variant_barcode_count_df = alt_variant_barcode_df.merge(c_dna_count_df, how = 'left', on='Barcodes')
    alt_variant_barcode_count_df.rename(columns={'Count': 'cDNA_Count'}, inplace=True)
    alt_variant_barcode_count_df = alt_variant_barcode_count_df.merge(p_dna_count_df, how = 'left', on='Barcodes')
    alt_variant_barcode_count_df.rename(columns={'Count': 'pDNA_Count'}, inplace=True)

    return ref_variant_barcode_count_df, alt_variant_barcode_count_df

def analyzeCounts(ref_count_df, alt_count_df):
    """
    calculate fold change and p-value
    :param ref_count_df: see createCombinedCountDataframe()
    :param alt_count_df: see createCombinedCountDataframe()
    :return:
    """
    # For reference variants, calculate normalized expression by barcode
    ref_count_df['Norm_Expr_ref'] = ref_count_df['cDNA_Count'] / ref_count_df['pDNA_Count']
    # group by variant and take average of normalized expression of each barcode
    ref_norm_expr_df = ref_count_df[['Variant_ID', 'Norm_Expr_ref']]

    # Do the same for alternate variants
    alt_count_df['Norm_Expr_alt'] = alt_count_df['cDNA_Count'] / alt_count_df['pDNA_Count']
    # group by variant and take average of normalized expression of each barcode
    alt_norm_expr_df = alt_count_df[['Variant_ID', 'Norm_Expr_alt']]

    # create list to hold pvalues
    pvalue_list = []
    # loop over list of variant_ids
    for i in ref_norm_expr_df['Variant_ID'].unique():
        # extract rows from ref and alt given a specific variant_id and calculate the mann whitney p value with those arrays. append to list
        pvalue_list.append(mannwhitneyu(ref_norm_expr_df.loc[ref_norm_expr_df['Variant_ID'] == i]['Norm_Expr_ref'],
                           alt_norm_expr_df.loc[alt_norm_expr_df['Variant_ID'] == i]['Norm_Expr_alt'])[1])

    # calculate mean count in ref and alt
    ref_norm_expr_df = ref_norm_expr_df.groupby(['Variant_ID']).mean()
    ref_norm_expr_df.reset_index(inplace=True)
    # see prev. comment
    alt_norm_expr_df = alt_norm_expr_df.groupby(['Variant_ID']).mean()
    alt_norm_expr_df.reset_index(inplace=True)

    # merge ref_norm_expr_df and alt_norm_expr_df after grouping by variant_id and taking mean
    variant_expression_df = ref_norm_expr_df.merge(alt_norm_expr_df, how = 'left', on='Variant_ID')
    # calculate log2 fold change
    variant_expression_df['log2_fc'] = variant_expression_df['Norm_Expr_alt'] - variant_expression_df['Norm_Expr_ref']
    # add p-value column
    variant_expression_df['p_value'] = pvalue_list

    return variant_expression_df[['Variant_ID', 'log2_fc', 'p_value']]

def writeOut(variant_expression_df, output_filename):
    """

    :param variant_expression_df:
    :param output_filename:
    :return:
    """
    print(variant_expression_df)
    print('\n')
    variant_expression_df.to_csv(output_filename, index=False, sep='\t')

def createLogger(filename, logger_name, logging_conf = None):
    """
    create logger in filemode append and format name-levelname-message with package/module __name__ (best practice from logger tutorial)
    :param filename: name of the file in which to log.
    :param logger_name: __name__ is recommended as a best practice in logger.config eg you can call this like so: createLogger(<your_filename>, __name__)
                     (__name__ is a special variable in python)
    :param logging_conf: path to logging configuration file
    :returns: an instance of the configured logger
    """
    # a config file is passed, load it
    if logging_conf:
        logging.config.fileConfig(logging_conf) # should include at least what is below
    # if it is not, configure as follows
    else:
        # create log for the year-month-day
        logging.basicConfig(
            filename='%s' % filename,
            filemode='w',
            format='%(name)s-%(levelname)s-%(asctime)s-%(message)s',
            datefmt='%I:%M:%S %p', # set 'datefmt' to hour-minute-second AM/PM
            level='DEBUG'
        )
    # return an instance of the configured logger
    return logging.getLogger(logger_name)


if __name__ == "__main__":
    main(sys.argv)