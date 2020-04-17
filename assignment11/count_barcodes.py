#! /usr/bin/env python3

"""
    count_barcodes.py filters out reads that do not perfectly match barcode and counts the total number of perfectly
    matching reads for each barcode
    usage: count_barcodes.py -b filtered_variant_to_barcode.tsv -f <cDNA.fq or pDNA.fq>
"""
import sys
import argparse
import numpy as np
import pandas as pd
import logging
import os
import subprocess


def main(argv):
    """ main method
    :param argv: cmd line arguments
    """
    # create logger
    global logger
    logger = createLogger('/home/chase/code/cmatkhan/genomics/assignment11/log/assignment11.log', __name__)
    # parse cmd line arguments
    args = parseArgs(argv)
    filtered_variant_to_barcode_path = args.filtered_variant_to_barcode
    fastq_file = args.fastq_file
    # create output path
    output_path = os.path.basename(fastq_file).replace('.fq', '_count.tsv')

    # create list of valid barcodes
    valid_barcode_list = createBarcodeList(filtered_variant_to_barcode_path)
    # count barcodes, print result, return dict to write
    barcode_count_results, unfiltered_read_count = countBarcodes(valid_barcode_list, fastq_file)
    # print result and write results
    reportResults(barcode_count_results, output_path, unfiltered_read_count)


def parseArgs(argv):
    parser = argparse.ArgumentParser(
        description="description here")
    parser.add_argument("-b", "--filtered_variant_to_barcode", required=True,
                        help='path to filtered_variant_to_barcode.txt')
    parser.add_argument("-f", "--fastq_file", required=True,
                        help='path to fastq file')
    args = parser.parse_args(argv[1:])
    return args


def createBarcodeList(variant_to_barcode_path):
    """
    read in filtered variant_to_barcode, return dict
    :param variant_to_barcode_path: path to variant_to_barcode.txt
    :return: dict of structure {variant_id_1: [[ref_barcodes...], [alt_barcodes]], variant_id_2: ...}
    """
    barcode_list = []

    with open(variant_to_barcode_path, 'r') as variant_to_barcode_file:
        line = variant_to_barcode_file.readline()  # skip the header

        if not line == 'Variant_ID\tREF_barcode\tALT_barcode\n':
            raise Exception('the variant_to_barcode dict does not have the proper header of:\n'
                            'Variant_ID\tREF_barcode\tALT_barcode\n'
                            'check the formatting of the file and try again')

        line = variant_to_barcode_file.readline().strip().split('\t')
        while line and not line == ['']:
            barcode_list.append(line[1].split(':'))
            barcode_list.append(line[2].split(':'))
            line = variant_to_barcode_file.readline().strip().split('\t')
    #logger.debug(barcode_list)
    # flatten list
    flat_barcode_list = [item for sublist in barcode_list for item in
                         sublist]  # credit https://stackoverflow.com/a/952952/9708266
    # drop any repeat barcodes, just in case
    flat_barcode_list = np.array(flat_barcode_list)
    flat_barcode_list = list(np.unique(flat_barcode_list))
    return flat_barcode_list


def countBarcodes(variant_barcode_list, fastq_file_path):
    """
    count barcodes in fastq_file
    :param variant_barcode_list. see createBarcodeList()
    :param fastq_file_path: path to fastq file
    :return: print out results, return dict in structure {}
    """
    barcode_count_dict = {}
    with open(fastq_file_path, 'r') as fastq_file:
        line = fastq_file.readline().strip()
        fastq_file_read_count = 0
        while line:
            if line.startswith('@SOLEXA6'):
                fastq_file_read_count = fastq_file_read_count + 1
                barcode = fastq_file.readline().strip()[14:23]
                logger.debug(barcode)
                if not len(barcode) == 9:
                    raise Exception('The extracted barcode is not 9 bases long. Check code')
                if barcode in variant_barcode_list:
                    try:
                        barcode_count_dict[barcode] += 1
                    except KeyError:
                        barcode_count_dict.setdefault(barcode, 1)
            line = fastq_file.readline().strip()
    return barcode_count_dict, fastq_file_read_count


def reportResults(barcode_count_dict, output_path, unfiltered_read_count):
    """
    report results from barcode_count_dict
    :param barcode_count_dict: see countBarcodes()
    :return: None -- print to stdout, write to file
    """
    barcode_count_df = pd.DataFrame.from_dict(barcode_count_dict, orient='index')
    barcode_count_df.reset_index(inplace=True)
    barcode_count_df.rename(columns={'index': 'Barcodes', 0: 'Count'}, inplace=True)
    barcode_count_df.to_csv(output_path, sep='\t', index=False)
    cmd = 'column -s\'\\t\' %s' % output_path
    os.system(cmd)
    print('\nPrior to filtering, there were %i reads\n'
          'The number of reads left after filtering is: %i\n' % (unfiltered_read_count,
                                                                 barcode_count_df['Count'].sum()))


def createLogger(filename, logger_name, logging_conf=None):
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
        logging.config.fileConfig(logging_conf)  # should include at least what is below
    # if it is not, configure as follows
    else:
        # create log for the year-month-day
        logging.basicConfig(
            filename='%s' % filename,
            filemode='w',
            format='%(name)s-%(levelname)s-%(asctime)s-%(message)s',
            datefmt='%I:%M:%S %p',  # set 'datefmt' to hour-minute-second AM/PM
            level='DEBUG'
        )
    # return an instance of the configured logger
    return logging.getLogger(logger_name)


if __name__ == "__main__":
    main(sys.argv)
