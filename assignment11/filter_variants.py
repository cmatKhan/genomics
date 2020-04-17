#! /usr/bin/env python3

"""
    counts barcodes assigned to each variant(consider alternative and reference alleles independently),
    and writes only variants that havepassed a filter of at least 6 barcodes to an output file called
    filtered_variant_to_barcode.txt.​ The output file should use the same format as
    variant_to_barcode.txt​.

    usage: filter_variants.py -i variant_to_barcode.txt
"""
import sys
import argparse
import logging
import numpy as np

def main(argv):
    """ main method
    :param argv: cmd line arguments
    """
    # create logger
    global logger
    logger = createLogger('./log/assignment11.log',__name__)

    # parse cmd line arguments
    args = parseArgs(argv)
    variants_to_barcode_path = args.variant_to_barcode

    # read in file, return list of lines as iterable object
    line_generator = createLineGenerator(variants_to_barcode_path)
    # parse lines into dict retaining only unique barcodes
    variant_dict = createVariantBarcodeDict(line_generator)
    # write to file, applying filter
    writeFilteredVariantToBarcode(variant_dict, 6, 'filtered_variant_to_barcode.tsv')

def parseArgs(argv):
    parser = argparse.ArgumentParser(
        description="description here")
    parser.add_argument("-i", "--variant_to_barcode", required=True,
                        help='path to variant_to_barcode.txt')

    args = parser.parse_args(argv[1:])
    return args

def createLineGenerator(variant_to_barcode_path):
    """
    read in variant_to_barcode, return a iter object
    :param variant_to_barcode_path:
    :return: iter object
    """
    with open(variant_to_barcode_path) as file:
        lines = file.readlines()

    return iter(lines)

# split line on tab, take only unique in cols 1 and 2. return [name, ref_list, alt_list]
def createVariantBarcodeDict(line_generator):
    """
    iterate through line_generator, create dictionary with structure ('variant':
    :param line_generator:
    :return: dict in structure {variant_id_!: [[UNIQUE ref barcodes], [UNIQUE alt barcodes]], variant_id_2: ...}
    """
    barcode_dict = {}
    header = next(line_generator).strip()
    if not header.strip() == 'Variant_ID\tREF_barcode\tALT_barcode':
        raise Exception('The variant_to_barcode file does not have the expected header.\n'
                        'Check your file. The header needs to be present and in format:\n'
                        'Variant_ID\tREF_barcode\tALT_barcode')

    for line in line_generator:
        # strip white space and split on \t
        line = line.strip().split('\t')
        barcode_dict.setdefault(line[0], [])
        if not len(barcode_dict[line[0]]) == 0:
            logger.warning('the barcode %s is not unique' % line[0])
        # split each column on ':'
        if len(line) < 3:
            logger.warning('Variant %s has less than 2 sets of barcodes. '
                           'Check the formatting and enter a dummy value for the missing values in the line.' % line)
        else:
            line[1] = np.array(line[1].split(':'))
            barcode_dict[line[0]].append(list(np.unique(line[1])))
            line[2] = np.array(line[2].split(':'))
            barcode_dict[line[0]].append(list(np.unique(line[2])))

    return barcode_dict

def writeFilteredVariantToBarcode(variant_barcode_dict, threshold, output_filename):
    """
    write variant_barcode_dict to file, filtering on threshold
    :param variant_barcode_dict: see return explanation in createVariantBarcodeDict()
    :param threshold: number of (unique) barcodes below which to discard variant (filters strictly less than)
    :param output_filename: path to output file
    :return: none -- write to file
    """
    with open(output_filename, 'w') as output_file:
        output_file.write('Variant_ID\tREF_barcode\tALT_barcode\n')
        for variant_id, ref_alt in variant_barcode_dict.items():
            if len(ref_alt[0]) > threshold and len(ref_alt[1]) > threshold:
                line = variant_id + '\t' + ':'.join(ref_alt[0]) + '\t' + ':'.join(ref_alt[1]) + '\n'
                output_file.write(line)

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