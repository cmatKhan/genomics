#! /usr/bin/env python3

"""
    count_barcodes.py filters out reads that do not perfectly match barcode and counts the total number of perfectly
    matching reads for each barcode
    usage: count_barcodes.py -b filtered_variant_to_barcode.tsv -f <cDNA.fq or pDNA.fq>
"""
import sys
import argparse
import logging


def main(argv):
    """ main method
    :param argv: cmd line arguments
    """
    # create logger
    global logger
    logger = createLogger('/home/chase/code/cmatkhan/genomics/assignment12/assignment12.log', __name__)
    # parse cmd line arguments


def parseArgs(argv):
    parser = argparse.ArgumentParser(
        description="description here")
    parser.add_argument("-b", "--filtered_variant_to_barcode", required=True,
                        help='path to filtered_variant_to_barcode.txt')
    parser.add_argument("-f", "--fastq_file", required=True,
                        help='path to fastq file')
    args = parser.parse_args(argv[1:])
    return args


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