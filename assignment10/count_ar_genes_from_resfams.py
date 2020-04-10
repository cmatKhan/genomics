#!/usr/bin/env python3
"""
    takes resfams_annotations.txt as input and calculates the num ber of antibiotic resistance genes identified
    in your contigs.

    usage: count_ar_genes_from_resfams.py -h <hmmscan_output>
"""

import sys
import argparse


def main(argv):
    """ main method
    :param argv: cmd line argument
    """
    # parse cmd line input
    args = parseArgs(argv)
    res_fam_output_path = args.res_fam_output

    # instantiate a counter
    count = 0

   # open the file
    with open(res_fam_output_path) as res_fam_file:
        # parse the file to a list of lines, loop over lines
        for line in res_fam_file.readlines():
            # if the line does not start with a #
            if not line.startswith('#'):
                # increment counter
                count = count + 1

    print('\nThere are %i genes identified in resfams_annotataions.txt\n' % count)

def parseArgs(argv):
    parser = argparse.ArgumentParser(
        description="count genes from resfam")
    parser.add_argument("-r", "--res_fam_output", required=True,
                        help="path to resfams_annotations.txt")

    args = parser.parse_args(argv[1:])
    return args


if __name__ == "__main__":
    main(sys.argv)