#!/usr/bin/env python3
"""
    takes blast_to_card.txt as input and calculates the number of antibiotic resistance genes identified in
    contigs using an identity threshold of >80% amino acid identity over >85% of the subject sequence.

    usage: count_ar_genes_from_blast.py -b <blast output>
"""

import sys
import argparse
import pandas as pd

def main(argv):
    """ main method
    :param argv: cmd line argument
    """
    # parse cmd line input
    args = parseArgs(argv)
    blastp_output = args.blast_to_card
    output_filename = args.output_filename

    # instantiate a counter
    count = 0

    # open the file
    with open(blastp_output) as annotation_file:
        # store the lines in the file as a list, iterate over them
        for line in annotation_file.readlines():
            line = line.strip()  # strip whitespace and newline characters
            line = line.split('\t')  # split line on tabs
            # extract percent identity, move to next line if less than 80%
            percent_identity = float(line[2])
            if percent_identity < 80:
                continue
            # extract identity_over_subject_sequence percent, move to next line if less than 85%
            match_length = line[3]
            sequence_length = line[12]
            identity_over_subject_sequence = float(match_length) / float(sequence_length)
            if identity_over_subject_sequence < .85:
                continue
            count = count + 1 # if the line passes, increment count
            # if the line passes, add line to output file
            with open(output_filename, 'a') as output_file:
                output_file.write('\t'.join(line))

    print('\nThe total number of antibiotic genes with identity greater than 80 percent\n'
          'over greater than 85 percent of the target sequence is: %i\n' % count)

def parseArgs(argv):
    parser = argparse.ArgumentParser(
        description="calculate number of antibiotic resistent genes using id threshold >80% aa identity over >85% subject sequence")
    parser.add_argument("-b", "--blast_to_card", required=True,
                        help="blast_to_card.txt output from blastp script in assignment")
    parser.add_argument('-o', '--output_filename', default='ar_genes_from_blast.txt',
                        help = '[Optional] optionally name the output file. By Default, ar_genes_from_blast.txt in $PWD')

    args = parser.parse_args(argv[1:])
    return args


if __name__ == "__main__":
    main(sys.argv)