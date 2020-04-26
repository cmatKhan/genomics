#!/usr/bin/env python3
"""
    Script that takes in amino acid sequence as input and outputs
    all possible DNA sequences that could encode this AA.
    Usage: python3 Polk.py -p <peptide sequence>
"""

import sys
import os
import argparse
import logging

def main(argv):
    """
    main method
    :param argv: cmd line args
    """
    # create global variable to store aa_to_codons dict (available to all methods)
    global aa_to_codons
    aa_to_codons = AminoAcidToCodonDict()

    # parse cmd line input
    args = parseArgs(argv)
    peptide_sequences = args.peptide_sequence
    melting_temps = args.melting_temp

    # output file
    output_file = 'polk_script_output.txt'

    # prompt user to rename or overwrite previous output, if exists
    if os.path.exists(output_file):
        response = input('do you wish to overwrite %s? If not, exist and rename. Otherwise, enter \'y\': ' %output_file)
        if response == 'y':
            cmd = 'rm polk_script_output.txt'
            exit_status = os.system(cmd)
            if exit_status == '1':
                raise Exception('Could not remove polk_script_output.txt')

    for sequence in peptide_sequences:
        sequence = sequence.upper()
        for melting_temp in melting_temps:
            with open(output_file, 'a') as file:
                melting_temp = float(melting_temp)
                file.write("Amino Acid Sequence: %s\nMelting Temperature: %.2f\nDNA Sequences:\n" %(sequence, melting_temp))
            checkCombinations("", sequence, melting_temp, output_file)
            file.close()
            with open(output_file, 'a') as file:
                file.write('\n')
            file.close()

    # print output to stdout
    cmd = 'cat ./%s' %output_file
    exit_status = os.system(cmd)
    if exit_status == '1':
        raise Exception('Could print polk_script_output.txt. It may be that it is not being printed to $PWD')

def parseArgs(argv):
    parser = argparse.ArgumentParser(
        description="")
    parser.add_argument("-p", "--peptide_sequence", required=True, nargs='+',
                        help='[REQUIRED] The peptide sequence of interest')
    parser.add_argument("-m", "--melting_temp", required=True, nargs='+',
                        help='[REQUIRED] The required melting temperature')
    args = parser.parse_args(argv[1:])
    return args

def AminoAcidToCodonDict():
    """
    create amino acide to codon dictionary
    :return: dictionary of structure { <amino acid symbol>: [list of codons] }
    """
    aa_to_codons = {}
    aa_to_codons[ "A" ] = ["GCA", "GCC", "GCG", "GCT" ]
    aa_to_codons[ "C" ] = ["TGC", "TGT" ]
    aa_to_codons[ "D" ] = ["GAC", "GAT" ]
    aa_to_codons[ "E" ] = ["GAA", "GAG" ]
    aa_to_codons[ "F" ] = ["TTC", "TTT" ]
    aa_to_codons[ "G" ] = ["GGA", "GGC", "GGG", "GGT" ]
    aa_to_codons[ "H" ] = ["CAC", "CAT" ]
    aa_to_codons[ "I" ] = ["ATA", "ATC", "ATT" ]
    aa_to_codons[ "K" ] = ["AAA", "AAG" ]
    aa_to_codons[ "L" ] = ["CTA", "CTC", "CTG", "CTT", "TTA", "TTG" ]
    aa_to_codons[ "M" ] = ["ATG" ]
    aa_to_codons[ "N" ] = ["AAC", "AAT" ]
    aa_to_codons[ "P" ] = ["CCA", "CCC", "CCG", "CCT" ]
    aa_to_codons[ "Q" ] = ["CAA", "CAG" ]
    aa_to_codons[ "R" ] = ["AGA", "AGG", "CGA", "CGC", "CGG", "CGT" ]
    aa_to_codons[ "S" ] = ["AGC", "AGT", "TCA", "TCC", "TCG", "TCT" ]
    aa_to_codons[ "T" ] = ["ACA", "ACC", "ACG", "ACT" ]
    aa_to_codons[ "V" ] = ["GTA", "GTC", "GTG", "GTT" ]
    aa_to_codons[ "W" ] = ["TGG" ]
    aa_to_codons[ "Y" ] = ["TAC", "TAT" ]
    aa_to_codons[ "*" ] = ["TAA", "TAG", "TGA" ]

    return aa_to_codons

def calculateMeltingTemp(dna_string):
    """

    :param dna_string:
    :return:
    """
    num_G = float(dna_string.count('G'))
    num_C = float(dna_string.count('C'))
    return 64.9 + (41.0 * ((num_G + num_C -16.4) / float(len(dna_string))) )

def checkCombinations(dna_string, aa_string, melting_temp_threshold, output_file):
    """
    Translate aa_string to all possible DNA strings (given the a dictionary of AA: [list of synonymous codons]
    :param dna_string:
    :param aa_string:
    :param restriction_sites:
    :return: a list of all possible DNA strings for the given aa_string
    """
    # restriction site list
    restriction_sites = ['CATATG', 'CTCGAG', 'TCGA', 'CTAG']
    # base case
    if (len(aa_string)==0):
        # discard DNA strings with restriction sites
        restriction_site = False
        for seq in restriction_sites:
            if seq in dna_string:
                restriction_site = True
        if not restriction_site:
            melting_temp = calculateMeltingTemp(dna_string)
            # do not print DNA strings if the melting temp threshold is not satisfied
            if melting_temp_threshold - .5 <= melting_temp <= melting_temp_threshold + .5:
                with open(output_file, 'a') as file:
                    file.write('%s\t%.2f\t' %(dna_string, melting_temp))

    # if the base case has not been reached...
    elif len(aa_string) != 0:

        # look at the first amino acid in aa_string
        current_AA = aa_string[0]

        # iterate over codons in the list corresponding to the dictionary entry in aa_to_codons for the particular amino acid in current_AA
        for single_codon in aa_to_codons[current_AA]:

            # add the (next) codon in the list of codons for this aa
            new_dna_string = dna_string + single_codon

            # call the function on the aa_string minus the first amino acid in the sequence (recursive step)
            checkCombinations(new_dna_string, aa_string[1:], melting_temp_threshold, output_file)

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
            level='WARNING'
        )
    # return an instance of the configured logger
    return logging.getLogger(logger_name)


if __name__ == "__main__":
    main(sys.argv)