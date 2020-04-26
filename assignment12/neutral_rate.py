#! /usr/bin/env python3

"""
    calculate the fraction of wobble positions that are conserved across all 4 species starting from a conserved start codon
    until the end of the sequence (which is a stop codon in all 4 species). Considers only third position wobbles.

    usage: neutral_rate.py -a PRE1_clustal_output.aln # note: the other flags are optional. enter neutral_rate.py -h to see description
"""
import sys
import argparse
import logging
import pandas as pd
from scipy import stats

def main(argv):
    """ main method
    :param argv: cmd line arguments
    """
    # create logger
    global logger
    logger = createLogger('assignment12.log', __name__)

    # parse cmd line arguments
    args = parseArgs(argv)
    alignment_file_path = args.alignment_file

    # parse apart the alignment file
    parsed_alignment_file = parseAlignmentFile(alignment_file_path)

    # determine rate of wobble position conserved across all 4 species starting from conserved start codon
    neutral_wobble_rate = calculateNeutralWobbleRate(parsed_alignment_file)
    print('\n 3rd position wobble rate: %f' %neutral_wobble_rate)

    # get number of positions that need to be conserved for the pval to be > .05
    significant_conservation_threshold = pvalGreaterThan05(neutral_wobble_rate, 10)
    print('\n %i or more out of 10 bp need to be conserved for the probability of this occuring to be <= .05\n' %significant_conservation_threshold)

    conserved_seq_df = examinePromoterAlignments(parsed_alignment_file, significant_conservation_threshold)

    printDataframeAndFasta(conserved_seq_df, 'conserved_promoter_seqs.fasta')



def parseArgs(argv):
    parser = argparse.ArgumentParser(
        description="Calculate wobble rate, search alignment file for putative promoter regions")
    parser.add_argument("-a", "--alignment_file", required=True,
                        help='[REQUIRED]path to .aln file')
    parser.add_argument("-m", "--min_bp_conserved",
                        help='[OPTIONAL] If this is not set, it will be calculated in the script. I suggest not setting this flag')
    parser.add_argument("-w", "--window_length", default=10,
                        help='[OPTIONAL] default is 10')
    args = parser.parse_args(argv[1:])
    return args

def parseAlignmentFile(alignment_file_path):
    """
    parse apart alignment file into dict structure {species: alignment_string} where alignment_string is either the aligned sequence or the *  ** track
    :param alignment_file_path:
    :return: a dict -- see structure in description above
    """
    # instantiate the dictionary
    alignment_dict = {'Skud': [], 'Smik': [], 'Scer': [], 'Sbay': [], 'alignment': []}

    # open the alignment file
    with open(alignment_file_path) as alignment_file:
        # read a line
        line = alignment_file.readline().strip()
        # loop over file while there are still lines to read in the file
        while line:
            if line.startswith('Skud'):
                line = line.strip()
                # enter the next 5 lines into the dictionary
                alignment_dict['Skud'].append(line[16:])
                line = alignment_file.readline().strip()
                alignment_dict['Smik'].append(line[16:])
                line = alignment_file.readline().strip()
                alignment_dict['Scer'].append(line[16:])
                line = alignment_file.readline().strip()
                alignment_dict['Sbay'].append(line[16:])
                line = alignment_file.readline()
                # cannot strip line b/c it would strip off final whitespace, which we need. -1 removes \n
                alignment_dict['alignment'].append(line[16:-1])
            line = alignment_file.readline()

        # instantiate a variable to hold the length of the first sequence once it has been joined
        aligned_seq_length = 0
        # join each list into a string
        for key in alignment_dict:
            alignment_dict[key] = ''.join(alignment_dict[key])
            # check that the strings in alignment_dict are all the same length
            if aligned_seq_length == 0:
                aligned_seq_length = len(alignment_dict[key])
            elif not len(alignment_dict[key]) == aligned_seq_length:
                # throw an error if any of the sequences are not the same length as the first
                raise Exception('The length of %s is not the same as the previous entry in alignment_dict. There is a mistake.')

        return alignment_dict


def calculateNeutralWobbleRate(alignment_dict):
    """
    calculate wobble rate
    :param alignment_dict: output of parseAlignmentFile(). See parseAlignmentFile() for details
    :return:
    """
    # mth_trp do not wobble
    mth_trp = ['TGG', 'ATG']
    # extract length of the aligned sequence
    aligned_seq_length = len(alignment_dict['alignment'])
    # instantiate variables to store counts
    conserved_wobble_position = 0
    unconserved_wobble_position = 0
    # dict to examine wobble positions
    synonymous_codon_dict = {'TT': [['T', 'C'], ['A', 'G']],
                             'CT': [['T', 'C', 'A', 'G']],
                             'AT': [['T', 'C', 'A']],
                             'GT': [['T', 'C', 'A', 'G']],
                             'TC': [['T', 'C', 'A', 'G']],
                             'CC': [['T', 'C', 'A', 'G']],
                             'AC': [['T', 'C', 'A', 'G']],
                             'GC': [['T', 'C', 'A', 'G']],
                             'TA': [['T', 'C'], ['A', 'G']],
                             'CA': [['T', 'C'], ['A', 'G']],
                             'AA': [['T', 'C'], ['A', 'G']],
                             'GA': [['T', 'C'], ['A', 'G']],
                             'TG': [['T', 'C'], ['A', 'G']],
                             'CG': [['T', 'C', 'A', 'G']],
                             'AG': [['T', 'C'], ['A', 'G']],
                             'GG': [['T', 'C', 'A', 'G']]}

    # 706 is the index of the aligned open reading frame
    for index in range(706, aligned_seq_length-3, 3):
        # look in the alignment track for ***
        if alignment_dict['alignment'][index : index + 3] == '***':
            logger.debug('conserved: %s' %alignment_dict['Skud'][index : index + 3])
            # if the codon is not mth or trp
            if not alignment_dict['Skud'][index : index + 3] in mth_trp:
                # increment conserved_wobble_position
                conserved_wobble_position += 1
        # if there is a wobble
        elif alignment_dict['alignment'][index : index + 3] == '** ':
            logger.debug([alignment_dict['Skud'][index : index + 3], alignment_dict['Smik'][index : index + 3],
                         alignment_dict['Scer'][index : index + 3], alignment_dict['Sbay'][index : index + 3]])
            # extract the first two nucleotides
            first_two_nucs = alignment_dict['Skud'][index : index + 2]
            # and the wobbler
            wobble_nucs = {alignment_dict['Skud'][index + 2], alignment_dict['Smik'][index + 2],
                           alignment_dict['Scer'][index + 2], alignment_dict['Sbay'][index + 2]}
            # look to see if the wobbler is a synonymous mutation
            for key in synonymous_codon_dict:
                if key == first_two_nucs:
                    for synon_wobbles in synonymous_codon_dict[key]:
                        if wobble_nucs.issubset(set(synon_wobbles)):
                            unconserved_wobble_position += 1

    return conserved_wobble_position / float(conserved_wobble_position + unconserved_wobble_position)

def pvalGreaterThan05(neutral_wobble_rate, n):
    """
    return the lowest value above which the probability is less than .05 given a total number of items from which to choose and the probability
    :param neutral_wobble_rate: calculated by calculateNeutralWobbleRate() above
    :param n: The total number of items from which to choose
    :return: The number of items to choose from n in order to satisfy p <= .05
    """
    # instantiate a variable with a dummy value to store the significant_threshold when it is found
    significant_threshold = -1
    # iterate over the possible values for k in B(k; n,p)
    for i in range(1,n+1,1):
        # calculate the probability. If it is <= .05, store in significant_threshold and stop iterating
        if stats.binom_test(i, n, neutral_wobble_rate, 'greater') <= .05: # 'greater' since we are only concerned with the right side
            significant_threshold = i
            break
    # raise an error if there is no value which equals or excees .05
    if significant_threshold == -1:
        raise Exception('At the given neutral_wobble_rate of %f, the p value of a full match is still above .05' %neutral_wobble_rate)

    return significant_threshold


def examinePromoterAlignments(alignment_dict, length_threshold):
    """
    Scan promoter sequences for conserved sequences of lengths >= 10
    :param alignment_dict:
    :param length_threshold:
    :return:
    """
    # set minimum acceptable length of window to call conserved
    minimum_match = 10

    # instantiate dataframe
    col_names = ['Scer_sequence', 'start_position', 'stop_position']
    conserved_seq_df = pd.DataFrame(columns=col_names)

    # a variable used to track the start position of the first window in a series that has 8/10 conserved positions
    start_index = -1
    # iterate from 0 to the aligned ORF (which is at index 706) minus the minimum match, which will be added to index to create the window
    for index in range(0,706-minimum_match, 1):
        # counter for how may * are in the window
        alignment_count = 0
        # create the window from index to minimum_match (in this case 10 bp) out
        alignment_window = alignment_dict['alignment'][index : index + minimum_match]
        # in the window, count the *
        for alignment_char in alignment_window:
            if alignment_char == '*':
                alignment_count += 1
        # if there are >= length_threshold (8 in this case) / 10 '*' and we have not identified the beginning of a possible promoter region yet
        if alignment_count >= length_threshold and start_index == -1:
            # set the start_index of the possible promoter region to the current index (see for loop)
            start_index = index
        # if the alignment_count for this window is less than the threshold and we are inside of a putative promoter region
        elif alignment_count < length_threshold and start_index != -1:
            # add the sequence to the dataframe
            row_dict = {'Scer_sequence': alignment_dict['Scer'][start_index : index + 9], # current index - 1 (move frame back 1) + 10 to get previous window
                        'start_position': start_index, 'stop_position': index + 9}
            conserved_seq_df = conserved_seq_df.append(row_dict, ignore_index=True)
            # reset the start_index
            start_index = -1

    return conserved_seq_df

def printDataframeAndFasta(df, fasta_file_path):
    """
    write S_cer_conserved.txt and format as a fasta for easy pasting into JASPAR
    :param df:
    :return:
    """

    df.to_csv('S_cer_conserved.csv', index=False)
    with open(fasta_file_path, 'w') as file:
        for row in df.itertuples():
            file.write('>%s_%s\n%s\n' %(row.start_position, row.stop_position, row.Scer_sequence))


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