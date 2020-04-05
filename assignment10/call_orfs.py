#!/usr/bin/env python3
"""
    scan all 6 reading frames of a given line in a fasta of contigs for ORFs, extract longest nonoverlapping ORFs from each reading frame in each contig
    input: fasta of contigs
    output:  in $PWD: all_orfs.fna -- longest non overlapping orfs in each contig
                      all_proteins.faa -- all_orfs sequences translated to amino acid sequence
    usage: call_orfs.py -f <contig_file>
"""
import sys
import argparse

def main(argv):
    """ main method
    :param argv: cmd line argument
    """
    # parse cmd line arguments
    args = parseArgs(argv)
    # extract longest orf from each contig
    longest_orf_dict = createLongestOrfDict(args.fasta_path)
    # translate longest_orf_dict sequences to amino acid sequences
    translated_longest_orf_dict = {fasta_header: translateSequenceToAminoAcid(nuc_seq)
                                   for fasta_header, nuc_seq in longest_orf_dict.items()}
    # write both to file
    writeFasta(longest_orf_dict, 'all_orfs.fna')
    writeFasta(translated_longest_orf_dict, 'all_proteins.faa')

def parseArgs(argv):
    parser = argparse.ArgumentParser(
        description="scan all 6 reading frames of a given line in a fasta of contigs for ORFs, extract longest ORF from each contig")
    parser.add_argument("-f", "--fasta_path", required=True,
                        help="[Required] Directory path of fastq files.\n")

    args = parser.parse_args(argv[1:])
    return args

def writeFasta(longest_orf_dict, filename):
    """ write longest_orf_dict to file in fasta format. header lines for each open reading frame will be:
        > longest_orf_dict[key] (see longestOpenReadingFrame for details)
    :param longest_orf_dict: see function longestOpenReadingFrame()
    :param filename: the name of the fasta file (if only basename.fasta given, output in $PWD. Give relative or absolute path if want to output elsewhere)
    :return: None. write to file. Any call to this function will overwrite a pre-existing file of the same name
    """
    with open(filename, 'w') as fasta_output:
        for key, sequence in longest_orf_dict.items():
            fasta_output.write('%s\n%s\n' % (key, sequence))

def translateSequenceToAminoAcid(orf_sequence):
    """ translate a sequence to amino acid. Assume sequence starts from first nucleotide in sequence
    :param orf_sequence: a sequence that starts with a start codon
    :return: translated sequence
    """
    # aa_dict stores {codon: amino_acid_symbol}
    aa_dict = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
               'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
               'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
               'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
               'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
               'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
               'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
               'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
               'TAT': 'Y', 'TAC': 'Y', 'TAA': 'STOP', 'TAG': 'STOP',
               'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
               'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
               'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
               'TGT': 'C', 'TGC': 'C', 'TGA': 'STOP', 'TGG': 'W',
               'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
               'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
               'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}
    # instantiate a list to hold translated amino acid sequence
    aa_seq_list = []
    # check that orf_sequence is a string and starts with a start codon
    if not isinstance(orf_sequence, str):
        raise TypeError('the item passed to translateSequenceToAminoAcid string is not a string. Try again.')
    if not orf_sequence.startswith('ATG'):
        raise Exception('sequence does not start with start codon')
    if not aa_dict[orf_sequence[-3:]] == 'STOP':
        raise Exception('sequence does not end with a stop codon. Are you sure this is an ORF?')
    # if the sequence passes the tests above, translate
    else:
        # iterate overs sequence in chunks of 3
        for index in range(0, len(orf_sequence) - 2, 3):
            # check to make sure there is not a stop codon before the end of the orf
            if aa_dict[orf_sequence[index:index + 3]] == 'STOP' and not index + 3 == len(orf_sequence):
                raise Exception('There is a stop codon in the middle of your ORF! Something is fishy here...')
            # append amino acid to list
            aa_seq_list.append(aa_dict[orf_sequence[index:index + 3]])
    aa_seq_list.pop(-1)  # pop the stop codon off the list
    return ''.join(aa_seq_list)  # turn list into string, return

def createLongestOrfDict(fasta_path):
    """ read fasta of contigs, extract longest reading frame from each. Return dict in form
        {longest_orf_fasta[key] : nucleotide_sequence} <-- see longestOpenReadingFrame for dict details
    :param fasta_path: path to input fasta of contigs to evaluate for longest ORF
    :return: dict of longest ORF from each contig
    """
    longest_nonoverlapping_orf_dict = {}

    with open(fasta_path, 'r') as fasta_file:
        line = fasta_file.readline()
        while line:
            if line.startswith('>'):
                # extract relevant contig data from fasta header
                orf_fasta_contig_description = line.split('_')[0:5]
                contig_header = "_".join(orf_fasta_contig_description)
                # extract sequence
                sequence = fasta_file.readline().strip()
                orfs = longestOpenReadingFrame(contig_header, sequence)
                # add orfs to longest_orf_dict
                longest_nonoverlapping_orf_dict.update(orfs)

            # read next line and continue while loop, if line exists
            line = fasta_file.readline()

    return longest_nonoverlapping_orf_dict

def longestOpenReadingFrame(contig_header, sequence):
    """ read sequence, evaluate all 6 reading frames, return longest ORF
    :param sequence: a nucleotide sequence
    :param contig_header: extracted from input fasta. see main method
    :return: non overlapping ORF as dictionary {[<contig_header_<orf_number>_<direction>_<length>] : ATACC..}
             where orf_number is arbitrary sequential number of orfs found in a given sequence
    """
    # instantiate dictionary for forward and reverse orfs
    orf_dict = {'forward': [], 'reverse': []}
    # instantiate a dict to hold the orfs in the structure described in the :return: statement in docstring
    fasta_orf_dict = {}

    # forward reading frames
    for index in range(0, 3):
        orf_dict['forward'].extend(findLongestOpenReadingFrame(sequence[index:]))

    # take reverse complement of the sequence
    reverse_complement = reverseComplement(sequence)
    for index in range(0, 3):
        orf_dict['reverse'].extend(findLongestOpenReadingFrame(reverse_complement[index:]))

    for direction, orfs_list in orf_dict.items():
        orf_number = 0
        for orf in orfs_list:
            if not orf == '':  # skip empty strings -- these are from reading frames that did not contain orfs above the threshold
                orf_number = orf_number + 1 # increment orf_number
                fasta_header = '%s_%s_%s_%s' % (contig_header, orf_number, direction, len(orf)) # create fasta header
                fasta_orf_dict.setdefault(fasta_header, orf) # set entry in dict

    return fasta_orf_dict

def findLongestOpenReadingFrame(sequence, **kwargs):
    """
    recursive method to return non overlapping ORFs (the longest ORFs with unique stop codons) on a given sequence. This does not account for reading frames.
    To account for reading frame, feed in string as follows: findOpenReadingFrames(sequence), findOpenReadingFrames(sequence[1:), findOpenReadingFrames(sequence[2:])
    Do the same with the reverse complement of the sequence
       alternative method: iterate over 1 nucleotide at a time. This has the advantage that it could be run in parallel. also easier to test -- could be wrapped in method to look at all 6
    :param sequence: a sequence to search for longest ORF
    :param kwargs: used to pass longest_orf and orf_list during recursion
    :return: the longest ORF in sequence (note: this does not account for reading frames -- starts at position 0)
    """
    threshold = 100 # threshold below which orfs will be filtered out
    # store stop codons
    stop_codons = ['TAA', 'TAG', 'TGA']
    # store longest_orf and orf_list if passed. else, initiate variables
    try:
        longest_orf = kwargs['longest_orf']
    except KeyError:
        longest_orf = 1
    try:
        orf_list = kwargs['orf_list']
    except KeyError:
        orf_list = ['']

    # discard nucleotides in chunks of 3 if a start codon has not yet been found
    while not sequence[0:3] == 'ATG' and len(sequence) > 3:
        sequence = sequence[3:]
    # if it is still possible to find a longer ORF in the sequence, look for it.
    if not len(sequence) < threshold:
        # scan for stop codons
        for i in range(3, len(sequence) - 2, 3):
            # if a stop codon is found
            if sequence[i : i+3] in stop_codons:
                orf = sequence[0 : i+3]
                if len(orf) > threshold:  # filter orfs less than threshold
                    orf_list.append(orf)
                    # if this orf is longer than longest_orf, update longest_orf
                    if len(orf) > longest_orf:
                        longest_orf = len(orf)-3 # do not count stop codon
                # recursive step
                return findLongestOpenReadingFrame(sequence[i+3:], longest_orf = longest_orf, orf_list = orf_list)
    # return statement for base case
    if len(orf_list) > 1:
        orf_list.pop(0)
    return orf_list

def reverseComplement(sequence):
    """
    :return: reverse complement
    :raises TypeError: if sequence is not a str
    """
    if not isinstance(sequence, str):
        raise TypeError('the item passed to reverseComplement string is not a string. Try again.')
    return reverse(complement(sequence))

def complement(sequence):
    """ return complement of sequence
    :param sequence: a nucleotide sequence
    :return: base pair complement of string, all upper case
    :raises TypeError: if sequence is not a string
    """
    if not isinstance(sequence, str):
        raise TypeError('the item passed to reverseComplement string is not a string. Try again.')
    # complement dictionary
    complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    # cast sequence to upper case
    sequence_uppercase = sequence.upper()
    # instantiate an empty list
    sequence_complement = []
    # take complement of sequence
    for base in sequence_uppercase:
        sequence_complement.append(complement_dict[base])  # collect as list b/c strings are immutable
    # turn list into string, return
    return ''.join(sequence_complement)

def reverse(sequence):
    """ reverse the string
    :param sequence: any string, but expected to be a string of nucleotide bases
    :return: a reversed string
    :raises TypeError: if sequence is not a string
    """
    if not isinstance(sequence, str):
        raise TypeError('the item passed to reverse string is not a string. Try again.')
    return sequence[::-1]


if __name__ == "__main__":
    main(sys.argv)
