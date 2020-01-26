#!/usr/bin/env python3
"""
Align sequencing reads to a set of genes.

Usage: map_sequence_starter.py <cDNA_fasta> <seq_reads_file>

Arguments:
    <cDNA_fasta>     = Path to cDNA fasta file.  (Required.)
    <seq_reads_file> = Path to sequencing reads file.  (Required.)
"""

import sys

def reverse_complement(sequence):
    # This function returns the reverse complement of a sequence
    # Args: a sequence of [ATGC]
    # Returns: Reverse complement. If non nucleotide passed, returns N

    # define complement dictionary
    comp_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    # reverse the string
    seq_str = str(sequence) # cast to string
    reverse_seq = seq_str[::-1]
    # define reverse_comp
    reverse_comp = ''
    # iterate over reverse_comp, return key value of comp_dict if exists, otherwise return N
    for nuc in reverse_seq:
        if comp_dict.setdefault(nuc):
            reverse_comp += comp_dict[nuc]
        else:
            reverse_comp +='N'
    # return the reverse complement
    return reverse_comp


def read_cDNA_file_to_dict(filename):
    # read in a cDNA file, return sequence: gene dictionary
    # Args: cDNA filepath
    # Returns: {sequence: gene}

    # initialize dictionary
    cDNA_dict = {}

    # open file
    with open(filename, 'r') as cDNA_file:
        # loop through lines file
        for line in cDNA_file:
            # add gene/sequence(upper case) to dictionary. Strip whitespace from both
            if line.startswith('>'):
                gene_name = line[1:].rstrip()
                seq = next(cDNA_file).rstrip().upper()
                cDNA_dict.setdefault(gene_name, seq)

    # return dictionary
    return cDNA_dict

def create_twenty_five_mer_dict(cDNA_dict):
    # This function creates a dictionary of 25mer sequences from the cDNA sequences.  The dictionary keys are the 25mer sequences and the values are the corresponding gene names.
    # Args: a dictionary of gene: sequence
    # Returns: a dictionary of 25-mer sequences : gene

    # initialize dictionary
    k_mer_dict = {}

    # loop through key, value of cDNA dictionary
    for key, value in cDNA_dict.items():
        # from 0 to twenty five minus the length of the sequence
        for i in range(len(value) - 25):
            # take 25 characters starting at i and add to dictionary with value as the current key of the cDNA_dict
            k_mer_dict.setdefault(value[i:i + 25], key)

    # return dictionary
    return k_mer_dict

##############main loop###################

def main(args):
    # check that the correct number of arguments were given
    if (len(sys.argv) != 3):
        sys.exit(__doc__)

    # parse command line
    cDNA_file = sys.argv[1]
    seq_reads_file = sys.argv[2]

    # read in cDNAs
    cDNA_dict = read_cDNA_file_to_dict(cDNA_file)
    print("Read in the cDNAs")

    # create dictionary of {25-mer substrings of a given sequence : gene} eg. {(ATCC..GG): geneA, (TCCG...GA): geneA...}
    twenty_five_mer_dict = create_twenty_five_mer_dict(cDNA_dict)
    print("Created dictionary")

    # create copy of cDNA_dict with all values set to zero
    count_dict = {}
    for key in cDNA_dict:
        count_dict[key] = 0

    # create filehandle for the sequenced reads
    fh = open(seq_reads_file, "r")

    # iterate through the lines of the sequenced reads file
    for line in fh:
        # strip whitespace and newline char
        line = line.rstrip()
        # transform string to upper case
        line = line.upper()
        # if the given read is an exact match to one of the key(25-mer sequences) in the twenty_five_mer dict
        if line in twenty_five_mer_dict:
            # and if the gene (value) of the twenty_five_mer_dict exists is in the count dict
            if twenty_five_mer_dict[line] in count_dict:
                # increment the value associated with that gene in the count dict
                count_dict[twenty_five_mer_dict[line]] = count_dict[twenty_five_mer_dict[line]] + 1
        else:
            # if the read sequence is not found in the keys of the twenty_five_mer_dict, check the reverse complement
            rc_seq = reverse_complement(line)
            if rc_seq in twenty_five_mer_dict:
                if twenty_five_mer_dict[rc_seq] in count_dict:
                    count_dict[twenty_five_mer_dict[rc_seq]] = count_dict[twenty_five_mer_dict[rc_seq]] + 1
    # close the reads filehandle
    fh.close()

    # print out results
    print("Name\tReads\tReads per BP")
    for name, count in sorted(count_dict.items()):
        print(name + "\t" + str(count) + "\t" + str(float(count)/len(cDNA_dict[name])))

# main function call
if __name__ == "__main__":
    # pass the cmd line arguments to the main function
    main(sys.argv)