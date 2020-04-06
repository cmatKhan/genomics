#!/usr/bin/env python3
"""
   compare orfs identified by call_orfs.py and MetaMark
   # identified in both, unique to call_orfs, unique to metagenemark

   usage = compare_orf_callers -ca all_proteins.faa -cn all_orfs.fna -ma mgm_all_proteins.fna -mn mgm_all_orfs.fna
"""
import sys
import argparse

def main(argv):
    """ main method
    :param argv: cmd line arguments
    """
    # parse cmd line arguments
    args = parseArgs(argv)
    call_orfs_aa = args.call_orfs_aa
    call_orfs_nt = args.call_orfs_nt
    mgm_aa = args.mgm_aa
    mgm_nt = args.mgm_nt

    # create set dictionaries structure {call_orfs.py: {set of sequences}, metagenemark: {set of sequences}
    aa_set_dict = {'call_orfs_aa': extractSequencesAsList(call_orfs_aa), 'mgm_aa': extractSequencesAsList(mgm_aa)}
    nucleotide_set_dict = {'call_orfs_nt': extractSequencesAsList(call_orfs_nt), 'mgm_nucleotide': extractSequencesAsList(mgm_nt)}

    # compare sets in each dictionary
    set_list = [aa_set_dict, nucleotide_set_dict]
    for set_dict in set_list:
        compareSets(set_dict)


def parseArgs(argv):
    parser = argparse.ArgumentParser(
        description="compare orfs identified by call_orfs.py and MetaMark")
    parser.add_argument("-ca", "--call_orfs_aa", required=True,
                        help="call_orfs.py protein sequence output")
    parser.add_argument("-cn" "call_orfs_nt", required=True,
                        help='call_orfs.py aa sequence output')
    parser.add_argument("-ma", "--mgm_aa", required=True,
                        help="metagenemark protein sequence output")
    parser.add_argument("-mn" "mgm_nt", required=True,
                        help='metagenemark nucleotide sequence output')

    args = parser.parse_args(argv[1:])
    return args

def extractSequencesAsList(fasta):
    """
    extract sequences from a fasta as a set
    :param fasta: a fasta format file
    :return: a set of sequences
    """
    sequence_set = set()
    with open(fasta, 'r') as fasta_file:
        for line in fasta_file.readline():
            if not line.startswith('>'):
                sequence_set.add(line.strip())  # strip \n and whitespace
    return sequence_set

def compareSets(set_dict):
    """
    compare two sets, return summary
    :param set_dict: a dictionary of two sets. key should be name/description of set. This will be printed to stdout for comparison
    :return: None. print summary to std out
    """
    key_list = list(set_dict.keys())
    if not len(key_list) == 2:
        raise Exception('compareSets can only compare two sets. Check that the set_dict has only two keys')

    key_1 = key_list[0]
    key_2 = key_list[1]

    intersection = len(set_dict[key_1].intersect(set_dict[key_2]))
    key_1_exclusive = len(set_dict[key_1].difference(set_dict[key_2]))
    key_2_exclusive = len(set_dict[key_2].difference(set_dict[key_1]))

    print("\nThe number of items in both %s and %s: %i" % (key_1, key_2, intersection))
    print("The number of items in exclusively in %s: %i" % (key_1, key_1_exclusive))
    print("The number of items in exclusively in %s: %i\n" % (key_2, key_2_exclusive))


if __name__ == "__main__":
    main(sys.argv)
