#!/usr/bin/env python3

"""
    Compare to scoring matrix, extract highest affinity binding site

    usage: highest_affinity_binding_site.py -m
"""
import argparse
import os
import sys

def main(argv):
    """
       main method
    """
    args = parseArgs(argv)

    # create scoring matrix from file
    scoring_matrix = create_scoring_matrix_from_file(args.matrix)

    # retrieve highest scoring sequence and score
    sequence, score = get_highest_score_sequence(scoring_matrix)

    print("The highest scoring sequence is {} with a score of {}".format(sequence, score))



def parseArgs(argv):
    """
        cmd line input
    """
    parser = argparse.ArgumentParser(description="This script generates sbatch script and submits sbatch job.")
    parser.add_argument("-m", "--matrix", required=True,
                        help="[Required] Path to scoring matrix. matricies must be in nucleotide x (sequence length) format.")

    # remove script call from list (this is passed as list of cmd line input with script at position 0)
    args = parser.parse_args(argv[1:])
    return args

def create_scoring_matrix_from_file(matrix_file):
    """
        Reads columnar data from file, returns a dictionary representation of a matrix with {nucleotide: score}
        representing the columns in the data
        Args: matrix_file
        Returns: a dictionary representation of a 4 x p matrix where the rows are the 4 nucleotide bases
        and p is the length of the sequence
    """
    file_data = [ x.split() for x in open(matrix_file).readlines() ]
    scoring_matrix = [ dict(A=float(a_score), C=float(c_score), G=float(g_score),
        T=float(t_score)) for a_score, c_score, g_score, t_score in
        zip(file_data[0], file_data[1], file_data[2], file_data[3])
    ]

    return scoring_matrix

def get_highest_score_sequence(matrix):
    """
        modified create_scoring_matrix_from_file() from scan_sequence.py
        based on https://stackoverflow.com/a/268285/9708266
        Args: a matrix with columns represented as dictionaries where the key is the rowname
        Returns: a string length p representing the sequence with the highest score across the n x p matrix
    """
    sequence_score_tuple = [[max(col, key=lambda key: col[key]), col[max(col, key=lambda key: col[key])]] for col in
                            matrix]
    sequence = ''.join([x[0] for x in sequence_score_tuple])
    score = sum([x[1] for x in sequence_score_tuple])

    return sequence, score #''.join([max(col, key=lambda key: col[key]) for col in matrix])


# call main method
if __name__ == "__main__":
    main(sys.argv)