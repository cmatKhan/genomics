#!/usr/bin/env python3
"""Scan a DNA sequence to find putative binding sites

Usage: python3 scan_sequence.py <scoring_matrix> <sequence_file>

Args:
    scoring_matrix = Path to scoring matrix. The rows of the matrix correspond 
       to A, C, G, and T, and the columns correspond to positions
    sequence_file = Path to DNA sequence file.
"""
import sys

# Helper dictionary for reverse complementation
reverse_comp = { 'A':'T', 'C':'G', 'G':'C', 'T':'A' }

###############################################################
# Begin functions
###############################################################

def create_scoring_matrix_from_file(matrix_file):
    """TODO: write function docstring
    """
    file_data = [ x.split() for x in open(matrix_file).readlines() ]
    scoring_matrix = [ dict(A=float(a_score), C=float(c_score), G=float(g_score), 
        T=float(t_score)) for a_score, c_score, g_score, t_score in 
        zip(file_data[0], file_data[1], file_data[2], file_data[3])
    ]

    return scoring_matrix


def read_sequence_from_file(file_name):
    """TODO: write function docstring
    """
    return open(file_name).readlines()[0].strip()



def score_with_matrix(subseq, matrix):
    """TODO: write function docstring
    """
    return sum([ score[ base ] for score, base in zip(matrix, subseq)])


def get_reverse_complement(in_string):
    """TODO: write function docstring
    """
    return (''.join([ reverse_comp[x] for x in in_string[::-1] ]))

###############################################################
# End functions
###############################################################

###############################################################
# Begin main script
###############################################################

# Check the correct number of command line arguments
if(len(sys.argv)!= 3):
    sys.exit(__doc__)

score_matrix_file = sys.argv[1]
sequence_file = sys.argv[2]

# TODO: explain what this code does
score_matrix = create_scoring_matrix_from_file(score_matrix_file)
motif_width = len(score_matrix)
search_sequence = read_sequence_from_file(sequence_file)

# Calculate the number of matrix 'windows' for calculating sequence scores
last_index = len(search_sequence) - motif_width + 1

# TODO: explain what this code does
forward_hit_list = [ (i, score_with_matrix(search_sequence[i:i+motif_width], score_matrix)) for i in range(last_index) ]

# TODO: explain what this code does
if len(forward_hit_list) == 0:
    print("No threshold-exceeding hits found in the forward direction!")
else:
    print("orientation\tposition\tscore")
    for hit in forward_hit_list:
        print("forward\t{position:d}\t{score:.2f}".format(position=hit[0], score=hit[1]))


