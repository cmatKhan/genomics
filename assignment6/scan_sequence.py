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
nucleotide_complement = { 'A':'T', 'C':'G', 'G':'C', 'T':'A' }

###############################################################
# Begin functions
###############################################################

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


def read_sequence_from_file(file_name):
    """
        Read the first line of file_name, strip whitespace at end of line, return
        Args: a file with one line (intended to be nucleotide string)
        Returns: the first line of the file contents as a string
    """
    return open(file_name).readlines()[0].strip()



def score_with_matrix(subseq, matrix):
    """
        the matrix is stored as a dictionary {row: column_value} st in column 1,
        the dictionary looks like {'A': score, 'C': score, 'G': score, 'C': score}.
        The zip function makes a tuple of ({row x column dictionary}, nucleotide in subseq)
        the sum() function returns the sum of the scores of the nucleotides in subseq across the columns of the matrix
        Args: a subsequence and a nucleotide x len(subseq) scoring matrix
        Returns: the score of the subsequence across the matrix
    """
    return sum([ score[ base ] for score, base in zip(matrix, subseq)])

def get_reverse_complement(in_string):
    """
        Starting from the end of in_string, append the complement and return the result
        Args: a nucleotide string
        Returns: reverse complement of the nucleotide string
    """
    return (''.join([ nucleotide_complement[x] for x in in_string[::-1] ]))

def filter_hit_list(hit_list, threshold):
    """
        filter the forward/reverse hit_lists by a given threshold
        Args: hit_list -- a list of tuples describing position, score and sequence. threshold at which to keep tuples
        Returns: a filtered list
    """

###############################################################
# End functions
###############################################################

###############################################################
# Begin main script
###############################################################

# Check the correct number of command line arguments
#if(len(sys.argv)!= 3):
#    sys.exit(__doc__)

#score_matrix_file = sys.argv[1]
#sequence_file = sys.argv[2]

polymerase_score_matrix = '/home/chase/code/cmatkhan/genomics/assignment6/polymerase_score_matrix.txt'
tf_score_matrix = '/home/chase/code/cmatkhan/genomics/assignment6/tf_score_matrix.txt'
promoter1 = '/home/chase/code/cmatkhan/genomics/assignment6/promoter1.txt'
promoter2 = '/home/chase/code/cmatkhan/genomics/assignment6/promoter2.txt'

score_matrix_file = polymerase_score_matrix
sequence_file = promoter1

# read a file containing a columnar representation of a matrix, return a dictionary representation of the same matrix
score_matrix = create_scoring_matrix_from_file(score_matrix_file)
# take length of sequence (matrix is 4 x p where p is the length of the sequence
# and the matrix values represent the score of one of the four nucleotides at each position in the sequence)
motif_width = len(score_matrix)
# read in the sequence (used later to search against the scoring matrix)
forward_search_sequence = read_sequence_from_file(sequence_file)
# take reverse complement
reverse_search_sequence = get_reverse_complement(forward_search_sequence)

# Calculate the number of matrix 'windows' for calculating sequence scores
last_index = len(forward_search_sequence) - motif_width + 1

# Using a substring of the inputted sequence, return a score derived from the scoring matrix. Scan the entire sequence moving one base at a time.
forward_hit_list = [ (i, score_with_matrix(forward_search_sequence[i:i + motif_width], score_matrix), forward_search_sequence[i:i + motif_width]) for i in range(last_index) ]
# do the same, but for the reverse complement
reverse_hit_list = [ (len(reverse_search_sequence) - i - motif_width, score_with_matrix(reverse_search_sequence[i:i+motif_width], score_matrix), reverse_search_sequence[i:i + motif_width]) for i in range(last_index) ]



# Print result to stdout. If either hit_list is 0, inform user
if len(forward_hit_list) == 0:
    print("No threshold-exceeding hits found in the forward direction!")
elif len(reverse_hit_list) == 0:
    print("No threshold-exceeding hits found in the reverse direction!")
else:
    print("orientation\tposition\tscore")
    for hit in forward_hit_list:
        print("forward\t{position:d}\t{score:.2f}\t{sequence:s}".format(position=hit[0], score=hit[1], sequence=hit[3]))



