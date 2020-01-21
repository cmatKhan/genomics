#!/usr/bin/env python3

""" 
nuc_count.py counts nucleotides in a fasta file

Usage: nuc_count.py <fasta>

<fasta> = path to a fasta file
""" 

# Import modules
import sys

def main(argv):
    # main method
    # Args: one cmd line input of a fastq file with a single entry
    # Output: Parts 2, 3 and 5

    fasta = sys.argv[1]
    # get nucleotide sequence from fastq with one chromosome, all upper case. Print nuc frequencies and sequence length
    # return dictionary with sequence and nuc counts
    seq_dict = getSeq(fasta)
    # print the frequency of the nucleotides (part 3)
    printNucFreq(seq_dict)
    # print the frequency of the dinucleotides (part 5)
    printDinucFreq(seq_dict)

# sys.arg is a list containing 2 elements: the script name and 1 command line argument
# Check that all the command line argument was given. If not, print the documentation and exit.
#if (len(sys.argv) != 2):
#    sys.exit(__doc__)

# Save the input arguments as variables
#fasta = sys.argv[1]

def getSeq(fasta):
    # read in sequence
    # Args: a fasta with one chr(cmd line input)
    # Return: a dictionary of the sequence and nuc counts
    # Output: Prints raw counts and sequence length

    # Initialize a nucleotide string
    nucleotides = ""

    # Load the fasta sequence
    # NOTE: this script assumes there is only *one* sequence in the fasta file
    # Open the fasta file
    with open(fasta) as f:
        # For each line in the file
        for line in f:
            # Skip lines starting with ">"
            if not line.startswith(">"):
                # Add each line to the nucleotide string
                nucleotides += line.rstrip()

    # Make the nucleotide string all capital letters
    nucleotides = nucleotides.upper()

    # Count the nucleotides and print output
    num_a = nucleotides.count('A')
    num_c = nucleotides.count('C')
    num_g = nucleotides.count('G')
    num_t = nucleotides.count('T')
    num_n = nucleotides.count('N')

    print("Raw Counts")
    print("A: ", num_a)
    print("C: ", num_c)
    print("G: ", num_g)
    print("T: ", num_t)
    print("N: ", num_n)
    print("Sequence length: {}" .format(len(nucleotides)))

    # create dictionary to store sequence and raw counts
    seq_dict = {'seq': nucleotides, 'A': num_a, 'T': num_t, 'G': num_g, 'C': num_c, 'N': num_n}
    return seq_dict

## Part 3
def printNucFreq(seq_dict):
    # print nucleotide frequencies (not N)
    # Args: seq_dict with keys seq, A, T, G, C, N
    # Return: None
    # Output: print relative frequencies of nucleotides

    chr20_total = len(seq_dict['seq'])

    # print frequencies -- N is included in the total calculation
    for nuc in ["A", "T", "G", "C"]:
        nuc_freq = seq_dict[nuc]/chr20_total
        print("\n The frequency of {} is {}".format(nuc, nuc_freq))
    print("\n")

## Part 5
def printDinucFreq(seq_dict):
    # calculate dinucleotide frequencies
    # Args: seq_dict with keys seq, A, T, G, C, N
    # Return: None
    # Output: prints dinucleotide frequencies

    seq = seq_dict['seq']
    seq_length = len(seq)

    # create dictionary to store key, value pairs {dinucleotide: total}
    dinuc_dict = {}

    for i in range(seq_length-1):
        # concat the ith and i+1 element of string
        dinuc = seq[i] + seq[i+1]
        # pass if N
        if 'N' in dinuc:
            pass
        else:
            # test if given dinuc is in dictrion
            if (dinuc_dict.setdefault(dinuc, False)):
                # if so, add one to the value associated with the given dinuc
                dinuc_dict[dinuc] = dinuc_dict[dinuc] + 1
            else:
                # if not, create a new key (the dinuc) in the dictionary and associate it with a value of 1
                dinuc_dict[dinuc] = 1

    # take sum of dinucleotides
    total_dinuc = sum(dinuc_dict.values())

    print("The Frequency of dinucleotides is: ")

    # I did look up how to sort the keys of a dictionary. Answer found here: https://stackoverflow.com/questions/16600174/return-output-of-dictionary-to-alphabetical-order
    for key, value in sorted(dinuc_dict.items()):
        dinuc_freq = value/total_dinuc
        print("The frequency of dinucleotide sequence {} in chromosome 20 is {}". format(key, dinuc_freq))

# run main method
if __name__ == '__main__':
	main(sys.argv)