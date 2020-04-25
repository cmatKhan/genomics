#!/usr/bin/env python3
"""Script that takes in amino acid sequence as input and outputs
all possible DNA sequences that could encode this AA.
Usage: python3 Polk.py <peptide sequence>"""

#Make sure # of arguments is correct
import sys
if len(sys.argv)!=2:
        sys.exit( __doc__)

#Standard table of codons - a dictionary of single-letter
#amino acid code to a list of codon choices

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

def check_combinations(dna_string, aa_string):
    ####Write a Doc string for this####
    # if this code is confusing to you, uncomment the print statement
    #print("Input DNA is:",dna_string,"Remaining AAs are:", aa_string, sep='\t')

    ####What's going on here####
    if (len(aa_string)==0):
        print(dna_string)

    ####What's going on here####
    else:

        ####What's going on here####
        current_AA = aa_string[0];

        ####What's going on here####
        for single_codon in aa_to_codons[current_AA]:

            ####What's going on here####
            new_dna_string = dna_string + single_codon

            ####What's going on here####
            check_combinations(new_dna_string,aa_string[1:])

# Main Script

check_combinations( "", sys.argv[1] )


