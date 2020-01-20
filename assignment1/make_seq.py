#!/usr/bin/env python3

"""make_seq.py prints a random sequence given a sequence length and nucleotide frequencies.  The random sequence will have the same nucleotide frequencies as the input nucleotide frequencies.

Usage: python3 make_seq.py <sequence_length> <a_freq> <c_freq> <g_freq> <t_freq>

<sequence_length> = Length of sequence to generate
<a_freq> = frequency of As
<c_freq> = frequency of Cs
<g_freq> = frequency of Gs
<t_freq> = frequency of Ts
"""

# Import modules 
import sys
import random

## Part 4: random nucleotide sequence

def getRandom(freq1, freq2):
	# generate a float btwn 0.0 and 1.0
	# Args: none
	# Return: a floating point number btwn 0.0 and 1.0

	# generate random
	random_dec = random.random()
	# do it again if the number is equal to either AT or GC total
	if random_dec == freq1 or random_dec == freq2:
		random_dec = getRandom(freq1, freq2)
	return random_dec

def chooseSides(wallflower, two_sets):
	# Choose a side based on which set a number belongs
	# Args: A value that doesn't know where it belongs
	# 		A dictionary with two key/values. These are the two possible sets that the wallflower may belong to
	# Return: A key of the dictionary two_sets representing the set to which the wallflower belongs
	# Constraint: there is a silly solution to comparing floats. couldn't think of anything better...need to look

	# python can't compare floats. turn the float into a big int and then run chooseSides
	wallflower_place = 'NN'
	if isinstance(wallflower, float):
		multiplier = 10000000
		wallflower = int(round((wallflower * multiplier), 7))
		for key, value in two_sets.items():
			two_sets[key] = round(value * multiplier, 7)
		wallflower_place = chooseSides(wallflower, two_sets)
	else:
		# iterate through the sets
		for key, value in two_sets.items():
			if wallflower < value:
				wallflower_place = key

	return wallflower_place

def createRandomSequence(sequence_len, a_freq, t_freq, g_freq, c_freq):
	# Initialize an empty string that nucleotides can be appended to
	random_nucleotide_seq = ''

	AT_freq = a_freq + t_freq
	GC_freq = g_freq + c_freq

	# dictionaries to store relative frequencies of bp to bp_totals
	bp_freq_dict = {"AT": AT_freq, "GC": GC_freq}
	# dictionaries to store relative frequencies of nucleotides to corresponding bp
	AT_freq_dict = {"A": a_freq / AT_freq, "T": t_freq / AT_freq}
	GC_freq_dict = {"G": g_freq / GC_freq, "C": c_freq / GC_freq}

	for i in range(sequence_length):
		# Generate a random decimal. If the random decimal is equal to AT or GC frequencies, it will be regenerated until it is not
		ran = getRandom(AT_freq, GC_freq)

		# based on relative frequencies of AT and GC, determine if my_ran is a AT or GC
		bp = chooseSides(ran, bp_freq_dict)

		# generate another random number
		if bp == "AT":
			another_ran = getRandom(AT_freq_dict["A"], AT_freq_dict["T"])
			nuc = chooseSides(another_ran, AT_freq_dict)
		else:
			another_ran = getRandom(GC_freq_dict["G"], GC_freq_dict["C"])
			nuc = chooseSides(another_ran, GC_freq_dict)

		# an array would be better. need numpy. not sure how else to append string in python (note to self to look)
		print(random_nucleotide_seq + ' and the next is... ' + nuc)
		random_nucleotide_seq = random_nucleotide_seq + nuc

	return random_nucleotide_seq

sequence_length = 100
a_freq = 0.28
c_freq = 0.28
g_freq = 0.21
t_freq = 0.23

# Check that frequencies add to 1. If not, exit the program
if (abs(a_freq + t_freq + c_freq + g_freq - 1) > 1e-4):
	sys.exit("ERROR: Nucleotide frequencies do not add up to 1!")

def main(sequence_length, a_freq, t_freq, g_freq, c_freq):
	#TODO: add error catching on chooseSides NN

	# sys.arg is a list containing 6 elements: the script name and 5 command line arguments
	# Check that all 5 command line arguments were given. If not, print the documentation and exit.
	#if (len(sys.argv) != 6):
	#	sys.exit("ERROR: incorrect number of arguments.\n" + __doc__)

	# Save the input arguments as variables
	# By default, the command line arguments are saved as strings. Convert them to numeric types.
	# sequence_length = int(sys.argv[1])
	#a_freq = float(sys.argv[2])
	#c_freq = float(sys.argv[3])
	#g_freq = float(sys.argv[4])
	#t_freq = float(sys.argv[5])

	random_sequence = createRandomSequence(sequence_length, a_freq, t_freq, g_freq, c_freq)
	print(random_sequence)

main(sequence_length, a_freq, t_freq, g_freq, c_freq)