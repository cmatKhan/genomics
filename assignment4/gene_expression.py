#!/usr/bin/env python3
"""
TODO: Write a summary of what this script does. 
Usage: python3 gene_expression.py raw_counts.txt
"""

import sys
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

#TODO: Print out the doc string and exit if the number of input parameters is not correct

#######################
##Part 0 -- Functions##
#######################

# CPM counts per million
# Convert raw counts to counts per million (cpm)
# Raw count x[i,j] (from sample i, gene j)
# Total counts N[i] (from sample i)
# cpm[i,j] = (10^6)*x[i,j]/N[i]
def counts_per_million(dictionary, list_of_samples):
    # Initialize output dictionary cpm_dict
    cpm_dict = {}
    # Find N, the list of each sample's library size
    N = library_sizes(dictionary, list_of_samples)
    # Calculate cpm one gene at a time using list comprehension
    for k,v in dictionary.items():
        # Note: do NOT use 10^6, which equals 12
        # Use either 10**6 or 1000000 for 1 million
        # k is the key (gene name)
        # v is the value (raw counts of gene k)
        cpm_dict[k] = [(10**6)*x/n for x,n in zip(v,N)]
    return(cpm_dict)
# End of CPM function #

# Library size
# Calculate library size of each sample (e.g. sum of RNA-seq counts)
def library_sizes(dictionary, list_of_samples):
    # Get the total number of samples
    num_samples = len(list_of_samples)
    # Initialize N, a list to hold the total counts from each sample
    N = []
    # For loop to iterate over each sample 
    for i in range(num_samples):
        # Append a new float zero value for each sample (goes to index i)
        N.append(0.0)
        # For loop to iterate over each value in our dictionary
        for v in dictionary.values():
            # Get the count from the i index of this gene and add it to the total for sample i
            N[i] += v[i]
    # Return the list containing each sample's library size
    return(N)
# End of lbirary sizes function #

# TODO: Fill in the code for the translate_dictionary function.
# Translate dictionary
# Function to translate a dictionary from {gene:[list of counts by sample]} to {sample:[list of counts by gene]}
# The new dictionary will have one key for each sample and the value of each key will be a list of counts associated with that sample
# This function is used in another function called upper_quartile_norm()
# Each comment below should correspond to one line of code
def translate_dictionary(dictionary, list_of_samples):
    # Initialize a new dictionary called translated_dictionary with empty curly brackets
    translated_dictionary = {}
    # Sort the keys of the input dictionary and save it in a list called 'genes'

    # Use a for loop to iterate over each sample (using a numbered index, not sample names)

        # Get the current sample name from list_of_samples

        # Initialize a new key in translated_dictionary using the current sample name. Let the value be an empty list.

        # Use a for loop to iterate over each gene in your list of genes

            # Find the RNA-seq count associated with this gene for this sample

            # Append the count to the list of counts for this sample in translated_dictionary

    # Return translated_dictionary
    return(translated_dictionary)
# End of translate dictionary function #

# TODO: Add comments to the code in the upper_quartile_norm function to explain what each line does
# Upper quartile normalization
# Compute the upper quartile normalization of raw counts
# Raw count x[i,j] (from sample i, gene j)
# D[i] value corresponding to 75th percentile of raw counts (from sample i)
# Mean of D, or the mean of all upper quartile means
# Upper quartile normalized count for sample i and gene j is (Mean of D)*x[i,j]/D[i]
def upper_quartile_norm(dictionary, list_of_samples):
    # Add comment: 
    upper_quartile_norm_dictionary = {}
    # Add comment:
    num_samples = len(list_of_samples)
    # Add comment: 
    translated_dict_for_D = translate_dictionary(dictionary, list_of_samples)
    # Add comment: 
    D = []
    # Add comment: 
    for i in range(num_samples):
        # Add comment: 
        sample_name = list_of_samples[i]
        # Add comment: 
        D.append(np.percentile(translated_dict_for_D[sample_name],75))
    # Add comment: 
    meanD = np.mean(D)
    # Add comment:
    for k,v in dictionary.items():
        # Add comment:
        upper_quartile_norm_dictionary[k] = [meanD*x/d for x,d in zip(v,D)]
    # Add comment:
    return(upper_quartile_norm_dictionary)
# End of upper quartile normalization function #

# TODO: Write a function to calculate Fisher's Linear Discriminant (add comments, too!) for all genes in your count dictionary. Call your function fishers_linear_discriminant. Remember that functions should be able to work using different data sets, so make sure that your function would work with data from a different experiment with different numbers of before and after samples. 
# Return the top ten genes according to the highest FLD values
# Input to this function should be a count dictionary and two lists letting the function know which index values go with each group
# You will need to iterate over each gene in the dictionary and calculate the FLD of each gene
# FLD of a gene = ((m1-m2)^2)/((s1)^2 + (s2)^2)
# m1 = mean of the first group
# m2 = mean of the second group
# s1 = standard deviation of the first group
# s2 = standard deviation of the second group


####################
##End of functions##
####################

# These first lines of code will get the data imported and in the right format for the rest of the homework
# You need to do the rest of the work starting from Part 1 -- Data filtering

#open the data file
data_file = open(sys.argv[1],"r")
#first line of data file contains sample names -- store in a list
sample_list = data_file.readline().strip().split()[1:]
#initialize raw count dictionary
raw_counts_dict = {}
#add each gene and expression values to dictionary
for line in data_file:
    line_list = line.strip().split() #split line into list at whitespace, strip removes leading and trailing whitespace
    gene = line_list[0] #name of gene is the first thing in the list
    expression_values = [int(float(v)) for v in line_list[1:]] #values of gene expression follow the gene name
    raw_counts_dict[gene] = expression_values #add keys and values to the dictionary {gene:expression}
#close the data file
data_file.close()

############################
##Part 1 -- Data filtering##
############################

# Filter out genes with zero expression in all samples

# Filter out genes with 20 or more samples with cpm < 1

################################
##Part 2 -- Data visualization##
################################

# Plot library sizes (save file as library_size.png)

################################
##Part 3 -- Data normalization##
################################

# Normalize count data left after filtering steps

# Plot normalized library sizes (save file as library_size_normalzied.png)

##############################
##Part 4 -- Data exploration##
##############################

# Calculate Fisher's Linear Discriminant for each gene basd on the normalized count data

# Pick a gene to explore further and plot the mean expression level for the before and after groups (save file as mean_expression.png)

######################################
##Extra credit -- Going the distance##
######################################

# Calculate Euclidean distance matrix of samples, output the most/least related samples

# Plot a dendrogram using Euclidean distance (save file as dendrogram.png)
