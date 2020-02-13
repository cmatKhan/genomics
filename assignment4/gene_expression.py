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
import pandas as pd

#TODO: Print out the doc string and exit if the number of input parameters is not correct
#if not len(sys.argv[0:]) == 2:
#    sys.exit("Please input raw_counts.txt to the script.")


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
    genes = sorted(dictionary.keys())
    # Use a for loop to iterate over each sample (using a numbered index, not sample names)
    for i in range(len(list_of_samples)):
        # Get the current sample name from list_of_samples
        sample_name = list_of_samples[i]
        # Initialize a new key in translated_dictionary using the current sample name. Let the value be an empty list.
        translated_dictionary.setdefault(sample_name, [])
        # Use a for loop to iterate over each gene in your list of genes
        for gene_name in genes:
            # Find the RNA-seq count associated with this gene for this sample
            gene_count = dictionary[gene_name][i]
            # Append the count to the list of counts for this sample in translated_dictionary
            translated_dictionary.setdefault(sample_name, []).append(gene_count) # it is not necessary to setdefault twice -- this line is sufficient for both line 83 and 77

    # Return translated_dictionary
    return translated_dictionary
# End of translate dictionary function #

# TODO: Add comments to the code in the upper_quartile_norm function to explain what each line does
# Upper quartile normalization
# Compute the upper quartile normalization of raw counts
# Raw count x[i,j] (from sample i, gene j)
# D[i] value corresponding to 75th percentile of raw counts (from sample i)
# Mean of D, or the mean of all upper quartile means
# Upper quartile normalized count for sample i and gene j is (Mean of D)*x[i,j]/D[i]
def upper_quartile_norm(dictionary, list_of_samples):
    # Add comment: instantiate dictionary to hold {input_dictionary[key]: upper_quartile_norm}
    upper_quartile_norm_dictionary = {}
    # Add comment: instantiate list with length corresponding to the number of samples
    num_samples = len(list_of_samples)
    # Add comment: use the function above to translate a dictionary from {gene:[list of counts by sample]} to {sample:[list of counts by gene]}
    translated_dict_for_D = translate_dictionary(dictionary, list_of_samples)
    # Add comment: instantiate a list
    D = []
    # Add comment: loop over the number of samples from 0 to len(num_samples)-1
    for i in range(num_samples):
        # Add comment: store the sample name in location i in variable sample_name
        sample_name = list_of_samples[i]
        # Add comment: add the 75 percentile of the sample from the translated_dict to list D
        D.append(np.percentile(translated_dict_for_D[sample_name],75))
    # Add comment: take mean of list D
    meanD = np.mean(D)
    # Add comment: loop over key, value of dictionary
    for k,v in dictionary.items():
        # Add comment: associate the upper quartile mean with each key from dictionary.keys()
        upper_quartile_norm_dictionary[k] = [meanD*x/d for x,d in zip(v,D)]
    # Add comment: function output statement
    return upper_quartile_norm_dictionary
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


#open the data file TODO: uncomment this sys.arg line before submitting!
data_file = open("/home/chase/code/cmatkhan/genomics/assignment4/raw_counts.txt")
#data_file = open(sys.argv[1],"r")
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
def filter_counts(dict, threshold):
    # filter dicttion by removing gene if all counts for the given gene are below the inputted threshold
    # Args: dictionary in format {gene: [count of gene in all samples]}
    # Returns: filtered dictionary

    # instantiate list to hold keys to remove from dict based on filtering below
    fltr_dict = dict.copy()
    pop_list = []
    # loop through keys, values of raw_counts_dict
    for k, v in dict.items():
        # assume for the time being that the gene has counts below threshold
        pop_it = True
        # loop over counts in all samples for a given gene
        for count in v:
            # if a count does not equal 0, turn pop_it to false and move onto the next gene
            if count > threshold:
                pop_it = False
                pass
        # if pop_it still true, add key to list
        if pop_it:
            pop_list.append(k)
    for gene in pop_list:
        fltr_dict.pop(gene)

    return fltr_dict

no_zero_count_dict = filter_counts(raw_counts_dict, 1)

# Filter out genes with 20 or more samples with cpm < 1
def fltr_cpm(dict, sample_threshold, threshold):
    # return dict that has been filtered for cpm with counts less than threshold in samples greater than sample_threshold
    # Args: a dict with {gene: [sample counts]}, a number of samples for which if counts are less than threshold to discard, a threshold for cpm
    # Return: a filtered dictionary

    # copy the inputted dict
    fltr_dict = dict.copy()
    cpm_fltr_dict = counts_per_million(fltr_dict, sample_list)
    # instantiate a list to hold the genes (keys) to pop from the dict
    pop_list = []

    for gene, sample_counts in cpm_fltr_dict.items():
        # a counter to count samples below a given cpm threshold
        sample_count_below_threshold = 0
        for count in sample_counts:
            if count < threshold:
                # increment sample_count_below_threshold if a given count for a given sample is below threshold
                sample_count_below_threshold += sample_count_below_threshold + 1
        # if the number of samples with cpm below the given threshold is greater than sample_threshold
        if sample_count_below_threshold >= sample_threshold:
            # add the gene (dict key) to the pop_list
            pop_list.append(gene)

    # pop the genes identified in the steps above from fltr_dict
    for gene in pop_list:
        fltr_dict.pop(gene)
    # return fltr_dict
    return fltr_dict

fltr_cpm_dict = fltr_cpm(no_zero_count_dict,20, 1)
#for k,v in fltr_cpm_dict.items():
#    fltr_cpm_dict[k] = [float(x) for x in v]


################################
##Part 2 -- Data visualization##
################################

# Plot library sizes (save file as library_size.png)
filtered_lib_sizes = library_sizes(fltr_cpm_dict,sample_list)

################################
##Part 3 -- Data normalization##
################################

# Normalize count data left after filtering steps
norm_fltr_dict = upper_quartile_norm(fltr_cpm_dict, sample_list)


# Plot normalized library sizes (save file as library_size_normalzied.png)
norm_lib_sizes = library_sizes(norm_fltr_dict, sample_list)

##############################
##Part 4 -- Data exploration##
##############################

def appendSummaryColumns(df, new_col_name, group_from, group_to, summary_stat):
    if summary_stat == "mean":
        df[new_col_name] = df.loc[:, group_from: group_to].apply(np.mean, axis=1)
    if summary_stat == 'sd':
        df[new_col_name] = df.loc[:, group_from: group_to].apply(np.std, axis=1)
    return df

# Calculate Fisher's Linear Discriminant for each gene basd on the normalized count data
def fldFormula(group1_mean, group1_sd, group2_mean, group2_sd):
    # calculates by formula fisher linear discriminant
    # Args: mean of group one, standard dev of group1, mean of group2, sd of group 2
    # Returns: calculation of FLD by formula (float)

    # cast denominator to a float for floating point division
    x = group1_mean - group1_mean**2 / float(group1_sd**2 + group2_sd**2)
    return x

def createFLDTable(norm_fltr_dict, list_of_samples):
    # calculate FLD between the before group (before_1 - 20) and after group (after_1 - 20)
    # Args: a normalized dictionary filtered by genes with zero count across samples and cpm < 1 in 20+ samples. NOTE: norm_fltr_dict must be in RAW COUNTS
    # Returns: a dataframe of gene by FLD

    gene_list = norm_fltr_dict.keys()
    trans_dict = translate_dictionary(norm_fltr_dict, sample_list)
    col_order = ["gene"]
    col_order.extend(list(trans_dict.keys()))
    sum_cols = ["Before_mean", "Before_sd", "After_mean", "After_sd"]
    col_order.extend(sum_cols)

    gene_by_sample_df = pd.DataFrame(trans_dict)
    gene_by_sample_df['gene'] = gene_list
    gene_by_sample_df = appendSummaryColumns(gene_by_sample_df, 'Before_mean', 'Before_1', 'Before_20', 'mean')
    gene_by_sample_df = appendSummaryColumns(gene_by_sample_df, 'After_mean', 'After_1', 'After_20', 'mean')
    gene_by_sample_df = appendSummaryColumns(gene_by_sample_df, 'Before_sd', 'Before_1', 'Before_20', 'sd')
    gene_by_sample_df = appendSummaryColumns(gene_by_sample_df, 'After_sd', 'After_1', 'After_20', 'sd')
    gene_by_sample_df = gene_by_sample_df[col_order]
    gene_by_sample_df.set_index('gene')

    fld_df = gene_by_sample_df.apply(lambda row: fldFormula(row.Before_mean, row.Before_sd, row.After_mean, row.After_sd), axis=1)

    gene_by_fld = pd.concat([gene_by_sample_df['gene'], fld_df], axis=1)
    gene_by_fld.columns = ["gene", "fld"]

    return gene_by_fld

fld_df = createFLDTable(norm_fltr_dict, sample_list)

fld_df
# Pick a gene to explore further and plot the mean expression level for the before and after groups (save file as mean_expression.png)

######################################
##Extra credit -- Going the distance##
######################################

# Calculate Euclidean distance matrix of samples, output the most/least related samples

# Plot a dendrogram using Euclidean distance (save file as dendrogram.png)
