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

if not len(sys.argv[0:]) == 2:
    sys.exit("Please input raw_counts.txt to the script.")


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
    return cpm_dict
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
        #translated_dictionary.setdefault(sample_name, [])
        # Use a for loop to iterate over each gene in your list of genes
        for gene_name in genes:
            # Find the RNA-seq count associated with this gene for this sample
            gene_count = dictionary[gene_name][i]
            # Append the count to the list of counts for this sample in translated_dictionary
            translated_dictionary.setdefault(sample_name, []).append(gene_count)

    # Return translated_dictionary
    return translated_dictionary
# End of translate dictionary function #

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

# FLD function in three parts. appendSummaryColumns was just a quicker way for me to put mean and sd columns on the dataframe in createFLDTable.
# fldFormula the FLD formula. In createFLDTable, I apply this formula over the mean and sd rows of the gene by sample dataframe

def appendSummaryColumns(df, new_col_name, group_from, group_to, summary_stat):
    # return a dataframe with new column. I should write this so that I simply pass the function in rather than hard coding mean and sd. If this comment remains, I didn't do it.
    # Args: a dataframe, the new column name, group_from and group_to are the columns over which to apply the summary stat, and summary_stat is either "mean" or "sd" TODO: pass function rather than hard code mean/sd
    # Return the dataframe with the appended columns

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
    num = (group1_mean - group2_mean)**2
    denom = float(group1_sd**2 + group2_sd**2)
    x = num / denom
    return x

def createFLDTable(norm_fltr_dict, print_full = False):
    # calculate FLD between the before group (before_1 - 20) and after group (after_1 - 20)
    # Args: a normalized dictionary filtered by genes with zero count across samples and cpm < 1 in 20+ samples. NOTE: norm_fltr_dict must be in RAW COUNTS
    # Returns: a dataframe of gene by FLD

    gene_list = sorted(norm_fltr_dict.keys())
    trans_dict = translate_dictionary(norm_fltr_dict, sample_list)
    col_order = ["gene"]
    col_order.extend(list(trans_dict.keys()))
    sum_cols = ["Before_mean", "Before_sd", "After_mean", "After_sd"]
    col_order.extend(sum_cols)

    gene_by_sample_mtrx = np.array([trans_dict[i] for i in sample_list]).transpose()
    gene_by_sample_df = pd.DataFrame(gene_by_sample_mtrx, columns=sample_list)
    gene_by_sample_df['gene'] = gene_list
    gene_by_sample_df = appendSummaryColumns(gene_by_sample_df, 'Before_mean', 'Before_1', 'Before_20', 'mean')
    gene_by_sample_df = appendSummaryColumns(gene_by_sample_df, 'After_mean', 'After_1', 'After_20', 'mean')
    gene_by_sample_df = appendSummaryColumns(gene_by_sample_df, 'Before_sd', 'Before_1', 'Before_20', 'sd')
    gene_by_sample_df = appendSummaryColumns(gene_by_sample_df, 'After_sd', 'After_1', 'After_20', 'sd')
    gene_by_sample_df = gene_by_sample_df[col_order]
    gene_by_sample_df.set_index('gene')
    if print_full:
        return gene_by_sample_df

    fld_df = gene_by_sample_df.apply(lambda row: fldFormula(row.Before_mean, row.Before_sd, row.After_mean, row.After_sd), axis=1)

    gene_by_sample_df['fld'] = fld_df

    gene_by_fld = pd.concat([gene_by_sample_df['gene'], fld_df], axis=1)
    gene_by_fld.columns = ["gene", "fld"]

    return gene_by_fld.sort_values(by=['fld'], ascending=False).iloc[:10,:]

####################
##End of functions##
####################


# open the data file
data_file = open(sys.argv[1],"r")
# first line of data file contains sample names -- store in a list
sample_list = data_file.readline().strip().split()[1:]
# initialize raw count dictionary
raw_counts_dict = {}
# add each gene and expression values to dictionary
for line in data_file:
    line_list = line.strip().split() #split line into list at whitespace, strip removes leading and trailing whitespace
    gene = line_list[0] #name of gene is the first thing in the list
    expression_values = [int(float(v)) for v in line_list[1:]] #values of gene expression follow the gene name
    raw_counts_dict[gene] = expression_values #add keys and values to the dictionary {gene:expression}
# close the data file
data_file.close()

############################
##Part 1 -- Data filtering##
############################

# Filter out genes with zero expression in all samples
def filterZeroCount(dict):
    # filter dictionary by removing gene if all counts for the given gene are below the inputted threshold
    # Args: dictionary in format {gene: [count of gene in all samples]}
    # Returns: filtered dictionary

    # instantiate list to hold keys to remove from dict based on filtering below
    fltr_dict = dict.copy()
    pop_list = []
    # loop through keys, values of raw_counts_dict
    for k, v in dict.items():
        # if the sum of the gene counts is zero, then add to the poplist
        if sum(v) == 0:
            pop_list.append(k)
    for gene in pop_list:
        fltr_dict.pop(gene)

    return fltr_dict

no_zero_count_dict = filterZeroCount(raw_counts_dict)

# Filter out genes with 20 or more samples with cpm < 1
def fltrCpm(dict, sample_threshold, threshold):
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
                sample_count_below_threshold = sample_count_below_threshold + 1
        # if the number of samples with cpm below the given threshold is greater than sample_threshold
        if sample_count_below_threshold >= sample_threshold:
            # add the gene (dict key) to the pop_list
            pop_list.append(gene)

    # pop the genes identified in the steps above from fltr_dict
    for gene in pop_list:
        fltr_dict.pop(gene)
    return fltr_dict

fltr_cpm_dict = fltrCpm(no_zero_count_dict, 20, 1)

range_list = [max(value) - min(value) for key, value in fltr_cpm_dict.items()]
print("\nThe minimum range of a gene by raw count is {}\nThe maximum range of a gene by raw count is {}".format(min(range_list), max(range_list)))


################################
##Part 2 -- Data visualization##
################################

# See R code. Plot called: library_size.png
filtered_lib_sizes = library_sizes(fltr_cpm_dict,sample_list)
# write out list to plot with R
with open('./filtered_lib_size.txt', 'w') as f:
    for item in my_list:
        f.write("%s\n" % item)


################################
##Part 3 -- Data normalization##
################################

# Normalize count data left after filtering steps
norm_fltr_dict = upper_quartile_norm(fltr_cpm_dict, sample_list)
with open('./norm_lib_size.txt', 'w') as f:
    for item in my_list:
        f.write("%s\n" % item)


# See R script. Plot called: library_size_normalzied.png)
norm_lib_sizes = library_sizes(norm_fltr_dict, sample_list)
print("\nThe minimum range by gene after normalization is {}\nThe maximum range by gene after normalization is {}".format(min(norm_lib_sizes), max(norm_lib_sizes)))

##############################
##Part 4 -- Data exploration##
##############################

fld_df = createFLDTable(norm_fltr_dict, sample_list)

print("\nThe top 10 differentially expressed genes by FLD are:\n")
print(fld_df)


# Pick a gene to explore further and plot the mean expression level for the before and after groups (save file as mean_expression.png)
norm_mean_df = createFLDTable(norm_fltr_dict, print_full=True)
gene_of_interest_df = norm_mean_df.loc['RAB30',:]
gene_of_interest_df.to_csv('./RAB30_exp.csv')

######################################
##Extra credit -- Going the distance##
######################################

# Calculate Euclidean distance matrix of samples, output the most/least related samples

# Plot a dendrogram using Euclidean distance (save file as dendrogram.png)
