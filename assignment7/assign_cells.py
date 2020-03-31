#!/usr/bin/env python3

"""
   Assign cells to one of three groups: control, treatment, non-determined

   usage:
   assign_cells.py -m call.matrix.csv
"""

import pandas as pd
import argparse
import sys

def main(argv):
    """
       main method
    """
    args = parseArgs(argv)

    control_tag = 'ATGTTGC'
    treatment_tag = 'GATTACA'

    # create scoring matrix from file
    cell_matrix = pd.read_csv(args.matrix)

    # assign cells to categories control, treatment, ambiguous
    cell_dict = assignCells(cell_matrix)

    # display number of cells in each group
    print_msg = '\nCells assigned unambiguously:\n\tcontrol: {}\n\ttreatment: {}\nNon-determined cells: {}'
    print(print_msg.format(len(cell_dict['control']), len(cell_dict['treatment']), len(cell_dict['ambiguous'])))

    # examine non-determined cells. If hamming distance of cellTag is <=1, assign to appropriate group in cell_dict
    cell_dict_with_ambiguous = assignAmbiguousCells(cell_dict, cell_matrix, control_tag, treatment_tag)

def parseArgs(argv):
    """
        cmd line input
    """
    parser = argparse.ArgumentParser(description = "This script assigns cells to one of three groups: control, treatment, non-determined")
    parser.add_argument("-m", "--matrix", required=True,
                        help="[Required] Path to cell matrix")

    # remove script call from list (this is passed as list of cmd line input with script at position 0)
    args = parser.parse_args(argv[1:])
    return args

def assignCells(cell_df):
    """
    parse binarized matrix of cell barcodes (rows_ to cellTags (columns)
    :param cell_df:
    :return: dictionary storing barcode assignments to control, treatment or ambiguous
    """
    # create dictionary with appropriate keys
    cell_dict = {'control':[], 'treatment':[], 'ambiguous': []}
    # loop through rows, assign barcode to appropriate key
    for index,row in cell_df.iterrows():
        if row['ATGTTGC'] == 1:
            cell_dict.setdefault('control').append(row['cell.barcode'])
        elif row['GATTACA'] == 1:
            cell_dict.setdefault('treatment').append(row['cell.barcode'])
        else:
            cell_dict.setdefault('ambiguous').append(row['cell.barcode'])

    return cell_dict

def assignAmbiguousCells(cell_dict, cell_matrix, control_tag, treatment_tag):
    """
    Examine ambiguous cell cellTag and assign to appropriate group if hamming distance between ambiguous tag
    and known group is <= 1
    :param cell_dict: dictionary with three keys, control, treatment, and ambiguous, each associated with a list
    of barcodes corresponding to the cells in the cell_matrix
    :param cell_matrix: rows = cell barcodes, columns = unique cell tags. values are binary and represent
    whether a given cellTag is present in a given cell barcode population above a certain threshold
    :return: input cell_dict plus the cells with cellTags with hamming distance <= 1 added to the appropriate key group
    in dict.
    """

    # counters to track how many ambiguous barcodes are re-assigned to treatment or control
    control_count = 0
    treatment_count = 0
    # instantiate list to store barcodes re-assigned to treatment/control
    reassigned_barcode_list = []

    # loop through cell barcodes in cell_dict['ambiguous']
    for barcode in cell_dict['ambiguous']:
        # get the cell tag associated with the barcode
        cell_tag = getCellTag(barcode, cell_matrix)
        # calculate hamming distance from both control and treatment
        hamming_distance_control = hammingDistance(cell_tag, control_tag)
        hamming_distance_treatment = hammingDistance(cell_tag, treatment_tag)
        # if both hamming distances are greater than 1, skip to next cell barcode
        if hamming_distance_control & hamming_distance_treatment > 1:
            pass
        # if both are not greater than 1 and equal, pass (this would be a problem)
        elif hamming_distance_control == hamming_distance_treatment:
            pass
        # if hamming distance to control is less than or equal to one, assign to control, increment count, add barcode to list to remove from cell_dict['ambiguous']
        elif hamming_distance_control <= 1:
            cell_dict.setdefault('control').append(barcode)
            control_count = control_count + 1
            reassigned_barcode_list.append(barcode)
        # only option left is that treatment is <=1. Assign it to cell_dict treatment group, increment count, add barcode to list to remove from cell_dict['ambiguous']
        else:
            cell_dict.setdefault('treatment').append(barcode)
            treatment_count = treatment_count + 1
            reassigned_barcode_list.append(barcode)

    # remove assigned barcodes from cell_dict['ambiguous']
    cell_dict['ambiguous'] = [barcode for barcode in cell_dict['ambiguous'] if barcode not in reassigned_barcode_list]

    print_msg = "\nThe number of cells with hamming distance less than or equal to 1 re-assigned to a group:" \
                "\n\tcontrol: {}" \
                "\n\ttreatment: {}" \
                "\n\nThe number of barcodes that remain non-determined is: {}"

    print(print_msg.format(control_count, treatment_count, len(cell_dict['ambiguous'])))

    return cell_dict

def getCellTag(barcode, cell_matrix):
    """
    Get the cellTag associated with a given barcode
    :param barcode: cell barcode derived originally from cell_matrix, but passed to this method from cell_dict
    :param cell_matrix: the original matrix of genes x cellTags with binarized values reflecting whether the celltag is present in a given cell population above a certain threshold
    :return: the cellTag associated with a given barcode
    """
    df_row = cell_matrix[cell_matrix['cell.barcode'] == barcode] == 1
    # get the max (the true value will be 1, all else will be 0) column heading. Return the 1st value of the object (which is the celltag)
    cell_tag = df_row.apply(lambda x: df_row.columns[x.argmax()], axis=1).values[0]

    return cell_tag

def hammingDistance(str1,str2):
    """
    If two strings have the same character in a given position, they are considered the same. If they have different characters,
    the hamming distance is incremented. For strings of different lengths, hamming distance is calculated across characters from 0 to len(shortest),
    and then the difference in lengths is added.
    :param str1: one of two input strings to compare
    :param str2: one of two input strings to compare
    :return: the hamming distance btwn str1 str2
    """

    # get length of shortest string
    if len(str1) <= len(str2):
        shortest_str = len(str1)
        length_diff = len(str2) - len(str1)
    else:
        shortest_str = len(str2)
        length_diff = len(str1) - len(str2)

    hamming_distance = 0

    for i in range(shortest_str):
        if str1[i] != str2[i]:
            hamming_distance = hamming_distance + 1

    hamming_distance = hamming_distance + length_diff

    return hamming_distance


# call main method
if __name__ == "__main__":
    main(sys.argv)