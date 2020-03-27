#!/usr/bin/env python3

"""
   count the different classes of genome variation (SNVs, small indels, as well as larger structural variants) for
   the individual NA12878. Plot the length distribution of small indels (output file: histogram_indels.png),
   large deletions (output file: histogram_deletions.png) and MEIs (output file: histogram_meis.png) observed in NA12878.

   usage: count_gv.py <SNV_indel VCF> <SV VCF>
"""
import argparse
import sys
import matplotlib.pyplot as plt


def main(argv):
    """
    main method
    """
    # parse cmd line arguments
    args = parseArgs(argv)

    # assign cmd line arguments to variables
    snv_file = args.small_genome_variations
    sv_file = args.structural_varations

    # list of categories of variants to quantify
    variant_categories = ['SNV', 'INDEL', 'DEL', 'DUP', 'INV', 'MEI', 'BND']

    # IndividualVariantObject for SNV file, individual NA12878
    snv_NA12878 = IndividualVariantObject('snv', snv_file, 'NA12878', variant_categories)
    snv_NA12878.quantifyVariants()
    # IndividualVariantObject for sv file, individual NA12878
    sv_NA12878 = IndividualVariantObject('sv', sv_file, 'NA12878', variant_categories)
    sv_NA12878.quantifyVariants()

    # combine snv_NA12878 and sv_NA12878 variant dicts and print
    combined_variant_count_dict = {}
    for variant_type in variant_categories:
        combined_variant_count_dict.setdefault(variant_type,
                                               snv_NA12878.variant_count_dict[variant_type][0] +
                                               sv_NA12878.variant_count_dict[variant_type][0])

    print("The following variations are present in NA12878\n{}\n".format(combined_variant_count_dict))
    # find proportion that are structural variants
    proportion_sv = (combined_variant_count_dict['DEL'] + combined_variant_count_dict['DUP'] +
                     combined_variant_count_dict['INV'] + combined_variant_count_dict['MEI'] +
                     combined_variant_count_dict['BND']) / float(sum(combined_variant_count_dict.values()))
    print("The proportion of genomic variants that are SVs is: %f" % proportion_sv)

    # do the same for the variant_length_distributions
    combined_length_dist_dict = {'INDEL': snv_NA12878.variant_length_dist['INDEL'],
                                 'DEL': sv_NA12878.variant_length_dist['DEL'],
                                 'MEI': sv_NA12878.variant_length_dist['MEI']}

    # plot the INDEL, DEL and MEI length distributions
    title = "%s_distribution" % 'INDEL'
    plt.hist(combined_length_dist_dict['INDEL'], bins=376)
    plt.xscale('log')
    plt.title("%s" % title, loc='left')
    plt.xlabel('length of variant')
    plt.ylabel('count')
    plt.savefig('%s.png' % 'histogram_indels')
    plt.close()

    title = "%s_distribution" % 'DEL'
    plt.hist(combined_length_dist_dict['DEL'], bins=[0, 10, 10 ** 2, 10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6])
    plt.xscale('log')
    plt.title("%s" % title, loc='left')
    plt.xlabel('length of variant')
    plt.ylabel('count')
    plt.savefig('%s.png' % 'histogram_deletions')
    plt.close()

    title = "%s_distribution" % 'MEI'
    plt.hist(combined_length_dist_dict['MEI'], bins=[0, 10, 10 ** 2, 10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6])
    plt.xscale('log')
    plt.title("%s" % title, loc='left')
    plt.xlabel('length of variant')
    plt.ylabel('count')
    plt.savefig('%s.png' % 'histogram_meis')
    plt.close()


def parseArgs(argv):
    """
    cmd line input
    """
    parser = argparse.ArgumentParser(
        description="This script counts the different classes of genome variation in the snv and sv files")
    parser.add_argument("-snv", "--small_genome_variations", required=True,
                        help="[Required] Path to SNV_indel.biallelic.vcf")
    parser.add_argument("-sv", "--structural_varations", required=True,
                        help="[Required] Path to sv.reclassed.filtered.vcf")

    # remove script call from list (this is passed as list of cmd line input with script at position 0)
    args = parser.parse_args(argv[1:])
    return args


#################
### classes #####
#################

class IndividualVariantObject:
    """

    """

    def __init__(self, vcf_type, vcf_path, individual, variant_categories, *args):
        """
        constructor
        :param vcf_type: either 'sv' or 'snv'
        :param vcf_path: path to vcf
        :param individual: eg NA12878
        :param variant_categories: a list of variant types eg ['SNV', 'INDEL', 'DEL', 'DUP', 'INV', 'MEI', 'BND']
        """
        # create class attributes
        self.variant_categories = variant_categories
        self.vcf_type = vcf_type
        self.vcf_path = vcf_path
        self.individual = individual

        self.variant_count_dict = {}
        # self.variant_length_dist = {'INDEL': {}, 'DEL': {}, 'MEI': {}}
        self.variant_length_dist = {'INDEL': [], 'DEL': [], 'MEI': []}
        self.individual_index = None
        self.ref_index = None
        self.alt_index = None
        self.info_index = None

        if args:
            self.parent_1 = args[0]
            self.parent_2 = args[1]
            try:
                self.threshold = args[2]
            except IndexError:
                pass

        self.setVariantCountDict()
        self.setVcfIndices()

    ## ~~ setter functions ~~ ##
    def setVariantCountDict(self):
        """
        set attribute and create variant_count_dict with structure {variant_type: [0]} where variant_type is an item from variant_categories
        """
        for variant_type in self.variant_categories:
            self.variant_count_dict.setdefault(variant_type, [0])

    def setVcfIndices(self):
        """

        :param vcf_format_line:
        :return:
        """
        if not (hasattr(self, 'vcf_path') and hasattr(self, 'vcf_type')):
            sys.exit('no vcf_path and/or file_type set. Cannot find indicies of attributes in vcf_file without both.')

        with open(self.vcf_path) as file:
            for line in file:
                # open only metadata lines
                if line.startswith('#CHROM'):
                    # extract index of individual of interest
                    setattr(self, 'individual_index', line.split('\t').index(self.individual))
                    # extract index of REF column
                    setattr(self, 'ref_index', line.split('\t').index('REF'))
                    # extract index of ALT column
                    setattr(self, 'alt_index', line.split('\t').index('ALT'))
                    # extract the info column (only used in sv)
                    if self.vcf_type == 'sv':
                        setattr(self, 'info_index', line.split('\t').index('INFO'))
                    # extract parent information if given (QuanitfyIndividualGenotype class)
                    elif hasattr(self, 'parent_1'):
                        setattr(self, 'parent_1_index', line.split('\t').index(self.parent_1))
                        setattr(self, 'parent_2_index', line.split('\t').index(self.parent_2))
                    # halt loop
                    break

    def quantifyVariants(self):  # TODO: make sure you combine these if there is more than one input!
        """
        iterate through each line of vcf and count
        :param args: a variable length list of IndividualVariantObjects
        :return: variant_count_dict with counts from NA12878 in vcf added
        """

        print("...counting {} in individual {}...".format(self.vcf_type, self.individual))
        with open(self.vcf_path) as file:
            for line in file:
                # skip metadata lines
                if not (line.startswith('#')):

                    # extract genotype and genotype metadata of individual of interest
                    genotype_only = self.extractGenotype(line, self.individual_index)
                    # if the individual of interest has the alternate (0/0 is no alternate)
                    if '1' in genotype_only:
                        # store line information
                        vcf_line_ref_alt = line.split('\t')[self.ref_index:self.alt_index + 1]
                        ref_alt_dict = self.createRefAltDict(vcf_line_ref_alt)
                        # if vcf_type is 'sv', the variant category will be found in info_index
                        if self.vcf_type == 'sv':
                            variant_category = self.extractVariantInfo(line, self.info_index)
                        # else, the vcf_type is 'snv' and we'll need to intuit the type of variant from the ref/alt column
                        else:
                            variant_category = self.classifySnvVariantType(ref_alt_dict)
                        # increment appropriate category in variant_count_dict
                        self.variant_count_dict[variant_category][0] = \
                            self.variant_count_dict[variant_category][0] + 1
                        # count lengths of variants for plotting

                        # hard coding in the variation categories for which to track length
                        length_dist_list = ['INDEL', 'DEL', 'MEI']
                        if variant_category in length_dist_list:
                            length_of_variant = self.getSvVariantLength(line, ref_alt_dict, self.info_index,
                                                                        variant_category)
                            length_of_variant = abs(int(length_of_variant))
                            self.variant_length_dist[variant_category].append(length_of_variant)
                            # try:
                            #     self.variant_length_dist[variant_category][length_of_variant] +=1
                            # except KeyError:
                            #     self.variant_length_dist[variant_category].setdefault(length_of_variant, 1)
        self.variant_count_dict['BND'][0] = int(self.variant_count_dict['BND'][0] / 2)
        print('complete')

    @staticmethod
    def classifySnvVariantType(ref_alt_dict):
        """
        classify SnvVariantType from the REF and ALT columns of a vcf file line
        :param ref_alt_dict: output of extractRefAltDict
        :return: the classification of the Snv variant type
        """

        if ref_alt_dict['ref'][1] == 1 and ref_alt_dict['alt'][1] == 1:
            return 'SNV'
        if ref_alt_dict['ref'][1] == 1 and ref_alt_dict['alt'][1] > 1:
            return 'INDEL'
        if ref_alt_dict['ref'][1] > 1 and ref_alt_dict['alt'][1] == 1:
            return 'INDEL'

    @staticmethod
    def extractVariantInfo(vcf_line, info_index):
        """

        """
        return vcf_line.split('\t')[info_index].split('=')[1].split(';')[0]

    @staticmethod
    def createRefAltDict(ref_alt_list):
        """
        set new attribute genotype_dict which contains a reference and alternate genotype and string lengths at a given location eg {'ref': ['A', 1], 'alt':['AT', 2]}
        :param ref_alt_list: a list of the values in the REF category and ALT category of the vcf line eg ['A', 'AT']
        :return: a dictionary of structure {'ref':[ref_genotype, ref_genotype_length], 'alt':[ditto]} where ref_genotype and
                 ref_genotype_length refer to the entries extracted from the vcf at REF and ALT and their string lengths
        """
        # vcf_genotype = vcf_line.split('\t')[ref_index:alt_index + 1]
        ref_alt_dict = dict(ref=[ref_alt_list[0], len(ref_alt_list[0])], alt=[ref_alt_list[1], len(ref_alt_list[1])])

        return ref_alt_dict

    @staticmethod
    def getSvVariantLength(sv_line, ref_alt_dict, info_index, variant_category):
        if variant_category == 'INDEL':
            ref_alt_lengths = [ref_alt_dict['ref'][1], ref_alt_dict['alt'][1]]
            return max(ref_alt_lengths)
        else:
            return sv_line.split('\t')[info_index].split(';')[2].split('=')[1]

    @staticmethod
    def extractGenotype(vcf_line, individual_index, genotype_only=True):
        """
        :param vcf_line: a line of vcf file
        :param individual_index: the index of the individual of interest eg NA12878
        :param genotype_only: if True, only returns genotype eg 1/1. If False, returns all metadata associated with the genotype
        :return: the value at individual_index in vcf_line
        """
        if genotype_only:
            return vcf_line.split('\t')[individual_index].split(':')[0]
        else:
            return vcf_line.split('\t')[individual_index]

    @staticmethod
    def extractVcfVariantInfo(vcf_line, info_index):
        """
        :param vcf_line: a non metadata line of the vcf
        :param info_index: the index of 'INFO'
        :return: the value at INFO in vcf_line
        """
        return vcf_line.split('\t')[info_index].split(';')[0].split('=')[1]

    @staticmethod
    def extractAlternateFeatureLength(vcf_type, ref_alt_dict):
        """
        get feature length of the alternate sequence
        :param: vcf_type: either snv or sv
        :param: ref_alt_dct: parsed from the vcf_line by extractRefAltDict
        :return: length of indel sequence
        """
        if vcf_type == 'snv':
            return max([ref_alt_dict['ref'][1], ref_alt_dict['alt'][1]])
        elif vcf_type == 'sv':
            pass
        else:
            sys.exit('vcf_tpye not recognized')

    @staticmethod
    def classifySnvVariantType(ref_alt_dict):
        """
        classify SnvVariantType from the REF and ALT columns of a vcf file line
        :param ref_alt_dict: output of extractRefAltDict
        :return: the classification of the Snv variant type
        """

        if ref_alt_dict['ref'][1] == 1 and ref_alt_dict['alt'][1] == 1:
            return 'SNV'
        if ref_alt_dict['ref'][1] == 1 and ref_alt_dict['alt'][1] > 1:
            return 'INDEL'
        if ref_alt_dict['ref'][1] > 1 and ref_alt_dict['alt'][1] == 1:
            return 'INDEL'

    @staticmethod
    def extractQualScore(line, individual_index):
        """
        extract the quality score for the genotype of a given individual
        :param line: a line from the vcf
        :param individual_index: the index of the individual to extract a quality score for
        :returns: a quality score for the genotype. If none exists in position 4 or the quality score is '.', return 0
        """
        if line.split('\t')[8].split(':')[3] == 'GQ':
            try:
                qual_score = line.split('\t')[individual_index].split(':')[3]
            except IndexError:
                qual_score = 0
            if qual_score == '.':
                qual_score = 0
            return qual_score


# call main method
if __name__ == "__main__":
    main(sys.argv)
