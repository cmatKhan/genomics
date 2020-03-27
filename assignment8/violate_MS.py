#!/usr/bin/env python3

"""
count the number of variants that clearly violate the rules of Mendelian segregation,
given the trio’s relationships to one another. With the -t flag, the user may add a threshold.
If any of the three genotypes fall below this threshold, they will be filtered out.
If -t is not passed, no filtering takes place.

usage: violate_MS.py <SNV_indel VCF> [GQ threshold*] *default: no thresholding based on genotype quality score

"""

import argparse
import sys

def main(argv):
    """
    main method
    """
    args = parseArgs(argv)

    # assign cmd line arguments to variables
    snv_file = args.small_genome_variations

    # list of categories of variants to quantify
    variant_categories = ['SNV', 'INDEL', 'DEL', 'DUP', 'INV', 'MEI', 'BND']

    if args.threshold is not None:
        snv_NA12878 = QuantifyIndividualGenotype('snv', snv_file, 'NA12878', variant_categories, 'NA12891', 'NA12892', int(args.threshold))
    else:
        snv_NA12878 = QuantifyIndividualGenotype('snv', snv_file, 'NA12878', variant_categories, 'NA12891', 'NA12892')

    snv_NA12878.evalulateMedelianSegregation()

    print(snv_NA12878.ms_genotype_violation_dict)
    print('The total number of violations is %i' % sum(snv_NA12878.ms_genotype_violation_dict.values()))

def parseArgs(argv):
    """
    cmd line input
    """
    parser = argparse.ArgumentParser(description = "This script counts the different classes of genome variation in the snv and sv files")
    parser.add_argument("-snv", "--small_genome_variations", required=True,
                        help="[Required] Path to SNV_indel.biallelic.vcf")
    parser.add_argument("-t", "--threshold",
                        help="[Optional] quality threshold. Variants with quality below threshold will not be counted")

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
                    setattr(self,'alt_index', line.split('\t').index('ALT'))
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
                            length_of_variant = self.getSvVariantLength(line, ref_alt_dict, self.info_index, variant_category)
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
        #vcf_genotype = vcf_line.split('\t')[ref_index:alt_index + 1]
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


class QuantifyIndividualGenotype(IndividualVariantObject):
    def quantifyGenotype(self):
        setattr(self, 'snv_indel_genotype_dict',
                {'homozygous_ref': 0, 'homozygous_alt': 0, 'heterozygous': 0, '1+ missing allele': 0})
        if self.vcf_type == 'sv':
            sys.exit('This function only quantifies SN variants. Please input SNV')
        print(
            "...counting SNV and INDEL genotypes in individual %s..." % self.individual)  # copied code -- no good. should be in method
        with open(self.vcf_path) as file:
            for line in file:
                # skip metadata lines
                if not (line.startswith('#')):
                    # extract genotype and genotype metadata of individual of interest
                    genotype_only = self.extractGenotype(line, self.individual_index)
                    if genotype_only == '0/0':
                        self.snv_indel_genotype_dict['homozygous_ref'] += 1
                    elif genotype_only == '1/1':
                        self.snv_indel_genotype_dict['homozygous_alt'] += 1
                    elif genotype_only == '0/1' or genotype_only == '1/0':
                        self.snv_indel_genotype_dict['heterozygous'] += 1
                    elif '/' in genotype_only:
                        self.snv_indel_genotype_dict['1+ missing allele'] += 1

    def evalulateMedelianSegregation(self):
        """
        A function which, inelegantly, attempts to evaluate violations of mendelian segretation
        """
        # set violation dict
        setattr(self, 'ms_genotype_violation_dict', {})
        # make sure the right file was passed
        if self.vcf_type == 'sv':
            sys.exit('This function only quantifies SN variants. Please input SNV')
        print(
            "...evaluating mendelian segregation in %s..." % self.individual)  # copied code -- no good. should be in method
        with open(self.vcf_path) as file:
            for line in file:
                # skip metadata lines and sex chromosome
                if not (line.startswith('#')) and line[0] is not '23':
                    # if threshold is set, filter
                    if hasattr(self, 'threshold'):
                        NA12878_qual = int(self.extractQualScore(line, self.individual_index))
                        NA12891_qual = int(self.extractQualScore(line, self.parent_1_index))
                        NA12892_qual = int(self.extractQualScore(line, self.parent_2_index))
                        # if any of three are less than threshold)
                        if NA12878_qual < self.threshold or NA12891_qual < self.threshold or NA12892_qual < self.threshold:
                            # skip to next line
                            continue
                    # extract genotypes of 12878 and parents
                    NA12878_genotype = self.extractGenotype(line, self.individual_index)
                    # split the parent genotypes on the / and count the alleles
                    NA12891_alleles = self.extractGenotype(line, self.parent_1_index).split('/')
                    NA12892_alleles = self.extractGenotype(line, self.parent_2_index).split('/')
                    # get possible offspring genotype given the parents
                    possible_offspring = self.punnetSquare(NA12891_alleles, NA12892_alleles)
                    # if the NA12878 genotype or its reverse is not in the punnetSquare list, it violates MS
                    if NA12878_genotype not in possible_offspring and NA12878_genotype[2]+'/'+NA12878_genotype[0] not in possible_offspring:
                        try:
                            self.ms_genotype_violation_dict[NA12878_genotype] = \
                                self.ms_genotype_violation_dict[NA12878_genotype] + 1
                        except KeyError:
                            self.ms_genotype_violation_dict.setdefault(NA12878_genotype, 1)

    @staticmethod
    def punnetSquare(geno_1, geno_2):
        """
        return the results of a 2x2 punnet square eg 0,0 x 1,1 = [0/1, 0/1, 0/1, 0/1]
        :param geno_1: the genotype of parent 1
        :param geno_2: the genotype of parent 2
        :returns: a list of possible offspring genotypes
        """
        possible_offspring = []
        for allele_parent_1 in geno_1:
            for allele_parent_2 in geno_2:
                possible_offspring.append('%s/%s' % (allele_parent_1, allele_parent_2))
        return possible_offspring


# call main method
if __name__ == "__main__":
    main(sys.argv)