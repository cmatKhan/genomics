from assignment8_tools.IndividualVariantObject import IndividualVariantObject
import sys
import re


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
            "...counting SNV and INDEL genotypes in individual %s..." % self.individual)  # copied code -- no good. should be in method
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
