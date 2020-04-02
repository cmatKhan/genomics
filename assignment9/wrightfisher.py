#!/usr/bin/env python3

"""
    estimate the approximate fitness of an allele through simulation. produce plots of one first iteration
    in each simulation in which mutant allele reaches fixation

    usage: wrightfisher.py -f <real_number> -g <int> -m <recessive/dominant>
           You may use -p <int> to change population_size. composition will aways initiate as population_size-1 homozygous wildtype,
           1 heterozygous mutant
    if -p is not passed, the program will produce the models requested for part 1 and part 2 using preset parameters. If p is passed,
    part1 and part2 models will be produced with the presets, and the user generated model will be simulated.
"""

import argparse
import sys
import random
from copy import deepcopy
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def main(argsv):
    """
    main method
    :param argsv: cmd line input. See function argParse
    """
    # parse cmd line argments
    args = parseArgs(argsv)
    fitness_multiplier = float(args.fitness_score)
    max_generations = int(args.max_generations)
    model = args.model
    if args.population_size:
        population_size = int(args.population_size)
    else:
        population_size = args.population_size

    print("Addressing part 1 using population size of 100, "
          "max generations 1000, model dominant, "
          "and fitness mutiplier 3")
    # instantiate WrightFisherModel with a dominant model,
    # fitness multiplier 3, population size 100 and composition 99 homozygous wildtype, 1 heterozygous mutant
    question_1_wfm = WrightFisherModel('dominant', 3)
    # run 100 simulations
    dominant_plot_dict = testFitnessScore(question_1_wfm, 1000, 'dominant')
    # plot
    writePlot(dominant_plot_dict, 100, 'Dominant', 'allelic_frequency_vs_generations_dominant.png')

    print("\nAddressing part 2 using population size of 100, "
          "max generations 1000, model recessive, "
          "and fitness multiplier 3")
    # instantiate WrightFisherModel with a recessive model,
    # fitness multiplier 3, population size 100 and composition 99 homozygous wildtype, 1 heterozygous mutant
    question_2_wfm = WrightFisherModel('recessive', 3)
    # run 100 simulations
    recessive_plot_dict = testFitnessScore(question_2_wfm, 1000, 'recessive')
    # plot
    writePlot(recessive_plot_dict, 100, 'Recessive', 'allelic_frequency_vs_generations_recessive.png')

    if not isinstance(population_size, int):
        sys.exit('\nParts 1 and 2 of the assignment are produced automtically with preset parameters '
                 '\nchosen according to the assignment. However, if you wish to produce a model based completely '
                 '\non user inputted parameters, you must enter all the parameters including population size')
    print("\nNow running the user input parameters")
    # instantiate user model
    user_wfm = WrightFisherModel(model, fitness_multiplier, population_size - 1, 0, 1, 0)
    # run 100 simulations
    user_plot_dict = testFitnessScore(user_wfm, max_generations, model)
    # plot
    writePlot(user_plot_dict, population_size, model, 'user_parameters_wfm.png')


def parseArgs(argv):
    """
    cmd line arguments
    :param argv: cmd line arguments. see add_arguments below
    :return: a keyword arguments parsed in main method
    """
    parser = argparse.ArgumentParser(description="Assignment 9")
    parser.add_argument("-f", "--fitness_score", required=True,
                        help="[Required] fitness score -- this is a multiplier of the chosen model")
    parser.add_argument("-g", "--max_generations", required=True,
                        help="[Required] maximum number of generations")
    parser.add_argument("-m", "--model", required=True,
                        help="[Required] Dominant or recessive")
    parser.add_argument("-p", "--population_size",
                        help="[Optional] population_size. Composition will always be only 1 heterozygous mutant. Default"
                             " is 99 homozygous wildtype, 1 heterozygous mutant")
    args = parser.parse_args(argv[1:])
    return args


def testFitnessScore(wfm, max_generations, model):
    """

    :param wfm: WrightFisherModel object. see classes below
    :param max_generations: total number of generations to generate from initial
    :param model: dominant or recessive (used in printout at end
    :return: plot_dict structure {'wt_allele_frequency': <int>, 'mutant_allele_frequency': <int>}
    """
    # number of simulations is fixed
    simulation_repetitions = 100
    # store number of times mutant allele fixes in 100 simulations
    number_mutant_fixation = 0
    # store a list initiated with 0 to track the first instance of mutant fixation. 0 is dropped before being returned in dictionary
    wt_allele_frequency = [0]
    mutant_allele_frequency = [0]
    homozygous_wildtype_frequency = [0]
    homozygous_mutant_frequency = [0]
    heterozygous_frequency = [0]
    # simulation_repetitions is decremented at the end of each iteration
    while simulation_repetitions > 0:
        # initially, fixation is false
        mutant_fixation = False
        wildtype_fixation = False
        # counter for number of generations. This is incremented after each iteration of the while loop below
        number_of_generations = 0
        # if the population reaches mutant fixation or max_generations, halt loop and move onto next iteration
        while mutant_fixation is False and wildtype_fixation is False and number_of_generations < max_generations:
            # mate the population to create the next generation
            wfm.matePopulation()
            # determine if wt or mutant is fixed, exit this iteration of the simulation if so
            mutant_fixation = wfm.determineMutantFixation()
            wildtype_fixation = wfm.determineWildtypeFixation()
            # store progress of allele until the first iteration at which it reaches fixation
            if mutant_allele_frequency[-1] < wfm.size * 2:
                mutant_allele_frequency.append(wfm.mutantAlleleFrequency())
                wt_allele_frequency.append(wfm.wildtypeAlleleFrequency())
                homozygous_wildtype_frequency.append(wfm.homozygousWildtypeFrequency())
                homozygous_mutant_frequency.append(wfm.homozygousMutantFrequency())
                heterozygous_frequency.append(wfm.heterozygousFrequency())
            # if the mutant reaches fixation,
            if mutant_fixation:
                # increment mutant_fixation
                number_mutant_fixation += 1
            # increment number_of_generations
            number_of_generations += 1
        # reset the population to initial state and start loop over if question_1_repeitions is not zero
        wfm.resetPopulation()
        # if mutant_allele_frequency less than wfm.size * 2, reset it also
        if mutant_allele_frequency[-1] < wfm.size * 2:
            wt_allele_frequency = [0]
            mutant_allele_frequency = [0]
            homozygous_wildtype_frequency = [0]
            homozygous_mutant_frequency = [0]
            heterozygous_frequency = [0]
        # decrement question_1_repitions
        simulation_repetitions = simulation_repetitions - 1

    print("fraction Bio5488 fixation in 100 simulations over %i "
          "generations using model \'%s\' with fitness score %4.2f: %.2f" % (max_generations, model,
                                                                             wfm.fitness_multiplier,
                                                                             number_mutant_fixation / 100))

    # return a dictionary of the wt_allele_frequency and mutant_allele_frequency in the first iteration in which mutant allele fixed
    return {'wt_allele_frequency': wt_allele_frequency[1:],
            'mutant_allele_frequency': mutant_allele_frequency[1:],
            'homozygous_wildtype_frequency': homozygous_wildtype_frequency[1:],
            'homozygous_mutant_frequency': homozygous_mutant_frequency[1:],
            'heterozygous_frequency': heterozygous_frequency[1:]}


def writePlot(plot_dict, population_size, model, filename):
    """
    plot values in plot_dict generated from function testFitnessScore above
    :param plot_dict: see return statement in docstring of testFitnessScore above
    :param population_size: population_size (present for question 1 and 2, otherwise cmd line argument)
    :param model: either dominant or recessive (present for question 1 and 2, otherwise cmd line argument)
    :param filename: filename of graph. written to present working directory
    """
    # check to ensure population_size dictionary contains information about alleles that reached fixation. print msg if not
    if not population_size * 2 == plot_dict['mutant_allele_frequency'][-1]:
        print('The mutant allele did not reach fixation in any of the simulations')
    else:
        # create plot
        plt.plot(plot_dict['mutant_allele_frequency'], 'm--', plot_dict['wt_allele_frequency'], 'k-',
                 plot_dict['homozygous_wildtype_frequency'], 'b-',
                 plot_dict['homozygous_mutant_frequency'], 'g--',
                 plot_dict['heterozygous_frequency'], 'y:')
        black_patch = mpatches.Patch(color='black', label='Wildtype Allele Frequency')
        magenta_patch = mpatches.Patch(color='magenta', label='Wildtype Allele Frequency')
        blue_patch = mpatches.Patch(color='blue', label='Homozygous Wildtype')
        green_patch = mpatches.Patch(color='green', label='Homozygous Mutant')
        yellow_patch = mpatches.Patch(color='yellow', label='Heterozygous')

        plt.legend(handles=[magenta_patch, black_patch, blue_patch, green_patch, yellow_patch], loc='best')
        plt.xlabel('Generations to fixations in first sample to reach Mutant Fixation')
        plt.ylabel('Allele Frequency')
        plt.title('%s Model' % model)
        plt.savefig(filename)
        plt.close()


class Individual:
    """
    A parent class for any individual, haploid or diploid
    Attributes:
        _possible_genotypes = a dictionary, eg {'homozygous_a': [0,0], 'homozygous_b': [1,1], 'heterzygous': self.heteroGenotype()}
        genotype =  a list, eg [0,1]. If haploid, still a list eg [0]
    Methods:
        donateAllele(self): randomly return an allele from genotype
    """
    _possible_genotypes = {}
    genotype = []

    def donateAllele(self):
        """
        randomly select an allele to return (handles either haploid or diploid individuals
        :return: one of the individual's allele. If diploid, randomly chooses.
        """
        # randomly generate a value that represents the position of an allele in self.genotype (a list)
        allele_index = random.randint(0, len(self.genotype) - 1)
        return self.genotype[allele_index]


class DiploidIndividual(Individual):
    """
    A child of Individual.
    Attributes:
        _possible_genotypes: set to {'homozygous_a': [0,0], 'homozygous_b': [1,1], 'heterzygous': self.heteroGenotype()}
        genotype: set based on parameter passed to constructor
    Methods:
        heteroGenotype(): static method, randomly generates either [0,1] or [1,0]
    """

    def __init__(self, genotype):
        """
        An object representing an individual in a population
        :param genotype: One of the following: [[0,0], [1,1], [0,1], [1,0]]
        """
        self._possible_genotypes = ['homozygous_wt', 'homozygous_mutant', 'heterozygous']
        if genotype not in self._possible_genotypes:
            raise Exception("%s is not a valid genotype_description. "
                            "Valid genotype_descriptions are: %s" % (genotype, self._possible_genotypes))
        if genotype == 'heterozygous':
            self.genotype = self.heteroGenotype()
        elif genotype == 'homozygous_wt':
            self.genotype = [0, 0]
        else:
            self.genotype = [1, 1]

    @staticmethod
    def heteroGenotype():
        """
        generate a heterozygous genotype randomly
        :return: a heterozygous genotype
        """
        first_allele = random.randint(0, 1)
        second_allele = 0 if first_allele is 1 else 1
        return [first_allele, second_allele]


class DiploidPopulation:
    def __init__(self, *args):
        """
        create a diploid population
        :param args: Optionally enter 4 ints that will construct a population with a given composition. must be in
                     order homozygous_wt, homozygous_mutant, heterozygous. If none are entered, or less than 4,
                     a population will be instantiated with 99 homogenous_wt and 1 heterzygous (randomly chosen order of alleles)
        """
        # create the various types of individuals in a diploid population
        self._homozygous_wt = DiploidIndividual('homozygous_wt')
        self._homozygous_mutant = DiploidIndividual('homozygous_mutant')
        self._heterozygous = DiploidIndividual('heterozygous')

        # set private attributes _template_generation and _initial_population
        self._template_generation = {'homozygous_wt': [self._homozygous_wt, 0],
                                     'homozygous_mutant': [self._homozygous_mutant, 0],
                                     'heterozygous': [self._heterozygous, 0]}

        # composition is a dictionary that gives the individual type (by genotype) and amount in the population
        self.composition = deepcopy(self._template_generation)
        self.composition['homozygous_wt'][1] = args[0]
        self.composition['homozygous_mutant'][1] = args[1]
        self.composition['heterozygous'][1] = args[2]

        # size stores the size of the population
        self.size = sum(args)


class Generation(DiploidPopulation):
    def __init__(self, *args):
        """
        generate a population of Individuals (see class Individual or one of its children)
        :param args: see DiploidPopulation
        """
        if len(args) == 4 and all(isinstance(items, int) for items in args):
            super(Generation, self).__init__(*args)
        else:
            print('...initializing population with 99 homozygous wt individuals and 1 heterozygous individual...')
            super(Generation, self).__init__(99, 0, 1)
        # set private attributes _template_generation and _initial_population
        self._initial_composition = deepcopy(self.composition)
        # initialize with equal fitness scores for all genotypes. This will be overwritten by wrightfisher model
        self.model = {'homozygous_wt': 1, 'homozygous_mutant': 1, 'heterozygous': 1}
        # a counter for the number of matings since initial (start at 0)
        self.number_generations_from_initial = 0

    def matePopulation(self):
        """
        mate the current population.
        update composition_dict and number_generations_from_initial
        """
        # used to translate from genotype to description
        genotype_dict = {(0, 0): 'homozygous_wt', (1, 1): 'homozygous_mutant', (0, 1): 'heterozygous',
                         (1, 0): 'heterozygous'}
        # list of possible genotypes
        genotype_list = ['homozygous_wt', 'homozygous_mutant', 'heterozygous']
        # instantiate list to store weights based on model and population composition
        relative_weights = []
        # make a copy of _template_generation
        new_generation_composition = deepcopy(self._template_generation.copy())

        # create list of relative weights
        for genotype in genotype_list:
            if not self.model[genotype] == 0:
                relative_weights.append(self.composition[genotype][1] * float(self.model[genotype]))
            else:
                relative_weights.append(self.composition[genotype][1])
        # create weights by adding weighted sum of each genotype and dividing by sum
        sum_of_weights = sum(relative_weights)
        weights = [weight / sum_of_weights for weight in relative_weights]
        # generate random sets of parents
        population_generator = self.randomlySelectIndividuals(genotype_list, weights, self.size * 2)

        # randomly select allele from sets of two chromosomes
        while True:
            try:
                # randomly select an allele from each parent
                allele_1 = self.composition[next(population_generator)][0].donateAllele()
                allele_2 = self.composition[next(population_generator)][0].donateAllele()
                # create child
                child_genotype = (allele_1, allele_2)
                # update new population count
                child_genotype_name = genotype_dict[child_genotype]
                new_generation_composition[child_genotype_name][1] += 1
            except StopIteration:
                break
        # update composition of generation
        self.composition = new_generation_composition
        # iterate number_of_generations_from_initial
        self.number_generations_from_initial += 1

    @staticmethod
    def randomlySelectIndividuals(population_composition_list, weights, size):
        """ Based on 'roulette wheel random selection
                https://www.tutorialspoint.com/genetic_algorithms/genetic_algorithms_parent_selection.html
                and implemented using python random.choices
        :returns: an iterable list (a generator)
        """
        return iter(random.choices(population=population_composition_list, weights=weights, k=size))

    def printHumanReadableCompositionDict(self):
        """
        print human readable composition dictionary
        """
        for individual_object, number in self.composition_dict.items():
            print("%s : %i" % (individual_object.genotype, number))


class WrightFisherModel(Generation):
    """
    initialize a population with either a dominant or recessive model used in mating
    Attributes:
        model: user chosen model, either dominant or recessive
        max_generations: total number of generations to generate
    """

    def __init__(self, model, fitness_multiplier, *args):
        """
        constructor for WrightFisherModelObject
        :param model: either recessive or dominant
        :param fitness_multiplier: an integer or float. Multiplies penetrance
        :param max_generations: total number of generations to generate
        :param args: see DiploidPopulation
        """
        # pass optional arguments to the parent (Generation) constructor
        super().__init__(*args[2:])
        # store models
        dominant_model = {'homozygous_wt': 0, 'homozygous_mutant': 1, 'heterozygous': 1}
        recessive_model = {'homozygous_wt': 0, 'homozygous_mutant': 1, 'heterozygous': 0}
        # set instance model based on user input
        model_dict = {'dominant': dominant_model, 'recessive': recessive_model}
        model = model.lower()
        if model not in model_dict.keys():
            sys.exit('the model is not a recognized model in object WrightFisherModel. Check input.')
        self.model = model_dict[model]
        self.fitness_multiplier = fitness_multiplier
        # multiply model by fitness_multiplier
        for genotype in self.model:
            self.model[genotype] *= fitness_multiplier

    def resetPopulation(self):
        """
        reset population to initial state. nothing stored from previous generation(s)
        """
        self.composition = deepcopy(self._initial_composition)
        self.number_generations_from_initial = 0

    def mutantAlleleFrequency(self):
        """
        :return: number of mutant alleles in current generation
        """
        return self.composition['homozygous_mutant'][1] * 2 + self.composition['heterozygous'][1]

    def wildtypeAlleleFrequency(self):
        """
        :return: number of wildtype alleles in current generation
        """
        return self.composition['homozygous_wt'][1] * 2 + self.composition['heterozygous'][1]

    def homozygousWildtypeFrequency(self):
        return self.composition['homozygous_wt'][1]

    def homozygousMutantFrequency(self):
        return self.composition['homozygous_mutant'][1]

    def heterozygousFrequency(self):
        return self.composition['heterozygous'][1]

    def determineMutantFixation(self):
        """
        return true if the homozygous_mutant represents 100% of the population
        :return: boolean representing whether the mutant allele fixed
        """
        return self.composition['homozygous_mutant'][1] == self.size

    def determineWildtypeFixation(self):
        """
        return true if the homozygous_wildtype represents 100% of the population
        :return: boolean representing whether the mutant allele fixed
        """
        return self.composition['homozygous_wt'][1] == self.size


if __name__ == "__main__":
    main(sys.argv)
