Question 1:
{Assumptions of Wright-Fisher}
- From the slides:
    2 allele system
    N diploid individuals in each generation
    2N gametes in each generation
    Random mating, no selection, no mutation
    discrete generations

Part 1
{Usage statement for dominant model}
- wrightfisher.py -f 3 -g 1000 -m recessive -p 100 # please read the usage statement at the top of wrightfisher for details
Question 2:
{Fitness value}
- 3
{Explanation}
- I expected to need a greater fitness value to overcome the population composition (allele frequency)
Part 2
{Usage statement for recessive model}
- wrightfisher.py -f 3 -g 1000 -m recessive -p 100 # please read the usage statement at the top of wrightfisher for details
Question 3:
{Fitness value}
- None that my simulation was able to handle.
{Explanation}
- From the slides, the strength of effect of drift vs selection is determined by allele frequency and population size. Without modulating those two factors, selection (fitness) is unlikely to overcome drift (especially at 90% in a recessive model of penetrance). The dominant model, on the other handle, allows for far greater effect of the fitness value as any presence of the mutant allele will make that genotype more likely to be selected.

Comments:
{Things that went wrong or you can not figure out}
- This assignment prompt is not clearly written. This was quite frustrating given that we are not able to come in for office hours.
