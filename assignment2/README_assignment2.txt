﻿name: chase mateusiak

USAGE:
{Please provide the exact command line arguments you used to run map_sequence_starter.py and generate your results}

/assignment2$ ./map_sequence_starter.py /home/assignments/assignment2/scott_mouse_cDNAs.fa /home/assignments/assignment2/mckinley_raw_reads.txt

-
Question 1:
{Output from running map_sequence_starter.py on the Scott cDNAs and McKinley sequencing reads}

Read in the cDNAs
Created dictionary
Name	Reads	Reads per BP
Abcg2	7	0.0035460992907801418
Atp1A1	369	0.1201171875
Cd34	7	0.006092254134029591
ChAT	0	0.0
Gap43	36	0.05263157894736842
Gfap	37	0.028615622583139984
Mbp	304	0.4037184594953519
Myod1	0	0.0
Olig2	6	0.006172839506172839
Tubb3	153	0.1130820399113082


############ question 2 numbers ##################

The probability that 10 of 20 highly enriched genes is a gene of brain function is 6.442043387911576e-06
p-val for this distribution is 7.150904021083785e-06

-- It is entirely possible that I am missing some critical information. However, without information about
the experimental conditions and a control, it is hard (for me at least -- I am not familiar with these
genes) to know if these numbers are outside of the range of what we'd expect.

{Question 1.2 Genes highly expressed, their function, and why?} -- hard to tell w/out control

Atp1A1 -- membrane protein. Na+K+ atpase pump. Associated with Charcot-Marie-tooth disease (axon degeneration disease)
Mbp    -- Involved in myelination of the nerves
Tubb3  -- highly enriched in neuronal cells. Microtubual binding protein

 -- This is a brain sample, and these are brain/neuron associated proteins. Presumably these are results from a single
 experimental group, which would lead me to believe that there was some sort of perturbation causing these genes to act unusually.
 If we assume that these are (significantly) highly expressed, then the neurons are expressing a protein pump in an unusual way.
 This is implicated in Charcot-Marie-Tooth disease, so possibly this is a marker of axon degeneration. I would be curious to know
 if the experimenters noted any coordination issues in their mice. Highly expressed Mbp seems to indicate increase myelination, another
 phenotypic observation that may be possible in these cells. Tubb3 may be overexpressing some structural components (although the
 information I found on Tubb3 says it is highly enriched in neuronal cells, so this may be a normal expression level).

{Question 1.3 Genes lowly or not expressed, their function, and why?}
ABCG2 -- extra and intracellular membrane transport. Implicated as a breast cancer resistance protein.
Cd34  -- transmembrane phospho-glyco protein. Expressed on blood vessels
Gfap  -- cytoskeletal protein
Myod1 -- muscle differentiation
Olig2 -- expressed in CNS. Also involved in motor neuron differentiation. Also oligodendrocyte differentiation

 -- These seem to have to do with transport and cell differentiation. It may indicate the stage of cell development,
 though it is also possible that the cells are diverting resources to other more critical needs.

Question 2:
{Yes or No} -- I'm not sure what this refers to
{Statistical test used and explanation}

    The binomial distribution. We can consider a random variable to represent either 1 for a brain function gene or 0
    for a non-brain function gene. Since this is a bernoulli trial over n iid experiments, we can use the binomial distribution
    to model the probability of a certain number of successes (drawing a brain function gene) over a certain number of
    trials.

    The probability is p(10; 20, 3/30) = .00000644

    The p-value is the probability of a result more extreme

    p = .0000071509, which passes our typical p < .05 to conclude that this result is statistically significant

{one-tailed or two-tailed}
    two tailed. The brain specific genes may be either over or under enriched.

{Show Me the, ehh, Calculation}
    (20 -C- 10)(3/30)^10(1-(3/30)^10 -- 20 choose 10, times the probability of a gene being annotated for brain function,
    times the probability of a gene not being annotated for brain function, both raised to the appropriate power.

{p-value}
    p-val of this probability is 7.150904021083785e-06. P(X=10) is therefore statistically significant at the .05 cutoff

Question 3:
To normalize the counts by gene length. Otherwise, we would not be able to meaningfully compare the number of reads mapping
to each transcript as the raw count depends on the gene length (the longer the gene, the more likely it is that there are more transcripts)
-
Question 4:

{Limitation 1}
    indels, sequencing errors, and variation are not accounted for (requires a strict match)
{Limitation 2}
    If gene sequences are have repeated 25 length sequences, the first gene with this sequence will receive all of the counts
{Limitation 3}
    Requires that all reads be fragments of one and only one gene
-
Comments:
{Things that went wrong or you can not figure out}
Re-phrase the question on highly/lowly expressed genes to something like, "Choose 8 genes and describe their function".
 -
Suggestions:
More math!
