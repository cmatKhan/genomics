USAGE:
{Please provide the exact command line arguments you used to run map_sequence_starter.py and generate your results}


-
Question 1:
{Output from running map_sequence_starter.py on the Scott cDNAs and McKinley sequencing reads}

Name	                       Reads	Reads per BP      # I fudged this and tabbed the column headings over rather than extract the gene names. My code currently does not extract the gene name
CCDS28177.1|Gap43|Mus musculus	36	0.05263157894736842
CCDS28323.1|Olig2|Mus musculus	6	0.006172839506172839
CCDS29373.1|Mbp|Mus musculus	304	0.4037184594953519
CCDS15640.1|Cd34|Mus musculus	7	0.006092254134029591
CCDS17683.1|Atp1A1|Mus musculus	369	0.1201171875
CCDS20195.1|Abcg2|Mus musculus	7	0.0035460992907801418
CCDS21277.1|Myod1|Mus musculus	0	0.0
CCDS22757.1|Tubb3|Mus musculus	153	0.1130820399113082
CCDS25507.1|Gfap|Mus musculus	37	0.028615622583139984
CCDS26920.1|ChAT|Mus musculus	0	0.0

gene CCDS25507.1|Gfap|Mus musculus is length 0.6550151975683891
gene CCDS22757.1|Tubb3|Mus musculus is length 0.6854103343465046
gene CCDS28323.1|Olig2|Mus musculus is length 0.49240121580547114
gene CCDS26920.1|ChAT|Mus musculus is length 0.9817629179331308
gene CCDS28177.1|Gap43|Mus musculus is length 0.3465045592705167
gene CCDS29373.1|Mbp|Mus musculus is length 0.38145896656534956
gene CCDS17683.1|Atp1A1|Mus musculus is length 1.5562310030395137
gene CCDS21277.1|Myod1|Mus musculus is length 0.4848024316109423
gene CCDS15640.1|Cd34|Mus musculus is length 0.5820668693009119
gene CCDS20195.1|Abcg2|Mus musculus is length 1.0

{Question 1.2 Genes highly expressed, their function, and why?}

{Question 1.3 Genes lowly or not expressed, their function, and why?}
-
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
    two tailed. The brain specific genes may be both over and under enriched.

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
 -
Suggestions:
More math! The quantitative portion of the lectures/lab-lectures has been pretty light so far.