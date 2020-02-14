Assignment 4 Due February 14, 2020 at 10am

Please provide the exact command line arguments you used to generate your results.
{How to run gene_expression.py}
./gene_expression.py raw_counts.txt

****PLEASE NOTE****
I work on these scripts on my local machine, which has python3 installed. Python2.7.17 is installed on the class server, however.
Surprisingly, this hasn't been an issue yet and I hadn't noticed, but it causes problems is in this case.
The script will run, but the output of the FLD function is quite different. My guess is it has
something to do with the floating point arithmetic. But, I didn't notice until Friday morning.

I would normally run a virtual environment to get around this sort of issue, but I hope that an otherwise working script
(there is output, it is just the wrong numbers) and a screenshot of the python3 output is enough proof that a virtual
environment would solve the issue. Please see the image "screenshot" for the python3 output.

Please see the R script to see the code I used to generate my graphs.

Question 1:
{How many genes are left after removing genes with zero expression in all samples?}
- 23259
Question 2:
{How many genes are left after removing genes where 20 or more samples have cpm < 1?}
- 13241
Question 3:
{What is the range of library sizes (min, max)?}
- The smallest range of a gene by raw count is 21
  The largest range of a gene by raw count is 6124081
Question 4:
{What is the range of library sizes (min, max) after normalization?}
- The minimum range by gene after normalization is 13755684.053136578
  The maximum range by gene after normalization is 22154262.742857095

Question 5:
{Compare the two library size bar charts you made. How did the distribution of library sizes change after normalization?}
- The graph compresses towards the average library size (samples with a large number of reads are scaled down, samples with a smaller number of reads are scaled up)

{Briefly discuss why it is important to normalize your RNA-seq data.}
- Some libraries are sequenced to a greater depth, and therefore have more reads for all genes, and vice versa.
  In order to compare the expression of a given gene in a given sample to another sample, it is necessary to choose
  a scale by which to account for the difference in library depth.

Question 6:
{What are the top ten differentially expressed genes according to your FLD analysis? (Copy and paste your function's output.)}
- The top 10 differentially expressed genes by FLD are:
           gene       fld
        AK000953    1.018488
        RAB30       0.988243
        DBNDD1      0.957612
        CTDSPL      0.915162
        GRB14       0.885167
        YTHDC2      0.851677
        MYLK4       0.842608
        ABCG1       0.746580
        FNDC5       0.739654
        BC010186    0.737283

{Do these genes make sense given the tissue and groups in the experiment?}
- RAB30 -- increases in diabetic mice (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3516129/)
- DBNDD1 -- implicated in a glycogen storage disease, also upregulated in the muscles of diabetic rats (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3931395/)
- CTDSPL -- downregulated in diabetic neuropathy (hopefully this is not the case for our friends in the study) (https://www.ncbi.nlm.nih.gov/pubmed/26063197)
- GRB14 -- upregulation leads to insulin resistence (https://www.ncbi.nlm.nih.gov/pubmed/15059968)

  Yes -- These genes do make sense given the experiment

Question 7: 
{Does your result point toward one gene with large effect or many genes with small effects?}
- Many with small effects
{Does RNA-seq expression data always give researchers a clear answer?}
- No. Not only is quality assessment of the data hard, but interpretation of good results is also challenging.
Question 8:
{How does the study design of this experiment relate to the assumptions made when studying gene expression data?}
- 20 biological replicates is quite a lot, so if there are actual changes in expression in the before/after exercise group, as long as they are large enough to overcome the biological and technical variability, we would hopefully be able to observe them.
  However, there is also a lot that can influence any individual over a year and it is not necessarily true that most other factors (diet, stress, etc) were held constant. This likely introduces good amount of variability that even 20 biological replicates may not
  overcome. Because we have so many replicates, we likely would be able to examine this variability and possibly adjust for it (e.g. outlier removal or something akin to 'batch effect' if we have other phenotype data that would allow for reasonable grouping).

Question 9:
{If you were going to spend time and money following up on one of these top ten genes, what would be your candidate and why? (There could be many correct answers.)}
- I would be interested in CTDSPL. There is enough interest in diabetes, and the possibility of losing a foot is bad enough, that maybe a grant would be funded.
  RAB30 may also be a good candidate since it is a highly conserved and we could follow up this experiment in a variety of models. It seems to be essential to the Golgi,
  so further study may lead to an interesting result in protein trafficking or other fundamental cell biology questions.

Extra credit 1:
{What do you expect to see?}
-
Extra credit 2:
{What did you actually see? If you did not find what you expected, what sorts of variation could account for this?}
-
Comments:
{Things that went wrong or you cannot figure out}
-
Suggestions:
{What programming and/or genomics topics should the TAs cover in the next class that would have made this assignment go smoother?}
-


