Assignment 5 Due Friday, February 21 at 10am

Part 1.0
{Command for running analyze_WGBS_methylation.py}
- ./analyze_WGBS_methylation.py​ -b data/BGM_WGBS.bed -o .

Question 1:
{What does DNA methylation look like across chromosome 21?}
- Bimodal with high likelihoods at 0 and 1.
  There are twice as many unmethylated CpGs than methylated. This seems to imply that the majority of CpGs on Chromosome 21 are unmethylated.
  I'm not sure if it should, but this was surprising to me. I'll be curious to know if this is close to the expected result.
  I was expecting more methylated (and repressed) regions than unmethylated. Possibly this means that regions of methylation regulation
  are more targeted. Alternatively, maybe there are a good number of CG dinucleotides in non regulatory regions and therefor their methylation
  status is inconsequential. I suspect this is not the case.

Question 2:
{What does the CpG coverage look like across chromosome 21?}
- CpG coverage is count data. Additionally, the mean and the standard deviation are both very close to equal at 26.
  This leads me to believe that a good model for CpG coverage across chromosome 21 is a poisson distribution.

Question 2.1:
{What fraction of the CpGs have 0X coverage?}
- The fraction of 0x CpG sites is: 0.08414904690309218

Part 1.1
{Command for creating a bed file with the average CpG methylation level in each CGI.}

bedtools intersect -a CGI.bed -b BGM_WGBS_CpG_methylation.bed -loj -sorted > data/intersect_CGI_BGM_WGBS_CpG.bed

bedtools groupby -i data/intersect_CGI_BGM_WGBS_CpG.bed -g 1,2,3,4 -c 9 -o mean > WGBS_CGI_methylation.bed

Part 1.2
{Command for plotting the distribution of average CGI methylation levels}
- analyze_CGI_methylation.py -b WGBS_CGI_methylation.bed -o .

Question 3:
{What does DNA methylation look like for CpGs in CGIs? How does it compare to all the CpGs on chromosome 21?}
- Bimodal with high likelihoods at 0 and 1.
  The DNA methylation of the CpGs is (maybe, depending on if I am even in the ballpark) reassuringly similar to that
  of the CpGs across chromosome 21.

Part 1.3.0
Gene promoters
{Command for generating the promoter bed file}
generate_promoters.py​​ -b data/refGene.bed -o . -bn refGene_promoter.bed

{Justification for promoter definition}

30 to 1000bp upstream of the exon closest to the 5' end

"Our observations that 68% of 40-bp core promoter fragments maintain basal promoter activity and that these fragments
contain much of the constraint observed in promoters emphasize the importance of the core promoter. However, the
deletion analyses we report also demonstrate that additional regulatory sequences are present throughout the extended promoter.
Successive removal of sequences in the -350- to -40-bp region of the promoters significantly reduces promoter activity
in the transient transfection assay, indicating that these regions contain positive regulatory elements.
In contrast, the region upstream of -350 tends to contain elements that negatively affect transcription initiation.
This trend was particularly striking within a few of the -1000- to -500-bp regions."

Cooper SJ, Trinklein ND, Anton ED, Nguyen L, Myers RM.
Comprehensive analysis of transcriptional promoter structure and function in 1% of the human genome.
Genome Res. 2006;16(1):1–10. doi:10.1101/gr.4222606

{Commands for generating promoter-CGI and non-promoter-CGI bed files}
- bedtools intersect -a refGene_promoter.bed -b data/CGI.bed > promoter_CGI.bed

- bedtools intersect -a refGene_promoter.bed -b data/CGI.bed -v > non_promoter_CGI.bed


{Justification for overlapping criteria}
- This is an inner join between a and b. I have tried this also with a left outer join, which retains all columns of a and enters
  NULL if there is no entry for b, but decided this was unnecessary.

 - -v reports rows in a that are not found in b

{Commands for calculating the average CpG methylation for each promoter-CGI and non-promoter-CGI}
bedtools intersect -a promoter_CGI.bed -b BGM_WGBS_CpG_methylation.bed​ -loj > data/intersect_promoter.bed
bedtools groupby -i data/intersect_promoter.bed -g 1,2,3,4 -c 9 -o mean > average_promoter_CGI_methylation.bed

bedtools intersect -a non_promoter_CGI.bed -b BGM_WGBS_CpG_methylation.bed​ -loj > data/intersect_non_promoter.bed
bedtools groupby -i data/intersect_non_promoter.bed -g 1,2,3,4 -c 9 -o mean > average_non_promoter_CGI_methylation.bed

{Commands for running analyze_CGI_methylation.py on average_promoter_CGI_methylation.bed and average_non_promoter_CGI_methylation.bed}
analyze_CGI_methylation.py -b ./average_promoter_CGI_methylation.bed -o .
analyze_CGI_methylation.py -b ./average_non_promoter_CGI_methylation.bed -o .

-
Question 4:
{How do the DNA methylation profiles of promoter-CGIs and non-promoter-CGIs differ?}
- The promoter CGIs are overwhelmingly unmethylated. The non-promoter CGIs are more uniformly distributed in their methylation.
  This may point to a mistake in my promoter defintion or calculations.

Part 1.3.1
{Commands for calculating CpG frequency for each promoter type}

bedtools getfasta -fi data/hg19_chr21.fa -bed promoter_CGI.bed > promoter_cgi.fasta
bedtools getfasta -fi hg19_chr21.fa -bed non_promoter_CGI.bed > non_promoter_cgi.fasta

./nuc_count_multisequence_fasta.py promoter_cgi.fasta
./nuc_count_multisequence_fasta.py non_promoter_cgi.fasta


{CpG frequencies for each promoter type}

promoter_cgi
CG:0.10756197050517728

non_promoter_cgi
CG:0.015648984051825063


Question 5:
{What is a possible biological explanation for the difference in CpG frequencies?  Interpret your results from parts 1.3.0 and 1.3.1: what are the “simple rules” for describing regulation by DNA methylation in promoters?}
- It is not surprising (in fact, reassuring that I'm doing something a little right) that there are higher rates of CpGs
  in promoter regions as they are functional units of regulation. Simple rules:
     1) promoter regions are (quite significantly, it seems) enriched for CpG
     2) Methylated CpG is more common in promoter regions

Part 2
{Commands to calculate CGI RPKM methylation scores}
perl bed_reads_RPKM.pl data/CGI.bed data/BGM_MeDIP.bed > MeDIP_CGI_rpkm.bed
perl bed_reads_RPKM.pl data/CGI.bed data/BGM_MRE.bed > MRE_CGI_rpkm.bed
perl bed_reads_RPKM.pl data/CGI.bed data/BGM_WGBS.bed > WGBS_CGI_rpkm.bed


### PLEASE NOTE: I was getting and error about cmd line input file length, so shortened the filenames

{Command to generate the correlation plots}
./compare_methylome_technologies.py -md MeDIP_CGI_rpkm.bed -mr MRE_CGI_rpkm.bed -bs WGBS_CGI_rpkm.bed -o .

{Correlations for each comparison}
(same -- see above. the correlations print to stdout)
{Justification for chosen correlation metric}
Pearson's correlation coefficient. I thought that there would be a positive correlation between counts in the given technologies.
However, either I am generating the data incorrectly, this is not the correct comparison statistic, or they don't correlate well.

-
Question 6:
{How do MeDIP-seq and methylation correlate? How do MRE-seq and methylation correlate? How do MeDIP-seq and MRE-seq correlate?}
- The below numbers and the graphs look like I'm be doing something wrong. Else, the correlation is low. The MRE vs WGBS seems to
  be a bit correlatd, but even that is weak.


The pearson R coefficient for MeDIP_Seq_RPKM_vs_MRE_Seq_RPKM
is (0.03108172430788713, 0.5539030120795345)
The pearson R coefficient for MeDIP_Seq_RPKM_vs_WGBS_Seq_RPKM
is (-0.25507994013728924, 7.872154733877305e-07)
The pearson R coefficient for MRE_Seq_RPKM_vs_WGBS_Seq_RPKM
is (0.39799148135390455, 2.6464007007930845e-15)

The pearson R coefficient for No_Outlier_MeDIP_Seq_RPKM_vs_MRE_Seq_RPKM_
is (-0.6046168054870494, 1.178327709101451e-37)
The pearson R coefficient for No_Outlier_MeDIP_Seq_RPKM_vs_WGBS_Seq_RPKM
is (-0.38917333749715266, 1.3054713382299519e-14)
The pearson R coefficient for No_Outlier_MRE_Seq_RPKM_vs_WGBS_Seq_RPKM
is (0.4036751575938664, 1.0665686790987536e-15)

Outliers
{Answers to outlier questions}
Outlier in the MRE seq:
chr21 	9825442 	9826296
This seems to be at the beginning of the chromosome.

{If applicable: correlations for each comparison}
See above
{If applicable: copy the updated figures to your submissions directory}
- Plots are in the submission directory

Comments:
{Things that went wrong or you can not figure out}
- This is a long assignment. Another weekend would have been very nice and allowed for more focus, checking, etc.

Suggestions:
{What programming and/or genomics topics should the TAs cover in the next class that would have made this assignment go smoother?}
- It would be great to talk a bit about generating test data and writing various types of tests and checks.

- I'd also really like to get a bit more hands on instruction on creating the sort of graphs and visualizations that are common in genomics, especially visualizations
  that aren't of the standard variety. This is not b/c the basic plots aren't important, but due to the fact that there is so much documentation on creating
  them already, whereas there is comparably less for the more specialized graphs.

- However, all that said, I'm writing this because there is space for it. I really appreciate these assignments. Thank you all for your time!

