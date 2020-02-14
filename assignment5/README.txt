Assignment 5 Due Friday, February 21 at 10am

Part 1.0
{Command for running analyze_WGBS_methylation.py}
{Copy your output files BGM_WGBS_CpG_methylation.bed, BGM_WGBS_methylation_distribution.png, and BGM_WGBS_CpG_coverage_distribution.png to your submissions directory}
-
Question 1:
{What does DNA methylation look like across chromosome 21?}
-
Question 2:
{What does the CpG coverage look like across chromosome 21?}
-
Question 2.1:
{What fraction of the CpGs have 0X coverage?}
-
Part 1.1
{Command for creating a bed file with the average CpG methylation level in each CGI.}
{Copy WGBS_CGI_methylation.bed to your submissions directory}
-
Part 1.2
{Command for plotting the distribution of average CGI methylation levels}
{Copy analyze_CGI_methylation.py and WGBS_CGI_methylation_distribution.png to your submissions directory}
-
Question 3:
{What does DNA methylation look like for CpGs in CGIs? How does it compare to all the CpGs on chromosome 21?}
-
Part 1.3.0
Gene promoters
{Command for generating the promoter bed file}
{Justification for promoter definition}
{Copy generate_promoters.py and refGene_promoters.bed to your submissions directory}
-
Promoter-CGI and non-promoter-CGI
{Commands for generating promoter-CGI and non-promoter-CGI bed files}
{Justification for overlapping criteria}
{Commands for calculating the average CpG methylation for each promoter-CGI and non-promoter-CGI}
{Commands for running analyze_CGI_methylation.py on average_promoter_CGI_methylation.bed and average_non_promoter_CGI_methylation.bed}
{Copy refGene_promoters.bed, promoter_CGI.bed, non_promoter_CGI.bed average_promoter_CGI_methylation.bed, average_non_promoter_CGI_methylation.bed, average_promoter_CGI_methylation.png and average_non_promoter_CGI_methylation.png to your submissions directory}
-
Question 4:
{How do the DNA methylation profiles of promoter-CGIs and non-promoter-CGIs differ?}
-
Part 1.3.1
{Commands for calculating CpG frequency for each promoter type}
{CpG frequencies for each promoter type}
-
Question 5:
{What is a possible biological explanation for the difference in CpG frequencies?  Interpret your results from parts 1.3.0 and 1.3.1: what are the “simple rules” for describing regulation by DNA methylation in promoters?}
-
Part 2
{Commands to calculate CGI RPKM methylation scores}
{Command to generate the correlation plots}
{Correlations for each comparison}
{Justification for chosen correlation metric}
{Copy compare_methylome_technologies.py, MeDIP_CGI_RPKM.bed, MRE_CGI_RPKM.bed MeDIP_CGI_RPKM_vs_MRE_CGI_RPKM.png, MeDIP_CGI_RPKM_vs_WGBS_CGI_methylation.png, and MRE_CGI_RPKM_vs_WGBS_CGI_methylation.png to your submissions directory}
-
Question 6:
{How do MeDIP-seq and methylation correlate? How do MRE-seq and methylation correlate? How do MeDIP-seq and MRE-seq correlate?}
-
Outliers
{Answers to outlier questions}
{If applicable: correlations for each comparison}
{If applicable: copy the updated figures to your submissions directory}
-
Comments:
{Things that went wrong or you can not figure out}
-
Suggestions:
{What programming and/or genomics topics should the TAs cover in the next class that would have made this assignment go smoother?}
-
Extra credit
{Commands for running bed_reads_RKPM.pl}
{Command for running analyze_H3K4me3_scores.py}
{Copy analyze_H3K4me3_scores.py, H3K4me3_RPKM_promoter_CGI.bed, H3K4me3_RPKM_non_promoter_CGI.bed, and H3K4me3_RPKM_promoter_CGI_and_H3K4me3_RPKM_non_promoter_CGI.png}
-
Question EC.1:
{How does the H3K4me3 signal differ in promoter-CGIs and non-promoter-CGIs?}
-
Question EC.2:
{What are some better alternatives to model MeDIP-seq data and MRE-seq data instead of using RPKM? Explain.}
-
Question EC.3:
{What would be a better way to compare H3K4me3 values instead of using boxplots? Explain.}
-


