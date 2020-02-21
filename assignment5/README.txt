Assignment 5 Due Friday, February 21 at 10am

Part 1.0
{Command for running analyze_WGBS_methylation.py}
- analyze_WGBS_methylation.py -b BGM_WGBS.bed -o .

Question 1:
{What does DNA methylation look like across chromosome 21?}
- Bimodal with high likelihoods at 0 and 1.
  There are twice as many unmethylated CpGs than methylated. This seems to imply that the majority of CpGs on Chromosome 21 are unmethylated.
  I'm not sure if it should, but this was surprising to me. I'll be curious to know if this is close to the expected result.
  If I were to try to interpret this, I would guess that there are many CpGs that are not in regulatory regions, and therefore less
  likely to be methylated. In order to test this, I would look at the rates of methylation in known regulatory regions vs those CpGs in
  areas that we may know are not regulatory.

Question 2:
{What does the CpG coverage look like across chromosome 21?}
- CpG coverage is count data. Additionally, the mean and the standard deviation are both very close to equal at 26.
  This leads me to believe that a good model for CpG coverage across chromosome 21 is a poisson distribution.

Question 2.1:
{What fraction of the CpGs have 0X coverage?}
- The fraction of 0x CpG sites is: 0.08414904690309218

Part 1.1
{Command for creating a bed file with the average CpG methylation level in each CGI.}

bedtools intersect -a CGI.bed -b BGM_WGBS_CpG_methylation.bed -loj -sorted > intersect_CGI_BGM_WGBS_CpG.bed

bedtools groupby -i intersect_CGI_BGM_WGBS_CpG.bed -g 1,2,3,4 -c 9 -o mean > WGBS_CGI_methylation.bed

Part 1.2
{Command for plotting the distribution of average CGI methylation levels}
- analyze_CGI_methylation.py -b /home/chase/code/cmatkhan/genomics/assignment5/WGBS_CGI_methylation.bed -o .

Question 3:
{What does DNA methylation look like for CpGs in CGIs? How does it compare to all the CpGs on chromosome 21?}
- Bimodal with high likelihoods at 0 and 1.
  The DNA methylation of the CpGs is (maybe, depending on if I am even in the ballpark) reassuringly similar to that
  of the CpGs across chromosome 21.

Part 1.3.0
Gene promoters
{Command for generating the promoter bed file}


{Justification for promoter definition}

30 to 1000bp upstream

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

{Copy generate_promoters.py and refGene_promoters.bed to your submissions directory}
-
Promoter-CGI and non-promoter-CGI
{Commands for generating promoter-CGI and non-promoter-CGI bed files}
bedtools intersect -a refGene_promoter.bed -b CGI.bed -loj > promoter_CGI.bed

bedtools intersect -a refGene_promoter.bed -b CGI.bed -v > non_promoter_CGI.bed

{Justification for overlapping criteria}
Left Outer Join returns all values in a which overlap with B. If no overlap is found, A is reported with value Null.
This returns all information in our reference which overlaps with CGI.bed, while also providing info on what does not.

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
bedtools getfasta -fi hg19_chr21.fa -bed promoter_CGI.bed > promoter_cgi.fasta
bedtools getfasta -fi hg19_chr21.fa -bed non_promoter_CGI.bed > non_promoter_cgi.fasta


./nuc_count_multisequence_fasta.py promoter_cgi.fasta
promoter_cgi
CG:0.02750335674884078

./nuc_count_multisequence_fasta.py non_promoter_cgi.fasta
non_promoter_cgi
CG:0.015648984051825063

Question 5:
{What is a possible biological explanation for the difference in CpG frequencies?  Interpret your results from parts 1.3.0 and 1.3.1: what are the “simple rules” for describing regulation by DNA methylation in promoters?}
- The higher rate of methylation in promoters signals that methylation is more likely to occur where it has a functional impact on transcription.

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
- It would be great to talk a bit about generating test data. It's something I know I should do before I start working on something and I hardly ever do.

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


