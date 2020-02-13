Part 1:
Question 1:
{Explanation}
- the nr database is a non-redundant -- similar sequences have been merged -- collection of GenBank and EMBL sequences.

We use the nr database as a comprehensive -- to the extent of our current knowledge -- search against all known protein sequences.

citation: http://arep.med.harvard.edu/seqanal/db.html
  
Question 2:
{Number of hits}
- 211
Question 3:
{Genus species}
Naumovozyma castellii
Score, % identity}
(with BLOSUM62 and default gap/extension)
- Max Score: 820 	Total Score: 820 	Query Coverage: 99% 	E-Val: 0.0	% identity: 51.87%
Question 4:
{Number of hits, explanation}
- Less by one. BLOSUM80's scores are calibrated to sequences with 80% homology and therefore more closely related in evolutionary time.

Question 5:
{Genus species}
- Torulaspora delbrueckii
Question 6:
{Score, % identity}
N. castellii with BLOSUM80:
- Max  Score: 746 	Total Score: 881 	Query Coverage: 85% 	E-val: 0.0	% identity: 69.14%
Question 7:
{Number of hits, explanation}
213. The increase in the extension penalty caused the number of expected hits to decrease, increasing our number of hits unlikely due to chance alone.

Question 8:
{Explanation}
N. castelli with BLOSUM62 and default gap/extension:
- Max Score: 820 	Total Score: 820 	Query Coverage: 99% 	E-Val: 0.0	% identity: 51.87%

N. castelli with BLOSUM62 and 7/2 gap/extension penalty:
- Max Score: 732 	Total Score: 732 	Query Coverage: 99% 	E-val: 0.0	% identity: 51.99%

- The affine gap penalty more severely deducts for the opening of a gap than the extension of a gap. This is done
  because "in nature, a series of k indels often comes as a single event rather than a series of k [single nucleotide indels]"*.

  By increasing the gap extension penalty from 1 to 2, we increased the number of short gaps in our final 'best' alignment.
  Since the gap opening penalty is still high at (minus) 7, this decreases our final max score from 820 to 732 even though the query
  coverage and % identity remain the same (actually, % identity increases very slightly).

  *cite: http://www.csbio.unc.edu/mcmillan/Comp555S16/Lecture14.html. This was also in our slides, but I found this set
  informative, also, and used the quoted sentence above more or less directly from slide 6

Question 9:
{Answer, Explanation}
- More time -- there will be more word fragments to compare
Question 10:
{}
- Quite.

Part 2:
Question 11
{Your command with specific file names}
bowtie2 -x chr22_idx/chr22_idx -U reads.fq -S reads_chr22.sam 2> alignment_report.txt

{Number of uniquely mapped reads}
    7115 (26.69%) aligned exactly 1 time

{Number of multi mapped reads}
    10126 (37.98%) aligned >1 times

{Number of unmapped reads}
    9420 (35.33%) aligned 0 times

-
Question 12
{What is enriched in this dataset} -- ratio of reads.fq / chr22.fq

dinucleotide CG
{Single nucleotide enrichment scores} <-- nothing very noteworthy
read_A / chr22_A : 1.0362139620338378
read_C / chr22_C : 1.0854251812551212
read_G / chr22_G : 1.1137595677932663
read_T / chr22_T : 0.7868097786252916


{Dinucleotide enrichment scores}

read_AA / chr22_AA : 1.1548495676422013
read_AC / chr22_AC : 0.9470807240212612
read_AG / chr22_AG : 0.7355022737268961
read_AT / chr22_AT : 1.3319437493655997
read_CA / chr22_CA : 0.8912113345409589
read_CC / chr22_CC : 1.0508820804164807
read_CG / chr22_CG : 3.8718498632509712   <-- this is the enriched dinucleotide
read_CT / chr22_CT : 0.7093604002470296
read_GA / chr22_GA : 1.313616652067337
read_GC / chr22_GC : 1.135835963884709
read_GG / chr22_GG : 1.1920619310024492
read_GT / chr22_GT : 0.7299651374274183
read_TA / chr22_TA : 0.6944987899479796
read_TC / chr22_TC : 1.1834792349695444
read_TG / chr22_TG : 0.8377493920020211
read_TT / chr22_TT : 0.4830063004242911

{Description of how enrichment scores were calculated}
{Assay, Explanation}
the enrichment for cg indicates that this assay was enriching for GC and was likely exploring gene regulation and methylation.
The assay that produced this data may have been whole genome bisulphate sequencing.

Extra Credit 1:
{Why BLASTn?}
- canis familiaris -- our best friends!
- Blastn translates protein sequences to nucleotide sequences prior to searching for homologous sequences while blastx
  translates nucleotide sequences into protein sequences prior to searching for homology.
Comments:
{Things that went wrong or you can not figure out}
- I found number 8 to be the most challenging question -- I am going to visit office hours again next week to make sure that I understand. Thank you for all your help!
Suggestions:
{What programming and/or genomics topics should the TAs cover in the next class that would have made this assignment go smoother?}
-