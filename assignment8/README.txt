Part I


Question 1:
{Exact command to run count_gv.py}
/count_gv.py -snv git staSNV_indel.biallelic.vcf -sv sv.reclassed.filtered.vcf
{Output count table}
- {'SNV': 4086192, 'INDEL': 479363, 'DEL': 1664, 'DUP': 530, 'INV': 95, 'MEI': 1200, 'BND': 2690}

Question 2:
{Proportion of genomic variants are SVs} 
- ~.135%
Question 3:
{Describe the spectrum of SVs for the individual NA12878}
- Inversions are quite uncommon, deletions, mobile element insertions and
- 'Indeterminant Break Ends' (unsurprisingly -- this seems to be a catchall category) are more common
Question 4:
{Describe the distributions observed in each of the histograms}
- The distribution of structural variant deletions is skewed right, meaning most deletions are in the shorter end
- of the range. However, there is a long tail, so longer deletions do exist in NA12878. The majority of SNV are
- small while MEIs seem to have a much tighter range than the other two. The Majority of MEIs seem to be much more similar
- in size.
{Speculate how the length distribution might differ if we limit the data to exonic indels?}
- I would hypothesize that exonic indels are, on average, shorter. Therefore, I would expect that there would be a more curtailed tail.


Part II


Question 5:
{Exact command to run quantify_genotype.py}
quantify_genotype.py -snv SNV_indel.biallelic.vcf
{Output count table}
- {'homozygous_ref': 3518113, 'homozygous_alt': 1643132, 'heterozygous': 2922423, '1+ missing allele': 95900}

Question 6:
{Does the difference in the number of homozygous alternate (or non-reference homozygous) and heterozygous SNVs and indels, make biological sense? Why, or why not}
- {'homozygous_ref': 3518113, 'homozygous_alt': 1643132, 'heterozygous': 2922423, '1+ missing allele': 95900}

Question 7:
{How many variants clearly violate the rules of Mendelian segregation? (Autosome only)}
- {'0/1': 103221, './.': 69046, '0/0': 73456, '1/1': 40791}
The total number of violations is 286514

Question 8:
{Describe four potential reasons that could explain the Mendelian violations.}
- sequencing error, de novo mutation (somatic or germ line), unexpected paternity (uh oh!),
"INDEL detection is difficult and error prone" (from slides), so possibly alignment error

Question 9:
{How many variants now violate mendelian segregation after filtering? (Autosome only)}
{'0/1': 26275, '1/1': 4728, '0/0': 4314}
The total number of violations is 35317

-
Comments:
{Things that went wrong or you can not figure out}
-
Suggestions:
{What programming and/or genomics topics should the TAs cover in the next class that would have made this assignment go smoother?}
