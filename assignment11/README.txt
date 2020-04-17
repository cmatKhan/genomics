Part 1:
{Command line argument you used for filter_variants.py}
- ./filter_variants.py -i data/variant_to_barcode.txt

Question 1:
{How many variants are not well represented by barcodes?}
- 37 barcodes in raw - 31 barcodes in filtered = 6 not well represented by barcodes

Question 2:
{Why do you want to know barcode counts in plasmid DNA?}
- In order to normalize the expression in the cDNA as a proportion of occurrences in the cDNA to the occurrences in the pDNA.

Question 3:
{Command line argument you used for count_barcodes.py}
 - ./count_barcodes.py -b filtered_variant_to_barcode.tsv -f data/pDNA.fq
   ./count_barcodes.py -b filtered_variant_to_barcode.tsv -f data/cDNA.fq

{How many reads in cDNA.fq are left after filtering?}
- The number of reads left after filtering is: 52073

Part 2:
{Command line argument you used for analyze_MPRA.py}
- ./analyze_mpra.py -c cDNA_count.tsv -p pDNA_count.tsv -f filtered_variant_to_barcode.tsv -e data/variant_eQTL_results.txt
Question 4:
  4.1 -                        Count_by_eQTL
Significant_eQTL_gene
ENSG00000120071                   17
ENSG00000159202                    6
ENSG00000184716                    7

  4.2 - The variants associated with ENSG00...120071 have the most significant p-values,
  but their effect is negative while the direction of effect of the eQTL is positive.
  In fact, all of my log2FC are in the opposite direction of the eQTL for the most part, which makes me think
  I've done something wrong.

Question 5:
{Possible reason?}
 - It is very hard to tell if my results are the result of a mistake on my part, or if this is a planned
   part of the assignment. It would be nice if, rather than hiding information like this in the prompt, it would be
   obvious. I don't understand how unclear questions benefit learning, nor how a check on accuracy that would prompt
   a student to re-check their work would be cheapening the learning experience.

   That said, it is hard for me to imagine that variants that are expected to have a positive effect (presumably
   less susceptibility to Parkinsons or relieving some of the symptoms) would across the board have a negative log2FC at
   a significant p-value. That points to a mistake in my code.

   It is also possible that while these variants are beneficial to alleviating traits of Parkinsons, they are otherwise
   detrimental to the cell. This may be an effect of the cell line and culturing.

{Way to improve in the study design?}

 - More careful checks on the code, quite possibly. I don't see my error, but it seems likely that it is there. More tests
   on my functions against data with known results would address this.

   There isn't any information on the timeframe in which the cells were collected, and that may have an effect on
   how many cells are in each group. It may be worth testing different times for collection.

   Also, since we're testing Parkinsons, it seems that it would be more appropriate to use a neuronal cell line, if one exists.
   At the very least, possibly not using a cancer cell line would be more appropriate.

   The docs on scipy for the mann whitney U test imply that there need to be more than items from which you are sampling
   to have an accurate p-value. If I understand correctly, what we are sampling is not cells, but counts. In that case,
   it may be that a different significance test is more appropriate and these p-values are misleading.

