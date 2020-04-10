Part 1:
Question 1:
{How many ORFs were predicted using the described method?}
- 73

- The number of items in both call_orfs_aa and mgm_aa: 8
  The number of items in exclusively in call_orfs_aa: 65
  The number of items in exclusively in mgm_aa: 11

  The number of items in both call_orfs_nt and mgm_nucleotide: 8
  The number of items in exclusively in call_orfs_nt: 65
  The number of items in exclusively in mgm_nucleotide: 11
Part 2:
Question 2:
{Explain what each of the flags for the script gmhmmp are doing. Would you have added or dropped any flags for your particular problem or use?}
- -a output protein sequences of predicted genes
  -d output nucleotide sequences of predicted genes
  -f output format G == gff
  -m file with gene finding parameters (config file)

  Adding the flags -A and -D output fasta files with protein sequences and nucleotide sequences respectively

Question 3:
{Number of ORFs from MetaGeneMark}
- 19 (I can't remember which predicted 20 and which 19, but the perl scripts predicted 1 different than using -A and -D with gmhmmp)

Part 3:
Question 4:
{Number of ORFs shared between MetaGeneMark and call_orfs.py output}
- 8
{Number of ORFs unique to call_orfs.py output}
- 65
{Number of ORFs unique to MetaGeneMark output}
- 11

{Interpretation}
- The search method in call_orfs is naive. Presumably MetaGeneMark uses known frequencies to better predict likely ORFs

Part 4:
Question 5:
{Explain briefly (in a few words) what each one of the parameters mean}
- -db  the database to search against
  -out output file
  -outfmt instructions to blastp on which columns to include in the output and in what order. In this case seq_id, %_identity,
                                                                                              match_length, number_mismatches, number_gaps,
                                                                                              start_alignment_index, end_alignment_index,
                                                                                              e_value, bitscore, sequence_length, sequence_title
Question 6:
{How many unique predicted antibiotic resistance genes were identified using BLAST against the CARD database? How many survive the filtering in your Python script?}
- wc -l blast_to_card.txt
  526 (unfiltered)

The total number of antibiotic genes with identity greater than 80 percent
over greater than 85 percent of the target sequence is: 3
-
Part 5:
Question 7:
{two uses of HMMER}
  Compare protein domains by comparing the alignment sof multiple sequence alignments as opposed to single sequence alignment searches for homology against a database
  (something blastp does not do)
  However, it can also be used for single sequence alignment against a database (in place of blastp)
{and how you would go about executing it}
- hmmsearch -- first, build a multi sequence alignment profile, then search that profile against a database of protein domain profiles
           this can also be done vice versa with hmmscan
    example command: hmmbuild <your multiple sequences>
                     hmmsearch <hmmbuild output> <protein domain database>
    for sequence alignment (multiple sequences):
                     jackhmmer <your sequences> <your database>

- In the case of this assignment, we could identify an antibiotic resistant sequences, build a profile and search for conserved domains
- Or, if we have a sequence that is not annotated, we might use HMMER (instead of blast) to search for homology.
     If you want more detail, I have a somewhat length jupyter notebook script which parses through some hmmer output.

Question 8:
{Total number of resistance genes annotated by Resfams}
- 7
{Comparison of the number of genes annotated by Resfams and the number of genes annotated by BLAST}
- Quite different, prior to filtering -- the HMM method seems to return much less.
Comments:
{Things that went wrong or you can not figure out}
-
