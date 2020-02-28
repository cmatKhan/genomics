 ----> reverse complement and scoring matrix

./highest_affinity_binding_site.py -m polymerase_score_matrix.txt -o .

- The highest scoring sequence is GCACGGCACG with a score of 62.0

./highest_affinity_binding_site.py -m tf_score_matrix.txt -o .
The highest scoring sequence is GCTGCGCACG with a score of 48.0

Question 1:
{Score for highest affinity binding site for the TF scoring matrix}
 - 48.0
{Sequence(s) corresponding to highest affinity binding site for the TF scoring matrix}
 - GCTGCGCACG
 * note: in position 5, there could be either a G or a C. The score is the same

{Score for highest affinity binding site for the polymerase scoring matrix}
 - 62.0
{Sequence(s) corresponding to highest affinity binding site for the polymerase scoring matrix}
 - GCACGGCACG

{Explain how you arrived at your answer.}
 - By taking the highest score at each position with a script and then doing a visual inspection of the matricies.
   See script 'highest_affinity_binding_site.py'.


{If you wrote a script, provide the command to run it.}

 - ./highest_affinity_binding_site.py -m <your_score_matrix.txt>

Question 2:
{Command to run scan_sequence.py on the TF and promoter1}
 - ./scan_sequence.py tf_score_matrix.txt promoter1.txt 40

{Output of command}

The results for tf_score_matrix.txt and promoter1.txt at threshold 40.0 are:

Forward and reverse with scores greater than or equal to 40.0
orientation     position        score   sequence
forward             68          45.00   GCTATGCACG
reverse             27          43.00   GCTGTGCAGG


{Command to run scan_sequence.py on the TF and promoter2}
 - ./scan_sequence.py tf_score_matrix.txt promoter2.txt 40

{Output of command}

The results for tf_score_matrix.txt and promoter2.txt at threshold 40.0 are:

Forward and reverse with scores greater than or equal to 40.0
orientation     position        score   sequence
forward             181            48.00   GCTGCGCACG


{Command to run scan_sequence.py on the polymerase and promoter1}
- ./scan_sequence.py polymerase_score_matrix.txt promoter1.txt 40

{Output of command}

The results for polymerase_score_matrix.txt and promoter1.txt at threshold 40.0 are:

Forward and reverse with scores greater than or equal to 40.0
orientation     position        score   sequence
forward             68          41.00   GCTATGCACG
forward             102         59.00   CCACGGCACG


{Command to run scan_sequence.py on the polymerase and promoter2}
 - ./scan_sequence.py polymerase_score_matrix.txt promoter2.txt 40

{Output of command}

The results for polymerase_score_matrix.txt and promoter2.txt at threshold 40.0 are:

Forward and reverse with scores greater than or equal to 40.0
orientation     position        score   sequence
forward             181         43.00   GCTGCGCACG
forward             186         62.00   GCACGGCACG
forward             202         42.50   TAACGCCACG
forward             386         41.50   GGAAAGCACG
forward             391         44.00   GCACGGTACC
reverse             28          42.00   GCTGGGCAAG


{Name of the promoter that you'd expect to be repressed by the TF}

Promoter 1

{Name of the promoter that you'd expect to be activated by the TF}

Promoter 2

{Explanation}

Promoter 1:
 - I'd expect the transcription factor at position 68 on promoter1 to block the polymerase from binding. This would
repress transcription.

The results for tf_score_matrix.txt and promoter1.txt at threshold 40.0 are:

Forward and reverse with scores greater than or equal to 40.0
orientation     position        score   sequence
forward             68          45.00   GCTATGCACG

The results for polymerase_score_matrix.txt and promoter1.txt at threshold 40.0 are:

Forward and reverse with scores greater than or equal to 40.0
orientation     position        score   sequence
forward             68          41.00   GCTATGCACG

Promoter 2:
- The polymerase has the highest score at position 186 and strong scores downstream of that while the TF scores a 48
at position 181. I'd hypothesize that the TF is functioning to enable the polymerase at promoter 2.

The results for tf_score_matrix.txt and promoter2.txt at threshold 40.0 are:

Forward and reverse with scores greater than or equal to 40.0
orientation     position        score   sequence
forward             181            48.00   GCTGCGCACG

The results for polymerase_score_matrix.txt and promoter2.txt at threshold 40.0 are:

Forward and reverse with scores greater than or equal to 40.0
orientation     position        score   sequence
forward             181         43.00   GCTGCGCACG
forward             186         62.00   GCACGGCACG
forward             202         42.50   TAACGCCACG
forward             386         41.50   GGAAAGCACG
forward             391         44.00   GCACGGTACC
reverse             28          42.00   GCTGGGCAAG

Comments:
{Things that went wrong or you can not figure out}
-
Suggestions:
{What programming and/or genomics topics should the TAs cover in the next class that would have made this assignment go smoother?}
