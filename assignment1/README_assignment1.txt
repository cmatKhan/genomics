Name: Chase Mateusiak

USAGE:

For nuc_count:

   $: ./nuc_count.py ~/Downloads/chr20.fna

For make_seq.py

   $: ./make_seq.py 1000000 .25 .20 .30 .25

OUTPUT:

nuc_count output

Question 1 ​Run ​nuc_count.py ​on ​chr20.fna​. How many times do each of the 4 nucleotidesoccur in chr20?

Raw Counts
A:  17867246
C:  13916133
G:  14094472
T:  18066406
N:  499910

Question 2 ​Run your modified ​nuc_count.py ​on ​chr20.fna​. What are the frequencies of the4 nucleotides on chr20?

Sequence length: 64444167
 The frequency of A is 0.27725156258129613
 The frequency of T is 0.2803419896792211
 The frequency of G is 0.21870826571472327
 The frequency of C is 0.21594092449049734

The Frequency of dinucleotides is:
The frequency of dinucleotide sequence AA in chromosome 20 is 0.08797102497291075
The frequency of dinucleotide sequence AC in chromosome 20 is 0.05043343436730618
The frequency of dinucleotide sequence AG in chromosome 20 is 0.07270115704687877
The frequency of dinucleotide sequence AT in chromosome 20 is 0.06831340740878825
The frequency of dinucleotide sequence CA in chromosome 20 is 0.07454682966552946
The frequency of dinucleotide sequence CC in chromosome 20 is 0.058469913941174435
The frequency of dinucleotide sequence CG in chromosome 20 is 0.012096130486833452
The frequency of dinucleotide sequence CT in chromosome 20 is 0.0725161519030164
The frequency of dinucleotide sequence GA in chromosome 20 is 0.061722345030746195
The frequency of dinucleotide sequence GC in chromosome 20 is 0.047939993526852985
The frequency of dinucleotide sequence GG in chromosome 20 is 0.05968420763563614
The frequency of dinucleotide sequence GT in chromosome 20 is 0.05107164737838171
The frequency of dinucleotide sequence TA in chromosome 20 is 0.05517868337891268
The frequency of dinucleotide sequence TC in chromosome 20 is 0.0607859031022188
The frequency of dinucleotide sequence TG in chromosome 20 is 0.07593666712498316
The frequency of dinucleotide sequence TT in chromosome 20 is 0.09063250302983065

Question 3 ​Run the modified ​nuc_count.py ​for both human chr20 and your generated ‘random_seq_1M.txt’
from part 4. Compare the two lists of frequencies. What are the differences? Can you provide a biological
explanation for these differences?

There is less G and C in the human chromosome than my randomly generated chromosome. This may be an error in my "random"
function or sampling variation as the difference is slight. However, it may also indicate that the G-C frequency on
chromosome 20 is less than that of the A-T frequency. I did a very quick search for papers and looked at one that said
the human GC content is ~ 40%, so it is not terribly unlikely that this is a reasonable observation.
(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6391780/#!po=22.7273)

** please note: I did notice that the frequencies of AT and GC nearly identically flipped from the chromosome 20 numbers
** and the random sequence. The script handles and prints by variable and key rather than numerical position in the input
** so even though I did change the order of the frequency output vs the input order, it is pretty easy to check that I
** don't (think) there is a typo there. However, my random sequence generation probably has a bias in it that I
** haven't figured out and fixed.

Similarly, dinucleotide frequencies were more uniform in the random sequence than in chromosome 20. The differences were
very small, but this could indicate that there is a higher frequency of repetition of certain dinucleotides in the actual
genome due to some reason other than chance. This may be accounted for by common genomic features such as start/stop codons,
transcription start sites, and others. Dinucleotides AA and TT seemed to be quite a bit more common (relatively) to the
other dinucleotides in chromosome 20, which points to there being regions in chromosome 20 that have
longer strings of A and T than in the random sequence. In the random sequence, there are more GG and CC dinucleotides.
GG and CC repetition seem to be less likely in the actual sequence. But, as noted above, this is probably an error
in my "random" generation.

make_seq output:
$: ./make_seq.py 1000000 0.28 0.22 0.22 0.28

# note the order of input is unchanged from the template script: seq_length, a_freq, c_freq, g_freq, t_freq

>random_sequence
CCGCCCCTGCGGACTGCATCATCTGCCAATTGCATCTCTTTGCTACAAATCACTGACTTACTAATCTGAACGAGGGGGGCGACCGCGATCGGATACGCAG
GGGGGCCAAGTCTCATGGGC...

Raw Counts
A:  220468
C:  279541
G:  280189
T:  219802
N:  0
Sequence length: 1000000

 The frequency of A is 0.220468

 The frequency of T is 0.219802

 The frequency of G is 0.280189

 The frequency of C is 0.279541


The Frequency of dinucleotides is:
The frequency of dinucleotide sequence AA in chromosome 20 is 0.0488040488040488
The frequency of dinucleotide sequence AC in chromosome 20 is 0.06136606136606137
The frequency of dinucleotide sequence AG in chromosome 20 is 0.0619000619000619
The frequency of dinucleotide sequence AT in chromosome 20 is 0.0483980483980484
The frequency of dinucleotide sequence CA in chromosome 20 is 0.06137406137406137
The frequency of dinucleotide sequence CC in chromosome 20 is 0.07809407809407809
The frequency of dinucleotide sequence CG in chromosome 20 is 0.0784920784920785
The frequency of dinucleotide sequence CT in chromosome 20 is 0.06158006158006158
The frequency of dinucleotide sequence GA in chromosome 20 is 0.06177606177606178
The frequency of dinucleotide sequence GC in chromosome 20 is 0.07867807867807868
The frequency of dinucleotide sequence GG in chromosome 20 is 0.07793907793907794
The frequency of dinucleotide sequence GT in chromosome 20 is 0.0617960617960618
The frequency of dinucleotide sequence TA in chromosome 20 is 0.048514048514048516
The frequency of dinucleotide sequence TC in chromosome 20 is 0.0614020614020614
The frequency of dinucleotide sequence TG in chromosome 20 is 0.06185806185806186
The frequency of dinucleotide sequence TT in chromosome 20 is 0.048028048028048026
