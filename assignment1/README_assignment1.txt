USAGE:

For nuc_count:

   $: ./nuc_count.py ~/Downloads/chr20.fna

For make_seq.py

   $: ./make_seq.py 1000 .28 .28 .21 .23

*** please note: make_seq.py will create a "random_sequence.fa" in same format as the chr20.fa in the directory from which
    the script is launched

OUTPUT:

nuc_count output

Raw Counts
A:  17867246
C:  13916133
G:  14094472
T:  18066406
N:  499910

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

make_seq output:

>random_sequence
CCGCTCCTCCGTGGCGTGGGTGCTCGTGCCCAAACAGCTCTATGTGGAGATCTTAAGGCTCTCTCAGTGTGCTTCTCGTGTCCTCGGTACCGTCGAAAAC
ACGAGAACCAGGATAACGCTCCCCACAGAGGGATGGTGACTAAGCTAGACCGCAGGTTGTAGAAGCGTAGCGGCTGCCTACGACCAGTACGGCGTCGGC
ACCGTTGTAGGCATCACGACCCCCGACGTTCCCGATACTCTGAGCAGATCCGCCGAATTCCATCGGGAGAATCCATACTTCCTCCTCCCCAGAGACGGT
ACTAGATCCACTCCAACGGCTAGTCAACCGTTAATCATGCGGTGGCCCCTCCACAGCTATAGAAATCAGGAACTAATGTCATGCCGTTGCCTCCAAAAC
GAATGCACCTGGCGTCCAGTATAGATCAGCACACCGGCACGTGAGGGTTTTAGCGACCGTCTGCTTTGCAACACCCCCCCACGCTCGCCAAGGGCCCCG
CGGCCTGTCGGTTGTATGGCGATTAGTGCGTGTTACAAGGATAGCGCAGCGGTGGTCACGCATTAAGCACAGCTCGGAAATGCCATGTCTCATTTTTCA
CTGGCGGCTCGTGATGTGACCTCTCAAGTGAATAGAGAGTACCATAACCTCCCGGGGAAATTAACGTATCGGATCGGTCTTATAAAGCTTGCATCCTAT
CATGCTTTGGCTGCACACGGTTGGAACTAGCTATGGGCTGTGCGGATTGCTCCCCTAGATCTGGCTGTTCTCGAGTACTTTTAACGTCTTAGCTTGCGA
GACGTTCATGGCTGTCTCCTCAACTCGCTATACATTTTCCCCAAGTTCTTCGGGCAGATGAATTCTGTGGTGAGTTACGCTAGGCGTTCTGGAACACAC
TCGGCTCTTGACCCCATGCGTTCAATTCGATCCGGTGTACCCGTATCAGCGGCTCCTTAGACTTCTAAGACTAAGGTCGCTAACAACTGGTCTGACCGG
GGCTATGTT

 $: ./nuc_count.py random_sequence.fa

Raw Counts
A:  216
C:  287
G:  253
T:  244
N:  0
Sequence length: 1000

 The frequency of A is 0.216

 The frequency of T is 0.244

 The frequency of G is 0.253

 The frequency of C is 0.287


The Frequency of dinucleotides is:
The frequency of dinucleotide sequence AA in chromosome 20 is 0.04804804804804805
The frequency of dinucleotide sequence AC in chromosome 20 is 0.06106106106106106
The frequency of dinucleotide sequence AG in chromosome 20 is 0.056056056056056056
The frequency of dinucleotide sequence AT in chromosome 20 is 0.05105105105105105
The frequency of dinucleotide sequence CA in chromosome 20 is 0.06106106106106106
The frequency of dinucleotide sequence CC in chromosome 20 is 0.07707707707707707
The frequency of dinucleotide sequence CG in chromosome 20 is 0.07107107107107107
The frequency of dinucleotide sequence CT in chromosome 20 is 0.07807807807807808
The frequency of dinucleotide sequence GA in chromosome 20 is 0.05405405405405406
The frequency of dinucleotide sequence GC in chromosome 20 is 0.07107107107107107
The frequency of dinucleotide sequence GG in chromosome 20 is 0.06306306306306306
The frequency of dinucleotide sequence GT in chromosome 20 is 0.06506506506506507
The frequency of dinucleotide sequence TA in chromosome 20 is 0.05305305305305305
The frequency of dinucleotide sequence TC in chromosome 20 is 0.07707707707707707
The frequency of dinucleotide sequence TG in chromosome 20 is 0.06306306306306306
The frequency of dinucleotide sequence TT in chromosome 20 is 0.05005005005005005
