# Notes on alignments
MAFFT as extremely slow for the filtered database
Database was filtered only to keep seqs below 450 and above 350. 

We decided to use FAMSA instead as it was much quicker as demostrated in the benchmarking.py.
When famsa alignment as trimmed all that was left for each seq was 6 amino acids. These corresponded with each sequences walker B motif.
This meant our current hidden markov model was being majority influenced by the presence of the walker b motif in all of these sequences.
This is only one of 4 important domains and highlights a serious limitation in our HMMs.

>YDI94551.1
VVVFDS
>YDI88709.1
IIFIDS
>YDD38469.1
IVFIDS
>YCR04191.1
IVFIDS
>YCQ78089.1
IIFIDS
>YCQ76618.1
IIFIDS
>YCP71262.1
LIIWDS
>sp|P04529.2|UVSX_BPT4
VVFIDS


Notes
make an alignments of only the sequences with hits for the scoring methods we will use going forward, dont include the redundant pfam domains.

