# Creating Curated Database
To do this is used the uniprot, ncbi and interpro databases. 
For uniprot and ncbi i used search queries for any uvsx-like/recombinase proteins.
For interpro i took the 2 interpro ids associated with uvsx and used an api call to extrct all proteins also containing these ids.
I saved all sequences to a fasta file depending on which database and search method they came from aswell as saving all the metadata I could get for each sequence.
I found that the accessions obtaining from the interpro search had complete overlap with the ncbi search query.
Because the ncbi api call returned more sequences I archived the interpro fasta as these sequences would all be found in the ncbi fasta.
The metadata however would be useful for me at somepoint

This returned roughly 25,000 sequences. However this was uncleaned data in which there may be many duplicate accessions and identical protein sequences from the same species.

## Cleaning data
To clean data I first sorted all accessions to only keep the unique ones. I had taxon ids for all my sequences and so could check what species each individual came from.
For all my unique accessions I grouped by sequence. All accessions with identical sequences would be grouped together. I then grouped by species using taxon id and disregarded any individual with the same sequence and from the same species as the first individual found of that species.
This meant i now had a clean database of non-redudant sequences. There wouldnt be any identical sequences found from the same species.
I made the decision of keeping identical sequences if they had came from different species as this may be useful for downstream analysis.
This left me with roughly 20,000 sequences.

## Filtering data
My database had been cleaned however there were still some very short sequences and very long sequences perhaps added due to poor annotation. 
To filter data i simply cut any sequences shorter than 341 amino acids and longer than 441 amino acids. 
UvsX t4 is a protein 391 amino acids long and so my reasoning was any protein significantly shorter or longer than this most likely arised due to some database annotation error.

# Gold standards
To guide our discovery process it was crucial to set a ground truth. This would obviously be majority guided by UvsX as it is the primary protein used in rpa
There were also some other orthologs that have been discovered which we could use. 7Z3M and 9GBG were 2 UvsX orthologs that had been discovered from a proprietary database and wet lab validated to confrim their thermal ranges
K35G/E198R was an engineered variant using directed evolution that displayed higher catalytic activity. This study also identified a few more engineered variants however using all of them in our analysis could cause bias as the sequences were entirely identical except for single amino acid changes. K35G/E198R was the best performing variant and so it was the one we chose to keep in the analysis
Crystal structures of all the gold standards are avaliable boosting exploratory power.
We had 4 gold standards UvsX t4, K35G/E198R, 7Z3M and 9GBG Which we added to our database

# Scoring System
To assess how likely each sequence in the curated database is to be a true UvsX ortholog for downstream thermal classification
and molecular dynamics analysis. We decided a scoring system would be best. This would be made up of hits from Hidden
Markov Models and presence of key domains using string matches.

## Hidden Markov Models
We created 5 Hidden Markov Models (hmms) to scan our database. These were created using hmmer3.
2 Hmms were built using the seed alignments of the domains annotated on pfam (PF00154 & PF21134).
Seed alignments were installed and `hmmbuild` was used to construct a profile.
The other 3 were custom hmms. 
One profile was seeded using an alignment of the 4 gold standards. 
The final 2 used the pfam hmms to search the gold standards, the regions that matched were extracted and aligned to seed the last two profiles that made up our 5 hmm profiles.
Interestingly PF21134 found no matches in the 2 novel discoveries (7Z3M and 9GBG) as this domain was in the c terminus, which is highly variable as it is involved in the polymerisation of recombinase proteins and so not conserved between species
This meant that these hmms would be much less useful and so it was removed from downstream analysis. 

This left only 3 hmms The two built from the PF00154 domain and the gold standard alignment hmm

## Motif checks
### ATP binding sites 
In UvsX there are 2 domains responsible for ATP binding The walker a and walker b motifs
#### Walker A
Found in T4 UvsX at positions 59-67. In the wider literature this is classified as a region that follows GxxxxGK[T/S].
Interestingly in UvsX the walker motif has a F in place of the signature G. GxxxFKS. The string we were searching for
was updated and we scored protein on whether this string was present Gxxxx[G/F]K[T/S]

#### Walker B
Found in T4 UvsX at position 138-143 it is characterized by a sequence of hhhh[D/E] with h being any hydrophobic residue
### DNA binding sites


# Phylogeny construction
To explore relationships between sequences in our database and their relation to our gold standards we wanted to
create phylogentic trees. 
To do this we first aligned the sequences in the filtered database. We ran a mafft alignment which took far too long and
so instead looked to use famsa. The famsa alignment was trimmed using trimal to remove low occupancy regions. 
Interestingly this trimmed alignment had only 6 residues for each sequence which corresponded to the walker b motif. 
This highlighted that the creation of the curated database was mainly influenced by the prescence of this motif. Highlighting 
a limitation in the computational annotation tools used by the uniprot and ncbi databases. 


# Predicting Optimal temperatures
We used seq2topt to predict the optimal temperatures. This is an ab-initio method that uses only the sequence to output
a thermal optima. When sandboxing it on our wet-lab validated gold standards we found it performed poorly,
classifying the deep sea ortholog found at 1000m below sea level at 0.4c, as having a much higher optimal temperature
than the ortholog discovered in a volcano. 


 