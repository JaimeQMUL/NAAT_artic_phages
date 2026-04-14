Seq2topt doesnt perform very well high mse and lack of psychrophile representtions makes it poor at estimating opt for my seqs.

Could look to create my own model
Found a website with 30,000 protein sequences with temps that we could train a model on.

Esm database has biomes sequences were collected in, could use this to help with model or in the literature search. 
https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2022_05/README.txt
downside is the license for the mgnify dataset


Tome
This repo has options for predicting optimal temperature
https://github.com/EngqvistLab/Tome


Enhancing ML prediction of topt
https://pmc.ncbi.nlm.nih.gov/articles/PMC11173260/#sec3-ijms-25-06252
This study removes conserved amino acids as they believe them to have little effect on thermostability
Which makes sense as regardless of thermostability conserved residues will remain, thermostability must dance around cnservation of keys functions
https://github.com/cyinyin/Tpho/tree/main


ProTherm
https://web.iitm.ac.in/bioinfo2/prothermdb/Organism.html
Has optimal temps of 12042 proteins can use this for predictions
