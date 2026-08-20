# Generating a3m alignments for each of the gold standards

# setting up directories
mkdir inputs outputs

# setting up inputs
# copying protein fastas
cp ../../../1-DatabaseCuration/gold_standards/P04529.fasta inputs
cp ../../../1-DatabaseCuration/gold_standards/7Z3M.fasta inputs
cp ../../../1-DatabaseCuration/gold_standards/9GBG.fasta inputs

# copying top hits database
cp ../../../3-InsilicoScreening/top_hits/sequences/top_hits.fasta inputs

cd inputs

# Create database of top hitsd
mmseqs createdb top_hits.fasta top_hits_DB

# create query dbs
mmseqs createdb P04529.fasta P04529DB
mmseqs createdb 7Z3M.fasta 7Z3MDB
mmseqs createdb 9GBG.fasta 9GBGDB

# Search database
mmseqs search P04529DB top_hits_DB P04529_result tmp -s 7.5 -e 1e-3
mmseqs search 7Z3MDB top_hits_DB 7Z3M_result tmp -s 7.5 -e 1e-3
mmseqs search 9GBGDB top_hits_DB 9GBG_result tmp -s 7.5 -e 1e-3

# Convert search to MSA
mmseqs result2msa P04529DB top_hits_DB P04529_result ../outputs/P04529.a3m
mmseqs result2msa 7Z3MDB top_hits_DB 7Z3M_result ../outputs/7Z3M.a3m
mmseqs result2msa 9GBGDB top_hits_DB 9GBG_result ../outputs/9GBG.a3m

cd ../outputs
#clean a3m of null characters '\000'
tr -d '\000' < P04529.a3m > tmp.a3m && mv tmp.a3m P04529.a3m
tr -d '\000' < 7Z3M.a3m > tmp.a3m && mv tmp.a3m 7Z3M.a3m
tr -d '\000' < 9GBG.a3m > tmp.a3m && mv tmp.a3m 9GBG.a3m





# potentially make sequence coverage plots similar to that seen in the colabfold alignments.