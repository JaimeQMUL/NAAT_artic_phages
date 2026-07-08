conda create -n colabfold -c conda-forge -c bioconda python=3.10 kalign2=2.04 hhsuite=3.3.0 mmseqs2 -y
conda init batch
conda activate colabfold
pip install colabfold[alphafold]

# making directories
mkdir sequences alignments
mkdir alignments/P04529 alignments/7Z3M alignments/9GBG

# copying gold standard sequences to this location
cp ../../uvsx/data/references/P04529.fasta sequences
cp ../../uvsx/data/references/7Z3M.fasta sequences
cp ../../uvsx/data/references/9GBG.fasta sequences

# Generating a3m
colabfold_batch sequences/P04529.fasta alignments/P04529/ --msa-only
colabfold_batch sequences/7Z3M.fasta alignments/7Z3M/ --msa-only
colabfold_batch sequences/9GBG.fasta alignments/9GBG --msa-only

# these alignments used the ref100 database, seemed to have poor coverage.
