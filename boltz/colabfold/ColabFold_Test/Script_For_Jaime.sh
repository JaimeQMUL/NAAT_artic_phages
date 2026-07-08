# Set up Environment

conda create -n colabfold -c conda-forge -c bioconda python=3.10 kalign2=2.04 hhsuite=3.3.0 mmseqs2 -y
conda init batch
conda activate colabfold
pip install colabfold[alphafold]

# Run colabfold_batch with --msa-only option to then generate the a3m file for downstream use. 
colabfold_batch ./test.fa output/ --msa-only

# within the output directory, you will find the a3m file for each sequence in the test.fa file. 
# You can use these a3m files for the boltz-2 analysis