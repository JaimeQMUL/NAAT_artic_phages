# Testing boltz using P04529
# Created YAMLs as specified in boltz documentation. Only contain a field for the sequence. for each of my gold standards

# Getting structural prediction of gold standards
# UvsX
boltz predict uvsx_input.yaml --out_dir uvsx_results --use_msa_server

# 7Z3M
boltz predict 7Z3M_input.yaml --out_dir 7Z3M_results --use_msa_server

# 9GBG
boltz predict 9GBG_input.yaml --out_dir 9GBG_results --use_msa_server


# To improve the prediction accuracy looked at trimming down the P04529 fasta based on pdb and swiss modelling'
cd uvsx_modelling/sequences
# pdb trim
awk 'BEGIN{seq=""}
     /^>/ {print; next}
     {seq=seq $0}
     END{
       print substr(seq,30,358-30+1)
     }' P04529.fasta > P04529_pdb.fasta

cat P04529_pdb.fasta | tail -n 1 >> ../inputs/pdb_input.yaml


# swiss trim
awk 'BEGIN{seq=""}
     /^>/ {print; next}
     {seq=seq $0}
     END{
       print substr(seq,31,342-31+1)
     }' P04529.fasta > P04529_swiss.fasta

cat P04529_swiss.fasta | tail -n 1 >> ../inputs/swiss_input.yaml

cd ..
boltz predict inputs/pdb_input.yaml --out_dir results --use_msa_server
boltz predict inputs/swiss_input.yaml --out_dir results --use_msa_server

