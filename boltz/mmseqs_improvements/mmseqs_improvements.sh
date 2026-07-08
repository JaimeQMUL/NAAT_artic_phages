# Running structural prediciton using the mmseqs alignments
# set up directories
mkdir inputs results
# creating yamls
cd inputs
touch P04529_mmseqs.yaml 7Z3M_mmseqs.yaml 9GBG_mmseqs.yaml
cd ..



for id in P04529 7Z3M 9GBG
# clean a3m headers
do
  msa="../../mmseqs/outputs/${id}.a3m"
  fasta="../../uvsx/data/references/${id}.fasta"
  seq=$(awk '!/^>/ {gsub(/\r/,""); printf "%s", $0}' "$fasta")
  {
    echo "sequences:"
    echo "  - protein:"
    echo "      id: A"
    echo "      sequence: $seq"
    echo "      msa: $msa"
  } > inputs/${id}_mmseqs.yaml
done



# Run boltz
boltz predict inputs/P04529_mmseqs.yaml --out_dir results
boltz predict inputs/7Z3M_mmseqs.yaml --out_dir results
boltz predict inputs/9GBG_mmseqs.yaml --out_dir results

# Align predictions to crystal
echo "===== P04529 =====" > results/usalign_results.txt
usalign results/boltz_results_P04529_mmseqs/predictions/P04529_mmseqs/P04529_mmseqs_model_0.cif ../../uvsx/data/structures/gold_standards/crystal/3IO5.pdb >> results/usalign_results.txt

echo "====== 7Z3M ======" >> results/usalign_results.txt
usalign results/boltz_results_7Z3M_mmseqs/predictions/7Z3M_mmseqs/7Z3M_mmseqs_model_0.cif ../../uvsx/data/structures/gold_standards/crystal/7Z3M.pdb >> results/usalign_results.txt

echo "====== 9GBG ======" >> results/usalign_results.txt
usalign results/boltz_results_9GBG_mmseqs/predictions/9GBG_mmseqs/9GBG_mmseqs_model_0.cif ../../uvsx/data/structures/gold_standards/crystal/9GBG.pdb >> results/usalign_results.txt


# unfortunately results see no improvement in terms of alignment accuracy between prediction and crystal