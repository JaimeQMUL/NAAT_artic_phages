# Using MSA to improve boltz confidence
# Created A3M MSA to better inform boltz about conserved motifs.

source /opt/anaconda3/etc/profile.d/conda.sh
conda activate cold
#
#boltz predict inputs/uvsx_input.yaml --out_dir results
#boltz predict inputs/7Z3M_input.yaml --out_dir results
#boltz predict inputs/9GBG_input.yaml --out_dir results

# Comparing structures
# Align predictions to crystal
echo "===== P04529 =====" > results/usalign_results.txt
usalign results/boltz_results_uvsx_input/predictions/uvsx_input/uvsx_input_model_0.cif ../../4-MolecularDynamics/structures/gold_standards/crystal/3IO5.pdb >> results/usalign_results.txt

echo "====== 7Z3M ======" >> results/usalign_results.txt
usalign results/boltz_results_7Z3M_input/predictions/7Z3M_input/7Z3M_input_model_0.cif ../../4-MolecularDynamics/structures/gold_standards/crystal/7Z3M.pdb >> results/usalign_results.txt

echo "====== 9GBG ======" >> results/usalign_results.txt
usalign results/boltz_results_9GBG_input/predictions/9GBG_input/9GBG_input_model_0.cif ../../4-MolecularDynamics/structures/gold_standards/crystal/9GBG.pdb >> results/usalign_results.txt



