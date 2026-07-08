# rename 3IO5.pdb to uvsx.pdb
mv ../../crystal/3IO5.pdb ../../crystal/uvsx.pdb

proteins=("uvsx" "7Z3M" "9GBG")

for p in "${proteins[@]}"
do
  outfile="${p}/structural_alignment.txt"
  ref="../../crystal/${p}.pdb"
  echo "========================== ${p} ==========================" > "$outfile"
  for cif in "${p}"/*.cif
  do
    usalign "$cif" "$ref" >> "$outfile"
  done
done

# rename it back
mv ../../crystal/uvsx.pdb ../../crystal/3IO5.pdb


# Converting best models to PDB
maxit -input uvsx/fold_uvsx_model_0.cif -output uvsx/uvsx.pdb -o 2
maxit -input 7Z3M/fold_7z3m_model_3.cif -output 7Z3M/7Z3M.pdb -o 2
maxit -input 9GBG/fold_9gbg_model_1.cif -output 9GBG/9GBG.pdb -o 2