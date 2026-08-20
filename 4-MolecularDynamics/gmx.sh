#!/bin/bash
#SBATCH --job-name=gromacs_md
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --time=240:00:00
#SBATCH --output=md_%j.out
#SBATCH --error=md_%j.err


set -e

# =========================================================
# INPUT ARGUMENTS
# =========================================================
PROTEIN=$1
CONDITION=$2
WORKDIR=$3
PROJECT_ROOT=$4

echo "======================================"
echo "Protein:   $PROTEIN"
echo "Condition: $CONDITION"
echo "Workdir:   $WORKDIR"
echo "Root:      $PROJECT_ROOT"
echo "======================================"

# =========================================================
# LOAD MODULES (IMPORTANT: INSIDE JOB)
# =========================================================
module load gromacs

# =========================================================
# MOVE INTO WORKDIR
# =========================================================
mkdir -p "$WORKDIR"
cd "$WORKDIR"


# =========================================================
# INPUT PDB
# =========================================================
PDB="${PROJECT_ROOT}/${PROTEIN}/${PROTEIN}.pdb"
if [ ! -f "$PDB" ]; then
    echo "ERROR: PDB not found at $PDB"
    exit 1
fi

# =========================================================
# CLEAN PDB (extra cleaning)
# =========================================================
grep -vE "HETATM|CONECT" "$PDB" > protein.pdb


grep "LYS" protein.pdb | head

# =========================================================
# TOPOLOGY GENERATION
# =========================================================
gmx pdb2gmx \
-f protein.pdb \
-o processed.gro \
-p topol.top \
-water tip3p \
-ff charmm27 \
-ignh

# =========================================================
# BOX + SOLVATION
# =========================================================
gmx editconf \
    -f processed.gro \
    -o boxed.gro \
    -c -d 1.0 -bt dodecahedron

gmx solvate \
    -cp boxed.gro \
    -cs spc216.gro \
    -o solvated.gro \
    -p topol.top

# =========================================================
# IONS
# =========================================================
cat > ions.mdp << EOF
integrator = steep
emtol = 1000
emstep = 0.01
nsteps = 500
EOF

gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr

printf "SOL\n" | gmx genion \
  -s ions.tpr \
  -o solv_ions.gro \
  -p topol.top \
  -pname K \
  -nname CL \
  -neutral \
  -conc 0.10

# =========================================================
# ENERGY MINIMISATION
# =========================================================
gmx grompp \
    -f "${PROJECT_ROOT}/${PROTEIN}/${CONDITION}/inputs/emin-charmm.mdp" \
    -c solv_ions.gro \
    -p topol.top \
    -o em.tpr \
    -maxwarn 1

gmx mdrun -deffnm em -ntomp 8

printf "Potential\n0\n" | gmx energy -f em.edr -o potential.xvg


# =========================================================
# NVT
# =========================================================
gmx grompp \
    -f "${PROJECT_ROOT}/${PROTEIN}/${CONDITION}/inputs/nvt-charmm.mdp" \
    -c em.gro \
    -r em.gro \
    -p topol.top \
    -o nvt.tpr

gmx mdrun -deffnm nvt -ntomp 8

echo "Temperature" | gmx energy -f nvt.edr -o temperature.xvg -b 20


# =========================================================
# NPT
# =========================================================
gmx grompp \
    -f "${PROJECT_ROOT}/${PROTEIN}/${CONDITION}/inputs/npt-charmm.mdp" \
    -c nvt.gro \
    -r nvt.gro \
    -t nvt.cpt \
    -p topol.top \
    -o npt.tpr

gmx mdrun -deffnm npt -ntomp 8

echo "Pressure" | gmx energy -f npt.edr -o pressure.xvg
echo "Density" | gmx energy -f npt.edr -o density.xvg


# =========================================================
# PRODUCTION MD
# =========================================================
gmx grompp \
    -f "${PROJECT_ROOT}/${PROTEIN}/${CONDITION}/inputs/md-charmm.mdp" \
    -c npt.gro \
    -t npt.cpt \
    -p topol.top \
    -o md.tpr

gmx mdrun -deffnm md -ntomp 8

# =========================================================
# ANALYSIS
# =========================================================
printf "1" | gmx trjconv -s md.tpr -f md.xtc -o md_whole.xtc -pbc whole
printf "1" | gmx trjconv -s md.tpr -f md_whole.xtc -o md_nojump.xtc -pbc nojump
printf "1\n1" | gmx trjconv -s md.tpr -f md_nojump.xtc -o md_center.xtc -center -pbc mol
printf "1\n" | gmx mindist -s md.tpr -f md_center.xtc -pi -od mindist.xvg

printf "4\n1\n" | gmx rms -s em.tpr -f md_center.xtc -o rmsd.xvg -tu ns
echo "1" | gmx gyrate -s md.tpr -f md_center.xtc -o gyrate.xvg
printf "1\n" | gmx rmsf -s md.tpr -f md_center.xtc -o rmsf.xvg -res


# =========================================================
# SAVE RESULTS
# =========================================================
RESULTS_DIR="${PROJECT_ROOT}/${PROTEIN}/results/${CONDITION}/results"

mkdir -p "$RESULTS_DIR"

cp *.xtc *.xvg *.gro *.edr *.log *.tpr \
   "$RESULTS_DIR/" 2>/dev/null || true

echo "======================================"
echo "DONE: ${PROTEIN} ${CONDITION}"
echo "======================================"