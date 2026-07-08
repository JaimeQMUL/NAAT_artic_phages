# Script to generate jobs for molecular dynamics
# To be used on apocrita

import subprocess
from pathlib import Path

proteins = ["3IO5", "9GBG"]

conditions = {
    "psychrophilic": 288,
    "mesophilic": 310,
    "thermophilic": 338
}

PROJECT_ROOT = str(Path("/data/home/bt24715/dissertation/gromacs/gold_standards/crystal").resolve())

# Looping over proteins to run MD
for protein in proteins:
    for condition in conditions.keys():
        workdir = Path(PROJECT_ROOT) / protein / condition
        subprocess.run([
            "sbatch",
            "job.sh",
            protein,
            condition,
            str(workdir),
            PROJECT_ROOT
        ])