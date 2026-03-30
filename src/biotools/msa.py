from pathlib import Path
import subprocess
import os

def RunMSA(protein_name):

    fasta_dir = Path(f"{protein_name}/data/references")  # change this
    combined_file = Path("combined.fasta")

    with combined_file.open("w") as out_f:
        for fasta_path in fasta_dir.glob("*.fasta"):
            with fasta_path.open() as f:
                out_f.write(f.read())
                if not f.read().endswith("\n"):
                    out_f.write("\n")



    input_fasta = "combined.fasta"
    aligned_fasta = f"{protein_name}/data/references/aligned.fasta"

    # "--auto" lets MAFFT pick the best algorithm for your dataset
    subprocess.run(
        ["/opt/anaconda3/envs/cold/bin/mafft", "--auto", input_fasta],
        stdout=open(aligned_fasta, "w"),
        env=os.environ  # passes your current PATH to the subprocess
    )

    #Delete combined fasta
    combined_file.unlink()