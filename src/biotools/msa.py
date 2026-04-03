from pathlib import Path
import subprocess
import os

def RunMSA(protein_name, files):

    fasta_dir = Path(f"{protein_name}/data/references")
    combined_file = Path("combined.fasta")

    with combined_file.open("w") as out_f:
        for file_name in files:
            fasta_path = fasta_dir / file_name  # use only specified files

            with fasta_path.open() as f:
                content = f.read()
                out_f.write(content)

                if not content.endswith("\n"):
                    out_f.write("\n")

    input_fasta = "combined.fasta"
    aligned_fasta = f"{protein_name}/data/references/aligned.fasta"

    subprocess.run(
        ["/opt/anaconda3/envs/cold/bin/mafft", "--auto", input_fasta],
        stdout=open(aligned_fasta, "w"),
        env=os.environ
    )

    # Delete combined fasta
    combined_file.unlink()