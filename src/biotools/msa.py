from pathlib import Path
import subprocess
import os

def RunMSA(save_dir, filename, files):

    fasta_dir = Path(save_dir)
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
    aligned_fasta = f"{save_dir}/data/alignments/{filename}.fasta"

    subprocess.run(
        ["/opt/anaconda3/envs/cold/bin/mafft", "--auto", input_fasta],
        stdout=open(aligned_fasta, "w"),
        env=os.environ
    )

    # Delete combined fasta
    combined_file.unlink()

