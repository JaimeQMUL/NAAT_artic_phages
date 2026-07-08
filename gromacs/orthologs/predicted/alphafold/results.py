# imports
from src.helpers.gromacs import *
import matplotlib.pyplot as plt
import numpy as np
from src.biotools.fasta_tools import *

# Mapping of accessions
mapping={'psychrophiles': ['YCB25778', 'UFK27161', 'YCQ78089'],
 'thermophiles': ['YP_874025', 'YP_004782219', 'AHZ95250'],
 'mesophiles': ['XCZ64124', 'WAX22921', 'UUV43823'],
 'hottest': ['YP_009791320', 'YP_004322308', 'XYR95885'],
 'coldest': ['XQO94789', 'YP_009151072', 'YP_007010273'],
 'random': ['YP_009213881', 'UVK61175', 'CAD5240622']}

accs=[item for sublist in mapping.values() for item in sublist]

for acc in accs:

    ########################################################################################################################
    # PLOT RMSD
    ########################################################################################################################


    files = {
        "Psychrophilic": f"{acc}/psychrophilic/results/rmsd.xvg",
        "Mesophilic": f"{acc}/mesophilic/results/rmsd.xvg",
        "Thermophilic": f"{acc}/thermophilic/results/rmsd.xvg",
    }

    colors = {
        "Psychrophilic": "blue",
        "Mesophilic": "orange",
        "Thermophilic": "red",
    }

    fig, ax = plt.subplots(figsize=(10, 4))

    for label, file in files.items():
        d = parse_xvg(file)

        x = d["time"]
        rmsd = d["columns"][:, 0]

        ax.plot(x, rmsd, linewidth=0.8, label=label, color=colors[label])

        print(f"\n{label}")
        print(f"Final RMSD:         {rmsd[-1]:.3f} nm ({rmsd[-1]*10:.1f} Å)")
        print(f"Mean RMSD:          {np.mean(rmsd):.3f} nm ({np.mean(rmsd)*10:.1f} Å)")
        print(f"Mean (last 25%):    {np.mean(rmsd[len(rmsd)*3//4:]):.3f} nm "
              f"({np.mean(rmsd[len(rmsd)*3//4:])*10:.1f} Å)")

    ax.set_title("Backbone RMSD vs Crystal Structure")
    ax.set_xlabel(d["xlabel"])
    ax.set_ylabel("RMSD (nm)")

    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'plots/{acc}_RMSD.png', dpi=300, bbox_inches='tight')



    ########################################################################################################################
    # PLOT RG
    ########################################################################################################################


    files = {
        "Psychrophilic": f"{acc}/psychrophilic/results/gyrate.xvg",
        "Mesophilic": f"{acc}/mesophilic/results/gyrate.xvg",
        "Thermophilic": f"{acc}/thermophilic/results/gyrate.xvg",
    }

    colors = {
        "Psychrophilic": "blue",
        "Mesophilic": "orange",
        "Thermophilic": "red",
    }


    fig, ax = plt.subplots(figsize=(10, 4))

    for label, file in files.items():
        d = parse_xvg(file)

        x = d["time"]
        rg = d["columns"][:, 0]

        ax.plot(x, rg, linewidth=0.8, label=label, color=colors[label])

        print(f"\n{label}")
        print(f"Mean Rg:            {np.mean(rg):.3f} nm")
        print(f"Std deviation:      {np.std(rg):.3f} nm")
        print(f"Mean Rg (last 25%): {np.mean(rg[len(rg)*3//4:]):.3f} nm")

    ax.set_title("Radius of Gyration")
    ax.set_xlabel(d["xlabel"])
    ax.set_ylabel("Rg (nm)")

    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'plots/{acc}_Rg.png', dpi=300, bbox_inches='tight')



    ########################################################################################################################
    # PLOT RMSF
    ########################################################################################################################

    files = {
        "Psychrophilic": f"{acc}/psychrophilic/results/rmsf.xvg",
        "Mesophilic": f"{acc}/mesophilic/results/rmsf.xvg",
        "Thermophilic": f"{acc}/thermophilic/results/rmsf.xvg",
    }

    colors = {
        "Psychrophilic": "blue",
        "Mesophilic": "orange",
        "Thermophilic": "red",
    }

    # Get motifs
    seq=ExtractSequence(acc, '../../../../uvsx/data/curated_database/top_hits.fasta')
    walker_a=r"[GA]....[GF]K[TS]"
    matches_a = [(m.start(),m.end(), m.group()) for m in re.finditer(walker_a, seq)]
    if matches_a:
        a_start=matches_a[0][0]
        a_end=matches_a[0][1]

    walker_b=r"[AVILMFWY]{4}DE?"
    matches_b = [(m.start(),m.end(), m.group()) for m in re.finditer(walker_b, seq)]
    if matches_b:
        b_start=matches_b[0][0]
        b_end=matches_b[0][1]


    fig, ax = plt.subplots(figsize=(10, 4))

    for label, file in files.items():
        d = parse_xvg(file)

        x = d["time"]
        rmsf = d["columns"][:, 0]

        ax.plot(x, rmsf, linewidth=0.8, label=label, color=colors[label])

        print(f"\n{label}")
        print(f"Mean Rmsf:            {np.mean(rmsf):.3f} nm")
        print(f"Std deviation:      {np.std(rmsf):.3f} nm")
        print(f"Mean Rmsf (last 25%): {np.mean(rmsf[len(rmsf)*3//4:]):.3f} nm")

    ax.set_title("Rmsf")
    ax.set_xlabel(d["xlabel"])
    ax.set_ylabel("Rmsf (nm)")


    # showing walker a and b
    ax.axvline(a_start, color="k", linestyle="--")
    ax.axvline(a_end, color="k", linestyle="--")

    ax.axvline(b_start, color="k", linestyle="--")
    ax.axvline(b_end, color="k", linestyle="--")

    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'plots/{acc}_RMSF.png', dpi=300, bbox_inches='tight')
