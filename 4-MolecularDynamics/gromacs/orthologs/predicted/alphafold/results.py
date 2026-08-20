# imports
from src.helpers.gromacs import *
import matplotlib.pyplot as plt
import numpy as np
from src.biotools.fasta_tools import *
import pandas as pd
import re
import os


# Mapping of accessions
mapping = {
    'psychrophiles': ['YCB25778', 'UFK27161', 'YCQ78089'],
    'thermophiles': ['YP_874025', 'YP_004782219', 'AHZ95250'],
    'mesophiles': ['XCZ64124', 'WAX22921', 'UUV43823'],
    'hottest': ['YP_009791320', 'YP_004322308', 'XYR95885'],
    'coldest': ['XQO94789', 'YP_009151072', 'YP_007010273'],
    'random': ['YP_009213881', 'UVK61175', 'CAD5240622']
}


# Flatten accession list
accs = [item for sublist in mapping.values() for item in sublist]


# Map each accession back to its original condition
accession_condition = {
    acc: condition
    for condition, accessions in mapping.items()
    for acc in accessions
}


# Create directory for plots if it doesn't already exist
os.makedirs("plots", exist_ok=True)


# Store all summary statistics here
summary = []


for acc in accs:

    # Dictionaries to store means for this protein
    rmsd_means = {}
    rg_means = {}
    rmsf_means = {}


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

        ax.plot(
            x,
            rmsd,
            linewidth=0.8,
            label=label,
            color=colors[label]
        )

        # Calculate statistics
        mean_rmsd = np.mean(rmsd)
        mean_rmsd_last25 = np.mean(rmsd[len(rmsd) * 3 // 4:])

        # Store mean RMSD
        rmsd_means[label] = mean_rmsd

        print(f"\n{acc} - {label}")
        print(f"Final RMSD:         {rmsd[-1]:.3f} nm ({rmsd[-1] * 10:.1f} Å)")
        print(f"Mean RMSD:          {mean_rmsd:.3f} nm ({mean_rmsd * 10:.1f} Å)")
        print(
            f"Mean (last 25%):    {mean_rmsd_last25:.3f} nm "
            f"({mean_rmsd_last25 * 10:.1f} Å)"
        )

    ax.set_title("Backbone RMSD vs Crystal Structure")
    ax.set_xlabel(d["xlabel"])
    ax.set_ylabel("RMSD (nm)")

    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(
        f"plots/{acc}_RMSD.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)


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

        ax.plot(
            x,
            rg,
            linewidth=0.8,
            label=label,
            color=colors[label]
        )

        # Calculate statistics
        mean_rg = np.mean(rg)
        std_rg = np.std(rg)
        mean_rg_last25 = np.mean(rg[len(rg) * 3 // 4:])

        # Store mean Rg
        rg_means[label] = mean_rg

        print(f"\n{acc} - {label}")
        print(f"Mean Rg:            {mean_rg:.3f} nm")
        print(f"Std deviation:      {std_rg:.3f} nm")
        print(f"Mean Rg (last 25%): {mean_rg_last25:.3f} nm")

    ax.set_title("Radius of Gyration")
    ax.set_xlabel(d["xlabel"])
    ax.set_ylabel("Rg (nm)")

    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(
        f"plots/{acc}_Rg.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)


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


    # Get sequence
    seq = ExtractSequence(
        acc,
        '../../../../../3-InsilicoScreening/top_hits/sequences/top_hits.fasta'
    )


    # Walker A motif
    walker_a = r"[GA]....[GF]K[TS]"
    matches_a = [
        (m.start(), m.end(), m.group())
        for m in re.finditer(walker_a, seq)
    ]

    if matches_a:
        a_start = matches_a[0][0]
        a_end = matches_a[0][1]
    else:
        a_start = None
        a_end = None


    # Walker B motif
    walker_b = r"[AVILMFWY]{4}DE?"
    matches_b = [
        (m.start(), m.end(), m.group())
        for m in re.finditer(walker_b, seq)
    ]

    if matches_b:
        b_start = matches_b[0][0]
        b_end = matches_b[0][1]
    else:
        b_start = None
        b_end = None


    fig, ax = plt.subplots(figsize=(10, 4))

    for label, file in files.items():

        d = parse_xvg(file)

        x = d["time"]
        rmsf = d["columns"][:, 0]

        ax.plot(
            x,
            rmsf,
            linewidth=0.8,
            label=label,
            color=colors[label]
        )

        # Calculate statistics
        mean_rmsf = np.mean(rmsf)
        std_rmsf = np.std(rmsf)
        mean_rmsf_last25 = np.mean(rmsf[len(rmsf) * 3 // 4:])

        # Store mean RMSF
        rmsf_means[label] = mean_rmsf

        print(f"\n{acc} - {label}")
        print(f"Mean RMSF:            {mean_rmsf:.3f} nm")
        print(f"Std deviation:        {std_rmsf:.3f} nm")
        print(f"Mean RMSF (last 25%): {mean_rmsf_last25:.3f} nm")


    ax.set_title("RMSF")
    ax.set_xlabel(d["xlabel"])
    ax.set_ylabel("RMSF (nm)")


    # Showing Walker A and B
    if a_start is not None:
        ax.axvline(
            a_start,
            color="k",
            linestyle="--"
        )

        ax.axvline(
            a_end,
            color="k",
            linestyle="--"
        )


    if b_start is not None:
        ax.axvline(
            b_start,
            color="k",
            linestyle="--"
        )

        ax.axvline(
            b_end,
            color="k",
            linestyle="--"
        )


    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(
        f"plots/{acc}_RMSF.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)


    ########################################################################################################################
    # STORE SUMMARY FOR THIS PROTEIN
    ########################################################################################################################

    summary.append({
        "accession": acc,
        "TMP/Seq2Topt_classification": accession_condition[acc],

        "mean_RMSD_psychrophilic": rmsd_means["Psychrophilic"],
        "mean_RMSD_mesophilic": rmsd_means["Mesophilic"],
        "mean_RMSD_thermophilic": rmsd_means["Thermophilic"],

        "mean_Rg_psychrophilic": rg_means["Psychrophilic"],
        "mean_Rg_mesophilic": rg_means["Mesophilic"],
        "mean_Rg_thermophilic": rg_means["Thermophilic"],

        "mean_RMSF_psychrophilic": rmsf_means["Psychrophilic"],
        "mean_RMSF_mesophilic": rmsf_means["Mesophilic"],
        "mean_RMSF_thermophilic": rmsf_means["Thermophilic"],
    })


############################################################################################################################
# SAVE SUMMARY
############################################################################################################################

summary_df = pd.DataFrame(summary)


# Save as CSV
summary_df.to_csv(
    "protein_summary.csv",
    index=False
)


print("\n" + "=" * 100)
print("SUMMARY")
print("=" * 100)

print(summary_df)

print("\nSummary saved to:")
print("  protein_summary.csv")