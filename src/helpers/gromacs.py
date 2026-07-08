
import numpy as np
import matplotlib.pyplot as plt

def parse_xvg(filepath):
    """
    Parse a GROMACS .xvg file, skipping comment (#) and metadata (@) lines.
    Returns a dict with 'time' (first column) and 'columns' (remaining columns),
    plus 'labels' extracted from @ axis and legend metadata.
    """
    data_lines = []
    xlabel, ylabel = "Time (ps)", "Value"
    legend_labels = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            elif line.startswith("@"):
                # Extract axis labels and legend entries
                if '"' in line:
                    label = line.split('"')[1]
                    if "xaxis" in line and "label" in line:
                        xlabel = label
                    elif "yaxis" in line and "label" in line:
                        ylabel = label
                    elif line.startswith("@ s") and "legend" in line:
                        legend_labels.append(label)
                continue
            else:
                vals = line.split()
                if vals:
                    data_lines.append([float(v) for v in vals])

    data = np.array(data_lines)
    return {
        "time": data[:, 0],
        "columns": data[:, 1:],
        "xlabel": xlabel,
        "ylabel": ylabel,
        "legend": legend_labels,
    }


def plot_xvg(filepath, title=None, figsize=(10, 4)):
    """
    Plot all data columns from an .xvg file with appropriate labels.
    """
    d = parse_xvg(filepath)
    fig, ax = plt.subplots(figsize=figsize)

    n_cols = d["columns"].shape[1]
    for i in range(n_cols):
        label = d["legend"][i] if i < len(d["legend"]) else f"Column {i+1}"
        ax.plot(d["time"], d["columns"][:, i], linewidth=0.8, label=label)

    ax.set_xlabel(d["xlabel"], fontsize=12)
    ax.set_ylabel(d["ylabel"], fontsize=12)
    if title:
        ax.set_title(title, fontsize=14)
    if n_cols > 1:
        ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
    return d


print("Helper functions loaded.")