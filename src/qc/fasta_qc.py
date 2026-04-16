# Quality control function for fasta files

def CountAlignmentGaps(file):
    gap_counts = []
    current_seq = ""

    with open(file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_seq:
                    gap_counts.append(current_seq.count('-'))
                    current_seq = ""
            else:
                current_seq += line.strip()

        # last sequence
        if current_seq:
            gap_counts.append(current_seq.count('-'))

    mean_gap = sum(gap_counts) / len(gap_counts)
    return mean_gap, gap_counts



import numpy as np

def AlignmentOccupancy(alignment_file):
    """
    Compute column-wise and mean occupancy from a multiple sequence alignment FASTA file.

    Parameters
    ----------
    alignment_file : str
        Path to aligned FASTA file

    Returns
    -------
    mean_occupancy : float
        Average occupancy across all columns
    occupancy_per_position : list of float
        Occupancy value for each alignment column
    """

    # Read sequences from FASTA
    sequences = []
    current_seq = ""

    with open(alignment_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_seq:
                    sequences.append(current_seq)
                    current_seq = ""
            else:
                current_seq += line

        if current_seq:
            sequences.append(current_seq)

    # Convert to numpy array (columns = alignment positions)
    arr = np.array([list(seq) for seq in sequences])

    # Compute occupancy per column
    occupancy_per_position = np.mean(arr != '-', axis=0)

    # Summary statistic
    mean_occupancy = np.mean(occupancy_per_position)

    return mean_occupancy, occupancy_per_position.tolist()