
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict

def PlotSequenceLengths(fasta_file):
    seq_lengths = []


    with open(fasta_file) as f:
        seq = ''

        for line in f:
            line = line.strip()

            if line.startswith('>'):
                if seq:
                    seq_lengths.append(len(seq))
                seq = ''
            else:
                seq += line

        if seq:
            seq_lengths.append(len(seq))
    bins=len(seq_lengths)//20

    fig, ax = plt.subplots()
    ax.hist(seq_lengths, bins)
    ax.set_xlabel("Sequence length")
    ax.set_ylabel("Frequency")
    ax.set_title("Distribution of sequence lengths")

    return fig, ax

def PlotRankCounts(counts, rank, top_n=10):
    df = pd.DataFrame(list(counts.items()), columns=['Virus', 'Count'])

    # Sort
    df = df.sort_values(by='Count', ascending=False)

    # Top N
    major = df.head(top_n)
    minor = df.iloc[top_n:]

    # Combine others
    other_count = minor['Count'].sum()
    if other_count > 0:
        major = pd.concat(
            [major, pd.DataFrame([['Other', other_count]], columns=['Virus', 'Count'])]
        )

    # Plot
    fig, ax = plt.subplots()
    major.plot(kind='bar', x='Virus', y='Count', ax=ax)

    ax.set_xlabel(rank)
    ax.set_ylabel('Count')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

    plt.tight_layout()

    return fig, ax




def PlotRanksSummary(lineages):
    collapse_map = {
        'subphylum': 'phylum',
        'suborder': 'order',
        'subfamily': 'family',
        'subgenus': 'genus'
    }

    rank_values = defaultdict(set)

    for entry in lineages:
        lineage = entry.get("lineage")
        if not lineage:
            continue

        for rank, value in lineage.items():
            new_rank = collapse_map.get(rank, rank)
            rank_values[new_rank].add(value)

    collapsed_counts = {rank: len(values) for rank, values in rank_values.items()}

    df = pd.DataFrame(list(collapsed_counts.items()), columns=['Rank', 'Count'])

    rank_order = [
        'acellular root', 'realm', 'kingdom', 'phylum', 'class', 'order',
        'family', 'genus', 'species'
    ]

    df['Rank'] = pd.Categorical(df['Rank'], categories=rank_order, ordered=True)
    df = df.sort_values('Rank')
    df = df[df['Rank'].notna()]

    sorted_counts = df['Count'].sort_values(ascending=False).values
    scale_value = sorted_counts[2] if len(sorted_counts) >= 3 else sorted_counts[-1]
    cap_value = scale_value * 2

    df_plot = df.copy()
    real_values = {}

    for i, row in df.iterrows():
        if row['Count'] > cap_value:
            real_values[row['Rank']] = row['Count']
            df_plot.loc[i, 'Count'] = cap_value

    # Create fig + ax
    fig, ax = plt.subplots()

    df_plot.plot(kind='bar', x='Rank', y='Count', legend=False, ax=ax)

    for i, row in df_plot.iterrows():
        value = row['Count']
        rank = row['Rank']

        # If capped, show real value (you already stored it)
        if rank in real_values:
            label = real_values[rank]
            y = cap_value
        else:
            label = value
            y = value

        ax.text(i, y, f"{label}",
                ha='center', va='bottom', fontsize=9)

    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_ylabel("Unique classifications")
    ax.set_xlabel('Ranks')

    plt.tight_layout()

    return fig, ax
