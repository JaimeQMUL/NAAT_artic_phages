from Bio import Phylo
import matplotlib.pyplot as plt
from collections import defaultdict
import random

def PlotPhyloGroups(tree_file, mapping, output):
    for key in mapping.keys():

        tree = Phylo.read(tree_file, "newick")

        keep = set(mapping[key])  # faster lookup

        for term in list(tree.get_terminals()):  # 👈 important: list()
            if term.name not in keep:
                tree.prune(term)

        fig = plt.figure(figsize=(20, 80))
        axes = fig.add_subplot(1, 1, 1)

        Phylo.draw(
            tree,
            axes=axes,
            branch_labels=lambda c: None
        )

        name=f'{key}_tree.png'
        filename = f'{output}/{name}'

        fig.savefig(filename, dpi=300, bbox_inches="tight")
        plt.close()



def PlotWholePhylogeny(tree_file, output_file):

    tree = Phylo.read(tree_file, "newick")
    fig = plt.figure(figsize=(20, 80))
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(
        tree,
        axes=axes,
        branch_labels=lambda c: None
    )

    fig.savefig(output_file, dpi=300)
    plt.close()


def PlotPhyloTreeGrouped(tree_file, mapping, output):
    # ----------------------------
    # 1. Load tree
    # ----------------------------
    tree = Phylo.read(tree_file, "newick")

    # ----------------------------
    # 2. Invert mapping: acc → [groups]
    # ----------------------------
    acc_to_groups = defaultdict(list)

    for group, accs in mapping.items():
        for acc in accs:
            acc_to_groups[acc].append(group)

    # ----------------------------
    # 3. Define colours (0–255 for Biopython)
    # ----------------------------
    def to_255(color):
        return tuple(int(255 * c) for c in color)

    cmap = plt.cm.tab20.colors

    groups = list(mapping.keys())

    group_to_color = {
        group: to_255(cmap[i % len(cmap)])
        for i, group in enumerate(groups)
    }

    MULTI_COLOR = (128, 0, 128)  # purple
    DEFAULT_COLOR = (0, 0, 0)  # black

    # ----------------------------
    # 4. Recursive colouring
    # ----------------------------
    def assign_clade_color(clade):

        # Leaf node
        if clade.is_terminal():
            name = clade.name
            groups_here = acc_to_groups.get(name, [])

            if len(groups_here) == 1:
                group = groups_here[0]
                clade.color = group_to_color.get(group, DEFAULT_COLOR)
                return group

            elif len(groups_here) > 1:
                clade.color = MULTI_COLOR
                return "MULTI"

            else:
                clade.color = DEFAULT_COLOR
                return None

        # Internal node
        child_groups = set()

        for child in clade.clades:
            g = assign_clade_color(child)
            if g:
                child_groups.add(g)

        # If all children share same group → colour branch
        if len(child_groups) == 1:
            group = next(iter(child_groups))
            clade.color = group_to_color.get(group, DEFAULT_COLOR)
            return group

        # Mixed clade → no colour
        clade.color = DEFAULT_COLOR
        return None

    assign_clade_color(tree.root)

    # ----------------------------
    # 5. Plot
    # ----------------------------
    fig = plt.figure(figsize=(20, 80))
    axes = fig.add_subplot(1, 1, 1)

    Phylo.draw(
        tree,
        axes=axes,
        label_func=lambda x: None  # remove labels
    )

    # ----------------------------
    # 6. Save
    # ----------------------------
    plt.savefig(output, dpi=300, bbox_inches="tight")
    plt.close()
