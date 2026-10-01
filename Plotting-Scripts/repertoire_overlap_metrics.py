#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.metrics.pairwise import cosine_similarity
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage

# --------------------
# Build clone matrix
# --------------------
def build_clone_matrix(
    df,
    normalize=False,
    min_cluster_size=1,
    binary=False
):

    required = {"ClusterID", "tissue", "nt_seq_count"}

    if not required.issubset(df.columns):
        raise ValueError(
            f"Missing required columns: {required}"
        )

    df = df.dropna(
        subset=["ClusterID", "tissue", "nt_seq_count"]
    )

    df["nt_seq_count"] = pd.to_numeric(
        df["nt_seq_count"],
        errors="coerce"
    )

    df = df.dropna(subset=["nt_seq_count"])

    # --------------------
    # Global cluster filter
    # --------------------
    cluster_sizes = (
        df.groupby("ClusterID")["nt_seq_count"]
        .sum()
    )

    keep_clusters = cluster_sizes[
        cluster_sizes >= min_cluster_size
    ].index

    df = df[
        df["ClusterID"].isin(keep_clusters)
    ]

    # --------------------
    # Tissue x clone matrix
    # --------------------
    clone_matrix = (
        df.groupby(["tissue", "ClusterID"])["nt_seq_count"]
        .sum()
        .unstack(fill_value=0)
    )

    # --------------------
    # Binary mode
    # --------------------
    if binary:
        clone_matrix = (
            clone_matrix > 0
        ).astype(int)

    # --------------------
    # Normalize tissues
    # --------------------
    if normalize:

        clone_matrix = clone_matrix.div(
            clone_matrix.sum(axis=1),
            axis=0
        ).fillna(0)

    return clone_matrix

# --------------------
# Cosine similarity
# --------------------
def compute_cosine_similarity_matrix(clone_matrix):

    sim = cosine_similarity(clone_matrix)

    sim_df = pd.DataFrame(
        sim,
        index=clone_matrix.index,
        columns=clone_matrix.index
    )

    return sim_df

# --------------------
# Bray-Curtis
# --------------------
def compute_bray_curtis_matrix(clone_matrix):

    dist = pdist(
        clone_matrix.values,
        metric="braycurtis"
    )

    dist_matrix = squareform(dist)

    bc_df = pd.DataFrame(
        dist_matrix,
        index=clone_matrix.index,
        columns=clone_matrix.index
    )

    return bc_df

# --------------------
# Morisita-Horn
# --------------------
def morisita_horn(x, y):

    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)

    Nx = x.sum()
    Ny = y.sum()

    if Nx == 0 or Ny == 0:
        return np.nan

    Dx = np.sum(x**2) / (Nx**2)
    Dy = np.sum(y**2) / (Ny**2)

    denominator = (Dx + Dy) * Nx * Ny

    if denominator == 0:
        return np.nan

    numerator = 2 * np.sum(x * y)

    return numerator / denominator

def compute_morisita_horn_matrix(clone_matrix):

    tissues = clone_matrix.index

    mh_df = pd.DataFrame(
        np.zeros((len(tissues), len(tissues))),
        index=tissues,
        columns=tissues
    )

    for i in tissues:
        for j in tissues:

            mh_df.loc[i, j] = morisita_horn(
                clone_matrix.loc[i],
                clone_matrix.loc[j]
            )

    return mh_df


# --------------------
# Shared lineage counts
# --------------------
def compute_shared_lineage_matrix(clone_matrix):
    """
    Count the number of lineages shared between every pair of tissues.

    Returns a tissue x tissue matrix where:
        diagonal = total lineages in tissue
        off-diagonal = shared lineages
    """

    presence = (clone_matrix > 0).astype(int)

    shared = pd.DataFrame(
        np.zeros((len(presence.index), len(presence.index)), dtype=int),
        index=presence.index,
        columns=presence.index
    )

    for tissue1 in presence.index:
        for tissue2 in presence.index:

            shared.loc[tissue1, tissue2] = (
                (presence.loc[tissue1] &
                 presence.loc[tissue2]).sum()
            )

    return shared

# --------------------
# Plot heatmap
# --------------------
def plot_heatmap(
    matrix_df,
    output_prefix,
    metric_name,
    cmap="viridis"
):

    linkage_matrix = linkage(
        matrix_df.values,
        method="average"
    )

    g = sns.clustermap(
        matrix_df,
        cmap=cmap,
        annot=True,
        fmt=".2f",
        linewidths=0.5,
        figsize=(10, 10),
        row_linkage=linkage_matrix,
        col_linkage=linkage_matrix
    )

    plt.title(metric_name)

    plt.savefig(
        f"{output_prefix}_{metric_name}.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.savefig(
        f"{output_prefix}_{metric_name}.pdf",
        bbox_inches="tight"
    )

    plt.savefig(
        f"{output_prefix}_{metric_name}.svg",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

# --------------------
# Main
# --------------------
def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "input_tsv"
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True
    )

    parser.add_argument(
        "--normalize",
        action="store_true"
    )

    parser.add_argument(
        "--binary",
        action="store_true"
    )

    parser.add_argument(
        "--min_cluster_size",
        type=int,
        default=1
    )

    args = parser.parse_args()

    # --------------------
    # Read input
    # --------------------
    df = pd.read_csv(
        args.input_tsv,
        sep="\t",
        low_memory=False
    )

    # --------------------
    # Build matrix
    # --------------------
    clone_matrix = build_clone_matrix(
        df,
        normalize=args.normalize,
        min_cluster_size=args.min_cluster_size,
        binary=args.binary
    )

    clone_matrix.to_csv(
        f"{args.output}_clone_matrix.tsv",
        sep="\t"
    )

    # --------------------
    # Shared lineage counts
    # --------------------
    shared_df = compute_shared_lineage_matrix(
        clone_matrix
    )

    shared_df.to_csv(
        f"{args.output}_shared_lineages.tsv",
        sep="\t"
    )

    plot_heatmap(
        shared_df,
        args.output,
        "SharedLineages",
        cmap="Blues"
    )

    # --------------------
    # Cosine
    # --------------------
    cosine_df = compute_cosine_similarity_matrix(
        clone_matrix
    )

    cosine_df.to_csv(
        f"{args.output}_cosine_similarity.tsv",
        sep="\t"
    )

    plot_heatmap(
        cosine_df,
        args.output,
        "CosineSimilarity"
    )

    # --------------------
    # Bray-Curtis
    # --------------------
    bray_df = compute_bray_curtis_matrix(
        clone_matrix
    )

    bray_df.to_csv(
        f"{args.output}_braycurtis.tsv",
        sep="\t"
    )

    plot_heatmap(
        bray_df,
        args.output,
        "BrayCurtisDissimilarity",
        cmap="magma_r"
    )

    # --------------------
    # Morisita-Horn
    # --------------------
    mh_df = compute_morisita_horn_matrix(
        clone_matrix
    )

    mh_df.to_csv(
        f"{args.output}_morisitahorn.tsv",
        sep="\t"
    )

    plot_heatmap(
        mh_df,
        args.output,
        "MorisitaHornSimilarity"
    )
 
# --------------------
# Run
# --------------------
if __name__ == "__main__":
    main()