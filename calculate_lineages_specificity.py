#!/usr/bin/env python3

import argparse
import pandas as pd
import os


def calculate_lineage_specificity(input_file, output_file):

    print(f"Reading: {input_file}")

    # Read TSV
    df = pd.read_csv(
        input_file,
        sep="\t",
        dtype={"ClusterID": str, "tissue": str}
    )

    # Check required columns
    required_columns = ["ClusterID", "tissue"]

    missing = [col for col in required_columns if col not in df.columns]

    if missing:
        raise ValueError(
            f"Missing required column(s): {', '.join(missing)}"
        )

    # Remove rows without a ClusterID or tissue
    df = df.dropna(subset=["ClusterID", "tissue"])

    # Remove accidental whitespace
    df["ClusterID"] = df["ClusterID"].astype(str).str.strip()
    df["tissue"] = df["tissue"].astype(str).str.strip()

    # Remove empty values after stripping
    df = df[
        (df["ClusterID"] != "") &
        (df["tissue"] != "")
    ]

    print(f"Total rows: {len(df):,}")
    print(f"Unique ClusterIDs: {df['ClusterID'].nunique():,}")
    print(f"Unique tissues: {df['tissue'].nunique():,}")

    # ---------------------------------------------------------
    # Determine which tissues contain each ClusterID
    #
    # Multiple rows for the same ClusterID/tissue are collapsed
    # here, so each lineage is counted only once per tissue.
    # ---------------------------------------------------------

    cluster_tissues = (
        df[["ClusterID", "tissue"]]
        .drop_duplicates()
        .groupby("ClusterID")["tissue"]
        .agg(set)
    )

    # ---------------------------------------------------------
    # Count total lineages and tissue-specific lineages
    # ---------------------------------------------------------

    results = []

    tissues = sorted(df["tissue"].unique())

    for tissue in tissues:

        # All unique ClusterIDs represented in this tissue
        tissue_clusters = set(
            df.loc[df["tissue"] == tissue, "ClusterID"]
        )

        total_lineages = len(tissue_clusters)

        # A lineage is tissue-specific if this is the ONLY
        # tissue in which that ClusterID occurs.
        tissue_specific_lineages = sum(
            1
            for cluster_id in tissue_clusters
            if cluster_tissues[cluster_id] == {tissue}
        )

        results.append({
            "tissue": tissue,
            "total_lineages": total_lineages,
            "tissue_specific_lineages": tissue_specific_lineages
        })

    # Create output dataframe
    results_df = pd.DataFrame(results)

    # Write CSV
    results_df.to_csv(output_file, index=False)

    print("\nResults:")
    print(results_df.to_string(index=False))

    print(f"\nOutput written to:")
    print(output_file)


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Calculate total and tissue-specific BCR lineages "
            "from a thresholded TSV file."
        )
    )

    parser.add_argument(
        "input_file",
        help="Input TSV containing ClusterID and tissue columns"
    )

    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output CSV filename"
    )

    args = parser.parse_args()

    # If no output name was provided, create one automatically
    if args.output is None:

        base = os.path.splitext(args.input_file)[0]

        args.output = base + "_lineage_counts.csv"

    calculate_lineage_specificity(
        args.input_file,
        args.output
    )


if __name__ == "__main__":
    main()