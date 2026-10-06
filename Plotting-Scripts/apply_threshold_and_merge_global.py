#!/usr/bin/env python3

import pandas as pd
import argparse
from pathlib import Path
import sys
from datetime import datetime

##########################################
# Log output to both terminal and file
##########################################

class Tee:

    def __init__(self, *files):
        self.files = files

    def write(self, message):
        for f in self.files:
            f.write(message)
            f.flush()

    def flush(self):
        for f in self.files:
            f.flush()

##########################################
# Parse isotype from c_call
##########################################

def extract_isotypes(c_call):

    if pd.isna(c_call):
        return []

    parts = str(c_call).upper().split(",")

    isos = set()

    for p in parts:
        if p.startswith("IGHG"):
            isos.add("IgG")
        elif p.startswith("IGHA"):
            isos.add("IgA")
        elif p.startswith("IGHM"):
            isos.add("IgM")

    return list(isos)


##########################################
# Find clustered files
##########################################

def find_clustered_files(parent):

    files = []

    for tissue_dir in parent.iterdir():

        if not tissue_dir.is_dir():
            continue

        for f in tissue_dir.glob("*_clustered.tsv"):

            files.append((tissue_dir.name, f))

    return files


##########################################
# Main
##########################################

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Apply a global ClusterID abundance threshold across "
            "all tissues, then split the surviving cells by isotype."
        )
    )

    parser.add_argument(
        "parent_directory",
        help="Parent directory containing tissue directories"
    )

    parser.add_argument(
        "-t",
        "--threshold",
        type=int,
        help="Minimum global number of cells required for a ClusterID"
    )

    parser.add_argument(
        "--threshold_file",
        default="recommended_global_threshold.txt",
        help="File containing the threshold if --threshold is not supplied"
    )

    parser.add_argument(
        "-o",
        "--outdir",
        default="thresholded_output",
        help="Output directory"
    )

    parser.add_argument(
        "--cluster_col",
        default="ClusterID",
        help="Column containing the cluster identifier"
    )

    args = parser.parse_args()

    parent = Path(args.parent_directory)
    outdir = Path(args.outdir)

    outdir.mkdir(exist_ok=True)

    ##########################################
    #   Create timestamped log
    ##########################################

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    log_file = outdir / f"thresholding_run_{timestamp}.log"

    log_handle = open(log_file, "w")

    sys.stdout = Tee(sys.__stdout__, log_handle)
    sys.stderr = Tee(sys.__stderr__, log_handle)

    print("=" * 60)
    print("THRESHOLDING RUN")
    print("=" * 60)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Parent directory: {parent}")
    print(f"Output directory: {outdir}")
    print(f"Log file: {log_file}")
    print()

    ##########################################
    # Get threshold
    ##########################################

    if args.threshold is not None:

        threshold = args.threshold

    else:

        threshold_path = parent / args.threshold_file

        if not threshold_path.exists():
            raise FileNotFoundError(
                f"Missing threshold file: {threshold_path}"
            )

        threshold = int(threshold_path.read_text().strip())

    print()
    print("=" * 60)
    print("GLOBAL CLUSTER THRESHOLDING")
    print("=" * 60)
    print(f"Using global threshold: {threshold}")
    print(f"Cluster column: {args.cluster_col}")
    print()


    ##########################################
    # Find files
    ##########################################

    files = find_clustered_files(parent)

    if not files:
        raise RuntimeError(
            f"No *_clustered.tsv files found under {parent}"
        )

    print(f"Found {len(files)} tissue files.")
    print()


    ##########################################
    # Read ALL tissues
    ##########################################

    all_data = []

    for tissue, filepath in files:

        print(f"Reading {tissue}: {filepath.name}")

        df = pd.read_csv(filepath, sep="\t")

        if args.cluster_col not in df.columns:
            raise KeyError(
                f"{args.cluster_col} not found in {filepath}"
            )

        if "c_call" not in df.columns:
            raise KeyError(
                f"'c_call' not found in {filepath}"
            )

        df["tissue"] = tissue

        # Determine isotypes now, but DO NOT threshold by isotype.
        df["isotypes"] = df["c_call"].apply(extract_isotypes)

        all_data.append(df)

        print(f"  Rows: {len(df):,}")

    ##########################################
    # Combine all tissues
    ##########################################

    print()
    print("=" * 60)
    print("COMBINING ALL TISSUES")
    print("=" * 60)

    all_df = pd.concat(
        all_data,
        ignore_index=True
    )

    print(f"Total rows: {len(all_df):,}")

    ##########################################
    # GLOBAL CLUSTER COUNTS
    ##########################################

    print()
    print("=" * 60)
    print("CALCULATING GLOBAL CLUSTER COUNTS")
    print("=" * 60)

    global_counts = all_df[args.cluster_col].value_counts()

    total_clusters = len(global_counts)

    print(f"Global unique clusters: {total_clusters:,}")

    ##########################################
    # Apply GLOBAL threshold
    ##########################################

    keep_clusters = global_counts[
        global_counts >= threshold
    ].index

    filtered_all = all_df[
        all_df[args.cluster_col].isin(keep_clusters)
    ].copy()

    retained_clusters = len(keep_clusters)

    print()
    print(f"Clusters before threshold: {total_clusters:,}")
    print(f"Clusters after threshold:  {retained_clusters:,}")

    print(
        f"Clusters retained: "
        f"{retained_clusters / total_clusters * 100:.1f}%"
    )

    print()

    print(f"Rows before threshold: {len(all_df):,}")
    print(f"Rows after threshold:  {len(filtered_all):,}")

    print(
        f"Rows retained: "
        f"{len(filtered_all) / len(all_df) * 100:.1f}%"
    )

    filtered_all.to_csv(
                outdir / f"combined_threshold{threshold}_all_isotypes.tsv",
                sep="\t",
                index=False
            )


    ##########################################
    # Save global cluster counts
    ##########################################

    cluster_summary = pd.DataFrame({
        args.cluster_col: global_counts.index,
        "global_cell_count": global_counts.values,
        "passes_threshold": global_counts.values >= threshold
    })

    cluster_summary.to_csv(
        outdir / "global_cluster_counts.tsv",
        sep="\t",
        index=False
    )

    print()
    print(
        f"Wrote: {outdir / 'global_cluster_counts.tsv'}"
    )


    ##########################################
    # Split surviving cells by isotype
    ##########################################

    print()
    print("=" * 60)
    print("SPLITTING SURVIVING CELLS BY ISOTYPE")
    print("=" * 60)

    combined_by_iso = {
        "IgG": [],
        "IgA": [],
        "IgM": []
    }

    tissue_summaries = []

    for tissue, filepath in files:

        print()
        print(f"Processing output for {tissue}")

        tissue_original = all_df[
            all_df["tissue"] == tissue
        ]

        tissue_filtered = filtered_all[
            filtered_all["tissue"] == tissue
        ]

        for iso in ["IgG", "IgA", "IgM"]:

            ##################################
            # Original rows for this isotype
            ##################################

            original_iso = tissue_original[
                tissue_original["isotypes"].apply(
                    lambda x: iso in x
                )
            ]

            ##################################
            # Globally thresholded rows
            ##################################

            filtered_iso = tissue_filtered[
                tissue_filtered["isotypes"].apply(
                    lambda x: iso in x
                )
            ].copy()

            ##################################
            # Cluster counts
            ##################################

            original_clusters = (
                original_iso[args.cluster_col].nunique()
            )

            filtered_clusters = (
                filtered_iso[args.cluster_col].nunique()
            )

            ##################################
            # Add isotype
            ##################################

            if not filtered_iso.empty:
                filtered_iso["isotype"] = iso

            ##################################
            # Save per-tissue output
            ##################################

            if not filtered_iso.empty:

                out_file = (
                    outdir
                    / f"{tissue}_{iso}_threshold{threshold}.tsv"
                )

                filtered_iso.to_csv(
                    out_file,
                    sep="\t",
                    index=False
                )

                combined_by_iso[iso].append(filtered_iso)

            ##################################
            # Summary
            ##################################

            original_rows = len(original_iso)
            filtered_rows = len(filtered_iso)

            row_retention = (
                filtered_rows / original_rows * 100
                if original_rows > 0
                else 0
            )

            cluster_retention = (
                filtered_clusters / original_clusters * 100
                if original_clusters > 0
                else 0
            )

            print(f"  {iso}:")
            print(
                f"    Rows: "
                f"{original_rows:,} → {filtered_rows:,} "
                f"({row_retention:.1f}% retained)"
            )

            print(
                f"    Clusters: "
                f"{original_clusters:,} → "
                f"{filtered_clusters:,} "
                f"({cluster_retention:.1f}% retained)"
            )

            tissue_summaries.append({
                "tissue": tissue,
                "isotype": iso,
                "rows_before": original_rows,
                "rows_after": filtered_rows,
                "row_retention_percent": row_retention,
                "clusters_before": original_clusters,
                "clusters_after": filtered_clusters,
                "cluster_retention_percent": cluster_retention,
                "global_threshold": threshold
            })


    ##########################################
    # Write tissue summary
    ##########################################

    tissue_summary_df = pd.DataFrame(
        tissue_summaries
    )

    tissue_summary_df.to_csv(
        outdir / "tissue_isotype_threshold_summary.tsv",
        sep="\t",
        index=False
    )

    print()
    print(
        f"Wrote: "
        f"{outdir / 'tissue_isotype_threshold_summary.tsv'}"
    )


    ##########################################
    # Write combined isotype outputs
    ##########################################

    print()
    print("=" * 60)
    print("WRITING COMBINED ISOTYPE OUTPUTS")
    print("=" * 60)

    for iso, dfs in combined_by_iso.items():

        if not dfs:
            print(f"No surviving {iso} rows.")
            continue

        combined = pd.concat(
            dfs,
            ignore_index=True
        )

        out_file = (
            outdir
            / f"combined_{iso}_threshold{threshold}.tsv"
        )

        combined.to_csv(
            out_file,
            sep="\t",
            index=False
        )

        print(
            f"{iso}: {len(combined):,} rows → {out_file}"
        )


    ##########################################
    # Global summary
    ##########################################

    global_summary = pd.DataFrame([{
        "global_threshold": threshold,

        "rows_before": len(all_df),
        "rows_after": len(filtered_all),

        "row_retention_percent":
            len(filtered_all) / len(all_df) * 100,

        "clusters_before": total_clusters,
        "clusters_after": retained_clusters,

        "cluster_retention_percent":
            retained_clusters / total_clusters * 100
    }])

    global_summary.to_csv(
        outdir / "global_threshold_summary.tsv",
        sep="\t",
        index=False
    )

    print()
    print(
        f"Wrote: "
        f"{outdir / 'global_threshold_summary.tsv'}"
    )


    ##########################################
    # Done
    ##########################################

    print()
    print("=" * 60)
    print("DONE")
    print("=" * 60)
    print()


if __name__ == "__main__":
    main()
