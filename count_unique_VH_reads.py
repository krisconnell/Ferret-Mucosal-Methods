#!/usr/bin/env python3
import argparse
import pandas as pd
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Count unique VH reads per clonotype from IgG/IgA cluster lists."
    )
    parser.add_argument(
        "--cluster_list", "-c",
        required=True,
        help="Path to text file listing clusters (e.g., IgG_clusters.txt)."
    )
    parser.add_argument(
        "--vh_database", "-v",
        required=True,
        help="Path to VH search database (TSV format). Must include a 'ClusterID' column."
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output TSV file path."
    )

    args = parser.parse_args()

    # --- Load cluster list ---
    try:
        with open(args.cluster_list, 'r') as f:
            clusters = [line.strip().split()[0] for line in f if line.strip()]
        clusters = list(set(clusters))  # remove duplicates
    except Exception as e:
        sys.exit(f"Error reading cluster list: {e}")

    print(f"Loaded {len(clusters)} unique cluster IDs from {args.cluster_list}")

    # --- Load VH search database ---
    try:
        vh_db = pd.read_csv(args.vh_database, sep='\t', dtype=str)
    except Exception as e:
        sys.exit(f"Error reading VH database: {e}")

    if "ClusterID" not in vh_db.columns:
        sys.exit("Error: VH database must contain a 'ClusterID' column.")

    print(f"Loaded VH database with {len(vh_db)} reads and {len(vh_db.columns)} columns.")

    # --- Filter VH DB for relevant clusters ---
    filtered_db = vh_db[vh_db["ClusterID"].isin(clusters)]

    # --- Count unique reads per cluster ---
    # Assuming each row is a unique read (so we just count rows per cluster)
    vh_counts = filtered_db.groupby("ClusterID").size().reset_index(name="Unique_VH_Reads")

    # --- Include clusters with zero reads ---
    all_clusters_df = pd.DataFrame({"ClusterID": clusters})
    result = pd.merge(all_clusters_df, vh_counts, on="ClusterID", how="left").fillna(0)
    result["Unique_VH_Reads"] = result["Unique_VH_Reads"].astype(int)

    # --- Save output ---
    result.to_csv(args.output, sep="\t", index=False)
    print(f"Output written to {args.output}")
    print(result.head())

if __name__ == "__main__":
    main()
