#for use with large clustered BCR databases that have tissue- and isotype-specific data. 
# extracts tissue-specific VH databases for downstream analysis 
#and also summarizes clusters 
#KMC 2/20/2026

#!/usr/bin/env python3

import os
import argparse
import pandas as pd
from collections import defaultdict


def main(clustered_file, output_dir, chunksize):

    os.makedirs(output_dir, exist_ok=True)

    written_headers = {}
    sequence_counts = defaultdict(int)
    clone_sets = defaultdict(set)

    print(f"Processing: {clustered_file}")

    for chunk in pd.read_csv(clustered_file, sep="\t", chunksize=chunksize):

        if "Tissue" not in chunk.columns:
            raise ValueError("Missing 'Tissue' column.")
        if "ClusterID" not in chunk.columns:
            raise ValueError("Missing 'ClusterID' column.")

        # Update summary stats
        for tissue, sub_df in chunk.groupby("Tissue"):

            sequence_counts[tissue] += len(sub_df)

            clone_sets[tissue].update(sub_df["ClusterID"].unique())

            # Write tissue-specific file
            out_file = os.path.join(output_dir, f"{tissue}_clustered.tsv")

            write_header = False
            if tissue not in written_headers:
                write_header = True
                written_headers[tissue] = True

            sub_df.to_csv(
                out_file,
                sep="\t",
                index=False,
                mode="a",
                header=write_header
            )

    print("Finished splitting. Now computing summaries...")

    # ----------------------------
    # Clone count summary
    # ----------------------------
    summary_rows = []
    for tissue in clone_sets:
        summary_rows.append({
            "Tissue": tissue,
            "Total_Sequences": sequence_counts[tissue],
            "Unique_Clones": len(clone_sets[tissue])
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(
        os.path.join(output_dir, "clone_count_summary.tsv"),
        sep="\t",
        index=False
    )

    # ----------------------------
    # Cross-tissue overlap matrix
    # ----------------------------
    tissues = sorted(clone_sets.keys())
    overlap_matrix = pd.DataFrame(index=tissues, columns=tissues)

    for t1 in tissues:
        for t2 in tissues:
            shared = len(clone_sets[t1].intersection(clone_sets[t2]))
            overlap_matrix.loc[t1, t2] = shared

    overlap_matrix.to_csv(
        os.path.join(output_dir, "clone_overlap_matrix.tsv"),
        sep="\t"
    )

    print("All summaries written.")
    print("Done.")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Split clustered BCR file by tissue and compute clone summaries."
    )

    parser.add_argument("clustered_file",
                        help="Clustered BCR file containing Tissue and ClusterID columns")

    parser.add_argument("-o", "--output",
                        default="clustered_by_tissue",
                        help="Output directory")

    parser.add_argument("--chunksize",
                        type=int,
                        default=100000,
                        help="Chunk size (default=100k)")

    args = parser.parse_args()

    main(args.clustered_file, args.output, args.chunksize)