#!/usr/bin/env python3

import os
import re
import argparse
import pandas as pd


def extract_tissue_and_isotype(filename):
    """
    Extract tissue and isotype from filename.
    Example:
    T8575nec_BAL_IgA_S10_L001_001.assembled_IgBLAST_filtered.tsv
    """

    iso_match = re.search(r'_(IgA|IgG|IgM)_', filename)
    isotype = iso_match.group(1) if iso_match else "Unknown"

    tissue_match = re.search(r'T\d+[A-Za-z]*_([^_]+)_(IgA|IgG|IgM)_', filename)
    parts = filename.split("_")
    #tissue = tissue_match.group(1) if tissue_match else "Unknown"
    tissue = parts[1]

    return tissue, isotype


def main(parent_dir, output_dir, chunksize):
    os.makedirs(output_dir, exist_ok=True)

    # Track whether we've written headers per tissue
    written_headers = {}

    for root, dirs, files in os.walk(parent_dir):
        for file in files:
            if file.endswith("_IgBLAST_filtered.tsv"):

                full_path = os.path.join(root, file)
                tissue, isotype = extract_tissue_and_isotype(file)

                print(f"\nProcessing: {file}")
                print(f"  Tissue: {tissue}")
                print(f"  Isotype: {isotype}")

                out_file = os.path.join(output_dir, f"{tissue}_combined_isotypes.tsv")

                # Stream file in chunks
                for chunk in pd.read_csv(full_path, sep="\t", chunksize=chunksize):

                    chunk["Tissue"] = tissue
                    chunk["Isotype"] = isotype

                    # Write header only once per tissue
                    write_header = False
                    if tissue not in written_headers:
                        write_header = True
                        written_headers[tissue] = True

                    chunk.to_csv(
                        out_file,
                        sep="\t",
                        index=False,
                        mode="a",
                        header=write_header
                    )

                print(f"  Appended to: {out_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Memory-safe combination of IgBLAST isotypes per tissue."
    )
    parser.add_argument("parent_dir", help="Parent directory containing IgBLAST directories")
    parser.add_argument("-o", "--output", default="combined_by_tissue",
                        help="Output directory")
    parser.add_argument("--chunksize", type=int, default=100000,
                        help="Rows per chunk (default=100k)")

    args = parser.parse_args()
    main(args.parent_dir, args.output, args.chunksize)