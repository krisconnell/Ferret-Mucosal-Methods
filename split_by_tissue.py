#!/usr/bin/env python3

import pandas as pd
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description="Split clustered TSV by tissue column.")
    parser.add_argument("input_tsv", help="Combined clustered TSV file")
    parser.add_argument("-o", "--outdir", default="tissue_split_clustered",
                        help="Output directory (default: tissue_split_clustered)")
    args = parser.parse_args()

    # Load file
    df = pd.read_csv(args.input_tsv, sep="\t")

    # Clean column names just in case
    df.columns = df.columns.str.strip()

    if "tissue" not in df.columns:
        raise ValueError("No 'tissue' column found in input TSV.")

    # Make output directory
    os.makedirs(args.outdir, exist_ok=True)

    # Split and write
    for tissue, sub_df in df.groupby("tissue"):
        safe_tissue = str(tissue).replace(" ", "_")
        outfile = os.path.join(args.outdir, f"{safe_tissue}_clustered.tsv")
        sub_df.to_csv(outfile, sep="\t", index=False)
        print(f"Wrote: {outfile} ({len(sub_df)} rows)")

if __name__ == "__main__":
    main()