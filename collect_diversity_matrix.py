#!/usr/bin/env python3

import argparse
import os
import glob
import pandas as pd


def main():

    parser = argparse.ArgumentParser(
        description="Collect diversity values across tissues for a specified clone threshold."
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Root directory containing tissue subdirectories."
    )

    parser.add_argument(
        "-t",
        "--threshold",
        required=True,
        type=int,
        help="Clone threshold to extract."
    )

    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output TSV file."
    )

    args = parser.parse_args()

    root = args.input
    threshold = args.threshold

    results = {}

    for tissue in sorted(os.listdir(root)):

        tissue_dir = os.path.join(root, tissue)

        if not os.path.isdir(tissue_dir):
            continue

        results[tissue] = {}

        for order in range(6):

            pattern = os.path.join(
                tissue_dir,
                f"*threshold{threshold}_r{order}_diversity*"
            )

            matches = glob.glob(pattern)

            if len(matches) == 0:
                results[tissue][f"Order{order}"] = None
                continue

            filename = matches[0]

            try:
                df = pd.read_csv(filename)

                if "diversity" not in df.columns:
                    raise ValueError(
                        f"'diversity' column not found in {filename}"
                    )

                diversity = df.iloc[0]["diversity"]

                results[tissue][f"Order{order}"] = diversity

            except Exception as e:
                print(f"Could not read {filename}: {e}")
                results[tissue][f"Order{order}"] = None

    matrix = pd.DataFrame.from_dict(results, orient="index")

    matrix.index.name = "Tissue"

    matrix = matrix[
        [f"Order{i}" for i in range(6)]
    ]

    if args.output is None:
        output = os.path.join(
            root,
            f"threshold{threshold}_diversity_matrix.tsv"
        )
    else:
        output = args.output

    matrix.to_csv(output, sep="\t")

    print(f"\nSaved matrix to:")
    print(output)

    print("\nMatrix:")
    print(matrix)


if __name__ == "__main__":
    main()