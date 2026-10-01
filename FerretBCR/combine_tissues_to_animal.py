#take in directory from combine_isotypes_tissuestreaming.py to make a donor-specific IgBLAST database

#!/usr/bin/env python3

#Example argument:
#python3 combine_tissues_to_animal.py tissue_combined \
#    -o T8575nec_combined_BCR.tsv

import os
import re
import argparse
import pandas as pd


def extract_donor_from_filename(filename):
    """
    Extract donor ID from filename if present.
    Example:
    T8575nec_BAL_IgA_... → T8575nec
    """
    match = re.search(r'(T\d+[A-Za-z]*)', filename)
    return match.group(1) if match else "Unknown"


def main(tissue_dir, output_file, chunksize):

    written_header = False
    donor_id = None

    for file in os.listdir(tissue_dir):

        if file.endswith("combined_isotypes.tsv"):

            full_path = os.path.join(tissue_dir, file)

            print(f"\nProcessing tissue file: {file}")

            # Extract donor from filename if possible
            if donor_id is None:
                donor_id = extract_donor_from_filename(file)

            # Stream read
            for chunk in pd.read_csv(full_path, sep="\t", chunksize=chunksize):

                # Ensure Tissue column exists (it should)
                if "Tissue" not in chunk.columns:
                    raise ValueError(f"Tissue column missing in {file}")

                # Add Donor column
                chunk["Donor"] = donor_id

                chunk.to_csv(
                    output_file,
                    sep="\t",
                    index=False,
                    mode="a",
                    header=not written_header
                )

                written_header = True

            print(f"  Appended to {output_file}")

    print("\nFinished combining all tissues.")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Combine tissue-level BCR files into one animal-level file."
    )

    parser.add_argument("tissue_dir",
                        help="Directory containing *_combined_isotypes.tsv files")

    parser.add_argument("-o", "--output",
                        required=True,
                        help="Output animal-level BCR file (e.g., T8575nec_combined_BCR.tsv)")

    parser.add_argument("--chunksize",
                        type=int,
                        default=100000,
                        help="Chunk size for streaming (default=100k)")

    args = parser.parse_args()

    main(args.tissue_dir, args.output, args.chunksize)