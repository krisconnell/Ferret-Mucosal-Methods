#!/usr/bin/env python3
#Use this after IgBLAST annotation to combine different isotype VH libraries into one library per tissue. 
#KMC 02/20/2026 
#Note I used ChatGPT to help write this code
#python3 
import os
import argparse

def extract_isotype(dirname):
    """
    Extract IgA / IgG / IgM etc from directory name.
    Example:
    T8575nec_ttBAL_IgG-ds.2828857a7bf845f6a52a6abc90f9f723
    """
    parts = dirname.split("_")
    for p in parts:
        if p.startswith("Ig"):
            return p.split("-")[0]  # IgG-ds -> IgG
    return "Unknown"

def extract_tissue(dirname):
    """
    Extract tissue from directory name.
    Example:
    T8575nec_ttBAL_IgG-ds -> tissue = ttBAL
    """
    parts = dirname.split("_")
    for p in parts:
        if p.startswith("tt") or p.startswith("LN") or p.startswith("BAL"):  # add more if needed
            return p
    return "Unknown"

def combine_files(input_root, output_file):

    if os.path.exists(output_file):
        print(f"Removing existing output file: {output_file}")
        os.remove(output_file)

    header_written = False
    reference_header = None

    with open(output_file, "w") as outfile:

        for dirname in sorted(os.listdir(input_root)):

            dirpath = os.path.join(input_root, dirname)

            if not os.path.isdir(dirpath):
                continue

            if "_Ig" not in dirname:
                continue

            isotype = extract_isotype(dirname)
            tissue = extract_tissue(dirname)

            for file in os.listdir(dirpath):

                if (
                    file.endswith(".tsv")
                    and "filtered" in file
                    and "clustered" not in file
                ):

                    filepath = os.path.join(dirpath, file)
                    print(f"Processing: {filepath}")

                    with open(filepath, "r") as infile:

                        header = infile.readline().rstrip("\n")

                        if not header_written:
                            reference_header = header
                            outfile.write(header + "\tisotype\ttissue\n")
                            header_written = True
                        else:
                            if header != reference_header:
                                raise ValueError(
                                    f"Header mismatch detected in {filepath}"
                                )

                        for line in infile:
                            outfile.write(
                                line.rstrip("\n") + f"\t{isotype}\t{tissue}\n"
                            )

    print("Done.")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_root", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    combine_files(args.input_root, args.output)