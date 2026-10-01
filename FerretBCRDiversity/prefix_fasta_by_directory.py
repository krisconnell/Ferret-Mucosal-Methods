#!/usr/bin/env python3

import argparse
from pathlib import Path

FASTA_EXTENSIONS = {".fasta", ".fa", ".fna"}

def process_fasta(fasta_path: Path, prefix: str, inplace: bool):
    """
    Add prefix to FASTA headers if not already present.
    """
    if not prefix.endswith("|"):
        prefix = prefix + "|"

    output_path = fasta_path if inplace else fasta_path.with_suffix(fasta_path.suffix + ".prefixed")

    with fasta_path.open("r") as fin, output_path.open("w") as fout:
        for line in fin:
            if line.startswith(">"):
                header = line[1:].rstrip("\n")
                if header.startswith(prefix):
                    fout.write(line)  # already prefixed
                else:
                    fout.write(f">{prefix}{header}\n")
            else:
                fout.write(line)

def main():
    parser = argparse.ArgumentParser(
        description="Prefix FASTA headers using subdirectory names"
    )
    parser.add_argument(
        "parent_dir",
        help="Parent directory containing tissue subdirectories"
    )
    parser.add_argument(
        "--inplace",
        action="store_true",
        help="Modify FASTA files in place (default: write .prefixed files)"
    )

    args = parser.parse_args()
    parent_dir = Path(args.parent_dir)

    if not parent_dir.is_dir():
        raise ValueError(f"{parent_dir} is not a directory")

    for subdir in sorted(p for p in parent_dir.iterdir() if p.is_dir()):
        prefix = subdir.name
        fasta_files = [
            f for f in subdir.iterdir()
            if f.suffix.lower() in FASTA_EXTENSIONS
        ]

        if not fasta_files:
            continue

        for fasta in fasta_files:
            process_fasta(fasta, prefix, args.inplace)
            print(f"[OK] {fasta}  → prefix='{prefix}|'")

if __name__ == "__main__":
    main()