#!/usr/bin/env python3

import argparse
import subprocess
import time
from pathlib import Path
from datetime import datetime

THRESHOLDS = [1, 5, 10] #edit this as needed

# Only accept prefixed FASTA files
VALID_SUFFIXES = (
    ".fasta.prefixed",
    ".fa.prefixed",
    ".fna.prefixed"
)

def run_step1(input_fasta, output_fasta):
    print(f"\n[Step 1] Collapsing clone instances -> {output_fasta}")
    subprocess.run([
        "python3", "collapse_clone_instances.py",
        "-i", str(input_fasta),
        "-o", str(output_fasta)
    ], check=True)

def run_step2(input_fasta, threshold, output_fasta):
    print(f"[Step 2] Threshold >= {threshold} -> {output_fasta}")
    subprocess.run([
        "python3", "filter_clones_by_instance_count.py",
        "-i", str(input_fasta),
        "-t", str(threshold),
        "-o", str(output_fasta)
    ], check=True)

def run_step3(input_fasta, output_fasta, sample_name):
    print(f"[Step 3] Final collapse -> {output_fasta}")
    subprocess.run([
        "python3", "collapse_filtered_clones.py",
        "-i", str(input_fasta),
        "-o", str(output_fasta),
        "-s", sample_name
    ], check=True)

def count_clones(fasta_file):
    clone_ids = set()
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                clone = line.split("|", 1)[1] if "|" in line else line.strip()
                clone_ids.add(clone)
    return len(clone_ids)

def main():
    parser = argparse.ArgumentParser(
        description="Run clone threshold pipeline per tissue directory (prefixed FASTA only)"
    )
    parser.add_argument(
        "parent_dir",
        help="Parent directory containing tissue subdirectories"
    )
    args = parser.parse_args()

    parent_dir = Path(args.parent_dir)
    if not parent_dir.is_dir():
        raise ValueError(f"{parent_dir} is not a directory")

    start_time = time.time()
    print("\n🚀 Starting directory-level clone threshold pipeline")
    print(f"Parent directory: {parent_dir}")
    print(f"Thresholds: {THRESHOLDS}")
    print("Input constraint: *.fasta.prefixed / *.fa.prefixed / *.fna.prefixed")

    for tissue_dir in sorted(d for d in parent_dir.iterdir() if d.is_dir()):
        tissue = tissue_dir.name

        prefixed_fastas = [
            f for f in tissue_dir.iterdir()
            if f.name.endswith(VALID_SUFFIXES)
        ]

        if not prefixed_fastas:
            print(f"\n⏭️  Skipping {tissue}: no prefixed FASTA found")
            continue

        if len(prefixed_fastas) > 1:
            raise RuntimeError(
                f"{tissue}: multiple prefixed FASTA files found:\n" +
                "\n".join(f"  - {f.name}" for f in prefixed_fastas)
            )

        input_fasta = prefixed_fastas[0]

        print(f"\n🧬 Processing tissue: {tissue}")
        print(f"Input FASTA: {input_fasta.name}")

        summary_lines = ["Threshold\tClonesKept"]

        step1_fasta = tissue_dir / f"{tissue}_step1_collapsed.fasta"
        run_step1(input_fasta, step1_fasta)

        for threshold in THRESHOLDS:
            print(f"\n--- {tissue}: threshold {threshold} ---")

            step2_fasta = tissue_dir / f"{tissue}_step2_thresh{threshold}.fasta"
            final_fasta = tissue_dir / f"{tissue}_threshold{threshold}.fasta"

            run_step2(step1_fasta, threshold, step2_fasta)
            run_step3(step2_fasta, final_fasta, tissue)

            n_clones = count_clones(final_fasta)
            summary_lines.append(f"{threshold}\t{n_clones}")

            print(f"✅ {tissue} threshold {threshold}: {n_clones} clones")

        summary_file = tissue_dir / f"{tissue}_threshold_summary.tsv"
        with open(summary_file, "w") as f:
            f.write("\n".join(summary_lines) + "\n")

        print(f"📁 Summary written to {summary_file}")

    elapsed = time.time() - start_time
    print(f"\n🕒 All tissues completed in {elapsed:.2f} seconds.")

if __name__ == "__main__":
    main()