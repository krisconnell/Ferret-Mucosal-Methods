#!/usr/bin/env python3

# how to run: 
# tmux new -s diversity
# cd /stor/work/Georgiou/kmc5996/Scripts/Bcelldiversity
# python3 run_diversity_parallel.py ../../Ferret_Mucosal_Immunity/.../T8575_tissue_clustered_files
# you can also run a dry run beforehand so you dont accidentally run a several-hour long process with no results 
#!/usr/bin/env python3

import argparse
import subprocess
from pathlib import Path
from datetime import datetime
import time
import re

# ----------------------------
# FASTA utilities
# ----------------------------

def get_total_reads(fasta_path: Path) -> int:
    """
    Sum read counts from FASTA headers of the form:
    >Tissue|something|READCOUNT
    """
    total = 0
    with fasta_path.open() as f:
        for line in f:
            if line.startswith(">"):
                parts = line.strip().split("|")
                if len(parts) >= 3:
                    try:
                        total += int(parts[-1])
                    except ValueError:
                        pass
    return total


def compute_I_value(total_reads: int, n_points: int) -> str:
    """
    Compute diversity -I argument as:
    start interval end
    """
    start = max(1, total_reads // (n_points * 5))
    end = total_reads
    interval = max(1, (end - start) // (n_points - 1))
    return f"{start} {interval} {end}"

def run_diversity(
    tissue,
    threshold,
    fasta_path,
    order,
    n_points,
    dry_run,
    log_file,
    stack_yaml
):

    total_reads = get_total_reads(fasta_path)

    I_value = compute_I_value(total_reads, n_points)

    out_prefix = fasta_path.parent / f"{tissue}_{threshold}_r{order}"

    diversity_file = out_prefix.with_name(out_prefix.name + "_diversity.csv")

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # ----------------------------
    # Skip completed
    # ----------------------------

    if diversity_file.exists():

        log_file.write(
            f"[{timestamp}] ⏭️ SKIP {tissue} {threshold} r={order} (diversity exists)\n"
        )
        log_file.flush()
        return

    cmd = [
        "stack",
        "--stack-yaml", str(stack_yaml),
        "exec", "diversity", "--",
        "-i", str(fasta_path),
        "-o", str(diversity_file),   # diversity output
        "-r", str(order),
        "-R", "5",
        "-S", "1",
        "-C", "3",
        "-a",
        "-I", I_value,
        "-d",

        "-l", tissue
    ]

    log_file.write(
        f"[{timestamp}] ▶ START {tissue} {threshold} r={order} | reads={total_reads}\n"
    )
    log_file.write(f"[{timestamp}] CMD: {' '.join(cmd)}\n")
    log_file.flush()

    if dry_run:
        log_file.write(f"[{timestamp}] 🟡 DRY RUN\n")
        log_file.flush()
        return

    start_time = time.time()

    try:
        subprocess.run(cmd, check=True)

        duration = round(time.time() - start_time, 2)

        log_file.write(
            f"[{timestamp}] ✅ DONE {tissue} {threshold} r={order} | {duration}s\n"
        )

    except subprocess.CalledProcessError:

        log_file.write(
            f"[{timestamp}] ❌ ERROR {tissue} {threshold} r={order}\n"
        )

    log_file.flush()


def main():
    parser = argparse.ArgumentParser(
        description="Run diversity on all threshold FASTA files per tissue directory"
    )
    parser.add_argument(
        "parent_dir",
        help="Parent directory containing tissue subdirectories"
    )
    parser.add_argument(
        "-n", "--n-subsamples",
        type=int,
        default=40,
        help="Number of rarefaction subsampling points (default: 40)"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing"
    )
    parser.add_argument(
        "--stack-yaml",
        default="/stor/work/Georgiou/kmc5996/Scripts/Bcelldiversity/diversity/stack.yaml",
        help="Path to diversity stack.yaml"
    )
    parser.add_argument(
        "--orders",
        default="0,1,2,3,4,5",
        help="Hill orders to compute (default: 0,1,2,3,4,5)"
    )

    args = parser.parse_args()
    parent_dir = Path(args.parent_dir).resolve()
    stack_yaml = Path(args.stack_yaml).resolve()

    orders = [int(x) for x in args.orders.split(",")]

    log_path = parent_dir / "diversity_run.log"

    print("\n🚀 Starting diversity analysis")
    print(f"Parent directory: {parent_dir}")
    print(f"Subsampling points: {args.n_subsamples}")
    print(f"Dry run: {args.dry_run}\n")
    print("Orders:", orders)

    

    with log_path.open("w") as log_file:
        log_file.write("# Diversity batch run log\n\n")

        for tissue_dir in sorted(p for p in parent_dir.iterdir() if p.is_dir()):
            tissue = tissue_dir.name

            fasta_files = sorted(tissue_dir.glob(f"{tissue}_step3_thresh3_collapsed.fasta"))

            if not fasta_files:
                continue

            threshold = "3"

            for fasta_path in fasta_files: 

                for order in orders:

                    run_diversity(
                        tissue,
                        threshold,
                        fasta_path,
                        order,
                        args.n_subsamples,
                        args.dry_run,
                        log_file,
                        stack_yaml
                    )
                
        log_file.write("\n🎉 All tissues / thresholds complete\n")

if __name__ == "__main__":
    main()