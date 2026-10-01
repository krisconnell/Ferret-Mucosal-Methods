#!/usr/bin/env python3
import argparse
from collections import defaultdict
from Bio import SeqIO

def get_clone_name_from_instance(clone_instance):
    # clone_instance like "4_515"
    return clone_instance.split('_')[0]

def main():
    parser = argparse.ArgumentParser(description="Filter clones by minimum number of clone instances")
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file (collapsed instances)')
    parser.add_argument('-t', '--threshold', type=int, required=True, help='Minimum clone instance count threshold')
    parser.add_argument('-o', '--output', required=True, help='Output filtered FASTA file')
    args = parser.parse_args()

    # Count how many clone instances per clone name
    clone_instance_counts = defaultdict(list)  # clone_name -> list of SeqRecords

    for rec in SeqIO.parse(args.input, "fasta"):
        parts = rec.id.split('|')
        if len(parts) != 3:
            raise ValueError(f"Unexpected header format: {rec.id}")
        sample_name, clone_instance, read_count = parts
        clone_name = get_clone_name_from_instance(clone_instance)
        clone_instance_counts[clone_name].append(rec)

    # Filter by clone instance count threshold
    filtered_recs = []
    for clone_name, recs in clone_instance_counts.items():
        if len(recs) >= args.threshold:
            filtered_recs.extend(recs)

    # Write filtered FASTA (headers unchanged)
    with open(args.output, "w") as out_handle:
        SeqIO.write(filtered_recs, out_handle, "fasta")

if __name__ == "__main__":
    main()
