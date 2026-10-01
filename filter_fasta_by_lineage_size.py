#!/usr/bin/env python3

import argparse
from collections import defaultdict
import os

def parse_fasta_by_lineage(input_file):
    lineage_seqs = defaultdict(list)
    with open(input_file, 'r') as f:
        current_header = ''
        current_seq = ''
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header and current_seq:
                    lineage_id = current_header.split('_')[2]
                    lineage_seqs[lineage_id].append((current_header, current_seq))
                current_header = line
                current_seq = ''
            else:
                current_seq += line
        if current_header and current_seq:
            lineage_id = current_header.split('_')[2]
            lineage_seqs[lineage_id].append((current_header, current_seq))
    return lineage_seqs

def write_filtered_fastas(lineage_seqs, thresholds, out_prefix):
    for t in thresholds:
        out_file = f"{out_prefix}_{t}plus.fasta"
        with open(out_file, 'w') as out:
            for lineage_id, seqs in lineage_seqs.items():
                if len(seqs) >= t:
                    for header, seq in seqs:
                        out.write(f"{header}\n{seq}\n")
        print(f"✅ Written: {out_file}")

def main():
    parser = argparse.ArgumentParser(description="Filter FASTA by lineage size thresholds.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-t", "--thresholds", nargs='+', type=int, required=True,
                        help="List of lineage size thresholds (e.g. 5 10 20)")
    parser.add_argument("-o", "--outprefix", default="filtered_lineages",
                        help="Prefix for output files (default: filtered_lineages)")

    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"❌ Input file not found: {args.input}")
        return

    print(f"🔍 Parsing FASTA file: {args.input}")
    lineage_seqs = parse_fasta_by_lineage(args.input)

    print(f"✏️ Writing filtered FASTA files for thresholds: {args.thresholds}")
    write_filtered_fastas(lineage_seqs, args.thresholds, args.outprefix)

if __name__ == "__main__":
    main()
