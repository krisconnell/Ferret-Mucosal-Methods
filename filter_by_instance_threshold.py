#!/usr/bin/env python3

import argparse
from collections import defaultdict
from Bio import SeqIO

def filter_clone_instances(input_fasta, output_fasta, threshold):
    clone_to_records = defaultdict(list)

    for record in SeqIO.parse(input_fasta, "fasta"):
        try:
            clone_part = record.id.split("_")[0]  # Clone4
            clone_id = clone_part.replace("Clone", "")  # 4
        except IndexError:
            print("Skipping malformed header: {}".format(record.id))
            continue

        clone_to_records[clone_id].append(record)

    kept_records = []
    for clone_id, records in clone_to_records.items():
        if len(records) >= threshold:
            kept_records.extend(records)

    SeqIO.write(kept_records, output_fasta, "fasta")
    print("✅ Kept {} clone instances (≥ {} per clone) in {}".format(
        len(kept_records), threshold, output_fasta))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter clone instances by minimum instance count threshold."
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file (from Step 1)")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA with filtered instances")
    parser.add_argument("-t", "--threshold", type=int, required=True, help="Minimum clone instance count")

    args = parser.parse_args()
    filter_clone_instances(args.input, args.output, args.threshold)
