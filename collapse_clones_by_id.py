#!/usr/bin/env python3

import argparse
from collections import defaultdict
from Bio import SeqIO

def collapse_by_clone_id(input_fasta, output_fasta):
    clone_to_total_copies = defaultdict(int)
    clone_to_best_seq = {}

    for record in SeqIO.parse(input_fasta, "fasta"):
        try:
            id_parts = record.id.split("_")
            clone_id = id_parts[0].replace("Clone", "")  # 4
            copy_count = int(id_parts[1])                # 515
        except (IndexError, ValueError):
            print("⚠️ Skipping malformed header: {}".format(record.id))
            continue

        clone_to_total_copies[clone_id] += copy_count

        if clone_id not in clone_to_best_seq:
            clone_to_best_seq[clone_id] = record
        else:
            if len(record.seq) > len(clone_to_best_seq[clone_id].seq):
                clone_to_best_seq[clone_id] = record

    with open(output_fasta, "w") as out_f:
        for clone_id, record in clone_to_best_seq.items():
            total_copies = clone_to_total_copies[clone_id]
            record.id = "Clone{}_Copies{}".format(clone_id, total_copies)
            record.description = ""
            SeqIO.write(record, out_f, "fasta")

    print("✅ Collapsed to {} clones in {}".format(len(clone_to_best_seq), output_fasta))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Collapse clone instances into one per clone ID by summing copy counts and keeping the longest sequence."
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA from thresholded clone instances")
    parser.add_argument("-o", "--output", required=True, help="Output collapsed FASTA")

    args = parser.parse_args()
    collapse_by_clone_id(args.input, args.output)
