#!/usr/bin/env python3
import argparse
from collections import defaultdict
from Bio import SeqIO

def get_clone_name_from_instance(clone_instance):
    return clone_instance.split('_')[0]

def main():
    parser = argparse.ArgumentParser(description="Collapse filtered clones by summing read counts")
    parser.add_argument('-i', '--input', required=True, help='Input filtered FASTA file (clone instances)')
    parser.add_argument('-o', '--output', required=True, help='Output collapsed FASTA file')
    parser.add_argument('-s', '--sample', required=True, help='Sample name')
    args = parser.parse_args()

    clones = defaultdict(lambda: {"seq_record": None, "total_count": 0})

    for rec in SeqIO.parse(args.input, "fasta"):
        parts = rec.id.split('|')
        if len(parts) != 3:
            raise ValueError(f"Unexpected header format: {rec.id}")
        sample_name, clone_instance, read_count_str = parts
        clone_name = get_clone_name_from_instance(clone_instance)
        try:
            read_count = int(read_count_str)
        except:
            raise ValueError(f"Cannot parse read count in {rec.id}")

        clones[clone_name]["total_count"] += read_count
        current_rec = clones[clone_name]["seq_record"]
        if current_rec is None or len(rec.seq) > len(current_rec.seq):
            clones[clone_name]["seq_record"] = rec

    with open(args.output, "w") as out_handle:
        for clone_name, data in clones.items():
            rec = data["seq_record"]
            rec.id = f"{args.sample}|{clone_name}|{data['total_count']}"
            rec.description = ""
            SeqIO.write(rec, out_handle, "fasta")

if __name__ == "__main__":
    main()
