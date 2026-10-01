#!/usr/bin/env python3
import argparse
from collections import defaultdict
from Bio import SeqIO

def parse_clone_instance_from_header(header):
    # Header: SampleName|full_clone_info
    # e.g. LN1|BCRseq_D8kLnhNo_4_515_ASLLRSGY_ISSDGSST_GFTFSNYY_...
    parts = header.split('|')
    if len(parts) != 2:
        raise ValueError(f"Unexpected header format: {header}")
    sample_name = parts[0]
    clone_info = parts[1]

    # Extract clone instance name and read count from clone_info
    # e.g. BCRseq_D8kLnhNo_4_515_... => clone_instance = "4_515"
    clone_parts = clone_info.split('_')
    if len(clone_parts) < 4:
        raise ValueError(f"Unexpected clone info format: {clone_info}")
    clone_name = clone_parts[2]
    read_count_str = clone_parts[3]
    clone_instance = f"{clone_name}_{read_count_str}"

    # Try convert read count to int
    try:
        read_count = int(read_count_str)
    except Exception as e:
        raise ValueError(f"Cannot parse read count as int from {read_count_str} in {header}")

    return sample_name, clone_instance, read_count

def main():
    parser = argparse.ArgumentParser(description="Collapse clone instances to longest sequences")
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output collapsed FASTA file')
    args = parser.parse_args()

    longest_seqs = dict()  # clone_instance -> (SeqRecord, read_count, sample_name)

    for rec in SeqIO.parse(args.input, "fasta"):
        sample_name, clone_instance, read_count = parse_clone_instance_from_header(rec.id)
        if clone_instance not in longest_seqs:
            longest_seqs[clone_instance] = (rec, read_count, sample_name)
        else:
            current_rec, current_count, current_sample = longest_seqs[clone_instance]
            if len(rec.seq) > len(current_rec.seq):
                longest_seqs[clone_instance] = (rec, read_count, sample_name)

    # Write output with header: >SampleName|CloneInstance|ReadCount
    with open(args.output, 'w') as out_handle:
        for clone_instance, (rec, read_count, sample_name) in longest_seqs.items():
            rec.id = f"{sample_name}|{clone_instance}|{read_count}"
            rec.description = ""
            SeqIO.write(rec, out_handle, "fasta")

if __name__ == "__main__":
    main()
