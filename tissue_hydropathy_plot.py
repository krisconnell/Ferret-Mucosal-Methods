#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import numpy as np

# Hydropathy function
def avg_hydro(aa_str):
    """Takes a string of amino acid symbols. Returns the average hydrophobicity.
    Ignores unknown amino acids (e.g., X)."""
    if len(aa_str) == 0:
        return np.nan
    hydros = {
        'I':4.5, 'V':4.2, 'L':3.8, 'F':2.8, 'C':2.5,
        'M':1.9, 'A':1.8, 'G':-0.4, 'T':-0.7, 'S':-0.8,
        'W':-0.9, 'Y':-1.3, 'P':-1.6, 'H':-3.2, 'E':-3.5,
        'Q':-3.5, 'D':-3.5, 'N':-3.5, 'K':-3.9, 'R':-4.5
    }
    known_aas = [aa for aa in aa_str if aa in hydros]
    if not known_aas:
        return np.nan
    total = sum(hydros[aa] for aa in known_aas)
    return total / len(known_aas)

def compute_avg_hydro_per_position(seqs):
    """Compute average hydrophobicity per amino acid position for a list of sequences."""
    if not seqs:
        return np.array([])
    
    max_len = max(len(s) for s in seqs)
    hydro_matrix = []

    for s in seqs:
        # Compute per-residue hydropathy, ignoring unknown AAs
        per_residue = [avg_hydro(aa) for aa in s]
        # Pad shorter sequences with np.nan
        per_residue += [np.nan] * (max_len - len(per_residue))
        hydro_matrix.append(per_residue)

    hydro_array = np.array(hydro_matrix, dtype=float)
    # Compute mean across sequences at each position, ignoring NaNs
    mean_per_pos = np.nanmean(hydro_array, axis=0)
    return mean_per_pos

def main(input_file, output_prefix):
    df = pd.read_csv(input_file, sep="\t", low_memory=False)

    regions = ['cdr1_aa','cdr2_aa','cdr3_aa','fwr1_aa','fwr2_aa','fwr3_aa','fwr4_aa']

    for region in regions:
        if region not in df.columns:
            print(f"Skipping region {region}: column not found in TSV.")
            continue

        # Group by tissue and cluster
        grouped = df.groupby(['tissue', 'ClusterID'])[region].apply(list).reset_index()
        
        # Compute average hydropathy per cluster
        tissue_data = {}
        for tissue in grouped['tissue'].unique():
            tissue_df = grouped[grouped['tissue']==tissue]
            cluster_hydro = []
            for seq_list in tissue_df[region]:
                # Flatten if needed (usually each cluster has multiple rows)
                all_seqs = [s for s in seq_list if isinstance(s, str)]
                if not all_seqs:
                    continue
                mean_hydro = compute_avg_hydro_per_position(all_seqs)
                cluster_hydro.append(mean_hydro)
            if cluster_hydro:
                # Average across clusters for tissue
                tissue_mean = np.nanmean(np.array([np.pad(h, (0,max(map(len, cluster_hydro))-len(h)), constant_values=np.nan)
                                                  for h in cluster_hydro]), axis=0)
                tissue_data[tissue] = tissue_mean

        # Plotting
        plt.figure(figsize=(10,6))
        for tissue, hydro_vals in tissue_data.items():
            positions = np.arange(1, len(hydro_vals)+1)
            plt.plot(positions, hydro_vals, label=tissue)
        plt.xlabel("Amino Acid Position")
        plt.ylabel("Average Hydropathy")
        plt.title(f"Average Hydropathy for {region.upper()}")
        plt.legend(title="Tissue")
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_{region}_hydropathy.png", dpi=300)
        plt.savefig(f"{output_prefix}_{region}_hydropathy.pdf")
        plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot average hydrophobicity per cluster/tissue")
    parser.add_argument("input_tsv", help="Input TSV file with ClusterID, tissue, and AA sequences for CDRs/FWRs")
    parser.add_argument("-o", "--output", default="hydropathy_output", help="Output file prefix")
    args = parser.parse_args()

    main(args.input_tsv, args.output)
