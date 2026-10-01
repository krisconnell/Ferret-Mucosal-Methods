#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import numpy as np

def main(input_file, output_prefix, plot_type, bar_style, bin_width, min_shm, max_shm):
    # Load TSV
    df = pd.read_csv(input_file, sep="\t")

    # Check required columns
    required_cols = {"ClusterID", "tissue", "v_identity", "nt_seq_count"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"Input file must contain columns: {required_cols}")

    # Compute %_V_SHM
    df["%_V_SHM"] = 100 - df["v_identity"]

    # Filter: ignore %_V_SHM = 0 and limit to specified range
    df = df[(df["%_V_SHM"] > 0) & (df["%_V_SHM"] >= min_shm) & (df["%_V_SHM"] <= max_shm)]

    # Weighted mean %_V_SHM per cluster within tissue
    df_cluster = (
        df.groupby(["ClusterID", "tissue"])
          .apply(lambda x: np.average(x["%_V_SHM"], weights=x["nt_seq_count"]))
          .reset_index(name="mean_%_V_SHM")
    )
    df_cluster["mean_%_V_SHM"] = df_cluster["mean_%_V_SHM"].round(1)

    # Compute tissue-level mean for vertical lines
    tissue_means = df_cluster.groupby("tissue")["mean_%_V_SHM"].mean().round(1)

    # Bin %_V_SHM within specified range
    bins = np.arange(min_shm, max_shm + bin_width, bin_width)
    df_cluster["SHM_bin"] = pd.cut(df_cluster["mean_%_V_SHM"], bins, include_lowest=True, right=False)
    df_cluster["SHM_bin_center"] = df_cluster["SHM_bin"].apply(lambda x: x.left + bin_width/2)

    # Count clusters per tissue and SHM bin
    counts = (
        df_cluster.groupby(["tissue", "SHM_bin_center"])["ClusterID"]
        .nunique()
        .reset_index(name="count")
    )

    # Normalize counts to % of clusters per tissue
    totals = counts.groupby("tissue")["count"].transform("sum")
    counts["percent"] = 100 * counts["count"] / totals

    # Pivot for plotting
    pivot_df = counts.pivot(index="SHM_bin_center", columns="tissue", values="percent").fillna(0)

    # --- PLOTTING ---
    fig, ax = plt.subplots(figsize=(10,6))

    if plot_type == "line":
        color_map = {}
        for tissue in pivot_df.columns:
            # Remove markers
            line, = ax.plot(pivot_df.index, pivot_df[tissue], label=tissue)
            color_map[tissue] = line.get_color()
        ax.set_xlabel("% V_SHM")
        ax.set_ylabel("% of Clusters")

        # Add vertical mean lines (no numeric labels)
        for tissue, mean_val in tissue_means.items():
            ax.axvline(x=mean_val, linestyle="--", color=color_map[tissue], alpha=0.8)

    elif plot_type == "bar":
        if bar_style == "group":
            pivot_df.plot(kind="bar", width=0.8, ax=ax)
        elif bar_style == "stack":
            pivot_df.plot(kind="bar", stacked=True, width=0.8, ax=ax)
        else:
            raise ValueError("Invalid bar style. Choose 'group' or 'stack'.")
        ax.set_xlabel("% V_SHM")
        ax.set_ylabel("% of Clusters")

        # Map categorical x positions
        x_positions = dict(zip(pivot_df.index, range(len(pivot_df.index))))
        handles, labels = ax.get_legend_handles_labels()
        color_map = {t: h.get_facecolor()[0:3] for t,h in zip(pivot_df.columns, handles)}

        for tissue, mean_val in tissue_means.items():
            xpos = x_positions.get(round(mean_val), None)
            if xpos is not None:
                ax.axvline(x=xpos, linestyle="--", color=color_map[tissue], alpha=0.8)

    else:
        raise ValueError("Invalid plot type. Choose 'line' or 'bar'.")

    ax.set_title(f"% V_SHM Distribution by Tissue ({plot_type.capitalize()} Plot)")
    ax.legend(title="Tissue", bbox_to_anchor=(1.05, 1), loc="upper left")
    ax.grid(True, axis="y", linestyle="--", alpha=0.7)
    plt.tight_layout()

    suffix = f"_{plot_type}_{bar_style}" if plot_type=="bar" else f"_{plot_type}"
    plt.savefig(f"{output_prefix}_V_SHM_distribution{suffix}.png", dpi=300)
    plt.savefig(f"{output_prefix}_V_SHM_distribution{suffix}.pdf")
    plt.close()

    # Save cluster-level SHM table
    df_cluster.to_csv(f"{output_prefix}_cluster_weighted_V_SHM.tsv", sep="\t", index=False)
    # Save frequency table
    counts.to_csv(f"{output_prefix}_V_SHM_frequencies.tsv", sep="\t", index=False)
    # Save tissue summary table
    stats = df_cluster.groupby("tissue")["mean_%_V_SHM"].agg(["mean", "std", "count"]).round(2)
    stats.to_csv(f"{output_prefix}_V_SHM_summary_stats.tsv", sep="\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot %V_SHM distributions by tissue (weighted per cluster)")
    parser.add_argument("input_tsv", help="Input TSV file with ClusterID, tissue, v_identity, nt_seq_count")
    parser.add_argument("-o", "--output", default="V_SHM_output", help="Output file prefix")
    parser.add_argument("--plot", choices=["line", "bar"], default="line", help="Plot type: line or bar")
    parser.add_argument("--barstyle", choices=["group", "stack"], default="group", help="Bar style if plot is 'bar'")
    parser.add_argument("--bin", type=float, default=1.0, help="Bin width for %_V_SHM (default 1%)")
    parser.add_argument("--min_shm", type=float, default=0.1, help="Minimum %_V_SHM to include (default 0.1)")
    parser.add_argument("--max_shm", type=float, default=30.0, help="Maximum %_V_SHM to include (default 30)")
    args = parser.parse_args()

    main(args.input_tsv, args.output, args.plot, args.barstyle, args.bin, args.min_shm, args.max_shm)
