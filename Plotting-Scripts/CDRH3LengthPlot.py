#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import argparse

def main(input_file, output_prefix, plot_type, bar_style):
    # Load TSV file
    df = pd.read_csv(input_file, sep="\t")

    # Check required columns
    required_cols = {"ClusterID", "tissue", "CDRH3_length"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"Input file must contain columns: {required_cols}")

    # Average CDRH3 length per (ClusterID, tissue), then round to nearest integer
    df_avg = (
        df.groupby(["ClusterID", "tissue"], as_index=False)["CDRH3_length"]
        .mean()
    )
    df_avg["CDRH3_length"] = df_avg["CDRH3_length"].round().astype(int)

    # Count clusters per tissue and averaged length
    counts = (
        df_avg.groupby(["tissue", "CDRH3_length"])["ClusterID"]
        .nunique()
        .reset_index(name="count")
    )

    # Normalize counts to percentage within each tissue
    totals = counts.groupby("tissue")["count"].transform("sum")
    counts["percent"] = 100 * counts["count"] / totals

    # Pivot for plotting
    pivot_df = counts.pivot(index="CDRH3_length", columns="tissue", values="percent").fillna(0)

    # Compute mean and SD of CDRH3 lengths per tissue (no rounding here)
    stats = (
        df_avg.groupby("tissue")["CDRH3_length"]
        .agg(["mean", "std", "count"])
        .round(2)   # keep two decimals for table
    )
    mean_lengths = stats["mean"].round(1)  # one decimal for plotting

    # --- PLOTTING ---
    fig, ax = plt.subplots(figsize=(10,6))

    if plot_type == "line":
        color_map = {}
        for tissue in pivot_df.columns:
            line, = ax.plot(pivot_df.index, pivot_df[tissue], label=tissue)
            color_map[tissue] = line.get_color()
        ax.set_ylabel("% of Clusters")

        # Add mean lines + labels
        y_max = pivot_df.max().max()
        for tissue, mean_val in mean_lengths.items():
            ax.axvline(x=mean_val, linestyle="--", color=color_map[tissue], alpha=0.8)
            

    elif plot_type == "bar":
        if bar_style == "group":
            pivot_df.plot(kind="bar", width=0.8, ax=ax)
        elif bar_style == "stack":
            pivot_df.plot(kind="bar", stacked=True, width=0.8, ax=ax)
        else:
            raise ValueError("Invalid bar style. Choose 'group' or 'stack'.")
        ax.set_ylabel("% of Clusters")

        # Map categorical x positions back to true lengths
        x_positions = dict(zip(pivot_df.index, range(len(pivot_df.index))))
        handles, labels = ax.get_legend_handles_labels()
        color_map = {t: h.get_facecolor()[0:3] for t,h in zip(pivot_df.columns, handles)}

        # Add mean lines + labels
        y_max = pivot_df.max().max()
        for tissue, mean_val in mean_lengths.items():
            xpos = x_positions.get(round(mean_val), None)
            if xpos is not None:
                ax.axvline(x=xpos, linestyle="--", color=color_map[tissue], alpha=0.8)
                ax.text(xpos, y_max*1.02, str(mean_val),
                        color=color_map[tissue], rotation=90,
                        ha="center", va="bottom", fontsize=8)

    else:
        raise ValueError("Invalid plot type. Choose 'line' or 'bar'.")

    ax.set_xlabel("CDRH3 Length (aa)")
    ax.set_title(f"CDRH3 Length Distribution by Tissue ({plot_type.capitalize()} Plot)")

    # Legend only for tissues (mean lines excluded)
    ax.legend(title="Tissue", bbox_to_anchor=(1.05, 1), loc="upper left")
    ax.grid(True, axis="y", linestyle="--", alpha=0.7)
    plt.tight_layout()

    # Save
    suffix = f"_{plot_type}_{bar_style}" if plot_type=="bar" else f"_{plot_type}"
    plt.savefig(f"{output_prefix}_CDRH3_length_distribution{suffix}.png", dpi=300)
    plt.savefig(f"{output_prefix}_CDRH3_length_distribution{suffix}.pdf")
    plt.close()

    # Save frequency table
    counts.to_csv(f"{output_prefix}_CDRH3_length_frequencies.tsv", sep="\t", index=False)
    # Save mean + SD summary
    stats.to_csv(f"{output_prefix}_CDRH3_length_stats.tsv", sep="\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot CDRH3 length distributions by tissue (cluster-averaged, rounded)")
    parser.add_argument("input_tsv", help="Input TSV file with ClusterID, tissue, CDRH3_length")
    parser.add_argument("-o", "--output", default="CDRH3_output", help="Output file prefix")
    parser.add_argument("--plot", choices=["line", "bar"], default="line", help="Type of plot: 'line' or 'bar' (default: line)")
    parser.add_argument("--barstyle", choices=["group", "stack"], default="group", help="Bar style if plot is 'bar': 'group' or 'stack' (default: group)")
    args = parser.parse_args()

    main(args.input_tsv, args.output, args.plot, args.barstyle)
