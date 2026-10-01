#!/usr/bin/env python3

import os
import re
import argparse
import pandas as pd
import matplotlib.pyplot as plt

# -------------------------
# Tissue color palette
# -------------------------

TISSUE_COLORS = {
    'PBMC': '#7C0002',
    'LN': '#E8192B',
    'Spleen': '#F46D43',
    'BAL': '#FDAE61',
    'Tonsil': '#66C2A5',
    'NasalTurb': '#3288BD',
    'NasalCell': "#47C234",
    'Oral': "#C28169",
    'Lung': '#5E4FA2'
}

# -------------------------
# Filename parser
# -------------------------

pattern = re.compile(
    r'(?P<tissue>[A-Za-z]+)_threshold(?P<threshold>\d+)_r(?P<order>\d+)_diversity\.csv'
)

def parse_filename(fname):

    m = pattern.search(fname)

    if not m:
        return None

    return (
        m.group("tissue"),
        int(m.group("threshold")),
        int(m.group("order"))
    )

# -------------------------
# Read diversity
# -------------------------

def read_diversity(path):

    df = pd.read_csv(path, sep=None, engine="python")

    if "diversity" not in df.columns:
        raise ValueError(f"No diversity column in {path}")

    return float(df["diversity"].iloc[0])

# -------------------------
# Collect dataset
# -------------------------

def collect_data(root):

    rows = []

    for dirpath, _, files in os.walk(root):

        for f in files:

            if not f.endswith("_diversity.csv"):
                continue

            parsed = parse_filename(f)

            if parsed is None:
                continue

            tissue, threshold, order = parsed

            # Only collect q=0–4
            if order > 4:
                continue

            full = os.path.join(dirpath, f)

            diversity = read_diversity(full)

            rows.append({
                "tissue": tissue,
                "threshold": threshold,
                "q": order,
                "diversity": diversity
            })

    return pd.DataFrame(rows)

# -------------------------
# Plot Hill curves
# -------------------------

def plot_hill_curves(df, outdir):

    os.makedirs(outdir, exist_ok=True)

    thresholds = sorted(df["threshold"].unique())

    for threshold in thresholds:

        sub = df[df["threshold"] == threshold]

        plt.figure(figsize=(7,5))

        for tissue, g in sub.groupby("tissue"):

            g = g.sort_values("q")

            color = TISSUE_COLORS.get(tissue, "black")

            plt.plot(
                g["q"],
                g["diversity"],
                marker="o",
                linewidth=2,
                label=tissue,
                color=color
            )

        plt.xlabel("Hill order (q)")
        plt.ylabel("Diversity")
        plt.title(f"Hill Diversity Curve (threshold {threshold})")

        plt.xticks([0,1,2,3,4])
        plt.yscale("log")

        plt.legend(frameon=False)

        plt.tight_layout()

        outfile = f"hill_curve_threshold{threshold}.png"

        plt.savefig(os.path.join(outdir, outfile), dpi=300)

        plt.close()

# -------------------------
# Main
# -------------------------

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "directory",
        help="Animal directory containing tissue folders"
    )

    parser.add_argument(
        "-o",
        "--outdir",
        default="hill_plots"
    )

    args = parser.parse_args()

    df = collect_data(args.directory)

    if df.empty:
        print("No diversity files found.")
        return

    os.makedirs(args.outdir, exist_ok=True)

    df.to_csv(
        os.path.join(args.outdir, "compiled_hill_diversities.csv"),
        index=False
    )

    plot_hill_curves(df, args.outdir)


if __name__ == "__main__":
    main()