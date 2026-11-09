#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# ============================================================
# Plot settings
# ============================================================

sns.set(style="whitegrid", context="talk")


TISSUE_COLORS = {
    "PBMC": "#7C0002",
    "LN": "#E8192B",
    "Spleen": "#F46D43",
    "BAL": "#FDAE61",
    "Tonsil": "#66C2A5",
    "NasalTurb": "#3288BD",
    "NasalCell": "#47C234",
    "Oral": "#C28169",
    "Lung": "#5E4FA2",
}


# ============================================================
# Parse threshold number
# ============================================================

def get_threshold(filename):
    """
    Extract threshold number from a filename.

    Examples
    --------
    BAL1_thresh3_diversity.csv -> 3
    LN1_thresh10_rarefaction_curve.csv -> 10
    """

    match = re.search(r"threshold(\d+)", filename)

    if match:
        return int(match.group(1))

    return None


# ============================================================
# Read a CSV or TSV file
# ============================================================

def read_table(file):
    """
    Read either comma-separated or tab-separated files.

    sep=None with engine='python' allows pandas to infer
    the delimiter.
    """

    return pd.read_csv(file, sep=None, engine="python")


# ============================================================
# Extract Hill diversity
# ============================================================

def extract_hill(div_file, sample, threshold):
    """
    Extract Hill diversity of order q=1 from a diversity file.

    Expected columns:
        order
        diversity

    Returns
    -------
    dict
        sample
        threshold
        Hill_q1
    """

    df = read_table(div_file)

    required = {"order", "diversity"}
    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            f"{div_file} is missing required columns: {sorted(missing)}"
        )

    # Find rows corresponding to Hill order q=1
    q1 = df.loc[df["order"] == 1, "diversity"]

    if len(q1) == 0:
        hill_q1 = np.nan

    else:
        # Normally there should only be one q=1 row.
        # Take the first one if there is more than one.
        hill_q1 = q1.iloc[0]

    return {
        "sample": sample,
        "threshold": threshold,
        "Hill_q1": hill_q1,
    }


# ============================================================
# Find corresponding diversity / rarefaction files
# ============================================================

def find_related_file(rare_file, replacement):
    """
    Given a rarefaction metrics filename, construct the
    corresponding file by replacing:

        rarefaction_metrics

    with the requested replacement.

    Example
    -------
    BAL1_thresh3_rarefaction_metrics.csv

    becomes

    BAL1_thresh3_diversity.csv
    """

    new_name = rare_file.name.replace(
        "rarefaction_metrics",
        replacement
    )

    candidate = rare_file.with_name(new_name)

    if candidate.exists():
        return candidate

    # If extension differs, try csv/tsv alternatives.
    stem = candidate.stem

    alternatives = [
        rare_file.with_name(stem + ".csv"),
        rare_file.with_name(stem + ".tsv"),
    ]

    for alternative in alternatives:
        if alternative.exists():
            return alternative

    return None


# ============================================================
# Load and normalize rarefaction curve
# ============================================================

def load_rarefaction_curve(curve_file, sample, threshold):
    """
    Load a rarefaction curve and normalize sampling depth
    and expected richness.

    Expected columns:
        subsample
        expected_richness

    Returns
    -------
    DataFrame
    """

    df = read_table(curve_file)

    required = {"subsample", "expected_richness"}
    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            f"{curve_file} is missing required columns: "
            f"{sorted(missing)}"
        )

    # Remove rows with missing values
    df = df.dropna(
        subset=["subsample", "expected_richness"]
    ).copy()

    if len(df) == 0:
        raise ValueError(
            f"{curve_file} contains no usable rarefaction data."
        )

    max_subsample = df["subsample"].max()
    max_richness = df["expected_richness"].max()

    if max_subsample == 0:
        raise ValueError(
            f"{curve_file}: maximum subsample is zero."
        )

    if max_richness == 0:
        raise ValueError(
            f"{curve_file}: maximum expected richness is zero."
        )

    df["subsample_norm"] = (
        df["subsample"] / max_subsample
    )

    df["richness_norm"] = (
        df["expected_richness"] / max_richness
    )

    df["sample"] = sample
    df["threshold"] = threshold

    return df


# ============================================================
# Process one sample/tissue
# ============================================================

def process_sample(sample_dir, outdir, metrics):
    """
    Process all thresholds found in one sample directory.

    For each threshold:
        - find rarefaction metrics file
        - find rarefaction curve
        - find diversity file
        - extract Hill q=1
        - normalize rarefaction curve

    Returns
    -------
    DataFrame or None
    """

    sample = sample_dir.name

    rare_files = sorted(
        sample_dir.glob("*rarefaction*")
    )

    if len(rare_files) == 0:
        print(f"⚠ No rarefaction metrics files found for {sample}")
        return None

    print(f"\nProcessing {sample}")

    curves = []

    for rare_file in rare_files:

        threshold = get_threshold(rare_file.name)

        if threshold is None:
            print(
                f"⚠ Could not determine threshold from "
                f"{rare_file.name}"
            )
            continue

        # ----------------------------------------------------
        # Find diversity file
        # ----------------------------------------------------

        div_file = find_related_file(
            rare_file,
            "diversity"
        )

        if div_file is None:
            print(
                f"⚠ Missing diversity file corresponding to "
                f"{rare_file.name}"
            )
            continue

        # ----------------------------------------------------
        # Find rarefaction curve
        # ----------------------------------------------------

        curve_file = find_related_file(
            rare_file,
            "curve"
        )

        if curve_file is None:
            print(
                f"⚠ Missing rarefaction curve corresponding to "
                f"{rare_file.name}"
            )
            continue

        # ----------------------------------------------------
        # Extract Hill q=1
        # ----------------------------------------------------

        try:

            hill_result = extract_hill(
                div_file,
                sample,
                threshold
            )

            metrics.append(hill_result)

            print(
                f"  threshold {threshold}: "
                f"Hill q=1 = {hill_result['Hill_q1']}"
            )

        except Exception as e:

            print(
                f"⚠ Could not extract Hill diversity from "
                f"{div_file}: {e}"
            )

            continue

        # ----------------------------------------------------
        # Load rarefaction curve
        # ----------------------------------------------------

        try:

            curve_df = load_rarefaction_curve(
                curve_file,
                sample,
                threshold
            )

            curves.append(curve_df)

        except Exception as e:

            print(
                f"⚠ Could not load curve {curve_file}: {e}"
            )

            continue

    # --------------------------------------------------------
    # Nothing usable
    # --------------------------------------------------------

    if len(curves) == 0:

        print(
            f"⚠ No valid rarefaction curves for {sample}"
        )

        return None

    # --------------------------------------------------------
    # Combine thresholds
    # --------------------------------------------------------

    curves = pd.concat(
        curves,
        ignore_index=True
    )

    # --------------------------------------------------------
    # Save normalized curve data
    # --------------------------------------------------------

    output_file = (
        outdir /
        f"{sample}_normalized_rarefaction.csv"
    )

    curves.to_csv(
        output_file,
        index=False
    )

    print(
        f"  Saved normalized curves: {output_file}"
    )

    # --------------------------------------------------------
    # Plot per-sample rarefaction curves
    # --------------------------------------------------------

    plot_sample_rarefaction(
        curves,
        sample,
        outdir
    )

    return curves


# ============================================================
# Plot rarefaction curves for one sample
# ============================================================

def plot_sample_rarefaction(curves, sample, outdir):
    """
    Plot normalized rarefaction curves for all thresholds
    within one sample.
    """

    plt.figure(figsize=(8, 6))

    sns.lineplot(
        data=curves,
        x="subsample_norm",
        y="richness_norm",
        hue="threshold",
        marker="o",
        palette="viridis"
    )

    plt.xlabel("Normalized sampling depth")
    plt.ylabel("Normalized expected richness")

    plt.title(
        f"{sample} Rarefaction Curves"
    )

    plt.ylim(0, 1.05)

    plt.tight_layout()

    output_file = (
        outdir /
        f"{sample}_rarefaction_curves.png"
    )

    plt.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print(
        f"  Saved rarefaction plot: {output_file}"
    )


# ============================================================
# Plot Hill q=1 diversity
# ============================================================

def plot_hill(df, outdir):
    """
    Plot Hill diversity q=1.

    Each sample is shown on the x-axis and threshold is
    represented by color.
    """

    if df.empty:
        print("⚠ No Hill diversity data to plot.")
        return

    plt.figure(figsize=(10, 6))

    sns.barplot(
        data=df,
        x="sample",
        y="Hill_q1",
        hue="threshold",
        palette="viridis",
        errorbar=None
    )

    plt.xlabel("Sample")
    plt.ylabel("Hill diversity q=1\n(effective number of clones)")

    plt.title(
        "Hill Diversity (q=1)"
    )

    plt.xticks(rotation=45)

    plt.tight_layout()

    output_file = (
        outdir /
        "hill_q1_diversity.png"
    )

    plt.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print(
        f"Saved Hill diversity plot: {output_file}"
    )


# ============================================================
# Plot Hill q=1 versus threshold
# ============================================================

def plot_hill_vs_threshold(df, outdir):
    """
    Plot Hill q=1 as a function of threshold.

    This is useful for seeing how the clone threshold
    affects estimated diversity.
    """

    if df.empty:
        return

    plt.figure(figsize=(9, 6))

    for sample, group in df.groupby("sample"):

        group = group.sort_values("threshold")

        # Try to identify the base tissue name.
        # For example BAL1 -> BAL.
        tissue = re.sub(r"\d+$", "", sample)

        color = TISSUE_COLORS.get(
            tissue,
            None
        )

        if color is not None:

            plt.plot(
                group["threshold"],
                group["Hill_q1"],
                marker="o",
                linewidth=2,
                label=sample,
                color=color
            )

        else:

            plt.plot(
                group["threshold"],
                group["Hill_q1"],
                marker="o",
                linewidth=2,
                label=sample
            )

    plt.xlabel("Clone threshold")
    plt.ylabel(
        "Hill diversity q=1\n(effective number of clones)"
    )

    plt.title(
        "Hill q=1 Diversity Across Clone Thresholds"
    )

    plt.legend(
        frameon=False,
        bbox_to_anchor=(1.02, 1),
        loc="upper left"
    )

    plt.tight_layout()

    output_file = (
        outdir /
        "hill_q1_vs_threshold.png"
    )

    plt.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print(
        f"Saved Hill-vs-threshold plot: {output_file}"
    )


# ============================================================
# Plot cross-sample normalized rarefaction
# ============================================================

def plot_cross_sample_rarefaction(
    all_curves,
    threshold,
    outdir
):
    """
    Compare samples at one specific clone threshold.

    Example:
        --curve-threshold 3

    will compare all samples using their threshold-3
    rarefaction curves.
    """

    if all_curves.empty:
        print(
            "⚠ No curve data available for "
            "cross-sample rarefaction."
        )
        return

    curves = all_curves[
        all_curves["threshold"] == threshold
    ].copy()

    if curves.empty:

        print(
            f"⚠ No rarefaction curves found for "
            f"threshold {threshold}"
        )

        return

    plt.figure(figsize=(9, 6))

    for sample, group in curves.groupby("sample"):

        group = group.sort_values("subsample_norm")

        # Remove numerical suffix from sample name
        # BAL1 -> BAL
        # LN2 -> LN
        # Lung1 -> Lung
        tissue = re.sub(r"\d+$", "", sample)

        color = TISSUE_COLORS.get(
            tissue,
            None
        )

        if color is not None:

            plt.plot(
                group["subsample_norm"],
                group["richness_norm"],
                linewidth=2,
                label=sample,
                color=color
            )

        else:

            plt.plot(
                group["subsample_norm"],
                group["richness_norm"],
                linewidth=2,
                label=sample
            )

    plt.xlabel("Normalized sampling depth")
    plt.ylabel("Normalized expected richness")

    plt.title(
        f"Normalized Rarefaction Curves "
        f"(Threshold {threshold})"
    )

    plt.ylim(0, 1.05)

    plt.legend(
        frameon=False,
        bbox_to_anchor=(1.02, 1),
        loc="upper left"
    )

    plt.tight_layout()

    output_file = (
        outdir /
        f"cross_sample_rarefaction_threshold{threshold}.png"
    )

    plt.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print(
        f"Saved cross-sample rarefaction plot: "
        f"{output_file}"
    )


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Summarize Hill q=1 diversity and rarefaction "
            "curves across tissue/sample directories."
        )
    )

    parser.add_argument(
        "parent_directory",
        help=(
            "Parent directory containing one directory "
            "per sample/tissue."
        )
    )

    parser.add_argument(
        "-o",
        "--output",
        default="diversity_summary",
        help=(
            "Output directory "
            "(default: diversity_summary)"
        )
    )

    parser.add_argument(
        "--curve-threshold",
        type=int,
        default=None,
        help=(
            "Threshold to use for the cross-sample "
            "rarefaction plot. For example: --curve-threshold 3"
        )
    )

    args = parser.parse_args()

    # --------------------------------------------------------
    # Directories
    # --------------------------------------------------------

    parent = Path(args.parent_directory)
    outdir = Path(args.output)

    if not parent.exists():
        raise FileNotFoundError(
            f"Parent directory does not exist: {parent}"
        )

    outdir.mkdir(
        parents=True,
        exist_ok=True
    )

    # --------------------------------------------------------
    # Storage
    # --------------------------------------------------------

    metrics = []
    all_curves = []

    # --------------------------------------------------------
    # Process each sample directory
    # --------------------------------------------------------

    sample_dirs = sorted(
        [
            d for d in parent.iterdir()
            if d.is_dir()
        ]
    )

    if len(sample_dirs) == 0:

        print(
            f"❌ No sample directories found in {parent}"
        )

        return

    for sample_dir in sample_dirs:

        curves = process_sample(
            sample_dir,
            outdir,
            metrics
        )

        if curves is not None:
            all_curves.append(curves)

    # --------------------------------------------------------
    # Check whether anything was processed
    # --------------------------------------------------------

    if len(metrics) == 0:

        print(
            "\n❌ No valid diversity data found."
        )

        return

    # --------------------------------------------------------
    # Hill diversity summary
    # --------------------------------------------------------

    df_metrics = pd.DataFrame(metrics)

    df_metrics = df_metrics.sort_values(
        ["sample", "threshold"]
    )

    hill_output = (
        outdir /
        "all_hill_diversity.csv"
    )

    df_metrics.to_csv(
        hill_output,
        index=False
    )

    print(
        f"\nSaved Hill diversity summary: {hill_output}"
    )

    # --------------------------------------------------------
    # Hill plots
    # --------------------------------------------------------

    plot_hill(
        df_metrics,
        outdir
    )

    plot_hill_vs_threshold(
        df_metrics,
        outdir
    )

    # --------------------------------------------------------
    # Combine all rarefaction curves
    # --------------------------------------------------------

    if len(all_curves) > 0:

        combined_curves = pd.concat(
            all_curves,
            ignore_index=True
        )

        combined_output = (
            outdir /
            "all_normalized_rarefaction_curves.csv"
        )

        combined_curves.to_csv(
            combined_output,
            index=False
        )

        print(
            f"Saved combined rarefaction data: "
            f"{combined_output}"
        )

        # ----------------------------------------------------
        # Cross-sample plot
        # ----------------------------------------------------

        if args.curve_threshold is not None:

            plot_cross_sample_rarefaction(
                combined_curves,
                args.curve_threshold,
                outdir
            )

    print(
        "\n✅ Diversity analysis complete."
    )


if __name__ == "__main__":
    main()
