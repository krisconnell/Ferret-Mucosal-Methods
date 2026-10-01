#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# ============================================================
# Plotting
# ============================================================

sns.set(
    style="whitegrid",
    context="talk"
)


# ============================================================
# Tissue colors
# ============================================================

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

    # Older tissue names, if present
    "Intestine": "#66C2A5",
    "Nasal": "#3288BD",
}


# ============================================================
# Tissue name cleanup
# ============================================================

def clean_tissue_name(tissue):
    """
    Clean tissue names for plotting.

    For example:

        BAL1       -> BAL
        LN1        -> LN
        Lung1      -> Lung
        PBMC1      -> PBMC

    This is useful when the directory itself is named BAL1,
    even though the "1" represents Hill order rather than
    tissue identity.

    Known current tissue names are handled explicitly.
    """

    tissue = str(tissue).strip()

    known_tissues = sorted(
        TISSUE_COLORS.keys(),
        key=len,
        reverse=True
    )

    for name in known_tissues:

        if tissue == name:
            return name

        if tissue == f"{name}1":
            return name

    # Generic fallback:
    # Remove a trailing "1" if present.
    if re.match(r".+1$", tissue):
        return tissue[:-1]

    return tissue


# ============================================================
# Filename parser
# ============================================================

def parse_filename(filename):
    """
    Supports the two filename conventions used in the analysis.

    OLD FORMAT
    ----------
    BAL1_thresh1_diversity.csv
    BAL1_thresh1_rarefaction_curve.csv
    BAL1_thresh1_rarefaction_metrics.csv

    NEW FORMAT
    ----------
    BAL_threshold1_rarefaction.csv

    Old format:
        BAL = tissue
        1   = Hill order
        thresh1 = clone threshold

    New format is assumed to be Hill order 1.
    """

    name = Path(filename).stem

    # --------------------------------------------------------
    # Old format
    # --------------------------------------------------------

    old_match = re.match(
        r"^(?P<tissue>.+?)(?P<order>\d+)_"
        r"thresh(?P<threshold>\d+)_"
        r"(?P<type>.+)$",
        name
    )

    if old_match:

        tissue = clean_tissue_name(
            old_match.group("tissue")
        )

        return {
            "tissue": tissue,
            "order": int(
                old_match.group("order")
            ),
            "threshold": int(
                old_match.group("threshold")
            ),
            "type": old_match.group("type"),
            "format": "old"
        }

    # --------------------------------------------------------
    # New format
    # --------------------------------------------------------

    new_match = re.match(
        r"^(?P<tissue>.+?)_"
        r"threshold(?P<threshold>\d+)_"
        r"(?P<type>.+)$",
        name
    )

    if new_match:

        tissue = clean_tissue_name(
            new_match.group("tissue")
        )

        return {
            "tissue": tissue,
            "order": 1,
            "threshold": int(
                new_match.group("threshold")
            ),
            "type": new_match.group("type"),
            "format": "new"
        }

    return None


# ============================================================
# Find rarefaction curve files
# ============================================================

def find_curve_files(tissue_dir):
    """
    Find only Hill-1 rarefaction curve files.

    OLD:
        BAL1_thresh1_rarefaction_curve.csv

    NEW:
        BAL_threshold1_rarefaction.csv

    Ignore:

        *_diversity.csv
        *_rarefaction_metrics.csv
        Hill orders other than 1
    """

    curve_files = []

    for file in tissue_dir.glob("*.csv"):

        info = parse_filename(
            file.name
        )

        if info is None:
            continue

        # Only Hill order 1
        if info["order"] != 1:
            continue

        # Old format
        if info["format"] == "old":

            if info["type"] == "rarefaction_curve":
                curve_files.append(file)

        # New format
        elif info["format"] == "new":

            if info["type"] == "curve":
                curve_files.append(file)

    return sorted(curve_files)


# ============================================================
# Final diversity from a rarefaction curve
# ============================================================

def get_final_diversity(group):
    """
    Return the final Hill-1 effective diversity.

    The final point is defined as the point with the largest
    subsample size.
    """

    group = group.sort_values(
        "subsample"
    )

    return group[
        "expected_richness"
    ].iloc[-1]


# ============================================================
# Calculate absolute diversity summary
# ============================================================

def calculate_threshold_summary(all_curves):
    """
    Calculate mean and median final Hill-1 effective diversity
    across tissues for each clone threshold.

    Returns:

        threshold
        mean_diversity
        median_diversity
        sd_diversity
        n_tissues
    """

    rows = []

    for threshold, threshold_group in (
        all_curves.groupby("threshold")
    ):

        values = []

        for tissue, tissue_group in (
            threshold_group.groupby("tissue")
        ):

            value = get_final_diversity(
                tissue_group
            )

            if pd.notna(value):
                values.append(value)

        if len(values) == 0:
            continue

        rows.append(
            {
                "threshold": int(threshold),
                "mean_diversity": np.mean(values),
                "median_diversity": np.median(values),
                "sd_diversity": (
                    np.std(values, ddof=1)
                    if len(values) > 1
                    else 0
                ),
                "n_tissues": len(values)
            }
        )

    if len(rows) == 0:

        return pd.DataFrame(
            columns=[
                "threshold",
                "mean_diversity",
                "median_diversity",
                "sd_diversity",
                "n_tissues"
            ]
        )

    return (
        pd.DataFrame(rows)
        .sort_values("threshold")
        .reset_index(drop=True)
    )


# ============================================================
# Calculate elbow distances
# ============================================================

def calculate_elbow_distances(all_curves):
    """
    Calculate the geometric distance of every threshold from
    the straight line connecting the first and last points.

    X and Y are scaled internally for the geometric calculation.

    The actual diversity values are NOT normalized in the
    resulting plot.
    """

    summary = calculate_threshold_summary(
        all_curves
    )

    if len(summary) < 3:

        summary["elbow_distance"] = np.nan

        return summary

    x = summary[
        "threshold"
    ].values.astype(float)

    y = summary[
        "mean_diversity"
    ].values.astype(float)

    x_range = x.max() - x.min()
    y_range = y.max() - y.min()

    if x_range == 0 or y_range == 0:

        summary["elbow_distance"] = 0.0

        return summary

    # Scale ONLY for elbow calculation
    x_scaled = (
        (x - x.min())
        / x_range
    )

    y_scaled = (
        (y - y.min())
        / y_range
    )

    start = np.array(
        [
            x_scaled[0],
            y_scaled[0]
        ]
    )

    end = np.array(
        [
            x_scaled[-1],
            y_scaled[-1]
        ]
    )

    line_vector = end - start

    line_length = np.linalg.norm(
        line_vector
    )

    if line_length == 0:

        summary["elbow_distance"] = 0.0

        return summary

    line_unit = (
        line_vector
        / line_length
    )

    distances = []

    for i in range(len(x_scaled)):

        point = np.array(
            [
                x_scaled[i],
                y_scaled[i]
            ]
        )

        vector = point - start

        projection = (
            np.dot(
                vector,
                line_unit
            )
            * line_unit
        )

        perpendicular = (
            vector
            - projection
        )

        distance = np.linalg.norm(
            perpendicular
        )

        distances.append(distance)

    summary[
        "elbow_distance"
    ] = distances

    return summary


# ============================================================
# Select global threshold
# ============================================================

def choose_global_threshold(all_curves):
    """
    Select one global clone threshold using the geometric elbow
    of mean absolute Hill-1 effective diversity.
    """

    summary = calculate_elbow_distances(
        all_curves
    )

    if len(summary) == 0:

        raise ValueError(
            "No valid threshold data available."
        )

    if len(summary) < 3:

        print(
            "\nWARNING:"
        )

        print(
            "Fewer than 3 thresholds are available."
        )

        print(
            "A geometric elbow cannot be reliably calculated."
        )

        selected = int(
            summary.iloc[0]["threshold"]
        )

        print(
            f"Using lowest available threshold: "
            f"{selected}"
        )

        return selected, summary

    distances = summary[
        "elbow_distance"
    ].values.copy()

    # Endpoints cannot be the elbow
    distances[0] = -np.inf
    distances[-1] = -np.inf

    elbow_index = np.argmax(
        distances
    )

    selected = int(
        summary.iloc[
            elbow_index
        ]["threshold"]
    )

    return selected, summary


# ============================================================
# Calculate tissue-specific diversity
# ============================================================

def calculate_tissue_threshold_values(
    all_curves
):
    """
    Calculate final Hill-1 effective diversity for every:

        tissue x threshold

    combination.
    """

    rows = []

    for (
        tissue,
        threshold
    ), group in all_curves.groupby(
        [
            "tissue",
            "threshold"
        ]
    ):

        final_diversity = (
            get_final_diversity(
                group
            )
        )

        rows.append(
            {
                "tissue": tissue,
                "threshold": int(threshold),
                "hill1_diversity": (
                    final_diversity
                )
            }
        )

    if len(rows) == 0:

        return pd.DataFrame(
            columns=[
                "tissue",
                "threshold",
                "hill1_diversity"
            ]
        )

    return (
        pd.DataFrame(rows)
        .sort_values(
            [
                "tissue",
                "threshold"
            ]
        )
        .reset_index(drop=True)
    )


# ============================================================
# Calculate diversity retention
# ============================================================

def calculate_diversity_retention(
    all_curves
):
    """
    Calculate Hill-1 diversity retention relative to threshold 1.

    For each tissue:

        retention =
            diversity_at_threshold
            /
            diversity_at_threshold_1
            * 100

    This is calculated independently within each tissue.
    """

    tissue_threshold = (
        calculate_tissue_threshold_values(
            all_curves
        )
    )

    if len(tissue_threshold) == 0:

        return tissue_threshold

    # --------------------------------------------------------
    # Get threshold-1 baseline
    # --------------------------------------------------------

    baseline = (
        tissue_threshold[
            tissue_threshold["threshold"] == 1
        ][
            [
                "tissue",
                "hill1_diversity"
            ]
        ]
        .rename(
            columns={
                "hill1_diversity":
                    "baseline_diversity"
            }
        )
    )

    # --------------------------------------------------------
    # Merge baseline
    # --------------------------------------------------------

    tissue_threshold = (
        tissue_threshold.merge(
            baseline,
            on="tissue",
            how="left"
        )
    )

    # --------------------------------------------------------
    # Calculate retention
    # --------------------------------------------------------

    tissue_threshold[
        "percent_retained"
    ] = (
        tissue_threshold[
            "hill1_diversity"
        ]
        /
        tissue_threshold[
            "baseline_diversity"
        ]
        * 100
    )

    return tissue_threshold


# ============================================================
# Plot main two-panel threshold figure
# ============================================================

def plot_threshold_analysis(
    all_curves,
    best_threshold,
    outdir
):
    """
    Create the main two-panel figure.

    PANEL A:
        Absolute Hill-1 effective diversity.

        Mean + median across tissues.

    PANEL B:
        Percent Hill-1 diversity retained relative to threshold 1.

        Individual tissue curves.
        Mean ± SD across tissues.

    The same global threshold is highlighted in both panels.
    """

    # --------------------------------------------------------
    # Absolute diversity
    # --------------------------------------------------------

    summary = calculate_elbow_distances(
        all_curves
    )

    summary.to_csv(
        outdir
        / "threshold_elbow_summary.csv",
        index=False
    )

    # --------------------------------------------------------
    # Retention
    # --------------------------------------------------------

    retention = calculate_diversity_retention(
        all_curves
    )

    retention.to_csv(
        outdir
        / "threshold_diversity_retention.csv",
        index=False
    )

    # --------------------------------------------------------
    # Calculate mean and SD retention
    # --------------------------------------------------------

    retention_summary = (
        retention
        .groupby("threshold")[
            "percent_retained"
        ]
        .agg(
            [
                "mean",
                "median",
                "std",
                "count"
            ]
        )
        .reset_index()
    )

    retention_summary = (
        retention_summary.rename(
            columns={
                "mean":
                    "mean_percent_retained",
                "median":
                    "median_percent_retained",
                "std":
                    "sd_percent_retained",
                "count":
                    "n_tissues"
            }
        )
    )

    retention_summary[
        "sd_percent_retained"
    ] = retention_summary[
        "sd_percent_retained"
    ].fillna(0)

    retention_summary.to_csv(
        outdir
        / "threshold_retention_summary.csv",
        index=False
    )

    # --------------------------------------------------------
    # Create figure
    # --------------------------------------------------------

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(17, 7)
    )

    # ========================================================
    # PANEL A
    # ========================================================

    ax = axes[0]

    # Mean
    ax.plot(
        summary["threshold"],
        summary["mean_diversity"],
        marker="o",
        linewidth=3,
        markersize=8,
        label=(
            "Mean Hill-1 effective diversity"
        )
    )

    # Median
    ax.plot(
        summary["threshold"],
        summary["median_diversity"],
        marker="s",
        linestyle="--",
        linewidth=2,
        markersize=6,
        alpha=0.75,
        label=(
            "Median Hill-1 effective diversity"
        )
    )

    # Global threshold
    ax.axvline(
        best_threshold,
        color="red",
        linestyle="--",
        linewidth=2.5,
        label=(
            f"Global threshold = "
            f"{best_threshold}"
        )
    )

    ax.set_xlabel(
        "Clone threshold"
    )

    ax.set_ylabel(
        "Hill-1 effective diversity"
    )

    ax.set_title(
        "A. Absolute diversity"
    )

    ax.legend(
        frameon=False,
        fontsize="small"
    )

    # ========================================================
    # PANEL B
    # ========================================================

    ax = axes[1]

    tissues = sorted(
        retention["tissue"]
        .dropna()
        .unique()
    )

    # --------------------------------------------------------
    # Individual tissues
    # --------------------------------------------------------

    for tissue in tissues:

        tissue_data = (
            retention[
                retention["tissue"] == tissue
            ]
            .sort_values("threshold")
        )

        color = TISSUE_COLORS.get(
            tissue,
            "gray"
        )

        ax.plot(
            tissue_data["threshold"],
            tissue_data["percent_retained"],
            marker="o",
            linewidth=2,
            markersize=5,
            color=color,
            label=tissue
        )

    # --------------------------------------------------------
    # Mean retention
    # --------------------------------------------------------

    ax.plot(
        retention_summary["threshold"],
        retention_summary[
            "mean_percent_retained"
        ],
        color="black",
        linewidth=3,
        linestyle="--",
        marker="s",
        markersize=5,
        label="Mean ± SD"
    )

    # --------------------------------------------------------
    # SD band
    # --------------------------------------------------------

    mean_retention = (
        retention_summary[
            "mean_percent_retained"
        ].values
    )

    sd_retention = (
        retention_summary[
            "sd_percent_retained"
        ].values
    )

    thresholds = (
        retention_summary[
            "threshold"
        ].values
    )

    lower = (
        mean_retention
        - sd_retention
    )

    upper = (
        mean_retention
        + sd_retention
    )

    ax.fill_between(
        thresholds,
        lower,
        upper,
        color="gray",
        alpha=0.15
    )

    # --------------------------------------------------------
    # Global threshold
    # --------------------------------------------------------

    ax.axvline(
        best_threshold,
        color="red",
        linestyle="--",
        linewidth=2.5,
        label=(
            f"Global threshold = "
            f"{best_threshold}"
        )
    )

    # --------------------------------------------------------
    # Threshold 1 = 100%
    # --------------------------------------------------------

    ax.axhline(
        100,
        color="black",
        linestyle=":",
        linewidth=1.5,
        alpha=0.6
    )

    ax.set_xlabel(
        "Clone threshold"
    )

    ax.set_ylabel(
        "Hill-1 diversity retained relative to threshold 1 (%)"
    )

    ax.set_title(
        "B. Within-tissue diversity retention"
    )

    ax.set_ylim(
        0,
        105
    )

    ax.legend(
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
        fontsize="small"
    )

    # --------------------------------------------------------
    # Figure title
    # --------------------------------------------------------

    fig.suptitle(
        "Clone Threshold Analysis",
        fontsize=22,
        y=1.02
    )

    plt.tight_layout()

    # --------------------------------------------------------
    # Save figure
    # --------------------------------------------------------

    output_file = (
        outdir
        / "clone_threshold_analysis.png"
    )

    plt.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print(
        "\nSaved main threshold figure:"
    )

    print(
        f"  {output_file}"
    )


# ============================================================
# Print threshold results
# ============================================================

def print_threshold_results(
    all_curves,
    best_threshold
):
    """
    Print a concise summary of diversity retention at the
    selected global threshold.
    """

    retention = calculate_diversity_retention(
        all_curves
    )

    selected = retention[
        retention["threshold"]
        == best_threshold
    ].copy()

    if len(selected) == 0:

        return

    print(
        "\n---------------------------------------"
    )

    print(
        "Diversity retained at selected "
        f"global threshold = {best_threshold}"
    )

    print(
        "---------------------------------------"
    )

    for _, row in (
        selected.sort_values("tissue")
        .iterrows()
    ):

        tissue = row["tissue"]

        percent = row[
            "percent_retained"
        ]

        diversity = row[
            "hill1_diversity"
        ]

        baseline = row[
            "baseline_diversity"
        ]

        if pd.isna(percent):

            print(
                f"{tissue:12s} "
                f"baseline={baseline:.1f} "
                f"threshold={diversity:.1f} "
                f"retention=NA"
            )

        else:

            print(
                f"{tissue:12s} "
                f"baseline={baseline:.1f} "
                f"threshold={diversity:.1f} "
                f"retention={percent:.1f}%"
            )

    valid = selected[
        selected["percent_retained"]
        .notna()
    ]

    if len(valid) > 0:

        mean_retention = (
            valid["percent_retained"]
            .mean()
        )

        median_retention = (
            valid["percent_retained"]
            .median()
        )

        print(
            "\nMean retention: "
            f"{mean_retention:.1f}%"
        )

        print(
            "Median retention: "
            f"{median_retention:.1f}%"
        )


# ============================================================
# Process one tissue directory
# ============================================================

def process_tissue(
    tissue_dir,
    outdir
):
    """
    Read all Hill-1 rarefaction curves for one tissue.
    """

    tissue = clean_tissue_name(
        tissue_dir.name
    )

    curve_files = find_curve_files(
        tissue_dir
    )

    if len(curve_files) == 0:

        print(
            f"\nWARNING: No Hill-1 rarefaction "
            f"curve files found for "
            f"{tissue_dir.name}"
        )

        return None

    print(
        f"\nProcessing {tissue}"
    )

    curves = []

    # --------------------------------------------------------
    # Read each curve
    # --------------------------------------------------------

    for curve_file in curve_files:

        info = parse_filename(
            curve_file.name
        )

        if info is None:
            continue

        threshold = info[
            "threshold"
        ]

        order = info[
            "order"
        ]

        print(
            f"   {curve_file.name}"
        )

        print(
            f"      format: "
            f"{info['format']}"
        )

        print(
            f"      Hill order: "
            f"{order}"
        )

        print(
            f"      threshold: "
            f"{threshold}"
        )

        # ----------------------------------------------------
        # Read CSV
        # ----------------------------------------------------

        try:

            df = pd.read_csv(
                curve_file
            )

        except Exception as e:

            print(
                f"   WARNING: Could not read "
                f"{curve_file.name}: {e}"
            )

            continue

        # ----------------------------------------------------
        # Required columns
        # ----------------------------------------------------

        required = [
            "subsample",
            "expected_richness"
        ]

        missing = [
            col
            for col in required
            if col not in df.columns
        ]

        if missing:

            print(
                f"   WARNING: Missing required "
                f"columns: {missing}"
            )

            print(
                f"      File: {curve_file.name}"
            )

            continue

        # ----------------------------------------------------
        # Numeric conversion
        # ----------------------------------------------------

        df[
            "subsample"
        ] = pd.to_numeric(
            df["subsample"],
            errors="coerce"
        )

        df[
            "expected_richness"
        ] = pd.to_numeric(
            df["expected_richness"],
            errors="coerce"
        )

        df = df.dropna(
            subset=[
                "subsample",
                "expected_richness"
            ]
        )

        if len(df) == 0:

            print(
                "   WARNING: No valid curve points."
            )

            continue

        # ----------------------------------------------------
        # Add metadata
        # ----------------------------------------------------

        df["threshold"] = (
            threshold
        )

        df["tissue"] = (
            tissue
        )

        df["hill_order"] = 1

        df["source_file"] = (
            curve_file.name
        )

        curves.append(
            df
        )

    # --------------------------------------------------------
    # Check data
    # --------------------------------------------------------

    if len(curves) == 0:

        print(
            f"\nWARNING: No valid Hill-1 "
            f"curves for {tissue}"
        )

        return None

    # --------------------------------------------------------
    # Combine
    # --------------------------------------------------------

    curves = pd.concat(
        curves,
        ignore_index=True
    )

    # --------------------------------------------------------
    # Save tissue curve data
    # --------------------------------------------------------

    curves.to_csv(
        outdir
        / f"{tissue}_hill1_rarefaction.csv",
        index=False
    )

    return curves


# ============================================================
# Tissue rarefaction plot
# ============================================================

def plot_rarefaction(
    curves,
    tissue,
    outdir,
    best_threshold
):
    """
    Plot Hill-1 rarefaction curves for one tissue.

    The selected global threshold is emphasized.
    """

    plt.figure(
        figsize=(8, 6)
    )

    thresholds = sorted(
        curves["threshold"]
        .dropna()
        .unique()
    )

    for threshold in thresholds:

        group = (
            curves[
                curves["threshold"]
                == threshold
            ]
            .sort_values("subsample")
        )

        if threshold == best_threshold:

            linewidth = 4
            alpha = 1.0

        else:

            linewidth = 2
            alpha = 0.65

        plt.plot(
            group["subsample"],
            group["expected_richness"],
            marker="o",
            markersize=3,
            linewidth=linewidth,
            alpha=alpha,
            label=f"Threshold {threshold}"
        )

    plt.title(
        f"{tissue} — Hill order 1"
    )

    plt.xlabel(
        "Subsample size"
    )

    plt.ylabel(
        "Hill-1 effective diversity"
    )

    plt.legend(
        title="Clone threshold",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
        fontsize="small"
    )

    plt.tight_layout(
        rect=[
            0,
            0,
            0.82,
            1
        ]
    )

    output_file = (
        outdir
        / f"{tissue}_hill1_rarefactioncurve.png"
    )

    plt.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()


# ============================================================
# Combined rarefaction plot
# ============================================================

def plot_combined_rarefaction(
    all_curves,
    outdir,
    best_threshold
):
    """
    Plot all tissue Hill-1 rarefaction curves.

    Tissue determines color.

    Selected global threshold is emphasized.
    """

    plt.figure(
        figsize=(11, 7)
    )

    for (
        tissue,
        threshold
    ), group in all_curves.groupby(
        [
            "tissue",
            "threshold"
        ]
    ):

        group = group.sort_values(
            "subsample"
        )

        color = TISSUE_COLORS.get(
            tissue,
            "gray"
        )

        if threshold == best_threshold:

            linewidth = 4
            alpha = 1.0

        else:

            linewidth = 1.5
            alpha = 0.35

        plt.plot(
            group["subsample"],
            group["expected_richness"],
            color=color,
            linewidth=linewidth,
            alpha=alpha
        )

    # --------------------------------------------------------
    # Tissue legend
    # --------------------------------------------------------

    handles = []

    present_tissues = (
        all_curves[
            "tissue"
        ]
        .dropna()
        .unique()
    )

    for tissue in TISSUE_COLORS:

        if tissue not in present_tissues:
            continue

        color = TISSUE_COLORS[
            tissue
        ]

        handles.append(
            plt.Line2D(
                [0],
                [0],
                color=color,
                lw=3,
                label=tissue
            )
        )

    # --------------------------------------------------------
    # Threshold legend
    # --------------------------------------------------------

    handles.append(
        plt.Line2D(
            [0],
            [0],
            color="black",
            lw=4,
            label=(
                f"Selected threshold = "
                f"{best_threshold}"
            )
        )
    )

    plt.xlabel(
        "Subsample size"
    )

    plt.ylabel(
        "Hill-1 effective diversity"
    )

    plt.title(
        "Hill-1 Rarefaction Across Tissues"
    )

    plt.legend(
        handles=handles,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False
    )

    plt.tight_layout(
        rect=[
            0,
            0,
            0.78,
            1
        ]
    )

    output_file = (
        outdir
        / "combined_hill1_rarefaction_curves.png"
    )

    plt.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Analyze Hill-1 rarefaction curves across "
            "clone thresholds and tissues. Select a "
            "single global clone threshold using a "
            "geometric elbow and calculate within-tissue "
            "diversity retention."
        )
    )

    parser.add_argument(
        "parent_directory",
        help=(
            "Parent directory containing one directory "
            "per tissue."
        )
    )

    parser.add_argument(
        "-o",
        "--output",
        default="rarefaction_output",
        help=(
            "Output directory. "
            "Default: rarefaction_output"
        )
    )

    args = parser.parse_args()

    parent = Path(
        args.parent_directory
    )

    outdir = Path(
        args.output
    )

    # --------------------------------------------------------
    # Validate input
    # --------------------------------------------------------

    if not parent.exists():

        print(
            f"ERROR: Directory does not exist: "
            f"{parent}"
        )

        return

    if not parent.is_dir():

        print(
            f"ERROR: Not a directory: "
            f"{parent}"
        )

        return

    # --------------------------------------------------------
    # Create output directory
    # --------------------------------------------------------

    outdir.mkdir(
        parents=True,
        exist_ok=True
    )

    # --------------------------------------------------------
    # Process tissues
    # --------------------------------------------------------

    all_curves = []

    tissue_curves = {}

    for tissue_dir in sorted(
        parent.iterdir()
    ):

        if not tissue_dir.is_dir():
            continue

        curves = process_tissue(
            tissue_dir,
            outdir
        )

        if curves is not None:

            all_curves.append(
                curves
            )

            tissue_curves[
                clean_tissue_name(
                    tissue_dir.name
                )
            ] = curves

    # --------------------------------------------------------
    # Check for data
    # --------------------------------------------------------

    if len(all_curves) == 0:

        print(
            "\nERROR: No valid Hill-1 "
            "rarefaction data found."
        )

        return

    # --------------------------------------------------------
    # Combine all tissues
    # --------------------------------------------------------

    all_curves = pd.concat(
        all_curves,
        ignore_index=True
    )

    all_curves.to_csv(
        outdir
        / "all_hill1_rarefaction_curves.csv",
        index=False
    )

    # --------------------------------------------------------
    # Select global threshold
    # --------------------------------------------------------

    best_threshold, threshold_summary = (
        choose_global_threshold(
            all_curves
        )
    )

    # --------------------------------------------------------
    # Save threshold summary
    # --------------------------------------------------------

    threshold_summary.to_csv(
        outdir
        / "threshold_elbow_summary.csv",
        index=False
    )

    # --------------------------------------------------------
    # Save selected threshold
    # --------------------------------------------------------

    threshold_file = (
        outdir
        / "recommended_global_threshold.txt"
    )

    with open(
        threshold_file,
        "w"
    ) as f:

        f.write(
            str(best_threshold)
        )

    # --------------------------------------------------------
    # Main two-panel figure
    # --------------------------------------------------------

    plot_threshold_analysis(
        all_curves,
        best_threshold,
        outdir
    )

    # --------------------------------------------------------
    # Print retention at selected threshold
    # --------------------------------------------------------

    print_threshold_results(
        all_curves,
        best_threshold
    )

    # --------------------------------------------------------
    # Tissue-specific rarefaction plots
    # --------------------------------------------------------

    for tissue, curves in (
        tissue_curves.items()
    ):

        plot_rarefaction(
            curves,
            tissue,
            outdir,
            best_threshold
        )

    # --------------------------------------------------------
    # Combined rarefaction plot
    # --------------------------------------------------------

    plot_combined_rarefaction(
        all_curves,
        outdir,
        best_threshold
    )

    # --------------------------------------------------------
    # Final summary
    # --------------------------------------------------------

    print(
        "\n======================================="
    )

    print(
        "Hill-1 threshold analysis complete"
    )

    print(
        "======================================="
    )

    print(
        f"\nGLOBAL CLONE THRESHOLD: "
        f"{best_threshold}"
    )

    print(
        f"\nOutput directory:"
        f"\n  {outdir}"
    )

    print(
        "\nFiles created:"
    )

    print(
        "  clone_threshold_analysis.png"
    )

    print(
        "  threshold_elbow_summary.csv"
    )

    print(
        "  threshold_diversity_retention.csv"
    )

    print(
        "  threshold_retention_summary.csv"
    )

    print(
        "  recommended_global_threshold.txt"
    )

    print(
        "  all_hill1_rarefaction_curves.csv"
    )

    print(
        "  combined_hill1_rarefaction_curves.png"
    )

    print(
        "  <tissue>_hill1_rarefaction.csv"
    )

    print(
        "  <tissue>_hill1_rarefactioncurve.png"
    )

    print(
        "\nAnalysis complete."
    )


# ============================================================
# Run
# ============================================================

if __name__ == "__main__":

    main()