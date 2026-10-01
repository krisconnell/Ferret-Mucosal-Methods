#!/usr/bin/env python3

import argparse
import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from matplotlib_venn import venn3


# --------------------------------------------------
# Define anatomical compartments
# --------------------------------------------------

REGIONS = {
    "URT": [
        "NasalCell",
        "Oral",
        "NasalTurb",
        "Tonsil"
    ],

    "LRT": [
        "LN",
        "Lung",
        "BAL"
    ],

    "Systemic": [
        "PBMC",
        "Spleen"
    ]
}


# --------------------------------------------------
# Build cluster sets per Tissue / isotype
# --------------------------------------------------

# python isotype_overlap.py \
#     clustered_database.tsv \
#     -o IsotypeOverlap \
#     --min_cluster_size 2 \
#     --plot_mode percent


def load_and_filter(df, min_cluster_size=1):

    required = {
        "ClusterID",
        "Tissue",
        "c_call",
        "nt_seq_count"
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            f"Missing required columns: {missing}"
        )

    df = df.dropna(
        subset=["ClusterID", "Tissue", "c_call"]
    )

    df["nt_seq_count"] = pd.to_numeric(
        df["nt_seq_count"],
        errors="coerce"
    ).fillna(0)

    # Global cluster abundance filter
    cluster_sizes = (
        df.groupby("ClusterID")["nt_seq_count"]
        .sum()
    )

    keep_clusters = cluster_sizes[
        cluster_sizes >= min_cluster_size
    ].index

    df = df[
        df["ClusterID"].isin(keep_clusters)
    ].copy()

    return df


# --------------------------------------------------
# Calculate Venn regions
# --------------------------------------------------

def calculate_overlap_regions(igm, igg, iga):

    igm = set(igm)
    igg = set(igg)
    iga = set(iga)

    all3 = igm & igg & iga

    igm_igg = (igm & igg) - all3
    igm_iga = (igm & iga) - all3
    igg_iga = (igg & iga) - all3

    igm_only = igm - igg - iga
    igg_only = igg - igm - iga
    iga_only = iga - igm - igg

    total_unique = len(
        igm | igg | iga
    )

    return {
        "IgM_only": len(igm_only),
        "IgG_only": len(igg_only),
        "IgA_only": len(iga_only),
        "IgM_IgG": len(igm_igg),
        "IgM_IgA": len(igm_iga),
        "IgG_IgA": len(igg_iga),
        "All_3": len(all3),
        "Total_unique": total_unique
    }


# --------------------------------------------------
# Make Venn diagram
# --------------------------------------------------

def make_venn(
    Tissue,
    counts,
    output_dir,
    plot_mode="counts"
):

    total = counts["Total_unique"]

    subsets_counts = (
        counts["IgM_only"],
        counts["IgG_only"],
        counts["IgM_IgG"],
        counts["IgA_only"],
        counts["IgM_IgA"],
        counts["IgG_IgA"],
        counts["All_3"]
    )

    plt.figure(figsize=(8, 8))

    venn = venn3(
        subsets=subsets_counts,
        set_labels=("IgM", "IgG", "IgA")
    )

    # --------------------------------------------------
    # Font sizes
    # --------------------------------------------------

    SET_LABEL_FONTSIZE = 24
    REGION_LABEL_FONTSIZE = 20
    TITLE_FONTSIZE = 22

    # Set labels: IgM, IgG, IgA
    for label in venn.set_labels:
        if label:
            label.set_fontsize(SET_LABEL_FONTSIZE)
            label.set_fontweight("bold")

    # --------------------------------------------------
    # Region labels
    # --------------------------------------------------

    labels = {
        "100": counts["IgM_only"],
        "010": counts["IgG_only"],
        "110": counts["IgM_IgG"],
        "001": counts["IgA_only"],
        "101": counts["IgM_IgA"],
        "011": counts["IgG_IgA"],
        "111": counts["All_3"]
    }

    for region, value in labels.items():

        label = venn.get_label_by_id(region)

        if label is not None:

            if plot_mode == "percent":

                pct = (
                    value / total * 100
                    if total > 0 else 0
                )

                label.set_text(
                    f"{pct:.1f}%"
                )

            else:

                label.set_text(
                    f"{value}"
                )

            label.set_fontsize(REGION_LABEL_FONTSIZE)
            label.set_fontweight("bold")

    # --------------------------------------------------
    # Title
    # --------------------------------------------------

    plt.title(
        f"{Tissue}\n({plot_mode})",
        fontsize=TITLE_FONTSIZE,
        fontweight="bold"
    )

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            output_dir,
            f"{Tissue}_venn_{plot_mode}.png"
        ),
        dpi=300,
        bbox_inches="tight"
    )

    plt.savefig(
        os.path.join(
            output_dir,
            f"{Tissue}_venn_{plot_mode}.pdf"
        ),
        bbox_inches="tight"
    )

    plt.close()
# --------------------------------------------------
# Analyze a tissue
# --------------------------------------------------

def analyze_Tissue(df_Tissue):

    igm = set(
        df_Tissue[
            df_Tissue["c_call"] == "IGHM*01"
        ]["ClusterID"]
    )

    igg = set(
        df_Tissue[
            df_Tissue["c_call"] == "IGHG1*02,IGHG2*01"
        ]["ClusterID"]
    )

    iga = set(
        df_Tissue[
            df_Tissue["c_call"] == "IGHA*01"
        ]["ClusterID"]
    )

    return calculate_overlap_regions(
        igm,
        igg,
        iga
    )


# --------------------------------------------------
# Analyze a regional compartment
# --------------------------------------------------

def analyze_region(df, region_tissues):

    region_df = df[
        df["Tissue"].isin(region_tissues)
    ]

    return analyze_Tissue(region_df)


# --------------------------------------------------
# Main
# --------------------------------------------------

def main():

    parser = argparse.ArgumentParser(
        description="Isotype overlap analysis"
    )

    parser.add_argument(
        "input_tsv"
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True
    )

    parser.add_argument(
        "--min_cluster_size",
        type=int,
        default=1,
        help="Global clone abundance threshold"
    )

    parser.add_argument(
        "--plot_mode",
        choices=["counts", "percent"],
        default="counts"
    )

    args = parser.parse_args()

    os.makedirs(
        args.output,
        exist_ok=True
    )

    df = pd.read_csv(
        args.input_tsv,
        sep="\t",
        low_memory=False
    )

    df = load_and_filter(
        df,
        args.min_cluster_size
    )

    summary_rows = []


    # ==================================================
    # 1. Original tissue-level analysis
    # ==================================================

    for Tissue in sorted(
        df["Tissue"].unique()
    ):

        Tissue_df = df[
            df["Tissue"] == Tissue
        ]

        counts = analyze_Tissue(
            Tissue_df
        )

        total_unique = counts["Total_unique"]

        shared_any = (
            counts["IgM_IgG"]
            + counts["IgM_IgA"]
            + counts["IgG_IgA"]
            + counts["All_3"]
        )

        pct_shared_any = (
            shared_any / total_unique * 100
            if total_unique > 0 else 0
        )

        summary_rows.append({

            "Tissue": Tissue,

            "Region": "Individual_Tissue",

            **counts,

            "Shared_with_other_isotype":
                shared_any,

            "Pct_shared_with_other_isotype":
                pct_shared_any

        })

        make_venn(
            Tissue,
            counts,
            args.output,
            args.plot_mode
        )


    # ==================================================
    # 2. Regional analysis
    # ==================================================

    for region, tissues in REGIONS.items():

        # Only analyze tissues that are actually present
        present_tissues = [
            tissue
            for tissue in tissues
            if tissue in df["Tissue"].unique()
        ]

        if not present_tissues:

            print(
                f"WARNING: No tissues found for "
                f"region {region}"
            )

            continue

        counts = analyze_region(
            df,
            present_tissues
        )

        total_unique = counts["Total_unique"]

        shared_any = (
            counts["IgM_IgG"]
            + counts["IgM_IgA"]
            + counts["IgG_IgA"]
            + counts["All_3"]
        )

        pct_shared_any = (
            shared_any / total_unique * 100
            if total_unique > 0 else 0
        )

        summary_rows.append({

            "Tissue": region,

            "Region": region,

            **counts,

            "Shared_with_other_isotype":
                shared_any,

            "Pct_shared_with_other_isotype":
                pct_shared_any

        })

        make_venn(
            region,
            counts,
            args.output,
            args.plot_mode
        )


    # ==================================================
    # 3. Save summary
    # ==================================================

    summary_df = pd.DataFrame(
        summary_rows
    )

    summary_df.to_csv(
        os.path.join(
            args.output,
            "isotype_overlap_summary.tsv"
        ),
        sep="\t",
        index=False
    )


    # ==================================================
    # 4. Percent table
    # ==================================================

    percent_df = summary_df.copy()

    region_cols = [
        "IgM_only",
        "IgG_only",
        "IgA_only",
        "IgM_IgG",
        "IgM_IgA",
        "IgG_IgA",
        "All_3",
        "Shared_with_other_isotype"
    ]

    for col in region_cols:

        pct_col = f"Pct_{col}"

        percent_df[pct_col] = np.where(
            percent_df["Total_unique"] > 0,
            percent_df[col]
            / percent_df["Total_unique"]
            * 100,
            0
        )

        percent_df[pct_col] = percent_df[pct_col].round(2)


    percent_df.to_csv(
        os.path.join(
            args.output,
            "isotype_overlap_percentages.tsv"
        ),
        sep="\t",
        index=False
    )
    
    # ==================================================
    # 5. Save regional-only summary
    # ==================================================

    regional_df = summary_df[
        summary_df["Region"].isin(
            REGIONS.keys()
        )
    ].copy()

    regional_df.to_csv(
        os.path.join(
            args.output,
            "regional_isotype_overlap_summary.tsv"
        ),
        sep="\t",
        index=False
    )


    print()
    print(
        f"Finished analysis for "
        f"{len(summary_df)} total entries"
    )

    print(
        f"  Individual tissues: "
        f"{len(df['Tissue'].unique())}"
    )

    print(
        f"  Regional compartments: "
        f"{len(REGIONS)}"
    )

    print(
        f"Results written to: {args.output}"
    )

    print()


if __name__ == "__main__":
    main()
