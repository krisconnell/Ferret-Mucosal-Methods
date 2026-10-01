#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import os
from collections import defaultdict
from itertools import combinations_with_replacement

import holoviews as hv
from holoviews import opts
from matplotlib.colors import to_rgb

hv.extension("matplotlib")

# --------------------
# Global color map
# --------------------
HEX_COLORS = {
    'PBMC': '#E8192B',
    'LN': '#DF6DFF',
    'Spleen': '#7C0002',
    'BAL': '#53006F',
    'Tonsil': '#FF8000',
    'NasalTurb': '#00C000',
    'NasalCell': "#C0C000",
    'Oral': "#55A0FB",
    'Lung': '#C700CF',
    'Intestine': '#8D8DFF'
}


def blend_colors(c1, c2):
    """Blend two hex colors evenly."""
    r1, g1, b1 = to_rgb(c1)
    r2, g2, b2 = to_rgb(c2)
    return ( (r1 + r2) / 2, (g1 + g2) / 2, (b1 + b2) / 2 )


def main():
    parser = argparse.ArgumentParser(description="Stable tissue-colored chord diagram")
    parser.add_argument("input_file", help="Input TSV file")
    parser.add_argument("-o", "--output", default=None, help="Output SVG filename")
    parser.add_argument("--read_counts", action="store_true",
                        help="Weight edges by nt_seq_count")
    parser.add_argument("--dpi", type=int, default=300)

    args = parser.parse_args()

    # --------------------
    # Load data
    # --------------------
    df = pd.read_csv(args.input_file, sep="\t")
    df["tissue"] = df["tissue"].astype(str)

    # --------------------
    # Build node list FIRST (deterministic)
    # --------------------
    tissues = sorted(df["tissue"].unique())
    node_index = {t: i for i, t in enumerate(tissues)}

    node_df = pd.DataFrame({
        "index": [node_index[t] for t in tissues],
        "name": tissues,
        "color": [HEX_COLORS.get(t, "#999999") for t in tissues]
    })

    # 🔒 CRITICAL: pandas index MUST be numeric and identical
    node_df.index = node_df["index"].astype(int)
    node_df = node_df.sort_index()

    # --------------------
    # Build edges
    # --------------------
    connections = defaultdict(float)

    grouped = df.groupby("ClusterID")

    if args.read_counts:
        for _, g in grouped:
            counts = g.groupby("tissue")["nt_seq_count"].sum()
            ts = counts.index.tolist()
            for a, b in combinations_with_replacement(ts, 2):
                w = counts[a] * counts[b]
                connections[(a, b)] += w
                if a != b:
                    connections[(b, a)] += w
    else:
        for _, g in grouped:
            ts = g["tissue"].unique()
            for a, b in combinations_with_replacement(ts, 2):
                connections[(a, b)] += 1
                if a != b:
                    connections[(b, a)] += 1

    edge_rows = []
    for (src, tgt), val in connections.items():
        edge_rows.append({
            "source_name": src,
            "target_name": tgt,
            "value": val
        })

    edge_df = pd.DataFrame(edge_rows)

    # Map to numeric indices ONLY
    edge_df["source"] = edge_df["source_name"].map(node_index).astype(int)
    edge_df["target"] = edge_df["target_name"].map(node_index).astype(int)

    # --------------------
    # Edge coloring (blend for shared)
    # --------------------

    def blend_hex(c1, c2):
        r1, g1, b1 = to_rgb(c1)
        r2, g2, b2 = to_rgb(c2)
        r, g, b = (r1+r2)/2, (g1+g2)/2, (b1+b2)/2
        return "#{:02x}{:02x}{:02x}".format(
            int(r*255), int(g*255), int(b*255)
        )


    edge_colors = []
    for _, r in edge_df.iterrows():
        c1 = HEX_COLORS.get(r["source_name"], "#999999")
        c2 = HEX_COLORS.get(r["target_name"], "#999999")
        if r["source_name"] == r["target_name"]:
            edge_colors.append(c1)          # hex string
        else:
            edge_colors.append(blend_hex(c1, c2))

    edge_df["edge_color"] = edge_colors

    # --------------------
    # Normalize linewidths
    # --------------------
    v = edge_df["value"].values
    vmin, vmax = v.min(), v.max()
    if vmax > vmin:
        edge_df["linewidth"] = 0.5 + 4.5 * (v - vmin) / (vmax - vmin)
    else:
        edge_df["linewidth"] = 2.0

    # 🔒 Strip to ONLY what HoloViews should see
    edge_df = edge_df[["source", "target", "value", "linewidth", "edge_color"]]
    edge_df = edge_df.reset_index(drop=True)

    # --------------------
    # HoloViews objects
    # --------------------
    nodes = hv.Dataset(
        node_df,
        kdims=[hv.Dimension("index", type=int)],
        vdims=["name", "color"]
    )

    def fix_node_labels(plot, element):
        """
        Make chord node labels upright and readable.
        """
        labels = plot.handles.get("labels", [])
        for lbl in labels:
            lbl.set_rotation(0)          # no rotation
            lbl.set_fontsize(20)         # requested size
            lbl.set_ha("center")
            lbl.set_va("center")

            # push label slightly outward from the circle
            x, y = lbl.get_position()
            lbl.set_position((x * 1.10, y * 1.10))

    chord = hv.Chord((edge_df, nodes)).opts(
        opts.Chord(
            labels="name",
            node_color="color",
            edge_color="edge_color",
            edge_linewidth="linewidth",
            colorbar=False,
            padding=0.03,      
            node_linewidth=1.2,
            node_edgecolors="white",
            show_frame=False,
            fig_inches=(16, 16),
            title="Tissue Sharing of BCR Lineages",
            hooks=[fix_node_labels]  
        )
    )

    # --------------------
    # Output
    # --------------------
    if args.output:
        out = args.output
    else:
        base = os.path.splitext(os.path.basename(args.input_file))[0]
        out = f"{base}_chord.svg"

    assert isinstance(edge_df["edge_color"].iloc[0], str)
    hv.save(chord, out, fmt="svg", dpi=args.dpi)
    print(f"Saved chord diagram → {out}")


if __name__ == "__main__":
    main()
