#!/usr/bin/env python3
"""
KS_allchr_allele_heatmap.py

For chr3, chr8, chr11, chr12: pair haplotypes by sample (h1/h2),
build a co-occurrence matrix of cluster assignments, and plot as heatmaps.
This reveals whether outliers are allele-specific or biallelic.

Uses the same sil-threshold + centroid-scan logic as KS_allchr_barplot.py
to ensure cluster labels match the barplot.

Usage:
  python3 scripts/KS_allchr_allele_heatmap.py \
    --assignments agent_results/allchr_structure.sequence_assignments.tsv \
    --bed /path/to/pangenome.BED \
    --sil-threshold 0.5 --centroid-sd 5
"""

import matplotlib
matplotlib.rcParams['svg.fonttype'] = 'none'

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from KS_allchr_dendrogram import EXCLUDED_CHROMS, CHROM_ORDER
from KS_allchr_barplot import apply_sil_and_centroid


def simplify_cluster_name(cluster):
    """chr3_Major -> Major, chr3_Outlier_5 -> Outlier_5"""
    parts = cluster.split("_", 1)
    return parts[1] if len(parts) > 1 else cluster


def build_allele_matrix(df, chrom):
    """
    For a single chromosome, pair h1/h2 by sample and build co-occurrence matrix.
    Returns: matrix (DataFrame), paired_samples, unpaired, pairs
    """
    chrom_df = df[df["chromosome"] == chrom].copy()

    chrom_df["sample_id"] = chrom_df["sequence"].str.extract(r"^(.+?)#[12]#")[0]
    chrom_df["hap_num"] = chrom_df["sequence"].str.extract(r"#([12])#")[0]
    chrom_df["short_cluster"] = chrom_df["cluster"].apply(simplify_cluster_name)

    h1 = chrom_df[chrom_df["hap_num"] == "1"].set_index("sample_id")["short_cluster"]
    h2 = chrom_df[chrom_df["hap_num"] == "2"].set_index("sample_id")["short_cluster"]

    paired_samples = h1.index.intersection(h2.index)
    unpaired = set(h1.index.symmetric_difference(h2.index))

    pairs = pd.DataFrame({"h1": h1.loc[paired_samples], "h2": h2.loc[paired_samples]})

    all_clusters = sorted(
        chrom_df["short_cluster"].unique(),
        key=lambda x: (0, "") if x == "Major" else (1, x),
    )

    matrix = pd.DataFrame(0, index=all_clusters, columns=all_clusters)
    for _, row in pairs.iterrows():
        matrix.loc[row["h1"], row["h2"]] += 1

    return matrix, paired_samples, unpaired, pairs


def plot_heatmaps(matrices, output_dir):
    """Plot 2x2 grid of heatmaps for 4 chromosomes."""
    chroms = ["chr3", "chr8", "chr11", "chr12"]

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle(
        "Allele-specific cluster assignment (h1 vs h2)",
        fontsize=14,
        fontweight="bold",
        y=0.98,
    )

    for idx, chrom in enumerate(chroms):
        ax = axes[idx // 2][idx % 2]
        mat = matrices[chrom]["matrix"]
        n_paired = matrices[chrom]["n_paired"]

        # Symmetric matrix: combine (A,B) and (B,A) for display
        sym = mat.values + mat.values.T
        np.fill_diagonal(sym, np.diag(mat.values))
        sym_df = pd.DataFrame(sym, index=mat.index, columns=mat.columns)

        vmax = sym_df.values.max()
        if vmax == 0:
            vmax = 1

        im = ax.imshow(
            sym_df.values,
            cmap="YlOrRd",
            aspect="auto",
            vmin=0,
            vmax=vmax,
        )

        ax.set_xticks(range(len(sym_df.columns)))
        ax.set_xticklabels(sym_df.columns, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(len(sym_df.index)))
        ax.set_yticklabels(sym_df.index, fontsize=8)

        for i in range(len(sym_df.index)):
            for j in range(len(sym_df.columns)):
                val = sym_df.values[i, j]
                if val > 0:
                    color = "white" if val > vmax * 0.6 else "black"
                    ax.text(
                        j, i, str(int(val)),
                        ha="center", va="center",
                        fontsize=9, fontweight="bold",
                        color=color,
                    )

        ax.set_title(f"{chrom}  (n={n_paired} paired samples)", fontsize=11)
        ax.set_xlabel("Haplotype 2 cluster")
        ax.set_ylabel("Haplotype 1 cluster")

        plt.colorbar(im, ax=ax, shrink=0.8, label="# samples")

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    svg_path = os.path.join(output_dir, "allchr_allele_heatmap.svg")
    png_path = os.path.join(output_dir, "allchr_allele_heatmap.png")
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {svg_path}")
    print(f"Saved: {png_path}")
    return svg_path, png_path


def print_summary(matrices, output_dir="agent_results"):
    """Print summary statistics and save to TSV files."""
    # 1) Summary table
    summary_rows = []
    for chrom in ["chr3", "chr8", "chr11", "chr12"]:
        info = matrices[chrom]
        pairs = info["pairs"]
        n_paired = info["n_paired"]
        n_unpaired = len(info["unpaired"])

        both_major = ((pairs["h1"] == "Major") & (pairs["h2"] == "Major")).sum()
        both_outlier = ((pairs["h1"] != "Major") & (pairs["h2"] != "Major")).sum()
        h1_only_outlier = ((pairs["h1"] != "Major") & (pairs["h2"] == "Major")).sum()
        h2_only_outlier = ((pairs["h1"] == "Major") & (pairs["h2"] != "Major")).sum()
        monoallelic = h1_only_outlier + h2_only_outlier

        both_out = pairs[(pairs["h1"] != "Major") & (pairs["h2"] != "Major")]
        same_cluster = (both_out["h1"] == both_out["h2"]).sum() if len(both_out) > 0 else 0
        diff_cluster = (both_out["h1"] != both_out["h2"]).sum() if len(both_out) > 0 else 0

        summary_rows.append({
            "chromosome": chrom,
            "n_paired": n_paired,
            "n_unpaired": n_unpaired,
            "both_major": both_major,
            "h1_outlier_h2_major": h1_only_outlier,
            "h1_major_h2_outlier": h2_only_outlier,
            "both_outlier": both_outlier,
            "monoallelic_outlier": monoallelic,
            "biallelic_outlier": both_outlier,
            "biallelic_same_cluster": same_cluster,
            "biallelic_diff_cluster": diff_cluster,
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(output_dir, "allchr_allele_summary.tsv")
    summary_df.to_csv(summary_path, sep="\t", index=False)
    print(f"Saved: {summary_path}")

    # 2) Per-sample pair assignments
    pair_rows = []
    for chrom in ["chr3", "chr8", "chr11", "chr12"]:
        pairs = matrices[chrom]["pairs"]
        for sample_id, row in pairs.iterrows():
            pair_rows.append({
                "chromosome": chrom,
                "sample": sample_id,
                "h1_cluster": row["h1"],
                "h2_cluster": row["h2"],
            })
    pairs_df = pd.DataFrame(pair_rows)
    pairs_path = os.path.join(output_dir, "allchr_allele_pairs.tsv")
    pairs_df.to_csv(pairs_path, sep="\t", index=False)
    print(f"Saved: {pairs_path}")

    # 3) Co-occurrence matrices (one sheet per chrom, long format)
    cooccur_rows = []
    for chrom in ["chr3", "chr8", "chr11", "chr12"]:
        mat = matrices[chrom]["matrix"]
        for h1_cl in mat.index:
            for h2_cl in mat.columns:
                count = mat.loc[h1_cl, h2_cl]
                if count > 0:
                    cooccur_rows.append({
                        "chromosome": chrom,
                        "h1_cluster": h1_cl,
                        "h2_cluster": h2_cl,
                        "count": count,
                    })
    cooccur_df = pd.DataFrame(cooccur_rows)
    cooccur_path = os.path.join(output_dir, "allchr_allele_cooccurrence.tsv")
    cooccur_df.to_csv(cooccur_path, sep="\t", index=False)
    print(f"Saved: {cooccur_path}")

    # Print to stdout as before
    for chrom in ["chr3", "chr8", "chr11", "chr12"]:
        info = matrices[chrom]
        mat = info["matrix"]
        pairs = info["pairs"]
        n_paired = info["n_paired"]

        print(f"\n{'='*60}")
        print(f"{chrom}: {n_paired} paired samples, {len(info['unpaired'])} unpaired")
        print(f"{'='*60}")

        both_major = ((pairs["h1"] == "Major") & (pairs["h2"] == "Major")).sum()
        both_outlier = ((pairs["h1"] != "Major") & (pairs["h2"] != "Major")).sum()
        h1_only_outlier = ((pairs["h1"] != "Major") & (pairs["h2"] == "Major")).sum()
        h2_only_outlier = ((pairs["h1"] == "Major") & (pairs["h2"] != "Major")).sum()

        print(f"  Both Major:           {both_major:>4}  ({100*both_major/n_paired:.1f}%)")
        print(f"  h1 Outlier, h2 Major: {h1_only_outlier:>4}  ({100*h1_only_outlier/n_paired:.1f}%)")
        print(f"  h1 Major, h2 Outlier: {h2_only_outlier:>4}  ({100*h2_only_outlier/n_paired:.1f}%)")
        print(f"  Both Outlier:         {both_outlier:>4}  ({100*both_outlier/n_paired:.1f}%)")

        one_outlier = h1_only_outlier + h2_only_outlier
        print(f"  --> Monoallelic outlier: {one_outlier} ({100*one_outlier/n_paired:.1f}%)")
        print(f"  --> Biallelic outlier:   {both_outlier} ({100*both_outlier/n_paired:.1f}%)")

        both_out = pairs[(pairs["h1"] != "Major") & (pairs["h2"] != "Major")]
        if len(both_out) > 0:
            same_cluster = (both_out["h1"] == both_out["h2"]).sum()
            diff_cluster = (both_out["h1"] != both_out["h2"]).sum()
            print(f"\n  Among biallelic outliers:")
            print(f"    Same cluster:      {same_cluster}")
            print(f"    Different cluster: {diff_cluster}")

        print(f"\n  Co-occurrence matrix (h1 rows x h2 cols):")
        print(mat.to_string(col_space=8))


def parse_args():
    parser = argparse.ArgumentParser(description="Allele-specific cluster heatmap")
    parser.add_argument(
        "--assignments", required=True,
        help="Path to cluster assignment TSV",
    )
    parser.add_argument(
        "--bed", required=True, nargs="+",
        help="BED file(s) with feature annotations",
    )
    parser.add_argument(
        "--output-dir", dest="output_dir",
        default="agent_results",
        help="Output directory",
    )
    parser.add_argument(
        "--exclude-features", dest="exclude_features",
        default="novel",
        help="Comma-separated features to exclude (default: novel)",
    )
    parser.add_argument(
        "--edge-mode", dest="edge_mode", default="directional",
        help="Edge counting mode (default: directional)",
    )
    parser.add_argument(
        "--matrix-type", dest="matrix_type",
        default="count_log1p_zscore",
        help="Matrix normalization (default: count_log1p_zscore)",
    )
    parser.add_argument(
        "--sil-threshold", dest="sil_threshold", type=float, default=0.0,
        help="Silhouette threshold for collapsing weak splits (default: 0.0)",
    )
    parser.add_argument(
        "--centroid-sd", dest="centroid_sd", type=float, default=3.0,
        help="SD threshold for centroid scan (default: 3.0)",
    )
    parser.add_argument(
        "--centroid-scan", dest="centroid_scan",
        action=argparse.BooleanOptionalAction, default=True,
        help="Enable stage 2 centroid scan (default: on)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    asgn = pd.read_csv(args.assignments, sep="\t")
    asgn = asgn[~asgn["chromosome"].isin(EXCLUDED_CHROMS)]
    print(f"Loaded {len(asgn)} haplotypes across {asgn['chromosome'].nunique()} chromosomes")

    # Apply same sil-threshold + centroid-scan as barplot
    asgn = apply_sil_and_centroid(asgn, args)

    matrices = {}
    for chrom in ["chr3", "chr8", "chr11", "chr12"]:
        mat, paired, unpaired, pairs = build_allele_matrix(asgn, chrom)
        matrices[chrom] = {
            "matrix": mat,
            "n_paired": len(paired),
            "unpaired": unpaired,
            "pairs": pairs,
        }

    print_summary(matrices, args.output_dir)
    plot_heatmaps(matrices, args.output_dir)


if __name__ == "__main__":
    main()
