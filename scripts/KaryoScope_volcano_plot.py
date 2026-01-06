#!/usr/bin/env python3
"""
KaryoScope_volcano_plot.py - Generate volcano plot from cluster analysis results.

Creates a volcano plot showing -log10(q-value) vs log2(odds ratio) for each cluster,
with optional labels from a curation file.
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adjustText import adjust_text


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate volcano plot from cluster analysis results")

    parser.add_argument("--cluster-analysis", dest="cluster_analysis", required=True,
                        help="Path to cluster_analysis.tsv file")
    parser.add_argument("--curation", dest="curation", default=None,
                        help="Path to curation file (TSV or Excel) with cluster labels")
    parser.add_argument("--output", dest="output", required=True,
                        help="Output path for SVG file")
    parser.add_argument("--label-column", dest="label_column", default="curated_annotation",
                        help="Column name containing labels in curation file (default: curated_annotation)")
    parser.add_argument("--size-scale", dest="size_scale", type=float, default=0.05,
                        help="Scaling factor for point sizes (default: 0.05)")
    parser.add_argument("--q-threshold", dest="q_threshold", type=float, default=0.05,
                        help="Q-value threshold for significance line (default: 0.05)")
    parser.add_argument("--figsize", dest="figsize", default="10,8",
                        help="Figure size as 'width,height' in inches (default: 10,8)")

    return parser.parse_args()


def load_cluster_analysis(path):
    """Load cluster analysis TSV file."""
    df = pd.read_csv(path, sep='\t')
    print(f"Loaded {len(df)} clusters from {path}")
    return df


def load_curation_labels(path, label_column):
    """Load curation labels from TSV or Excel file."""
    if path is None:
        return {}

    if path.endswith('.xlsx') or path.endswith('.xls'):
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path, sep='\t')

    labels = {}
    if 'cluster_id' in df.columns and label_column in df.columns:
        for _, row in df.iterrows():
            if pd.notna(row[label_column]) and str(row[label_column]).strip():
                labels[int(row['cluster_id'])] = str(row[label_column])

    print(f"Loaded {len(labels)} labels from {path}")
    return labels


def create_volcano_plot(df, labels, args):
    """Create volcano plot with labeled points."""

    # Parse figure size
    figsize = tuple(float(x) for x in args.figsize.split(','))
    fig, ax = plt.subplots(figsize=figsize)

    # Prepare data
    # Handle odds_ratio = 0 (avoid log2(0))
    df = df.copy()
    df['log2_or'] = np.log2(df['odds_ratio'].replace(0, 1e-10))
    df['-log10_q'] = -np.log10(df['q_value'].replace(0, 1e-300))

    # Define colors for enrichment categories
    color_map = {
        'Normal-enriched': '#3b82f6',  # blue
        'Tumor-enriched': '#ef4444',   # red
        'mixed': '#9ca3af'             # gray
    }

    # Plot points by enrichment category
    for enrichment, color in color_map.items():
        mask = df['enrichment'] == enrichment
        subset = df[mask]
        if len(subset) > 0:
            sizes = subset['size'] * args.size_scale
            sizes = sizes.clip(lower=20, upper=500)  # reasonable size range
            ax.scatter(
                subset['log2_or'],
                subset['-log10_q'],
                c=color,
                s=sizes,
                alpha=0.6,
                edgecolors='white',
                linewidths=0.5,
                label=f"{enrichment} (n={len(subset)})"
            )

    # Add significance threshold line
    if args.q_threshold > 0:
        threshold_y = -np.log10(args.q_threshold)
        ax.axhline(threshold_y, color='black', linestyle='--', alpha=0.5, linewidth=1)
        ax.text(ax.get_xlim()[1], threshold_y, f' q={args.q_threshold}',
                va='center', ha='left', fontsize=8, alpha=0.7)

    # Add vertical line at x=0
    ax.axvline(0, color='black', linestyle='-', alpha=0.2, linewidth=1)

    # Add labels for curated clusters
    texts = []
    for cluster_id, label in labels.items():
        row = df[df['cluster_id'] == cluster_id]
        if len(row) > 0:
            x = row['log2_or'].iloc[0]
            y = row['-log10_q'].iloc[0]
            texts.append(ax.text(x, y, label, fontsize=7, ha='left', va='bottom'))

    # Adjust text positions to avoid overlap
    if texts:
        adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', alpha=0.5, lw=0.5))

    # Labels and title
    ax.set_xlabel('log₂(Odds Ratio)', fontsize=11)
    ax.set_ylabel('-log₁₀(q-value)', fontsize=11)
    ax.set_title('Cluster Enrichment Volcano Plot', fontsize=13)

    # Legend
    ax.legend(loc='upper left', framealpha=0.9, fontsize=9)

    # Grid
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)

    # Tight layout
    plt.tight_layout()

    return fig


def main():
    args = parse_args()

    # Load data
    df = load_cluster_analysis(args.cluster_analysis)
    labels = load_curation_labels(args.curation, args.label_column)

    # Create plot
    fig = create_volcano_plot(df, labels, args)

    # Save
    output_path = args.output
    if not output_path.endswith('.svg'):
        output_path += '.svg'

    fig.savefig(output_path, format='svg', bbox_inches='tight')
    print(f"Saved volcano plot to {output_path}")

    # Also save PNG for quick preview
    png_path = output_path.replace('.svg', '.png')
    fig.savefig(png_path, format='png', dpi=150, bbox_inches='tight')
    print(f"Saved PNG preview to {png_path}")

    plt.close(fig)


if __name__ == '__main__':
    main()
