#!/usr/bin/env python3
"""
KaryoScope_enrichment_bubbles.py - Vertical bubble plot showing cluster enrichment.

Creates a vertical bubble plot with clusters as rows and groups as columns,
showing enrichment/depletion with bubble size based on significance.
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Keep text editable in output files (not converted to paths)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate vertical enrichment bubble plot for curated clusters")

    parser.add_argument("--cluster-analysis", dest="cluster_analysis", required=True,
                        help="Path to cluster_analysis.tsv file")
    parser.add_argument("--curation", dest="curation", required=True,
                        help="Path to curation file (TSV or Excel) with cluster labels")
    parser.add_argument("--output", dest="output", required=True,
                        help="Output path for output file")
    parser.add_argument("--label-column", dest="label_column", default="curated_annotation",
                        help="Column name containing labels in curation file (default: curated_annotation)")
    parser.add_argument("--dark-mode", dest="dark_mode", action="store_true",
                        help="Use dark background")
    parser.add_argument("--figsize", dest="figsize", default="6,12",
                        help="Figure size as 'width,height' in inches (default: 6,12)")

    return parser.parse_args()


def load_cluster_analysis(path):
    """Load cluster analysis TSV file."""
    df = pd.read_csv(path, sep='\t')
    print(f"Loaded {len(df)} clusters from {path}")
    return df


def load_curation_labels(path, label_column):
    """Load curation labels from TSV or Excel file."""
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


def create_vertical_bubble_plot(df, labels, args):
    """Create vertical bubble plot with clusters as rows, groups as columns."""

    # Filter to only curated clusters
    curated_ids = list(labels.keys())
    df_curated = df[df['cluster_id'].isin(curated_ids)].copy()

    if len(df_curated) == 0:
        print("No curated clusters found in data!")
        return None

    # Add labels
    df_curated['label'] = df_curated['cluster_id'].map(labels)

    # Find group columns (Normal_pct, Tumor_pct, etc.)
    pct_cols = [c for c in df_curated.columns if c.endswith('_pct')]
    groups = [c.replace('_pct', '') for c in pct_cols]

    if len(groups) == 0:
        print("No group percentage columns found!")
        return None

    print(f"Found groups: {groups}")

    # Sort clusters by enrichment pattern
    df_curated = df_curated.sort_values(['enrichment', 'q_value'], ascending=[True, True])

    n_clusters = len(df_curated)
    n_groups = len(groups)

    # Set up colors
    if args.dark_mode:
        plt.style.use('dark_background')
        bg_color = 'black'
        text_color = 'white'
        grid_color = '#444444'
    else:
        bg_color = 'white'
        text_color = 'black'
        grid_color = '#cccccc'

    # Parse figure size
    figsize = tuple(float(x) for x in args.figsize.split(','))
    fig, ax = plt.subplots(figsize=figsize)

    if args.dark_mode:
        fig.patch.set_facecolor(bg_color)
        ax.set_facecolor(bg_color)

    # Calculate expected percentages (from total counts)
    # Assuming Normal and Tumor are the main groups
    expected_pcts = {}
    for g in groups:
        count_col = f'{g}_count'
        if count_col in df.columns:
            total_in_group = df[count_col].sum()
            total_all = sum(df[f'{grp}_count'].sum() for grp in groups if f'{grp}_count' in df.columns)
            expected_pcts[g] = (total_in_group / total_all * 100) if total_all > 0 else 50
        else:
            expected_pcts[g] = 100 / len(groups)

    print(f"Expected percentages: {expected_pcts}")

    # Create bubble plot: rows = clusters, columns = groups
    for i, (_, row) in enumerate(df_curated.iterrows()):
        cluster_qval = row.get('q_value', 1.0)
        enrichment = row.get('enrichment', 'mixed')

        for j, group in enumerate(groups):
            observed_pct = row.get(f'{group}_pct', 0)
            expected_pct = expected_pcts.get(group, 50)

            # Compute odds ratio as observed/expected
            if expected_pct > 0:
                odds = observed_pct / expected_pct
            else:
                odds = 1.0

            # Size based on significance (for this cluster overall)
            neg_log_p = -np.log10(cluster_qval) if cluster_qval > 0 else 10
            neg_log_p = min(neg_log_p, 10)  # Cap at 10
            size = neg_log_p * 80 + 20

            # Color by odds ratio: red = enriched (>1), blue = depleted (<1)
            if odds > 1:
                intensity = min(1.0, np.log2(odds) / 2)
                color = plt.cm.Reds(0.3 + intensity * 0.7)
            else:
                intensity = min(1.0, -np.log2(max(odds, 0.01)) / 2)
                color = plt.cm.Blues(0.3 + intensity * 0.7)

            # Significant edge
            if cluster_qval < 0.05:
                edgecolor = 'white' if args.dark_mode else 'black'
                linewidth = 1.5
            else:
                edgecolor = 'gray'
                linewidth = 0.5

            ax.scatter(j, i, s=size, c=[color], edgecolors=edgecolor, linewidths=linewidth)

    # Labels
    ax.set_xticks(range(n_groups))
    ax.set_xticklabels(groups, fontsize=11, color=text_color)
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')

    ax.set_yticks(range(n_clusters))
    ax.set_yticklabels(df_curated['label'].values, fontsize=9, color=text_color)

    ax.set_xlabel('Group', fontsize=12, color=text_color)
    ax.set_ylabel('Cluster', fontsize=12, color=text_color)

    ax.set_xlim(-0.5, n_groups - 0.5)
    ax.set_ylim(-0.5, n_clusters - 0.5)

    # Invert y-axis so first cluster is at top
    ax.invert_yaxis()

    # Grid
    ax.grid(True, alpha=0.3, color=grid_color)

    # Add legend for size
    legend_elements = []
    for pval, label in [(0.05, 'p=0.05'), (0.01, 'p=0.01'), (0.001, 'p=0.001')]:
        neg_log = -np.log10(pval)
        size = neg_log * 80 + 20
        legend_elements.append(plt.scatter([], [], s=size, c='gray', alpha=0.7, label=label))

    leg1 = ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1),
                     title='Significance', fontsize=8, title_fontsize=9)
    if args.dark_mode:
        leg1.get_frame().set_facecolor('#222222')

    # Add color legend
    from matplotlib.patches import Patch
    color_legend = [
        Patch(facecolor=plt.cm.Reds(0.7), label='Enriched'),
        Patch(facecolor=plt.cm.Blues(0.7), label='Depleted'),
    ]
    leg2 = ax.legend(handles=color_legend, loc='lower left', bbox_to_anchor=(1.02, 0),
                     title='Direction', fontsize=8, title_fontsize=9)
    if args.dark_mode:
        leg2.get_frame().set_facecolor('#222222')
    ax.add_artist(leg1)

    # Adjust spines
    for spine in ax.spines.values():
        spine.set_color(text_color if args.dark_mode else 'black')

    ax.tick_params(colors=text_color)

    plt.tight_layout()
    return fig


def main():
    args = parse_args()

    # Load data
    df = load_cluster_analysis(args.cluster_analysis)
    labels = load_curation_labels(args.curation, args.label_column)

    # Create plot
    fig = create_vertical_bubble_plot(df, labels, args)

    if fig is None:
        sys.exit(1)

    # Save
    output_path = args.output
    if not output_path.endswith('.svg') and not output_path.endswith('.pdf') and not output_path.endswith('.png'):
        output_path += '.svg'

    fig.savefig(output_path, format=output_path.split('.')[-1], bbox_inches='tight',
                facecolor=fig.get_facecolor())
    print(f"Saved bubble plot to {output_path}")

    # Also save PNG for quick preview
    png_path = output_path.rsplit('.', 1)[0] + '.png'
    fig.savefig(png_path, format='png', dpi=150, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    print(f"Saved PNG preview to {png_path}")

    plt.close(fig)


if __name__ == '__main__':
    main()
