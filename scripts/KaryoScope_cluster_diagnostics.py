#!/usr/bin/env python3
"""
KaryoScope Cluster Diagnostics

Creates diagnostic plots comparing clusters across various metrics from
annotated sequence data.

Usage:
  python KaryoScope_cluster_diagnostics.py \
    --annotated analysis.read_assignments.annotated.tsv \
    --output-prefix analysis_diagnostics

Generates:
  - {prefix}.cluster_metrics.pdf: Box/violin plots of metrics by cluster
  - {prefix}.cluster_composition.pdf: Sample/group composition per cluster
  - {prefix}.cluster_summary.tsv: Summary statistics per cluster
"""

import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats


def load_annotated_data(filepath):
    """Load annotated sequence data."""
    print(f"Loading annotated data: {filepath}")
    df = pd.read_csv(filepath, sep='\t')
    print(f"  Total sequences: {len(df)}")
    print(f"  Columns: {', '.join(df.columns)}")
    print(f"  Clusters: {df['cluster'].nunique()}")
    return df


def plot_metric_by_cluster(df, metric, ax, title=None, ylabel=None):
    """Create a box plot of a metric by cluster with Kruskal-Wallis test."""
    if metric not in df.columns:
        ax.text(0.5, 0.5, f"Column '{metric}' not found", ha='center', va='center')
        ax.set_title(title or metric)
        return

    # Sort clusters by median value
    cluster_order = df.groupby('cluster')[metric].median().sort_values(ascending=False).index

    sns.boxplot(data=df, x='cluster', y=metric, order=cluster_order, ax=ax,
                palette='viridis', showfliers=False)
    sns.stripplot(data=df, x='cluster', y=metric, order=cluster_order, ax=ax,
                  color='black', alpha=0.3, size=2)

    # Kruskal-Wallis test
    groups = [group[metric].dropna().values for name, group in df.groupby('cluster')]
    groups = [g for g in groups if len(g) > 0]
    if len(groups) >= 2:
        stat, p_value = stats.kruskal(*groups)
        p_str = f"p={p_value:.2e}" if p_value < 0.001 else f"p={p_value:.3f}"
        title_text = f"{title or metric}\n(Kruskal-Wallis {p_str})"
    else:
        title_text = title or metric

    ax.set_title(title_text)
    ax.set_ylabel(ylabel or metric)
    ax.set_xlabel('Cluster')
    ax.tick_params(axis='x', rotation=45)


def plot_metric_by_enrichment(df, metric, ax, title=None, ylabel=None):
    """Create a box plot of a metric by enrichment category with Kruskal-Wallis test."""
    if metric not in df.columns or 'enrichment' not in df.columns:
        ax.text(0.5, 0.5, f"Required columns not found", ha='center', va='center')
        ax.set_title(title or metric)
        return

    sns.boxplot(data=df, x='enrichment', y=metric, ax=ax, palette='Set2', showfliers=False)
    sns.stripplot(data=df, x='enrichment', y=metric, ax=ax, color='black', alpha=0.3, size=2)

    # Kruskal-Wallis test
    groups = [group[metric].dropna().values for name, group in df.groupby('enrichment')]
    groups = [g for g in groups if len(g) > 0]
    if len(groups) >= 2:
        stat, p_value = stats.kruskal(*groups)
        p_str = f"p={p_value:.2e}" if p_value < 0.001 else f"p={p_value:.3f}"
        title_text = f"{title or metric}\n(Kruskal-Wallis {p_str})"
    else:
        title_text = title or metric

    ax.set_title(title_text)
    ax.set_ylabel(ylabel or metric)
    ax.set_xlabel('Enrichment')
    ax.tick_params(axis='x', rotation=45)


def plot_metric_by_sample(df, metric, ax, title=None, ylabel=None):
    """Create a box plot of a metric by sample with Kruskal-Wallis test."""
    if metric not in df.columns or 'sample' not in df.columns:
        ax.text(0.5, 0.5, f"Required columns not found", ha='center', va='center')
        ax.set_title(title or metric)
        return

    # Sort samples by median value
    sample_order = df.groupby('sample')[metric].median().sort_values(ascending=False).index

    sns.boxplot(data=df, x='sample', y=metric, order=sample_order, ax=ax,
                palette='tab10', showfliers=False)
    sns.stripplot(data=df, x='sample', y=metric, order=sample_order, ax=ax,
                  color='black', alpha=0.3, size=2)

    # Kruskal-Wallis test
    groups = [group[metric].dropna().values for name, group in df.groupby('sample')]
    groups = [g for g in groups if len(g) > 0]
    if len(groups) >= 2:
        stat, p_value = stats.kruskal(*groups)
        p_str = f"p={p_value:.2e}" if p_value < 0.001 else f"p={p_value:.3f}"
        title_text = f"{title or metric}\n(Kruskal-Wallis {p_str})"
    else:
        title_text = title or metric

    ax.set_title(title_text)
    ax.set_ylabel(ylabel or metric)
    ax.set_xlabel('Sample')
    ax.tick_params(axis='x', rotation=45)


def plot_metric_by_cluster_and_sample(df, metric, ax, title=None, ylabel=None):
    """Create a box plot of a metric by cluster, colored by sample."""
    if metric not in df.columns or 'sample' not in df.columns:
        ax.text(0.5, 0.5, f"Required columns not found", ha='center', va='center')
        ax.set_title(title or metric)
        return

    # Sort clusters by median value
    cluster_order = df.groupby('cluster')[metric].median().sort_values(ascending=False).index

    sns.boxplot(data=df, x='cluster', y=metric, hue='sample', order=cluster_order, ax=ax,
                palette='tab10', showfliers=False)

    # Kruskal-Wallis test (across clusters, ignoring sample)
    groups = [group[metric].dropna().values for name, group in df.groupby('cluster')]
    groups = [g for g in groups if len(g) > 0]
    if len(groups) >= 2:
        stat, p_value = stats.kruskal(*groups)
        p_str = f"p={p_value:.2e}" if p_value < 0.001 else f"p={p_value:.3f}"
        title_text = f"{title or metric}\n(Kruskal-Wallis {p_str})"
    else:
        title_text = title or metric

    ax.set_title(title_text)
    ax.set_ylabel(ylabel or metric)
    ax.set_xlabel('Cluster')
    ax.legend(title='Sample', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize='small')
    ax.tick_params(axis='x', rotation=45)


def plot_cluster_composition(df, group_col, ax, title):
    """Create a stacked bar chart of cluster composition."""
    if group_col not in df.columns:
        ax.text(0.5, 0.5, f"Column '{group_col}' not found", ha='center', va='center')
        ax.set_title(title)
        return

    # Calculate composition
    composition = df.groupby(['cluster', group_col]).size().unstack(fill_value=0)
    composition_pct = composition.div(composition.sum(axis=1), axis=0) * 100

    # Sort clusters by size
    cluster_order = df['cluster'].value_counts().index

    composition_pct = composition_pct.reindex(cluster_order)
    composition_pct.plot(kind='bar', stacked=True, ax=ax, colormap='tab10')

    ax.set_title(title)
    ax.set_ylabel('Percentage')
    ax.set_xlabel('Cluster')
    ax.legend(title=group_col, bbox_to_anchor=(1.02, 1), loc='upper left')
    ax.tick_params(axis='x', rotation=45)


def plot_cluster_sizes(df, ax):
    """Create a bar chart of cluster sizes."""
    cluster_sizes = df['cluster'].value_counts().sort_index()

    colors = []
    if 'enrichment' in df.columns:
        enrichment_map = df.groupby('cluster')['enrichment'].first()
        color_map = {'mixed': '#999999'}
        palette = sns.color_palette('Set2', n_colors=10)
        enrichments = df['enrichment'].unique()
        for i, e in enumerate(enrichments):
            if e != 'mixed' and e not in color_map:
                color_map[e] = palette[i % len(palette)]
        colors = [color_map.get(enrichment_map.get(c, 'mixed'), '#999999') for c in cluster_sizes.index]
    else:
        colors = 'steelblue'

    cluster_sizes.plot(kind='bar', ax=ax, color=colors)
    ax.set_title('Cluster Sizes')
    ax.set_ylabel('Number of sequences')
    ax.set_xlabel('Cluster')
    ax.tick_params(axis='x', rotation=45)


def plot_correlation_heatmap(df, metrics, ax):
    """Create a correlation heatmap of numeric metrics."""
    available_metrics = [m for m in metrics if m in df.columns]
    if len(available_metrics) < 2:
        ax.text(0.5, 0.5, "Not enough metrics for correlation", ha='center', va='center')
        return

    corr = df[available_metrics].corr()
    sns.heatmap(corr, annot=True, fmt='.2f', cmap='coolwarm', center=0, ax=ax,
                square=True, linewidths=0.5)
    ax.set_title('Metric Correlations')


def compute_cluster_summary(df):
    """Compute summary statistics per cluster."""
    numeric_cols = ['read_length', 'centroid_distance', 'mapq', 'de', 'align_len', 'align_fraction',
                    'primary_mapq', 'primary_de', 'primary_align_len', 'primary_align_fraction',
                    'total_align_len', 'total_align_fraction', 'n_alignments', 'n_secondary', 'n_supplementary']
    available_cols = [c for c in numeric_cols if c in df.columns]

    summary_rows = []
    for cluster in sorted(df['cluster'].unique()):
        cluster_df = df[df['cluster'] == cluster]
        row = {
            'cluster': cluster,
            'n_sequences': len(cluster_df),
        }

        # Add group/sample composition
        if 'group' in df.columns:
            for group in df['group'].unique():
                row[f'{group}_count'] = (cluster_df['group'] == group).sum()
                row[f'{group}_pct'] = (cluster_df['group'] == group).mean() * 100

        # Add enrichment if available
        if 'enrichment' in df.columns:
            row['enrichment'] = cluster_df['enrichment'].iloc[0] if len(cluster_df) > 0 else 'unknown'

        # Add numeric summaries
        for col in available_cols:
            row[f'{col}_mean'] = cluster_df[col].mean()
            row[f'{col}_median'] = cluster_df[col].median()
            row[f'{col}_std'] = cluster_df[col].std()

        summary_rows.append(row)

    return pd.DataFrame(summary_rows)


def main():
    parser = argparse.ArgumentParser(
        description="Generate diagnostic plots for cluster analysis",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--annotated", required=True,
                        help="Annotated sequences TSV file (from KaryoScope_annotate_sequences.py)")
    parser.add_argument("--cluster-analysis", dest="cluster_analysis", default=None,
                        help="Cluster analysis TSV file (optional, for enrichment info)")
    parser.add_argument("--output-prefix", dest="output_prefix", required=True,
                        help="Output prefix for generated files")
    parser.add_argument("--dark-mode", dest="dark_mode", action="store_true", default=False,
                        help="Output plots with dark background (default: False)")

    args = parser.parse_args()

    # Set up plot style
    if args.dark_mode:
        plt.style.use('dark_background')
        sns.set_style("darkgrid")
    else:
        sns.set_style("whitegrid")

    print("=" * 60)
    print("KaryoScope Cluster Diagnostics")
    print("=" * 60)

    # Load data
    df = load_annotated_data(args.annotated)

    # Optionally load enrichment info
    if args.cluster_analysis and os.path.exists(args.cluster_analysis):
        print(f"\nLoading cluster analysis: {args.cluster_analysis}")
        ca = pd.read_csv(args.cluster_analysis, sep='\t')
        enrichment_map = dict(zip(ca['cluster_id'], ca['enrichment']))
        df['enrichment'] = df['cluster'].map(enrichment_map).fillna('unknown')
        print(f"  Added enrichment labels")

    # Infer cluster analysis file from annotated file path if not provided
    if 'enrichment' not in df.columns:
        # Try to find cluster_analysis.tsv based on naming convention
        annotated_path = args.annotated
        if '.read_assignments.annotated.tsv' in annotated_path:
            ca_path = annotated_path.replace('.read_assignments.annotated.tsv', '.cluster_analysis.tsv')
            if os.path.exists(ca_path):
                print(f"\nAuto-discovered cluster analysis: {ca_path}")
                ca = pd.read_csv(ca_path, sep='\t')
                enrichment_map = dict(zip(ca['cluster_id'], ca['enrichment']))
                df['enrichment'] = df['cluster'].map(enrichment_map).fillna('unknown')
                print(f"  Added enrichment labels")

    # Set up plotting style
    plt.rcParams['figure.dpi'] = 150

    # === Page 1: Metric distributions by cluster ===
    print("\nGenerating metric plots...")
    metrics_pdf = f"{args.output_prefix}.cluster_metrics.pdf"

    with PdfPages(metrics_pdf) as pdf:
        # Page 1: Key metrics by cluster (3x2 grid)
        fig, axes = plt.subplots(3, 2, figsize=(14, 14))
        fig.suptitle('Sequence Metrics by Cluster', fontsize=14, fontweight='bold')

        plot_metric_by_cluster(df, 'read_length', axes[0, 0],
                               title='Read Length by Cluster', ylabel='Read Length (bp)')
        plot_metric_by_cluster(df, 'align_fraction', axes[0, 1],
                               title='Alignment Fraction by Cluster', ylabel='Alignment Fraction')
        plot_metric_by_cluster(df, 'de', axes[1, 0],
                               title='Divergence (Error Rate) by Cluster', ylabel='Divergence')
        plot_metric_by_cluster(df, 'centroid_distance', axes[1, 1],
                               title='Centroid Distance by Cluster', ylabel='Centroid Distance')
        plot_metric_by_cluster(df, 'mapq', axes[2, 0],
                               title='Mapping Quality by Cluster', ylabel='MAPQ')
        plot_metric_by_cluster(df, 'align_len', axes[2, 1],
                               title='Alignment Length by Cluster', ylabel='Alignment Length (bp)')

        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # Page 2: Metrics by cluster, colored by sample (3x2 grid)
        if 'sample' in df.columns and df['sample'].nunique() > 1:
            fig, axes = plt.subplots(3, 2, figsize=(16, 14))
            fig.suptitle('Sequence Metrics by Cluster (Colored by Sample)', fontsize=14, fontweight='bold')

            plot_metric_by_cluster_and_sample(df, 'read_length', axes[0, 0],
                                              title='Read Length by Cluster', ylabel='Read Length (bp)')
            plot_metric_by_cluster_and_sample(df, 'align_fraction', axes[0, 1],
                                              title='Alignment Fraction by Cluster', ylabel='Alignment Fraction')
            plot_metric_by_cluster_and_sample(df, 'de', axes[1, 0],
                                              title='Divergence by Cluster', ylabel='Divergence')
            plot_metric_by_cluster_and_sample(df, 'centroid_distance', axes[1, 1],
                                              title='Centroid Distance by Cluster', ylabel='Centroid Distance')
            plot_metric_by_cluster_and_sample(df, 'mapq', axes[2, 0],
                                              title='Mapping Quality by Cluster', ylabel='MAPQ')
            plot_metric_by_cluster_and_sample(df, 'align_len', axes[2, 1],
                                              title='Alignment Length by Cluster', ylabel='Alignment Length (bp)')

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

        # Page 3: Metrics by enrichment (if available, 3x2 grid)
        if 'enrichment' in df.columns:
            fig, axes = plt.subplots(3, 2, figsize=(12, 14))
            fig.suptitle('Sequence Metrics by Enrichment Category', fontsize=14, fontweight='bold')

            plot_metric_by_enrichment(df, 'read_length', axes[0, 0],
                                      title='Read Length by Enrichment', ylabel='Read Length (bp)')
            plot_metric_by_enrichment(df, 'align_fraction', axes[0, 1],
                                      title='Alignment Fraction by Enrichment', ylabel='Alignment Fraction')
            plot_metric_by_enrichment(df, 'de', axes[1, 0],
                                      title='Divergence by Enrichment', ylabel='Divergence')
            plot_metric_by_enrichment(df, 'centroid_distance', axes[1, 1],
                                      title='Centroid Distance by Enrichment', ylabel='Centroid Distance')
            plot_metric_by_enrichment(df, 'mapq', axes[2, 0],
                                      title='Mapping Quality by Enrichment', ylabel='MAPQ')
            plot_metric_by_enrichment(df, 'align_len', axes[2, 1],
                                      title='Alignment Length by Enrichment', ylabel='Alignment Length (bp)')

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

        # Page 4: Metrics by sample (batch effect detection, 3x2 grid)
        if 'sample' in df.columns and df['sample'].nunique() > 1:
            fig, axes = plt.subplots(3, 2, figsize=(12, 14))
            fig.suptitle('Sequence Metrics by Sample (Batch Effect Check)', fontsize=14, fontweight='bold')

            plot_metric_by_sample(df, 'read_length', axes[0, 0],
                                  title='Read Length by Sample', ylabel='Read Length (bp)')
            plot_metric_by_sample(df, 'align_fraction', axes[0, 1],
                                  title='Alignment Fraction by Sample', ylabel='Alignment Fraction')
            plot_metric_by_sample(df, 'de', axes[1, 0],
                                  title='Divergence by Sample', ylabel='Divergence')
            plot_metric_by_sample(df, 'centroid_distance', axes[1, 1],
                                  title='Centroid Distance by Sample', ylabel='Centroid Distance')
            plot_metric_by_sample(df, 'mapq', axes[2, 0],
                                  title='Mapping Quality by Sample', ylabel='MAPQ')
            plot_metric_by_sample(df, 'align_len', axes[2, 1],
                                  title='Alignment Length by Sample', ylabel='Alignment Length (bp)')

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

        # Page 5: Alignment statistics by cluster (if available)
        if 'n_alignments' in df.columns:
            fig, axes = plt.subplots(2, 3, figsize=(16, 10))
            fig.suptitle('Alignment Statistics by Cluster', fontsize=14, fontweight='bold')

            plot_metric_by_cluster(df, 'n_alignments', axes[0, 0],
                                   title='Number of Alignments', ylabel='Alignments per Read')
            plot_metric_by_cluster(df, 'n_secondary', axes[0, 1],
                                   title='Secondary Alignments', ylabel='Secondary Alignments')
            plot_metric_by_cluster(df, 'n_supplementary', axes[0, 2],
                                   title='Supplementary Alignments', ylabel='Supplementary Alignments')
            plot_metric_by_cluster(df, 'primary_align_fraction', axes[1, 0],
                                   title='Primary Alignment Fraction', ylabel='Fraction')
            plot_metric_by_cluster(df, 'total_align_fraction', axes[1, 1],
                                   title='Total Alignment Fraction', ylabel='Fraction')
            # Difference between total and primary (contribution of supplementary)
            if 'total_align_fraction' in df.columns and 'primary_align_fraction' in df.columns:
                df['supplementary_contribution'] = df['total_align_fraction'] - df['primary_align_fraction']
                plot_metric_by_cluster(df, 'supplementary_contribution', axes[1, 2],
                                       title='Supplementary Contribution', ylabel='Additional Fraction')

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

        # Page 6: Alignment statistics by sample (batch effects)
        if 'n_alignments' in df.columns and 'sample' in df.columns and df['sample'].nunique() > 1:
            fig, axes = plt.subplots(2, 3, figsize=(16, 10))
            fig.suptitle('Alignment Statistics by Sample (Batch Effect Check)', fontsize=14, fontweight='bold')

            plot_metric_by_sample(df, 'n_alignments', axes[0, 0],
                                  title='Number of Alignments', ylabel='Alignments per Read')
            plot_metric_by_sample(df, 'n_secondary', axes[0, 1],
                                  title='Secondary Alignments', ylabel='Secondary Alignments')
            plot_metric_by_sample(df, 'n_supplementary', axes[0, 2],
                                  title='Supplementary Alignments', ylabel='Supplementary Alignments')
            plot_metric_by_sample(df, 'primary_align_fraction', axes[1, 0],
                                  title='Primary Alignment Fraction', ylabel='Fraction')
            plot_metric_by_sample(df, 'total_align_fraction', axes[1, 1],
                                  title='Total Alignment Fraction', ylabel='Fraction')
            if 'supplementary_contribution' in df.columns:
                plot_metric_by_sample(df, 'supplementary_contribution', axes[1, 2],
                                      title='Supplementary Contribution', ylabel='Additional Fraction')

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

        # Page 7: Correlation heatmap and cluster sizes
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        metrics = ['read_length', 'centroid_distance', 'primary_mapq', 'primary_de',
                   'primary_align_fraction', 'total_align_fraction', 'n_alignments',
                   'n_secondary', 'n_supplementary']
        plot_correlation_heatmap(df, metrics, axes[0])
        plot_cluster_sizes(df, axes[1])

        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    print(f"  Saved: {metrics_pdf}")

    # === Page 2: Composition plots ===
    print("\nGenerating composition plots...")
    composition_pdf = f"{args.output_prefix}.cluster_composition.pdf"

    with PdfPages(composition_pdf) as pdf:
        # Group composition
        if 'group' in df.columns:
            fig, ax = plt.subplots(figsize=(12, 6))
            plot_cluster_composition(df, 'group', ax, 'Group Composition by Cluster')
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

        # Sample composition
        if 'sample' in df.columns:
            fig, ax = plt.subplots(figsize=(14, 6))
            plot_cluster_composition(df, 'sample', ax, 'Sample Composition by Cluster')
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

        # Sequencing approach composition
        if 'sequencing_approach' in df.columns:
            fig, ax = plt.subplots(figsize=(12, 6))
            plot_cluster_composition(df, 'sequencing_approach', ax, 'Sequencing Approach by Cluster')
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

    print(f"  Saved: {composition_pdf}")

    # === Summary table ===
    print("\nGenerating summary table...")
    summary = compute_cluster_summary(df)
    summary_file = f"{args.output_prefix}.cluster_summary.tsv"
    summary.to_csv(summary_file, sep='\t', index=False)
    print(f"  Saved: {summary_file}")

    print(f"\n{'=' * 60}")
    print("Summary")
    print("=" * 60)
    print(f"Total sequences: {len(df)}")
    print(f"Clusters: {df['cluster'].nunique()}")
    if 'enrichment' in df.columns:
        print(f"Enrichment categories: {df['enrichment'].value_counts().to_dict()}")
    print(f"\nOutput files:")
    print(f"  - {metrics_pdf}")
    print(f"  - {composition_pdf}")
    print(f"  - {summary_file}")


if __name__ == "__main__":
    main()
