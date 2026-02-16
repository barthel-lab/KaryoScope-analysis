#!/usr/bin/env python3
"""
KaryoScope Cluster Diagnostics

Creates diagnostic plots comparing clusters across various metrics from
annotated sequence data. Supports both exploratory analysis and publication-
quality figure generation.

Usage:
  # Exploratory diagnostics (default)
  python KaryoScope_cluster_diagnostics.py \
    --annotated analysis.read_assignments.annotated.tsv \
    --output-prefix analysis_diagnostics

  # Publication-quality figures for specific clusters
  python KaryoScope_cluster_diagnostics.py \
    --annotated analysis.annotated.tsv \
    --cluster-analysis analysis.cluster_analysis.tsv \
    --output-prefix figures \
    --pub-comparison --clusters 82,88 --cluster-labels "82:α-sat,88:rDNA"

  # Summary heatmap of all clusters
  python KaryoScope_cluster_diagnostics.py \
    --annotated analysis.annotated.tsv \
    --cluster-analysis analysis.cluster_analysis.tsv \
    --output-prefix figures \
    --pub-heatmap

Generates:
  Exploratory (default):
    - {prefix}.cluster_metrics.pdf: Box/violin plots of metrics by cluster
    - {prefix}.cluster_composition.pdf: Sample/group composition per cluster
    - {prefix}.cluster_summary.tsv: Summary statistics per cluster

  Publication (--pub-*):
    - {prefix}.comparison.pdf: Focused comparison with statistical annotations
    - {prefix}.counts.pdf: Count metric distributions
    - {prefix}.heatmap.pdf: Z-scored summary heatmap
    - {prefix}.statistics.tsv: Detailed statistical comparisons
"""

import argparse
import os
import sys
from typing import List, Dict, Tuple, Optional

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats

# Global plot style variables (set by main based on --dark-mode)
STRIP_COLOR = 'black'
STRIP_ALPHA = 0.3

# =============================================================================
# Publication-quality style settings
# =============================================================================

PUBLICATION_RCPARAMS = {
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9,
    'figure.titlesize': 14,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
    'axes.linewidth': 1.0,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'lines.linewidth': 1.5,
    'patch.linewidth': 0.5,
    'pdf.fonttype': 42,      # TrueType fonts (editable text)
    'svg.fonttype': 'none',  # Keep text as text elements
}

ENRICHMENT_COLORS = {
    'E6E7': '#c0392b',       # Dark red
    'primary': '#2980b9',    # Dark blue
    'other': '#7f8c8d',      # Gray
    'mixed': '#95a5a6',      # Light gray
}

ENRICHMENT_COLORS_DARK = {
    'E6E7': '#e74c3c',       # Bright red
    'primary': '#3498db',    # Bright blue
    'other': '#95a5a6',      # Gray
    'mixed': '#7f8c8d',      # Dark gray
}


def format_pvalue(p: float) -> str:
    """Format p-value for display."""
    if p < 0.0001:
        return "p < 0.0001"
    elif p < 0.001:
        return f"p = {p:.4f}"
    elif p < 0.01:
        return f"p = {p:.3f}"
    else:
        return f"p = {p:.2f}"


def format_pvalue_stars(p: float) -> str:
    """Format p-value as significance stars."""
    if p < 0.0001:
        return "****"
    elif p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"


def compute_effect_size(group1: np.ndarray, group2: np.ndarray) -> Tuple[float, str]:
    """
    Compute rank-biserial correlation (effect size for Mann-Whitney U).
    Returns (effect_size, interpretation).
    """
    if len(group1) == 0 or len(group2) == 0:
        return 0.0, "n/a"
    try:
        stat, _ = stats.mannwhitneyu(group1, group2, alternative='two-sided')
    except ValueError:
        # Can happen with identical values or other edge cases
        return 0.0, "n/a"
    n1, n2 = len(group1), len(group2)
    r = 1 - (2 * stat) / (n1 * n2)

    abs_r = abs(r)
    if abs_r < 0.1:
        interp = "negligible"
    elif abs_r < 0.3:
        interp = "small"
    elif abs_r < 0.5:
        interp = "medium"
    else:
        interp = "large"

    return r, interp


def format_metric_label(metric: str) -> str:
    """Format metric name for y-axis label."""
    labels = {
        'read_length': 'Read Length (bp)',
        'primary_mapq': 'Mapping Quality',
        'primary_de': 'Divergence Rate',
        'primary_align_fraction': 'Primary Alignment\nFraction',
        'total_align_fraction': 'Total Alignment\nFraction',
        'n_alignments': 'Alignments (n)',
        'n_secondary': 'Secondary\nAlignments (n)',
        'n_supplementary': 'Supplementary\nAlignments (n)',
        'centroid_distance': 'Centroid Distance',
    }
    return labels.get(metric, metric.replace('_', ' ').title())


def format_metric_title(metric: str) -> str:
    """Format metric name for panel title."""
    titles = {
        'read_length': 'Read Length',
        'primary_mapq': 'Mapping Quality (MAPQ)',
        'primary_de': 'Sequence Divergence',
        'primary_align_fraction': 'Primary Alignment Fraction',
        'total_align_fraction': 'Total Alignment Fraction',
        'n_alignments': 'Total Alignments',
        'n_secondary': 'Secondary Alignments',
        'n_supplementary': 'Supplementary Alignments',
        'centroid_distance': 'Cluster Centroid Distance',
    }
    return titles.get(metric, metric.replace('_', ' ').title())


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
                  color=STRIP_COLOR, alpha=STRIP_ALPHA, size=2)

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
    sns.stripplot(data=df, x='enrichment', y=metric, ax=ax, color=STRIP_COLOR, alpha=STRIP_ALPHA, size=2)

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
                  color=STRIP_COLOR, alpha=STRIP_ALPHA, size=2)

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


def plot_significant_vs_nonsig(df, significant_clusters, nonsig_clusters, metrics,
                                enrichment_map=None, q_value_map=None, dark_mode=False):
    """
    Create a multi-box comparison plot: each significant cluster + pooled non-sig.
    Uses boxplots for continuous metrics, stacked bars for count metrics (n_*).
    Returns a matplotlib figure.
    """
    # Separate count metrics from continuous metrics
    count_metrics = [m for m in metrics if m.startswith('n_') and m in df.columns]
    continuous_metrics = [m for m in metrics if not m.startswith('n_') and m in df.columns]

    all_metrics = continuous_metrics + count_metrics
    n_metrics = len(all_metrics)

    if n_metrics == 0:
        return None

    # Create figure - wider to accommodate many clusters
    n_cols = 2
    n_rows = (n_metrics + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 5 * n_rows))
    axes = axes.flatten() if n_metrics > 1 else [axes]

    # Prepare data
    sig_df = df[df['cluster'].isin(significant_clusters)]
    nonsig_df = df[df['cluster'].isin(nonsig_clusters)]
    n_sig = len(significant_clusters)

    # Colors based on enrichment direction
    # E6E7-enriched: red, primary-enriched: blue, non-sig: gray
    labels = []
    colors = []
    for cid in significant_clusters:
        labels.append(str(cid))
        if enrichment_map:
            enr = enrichment_map.get(cid, '')
            if 'E6E7' in enr:
                colors.append('#e74c3c')  # red
            elif 'primary' in enr:
                colors.append('#3498db')  # blue
            else:
                colors.append('#9b59b6')  # purple for unknown
        else:
            colors.append('#e74c3c')  # default red

    labels.append('NS')  # Non-significant
    colors.append('#999999')  # gray

    fig.suptitle(f"Significant Clusters (n={n_sig}) vs Non-significant (pooled, n={len(nonsig_clusters)})",
                 fontsize=14, fontweight='bold')

    for i, metric in enumerate(all_metrics):
        ax = axes[i]
        is_count_metric = metric.startswith('n_')

        if is_count_metric:
            # Stacked bar plot for count data - show proportions of count values
            # Determine count value categories (cap at reasonable max for visualization)
            all_vals = df[metric].dropna()
            max_val = min(int(all_vals.max()), 10)  # Cap at 10 for readability
            count_categories = list(range(max_val + 1)) + [f'{max_val + 1}+'] if max_val < all_vals.max() else list(range(max_val + 1))
            n_cats = len(count_categories)

            # Create a colormap for count categories
            count_cmap = plt.cm.viridis(np.linspace(0.2, 0.9, n_cats))

            x = np.arange(len(labels))
            bar_width = 0.8

            # Compute proportions for each cluster
            bottom = np.zeros(len(labels))

            for cat_idx, cat_val in enumerate(count_categories):
                proportions = []
                for cid in significant_clusters:
                    cluster_vals = df[df['cluster'] == cid][metric].dropna()
                    if len(cluster_vals) > 0:
                        if isinstance(cat_val, str):  # "N+" category
                            prop = (cluster_vals > max_val).mean()
                        else:
                            prop = (cluster_vals == cat_val).mean()
                        proportions.append(prop)
                    else:
                        proportions.append(0)

                # Non-significant pooled
                if len(nonsig_df) > 0:
                    nonsig_vals = nonsig_df[metric].dropna()
                    if isinstance(cat_val, str):
                        prop = (nonsig_vals > max_val).mean()
                    else:
                        prop = (nonsig_vals == cat_val).mean()
                    proportions.append(prop)
                else:
                    proportions.append(0)

                proportions = np.array(proportions)
                ax.bar(x, proportions, bar_width, bottom=bottom, color=count_cmap[cat_idx],
                       label=str(cat_val), edgecolor='white', linewidth=0.3)
                bottom += proportions

            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=90, fontsize=6, ha='center')
            ax.set_title(f"{metric} (proportion by count)")
            ax.set_ylabel('Proportion')
            ax.set_ylim(0, 1)

            # Add legend for count categories (only on first count metric)
            if metric == count_metrics[0]:
                ax.legend(title='Count', loc='upper right', fontsize=5, title_fontsize=6,
                         ncol=min(4, n_cats), framealpha=0.7)

        else:
            # Box plot for continuous data
            box_data = []
            valid_labels = []
            valid_colors = []

            for j, cid in enumerate(significant_clusters):
                cluster_vals = df[df['cluster'] == cid][metric].dropna()
                if len(cluster_vals) > 0:
                    box_data.append(cluster_vals.values)
                    valid_labels.append(labels[j])
                    valid_colors.append(colors[j])

            nonsig_vals = nonsig_df[metric].dropna()
            if len(nonsig_vals) > 0:
                box_data.append(nonsig_vals.values)
                valid_labels.append('NS')
                valid_colors.append('#999999')

            if len(box_data) == 0:
                ax.text(0.5, 0.5, f"No data for {metric}", ha='center', va='center')
                continue

            bp = ax.boxplot(box_data, patch_artist=True, showfliers=False)

            for patch, color in zip(bp['boxes'], valid_colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)

            if dark_mode:
                for element in ['whiskers', 'caps', 'medians']:
                    for item in bp[element]:
                        item.set_color('white')

            ax.set_xticklabels(valid_labels, rotation=90, fontsize=6, ha='center')
            ax.set_title(metric)
            ax.set_ylabel(metric)

    # Hide unused axes
    for i in range(n_metrics, len(axes)):
        axes[i].set_visible(False)

    plt.tight_layout()
    return fig


def plot_cluster_vs_others(df, cluster_ids, metrics, group_label=None, enrichment_label=None,
                           dark_mode=False, background_clusters=None):
    """
    Create a comparison plot of one or more clusters vs others.
    cluster_ids can be a single int or a list of ints.
    background_clusters: if provided, compare against only these clusters (not all others)
    Returns a matplotlib figure.
    """
    # Handle single cluster or list
    if isinstance(cluster_ids, int):
        cluster_ids = [cluster_ids]

    # Split data
    cluster_df = df[df['cluster'].isin(cluster_ids)].copy()

    if background_clusters is not None:
        others_df = df[df['cluster'].isin(background_clusters)].copy()
        others_label = 'Non-sig'
    else:
        others_df = df[~df['cluster'].isin(cluster_ids)].copy()
        others_label = 'Others'

    if group_label is None:
        if len(cluster_ids) == 1:
            group_label = f'Cluster {cluster_ids[0]}'
        else:
            group_label = f'Clusters ({len(cluster_ids)})'

    cluster_df['group_label'] = group_label
    others_df['group_label'] = 'All Others'

    combined = pd.concat([cluster_df, others_df])

    # Filter to available metrics
    available_metrics = [m for m in metrics if m in df.columns]
    n_metrics = len(available_metrics)

    if n_metrics == 0:
        return None

    # Create figure
    n_cols = 3
    n_rows = (n_metrics + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, 4 * n_rows))
    axes = axes.flatten() if n_metrics > 1 else [axes]

    # Title
    if len(cluster_ids) == 1:
        title = f"Cluster {cluster_ids[0]} (n={len(cluster_df)}) vs {others_label} (n={len(others_df)})"
        box_label = f'C{cluster_ids[0]}'
    else:
        title = f"{group_label} (n={len(cluster_df)}) vs {others_label} (n={len(others_df)})"
        box_label = 'Selected'
    if enrichment_label:
        title += f" [{enrichment_label}]"
    fig.suptitle(title, fontsize=14, fontweight='bold')

    # Box colors
    box_colors = ['#e74c3c', '#3498db']  # red for cluster, blue for others

    for i, metric in enumerate(available_metrics):
        ax = axes[i]

        # Box plot
        bp = ax.boxplot(
            [cluster_df[metric].dropna(), others_df[metric].dropna()],
            labels=[box_label, others_label],
            patch_artist=True,
            showfliers=False
        )

        # Color boxes
        for patch, color in zip(bp['boxes'], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        # Style whiskers/caps/medians for dark mode
        if dark_mode:
            for element in ['whiskers', 'caps', 'medians']:
                for item in bp[element]:
                    item.set_color('white')

        # Mann-Whitney U test
        cluster_vals = cluster_df[metric].dropna()
        others_vals = others_df[metric].dropna()
        if len(cluster_vals) > 0 and len(others_vals) > 0:
            try:
                stat, p_val = stats.mannwhitneyu(cluster_vals, others_vals, alternative='two-sided')
                p_str = f"p={p_val:.2e}" if p_val < 0.001 else f"p={p_val:.3f}"

                # Effect size (rank-biserial correlation)
                n1, n2 = len(cluster_vals), len(others_vals)
                effect_size = 1 - (2 * stat) / (n1 * n2)

                ax.set_title(f"{metric}\n({p_str}, r={effect_size:.2f})")
            except (ValueError, ZeroDivisionError):
                ax.set_title(metric)
        else:
            ax.set_title(metric)

        ax.set_ylabel(metric)

    # Hide unused axes
    for i in range(n_metrics, len(axes)):
        axes[i].set_visible(False)

    plt.tight_layout()
    return fig


# =============================================================================
# Publication-quality figure functions
# =============================================================================

def pub_comparison_figure(
    df: pd.DataFrame,
    target_clusters: List[int],
    metrics: List[str],
    enrichment_map: Dict[int, str] = None,
    cluster_labels: Dict[int, str] = None,
    dark_mode: bool = False,
    show_points: bool = True,
    figsize: Tuple[float, float] = None,
) -> plt.Figure:
    """
    Create a publication-quality comparison figure for specific clusters vs pooled others.

    Features:
    - Statistical annotations with significance brackets
    - Effect sizes displayed
    - Clean typography (12pt+ fonts)
    - Proper axis labels and titles
    - Sample sizes in labels
    """
    # Apply publication style
    plt.rcParams.update(PUBLICATION_RCPARAMS)

    if dark_mode:
        plt.rcParams.update({
            'axes.facecolor': '#1a1a1a',
            'figure.facecolor': '#1a1a1a',
            'text.color': 'white',
            'axes.labelcolor': 'white',
            'xtick.color': 'white',
            'ytick.color': 'white',
            'axes.edgecolor': 'white',
        })

    metrics = [m for m in metrics if m in df.columns]
    n_metrics = len(metrics)

    if n_metrics == 0:
        raise ValueError("No valid metrics found")

    # Layout
    if n_metrics <= 3:
        n_cols, n_rows = n_metrics, 1
    elif n_metrics <= 6:
        n_cols, n_rows = 3, 2
    else:
        n_cols = 4
        n_rows = (n_metrics + 3) // 4

    if figsize is None:
        figsize = (4.5 * n_cols, 4.5 * n_rows)

    # Data prep
    target_df = df[df['cluster'].isin(target_clusters)]
    other_df = df[~df['cluster'].isin(target_clusters)]

    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    if n_metrics == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    colors = ENRICHMENT_COLORS_DARK if dark_mode else ENRICHMENT_COLORS

    # Determine colors per cluster
    cluster_colors = []
    for cid in target_clusters:
        if enrichment_map:
            enr = enrichment_map.get(cid, '')
            if 'E6E7' in enr:
                cluster_colors.append(colors['E6E7'])
            elif 'primary' in enr:
                cluster_colors.append(colors['primary'])
            else:
                cluster_colors.append(colors['mixed'])
        else:
            cluster_colors.append(colors['E6E7'])
    cluster_colors.append(colors['other'])

    # Labels with sample sizes
    labels = []
    for cid in target_clusters:
        n = len(df[df['cluster'] == cid])
        if cluster_labels and cid in cluster_labels:
            labels.append(f"{cluster_labels[cid]}\n(n={n})")
        else:
            labels.append(f"Cluster {cid}\n(n={n})")
    labels.append(f"Others\n(n={len(other_df)})")

    # Plot each metric
    for i, metric in enumerate(metrics):
        ax = axes[i]

        # Build box data
        box_data = []
        for cid in target_clusters:
            vals = df[df['cluster'] == cid][metric].dropna().values
            box_data.append(vals)
        box_data.append(other_df[metric].dropna().values)

        # Create boxplot
        bp = ax.boxplot(box_data, patch_artist=True, widths=0.6, showfliers=False,
                       medianprops={'color': 'white' if dark_mode else 'black', 'linewidth': 2})

        for patch, color in zip(bp['boxes'], cluster_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
            patch.set_edgecolor('white' if dark_mode else 'black')
            patch.set_linewidth(1)

        whisker_color = 'white' if dark_mode else 'black'
        for element in ['whiskers', 'caps']:
            for item in bp[element]:
                item.set_color(whisker_color)
                item.set_linewidth(1)

        # Overlay points (behind boxes, smaller and more transparent)
        if show_points:
            for j, data in enumerate(box_data):
                if len(data) > 0:
                    jitter = np.random.uniform(-0.15, 0.15, size=len(data))
                    point_color = 'white' if dark_mode else 'black'
                    ax.scatter(np.full_like(data, j + 1, dtype=float) + jitter, data,
                              alpha=0.15, s=8, color=point_color, zorder=1, edgecolors='none')

        # Statistical annotations
        y_max = max([np.percentile(d, 95) if len(d) > 0 else 0 for d in box_data])
        y_min = min([np.percentile(d, 5) if len(d) > 0 else 0 for d in box_data])
        y_range = y_max - y_min if (y_max - y_min) != 0 else 1

        annotation_y = y_max + 0.08 * y_range

        for j, cid in enumerate(target_clusters):
            target_vals = box_data[j]
            other_vals = box_data[-1]

            if len(target_vals) > 0 and len(other_vals) > 0:
                stat, p_val = stats.mannwhitneyu(target_vals, other_vals, alternative='two-sided')
                r, _ = compute_effect_size(target_vals, other_vals)
                stars = format_pvalue_stars(p_val)

                bracket_y = annotation_y + j * 0.14 * y_range

                # Draw bracket
                ax.plot([j + 1, len(box_data)], [bracket_y, bracket_y],
                       color=whisker_color, linewidth=0.8, clip_on=False)
                ax.plot([j + 1, j + 1], [bracket_y - 0.025 * y_range, bracket_y],
                       color=whisker_color, linewidth=0.8, clip_on=False)
                ax.plot([len(box_data), len(box_data)], [bracket_y - 0.025 * y_range, bracket_y],
                       color=whisker_color, linewidth=0.8, clip_on=False)

                # Annotation
                mid_x = (j + 1 + len(box_data)) / 2
                if stars != "ns":
                    ax.text(mid_x, bracket_y + 0.015 * y_range, stars,
                           ha='center', va='bottom', fontsize=11, fontweight='bold')
                else:
                    ax.text(mid_x, bracket_y + 0.015 * y_range, "ns",
                           ha='center', va='bottom', fontsize=9, style='italic', color='gray')

        # Formatting
        ax.set_xticklabels(labels, fontsize=10)
        ax.set_ylabel(format_metric_label(metric), fontsize=11)
        ax.set_title(format_metric_title(metric), fontsize=12, fontweight='bold', pad=10)

        # Extend y-axis for annotations
        ax.set_ylim(y_min - 0.05 * y_range,
                   annotation_y + (len(target_clusters) + 0.8) * 0.14 * y_range)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    for i in range(n_metrics, len(axes)):
        axes[i].set_visible(False)

    # Figure title
    cluster_str = ", ".join([cluster_labels.get(c, str(c)) if cluster_labels else str(c)
                             for c in target_clusters])
    fig.suptitle(f"Comparison: {cluster_str} vs Others", fontsize=14, fontweight='bold', y=1.02)

    # Test description
    fig.text(0.5, -0.02,
             "Mann-Whitney U test with rank-biserial effect size; * p<0.05, ** p<0.01, *** p<0.001, **** p<0.0001",
             ha='center', fontsize=9, style='italic')

    plt.tight_layout()
    return fig


def pub_count_distribution_figure(
    df: pd.DataFrame,
    target_clusters: List[int],
    count_metrics: List[str] = None,
    cluster_labels: Dict[int, str] = None,
    dark_mode: bool = False,
    figsize: Tuple[float, float] = None,
) -> plt.Figure:
    """
    Create publication-quality stacked bar chart for count metrics.
    Shows proportion of reads at each count value with Chi-square test.
    """
    plt.rcParams.update(PUBLICATION_RCPARAMS)

    if dark_mode:
        plt.rcParams.update({
            'axes.facecolor': '#1a1a1a',
            'figure.facecolor': '#1a1a1a',
            'text.color': 'white',
            'axes.labelcolor': 'white',
            'xtick.color': 'white',
            'ytick.color': 'white',
            'axes.edgecolor': 'white',
        })

    if count_metrics is None:
        count_metrics = ['n_alignments', 'n_secondary', 'n_supplementary']
    count_metrics = [m for m in count_metrics if m in df.columns]

    if not count_metrics:
        return None

    n_metrics = len(count_metrics)
    if figsize is None:
        figsize = (4.5 * n_metrics, 5.5)  # Slightly taller for stats annotation

    fig, axes = plt.subplots(1, n_metrics, figsize=figsize)
    if n_metrics == 1:
        axes = [axes]

    other_df = df[~df['cluster'].isin(target_clusters)]

    # Labels with sample sizes
    labels = []
    for cid in target_clusters:
        n = len(df[df['cluster'] == cid])
        if cluster_labels and cid in cluster_labels:
            labels.append(f"{cluster_labels[cid]}\n(n={n})")
        else:
            labels.append(f"C{cid}\n(n={n})")
    labels.append(f"Others\n(n={len(other_df)})")

    legend_handles = None
    legend_labels = None

    for ax_idx, (ax, metric) in enumerate(zip(axes, count_metrics)):
        all_vals = df[metric].dropna()
        max_val = min(int(all_vals.quantile(0.95)), 10)
        categories = list(range(max_val + 1))
        if all_vals.max() > max_val:
            categories.append(f"{max_val + 1}+")

        # Compute counts (not proportions) for chi-square test
        count_data = []
        proportions = []

        for cid in target_clusters:
            cluster_vals = df[df['cluster'] == cid][metric].dropna()
            counts = []
            props = []
            for cat in categories:
                if isinstance(cat, str):
                    c = (cluster_vals > max_val).sum()
                    p = (cluster_vals > max_val).mean() if len(cluster_vals) > 0 else 0
                else:
                    c = (cluster_vals == cat).sum()
                    p = (cluster_vals == cat).mean() if len(cluster_vals) > 0 else 0
                counts.append(c)
                props.append(p)
            count_data.append(counts)
            proportions.append(props)

        other_vals = other_df[metric].dropna()
        counts = []
        props = []
        for cat in categories:
            if isinstance(cat, str):
                c = (other_vals > max_val).sum()
                p = (other_vals > max_val).mean() if len(other_vals) > 0 else 0
            else:
                c = (other_vals == cat).sum()
                p = (other_vals == cat).mean() if len(other_vals) > 0 else 0
            counts.append(c)
            props.append(p)
        count_data.append(counts)
        proportions.append(props)

        # Stacked bars
        x = np.arange(len(labels))
        bar_width = 0.65
        cmap = plt.cm.viridis(np.linspace(0.15, 0.85, len(categories)))

        bottom = np.zeros(len(labels))
        bars_for_legend = []
        for cat_idx, cat in enumerate(categories):
            heights = [proportions[i][cat_idx] for i in range(len(labels))]
            bar = ax.bar(x, heights, bar_width, bottom=bottom, color=cmap[cat_idx],
                        label=str(cat), edgecolor='white', linewidth=0.5)
            bars_for_legend.append(bar[0])
            bottom += heights

        # Save legend info from first plot
        if ax_idx == 0:
            legend_handles = bars_for_legend
            legend_labels = [str(c) for c in categories]

        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=10)
        ax.set_ylabel('Proportion of Reads', fontsize=11)
        ax.set_ylim(0, 1.15)  # Extra space for p-value annotation

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Chi-square test: each target cluster vs others
        p_values = []
        for i, cid in enumerate(target_clusters):
            # Create contingency table: cluster vs others
            obs_cluster = np.array(count_data[i])
            obs_others = np.array(count_data[-1])

            # Only test if we have enough counts
            if obs_cluster.sum() >= 5 and obs_others.sum() >= 5:
                contingency = np.array([obs_cluster, obs_others])
                # Remove zero columns to avoid issues
                nonzero_cols = contingency.sum(axis=0) > 0
                if nonzero_cols.sum() >= 2:
                    contingency = contingency[:, nonzero_cols]
                    try:
                        chi2, p_val, dof, expected = stats.chi2_contingency(contingency)
                        p_values.append(p_val)
                    except (ValueError, ZeroDivisionError):
                        p_values.append(1.0)
                else:
                    p_values.append(1.0)
            else:
                p_values.append(1.0)

        # Add p-value annotations above each target cluster bar
        for i, (cid, p_val) in enumerate(zip(target_clusters, p_values)):
            stars = format_pvalue_stars(p_val)
            if stars != "ns":
                ax.text(i, 1.02, stars, ha='center', va='bottom', fontsize=10, fontweight='bold')
            else:
                ax.text(i, 1.02, "ns", ha='center', va='bottom', fontsize=8, style='italic',
                       color='gray' if not dark_mode else '#888888')

        # Title with test name
        ax.set_title(f"{format_metric_title(metric)}", fontsize=12, fontweight='bold', pad=10)

    # Add shared legend below plots
    fig.legend(legend_handles, legend_labels, title='Count', loc='lower center',
              bbox_to_anchor=(0.5, -0.02), ncol=len(legend_labels), fontsize=9,
              title_fontsize=10, frameon=False)

    # Add test description
    fig.text(0.5, -0.08, "Chi-square test vs Others; * p<0.05, ** p<0.01, *** p<0.001, **** p<0.0001",
             ha='center', fontsize=9, style='italic')

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.18)  # Make room for legend
    return fig


def pub_summary_heatmap(
    df: pd.DataFrame,
    metrics: List[str],
    enrichment_map: Dict[int, str] = None,
    q_value_map: Dict[int, float] = None,
    cluster_order: str = 'enrichment',
    show_significance: bool = True,
    dark_mode: bool = False,
    figsize: Tuple[float, float] = None,
    min_cluster_size: int = 5,
) -> plt.Figure:
    """
    Create publication-quality heatmap showing z-scored metrics for all clusters.

    Clusters are colored by enrichment direction in row labels.
    Significance markers (* ** ***) indicate FDR-significant clusters.
    """
    plt.rcParams.update(PUBLICATION_RCPARAMS)

    metrics = [m for m in metrics if m in df.columns]
    clusters = sorted(df['cluster'].unique())

    # Filter by size
    cluster_sizes = df.groupby('cluster').size()
    clusters = [c for c in clusters if cluster_sizes[c] >= min_cluster_size]

    # Build matrix
    data = []
    for cid in clusters:
        cluster_df = df[df['cluster'] == cid]
        row = [cluster_df[m].median() for m in metrics]
        data.append(row)

    matrix = pd.DataFrame(data, index=clusters, columns=metrics)
    # Handle constant columns (std=0) by replacing with 1 to avoid NaN
    std = matrix.std()
    std[std == 0] = 1.0
    z_matrix = (matrix - matrix.mean()) / std

    # Order clusters
    if cluster_order == 'enrichment' and enrichment_map:
        def sort_key(cid):
            enr = enrichment_map.get(cid, 'zzz')
            return (0 if 'E6E7' in enr else 1 if 'primary' in enr else 2, -cluster_sizes[cid])
        clusters = sorted(clusters, key=sort_key)
    elif cluster_order == 'size':
        clusters = sorted(clusters, key=lambda c: -cluster_sizes[c])
    elif cluster_order == 'hierarchical':
        from scipy.cluster.hierarchy import linkage, leaves_list
        if len(clusters) > 2:
            Z = linkage(z_matrix.values, method='ward')
            order = leaves_list(Z)
            clusters = [clusters[i] for i in order]

    z_matrix = z_matrix.loc[clusters]

    # Figure size
    if figsize is None:
        figsize = (3 + len(metrics) * 0.9, max(4, 1 + len(clusters) * 0.22))

    fig, ax = plt.subplots(figsize=figsize)

    cmap = 'RdBu_r'
    im = ax.imshow(z_matrix.values, aspect='auto', cmap=cmap, vmin=-3, vmax=3)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label('Z-score (median)', fontsize=10)
    cbar.ax.tick_params(labelsize=9)

    # Y-axis (cluster labels)
    y_labels = []
    y_colors = []
    for cid in clusters:
        label = str(cid)
        color = 'black'

        if enrichment_map:
            enr = enrichment_map.get(cid, '')
            if 'E6E7' in enr:
                color = ENRICHMENT_COLORS['E6E7']
            elif 'primary' in enr:
                color = ENRICHMENT_COLORS['primary']
            else:
                color = 'gray'

        if show_significance and q_value_map:
            q = q_value_map.get(cid, 1.0)
            if q < 0.001:
                label += " ***"
            elif q < 0.01:
                label += " **"
            elif q < 0.05:
                label += " *"

        y_labels.append(label)
        y_colors.append(color)

    ax.set_yticks(range(len(clusters)))
    ax.set_yticklabels(y_labels, fontsize=8)
    for i, (label, color) in enumerate(zip(ax.get_yticklabels(), y_colors)):
        label.set_color(color)
        label.set_fontweight('bold' if q_value_map and q_value_map.get(clusters[i], 1) < 0.05 else 'normal')

    # X-axis (metrics)
    x_labels = [format_metric_title(m) for m in metrics]
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=10)

    ax.set_title(f"Cluster Metric Summary (n={len(clusters)} clusters, size ≥ {min_cluster_size})",
                 fontsize=12, fontweight='bold', pad=15)

    # Legend for enrichment colors
    if enrichment_map:
        legend_elements = [
            mpatches.Patch(facecolor=ENRICHMENT_COLORS['E6E7'], label='E6E7-enriched'),
            mpatches.Patch(facecolor=ENRICHMENT_COLORS['primary'], label='Primary-enriched'),
            mpatches.Patch(facecolor='gray', label='Mixed/NS'),
        ]
        ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 0.5),
                 fontsize=8, frameon=False, title='Enrichment', title_fontsize=9)

    plt.tight_layout()
    return fig


def pub_statistics_table(
    df: pd.DataFrame,
    target_clusters: List[int],
    metrics: List[str],
    output_path: str,
) -> pd.DataFrame:
    """
    Generate detailed statistics table for publication supplementary material.
    """
    other_df = df[~df['cluster'].isin(target_clusters)]

    rows = []
    for metric in metrics:
        if metric not in df.columns:
            continue

        for cid in target_clusters:
            target_vals = df[df['cluster'] == cid][metric].dropna()
            other_vals = other_df[metric].dropna()

            if len(target_vals) == 0 or len(other_vals) == 0:
                continue

            stat, p_val = stats.mannwhitneyu(target_vals, other_vals, alternative='two-sided')
            r, effect_interp = compute_effect_size(target_vals.values, other_vals.values)

            rows.append({
                'metric': metric,
                'cluster': cid,
                'cluster_n': len(target_vals),
                'cluster_median': round(target_vals.median(), 4),
                'cluster_mean': round(target_vals.mean(), 4),
                'cluster_std': round(target_vals.std(), 4),
                'cluster_q25': round(target_vals.quantile(0.25), 4),
                'cluster_q75': round(target_vals.quantile(0.75), 4),
                'others_n': len(other_vals),
                'others_median': round(other_vals.median(), 4),
                'others_mean': round(other_vals.mean(), 4),
                'others_std': round(other_vals.std(), 4),
                'U_statistic': stat,
                'p_value': p_val,
                'p_value_formatted': format_pvalue(p_val),
                'significance': format_pvalue_stars(p_val),
                'effect_size_r': round(r, 3),
                'effect_interpretation': effect_interp,
            })

    result_df = pd.DataFrame(rows)
    result_df.to_csv(output_path, sep='\t', index=False)
    return result_df


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
    parser.add_argument("--format", dest="output_format", default="pdf",
                        choices=["pdf", "svg", "png"],
                        help="Output format for publication figures (default: pdf)")

    # Exploratory comparison options
    parser.add_argument("--compare-clusters", dest="compare_clusters", default=None,
                        help="Comma-separated cluster IDs to compare vs all others (e.g., '82,88')")
    parser.add_argument("--compare-significant", dest="compare_significant", action="store_true", default=False,
                        help="Compare each FDR-significant cluster vs all others (requires --cluster-analysis)")

    # Publication-quality figure options
    pub_group = parser.add_argument_group('Publication figures',
                                          'Generate publication-quality figures with statistical annotations')
    pub_group.add_argument("--pub-comparison", dest="pub_comparison", action="store_true", default=False,
                           help="Generate publication-quality comparison figure (requires --clusters)")
    pub_group.add_argument("--pub-counts", dest="pub_counts", action="store_true", default=False,
                           help="Generate publication-quality count distribution figure (requires --clusters)")
    pub_group.add_argument("--pub-heatmap", dest="pub_heatmap", action="store_true", default=False,
                           help="Generate publication-quality summary heatmap of all clusters")
    pub_group.add_argument("--pub-stats", dest="pub_stats", action="store_true", default=False,
                           help="Generate detailed statistics table (requires --clusters)")
    pub_group.add_argument("--pub-all", dest="pub_all", action="store_true", default=False,
                           help="Generate all publication figures")
    pub_group.add_argument("--clusters", dest="clusters", default=None,
                           help="Comma-separated cluster IDs for publication figures (e.g., '82,88')")
    pub_group.add_argument("--cluster-labels", dest="cluster_labels", default=None,
                           help="Custom labels for clusters (e.g., '82:α-sat,88:rDNA')")
    pub_group.add_argument("--heatmap-order", dest="heatmap_order", default="enrichment",
                           choices=["enrichment", "size", "hierarchical"],
                           help="How to order clusters in heatmap (default: enrichment)")
    pub_group.add_argument("--min-cluster-size", dest="min_cluster_size", type=int, default=5,
                           help="Minimum cluster size for heatmap (default: 5)")
    pub_group.add_argument("--no-exploratory", dest="no_exploratory", action="store_true", default=False,
                           help="Skip exploratory diagnostic plots (only generate publication figures)")

    args = parser.parse_args()

    # Set up plot style
    global STRIP_COLOR, STRIP_ALPHA
    if args.dark_mode:
        plt.style.use('dark_background')
        sns.set_style("dark", {
            'axes.facecolor': 'black',
            'figure.facecolor': 'black',
            'axes.edgecolor': 'white',
            'axes.labelcolor': 'white',
            'xtick.color': 'white',
            'ytick.color': 'white',
            'text.color': 'white',
            'grid.color': '#555555'
        })
        # Override stripplot color for visibility
        STRIP_COLOR = 'white'
        STRIP_ALPHA = 0.5
    else:
        sns.set_style("whitegrid")
        STRIP_COLOR = 'black'
        STRIP_ALPHA = 0.3

    # Re-apply fonttype settings (style.use resets them)
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['svg.fonttype'] = 'none'

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

    # Track output files
    output_files = []

    # =========================================================================
    # Exploratory diagnostic plots
    # =========================================================================

    # === Page 1: Metric distributions by cluster ===
    print("\nGenerating metric plots...")
    metrics_pdf = f"{args.output_prefix}.cluster_metrics.pdf"

    with PdfPages(metrics_pdf) as pdf:
        # Page 1: Key metrics by cluster (3x2 grid)
        fig, axes = plt.subplots(3, 2, figsize=(14, 14))
        fig.suptitle('Sequence Metrics by Cluster', fontsize=14, fontweight='bold')

        plot_metric_by_cluster(df, 'read_length', axes[0, 0],
                               title='Read Length by Cluster', ylabel='Read Length (bp)')
        plot_metric_by_cluster(df, 'primary_align_fraction', axes[0, 1],
                               title='Primary Alignment Fraction by Cluster', ylabel='Alignment Fraction')
        plot_metric_by_cluster(df, 'primary_de', axes[1, 0],
                               title='Divergence (Error Rate) by Cluster', ylabel='Divergence')
        plot_metric_by_cluster(df, 'centroid_distance', axes[1, 1],
                               title='Centroid Distance by Cluster', ylabel='Centroid Distance')
        plot_metric_by_cluster(df, 'primary_mapq', axes[2, 0],
                               title='Mapping Quality by Cluster', ylabel='MAPQ')
        plot_metric_by_cluster(df, 'primary_align_len', axes[2, 1],
                               title='Primary Alignment Length by Cluster', ylabel='Alignment Length (bp)')

        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # Page 2: Metrics by cluster, colored by sample (3x2 grid)
        if 'sample' in df.columns and df['sample'].nunique() > 1:
            fig, axes = plt.subplots(3, 2, figsize=(16, 14))
            fig.suptitle('Sequence Metrics by Cluster (Colored by Sample)', fontsize=14, fontweight='bold')

            plot_metric_by_cluster_and_sample(df, 'read_length', axes[0, 0],
                                              title='Read Length by Cluster', ylabel='Read Length (bp)')
            plot_metric_by_cluster_and_sample(df, 'primary_align_fraction', axes[0, 1],
                                              title='Primary Alignment Fraction by Cluster', ylabel='Alignment Fraction')
            plot_metric_by_cluster_and_sample(df, 'primary_de', axes[1, 0],
                                              title='Divergence by Cluster', ylabel='Divergence')
            plot_metric_by_cluster_and_sample(df, 'centroid_distance', axes[1, 1],
                                              title='Centroid Distance by Cluster', ylabel='Centroid Distance')
            plot_metric_by_cluster_and_sample(df, 'primary_mapq', axes[2, 0],
                                              title='Mapping Quality by Cluster', ylabel='MAPQ')
            plot_metric_by_cluster_and_sample(df, 'primary_align_len', axes[2, 1],
                                              title='Primary Alignment Length by Cluster', ylabel='Alignment Length (bp)')

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

        # Page 3: Metrics by enrichment (if available, 3x2 grid)
        if 'enrichment' in df.columns:
            fig, axes = plt.subplots(3, 2, figsize=(12, 14))
            fig.suptitle('Sequence Metrics by Enrichment Category', fontsize=14, fontweight='bold')

            plot_metric_by_enrichment(df, 'read_length', axes[0, 0],
                                      title='Read Length by Enrichment', ylabel='Read Length (bp)')
            plot_metric_by_enrichment(df, 'primary_align_fraction', axes[0, 1],
                                      title='Primary Alignment Fraction by Enrichment', ylabel='Alignment Fraction')
            plot_metric_by_enrichment(df, 'primary_de', axes[1, 0],
                                      title='Divergence by Enrichment', ylabel='Divergence')
            plot_metric_by_enrichment(df, 'centroid_distance', axes[1, 1],
                                      title='Centroid Distance by Enrichment', ylabel='Centroid Distance')
            plot_metric_by_enrichment(df, 'primary_mapq', axes[2, 0],
                                      title='Mapping Quality by Enrichment', ylabel='MAPQ')
            plot_metric_by_enrichment(df, 'primary_align_len', axes[2, 1],
                                      title='Primary Alignment Length by Enrichment', ylabel='Alignment Length (bp)')

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

        # Page 4: Metrics by sample (batch effect detection, 3x2 grid)
        if 'sample' in df.columns and df['sample'].nunique() > 1:
            fig, axes = plt.subplots(3, 2, figsize=(12, 14))
            fig.suptitle('Sequence Metrics by Sample (Batch Effect Check)', fontsize=14, fontweight='bold')

            plot_metric_by_sample(df, 'read_length', axes[0, 0],
                                  title='Read Length by Sample', ylabel='Read Length (bp)')
            plot_metric_by_sample(df, 'primary_align_fraction', axes[0, 1],
                                  title='Primary Alignment Fraction by Sample', ylabel='Alignment Fraction')
            plot_metric_by_sample(df, 'primary_de', axes[1, 0],
                                  title='Divergence by Sample', ylabel='Divergence')
            plot_metric_by_sample(df, 'centroid_distance', axes[1, 1],
                                  title='Centroid Distance by Sample', ylabel='Centroid Distance')
            plot_metric_by_sample(df, 'primary_mapq', axes[2, 0],
                                  title='Mapping Quality by Sample', ylabel='MAPQ')
            plot_metric_by_sample(df, 'primary_align_len', axes[2, 1],
                                  title='Primary Alignment Length by Sample', ylabel='Alignment Length (bp)')

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

    # === Cluster comparison plots ===
    comparison_metrics = [
        'read_length', 'centroid_distance', 'primary_mapq', 'primary_de',
        'primary_align_len', 'primary_align_fraction', 'total_align_fraction',
        'n_alignments', 'n_secondary', 'n_supplementary'
    ]

    output_files.extend([metrics_pdf, composition_pdf, summary_file])

    # Get enrichment info if available
    enrichment_map = {}
    q_value_map = {}
    if 'enrichment' in df.columns:
        # Try to get q-values from cluster analysis
        if args.cluster_analysis and os.path.exists(args.cluster_analysis):
            ca = pd.read_csv(args.cluster_analysis, sep='\t')
            if 'q_value' in ca.columns:
                q_value_map = dict(zip(ca['cluster_id'], ca['q_value']))
            enrichment_map = dict(zip(ca['cluster_id'], ca['enrichment']))

    # Compare user-specified clusters (same stacked style as --compare-significant)
    if args.compare_clusters:
        cluster_ids = [int(c.strip()) for c in args.compare_clusters.split(',')]
        # Filter to clusters that exist
        cluster_ids = [c for c in cluster_ids if c in df['cluster'].values]
        if not cluster_ids:
            print(f"\nWarning: None of the specified clusters found")
        else:
            # All other clusters become background
            other_clusters = [c for c in df['cluster'].unique() if c not in cluster_ids]

            print(f"\nGenerating comparison plot for clusters {cluster_ids} vs {len(other_clusters)} others")

            compare_pdf = f"{args.output_prefix}.cluster_comparison.pdf"
            with PdfPages(compare_pdf) as pdf:
                fig = plot_significant_vs_nonsig(
                    df, cluster_ids, other_clusters, comparison_metrics,
                    enrichment_map=enrichment_map, q_value_map=q_value_map,
                    dark_mode=args.dark_mode
                )
                if fig:
                    pdf.savefig(fig, bbox_inches='tight')
                    plt.close(fig)

            print(f"  Saved: {compare_pdf}")
            output_files.append(compare_pdf)

    # Compare all FDR-significant clusters vs pooled non-significant
    if args.compare_significant:
        if 'enrichment' not in df.columns:
            print("\nWarning: --compare-significant requires --cluster-analysis with enrichment info")
        elif not q_value_map:
            print("\nWarning: --compare-significant requires cluster_analysis with q_value column")
        else:
            # Find significant clusters (q < 0.05 and not 'mixed')
            significant_clusters = [
                cid for cid, enr in enrichment_map.items()
                if enr != 'mixed' and q_value_map.get(cid, 1.0) < 0.05
            ]
            significant_clusters = sorted(significant_clusters)

            # Find non-significant clusters (mixed or q >= 0.05) as background
            nonsig_clusters = [
                cid for cid in df['cluster'].unique()
                if cid not in significant_clusters
            ]

            if significant_clusters:
                print(f"\nGenerating comparison plot: {len(significant_clusters)} significant vs {len(nonsig_clusters)} non-significant clusters")

                sig_pdf = f"{args.output_prefix}.significant_clusters_comparison.pdf"
                with PdfPages(sig_pdf) as pdf:
                    fig = plot_significant_vs_nonsig(
                        df, significant_clusters, nonsig_clusters, comparison_metrics,
                        enrichment_map=enrichment_map, q_value_map=q_value_map,
                        dark_mode=args.dark_mode
                    )
                    if fig:
                        pdf.savefig(fig, bbox_inches='tight')
                        plt.close(fig)

                print(f"  Saved: {sig_pdf}")
                output_files.append(sig_pdf)
            else:
                print("\nNo FDR-significant clusters found (q < 0.05)")

    # =========================================================================
    # Publication-quality figures
    # =========================================================================

    pub_any = args.pub_all or args.pub_comparison or args.pub_counts or args.pub_heatmap or args.pub_stats

    if pub_any:
        print(f"\n{'=' * 60}")
        print("Publication Figures")
        print("=" * 60)

        # Parse target clusters
        target_clusters = []
        if args.clusters:
            target_clusters = [int(c.strip()) for c in args.clusters.split(',')]
            target_clusters = [c for c in target_clusters if c in df['cluster'].values]

        # Parse cluster labels
        cluster_labels = {}
        if args.cluster_labels:
            for item in args.cluster_labels.split(','):
                if ':' in item:
                    cid, label = item.split(':', 1)
                    cluster_labels[int(cid.strip())] = label.strip()

        # Metrics for publication figures
        continuous_metrics = ['read_length', 'primary_mapq', 'primary_de',
                             'primary_align_fraction', 'total_align_fraction']
        count_metrics = ['n_alignments', 'n_secondary', 'n_supplementary']
        all_pub_metrics = continuous_metrics + count_metrics

        fmt = args.output_format

        # Comparison figure
        if (args.pub_all or args.pub_comparison) and target_clusters:
            print(f"\nGenerating publication comparison figure...")
            try:
                fig = pub_comparison_figure(
                    df, target_clusters, continuous_metrics,
                    enrichment_map=enrichment_map,
                    cluster_labels=cluster_labels,
                    dark_mode=args.dark_mode,
                )
                outpath = f"{args.output_prefix}.pub_comparison.{fmt}"
                fig.savefig(outpath, dpi=300, bbox_inches='tight')
                plt.close(fig)
                print(f"  Saved: {outpath}")
                output_files.append(outpath)
            except Exception as e:
                print(f"  Error: {e}")

        # Count distribution figure
        if (args.pub_all or args.pub_counts) and target_clusters:
            print(f"\nGenerating publication count distribution figure...")
            try:
                fig = pub_count_distribution_figure(
                    df, target_clusters, count_metrics,
                    cluster_labels=cluster_labels,
                    dark_mode=args.dark_mode,
                )
                if fig:
                    outpath = f"{args.output_prefix}.pub_counts.{fmt}"
                    fig.savefig(outpath, dpi=300, bbox_inches='tight')
                    plt.close(fig)
                    print(f"  Saved: {outpath}")
                    output_files.append(outpath)
            except Exception as e:
                print(f"  Error: {e}")

        # Summary heatmap
        if args.pub_all or args.pub_heatmap:
            print(f"\nGenerating publication summary heatmap...")
            try:
                fig = pub_summary_heatmap(
                    df, all_pub_metrics,
                    enrichment_map=enrichment_map,
                    q_value_map=q_value_map,
                    cluster_order=args.heatmap_order,
                    dark_mode=args.dark_mode,
                    min_cluster_size=args.min_cluster_size,
                )
                outpath = f"{args.output_prefix}.pub_heatmap.{fmt}"
                fig.savefig(outpath, dpi=300, bbox_inches='tight')
                plt.close(fig)
                print(f"  Saved: {outpath}")
                output_files.append(outpath)
            except Exception as e:
                print(f"  Error: {e}")

        # Statistics table
        if (args.pub_all or args.pub_stats) and target_clusters:
            print(f"\nGenerating publication statistics table...")
            try:
                outpath = f"{args.output_prefix}.pub_statistics.tsv"
                pub_statistics_table(df, target_clusters, all_pub_metrics, outpath)
                print(f"  Saved: {outpath}")
                output_files.append(outpath)
            except Exception as e:
                print(f"  Error: {e}")

    print(f"\n{'=' * 60}")
    print("Summary")
    print("=" * 60)
    print(f"Total sequences: {len(df)}")
    print(f"Clusters: {df['cluster'].nunique()}")
    if 'enrichment' in df.columns:
        print(f"Enrichment categories: {df['enrichment'].value_counts().to_dict()}")
    print(f"\nOutput files:")
    for f in output_files:
        print(f"  - {f}")


if __name__ == "__main__":
    main()
