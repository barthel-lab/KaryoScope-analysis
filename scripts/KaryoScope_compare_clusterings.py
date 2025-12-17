#!/usr/bin/env python3
"""
KaryoScope Clustering Comparison

Compares clustering results from two different featuresets (e.g., region vs repeat)
to assess concordance and identify complementary signals.

Usage:
  python KaryoScope_compare_clusterings.py \
    --clustering1 analysis1.read_assignments.tsv \
    --clustering2 analysis2.read_assignments.tsv \
    --label1 "region" \
    --label2 "repeat" \
    --output-prefix comparison

Generates:
  - {prefix}.comparison_report.txt: Text summary of comparison
  - {prefix}.comparison_plots.pdf: Visualization of clustering concordance
  - {prefix}.comparison_matrix.tsv: Cluster-to-cluster mapping matrix
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


def load_clustering(assignments_file, analysis_file=None):
    """Load read assignments and optionally cluster analysis for enrichment labels."""
    df = pd.read_csv(assignments_file, sep='\t')

    # Try to auto-discover cluster analysis file
    if analysis_file is None:
        analysis_file = assignments_file.replace('.read_assignments.tsv', '.cluster_analysis.tsv')

    enrichment_map = {}
    if os.path.exists(analysis_file):
        ca = pd.read_csv(analysis_file, sep='\t')
        enrichment_map = dict(zip(ca['cluster_id'], ca['enrichment']))

    df['enrichment'] = df['cluster'].map(enrichment_map).fillna('unknown')

    return df, enrichment_map


def compute_cluster_overlap_matrix(merged, col1, col2):
    """Compute overlap matrix between two clusterings."""
    return pd.crosstab(merged[col1], merged[col2])


def compute_enrichment_concordance(merged, enrich1, enrich2):
    """Compute concordance between enrichment labels."""
    return pd.crosstab(merged[enrich1], merged[enrich2], margins=True)


def compute_adjusted_rand_index(labels1, labels2):
    """Compute Adjusted Rand Index between two clusterings."""
    from sklearn.metrics import adjusted_rand_index
    return adjusted_rand_index(labels1, labels2)


def compute_normalized_mutual_info(labels1, labels2):
    """Compute Normalized Mutual Information between two clusterings."""
    from sklearn.metrics import normalized_mutual_info_score
    return normalized_mutual_info_score(labels1, labels2)


def plot_enrichment_sankey(merged, enrich1, enrich2, label1, label2, ax):
    """Create a heatmap showing flow between enrichment categories."""
    ct = pd.crosstab(merged[enrich1], merged[enrich2])

    # Normalize by row (what fraction of each label1 category goes to each label2 category)
    ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100

    sns.heatmap(ct_pct, annot=True, fmt='.1f', cmap='Blues', ax=ax,
                cbar_kws={'label': '% of reads'})
    ax.set_xlabel(f'{label2} enrichment')
    ax.set_ylabel(f'{label1} enrichment')
    ax.set_title(f'Enrichment Flow: {label1} → {label2}\n(row-normalized percentages)')


def plot_cluster_size_comparison(df1, df2, label1, label2, ax):
    """Compare cluster size distributions."""
    sizes1 = df1.groupby('cluster').size()
    sizes2 = df2.groupby('cluster').size()

    ax.hist(sizes1, bins=20, alpha=0.5, label=f'{label1} (n={len(sizes1)})', color='blue')
    ax.hist(sizes2, bins=20, alpha=0.5, label=f'{label2} (n={len(sizes2)})', color='orange')
    ax.set_xlabel('Cluster size')
    ax.set_ylabel('Number of clusters')
    ax.set_title('Cluster Size Distributions')
    ax.legend()


def plot_top_cluster_mappings(merged, col1, col2, enrich1_map, enrich2_map, label1, label2, ax):
    """Show top cluster-to-cluster mappings."""
    cluster_map = merged.groupby([col1, col2]).size().reset_index(name='count')
    cluster_map = cluster_map.sort_values('count', ascending=True).tail(15)

    # Create labels with enrichment info
    labels = []
    for _, row in cluster_map.iterrows():
        e1 = enrich1_map.get(row[col1], 'unk')[:8]
        e2 = enrich2_map.get(row[col2], 'unk')[:8]
        labels.append(f"{label1[:3]}{row[col1]}({e1}) → {label2[:3]}{row[col2]}({e2})")

    colors = []
    for _, row in cluster_map.iterrows():
        e1 = enrich1_map.get(row[col1], 'mixed')
        e2 = enrich2_map.get(row[col2], 'mixed')
        if 'Post' in e1 and 'Post' in e2:
            colors.append('#377EB8')  # Blue - both Post
        elif 'Pre' in e1 and 'Pre' in e2:
            colors.append('#E41A1C')  # Red - both Pre
        elif e1 == e2:
            colors.append('#4DAF4A')  # Green - same
        else:
            colors.append('#999999')  # Gray - different

    ax.barh(range(len(labels)), cluster_map['count'], color=colors)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel('Number of shared reads')
    ax.set_title(f'Top Cluster Mappings: {label1} → {label2}')


def plot_concordance_by_cluster(merged, col1, enrich1, enrich2, label1, label2, ax):
    """For each cluster in clustering1, show enrichment distribution in clustering2."""
    # Get clusters sorted by size
    cluster_sizes = merged.groupby(col1).size().sort_values(ascending=False)
    top_clusters = cluster_sizes.head(15).index

    data = []
    for cluster in top_clusters:
        cluster_data = merged[merged[col1] == cluster]
        enrich_counts = cluster_data[enrich2].value_counts()
        total = len(cluster_data)
        for enrich, count in enrich_counts.items():
            data.append({
                'cluster': cluster,
                'enrichment': enrich,
                'count': count,
                'pct': 100 * count / total
            })

    plot_df = pd.DataFrame(data)
    pivot = plot_df.pivot(index='cluster', columns='enrichment', values='pct').fillna(0)
    pivot = pivot.reindex(top_clusters)

    pivot.plot(kind='barh', stacked=True, ax=ax, colormap='Set2')
    ax.set_xlabel(f'% in {label2} enrichment category')
    ax.set_ylabel(f'{label1} cluster')
    ax.set_title(f'{label1} Clusters → {label2} Enrichment Distribution')
    ax.legend(title=f'{label2} enrichment', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)


def generate_report(merged, df1, df2, enrich1_map, enrich2_map, label1, label2, output_file):
    """Generate text report of comparison findings."""
    with open(output_file, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write(f"KARYOSCOPE CLUSTERING COMPARISON: {label1.upper()} vs {label2.upper()}\n")
        f.write("=" * 70 + "\n\n")

        # Basic stats
        f.write("--- BASIC STATISTICS ---\n\n")
        f.write(f"Total reads compared: {len(merged)}\n")
        f.write(f"{label1} clusters: {merged[f'{label1}_cluster'].nunique()}\n")
        f.write(f"{label2} clusters: {merged[f'{label2}_cluster'].nunique()}\n\n")

        # Enrichment counts
        f.write(f"{label1} enrichment distribution:\n")
        for enrich, count in merged[f'{label1}_enrich'].value_counts().items():
            f.write(f"  {enrich}: {count} reads ({100*count/len(merged):.1f}%)\n")
        f.write(f"\n{label2} enrichment distribution:\n")
        for enrich, count in merged[f'{label2}_enrich'].value_counts().items():
            f.write(f"  {enrich}: {count} reads ({100*count/len(merged):.1f}%)\n")

        # Concordance metrics
        f.write("\n--- CONCORDANCE METRICS ---\n\n")

        # Enrichment agreement
        same_enrich = (merged[f'{label1}_enrich'] == merged[f'{label2}_enrich']).sum()
        f.write(f"Enrichment label agreement: {same_enrich} / {len(merged)} ({100*same_enrich/len(merged):.1f}%)\n")

        # Try to compute ARI and NMI
        try:
            ari = compute_adjusted_rand_index(
                merged[f'{label1}_cluster'].values,
                merged[f'{label2}_cluster'].values
            )
            f.write(f"Adjusted Rand Index: {ari:.4f}\n")
        except ImportError:
            f.write("Adjusted Rand Index: (sklearn not available)\n")

        try:
            nmi = compute_normalized_mutual_info(
                merged[f'{label1}_cluster'].values,
                merged[f'{label2}_cluster'].values
            )
            f.write(f"Normalized Mutual Information: {nmi:.4f}\n")
        except ImportError:
            f.write("Normalized Mutual Information: (sklearn not available)\n")

        # Enrichment cross-tabulation
        f.write("\n--- ENRICHMENT CROSS-TABULATION ---\n\n")
        ct = pd.crosstab(merged[f'{label1}_enrich'], merged[f'{label2}_enrich'], margins=True)
        f.write(ct.to_string() + "\n")

        # Flow analysis
        f.write("\n--- ENRICHMENT FLOW ANALYSIS ---\n\n")

        for enrich_cat in ['SW26_Post-enriched', 'SW26_Pre-enriched', 'mixed']:
            subset = merged[merged[f'{label1}_enrich'] == enrich_cat]
            if len(subset) == 0:
                continue
            f.write(f"Reads in {label1} '{enrich_cat}' clusters (n={len(subset)}):\n")
            for enrich2, count in subset[f'{label2}_enrich'].value_counts().items():
                f.write(f"  → {label2} '{enrich2}': {count} ({100*count/len(subset):.1f}%)\n")
            f.write("\n")

        # Top cluster mappings
        f.write("--- TOP CLUSTER MAPPINGS ---\n\n")
        cluster_map = merged.groupby([f'{label1}_cluster', f'{label2}_cluster']).size().reset_index(name='count')
        cluster_map = cluster_map.sort_values('count', ascending=False).head(20)

        f.write(f"{'Cluster1':>10} {'Enrich1':>15} {'Cluster2':>10} {'Enrich2':>15} {'Reads':>8}\n")
        f.write("-" * 60 + "\n")
        for _, row in cluster_map.iterrows():
            e1 = enrich1_map.get(row[f'{label1}_cluster'], 'unknown')
            e2 = enrich2_map.get(row[f'{label2}_cluster'], 'unknown')
            f.write(f"{row[f'{label1}_cluster']:>10} {e1:>15} {row[f'{label2}_cluster']:>10} {e2:>15} {row['count']:>8}\n")

        # Interpretation
        f.write("\n--- INTERPRETATION ---\n\n")

        if same_enrich / len(merged) < 0.5:
            f.write("LOW CONCORDANCE: The two featuresets capture different signals.\n")
            f.write("  - Reads that cluster together by one featureset often don't\n")
            f.write("    cluster together by the other.\n")
            f.write("  - Consider combining featuresets for richer analysis.\n")
        elif same_enrich / len(merged) < 0.75:
            f.write("MODERATE CONCORDANCE: The featuresets show partial overlap.\n")
            f.write("  - Some biological signal is shared between featuresets.\n")
            f.write("  - Each featureset may capture unique aspects.\n")
        else:
            f.write("HIGH CONCORDANCE: The featuresets capture similar signals.\n")
            f.write("  - Clustering results are largely consistent.\n")
            f.write("  - Either featureset may be sufficient for analysis.\n")

        f.write("\n" + "=" * 70 + "\n")
        f.write("END OF REPORT\n")
        f.write("=" * 70 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Compare clustering results from two different featuresets",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--clustering1", required=True,
                        help="First read_assignments.tsv file")
    parser.add_argument("--clustering2", required=True,
                        help="Second read_assignments.tsv file")
    parser.add_argument("--label1", default="clustering1",
                        help="Label for first clustering (e.g., 'region')")
    parser.add_argument("--label2", default="clustering2",
                        help="Label for second clustering (e.g., 'repeat')")
    parser.add_argument("--output-prefix", dest="output_prefix", required=True,
                        help="Output prefix for generated files")

    args = parser.parse_args()

    print("=" * 60)
    print("KaryoScope Clustering Comparison")
    print("=" * 60)

    # Load clusterings
    print(f"\nLoading {args.label1} clustering: {args.clustering1}")
    df1, enrich1_map = load_clustering(args.clustering1)
    print(f"  Reads: {len(df1)}, Clusters: {df1['cluster'].nunique()}")

    print(f"\nLoading {args.label2} clustering: {args.clustering2}")
    df2, enrich2_map = load_clustering(args.clustering2)
    print(f"  Reads: {len(df2)}, Clusters: {df2['cluster'].nunique()}")

    # Prepare merged dataframe
    df1_subset = df1[['read', 'cluster', 'enrichment']].copy()
    df1_subset.columns = ['read', f'{args.label1}_cluster', f'{args.label1}_enrich']

    df2_subset = df2[['read', 'cluster', 'enrichment']].copy()
    df2_subset.columns = ['read', f'{args.label2}_cluster', f'{args.label2}_enrich']

    merged = df1_subset.merge(df2_subset, on='read')
    print(f"\nMerged reads: {len(merged)}")

    if len(merged) == 0:
        print("ERROR: No overlapping reads between clusterings!")
        sys.exit(1)

    # Generate report
    print("\nGenerating comparison report...")
    report_file = f"{args.output_prefix}.comparison_report.txt"
    generate_report(merged, df1, df2, enrich1_map, enrich2_map,
                    args.label1, args.label2, report_file)
    print(f"  Saved: {report_file}")

    # Generate plots
    print("\nGenerating comparison plots...")
    plots_file = f"{args.output_prefix}.comparison_plots.pdf"

    sns.set_style("whitegrid")
    plt.rcParams['figure.dpi'] = 150

    with PdfPages(plots_file) as pdf:
        # Page 1: Enrichment flow and cluster sizes
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        fig.suptitle(f'Clustering Comparison: {args.label1} vs {args.label2}',
                     fontsize=14, fontweight='bold')

        plot_enrichment_sankey(merged, f'{args.label1}_enrich', f'{args.label2}_enrich',
                               args.label1, args.label2, axes[0])
        plot_cluster_size_comparison(df1, df2, args.label1, args.label2, axes[1])

        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # Page 2: Top cluster mappings
        fig, ax = plt.subplots(figsize=(12, 8))
        plot_top_cluster_mappings(merged, f'{args.label1}_cluster', f'{args.label2}_cluster',
                                  enrich1_map, enrich2_map, args.label1, args.label2, ax)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # Page 3: Concordance by cluster (label1 → label2)
        fig, ax = plt.subplots(figsize=(12, 8))
        plot_concordance_by_cluster(merged, f'{args.label1}_cluster',
                                    f'{args.label1}_enrich', f'{args.label2}_enrich',
                                    args.label1, args.label2, ax)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # Page 4: Concordance by cluster (label2 → label1)
        fig, ax = plt.subplots(figsize=(12, 8))
        plot_concordance_by_cluster(merged, f'{args.label2}_cluster',
                                    f'{args.label2}_enrich', f'{args.label1}_enrich',
                                    args.label2, args.label1, ax)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    print(f"  Saved: {plots_file}")

    # Save cluster mapping matrix
    print("\nGenerating cluster mapping matrix...")
    matrix_file = f"{args.output_prefix}.comparison_matrix.tsv"
    overlap_matrix = compute_cluster_overlap_matrix(
        merged, f'{args.label1}_cluster', f'{args.label2}_cluster'
    )
    overlap_matrix.to_csv(matrix_file, sep='\t')
    print(f"  Saved: {matrix_file}")

    # Print summary
    print(f"\n{'=' * 60}")
    print("Summary")
    print("=" * 60)

    same_enrich = (merged[f'{args.label1}_enrich'] == merged[f'{args.label2}_enrich']).sum()
    print(f"Enrichment agreement: {same_enrich} / {len(merged)} ({100*same_enrich/len(merged):.1f}%)")

    try:
        ari = compute_adjusted_rand_index(
            merged[f'{args.label1}_cluster'].values,
            merged[f'{args.label2}_cluster'].values
        )
        print(f"Adjusted Rand Index: {ari:.4f}")
    except ImportError:
        pass

    print(f"\nOutput files:")
    print(f"  - {report_file}")
    print(f"  - {plots_file}")
    print(f"  - {matrix_file}")


if __name__ == "__main__":
    main()
