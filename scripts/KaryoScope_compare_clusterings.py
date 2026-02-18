#!/usr/bin/env python3
"""
KaryoScope Clustering Comparison

Compares clustering results from two different featuresets (e.g., region vs repeat)
to assess concordance and identify complementary signals.

Usage:
  python KaryoScope_compare_clusterings.py \
    --clustering1 analysis1.sequence_assignments.tsv \
    --clustering2 analysis2.sequence_assignments.tsv \
    --label1 "region" \
    --label2 "repeat" \
    --output-prefix comparison

Generates:
  - {prefix}.comparison_report.txt: Text summary of comparison
  - {prefix}.comparison_plots.pdf: Visualization of clustering concordance
  - {prefix}.comparison_matrix.tsv: Cluster-to-cluster mapping matrix
  - {prefix}.comparison_labels.tsv: Auto-label cross-tabulation
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

# Keep text editable in output files (not converted to paths)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'


def load_clustering(assignments_file, analysis_file=None):
    """Load read assignments and optionally cluster analysis for enrichment labels."""
    df = pd.read_csv(assignments_file, sep='\t')

    # Normalize column name: newer files may use 'sequence' instead of 'read'
    if 'sequence' in df.columns and 'read' not in df.columns:
        df.rename(columns={'sequence': 'read'}, inplace=True)

    # Try to auto-discover cluster analysis file
    if analysis_file is None:
        analysis_file = assignments_file.replace('.sequence_assignments.tsv', '.cluster_analysis.tsv')

    enrichment_map = {}
    if os.path.exists(analysis_file):
        ca = pd.read_csv(analysis_file, sep='\t')
        enrichment_map = dict(zip(ca['cluster_id'], ca['enrichment']))

    df['enrichment'] = df['cluster'].map(enrichment_map).fillna('unknown')

    # Try to auto-discover cluster annotations file for auto-labels
    annotations_file = assignments_file.replace('.sequence_assignments.tsv', '.cluster_annotations.tsv')
    cluster_name_map = {}
    if os.path.exists(annotations_file):
        ann = pd.read_csv(annotations_file, sep='\t')
        if 'cluster_id' in ann.columns and 'cluster_name' in ann.columns:
            cluster_name_map = dict(zip(ann['cluster_id'], ann['cluster_name']))

    df['cluster_name'] = df['cluster'].map(cluster_name_map).fillna('(unlabeled)')
    df.loc[df['cluster_name'] == '', 'cluster_name'] = '(unlabeled)'

    return df, enrichment_map, cluster_name_map


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


def compute_label_retention(merged, label_col_a, label_col_b):
    """Forward mapping: for each label in A, what % of reads get the same label in B."""
    rows = []
    for label_a, grp in merged.groupby(label_col_a):
        n_reads = len(grp)
        n_clusters = grp[label_col_a.replace('_label', '_cluster')].nunique() if label_col_a.replace('_label', '_cluster') in merged.columns else 0
        same = (grp[label_col_b] == label_a).sum()
        same_pct = 100 * same / n_reads if n_reads > 0 else 0
        # Top destinations
        dest_counts = grp[label_col_b].value_counts()
        top_dests = '; '.join(f"{lbl}({100*cnt/n_reads:.0f}%)" for lbl, cnt in dest_counts.head(3).items())
        rows.append({
            'label': label_a,
            'n_reads': n_reads,
            'n_clusters': n_clusters,
            'same_label_pct': round(same_pct, 1),
            'top_destinations': top_dests,
        })
    return pd.DataFrame(rows).sort_values('n_reads', ascending=False)


def compute_label_purity(merged, label_col_a, label_col_b):
    """Reverse mapping: for each label in B, what % came from the same label in A."""
    rows = []
    for label_b, grp in merged.groupby(label_col_b):
        n_reads = len(grp)
        n_clusters = grp[label_col_b.replace('_label', '_cluster')].nunique() if label_col_b.replace('_label', '_cluster') in merged.columns else 0
        same = (grp[label_col_a] == label_b).sum()
        same_pct = 100 * same / n_reads if n_reads > 0 else 0
        # Top sources
        src_counts = grp[label_col_a].value_counts()
        top_srcs = '; '.join(f"{lbl}({100*cnt/n_reads:.0f}%)" for lbl, cnt in src_counts.head(3).items())
        rows.append({
            'label': label_b,
            'n_reads': n_reads,
            'n_clusters': n_clusters,
            'same_label_pct': round(same_pct, 1),
            'top_sources': top_srcs,
        })
    return pd.DataFrame(rows).sort_values('n_reads', ascending=False)


def plot_enrichment_sankey(merged, enrich1, enrich2, label1, label2, ax, dark_mode=False):
    """Create a heatmap showing flow between enrichment categories."""
    ct = pd.crosstab(merged[enrich1], merged[enrich2])

    # Normalize by row (what fraction of each label1 category goes to each label2 category)
    ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100

    # Use different colormap for dark mode
    cmap = 'YlGnBu' if dark_mode else 'Blues'
    annot_color = 'white' if dark_mode else 'black'

    sns.heatmap(ct_pct, annot=True, fmt='.1f', cmap=cmap, ax=ax,
                cbar_kws={'label': '% of reads'},
                annot_kws={'color': annot_color})
    ax.set_xlabel(f'{label2} enrichment')
    ax.set_ylabel(f'{label1} enrichment')
    ax.set_title(f'Enrichment Flow: {label1} → {label2}\n(row-normalized percentages)')


def plot_cluster_size_comparison(df1, df2, label1, label2, ax, colors=None):
    """Compare cluster size distributions."""
    if colors is None:
        colors = ['blue', 'orange']

    sizes1 = df1.groupby('cluster').size()
    sizes2 = df2.groupby('cluster').size()

    ax.hist(sizes1, bins=20, alpha=0.6, label=f'{label1} (n={len(sizes1)})', color=colors[0])
    ax.hist(sizes2, bins=20, alpha=0.6, label=f'{label2} (n={len(sizes2)})', color=colors[1])
    ax.set_xlabel('Cluster size')
    ax.set_ylabel('Number of clusters')
    ax.set_title('Cluster Size Distributions')
    ax.legend()


def plot_top_cluster_mappings(merged, col1, col2, enrich1_map, enrich2_map, label1, label2, ax,
                              dark_mode=False):
    """Show top cluster-to-cluster mappings."""
    cluster_map = merged.groupby([col1, col2]).size().reset_index(name='count')
    cluster_map = cluster_map.sort_values('count', ascending=True).tail(15)

    # Create labels with enrichment info
    labels = []
    for _, row in cluster_map.iterrows():
        e1 = enrich1_map.get(row[col1], 'unk')[:8]
        e2 = enrich2_map.get(row[col2], 'unk')[:8]
        labels.append(f"{label1[:3]}{row[col1]}({e1}) → {label2[:3]}{row[col2]}({e2})")

    # Adjust colors for dark mode visibility
    colors = []
    for _, row in cluster_map.iterrows():
        e1 = enrich1_map.get(row[col1], 'mixed')
        e2 = enrich2_map.get(row[col2], 'mixed')
        if 'Post' in e1 and 'Post' in e2:
            colors.append('#5dade2' if dark_mode else '#377EB8')  # Blue - both Post
        elif 'Pre' in e1 and 'Pre' in e2:
            colors.append('#ec7063' if dark_mode else '#E41A1C')  # Red - both Pre
        elif e1 == e2:
            colors.append('#58d68d' if dark_mode else '#4DAF4A')  # Green - same
        else:
            colors.append('#aab7b8' if dark_mode else '#999999')  # Gray - different

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


def generate_report(merged, df1, df2, enrich1_map, enrich2_map, label1, label2, output_file,
                     name1_map=None, name2_map=None):
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

        for enrich_cat in merged[f'{label1}_enrich'].value_counts().index:
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

        # Auto-label analysis (if labels available)
        label1_col = f'{label1}_label'
        label2_col = f'{label2}_label'
        if label1_col in merged.columns and label2_col in merged.columns:
            has_labels_1 = (merged[label1_col] != '(unlabeled)').any()
            has_labels_2 = (merged[label2_col] != '(unlabeled)').any()

            if has_labels_1 or has_labels_2:
                f.write("\n--- AUTO-LABEL RETENTION (forward: {l1} → {l2}) ---\n\n".format(l1=label1, l2=label2))
                f.write("For each {l1} label, what % of reads receive the same label in {l2}?\n\n".format(l1=label1, l2=label2))
                retention = compute_label_retention(merged, label1_col, label2_col)
                f.write(f"{'Label':<35} {'Reads':>7} {'Clust':>6} {'Same%':>7}  Top destinations\n")
                f.write("-" * 100 + "\n")
                for _, row in retention.iterrows():
                    f.write(f"{row['label']:<35} {row['n_reads']:>7} {row['n_clusters']:>6} {row['same_label_pct']:>6.1f}%  {row['top_destinations']}\n")

                f.write("\n--- AUTO-LABEL PURITY (reverse: {l1} ← {l2}) ---\n\n".format(l1=label1, l2=label2))
                f.write("For each {l2} label, what % of reads came from the same label in {l1}?\n\n".format(l1=label1, l2=label2))
                purity = compute_label_purity(merged, label1_col, label2_col)
                f.write(f"{'Label':<35} {'Reads':>7} {'Clust':>6} {'Same%':>7}  Top sources\n")
                f.write("-" * 100 + "\n")
                for _, row in purity.iterrows():
                    f.write(f"{row['label']:<35} {row['n_reads']:>7} {row['n_clusters']:>6} {row['same_label_pct']:>6.1f}%  {row['top_sources']}\n")

                f.write("\n--- AUTO-LABEL CROSS-TABULATION ---\n\n")
                label_ct = pd.crosstab(merged[label1_col], merged[label2_col], margins=True)
                f.write(label_ct.to_string() + "\n")

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
                        help="First sequence_assignments.tsv file")
    parser.add_argument("--clustering2", required=True,
                        help="Second sequence_assignments.tsv file")
    parser.add_argument("--label1", default="clustering1",
                        help="Label for first clustering (e.g., 'region')")
    parser.add_argument("--label2", default="clustering2",
                        help="Label for second clustering (e.g., 'repeat')")
    parser.add_argument("--output-prefix", dest="output_prefix", required=True,
                        help="Output prefix for generated files")
    parser.add_argument("--dark-mode", dest="dark_mode", action="store_true", default=False,
                        help="Use dark background for plots (default: False)")

    args = parser.parse_args()

    print("=" * 60)
    print("KaryoScope Clustering Comparison")
    print("=" * 60)

    # Load clusterings
    print(f"\nLoading {args.label1} clustering: {args.clustering1}")
    df1, enrich1_map, name1_map = load_clustering(args.clustering1)
    print(f"  Reads: {len(df1)}, Clusters: {df1['cluster'].nunique()}")
    if name1_map:
        print(f"  Auto-labels: {len([v for v in name1_map.values() if v])} clusters labeled")

    print(f"\nLoading {args.label2} clustering: {args.clustering2}")
    df2, enrich2_map, name2_map = load_clustering(args.clustering2)
    print(f"  Reads: {len(df2)}, Clusters: {df2['cluster'].nunique()}")
    if name2_map:
        print(f"  Auto-labels: {len([v for v in name2_map.values() if v])} clusters labeled")

    # Prepare merged dataframe
    df1_subset = df1[['read', 'cluster', 'enrichment', 'cluster_name']].copy()
    df1_subset.columns = ['read', f'{args.label1}_cluster', f'{args.label1}_enrich', f'{args.label1}_label']

    df2_subset = df2[['read', 'cluster', 'enrichment', 'cluster_name']].copy()
    df2_subset.columns = ['read', f'{args.label2}_cluster', f'{args.label2}_enrich', f'{args.label2}_label']

    merged = df1_subset.merge(df2_subset, on='read')
    print(f"\nMerged reads: {len(merged)}")

    if len(merged) == 0:
        print("ERROR: No overlapping reads between clusterings!")
        sys.exit(1)

    # Generate report
    print("\nGenerating comparison report...")
    report_file = f"{args.output_prefix}.comparison_report.txt"
    generate_report(merged, df1, df2, enrich1_map, enrich2_map,
                    args.label1, args.label2, report_file,
                    name1_map=name1_map, name2_map=name2_map)
    print(f"  Saved: {report_file}")

    # Generate plots
    print("\nGenerating comparison plots...")
    plots_file = f"{args.output_prefix}.comparison_plots.pdf"

    # Set up plot style based on dark mode
    if args.dark_mode:
        plt.style.use('dark_background')
        sns.set_style("darkgrid", {"axes.facecolor": "#1a1a1a", "grid.color": "#404040"})
        PLOT_BG_COLOR = '#1a1a1a'
        TEXT_COLOR = 'white'
        HIST_COLORS = ['#5dade2', '#f5b041']  # Light blue and orange for dark mode
    else:
        sns.set_style("whitegrid")
        PLOT_BG_COLOR = 'white'
        TEXT_COLOR = 'black'
        HIST_COLORS = ['blue', 'orange']
    plt.rcParams['figure.dpi'] = 150
    # Re-apply fonttype settings (style.use resets them)
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['svg.fonttype'] = 'none'

    with PdfPages(plots_file) as pdf:
        # Page 1: Enrichment flow and cluster sizes
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        fig.patch.set_facecolor(PLOT_BG_COLOR)
        fig.suptitle(f'Clustering Comparison: {args.label1} vs {args.label2}',
                     fontsize=14, fontweight='bold', color=TEXT_COLOR)

        plot_enrichment_sankey(merged, f'{args.label1}_enrich', f'{args.label2}_enrich',
                               args.label1, args.label2, axes[0], dark_mode=args.dark_mode)
        plot_cluster_size_comparison(df1, df2, args.label1, args.label2, axes[1],
                                     colors=HIST_COLORS)

        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight', facecolor=PLOT_BG_COLOR)
        plt.close(fig)

        # Page 2: Top cluster mappings
        fig, ax = plt.subplots(figsize=(12, 8))
        fig.patch.set_facecolor(PLOT_BG_COLOR)
        plot_top_cluster_mappings(merged, f'{args.label1}_cluster', f'{args.label2}_cluster',
                                  enrich1_map, enrich2_map, args.label1, args.label2, ax,
                                  dark_mode=args.dark_mode)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight', facecolor=PLOT_BG_COLOR)
        plt.close(fig)

        # Page 3: Concordance by cluster (label1 → label2)
        fig, ax = plt.subplots(figsize=(12, 8))
        fig.patch.set_facecolor(PLOT_BG_COLOR)
        plot_concordance_by_cluster(merged, f'{args.label1}_cluster',
                                    f'{args.label1}_enrich', f'{args.label2}_enrich',
                                    args.label1, args.label2, ax)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight', facecolor=PLOT_BG_COLOR)
        plt.close(fig)

        # Page 4: Concordance by cluster (label2 → label1)
        fig, ax = plt.subplots(figsize=(12, 8))
        fig.patch.set_facecolor(PLOT_BG_COLOR)
        plot_concordance_by_cluster(merged, f'{args.label2}_cluster',
                                    f'{args.label2}_enrich', f'{args.label1}_enrich',
                                    args.label2, args.label1, ax)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight', facecolor=PLOT_BG_COLOR)
        plt.close(fig)

        # Page 5: Auto-label flow heatmap (if labels available)
        label1_col = f'{args.label1}_label'
        label2_col = f'{args.label2}_label'
        if label1_col in merged.columns and label2_col in merged.columns:
            has_labels = ((merged[label1_col] != '(unlabeled)').any() or
                          (merged[label2_col] != '(unlabeled)').any())
            if has_labels:
                # Forward: row = label1, col = label2, row-normalized
                ct = pd.crosstab(merged[label1_col], merged[label2_col])
                ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100

                # Size figure based on number of labels
                n_rows, n_cols = ct_pct.shape
                fig_w = max(10, n_cols * 1.2 + 4)
                fig_h = max(6, n_rows * 0.6 + 3)
                fig, ax = plt.subplots(figsize=(fig_w, fig_h))
                fig.patch.set_facecolor(PLOT_BG_COLOR)

                cmap = 'YlGnBu' if args.dark_mode else 'Blues'
                # Show percentages, suppress <1% for readability
                annot_matrix = ct_pct.copy()
                sns.heatmap(ct_pct, annot=True, fmt='.0f', cmap=cmap, ax=ax,
                            cbar_kws={'label': '% of reads (row-normalized)'},
                            linewidths=0.5, vmin=0, vmax=100,
                            mask=(ct_pct < 0.5))
                ax.set_xlabel(f'{args.label2} auto-label')
                ax.set_ylabel(f'{args.label1} auto-label')
                ax.set_title(f'Auto-label Flow: {args.label1} → {args.label2}\n(row-normalized %, cells <1% hidden)')
                plt.tight_layout()
                pdf.savefig(fig, bbox_inches='tight', facecolor=PLOT_BG_COLOR)
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

    # Save auto-label cross-tabulation
    label1_col = f'{args.label1}_label'
    label2_col = f'{args.label2}_label'
    if label1_col in merged.columns and label2_col in merged.columns:
        labels_file = f"{args.output_prefix}.comparison_labels.tsv"
        label_xtab = merged.groupby([label1_col, label2_col]).size().reset_index(name='n_reads')
        # Add pct_of_a: what fraction of label_a goes to this label_b
        a_totals = merged[label1_col].value_counts()
        label_xtab['pct_of_a'] = label_xtab.apply(
            lambda r: round(100 * r['n_reads'] / a_totals[r[label1_col]], 1), axis=1)
        # Add pct_of_b: what fraction of label_b comes from this label_a
        b_totals = merged[label2_col].value_counts()
        label_xtab['pct_of_b'] = label_xtab.apply(
            lambda r: round(100 * r['n_reads'] / b_totals[r[label2_col]], 1), axis=1)
        label_xtab.columns = ['label_a', 'label_b', 'n_reads', 'pct_of_a', 'pct_of_b']
        label_xtab = label_xtab.sort_values('n_reads', ascending=False)
        label_xtab.to_csv(labels_file, sep='\t', index=False)
        print(f"  Saved: {labels_file}")

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
    if label1_col in merged.columns and label2_col in merged.columns:
        print(f"  - {labels_file}")


if __name__ == "__main__":
    main()
