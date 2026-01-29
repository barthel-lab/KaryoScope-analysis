#!/usr/bin/env python3
"""
KaryoScope Cluster Annotation

Summarizes the dominant features for each cluster based on BED file annotations.
Annotates each featureset layer separately (e.g., region, subtelomeric, chromosome).

Usage:
  python KaryoScope_cluster_annotate.py \
    --prefix analysis_output_prefix \
    --bed-dir results \
    --output cluster_annotations.tsv

  # With specific featuresets:
  python KaryoScope_cluster_annotate.py \
    --prefix analysis_output_prefix \
    --bed-dir results \
    --featuresets region,subtelomeric,chromosome \
    --output cluster_annotations.tsv

The script automatically finds these files from the prefix:
  - {prefix}.sequence_assignments.tsv
  - {prefix}.cluster_analysis.tsv
"""

import argparse
import fnmatch
import gzip
import math
import os
import sys

import pandas as pd

# Constants for log2 fold change capping (avoids Excel Inf display issues)
LOG2FC_MIN = -10
LOG2FC_MAX = 10


def _calculate_log2_fc(odds):
    """Calculate clamped log2 fold change from odds ratio."""
    if odds is None:
        return None
    if odds == 0:
        return LOG2FC_MIN
    if math.isinf(odds):
        return LOG2FC_MAX
    if odds > 0:
        log2_fc = math.log2(odds)
        return round(max(LOG2FC_MIN, min(LOG2FC_MAX, log2_fc)), 2)
    return None


def load_bed_file(filepath):
    """Load a BED file (optionally gzipped) into a DataFrame."""
    open_func = gzip.open if filepath.endswith('.gz') else open
    mode = 'rt' if filepath.endswith('.gz') else 'r'

    records = []
    with open_func(filepath, mode) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                read = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                feature = parts[3]
                length = end - start
                records.append({
                    'read': read,
                    'start': start,
                    'end': end,
                    'feature': feature,
                    'length': length
                })
    return pd.DataFrame(records)


def find_featureset_beds(bed_dir, samples, featuresets, database="KS_human_CHM13", smoothness="smoothed"):
    """Find BED files for each featureset for each sample."""
    beds_by_featureset = {fs: [] for fs in featuresets}

    for sample in samples:
        base_path = f"{bed_dir}/{sample}/telogator/1/KaryoScope/{database}"

        for fs in featuresets:
            # Try different naming patterns
            patterns = [
                f"{base_path}/{sample}.telogator.1.{database}.{fs}.{smoothness}.KaryoScope.bed",
                f"{base_path}/{sample}.telogator.1.{database}.{fs}.{smoothness}.bed",
            ]

            for pattern in patterns:
                if os.path.exists(pattern):
                    beds_by_featureset[fs].append(pattern)
                    break
                elif os.path.exists(pattern + '.gz'):
                    beds_by_featureset[fs].append(pattern + '.gz')
                    break

    return beds_by_featureset


def matches_any_pattern(feature, patterns):
    """Check if feature matches any of the exclude patterns (supports wildcards)."""
    for pattern in patterns:
        if fnmatch.fnmatch(feature, pattern):
            return True
    return False


def summarize_featureset(cluster_reads, bed_df, top_n=3, exclude_patterns=None):
    """Summarize features for a cluster from a single featureset."""
    cluster_bed = bed_df[bed_df['read'].isin(cluster_reads)]

    if len(cluster_bed) == 0:
        return ''

    # Count by feature (weighted by length)
    feature_bp = cluster_bed.groupby('feature')['length'].sum()

    if len(feature_bp) == 0:
        return ''

    # Calculate total BEFORE filtering (so percentages reflect true proportion)
    total_bp = feature_bp.sum()

    # Filter out excluded features for display only
    if exclude_patterns:
        feature_bp = feature_bp[~feature_bp.index.map(lambda f: matches_any_pattern(f, exclude_patterns))]

    if len(feature_bp) == 0:
        return ''

    # Top features by coverage (percentages relative to total including excluded)
    top_features = feature_bp.nlargest(top_n)
    top_str = '; '.join([f"{f}({100*v/total_bp:.1f}%)" for f, v in top_features.items()])

    return top_str


def main():
    parser = argparse.ArgumentParser(
        description="Annotate clusters with dominant features per featureset",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--prefix", required=True,
                        help="Analysis prefix (auto-finds {prefix}.sequence_assignments.tsv, {prefix}.cluster_analysis.tsv)")
    parser.add_argument("--bed-dir", dest="bed_dir", required=True,
                        help="Base directory containing sample BED files")
    parser.add_argument("--featuresets", default="region,subtelomeric,chromosome,acrocentric,repeat,gene",
                        help="Comma-separated featuresets to annotate (default: region,subtelomeric,chromosome,acrocentric,repeat,gene)")
    parser.add_argument("--database", default="KS_human_CHM13",
                        help="Database name (default: KS_human_CHM13)")
    parser.add_argument("--smoothness", default="smoothed",
                        choices=["smoothed", "presmoothed"],
                        help="BED file smoothness (default: smoothed)")
    parser.add_argument("--output", "-o", required=True,
                        help="Output TSV file")
    parser.add_argument("--top-n", dest="top_n", type=int, default=3,
                        help="Number of top features per featureset (default: 3)")
    parser.add_argument("--clusters",
                        help="Comma-separated cluster IDs to analyze (default: all)")
    parser.add_argument("--min-size", dest="min_size", type=int, default=1,
                        help="Minimum cluster size to include (default: 1)")
    parser.add_argument("--exclude-features", dest="exclude_features",
                        default="*multigroup*,*_arm,nonsubtelomeric,nonacrocentric,nonrepeat,categorized,canonical_telomere*",
                        help="Comma-separated features to exclude, supports wildcards (default: '*multigroup*,*_arm,nonsubtelomeric,nonacrocentric,nonrepeat,categorized,canonical_telomere*')")

    args = parser.parse_args()

    print("=" * 60)
    print("KaryoScope Cluster Annotation")
    print("=" * 60)

    # Derive file paths from prefix
    read_assignments_file = f"{args.prefix}.sequence_assignments.tsv"
    cluster_analysis_file = f"{args.prefix}.cluster_analysis.tsv"

    print(f"\nPrefix: {args.prefix}")

    # Load read assignments
    if not os.path.exists(read_assignments_file):
        print(f"ERROR: Read assignments file not found: {read_assignments_file}")
        sys.exit(1)

    print(f"\nLoading read assignments: {read_assignments_file}")
    assignments = pd.read_csv(read_assignments_file, sep='\t')
    print(f"  Total reads: {len(assignments)}")
    print(f"  Total clusters: {assignments['cluster'].nunique()}")

    # Get samples
    samples = assignments['sample'].unique().tolist()
    print(f"  Samples: {len(samples)}")

    # Load cluster analysis
    cluster_info = {}
    group_names = []
    per_sample_mode = False
    if os.path.exists(cluster_analysis_file):
        print(f"\nLoading cluster analysis: {cluster_analysis_file}")
        ca = pd.read_csv(cluster_analysis_file, sep='\t')

        # Dynamically detect group names from columns ending in _count (excluding cluster_id, size, etc.)
        count_cols = [c for c in ca.columns if c.endswith('_count')]
        group_names = [c.replace('_count', '') for c in count_cols]
        if group_names:
            print(f"  Detected groups: {group_names}")

        # Detect per-sample mode by checking for sample-specific pval columns
        sample_pval_cols = [c for c in ca.columns if c.endswith('_pval')]
        per_sample_mode = len(sample_pval_cols) > 0
        if per_sample_mode:
            print(f"  Detected per-sample comparison mode")

        for _, row in ca.iterrows():
            info = {
                'enrichment': row.get('enrichment', 'unknown'),
                'p_value': row.get('p_value', None),
                'q_value': row.get('q_value', None),
                'odds_ratio': row.get('odds_ratio', None),
            }
            # Add group-specific counts and percentages
            for grp in group_names:
                info[f'{grp}_count'] = row.get(f'{grp}_count', None)
                info[f'{grp}_pct'] = row.get(f'{grp}_pct', None)
                # Add per-sample p-values and odds ratios if available
                if per_sample_mode:
                    info[f'{grp}_pval'] = row.get(f'{grp}_pval', None)
                    info[f'{grp}_odds'] = row.get(f'{grp}_odds', None)

            cluster_info[row['cluster_id']] = info
    else:
        print(f"\nWARNING: Cluster analysis file not found: {cluster_analysis_file}")

    # Parse exclude patterns
    exclude_patterns = []
    if args.exclude_features:
        exclude_patterns = [p.strip() for p in args.exclude_features.split(',') if p.strip()]
        print(f"\nExcluding features matching: {exclude_patterns}")

    # Parse featuresets
    featuresets = [fs.strip() for fs in args.featuresets.split(',')]
    print(f"\nFeaturesets to annotate: {featuresets}")

    # Find BED files for each featureset
    print(f"\nFinding BED files in: {args.bed_dir}")
    beds_by_featureset = find_featureset_beds(
        args.bed_dir, samples, featuresets, args.database, args.smoothness
    )

    # Load BED data for each featureset
    bed_data = {}
    for fs in featuresets:
        bed_files = beds_by_featureset[fs]
        if not bed_files:
            print(f"  WARNING: No BED files found for featureset '{fs}'")
            continue

        print(f"\n  Loading {fs}: {len(bed_files)} files")
        all_data = []
        for bf in bed_files:
            df = load_bed_file(bf)
            all_data.append(df)

        bed_data[fs] = pd.concat(all_data, ignore_index=True)
        print(f"    Total records: {len(bed_data[fs])}")

    if not bed_data:
        print("ERROR: No BED data loaded")
        sys.exit(1)

    # Determine clusters to analyze
    clusters = sorted(assignments['cluster'].unique())
    if args.clusters:
        clusters = [int(c) for c in args.clusters.split(',')]
        print(f"\nFiltering to {len(clusters)} specified clusters")
    else:
        print(f"\nAnalyzing all {len(clusters)} clusters")

    # Summarize each cluster
    print("\nAnnotating clusters...")
    results = []

    for cluster_id in clusters:
        cluster_reads = set(assignments[assignments['cluster'] == cluster_id]['sequence'])

        if len(cluster_reads) < args.min_size:
            continue

        # Basic info
        row = {
            'cluster_id': cluster_id,
            'size': len(cluster_reads),
        }

        # Add info from cluster analysis
        if cluster_id in cluster_info:
            info = cluster_info[cluster_id]
            enrichment = info.get('enrichment', 'unknown')
            row['enrichment'] = enrichment

            # Add group-specific counts and percentages dynamically
            for grp in group_names:
                count_val = info.get(f'{grp}_count')
                pct_val = info.get(f'{grp}_pct')
                row[f'{grp}_count'] = count_val
                row[f'{grp}_pct'] = round(pct_val, 1) if pct_val is not None else None

            # Handle per-sample mode: extract enriched sample's stats
            if per_sample_mode and enrichment and enrichment not in ('mixed', 'unknown'):
                # Extract sample name from enrichment (e.g., "BJ-enriched" -> "BJ")
                enriched_sample = enrichment.replace('-enriched', '')

                # Get the enriched sample's p-value and odds ratio
                sample_pval = info.get(f'{enriched_sample}_pval')
                sample_odds = info.get(f'{enriched_sample}_odds')

                # Use sample-specific values
                row['enriched_sample'] = enriched_sample
                row['enriched_pval'] = f"{sample_pval:.4e}" if sample_pval is not None else None

                # Calculate log2_fc from sample-specific odds ratio
                row['log2_fc'] = _calculate_log2_fc(sample_odds)
                row['odds_ratio'] = round(sample_odds, 2) if sample_odds is not None else None
            else:
                # Fallback to global q_value and odds_ratio (for group comparison mode)
                q = info.get('q_value')
                row['q_value'] = f"{q:.4e}" if q is not None else None
                odds = info.get('odds_ratio')
                row['odds_ratio'] = round(odds, 2) if odds is not None else None
                row['log2_fc'] = _calculate_log2_fc(odds)

        # Annotate each featureset
        for fs in featuresets:
            if fs in bed_data:
                row[f'{fs}_top'] = summarize_featureset(cluster_reads, bed_data[fs], args.top_n, exclude_patterns)

        results.append(row)

    # Create output DataFrame
    result_df = pd.DataFrame(results)

    # Sort by enrichment group, then by log2_fc (descending) within each group
    if 'enrichment' in result_df.columns and 'log2_fc' in result_df.columns:
        # Sort by enrichment first (alphabetically), then by log2_fc descending
        result_df = result_df.sort_values(
            ['enrichment', 'log2_fc'],
            ascending=[True, False],
            na_position='last'
        )
    elif 'log2_fc' in result_df.columns:
        # Sort so most enriched in first group are at top (highest log2_fc)
        result_df = result_df.sort_values('log2_fc', ascending=False, na_position='last')
    elif group_names:
        first_grp_pct = f'{group_names[0]}_pct'
        if first_grp_pct in result_df.columns:
            result_df = result_df.sort_values(first_grp_pct, ascending=False)

    # Save output
    result_df.to_csv(args.output, sep='\t', index=False)
    print(f"\nSaved cluster annotations to: {args.output}")

    # Print summary
    print(f"\n{'=' * 60}")
    print("Summary")
    print("=" * 60)
    print(f"Clusters annotated: {len(result_df)}")

    if 'enrichment' in result_df.columns:
        print("\nBy enrichment:")
        for enrich in result_df['enrichment'].unique():
            count = (result_df['enrichment'] == enrich).sum()
            print(f"  {enrich}: {count}")

    # Show 100% clusters for each group
    for grp in group_names:
        pct_col = f'{grp}_pct'
        if pct_col in result_df.columns:
            n_100pct = (result_df[pct_col] == 100).sum()
            if n_100pct > 0:
                print(f"\n100% {grp} clusters: {n_100pct}")

    # Show fold change statistics for per-sample mode
    if per_sample_mode and 'log2_fc' in result_df.columns:
        print("\nFold change statistics (log2):")
        fc_data = result_df[result_df['log2_fc'].notna()]['log2_fc']
        if len(fc_data) > 0:
            print(f"  Mean: {fc_data.mean():.2f}")
            print(f"  Median: {fc_data.median():.2f}")
            print(f"  Range: [{fc_data.min():.2f}, {fc_data.max():.2f}]")

        # Show per-enrichment group stats
        print("\nBy enrichment group:")
        for enrich in sorted(result_df['enrichment'].unique()):
            subset = result_df[result_df['enrichment'] == enrich]
            if 'log2_fc' in subset.columns:
                fc_vals = subset['log2_fc'].dropna()
                if len(fc_vals) > 0:
                    print(f"  {enrich}: n={len(subset)}, median log2FC={fc_vals.median():.2f}")



if __name__ == "__main__":
    main()
