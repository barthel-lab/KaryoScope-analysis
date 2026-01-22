#!/usr/bin/env python3
"""
KaryoScope Representative Read Selection

Selects the best representative reads for each cluster based on:
1. Feature matching (reads must contain cluster's defining features)
2. Read length (prefer longer reads for visualization)

Outputs a TSV file that can be used with --reads-file in KaryoScope_cluster_plot.py

Usage:
  python KaryoScope_select_representatives.py \
    --cluster-analysis tmp/IDH_astro.cluster_analysis.tsv \
    --read-assignments tmp/IDH_astro.read_assignments.tsv \
    --cluster-labels tmp/IDH_astro.cluster_annotations-curated.xlsx \
    --bed-prefix results \
    --n-per-cluster 5 \
    --output tmp/IDH_astro.representative_reads.tsv
"""

import argparse
import gzip
import os
import re
from collections import defaultdict

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Select representative reads for cluster visualization")
    parser.add_argument("--cluster-analysis", required=True, help="Path to cluster_analysis.tsv")
    parser.add_argument("--read-assignments", required=True, help="Path to read_assignments.tsv")
    parser.add_argument("--cluster-labels", required=True, help="Path to cluster labels file with _top columns (TSV or Excel)")
    parser.add_argument("--bed-prefix", required=True, help="Base directory for BED files")
    parser.add_argument("--database", default="KS_human_CHM13", help="Database name")
    parser.add_argument("--smoothness", default="smoothed", help="Smoothness level")
    parser.add_argument("--n-per-cluster", type=int, default=5, help="Number of representatives per cluster")
    parser.add_argument("--feature-min-pct", type=float, default=10.0, help="Minimum feature %% to be characteristic")
    parser.add_argument("--clusters", help="Comma-separated cluster IDs to process (default: all)")
    parser.add_argument("--preferred-min-length", type=int, default=20000, help="Minimum preferred read length in bp (default: 20000)")
    parser.add_argument("--preferred-max-length", type=int, default=30000, help="Maximum preferred read length in bp (default: 30000)")
    parser.add_argument("--output", required=True, help="Output TSV file")
    return parser.parse_args()


def parse_top_features(cluster_row, min_pct=10.0):
    """Parse features with >= min_pct% coverage from _top columns."""
    features = []
    for col in ['region_top', 'subtelomeric_top', 'repeat_top']:
        top_string = cluster_row.get(col, '')
        if pd.isna(top_string) or not top_string:
            continue
        featureset = col.replace('_top', '')
        for part in str(top_string).split(';'):
            match = re.match(r'\s*(\w+)\s*\(([\d.]+)%\)', part.strip())
            if match:
                feat_name, pct = match.group(1), float(match.group(2))
                if pct >= min_pct:
                    features.append((feat_name, pct, featureset))
    return features


def load_all_features(samples, read_ids, bed_prefix, featuresets, smoothness, database):
    """Load features for all reads from all samples into memory.

    Returns:
        dict: {(read_id, sample): {featureset: {feature: bp}}}
    """
    # Build set of reads we need
    reads_needed = set(read_ids)

    # Cache: (read_id, sample) -> {featureset: {feature: bp}}
    feature_cache = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    print(f"  Pre-loading features for {len(reads_needed)} reads from {len(samples)} samples...")

    for sample in samples:
        for featureset in featuresets:
            bed_patterns = [
                f"{bed_prefix}/{sample}/telogator/1/KaryoScope/{database}/{sample}.telogator.1.{database}.{featureset}.{smoothness}.features.bed.gz",
                f"{bed_prefix}/{sample}/telogator/1/KaryoScope/{database}/{sample}.telogator.1.{database}.{featureset}.{smoothness}.features.bed",
            ]
            for bed_path in bed_patterns:
                if os.path.exists(bed_path):
                    try:
                        open_func = gzip.open if bed_path.endswith('.gz') else open
                        mode = 'rt' if bed_path.endswith('.gz') else 'r'
                        with open_func(bed_path, mode) as f:
                            for line in f:
                                parts = line.strip().split('\t')
                                if len(parts) >= 4:
                                    read_id = parts[0]
                                    if read_id in reads_needed:
                                        start, end = int(parts[1]), int(parts[2])
                                        feature = parts[3]
                                        feature_cache[(read_id, sample)][featureset][feature] += (end - start)
                    except Exception:
                        pass
                    break

    print(f"  Loaded features for {len(feature_cache)} read-sample pairs")
    return feature_cache


def check_read_features_cached(read_id, sample, cluster_features, feature_cache):
    """Check if a read has the cluster's defining features (using cache)."""
    if not cluster_features:
        return False, 0

    read_features = feature_cache.get((read_id, sample), {})

    # Check for top feature
    top_feat, top_pct, top_fs = cluster_features[0]
    has_top = read_features.get(top_fs, {}).get(top_feat, 0) > 0

    # Count matching features
    matches = sum(1 for f, p, fs in cluster_features if read_features.get(fs, {}).get(f, 0) > 0)

    return has_top, matches


def select_representatives(cluster_data, cluster_features, n_reps, feature_cache,
                           preferred_min_length=20000, preferred_max_length=30000,
                           length_column='read_span'):
    """Select best representative reads for a cluster.

    Priority: reads in preferred length range that have the defining features.
    For reads outside preferred range, prefer those closest to the range boundaries.
    Uses pre-loaded feature cache for fast lookup.

    Args:
        cluster_data: DataFrame with cluster reads
        cluster_features: List of (feature, pct, featureset) tuples
        n_reps: Number of representatives to select
        feature_cache: Pre-loaded feature cache
        preferred_min_length: Minimum preferred read length (default: 20000)
        preferred_max_length: Maximum preferred read length (default: 30000)
        length_column: Column to use for length filtering (default: 'read_span')
    """
    # Fall back to read_length if read_span not available
    if length_column not in cluster_data.columns:
        length_column = 'read_length'

    preferred_mid = (preferred_min_length + preferred_max_length) / 2

    # For preferred range: sort by closeness to midpoint
    in_range = cluster_data[(cluster_data[length_column] >= preferred_min_length) &
                            (cluster_data[length_column] <= preferred_max_length)].copy()
    in_range['_dist'] = abs(in_range[length_column] - preferred_mid)
    in_range = in_range.sort_values('_dist')

    # For below range: sort by closeness to min (prefer longer)
    below_range = cluster_data[(cluster_data[length_column] >= 10000) &
                               (cluster_data[length_column] < preferred_min_length)].copy()
    below_range['_dist'] = preferred_min_length - below_range[length_column]
    below_range = below_range.sort_values('_dist')

    # For above range: sort by closeness to max (prefer shorter)
    above_range = cluster_data[cluster_data[length_column] > preferred_max_length].copy()
    above_range['_dist'] = above_range[length_column] - preferred_max_length
    above_range = above_range.sort_values('_dist')

    # Very short reads last
    very_short = cluster_data[cluster_data[length_column] < 10000].copy()
    very_short = very_short.sort_values(length_column, ascending=False)

    # Tier order: preferred range, then closest to range boundaries
    tiers = [in_range, below_range, above_range, very_short]

    selected = []

    for tier_df in tiers:
        if len(selected) >= n_reps:
            break
        for _, row in tier_df.iterrows():
            if len(selected) >= n_reps:
                break

            has_top, n_matches = check_read_features_cached(
                row['read'], row['sample'], cluster_features, feature_cache
            )

            if has_top or n_matches > 0:
                selected.append({
                    'read': row['read'],
                    'sample': row['sample'],
                    'read_length': row.get('read_length', row.get(length_column, 0)),
                    'read_span': row.get('read_span', row.get(length_column, 0)),
                    'centroid_distance': row['centroid_distance'],
                    'has_top_feature': has_top,
                    'n_feature_matches': n_matches
                })

    # If not enough found, fall back to reads closest to preferred range
    if len(selected) < n_reps:
        # Combine all tiers and try without feature requirement
        all_reads = pd.concat(tiers, ignore_index=True)
        for _, row in all_reads.iterrows():
            if len(selected) >= n_reps:
                break
            if row['read'] not in [s['read'] for s in selected]:
                selected.append({
                    'read': row['read'],
                    'sample': row['sample'],
                    'read_length': row.get('read_length', row.get(length_column, 0)),
                    'read_span': row.get('read_span', row.get(length_column, 0)),
                    'centroid_distance': row['centroid_distance'],
                    'has_top_feature': False,
                    'n_feature_matches': 0
                })

    return selected


def normalize_by_rank(cluster_reps, n_per_cluster):
    """Reorder reads so that index N reads have similar lengths across clusters.

    For each rank position (1, 2, 3...), assign the read whose length best fits
    the median length for that position across all clusters.
    Uses read_span (full coordinate range) for length comparisons.
    """
    # First, compute target lengths for each rank
    # Use the median length at each rank across all clusters
    # Prefer read_span (full length) over read_length (filtered)
    all_lengths = []
    for cluster_id, reps in cluster_reps.items():
        lengths = sorted([r.get('read_span', r.get('read_length', 0)) for r in reps], reverse=True)
        all_lengths.append(lengths)

    # Compute median length for each rank position
    target_lengths = []
    for rank in range(n_per_cluster):
        lengths_at_rank = [ls[rank] if rank < len(ls) else 0 for ls in all_lengths]
        lengths_at_rank = [l for l in lengths_at_rank if l > 0]
        if lengths_at_rank:
            target_lengths.append(sorted(lengths_at_rank)[len(lengths_at_rank) // 2])
        else:
            target_lengths.append(0)

    print(f"  Target lengths by rank: {[f'{l:,}bp' for l in target_lengths]}")

    # For each cluster, assign reads to ranks based on closest length match
    normalized = {}
    for cluster_id, reps in cluster_reps.items():
        available = list(reps)
        assigned = [None] * n_per_cluster

        for rank in range(n_per_cluster):
            if not available:
                break
            target = target_lengths[rank]
            # Find read with length closest to target (using read_span)
            best_idx = min(range(len(available)),
                          key=lambda i: abs(available[i].get('read_span', available[i].get('read_length', 0)) - target))
            assigned[rank] = available.pop(best_idx)

        normalized[cluster_id] = [r for r in assigned if r is not None]

    return normalized


def main():
    args = parse_args()

    print("Loading data...")

    # Load read assignments
    reads_df = pd.read_csv(args.read_assignments, sep='\t')
    print(f"  Loaded {len(reads_df)} reads from {args.read_assignments}")

    # Load cluster labels (for feature info)
    if args.cluster_labels.endswith('.xlsx') or args.cluster_labels.endswith('.xls'):
        labels_df = pd.read_excel(args.cluster_labels)
    else:
        labels_df = pd.read_csv(args.cluster_labels, sep='\t')
    print(f"  Loaded {len(labels_df)} clusters from {args.cluster_labels}")

    # Check for _top columns
    top_cols = [c for c in labels_df.columns if c.endswith('_top')]
    if not top_cols:
        print("  ERROR: No _top columns found in cluster_labels file")
        return
    print(f"  Found {len(top_cols)} feature columns: {top_cols}")

    # Determine clusters to process
    if args.clusters:
        cluster_ids = [int(c.strip()) for c in args.clusters.split(',')]
    else:
        cluster_ids = sorted(reads_df['cluster'].unique())
    print(f"  Processing {len(cluster_ids)} clusters")

    # Filter reads to only those in clusters we're processing
    cluster_reads_df = reads_df[reads_df['cluster'].isin(cluster_ids)]
    print(f"  Total reads in selected clusters: {len(cluster_reads_df)}")

    # Pre-load all features into memory (fast approach)
    print("\nPre-loading feature data...")
    samples = cluster_reads_df['sample'].unique()
    read_ids = cluster_reads_df['read'].unique()
    featuresets = ['region', 'subtelomeric', 'repeat']

    feature_cache = load_all_features(
        samples, read_ids, args.bed_prefix, featuresets,
        args.smoothness, args.database
    )

    # Select representatives for each cluster
    print("\nSelecting representatives...")
    cluster_reps = {}

    for i, cluster_id in enumerate(cluster_ids):
        cluster_data = reads_df[reads_df['cluster'] == cluster_id]
        if cluster_data.empty:
            continue

        # Get cluster features
        cluster_row = labels_df[labels_df['cluster_id'] == cluster_id]
        if cluster_row.empty:
            cluster_features = []
        else:
            cluster_features = parse_top_features(cluster_row.iloc[0].to_dict(), args.feature_min_pct)

        reps = select_representatives(
            cluster_data, cluster_features, args.n_per_cluster, feature_cache,
            preferred_min_length=args.preferred_min_length,
            preferred_max_length=args.preferred_max_length
        )

        cluster_reps[cluster_id] = reps

        # Progress
        if (i + 1) % 10 == 0 or i == len(cluster_ids) - 1:
            print(f"  Processed {i + 1}/{len(cluster_ids)} clusters")

    # Normalize by rank (similar lengths at same index)
    print("\nNormalizing read lengths by rank...")
    normalized_reps = normalize_by_rank(cluster_reps, args.n_per_cluster)

    # Build output
    print("\nWriting output...")
    rows = []
    for cluster_id in cluster_ids:
        if cluster_id not in normalized_reps:
            continue
        for rank, rep in enumerate(normalized_reps[cluster_id], 1):
            rows.append({
                'cluster_id': cluster_id,
                'rank': rank,
                'read': rep['read'],
                'sample': rep['sample'],
                'read_length': rep.get('read_length', 0),
                'read_span': rep.get('read_span', rep.get('read_length', 0)),
                'centroid_distance': rep['centroid_distance'],
                'has_top_feature': rep['has_top_feature'],
                'n_feature_matches': rep['n_feature_matches']
            })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(args.output, sep='\t', index=False)
    print(f"  Saved {len(out_df)} representatives to {args.output}")

    # Also write a simple reads file for --reads-file
    reads_file = args.output.replace('.tsv', '.reads.txt')
    with open(reads_file, 'w') as f:
        for row in rows:
            f.write(row['read'] + '\n')
    print(f"  Saved read IDs to {reads_file}")

    # Summary stats
    print("\nSummary:")
    print(f"  Clusters: {len(normalized_reps)}")
    print(f"  Total reads: {len(rows)}")

    has_features = sum(1 for r in rows if r['has_top_feature'])
    print(f"  Reads with top feature: {has_features} ({has_features/len(rows)*100:.1f}%)")

    # Use read_span for length statistics (actual read length)
    spans = [r['read_span'] for r in rows]
    print(f"  Read span range: {min(spans):,} - {max(spans):,} bp")


if __name__ == '__main__':
    main()
