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
  - {prefix}.read_assignments.tsv
  - {prefix}.cluster_analysis.tsv
"""

import argparse
import fnmatch
import gzip
import os
import re
import sys
from collections import Counter, defaultdict

import pandas as pd


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
        return '', 0

    # Count by feature (weighted by length)
    feature_bp = cluster_bed.groupby('feature')['length'].sum()

    # Filter out excluded features
    if exclude_patterns:
        feature_bp = feature_bp[~feature_bp.index.map(lambda f: matches_any_pattern(f, exclude_patterns))]

    if len(feature_bp) == 0:
        return '', 0

    total_bp = feature_bp.sum()

    # Top features by coverage
    top_features = feature_bp.nlargest(top_n)
    top_str = '; '.join([f"{f}({100*v/total_bp:.1f}%)" for f, v in top_features.items()])

    return top_str, len(feature_bp)


def categorize_feature(feature):
    """Categorize a feature into high-level categories."""
    feature_lower = feature.lower()

    # Satellites
    if 'hsat' in feature_lower:
        return 'satellite', 'HSat'
    if 'asat' in feature_lower or 'a-sat' in feature_lower or 'alpha' in feature_lower or feature_lower in ['active', 'inactive', 'divergent', 'monomeric']:
        return 'satellite', 'aSat'
    if 'bsat' in feature_lower:
        return 'satellite', 'BSat'
    if 'gsat' in feature_lower:
        return 'satellite', 'GSat'

    # rDNA
    if 'rdna' in feature_lower or 'rrna' in feature_lower:
        return 'rDNA', 'rDNA'

    # Chromosome
    if feature.startswith('chr') and '_' not in feature:
        return 'chromosome', feature

    return 'other', feature


def get_satellite_summary(cluster_reads, bed_df):
    """Get satellite and rDNA summary for a cluster."""
    cluster_bed = bed_df[bed_df['read'].isin(cluster_reads)]

    if len(cluster_bed) == 0:
        return '', False

    feature_bp = cluster_bed.groupby('feature')['length'].sum()
    total_bp = feature_bp.sum()

    satellites = Counter()
    has_rDNA = False

    for feature, bp in feature_bp.items():
        cat, subcat = categorize_feature(feature)
        if cat == 'satellite':
            satellites[subcat] += bp
        elif cat == 'rDNA':
            has_rDNA = True

    sat_str = '; '.join([f"{s}({100*v/total_bp:.1f}%)" for s, v in satellites.most_common() if v/total_bp > 0.001]) if satellites else ''

    return sat_str, has_rDNA


def main():
    parser = argparse.ArgumentParser(
        description="Annotate clusters with dominant features per featureset",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--prefix", required=True,
                        help="Analysis prefix (auto-finds {prefix}.read_assignments.tsv, {prefix}.cluster_analysis.tsv)")
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
                        default="*multigroup*,*_arm,nonsubtelomeric,nonacrocentric,nonrepeat,categorized",
                        help="Comma-separated features to exclude, supports wildcards (default: '*multigroup*,*_arm,nonsubtelomeric,nonacrocentric,nonrepeat,categorized')")

    args = parser.parse_args()

    print("=" * 60)
    print("KaryoScope Cluster Annotation")
    print("=" * 60)

    # Derive file paths from prefix
    read_assignments_file = f"{args.prefix}.read_assignments.tsv"
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
    if os.path.exists(cluster_analysis_file):
        print(f"\nLoading cluster analysis: {cluster_analysis_file}")
        ca = pd.read_csv(cluster_analysis_file, sep='\t')
        for _, row in ca.iterrows():
            cluster_info[row['cluster_id']] = {
                'enrichment': row.get('enrichment', 'unknown'),
                'p_value': row.get('p_value', None),
                'q_value': row.get('q_value', None),
                'Tumor_count': row.get('Tumor_count', None),
                'Tumor_pct': row.get('Tumor_pct', None),
                'Normal_count': row.get('Normal_count', None),
                'Normal_pct': row.get('Normal_pct', None),
            }
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
        cluster_reads = set(assignments[assignments['cluster'] == cluster_id]['read'].tolist())

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
            row['enrichment'] = info.get('enrichment', 'unknown')
            row['Tumor_count'] = info.get('Tumor_count')
            row['Tumor_pct'] = round(info.get('Tumor_pct'), 1) if info.get('Tumor_pct') is not None else None
            row['Normal_count'] = info.get('Normal_count')
            row['Normal_pct'] = round(info.get('Normal_pct'), 1) if info.get('Normal_pct') is not None else None
            row['q_value'] = info.get('q_value')

        # Annotate each featureset
        for fs in featuresets:
            if fs in bed_data:
                top_str, _ = summarize_featureset(cluster_reads, bed_data[fs], args.top_n, exclude_patterns)
                row[f'{fs}_top'] = top_str

        # Get satellite summary (from any available featureset)
        for fs in ['region', 'centromeric', 'repeat']:
            if fs in bed_data:
                sat_str, _ = get_satellite_summary(cluster_reads, bed_data[fs])
                row['satellites'] = sat_str
                break

        results.append(row)

    # Create output DataFrame
    result_df = pd.DataFrame(results)

    # Sort by tumor percentage (descending)
    if 'Tumor_pct' in result_df.columns:
        result_df = result_df.sort_values('Tumor_pct', ascending=False)

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

    if 'Tumor_pct' in result_df.columns:
        n_100pct = (result_df['Tumor_pct'] == 100).sum()
        print(f"\n100% Tumor clusters: {n_100pct}")

    # Print clusters with satellites
    if 'satellites' in result_df.columns:
        sat_clusters = result_df[result_df['satellites'].notna() & (result_df['satellites'] != '')]
        if len(sat_clusters) > 0:
            print(f"\nClusters with satellite signal: {len(sat_clusters)}")
            for _, row in sat_clusters.head(5).iterrows():
                print(f"  Cluster {row['cluster_id']}: {row.get('satellites', '')}")


if __name__ == "__main__":
    main()
