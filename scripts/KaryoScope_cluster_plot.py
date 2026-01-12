#!/usr/bin/env python3
"""
KaryoScope Cluster Representative Plotting

Plots representative reads from each cluster with sample and cluster annotations.
Designed to work with outputs from KaryoScope_cluster_analysis.py.

Usage with pre-selected representative reads (recommended):
  # Step 1: Select representative reads using KaryoScope_select_representatives.py
  python KaryoScope_select_representatives.py \
    --cluster-analysis tmp/analysis.cluster_analysis.tsv \
    --read-assignments tmp/analysis.read_assignments.tsv \
    --cluster-labels tmp/analysis.cluster_annotations.xlsx \
    --bed-prefix results \
    --n-per-cluster 5 \
    --output tmp/analysis.representative_reads.tsv

  # Step 2: Plot with pre-selected reads
  python KaryoScope_cluster_plot.py \
    --cluster-analysis-prefix tmp/analysis \
    --input-bed-prefix results \
    --colors resources/KS_human_CHM13 \
    --featuresets repeat,region \
    --reads-file tmp/analysis.representative_reads.reads.txt \
    --output cluster_representatives.svg

Usage with auto-discovery (all reads):
  python KaryoScope_cluster_plot.py \
    --cluster-analysis-prefix tmp/NHA_repeat_region_composite \
    --input-bed-prefix results \
    --database KS_human_CHM13 \
    --colors resources/KS_human_CHM13 \
    --featuresets repeat,region \
    --smoothness smoothed \
    --output cluster_representatives.svg
"""

import argparse
import glob
import gzip
import os
import sys
from collections import defaultdict, OrderedDict
from math import floor

# Capture original command line for logging
_original_command = ' '.join(sys.argv)

import drawsvg as draw
import matplotlib
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd


# =============================================================================
# Command Line Arguments
# =============================================================================

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate KaryoScope SVG for cluster representative reads.",
        formatter_class=argparse.RawTextHelpFormatter)

    # Input/Output
    parser.add_argument("--cluster-analysis-prefix", dest="cluster_prefix", required=True,
                        help="Prefix from cluster_analysis.py outputs (auto-discovers .read_assignments.tsv, .feature_matrix.npz, etc.)")
    parser.add_argument("--output", required=True,
                        help="Output SVG file path")

    # Data sources
    parser.add_argument("--bed", dest="bed_files", nargs='+',
                        help="Full paths to BED files. If not provided, uses --input-bed-prefix to auto-discover.")
    parser.add_argument("--input-bed-prefix", dest="input_bed_prefix",
                        help="Base directory for auto-discovery of BED files (e.g., 'results'). "
                             "Structure: {prefix}/{sample}/telogator/1/KaryoScope/{database}/")
    parser.add_argument("--database", dest="database",
                        help="Database name (e.g., KS_human_CHM13). Auto-detected from --bed paths if not provided.")
    parser.add_argument("--colors", dest="colors_dir", required=True,
                        help="Full path to colors database directory (contains {database}.{featureset}.colors.txt files)")
    parser.add_argument("--featuresets", default="chromosome,subtelomeric,region",
                        help="Comma-separated list of feature sets to plot (default: chromosome,subtelomeric,region)")
    parser.add_argument("--custom-beds", dest="custom_beds", nargs='+', metavar="NAME:PATH",
                        help="Custom BED files to add as feature tracks. Format: 'featureset_name:/path/to/file.bed'. "
                             "These are added after standard featuresets. Requires matching color file in --colors dir.")
    parser.add_argument("--density-featuresets", dest="density_featuresets", default=None,
                        help="Comma-separated list of featuresets to render as density tracks (e.g., 'fiberseq_m6A,fiberseq_5mC'). "
                             "Small features in these tracks are binned and colored by density level.")
    parser.add_argument("--density-bin-size", dest="density_bin_size", type=int, default=300,
                        help="Bin size in bp for density computation (default: 300)")
    parser.add_argument("--density-line-plot", dest="density_line_plot", default=None,
                        help="Combine multiple featuresets into one track as overlaid density line plots. "
                             "Format: 'fiberseq_m6A:fiberseq_5mC'. Uses --density-bin-size for binning. "
                             "Lines are colored by the first color in each featureset's color file.")
    parser.add_argument("--rect-plot", dest="rect_plot", default=None,
                        help="Combine multiple featuresets into one track as stacked rectangles (exact calls). "
                             "Format: 'fiberseq_FIRE:fiberseq_LINKER'. Unlike density plots, this shows "
                             "the exact feature regions as colored rectangles.")

    # Fiberseq-specific options
    parser.add_argument("--fiberseq", dest="fiberseq_dir", default=None,
                        help="Directory containing fiberseq BED files. Auto-discovers files matching "
                             "*.FIRE.bed, *.LINKER.bed, *.m6A.bed, *.5mC.bed patterns. "
                             "Sets up FIRE_LINKER as a combined feature track and m6A/5mC as density lines.")

    # Display options
    parser.add_argument("--background", dest="background_color", default="black",
                        choices=["white", "black"],
                        help="Background color for the SVG (default: black)")
    parser.add_argument("--bar-width", dest="bar_width", type=int, default=8,
                        help="Width of each feature bar in pixels (default: 8)")
    parser.add_argument("--bar-spacing", dest="bar_spacing", type=int, default=0,
                        help="Spacing between bars within a read group (default: 0)")
    parser.add_argument("--read-spacing", dest="read_spacing", type=int, default=12,
                        help="Spacing between read groups (default: 12)")
    parser.add_argument("--cluster-spacing", dest="cluster_spacing", type=int, default=30,
                        help="Spacing between clusters (default: 30)")
    parser.add_argument("--ratio", type=float, default=1/300,
                        help="Ratio for scaling bp to pixels (default: 1/300)")
    parser.add_argument("--smoothness", default="smoothed",
                        help="Smoothness level (default: smoothed)")

    # Read selection (filtering is done by KaryoScope_select_representatives.py)
    parser.add_argument("--reads-file", dest="reads_file", default=None,
                        help="File containing read names to include (one per line). "
                             "Only these reads will be plotted. Use KaryoScope_select_representatives.py "
                             "to pre-select representative reads.")

    # Mode options
    parser.add_argument("--show-dendrogram", dest="show_dendrogram", action="store_true",
                        help="Show hierarchical clustering dendrogram to the left of feature plots")
    parser.add_argument("--hide-brackets", dest="hide_brackets", action="store_true",
                        help="Hide cluster brackets and labels (cleaner dendrogram view)")
    parser.add_argument("--no-reorder", dest="no_reorder", action="store_true",
                        help="Disable dendrogram reordering - keep reads grouped by cluster")
    parser.add_argument("--hide-dendrogram", dest="hide_dendrogram", action="store_true",
                        help="Completely hide the dendrogram (sets dendrogram height to 0)")
    parser.add_argument("--full-dendrogram", dest="full_dendrogram", action="store_true",
                        help="Show complete hierarchical tree down to individual reads/taxa, "
                             "instead of cluster-level dendrogram. Computes linkage from adjacency matrix.")
    parser.add_argument("--target-width", dest="target_width", type=int, default=None,
                        help="Target image width in pixels (auto-calculates ratio to fit)")
    parser.add_argument("--target-height", dest="target_height", type=int, default=None,
                        help="Target image height in pixels (auto-calculates ratio to fit)")
    parser.add_argument("--dendro-cluster-gap", dest="dendro_cluster_gap", type=int, default=0,
                        help="Extra gap (pixels) between cluster groups in full dendrogram mode (default: 0)")
    parser.add_argument("--cluster-labels", dest="cluster_labels", default=None,
                        help="TSV or Excel file with custom cluster labels. "
                             "Must have 'cluster_id' and label column (default column: 'curated_annotation')")
    parser.add_argument("--label-column", dest="label_column", default="curated_annotation",
                        help="Column name for custom labels in --cluster-labels file (default: curated_annotation)")
    parser.add_argument("--vertical", dest="vertical", action="store_true",
                        help="Rotate plot 90 degrees (dendrogram on left, reads vertical)")
    parser.add_argument("--show-matrix", dest="show_matrix", action="store_true",
                        help="Show sample × cluster read count matrix (vertical mode only)")
    parser.add_argument("--column-tracks", dest="column_tracks", action="store_true",
                        help="Display featuresets as separate columns instead of stacked rows. "
                             "In vertical mode: each featureset gets its own column area. "
                             "In horizontal mode: each featureset gets its own row area.")
    parser.add_argument("--n-per-cluster", dest="max_reps", type=int, default=None,
                        help="Maximum number of sequences per cluster (optional, for fallback selection)")
    parser.add_argument("--curated-reps", dest="curated_reps", default=None,
                        help="TSV file with curated representative selection. Must have 'cluster_id' and "
                             "'curated_rep_i' columns. curated_rep_i indicates which rank (1-based) to plot "
                             "for each cluster. If not specified, plots rank 1 for each cluster.")
    parser.add_argument("--show-read-indices", dest="show_read_indices", action="store_true",
                        help="Show read index labels (1, 2, 3, ...) next to each read (default: hidden)")
    parser.add_argument("--show-threshold", dest="show_threshold", action="store_true",
                        help="Visualize the structural distance threshold on the dendrogram")
    parser.add_argument("--structural-threshold", "--st", dest="structural_threshold", type=float, default=0.25,
                        help="Threshold for structural outlier clustering (default: 0.25)")
    parser.add_argument("--priority-samples", dest="priority_samples", default=None,
                        help="Comma-separated list of sample names to prioritize as representatives for clusters.")
    parser.add_argument("--log-file", dest="log_file",
                        action=argparse.BooleanOptionalAction, default=True,
                        help="Save console output to {output}.log (default: True)")

    return parser.parse_args()


# =============================================================================
# Helper Functions: Data Loading
# =============================================================================

def load_sample_metadata(metadata_file):
    """Load sample metadata from TSV file.

    Returns:
        tuple: (sample_to_group, sample_colors, group_colors)
    """
    sample_to_group = {}
    sample_colors = {}
    group_colors = {}

    if metadata_file and os.path.exists(metadata_file):
        try:
            meta_df = pd.read_csv(metadata_file, sep='\t')
            for _, row in meta_df.iterrows():
                sample = row['sample']
                group = row.get('group', sample)
                sample_to_group[sample] = group
                if 'color' in meta_df.columns and pd.notna(row.get('color')):
                    sample_colors[sample] = row['color']
                    # Also store as group color
                    if group not in group_colors:
                        group_colors[group] = row['color']
            print(f"  Loaded sample metadata: {len(meta_df)} samples")
        except Exception as e:
            print(f"  Warning: Could not load sample metadata: {e}")

    return sample_to_group, sample_colors, group_colors


def load_cluster_labels(labels_file, label_column="curated_annotation"):
    """Load custom cluster labels from TSV or Excel file.

    Args:
        labels_file: Path to TSV or Excel file with cluster_id and label columns
        label_column: Column name containing the labels

    Returns:
        dict: cluster_id -> label mapping
    """
    labels = {}

    if labels_file is None or not os.path.exists(labels_file):
        return labels

    try:
        if labels_file.endswith('.xlsx') or labels_file.endswith('.xls'):
            df = pd.read_excel(labels_file)
        else:
            df = pd.read_csv(labels_file, sep='\t')

        if 'cluster_id' in df.columns and label_column in df.columns:
            for _, row in df.iterrows():
                if pd.notna(row[label_column]) and str(row[label_column]).strip():
                    labels[int(row['cluster_id'])] = str(row[label_column])
            print(f"  Loaded {len(labels)} custom cluster labels from {labels_file}")
        else:
            print(f"  Warning: Could not find 'cluster_id' or '{label_column}' columns in {labels_file}")
    except Exception as e:
        print(f"  Warning: Could not load cluster labels: {e}")

    return labels


def load_cluster_analysis(cluster_analysis_file):
    """Load cluster analysis results to get enrichment info and cluster order.

    Returns:
        tuple: (cluster_enrichments dict, cluster_order list, cluster_stats dict, cluster_df DataFrame)
            - cluster_enrichments: cluster_id -> enrichment label
            - cluster_order: list of cluster_ids sorted by enrichment tier then p-value:
                Tier 0: 100% enriched (perfect)
                Tier 1: 80%+ enriched (strong)
                Tier 2: all others
            - cluster_stats: cluster_id -> {'odds_ratio': float, 'size': int, 'q_value': float}
            - cluster_df: DataFrame with full cluster analysis (for feature info)
    """
    cluster_enrichments = {}
    cluster_order = []
    cluster_stats = {}
    cluster_df = None

    if cluster_analysis_file and os.path.exists(cluster_analysis_file):
        try:
            df = pd.read_csv(cluster_analysis_file, sep='\t')
            cluster_df = df.copy()  # Keep original for feature info

            # Find percentage columns (they end with _pct)
            pct_cols = [c for c in df.columns if c.endswith('_pct')]

            # Determine enrichment tier for each cluster
            def get_enrichment_tier(row):
                if not pct_cols:
                    return 2  # No pct columns, treat as "other"
                max_pct = max(row[col] for col in pct_cols)
                if max_pct == 100.0:
                    return 0  # Perfect: 100% enriched
                elif max_pct >= 80.0:
                    return 1  # Strong: 80%+ enriched
                else:
                    return 2  # Other

            df['enrichment_tier'] = df.apply(get_enrichment_tier, axis=1)

            # Sort by: enrichment tier (lower = better), then by p-value
            df = df.sort_values(['enrichment_tier', 'p_value'], ascending=[True, True])

            for _, row in df.iterrows():
                cluster_id = row['cluster_id']
                cluster_enrichments[cluster_id] = row['enrichment']
                cluster_order.append(cluster_id)
                # Store stats for bubble plots
                cluster_stats[cluster_id] = {
                    'odds_ratio': row.get('odds_ratio', 1.0),
                    'size': row.get('size', 0),
                    'q_value': row.get('q_value', 1.0)
                }
            print(f"  Loaded cluster analysis: {len(df)} clusters")
        except Exception as e:
            print(f"  Warning: Could not load cluster analysis: {e}")

    return cluster_enrichments, cluster_order, cluster_stats, cluster_df


# =============================================================================
# Helper Functions: Representative Read Selection Strategy
# =============================================================================

def parse_top_features(cluster_row, min_pct=10.0):
    """Parse top features from cluster analysis row.

    Considers: region_top, subtelomeric_top, repeat_top (in priority order)
    Only includes features with >= min_pct% coverage.

    Args:
        cluster_row: Row from cluster_analysis.tsv (dict-like)
        min_pct: Minimum percentage threshold for a feature to be considered characteristic

    Returns:
        list of (feature_name, coverage_pct, featureset) tuples
    """
    import re
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


def get_read_features(read_id, sample, bed_prefix, featuresets=['region', 'subtelomeric', 'repeat'], smoothness='smoothed', database='KS_human_CHM13'):
    """Load features for a single read from BED files.

    Args:
        read_id: Read identifier
        sample: Sample name
        bed_prefix: Base directory for BED files
        featuresets: List of featuresets to check
        smoothness: 'smoothed' or 'presmoothed'
        database: Database name

    Returns:
        dict: {featureset: {feature_name: total_bp}} mapping feature names to total base pairs
    """
    read_features = {}
    for featureset in featuresets:
        read_features[featureset] = defaultdict(int)
        # Try both .gz and uncompressed
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
                            if len(parts) >= 3 and parts[0] == read_id:
                                start, end = int(parts[1]), int(parts[2])
                                feature = parts[3] if len(parts) > 3 else 'unknown'
                                read_features[featureset][feature] += (end - start)
                except Exception:
                    pass
                break
    return read_features


def score_read_features(read_id, sample, cluster_features, bed_prefix, smoothness='smoothed', database='KS_human_CHM13'):
    """Score how well a read matches the cluster's characteristic features.

    Args:
        read_id: Read identifier
        sample: Sample name
        cluster_features: List of (feature_name, pct, featureset) from parse_top_features()
        bed_prefix: Base directory for BED files
        smoothness: 'smoothed' or 'presmoothed'
        database: Database name

    Returns:
        dict: {
            'feature_score': fraction of top cluster features present in read (0-1),
            'top_feature_bp': bp coverage of the #1 cluster feature,
            'has_top_feature': bool
        }
    """
    if not cluster_features:
        return {'feature_score': 0.0, 'top_feature_bp': 0, 'has_top_feature': False}

    # Get unique featuresets to check
    featuresets = list(set(f[2] for f in cluster_features))
    read_features = get_read_features(read_id, sample, bed_prefix, featuresets, smoothness, database)

    # Count how many cluster features are present in the read
    features_found = 0
    top_feature_bp = 0

    for i, (feat_name, pct, featureset) in enumerate(cluster_features):
        bp = read_features.get(featureset, {}).get(feat_name, 0)
        if bp > 0:
            features_found += 1
            if i == 0:  # Top feature
                top_feature_bp = bp

    feature_score = features_found / len(cluster_features) if cluster_features else 0.0
    has_top_feature = top_feature_bp > 0

    return {
        'feature_score': feature_score,
        'top_feature_bp': top_feature_bp,
        'has_top_feature': has_top_feature
    }


def compute_balanced_score(feature_score, read_length, centroid_distance, max_length, max_distance):
    """Compute composite score for balanced representative selection.

    Weights (per user specification):
      - feature_match: 0.5 (must contain cluster features at >10% coverage)
      - length_norm: 0.4 (prefer longer reads for visualization)
      - distance_inv: 0.1 (minor consideration for centroid proximity)

    Args:
        feature_score: Fraction of top cluster features present (0-1)
        read_length: Read length in bp
        centroid_distance: Distance from cluster centroid
        max_length: Maximum read length in cluster (for normalization)
        max_distance: Maximum centroid distance in cluster (for normalization)

    Returns:
        float: Composite score (higher is better)
    """
    length_norm = read_length / max_length if max_length > 0 else 0
    distance_inv = 1 - (centroid_distance / max_distance) if max_distance > 0 else 0

    return 0.5 * feature_score + 0.4 * length_norm + 0.1 * distance_inv


def select_fallback_read(cluster_df):
    """Fallback selection: longest read among those closest to centroid.

    Used when no reads contain the cluster's characteristic features.

    Args:
        cluster_df: DataFrame of reads in the cluster

    Returns:
        Series: The selected read row
    """
    # Get top 25% by centroid distance, then select longest
    n_candidates = max(1, len(cluster_df) // 4)
    candidates = cluster_df.nsmallest(n_candidates, 'centroid_distance')
    return candidates.loc[candidates['read_length'].idxmax()]


def _select_by_centroid(cluster_data, max_reps):
    """Original selection: by rank (centroid distance) with proportional sampling.

    Args:
        cluster_data: DataFrame of reads in the cluster
        max_reps: Maximum number of representatives to select

    Returns:
        list of (read_id, sample) tuples
    """
    sample_counts = cluster_data['sample'].value_counts()
    total_reads = len(cluster_data)

    selected_reads = []
    remaining_slots = max_reps

    for sample in sample_counts.index:
        proportion = sample_counts[sample] / total_reads
        n_select = max(1, round(proportion * max_reps))
        n_select = min(n_select, remaining_slots, sample_counts[sample])

        if n_select > 0:
            sample_reads = cluster_data[cluster_data['sample'] == sample].sort_values('rank')
            selected = list(zip(sample_reads['read'].iloc[:n_select],
                               sample_reads['sample'].iloc[:n_select]))
            selected_reads.extend(selected)
            remaining_slots -= n_select

        if remaining_slots <= 0:
            break

    return selected_reads


def _select_by_strategy(cluster_data, cluster_features, max_reps, strategy, bed_prefix, smoothness, database):
    """Select representative reads using feature-based scoring.

    Fast approach: prioritize long reads (16-20kb), check them one by one until
    we find one that has the cluster's defining features. Stops early once enough
    matches are found.

    Args:
        cluster_data: DataFrame of reads in the cluster
        cluster_features: List of (feature_name, pct, featureset) tuples
        max_reps: Maximum number of representatives to select
        strategy: 'feature-match', 'longest', or 'balanced'
        bed_prefix: Base directory for BED files
        smoothness: 'smoothed' or 'presmoothed'
        database: Database name

    Returns:
        list of (read_id, sample) tuples, or empty list if no reads match features
    """
    if not cluster_features:
        # No features to match, return empty to trigger fallback
        return []

    # Sort reads by length (descending) - prefer longer reads
    sorted_data = cluster_data.sort_values('read_length', ascending=False)

    # Prioritize reads by length tiers (check longest first)
    # Tier 1: 16-20kb (ideal visualization length)
    # Tier 2: 10-16kb (good length)
    # Tier 3: <10kb (shorter reads)
    ideal_reads = sorted_data[(sorted_data['read_length'] >= 16000) & (sorted_data['read_length'] <= 20000)]
    good_reads = sorted_data[(sorted_data['read_length'] >= 10000) & (sorted_data['read_length'] < 16000)]
    other_reads = sorted_data[sorted_data['read_length'] < 10000]

    selected = []

    # Check reads in priority order until we find enough matches
    for candidate_df in [ideal_reads, good_reads, other_reads]:
        if len(selected) >= max_reps:
            break

        for _, row in candidate_df.iterrows():
            if len(selected) >= max_reps:
                break

            read_id = row['read']
            sample = row['sample']

            # Check if this read has the cluster's defining features
            score_info = score_read_features(
                read_id, sample, cluster_features,
                bed_prefix, smoothness, database
            )

            # Accept if read has the top feature or any matching features
            if score_info['has_top_feature'] or score_info['feature_score'] > 0:
                selected.append((read_id, sample))

    return selected


def _select_fallback(cluster_data, max_reps):
    """Fallback selection: longest reads among those closest to centroid.

    Args:
        cluster_data: DataFrame of reads in the cluster
        max_reps: Maximum number of representatives to select

    Returns:
        list of (read_id, sample) tuples
    """
    # Get top 25% by centroid distance
    n_candidates = max(1, len(cluster_data) // 4)
    candidates = cluster_data.nsmallest(n_candidates, 'centroid_distance')

    # Sort by length (descending) and select top max_reps
    candidates = candidates.sort_values('read_length', ascending=False)
    selected = candidates.head(max_reps)

    return list(zip(selected['read'], selected['sample']))


def load_curated_representatives(curated_reps_file, cluster_labels_file=None):
    """Load curated representative selection from TSV file.

    Args:
        curated_reps_file: Path to TSV file with cluster_id, rank, read columns.
        cluster_labels_file: Optional path to TSV/Excel file with cluster_id and curated_rep_i columns.
                            If provided, curated_rep_i is read from here instead of curated_reps_file.

    Returns:
        set: Read IDs to include (one per cluster based on curated_rep_i or rank=1)
    """
    df = pd.read_csv(curated_reps_file, sep='\t')

    # Check required columns
    if 'cluster_id' not in df.columns or 'rank' not in df.columns or 'read' not in df.columns:
        print(f"  Warning: curated_reps file missing required columns (cluster_id, rank, read)")
        return set()

    # Load curated_rep_i from cluster_labels_file if provided
    curated_rep_map = {}
    if cluster_labels_file and os.path.exists(cluster_labels_file):
        try:
            if cluster_labels_file.endswith('.xlsx') or cluster_labels_file.endswith('.xls'):
                labels_df = pd.read_excel(cluster_labels_file)
            else:
                labels_df = pd.read_csv(cluster_labels_file, sep='\t')

            if 'cluster_id' in labels_df.columns and 'curated_rep_i' in labels_df.columns:
                for _, row in labels_df.iterrows():
                    if pd.notna(row.get('curated_rep_i')):
                        curated_rep_map[row['cluster_id']] = int(row['curated_rep_i'])
                print(f"  Loaded {len(curated_rep_map)} curated_rep_i values from {cluster_labels_file}")
        except Exception as e:
            print(f"  Warning: Could not load curated_rep_i from cluster_labels: {e}")

    # Fallback: check if curated_rep_i exists in curated_reps_file itself
    if not curated_rep_map and 'curated_rep_i' in df.columns:
        for cluster_id in df['cluster_id'].unique():
            cluster_df = df[df['cluster_id'] == cluster_id]
            curated_vals = cluster_df['curated_rep_i'].dropna()
            if len(curated_vals) > 0:
                curated_rep_map[cluster_id] = int(curated_vals.iloc[0])

    selected_reads = set()

    for cluster_id in df['cluster_id'].unique():
        cluster_df = df[df['cluster_id'] == cluster_id]
        target_rank = curated_rep_map.get(cluster_id, 1)

        # Select read with matching rank
        target_row = cluster_df[cluster_df['rank'] == target_rank]
        if len(target_row) > 0:
            selected_reads.add(target_row.iloc[0]['read'])
        elif len(cluster_df) > 0:
            # Fallback to first available read
            selected_reads.add(cluster_df.iloc[0]['read'])

    print(f"  Selected {len(selected_reads)} curated representatives from {curated_reps_file}")

    return selected_reads


def load_representative_reads(reps_file, cluster_enrichments=None, cluster_order=None, max_reps=None, reads_file=None, curated_reps_file=None, cluster_labels_file=None):
    """Load read assignments from TSV file.

    Representative read selection should be done beforehand using KaryoScope_select_representatives.py.
    This function loads the selected reads and groups them by cluster.

    Args:
        reps_file: Path to read_assignments.tsv (all reads with cluster assignments and stats)
        cluster_enrichments: Dict of cluster_id -> enrichment label from cluster_analysis.tsv
        cluster_order: List of cluster_ids in priority order (for ordering output)
        max_reps: Maximum representatives per cluster (optional fallback if no reads_file)
        reads_file: Path to file with read names to include (one per line)
        curated_reps_file: Path to TSV with cluster_id, rank, read columns for selecting specific reads
        cluster_labels_file: Path to TSV/Excel with curated_rep_i column (optional)

    Returns:
        tuple: (cluster_reads OrderedDict, unique_enrichments set)
    """
    print(f"\nLoading read assignments from: {reps_file}")
    reps_df = pd.read_csv(reps_file, sep='\t')
    print(f"  Total reads: {len(reps_df)}")

    # Load rank info from curated_reps_file or infer from reads_file
    read_ranks = {}
    rank_source = None
    if curated_reps_file and os.path.exists(curated_reps_file):
        rank_source = curated_reps_file
    elif reads_file:
        # Check if there's a .tsv version with rank info
        tsv_path = reads_file.replace('.reads.txt', '.tsv')
        if os.path.exists(tsv_path):
            rank_source = tsv_path

    if rank_source:
        try:
            rank_df = pd.read_csv(rank_source, sep='\t')
            if 'read' in rank_df.columns and 'rank' in rank_df.columns:
                for _, row in rank_df.iterrows():
                    read_ranks[row['read']] = row['rank']
        except Exception:
            pass  # Silently ignore if rank loading fails

    # Load curated representatives if provided (takes precedence over reads_file)
    if curated_reps_file and os.path.exists(curated_reps_file):
        allowed_reads = load_curated_representatives(curated_reps_file, cluster_labels_file)
        if allowed_reads:
            reps_df = reps_df[reps_df['read'].isin(allowed_reads)]
            print(f"  After curated filter: {len(reps_df)} reads")
    # Otherwise filter by reads file if provided
    elif reads_file:
        if not os.path.exists(reads_file):
            print(f"  Warning: reads_file not found: {reads_file}")
        else:
            try:
                with open(reads_file, 'r') as f:
                    allowed_reads = set(line.strip() for line in f if line.strip())
                if allowed_reads:
                    reps_df = reps_df[reps_df['read'].isin(allowed_reads)]
                    print(f"  After reads filter: {len(reps_df)} reads (from {len(allowed_reads)} in file)")
                else:
                    print(f"  Warning: reads_file is empty, no filtering applied")
            except Exception as e:
                print(f"  Warning: Could not read reads_file: {e}")

    # Merge enrichment info from cluster_analysis.tsv
    if cluster_enrichments:
        reps_df['enrichment'] = reps_df['cluster'].map(cluster_enrichments)
        reps_df['enrichment'] = reps_df['enrichment'].fillna('unknown')
    else:
        # Fallback: use group as enrichment proxy
        reps_df['enrichment'] = reps_df['group'].apply(lambda x: f"{x}-enriched" if pd.notna(x) else 'unknown')

    # Get unique enrichment labels from data
    unique_enrichments = set(reps_df['enrichment'].unique())
    print(f"  Available enrichment categories: {sorted(unique_enrichments)}")

    # Get unique clusters - use priority order from cluster_analysis.tsv if available
    available_clusters = set(reps_df['cluster'].unique())
    if cluster_order:
        clusters_to_plot = [c for c in cluster_order if c in available_clusters]
    else:
        clusters_to_plot = list(reps_df['cluster'].unique())
    print(f"  Clusters to plot: {len(clusters_to_plot)}")

    # Group reads by cluster
    cluster_reads = OrderedDict()

    for cluster_id in clusters_to_plot:
        cluster_data = reps_df[reps_df['cluster'] == cluster_id].copy()
        if cluster_data.empty:
            continue
        enrichment = cluster_data['enrichment'].iloc[0]

        # Optionally limit to max_reps per cluster (simple centroid-based fallback)
        if max_reps is not None and len(cluster_data) > max_reps:
            selected_reads = _select_by_centroid(cluster_data, max_reps)
        else:
            selected_reads = list(zip(cluster_data['read'], cluster_data['sample']))

        # Sort by rank if rank info is available
        if read_ranks:
            selected_reads = sorted(selected_reads, key=lambda x: read_ranks.get(x[0], 999))

        cluster_reads[cluster_id] = {
            'enrichment': enrichment,
            'reads': selected_reads
        }

    if max_reps is not None:
        print(f"  Limited to max {max_reps} representatives per cluster")

    return cluster_reads, unique_enrichments


def load_feature_matrix(matrix_file):
    """Load feature matrix from NPZ file.

    Returns:
        dict or None: Feature matrix data
    """
    if matrix_file and os.path.exists(matrix_file):
        try:
            data = np.load(matrix_file, allow_pickle=True)
            print(f"Loaded feature matrix from: {matrix_file}")
            return data
        except Exception as e:
            print(f"Warning: Could not load feature matrix: {e}")
    return None


def parse_bed_paths(bed_files):
    """Parse BED file paths to extract sample names, directories, and database.

    Args:
        bed_files: List of full paths to BED files

    Returns:
        tuple: (sample_bed_paths dict, database name)
            - sample_bed_paths: sample_name -> bed_directory
            - database: extracted database name (e.g., KS_human_CHM13)
    """
    print(f"\nParsing BED file paths...")
    sample_bed_paths = {}
    database = None

    for bed_path in bed_files:
        if not os.path.exists(bed_path):
            sys.stderr.write(f"Error: BED file not found: {bed_path}\n")
            sys.exit(1)

        bed_dir = os.path.dirname(bed_path)
        filename = os.path.basename(bed_path)

        # Expected format: {sample}.telogator.1.{database}.{featureset}.{smoothness}.KaryoScope.bed
        # or: {sample}.telogator.1.{database}.{featureset}.{smoothness}.features.bed.gz
        parts = filename.split('.')
        if len(parts) >= 4:
            sample_name = parts[0]
            # Database is typically the 4th part (index 3)
            db_name = parts[3]

            if sample_name not in sample_bed_paths:
                sample_bed_paths[sample_name] = bed_dir
                print(f"  {sample_name} -> {bed_dir}")

            # Use first database found, verify consistency
            if database is None:
                database = db_name
            elif database != db_name:
                sys.stderr.write(f"Warning: Inconsistent database names: {database} vs {db_name}\n")

    if database:
        print(f"  Database: {database}")

    return sample_bed_paths, database


def load_color_files(colors_dir, database, featuresets):
    """Load color mappings for featuresets.

    Color files contain features with _specific suffix (e.g., active_specific).
    For smoothed BED files, features don't have the suffix (e.g., active).
    This function creates mappings for both versions.

    Returns:
        tuple: (featureset_colors, featureset_color_order)
    """
    print(f"\nLoading color files...")
    featureset_colors = {}
    featureset_color_order = {}

    for fs in featuresets:
        colors_pattern = f"{database}.{fs}.colors.txt"
        colors_path = os.path.join(colors_dir, colors_pattern)

        if not os.path.exists(colors_path):
            sys.stderr.write(f"Error: Colors file not found: {colors_path}\n")
            sys.exit(1)

        featureset_colors[fs] = {0: ("#ffffff", 1.0)}
        featureset_color_order[fs] = []

        with open(colors_path, "r") as f:
            for i, line in enumerate(f):
                parts = line.strip().split()
                if len(parts) >= 2:
                    feature = parts[0]
                    color = parts[1]
                    # Skip header line (feature/color)
                    if i == 0 and feature.lower() == 'feature':
                        continue
                    featureset_colors[fs][feature] = (color, 1.0)
                    featureset_color_order[fs].append(feature)

                    # Also add mapping without _specific suffix for smoothed BED files
                    if feature.endswith('_specific'):
                        base_feature = feature[:-9]  # Remove '_specific'
                        featureset_colors[fs][base_feature] = (color, 1.0)

        print(f"  {fs}: {len(featureset_color_order[fs])} colors")

    return featureset_colors, featureset_color_order


def load_bed_data(sample_bed_paths, database, featuresets, smoothness, reads_needed):
    """Load BED data for specified reads.

    Returns:
        dict: read_data[read][featureset] = list of features
    """
    print(f"\nLoading BED data for representative reads...")
    print(f"  Reads to load: {len(reads_needed)}")

    read_data = defaultdict(lambda: defaultdict(list))

    for fs in featuresets:
        for sample_name, bed_dir in sample_bed_paths.items():
            bed_pattern = f"{sample_name}.telogator.1.{database}.{fs}.{smoothness}.features.bed.gz"
            bed_path = os.path.join(bed_dir, bed_pattern)
            if not os.path.exists(bed_path):
                # Restore lenient matching for manual inputs or non-standard names
                # Pattern looks for {sample} and {featureset} in the same directory
                import glob
                candidates = glob.glob(os.path.join(bed_dir, f"*{sample_name}*{fs}*.bed*"))
                if candidates:
                    # Filter for database and smoothness if possible to be more specific
                    best_candidates = [c for c in candidates if database in c and smoothness in c]
                    bed_path = best_candidates[0] if best_candidates else candidates[0]
                else:
                    continue

            open_func = gzip.open if bed_path.endswith(".gz") else open
            mode = "rt" if bed_path.endswith(".gz") else "r"

            with open_func(bed_path, mode) as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) < 4:
                        continue  # Skip malformed BED lines
                    scaffold, start, stop, feature = parts[:4]
                    try:
                        start, stop = int(start), int(stop)
                    except ValueError:
                        continue  # Skip lines with non-integer coordinates

                    if scaffold in reads_needed:
                        read_data[scaffold][fs].append({
                            'start': start,
                            'stop': stop,
                            'feature': feature
                        })

    print(f"  Loaded data for {len(read_data)} reads")
    return read_data


def load_custom_bed_files(custom_bed_files, reads_needed, read_data=None):
    """Load custom BED files into read_data structure.

    Args:
        custom_bed_files: Dict of featureset_name -> path
        reads_needed: Set of read names to load
        read_data: Existing read_data dict to update (or None to create new)

    Returns:
        dict: read_data[read][featureset] = list of features
    """
    if read_data is None:
        read_data = defaultdict(lambda: defaultdict(list))

    if not custom_bed_files:
        return read_data

    print(f"\nLoading custom BED files...")
    for fs_name, bed_path in custom_bed_files.items():
        print(f"  {fs_name}: {bed_path}")
        try:
            with open(bed_path) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) < 4:
                        continue
                    read_name, start, end, feature = parts[0], int(parts[1]), int(parts[2]), parts[3]
                    if read_name not in reads_needed:
                        continue
                    if read_name not in read_data:
                        read_data[read_name] = {}
                    if fs_name not in read_data[read_name]:
                        read_data[read_name][fs_name] = []
                    read_data[read_name][fs_name].append({
                        'start': start,
                        'stop': end,
                        'feature': feature
                    })
        except Exception as e:
            print(f"  Error loading {bed_path}: {e}")

    return read_data


# =============================================================================
# Helper Functions: Enrichment Handling
# =============================================================================

def get_enrichment_colors(group_colors, unique_enrichments, sample_colors=None):
    """Generate enrichment colors based on group or sample colors.

    Args:
        group_colors: Dict of group -> color from metadata
        unique_enrichments: Set of enrichment labels from data
        sample_colors: Dict of sample -> color (for per-sample mode)

    Returns:
        dict: enrichment_label -> color
    """
    enrichment_colors = {'mixed': '#999999'}

    for enrich in unique_enrichments:
        if enrich == 'mixed':
            continue

        # Extract group/sample name from enrichment label (e.g., "post-enriched" -> "post")
        name = enrich.replace('-enriched', '')

        # Try to find matching color (case-insensitive)
        color_found = False

        # First check sample_colors (for per-sample mode)
        if sample_colors:
            for s, c in sample_colors.items():
                if s.lower() == name.lower():
                    enrichment_colors[enrich] = c
                    color_found = True
                    break

        # Then check group_colors
        if not color_found:
            for g, c in group_colors.items():
                if g.lower() == name.lower():
                    enrichment_colors[enrich] = c
                    color_found = True
                    break

        # Default colors if not found
        if not color_found:
            if 'post' in enrich.lower():
                enrichment_colors[enrich] = '#E41A1C'
            elif 'pre' in enrich.lower():
                enrichment_colors[enrich] = '#377EB8'
            else:
                enrichment_colors[enrich] = '#666666'

    return enrichment_colors


def get_cluster_colors(unique_clusters):
    """Generate colors for clusters.

    Returns:
        dict: cluster_id -> color
    """
    cluster_cmap = matplotlib.colormaps.get_cmap('tab20')
    cluster_colors = {}

    for i, cid in enumerate(sorted(unique_clusters)):
        color_idx = i % 20
        cluster_colors[cid] = mcolors.rgb2hex(cluster_cmap(color_idx))

    return cluster_colors


def get_primary_color(fs, featureset_colors, featureset_color_order, default="#FFFFFF"):
    """Get the first non-unannotated color for a featureset.

    Args:
        fs: Featureset name
        featureset_colors: Dict of featureset -> feature -> (color, opacity)
        featureset_color_order: Dict of featureset -> list of feature names
        default: Default color if none found

    Returns:
        str: Hex color code
    """
    if fs in featureset_colors and featureset_colors[fs]:
        for feat_name in featureset_color_order.get(fs, []):
            if feat_name != 'unannotated':
                return featureset_colors[fs][feat_name][0]
    return default


def generate_sample_colors(samples, existing_colors=None):
    """Generate colors for samples.

    Args:
        samples: List of sample names
        existing_colors: Existing sample -> color mapping

    Returns:
        dict: sample -> color
    """
    if existing_colors is None:
        existing_colors = {}

    sample_colors = existing_colors.copy()
    samples_needing_colors = [s for s in samples if s not in sample_colors]

    if samples_needing_colors:
        n_samples = len(samples_needing_colors)
        tab10 = matplotlib.colormaps.get_cmap('tab10')

        for i, sample in enumerate(sorted(samples_needing_colors)):
            if n_samples == 2:
                sample_colors[sample] = '#377EB8' if i == 0 else '#E41A1C'
            else:
                sample_colors[sample] = mcolors.rgb2hex(tab10(i % 10))

    return sample_colors


# =============================================================================
# Helper Functions: Drawing Components
# =============================================================================

def compute_cluster_dendrogram_order(feature_matrix_data, cluster_reads):
    """Order clusters using pre-computed cluster-level linkage from cluster_analysis.py.

    Uses the cluster-level dendrogram (based on cluster centroid distances) to order
    clusters, keeping all reads within each cluster grouped together.

    Returns:
        tuple: (reordered_cluster_reads, cluster_dendro_data, read_to_original_cluster, read_to_original_enrichment)
               cluster_dendro_data contains 'linkage' and 'cluster_order' for drawing cluster dendrogram
    """
    from scipy.cluster.hierarchy import leaves_list, optimal_leaf_ordering
    from scipy.spatial.distance import squareform

    read_to_original_cluster = {}
    read_to_original_enrichment = {}

    # Build read mappings
    for cluster_id, data in cluster_reads.items():
        for read, sample in data['reads']:
            read_to_original_cluster[read] = cluster_id
            read_to_original_enrichment[read] = data['enrichment']

    # Check if cluster-level linkage is available
    if 'cluster_linkage' not in feature_matrix_data:
        print("  No cluster-level linkage found in feature matrix")
        return cluster_reads, None, read_to_original_cluster, read_to_original_enrichment

    try:
        cluster_linkage = feature_matrix_data['cluster_linkage']
        cluster_ids_ordered = list(feature_matrix_data['cluster_ids_ordered'])

        # Convert to int if needed
        cluster_ids_ordered = [int(c) for c in cluster_ids_ordered]

        # Get displayed cluster IDs
        displayed_cluster_ids = set(cluster_reads.keys())

        # Find which clusters from the original linkage are being displayed
        # Map original cluster indices to displayed ones
        original_to_displayed_idx = {}
        displayed_clusters_in_order = []
        for i, cid in enumerate(cluster_ids_ordered):
            if cid in displayed_cluster_ids:
                original_to_displayed_idx[i] = len(displayed_clusters_in_order)
                displayed_clusters_in_order.append(cid)

        if len(displayed_clusters_in_order) <= 1:
            print(f"  Only {len(displayed_clusters_in_order)} cluster(s) to display, skipping dendrogram")
            return cluster_reads, None, read_to_original_cluster, read_to_original_enrichment

        # Extract subset linkage for displayed clusters
        # We need to recompute linkage on the subset of cluster centroids
        cluster_centroids = feature_matrix_data['cluster_centroids']

        # Get indices of displayed clusters in original centroid array
        displayed_centroid_indices = [cluster_ids_ordered.index(cid) for cid in displayed_clusters_in_order]
        subset_centroids = cluster_centroids[displayed_centroid_indices]

        # Use Bio.Phylo to extract subtree from the full dendrogram
        # This preserves the original tree structure for displayed clusters
        from scipy.cluster.hierarchy import to_tree
        from io import StringIO
        from Bio import Phylo

        def parse_terminal_name(name):
            """Parse cluster ID from terminal name with validation."""
            if name is None:
                raise ValueError("Terminal name is None")
            if not name.startswith('n'):
                raise ValueError(f"Invalid terminal name format: {name}")
            try:
                return int(name[1:])
            except ValueError:
                raise ValueError(f"Cannot parse cluster ID from: {name}")

        def linkage_to_newick(linkage_matrix, labels):
            """Convert scipy linkage matrix to Newick format with proper branch lengths."""
            tree = to_tree(linkage_matrix)

            def to_newick(node):
                if node.is_leaf():
                    return f"n{labels[node.id]}"  # Prefix with 'n' to ensure valid names
                else:
                    left_child = node.get_left()
                    right_child = node.get_right()
                    # Branch length = parent height - child height
                    # For leaves, child height is 0
                    left_height = 0 if left_child.is_leaf() else left_child.dist
                    right_height = 0 if right_child.is_leaf() else right_child.dist
                    left_branch = node.dist - left_height
                    right_branch = node.dist - right_height
                    left = to_newick(left_child)
                    right = to_newick(right_child)
                    return f'({left}:{left_branch:.6f},{right}:{right_branch:.6f})'

            return to_newick(tree) + ';'

        def phylo_tree_to_linkage(tree, label_to_idx):
            """Convert Bio.Phylo tree back to scipy linkage matrix format."""
            n = len(label_to_idx)
            linkage_rows = []
            node_counter = [n]  # Next internal node ID

            def get_node_info(clade):
                """Recursively build linkage. Returns (node_id, height, count)."""
                if clade.is_terminal():
                    # Leaf node - extract cluster ID from name (remove 'n' prefix)
                    cluster_id = parse_terminal_name(clade.name)
                    if cluster_id not in label_to_idx:
                        raise ValueError(f"Unknown cluster ID: {cluster_id}")
                    return label_to_idx[cluster_id], 0, 1

                # Internal node - process children
                children = clade.clades
                if len(children) == 0:
                    raise ValueError("Internal node has no children")

                if len(children) != 2:
                    # Handle non-binary nodes by sequential merging
                    left_id, left_h, left_c = get_node_info(children[0])
                    for child in children[1:]:
                        right_id, right_h, right_c = get_node_info(child)
                        height = max(left_h, right_h) + (clade.branch_length or 0.1)
                        linkage_rows.append([left_id, right_id, height, left_c + right_c])
                        left_id = node_counter[0]
                        node_counter[0] += 1
                        left_h = height
                        left_c = left_c + right_c
                    return left_id, left_h, left_c

                left_id, left_h, left_c = get_node_info(children[0])
                right_id, right_h, right_c = get_node_info(children[1])

                # Height is the max child height plus branch length
                height = max(left_h, right_h) + (clade.branch_length or 0.1)
                linkage_rows.append([left_id, right_id, height, left_c + right_c])

                new_id = node_counter[0]
                node_counter[0] += 1
                return new_id, height, left_c + right_c

            get_node_info(tree.root)
            return np.array(linkage_rows) if linkage_rows else None

        original_full_linkage = feature_matrix_data['cluster_linkage']

        if len(displayed_clusters_in_order) > 1:
            # Convert full linkage to Newick and load with Bio.Phylo
            newick_str = linkage_to_newick(original_full_linkage, cluster_ids_ordered)
            full_tree = Phylo.read(StringIO(newick_str), 'newick')

            # Get terminal names to keep (with 'n' prefix)
            terminals_to_keep = {f"n{cid}" for cid in displayed_cluster_ids}

            # Prune tree to only displayed clusters
            # Remove terminals not in our set
            all_terminals = list(full_tree.get_terminals())
            for terminal in all_terminals:
                if terminal.name not in terminals_to_keep:
                    full_tree.prune(terminal)

            # Collapse internal nodes with single children
            def collapse_single_children(clade):
                if clade.is_terminal():
                    return
                # Recursively process children first
                for child in list(clade.clades):
                    collapse_single_children(child)
                # If this node has only one child, bypass it
                while len(clade.clades) == 1:
                    child = clade.clades[0]
                    # Accumulate branch lengths (handle None values)
                    clade.branch_length = (clade.branch_length or 0) + (child.branch_length or 0)
                    clade.clades = list(child.clades)  # Explicit copy

            collapse_single_children(full_tree.root)

            # Check for degenerate tree after collapse
            if not full_tree.root.clades and not full_tree.root.is_terminal():
                print("  Warning: Tree collapsed to empty root, falling back to single cluster")
                reordered_cluster_ids = displayed_clusters_in_order
                optimized_linkage = None
            else:
                # Get leaf order from pruned tree (left-to-right traversal)
                terminals = list(full_tree.get_terminals())
                reordered_cluster_ids = []
                for t in terminals:
                    if t.name is not None:
                        reordered_cluster_ids.append(parse_terminal_name(t.name))

                print(f"  Extracted subtree with {len(reordered_cluster_ids)} clusters using Bio.Phylo")

                # Convert pruned tree back to linkage format
                label_to_idx = {cid: i for i, cid in enumerate(reordered_cluster_ids)}
                optimized_linkage = phylo_tree_to_linkage(full_tree, label_to_idx)

            # Identity mapping since we built linkage in display order
            leaf_order = list(range(len(reordered_cluster_ids)))
        else:
            reordered_cluster_ids = displayed_clusters_in_order
            optimized_linkage = None
            print(f"  Single cluster, no reordering needed")

        # Rebuild cluster_reads in dendrogram order
        cluster_reads_reordered = OrderedDict()
        for cid in reordered_cluster_ids:
            cluster_reads_reordered[cid] = cluster_reads[cid]

        # Store data for drawing cluster dendrogram
        # leaf_order maps: display_position -> original_linkage_index
        # So leaf_order[0] is the linkage index of the leftmost cluster
        cluster_dendro_data = {
            'linkage': optimized_linkage,
            'cluster_order': reordered_cluster_ids,
            'n_clusters': len(reordered_cluster_ids),
            'leaf_order': leaf_order  # Maps display position to linkage index
        }

        return cluster_reads_reordered, cluster_dendro_data, read_to_original_cluster, read_to_original_enrichment

    except Exception as e:
        print(f"  Warning: Could not compute cluster dendrogram order: {e}")
        import traceback
        traceback.print_exc()
        return cluster_reads, None, read_to_original_cluster, read_to_original_enrichment


def compute_full_dendrogram(feature_matrix_data, displayed_reads):
    """Compute full hierarchical dendrogram down to individual reads.

    Args:
        feature_matrix_data: Dictionary from feature_matrix.npz containing adj_matrix and read_names
        displayed_reads: List of read names to include in dendrogram

    Returns:
        full_dendro_data: Dictionary with linkage matrix and read ordering, or None if failed
    """
    from scipy.cluster.hierarchy import linkage, leaves_list, optimal_leaf_ordering
    from scipy.spatial.distance import pdist

    try:
        adj_matrix = feature_matrix_data['adj_matrix']
        read_names = list(feature_matrix_data['read_names'])

        # Filter to displayed reads
        displayed_set = set(displayed_reads)
        indices = [i for i, name in enumerate(read_names) if name in displayed_set]

        if len(indices) < 2:
            print(f"  Not enough reads for full dendrogram ({len(indices)} reads)")
            return None

        subset_matrix = adj_matrix[indices]
        subset_names = [read_names[i] for i in indices]

        # Compute pairwise distances and linkage
        dist = pdist(subset_matrix)
        Z = linkage(dist, method='ward')

        # Optimize leaf ordering for cleaner visualization
        try:
            Z_optimized = optimal_leaf_ordering(Z, dist)
        except:
            Z_optimized = Z

        leaf_order = leaves_list(Z_optimized)
        ordered_reads = [subset_names[i] for i in leaf_order]

        full_dendro_data = {
            'linkage': Z_optimized,
            'read_order': ordered_reads,
            'read_names': subset_names,
            'n_reads': len(subset_names),
            'leaf_order': leaf_order
        }

        print(f"  Computed full dendrogram for {len(subset_names)} reads")
        return full_dendro_data

    except Exception as e:
        print(f"  Warning: Could not compute full dendrogram: {e}")
        import traceback
        traceback.print_exc()
        return None


def draw_full_dendrogram(d, full_dendro_data, read_y_positions, read_names_displayed,
                         left_margin, dendrogram_width, background_color):
    """Draw complete hierarchical dendrogram showing all reads as leaves.

    Args:
        d: drawsvg Drawing object
        full_dendro_data: Dictionary from compute_full_dendrogram
        read_y_positions: Dict mapping read names to their Y positions
        read_names_displayed: List of read names in display order
        left_margin: X position where dendrogram ends (right edge)
        dendrogram_width: Width allocated for dendrogram
        background_color: 'white' or 'black'
    """
    if full_dendro_data is None:
        return

    linkage_matrix = full_dendro_data['linkage']
    read_names = full_dendro_data['read_names']
    n_reads = len(read_names)

    line_color = '#FFFFFF' if background_color == 'black' else '#333333'

    # Build mapping from read name to Y position
    name_to_y = {}
    for name in read_names:
        if name in read_y_positions:
            name_to_y[name] = read_y_positions[name]

    if not name_to_y:
        print("  Warning: No read positions found for full dendrogram")
        return

    # Compute node positions - leaves at their read Y positions
    # Internal nodes at average Y of children
    node_y = {}  # node_id -> y_position
    node_x = {}  # node_id -> x_position (based on merge distance)

    # Leaves (indices 0 to n-1)
    for i, name in enumerate(read_names):
        if name in name_to_y:
            node_y[i] = name_to_y[name]
            node_x[i] = left_margin  # Leaves at right edge

    # Get max distance for scaling
    max_dist = linkage_matrix[:, 2].max() if len(linkage_matrix) > 0 else 1.0

    # Process internal nodes (indices n to 2n-2)
    for i, (idx1, idx2, dist, count) in enumerate(linkage_matrix):
        idx1, idx2 = int(idx1), int(idx2)
        new_node_id = n_reads + i

        # Y position is average of children
        if idx1 in node_y and idx2 in node_y:
            node_y[new_node_id] = (node_y[idx1] + node_y[idx2]) / 2
        elif idx1 in node_y:
            node_y[new_node_id] = node_y[idx1]
        elif idx2 in node_y:
            node_y[new_node_id] = node_y[idx2]
        else:
            continue

        # X position based on distance (scaled to dendrogram width)
        # Root at left, leaves at right
        node_x[new_node_id] = left_margin - (dist / max_dist) * (dendrogram_width - 20)

    # Draw branches
    for i, (idx1, idx2, dist, count) in enumerate(linkage_matrix):
        idx1, idx2 = int(idx1), int(idx2)
        new_node_id = n_reads + i

        if new_node_id not in node_y or new_node_id not in node_x:
            continue
        if idx1 not in node_y or idx2 not in node_y:
            continue

        # Horizontal line at merge height
        x_merge = node_x[new_node_id]
        y_merge = node_y[new_node_id]

        # Vertical lines to children
        y1, y2 = node_y[idx1], node_y[idx2]
        x1 = node_x.get(idx1, left_margin)
        x2 = node_x.get(idx2, left_margin)

        # Draw: vertical from child1 to merge height, horizontal across, vertical to child2
        d.append(draw.Line(x1, y1, x_merge, y1, stroke=line_color, stroke_width=1))
        d.append(draw.Line(x_merge, y1, x_merge, y2, stroke=line_color, stroke_width=1))
        d.append(draw.Line(x_merge, y2, x2, y2, stroke=line_color, stroke_width=1))

    print(f"  Drew full dendrogram for {n_reads} reads ({len(linkage_matrix)} branches)")


def draw_full_dendrogram_header(d, full_dendro_data, read_x_positions, group_width,
                                 top_margin, dendrogram_height, background_color):
    """Draw complete hierarchical dendrogram as header showing all reads as leaves.

    Args:
        d: drawsvg Drawing object
        full_dendro_data: Dictionary from compute_full_dendrogram
        read_x_positions: Dict mapping read names to their X positions
        group_width: Width of each read column
        top_margin: Y position of dendrogram bottom (where leaves connect)
        dendrogram_height: Height allocated for dendrogram
        background_color: 'white' or 'black'
    """
    if full_dendro_data is None:
        return

    linkage_matrix = full_dendro_data['linkage']
    read_names = full_dendro_data['read_names']
    n_reads = len(read_names)

    line_color = '#FFFFFF' if background_color == 'black' else '#333333'

    # Build mapping from read name to X position (center of read column)
    name_to_x = {}
    for name in read_names:
        if name in read_x_positions:
            name_to_x[name] = read_x_positions[name] + group_width / 2

    if not name_to_x:
        print("  Warning: No read positions found for full dendrogram header")
        return

    # Compute node positions
    # Leaves at their read X positions, Y at bottom of dendrogram
    # Internal nodes at average X of children, Y based on merge distance
    node_x = {}  # node_id -> x_position
    node_y = {}  # node_id -> y_position (merge height)

    # Get max distance for scaling
    max_dist = linkage_matrix[:, 2].max() if len(linkage_matrix) > 0 else 1.0

    # Y coordinates: leaves at top_margin, root at top_margin - dendrogram_height
    y_bottom = top_margin - 10  # Leaves (bottom of dendrogram, just above feature bars)
    y_top = top_margin - dendrogram_height + 10  # Root (top of dendrogram)

    # Leaves (indices 0 to n-1)
    for i, name in enumerate(read_names):
        if name in name_to_x:
            node_x[i] = name_to_x[name]
            node_y[i] = y_bottom  # Leaves at bottom

    # Process internal nodes (indices n to 2n-2)
    for i, (idx1, idx2, dist, count) in enumerate(linkage_matrix):
        idx1, idx2 = int(idx1), int(idx2)
        new_node_id = n_reads + i

        # X position is average of children
        if idx1 in node_x and idx2 in node_x:
            node_x[new_node_id] = (node_x[idx1] + node_x[idx2]) / 2
        elif idx1 in node_x:
            node_x[new_node_id] = node_x[idx1]
        elif idx2 in node_x:
            node_x[new_node_id] = node_x[idx2]
        else:
            continue

        # Y position based on distance (scaled to dendrogram height)
        # Higher distance = higher in the tree (lower Y value)
        node_y[new_node_id] = y_bottom - (dist / max_dist) * (y_bottom - y_top)

    # Draw branches
    for i, (idx1, idx2, dist, count) in enumerate(linkage_matrix):
        idx1, idx2 = int(idx1), int(idx2)
        new_node_id = n_reads + i

        if new_node_id not in node_x or new_node_id not in node_y:
            continue
        if idx1 not in node_x or idx2 not in node_x:
            continue

        # Merge point
        x_merge = node_x[new_node_id]
        y_merge = node_y[new_node_id]

        # Child positions
        x1, x2 = node_x[idx1], node_x[idx2]
        y1 = node_y.get(idx1, y_bottom)
        y2 = node_y.get(idx2, y_bottom)

        # Draw: vertical from child1 up to merge height, horizontal across, vertical down to child2
        d.append(draw.Line(x1, y1, x1, y_merge, stroke=line_color, stroke_width=1))
        d.append(draw.Line(x1, y_merge, x2, y_merge, stroke=line_color, stroke_width=1))
        d.append(draw.Line(x2, y_merge, x2, y2, stroke=line_color, stroke_width=1))

    print(f"  Drew full dendrogram header for {n_reads} reads ({len(linkage_matrix)} branches)")


def draw_dendrogram(d, dendro_data, read_x_positions, displayed_reads,
                    group_width, top_margin, dendrogram_height, background_color):
    """Draw dendrogram manually using linkage matrix with original distances.

    The linkage matrix Z has rows [idx1, idx2, distance, count]:
    - idx1, idx2: indices of clusters being merged (< n are leaves, >= n are internal nodes)
    - distance: the merge distance (height of the horizontal line)
    - count: number of leaves in the new cluster

    We recursively traverse the tree to get proper leaf positions for each subtree,
    which ensures branches don't cross.
    """
    n_leaves = len(displayed_reads)
    linkage_matrix = dendro_data['linkage']
    leaf_order = dendro_data['leaf_order']

    # Base Y position (bottom of dendrogram, where leaves attach)
    dendro_base_y = top_margin + 20

    # Find max distance in linkage for scaling
    max_distance = linkage_matrix[:, 2].max() if len(linkage_matrix) > 0 else 1
    max_distance = max(max_distance, 1)

    def distance_to_y(dist):
        """Convert distance to Y pixel coordinate (higher distance = higher up)."""
        return dendro_base_y - (dist / max_distance) * (dendrogram_height - 15)

    line_color = '#AAAAAA' if background_color == 'black' else '#444444'

    # Build a mapping from original leaf index to display position
    # leaf_order tells us: position i in display has original leaf leaf_order[i]
    # We need: original leaf j is at display position where leaf_order[pos] == j
    original_to_display_pos = {orig_idx: pos for pos, orig_idx in enumerate(leaf_order)}

    # Get x position for each display position
    display_pos_x = {}
    for pos, read in enumerate(displayed_reads):
        if read in read_x_positions:
            display_pos_x[pos] = read_x_positions[read] + group_width / 2

    # Recursively compute the x-center and leaf range for each node
    # This ensures we get the correct x positions based on actual leaf positions
    node_x_center = {}
    node_height = {}
    node_leaves = {}  # Store which leaves are under each node

    def get_node_info(node_idx):
        """Recursively get x-center, height, and leaves for a node."""
        if node_idx in node_x_center:
            return node_x_center[node_idx], node_height[node_idx], node_leaves[node_idx]

        if node_idx < n_leaves:
            # Leaf node - use display position
            display_pos = original_to_display_pos.get(node_idx)
            if display_pos is not None and display_pos in display_pos_x:
                x = display_pos_x[display_pos]
                node_x_center[node_idx] = x
                node_height[node_idx] = 0
                node_leaves[node_idx] = [display_pos]
                return x, 0, [display_pos]
            return None, 0, []

        # Internal node - get from linkage matrix
        row_idx = node_idx - n_leaves
        if row_idx >= len(linkage_matrix):
            return None, 0, []

        idx1, idx2, dist, count = linkage_matrix[row_idx]
        idx1, idx2 = int(idx1), int(idx2)

        x1, h1, leaves1 = get_node_info(idx1)
        x2, h2, leaves2 = get_node_info(idx2)

        if x1 is None or x2 is None:
            return None, dist, []

        # Combine leaves and compute center as mean of all leaf positions
        all_leaves = leaves1 + leaves2
        x_center = sum(display_pos_x[pos] for pos in all_leaves) / len(all_leaves)

        node_x_center[node_idx] = x_center
        node_height[node_idx] = dist
        node_leaves[node_idx] = all_leaves

        return x_center, dist, all_leaves

    # First pass: compute all node info
    root_idx = n_leaves + len(linkage_matrix) - 1
    get_node_info(root_idx)

    # Second pass: draw the dendrogram
    for row_idx, (idx1, idx2, dist, count) in enumerate(linkage_matrix):
        idx1, idx2 = int(idx1), int(idx2)

        x1 = node_x_center.get(idx1)
        x2 = node_x_center.get(idx2)
        h1 = node_height.get(idx1, 0)
        h2 = node_height.get(idx2, 0)

        if x1 is None or x2 is None:
            continue

        # Y coordinates
        y_bottom_left = distance_to_y(h1)
        y_bottom_right = distance_to_y(h2)
        y_top = distance_to_y(dist)

        # Draw the U-shape: left vertical, horizontal, right vertical
        d.append(draw.Line(x1, y_bottom_left, x1, y_top, stroke=line_color, stroke_width=1.5))
        d.append(draw.Line(x1, y_top, x2, y_top, stroke=line_color, stroke_width=1.5))
        d.append(draw.Line(x2, y_bottom_right, x2, y_top, stroke=line_color, stroke_width=1.5))

    n_branches = len(linkage_matrix)
    print(f"  Drew dendrogram for {n_leaves} reads ({n_branches} branches)")


def draw_cluster_dendrogram(d, cluster_dendro_data, cluster_x_start, cluster_x_end,
                            top_margin, dendrogram_height, background_color):
    """Draw cluster-level dendrogram using linkage matrix based on cluster centroids.

    This draws a dendrogram showing how clusters relate to each other,
    with branches connecting to the center of each cluster's read group.

    Key insight: The linkage matrix indices refer to the ORIGINAL order before
    optimal_leaf_ordering. We need to map these to DISPLAY positions using leaf_order.
    leaf_order[display_pos] = original_linkage_index
    """
    if cluster_dendro_data is None or cluster_dendro_data.get('linkage') is None:
        return

    linkage_matrix = cluster_dendro_data['linkage']
    cluster_order = cluster_dendro_data['cluster_order']  # Reordered cluster IDs for display
    leaf_order = cluster_dendro_data.get('leaf_order')  # Maps display_pos -> linkage_index
    n_clusters = len(cluster_order)

    if n_clusters <= 1 or len(linkage_matrix) == 0:
        return

    # Base Y position (bottom of dendrogram, where cluster leaves attach)
    dendro_base_y = top_margin + 20

    # Find max distance in linkage for scaling
    max_distance = linkage_matrix[:, 2].max() if len(linkage_matrix) > 0 else 1
    max_distance = max(max_distance, 1)

    def distance_to_y(dist):
        """Convert distance to Y pixel coordinate (higher distance = higher up)."""
        return dendro_base_y - (dist / max_distance) * (dendrogram_height - 15)

    # Thin white lines
    line_color = '#FFFFFF'

    # Build mapping from linkage index to display x-center
    # cluster_order[i] is the cluster_id at display position i
    # leaf_order[i] is the linkage index that should be at display position i
    # So: linkage_index_to_display_pos[leaf_order[i]] = i
    linkage_idx_to_display_pos = {}
    if leaf_order is not None:
        for display_pos, linkage_idx in enumerate(leaf_order):
            linkage_idx_to_display_pos[linkage_idx] = display_pos
    else:
        # No reordering - identity mapping
        for i in range(n_clusters):
            linkage_idx_to_display_pos[i] = i

    # Compute x center for each display position
    display_pos_x_center = {}
    for display_pos, cid in enumerate(cluster_order):
        x_start = cluster_x_start.get(cid, 0)
        x_end = cluster_x_end.get(cid, 0)
        display_pos_x_center[display_pos] = (x_start + x_end) / 2

    # Recursively compute the x-center for each node based on ALL leaves underneath
    # Node indices in linkage: 0 to n-1 are leaves, n to 2n-2 are internal nodes
    node_x_center = {}
    node_height = {}
    node_display_positions = {}  # Track which display positions are under each node

    def get_node_info(node_idx):
        """Recursively get x-center, height, and display positions for a node."""
        if node_idx in node_x_center:
            return node_x_center[node_idx], node_height[node_idx], node_display_positions[node_idx]

        if node_idx < n_clusters:
            # Leaf node - map linkage index to display position
            display_pos = linkage_idx_to_display_pos.get(node_idx)
            if display_pos is None:
                return None, 0, []
            x = display_pos_x_center.get(display_pos, 0)
            node_x_center[node_idx] = x
            node_height[node_idx] = 0
            node_display_positions[node_idx] = [display_pos]
            return x, 0, [display_pos]

        # Internal node - get from linkage matrix
        row_idx = node_idx - n_clusters
        if row_idx >= len(linkage_matrix):
            return None, 0, []

        idx1, idx2, dist, count = linkage_matrix[row_idx]
        idx1, idx2 = int(idx1), int(idx2)

        x1, h1, positions1 = get_node_info(idx1)
        x2, h2, positions2 = get_node_info(idx2)

        if x1 is None or x2 is None:
            return None, dist, []

        # Combine all display positions and compute center as mean of their x-positions
        all_positions = positions1 + positions2
        x_center = sum(display_pos_x_center[pos] for pos in all_positions) / len(all_positions)

        node_x_center[node_idx] = x_center
        node_height[node_idx] = dist
        node_display_positions[node_idx] = all_positions

        return x_center, dist, all_positions

    # First pass: compute all node info
    root_idx = n_clusters + len(linkage_matrix) - 1
    get_node_info(root_idx)

    # Second pass: draw the dendrogram
    # Draw vertical lines from each child up to the merge height,
    # then a horizontal line connecting them at the merge height
    for row_idx, (idx1, idx2, dist, count) in enumerate(linkage_matrix):
        idx1, idx2 = int(idx1), int(idx2)

        x1 = node_x_center.get(idx1)
        x2 = node_x_center.get(idx2)
        h1 = node_height.get(idx1, 0)
        h2 = node_height.get(idx2, 0)

        if x1 is None or x2 is None:
            continue

        # Y coordinates
        y_bottom_left = distance_to_y(h1)
        y_bottom_right = distance_to_y(h2)
        y_top = distance_to_y(dist)

        # Draw vertical lines from each child's x-position up to merge height
        # Then horizontal line connecting at merge height
        d.append(draw.Line(x1, y_bottom_left, x1, y_top, stroke=line_color, stroke_width=1))
        d.append(draw.Line(x2, y_bottom_right, x2, y_top, stroke=line_color, stroke_width=1))
        d.append(draw.Line(x1, y_top, x2, y_top, stroke=line_color, stroke_width=1))

    n_branches = len(linkage_matrix)
    print(f"  Drew cluster dendrogram for {n_clusters} clusters ({n_branches} branches)")


def draw_cluster_brackets(d, cluster_reads, cluster_x_start, cluster_x_end,
                          enrichment_colors, read_heights, label_height, text_color,
                          cluster_labels=None):
    """Draw cluster brackets and labels below the feature bars (inverted, pointing up).

    Args:
        read_heights: Dict of read -> (min_y, max_y, x_start, total_width) for y-positioning
        label_height: Height of rotated featureset labels below bars
        cluster_labels: Optional dict of cluster_id -> custom label string
    """
    if cluster_labels is None:
        cluster_labels = {}

    for cluster_id, data in cluster_reads.items():
        if cluster_id == 'all':  # Skip when in dendrogram mode
            continue

        x_start = cluster_x_start[cluster_id]
        x_end = cluster_x_end[cluster_id]
        enrichment = data['enrichment']
        color = enrichment_colors.get(enrichment, '#999999')

        # Find the max_y of the longest read in this cluster
        cluster_max_y = 0
        for read, sample in data['reads']:
            if read in read_heights:
                _, max_y, _, _ = read_heights[read]
                cluster_max_y = max(cluster_max_y, max_y)

        # Position bracket below this cluster's feature labels
        # Labels start at max_y + 5 and extend label_height pixels down
        bracket_y = cluster_max_y + 5 + label_height

        # Draw inverted bracket (horizontal line at bottom, vertical lines pointing up)
        d.append(draw.Line(x_start, bracket_y, x_end, bracket_y, stroke=color, stroke_width=3))
        d.append(draw.Line(x_start, bracket_y, x_start, bracket_y - 8, stroke=color, stroke_width=2))
        d.append(draw.Line(x_end, bracket_y, x_end, bracket_y - 8, stroke=color, stroke_width=2))

        # Labels below the bracket
        label_x = (x_start + x_end) / 2

        # Use custom label if available, otherwise default to cluster ID
        if cluster_id in cluster_labels:
            label_text = cluster_labels[cluster_id]
            label_font_size = 9  # slightly smaller for longer labels
        else:
            label_text = f"c{cluster_id}"
            label_font_size = 10

        d.append(draw.Text(
            label_text,
            font_size=label_font_size, x=label_x, y=bracket_y + 15,
            fill=color, font_family='sans-serif',
            text_anchor='middle', font_weight='bold'
        ))

        d.append(draw.Text(
            enrichment,
            font_size=8, x=label_x, y=bracket_y + 27,
            fill=color, font_family='sans-serif', text_anchor='middle'
        ))


# =============================================================================
# Vertical Mode Drawing Functions
# =============================================================================

def draw_cluster_dendrogram_vertical(d, cluster_dendro_data, cluster_y_start, cluster_y_end,
                                     left_margin, dendrogram_width, background_color):
    """Draw cluster-level dendrogram on the LEFT side for vertical mode.

    Clusters are arranged vertically, dendrogram branches extend horizontally from left.

    Key insight: The linkage matrix indices refer to the ORIGINAL order before
    optimal_leaf_ordering. We need to map these to DISPLAY positions using leaf_order.
    """
    if cluster_dendro_data is None or cluster_dendro_data.get('linkage') is None:
        return

    linkage_matrix = cluster_dendro_data['linkage']
    cluster_order = cluster_dendro_data['cluster_order']
    leaf_order = cluster_dendro_data.get('leaf_order')
    n_clusters = len(cluster_order)

    if n_clusters <= 1 or len(linkage_matrix) == 0:
        return

    # Map cluster display position to y-center
    display_pos_y_center = {}
    for display_idx, cluster_id in enumerate(cluster_order):
        y_start = cluster_y_start.get(cluster_id, 0)
        y_end = cluster_y_end.get(cluster_id, 0)
        y_center = (y_start + y_end) / 2
        display_pos_y_center[display_idx] = y_center

    # Build mapping from linkage index to display position
    # leaf_order[i] is the linkage index at display position i
    linkage_idx_to_display_pos = {}
    if leaf_order is not None:
        for display_pos, linkage_idx in enumerate(leaf_order):
            linkage_idx_to_display_pos[linkage_idx] = display_pos
    else:
        for i in range(n_clusters):
            linkage_idx_to_display_pos[i] = i

    max_distance = linkage_matrix[:, 2].max() if len(linkage_matrix) > 0 else 1
    max_distance = max(max_distance, 1)

    # Dendrogram base X (rightmost point, where leaves attach)
    # The dendrogram occupies x=50 to x=50+dendrogram_width
    dendro_base_x = 50 + dendrogram_width

    def distance_to_x(dist):
        """Convert distance to X pixel coordinate (higher distance = further left)."""
        return dendro_base_x - (dist / max_distance) * (dendrogram_width - 15)

    # Thin white lines
    line_color = '#FFFFFF'

    # Track node positions
    node_y_center = {}
    node_height = {}
    node_display_positions = {}

    def get_node_info(node_idx):
        """Recursively get y-center, height, and display positions for a node."""
        if node_idx in node_y_center:
            return node_y_center[node_idx], node_height[node_idx], node_display_positions[node_idx]

        if node_idx < n_clusters:
            # Leaf node - map linkage index to display position
            display_pos = linkage_idx_to_display_pos.get(node_idx)
            if display_pos is None:
                return None, 0, []
            y = display_pos_y_center.get(display_pos, 0)
            node_y_center[node_idx] = y
            node_height[node_idx] = 0
            node_display_positions[node_idx] = [display_pos]
            return y, 0, [display_pos]

        # Internal node
        row_idx = node_idx - n_clusters
        if row_idx >= len(linkage_matrix):
            return None, 0, []

        idx1, idx2, dist, count = linkage_matrix[row_idx]
        idx1, idx2 = int(idx1), int(idx2)

        y1, h1, positions1 = get_node_info(idx1)
        y2, h2, positions2 = get_node_info(idx2)

        if y1 is None or y2 is None:
            return None, dist, []

        all_positions = positions1 + positions2
        # Internal node y-center should be midpoint between its two children
        y_center = (y1 + y2) / 2

        node_y_center[node_idx] = y_center
        node_height[node_idx] = dist
        node_display_positions[node_idx] = all_positions

        return y_center, dist, all_positions

    # First pass: compute all node info
    root_idx = n_clusters + len(linkage_matrix) - 1
    get_node_info(root_idx)

    # Draw the dendrogram (horizontal lines from each child, vertical line connecting)
    for row_idx, (idx1, idx2, dist, count) in enumerate(linkage_matrix):
        idx1, idx2 = int(idx1), int(idx2)

        y1 = node_y_center.get(idx1)
        y2 = node_y_center.get(idx2)
        h1 = node_height.get(idx1, 0)
        h2 = node_height.get(idx2, 0)

        if y1 is None or y2 is None:
            continue

        x_right_top = distance_to_x(h1)
        x_right_bottom = distance_to_x(h2)
        x_left = distance_to_x(dist)

        # Horizontal lines from each child to merge point
        d.append(draw.Line(x_right_top, y1, x_left, y1, stroke=line_color, stroke_width=1))
        d.append(draw.Line(x_right_bottom, y2, x_left, y2, stroke=line_color, stroke_width=1))
        # Vertical line connecting at merge point
        d.append(draw.Line(x_left, y1, x_left, y2, stroke=line_color, stroke_width=1))

    print(f"  Drew vertical dendrogram for {n_clusters} clusters")


def draw_feature_bars_vertical(d, drawing_data, featuresets, bar_width, read_positions,
                               num_featuresets, bar_spacing=2, background_color='black'):
    """Draw feature bars VERTICALLY (top to bottom) for vertical mode.

    In vertical mode:
    - Reads are stacked vertically (y varies per read)
    - Each read has horizontal feature bars extending to the right
    - Features are drawn as horizontal rectangles at y positions based on scaled coordinates
    - Multiple featuresets are stacked vertically within each read's bar area
    """
    stroke_width = 0.5

    # First, draw all feature rectangles
    for read, (x_start, y_start, bar_length) in read_positions.items():
        if read not in drawing_data:
            continue

        for fs_idx, fs in enumerate(featuresets):
            y_offset = fs_idx * (bar_width + bar_spacing)  # Stack featuresets vertically

            for feat in drawing_data[read].get(fs, []):
                x = x_start + feat['scaled_start']
                width = max(2, feat['scaled_stop'] - feat['scaled_start'])  # Minimum 2px width

                d.append(draw.Rectangle(
                    x, y_start + y_offset,
                    width, bar_width,
                    fill=feat['color'],
                    fill_opacity=feat.get('fill_opacity', 1.0)
                ))

    # Draw borders around each read's bar area
    for read, (x_start, y_start, bar_length) in read_positions.items():
        if read not in drawing_data:
            continue

        total_height = num_featuresets * bar_width + (num_featuresets - 1) * bar_spacing

        # Outer border
        d.append(draw.Rectangle(
            x_start, y_start,
            bar_length, total_height,
            fill='none',
            stroke='black',
            stroke_width=stroke_width
        ))

        # Horizontal lines between featureset bars
        for i in range(1, num_featuresets):
            line_y = y_start + i * (bar_width + bar_spacing) - bar_spacing / 2
            d.append(draw.Line(
                x_start, line_y, x_start + bar_length, line_y,
                stroke='black',
                stroke_width=stroke_width
            ))


def draw_scale_bar(d, x_start, y_pos, ratio, text_color='white'):
    """Draw a scale bar showing read length scale.

    Args:
        d: Drawing object
        x_start: X position for scale bar start
        y_pos: Y position for scale bar
        ratio: Pixels per base pair (used for scaling)
        text_color: Color for text and bar
    """
    # Choose a nice round scale bar length based on what would fit
    # Common choices: 1kb, 2kb, 5kb, 10kb, 20kb
    scale_options = [1000, 2000, 5000, 10000, 20000]

    # Find a scale bar that's reasonably sized (50-150 pixels wide)
    chosen_bp = 5000  # Default 5kb
    for bp in scale_options:
        bar_width_px = bp * ratio
        if 50 <= bar_width_px <= 150:
            chosen_bp = bp
            break

    bar_width_px = chosen_bp * ratio

    # Format label
    if chosen_bp >= 1000:
        label = f"{chosen_bp // 1000} kb"
    else:
        label = f"{chosen_bp} bp"

    # Draw scale bar line (same stroke_width=1 as dendrogram)
    d.append(draw.Line(x_start, y_pos, x_start + bar_width_px, y_pos,
                       stroke=text_color, stroke_width=1))

    # Draw end ticks
    tick_height = 4
    d.append(draw.Line(x_start, y_pos - tick_height/2, x_start, y_pos + tick_height/2,
                       stroke=text_color, stroke_width=1))
    d.append(draw.Line(x_start + bar_width_px, y_pos - tick_height/2, x_start + bar_width_px, y_pos + tick_height/2,
                       stroke=text_color, stroke_width=1))

    # Draw label centered above
    d.append(draw.Text(
        label,
        font_size=8, x=x_start + bar_width_px / 2, y=y_pos - 5,
        fill=text_color, font_family='sans-serif', text_anchor='middle'
    ))


def draw_feature_bars_column_mode(d, drawing_data, featuresets, bar_width, read_y_positions,
                                   x_start, max_bar_length, column_spacing=10, background_color='black'):
    """Draw feature bars with each featureset in its own column.

    In column mode:
    - Each featureset gets its own column
    - Within each column, reads are stacked vertically (one bar per read)
    - Features extend horizontally along the read length axis
    - Columns are separated by column_spacing

    Args:
        d: Drawing object
        drawing_data: Dict of read -> featureset -> feature list
        featuresets: List of featuresets to draw
        bar_width: Height of each feature bar
        read_y_positions: Dict of read -> y position
        x_start: Starting x position for first column
        max_bar_length: Maximum bar length (for consistent column width)
        column_spacing: Gap between columns
        background_color: Background color
    """
    stroke_width = 0.5
    num_fs = len(featuresets)

    # Draw each featureset as a separate column
    for fs_idx, fs in enumerate(featuresets):
        # Column x position
        col_x = x_start + fs_idx * (max_bar_length + column_spacing)

        # Draw features for each read in this column
        for read, y_pos in read_y_positions.items():
            if read not in drawing_data:
                continue

            for feat in drawing_data[read].get(fs, []):
                x = col_x + feat['scaled_start']
                width = max(2, feat['scaled_stop'] - feat['scaled_start'])

                d.append(draw.Rectangle(
                    x, y_pos,
                    width, bar_width,
                    fill=feat['color'],
                    fill_opacity=feat.get('fill_opacity', 1.0)
                ))

        # Draw border around this column's bars for each read
        for read, y_pos in read_y_positions.items():
            if read not in drawing_data:
                continue

            # Get read's bar length from drawing data
            read_bar_length = 0
            for fs_check in featuresets:
                for feat in drawing_data[read].get(fs_check, []):
                    read_bar_length = max(read_bar_length, feat['scaled_stop'])

            d.append(draw.Rectangle(
                col_x, y_pos,
                read_bar_length, bar_width,
                fill='none',
                stroke='black',
                stroke_width=stroke_width
            ))

    # Calculate total width used
    total_width = num_fs * (max_bar_length + column_spacing) - column_spacing
    return total_width


def draw_sample_matrix(d, cluster_ids, cluster_y_start, cluster_y_end, sample_metadata,
                       read_assignments_file, x_start, cell_width, cell_height, text_color,
                       background_color='black'):
    """Draw sample × cluster read count matrix for vertical mode.

    Uses FULL cluster read counts from read_assignments file, not just representatives.

    Args:
        d: Drawing object
        cluster_ids: List of cluster IDs to include in matrix
        cluster_y_start: Dict of cluster_id -> y start position
        cluster_y_end: Dict of cluster_id -> y end position
        sample_metadata: DataFrame with 'sample' and 'group' columns
        read_assignments_file: Path to read_assignments.tsv with all reads
        x_start: X position for matrix start
        cell_width: Width of each cell (sample column)
        cell_height: Height of each cell (matches bar_width typically)
        text_color: Color for text
        background_color: Background color ('black' or 'white')

    Returns:
        float: Total width of the matrix (for layout calculation)
    """
    # Load full read assignments to get complete cluster counts
    full_df = pd.read_csv(read_assignments_file, sep='\t')

    # Get all samples from metadata, ordered by group
    if sample_metadata is not None and 'group' in sample_metadata.columns:
        # Sort samples by group (Normal first, then Tumor/other)
        sorted_metadata = sample_metadata.sort_values(
            by='group',
            key=lambda x: x.map({'Normal': 0, 'Control': 0}).fillna(1)
        )
        all_samples = sorted_metadata['sample'].tolist()
        sample_groups = sorted_metadata.set_index('sample')['group'].to_dict()
    else:
        # Extract all unique samples from read assignments
        all_samples = sorted(full_df['sample'].unique().tolist())
        sample_groups = {}

    n_samples = len(all_samples)

    # Calculate x positions with gaps between groups
    group_gap = 5  # Gap between Normal and Tumor groups
    sample_x_positions = {}
    current_x = 0
    prev_group = None
    for sample in all_samples:
        group = sample_groups.get(sample)
        if prev_group is not None and group != prev_group:
            current_x += group_gap  # Add gap when group changes
        sample_x_positions[sample] = current_x
        current_x += cell_width
        prev_group = group

    total_matrix_width = current_x  # Total width including gaps

    # Compute FULL read counts per cluster × sample from all assignments
    cluster_sample_counts = {}
    max_count = 1  # Avoid division by zero

    for cluster_id in cluster_ids:
        cluster_df = full_df[full_df['cluster'] == cluster_id]
        counts = {sample: 0 for sample in all_samples}
        for sample in cluster_df['sample']:
            if sample in counts:
                counts[sample] += 1
        cluster_sample_counts[cluster_id] = counts
        max_count = max(max_count, max(counts.values()) if counts.values() else 1)

    # Cluster samples within each group by their count profiles
    from scipy.cluster.hierarchy import linkage, leaves_list

    # Group samples by their group
    group_to_samples = {}
    for sample in all_samples:
        group = sample_groups.get(sample, 'Unknown')
        if group not in group_to_samples:
            group_to_samples[group] = []
        group_to_samples[group].append(sample)

    # Store linkage matrices and sample orders for dendrograms
    # group_name -> (Z, original_samples, ordered_samples)
    group_linkages = {}

    # Cluster samples within each group
    clustered_samples = []
    for group_name in ['Normal', 'Control']:  # Normal/Control first
        if group_name in group_to_samples:
            group_samples = group_to_samples.pop(group_name)
            original_samples = group_samples.copy()  # Keep original order for linkage indices
            if len(group_samples) > 2:
                # Build count matrix: samples × clusters
                count_matrix = np.array([
                    [cluster_sample_counts[cid].get(s, 0) for cid in cluster_ids]
                    for s in group_samples
                ])
                # Cluster using log-transformed counts (but display absolute)
                if count_matrix.sum() > 0:  # Only cluster if there's data
                    log_matrix = np.log1p(count_matrix)  # log(1 + count)
                    Z = linkage(log_matrix, method='ward')
                    order = leaves_list(Z)
                    ordered_samples = [group_samples[i] for i in order]
                    group_linkages[group_name] = (Z, original_samples, ordered_samples)
                    group_samples = ordered_samples
            clustered_samples.extend(group_samples)

    # Add remaining groups (Tumor, etc.) - also clustered
    for group_name in sorted(group_to_samples.keys()):
        group_samples = group_to_samples[group_name]
        original_samples = group_samples.copy()
        if len(group_samples) > 2:
            count_matrix = np.array([
                [cluster_sample_counts[cid].get(s, 0) for cid in cluster_ids]
                for s in group_samples
            ])
            # Cluster using log-transformed counts (but display absolute)
            if count_matrix.sum() > 0:
                log_matrix = np.log1p(count_matrix)  # log(1 + count)
                Z = linkage(log_matrix, method='ward')
                order = leaves_list(Z)
                ordered_samples = [group_samples[i] for i in order]
                group_linkages[group_name] = (Z, original_samples, ordered_samples)
                group_samples = ordered_samples
        clustered_samples.extend(group_samples)

    all_samples = clustered_samples

    # Recalculate x positions with gaps between groups
    sample_x_positions = {}
    current_x = 0
    prev_group = None
    for sample in all_samples:
        group = sample_groups.get(sample)
        if prev_group is not None and group != prev_group:
            current_x += group_gap  # Add gap when group changes
        sample_x_positions[sample] = current_x
        current_x += cell_width
        prev_group = group

    total_matrix_width = current_x

    # Draw column headers (rotated sample names)
    header_y = min(cluster_y_start.values()) - 5 if cluster_y_start else 30
    for sample in all_samples:
        x = x_start + sample_x_positions[sample] + cell_width / 2

        # Sample name, rotated 90 degrees
        d.append(draw.Text(
            sample, font_size=7, x=x, y=header_y,
            fill=text_color, font_family='sans-serif',
            transform=f"rotate(-90 {x} {header_y})",
            text_anchor='start'
        ))

    # Draw cells for each cluster
    for cluster_id in cluster_y_start:
        if cluster_id not in cluster_sample_counts:
            continue

        y_mid = (cluster_y_start[cluster_id] + cluster_y_end[cluster_id]) / 2
        y = y_mid - cell_height / 2

        for sample in all_samples:
            x = x_start + sample_x_positions[sample]
            count = cluster_sample_counts[cluster_id].get(sample, 0)

            # Heatmap color intensity based on count (log-scaled)
            # Use 2-color sequential scale: dark -> yellow (for black background)
            # Log scaling helps visualize low values when max is high
            import math
            if count == 0:
                intensity = 0
            elif max_count > 1:
                # log1p scaling: log(1+count) / log(1+max_count)
                intensity = math.log1p(count) / math.log1p(max_count)
            else:
                intensity = 1.0

            if count == 0:
                # Zero counts: dark background
                fill_color = '#1a1a1a' if background_color == 'black' else '#f0f0f0'
            else:
                # Sequential 2-color: dark gray -> yellow
                # Interpolate from #333333 (dark) to #ffff00 (yellow)
                r = int(51 + intensity * (255 - 51))  # 51 -> 255
                g = int(51 + intensity * (255 - 51))  # 51 -> 255
                b = int(51 * (1 - intensity))          # 51 -> 0
                fill_color = f'#{r:02x}{g:02x}{b:02x}'

            # Draw cell with white border
            d.append(draw.Rectangle(
                x, y,
                cell_width, cell_height,
                fill=fill_color,
                stroke='#FFFFFF',
                stroke_width=0.5
            ))

            # Draw count text - always show the number
            if count > 0:
                # Text color contrasts with cell - white for dark cells, black for bright
                count_text_color = '#000000' if intensity > 0.4 else '#ffffff'
                font_size = 5 if cell_width >= 10 else 4
                d.append(draw.Text(
                    str(count), font_size=font_size, x=x + cell_width / 2, y=y + cell_height / 2 + 2,
                    fill=count_text_color, font_family='sans-serif',
                    text_anchor='middle'
                ))
            else:
                # Always show "0" for zero counts - lighter color for visibility
                d.append(draw.Text(
                    '0', font_size=4, x=x + cell_width / 2, y=y + cell_height / 2 + 1.5,
                    fill='#666666', font_family='sans-serif',
                    text_anchor='middle'
                ))

    # Return data needed for bar plot and dendrogram
    return {
        'width': total_matrix_width,
        'sample_x_positions': sample_x_positions,
        'cluster_sample_counts': cluster_sample_counts,
        'all_samples': all_samples,
        'sample_groups': sample_groups,
        'group_linkages': group_linkages,
        'cell_width': cell_width,
        'max_count': max_count
    }


def draw_matrix_legend(d, x_start, y_start, max_count, text_color='white', background_color='black'):
    """Draw a legend for the matrix color scale (log-scaled).

    Args:
        d: Drawing object
        x_start: X position for legend start
        y_start: Y position for legend
        max_count: Maximum count value in the matrix
        text_color: Color for text
        background_color: Background color ('black' or 'white')
    """
    import math

    legend_width = 100
    legend_height = 10
    n_steps = 20

    # Title
    d.append(draw.Text(
        "Read Count (log scale)",
        font_size=9, x=x_start, y=y_start,
        fill=text_color, font_family='sans-serif', font_weight='bold'
    ))

    bar_y = y_start + 12

    # Draw gradient bar (using same log scaling as matrix)
    step_width = legend_width / n_steps
    for i in range(n_steps):
        intensity = i / (n_steps - 1)
        # Same color calculation as matrix cells
        r = int(51 + intensity * (255 - 51))
        g = int(51 + intensity * (255 - 51))
        b = int(51 * (1 - intensity))
        fill_color = f'#{r:02x}{g:02x}{b:02x}'

        d.append(draw.Rectangle(
            x_start + i * step_width, bar_y,
            step_width + 0.5, legend_height,  # +0.5 to avoid gaps
            fill=fill_color, stroke='none'
        ))

    # Border around gradient bar
    d.append(draw.Rectangle(
        x_start, bar_y, legend_width, legend_height,
        fill='none', stroke=text_color, stroke_width=0.5
    ))

    # Labels - show a few tick marks on log scale
    label_y = bar_y + legend_height + 10
    d.append(draw.Text(
        "0", font_size=7, x=x_start, y=label_y,
        fill=text_color, font_family='sans-serif', text_anchor='start'
    ))

    # Add intermediate tick at log midpoint
    if max_count > 10:
        mid_val = int(math.sqrt(max_count))  # Geometric midpoint
        mid_intensity = math.log1p(mid_val) / math.log1p(max_count)
        mid_x = x_start + mid_intensity * legend_width
        d.append(draw.Text(
            str(mid_val), font_size=7, x=mid_x, y=label_y,
            fill=text_color, font_family='sans-serif', text_anchor='middle'
        ))

    d.append(draw.Text(
        str(max_count), font_size=7, x=x_start + legend_width, y=label_y,
        fill=text_color, font_family='sans-serif', text_anchor='end'
    ))


def draw_sample_bar_plot(d, matrix_data, cluster_ids, cluster_enrichments, x_start, y_start,
                          cell_width, bar_height, text_color, background_color='black'):
    """Draw stacked vertical bar plot showing reads per sample by enrichment type.

    Args:
        d: Drawing object
        matrix_data: Dict returned from draw_sample_matrix containing sample info
        cluster_ids: List of cluster IDs included in the matrix
        cluster_enrichments: Dict of cluster_id -> enrichment type
        x_start: X position for bar plot start (same as matrix)
        y_start: Y position for bar plot top
        cell_width: Width of each bar (matches matrix cell width)
        bar_height: Maximum height of bars
        text_color: Color for text labels
        background_color: Background color ('black' or 'white')
    """
    sample_x_positions = matrix_data['sample_x_positions']
    cluster_sample_counts = matrix_data['cluster_sample_counts']
    all_samples = matrix_data['all_samples']

    # Colors for enrichment types
    enrichment_colors = {
        'Normal-enriched': '#3b82f6',  # Blue
        'Tumor-enriched': '#ef4444',   # Red
        'mixed': '#9ca3af'             # Gray
    }

    # Compute reads per sample by enrichment type
    sample_enrichment_counts = {sample: {'Normal-enriched': 0, 'Tumor-enriched': 0, 'mixed': 0}
                                 for sample in all_samples}
    sample_totals = {sample: 0 for sample in all_samples}

    for cid in cluster_ids:
        enrichment = cluster_enrichments.get(cid, 'mixed')
        for sample in all_samples:
            count = cluster_sample_counts.get(cid, {}).get(sample, 0)
            sample_enrichment_counts[sample][enrichment] += count
            sample_totals[sample] += count

    max_total = max(sample_totals.values()) if sample_totals.values() else 1

    # Draw stacked bars growing downward from y_start
    # Use thinner bars (max 8px) to match row barplot thickness
    thin_bar_width = min(cell_width - 2, 8)

    for sample in all_samples:
        # Center the thin bar within the cell
        x = x_start + sample_x_positions[sample] + (cell_width - thin_bar_width) / 2
        current_y = y_start

        # Stack order: Normal-enriched, mixed, Tumor-enriched (bottom to top visually = top to bottom in y)
        for enrichment in ['Normal-enriched', 'mixed', 'Tumor-enriched']:
            count = sample_enrichment_counts[sample][enrichment]
            if count > 0:
                bar_len = (count / max_total) * bar_height
                color = enrichment_colors.get(enrichment, '#888888')

                d.append(draw.Rectangle(
                    x, current_y,
                    thin_bar_width, bar_len,
                    fill=color,
                    stroke='none'
                ))
                current_y += bar_len

    # Draw axis line on the left
    axis_x = x_start - 3
    d.append(draw.Line(
        axis_x, y_start, axis_x, y_start + bar_height,
        stroke=text_color, stroke_width=1
    ))

    # Add tick marks and labels for axis
    tick_values = [0, max_total // 2, max_total]
    for val in tick_values:
        tick_y = y_start + (val / max_total) * bar_height if max_total > 0 else y_start
        # Tick mark
        d.append(draw.Line(
            axis_x - 3, tick_y, axis_x, tick_y,
            stroke=text_color, stroke_width=1
        ))
        # Label
        d.append(draw.Text(
            str(val), font_size=6,
            x=axis_x - 5, y=tick_y + 2,
            fill=text_color, font_family='sans-serif',
            text_anchor='end'
        ))


def draw_cluster_bar_plot(d, matrix_data, cluster_ids, cluster_y_start, cluster_y_end,
                          cluster_enrichments, x_start, bar_max_width, text_color, background_color='black'):
    """Draw horizontal stacked bar plot showing tumor vs normal reads per cluster (row sums).

    Args:
        d: Drawing object
        matrix_data: Dict returned from draw_sample_matrix containing sample info
        cluster_ids: List of cluster IDs
        cluster_y_start: Dict of cluster_id -> y_start position
        cluster_y_end: Dict of cluster_id -> y_end position
        cluster_enrichments: Dict of cluster_id -> enrichment type (unused, kept for API compat)
        x_start: X position for bar plot start (right edge of matrix)
        bar_max_width: Maximum width of bars
        text_color: Color for text labels
        background_color: Background color ('black' or 'white')
    """
    cluster_sample_counts = matrix_data['cluster_sample_counts']
    all_samples = matrix_data['all_samples']
    sample_groups = matrix_data.get('sample_groups', {})

    # Colors for sample groups
    group_colors = {
        'Normal': '#3b82f6',  # Blue
        'Tumor': '#ef4444',   # Red
    }

    # Separate samples by group
    normal_samples = [s for s in all_samples if sample_groups.get(s) == 'Normal']
    tumor_samples = [s for s in all_samples if sample_groups.get(s) == 'Tumor']

    # Compute reads per cluster by group
    cluster_group_counts = {}
    cluster_totals = {}
    for cid in cluster_ids:
        normal_count = sum(cluster_sample_counts.get(cid, {}).get(sample, 0) for sample in normal_samples)
        tumor_count = sum(cluster_sample_counts.get(cid, {}).get(sample, 0) for sample in tumor_samples)
        cluster_group_counts[cid] = {'Normal': normal_count, 'Tumor': tumor_count}
        cluster_totals[cid] = normal_count + tumor_count

    max_total = max(cluster_totals.values()) if cluster_totals.values() else 1

    # Draw horizontal stacked bars for each cluster
    for cid in cluster_ids:
        if cid not in cluster_y_start:
            continue

        y_start = cluster_y_start[cid]
        y_end = cluster_y_end[cid]
        y_center = (y_start + y_end) / 2
        bar_height = min(y_end - y_start - 2, 8)  # Bar height, max 8px

        current_x = x_start
        # Stack order: Normal first (blue), then Tumor (red)
        for group in ['Normal', 'Tumor']:
            count = cluster_group_counts.get(cid, {}).get(group, 0)
            if count > 0:
                bar_width = (count / max_total) * bar_max_width
                color = group_colors.get(group, '#888888')

                d.append(draw.Rectangle(
                    current_x, y_center - bar_height / 2,
                    bar_width, bar_height,
                    fill=color,
                    stroke='none'
                ))
                current_x += bar_width

    # Draw axis line on top
    axis_y = min(cluster_y_start.values()) - 5
    d.append(draw.Line(
        x_start, axis_y, x_start + bar_max_width, axis_y,
        stroke=text_color, stroke_width=1
    ))

    # Add tick marks and labels for axis
    tick_values = [0, max_total // 2, max_total]
    for val in tick_values:
        tick_x = x_start + (val / max_total) * bar_max_width if max_total > 0 else x_start
        # Tick mark
        d.append(draw.Line(
            tick_x, axis_y - 3, tick_x, axis_y,
            stroke=text_color, stroke_width=1
        ))
        # Label
        d.append(draw.Text(
            str(val), font_size=6,
            x=tick_x, y=axis_y - 5,
            fill=text_color, font_family='sans-serif',
            text_anchor='middle'
        ))


def draw_sample_dendrogram(d, matrix_data, x_start, y_bottom, dendro_height, line_color='#FFFFFF'):
    """Draw horizontal dendrogram above sample columns for each group.

    Args:
        d: Drawing object
        matrix_data: Dict returned from draw_sample_matrix
        x_start: X position where matrix starts
        y_bottom: Y position of dendrogram bottom (top of matrix headers)
        dendro_height: Height of dendrogram area
        line_color: Color for dendrogram lines
    """
    group_linkages = matrix_data.get('group_linkages', {})
    sample_x_positions = matrix_data['sample_x_positions']
    cell_width = matrix_data['cell_width']

    if not group_linkages:
        return

    # Draw dendrogram for each group
    for group_name, (Z, original_samples, ordered_samples) in group_linkages.items():
        n = len(original_samples)
        if n < 2:
            continue

        # Build sample to x-center mapping using DISPLAY positions
        sample_to_x = {}
        for sample in ordered_samples:
            sample_to_x[sample] = x_start + sample_x_positions[sample] + cell_width / 2

        # Normalize heights to fit in dendro_height
        max_height = Z[:, 2].max() if len(Z) > 0 else 1
        height_scale = (dendro_height - 5) / max_height if max_height > 0 else 1

        # Track node positions: node_id -> x_center
        # Leaf nodes are 0 to n-1 (indices into ORIGINAL sample order)
        # Internal nodes are n to 2n-2
        node_x = {}
        node_y = {}  # y position of the node (bottom of its subtree connection)

        # Initialize leaf nodes - linkage index i refers to original_samples[i]
        # but we need to look up the x position from the reordered display
        for i, sample in enumerate(original_samples):
            node_x[i] = sample_to_x[sample]  # Use display position of this sample
            node_y[i] = y_bottom

        # Process each merge in the linkage matrix
        for i, (left, right, height, _) in enumerate(Z):
            left, right = int(left), int(right)
            new_node = n + i

            # X position is midpoint of children
            x_left = node_x[left]
            x_right = node_x[right]
            x_mid = (x_left + x_right) / 2
            node_x[new_node] = x_mid

            # Y position based on height (growing upward from y_bottom)
            y_node = y_bottom - height * height_scale
            node_y[new_node] = y_node

            # Draw horizontal line connecting children
            d.append(draw.Line(
                x_left, y_node, x_right, y_node,
                stroke=line_color, stroke_width=1
            ))

            # Draw vertical lines down to children
            d.append(draw.Line(
                x_left, y_node, x_left, node_y[left],
                stroke=line_color, stroke_width=1
            ))
            d.append(draw.Line(
                x_right, y_node, x_right, node_y[right],
                stroke=line_color, stroke_width=1
            ))


def draw_cluster_labels_vertical(d, cluster_y_start, cluster_y_end, x_start, text_color,
                                  cluster_labels=None, enrichment_colors=None, cluster_enrichments=None):
    """Draw cluster labels on the RIGHT side for vertical mode (names only, no brackets).

    Args:
        d: Drawing object
        cluster_y_start: Dict of cluster_id -> y start position
        cluster_y_end: Dict of cluster_id -> y end position
        x_start: X position for labels
        text_color: Default text color
        cluster_labels: Dict of cluster_id -> custom label text
        enrichment_colors: Optional dict of enrichment -> color for coloring labels
        cluster_enrichments: Optional dict of cluster_id -> enrichment for coloring
    """
    if cluster_labels is None:
        cluster_labels = {}

    for cluster_id in cluster_y_start:
        y_start = cluster_y_start[cluster_id]
        y_end = cluster_y_end[cluster_id]
        label_y = (y_start + y_end) / 2

        # Get label text - include cluster ID in parenthesis for named clusters
        if cluster_id in cluster_labels and cluster_labels[cluster_id]:
            label_text = f"{cluster_labels[cluster_id]} ({cluster_id})"
        else:
            label_text = f"Cluster {cluster_id}"

        # Get color based on enrichment if available
        if enrichment_colors and cluster_enrichments and cluster_id in cluster_enrichments:
            enrichment = cluster_enrichments[cluster_id]
            color = enrichment_colors.get(enrichment, text_color)
        else:
            color = text_color

        d.append(draw.Text(
            label_text,
            font_size=9, x=x_start, y=label_y + 3,
            fill=color, font_family='sans-serif',
            text_anchor='start', font_weight='bold'
        ))


def draw_read_index_labels(d, cluster_reads, read_y_positions, x_position, bar_width, text_color='white'):
    """Draw index labels (1, 2, 3, ...) for each read within each cluster.

    Args:
        d: Drawing object
        cluster_reads: OrderedDict of cluster_id -> {'enrichment': str, 'reads': list of (read, sample)}
        read_y_positions: Dict of read -> y position
        x_position: X position for the labels (right side of feature bars)
        bar_width: Height of each read bar (for vertical centering)
        text_color: Color for the labels
    """
    for cluster_id, data in cluster_reads.items():
        if cluster_id == 'all':
            continue
        for idx, (read, sample) in enumerate(data['reads'], 1):
            if read not in read_y_positions:
                continue
            y_pos = read_y_positions[read]
            # Center the label vertically on the read bar
            label_y = y_pos + bar_width / 2 + 3
            d.append(draw.Text(
                str(idx),
                font_size=7, x=x_position, y=label_y,
                fill=text_color, font_family='sans-serif',
                text_anchor='start'
            ))


def draw_enrichment_bubbles(d, cluster_y_start, cluster_y_end, x_center, cluster_stats,
                            max_radius=8, min_radius=2):
    """Draw enrichment bubbles next to tree tips in vertical mode.

    Args:
        d: Drawing object
        cluster_y_start: Dict of cluster_id -> y start position
        cluster_y_end: Dict of cluster_id -> y end position
        x_center: X position for bubble centers
        cluster_stats: Dict of cluster_id -> {'odds_ratio': float, 'size': int, 'q_value': float}
        max_radius: Maximum bubble radius in pixels
        min_radius: Minimum bubble radius in pixels

    Bubble encoding:
        - Color: log2(odds_ratio) mapped to diverging colormap (blue=Normal, red=Tumor)
        - Size: cluster size (number of reads)
        - Alpha: -log10(q_value) mapped to opacity (more significant = more opaque)
    """
    import math

    if not cluster_stats:
        return

    # Get ranges for normalization
    sizes = [s['size'] for s in cluster_stats.values() if s['size'] > 0]
    if not sizes:
        return

    min_size = min(sizes)
    max_size = max(sizes)
    size_range = max_size - min_size if max_size > min_size else 1

    # Max log2(OR) for color scaling (cap at 4 for visualization)
    max_log2_or = 4.0

    for cluster_id in cluster_y_start:
        if cluster_id not in cluster_stats:
            continue

        stats = cluster_stats[cluster_id]
        odds_ratio = stats.get('odds_ratio', 1.0)
        size = stats.get('size', 0)
        q_value = stats.get('q_value', 1.0)

        if size == 0:
            continue

        # Y position (center of cluster)
        y_start = cluster_y_start[cluster_id]
        y_end = cluster_y_end[cluster_id]
        y_center = (y_start + y_end) / 2

        # Size -> radius (linear scaling)
        size_norm = (size - min_size) / size_range if size_range > 0 else 0.5
        radius = min_radius + size_norm * (max_radius - min_radius)

        # log2(OR) -> color (diverging: red for OR<1 Tumor, blue for OR>1 Normal)
        log2_or = math.log2(odds_ratio) if odds_ratio > 0 else 0
        # Clamp to [-max_log2_or, max_log2_or]
        log2_or_clamped = max(-max_log2_or, min(max_log2_or, log2_or))
        # Normalize to [-1, 1]
        color_norm = log2_or_clamped / max_log2_or

        # Diverging colormap: red (Tumor-enriched, OR<1) to blue (Normal-enriched, OR>1)
        if color_norm < 0:
            # Red side (Tumor-enriched, OR<1): interpolate white to red
            intensity = abs(color_norm)
            r = 255
            g = int(255 * (1 - intensity))
            b = int(255 * (1 - intensity))
        else:
            # Blue side (Normal-enriched, OR>1): interpolate white to blue
            intensity = color_norm
            r = int(255 * (1 - intensity))
            g = int(255 * (1 - intensity))
            b = 255

        color = f'rgb({r},{g},{b})'

        # q_value -> alpha (more significant = more opaque)
        # Use -log10(q_value), capped at 10 (q=1e-10)
        if q_value > 0:
            neg_log_q = -math.log10(q_value)
        else:
            neg_log_q = 10  # Very significant
        neg_log_q = min(10, max(0, neg_log_q))
        # Map to alpha: 0.3 (q=1) to 1.0 (q<=1e-10)
        alpha = 0.3 + 0.7 * (neg_log_q / 10)

        # Draw the bubble
        d.append(draw.Circle(
            x_center, y_center, radius,
            fill=color, fill_opacity=alpha,
            stroke='white', stroke_width=0.5
        ))


def draw_bubble_legend(d, x_start, y_start, cluster_stats, text_color='white', max_radius=8, min_radius=2):
    """Draw a legend explaining the enrichment bubble encoding.

    Args:
        d: Drawing object
        x_start: X position to start the legend
        y_start: Y position for legend
        cluster_stats: Dict with cluster stats (for computing size range)
        text_color: Color for legend text
        max_radius: Maximum bubble radius
        min_radius: Minimum bubble radius
    """
    import math

    if not cluster_stats:
        return

    # Get size range for legend
    sizes = [s['size'] for s in cluster_stats.values() if s['size'] > 0]
    if not sizes:
        return
    min_size = min(sizes)
    max_size = max(sizes)

    font_size = 9
    legend_y = y_start

    # Title
    d.append(draw.Text(
        "Enrichment Bubble Legend",
        font_size=10, x=x_start, y=legend_y,
        fill=text_color, font_family='sans-serif', font_weight='bold'
    ))
    legend_y += 18

    # --- Color legend (log2 FC) ---
    d.append(draw.Text(
        "Color: log₂(Odds Ratio)",
        font_size=font_size, x=x_start, y=legend_y,
        fill=text_color, font_family='sans-serif'
    ))
    legend_y += 14

    # Draw color gradient with labels
    gradient_width = 120
    gradient_height = 10
    n_steps = 20
    step_width = gradient_width / n_steps

    for i in range(n_steps):
        # Map i to [-1, 1]
        color_norm = (i / (n_steps - 1)) * 2 - 1
        # Red (Tumor, OR<1) on left, Blue (Normal, OR>1) on right
        if color_norm < 0:
            # Red side (Tumor)
            intensity = abs(color_norm)
            r = 255
            g = int(255 * (1 - intensity))
            b = int(255 * (1 - intensity))
        else:
            # Blue side (Normal)
            intensity = color_norm
            r = int(255 * (1 - intensity))
            g = int(255 * (1 - intensity))
            b = 255
        color = f'rgb({r},{g},{b})'

        d.append(draw.Rectangle(
            x_start + i * step_width, legend_y,
            step_width + 1, gradient_height,
            fill=color, stroke='none'
        ))

    # Gradient labels
    legend_y += gradient_height + 12
    d.append(draw.Text("-4", font_size=8, x=x_start, y=legend_y,
                       fill=text_color, font_family='sans-serif', text_anchor='start'))
    d.append(draw.Text("0", font_size=8, x=x_start + gradient_width/2, y=legend_y,
                       fill=text_color, font_family='sans-serif', text_anchor='middle'))
    d.append(draw.Text("+4", font_size=8, x=x_start + gradient_width, y=legend_y,
                       fill=text_color, font_family='sans-serif', text_anchor='end'))
    d.append(draw.Text("(Tumor)", font_size=7, x=x_start, y=legend_y + 10,
                       fill='#ff6666', font_family='sans-serif', text_anchor='start'))
    d.append(draw.Text("(Normal)", font_size=7, x=x_start + gradient_width, y=legend_y + 10,
                       fill='#6666ff', font_family='sans-serif', text_anchor='end'))

    # --- Size legend ---
    size_x = x_start + gradient_width + 40
    size_y = y_start + 18

    d.append(draw.Text(
        "Size: # reads",
        font_size=font_size, x=size_x, y=size_y,
        fill=text_color, font_family='sans-serif'
    ))
    size_y += 16

    # Create sensible size breakpoints based on data range
    size_range = max_size - min_size
    if size_range > 0:
        # Pick 3 representative sizes: small, medium, large
        size_values = [min_size, int(min_size + size_range * 0.5), max_size]
        # Remove duplicates and sort
        size_values = sorted(set(size_values))
    else:
        size_values = [min_size]

    # Draw example bubbles for each size
    bubble_x_offset = 0
    for size_val in size_values:
        size_norm = (size_val - min_size) / size_range if size_range > 0 else 0.5
        radius = min_radius + size_norm * (max_radius - min_radius)
        d.append(draw.Circle(size_x + bubble_x_offset + radius, size_y + 5, radius,
                             fill='white', fill_opacity=0.8, stroke='gray', stroke_width=0.5))
        d.append(draw.Text(f"{size_val}", font_size=7, x=size_x + bubble_x_offset + radius, y=size_y + 18,
                           fill=text_color, font_family='sans-serif', text_anchor='middle'))
        bubble_x_offset += radius * 2 + 15

    # --- Alpha/FDR legend ---
    alpha_x = size_x + bubble_x_offset + 20
    alpha_y = y_start + 18

    d.append(draw.Text(
        "Opacity: FDR",
        font_size=font_size, x=alpha_x, y=alpha_y,
        fill=text_color, font_family='sans-serif'
    ))
    alpha_y += 16

    # Draw example bubbles for alpha
    for i, (alpha, label) in enumerate([(0.3, "1"), (0.65, "0.01"), (1.0, "<1e-10")]):
        cx = alpha_x + i * 40
        d.append(draw.Circle(cx + 5, alpha_y + 5, 5,
                             fill='white', fill_opacity=alpha, stroke='gray', stroke_width=0.5))
        d.append(draw.Text(label, font_size=7, x=cx + 5, y=alpha_y + 18,
                           fill=text_color, font_family='sans-serif', text_anchor='middle'))


def draw_enrichment_text_legend(d, x_start, y_start, enrichment_colors, text_color='white'):
    """Draw a legend explaining the text color encoding for cluster labels.

    Args:
        d: Drawing object
        x_start: X position to start the legend
        y_start: Y position for legend
        enrichment_colors: Dict of enrichment -> color
        text_color: Color for legend text
    """
    if not enrichment_colors:
        return

    font_size = 9
    legend_y = y_start

    # Title
    d.append(draw.Text(
        "Cluster Label Colors",
        font_size=10, x=x_start, y=legend_y,
        fill=text_color, font_family='sans-serif', font_weight='bold'
    ))
    legend_y += 16

    # Draw color legend entries horizontally
    entry_spacing = 120
    x_offset = 0

    for enrichment, color in sorted(enrichment_colors.items()):
        # Draw colored square
        d.append(draw.Rectangle(
            x_start + x_offset, legend_y - 7,
            10, 10,
            fill=color, stroke='none'
        ))
        # Draw label
        d.append(draw.Text(
            enrichment,
            font_size=font_size, x=x_start + x_offset + 14, y=legend_y,
            fill=text_color, font_family='sans-serif'
        ))
        x_offset += entry_spacing


def draw_cluster_brackets_vertical(d, cluster_reads, cluster_y_start, cluster_y_end,
                                   enrichment_colors, read_positions, label_width, text_color,
                                   right_margin, cluster_labels=None):
    """Draw cluster brackets on the RIGHT side for vertical mode (with brackets and enrichment)."""
    if cluster_labels is None:
        cluster_labels = {}

    for cluster_id, data in cluster_reads.items():
        if cluster_id == 'all':
            continue

        y_start = cluster_y_start[cluster_id]
        y_end = cluster_y_end[cluster_id]
        enrichment = data['enrichment']
        color = enrichment_colors.get(enrichment, '#999999')

        # Find the rightmost x of reads in this cluster
        cluster_max_x = 0
        for read, sample in data['reads']:
            if read in read_positions:
                x_start, _, _ = read_positions[read]
                cluster_max_x = max(cluster_max_x, x_start + label_width)

        bracket_x = cluster_max_x + 10

        # Draw bracket (vertical line on right, horizontal lines pointing left)
        d.append(draw.Line(bracket_x, y_start, bracket_x, y_end, stroke=color, stroke_width=3))
        d.append(draw.Line(bracket_x, y_start, bracket_x - 8, y_start, stroke=color, stroke_width=2))
        d.append(draw.Line(bracket_x, y_end, bracket_x - 8, y_end, stroke=color, stroke_width=2))

        # Labels to the right of bracket
        label_y = (y_start + y_end) / 2

        if cluster_id in cluster_labels:
            label_text = cluster_labels[cluster_id]
            label_font_size = 9
        else:
            label_text = f"c{cluster_id}"
            label_font_size = 10

        d.append(draw.Text(
            label_text,
            font_size=label_font_size, x=bracket_x + 5, y=label_y - 8,
            fill=color, font_family='sans-serif',
            text_anchor='start', font_weight='bold'
        ))

        d.append(draw.Text(
            enrichment,
            font_size=8, x=bracket_x + 5, y=label_y + 8,
            fill=color, font_family='sans-serif', text_anchor='start'
        ))


def draw_annotation_bars(d, cluster_reads, read_x_positions, read_to_original_cluster,
                         read_to_original_enrichment, sample_colors, cluster_colors,
                         enrichment_colors, group_width, top_margin, left_margin, text_color):
    """Draw cluster, enrichment, and sample annotation bars directly below dendrogram tips with labels."""
    # Annotation bars start at dendrogram base (top_margin + 22, just below dendrogram tips at +20)
    annot_start_y = top_margin + 22

    # Draw labels on the left side
    label_x = left_margin - 10

    # Cluster label
    d.append(draw.Text(
        "Cluster", font_size=8, x=label_x, y=annot_start_y + 4,
        fill=text_color, font_family='sans-serif',
        text_anchor='end'
    ))

    # Sample label
    d.append(draw.Text(
        "Sample", font_size=8, x=label_x, y=annot_start_y + 20,
        fill=text_color, font_family='sans-serif',
        text_anchor='end'
    ))

    # Draw annotation bars for each read
    for cluster_id, data in cluster_reads.items():
        for read, sample in data['reads']:
            if read not in read_x_positions:
                continue

            base_x = read_x_positions[read]
            sample_color = sample_colors.get(sample, '#999999')
            orig_cluster = read_to_original_cluster.get(read, 'unknown')
            orig_cluster_color = cluster_colors.get(orig_cluster, '#666666')
            orig_enrichment = read_to_original_enrichment.get(read, 'mixed')
            enrichment_color = enrichment_colors.get(orig_enrichment, '#999999')

            # Cluster indicator bar (top, 8px height)
            d.append(draw.Rectangle(base_x, annot_start_y, group_width, 8, fill=orig_cluster_color))

            # Enrichment indicator bar (middle, thin 3px height)
            d.append(draw.Rectangle(base_x, annot_start_y + 10, group_width, 3, fill=enrichment_color))

            # Sample indicator bar (bottom, 8px height)
            d.append(draw.Rectangle(base_x, annot_start_y + 15, group_width, 8, fill=sample_color))


def draw_feature_bars(d, drawing_data, featuresets, bar_width, read_heights, num_featuresets,
                      density_line_data=None, rect_plot_data=None, background_color="black"):
    """Draw feature rectangles with borders between featuresets and around outer edge.

    Args:
        d: Drawing object
        drawing_data: Feature data per read
        featuresets: List of feature sets
        bar_width: Width of each bar
        read_heights: Dict of read -> (min_y, max_y, x_start, total_width) for borders
        num_featuresets: Number of feature sets
        density_line_data: Optional dict of read -> featureset -> {points, color, base_x} for line plots
        rect_plot_data: Optional dict of read -> featureset -> list of {y, height, color, base_x} for rect plots
        background_color: Background color for determining line plot background fill
    """
    stroke_width = 0.5

    # First draw all feature rectangles (no individual borders)
    for read in drawing_data:
        for fs in featuresets:
            for rect in drawing_data[read][fs]:
                if rect["height"] > 0 and rect["fill"] != "none":
                    d.append(draw.Rectangle(
                        rect["x"], rect["y"],
                        bar_width, rect["height"],
                        fill=rect["fill"],
                        fill_opacity=rect["fill_opacity"]
                    ))

    # Draw density line plots if present
    if density_line_data:
        for read in density_line_data:
            if read not in read_heights:
                continue
            min_y, max_y, x_start, total_width = read_heights[read]

            # Draw background for line plot area (dark gray)
            for line_fs, line_data in density_line_data[read].items():
                base_x = line_data['base_x']
                # Draw background rectangle for the line plot area
                d.append(draw.Rectangle(
                    base_x, min_y,
                    bar_width, max_y - min_y,
                    fill="#1a1a1a" if background_color == "black" else "#f0f0f0",
                    fill_opacity=1.0
                ))

            # Draw each line on top
            for line_fs, line_data in density_line_data[read].items():
                points = line_data['points']
                color = line_data['color']
                base_x = line_data['base_x']

                if len(points) >= 2:
                    # Create polyline path
                    # Build points string for polyline
                    points_str = " ".join(f"{x},{y}" for x, y in points)

                    # Add filled area under the line (to base_x)
                    # Create path: start at bottom-left, go up the line, then back down
                    path_d = f"M {base_x},{points[0][1]} "
                    for x, y in points:
                        path_d += f"L {x},{y} "
                    path_d += f"L {base_x},{points[-1][1]} Z"

                    d.append(draw.Path(
                        d=path_d,
                        fill=color,
                        fill_opacity=0.3,
                        stroke='none'
                    ))

                    # Draw the line itself
                    d.append(draw.Lines(
                        *[coord for point in points for coord in point],
                        stroke=color,
                        stroke_width=1.0,
                        fill='none'
                    ))

    # Draw rect plot rectangles if present (e.g., FIRE/Linker exact calls)
    if rect_plot_data:
        for read in rect_plot_data:
            if read not in read_heights:
                continue
            min_y, max_y, x_start, total_width = read_heights[read]

            # Draw background for rect plot area (dark gray, same as density line)
            first_rect = None
            for rect_fs in rect_plot_data[read]:
                if rect_plot_data[read][rect_fs]:
                    first_rect = rect_plot_data[read][rect_fs][0]
                    break
            if first_rect:
                d.append(draw.Rectangle(
                    first_rect['base_x'], min_y,
                    bar_width, max_y - min_y,
                    fill="#1a1a1a" if background_color == "black" else "#f0f0f0",
                    fill_opacity=1.0
                ))

            # Draw each feature rectangle
            for rect_fs, rects in rect_plot_data[read].items():
                for rect in rects:
                    rect_height = max(rect['height'], 2)  # Minimum 2px height
                    d.append(draw.Rectangle(
                        rect['base_x'], rect['y'],
                        bar_width, rect_height,
                        fill=rect['color'],
                        fill_opacity=0.85
                    ))

    # Draw borders: outer border + vertical lines between featuresets
    for read, (min_y, max_y, x_start, total_width) in read_heights.items():
        # Outer border
        d.append(draw.Rectangle(
            x_start, min_y,
            total_width, max_y - min_y,
            fill='none',
            stroke='black',
            stroke_width=stroke_width
        ))

        # Vertical lines between featureset bars
        for i in range(1, num_featuresets):
            line_x = x_start + i * bar_width
            d.append(draw.Line(
                line_x, min_y, line_x, max_y,
                stroke='black',
                stroke_width=stroke_width
            ))


def abbreviate_read_name(read_name, max_len=12):
    """Abbreviate read name to a unique, short identifier.

    Handles different read name formats:
    - PacBio HiFi: m84132_240112_213928_s2/201131976/ccs -> 201131976
    - ONT: uuid format -> first 8 chars
    - Generic: first max_len chars
    """
    # PacBio format: movie/zmw/ccs or movie/zmw/subread
    if '/' in read_name:
        parts = read_name.split('/')
        if len(parts) >= 2:
            # Return the ZMW number (second part)
            return parts[1][:max_len]

    # Default: first max_len characters
    return read_name[:max_len]


def draw_read_labels(d, cluster_reads, read_x_positions, group_width, top_margin, text_color):
    """Draw read ID labels above the annotation bars, rotated 90 degrees."""
    for cluster_id, data in cluster_reads.items():
        for read, sample in data['reads']:
            if read not in read_x_positions:
                continue

            base_x = read_x_positions[read]
            label_x = base_x + group_width / 2
            label_y = top_margin - 5  # Position above annotation bars
            short_id = abbreviate_read_name(read)

            d.append(draw.Text(
                short_id, font_size=5, x=label_x, y=label_y,
                fill=text_color, font_family='monospace',
                transform=f"rotate(-90 {label_x} {label_y})",
                text_anchor='start'
            ))


def draw_top_legends(d, sample_colors, cluster_colors, read_to_original_cluster,
                     read_to_original_enrichment, enrichment_colors, legend_x, legend_y, text_color):
    """Draw sample, cluster, and enrichment legends stacked vertically."""
    row_height = 20

    # --- Sample legend (row 1) ---
    current_y = legend_y
    d.append(draw.Text(
        "Sample:", font_size=10, x=legend_x, y=current_y,
        fill=text_color, font_family='sans-serif', font_weight='bold'
    ))

    # Calculate spacing based on sample name lengths (font_size=9, ~5.5px per char)
    item_x = legend_x + 60
    for sample, color in sample_colors.items():
        d.append(draw.Rectangle(item_x, current_y - 8, 12, 12, fill=color))
        d.append(draw.Text(
            sample, font_size=9, x=item_x + 16, y=current_y,
            fill=text_color, font_family='sans-serif'
        ))
        # Move to next position: box(12) + gap(4) + text width + padding(15)
        item_x += 12 + 4 + len(sample) * 5.5 + 15

    # --- Cluster legend (row 2) ---
    current_y += row_height
    d.append(draw.Text(
        "Cluster:", font_size=10, x=legend_x, y=current_y,
        fill=text_color, font_family='sans-serif', font_weight='bold'
    ))

    # Group clusters by enrichment
    clusters_by_enrichment = defaultdict(list)
    for read, cid in read_to_original_cluster.items():
        enrich = read_to_original_enrichment.get(read, 'mixed')
        if cid not in clusters_by_enrichment[enrich]:
            clusters_by_enrichment[enrich].append(cid)

    # Sort clusters within each enrichment
    for enrich in clusters_by_enrichment:
        clusters_by_enrichment[enrich] = sorted(set(clusters_by_enrichment[enrich]))

    # Draw cluster items
    cluster_legend_x = legend_x + 60

    # Sort enrichment types for consistent ordering
    enrichment_order = sorted(clusters_by_enrichment.keys(),
                             key=lambda x: (x == 'mixed', x))  # mixed last

    for enrich_type in enrichment_order:
        for cid in clusters_by_enrichment[enrich_type]:
            color = cluster_colors.get(cid, '#666666')
            d.append(draw.Rectangle(cluster_legend_x, current_y - 8, 12, 12, fill=color))

            # Enrichment indicator
            enrich_color = enrichment_colors.get(enrich_type, '#999999')
            d.append(draw.Rectangle(cluster_legend_x, current_y + 5, 12, 3, fill=enrich_color))

            d.append(draw.Text(
                f"C{cid}", font_size=8, x=cluster_legend_x + 15, y=current_y,
                fill=text_color, font_family='sans-serif'
            ))
            cluster_legend_x += 45

    # --- Enrichment legend (row 3) ---
    current_y += row_height
    d.append(draw.Text(
        "Enrichment:", font_size=10, x=legend_x, y=current_y,
        fill=text_color, font_family='sans-serif', font_weight='bold'
    ))

    for i, (enrich, color) in enumerate(sorted(enrichment_colors.items())):
        item_x = legend_x + 75 + i * 100
        d.append(draw.Rectangle(item_x, current_y - 8, 12, 12, fill=color))

        # Clean up label for display
        short_enrich = enrich.replace('-enriched', '')
        d.append(draw.Text(
            short_enrich, font_size=9, x=item_x + 16, y=current_y,
            fill=text_color, font_family='sans-serif'
        ))


def draw_density_line_legend(d, density_line_colors, rect_plot_colors, legend_x, legend_y, text_color):
    """Draw legend for density line plot and rect plot tracks.

    Args:
        d: Drawing object
        density_line_colors: Dict of featureset -> color for density lines
        rect_plot_colors: Dict of featureset -> color for rect plots (or None if not shown)
        legend_x: X position for legend
        legend_y: Y position for legend
        text_color: Text color
    """
    if not density_line_colors and not rect_plot_colors:
        return

    current_x = legend_x

    # Density line (m6A/5mC) legend
    if density_line_colors:
        d.append(draw.Text(
            "m6A/5mC:", font_size=10, x=current_x, y=legend_y,
            fill=text_color, font_family='sans-serif', font_weight='bold'
        ))
        current_x += 60

        # Display name mapping
        display_names = {
            'fiberseq_m6A': 'm6A',
            'fiberseq_5mC': '5mC',
        }

        for fs, color in density_line_colors.items():
            # Draw a small line sample
            d.append(draw.Lines(
                current_x, legend_y - 3,
                current_x + 20, legend_y - 3,
                stroke=color, stroke_width=2, fill='none'
            ))
            # Draw filled area sample under line
            d.append(draw.Rectangle(
                current_x, legend_y - 3, 20, 6,
                fill=color, fill_opacity=0.3
            ))
            display_name = display_names.get(fs, fs.replace('fiberseq_', ''))
            d.append(draw.Text(
                display_name, font_size=9, x=current_x + 25, y=legend_y,
                fill=text_color, font_family='sans-serif'
            ))
            current_x += 25 + len(display_name) * 6 + 20

        current_x += 20  # Extra spacing before rect plot legend

    # Rect plot (FIRE/Linker) legend
    if rect_plot_colors:
        d.append(draw.Text(
            "FIRE/Linker:", font_size=10, x=current_x, y=legend_y,
            fill=text_color, font_family='sans-serif', font_weight='bold'
        ))
        current_x += 85

        # Display name mapping for rect plot features
        rect_display_names = {
            'fiberseq_FIRE': 'FIRE',
            'fiberseq_LINKER': 'Linker',
        }

        for fs, color in rect_plot_colors.items():
            d.append(draw.Rectangle(
                current_x, legend_y - 6, 20, 10,
                fill=color, fill_opacity=0.8,
                stroke=color, stroke_width=1
            ))
            display_name = rect_display_names.get(fs, fs.replace('fiberseq_', ''))
            d.append(draw.Text(
                display_name, font_size=9, x=current_x + 25, y=legend_y,
                fill=text_color, font_family='sans-serif'
            ))
            current_x += 25 + len(display_name) * 6 + 20


def _find_common_label(feature_names):
    """
    Find a common/generic label for multiple features that share a color.

    Examples:
        - ['multigroup1', 'multigroup2', 'multigroup3'] → 'multigroup'
        - ['centromeric_multigroup1', 'hsat_multigroup1'] → 'multigroup'
        - ['p_arm_specific', 'q_arm_specific'] → 'arm'
        - ['chr1_specific', 'chr2_specific'] → 'chromosome'
    """
    import re

    if not feature_names:
        return 'unknown'

    # Clean up names: remove _specific suffix and replace _ with space
    def clean_name(f):
        name = f.replace('_specific', '').replace('_', ' ')
        # Strip trailing numbers from multigroup names
        if 'multigroup' in name.lower():
            name = re.sub(r'\s*\d+$', '', name).rstrip()
        return name

    cleaned = [clean_name(f) for f in feature_names]

    # Remove duplicates while preserving order
    seen = set()
    unique_cleaned = []
    for name in cleaned:
        if name not in seen:
            seen.add(name)
            unique_cleaned.append(name)

    # Check if all unique features are just "multigroup" variants
    if all('multigroup' in name.lower() for name in unique_cleaned):
        return 'multigroup'

    # Check if all features contain "arm" - use generic "arm"
    if all('arm' in name.lower() for name in unique_cleaned):
        return 'arm'

    # If only one unique name after cleaning, return it
    if len(unique_cleaned) == 1:
        return unique_cleaned[0]

    # Check for numbered suffixes with same base
    base_names = set()
    for name in unique_cleaned:
        base = re.sub(r'\d+$', '', name).strip()
        base_names.add(base)

    if len(base_names) == 1:
        return base_names.pop()

    # Find longest common prefix
    if unique_cleaned:
        prefix = unique_cleaned[0]
        for name in unique_cleaned[1:]:
            while not name.startswith(prefix) and prefix:
                prefix = prefix[:-1]
        if len(prefix) >= 3:  # Minimum meaningful prefix
            return prefix.rstrip()

    # Fallback: use first cleaned name
    return unique_cleaned[0] if unique_cleaned else 'unknown'


def draw_color_legends(d, featuresets, featureset_colors, featureset_color_order,
                       fs_display_names, color_legend_y_start, left_margin, text_color,
                       displayed_features=None):
    """Draw featureset color legends at bottom.

    Args:
        d: Drawing object
        featuresets: List of featuresets to include
        featureset_colors: Dict of fs -> feature -> (color, opacity)
        featureset_color_order: Dict of fs -> list of feature names (ordering)
        fs_display_names: Dict of fs -> display name
        color_legend_y_start: Y position to start legend
        left_margin: Left margin for positioning
        text_color: Color for text
        displayed_features: Optional dict of fs -> set of actually displayed features.
                           If provided, only these features are shown in the legend,
                           grouped by color with common labels.
    """
    import re

    color_box_size = 10
    color_text_offset = 14
    colors_per_column = 12  # More items per column = fewer columns
    item_width = 120  # Narrower columns
    row_height = 14  # Tighter row spacing

    # Build legend items for each featureset
    # If displayed_features is provided, filter and group by color
    legend_items = {}  # fs -> list of (color, display_name)

    for fs in featuresets:
        fs_colors = featureset_colors.get(fs, {})

        if displayed_features is not None and fs in displayed_features:
            # Filter to only displayed features, group by color
            fs_displayed = displayed_features.get(fs, set())
            if not fs_displayed:
                legend_items[fs] = []
                continue

            # Group features by color
            color_to_features = {}
            for feature_name in fs_displayed:
                color_info = fs_colors.get(feature_name)
                if color_info is None:
                    color_hex = '#808080'  # Default gray for unknown features
                elif isinstance(color_info, tuple):
                    color_hex = color_info[0]
                else:
                    color_hex = color_info
                if color_hex not in color_to_features:
                    color_to_features[color_hex] = []
                color_to_features[color_hex].append(feature_name)

            # Sort by feature name for consistency
            sorted_items = sorted(color_to_features.items(), key=lambda x: x[1][0] if x[1] else '')

            items = []
            for color_hex, feature_names in sorted_items:
                if len(feature_names) > 1:
                    display_name = _find_common_label(feature_names)
                else:
                    # Single feature - clean up and strip trailing numbers from multigroup names
                    display_name = feature_names[0].replace('_', ' ').replace(' specific', '')
                    if 'multigroup' in display_name.lower():
                        display_name = re.sub(r'\d+$', '', display_name).rstrip()
                items.append((color_hex, display_name))
            legend_items[fs] = items
        else:
            # No filtering - use all features in order
            items = []
            for feature_name in featureset_color_order.get(fs, []):
                color, opacity = fs_colors.get(feature_name, ("#ffffff", 1.0))
                items.append((color, feature_name))
            legend_items[fs] = items

    def get_featureset_width(fs):
        num_items = len(legend_items.get(fs, []))
        if num_items == 0:
            return 0
        num_cols = (num_items + colors_per_column - 1) // colors_per_column
        return max(num_cols * item_width, 100)

    # Calculate x positions for each featureset
    featureset_legend_x = {}
    current_legend_x = left_margin
    for fs in featuresets:
        if not legend_items.get(fs):
            continue
        featureset_legend_x[fs] = current_legend_x
        current_legend_x += get_featureset_width(fs) + 15  # Reduced spacing between sections

    # Draw legends
    for fs in featuresets:
        items = legend_items.get(fs, [])
        if not items:
            continue

        section_x = featureset_legend_x[fs]
        display_name = fs_display_names.get(fs, fs)

        d.append(draw.Text(
            display_name, font_size=9, x=section_x, y=color_legend_y_start,
            fill=text_color, font_family='sans-serif', font_weight='bold'
        ))

        for i, (color, feature_label) in enumerate(items):
            row = i % colors_per_column
            col = i // colors_per_column

            item_x = section_x + col * item_width
            item_y = color_legend_y_start + 18 + row * row_height

            d.append(draw.Rectangle(
                item_x, item_y - 7, color_box_size, color_box_size,
                fill=color, stroke=text_color, stroke_width=0.5
            ))

            d.append(draw.Text(
                feature_label, font_size=7, x=item_x + color_text_offset, y=item_y,
                fill=text_color, font_family='sans-serif'
            ))


def draw_color_legends_vertical(d, featuresets, featureset_colors, featureset_color_order,
                                fs_display_names, legend_x, legend_y_start, text_color,
                                displayed_features=None):
    """Draw featureset color legends vertically on the right side.

    Args:
        d: Drawing object
        featuresets: List of featuresets to include
        featureset_colors: Dict of fs -> feature -> (color, opacity)
        featureset_color_order: Dict of fs -> list of feature names (ordering)
        fs_display_names: Dict of fs -> display name
        legend_x: X position for legend (right side)
        legend_y_start: Y position to start legend
        text_color: Color for text
        displayed_features: Optional dict of fs -> set of actually displayed features.

    Returns:
        float: Total height used by the legend
    """
    import re

    swatch_size = 8
    item_height = 11
    section_gap = 15  # Gap between featureset sections

    # Build legend items for each featureset (same logic as horizontal version)
    legend_items = {}

    for fs in featuresets:
        fs_colors = featureset_colors.get(fs, {})

        if displayed_features is not None and fs in displayed_features:
            fs_displayed = displayed_features.get(fs, set())
            if not fs_displayed:
                legend_items[fs] = []
                continue

            # Group features by color
            color_to_features = {}
            for feature_name in fs_displayed:
                color_info = fs_colors.get(feature_name)
                if color_info is None:
                    color_hex = '#808080'
                elif isinstance(color_info, tuple):
                    color_hex = color_info[0]
                else:
                    color_hex = color_info
                if color_hex not in color_to_features:
                    color_to_features[color_hex] = []
                color_to_features[color_hex].append(feature_name)

            # Sort by feature name for consistency
            sorted_items = sorted(color_to_features.items(), key=lambda x: x[1][0] if x[1] else '')

            items = []
            for color_hex, feature_names in sorted_items:
                if len(feature_names) > 1:
                    display_name = _find_common_label(feature_names)
                else:
                    display_name = feature_names[0].replace('_', ' ').replace(' specific', '')
                    if 'multigroup' in display_name.lower():
                        display_name = re.sub(r'\d+$', '', display_name).rstrip()
                items.append((color_hex, display_name))
            legend_items[fs] = items
        else:
            items = []
            for feature_name in featureset_color_order.get(fs, []):
                color, opacity = fs_colors.get(feature_name, ("#ffffff", 1.0))
                items.append((color, feature_name))
            legend_items[fs] = items

    # Draw legends vertically
    current_y = legend_y_start
    first_section = True

    for fs in featuresets:
        items = legend_items.get(fs, [])
        if not items:
            continue

        # Add gap before section (except first)
        if not first_section:
            current_y += section_gap
        first_section = False

        # Section header
        display_name = fs_display_names.get(fs, fs.title())
        d.append(draw.Text(
            display_name, font_size=7, x=legend_x, y=current_y,
            fill=text_color, font_family='sans-serif',
            text_anchor='start', font_weight='bold'
        ))
        current_y += 4  # Small gap after header

        # Draw each item
        for color_hex, feature_label in items:
            current_y += item_height

            # Color swatch
            d.append(draw.Rectangle(
                legend_x + 3, current_y - swatch_size, swatch_size, swatch_size,
                fill=color_hex
            ))

            # Label
            d.append(draw.Text(
                feature_label, font_size=6, x=legend_x + swatch_size + 6,
                y=current_y - swatch_size/2 + 1,
                fill=text_color, font_family='sans-serif',
                text_anchor='start', dominant_baseline='middle'
            ))

    return current_y - legend_y_start  # Return total height used


# =============================================================================
# Main Script
# =============================================================================



class TeeLogger:
    """Write to both stdout and a log file."""
    def __init__(self, log_path):
        self.terminal = sys.stdout
        self.log = open(log_path, 'w')

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()

    def flush(self):
        self.terminal.flush()
        self.log.flush()

    def close(self):
        self.log.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False


def compute_density_features(features, bin_size, read_length, feature_name):
    """Convert individual features to density bins.

    For features with many small marks (e.g., m6A, 5mC), this bins them
    and returns density levels (none/low/medium/high) per bin.
    """
    if not features or read_length <= 0:
        return []

    # Create bins
    n_bins = (read_length + bin_size - 1) // bin_size
    bin_counts = [0] * n_bins

    # Count features per bin (count each feature once per bin it overlaps)
    for feat in features:
        start_bin = feat['start'] // bin_size
        end_bin = min((feat['stop'] - 1) // bin_size, n_bins - 1)
        if end_bin < 0:
            continue
        for b in range(max(0, start_bin), end_bin + 1):
            bin_counts[b] += 1

    # Convert counts to density labels
    # Thresholds for 300bp bin with single-base features:
    # - none: 0 marks
    # - low: 1-5 marks (~1-2% density)
    # - medium: 6-15 marks (~2-5% density)
    # - high: 16+ marks (>5% density)
    density_features = []
    for b in range(n_bins):
        count = bin_counts[b]
        if count == 0:
            level = 'none'
        elif count <= 5:
            level = 'low'
        elif count <= 15:
            level = 'medium'
        else:
            level = 'high'

        density_features.append({
            'start': b * bin_size,
            'stop': min((b + 1) * bin_size, read_length),
            'feature': f"{feature_name}_{level}"
        })

    return density_features


def compute_density_line(features, bin_size, read_length, max_density=50):
    """Compute density as percentage per bin for line plot rendering.

    Args:
        features: List of feature dicts with 'start' and 'stop' keys
        bin_size: Size of each bin in bp
        read_length: Total length of the read
        max_density: Number of marks considered 100% density (default: 50)

    Returns:
        List of (bin_center_bp, density_pct) tuples for polyline rendering
    """
    if not features or read_length <= 0:
        return []

    n_bins = (read_length + bin_size - 1) // bin_size
    bin_counts = [0] * n_bins

    # Count features per bin
    for feat in features:
        if feat.get('feature') == 'unannotated':
            continue
        start_bin = feat['start'] // bin_size
        end_bin = min((feat['stop'] - 1) // bin_size, n_bins - 1)
        if end_bin < 0:
            continue
        for b in range(max(0, start_bin), end_bin + 1):
            bin_counts[b] += 1

    # Convert to density percentages
    density_points = []
    for b in range(n_bins):
        bin_center = b * bin_size + bin_size // 2
        density_pct = min(100.0, (bin_counts[b] / max_density) * 100.0)
        density_points.append((bin_center, density_pct))

    return density_points

# =============================================================================
# Structural Analysis Plotting (Additive)
# =============================================================================

def draw_mini_dendrogram(d, linkage_matrix, x_base, y_base, width, height, row_y_centers, color="white", threshold=None, show_threshold=False):
    """Draw a small dendrogram to the left of feature bars."""
    n_leaves = len(row_y_centers)
    if n_leaves < 2:
        return
        
    max_dist = linkage_matrix[-1, 2] if len(linkage_matrix) > 0 else 1.0
    if max_dist == 0: max_dist = 1.0
    
    node_x = {}
    node_y = {}
    
    for i in range(n_leaves):
        node_x[i] = x_base + width
        node_y[i] = row_y_centers[i]
        
    for i, (idx1, idx2, dist, count) in enumerate(linkage_matrix):
        idx1, idx2 = int(idx1), int(idx2)
        node_idx = n_leaves + i
        y1, y2 = node_y[idx1], node_y[idx2]
        x_merge = x_base + width - (dist / max_dist) * width
        
        d.append(draw.Line(node_x[idx1], y1, x_merge, y1, stroke=color, stroke_width=2.5))
        d.append(draw.Line(node_x[idx2], y2, x_merge, y2, stroke=color, stroke_width=2.5))
        d.append(draw.Line(x_merge, y1, x_merge, y2, stroke=color, stroke_width=2.5))
        
        node_x[node_idx] = x_merge
        node_y[node_idx] = (y1 + y2) / 2

    # Draw threshold line
    if show_threshold and threshold is not None:
        x_thresh = x_base + width - (threshold / max_dist) * width
        if x_thresh >= x_base and x_thresh <= x_base + width:
            d.append(draw.Line(x_thresh, row_y_centers[0] - 20, x_thresh, row_y_centers[-1] + 20,
                              stroke="#FF4444", stroke_width=2, stroke_dasharray="5,5"))


def plot_structural_mode(args, matrix_data):
    """Generate multi-panel structural plot (separate files per chromosome)."""
    print("structural mode plotting...")
    
    import scipy.cluster.hierarchy as sch
    from scipy.spatial.distance import squareform
    
    # Load assignments
    reps_file = f"{args.cluster_prefix}.read_assignments.tsv"
    if not os.path.exists(reps_file):
        print(f"Error: {reps_file} not found.")
        sys.exit(1)
        
    reps_df = pd.read_csv(reps_file, sep='\t')
    if 'chromosome' not in reps_df.columns:
        print("Error: 'chromosome' column missing. Not a structural analysis output?")
        sys.exit(1)
        
    chromosomes = sorted(reps_df['chromosome'].unique(), key=lambda x: (len(x), x))
    print(f"Plotting {len(chromosomes)} chromosomes")
    
    # BED Data Source
    if args.bed_files:
        sample_bed_paths, database = parse_bed_paths(args.bed_files)
    elif args.input_bed_prefix:
        unique_samples = reps_df['sample'].unique()
        sample_bed_paths = {}
        database = args.database if args.database else "unknown_db"
        for s in unique_samples:
            path = os.path.join(args.input_bed_prefix, s, 'telogator', '1', 'KaryoScope', database)
            sample_bed_paths[s] = path
    else:
        sample_bed_paths = {}
        database = args.database

    featuresets = args.featuresets.split(',')
    featureset_colors, _ = load_color_files(args.colors_dir, database, featuresets)
    
    panel_width = 1200
    margin_x = 40
    margin_y = 60
    
    out_base = args.output[:-4] if args.output.endswith('.svg') else args.output
    
    chrom_drawings = [] # Store (drawing, height, width)
    padding = 60
    
    # Priority samples list
    priority_samples = []
    if getattr(args, 'priority_samples', None):
        priority_samples = [s.strip() for s in args.priority_samples.split(',')]
        
    for chrom in chromosomes:
        chrom_df = reps_df[reps_df['chromosome'] == chrom]
        unique_cids = chrom_df['cluster'].unique()
        
        # Helper to pick the best representative for a given cluster dataframe
        def pick_cluster_rep(cluster_subset, c_type):
            # First, check for priority samples
            for ps in priority_samples:
                match = cluster_subset[cluster_subset['sample'] == ps]
                if not match.empty: return match.iloc[0]
                match = cluster_subset[cluster_subset['read'].str.startswith(ps)]
                if not match.empty: return match.iloc[0]
                
            if c_type == 'Major':
                return cluster_subset.iloc[0]
                
            # For outliers, pick the one with highest raw_divergence
            return cluster_subset.sort_values('raw_divergence', ascending=False).iloc[0]

        cluster_reps = []
        for cid in unique_cids:
            c_data = chrom_df[chrom_df['cluster'] == cid]
            if c_data.empty: continue
            
            c_type = c_data.iloc[0]['cluster_type']
            rep = pick_cluster_rep(c_data, c_type)
            
            cluster_reps.append({
                'read': rep['read'],
                'cluster': cid,
                'type': c_type,
                'sample': rep['sample'] if 'sample' in rep else 'unknown',
                'raw_div': rep['raw_divergence'] if 'raw_divergence' in rep else 0.0,
                'norm_div': rep['norm_divergence'] if 'norm_divergence' in rep else 0.0
            })
            
        # Sort reps: Major first, then Outliers by raw_div descending
        cluster_reps.sort(key=lambda x: (0 if x['type'] == 'Major' else 1, -x['raw_div']))
        
        selected_reads = cluster_reps
                    
        if not selected_reads: continue

        if not getattr(args, 'show_dendrogram', False):
            selected_reads.sort(key=lambda x: (0 if x['type'] == 'Major' else 1))

        fs_height = 14
        row_spacing = 40
        canvas_height = 120 + len(selected_reads) * (len(featuresets) * fs_height + row_spacing) + 120
        canvas_width = panel_width + 2 * margin_x
        
        d = draw.Drawing(canvas_width, canvas_height, displayInline=False)
        bg_color = args.background_color
        text_color = "white" if bg_color == "black" else "black"
        d.append(draw.Rectangle(0, 0, canvas_width, canvas_height, fill=bg_color))
        d.append(draw.Text(f"KaryoScope: {chrom} Structural Analysis", 24, canvas_width/2, 40, 
                          fill=text_color, font_weight='bold', text_anchor='middle'))
        
        panel_bg = "#111111" if bg_color == "black" else "#F5F5F5"
        d.append(draw.Rectangle(margin_x - 5, 75, panel_width + 10, canvas_height - 180, fill=panel_bg, rx=10, ry=10))

        reads_needed = set(r['read'] for r in selected_reads)
        read_bed_data = load_bed_data(sample_bed_paths, database, featuresets, args.smoothness, reads_needed)
        
        # Local clustering
        local_Z = None
        if getattr(args, 'show_dendrogram', False) and len(selected_reads) > 1:
            unique_f = sorted(list(set(f['feature'] for r in read_bed_data for fs in read_bed_data[r] for f in read_bed_data[r][fs])))
            f_map = {f: chr(j+200) for j, f in enumerate(unique_f)}
            encoded = []
            for r_obj in selected_reads:
                r = r_obj['read']
                if r in read_bed_data:
                    fs0 = featuresets[0]
                    feats = sorted(read_bed_data[r].get(fs0, []), key=lambda x: x['start'])
                    encoded.append("".join([f_map.get(f['feature'], '?') for f in feats]))
                else: encoded.append("")

            dm = np.zeros((len(encoded), len(encoded)))
            def _lev(s1, s2):
                if len(s1) < len(s2): return _lev(s2, s1)
                if not s2: return len(s1)
                p = range(len(s2) + 1)
                for c in s1:
                    cur = [p[0]+1]
                    for j, c2 in enumerate(s2): cur.append(min(p[j+1]+1, cur[j]+1, p[j]+(c!=c2)))
                    p = cur
                return p[-1]

            for j1 in range(len(encoded)):
                for j2 in range(j1+1, len(encoded)):
                    d_val = _lev(encoded[j1], encoded[j2])
                    dm[j1, j2] = dm[j2, j1] = d_val / max(len(encoded[j1]), len(encoded[j2]), 1)
            
            local_Z = sch.linkage(squareform(dm), method='ward')
            leaf_order = sch.leaves_list(local_Z)
            selected_reads = [selected_reads[idx] for idx in leaf_order]

        max_len = 0
        for r_obj in selected_reads:
            r = r_obj['read']
            if r in read_bed_data:
                for fs in read_bed_data[r]:
                    for feat in read_bed_data[r][fs]: max_len = max(max_len, feat['stop'])
        if max_len == 0: max_len = 10000 
        
        show_d = getattr(args, 'show_dendrogram', False)
        dendro_w = 250 if show_d else 20
        label_w = 400 # Slightly wider for detailed labels
        bars_x_start = margin_x + (dendro_w if show_d else 0) + label_w + 20
        ratio = (panel_width - (dendro_w if show_d else 0) - label_w - 60) / max_len
        
        row_y_centers = []
        for j in range(len(selected_reads)):
            row_y_centers.append(110 + j * (len(featuresets)*fs_height + row_spacing) + (len(featuresets)*fs_height)/2)

        if show_d and local_Z is not None:
            draw_mini_dendrogram(d, local_Z, margin_x + 10, 110, dendro_w - 20, 0, row_y_centers, 
                                 color=text_color, threshold=getattr(args, 'structural_threshold', 0.25), 
                                 show_threshold=getattr(args, 'show_threshold', False))
            
        for j, r_obj in enumerate(selected_reads):
            read, ry = r_obj['read'], row_y_centers[j] - (len(featuresets)*fs_height)/2
            l_color = "#888888" if r_obj['type'] == "Major" else "#FF4444"
            
            # Detailed label
            sample_val = r_obj.get('sample', 'unknown')
            if sample_val == 'pangenome' and '#' in read:
                sample_val = read.split('#')[0]
                
            clean_cid = r_obj['cluster'].replace(f"{chrom}_", "")
            raw_div = r_obj.get('raw_div', 0)
            norm_div = r_obj.get('norm_div', 0)
            
            label_text = f"[{sample_val}] {clean_cid}"
            if r_obj['type'] != "Major":
                label_text += f" (raw:{int(raw_div)}, norm:{norm_div:.2f})"
                
            display_name = label_text if len(label_text) <= 65 else label_text[:35] + "..." + label_text[-25:]
            d.append(draw.Text(display_name, 12, margin_x + (dendro_w if show_d else 10), ry + (len(featuresets)*fs_height)/2 + 4, fill=l_color, font_family='monospace'))
            if read in read_bed_data:
                for fs_idx, fs in enumerate(featuresets):
                    for feat in read_bed_data[read].get(fs, []):
                        w = max((feat['stop'] - feat['start']) * ratio, 2.5)
                        x = bars_x_start + (feat['start'] * ratio)
                        color, op = featureset_colors[fs].get(feat['feature'], ("#ffffff", 1.0))
                        d.append(draw.Rectangle(x, ry + (fs_idx * fs_height), w, fs_height, fill=color, fill_opacity=op))

        legend_y = canvas_height - 60
        d.append(draw.Text("Legend:", 14, margin_x, legend_y, fill=text_color, font_weight='bold'))
        d.append(draw.Text("Major (Dominant)", 12, margin_x + 100, legend_y, fill="#888888"))
        d.append(draw.Text("Outlier (Variant)", 12, margin_x + 300, legend_y, fill="#FF4444"))
        if featuresets: d.append(draw.Text(f"Tracks: {', '.join(featuresets)}", 12, margin_x + 500, legend_y, fill=text_color))
        d.save_svg(f"{out_base}.{chrom}.svg")
        
        # Store for combined grid
        chrom_drawings.append((d, canvas_height, canvas_width))

    if chrom_drawings:
        print(f"Generating combined structural grid plot: {out_base}.all_chromosomes.svg")
        n_cols = 5
        
        # Calculate row heights
        row_heights = []
        for i in range(0, len(chrom_drawings), n_cols):
            batch = chrom_drawings[i:i+n_cols]
            row_heights.append(max(h for d, h, w in batch))
            
        total_width = n_cols * (canvas_width + padding) + padding
        total_height = sum(row_heights) + (len(row_heights) + 1) * padding
        
        all_d = draw.Drawing(total_width, total_height, displayInline=False)
        all_d.append(draw.Rectangle(0, 0, total_width, total_height, fill=args.background_color))
        
        current_y = padding
        for r_idx in range(len(row_heights)):
            current_x = padding
            for c_idx in range(n_cols):
                idx = r_idx * n_cols + c_idx
                if idx < len(chrom_drawings):
                    d_obj, d_h, d_w = chrom_drawings[idx]
                    g = draw.Group(transform=f"translate({current_x},{current_y})")
                    for elem in d_obj.elements:
                        g.append(elem)
                    all_d.append(g)
                    current_x += d_w + padding
            current_y += row_heights[r_idx] + padding
            
        all_d.save_svg(f"{out_base}.all_chromosomes.svg")
        
    print("Finished generating structural plots.")
    sys.exit(0)


def main():
    args = parse_args()

    # --- Set up logging ---
    if args.log_file:
        # Derive log path from output file (replace .svg with .log)
        if args.output.endswith('.svg'):
            log_path = args.output[:-4] + '.log'
        else:
            log_path = args.output + '.log'
        sys.stdout = TeeLogger(log_path)

    # --- Auto-discover files from prefix ---
    prefix = args.cluster_prefix
    representatives_file = f"{prefix}.read_assignments.tsv"
    feature_matrix_file = f"{prefix}.feature_matrix.npz"
    sample_metadata_file = f"{prefix}.sample_metadata.tsv"
    cluster_analysis_file = f"{prefix}.cluster_analysis.tsv"

    # Check for structural mode
    if os.path.exists(feature_matrix_file):
        is_structure = False
        fm = None
        try:
             fm = np.load(feature_matrix_file, allow_pickle=True)
             if 'mode' in fm and str(fm['mode']) == 'structure':
                 is_structure = True
        except Exception:
             pass
             
        if is_structure:
             plot_structural_mode(args, fm)
             sys.exit(0)

    # Verify required files exist
    if not os.path.exists(representatives_file):
        sys.stderr.write(f"Error: Read assignments file not found: {representatives_file}\n")
        sys.exit(1)

    # --- Parse BED paths to get sample directories and database ---
    if args.bed_files:
        # Use explicit BED file paths
        sample_bed_paths, database = parse_bed_paths(args.bed_files)
        if not database:
            sys.stderr.write("Error: Could not determine database from BED file paths\n")
            sys.exit(1)
    elif args.input_bed_prefix:
        # Auto-discover from sample metadata + input prefix
        if not args.database:
            sys.stderr.write("Error: --database is required when using --input-bed-prefix\n")
            sys.exit(1)
        database = args.database

        # Load sample names from metadata
        if not os.path.exists(sample_metadata_file):
            sys.stderr.write(f"Error: Sample metadata file not found: {sample_metadata_file}\n")
            sys.exit(1)

        meta_df = pd.read_csv(sample_metadata_file, sep='\t')
        sample_names = meta_df['sample'].tolist()

        # Build sample_bed_paths from input prefix
        sample_bed_paths = {}
        print(f"\nAuto-discovering BED paths from --input-bed-prefix...")
        for sample in sample_names:
            bed_dir = os.path.join(args.input_bed_prefix, sample, 'telogator', '1', 'KaryoScope', database)
            if os.path.exists(bed_dir):
                sample_bed_paths[sample] = bed_dir
                print(f"  {sample} -> {bed_dir}")
            else:
                print(f"  Warning: Directory not found for {sample}: {bed_dir}")

        if not sample_bed_paths:
            sys.stderr.write("Error: No valid BED directories found\n")
            sys.exit(1)
    else:
        sys.stderr.write("Error: Either --bed or --input-bed-prefix is required\n")
        sys.exit(1)

    # Override database if explicitly provided
    if args.database:
        database = args.database

    # --- Setup ---
    background_color = args.background_color
    text_color = "#000000" if background_color == "white" else "#FFFFFF"
    featuresets = [f.strip() for f in args.featuresets.split(",")]

    # Parse custom BED files
    custom_bed_files = {}  # featureset_name -> path
    if args.custom_beds:
        for item in args.custom_beds:
            if ':' not in item:
                print(f"Error: Custom BED format should be 'name:path', got '{item}'")
                sys.exit(1)
            name, path = item.split(':', 1)
            custom_bed_files[name] = path
            if name not in featuresets:
                featuresets.append(name)

    # Process --fiberseq directory (auto-discover fiberseq BED files)
    fiberseq_files = {}  # feature_type -> path
    if args.fiberseq_dir:
        fiberseq_dir = args.fiberseq_dir
        print(f"\nAuto-discovering fiberseq BED files in: {fiberseq_dir}")

        # Look for each fiberseq feature type
        for feature_type in ['FIRE', 'LINKER', 'm6A', '5mC']:
            pattern = os.path.join(fiberseq_dir, f"*.{feature_type}.bed")
            matches = glob.glob(pattern)
            if matches:
                fiberseq_files[feature_type] = matches[0]
                print(f"  Found {feature_type}: {os.path.basename(matches[0])}")

        # Create combined FIRE_LINKER file if both exist
        if 'FIRE' in fiberseq_files and 'LINKER' in fiberseq_files:
            combined_path = os.path.join(fiberseq_dir, os.path.basename(fiberseq_files['FIRE']).replace('.FIRE.bed', '.FIRE_LINKER.bed'))
            if not os.path.exists(combined_path):
                print(f"  Creating combined FIRE_LINKER file...")
                with open(combined_path, 'w') as out:
                    for bed_file in [fiberseq_files['FIRE'], fiberseq_files['LINKER']]:
                        with open(bed_file) as f:
                            for line in f:
                                out.write(line)
            custom_bed_files['fiberseq_FIRE_LINKER'] = combined_path
            if 'fiberseq_FIRE_LINKER' not in featuresets:
                featuresets.append('fiberseq_FIRE_LINKER')
            print(f"  Added fiberseq_FIRE_LINKER track")

        # Add m6A and 5mC as custom beds for density line plot
        if 'm6A' in fiberseq_files:
            custom_bed_files['fiberseq_m6A'] = fiberseq_files['m6A']
        if '5mC' in fiberseq_files:
            custom_bed_files['fiberseq_5mC'] = fiberseq_files['5mC']

        # Auto-set density_line_plot if m6A and 5mC are available and not already set
        if 'm6A' in fiberseq_files and '5mC' in fiberseq_files and not args.density_line_plot:
            args.density_line_plot = 'fiberseq_m6A:fiberseq_5mC'
            print(f"  Auto-configured density line plot for m6A/5mC")

    num_featuresets = len(featuresets)

    # Parse density featuresets
    density_featuresets = set()
    if args.density_featuresets:
        density_featuresets = set(f.strip() for f in args.density_featuresets.split(','))
        print(f"Density featuresets: {density_featuresets}")
        print(f"Density bin size: {args.density_bin_size} bp")

    # Parse density line plot featuresets
    density_line_plot_featuresets = []
    density_line_plot_name = None
    if args.density_line_plot:
        density_line_plot_featuresets = [f.strip() for f in args.density_line_plot.split(':')]
        density_line_plot_name = "density_line"
        # Remove individual featuresets from main list and add combined track
        featuresets = [fs for fs in featuresets if fs not in density_line_plot_featuresets]
        featuresets.append(density_line_plot_name)
        print(f"Density line plot featuresets: {density_line_plot_featuresets}")
        print(f"Density bin size: {args.density_bin_size} bp")

    # Parse rect plot featuresets (exact feature rectangles)
    rect_plot_featuresets = []
    rect_plot_name = None
    if args.rect_plot:
        rect_plot_featuresets = [f.strip() for f in args.rect_plot.split(':')]
        rect_plot_name = "rect_plot"
        # Remove individual featuresets from main list and add combined track
        featuresets = [fs for fs in featuresets if fs not in rect_plot_featuresets]
        featuresets.append(rect_plot_name)
        print(f"Rect plot featuresets: {rect_plot_featuresets}")

    num_featuresets = len(featuresets)

    # Featureset display names
    fs_display_names = {
        "chromosome": "Chromosome",
        "subtelomeric": "Subtelomere",
        "region": "Satellite",
        "acrocentric": "Acrocentric",
        "repeat": "Interspersed repeat",
        "fiberseq": "Fiberseq",
        "fiberseq_FIRE": "FIRE",
        "fiberseq_m6A": "m6A",
        "fiberseq_5mC": "5mC",
        "fiberseq_LINKER": "Linker",
        "fiberseq_FIRE_LINKER": "FIRE/Linker",
        "density_line": "m6A/5mC",
        "rect_plot": "FIRE/Linker"
    }

    max_reps = args.max_reps

    print(f"Feature sets to plot: {featuresets}")
    print(f"Background color: {background_color}")

    # --- Load data ---
    # Load sample metadata
    sample_to_group, sample_colors, group_colors = load_sample_metadata(sample_metadata_file)

    # Load cluster analysis to get enrichment info, cluster priority order, stats for bubbles, and full df for features
    cluster_enrichments, cluster_order, cluster_stats, cluster_analysis_df = load_cluster_analysis(cluster_analysis_file)

    # Load custom cluster labels if provided
    cluster_labels = load_cluster_labels(args.cluster_labels, args.label_column)

    # Load read assignments (filtering via --curated-reps or --reads-file)
    cluster_reads, unique_enrichments = load_representative_reads(
        representatives_file,
        cluster_enrichments=cluster_enrichments,
        cluster_order=cluster_order,
        max_reps=max_reps,
        reads_file=args.reads_file,
        curated_reps_file=args.curated_reps,
        cluster_labels_file=args.cluster_labels
    )

    # Load feature matrix
    feature_matrix_data = load_feature_matrix(feature_matrix_file)

    # Load color files - include density line plot and rect plot featuresets for color lookup
    # Exclude placeholder tracks from color loading (they have no color files)
    placeholder_tracks = {density_line_plot_name, rect_plot_name} - {None}
    featuresets_for_colors = list(set(
        [fs for fs in featuresets if fs not in placeholder_tracks] +
        density_line_plot_featuresets +
        rect_plot_featuresets
    ))

    featureset_colors, featureset_color_order = load_color_files(
        args.colors_dir, database, featuresets_for_colors
    )

    # Add empty placeholders for combined tracks (they don't have individual color files)
    for track_name in placeholder_tracks:
        featureset_colors[track_name] = {}
        featureset_color_order[track_name] = []

    # Get all reads we need
    all_reads_needed = set()
    read_to_sample = {}
    for cluster_id, data in cluster_reads.items():
        for read, sample in data['reads']:
            all_reads_needed.add(read)
            read_to_sample[read] = sample

    # Load BED data
    # Filter out custom featuresets for standard loading
    standard_featuresets = [fs for fs in featuresets if fs not in custom_bed_files]
    read_data = load_bed_data(
        sample_bed_paths, database, standard_featuresets, args.smoothness, all_reads_needed
    )

    # Load custom BED files
    read_data = load_custom_bed_files(custom_bed_files, all_reads_needed, read_data)

    # --- Generate colors ---
    # Get all unique samples
    all_samples = sorted(set(sample for data in cluster_reads.values() for _, sample in data['reads']))
    sample_colors = generate_sample_colors(all_samples, sample_colors)

    # Generate enrichment colors from group/sample colors
    enrichment_colors = get_enrichment_colors(group_colors, unique_enrichments, sample_colors)

    # --- Compute cluster-level dendrogram order if feature matrix provided ---
    cluster_dendro_data = None
    read_to_original_cluster = {}
    read_to_original_enrichment = {}

    if feature_matrix_data is not None and not args.no_reorder:
        cluster_reads, cluster_dendro_data, read_to_original_cluster, read_to_original_enrichment = \
            compute_cluster_dendrogram_order(feature_matrix_data, cluster_reads)
    elif args.no_reorder:
        print("  Dendrogram reordering disabled (--no-reorder)")

    # Hide dendrogram if requested (but still use ordering)
    if args.hide_dendrogram:
        cluster_dendro_data = None
        print("  Dendrogram hidden (--hide-dendrogram)")

    # Compute full dendrogram if requested
    full_dendro_data = None
    if getattr(args, 'full_dendrogram', False) and feature_matrix_data is not None:
        # Collect all displayed reads
        all_displayed_reads = []
        for cluster_id, data in cluster_reads.items():
            for read, sample in data['reads']:
                all_displayed_reads.append(read)
        full_dendro_data = compute_full_dendrogram(feature_matrix_data, all_displayed_reads)

    # ==========================================================================
    # VERTICAL MODE - Separate drawing path
    # ==========================================================================
    if args.vertical:
        # Build mappings if not already done
        if not read_to_original_cluster:
            for cluster_id, data in cluster_reads.items():
                for read, sample in data['reads']:
                    read_to_original_cluster[read] = cluster_id
                    read_to_original_enrichment[read] = data['enrichment']

        unique_clusters = set(read_to_original_cluster.values())
        cluster_colors = get_cluster_colors(unique_clusters)

        # Vertical layout parameters
        dendrogram_width = 100 if cluster_dendro_data is not None else 0
        bubble_radius = 8
        dendrogram_to_bubble_gap = 1  # Gap from dendrogram tip to bubble left edge
        bubble_to_bars_gap = 12  # Gap from bubble right edge to feature bars
        bubble_space = dendrogram_to_bubble_gap + bubble_radius * 2 + bubble_to_bars_gap
        left_margin = 50 + dendrogram_width + bubble_space
        sample_dendro_height = 40 if args.show_matrix else 0  # Space for sample dendrogram
        top_margin = 80 + sample_dendro_height  # More space for rotated sample headers + dendrogram
        bar_width = args.bar_width
        num_fs = len(featuresets)
        column_spacing = 10  # Space between columns in column mode

        # In column mode, each read only needs one bar height
        # In row mode (default), each read has all featuresets stacked
        if args.column_tracks:
            group_height = bar_width
        else:
            group_height = bar_width * num_fs + args.bar_spacing * (num_fs - 1)

        # Use compact spacing when matrix is shown or column-tracks mode
        if args.show_matrix or args.column_tracks:
            row_spacing = 2  # Minimal spacing between reads
            cluster_gap = 4  # Small gap between clusters
        else:
            row_spacing = args.read_spacing
            cluster_gap = args.cluster_spacing

        # Calculate y positions (reads stacked vertically)
        read_y_positions = {}
        cluster_y_start = {}
        cluster_y_end = {}
        current_y = top_margin

        # If full dendrogram is available, use its leaf order instead of cluster order
        if full_dendro_data is not None and 'read_order' in full_dendro_data:
            # Build reverse mapping: read -> cluster_id
            read_to_cluster = {}
            for cluster_id, data in cluster_reads.items():
                for read, sample in data['reads']:
                    read_to_cluster[read] = cluster_id

            # Use dendrogram leaf order with optional cluster gap
            ordered_reads = full_dendro_data['read_order']
            dendro_cluster_gap = getattr(args, 'dendro_cluster_gap', 0)
            prev_cluster = None
            for read in ordered_reads:
                curr_cluster = read_to_cluster.get(read)
                # Add extra gap when crossing cluster boundary
                if dendro_cluster_gap > 0 and prev_cluster is not None and curr_cluster != prev_cluster:
                    current_y += dendro_cluster_gap
                read_y_positions[read] = current_y
                current_y += group_height + row_spacing
                prev_cluster = curr_cluster
            # Set cluster boundaries to span all reads (single group)
            cluster_y_start[1] = top_margin
            cluster_y_end[1] = current_y - row_spacing
        else:
            # Original cluster-based ordering
            for cluster_id, data in cluster_reads.items():
                cluster_y_start[cluster_id] = current_y
                for read, sample in data['reads']:
                    read_y_positions[read] = current_y
                    current_y += group_height + row_spacing
                cluster_y_end[cluster_id] = current_y - row_spacing
                current_y += cluster_gap

        # Calculate scaffold lengths and max length for uniform bar width
        max_read_length = 0
        scaffold_lengths = {}
        scaffold_min_starts = {}

        for read in read_data:
            for fs in read_data[read]:
                for feat in read_data[read][fs]:
                    scaffold_lengths[read] = max(scaffold_lengths.get(read, 0), feat['stop'])
                    if read not in scaffold_min_starts:
                        scaffold_min_starts[read] = feat['start']
                    else:
                        scaffold_min_starts[read] = min(scaffold_min_starts[read], feat['start'])
            if read in scaffold_lengths:
                max_read_length = max(max_read_length, scaffold_lengths[read] - scaffold_min_starts.get(read, 0))

        # Calculate ratio - either from args or auto-calculated from target dimensions
        ratio = args.ratio
        if args.target_width is not None and max_read_length > 0:
            # For vertical mode, target_width controls bar length
            # Estimate margins: left_margin + label_width + right_legend_width ≈ 400
            estimated_margins = 400
            target_bar_length = args.target_width - estimated_margins
            if target_bar_length > 0:
                ratio = target_bar_length / max_read_length
                print(f"  Auto-calculated ratio from target_width: {ratio:.6f}")
        max_bar_length = floor(max_read_length * ratio)

        # Build drawing data for vertical mode
        drawing_data_vertical = defaultdict(lambda: defaultdict(list))
        uncolored_features = defaultdict(set)
        displayed_features = defaultdict(set)  # Track actually displayed features for legend filtering

        for read in read_data:
            if read not in read_y_positions:
                continue

            base_y = read_y_positions[read]
            scaffold_min_start = scaffold_min_starts.get(read, 0)

            for fs in featuresets:
                features = read_data[read].get(fs, [])

                for feat in features:
                    final_start = feat['start'] - scaffold_min_start
                    final_stop = feat['stop'] - scaffold_min_start
                    scaled_start = floor(final_start * ratio)
                    scaled_stop = floor(final_stop * ratio)

                    feature_name = feat.get('feature', 'unknown')
                    displayed_features[fs].add(feature_name)  # Track displayed feature
                    color_info = featureset_colors.get(fs, {}).get(feature_name)
                    if color_info is None:
                        color = '#444444'
                        fill_opacity = 1.0
                        uncolored_features[fs].add(feature_name)
                    else:
                        color, fill_opacity = color_info

                    drawing_data_vertical[read][fs].append({
                        'scaled_start': scaled_start,
                        'scaled_stop': scaled_stop,
                        'color': color,
                        'fill_opacity': fill_opacity
                    })

        # Calculate read positions dict for vertical drawing: (x_start, y_start, bar_length)
        read_positions_vertical = {}
        for read in read_y_positions:
            y_pos = read_y_positions[read]
            read_length = scaffold_lengths.get(read, 0) - scaffold_min_starts.get(read, 0)
            bar_length = floor(read_length * ratio)
            read_positions_vertical[read] = (left_margin, y_pos, bar_length)

        # Calculate feature bar width (changes in column mode)
        if args.column_tracks:
            # In column mode: total width = num_featuresets * (max_bar_length + spacing)
            feature_bars_width = num_fs * (max_bar_length + column_spacing) - column_spacing
        else:
            # In row mode: total width = max_bar_length
            feature_bars_width = max_bar_length

        # Calculate matrix parameters if enabled
        matrix_width = 0
        matrix_x_start = left_margin + feature_bars_width + 15  # After feature bars
        # Square cells: size matches the row height (group_height)
        cell_size = group_height + row_spacing  # Square cells matching row spacing
        cell_width = cell_size
        cell_height = cell_size

        if args.show_matrix:
            # Load sample metadata for matrix
            meta_df = pd.read_csv(sample_metadata_file, sep='\t')
            n_samples = len(meta_df)
            matrix_width = n_samples * cell_width + 15  # Small padding
            print(f"  Matrix enabled: {n_samples} samples × {len(cluster_y_start)} clusters (cell size: {cell_size}px)")

        # Image dimensions
        # Reduce label width and skip bubble legend when using full dendrogram
        if full_dendro_data is not None:
            label_width = 120  # Smaller space for read labels
            bubble_legend_height = 30  # Minimal bottom margin
        else:
            label_width = 200  # Space for cluster labels
            bubble_legend_height = 150  # Space for bubble legend + enrichment text legend at bottom
        right_legend_width = 120  # Space for vertical color legend on right
        bar_plot_height = 100 if args.show_matrix else 0  # Space for bar plot below matrix
        row_bar_width = 60 if args.show_matrix else 0  # Space for row barplot to right of matrix
        image_width = left_margin + feature_bars_width + 20 + matrix_width + row_bar_width + label_width + right_legend_width
        image_height = current_y + 50 + bar_plot_height + bubble_legend_height

        if args.column_tracks:
            print(f"\nVertical mode (column tracks): {num_fs} columns × {image_width} x {image_height}")
        else:
            print(f"\nVertical mode image dimensions: {image_width} x {image_height}")

        # Create drawing
        d = draw.Drawing(image_width, image_height)
        d.append(draw.Rectangle(0, 0, image_width, image_height, fill=background_color))

        # Draw vertical dendrogram on left
        if full_dendro_data is not None:
            # Full dendrogram showing all individual reads
            read_names_displayed = [read for cluster_id, data in cluster_reads.items()
                                    for read, sample in data['reads']]
            draw_full_dendrogram(d, full_dendro_data, read_y_positions, read_names_displayed,
                                 left_margin, dendrogram_width, background_color)
        elif cluster_dendro_data is not None:
            draw_cluster_dendrogram_vertical(d, cluster_dendro_data, cluster_y_start, cluster_y_end,
                                             left_margin, dendrogram_width, background_color)

        # Draw scale bar above first featureset column
        if args.column_tracks:
            scale_bar_y = min(cluster_y_start.values()) - 15  # Position above first cluster
            draw_scale_bar(d, left_margin, scale_bar_y, ratio, text_color)

        # Draw feature bars
        if args.column_tracks:
            # Column mode: each featureset in its own column
            draw_feature_bars_column_mode(d, drawing_data_vertical, featuresets, bar_width,
                                          read_y_positions, left_margin, max_bar_length,
                                          column_spacing, background_color)
        else:
            # Row mode (default): featuresets stacked as rows within each read
            draw_feature_bars_vertical(d, drawing_data_vertical, featuresets, bar_width,
                                       read_positions_vertical, num_fs, args.bar_spacing, background_color)

        # Draw read index labels (1, 2, 3, ...) after feature bars (if enabled)
        if args.show_read_indices:
            read_index_x = left_margin + feature_bars_width + 5
            draw_read_index_labels(d, cluster_reads, read_y_positions, read_index_x, bar_width, text_color)

        # Draw read name labels with cluster colors when using full dendrogram
        if full_dendro_data is not None and 'read_order' in full_dendro_data:
            # Build reverse mapping: read -> cluster_id
            read_to_cluster = {}
            for cluster_id, data in cluster_reads.items():
                for read, sample in data['reads']:
                    read_to_cluster[read] = cluster_id

            # Define cluster colors (cycle through a palette)
            cluster_palette = ['#40D392', '#60A5FA', '#F07167', '#FBBF24', '#C4A9E8',
                              '#10B981', '#3B82F6', '#EF4444', '#F59E0B', '#8B5CF6']
            cluster_ids_sorted = sorted(set(read_to_cluster.values()))
            cluster_color_map = {cid: cluster_palette[i % len(cluster_palette)]
                                for i, cid in enumerate(cluster_ids_sorted)}

            label_x = left_margin + feature_bars_width + 8
            for read in full_dendro_data['read_order']:
                if read not in read_y_positions:
                    continue
                y_pos = read_y_positions[read] + group_height / 2
                cluster_id = read_to_cluster.get(read, 0)
                cluster_color = cluster_color_map.get(cluster_id, text_color)

                # Draw small cluster color indicator
                d.append(draw.Circle(label_x + 4, y_pos, 3, fill=cluster_color))

                # Draw read name (abbreviated)
                short_name = abbreviate_read_name(read)
                d.append(draw.Text(
                    short_name, font_size=8, x=label_x + 12, y=y_pos + 3,
                    fill=text_color, font_family='monospace', text_anchor='start'
                ))

        # Draw sample matrix if enabled
        matrix_data = None
        if args.show_matrix:
            cluster_ids = list(cluster_y_start.keys())
            matrix_data = draw_sample_matrix(d, cluster_ids, cluster_y_start, cluster_y_end, meta_df,
                              representatives_file, matrix_x_start, cell_width, cell_height,
                              text_color, background_color)

            # Draw sample dendrogram just above sample labels
            header_y = min(cluster_y_start.values()) - 5
            dendro_bottom = header_y - 30  # Above sample name labels with spacing
            draw_sample_dendrogram(d, matrix_data, matrix_x_start, dendro_bottom, sample_dendro_height)

            # Draw bar plot below matrix (column sums)
            bar_plot_y = max(cluster_y_end.values()) + cell_height / 2 + 5
            draw_sample_bar_plot(d, matrix_data, cluster_ids, cluster_enrichments, matrix_x_start, bar_plot_y,
                                cell_width, 40, text_color, background_color)

            # Draw row bar plot to right of matrix (row sums)
            row_bar_x_start = matrix_x_start + matrix_width + 5
            draw_cluster_bar_plot(d, matrix_data, cluster_ids, cluster_y_start, cluster_y_end,
                                 cluster_enrichments, row_bar_x_start, row_bar_width - 10,
                                 text_color, background_color)

        # Draw enrichment bubbles (to the RIGHT of dendrogram tips, before feature bars)
        # Dendrogram tips are at 50 + dendrogram_width (consistent with draw_cluster_dendrogram_vertical)
        dendro_tip_x = 50 + dendrogram_width
        bubble_x = dendro_tip_x + dendrogram_to_bubble_gap + bubble_radius

        # Skip enrichment bubbles and cluster labels when using full dendrogram (individual read view)
        if full_dendro_data is None:
            # Draw faint gray connecting lines from dendrogram tips to bubble centers (BEFORE bubbles so bubbles overlay)
            connector_color = '#555555'  # Faint gray, lighter than white tree
            for cluster_id in cluster_y_start:
                y_start = cluster_y_start[cluster_id]
                y_end = cluster_y_end[cluster_id]
                y_center = (y_start + y_end) / 2
                # Line from dendrogram tip to bubble center (bubble will be drawn on top)
                d.append(draw.Line(
                    dendro_tip_x, y_center,
                    bubble_x, y_center,
                    stroke=connector_color, stroke_width=0.5
                ))

            # Draw bubbles on top of connector lines
            draw_enrichment_bubbles(d, cluster_y_start, cluster_y_end, bubble_x, cluster_stats,
                                    max_radius=bubble_radius, min_radius=2)

            # Draw cluster labels on right (after feature bars and read indices)
            if not args.hide_brackets:
                if args.show_matrix:
                    label_x = matrix_x_start + matrix_width + row_bar_width + 10
                else:
                    label_x = left_margin + feature_bars_width + 20  # Extra space for read index labels
                # Build cluster_enrichments dict for coloring labels
                cluster_enrichments_dict = {cid: data['enrichment'] for cid, data in cluster_reads.items() if cid != 'all'}
                draw_cluster_labels_vertical(d, cluster_y_start, cluster_y_end, label_x, text_color,
                                             cluster_labels=cluster_labels,
                                             enrichment_colors=enrichment_colors,
                                             cluster_enrichments=cluster_enrichments_dict)

            # Draw bubble legend at bottom
            legend_y = max(cluster_y_end.values()) + bar_plot_height + 30
            draw_bubble_legend(d, left_margin, legend_y, cluster_stats, text_color,
                              max_radius=bubble_radius, min_radius=2)

            # Draw enrichment text color legend below bubble legend
            enrichment_legend_y = legend_y + 75  # Below bubble legend (bubble legend is ~65px tall)
            draw_enrichment_text_legend(d, left_margin, enrichment_legend_y, enrichment_colors, text_color)

        # Draw matrix color legend if matrix is enabled (to the right of enrichment legend)
        if args.show_matrix and matrix_data:
            matrix_legend_x = left_margin + 400  # Position to the right of enrichment legend
            draw_matrix_legend(d, matrix_legend_x, enrichment_legend_y, matrix_data['max_count'],
                              text_color, background_color)

        # Draw featureset color legends vertically on the right side
        color_legend_x = image_width - right_legend_width + 10
        color_legend_y_start = top_margin
        draw_color_legends_vertical(d, featuresets, featureset_colors, featureset_color_order,
                                    fs_display_names, color_legend_x, color_legend_y_start, text_color,
                                    displayed_features=displayed_features)

        # Save vertical plot
        d.save_svg(args.output)

        # Report warnings
        for fs in featuresets:
            if uncolored_features[fs]:
                sys.stderr.write(f"Warning: {fs} - features not in colors file:\n")
                for feature in sorted(list(uncolored_features[fs])):
                    sys.stderr.write(f"  - {feature}\n")

        print(f"\n--- Summary ---")
        original_cluster_count = len(set(read_to_original_cluster.values())) if read_to_original_cluster else len(cluster_reads)
        print(f"Clusters plotted: {original_cluster_count}")
        total_reads = sum(len(data['reads']) for data in cluster_reads.values())
        print(f"Total reads plotted: {total_reads}")
        print(f"\n✅ Saved to {args.output}")

        # Print parameters
        params = [
            ("cluster-analysis-prefix", args.cluster_prefix),
            ("output", args.output),
            ("featuresets", ','.join(featuresets)),
            ("smoothness", args.smoothness),
            ("background", background_color),
            ("reads-file", args.reads_file if args.reads_file else "None"),
            ("curated-reps", args.curated_reps if args.curated_reps else "None"),
            ("show-read-indices", args.show_read_indices),
            ("hide-brackets", args.hide_brackets),
            ("hide-dendrogram", args.hide_dendrogram),
            ("cluster-labels", args.cluster_labels if args.cluster_labels else "None"),
            ("vertical", args.vertical),
            ("show-matrix", args.show_matrix),
            ("column-tracks", args.column_tracks),
            ("n-per-cluster", args.max_reps if args.max_reps else "None"),
        ]
        print(f"\n{'='*60}")
        print("Parameters")
        print(f"{'='*60}")
        print(f"{'Parameter':<25} {'Value':<35}")
        print(f"{'-'*25} {'-'*35}")
        for param, value in params:
            print(f"{param:<25} {str(value):<35}")

        return  # Exit main after vertical mode

    # ==========================================================================
    # HORIZONTAL MODE (original code continues below)
    # ==========================================================================

    # If not reordered, build mappings
    if not read_to_original_cluster:
        for cluster_id, data in cluster_reads.items():
            for read, sample in data['reads']:
                read_to_original_cluster[read] = cluster_id
                read_to_original_enrichment[read] = data['enrichment']

    # Generate cluster colors
    unique_clusters = set(read_to_original_cluster.values())
    cluster_colors = get_cluster_colors(unique_clusters)

    # --- Calculate positions ---
    group_width = (args.bar_width * num_featuresets) + (args.bar_spacing * (num_featuresets - 1))
    left_margin = 150
    # Show dendrogram space if we have cluster or full dendrogram data
    has_dendrogram = cluster_dendro_data is not None or full_dendro_data is not None
    dendrogram_height = 100 if has_dendrogram else 0
    bracket_height = 0 if args.hide_brackets else 50
    top_margin = 100 + dendrogram_height + bracket_height

    # Calculate x positions
    read_x_positions = {}
    cluster_x_start = {}
    cluster_x_end = {}
    current_x = left_margin

    # If full dendrogram is available, use its leaf order instead of cluster order
    if full_dendro_data is not None and 'read_order' in full_dendro_data:
        # Build reverse mapping: read -> cluster_id
        read_to_cluster = {}
        for cluster_id, data in cluster_reads.items():
            for read, sample in data['reads']:
                read_to_cluster[read] = cluster_id

        # Use dendrogram leaf order with optional cluster gap
        ordered_reads = full_dendro_data['read_order']
        dendro_cluster_gap = getattr(args, 'dendro_cluster_gap', 0)
        prev_cluster = None
        for read in ordered_reads:
            curr_cluster = read_to_cluster.get(read)
            # Add extra gap when crossing cluster boundary
            if dendro_cluster_gap > 0 and prev_cluster is not None and curr_cluster != prev_cluster:
                current_x += dendro_cluster_gap
            read_x_positions[read] = current_x
            current_x += group_width + args.read_spacing
            prev_cluster = curr_cluster
        # Set cluster boundaries to span all reads (single group)
        cluster_x_start[1] = left_margin
        cluster_x_end[1] = current_x - args.read_spacing
    else:
        # Original cluster-based ordering
        for cluster_id, data in cluster_reads.items():
            cluster_x_start[cluster_id] = current_x
            for read, sample in data['reads']:
                read_x_positions[read] = current_x
                current_x += group_width + args.read_spacing
            cluster_x_end[cluster_id] = current_x - args.read_spacing
            current_x += args.cluster_spacing

    # --- Calculate scaffold lengths ---
    scaffold_lengths = {}
    scaffold_min_starts = {}

    for read in read_data:
        for fs in read_data[read]:
            for feat in read_data[read][fs]:
                scaffold_lengths[read] = max(scaffold_lengths.get(read, 0), feat['stop'])
                if read not in scaffold_min_starts:
                    scaffold_min_starts[read] = feat['start']
                else:
                    scaffold_min_starts[read] = min(scaffold_min_starts[read], feat['start'])

    # --- Calculate drawing data ---
    # Calculate ratio - either from args or auto-calculated from target dimensions
    ratio = args.ratio

    # For horizontal mode, target_height controls bar length (vertical dimension)
    if args.target_height is not None:
        # Find max read length for scaling
        max_read_len = 0
        for read in scaffold_lengths:
            read_len = scaffold_lengths[read] - scaffold_min_starts.get(read, 0)
            max_read_len = max(max_read_len, read_len)

        if max_read_len > 0:
            # Estimate margins: top_margin (~200) + 50 + legend_margin (~100) ≈ 350
            # (legend is smaller in full_dendro mode)
            estimated_margins = 350
            target_bar_height = args.target_height - estimated_margins
            if target_bar_height > 50:  # Minimum bar height
                ratio = target_bar_height / max_read_len
                print(f"  Auto-calculated ratio from target_height: {ratio:.6f}")

    drawing_data = defaultdict(lambda: defaultdict(list))
    density_line_data = defaultdict(lambda: defaultdict(list))  # read -> featureset -> [(y, x_offset)]
    rect_plot_data = defaultdict(lambda: defaultdict(list))  # read -> featureset -> [{y, height, color}]
    uncolored_features = defaultdict(set)
    displayed_features = defaultdict(set)  # Track actually displayed features for legend filtering

    # Get colors for density line plot featuresets
    density_line_colors = {
        fs: get_primary_color(fs, featureset_colors, featureset_color_order)
        for fs in density_line_plot_featuresets
    }
    if density_line_colors:
        print(f"Density line colors: {density_line_colors}")

    # Get colors for rect plot featuresets - use explicit colors for known features
    rect_plot_explicit_colors = {
        'fiberseq_FIRE': '#FF4500',    # Orange-red for FIRE
        'fiberseq_LINKER': '#00CC66',  # Green for Linker
    }
    rect_plot_colors = {
        fs: rect_plot_explicit_colors.get(fs,
            get_primary_color(fs, featureset_colors, featureset_color_order))
        for fs in rect_plot_featuresets
    }
    if rect_plot_colors:
        print(f"Rect plot colors: {rect_plot_colors}")

    for read in read_data:
        if read not in read_x_positions:
            continue

        base_x = read_x_positions[read]
        scaffold_min_start = scaffold_min_starts.get(read, 0)
        read_length = scaffold_lengths.get(read, 0)

        for fs_idx, fs in enumerate(featuresets):
            x_offset = fs_idx * (args.bar_width + args.bar_spacing)

            # Handle density line plot track
            if fs == density_line_plot_name and density_line_plot_featuresets:
                # Compute density lines for each featureset in the combined track
                for line_fs in density_line_plot_featuresets:
                    features = read_data[read].get(line_fs, [])
                    if not features:
                        continue

                    # Compute density percentages
                    density_points = compute_density_line(
                        features, args.density_bin_size, read_length, max_density=50
                    )

                    # Convert to drawing coordinates
                    line_points = []
                    for bp_pos, density_pct in density_points:
                        # Y position based on bp position
                        y_pos = top_margin + 50 + floor((bp_pos - scaffold_min_start) * ratio)
                        # X offset based on density (0% = left edge, 100% = right edge of bar)
                        x_pos = base_x + x_offset + (density_pct / 100.0) * args.bar_width
                        line_points.append((x_pos, y_pos))

                    if line_points:
                        density_line_data[read][line_fs] = {
                            'points': line_points,
                            'color': density_line_colors.get(line_fs, "#FFFFFF"),
                            'base_x': base_x + x_offset
                        }

                # Add a placeholder rectangle for read_heights calculation
                if read_length > 0:
                    start_y = top_margin + 50
                    stop_y = top_margin + 50 + floor((read_length - scaffold_min_start) * ratio)
                    drawing_data[read][fs].append({
                        "x": base_x + x_offset,
                        "y": start_y,
                        "height": stop_y - start_y,
                        "fill": "none",
                        "fill_opacity": 0
                    })
                continue

            # Handle rect plot track (exact feature rectangles)
            if fs == rect_plot_name and rect_plot_featuresets:
                # Collect rectangles for each featureset in the combined track
                for rect_fs in rect_plot_featuresets:
                    features = read_data[read].get(rect_fs, [])
                    if not features:
                        continue

                    for feat in features:
                        if feat.get('feature') == 'unannotated':
                            continue
                        final_start = feat['start'] - scaffold_min_start
                        final_stop = feat['stop'] - scaffold_min_start
                        start_y = top_margin + 50 + floor(final_start * ratio)
                        stop_y = top_margin + 50 + floor(final_stop * ratio)
                        rect_plot_data[read][rect_fs].append({
                            'y': start_y,
                            'height': max(stop_y - start_y, 2),  # Minimum 2px height
                            'color': rect_plot_colors.get(rect_fs, "#FFFFFF"),
                            'base_x': base_x + x_offset
                        })

                # Add a placeholder rectangle for read_heights calculation
                if read_length > 0:
                    start_y = top_margin + 50
                    stop_y = top_margin + 50 + floor((read_length - scaffold_min_start) * ratio)
                    drawing_data[read][fs].append({
                        "x": base_x + x_offset,
                        "y": start_y,
                        "height": stop_y - start_y,
                        "fill": "none",
                        "fill_opacity": 0
                    })
                continue

            # Get features for this read/featureset
            features = read_data[read].get(fs, [])

            # Apply density computation if this featureset is marked for density
            if fs in density_featuresets and features:
                # Filter out unannotated features - we only count actual marks for density
                density_features = [f for f in features if f['feature'] != 'unannotated']
                # Extract base feature name from featureset (e.g., "fiberseq_m6A" -> "m6A")
                base_feature = fs.split('_')[-1] if '_' in fs else fs
                features = compute_density_features(
                    density_features, args.density_bin_size, read_length, base_feature
                )

            for feat in features:
                final_start = feat['start'] - scaffold_min_start
                final_stop = feat['stop'] - scaffold_min_start

                # Feature bars start after annotation bars (annot at +22, height ~25, so start at +50)
                start_y = top_margin + 50 + floor(final_start * ratio)
                stop_y = top_margin + 50 + floor(final_stop * ratio)

                feature_name = feat['feature']
                displayed_features[fs].add(feature_name)  # Track displayed feature
                color, fill_opacity = featureset_colors[fs].get(feature_name, ("#ffffff", 1.0))
                if feature_name not in featureset_colors[fs]:
                    uncolored_features[fs].add(feature_name)

                drawing_data[read][fs].append({
                    "x": base_x + x_offset,
                    "y": start_y,
                    "height": stop_y - start_y,
                    "fill": color,
                    "fill_opacity": fill_opacity
                })

    # --- Calculate image dimensions and read heights for borders ---
    max_stop_y = 0
    read_heights = {}  # read -> (min_y, max_y, x_start, total_width)

    for read in drawing_data:
        min_y = float('inf')
        max_y = 0
        min_x = float('inf')
        max_x = 0
        for fs in drawing_data[read]:
            for rect in drawing_data[read][fs]:
                if rect["height"] > 0:
                    min_y = min(min_y, rect["y"])
                    max_y = max(max_y, rect["y"] + rect["height"])
                    min_x = min(min_x, rect["x"])
                    max_x = max(max_x, rect["x"] + args.bar_width)
        if max_y > 0:
            read_heights[read] = (min_y, max_y, min_x, max_x - min_x)
            max_stop_y = max(max_stop_y, max_y)

    # Calculate required width for color legends (must match draw_color_legends)
    colors_per_column = 12
    item_width = 120
    legend_spacing = 15

    def get_legend_width(fs):
        num_items = len(featureset_color_order[fs])
        num_cols = (num_items + colors_per_column - 1) // colors_per_column
        return max(num_cols * item_width, 100)

    total_legend_width = left_margin + sum(get_legend_width(fs) + legend_spacing for fs in featuresets)

    # Image width is max of data width and legend width
    image_width = max(current_x + 50, total_legend_width + 50)
    # Reduce legend margin in full dendrogram mode (fewer legends to show)
    legend_bottom_margin = 100 if full_dendro_data is not None else 350
    image_height = max_stop_y + 50 + legend_bottom_margin

    print(f"\nImage dimensions: {image_width} x {image_height}")

    # --- Create drawing ---
    d = draw.Drawing(image_width, image_height)
    d.append(draw.Rectangle(0, 0, image_width, image_height, fill=background_color))

    # --- Draw components ---
    # Dendrogram header - either full (all reads) or cluster-level
    if full_dendro_data is not None:
        # Draw full dendrogram showing all individual reads
        draw_full_dendrogram_header(d, full_dendro_data, read_x_positions, group_width,
                                    top_margin, dendrogram_height, background_color)
    elif cluster_dendro_data is not None:
        draw_cluster_dendrogram(d, cluster_dendro_data, cluster_x_start, cluster_x_end,
                                top_margin, dendrogram_height, background_color)

    # Cluster brackets - positioned below feature labels per cluster
    # Calculate label height: longest featureset name × font_size (~4.5px per char for font_size=6)
    longest_label = max((len(fs_display_names.get(fs, fs)) for fs in featuresets), default=10)
    label_height = longest_label * 4.5  # font_size=6, tight fit
    # Skip cluster brackets when using full dendrogram (individual read view)
    if not args.hide_brackets and full_dendro_data is None:
        draw_cluster_brackets(d, cluster_reads, cluster_x_start, cluster_x_end,
                             enrichment_colors, read_heights, label_height, text_color,
                             cluster_labels=cluster_labels)

    # Annotation bars
    draw_annotation_bars(d, cluster_reads, read_x_positions, read_to_original_cluster,
                        read_to_original_enrichment, sample_colors, cluster_colors,
                        enrichment_colors, group_width, top_margin, left_margin, text_color)

    # Feature bars
    draw_feature_bars(d, drawing_data, featuresets, args.bar_width, read_heights, num_featuresets,
                      density_line_data=density_line_data,
                      rect_plot_data=rect_plot_data if rect_plot_featuresets else None,
                      background_color=background_color)

    # Read labels
    draw_read_labels(d, cluster_reads, read_x_positions, group_width, top_margin, text_color)

    # Featureset labels below each read's bars
    for read, (min_y, max_y, x_start, total_width) in read_heights.items():
        if read not in read_x_positions:
            continue

        base_x = read_x_positions[read]
        label_base_y = max_y + 5  # Just below this read's bars

        for fs_idx, fs in enumerate(featuresets):
            x_offset = fs_idx * (args.bar_width + args.bar_spacing)
            display_name = fs_display_names.get(fs, fs)
            label_x = base_x + x_offset + args.bar_width / 2

            d.append(draw.Text(
                display_name, font_size=6, x=label_x, y=label_base_y,
                fill=text_color, font_family='sans-serif',
                text_anchor='start', dominant_baseline='middle',
                transform=f"rotate(90 {label_x} {label_base_y})"
            ))

    # --- Draw legends ---
    legend_y = 20
    legend_x = left_margin

    # Sample, Cluster, Enrichment legends (stacked vertically)
    draw_top_legends(d, sample_colors, cluster_colors, read_to_original_cluster,
                     read_to_original_enrichment, enrichment_colors,
                     legend_x, legend_y, text_color)

    # Density line / Rect plot legend (row 4, if applicable)
    if density_line_colors or rect_plot_colors:
        density_legend_y = legend_y + 60  # After sample, cluster, enrichment rows
        draw_density_line_legend(d, density_line_colors,
                                 rect_plot_colors if rect_plot_featuresets else None,
                                 legend_x, density_legend_y, text_color)

    # Color legends at bottom (exclude density_line since it has its own legend)
    color_legend_y_start = max_stop_y + 130
    bottom_legend_featuresets = [fs for fs in featuresets if fs != density_line_plot_name]
    draw_color_legends(d, bottom_legend_featuresets, featureset_colors, featureset_color_order,
                      fs_display_names, color_legend_y_start, left_margin, text_color,
                      displayed_features=displayed_features)

    # --- Save ---
    d.save_svg(args.output)

    # --- Report warnings ---
    for fs in featuresets:
        if uncolored_features[fs]:
            sys.stderr.write(f"Warning: {fs} - features not in colors file:\n")
            for feature in sorted(list(uncolored_features[fs])):
                sys.stderr.write(f"  - {feature}\n")

    print(f"\n--- Summary ---")
    # Use original cluster count (before dendrogram reordering merges them)
    original_cluster_count = len(set(read_to_original_cluster.values())) if read_to_original_cluster else len(cluster_reads)
    print(f"Clusters plotted: {original_cluster_count}")
    total_reads = sum(len(data['reads']) for data in cluster_reads.values())
    print(f"Total reads plotted: {total_reads}")
    print(f"\n✅ Saved to {args.output}")

    # --- Print parameters table ---
    print("\n" + "=" * 60)
    print("Parameters")
    print("=" * 60)
    params = [
        ("cluster-analysis-prefix", args.cluster_prefix),
        ("output", args.output),
        ("bed", f"{len(args.bed_files)} file(s)" if args.bed_files else "auto-discovered"),
        ("input-bed-prefix", args.input_bed_prefix if args.input_bed_prefix else "N/A"),
        ("database", database),
        ("colors", args.colors_dir),
        ("featuresets", args.featuresets),
        ("background", args.background_color),
        ("bar-width", args.bar_width),
        ("bar-spacing", args.bar_spacing),
        ("read-spacing", args.read_spacing),
        ("cluster-spacing", args.cluster_spacing),
        ("ratio", args.ratio),
        ("smoothness", args.smoothness),
        ("reads-file", args.reads_file if args.reads_file else "None"),
        ("hide-brackets", args.hide_brackets),
        ("hide-dendrogram", args.hide_dendrogram),
        ("no-reorder", args.no_reorder),
        ("cluster-labels", args.cluster_labels if args.cluster_labels else "None"),
        ("vertical", args.vertical),
        ("n-per-cluster", args.max_reps if args.max_reps else "None"),
        ("log-file", args.log_file),
    ]
    print(f"{'Parameter':<25} {'Value':<35}")
    print(f"{'-' * 25} {'-' * 35}")
    for param, value in params:
        print(f"{param:<25} {str(value):<35}")

    # --- Print command ---
    print("\n" + "=" * 60)
    print("Command")
    print("=" * 60)
    print(_original_command)


if __name__ == "__main__":
    main()
