# KaryoScope Cluster Analysis
# Analyzes hierarchical clustering results to identify biologically interesting clusters
# with read assignments sorted by centroid distance for visualization
#
# Usage:
# python KaryoScope_cluster_analysis.py \
#   --bed /path/to/sample1.bed.gz /path/to/sample2.bed.gz \
#   --sample-metadata samples.tsv \
#   --output-prefix analysis_output \
#   --n-clusters 10
#
# Sample metadata file format (TSV):
#   sample      group   color
#   SW26_Pre    pre     #377EB8
#   SW26_Post   post    #E41A1C
#
# Comparison modes:
#   two-group: Fisher's exact test between control and treatment groups
#   multi-group: Chi-square test across all groups
#   per-sample: Each sample vs all others (no groups required)

import argparse
import gzip
import numpy as np
import pandas as pd
from collections import defaultdict, Counter
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list, fcluster, cut_tree
from scipy.spatial.distance import pdist, squareform
from scipy.stats import fisher_exact, chi2_contingency
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
matplotlib.use('Agg')

# --- Argument parsing ---
parser = argparse.ArgumentParser(
    description="Analyze KaryoScope clustering to identify sample-enriched clusters.",
    formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument("--bed", required=True, nargs='+',
                    help="Path to input BED file(s) (can be gzipped). Multiple files will be concatenated.")
parser.add_argument("--bed2", dest="bed2", nargs='+', default=None,
                    help="Optional second BED file(s) to merge with --bed by position overlay.\n"
                         "Must have same number of files as --bed, in corresponding order.\n"
                         "Creates combined feature labels (e.g., 'chr1:p_arm' from chromosome and region).")
parser.add_argument("--output-prefix", dest="output_prefix", required=True,
                    help="Prefix for output files")
parser.add_argument("--sample-metadata", dest="sample_metadata", default=None,
                    help="TSV file with sample metadata (columns: sample, group, color).\n"
                         "If not provided, each sample becomes its own group.\n"
                         "Use --control-group to specify the reference group.")
parser.add_argument("--comparison-mode", dest="comparison_mode", default="two-group",
                    choices=["two-group", "multi-group", "per-sample"],
                    help="Comparison mode for enrichment testing:\n"
                         "  two-group: Fisher's exact test between control and treatment\n"
                         "  multi-group: Chi-square test across all groups\n"
                         "  per-sample: Each sample vs all others (default: two-group)")
parser.add_argument("--control-group", dest="control_group", default=None,
                    help="Name of control group for two-group comparison (default: auto-detect)")
parser.add_argument("--n-clusters", dest="n_clusters", type=int, default=None,
                    help="Number of clusters to cut tree into (default: auto-determine)")
parser.add_argument("--min-k", dest="min_k", type=int, default=20,
                    help="Minimum number of clusters to test during auto-detection (default: 20)")
parser.add_argument("--max-k", dest="max_k", type=int, default=100,
                    help="Maximum number of clusters to test during auto-detection (default: 100)")
parser.add_argument("--min-cluster-size", dest="min_cluster_size", type=int, default=10,
                    help="Minimum cluster size to consider (default: 10)")
parser.add_argument("--min-read-length", dest="min_read_length", type=int, default=10000,
                    help="Minimum read length in bp to include (default: 10000)")
parser.add_argument("--small-feature-quantile", dest="small_quantile", type=float, default=0.1,
                    help="Quantile threshold for small feature collapsing (default: 0.1)")
parser.add_argument("--linkage-method", dest="linkage_method", default="ward",
                    help="Linkage method for hierarchical clustering (default: ward)")
parser.add_argument("--matrix-type", dest="matrix_type", default="length_weighted",
                    choices=["binary", "count", "length_weighted"],
                    help="Type of adjacency matrix:\n"
                         "  binary: 0/1 for presence/absence of transitions\n"
                         "  count: count of each transition\n"
                         "  length_weighted: transitions weighted by feature length (default: length_weighted)")
parser.add_argument("--edges", dest="edge_mode", default="bidirectional",
                    choices=["directional", "bidirectional", "symmetric"],
                    help="Edge counting mode:\n"
                         "  directional: standard A->B edge counting\n"
                         "  bidirectional: A->B and B->A are both counted separately\n"
                         "  symmetric: edges are sorted alphabetically, A->B and B->A both count as A->B (default: bidirectional)")
parser.add_argument("--abundance", dest="include_abundance",
                    action=argparse.BooleanOptionalAction, default=True,
                    help="Include feature abundance dimensions (default: True)")
parser.add_argument("--umap", dest="plot_umap",
                    action=argparse.BooleanOptionalAction, default=True,
                    help="Generate UMAP visualization (default: True, requires umap-learn)")
parser.add_argument("--circular-dendrogram", dest="plot_circular_dendrogram",
                    action=argparse.BooleanOptionalAction, default=True,
                    help="Generate circular dendrogram visualization (default: True)")
parser.add_argument("--umap-neighbors", dest="umap_neighbors", type=int, default=25,
                    help="UMAP n_neighbors parameter (default: 25)")
parser.add_argument("--umap-min-dist", dest="umap_min_dist", type=float, default=0.2,
                    help="UMAP min_dist parameter (default: 0.2)")
parser.add_argument("--perfect-threshold", dest="perfect_threshold", type=float, default=0.95,
                    help="Threshold for perfect enrichment (default: 0.95 = 95%%)")
parser.add_argument("--strong-threshold", dest="strong_threshold", type=float, default=0.80,
                    help="Threshold for strong enrichment (default: 0.80 = 80%%)")
parser.add_argument("--davies-bouldin", dest="compute_davies_bouldin",
                    action=argparse.BooleanOptionalAction, default=False,
                    help="Compute Davies-Bouldin index (not recommended, favors high k) (default: False)")
parser.add_argument("--early-stopping", dest="early_stopping", type=int, default=25,
                    help="Stop k search if no improvement for N iterations (0 to disable) (default: 25)")
parser.add_argument("--nested", dest="nested",
                    action=argparse.BooleanOptionalAction, default=False,
                    help="Hierarchical testing: first test groups, then test samples within enriched groups (default: False)")
parser.add_argument("--also-test-samples", dest="also_test_samples",
                    action=argparse.BooleanOptionalAction, default=False,
                    help="Run per-sample Fisher's tests in addition to group-level tests (default: False)")
parser.add_argument("--stratified", dest="stratified",
                    action=argparse.BooleanOptionalAction, default=False,
                    help="Report within-group sample breakdown and variance metrics (default: False)")

args = parser.parse_args()

# --- Helper functions ---
def load_bed_file(filepath, sample_label=None):
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
                records.append({
                    'read': read,
                    'start': start,
                    'end': end,
                    'feature': feature,
                    'sample': sample_label
                })
    return pd.DataFrame(records)

def extract_sample_label(filepath):
    """Extract sample label from filepath."""
    import os
    basename = os.path.basename(filepath)
    parts = basename.split('.')
    if parts:
        return parts[0]
    return basename


def merge_beds_by_position(df1, df2, sample_label=None, sep=":"):
    """Merge two BED DataFrames by position overlay, creating combined feature labels.

    Uses pyranges for fast interval intersection. Falls back to numpy if unavailable.

    Args:
        df1: DataFrame with columns [read, start, end, feature] from first BED file
        df2: DataFrame with columns [read, start, end, feature] from second BED file
        sample_label: Optional sample label to add to the merged records
        sep: Separator for combined feature labels (default: ":")

    Returns:
        DataFrame with merged intervals and combined feature labels

    Example:
        df1 (chromosome):  read1  0-100  chr1_specific
        df2 (region):      read1  0-50   p_arm_specific
                           read1  50-100 arm_multigroup1
        Result:            read1  0-50   chr1_specific:p_arm_specific
                           read1  50-100 chr1_specific:arm_multigroup1
    """
    try:
        import pyranges as pr
        return _merge_beds_pyranges(df1, df2, sample_label, sep)
    except ImportError:
        print("  Warning: pyranges not available, using numpy method")
        return _merge_beds_pandas(df1, df2, sample_label, sep)


def _merge_beds_pyranges(df1, df2, sample_label=None, sep=":"):
    """Fast merge using pyranges join/intersect."""
    import pyranges as pr

    # Filter to common reads first
    common_reads = set(df1['read'].unique()) & set(df2['read'].unique())
    if not common_reads:
        print(f"  Warning: No common reads between BED files for sample {sample_label}")
        return pd.DataFrame(columns=['read', 'start', 'end', 'feature', 'sample'])

    df1_f = df1[df1['read'].isin(common_reads)][['read', 'start', 'end', 'feature']].copy()
    df2_f = df2[df2['read'].isin(common_reads)][['read', 'start', 'end', 'feature']].copy()

    # Rename columns for pyranges (Chromosome, Start, End)
    df1_f.columns = ['Chromosome', 'Start', 'End', 'Feature1']
    df2_f.columns = ['Chromosome', 'Start', 'End', 'Feature2']

    # Create PyRanges objects
    pr1 = pr.PyRanges(df1_f)
    pr2 = pr.PyRanges(df2_f)

    # Join/intersect the two ranges
    # This returns overlapping intervals with features from both
    joined = pr1.join(pr2)

    if len(joined) == 0:
        return pd.DataFrame(columns=['read', 'start', 'end', 'feature', 'sample'])

    # Convert back to DataFrame
    result_df = joined.df

    # The result has columns: Chromosome, Start, End, Feature1, Start_b, End_b, Feature2
    # Calculate actual overlap region
    result_df['overlap_start'] = result_df[['Start', 'Start_b']].max(axis=1)
    result_df['overlap_end'] = result_df[['End', 'End_b']].min(axis=1)

    # Keep only actual overlaps
    result_df = result_df[result_df['overlap_end'] > result_df['overlap_start']]

    # Create combined feature labels
    result_df['feature'] = result_df['Feature1'] + sep + result_df['Feature2']

    # Build final result
    result = pd.DataFrame({
        'read': result_df['Chromosome'],
        'start': result_df['overlap_start'].astype(int),
        'end': result_df['overlap_end'].astype(int),
        'feature': result_df['feature'],
        'sample': sample_label
    })

    return result


def _merge_beds_pandas(df1, df2, sample_label=None, sep=":"):
    """Efficient merge using numpy vectorized operations per-read."""
    # Get common reads
    common_reads = set(df1['read'].unique()) & set(df2['read'].unique())

    if not common_reads:
        print(f"  Warning: No common reads between BED files for sample {sample_label}")
        return pd.DataFrame(columns=['read', 'start', 'end', 'feature', 'sample'])

    # Filter and index by read
    df1_f = df1[df1['read'].isin(common_reads)].copy()
    df2_f = df2[df2['read'].isin(common_reads)].copy()

    # Group by read for efficient per-read processing
    grouped1 = {read: grp[['start', 'end', 'feature']].values for read, grp in df1_f.groupby('read')}
    grouped2 = {read: grp[['start', 'end', 'feature']].values for read, grp in df2_f.groupby('read')}

    merged_records = []

    for read in common_reads:
        intervals1 = grouped1.get(read)
        intervals2 = grouped2.get(read)

        if intervals1 is None or intervals2 is None or len(intervals1) == 0 or len(intervals2) == 0:
            continue

        # Extract starts, ends, features as numpy arrays
        starts1 = intervals1[:, 0].astype(int)
        ends1 = intervals1[:, 1].astype(int)
        feats1 = intervals1[:, 2]

        starts2 = intervals2[:, 0].astype(int)
        ends2 = intervals2[:, 1].astype(int)
        feats2 = intervals2[:, 2]

        # Collect all breakpoints
        breakpoints = np.unique(np.concatenate([starts1, ends1, starts2, ends2]))

        # For each interval between breakpoints, find covering features
        for i in range(len(breakpoints) - 1):
            bp_start = breakpoints[i]
            bp_end = breakpoints[i + 1]

            if bp_end <= bp_start:
                continue

            midpoint = (bp_start + bp_end) / 2

            # Find feature from df1 covering midpoint (vectorized)
            mask1 = (starts1 <= midpoint) & (midpoint < ends1)
            if not mask1.any():
                continue
            feat1 = feats1[mask1][0]

            # Find feature from df2 covering midpoint (vectorized)
            mask2 = (starts2 <= midpoint) & (midpoint < ends2)
            if not mask2.any():
                continue
            feat2 = feats2[mask2][0]

            merged_records.append({
                'read': read,
                'start': int(bp_start),
                'end': int(bp_end),
                'feature': f"{feat1}{sep}{feat2}",
                'sample': sample_label
            })

    return pd.DataFrame(merged_records)

def collapse_small_clusters(df, small_len_thresh, small_gap_thresh):
    """Collapse small adjacent features within each read.

    Original version - slower but proven to work correctly.
    Only collapses consecutive small features that are close together.
    """
    df = df.sort_values(['read', 'start']).reset_index(drop=True)

    def collapse_read_features(group):
        """Collapse small adjacent features within a single read."""
        group = group.sort_values('start').reset_index(drop=True)

        if len(group) <= 1:
            return group

        collapsed_rows = []
        i = 0

        while i < len(group):
            row = group.iloc[i]

            # If this feature is not small, keep it as-is
            if row['length'] > small_len_thresh:
                collapsed_rows.append({
                    'read': row['read'],
                    'start': row['start'],
                    'end': row['end'],
                    'feature': row['feature'],
                    'length': row['length']
                })
                i += 1
                continue

            # This feature is small - try to merge with following small features
            merge_start = row['start']
            merge_end = row['end']
            features_in_merge = [row['feature']]

            j = i + 1
            while j < len(group):
                next_row = group.iloc[j]
                gap = next_row['start'] - merge_end

                # Stop merging if gap too big or next feature not small
                if gap > small_gap_thresh or next_row['length'] > small_len_thresh:
                    break

                # Merge this feature
                merge_end = next_row['end']
                features_in_merge.append(next_row['feature'])
                j += 1

            # Use mode (most common) feature for the merged region
            feature_counts = Counter(features_in_merge)
            mode_feature = feature_counts.most_common(1)[0][0]

            collapsed_rows.append({
                'read': row['read'],
                'start': merge_start,
                'end': merge_end,
                'feature': mode_feature,
                'length': merge_end - merge_start
            })

            i = j

        return pd.DataFrame(collapsed_rows)

    # Apply to each read
    result = df.groupby('read', group_keys=False).apply(collapse_read_features)
    return result.reset_index(drop=True)

def get_edges(features, edge_mode="directional"):
    """Get edges from a list of features based on edge counting mode.

    Args:
        features: List of feature names in order
        edge_mode: One of:
            - "directional": standard A->B edge counting
            - "bidirectional": A->B and B->A are both counted separately
            - "symmetric": edges sorted alphabetically, A->B and B->A both count as A->B

    Returns:
        List of (from, to) tuples
    """
    if len(features) <= 1:
        return []
    edges = []
    for i in range(len(features) - 1):
        from_feat, to_feat = features[i], features[i + 1]

        if edge_mode == "directional":
            edges.append((from_feat, to_feat))
        elif edge_mode == "bidirectional":
            edges.append((from_feat, to_feat))
            edges.append((to_feat, from_feat))
        elif edge_mode == "symmetric":
            # Sort alphabetically so A->B and B->A both become the same edge
            sorted_pair = tuple(sorted([from_feat, to_feat]))
            edges.append(sorted_pair)
    return edges

def load_sample_metadata(metadata_file, sample_labels):
    """Load sample metadata from TSV file or auto-generate defaults.

    When a metadata file is provided:
      - Uses group names exactly as specified in the file
      - Colors are optional

    When no metadata file:
      - Each sample becomes its own group (named after the sample)
      - This allows --control-group to specify the reference
    """
    if metadata_file:
        # Load from file
        meta_df = pd.read_csv(metadata_file, sep='\t')

        # Validate required columns
        if 'sample' not in meta_df.columns:
            raise ValueError("Sample metadata file must have 'sample' column")

        # Build mappings
        sample_to_group = {}
        sample_to_color = {}

        for _, row in meta_df.iterrows():
            sample = row['sample']
            # Use group if provided, otherwise use sample name as group
            sample_to_group[sample] = row.get('group', sample) if 'group' in meta_df.columns else sample
            if 'color' in meta_df.columns and pd.notna(row.get('color')):
                sample_to_color[sample] = row['color']

        # Check all samples are covered
        missing = set(sample_labels) - set(sample_to_group.keys())
        if missing:
            print(f"  Warning: Samples not in metadata file: {missing}")
            # Auto-assign missing samples to their own groups
            for s in missing:
                sample_to_group[s] = s

        return sample_to_group, sample_to_color

    else:
        # No metadata file: each sample is its own group
        # This is clean and doesn't rely on text matching
        sample_to_group = {sample: sample for sample in sample_labels}
        return sample_to_group, {}


def generate_group_colors(groups, existing_colors=None):
    """Generate colors for groups, using existing colors where provided.

    Color assignment logic:
    1. Use explicitly provided colors from metadata file
    2. For exactly 2 groups: use blue (#377EB8) for first, red (#E41A1C) for second
    3. For more groups: use tab10 colormap
    """
    if existing_colors is None:
        existing_colors = {}

    sorted_groups = sorted(set(groups))
    n_groups = len(sorted_groups)

    # Color palette for multiple groups
    tab10 = matplotlib.colormaps.get_cmap('tab10')

    group_colors = {}

    for i, group in enumerate(sorted_groups):
        if group in existing_colors:
            # Use explicitly provided color
            group_colors[group] = existing_colors[group]
        elif n_groups == 2:
            # Two groups: blue for first (alphabetically), red for second
            group_colors[group] = '#377EB8' if i == 0 else '#E41A1C'
        else:
            # Multiple groups: use colormap
            group_colors[group] = mcolors.rgb2hex(tab10(i % 10))

    return group_colors


def calculate_enrichment_two_group(cluster_samples, sample_to_group, control_group, group_totals):
    """Calculate enrichment using Fisher's exact test (two-group comparison)."""
    # Get treatment group (the other group)
    groups = list(group_totals.keys())
    treatment_group = [g for g in groups if g != control_group][0] if len(groups) == 2 else None

    if treatment_group is None:
        # If more than 2 groups, pick the most enriched non-control group
        group_counts = Counter(sample_to_group.get(s, s) for s in cluster_samples)
        non_control = {g: c for g, c in group_counts.items() if g != control_group}
        treatment_group = max(non_control.keys(), key=lambda g: non_control[g]) if non_control else control_group

    # Count samples in cluster by group
    control_in = sum(1 for s in cluster_samples if sample_to_group.get(s, s) == control_group)
    treatment_in = sum(1 for s in cluster_samples if sample_to_group.get(s, s) == treatment_group)
    control_out = group_totals.get(control_group, 0) - control_in
    treatment_out = group_totals.get(treatment_group, 0) - treatment_in

    # Fisher's exact test
    contingency = [[control_in, treatment_in], [control_out, treatment_out]]
    odds_ratio, p_value = fisher_exact(contingency)

    # Calculate percentages
    total_in = control_in + treatment_in
    control_pct = (control_in / total_in * 100) if total_in > 0 else 0
    treatment_pct = (treatment_in / total_in * 100) if total_in > 0 else 0

    # Determine enrichment direction
    if p_value < 0.05:
        if control_pct > treatment_pct:
            enrichment = f"{control_group}-enriched"
        elif treatment_pct > control_pct:
            enrichment = f"{treatment_group}-enriched"
        else:
            enrichment = "mixed"
    else:
        enrichment = "mixed"

    # Build group counts dict
    group_counts = {control_group: control_in, treatment_group: treatment_in}
    group_pcts = {control_group: control_pct, treatment_group: treatment_pct}

    return {
        'group_counts': group_counts,
        'group_pcts': group_pcts,
        'odds_ratio': odds_ratio,
        'p_value': p_value,
        'enrichment': enrichment,
        'dominant_group': control_group if control_pct > treatment_pct else treatment_group
    }


def calculate_enrichment_multi_group(cluster_samples, sample_to_group, group_totals):
    """Calculate enrichment using chi-square test (multi-group comparison)."""
    groups = list(group_totals.keys())

    # Count samples in cluster by group
    group_counts = Counter(sample_to_group.get(s, s) for s in cluster_samples)

    # Build contingency table: [in_cluster, out_cluster] for each group
    observed_in = [group_counts.get(g, 0) for g in groups]
    observed_out = [group_totals[g] - group_counts.get(g, 0) for g in groups]

    contingency = [observed_in, observed_out]

    # Chi-square test (handle edge cases)
    try:
        if sum(observed_in) > 0 and len(groups) > 1:
            chi2, p_value, dof, expected = chi2_contingency(contingency)
        else:
            p_value = 1.0
    except ValueError:
        p_value = 1.0

    # Calculate percentages
    total_in = sum(observed_in)
    group_pcts = {g: (group_counts.get(g, 0) / total_in * 100) if total_in > 0 else 0 for g in groups}

    # Determine dominant group
    dominant_group = max(groups, key=lambda g: group_pcts[g])

    # Enrichment label
    if p_value < 0.05:
        enrichment = f"{dominant_group}-enriched"
    else:
        enrichment = "mixed"

    return {
        'group_counts': {g: group_counts.get(g, 0) for g in groups},
        'group_pcts': group_pcts,
        'odds_ratio': np.nan,  # Not applicable for multi-group
        'p_value': p_value,
        'enrichment': enrichment,
        'dominant_group': dominant_group
    }


def calculate_enrichment_per_sample(cluster_samples, sample_totals):
    """Calculate enrichment per sample (each sample vs all others)."""
    samples = list(sample_totals.keys())

    # Count samples in cluster
    sample_counts = Counter(cluster_samples)

    # Calculate percentages
    total_in = len(cluster_samples)
    sample_pcts = {s: (sample_counts.get(s, 0) / total_in * 100) if total_in > 0 else 0 for s in samples}

    # Fisher's exact for each sample vs rest
    p_values = {}
    for sample in samples:
        in_cluster = sample_counts.get(sample, 0)
        out_cluster = sample_totals[sample] - in_cluster
        other_in = total_in - in_cluster
        other_out = sum(sample_totals.values()) - sample_totals[sample] - other_in

        _, p_val = fisher_exact([[in_cluster, other_in], [out_cluster, other_out]])
        p_values[sample] = p_val

    # Find most significant sample
    min_p_sample = min(p_values.keys(), key=lambda s: p_values[s])
    min_p = p_values[min_p_sample]

    # Enrichment label
    if min_p < 0.05:
        enrichment = f"{min_p_sample}-enriched"
    else:
        enrichment = "mixed"

    return {
        'group_counts': {s: sample_counts.get(s, 0) for s in samples},
        'group_pcts': sample_pcts,
        'odds_ratio': np.nan,
        'p_value': min_p,
        'enrichment': enrichment,
        'dominant_group': max(samples, key=lambda s: sample_pcts[s])
    }


def calculate_nested_within_group(cluster_samples, enriched_group, sample_to_group, sample_totals):
    """Calculate within-group sample enrichment for nested/hierarchical testing.

    Tests each sample within the enriched group to see if it drives the enrichment.
    """
    # Get samples belonging to the enriched group
    group_samples = [s for s in sample_totals.keys() if sample_to_group.get(s, s) == enriched_group]

    # Count samples in cluster
    sample_counts = Counter(cluster_samples)
    cluster_group_samples = [s for s in cluster_samples if sample_to_group.get(s, s) == enriched_group]
    total_in_group_cluster = len(cluster_group_samples)

    within_group_results = {}
    for sample in group_samples:
        in_cluster = sample_counts.get(sample, 0)
        out_cluster = sample_totals[sample] - in_cluster
        other_in = total_in_group_cluster - in_cluster
        # Other samples in this group, not in cluster
        other_group_total = sum(sample_totals[s] for s in group_samples if s != sample)
        other_out = other_group_total - other_in

        if other_group_total > 0:
            _, p_val = fisher_exact([[in_cluster, other_in], [out_cluster, other_out]])
        else:
            p_val = 1.0

        pct_of_group = (in_cluster / total_in_group_cluster * 100) if total_in_group_cluster > 0 else 0

        within_group_results[sample] = {
            'count': in_cluster,
            'pct_of_group': pct_of_group,
            'p_value': p_val,
            'enriched': p_val < 0.05 and pct_of_group > (100 / len(group_samples)),
            'depleted': p_val < 0.05 and pct_of_group < (100 / len(group_samples))
        }

    return within_group_results


def calculate_all_sample_tests(cluster_samples, sample_totals):
    """Calculate per-sample Fisher's tests for all samples (Option B: --also-test-samples)."""
    sample_counts = Counter(cluster_samples)
    total_in = len(cluster_samples)
    total_all = sum(sample_totals.values())

    sample_results = {}
    for sample, sample_total in sample_totals.items():
        in_cluster = sample_counts.get(sample, 0)
        out_cluster = sample_total - in_cluster
        other_in = total_in - in_cluster
        other_out = total_all - sample_total - other_in

        _, p_val = fisher_exact([[in_cluster, other_in], [out_cluster, other_out]])
        pct = (in_cluster / total_in * 100) if total_in > 0 else 0
        expected_pct = (sample_total / total_all * 100) if total_all > 0 else 0

        sample_results[sample] = {
            'count': in_cluster,
            'pct': pct,
            'expected_pct': expected_pct,
            'p_value': p_val,
            'enriched': p_val < 0.05 and pct > expected_pct,
            'depleted': p_val < 0.05 and pct < expected_pct
        }

    return sample_results


def calculate_stratified_variance(cluster_samples, sample_to_group, sample_totals, group_totals):
    """Calculate within-group variance metrics (Option C: --stratified)."""
    sample_counts = Counter(cluster_samples)

    # Group samples by their group
    groups = set(sample_to_group.values())

    stratified_results = {}
    for group in groups:
        group_sample_names = [s for s, g in sample_to_group.items() if g == group]
        if not group_sample_names:
            continue

        # Get counts for each sample in this group
        counts = [sample_counts.get(s, 0) for s in group_sample_names]

        # Calculate statistics
        mean_count = np.mean(counts)
        std_count = np.std(counts)
        cv = (std_count / mean_count) if mean_count > 0 else 0  # Coefficient of variation

        # Normalized counts (as % of each sample's total reads)
        normalized = []
        for s in group_sample_names:
            if sample_totals[s] > 0:
                normalized.append(sample_counts.get(s, 0) / sample_totals[s] * 100)
            else:
                normalized.append(0)

        mean_normalized = np.mean(normalized)
        std_normalized = np.std(normalized)
        cv_normalized = (std_normalized / mean_normalized) if mean_normalized > 0 else 0

        stratified_results[group] = {
            'samples': group_sample_names,
            'counts': dict(zip(group_sample_names, counts)),
            'mean_count': mean_count,
            'std_count': std_count,
            'cv': cv,
            'normalized_pcts': dict(zip(group_sample_names, normalized)),
            'mean_normalized': mean_normalized,
            'std_normalized': std_normalized,
            'cv_normalized': cv_normalized,
            'consistent': cv_normalized < 0.5,  # Flag if CV > 50%
            'n_samples': len(group_sample_names)
        }

    return stratified_results


# --- Load and process data ---
print("=" * 60)
print("KaryoScope Cluster Analysis")
print("=" * 60)

print(f"\nLoading BED file(s)...")
dfs = []
sample_labels = []

# Check if we're merging two featuresets
if args.bed2:
    if len(args.bed) != len(args.bed2):
        raise ValueError(f"--bed and --bed2 must have same number of files "
                        f"({len(args.bed)} vs {len(args.bed2)})")
    print(f"  Merging two featuresets by position...")

    for bed_file1, bed_file2 in zip(args.bed, args.bed2):
        sample_label = extract_sample_label(bed_file1)
        sample_labels.append(sample_label)
        print(f"  - Merging: {bed_file1}")
        print(f"           + {bed_file2}")
        print(f"             (sample: {sample_label})")

        # Load both BED files
        df1 = load_bed_file(bed_file1, sample_label)
        df2 = load_bed_file(bed_file2, sample_label)

        # Merge by position
        merged_df = merge_beds_by_position(df1, df2, sample_label)
        print(f"             → {len(merged_df):,} merged intervals")
        dfs.append(merged_df)
else:
    # Standard single featureset loading
    for bed_file in args.bed:
        sample_label = extract_sample_label(bed_file)
        sample_labels.append(sample_label)
        print(f"  - {bed_file} (sample: {sample_label})")
        df = load_bed_file(bed_file, sample_label)
        dfs.append(df)

in_data = pd.concat(dfs, ignore_index=True)
read_to_sample = in_data.groupby('read')['sample'].first().to_dict()

print(f"\nTotal records: {len(in_data):,}")
print(f"Unique reads: {in_data['read'].nunique():,}")

# Calculate read lengths
in_data['length'] = in_data['end'] - in_data['start']
read_lengths = in_data.groupby('read').agg({'start': 'min', 'end': 'max'})
read_lengths['read_length'] = read_lengths['end'] - read_lengths['start']
read_length_dict = read_lengths['read_length'].to_dict()

# --- Filter by minimum read length ---
if args.min_read_length > 0:
    print(f"\n--- Filtering reads by length ---")
    reads_before = in_data['read'].nunique()
    valid_reads = set(r for r, l in read_length_dict.items() if l >= args.min_read_length)
    in_data = in_data[in_data['read'].isin(valid_reads)]
    read_to_sample = {r: s for r, s in read_to_sample.items() if r in valid_reads}
    read_length_dict = {r: l for r, l in read_length_dict.items() if r in valid_reads}
    reads_after = in_data['read'].nunique()
    print(f"  Minimum read length: {args.min_read_length:,} bp")
    print(f"  Reads before filter: {reads_before:,}")
    print(f"  Reads after filter: {reads_after:,}")
    print(f"  Reads removed: {reads_before - reads_after:,}")

# --- Load sample metadata ---
print(f"\n--- Loading sample metadata ---")
sample_to_group, sample_colors_from_meta = load_sample_metadata(args.sample_metadata, sample_labels)

# Determine unique groups
all_groups = sorted(set(sample_to_group.values()))
print(f"  Groups found: {', '.join(all_groups)}")
print(f"  Comparison mode: {args.comparison_mode}")

# Map reads to groups
read_to_group = {r: sample_to_group.get(s, s) for r, s in read_to_sample.items()}

# Count samples by group
group_totals = Counter(read_to_group.values())
sample_totals = Counter(read_to_sample.values())

print(f"\n  Group counts:")
for group, count in sorted(group_totals.items()):
    print(f"    {group}: {count:,} reads")

# Determine control group for two-group comparison
control_group = args.control_group
if args.comparison_mode == "two-group":
    if control_group is None:
        # Default: alphabetically first group is the control/reference
        control_group = all_groups[0]
        if len(all_groups) > 2:
            print(f"  Note: {len(all_groups)} groups found, using '{control_group}' as reference.")
            print(f"        Use --control-group to specify a different reference group.")
    print(f"  Reference group: {control_group}")

# Convert sample colors to group colors (use first sample's color for each group)
group_colors_from_meta = {}
for sample, color in sample_colors_from_meta.items():
    group = sample_to_group.get(sample, sample)
    if group not in group_colors_from_meta:
        group_colors_from_meta[group] = color

# Generate colors for groups
group_colors = generate_group_colors(all_groups, group_colors_from_meta)

# Generate colors for samples (derive from group colors or use custom)
sample_colors = {}
for sample in sample_labels:
    if sample in sample_colors_from_meta:
        sample_colors[sample] = sample_colors_from_meta[sample]
    else:
        # Use group color
        group = sample_to_group.get(sample, sample)
        sample_colors[sample] = group_colors.get(group, '#999999')

# --- Collapse small features ---
small_len_thresh = in_data['length'].quantile(args.small_quantile)
small_gap_thresh = small_len_thresh / 2

print(f"\n--- Collapsing small features ---")
collapsed_df = collapse_small_clusters(in_data, small_len_thresh, small_gap_thresh)
print(f"  Features before collapse: {len(in_data):,}")
print(f"  Features after collapse: {len(collapsed_df):,}")

# --- Build adjacency matrix ---
edge_mode_str = f" [{args.edge_mode}]"
print(f"\n--- Building adjacency matrix ({args.matrix_type}{edge_mode_str}) ---")

# Get feature data per read (with lengths for weighting)
read_feature_data = collapsed_df.groupby('read').apply(
    lambda x: list(zip(x['feature'], x['length']))
).to_dict()

read_features = {r: [f[0] for f in data] for r, data in read_feature_data.items()}
read_feature_lengths = {r: {f: l for f, l in data} for r, data in read_feature_data.items()}

# Get edges with weights (length of source feature)
def get_weighted_edges(features_with_lengths, edge_mode="directional"):
    """Get edges with weights based on edge counting mode.

    Args:
        features_with_lengths: List of (feature_name, length) tuples in order
        edge_mode: One of:
            - "directional": standard A->B edge counting, weight = source length
            - "bidirectional": A->B and B->A counted separately, reverse uses avg weight
            - "symmetric": edges sorted alphabetically, weight = avg of both features

    Returns:
        List of (from, to, weight) tuples
    """
    if len(features_with_lengths) <= 1:
        return []
    edges = []
    for i in range(len(features_with_lengths) - 1):
        from_feat, from_len = features_with_lengths[i]
        to_feat, to_len = features_with_lengths[i + 1]
        avg_len = (from_len + to_len) / 2

        if edge_mode == "directional":
            edges.append((from_feat, to_feat, from_len))
        elif edge_mode == "bidirectional":
            edges.append((from_feat, to_feat, from_len))
            edges.append((to_feat, from_feat, avg_len))
        elif edge_mode == "symmetric":
            # Sort alphabetically so A->B and B->A both become the same edge
            sorted_pair = tuple(sorted([from_feat, to_feat]))
            edges.append((sorted_pair[0], sorted_pair[1], avg_len))
    return edges

read_edges = {}
read_weighted_edges = {}
for read_name, features in read_features.items():
    read_edges[read_name] = get_edges(features, edge_mode=args.edge_mode)
    read_weighted_edges[read_name] = get_weighted_edges(read_feature_data[read_name], edge_mode=args.edge_mode)

all_features = sorted(collapsed_df['feature'].unique())
all_pairs = []
if args.edge_mode == "symmetric":
    # For symmetric mode, only include A->B where A < B (alphabetically)
    # Both A->B and B->A edges in the data will map to this single column
    for i, f1 in enumerate(all_features):
        for f2 in all_features[i+1:]:
            all_pairs.append(f"{f1}->{f2}")
else:
    # For directional and bidirectional, include all A->B pairs
    for f1 in all_features:
        for f2 in all_features:
            if f1 != f2:
                all_pairs.append(f"{f1}->{f2}")

read_names = sorted(read_features.keys())
pair_to_idx = {pair: i for i, pair in enumerate(all_pairs)}

# Build transition matrix based on matrix type
if args.matrix_type == "binary":
    adj_matrix = np.zeros((len(read_names), len(all_pairs)), dtype=np.float32)
    for i, read_name in enumerate(read_names):
        edges = read_edges[read_name]
        for from_feat, to_feat in edges:
            pair_name = f"{from_feat}->{to_feat}"
            if pair_name in pair_to_idx:
                adj_matrix[i, pair_to_idx[pair_name]] = 1

elif args.matrix_type == "count":
    adj_matrix = np.zeros((len(read_names), len(all_pairs)), dtype=np.float32)
    for i, read_name in enumerate(read_names):
        edges = read_edges[read_name]
        for from_feat, to_feat in edges:
            pair_name = f"{from_feat}->{to_feat}"
            if pair_name in pair_to_idx:
                adj_matrix[i, pair_to_idx[pair_name]] += 1

elif args.matrix_type == "length_weighted":
    adj_matrix = np.zeros((len(read_names), len(all_pairs)), dtype=np.float32)
    for i, read_name in enumerate(read_names):
        edges = read_weighted_edges[read_name]
        read_len = read_length_dict.get(read_name, 1)
        for from_feat, to_feat, weight in edges:
            pair_name = f"{from_feat}->{to_feat}"
            if pair_name in pair_to_idx:
                # Normalize weight by read length
                adj_matrix[i, pair_to_idx[pair_name]] += weight / read_len

print(f"Transition matrix shape: {adj_matrix.shape}")
print(f"Non-zero entries: {np.count_nonzero(adj_matrix):,}")

# --- Add feature abundance if requested ---
if args.include_abundance:
    print(f"\n--- Adding feature abundance dimensions ---")

    # Calculate feature proportions per read
    abundance_matrix = np.zeros((len(read_names), len(all_features)), dtype=np.float32)
    feature_to_idx = {f: i for i, f in enumerate(all_features)}

    for i, read_name in enumerate(read_names):
        read_len = read_length_dict.get(read_name, 1)
        # Sum lengths per feature
        feature_lengths = defaultdict(float)
        for feat, length in read_feature_data[read_name]:
            feature_lengths[feat] += length

        # Convert to proportions
        for feat, total_len in feature_lengths.items():
            if feat in feature_to_idx:
                abundance_matrix[i, feature_to_idx[feat]] = total_len / read_len

    # Concatenate transition matrix with abundance matrix
    adj_matrix = np.hstack([adj_matrix, abundance_matrix])
    print(f"Feature abundance dimensions: {len(all_features)}")
    print(f"Combined matrix shape: {adj_matrix.shape}")

print(f"Final matrix shape: {adj_matrix.shape}")

# --- Hierarchical clustering ---
print(f"\n--- Performing hierarchical clustering ---")
dist_matrix = pdist(adj_matrix, metric='euclidean')
linkage_matrix = linkage(dist_matrix, method=args.linkage_method)

# --- Determine optimal number of clusters ---
print(f"\n--- Analyzing cluster structure ---")

# Pre-compute arrays for faster k-optimization
read_names_arr = np.array(read_names)
read_groups_arr = np.array([read_to_group[r] for r in read_names])

# Fast enrichment check for k-optimization (returns enrichment type and max purity)
def fast_enrichment_check(cluster_mask, read_groups, group_totals_dict, control_grp, all_grps):
    """Fast enrichment calculation for k-optimization loop.

    Returns: (enrichment_label, max_purity) where max_purity is the highest group percentage.
    """
    cluster_groups = read_groups[cluster_mask]
    group_counts = Counter(cluster_groups)
    total_in = len(cluster_groups)
    if total_in == 0:
        return "Mixed", 0.0

    # Calculate max purity (highest percentage from any single group)
    max_purity = max(group_counts.get(g, 0) / total_in for g in all_grps)

    # For two groups, use simple proportion comparison
    if len(all_grps) == 2:
        ctrl_count = group_counts.get(control_grp, 0)
        ctrl_pct = ctrl_count / total_in
        other_grp = [g for g in all_grps if g != control_grp][0]

        # Quick significance check using binomial proportion
        expected_ctrl = group_totals_dict[control_grp] / sum(group_totals_dict.values())
        if ctrl_pct > expected_ctrl + 0.15:  # >15% enrichment threshold for speed
            return f"{control_grp}-enriched", max_purity
        elif ctrl_pct < expected_ctrl - 0.15:
            return f"{other_grp}-enriched", max_purity
    return "Mixed", max_purity

# Calculate clustering metrics for different k values
if args.n_clusters is None:
    # Try different numbers of clusters and evaluate
    # Upper bound: min of max_k or n_reads/10 (need at least 10 reads per cluster on average)
    max_k_limit = min(args.max_k + 1, len(read_names) // 10)
    k_range = list(range(args.min_k, max_k_limit))
    n_samples = len(read_names)

    print(f"Testing cluster counts from {min(k_range)} to {max(k_range)} ({len(k_range)} values)...")

    # Sequential evaluation with early stopping
    cluster_stats = []
    best_composite = -1
    no_improvement_count = 0

    for k in k_range:
        labels = fcluster(linkage_matrix, k, criterion='maxclust')
        labels_arr = np.array(labels)

        # Calculate standard clustering metrics
        if k > 1:
            if n_samples > 2000:
                silhouette = silhouette_score(adj_matrix, labels, sample_size=min(2000, n_samples))
            else:
                silhouette = silhouette_score(adj_matrix, labels)
            calinski_harabasz = calinski_harabasz_score(adj_matrix, labels)
            if args.compute_davies_bouldin:
                davies_bouldin = davies_bouldin_score(adj_matrix, labels)
            else:
                davies_bouldin = np.nan
        else:
            silhouette = 0
            calinski_harabasz = 0
            davies_bouldin = np.nan

        # Count clusters meeting minimum size
        cluster_sizes = Counter(labels)
        valid_clusters = sum(1 for size in cluster_sizes.values() if size >= args.min_cluster_size)

        # Calculate enrichment diversity with strength categories
        enrichments = []
        purities = []
        perfect_enriched = 0
        strong_enriched = 0
        any_enriched = 0

        for cluster_id in range(1, k + 1):
            cluster_mask = labels_arr == cluster_id
            if cluster_mask.sum() >= args.min_cluster_size:
                enrich, purity = fast_enrichment_check(cluster_mask, read_groups_arr,
                                                       group_totals, control_group, all_groups)
                enrichments.append(enrich)
                purities.append(purity)

                if enrich != "Mixed":
                    any_enriched += 1
                    if purity >= args.perfect_threshold:
                        perfect_enriched += 1
                    if purity >= args.strong_threshold:
                        strong_enriched += 1

        # Count enrichments by group
        group_enriched_counts = {g: sum(1 for e in enrichments if e == f"{g}-enriched") for g in all_groups}
        mixed = enrichments.count("Mixed")

        # Calculate enrichment ratios
        enriched_ratio = any_enriched / valid_clusters if valid_clusters > 0 else 0
        perfect_ratio = perfect_enriched / valid_clusters if valid_clusters > 0 else 0
        strong_ratio = strong_enriched / valid_clusters if valid_clusters > 0 else 0

        # Composite score
        silhouette_norm = (silhouette + 1) / 2
        composite_score = (0.3 * silhouette_norm +
                          0.3 * enriched_ratio +
                          0.2 * strong_ratio +
                          0.2 * perfect_ratio)

        stats = {
            'k': k,
            'silhouette': silhouette,
            'calinski_harabasz': calinski_harabasz,
            'davies_bouldin': davies_bouldin,
            'valid_clusters': valid_clusters,
            'any_enriched': any_enriched,
            'strong_enriched': strong_enriched,
            'perfect_enriched': perfect_enriched,
            'enriched_ratio': enriched_ratio,
            'strong_ratio': strong_ratio,
            'perfect_ratio': perfect_ratio,
            'composite_score': composite_score,
            **{f'{g}_enriched': group_enriched_counts.get(g, 0) for g in all_groups},
            'mixed': mixed
        }
        cluster_stats.append(stats)

        # Early stopping check
        if stats['composite_score'] > best_composite:
            best_composite = stats['composite_score']
            no_improvement_count = 0
        else:
            no_improvement_count += 1

        if args.early_stopping > 0 and no_improvement_count >= args.early_stopping:
            print(f"  Early stopping at k={k} (no improvement for {args.early_stopping} iterations)")
            break

    # Save cluster analysis
    stats_df = pd.DataFrame(cluster_stats)
    stats_file = f"{args.output_prefix}.cluster_k_analysis.tsv"
    stats_df.to_csv(stats_file, sep='\t', index=False)
    print(f"  Saved k analysis to: {stats_file}")

    # Generate k-selection diagnostic plot (2x3 grid)
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    # 1. Silhouette score (higher is better)
    ax1 = axes[0, 0]
    ax1.plot(stats_df['k'], stats_df['silhouette'], 'b-o', markersize=3)
    ax1.set_xlabel('Number of clusters (k)')
    ax1.set_ylabel('Silhouette Score')
    ax1.set_title('Silhouette Score (higher = better)')
    best_silhouette_k = int(stats_df.loc[stats_df['silhouette'].idxmax(), 'k'])
    ax1.axvline(x=best_silhouette_k, color='g', linestyle='--', alpha=0.5, label=f'Best k={best_silhouette_k}')
    ax1.legend()

    # 2. Calinski-Harabasz index (higher is better)
    ax2 = axes[0, 1]
    ax2.plot(stats_df['k'], stats_df['calinski_harabasz'], 'b-o', markersize=3)
    ax2.set_xlabel('Number of clusters (k)')
    ax2.set_ylabel('Calinski-Harabasz Index')
    ax2.set_title('Calinski-Harabasz (higher = better)')
    best_ch_k = int(stats_df.loc[stats_df['calinski_harabasz'].idxmax(), 'k'])
    ax2.axvline(x=best_ch_k, color='g', linestyle='--', alpha=0.5, label=f'Best k={best_ch_k}')
    ax2.legend()

    # 3. Enrichment ratios
    ax3 = axes[0, 2]
    ax3.plot(stats_df['k'], stats_df['enriched_ratio'], '-o', markersize=3, label='Any enriched', color='blue')
    ax3.plot(stats_df['k'], stats_df['strong_ratio'], '-o', markersize=3, label=f'Strong (>={int(args.strong_threshold*100)}%)', color='orange')
    ax3.plot(stats_df['k'], stats_df['perfect_ratio'], '-o', markersize=3, label=f'Perfect (>={int(args.perfect_threshold*100)}%)', color='red')
    ax3.set_xlabel('Number of clusters (k)')
    ax3.set_ylabel('Ratio (enriched / valid clusters)')
    ax3.set_title('Enrichment Ratios')
    ax3.set_ylim(0, 1)
    ax3.legend()

    # 4. Enrichment counts (absolute)
    ax4 = axes[1, 0]
    ax4.plot(stats_df['k'], stats_df['any_enriched'], '-o', markersize=3, label='Any enriched', color='blue')
    ax4.plot(stats_df['k'], stats_df['strong_enriched'], '-o', markersize=3, label=f'Strong (>={int(args.strong_threshold*100)}%)', color='orange')
    ax4.plot(stats_df['k'], stats_df['perfect_enriched'], '-o', markersize=3, label=f'Perfect (>={int(args.perfect_threshold*100)}%)', color='red')
    ax4.set_xlabel('Number of clusters (k)')
    ax4.set_ylabel('Number of clusters')
    ax4.set_title('Enrichment Counts (absolute)')
    ax4.legend()

    # 5. Enrichment by group
    ax5 = axes[1, 1]
    for g in all_groups:
        col = f'{g}_enriched'
        if col in stats_df.columns:
            ax5.plot(stats_df['k'], stats_df[col], '-o', markersize=3, label=f'{g}-enriched', color=group_colors.get(g, None))
    ax5.plot(stats_df['k'], stats_df['mixed'], '-o', markersize=3, label='Mixed', color='gray')
    ax5.set_xlabel('Number of clusters (k)')
    ax5.set_ylabel('Number of clusters')
    ax5.set_title('Enrichment by Group')
    ax5.legend()

    # 6. Composite score (our recommended metric)
    ax6 = axes[1, 2]
    ax6.plot(stats_df['k'], stats_df['composite_score'], 'b-o', markersize=3)
    ax6.set_xlabel('Number of clusters (k)')
    ax6.set_ylabel('Composite Score')
    ax6.set_title('Composite Score (silhouette + enrichment)')
    best_composite_k = int(stats_df.loc[stats_df['composite_score'].idxmax(), 'k'])
    ax6.axvline(x=best_composite_k, color='g', linestyle='--', alpha=0.5, label=f'Best k={best_composite_k}')
    ax6.axhline(y=stats_df['composite_score'].max(), color='r', linestyle='--', alpha=0.3)
    ax6.legend()

    plt.tight_layout()
    k_plot_file = f"{args.output_prefix}.k_selection.pdf"
    plt.savefig(k_plot_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved k-selection plot to: {k_plot_file}")

    # Find best k by each metric
    print(f"\n  Optimal k by metric:")
    print(f"    Silhouette:        k={best_silhouette_k} (score={stats_df['silhouette'].max():.4f})")
    print(f"    Calinski-Harabasz: k={best_ch_k} (score={stats_df['calinski_harabasz'].max():.1f})")
    if args.compute_davies_bouldin:
        best_db_k = int(stats_df.loc[stats_df['davies_bouldin'].idxmin(), 'k'])
        print(f"    Davies-Bouldin:    k={best_db_k} (score={stats_df['davies_bouldin'].min():.4f})")
    print(f"    Composite:         k={best_composite_k} (score={stats_df['composite_score'].max():.4f})")

    # Report enrichment stats at best composite k
    best_row = stats_df[stats_df['k'] == best_composite_k].iloc[0]
    print(f"\n  At k={best_composite_k}:")
    print(f"    Valid clusters:    {int(best_row['valid_clusters'])}")
    print(f"    Any enriched:      {int(best_row['any_enriched'])} ({best_row['enriched_ratio']*100:.1f}%)")
    print(f"    Strong enriched:   {int(best_row['strong_enriched'])} ({best_row['strong_ratio']*100:.1f}%)")
    print(f"    Perfect enriched:  {int(best_row['perfect_enriched'])} ({best_row['perfect_ratio']*100:.1f}%)")

    print(f"\n  Using k={best_composite_k} (based on composite score)")

    n_clusters = best_composite_k
else:
    n_clusters = args.n_clusters
    print(f"  Using specified k: {n_clusters}")

# --- Cut tree and analyze clusters ---
print(f"\n--- Cutting tree into {n_clusters} clusters ---")
cluster_labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')

# Analyze each cluster
cluster_analysis = []
cluster_reads_dict = {}
read_centroid_distances = {}  # Store centroid distance for each read

# Pre-compute arrays for faster cluster analysis
cluster_labels_arr = np.array(cluster_labels)
read_samples_arr = np.array([read_to_sample[r] for r in read_names])

for cluster_id in range(1, n_clusters + 1):
    # Use numpy boolean indexing (faster than list comprehension)
    cluster_mask = cluster_labels_arr == cluster_id
    if cluster_mask.sum() < args.min_cluster_size:
        continue

    cluster_indices = np.where(cluster_mask)[0]
    cluster_reads = read_names_arr[cluster_mask].tolist()
    cluster_samples = read_samples_arr[cluster_mask].tolist()

    cluster_reads_dict[cluster_id] = cluster_reads

    # Calculate enrichment using appropriate method
    if args.comparison_mode == "two-group":
        stats = calculate_enrichment_two_group(cluster_samples, sample_to_group, control_group, group_totals)
    elif args.comparison_mode == "multi-group":
        stats = calculate_enrichment_multi_group(cluster_samples, sample_to_group, group_totals)
    else:  # per-sample
        stats = calculate_enrichment_per_sample(cluster_samples, sample_totals)

    # Find centroid (read closest to cluster mean)
    cluster_matrix = adj_matrix[cluster_indices]
    centroid = cluster_matrix.mean(axis=0)
    distances_to_centroid = np.linalg.norm(cluster_matrix - centroid, axis=1)
    centroid_idx = np.argmin(distances_to_centroid)
    centroid_read = cluster_reads[centroid_idx]

    # Store centroid distances for all reads in this cluster
    for i, read in enumerate(cluster_reads):
        read_centroid_distances[read] = distances_to_centroid[i]

    # Build cluster analysis record
    cluster_record = {
        'cluster_id': cluster_id,
        'size': len(cluster_reads),
        'odds_ratio': stats['odds_ratio'],
        'p_value': stats['p_value'],
        'enrichment': stats['enrichment'],
        'centroid_read': centroid_read,
        'centroid_sample': read_to_sample[centroid_read],
        'centroid_group': sample_to_group.get(read_to_sample[centroid_read], read_to_sample[centroid_read])
    }

    # Add group-specific columns
    for g in all_groups:
        cluster_record[f'{g}_count'] = stats['group_counts'].get(g, 0)
        cluster_record[f'{g}_pct'] = stats['group_pcts'].get(g, 0)

    # Option A: Nested/hierarchical testing within enriched group
    if args.nested and stats['enrichment'] != 'mixed' and args.comparison_mode != 'per-sample':
        enriched_group = stats['enrichment'].replace('-enriched', '')
        nested_results = calculate_nested_within_group(
            cluster_samples, enriched_group, sample_to_group, sample_totals
        )
        cluster_record['nested_results'] = nested_results
        # Add summary columns
        enriched_samples = [s for s, r in nested_results.items() if r['enriched']]
        depleted_samples = [s for s, r in nested_results.items() if r['depleted']]
        cluster_record['nested_enriched_samples'] = ','.join(enriched_samples) if enriched_samples else ''
        cluster_record['nested_depleted_samples'] = ','.join(depleted_samples) if depleted_samples else ''
        cluster_record['nested_n_enriched'] = len(enriched_samples)
        cluster_record['nested_n_depleted'] = len(depleted_samples)

    # Option B: Per-sample tests for all samples
    if args.also_test_samples and args.comparison_mode != 'per-sample':
        sample_test_results = calculate_all_sample_tests(cluster_samples, sample_totals)
        cluster_record['sample_test_results'] = sample_test_results
        # Add summary columns for each sample
        for sample, result in sample_test_results.items():
            cluster_record[f'sample_{sample}_count'] = result['count']
            cluster_record[f'sample_{sample}_pct'] = result['pct']
            cluster_record[f'sample_{sample}_pval'] = result['p_value']
            cluster_record[f'sample_{sample}_enriched'] = result['enriched']

    # Option C: Stratified variance metrics
    if args.stratified and args.comparison_mode != 'per-sample':
        stratified_results = calculate_stratified_variance(
            cluster_samples, sample_to_group, sample_totals, group_totals
        )
        cluster_record['stratified_results'] = stratified_results
        # Add summary columns per group
        for g, result in stratified_results.items():
            cluster_record[f'{g}_cv'] = result['cv_normalized']
            cluster_record[f'{g}_consistent'] = result['consistent']
            cluster_record[f'{g}_sample_counts'] = ','.join(f"{s}:{c}" for s, c in result['counts'].items())

    cluster_analysis.append(cluster_record)

# Sort by enrichment type and p-value
cluster_df = pd.DataFrame(cluster_analysis)
cluster_df = cluster_df.sort_values(['enrichment', 'p_value'])

# --- Output results ---
print(f"\n--- Cluster Summary ---")

# Build dynamic header based on groups
header_parts = ['Cluster', 'Size']
for g in all_groups:
    header_parts.extend([g[:6], f'{g[:4]}%'])  # Truncate long group names
header_parts.extend(['P-value', 'Enrichment', 'Centroid'])
header_fmt = '{:<8} {:<6} ' + ' '.join(['{:<6}'] * (len(all_groups) * 2)) + ' {:<10} {:<20} {:<10}'
print(header_fmt.format(*header_parts))
print("-" * (95 + 12 * len(all_groups)))

for _, row in cluster_df.iterrows():
    centroid_group = row.get('centroid_group', sample_to_group.get(row['centroid_sample'], row['centroid_sample']))
    # Flag if centroid disagrees with enrichment
    flag = ""
    expected_group = row['enrichment'].replace('-enriched', '') if '-enriched' in row['enrichment'] else None
    if expected_group and centroid_group != expected_group:
        flag = " ⚠️"

    row_values = [row['cluster_id'], row['size']]
    for g in all_groups:
        row_values.append(int(row.get(f'{g}_count', 0)))
        row_values.append(f"{row.get(f'{g}_pct', 0):.1f}")
    row_values.extend([f"{row['p_value']:.2e}", row['enrichment'][:18], f"{centroid_group}{flag}"])

    row_fmt = '{:<8} {:<6} ' + ' '.join(['{:<6}'] * (len(all_groups) * 2)) + ' {:<10} {:<20} {:<10}'
    print(row_fmt.format(*row_values))

# Save cluster analysis
analysis_file = f"{args.output_prefix}.cluster_analysis.tsv"
# Drop columns containing nested dictionaries (can't be serialized to TSV)
cols_to_drop = [c for c in cluster_df.columns if c in ['nested_results', 'sample_test_results', 'stratified_results']]
cluster_df_export = cluster_df.drop(columns=cols_to_drop, errors='ignore')
cluster_df_export.to_csv(analysis_file, sep='\t', index=False)
print(f"\n  Saved cluster analysis to: {analysis_file}")

# Save read assignments with stats (sorted by cluster, then centroid distance)
# Build list of read records
read_records = []
for i, read in enumerate(read_names):
    cluster = cluster_labels[i]
    sample = read_to_sample[read]
    group = sample_to_group.get(sample, sample)
    centroid_dist = read_centroid_distances.get(read, np.nan)
    read_len = read_length_dict.get(read, 0)
    read_records.append({
        'read': read,
        'cluster': cluster,
        'sample': sample,
        'group': group,
        'centroid_distance': centroid_dist,
        'read_length': read_len
    })

# Create DataFrame and sort by cluster, then centroid_distance
assignments = pd.DataFrame(read_records)
assignments = assignments.sort_values(['cluster', 'centroid_distance'], ascending=[True, True])

# Add rank within cluster (1 = closest to centroid)
assignments['rank'] = assignments.groupby('cluster').cumcount() + 1

# Reorder columns
assignments = assignments[['read', 'cluster', 'sample', 'group', 'centroid_distance', 'read_length', 'rank']]

assignments_file = f"{args.output_prefix}.read_assignments.tsv"
assignments.to_csv(assignments_file, sep='\t', index=False)
print(f"  Saved read assignments to: {assignments_file}")

# Save feature matrix data for cluster_plot.py (for computing subset dendrograms)
matrix_file = f"{args.output_prefix}.feature_matrix.npz"
np.savez(matrix_file,
         adj_matrix=adj_matrix,
         read_names=np.array(read_names),
         cluster_labels=cluster_labels,
         linkage_method=args.linkage_method)
print(f"  Saved feature matrix to: {matrix_file}")

# Save sample metadata for cluster_plot.py
metadata_out_file = f"{args.output_prefix}.sample_metadata.tsv"
meta_records = []
for sample in sample_labels:
    group = sample_to_group.get(sample, sample)
    color = sample_colors.get(sample, '#999999')
    meta_records.append({'sample': sample, 'group': group, 'color': color})
meta_out_df = pd.DataFrame(meta_records)
meta_out_df.to_csv(metadata_out_file, sep='\t', index=False)
print(f"  Saved sample metadata to: {metadata_out_file}")

# --- Generate visualization ---
print(f"\n--- Generating cluster visualization ---")

fig, axes = plt.subplots(2, 2, figsize=(16, 14))
from scipy.cluster.hierarchy import set_link_color_palette
from matplotlib.patches import Patch

# Color palette for clusters
colors = plt.cm.tab20(np.linspace(0, 1, n_clusters))
set_link_color_palette([matplotlib.colors.rgb2hex(c) for c in colors])

# Create cluster color map
cluster_color_map = {c: matplotlib.colors.rgb2hex(colors[i % len(colors)]) for i, c in enumerate(sorted(set(cluster_labels)))}

# Create cluster to enrichment mapping
cluster_to_enrichment = dict(zip(cluster_df['cluster_id'], cluster_df['enrichment']))

# 1. Linear Dendrogram with cluster coloring
ax1 = axes[0, 0]
dend = dendrogram(
    linkage_matrix,
    ax=ax1,
    no_labels=True,
    color_threshold=linkage_matrix[-(n_clusters-1), 2] if n_clusters > 1 else 0
)
ax1.set_title(f"Dendrogram (k={n_clusters} clusters)")
ax1.set_ylabel("Distance")
ax1.axhline(y=linkage_matrix[-(n_clusters-1), 2] if n_clusters > 1 else 0,
            color='red', linestyle='--', alpha=0.5, label='Cut height')
ax1.legend()

# 2. Cluster size distribution
ax2 = axes[0, 1]
cluster_sizes = cluster_df['size'].values
cluster_ids = cluster_df['cluster_id'].values
enrichments = cluster_df['enrichment'].values

# Map enrichments to colors using group_colors
def get_enrichment_color(e):
    for g, color in group_colors.items():
        if e == f'{g}-enriched':
            return color
    return '#999999'  # Mixed

colors_bar = [get_enrichment_color(e) for e in enrichments]
bars = ax2.bar(range(len(cluster_sizes)), cluster_sizes, color=colors_bar)
ax2.set_xlabel("Cluster")
ax2.set_ylabel("Size (reads)")
ax2.set_title("Cluster Sizes by Enrichment")
ax2.set_xticks(range(len(cluster_sizes)))
ax2.set_xticklabels([f"C{c}" for c in cluster_ids], rotation=45)

# Legend - dynamic based on groups
legend_patches = [Patch(facecolor=c, label=f'{g}-enriched') for g, c in group_colors.items()]
legend_patches.append(Patch(facecolor='#999999', label='Mixed'))
ax2.legend(handles=legend_patches, loc='upper right')

# 3. Group composition per cluster
ax3 = axes[1, 0]
x = np.arange(len(cluster_df))
n_groups = len(all_groups)
width = 0.8 / n_groups

for i, g in enumerate(all_groups):
    offset = (i - n_groups / 2 + 0.5) * width
    pct_col = f'{g}_pct'
    if pct_col in cluster_df.columns:
        ax3.bar(x + offset, cluster_df[pct_col], width, label=g, color=group_colors[g])

ax3.set_xlabel("Cluster")
ax3.set_ylabel("Percentage")
ax3.set_title("Group Composition per Cluster")
ax3.set_xticks(x)
ax3.set_xticklabels([f"C{c}" for c in cluster_df['cluster_id']], rotation=45)
ax3.legend()
if n_groups == 2:
    ax3.axhline(y=50, color='gray', linestyle='--', alpha=0.5)

# 4. P-value distribution
ax4 = axes[1, 1]
pvals = cluster_df['p_value'].values
colors_pval = [get_enrichment_color(e) for e in cluster_df['enrichment']]
ax4.bar(range(len(pvals)), -np.log10(pvals + 1e-300), color=colors_pval)
ax4.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='p=0.05')
ax4.axhline(y=-np.log10(0.01), color='orange', linestyle='--', alpha=0.5, label='p=0.01')
ax4.set_xlabel("Cluster")
ax4.set_ylabel("-log10(p-value)")
ax4.set_title("Enrichment Significance")
ax4.set_xticks(range(len(pvals)))
ax4.set_xticklabels([f"C{c}" for c in cluster_df['cluster_id']], rotation=45)
ax4.legend()

plt.tight_layout()
plot_file = f"{args.output_prefix}.cluster_analysis.pdf"
plt.savefig(plot_file, dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved cluster visualization to: {plot_file}")

# --- Circular Dendrogram (separate output) ---
if args.plot_circular_dendrogram:
    print(f"\n--- Generating circular dendrogram ---")

    # Get dendrogram structure
    dend_data = dendrogram(linkage_matrix, no_plot=True)
    dend_colored = dendrogram(
        linkage_matrix,
        no_plot=True,
        color_threshold=linkage_matrix[-(n_clusters-1), 2] if n_clusters > 1 else 0
    )

    n_leaves = len(dend_data['leaves'])
    max_dist = max(max(y) for y in dend_data['dcoord'])
    leaf_order = dend_data['leaves']

    # Prepare annotation data for each leaf (in dendrogram order)
    leaf_samples = [read_to_sample[read_names[i]] for i in leaf_order]
    leaf_clusters = [cluster_labels[i] for i in leaf_order]
    leaf_enrichments = [cluster_to_enrichment.get(c, 'Mixed') for c in leaf_clusters]

    # Color mappings for annotations
    enrichment_colors = {f'{g}-enriched': c for g, c in group_colors.items()}
    enrichment_colors['Mixed'] = '#CCCCCC'

    # Create figure with circular dendrogram
    fig_circ = plt.figure(figsize=(16, 14))
    ax_circ = fig_circ.add_subplot(111, polar=True)

    # Plot each link in polar coordinates
    for xcoord, ycoord, color in zip(dend_colored['icoord'], dend_colored['dcoord'], dend_colored['color_list']):
        # Convert x coordinates to angles (0 to 2*pi)
        theta = [2 * np.pi * x / (n_leaves * 10) for x in xcoord]
        # Convert y (distance) to radius (invert so root is at center)
        r = [max_dist - y + max_dist * 0.1 for y in ycoord]
        ax_circ.plot(theta, r, color=color, linewidth=0.5)

    # Calculate theta positions for leaves
    theta_leaves = np.linspace(0, 2 * np.pi, n_leaves, endpoint=False)

    # Ring 1 (innermost): Sample of origin
    ring1_bottom = max_dist * 1.12
    ring1_height = max_dist * 0.06
    for theta, sample in zip(theta_leaves, leaf_samples):
        ax_circ.bar(theta, ring1_height, width=2 * np.pi / n_leaves, bottom=ring1_bottom,
                    color=sample_colors.get(sample, '#999999'), alpha=0.9, edgecolor='none')

    # Ring 2 (middle): Cluster number
    ring2_bottom = max_dist * 1.20
    ring2_height = max_dist * 0.06
    for theta, cluster in zip(theta_leaves, leaf_clusters):
        ax_circ.bar(theta, ring2_height, width=2 * np.pi / n_leaves, bottom=ring2_bottom,
                    color=cluster_color_map.get(cluster, '#999999'), alpha=0.9, edgecolor='none')

    # Ring 3 (outermost): Cluster enrichment
    ring3_bottom = max_dist * 1.28
    ring3_height = max_dist * 0.06
    for theta, enrich in zip(theta_leaves, leaf_enrichments):
        ax_circ.bar(theta, ring3_height, width=2 * np.pi / n_leaves, bottom=ring3_bottom,
                    color=enrichment_colors.get(enrich, '#CCCCCC'), alpha=0.9, edgecolor='none')

    ax_circ.set_title(f"Circular Dendrogram (k={n_clusters} clusters)\nRings: Sample | Cluster | Enrichment", pad=30, fontsize=12)
    ax_circ.set_yticklabels([])
    ax_circ.set_xticklabels([])
    ax_circ.grid(False)

    # Add legends
    # Sample legend (inner ring)
    sample_patches = [Patch(facecolor=sample_colors[s], label=s) for s in sorted(sample_colors.keys())]
    leg1 = ax_circ.legend(handles=sample_patches, loc='upper left', bbox_to_anchor=(1.05, 1.0), title='Sample (inner)')

    # Cluster legend (middle ring)
    cluster_patches = [Patch(facecolor=cluster_color_map[c], label=f'C{c}') for c in sorted(cluster_color_map.keys())]
    leg2 = ax_circ.legend(handles=cluster_patches, loc='upper left', bbox_to_anchor=(1.05, 0.75),
                          title='Cluster (middle)', ncol=2, fontsize=7)
    ax_circ.add_artist(leg1)

    # Enrichment legend (outer ring)
    enrich_patches = [Patch(facecolor=enrichment_colors[e], label=e) for e in sorted(enrichment_colors.keys())]
    leg3 = ax_circ.legend(handles=enrich_patches, loc='upper left', bbox_to_anchor=(1.05, 0.35), title='Enrichment (outer)')
    ax_circ.add_artist(leg2)

    plt.tight_layout()
    circ_dend_file = f"{args.output_prefix}.circular_dendrogram.pdf"
    plt.savefig(circ_dend_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved circular dendrogram to: {circ_dend_file}")

# --- UMAP Visualization ---
if args.plot_umap and len(read_names) > 10:
    print(f"\n--- Generating UMAP visualization ---")
    try:
        import umap

        # Fit UMAP
        reducer = umap.UMAP(
            n_neighbors=min(args.umap_neighbors, len(read_names) - 1),
            min_dist=args.umap_min_dist,
            metric='euclidean',
            random_state=42
        )
        embedding = reducer.fit_transform(adj_matrix)
        print(f"  UMAP parameters: n_neighbors={args.umap_neighbors}, min_dist={args.umap_min_dist}")

        # Create read-to-cluster mapping
        read_to_cluster = dict(zip(read_names, cluster_labels))
        cluster_to_enrichment = dict(zip(cluster_df['cluster_id'], cluster_df['enrichment']))

        # Prepare data for plotting
        sample_list = [read_to_sample[r] for r in read_names]
        group_list = [sample_to_group.get(s, s) for s in sample_list]
        cluster_list = [read_to_cluster[r] for r in read_names]

        # Color mappings
        enrichment_colors = {f'{g}-enriched': c for g, c in group_colors.items()}
        enrichment_colors['Mixed'] = '#CCCCCC'

        # Generate distinct colors for clusters
        n_unique_clusters = len(set(cluster_list))
        cluster_cmap = plt.cm.tab20(np.linspace(0, 1, max(20, n_unique_clusters)))
        cluster_color_map = {c: matplotlib.colors.rgb2hex(cluster_cmap[i % 20]) for i, c in enumerate(sorted(set(cluster_list)))}

        # Create 2x2 UMAP plot (group, enrichment, cluster numbers, cluster colors)
        fig, axes = plt.subplots(2, 2, figsize=(16, 14))

        # 1. Top-left: Colored by group (from metadata)
        point_colors_group = [group_colors.get(g, '#999999') for g in group_list]
        axes[0, 0].scatter(embedding[:, 0], embedding[:, 1], c=point_colors_group, s=15, alpha=0.6)
        axes[0, 0].set_title("UMAP - Colored by Group")
        axes[0, 0].set_xlabel("UMAP 1")
        axes[0, 0].set_ylabel("UMAP 2")
        group_patches = [Patch(facecolor=group_colors[g], label=g) for g in sorted(group_colors.keys())]
        axes[0, 0].legend(handles=group_patches, loc='upper right')

        # 2. Top-right: Colored by enrichment
        point_colors_enrich = []
        for r in read_names:
            cluster = read_to_cluster[r]
            enrich = cluster_to_enrichment.get(cluster, 'Mixed')
            point_colors_enrich.append(enrichment_colors.get(enrich, '#CCCCCC'))
        axes[0, 1].scatter(embedding[:, 0], embedding[:, 1], c=point_colors_enrich, s=15, alpha=0.6)
        axes[0, 1].set_title("UMAP - Colored by Cluster Enrichment")
        axes[0, 1].set_xlabel("UMAP 1")
        axes[0, 1].set_ylabel("UMAP 2")
        enrich_patches = [Patch(facecolor=c, label=e) for e, c in enrichment_colors.items()]
        axes[0, 1].legend(handles=enrich_patches, loc='upper right')

        # 3. Bottom-left: Colored by cluster with cluster number labels
        point_colors_cluster = [cluster_color_map[c] for c in cluster_list]
        axes[1, 0].scatter(embedding[:, 0], embedding[:, 1], c=point_colors_cluster, s=15, alpha=0.6)
        axes[1, 0].set_title("UMAP - Colored by Cluster")
        axes[1, 0].set_xlabel("UMAP 1")
        axes[1, 0].set_ylabel("UMAP 2")

        # Add cluster number labels at cluster centroids
        for cluster_id in sorted(set(cluster_list)):
            mask = [c == cluster_id for c in cluster_list]
            cluster_points = embedding[mask]
            centroid_x = np.mean(cluster_points[:, 0])
            centroid_y = np.mean(cluster_points[:, 1])
            axes[1, 0].annotate(str(cluster_id), (centroid_x, centroid_y),
                               fontsize=8, fontweight='bold', ha='center', va='center',
                               bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7, edgecolor='gray'))

        # 4. Bottom-right: Cluster numbers only (no points, just labels for clarity)
        axes[1, 1].scatter(embedding[:, 0], embedding[:, 1], c=point_colors_cluster, s=10, alpha=0.3)
        axes[1, 1].set_title("UMAP - Cluster Labels")
        axes[1, 1].set_xlabel("UMAP 1")
        axes[1, 1].set_ylabel("UMAP 2")

        # Add cluster labels with enrichment info
        for cluster_id in sorted(set(cluster_list)):
            mask = [c == cluster_id for c in cluster_list]
            cluster_points = embedding[mask]
            centroid_x = np.mean(cluster_points[:, 0])
            centroid_y = np.mean(cluster_points[:, 1])
            enrich = cluster_to_enrichment.get(cluster_id, 'Mixed')
            enrich_short = enrich.replace('-enriched', '').replace('Mixed', 'M')[:4]
            label = f"C{cluster_id}\n({enrich_short})"
            axes[1, 1].annotate(label, (centroid_x, centroid_y),
                               fontsize=7, ha='center', va='center',
                               bbox=dict(boxstyle='round,pad=0.2', facecolor=enrichment_colors.get(enrich, '#CCCCCC'),
                                        alpha=0.8, edgecolor='gray'))

        plt.tight_layout()
        umap_file = f"{args.output_prefix}.umap.pdf"
        plt.savefig(umap_file, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved UMAP plot: {umap_file}")

        # Save UMAP coordinates
        umap_coords = pd.DataFrame({
            'read': read_names,
            'umap_1': embedding[:, 0],
            'umap_2': embedding[:, 1],
            'sample': sample_list,
            'cluster': cluster_labels,
            'enrichment': [cluster_to_enrichment.get(c, 'Mixed') for c in cluster_labels]
        })
        umap_coords_file = f"{args.output_prefix}.umap_coordinates.tsv"
        umap_coords.to_csv(umap_coords_file, sep='\t', index=False)
        print(f"  Saved UMAP coordinates: {umap_coords_file}")

        # Try to generate interactive Plotly version
        try:
            import plotly.graph_objects as go
            from plotly.subplots import make_subplots

            # Create figure with dropdown menu for different colorings
            fig = go.Figure()

            # Prepare data arrays
            embedding_x = embedding[:, 0]
            embedding_y = embedding[:, 1]
            enrichment_list = [cluster_to_enrichment.get(c, 'Mixed') for c in cluster_list]

            # Hover text (same for all views)
            hover_text = [f"Read: {r[:20]}...<br>Group: {g}<br>Cluster: {c}<br>Enrichment: {e}"
                         for r, g, c, e in zip(read_names, group_list, cluster_list, enrichment_list)]

            # 1. Add traces for GROUP coloring (default view)
            for group in sorted(set(group_list)):
                mask = np.array([g == group for g in group_list])
                fig.add_trace(go.Scatter(
                    x=embedding_x[mask],
                    y=embedding_y[mask],
                    mode='markers',
                    name=group,
                    text=[hover_text[i] for i in range(len(mask)) if mask[i]],
                    hoverinfo='text',
                    marker=dict(size=6, color=group_colors.get(group, '#999999'), opacity=0.7),
                    visible=True
                ))

            n_group_traces = len(set(group_list))

            # 2. Add traces for CLUSTER coloring
            for cluster_id in sorted(set(cluster_list)):
                mask = np.array([c == cluster_id for c in cluster_list])
                enrich = cluster_to_enrichment.get(cluster_id, 'Mixed')
                fig.add_trace(go.Scatter(
                    x=embedding_x[mask],
                    y=embedding_y[mask],
                    mode='markers',
                    name=f"C{cluster_id} ({enrich[:4]})",
                    text=[hover_text[i] for i in range(len(mask)) if mask[i]],
                    hoverinfo='text',
                    marker=dict(size=6, color=cluster_color_map.get(cluster_id, '#999999'), opacity=0.7),
                    visible=False
                ))

            n_cluster_traces = len(set(cluster_list))

            # 3. Add traces for ENRICHMENT coloring
            for enrich in sorted(set(enrichment_list)):
                mask = np.array([e == enrich for e in enrichment_list])
                fig.add_trace(go.Scatter(
                    x=embedding_x[mask],
                    y=embedding_y[mask],
                    mode='markers',
                    name=enrich,
                    text=[hover_text[i] for i in range(len(mask)) if mask[i]],
                    hoverinfo='text',
                    marker=dict(size=6, color=enrichment_colors.get(enrich, '#CCCCCC'), opacity=0.7),
                    visible=False
                ))

            n_enrichment_traces = len(set(enrichment_list))

            # Create visibility arrays for dropdown
            total_traces = n_group_traces + n_cluster_traces + n_enrichment_traces

            vis_group = [True] * n_group_traces + [False] * n_cluster_traces + [False] * n_enrichment_traces
            vis_cluster = [False] * n_group_traces + [True] * n_cluster_traces + [False] * n_enrichment_traces
            vis_enrichment = [False] * n_group_traces + [False] * n_cluster_traces + [True] * n_enrichment_traces

            # Add dropdown menu
            fig.update_layout(
                updatemenus=[
                    dict(
                        active=0,
                        buttons=[
                            dict(label="Color by Group",
                                 method="update",
                                 args=[{"visible": vis_group},
                                       {"title": f"UMAP - Colored by Group ({len(read_names):,} reads)"}]),
                            dict(label="Color by Cluster",
                                 method="update",
                                 args=[{"visible": vis_cluster},
                                       {"title": f"UMAP - Colored by Cluster ({len(read_names):,} reads)"}]),
                            dict(label="Color by Enrichment",
                                 method="update",
                                 args=[{"visible": vis_enrichment},
                                       {"title": f"UMAP - Colored by Enrichment ({len(read_names):,} reads)"}]),
                        ],
                        direction="down",
                        showactive=True,
                        x=0.0,
                        xanchor="left",
                        y=1.15,
                        yanchor="top"
                    )
                ],
                title=f"UMAP - Colored by Group ({len(read_names):,} reads)",
                xaxis_title="UMAP 1",
                yaxis_title="UMAP 2",
                hovermode='closest',
                width=1000,
                height=750,
                legend=dict(yanchor="top", y=0.99, xanchor="left", x=1.02)
            )

            umap_html_file = f"{args.output_prefix}.umap.html"
            fig.write_html(umap_html_file)
            print(f"  Saved interactive UMAP: {umap_html_file}")

        except ImportError:
            print("  Note: Install plotly for interactive UMAP (pip install plotly)")

    except ImportError:
        print("  Warning: umap-learn not installed. Skipping UMAP plot.")
        print("  Install with: pip install umap-learn")

# --- Summary ---
print(f"\n" + "=" * 60)
print("Summary")
print("=" * 60)
print(f"Total reads: {len(read_names):,}")
print(f"Number of clusters: {n_clusters}")
print(f"Valid clusters (size >= {args.min_cluster_size}): {len(cluster_df)}")
for g in all_groups:
    enriched_label = f'{g}-enriched'
    count = sum(cluster_df['enrichment'] == enriched_label)
    print(f"  - {enriched_label}: {count}")
print(f"  - Mixed: {sum(cluster_df['enrichment'] == 'Mixed')}")
print(f"\nOutput files:")
print(f"  - {analysis_file}")
print(f"  - {assignments_file}")
print(f"  - {matrix_file}")
print(f"  - {plot_file}")
print(f"\nNext step: Use read_assignments.tsv with KaryoScope_cluster_plot.py")
print(f"to visualize reads from each cluster (sorted by centroid distance).")
print(f"Use --feature-matrix with the .feature_matrix.npz file for dendrogram header.")
