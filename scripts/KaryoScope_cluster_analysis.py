# KaryoScope Cluster Analysis
# Analyzes hierarchical clustering results to identify biologically interesting clusters
# and select representative reads for visualization
#
# Usage:
# python KaryoScope_cluster_analysis.py \
#   --bed /path/to/sample1.bed.gz /path/to/sample2.bed.gz \
#   --sample-metadata samples.tsv \
#   --output-prefix analysis_output \
#   --n-clusters 10
#
# Sample metadata file format (TSV):
#   sample    group    color
#   SW26_Pre  control  #377EB8
#   SW26_Post treatment #E41A1C
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
parser.add_argument("--output-prefix", dest="output_prefix", required=True,
                    help="Prefix for output files")
parser.add_argument("--sample-metadata", dest="sample_metadata", default=None,
                    help="TSV file with sample metadata (columns: sample, group, color).\n"
                         "If not provided, groups are auto-inferred from sample names:\n"
                         "  - Names containing 'Pre' -> 'control' group\n"
                         "  - Names containing 'Post' -> 'treatment' group\n"
                         "  - Otherwise -> sample name becomes group name")
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
parser.add_argument("--min-cluster-size", dest="min_cluster_size", type=int, default=10,
                    help="Minimum cluster size to consider (default: 10)")
parser.add_argument("--min-read-length", dest="min_read_length", type=int, default=0,
                    help="Minimum read length in bp to include (default: 0, no filter)")
parser.add_argument("--n-representatives", dest="n_reps", type=int, default=5,
                    help="Number of representative reads per cluster (default: 5)")
parser.add_argument("--small-feature-quantile", dest="small_quantile", type=float, default=0.1,
                    help="Quantile threshold for small feature collapsing (default: 0.1)")
parser.add_argument("--linkage-method", dest="linkage_method", default="ward",
                    help="Linkage method for hierarchical clustering (default: ward)")
parser.add_argument("--matrix-type", dest="matrix_type", default="binary",
                    choices=["binary", "count", "length_weighted"],
                    help="Type of adjacency matrix:\n"
                         "  binary: 0/1 for presence/absence of transitions\n"
                         "  count: count of each transition\n"
                         "  length_weighted: transitions weighted by feature length (default: binary)")
parser.add_argument("--include-abundance", dest="include_abundance", action="store_true",
                    help="Include feature abundance (proportion of read) as additional dimensions")
parser.add_argument("--plot-umap", dest="plot_umap", action="store_true",
                    help="Generate UMAP visualization (requires umap-learn)")
parser.add_argument("--umap-neighbors", dest="umap_neighbors", type=int, default=25,
                    help="UMAP n_neighbors parameter (default: 25)")
parser.add_argument("--umap-min-dist", dest="umap_min_dist", type=float, default=0.2,
                    help="UMAP min_dist parameter (default: 0.2)")

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

def get_edges(features):
    """Get directed edges from a list of features."""
    if len(features) <= 1:
        return []
    edges = []
    for i in range(len(features) - 1):
        edges.append((features[i], features[i + 1]))
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
            enrichment = "Mixed"
    else:
        enrichment = "Mixed"

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
        enrichment = "Mixed"

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
        enrichment = "Mixed"

    return {
        'group_counts': {s: sample_counts.get(s, 0) for s in samples},
        'group_pcts': sample_pcts,
        'odds_ratio': np.nan,
        'p_value': min_p,
        'enrichment': enrichment,
        'dominant_group': max(samples, key=lambda s: sample_pcts[s])
    }

# --- Load and process data ---
print("=" * 60)
print("KaryoScope Cluster Analysis")
print("=" * 60)

print(f"\nLoading BED file(s)...")
dfs = []
sample_labels = []
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

# Generate colors for groups
group_colors = generate_group_colors(all_groups, sample_colors_from_meta)

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
print(f"\n--- Building adjacency matrix ({args.matrix_type}) ---")

# Get feature data per read (with lengths for weighting)
read_feature_data = collapsed_df.groupby('read').apply(
    lambda x: list(zip(x['feature'], x['length']))
).to_dict()

read_features = {r: [f[0] for f in data] for r, data in read_feature_data.items()}
read_feature_lengths = {r: {f: l for f, l in data} for r, data in read_feature_data.items()}

# Get edges with weights (length of source feature)
def get_weighted_edges(features_with_lengths):
    """Get directed edges with weights (length of source feature)."""
    if len(features_with_lengths) <= 1:
        return []
    edges = []
    for i in range(len(features_with_lengths) - 1):
        from_feat, from_len = features_with_lengths[i]
        to_feat, _ = features_with_lengths[i + 1]
        edges.append((from_feat, to_feat, from_len))
    return edges

read_edges = {}
read_weighted_edges = {}
for read_name, features in read_features.items():
    read_edges[read_name] = get_edges(features)
    read_weighted_edges[read_name] = get_weighted_edges(read_feature_data[read_name])

all_features = sorted(collapsed_df['feature'].unique())
all_pairs = []
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

# Fast enrichment check for k-optimization (simplified - just needs direction)
def fast_enrichment_check(cluster_mask, read_groups, group_totals, control_grp, all_grps):
    """Fast enrichment calculation for k-optimization loop."""
    cluster_groups = read_groups[cluster_mask]
    group_counts = Counter(cluster_groups)
    total_in = len(cluster_groups)
    if total_in == 0:
        return "Mixed"

    # For two groups, use simple proportion comparison
    if len(all_grps) == 2:
        ctrl_count = group_counts.get(control_grp, 0)
        ctrl_pct = ctrl_count / total_in
        other_grp = [g for g in all_grps if g != control_grp][0]
        other_count = group_counts.get(other_grp, 0)

        # Quick significance check using binomial proportion
        expected_ctrl = group_totals[control_grp] / sum(group_totals.values())
        if ctrl_pct > expected_ctrl + 0.15:  # >15% enrichment threshold for speed
            return f"{control_grp}-enriched"
        elif ctrl_pct < expected_ctrl - 0.15:
            return f"{other_grp}-enriched"
    return "Mixed"

# Calculate within-cluster variance for different k values
if args.n_clusters is None:
    # Try different numbers of clusters and evaluate
    k_range = range(5, min(51, len(read_names) // 10))

    print(f"Testing cluster counts from {min(k_range)} to {max(k_range)}...")

    best_k = 10
    best_score = -np.inf

    cluster_stats = []
    for k in k_range:
        labels = fcluster(linkage_matrix, k, criterion='maxclust')
        labels_arr = np.array(labels)

        # Count clusters meeting minimum size
        cluster_sizes = Counter(labels)
        valid_clusters = sum(1 for size in cluster_sizes.values() if size >= args.min_cluster_size)

        # Calculate enrichment diversity (using fast method)
        enrichments = []
        for cluster_id in range(1, k + 1):
            cluster_mask = labels_arr == cluster_id
            if cluster_mask.sum() >= args.min_cluster_size:
                enrich = fast_enrichment_check(cluster_mask, read_groups_arr, group_totals, control_group, all_groups)
                enrichments.append(enrich)

        # Score: maximize variety of enrichments while having reasonable cluster count
        group_enriched_counts = {g: sum(1 for e in enrichments if e == f"{g}-enriched") for g in all_groups}
        mixed = enrichments.count("Mixed")

        # We want some of each type - score based on having at least some enriched clusters per group
        score = sum(min(c, 2) for c in group_enriched_counts.values()) + min(mixed, 3) + valid_clusters * 0.1

        cluster_stats.append({
            'k': k,
            'valid_clusters': valid_clusters,
            **{f'{g}_enriched': group_enriched_counts.get(g, 0) for g in all_groups},
            'mixed': mixed,
            'score': score
        })

        if score > best_score:
            best_score = score
            best_k = k

    # Save cluster analysis
    stats_df = pd.DataFrame(cluster_stats)
    stats_file = f"{args.output_prefix}.cluster_k_analysis.tsv"
    stats_df.to_csv(stats_file, sep='\t', index=False)
    print(f"  Saved k analysis to: {stats_file}")

    n_clusters = best_k
    print(f"\n  Optimal k: {n_clusters} (score: {best_score:.2f})")
else:
    n_clusters = args.n_clusters
    print(f"  Using specified k: {n_clusters}")

# --- Cut tree and analyze clusters ---
print(f"\n--- Cutting tree into {n_clusters} clusters ---")
cluster_labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')

# Analyze each cluster
cluster_analysis = []
cluster_reads_dict = {}

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

    # Get representative reads - balanced by group but prioritizing centroid proximity
    # Pre-compute cluster read groups for efficiency
    cluster_read_groups = [read_to_group[r] for r in cluster_reads]

    group_reads = {}
    for g in all_groups:
        # Build list of (index, read, distance) for this group
        group_reads[g] = [(i, cluster_reads[i], distances_to_centroid[i])
                         for i in range(len(cluster_reads))
                         if cluster_read_groups[i] == g]
        # Sort by distance to centroid
        group_reads[g].sort(key=lambda x: x[2])

    # Calculate proportional representation for each group
    group_n_reps = {}
    total_pct = sum(stats['group_pcts'].values())
    remaining_reps = args.n_reps

    for g in all_groups:
        pct = stats['group_pcts'].get(g, 0)
        n_reps = round(args.n_reps * pct / 100) if total_pct > 0 else 0
        n_reps = min(n_reps, len(group_reads[g]), remaining_reps)
        group_n_reps[g] = n_reps
        remaining_reps -= n_reps

    # Fill remaining slots from groups with available reads
    while remaining_reps > 0:
        filled = False
        for g in sorted(all_groups, key=lambda x: stats['group_pcts'].get(x, 0), reverse=True):
            if len(group_reads[g]) > group_n_reps[g] and remaining_reps > 0:
                group_n_reps[g] += 1
                remaining_reps -= 1
                filled = True
        if not filled:
            break

    # Select best (closest to centroid) from each group
    representative_reads = []
    avg_dists = {}
    for g in all_groups:
        n_reps = group_n_reps[g]
        selected = group_reads[g][:n_reps]
        representative_reads.extend([r[1] for r in selected])
        avg_dists[g] = np.mean([r[2] for r in selected]) if n_reps > 0 else np.nan

    # Build cluster analysis record
    cluster_record = {
        'cluster_id': cluster_id,
        'size': len(cluster_reads),
        'odds_ratio': stats['odds_ratio'],
        'p_value': stats['p_value'],
        'enrichment': stats['enrichment'],
        'centroid_read': centroid_read,
        'centroid_sample': read_to_sample[centroid_read],
        'centroid_group': sample_to_group.get(read_to_sample[centroid_read], read_to_sample[centroid_read]),
        'representative_reads': ','.join(representative_reads)
    }

    # Add group-specific columns
    for g in all_groups:
        cluster_record[f'{g}_count'] = stats['group_counts'].get(g, 0)
        cluster_record[f'{g}_pct'] = stats['group_pcts'].get(g, 0)
        cluster_record[f'avg_dist_{g}'] = avg_dists.get(g, np.nan)

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
cluster_df.to_csv(analysis_file, sep='\t', index=False)
print(f"\n  Saved cluster analysis to: {analysis_file}")

# Save cluster assignments for all reads
assignments = pd.DataFrame({
    'read': read_names,
    'cluster': cluster_labels,
    'sample': [read_to_sample[r] for r in read_names]
})
assignments_file = f"{args.output_prefix}.cluster_assignments.tsv"
assignments.to_csv(assignments_file, sep='\t', index=False)
print(f"  Saved cluster assignments to: {assignments_file}")

# Save representative reads per cluster (for plotting)
reps_file = f"{args.output_prefix}.representative_reads.tsv"
with open(reps_file, 'w') as f:
    f.write("cluster_id\tenrichment\tread\tsample\n")
    for _, row in cluster_df.iterrows():
        for read in row['representative_reads'].split(','):
            sample = read_to_sample[read]
            f.write(f"{row['cluster_id']}\t{row['enrichment']}\t{read}\t{sample}\n")
print(f"  Saved representative reads to: {reps_file}")

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

# 1. Dendrogram with cluster coloring
ax1 = axes[0, 0]
from scipy.cluster.hierarchy import set_link_color_palette
# Color palette for clusters
colors = plt.cm.tab20(np.linspace(0, 1, n_clusters))
set_link_color_palette([matplotlib.colors.rgb2hex(c) for c in colors])

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
from matplotlib.patches import Patch
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

        # Create dual UMAP plot (sample + enrichment coloring)
        fig, axes = plt.subplots(1, 2, figsize=(16, 7))

        # Left: Colored by sample (use sample_colors generated earlier)
        sample_list = [read_to_sample[r] for r in read_names]
        point_colors_sample = [sample_colors.get(s, '#999999') for s in sample_list]

        axes[0].scatter(embedding[:, 0], embedding[:, 1], c=point_colors_sample, s=15, alpha=0.6)
        axes[0].set_title("UMAP - Colored by Sample")
        axes[0].set_xlabel("UMAP 1")
        axes[0].set_ylabel("UMAP 2")

        # Add sample legend
        sample_patches = [Patch(facecolor=sample_colors[s], label=s) for s in sorted(sample_colors.keys())]
        axes[0].legend(handles=sample_patches, loc='upper right')

        # Right: Colored by enrichment (use group_colors for dynamic mapping)
        enrichment_colors = {f'{g}-enriched': c for g, c in group_colors.items()}
        enrichment_colors['Mixed'] = '#CCCCCC'
        point_colors_enrich = []
        for r in read_names:
            cluster = read_to_cluster[r]
            enrich = cluster_to_enrichment.get(cluster, 'Mixed')
            point_colors_enrich.append(enrichment_colors.get(enrich, '#CCCCCC'))

        axes[1].scatter(embedding[:, 0], embedding[:, 1], c=point_colors_enrich, s=15, alpha=0.6)
        axes[1].set_title("UMAP - Colored by Cluster Enrichment")
        axes[1].set_xlabel("UMAP 1")
        axes[1].set_ylabel("UMAP 2")

        # Add enrichment legend
        enrich_patches = [Patch(facecolor=c, label=e) for e, c in enrichment_colors.items()]
        axes[1].legend(handles=enrich_patches, loc='upper right')

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

            fig = go.Figure()

            # Add traces by sample
            for sample in sorted(set(sample_list)):
                mask = [s == sample for s in sample_list]
                sample_embedding = embedding[mask]
                sample_reads = [r for r, m in zip(read_names, mask) if m]
                sample_clusters = [read_to_cluster[r] for r in sample_reads]
                sample_enrichments = [cluster_to_enrichment.get(c, 'Mixed') for c in sample_clusters]

                hover_text = [f"Read: {r[:12]}...<br>Cluster: {c}<br>Enrichment: {e}"
                             for r, c, e in zip(sample_reads, sample_clusters, sample_enrichments)]

                fig.add_trace(go.Scatter(
                    x=sample_embedding[:, 0],
                    y=sample_embedding[:, 1],
                    mode='markers',
                    name=sample,
                    text=hover_text,
                    hoverinfo='text',
                    marker=dict(
                        size=6,
                        color=sample_colors.get(sample, '#999999'),
                        opacity=0.7
                    )
                ))

            fig.update_layout(
                title=f"UMAP - Interactive ({len(read_names):,} reads)",
                xaxis_title="UMAP 1",
                yaxis_title="UMAP 2",
                hovermode='closest',
                width=900,
                height=700
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
print(f"  - {reps_file}")
print(f"  - {matrix_file}")
print(f"  - {plot_file}")
print(f"\nNext step: Use representative_reads.tsv with KaryoScope_cluster_plot.py")
print(f"to visualize representative reads from each cluster.")
print(f"Use --feature-matrix with the .feature_matrix.npz file for dendrogram header.")
