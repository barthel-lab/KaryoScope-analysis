# KaryoScope Cluster Analysis
# Analyzes hierarchical clustering results to identify biologically interesting clusters
# with sequence assignments sorted by centroid distance for visualization
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
#   per-sample: Each sample vs all others (no groups required)

import argparse
import sys
import gzip
import fnmatch
import numpy as np
import pandas as pd

# Capture original command line for logging
_original_command = ' '.join(sys.argv)
from collections import defaultdict, Counter
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist, squareform
from scipy.stats import fisher_exact, false_discovery_control
from sklearn.metrics import silhouette_score, calinski_harabasz_score
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
matplotlib.use('Agg')

# Keep text editable in output files (not converted to paths)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

# --- Argument parsing ---
parser = argparse.ArgumentParser(
    description="Analyze KaryoScope clustering to identify sample-enriched clusters.",
    formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument("--bed", required=True, nargs='+',
                    help="Path to input BED file(s) (can be gzipped). Multiple files will be concatenated.\n"
                         "To merge multiple featuresets, use KaryoScope_merge_beds.py first.")
parser.add_argument("--output-prefix", dest="output_prefix", required=True,
                    help="Prefix for output files")
parser.add_argument("--sample-metadata", dest="sample_metadata", default=None,
                    help="TSV file with sample metadata (columns: sample, group, color).\n"
                         "If not provided, each sample becomes its own group.\n"
                         "Use --control-group to specify the reference group.")
parser.add_argument("--comparison-mode", dest="comparison_mode", default="two-group",
                    choices=["two-group", "per-sample"],
                    help="Comparison mode for enrichment testing:\n"
                         "  two-group: Fisher's exact test between control and treatment\n"
                         "  per-sample: Each sample vs all others (default: two-group)")
parser.add_argument("--control-group", dest="control_group", default=None,
                    help="Name of control group for two-group comparison (default: auto-detect)")
parser.add_argument("--n-clusters", dest="n_clusters", type=int, default=None,
                    help="Number of clusters to cut tree into (default: auto-determine)")
parser.add_argument("--min-k", dest="min_k", type=int, default=40,
                    help="Minimum number of clusters to test during auto-detection (default: 40)")
parser.add_argument("--max-k", dest="max_k", type=int, default=300,
                    help="Maximum number of clusters to test during auto-detection (default: 300)")
parser.add_argument("--k-selection", dest="k_selection", default="composite-knee",
                    choices=["composite", "silhouette", "calinski", "composite-knee"],
                    help="Metric for selecting optimal k:\n"
                         "  composite: max weighted combination of silhouette + enrichment\n"
                         "  silhouette: cluster cohesion (favors fewer, tighter clusters)\n"
                         "  calinski: Calinski-Harabasz index\n"
                         "  composite-knee: knee/elbow of composite score curve (default, diminishing returns)")
parser.add_argument("--min-cluster-size", dest="min_cluster_size", type=int, default=3,
                    help="Minimum cluster size to consider (default: 3)")
parser.add_argument("--min-sequence-length", dest="min_sequence_length", type=int, default=10000,
                    help="Minimum sequence length in bp to include (default: 10000)")
parser.add_argument("--max-sequence-length", dest="max_sequence_length", type=int, default=50000,
                    help="Maximum sequence length in bp to include (default: 50000)")
parser.add_argument("--sequence-list", dest="sequence_list", default=None,
                    help="File with sequence names to include (one per line). Filters BED to only these sequences.")
parser.add_argument("--exclude-features", dest="exclude_features", default="novel,canonical_telomere*",
                    help="Comma-separated list of features to exclude, supports wildcards (* and ?) (default: 'novel,canonical_telomere*')")
parser.add_argument("--linkage-method", dest="linkage_method", default="ward",
                    help="Linkage method for hierarchical clustering (default: ward)")
parser.add_argument("--matrix-type", dest="matrix_type", default="count_log1p_zscore_blockweight",
                    choices=["binary", "count", "length_weighted", "count_log1p",
                             "count_log1p_zscore", "count_log1p_zscore_blockweight"],
                    help="Type of adjacency matrix (controls both edges and abundance):\n"
                         "  binary: 0/1 presence/absence\n"
                         "  count: raw transition counts / raw feature bp totals\n"
                         "  length_weighted: feature-length proportions\n"
                         "  count_log1p: log(count+1), compresses dynamic range for rare features\n"
                         "  count_log1p_zscore: log(count+1) with per-column z-score normalization\n"
                         "  count_log1p_zscore_blockweight: log(count+1) with z-score + block reweighting for equal edge/abundance contribution (default)")
parser.add_argument("--edges", dest="edge_mode", default="symmetric",
                    choices=["directional", "symmetric"],
                    help="Edge counting mode:\n"
                         "  directional: standard A->B edge counting\n"
                         "  symmetric: edges are sorted alphabetically, A->B and B->A both count as A->B (default: symmetric)")
parser.add_argument("--matrix-mode", dest="matrix_mode", default="combined",
                    choices=["layered", "combined"],
                    help="Matrix building mode for merged featuresets:\n"
                         "  layered: split colon-separated features into layers, build matrices per layer\n"
                         "  combined: treat merged features as atomic (default, higher dimensionality)")
parser.add_argument("--include-edges", dest="include_edges",
                    action=argparse.BooleanOptionalAction, default=True,
                    help="Include edge (transition) dimensions (default: True)")
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
parser.add_argument("--umap-html", dest="umap_html",
                    action=argparse.BooleanOptionalAction, default=False,
                    help="Generate interactive HTML UMAP with Plotly (default: False, requires plotly)")
parser.add_argument("--perfect-threshold", dest="perfect_threshold", type=float, default=1.0,
                    help="Threshold for perfect enrichment (default: 1.0 = 100%%)")
parser.add_argument("--strong-threshold", dest="strong_threshold", type=float, default=0.80,
                    help="Threshold for strong enrichment (default: 0.80 = 80%%)")
parser.add_argument("--early-stopping", dest="early_stopping", type=int, default=150,
                    help="Stop k search if no improvement for N iterations (0 to disable) (default: 150)")
parser.add_argument("--silhouette-sample-size", dest="silhouette_sample_size", type=int, default=2000,
                    help="Sample size for silhouette score calculation (default: 2000)")
parser.add_argument("--reduce-dims", dest="reduce_dims", type=int, default=500,
                    help="Reduce matrix to N dimensions using truncated SVD before clustering.\n"
                         "Recommended for merged BED files which can create very high-dimensional matrices.\n"
                         "Set to 0 to disable reduction (default: 500)")
parser.add_argument("--background", dest="background", default="white",
                    choices=["white", "black", "both"],
                    help="Background color for diagnostic plots:\n"
                         "  white: white background (default)\n"
                         "  black: dark background\n"
                         "  both: generate both versions (dark files have _dark suffix)")
parser.add_argument("--log-file", dest="log_file",
                    action=argparse.BooleanOptionalAction, default=True,
                    help="Save console output to {output_prefix}.log (default: True)")
parser.add_argument("--fdr-threshold", dest="fdr_threshold", type=float, default=0.05,
                    help="FDR q-value threshold for calling enrichment (default: 0.05)")
parser.add_argument("--fdr-method", dest="fdr_method", default="bh",
                    choices=["bh", "by"],
                    help="FDR correction method:\n"
                         "  bh: Benjamini-Hochberg (default, assumes independence or PRDS)\n"
                         "  by: Benjamini-Yekutieli (more conservative, valid for any dependency)")
parser.add_argument("--analysis-mode", dest="analysis_mode", default="enrichment",
                    choices=["enrichment", "structure"],
                    help="Analysis mode:\n"
                         "  enrichment: Cluster all sequences and test for group enrichment (original mode)\n"
                         "  structure: Cluster per-chromosome to identify structural outliers (default: enrichment)")
parser.add_argument("--structural-threshold", "--st", dest="structural_threshold", type=float, default=0.25,
                    help="Distance threshold for structural outlier clustering (default: 0.25)")

args = parser.parse_args()

# --- Set up logging ---
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

if not args.include_edges and not args.include_abundance:
    parser.error("At least one of --edges or --abundance must be enabled")

if args.log_file:
    log_path = f"{args.output_prefix}.log"
    sys.stdout = TeeLogger(log_path)


# =============================================================================
# Structural Analysis Mode (Additive)
# =============================================================================

def run_structure_mode(in_data, args):
    """Run the structure-based analysis workflow (per-chromosome clustering)."""
    print("=" * 60)
    print("Mode: Structural Analysis (Per-Chromosome)")
    print("=" * 60)
    
    # Identify chromosomes
    chromosomes = sorted(in_data['chromosome'].unique())
    chromosomes = [c for c in chromosomes if c != 'unknown'] # Filter unknown
    
    print(f"Found {len(chromosomes)} chromosomes: {', '.join(chromosomes)}")
    
    all_cluster_assignments = []
    chromosome_linkages = {}
    
    # Standard hierarchical clustering imports
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    
    # Simple pure-python Levenshtein implementation
    def levenshtein_distance(s1, s2):
        if len(s1) < len(s2):
            return levenshtein_distance(s2, s1)
        if len(s2) == 0:
            return len(s1)
        previous_row = range(len(s2) + 1)
        for i, c1 in enumerate(s1):
            current_row = [i + 1]
            for j, c2 in enumerate(s2):
                insertions = previous_row[j + 1] + 1
                deletions = current_row[j] + 1
                substitutions = previous_row[j] + (c1 != c2)
                current_row.append(min(insertions, deletions, substitutions))
            previous_row = current_row
        return previous_row[-1]

    for chrom in chromosomes:
        print(f"\n--- Processing {chrom} ---")
        chrom_data = in_data[in_data['chromosome'] == chrom].copy()
        
        # Pre-calculate sequence -> sample map for speed
        read_sample_map = dict(zip(chrom_data['sequence'], chrom_data['sample']))

        # Create "sequence strings" for distance calculation
        seq_strings = {}
        for read, group in chrom_data.groupby('sequence'):
            # Sort by coordinate
            sorted_group = group.sort_values('start')
            # Create feature string
            feat_list = sorted_group['feature'].tolist()
            seq_strings[read] = feat_list
            
        seq_names = sorted(seq_strings.keys())
        n_reads = len(seq_names)
        
        if n_reads < 2:
            print(f"  Skipping {chrom} (only {n_reads} sample)")
            for r in seq_names:
                all_cluster_assignments.append({
                    'read': r, 'chromosome': chrom, 
                    'cluster_id': f"{chrom}_Major", 'cluster_type': 'Major',
                    'enrichment': 'Major', 
                    'sample': read_sample_map[r]
                })
            continue

        # Optimization: Deduplicate sequences based on feature content
        unique_feats = sorted(list(set(f for r in seq_names for f in seq_strings[r])))
        feat_map = {f: chr(i+200) for i, f in enumerate(unique_feats)} 
        
        unique_structures = {} 
        for r in seq_names:
            encoded = "".join([feat_map[f] for f in seq_strings[r]])
            if encoded not in unique_structures:
                unique_structures[encoded] = []
            unique_structures[encoded].append(r)
            
        unique_structure_list = sorted(unique_structures.keys())
        n_unique = len(unique_structure_list)
        print(f"  {n_reads} sequences -> {n_unique} unique structures")
        
        if n_unique > 1:
            print(f"  Computing distances for {n_unique} unique structures...")
            dist_matrix = np.zeros((n_unique, n_unique))
            raw_dist_matrix = np.zeros((n_unique, n_unique))
            for i in range(n_unique):
                for j in range(i + 1, n_unique):
                    raw_d = levenshtein_distance(unique_structure_list[i], unique_structure_list[j])
                    max_len = max(len(unique_structure_list[i]), len(unique_structure_list[j]))
                    raw_dist_matrix[i, j] = raw_dist_matrix[j, i] = raw_d
                    if max_len > 0:
                        dist_matrix[i, j] = dist_matrix[j, i] = raw_d / max_len

            # Clustering (still using normalized distance for Ward linkage)
            Z = linkage(squareform(dist_matrix), method='ward')
            chromosome_linkages[chrom] = Z
            
            # Use user-defined threshold
            threshold = args.structural_threshold
            clusters = fcluster(Z, t=threshold, criterion='distance')
            
            # Report diagnostics
            d_vals = dist_matrix[np.triu_indices(n_unique, k=1)]
            print(f"  Distance Stats (n={len(d_vals)}): min={np.min(d_vals):.3f}, med={np.median(d_vals):.3f}, max={np.max(d_vals):.3f}")
            print(f"  Clustering with threshold={threshold} -> {len(set(clusters))} clusters")
        else:
            clusters = [1] * n_unique
            chromosome_linkages[chrom] = None

        # Determine Major Cluster
        cluster_read_counts = defaultdict(int)
        structure_to_cluster = {s: c for s, c in zip(unique_structure_list, clusters)}
        
        # Track which structure is most abundant in the major cluster
        major_cluster_struct_counts = defaultdict(int)
        
        for s, reads in unique_structures.items():
            cid = structure_to_cluster[s]
            cluster_read_counts[cid] += len(reads)
            
        major_cluster_id = max(cluster_read_counts, key=cluster_read_counts.get)
        
        # Find the most abundant structure in the Major cluster to use as consensus
        major_structs = []
        for s, reads in unique_structures.items():
            if structure_to_cluster[s] == major_cluster_id:
                major_structs.append((s, len(reads)))
        
        # Sort by count desc
        major_structs.sort(key=lambda x: (-x[1], x[0]))
        major_consensus_s = major_structs[0][0]
        major_consensus_idx = unique_structure_list.index(major_consensus_s)
        
        # Pre-calculate lengths per feature for the length-weighted metric
        read_feature_lengths = chrom_data.groupby(['read', 'feature'])['length'].sum().unstack(fill_value=0)
        consensus_read = unique_structures[major_consensus_s][0]
        consensus_lengths = read_feature_lengths.loc[consensus_read]
        consensus_feat_set = set(major_consensus_s)

        # Assign with divergence score
        for s, reads in unique_structures.items():
            cid = structure_to_cluster[s]
            is_major = (cid == major_cluster_id)
            c_type = "Major" if is_major else "Outlier"
            full_cluster_id = f"{chrom}_{c_type}"
            if not is_major:
                full_cluster_id += f"_{cid}"
            
            # Relative to the major consensus structure
            s_idx = unique_structure_list.index(s)
            
            if n_unique > 1:
                norm_div = dist_matrix[s_idx, major_consensus_idx]
                raw_div = raw_dist_matrix[s_idx, major_consensus_idx]
            else:
                norm_div = 0.0
                raw_div = 0.0
            
            # Binary Metric: Symmetric difference of feature types present
            read_feat_set = set(s)
            binary_div = len(consensus_feat_set ^ read_feat_set)
                
            for r in reads:
                # Length-Weighted Metric: Sum of absolute differences in bp per feature type
                r_lengths = read_feature_lengths.loc[r]
                # Ensure we handle features present in either the read or the consensus
                all_feats = sorted(list(set(consensus_lengths.index) | set(r_lengths.index)))
                l_div = 0.0
                for f in all_feats:
                    v_cons = consensus_lengths.get(f, 0.0)
                    v_read = r_lengths.get(f, 0.0)
                    l_div += abs(v_cons - v_read)

                all_cluster_assignments.append({
                    'read': r, 'chromosome': chrom, 
                    'cluster': full_cluster_id, 'cluster_type': c_type,
                    'norm_divergence': norm_div,
                    'raw_divergence': raw_div,
                    'binary_divergence': float(binary_div),
                    'length_weighted_divergence': float(l_div),
                    'enrichment': c_type, 'sample': read_sample_map[r]
                })

    # Save
    matrix_file = f"{args.output_prefix}.feature_matrix.npz"
    assignments_file = f"{args.output_prefix}.sequence_assignments.tsv"
    assignment_df = pd.DataFrame(all_cluster_assignments)
    if 'group' not in assignment_df.columns:
        assignment_df['group'] = assignment_df['sample']
    assignment_df.to_csv(assignments_file, sep='\t', index=False)
    
    np.savez_compressed(
        matrix_file, mode='structure', chromosomes=chromosomes,
        linkages=chromosome_linkages, info="Per-chromosome structural analysis"
    )
    print(f"\nSaved structural assignments: {assignments_file}")
    sys.exit(0)

# --- Plot style management ---
def get_plot_style(bg_mode):
    """Return style configuration for given background mode."""
    if bg_mode == 'black':
        return {
            'bg_color': '#1a1a1a',
            'text_color': 'white',
            'grid_color': '#404040',
            'annotation_bg': '#333333',
            'annotation_edge': 'white',
            'edge_color': 'white',
            'style': 'dark_background'
        }
    else:  # white
        return {
            'bg_color': 'white',
            'text_color': 'black',
            'grid_color': '#cccccc',
            'annotation_bg': 'white',
            'annotation_edge': 'gray',
            'edge_color': 'black',
            'style': 'default'
        }

def apply_plot_style(bg_mode):
    """Apply matplotlib style for given background mode."""
    style = get_plot_style(bg_mode)
    if bg_mode == 'black':
        plt.style.use('dark_background')
    else:
        plt.style.use('default')
    # Re-apply fonttype settings (style.use resets them)
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['svg.fonttype'] = 'none'
    return style

def get_backgrounds_to_generate():
    """Return list of (bg_mode, suffix) tuples based on --background argument."""
    if args.background == 'both':
        return [('white', ''), ('black', '_dark')]
    elif args.background == 'black':
        return [('black', '')]
    else:
        return [('white', '')]

# Initialize with first background mode (for any early plots)
_bg_modes = get_backgrounds_to_generate()
_current_style = apply_plot_style(_bg_modes[0][0])
PLOT_BG_COLOR = _current_style['bg_color']
PLOT_TEXT_COLOR = _current_style['text_color']
PLOT_GRID_COLOR = _current_style['grid_color']

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
                # Try to extract chromosome from 5th column if available (standard KaryoScope output)
                chrom = parts[4] if len(parts) >= 5 else "unknown"
                
                records.append({
                    'sequence': read,
                    'start': start,
                    'end': end,
                    'feature': feature,
                    'chromosome': chrom,
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
        elif edge_mode == "symmetric":
            # Sort alphabetically so A->B and B->A both become the same edge
            sorted_pair = tuple(sorted([from_feat, to_feat]))
            edges.append(sorted_pair)
    return edges


def detect_feature_layers(features):
    """Detect if features are merged (colon-separated) and count layers.

    Args:
        features: List of feature names

    Returns:
        n_layers: Number of layers (1 if not merged, >1 if merged)
    """
    # Sample some features to detect layer count
    sample_features = list(features)[:100]
    layer_counts = [len(f.split(':')) for f in sample_features]

    # Use mode of layer counts (most features should have same structure)
    if not layer_counts:
        return 1
    n_layers = max(set(layer_counts), key=layer_counts.count)
    return n_layers


def split_feature_into_layers(feature, n_layers):
    """Split a merged feature into its component layers.

    Args:
        feature: Feature string (e.g., "chr7:p_arm:nonsubtelomeric")
        n_layers: Expected number of layers

    Returns:
        List of layer components (e.g., ["chr7", "p_arm", "nonsubtelomeric"])
    """
    parts = feature.split(':')
    if len(parts) == n_layers:
        return parts
    elif len(parts) < n_layers:
        # Pad with 'unknown' if fewer parts than expected
        return parts + ['unknown'] * (n_layers - len(parts))
    else:
        # If more parts than expected, join extra parts into last layer
        return parts[:n_layers-1] + [':'.join(parts[n_layers-1:])]


def get_layer_features(features_with_lengths, layer_idx, n_layers):
    """Extract features for a specific layer from merged feature data.

    Args:
        features_with_lengths: List of (feature, length) tuples
        layer_idx: Which layer to extract (0-indexed)
        n_layers: Total number of layers

    Returns:
        List of (layer_feature, length) tuples
    """
    result = []
    for feat, length in features_with_lengths:
        parts = split_feature_into_layers(feat, n_layers)
        result.append((parts[layer_idx], length))
    return result

def load_sample_metadata(metadata_file, sample_labels):
    """Load sample metadata from TSV file or auto-generate defaults.

    When a metadata file is provided:
      - Uses group names exactly as specified in the file
      - Colors are optional (sample colors in 'color' column)
      - Group colors can be specified in 'group_color' column

    When no metadata file:
      - Each sample becomes its own group (named after the sample)
      - This allows --control-group to specify the reference

    Returns: (sample_to_group, sample_to_color, group_to_color)
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
        group_to_color = {}

        for _, row in meta_df.iterrows():
            sample = row['sample']
            # Use group if provided, otherwise use sample name as group
            group = row.get('group', sample) if 'group' in meta_df.columns else sample
            sample_to_group[sample] = group
            if 'color' in meta_df.columns and pd.notna(row.get('color')):
                sample_to_color[sample] = row['color']
            # Parse explicit group colors if provided
            if 'group_color' in meta_df.columns and pd.notna(row.get('group_color')):
                group_to_color[group] = row['group_color']

        # Check all samples are covered
        missing = set(sample_labels) - set(sample_to_group.keys())
        if missing:
            print(f"  Warning: Samples not in metadata file: {missing}")
            # Auto-assign missing samples to their own groups
            for s in missing:
                sample_to_group[s] = s

        return sample_to_group, sample_to_color, group_to_color

    else:
        # No metadata file: each sample is its own group
        # This is clean and doesn't rely on text matching
        sample_to_group = {sample: sample for sample in sample_labels}
        return sample_to_group, {}, {}


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

    # Handle edge case: single group or empty groups
    if len(groups) < 2:
        return {
            'group_counts': {g: len(cluster_samples) for g in groups} if groups else {},
            'group_pcts': {g: 100.0 for g in groups} if groups else {},
            'odds_ratio': np.nan,
            'p_value': 1.0,
            'enrichment': 'mixed',
            'dominant_group': groups[0] if groups else None
        }

    non_control_groups = [g for g in groups if g != control_group]
    treatment_group = non_control_groups[0] if len(groups) == 2 and non_control_groups else None

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

    # Determine enrichment direction using odds ratio (accounts for unequal group sizes)
    # odds_ratio > 1 means control is over-represented relative to baseline
    # odds_ratio < 1 means treatment is over-represented relative to baseline
    if p_value < 0.05:
        if odds_ratio > 1:
            enrichment = f"{control_group}-enriched"
        elif odds_ratio < 1:
            enrichment = f"{treatment_group}-enriched"
        else:
            enrichment = "mixed"
    else:
        enrichment = "mixed"

    # Build group counts dict
    group_counts = {control_group: control_in, treatment_group: treatment_in}
    group_pcts = {control_group: control_pct, treatment_group: treatment_pct}

    # Dominant group based on odds ratio (which group is over-represented)
    dominant_group = control_group if odds_ratio > 1 else treatment_group

    return {
        'group_counts': group_counts,
        'group_pcts': group_pcts,
        'odds_ratio': odds_ratio,
        'p_value': p_value,
        'enrichment': enrichment,
        'dominant_group': dominant_group
    }


def calculate_enrichment_per_sample(cluster_samples, sample_totals):
    """Calculate enrichment per sample (each sample vs all others)."""
    samples = list(sample_totals.keys())

    # Handle edge case: no samples
    if not samples:
        return {
            'group_counts': {},
            'group_pcts': {},
            'odds_ratio': np.nan,
            'p_value': 1.0,
            'enrichment': 'mixed',
            'dominant_group': None,
            'all_p_values': {},
            'all_odds_ratios': {}
        }

    # Count samples in cluster
    sample_counts = Counter(cluster_samples)

    # Calculate percentages
    total_in = len(cluster_samples)
    sample_pcts = {s: (sample_counts.get(s, 0) / total_in * 100) if total_in > 0 else 0 for s in samples}

    # Fisher's exact for each sample vs rest
    p_values = {}
    odds_ratios = {}
    for sample in samples:
        in_cluster = sample_counts.get(sample, 0)
        out_cluster = sample_totals[sample] - in_cluster
        other_in = total_in - in_cluster
        other_out = sum(sample_totals.values()) - sample_totals[sample] - other_in

        odds, p_val = fisher_exact([[in_cluster, other_in], [out_cluster, other_out]], alternative='greater')
        p_values[sample] = p_val
        odds_ratios[sample] = odds

    # Find most significant ENRICHED sample (p < 0.05 AND odds > 1)
    # Bug fix: previously picked min p-value without checking if enriched or depleted
    enriched_samples = {s: p for s, p in p_values.items() if p < 0.05 and odds_ratios[s] > 1}

    if enriched_samples:
        # Pick the most significant enriched sample
        min_p_sample = min(enriched_samples.keys(), key=lambda s: enriched_samples[s])
        min_p = enriched_samples[min_p_sample]
        enrichment = f"{min_p_sample}-enriched"
    elif p_values:
        # No significant enrichment - find min p-value for reporting
        min_p_sample = min(p_values.keys(), key=lambda s: p_values[s])
        min_p = p_values[min_p_sample]
        enrichment = "mixed"
    else:
        min_p_sample = None
        min_p = 1.0
        enrichment = "mixed"

    # Determine dominant sample (guard against empty samples)
    dominant_sample = max(samples, key=lambda s: sample_pcts[s]) if samples else None

    return {
        'group_counts': {s: sample_counts.get(s, 0) for s in samples},
        'group_pcts': sample_pcts,
        'odds_ratio': np.nan,
        'p_value': min_p,
        'enrichment': enrichment,
        'dominant_group': dominant_sample,
        'all_p_values': p_values,
        'all_odds_ratios': odds_ratios
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
seq_to_sample = in_data.groupby('sequence')['sample'].first().to_dict()

print(f"\nTotal records: {len(in_data):,}")
print(f"Unique sequences: {in_data['sequence'].nunique():,}")

# --- Filter by sequence list (if provided) ---
if args.sequence_list:
    print(f"\n--- Filtering by sequence list ---")
    print(f"  Sequence list file: {args.sequence_list}")
    with open(args.sequence_list) as f:
        include_sequences = set(line.strip() for line in f if line.strip())
    print(f"  Sequences in list: {len(include_sequences):,}")
    sequences_before = in_data['sequence'].nunique()
    in_data = in_data[in_data['sequence'].isin(include_sequences)]
    seq_to_sample = {r: s for r, s in seq_to_sample.items() if r in include_sequences}
    sequences_after = in_data['sequence'].nunique()
    print(f"  Sequences before filter: {sequences_before:,}")
    print(f"  Sequences after filter: {sequences_after:,}")
    print(f"  Sequences from list found in BED: {sequences_after:,}")

# Calculate feature lengths
in_data['length'] = in_data['end'] - in_data['start']

# Calculate read span (full coordinate range) BEFORE filtering
# This is the actual read length including all features
seq_span_dict = in_data.groupby('sequence').apply(
    lambda x: x['end'].max() - x['start'].min()
).to_dict()

# --- Filter sequences by read span ---
# Use the full coordinate range (read span) BEFORE excluding features,
# so reads with large excluded regions aren't incorrectly dropped
if args.min_sequence_length > 0 or args.max_sequence_length is not None:
    print(f"\n--- Filtering sequences by read span ---")
    sequences_before = in_data['sequence'].nunique()

    # Apply min and max filters using read span
    valid_sequences = set(r for r, l in seq_span_dict.items()
                     if l >= args.min_sequence_length and l <= args.max_sequence_length)
    print(f"  Length range: {args.min_sequence_length:,} - {args.max_sequence_length:,} bp")

    in_data = in_data[in_data['sequence'].isin(valid_sequences)]
    seq_to_sample = {r: s for r, s in seq_to_sample.items() if r in valid_sequences}
    seq_span_dict = {r: l for r, l in seq_span_dict.items() if r in valid_sequences}
    sequences_after = in_data['sequence'].nunique()
    print(f"  Sequences before filter: {sequences_before:,}")
    print(f"  Sequences after filter: {sequences_after:,}")
    print(f"  Sequences removed: {sequences_before - sequences_after:,}")

# --- Filter excluded features AFTER read span filtering ---
if args.exclude_features:
    exclude_patterns = [f.strip() for f in args.exclude_features.split(',')]
    print(f"\n--- Filtering excluded features ---")
    print(f"  Excluding (component match): {sorted(exclude_patterns)}")
    before_count = len(in_data)
    before_sequences = in_data['sequence'].nunique()

    # Component matching: filter if any colon-separated part matches exclude patterns (supports wildcards)
    def has_excluded_component(feature):
        components = feature.split(':')
        for comp in components:
            for pattern in exclude_patterns:
                if fnmatch.fnmatch(comp, pattern):
                    return True
        return False

    mask = in_data['feature'].apply(has_excluded_component)
    in_data = in_data[~mask]

    # Remove reads that have no remaining features after exclusion
    sequences_with_features = set(in_data['sequence'].unique())
    seq_to_sample = {r: s for r, s in seq_to_sample.items() if r in sequences_with_features}
    seq_span_dict = {r: l for r, l in seq_span_dict.items() if r in sequences_with_features}

    print(f"  Records before filter: {before_count:,}")
    print(f"  Records after filter: {len(in_data):,}")
    print(f"  Sequences before filter: {before_sequences:,}")
    print(f"  Sequences after filter: {len(sequences_with_features):,}")

# Calculate annotated length per read (sum of remaining feature lengths)
# This is the total annotated sequence, not the span from start to end
seq_lengths = in_data.groupby('sequence')['length'].sum()
seq_length_dict = seq_lengths.to_dict()

# --- Branching Point ---
if args.analysis_mode == 'structure':
    run_structure_mode(in_data, args)
    # The function will sys.exit(0), so we won't reach here
    sys.exit(0)

# --- Load sample metadata ---
print(f"\n--- Loading sample metadata ---")
sample_to_group, sample_colors_from_meta, explicit_group_colors = load_sample_metadata(args.sample_metadata, sample_labels)

# Determine unique groups
all_groups = sorted(set(sample_to_group.values()))
print(f"  Groups found: {', '.join(all_groups)}")
print(f"  Comparison mode: {args.comparison_mode}")

# Map reads to groups
seq_to_group = {r: sample_to_group.get(s, s) for r, s in seq_to_sample.items()}

# Count samples by group
group_totals = Counter(seq_to_group.values())
sample_totals = Counter(seq_to_sample.values())

print(f"\n  Group counts:")
for group, count in sorted(group_totals.items()):
    print(f"    {group}: {count:,} reads")

# Determine control group for group-level comparisons
# (used in two-group mode and for group-level stats in per-sample mode)
control_group = args.control_group
if control_group is None and len(all_groups) >= 2:
    # Default: alphabetically first group is the control/reference
    control_group = all_groups[0]
if args.comparison_mode == "two-group":
    if len(all_groups) > 2:
        print(f"  Note: {len(all_groups)} groups found, using '{control_group}' as reference.")
        print(f"        Use --control-group to specify a different reference group.")
    print(f"  Reference group: {control_group}")

# Build group colors: explicit > derived from first sample > auto-generated
# Start with colors derived from first sample per group
group_colors_from_meta = {}
for sample, color in sample_colors_from_meta.items():
    group = sample_to_group.get(sample, sample)
    if group not in group_colors_from_meta:
        group_colors_from_meta[group] = color
# Override with explicit group colors if provided
group_colors_from_meta.update(explicit_group_colors)
if explicit_group_colors:
    print(f"  Explicit group colors: {explicit_group_colors}")

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

# --- Build adjacency matrix ---
edge_mode_str = f" [{args.edge_mode}]"
print(f"\n--- Building adjacency matrix ({args.matrix_type}{edge_mode_str}) ---")

# Get feature data per read (with lengths for weighting)
seq_feature_data = in_data.groupby('sequence').apply(
    lambda x: list(zip(x['feature'], x['length']))
).to_dict()

seq_names = sorted(seq_feature_data.keys())

# Detect if features are merged (colon-separated) and count layers
all_features = sorted(in_data['feature'].unique())
n_layers = detect_feature_layers(all_features)
print(f"  Detected {n_layers} feature layer(s)")

# Get edges with weights (length of source feature)
def get_weighted_edges(features_with_lengths, edge_mode="directional"):
    """Get edges with weights based on edge counting mode."""
    if len(features_with_lengths) <= 1:
        return []
    edges = []
    for i in range(len(features_with_lengths) - 1):
        from_feat, from_len = features_with_lengths[i]
        to_feat, to_len = features_with_lengths[i + 1]

        if edge_mode == "directional":
            edges.append((from_feat, to_feat, from_len))
        elif edge_mode == "symmetric":
            avg_len = (from_len + to_len) / 2
            sorted_pair = tuple(sorted([from_feat, to_feat]))
            edges.append((sorted_pair[0], sorted_pair[1], avg_len))
    return edges


def _build_matrix_from_features(seq_names, read_feature_data, seq_length_dict, feature_list,
                                edge_mode, matrix_type, include_abundance, include_edges=True):
    """Shared matrix builder used by both build_layer_matrix and build_combined_matrix.

    Args:
        seq_names: Sorted list of read names (row order).
        read_feature_data: dict of read_name -> [(feature, length), ...].
        seq_length_dict: dict of read_name -> total read length.
        feature_list: Sorted list of unique features (column vocabulary).
        edge_mode: "directional" or "symmetric".
        matrix_type: One of "binary", "count", "count_log1p", "count_log1p_zscore",
                     "count_log1p_zscore_blockweight", "length_weighted".
        include_abundance: Whether to build the abundance block.
        include_edges: Whether to build the edge block.

    Returns:
        matrix: Combined (edge + abundance) or single-block matrix.
        n_edge_cols: Number of edge columns.
        n_abundance_cols: Number of abundance columns.
    """
    # Build edge matrix if requested
    if include_edges:
        all_pairs = []
        if edge_mode == "symmetric":
            for i, f1 in enumerate(feature_list):
                for f2 in feature_list[i+1:]:
                    all_pairs.append(f"{f1}->{f2}")
        else:
            for f1 in feature_list:
                for f2 in feature_list:
                    if f1 != f2:
                        all_pairs.append(f"{f1}->{f2}")

        pair_to_idx = {pair: i for i, pair in enumerate(all_pairs)}
        edge_matrix = np.zeros((len(seq_names), len(all_pairs)), dtype=np.float32)

        for i, read_name in enumerate(seq_names):
            feat_data = read_feature_data[read_name]
            read_len = seq_length_dict.get(read_name, 1)

            if matrix_type == "binary":
                edges = get_edges([f for f, _ in feat_data], edge_mode=edge_mode)
                for from_feat, to_feat in edges:
                    pair_name = f"{from_feat}->{to_feat}"
                    if pair_name in pair_to_idx:
                        edge_matrix[i, pair_to_idx[pair_name]] = 1

            elif matrix_type in ("count", "count_log1p", "count_log1p_zscore", "count_log1p_zscore_blockweight"):
                edges = get_edges([f for f, _ in feat_data], edge_mode=edge_mode)
                for from_feat, to_feat in edges:
                    pair_name = f"{from_feat}->{to_feat}"
                    if pair_name in pair_to_idx:
                        edge_matrix[i, pair_to_idx[pair_name]] += 1

            elif matrix_type == "length_weighted":
                edges = get_weighted_edges(feat_data, edge_mode=edge_mode)
                for from_feat, to_feat, weight in edges:
                    pair_name = f"{from_feat}->{to_feat}"
                    if pair_name in pair_to_idx:
                        edge_matrix[i, pair_to_idx[pair_name]] += weight / read_len

        if matrix_type in ("count_log1p", "count_log1p_zscore", "count_log1p_zscore_blockweight"):
            edge_matrix = np.log1p(edge_matrix)

        n_edge_cols = len(all_pairs)
    else:
        edge_matrix = None
        n_edge_cols = 0

    # Build abundance matrix if requested
    if include_abundance:
        feature_to_idx = {f: i for i, f in enumerate(feature_list)}
        abundance_matrix = np.zeros((len(seq_names), len(feature_list)), dtype=np.float32)

        for i, read_name in enumerate(seq_names):
            read_len = seq_length_dict.get(read_name, 1)
            feat_data = read_feature_data[read_name]

            if matrix_type == "binary":
                for feat, _length in feat_data:
                    if feat in feature_to_idx:
                        abundance_matrix[i, feature_to_idx[feat]] = 1

            elif matrix_type in ("count", "count_log1p", "count_log1p_zscore", "count_log1p_zscore_blockweight"):
                # Raw total bp per feature (length-based, not occurrence-based)
                feature_lengths = defaultdict(float)
                for feat, length in feat_data:
                    feature_lengths[feat] += length
                for feat, total_len in feature_lengths.items():
                    if feat in feature_to_idx:
                        abundance_matrix[i, feature_to_idx[feat]] = total_len

            else:  # length_weighted
                feature_lengths = defaultdict(float)
                for feat, length in feat_data:
                    feature_lengths[feat] += length
                for feat, total_len in feature_lengths.items():
                    if feat in feature_to_idx:
                        abundance_matrix[i, feature_to_idx[feat]] = total_len / read_len

        if matrix_type in ("count_log1p", "count_log1p_zscore", "count_log1p_zscore_blockweight"):
            abundance_matrix = np.log1p(abundance_matrix)

        if matrix_type in ("count_log1p_zscore", "count_log1p_zscore_blockweight"):
            if edge_matrix is not None:
                scaler_e = StandardScaler()
                edge_matrix = scaler_e.fit_transform(edge_matrix)

            scaler_a = StandardScaler()
            abundance_matrix = scaler_a.fit_transform(abundance_matrix)

            print(f"  Applied z-score normalization to all columns")

            if matrix_type == "count_log1p_zscore_blockweight":
                if edge_matrix is not None:
                    n_edge = edge_matrix.shape[1]
                    n_abund = abundance_matrix.shape[1]
                    if n_abund > 0:
                        scale_factor = np.sqrt(n_edge / n_abund)
                        abundance_matrix *= scale_factor
                        print(f"  Applied block reweighting: abundance scaled by {scale_factor:.2f}x for equal variance contribution")
                else:
                    print(f"  Skipping block reweighting (only one feature block present)")

        if edge_matrix is not None:
            matrix = np.hstack([edge_matrix, abundance_matrix])
        else:
            matrix = abundance_matrix
        n_abundance_cols = len(feature_list)
    elif include_edges:
        if matrix_type in ("count_log1p_zscore", "count_log1p_zscore_blockweight"):
            scaler_e = StandardScaler()
            edge_matrix = scaler_e.fit_transform(edge_matrix)
            print(f"  Applied z-score normalization to all columns")
            if matrix_type == "count_log1p_zscore_blockweight":
                print(f"  Skipping block reweighting (only one feature block present)")

        matrix = edge_matrix
        n_abundance_cols = 0
    else:
        raise ValueError("At least one of include_edges or include_abundance must be True")

    return matrix, n_edge_cols, n_abundance_cols


def build_layer_matrix(seq_names, seq_feature_data, seq_length_dict, layer_idx, n_layers,
                       edge_mode, matrix_type, include_abundance, include_edges=True):
    """Build edge and abundance matrices for a single layer.

    Returns:
        matrix: Combined edge + abundance matrix for this layer
        n_edge_cols: Number of edge columns
        n_abundance_cols: Number of abundance columns
        layer_features: List of unique features in this layer
    """
    # Extract layer-specific features for each read
    layer_read_data = {}
    for read_name, feat_data in seq_feature_data.items():
        layer_read_data[read_name] = get_layer_features(feat_data, layer_idx, n_layers)

    # Get unique features for this layer
    layer_features_set = set()
    for feat_data in layer_read_data.values():
        for feat, _ in feat_data:
            layer_features_set.add(feat)
    layer_features = sorted(layer_features_set)

    matrix, n_edge_cols, n_abundance_cols = _build_matrix_from_features(
        seq_names, layer_read_data, seq_length_dict, layer_features,
        edge_mode, matrix_type, include_abundance, include_edges)

    return matrix, n_edge_cols, n_abundance_cols, layer_features


def build_combined_matrix(seq_names, seq_feature_data, seq_length_dict, all_features,
                          edge_mode, matrix_type, include_abundance, include_edges=True):
    """Build edge and abundance matrices treating features as atomic (original approach).

    This treats merged features like "chr7:p_arm:nonsubtelomeric" as single atomic features,
    rather than splitting them into layers. Results in higher dimensionality but preserves
    the original feature combinations.

    Returns:
        matrix: Combined edge + abundance matrix
        n_edge_cols: Number of edge columns
        n_abundance_cols: Number of abundance columns
    """
    return _build_matrix_from_features(
        seq_names, seq_feature_data, seq_length_dict, all_features,
        edge_mode, matrix_type, include_abundance, include_edges)


# Build matrix based on selected mode
if args.matrix_mode == "layered" and n_layers > 1:
    # Layered mode: split colon-separated features into layers, build matrices per layer
    print(f"  Using layered matrix mode ({n_layers} layers)")

    all_layer_matrices = []
    total_edge_cols = 0
    total_abundance_cols = 0

    for layer_idx in range(n_layers):
        layer_matrix, n_edge, n_abund, layer_feats = build_layer_matrix(
            seq_names, seq_feature_data, seq_length_dict, layer_idx, n_layers,
            args.edge_mode, args.matrix_type, args.include_abundance, args.include_edges
        )
        all_layer_matrices.append(layer_matrix)
        total_edge_cols += n_edge
        total_abundance_cols += n_abund
        print(f"  Layer {layer_idx + 1}: {len(layer_feats)} features, {n_edge} edge cols, {n_abund} abundance cols")

    # Concatenate all layer matrices
    adj_matrix = np.hstack(all_layer_matrices)
    svd_feature_name_list = None  # Feature names not tracked for layered mode
    print("  Note: SVD feature names not available in layered mode")

else:
    # Combined mode: treat merged features as atomic (original approach)
    if n_layers > 1:
        print(f"  Using combined matrix mode (treating {n_layers}-layer features as atomic)")
    else:
        print(f"  Using combined matrix mode (single featureset)")

    adj_matrix, total_edge_cols, total_abundance_cols = build_combined_matrix(
        seq_names, seq_feature_data, seq_length_dict, all_features,
        args.edge_mode, args.matrix_type, args.include_abundance, args.include_edges
    )

    # Build feature name list matching exact column order of combined matrix
    svd_feature_name_list = []
    if args.include_edges:
        if args.edge_mode == "symmetric":
            for i, f1 in enumerate(all_features):
                for f2 in all_features[i+1:]:
                    svd_feature_name_list.append(f"edge:{f1}->{f2}")
        else:
            for f1 in all_features:
                for f2 in all_features:
                    if f1 != f2:
                        svd_feature_name_list.append(f"edge:{f1}->{f2}")
    if args.include_abundance:
        for f in all_features:
            svd_feature_name_list.append(f"abundance:{f}")

    print(f"  Unique features: {len(all_features)}")
    if args.include_edges:
        print(f"  Edge columns: {total_edge_cols}")
    if args.include_abundance:
        print(f"  Abundance columns: {total_abundance_cols}")

if args.include_edges:
    print(f"\nTotal edge dimensions: {total_edge_cols}")
if args.include_abundance:
    print(f"Total abundance dimensions: {total_abundance_cols}")
print(f"Final matrix shape: {adj_matrix.shape}")
n_nonzero = np.count_nonzero(adj_matrix)
n_total = adj_matrix.shape[0] * adj_matrix.shape[1]
sparsity = 1 - (n_nonzero / n_total)
print(f"Non-zero entries: {n_nonzero:,} ({100*(1-sparsity):.2f}% dense, {100*sparsity:.2f}% sparse)")

# Report feature usage statistics
col_usage = (adj_matrix != 0).sum(axis=0)
active_cols = (col_usage > 0).sum()
print(f"Active features: {active_cols:,} / {adj_matrix.shape[1]:,} ({100*active_cols/adj_matrix.shape[1]:.1f}%)")
if active_cols > 0:
    print(f"  Reads per feature: min={col_usage[col_usage > 0].min()}, median={int(np.median(col_usage[col_usage > 0]))}, max={col_usage.max()}")

# --- Optional dimensionality reduction ---
adj_matrix_full = adj_matrix  # Keep full matrix for reference
svd_components_export = None
svd_explained_variance_ratio_export = None
if args.reduce_dims and args.reduce_dims > 0:
    n_samples, n_features = adj_matrix.shape
    n_components = min(args.reduce_dims, n_samples - 1, n_features)

    if n_features > n_components:
        print(f"\n--- Reducing dimensions with truncated SVD ---")
        print(f"  Original dimensions: {n_features}")
        print(f"  Target dimensions: {n_components}")

        svd = TruncatedSVD(n_components=n_components, random_state=42)
        adj_matrix = svd.fit_transform(adj_matrix)

        # Report explained variance
        explained_var = svd.explained_variance_ratio_.sum()
        print(f"  Total explained variance: {explained_var:.1%}")

        # Report variance by component ranges
        cumvar = np.cumsum(svd.explained_variance_ratio_)
        var_50 = np.searchsorted(cumvar, 0.5) + 1
        var_90 = np.searchsorted(cumvar, 0.9) + 1
        var_95 = np.searchsorted(cumvar, 0.95) + 1
        print(f"  Components for 50% variance: {var_50}")
        print(f"  Components for 90% variance: {var_90}")
        print(f"  Components for 95% variance: {var_95}")

        # Report top singular values
        top_k = min(5, len(svd.singular_values_))
        top_sv = svd.singular_values_[:top_k]
        top_var = svd.explained_variance_ratio_[:top_k] * 100
        print(f"  Top {top_k} singular values: {', '.join(f'{v:.1f}' for v in top_sv)}")
        print(f"  Top {top_k} variance %: {', '.join(f'{v:.1f}%' for v in top_var)}")

        print(f"  Reduced matrix shape: {adj_matrix.shape}")

        # Stash SVD data for NPZ export
        svd_components_export = svd.components_
        svd_explained_variance_ratio_export = svd.explained_variance_ratio_

        # Generate SVD scree plot
        var_ratio = svd.explained_variance_ratio_
        cumvar_full = np.cumsum(var_ratio)
        n_plot = len(var_ratio)

        for bg_mode, suffix in get_backgrounds_to_generate():
            style = apply_plot_style(bg_mode)
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
            fig.patch.set_facecolor(style['bg_color'])
            ax1.set_facecolor(style['bg_color'])
            ax2.set_facecolor(style['bg_color'])

            # Left: individual explained variance (log scale for visibility)
            ax1.bar(range(1, n_plot + 1), var_ratio * 100, color='#60A5FA', alpha=0.7, width=1.0)
            ax1.set_xlabel('Component')
            ax1.set_ylabel('Explained Variance (%)')
            ax1.set_title('Explained Variance per Component')
            ax1.set_xlim(0, min(n_plot + 1, 100))  # Zoom to first 100
            if n_plot > 100:
                ax1.annotate(f'({n_plot} total components)', xy=(0.95, 0.95),
                           xycoords='axes fraction', ha='right', va='top',
                           color=style['text_color'], fontsize=9)

            # Right: cumulative variance with threshold lines
            ax2.plot(range(1, n_plot + 1), cumvar_full * 100, color='#40D392', linewidth=2)
            for thresh, label_text in [(50, '50%'), (90, '90%'), (95, '95%'), (99, '99%')]:
                idx = np.searchsorted(cumvar_full, thresh / 100)
                if idx < n_plot:
                    k_val = idx + 1
                    ax2.axhline(y=thresh, color=style['grid_color'], linestyle='--', alpha=0.5)
                    ax2.axvline(x=k_val, color='#F07167', linestyle='--', alpha=0.4)
                    ax2.annotate(f'{label_text} at k={k_val}', xy=(k_val, thresh),
                               xytext=(10, -5), textcoords='offset points',
                               color=style['text_color'], fontsize=8)
            ax2.set_xlabel('Number of Components')
            ax2.set_ylabel('Cumulative Explained Variance (%)')
            ax2.set_title('Cumulative Explained Variance')
            ax2.set_ylim(0, 105)
            ax2.set_xlim(0, min(n_plot + 1, 100))

            plt.tight_layout()
            scree_file = f"{args.output_prefix}{suffix}.svd_scree.pdf"
            plt.savefig(scree_file, dpi=150, bbox_inches='tight', facecolor=style['bg_color'])
            plt.close()
            print(f"  Saved SVD scree plot to: {scree_file}")
    else:
        print(f"\n--- Skipping dimensionality reduction (already at {n_features} dims) ---")

# --- Hierarchical clustering ---
print(f"\n--- Performing hierarchical clustering ---")
dist_matrix = pdist(adj_matrix, metric='euclidean')
linkage_matrix = linkage(dist_matrix, method=args.linkage_method)

# --- Determine optimal number of clusters ---
print(f"\n--- Analyzing cluster structure ---")

# Pre-compute arrays for faster k-optimization
seq_names_arr = np.array(seq_names)
read_groups_arr = np.array([seq_to_group[r] for r in seq_names])

# Fast enrichment check for k-optimization (returns enrichment type and max purity)
def fast_enrichment_check(cluster_mask, read_groups, group_totals_dict, control_grp, all_grps):
    """Fast enrichment calculation for k-optimization loop.

    Returns: (enrichment_label, max_purity) where max_purity is the highest group percentage.

    Uses odds ratio to determine enrichment, which works correctly even when
    one group dominates the dataset (e.g., 90% tumor vs 10% normal).
    """
    cluster_groups = read_groups[cluster_mask]
    group_counts = Counter(cluster_groups)
    total_in = len(cluster_groups)
    if total_in == 0 or not all_grps:
        return "mixed", 0.0

    # Calculate max purity (highest percentage from any single group)
    max_purity = max((group_counts.get(g, 0) / total_in for g in all_grps), default=0.0)

    # Calculate expected proportions
    total_sequences = sum(group_totals_dict.values())

    # Use odds ratio to detect enrichment (works with skewed group sizes)
    # Odds ratio > 1.5 indicates meaningful enrichment
    best_enrichment = "mixed"
    best_odds = 1.0

    for grp in all_grps:
        in_cluster = group_counts.get(grp, 0)
        out_cluster = group_totals_dict[grp] - in_cluster
        other_in = total_in - in_cluster
        other_out = total_sequences - group_totals_dict[grp] - other_in

        # Calculate odds ratio with small pseudocount to avoid division by zero
        odds_ratio = ((in_cluster + 0.5) * (other_out + 0.5)) / ((out_cluster + 0.5) * (other_in + 0.5))

        if odds_ratio > best_odds and odds_ratio > 1.5:
            best_odds = odds_ratio
            best_enrichment = f"{grp}-enriched"

    return best_enrichment, max_purity

# Calculate clustering metrics for different k values
if args.n_clusters is None:
    # Try different numbers of clusters and evaluate
    # Upper bound: min of max_k or n_reads/10 (need at least 10 reads per cluster on average)
    max_k_limit = min(args.max_k + 1, len(seq_names) // 10)
    k_range = list(range(args.min_k, max_k_limit))
    n_samples = len(seq_names)

    if not k_range:
        print(f"ERROR: Empty k-range (min_k={args.min_k} >= max_k_limit={max_k_limit})")
        print(f"  Try lowering --min-k or ensuring you have enough reads (need at least 10x max_k reads)")
        sys.exit(1)

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
            if n_samples > args.silhouette_sample_size:
                silhouette = silhouette_score(adj_matrix, labels, sample_size=args.silhouette_sample_size, random_state=42)
                cosine_silhouette = silhouette_score(adj_matrix, labels, metric='cosine', sample_size=args.silhouette_sample_size, random_state=42)
            else:
                silhouette = silhouette_score(adj_matrix, labels, random_state=42)
                cosine_silhouette = silhouette_score(adj_matrix, labels, metric='cosine', random_state=42)
            calinski_harabasz = calinski_harabasz_score(adj_matrix, labels)
        else:
            silhouette = 0
            cosine_silhouette = 0
            calinski_harabasz = 0

        # Count clusters meeting minimum size
        cluster_sizes = Counter(labels)
        valid_clusters = sum(1 for size in cluster_sizes.values() if size >= args.min_cluster_size)

        # Calculate enrichment diversity with strength categories
        enrichments = []
        purities = []
        perfect_enriched = 0
        strong_enriched = 0
        any_enriched = 0
        # Track reads in enriched clusters
        perfect_sequences = 0
        strong_sequences = 0
        any_enriched_sequences = 0

        for cluster_id in range(1, k + 1):
            cluster_mask = labels_arr == cluster_id
            cluster_size = cluster_mask.sum()
            if cluster_size >= args.min_cluster_size:
                enrich, purity = fast_enrichment_check(cluster_mask, read_groups_arr,
                                                       group_totals, control_group, all_groups)
                enrichments.append(enrich)
                purities.append(purity)

                if enrich != "mixed":
                    any_enriched += 1
                    any_enriched_sequences += cluster_size
                    if purity >= args.perfect_threshold:
                        perfect_enriched += 1
                        perfect_sequences += cluster_size
                    if purity >= args.strong_threshold:
                        strong_enriched += 1
                        strong_sequences += cluster_size

        # Count enrichments by group
        group_enriched_counts = {g: sum(1 for e in enrichments if e == f"{g}-enriched") for g in all_groups}
        mixed = enrichments.count("mixed")

        # Calculate enrichment ratios
        enriched_ratio = any_enriched / valid_clusters if valid_clusters > 0 else 0
        perfect_ratio = perfect_enriched / valid_clusters if valid_clusters > 0 else 0
        strong_ratio = strong_enriched / valid_clusters if valid_clusters > 0 else 0

        # Calculate read percentages
        total_sequences = len(seq_names)
        perfect_sequences_pct = (perfect_sequences / total_sequences * 100) if total_sequences > 0 else 0
        strong_sequences_pct = (strong_sequences / total_sequences * 100) if total_sequences > 0 else 0
        any_enriched_sequences_pct = (any_enriched_sequences / total_sequences * 100) if total_sequences > 0 else 0

        # Composite: 50% cluster quality (silhouette), 10% any enrichment, 40% perfect purity
        # Weights favor biologically pure clusters over moderate enrichment
        silhouette_norm = (silhouette + 1) / 2  # normalize from [-1,1] to [0,1]
        composite_score = (0.5 * silhouette_norm +
                          0.1 * enriched_ratio +
                          0.4 * perfect_ratio)

        stats = {
            'k': k,
            'silhouette': silhouette,
            'cosine_silhouette': cosine_silhouette,
            'calinski_harabasz': calinski_harabasz,
            'valid_clusters': valid_clusters,
            'any_enriched': any_enriched,
            'strong_enriched': strong_enriched,
            'perfect_enriched': perfect_enriched,
            'enriched_ratio': enriched_ratio,
            'strong_ratio': strong_ratio,
            'perfect_ratio': perfect_ratio,
            'any_enriched_sequences': any_enriched_sequences,
            'strong_sequences': strong_sequences,
            'perfect_sequences': perfect_sequences,
            'any_enriched_sequences_pct': any_enriched_sequences_pct,
            'strong_sequences_pct': strong_sequences_pct,
            'perfect_sequences_pct': perfect_sequences_pct,
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

    # Generate k-selection diagnostic plot (2x4 grid)
    for bg_mode, suffix in get_backgrounds_to_generate():
        style = apply_plot_style(bg_mode)
        fig, axes = plt.subplots(2, 4, figsize=(20, 10))
        fig.patch.set_facecolor(style['bg_color'])
        for ax in axes.flat:
            ax.set_facecolor(style['bg_color'])

        # 1. Silhouette scores (higher is better)
        ax1 = axes[0, 0]
        ax1.plot(stats_df['k'], stats_df['silhouette'], '-o', markersize=3, color='#60A5FA', label='Euclidean')
        ax1.plot(stats_df['k'], stats_df['cosine_silhouette'], '-o', markersize=3, color='#F07167', label='Cosine')
        ax1.set_xlabel('Number of clusters (k)')
        ax1.set_ylabel('Silhouette Score')
        ax1.set_title('Silhouette Score (higher = better)')
        best_silhouette_k = int(stats_df.loc[stats_df['silhouette'].idxmax(), 'k'])
        best_cosine_k = int(stats_df.loc[stats_df['cosine_silhouette'].idxmax(), 'k'])
        ax1.axvline(x=best_silhouette_k, color='#60A5FA', linestyle='--', alpha=0.5, label=f'Best Euclidean k={best_silhouette_k}')
        ax1.axvline(x=best_cosine_k, color='#F07167', linestyle='--', alpha=0.5, label=f'Best Cosine k={best_cosine_k}')
        ax1.legend(fontsize=7)

        # 2. Calinski-Harabasz index (higher is better)
        ax2 = axes[0, 1]
        ax2.plot(stats_df['k'], stats_df['calinski_harabasz'], '-o', markersize=3, color='#60A5FA')
        ax2.set_xlabel('Number of clusters (k)')
        ax2.set_ylabel('Calinski-Harabasz Index')
        ax2.set_title('Calinski-Harabasz (higher = better)')
        best_ch_k = int(stats_df.loc[stats_df['calinski_harabasz'].idxmax(), 'k'])
        ax2.axvline(x=best_ch_k, color='#40D392', linestyle='--', alpha=0.5, label=f'Best k={best_ch_k}')
        ax2.legend()

        # 3. Enrichment ratios
        ax3 = axes[0, 2]
        ax3.plot(stats_df['k'], stats_df['enriched_ratio'], '-o', markersize=3, label='Any enriched', color='#60A5FA')
        ax3.plot(stats_df['k'], stats_df['strong_ratio'], '-o', markersize=3, label=f'Strong (>={int(args.strong_threshold*100)}%)', color='#FBBF24')
        ax3.plot(stats_df['k'], stats_df['perfect_ratio'], '-o', markersize=3, label=f'Perfect (>={int(args.perfect_threshold*100)}%)', color='#F07167')
        ax3.set_xlabel('Number of clusters (k)')
        ax3.set_ylabel('Ratio (enriched / valid clusters)')
        ax3.set_title('Enrichment Ratios')
        ax3.set_ylim(0, 1)
        ax3.legend()

        # 4. Enrichment counts (absolute)
        ax4 = axes[0, 3]
        ax4.plot(stats_df['k'], stats_df['any_enriched'], '-o', markersize=3, label='Any enriched', color='#60A5FA')
        ax4.plot(stats_df['k'], stats_df['strong_enriched'], '-o', markersize=3, label=f'Strong (>={int(args.strong_threshold*100)}%)', color='#FBBF24')
        ax4.plot(stats_df['k'], stats_df['perfect_enriched'], '-o', markersize=3, label=f'Perfect (>={int(args.perfect_threshold*100)}%)', color='#F07167')
        ax4.set_xlabel('Number of clusters (k)')
        ax4.set_ylabel('Number of clusters')
        ax4.set_title('Enrichment Counts (absolute)')
        ax4.legend()

        # 5. Enrichment by group
        ax5 = axes[1, 0]
        for g in all_groups:
            col = f'{g}_enriched'
            if col in stats_df.columns:
                ax5.plot(stats_df['k'], stats_df[col], '-o', markersize=3, label=f'{g}-enriched', color=group_colors.get(g, None))
        ax5.plot(stats_df['k'], stats_df['mixed'], '-o', markersize=3, label='Mixed', color='#545454')
        ax5.set_xlabel('Number of clusters (k)')
        ax5.set_ylabel('Number of clusters')
        ax5.set_title('Enrichment by Group')
        ax5.legend()

        # 6. Reads in enriched clusters (absolute counts) - log scale
        ax6 = axes[1, 1]
        ax6.plot(stats_df['k'], stats_df['any_enriched_sequences'], '-o', markersize=3, label='Any enriched', color='#60A5FA')
        ax6.plot(stats_df['k'], stats_df['strong_sequences'], '-o', markersize=3, label=f'Strong (>={int(args.strong_threshold*100)}%)', color='#FBBF24')
        ax6.plot(stats_df['k'], stats_df['perfect_sequences'], '-o', markersize=3, label=f'Perfect (>={int(args.perfect_threshold*100)}%)', color='#F07167')
        ax6.set_xlabel('Number of clusters (k)')
        ax6.set_ylabel('Number of reads')
        ax6.set_title('Reads in Enriched Clusters (count)')
        ax6.legend()
        ax6.set_yscale('log')

        # 7. Reads in enriched clusters (percentages)
        ax7 = axes[1, 2]
        ax7.plot(stats_df['k'], stats_df['any_enriched_sequences_pct'], '-o', markersize=3, label='Any enriched', color='#60A5FA')
        ax7.plot(stats_df['k'], stats_df['strong_sequences_pct'], '-o', markersize=3, label=f'Strong (>={int(args.strong_threshold*100)}%)', color='#FBBF24')
        ax7.plot(stats_df['k'], stats_df['perfect_sequences_pct'], '-o', markersize=3, label=f'Perfect (>={int(args.perfect_threshold*100)}%)', color='#F07167')
        ax7.set_xlabel('Number of clusters (k)')
        ax7.set_ylabel('Percent of reads')
        ax7.set_title('Reads in Enriched Clusters (%)')
        ax7.set_ylim(0, 100)
        ax7.legend()

        # 8. Composite score (our recommended metric) - last position
        ax8 = axes[1, 3]
        ax8.plot(stats_df['k'], stats_df['composite_score'], '-o', markersize=3, color='#60A5FA')
        ax8.set_xlabel('Number of clusters (k)')
        ax8.set_ylabel('Composite Score')
        ax8.set_title('Composite Score (silhouette + enrichment)')
        best_composite_k = int(stats_df.loc[stats_df['composite_score'].idxmax(), 'k'])
        ax8.axvline(x=best_composite_k, color='#40D392', linestyle='--', alpha=0.5, label=f'Best k={best_composite_k}')
        ax8.axhline(y=stats_df['composite_score'].max(), color='#F07167', linestyle='--', alpha=0.3)
        ax8.legend()

        plt.tight_layout()
        k_plot_file = f"{args.output_prefix}{suffix}.k_selection.pdf"
        plt.savefig(k_plot_file, dpi=150, bbox_inches='tight', facecolor=style['bg_color'])
        plt.close()
        print(f"  Saved k-selection plot to: {k_plot_file}")

    # Calculate knee point (diminishing returns)
    # Use FIXED normalization bounds to make knee detection independent of max-k
    k_values = stats_df['k'].values
    scores = stats_df['composite_score'].values

    # Kneedle-style knee detection: normalize to observed range, find max distance from diagonal
    # Note: knee point may vary slightly with max-k due to normalization; diagnostic plot shows this
    k_min_obs, k_max_obs = k_values.min(), k_values.max()
    score_min, score_max = scores.min(), scores.max()

    # Normalize both to [0, 1] using observed range
    k_norm = (k_values - k_min_obs) / (k_max_obs - k_min_obs) if k_max_obs > k_min_obs else np.zeros_like(k_values)
    score_norm = (scores - score_min) / (score_max - score_min) if score_max > score_min else np.zeros_like(scores)

    # Knee = point with maximum perpendicular distance from diagonal
    # This is equivalent to finding where (score_norm - k_norm) is maximized
    knee_distance = score_norm - k_norm

    # Apply smoothing to reduce noise (window = 20% of k range for adaptive stability)
    k_range = k_max_obs - k_min_obs
    window_size = max(3, int(k_range * 0.2))  # At least 3 for meaningful smoothing
    knee_smooth = pd.Series(knee_distance).rolling(window_size, center=True, min_periods=1).mean().values
    best_knee_idx = np.argmax(knee_smooth)
    best_knee_k = int(k_values[best_knee_idx])

    # Also compute raw (unsmoothed) knee for comparison
    best_knee_raw_idx = np.argmax(knee_distance)
    best_knee_raw_k = int(k_values[best_knee_raw_idx])

    # Generate knee diagnostic plot
    for bg_mode, suffix in get_backgrounds_to_generate():
        style = apply_plot_style(bg_mode)
        fig_knee, axes_knee = plt.subplots(2, 2, figsize=(12, 10))
        fig_knee.patch.set_facecolor(style['bg_color'])
        for ax in axes_knee.flat:
            ax.set_facecolor(style['bg_color'])

        # 1. Raw composite score
        ax1 = axes_knee[0, 0]
        ax1.plot(k_values, scores, 'b-', linewidth=1.5)
        ax1.axvline(x=best_knee_k, color='r', linestyle='--', label=f'Knee k={best_knee_k}')
        ax1.set_xlabel('Number of clusters (k)')
        ax1.set_ylabel('Composite Score')
        ax1.set_title('1. Raw Composite Score')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # 2. Normalized values with diagonal
        ax2 = axes_knee[0, 1]
        ax2.plot(k_norm, score_norm, 'b-', linewidth=1.5, label='Normalized curve')
        ax2.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Diagonal (no improvement)')
        knee_k_norm = k_norm[best_knee_idx]
        knee_score_norm = score_norm[best_knee_idx]
        ax2.scatter([knee_k_norm], [knee_score_norm], color='r', s=100, zorder=5, label=f'Knee k={best_knee_k}')
        ax2.plot([knee_k_norm, knee_k_norm], [knee_k_norm, knee_score_norm], 'r-', linewidth=2, alpha=0.7)
        ax2.set_xlabel(f'Normalized k (range: {int(k_min_obs)}-{int(k_max_obs)})')
        ax2.set_ylabel('Normalized Score')
        ax2.set_title('2. Kneedle: Observed Range Normalization')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_xlim(-0.05, 1.05)
        ax2.set_ylim(-0.05, 1.05)

        # 3. Composite-knee distance (raw and smoothed)
        ax3 = axes_knee[1, 0]
        ax3.plot(k_values, knee_distance, 'b-', alpha=0.5, label='Raw distance')
        ax3.plot(k_values, knee_smooth, 'b-', linewidth=2, label=f'Smoothed (window={window_size})')
        ax3.axvline(x=best_knee_k, color='r', linestyle='--', label=f'Smoothed k={best_knee_k}')
        ax3.axvline(x=best_knee_raw_k, color='orange', linestyle=':', linewidth=2, label=f'Raw k={best_knee_raw_k}')
        ax3.axhline(y=0, color='k', linestyle='-', alpha=0.3)
        ax3.set_xlabel('Number of clusters (k)')
        ax3.set_ylabel('Distance from Diagonal')
        ax3.set_title('3. Composite-Knee Distance')
        ax3.legend()
        ax3.grid(True, alpha=0.3)

        # 4. Enrichment vs k with composite-knee point
        ax4 = axes_knee[1, 1]
        ax4.plot(k_values, stats_df['enriched_ratio'].values * 100, 'g-', label='Any enriched %')
        ax4.plot(k_values, stats_df['strong_ratio'].values * 100, 'orange', label='Strong (≥80%) %')
        ax4.axvline(x=best_knee_raw_k, color='orange', linestyle=':', linewidth=2, label=f'Composite-knee k={best_knee_raw_k}')
        ax4.set_xlabel('Number of clusters (k)')
        ax4.set_ylabel('Enrichment Ratio (%)')
        ax4.set_title('4. Enrichment at Composite-Knee')
        ax4.legend()
        ax4.grid(True, alpha=0.3)

        plt.tight_layout()
        knee_plot_file = f"{args.output_prefix}{suffix}.composite_knee_diagnostic.pdf"
        plt.savefig(knee_plot_file, dpi=150, bbox_inches='tight', facecolor=style['bg_color'])
        plt.close()
        print(f"  Saved composite-knee diagnostic plot to: {knee_plot_file}")

    # Find best k by each metric
    print(f"\n  Optimal k by metric:")
    print(f"    Silhouette:        k={best_silhouette_k} (score={stats_df['silhouette'].max():.4f})")
    print(f"    Cosine Silhouette: k={best_cosine_k} (score={stats_df['cosine_silhouette'].max():.4f})")
    print(f"    Calinski-Harabasz: k={best_ch_k} (score={stats_df['calinski_harabasz'].max():.1f})")
    print(f"    Composite:         k={best_composite_k} (score={stats_df['composite_score'].max():.4f})")
    print(f"    Composite-knee:    k={best_knee_raw_k} (raw), k={best_knee_k} (smoothed, window={window_size})")

    # Select k based on chosen metric
    if args.k_selection == "silhouette":
        selected_k = best_silhouette_k
        selection_metric = "silhouette"
    elif args.k_selection == "calinski":
        selected_k = best_ch_k
        selection_metric = "Calinski-Harabasz"
    elif args.k_selection == "composite-knee":
        selected_k = best_knee_k  # Use smoothed for stability
        selection_metric = "composite-knee (diminishing returns)"
    else:  # composite (default)
        selected_k = best_composite_k
        selection_metric = "composite score"

    # Report enrichment stats at selected k
    best_row = stats_df[stats_df['k'] == selected_k].iloc[0]
    print(f"\n  At k={selected_k}:")
    print(f"    Valid clusters:    {int(best_row['valid_clusters'])}")
    print(f"    Any enriched:      {int(best_row['any_enriched'])} ({best_row['enriched_ratio']*100:.1f}%)")
    print(f"    Strong enriched:   {int(best_row['strong_enriched'])} ({best_row['strong_ratio']*100:.1f}%)")
    print(f"    Perfect enriched:  {int(best_row['perfect_enriched'])} ({best_row['perfect_ratio']*100:.1f}%)")

    print(f"\n  Using k={selected_k} (based on {selection_metric})")

    n_clusters = selected_k
else:
    n_clusters = args.n_clusters
    print(f"  Using specified k: {n_clusters}")

# --- Cut tree and analyze clusters ---
print(f"\n--- Cutting tree into {n_clusters} clusters ---")
cluster_labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')

# Analyze each cluster
cluster_analysis = []
cluster_reads_dict = {}
seq_centroid_distances = {}  # Store centroid distance for each read

# Pre-compute arrays for faster cluster analysis
cluster_labels_arr = np.array(cluster_labels)
read_samples_arr = np.array([seq_to_sample[r] for r in seq_names])

for cluster_id in range(1, n_clusters + 1):
    # Use numpy boolean indexing (faster than list comprehension)
    cluster_mask = cluster_labels_arr == cluster_id
    if cluster_mask.sum() < args.min_cluster_size:
        continue

    cluster_indices = np.where(cluster_mask)[0]
    cluster_reads = seq_names_arr[cluster_mask].tolist()
    cluster_samples = read_samples_arr[cluster_mask].tolist()

    cluster_reads_dict[cluster_id] = cluster_reads

    # Calculate enrichment using appropriate method
    if args.comparison_mode == "two-group":
        stats = calculate_enrichment_two_group(cluster_samples, sample_to_group, control_group, group_totals)
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
        seq_centroid_distances[read] = distances_to_centroid[i]

    # Build cluster analysis record
    cluster_record = {
        'cluster_id': cluster_id,
        'size': len(cluster_reads),
        'odds_ratio': stats['odds_ratio'],
        'p_value': stats['p_value'],
        'enrichment': stats['enrichment'],
        'centroid_read': centroid_read,
        'centroid_sample': seq_to_sample[centroid_read],
        'centroid_group': sample_to_group.get(seq_to_sample[centroid_read], seq_to_sample[centroid_read])
    }

    # Add group/sample-specific columns
    if args.comparison_mode == 'per-sample':
        # In per-sample mode, stats contains sample-level data
        for s in sample_labels:
            cluster_record[f'{s}_count'] = stats['group_counts'].get(s, 0)
            cluster_record[f'{s}_pct'] = stats['group_pcts'].get(s, 0)
            cluster_record[f'{s}_pval'] = stats['all_p_values'].get(s, 1.0)
            cluster_record[f'{s}_odds'] = stats['all_odds_ratios'].get(s, 1.0)

        # Also compute group-level stats for visualization
        cluster_groups = [sample_to_group.get(s, s) for s in cluster_samples]
        group_counts_cluster = Counter(cluster_groups)
        total_in_cluster = len(cluster_samples)
        for g in all_groups:
            cluster_record[f'{g}_count'] = group_counts_cluster.get(g, 0)
            cluster_record[f'{g}_pct'] = (group_counts_cluster.get(g, 0) / total_in_cluster * 100) if total_in_cluster > 0 else 0

        # Also run group-level Fisher's test for comparison
        if len(all_groups) == 2:
            group_stats = calculate_enrichment_two_group(cluster_samples, sample_to_group, control_group, group_totals)
            cluster_record['group_enrichment'] = group_stats['enrichment']
            cluster_record['group_p_value'] = group_stats['p_value']
            cluster_record['group_odds_ratio'] = group_stats['odds_ratio']
    else:
        # In group modes, stats contains group-level data
        for g in all_groups:
            cluster_record[f'{g}_count'] = stats['group_counts'].get(g, 0)
            cluster_record[f'{g}_pct'] = stats['group_pcts'].get(g, 0)

    cluster_analysis.append(cluster_record)

# --- Compute cluster-level dendrogram ---
# This creates a dendrogram showing how clusters relate to each other
# based on their centroid feature vectors
print(f"\n--- Computing cluster-level dendrogram ---")

from scipy.cluster.hierarchy import linkage as scipy_linkage
from scipy.spatial.distance import pdist

# Get unique cluster IDs (excluding noise cluster -1)
cluster_ids_ordered = sorted([c for c in set(cluster_labels) if c != -1])

# Compute cluster centroids efficiently
# Pre-compute cluster masks to avoid repeated boolean operations
if len(cluster_ids_ordered) > 0:
    n_clusters_valid = len(cluster_ids_ordered)
    n_features = adj_matrix.shape[1]

    # Pre-compute indices for each cluster (faster than repeated boolean masking)
    cluster_indices = {cid: np.where(cluster_labels == cid)[0] for cid in cluster_ids_ordered}

    # Compute centroids using pre-computed indices
    cluster_centroids = np.zeros((n_clusters_valid, n_features), dtype=adj_matrix.dtype)
    for i, cid in enumerate(cluster_ids_ordered):
        indices = cluster_indices[cid]
        cluster_centroids[i] = adj_matrix[indices].mean(axis=0)

    print(f"  Computed centroids for {len(cluster_ids_ordered)} clusters")
else:
    cluster_centroids = np.array([])
    print(f"  No valid clusters for centroid computation")

# Compute pairwise distances between cluster centroids
if len(cluster_centroids) > 1:
    cluster_distances = pdist(cluster_centroids, metric='euclidean')
    # Compute hierarchical linkage on clusters
    cluster_linkage = scipy_linkage(cluster_distances, method=args.linkage_method)
    print(f"  Computed cluster-level linkage using '{args.linkage_method}' method")
else:
    cluster_linkage = None
    print(f"  Only 1 cluster, skipping linkage computation")

# Create DataFrame and apply FDR correction
cluster_df = pd.DataFrame(cluster_analysis)

# Apply FDR correction to p-values
print(f"\n--- Applying FDR Correction (Benjamini-{'Hochberg' if args.fdr_method == 'bh' else 'Yekutieli'}) ---")
raw_pvals = cluster_df['p_value'].values
q_values = false_discovery_control(raw_pvals, method=args.fdr_method)
cluster_df['q_value'] = q_values

# Store original enrichment labels (based on raw p < 0.05)
cluster_df['enrichment_raw'] = cluster_df['enrichment']

# Update enrichment labels based on FDR-corrected q-values
def update_enrichment_label(row):
    """Update enrichment label based on q-value threshold."""
    if row['q_value'] >= args.fdr_threshold:
        return 'mixed'
    # If q-value passes, keep the original enrichment direction
    if '-enriched' in row['enrichment_raw']:
        return row['enrichment_raw']
    return 'mixed'

cluster_df['enrichment'] = cluster_df.apply(update_enrichment_label, axis=1)

# Report FDR correction impact
n_raw_sig = (cluster_df['enrichment_raw'] != 'mixed').sum()
n_fdr_sig = (cluster_df['enrichment'] != 'mixed').sum()
print(f"  Raw p < 0.05: {n_raw_sig} enriched clusters")
print(f"  FDR q < {args.fdr_threshold}: {n_fdr_sig} enriched clusters")
if n_raw_sig > n_fdr_sig:
    n_lost = n_raw_sig - n_fdr_sig
    print(f"  {n_lost} cluster(s) lost significance after FDR correction")

# Sort by enrichment type and q-value
cluster_df = cluster_df.sort_values(['enrichment', 'q_value'])

# --- Output results ---
print(f"\n--- Cluster Summary ---")

# Build dynamic header based on groups or samples (per-sample mode)
if args.comparison_mode == 'per-sample':
    # Use sample names for columns
    col_items = sample_labels
else:
    # Use group names for columns
    col_items = all_groups

header_parts = ['Cluster', 'Size']
for item in col_items:
    header_parts.extend([item[:6], f'{item[:4]}%'])  # Truncate long names
header_parts.extend(['P-value', 'Q-value', 'Enrichment', 'Centroid'])
header_fmt = '{:<8} {:<6} ' + ' '.join(['{:<6}'] * (len(col_items) * 2)) + ' {:<10} {:<10} {:<20} {:<10}'
print(header_fmt.format(*header_parts))
print("-" * (105 + 12 * len(col_items)))

for _, row in cluster_df.iterrows():
    if args.comparison_mode == 'per-sample':
        centroid_label = row['centroid_sample']
    else:
        centroid_label = row.get('centroid_group', sample_to_group.get(row['centroid_sample'], row['centroid_sample']))

    # Flag if centroid disagrees with enrichment
    flag = ""
    expected_item = row['enrichment'].replace('-enriched', '') if '-enriched' in row['enrichment'] else None
    if expected_item and centroid_label != expected_item:
        flag = " ⚠️"

    row_values = [row['cluster_id'], row['size']]
    for item in col_items:
        row_values.append(int(row.get(f'{item}_count', 0)))
        row_values.append(f"{row.get(f'{item}_pct', 0):.1f}")
    row_values.extend([f"{row['p_value']:.2e}", f"{row['q_value']:.2e}", row['enrichment'][:18], f"{centroid_label}{flag}"])

    row_fmt = '{:<8} {:<6} ' + ' '.join(['{:<6}'] * (len(col_items) * 2)) + ' {:<10} {:<10} {:<20} {:<10}'
    print(row_fmt.format(*row_values))

# Save cluster analysis
analysis_file = f"{args.output_prefix}.cluster_analysis.tsv"
cluster_df.to_csv(analysis_file, sep='\t', index=False)
print(f"\n  Saved cluster analysis to: {analysis_file}")

# Save read assignments with stats (sorted by cluster, then centroid distance)
# Build list of read records
seq_records = []
for i, read in enumerate(seq_names):
    cluster = cluster_labels[i]
    sample = seq_to_sample[read]
    group = sample_to_group.get(sample, sample)
    centroid_dist = seq_centroid_distances.get(read, np.nan)
    read_len = seq_length_dict.get(read, 0)
    read_span = seq_span_dict.get(read, 0)
    seq_records.append({
        'read': read,
        'cluster': cluster,
        'sample': sample,
        'group': group,
        'centroid_distance': centroid_dist,
        'read_length': read_len,
        'read_span': read_span
    })

# Create DataFrame and sort by cluster, then centroid_distance
assignments = pd.DataFrame(seq_records)
assignments = assignments.sort_values(['cluster', 'centroid_distance'], ascending=[True, True])

# Add rank within cluster (1 = closest to centroid)
assignments['rank'] = assignments.groupby('cluster').cumcount() + 1

# Reorder columns
# read_length = annotated length after excluding features (used for clustering)
# read_span = full coordinate range (actual read length for visualization)
assignments = assignments[['read', 'cluster', 'sample', 'group', 'centroid_distance', 'read_length', 'read_span', 'rank']]

assignments_file = f"{args.output_prefix}.sequence_assignments.tsv"
assignments.to_csv(assignments_file, sep='\t', index=False)
print(f"  Saved read assignments to: {assignments_file}")

# Save feature matrix data for cluster_plot.py
# Includes cluster-level linkage for drawing cluster dendrograms
matrix_file = f"{args.output_prefix}.feature_matrix.npz"
save_dict = {
    'adj_matrix': adj_matrix,
    'seq_names': np.array(seq_names),
    'cluster_labels': cluster_labels,
    'linkage_method': args.linkage_method,
    'cluster_ids_ordered': np.array(cluster_ids_ordered),
    'cluster_centroids': cluster_centroids,
}
if cluster_linkage is not None:
    save_dict['cluster_linkage'] = cluster_linkage
if svd_components_export is not None:
    save_dict['svd_components'] = svd_components_export
    save_dict['svd_explained_variance_ratio'] = svd_explained_variance_ratio_export
    if svd_feature_name_list is not None:
        save_dict['svd_feature_names'] = np.array(svd_feature_name_list)
        print(f"  Including SVD data: {svd_components_export.shape[0]} components x {len(svd_feature_name_list)} features")
    else:
        print(f"  Including SVD components but not feature names (layered mode)")
else:
    print(f"  SVD data not included (dimensionality reduction not used)")
np.savez(matrix_file, **save_dict)
print(f"  Saved feature matrix to: {matrix_file}")
if cluster_linkage is not None:
    print(f"    Includes cluster-level linkage ({len(cluster_ids_ordered)} clusters)")

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

from scipy.cluster.hierarchy import set_link_color_palette
from matplotlib.patches import Patch

# Color palette for clusters - use color + hatch combinations for many clusters
base_colors = list(plt.cm.tab20.colors)  # 20 distinct colors
hatches = ['', '///', '\\\\\\', '...', 'xxx', '+++', 'ooo', '---']  # 8 hatch patterns

def get_cluster_style(cluster_idx):
    """Get color and hatch pattern for a cluster index."""
    color_idx = cluster_idx % len(base_colors)
    hatch_idx = cluster_idx // len(base_colors)
    color = matplotlib.colors.rgb2hex(base_colors[color_idx])
    hatch = hatches[hatch_idx % len(hatches)]
    return color, hatch

# Create cluster color and hatch maps
cluster_color_map = {}
cluster_hatch_map = {}
for i, c in enumerate(sorted(set(cluster_labels))):
    color, hatch = get_cluster_style(i)
    cluster_color_map[c] = color
    cluster_hatch_map[c] = hatch

# For dendrogram coloring, use just colors (hatches not supported in dendrogram)
colors = plt.cm.tab20(np.linspace(0, 1, min(n_clusters, 20)))
set_link_color_palette([matplotlib.colors.rgb2hex(c) for c in colors])

# Create cluster to enrichment mapping
cluster_to_enrichment = dict(zip(cluster_df['cluster_id'], cluster_df['enrichment']))

for bg_mode, suffix in get_backgrounds_to_generate():
    style = apply_plot_style(bg_mode)

    # Use 3x2 layout for per-sample mode (to show both sample and group composition)
    if args.comparison_mode == 'per-sample':
        fig, axes = plt.subplots(3, 2, figsize=(16, 18))
    else:
        fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    fig.patch.set_facecolor(style['bg_color'])
    for ax in axes.flat:
        ax.set_facecolor(style['bg_color'])

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

    # Map enrichments to colors - use sample_colors in per-sample mode, group_colors otherwise
    def get_enrichment_color(e):
        if args.comparison_mode == 'per-sample':
            # Check sample colors first
            for s, color in sample_colors.items():
                if e == f'{s}-enriched':
                    return color
        # Then check group colors
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

    # Legend - dynamic based on groups or samples
    if args.comparison_mode == 'per-sample':
        legend_patches = [Patch(facecolor=c, label=f'{s}-enriched') for s, c in sample_colors.items()]
    else:
        legend_patches = [Patch(facecolor=c, label=f'{g}-enriched') for g, c in group_colors.items()]
    legend_patches.append(Patch(facecolor='#999999', label='Mixed'))
    ax2.legend(handles=legend_patches, loc='upper right', fontsize=7)

    # 3. Group composition per cluster
    ax3 = axes[1, 0]
    x = np.arange(len(cluster_df))

    # Use sample labels/colors in per-sample mode, group labels/colors otherwise
    if args.comparison_mode == 'per-sample':
        comp_items = sample_labels
        comp_colors = sample_colors
        comp_title = "Sample Composition per Cluster"
    else:
        comp_items = all_groups
        comp_colors = group_colors
        comp_title = "Group Composition per Cluster"

    n_comp = len(comp_items)
    width = 0.8 / n_comp

    for i, item in enumerate(comp_items):
        offset = (i - n_comp / 2 + 0.5) * width
        pct_col = f'{item}_pct'
        if pct_col in cluster_df.columns:
            ax3.bar(x + offset, cluster_df[pct_col], width, label=item, color=comp_colors.get(item, '#999999'))

    ax3.set_xlabel("Cluster")
    ax3.set_ylabel("Percentage")
    ax3.set_title(comp_title)
    ax3.set_xticks(x)
    ax3.set_xticklabels([f"C{c}" for c in cluster_df['cluster_id']], rotation=45)
    ax3.legend()
    if n_comp == 2:
        ax3.axhline(y=50, color='gray', linestyle='--', alpha=0.5)

    # 4. Q-value (FDR-corrected) distribution
    ax4 = axes[1, 1]
    qvals = cluster_df['q_value'].values
    colors_qval = [get_enrichment_color(e) for e in cluster_df['enrichment']]
    ax4.bar(range(len(qvals)), -np.log10(qvals + 1e-300), color=colors_qval)
    ax4.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.7, label='q=0.05')
    ax4.axhline(y=-np.log10(0.10), color='orange', linestyle='--', alpha=0.5, label='q=0.10')
    if args.fdr_threshold not in (0.05, 0.10):
        ax4.axhline(y=-np.log10(args.fdr_threshold), color='blue', linestyle='--', alpha=0.7, label=f'q={args.fdr_threshold}')
    ax4.set_xlabel("Cluster")
    ax4.set_ylabel("-log10(q-value)")
    ax4.set_title("FDR-Corrected Significance (per-sample)" if args.comparison_mode == 'per-sample' else "FDR-Corrected Significance")
    ax4.set_xticks(range(len(qvals)))
    ax4.set_xticklabels([f"C{c}" for c in cluster_df['cluster_id']], rotation=45)
    ax4.legend()

    # 5 & 6. Additional panels for per-sample mode: Group composition and group-level p-values
    if args.comparison_mode == 'per-sample':
        # 5. Group composition per cluster
        ax5 = axes[2, 0]
        x = np.arange(len(cluster_df))
        n_groups = len(all_groups)
        width = 0.8 / n_groups

        for i, g in enumerate(all_groups):
            offset = (i - n_groups / 2 + 0.5) * width
            pct_col = f'{g}_pct'
            if pct_col in cluster_df.columns:
                ax5.bar(x + offset, cluster_df[pct_col], width, label=g, color=group_colors.get(g, '#999999'))

        ax5.set_xlabel("Cluster")
        ax5.set_ylabel("Percentage")
        ax5.set_title("Group Composition per Cluster")
        ax5.set_xticks(x)
        ax5.set_xticklabels([f"C{c}" for c in cluster_df['cluster_id']], rotation=45)
        ax5.legend()
        if n_groups == 2:
            ax5.axhline(y=50, color='gray', linestyle='--', alpha=0.5)

        # 6. Group-level p-values (if available)
        ax6 = axes[2, 1]
        if 'group_p_value' in cluster_df.columns:
            group_pvals = cluster_df['group_p_value'].values
            group_enrichments = cluster_df['group_enrichment'].values

            def get_group_enrichment_color(e):
                for grp, clr in group_colors.items():
                    if e == f'{grp}-enriched':
                        return clr
                return '#999999'

            colors_group_pval = [get_group_enrichment_color(e) for e in group_enrichments]
            ax6.bar(range(len(group_pvals)), -np.log10(group_pvals + 1e-300), color=colors_group_pval)
            ax6.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='p=0.05')
            ax6.axhline(y=-np.log10(0.01), color='orange', linestyle='--', alpha=0.5, label='p=0.01')
            ax6.set_xlabel("Cluster")
            ax6.set_ylabel("-log10(p-value)")
            ax6.set_title("Enrichment Significance (group-level)")
            ax6.set_xticks(range(len(group_pvals)))
            ax6.set_xticklabels([f"C{c}" for c in cluster_df['cluster_id']], rotation=45)

            # Legend for group enrichment
            group_legend_patches = [Patch(facecolor=c, label=f'{g}-enriched') for g, c in group_colors.items()]
            group_legend_patches.append(Patch(facecolor='#999999', label='Mixed'))
            ax6.legend(handles=group_legend_patches, loc='upper right', fontsize=7)
        else:
            ax6.text(0.5, 0.5, "Group-level enrichment\nnot available\n(requires 2 groups)",
                     ha='center', va='center', fontsize=12, transform=ax6.transAxes)
            ax6.set_axis_off()

    plt.tight_layout()
    plot_file = f"{args.output_prefix}{suffix}.cluster_analysis.pdf"
    plt.savefig(plot_file, dpi=150, bbox_inches='tight', facecolor=style['bg_color'])
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
    leaf_samples = [seq_to_sample[seq_names[i]] for i in leaf_order]
    leaf_clusters = [cluster_labels[i] for i in leaf_order]
    leaf_enrichments = [cluster_to_enrichment.get(c, 'mixed') for c in leaf_clusters]

    # Color mappings for annotations
    if args.comparison_mode == 'per-sample':
        enrichment_colors = {f'{s}-enriched': c for s, c in sample_colors.items()}
    else:
        enrichment_colors = {f'{g}-enriched': c for g, c in group_colors.items()}
    enrichment_colors['mixed'] = '#CCCCCC'

    # Generate circular dendrogram for each background mode
    for bg_mode, suffix in get_backgrounds_to_generate():
        style = apply_plot_style(bg_mode)

        # Create figure with circular dendrogram
        fig_circ = plt.figure(figsize=(16, 14))
        fig_circ.patch.set_facecolor(style['bg_color'])
        ax_circ = fig_circ.add_subplot(111, polar=True)
        ax_circ.set_facecolor(style['bg_color'])

        # Plot each link in polar coordinates
        # scipy dendrogram icoord/dcoord are 4-point U-shapes: [x1, x1, x2, x2], [y1, y_merge, y_merge, y2]
        # In polar coords, we need: radial lines (constant theta) and arcs (constant radius)
        for xcoord, ycoord, color in zip(dend_colored['icoord'], dend_colored['dcoord'], dend_colored['color_list']):
            # Convert to polar coordinates
            theta1 = 2 * np.pi * xcoord[0] / (n_leaves * 10)
            theta2 = 2 * np.pi * xcoord[3] / (n_leaves * 10)
            r_leaf1 = max_dist - ycoord[0] + max_dist * 0.1  # bottom of left branch
            r_leaf2 = max_dist - ycoord[3] + max_dist * 0.1  # bottom of right branch
            r_merge = max_dist - ycoord[1] + max_dist * 0.1  # merge height (top of U)

            # Draw left vertical (radial) line
            ax_circ.plot([theta1, theta1], [r_leaf1, r_merge], color=color, linewidth=0.5)

            # Draw right vertical (radial) line
            ax_circ.plot([theta2, theta2], [r_leaf2, r_merge], color=color, linewidth=0.5)

            # Draw horizontal arc at merge height
            # Need to interpolate between theta1 and theta2 at constant radius
            n_arc_points = max(10, int(abs(theta2 - theta1) * 20))  # more points for larger arcs
            theta_arc = np.linspace(theta1, theta2, n_arc_points)
            r_arc = np.full(n_arc_points, r_merge)
            ax_circ.plot(theta_arc, r_arc, color=color, linewidth=0.5)

        # Calculate theta positions for leaves - must match dendrogram icoord positions
        # scipy dendrogram places leaves at x = 5, 15, 25, ... (i.e., 10*i + 5 for leaf i)
        # We convert these to angles the same way as the dendrogram branches
        theta_leaves = np.array([2 * np.pi * (10 * i + 5) / (n_leaves * 10) for i in range(n_leaves)])

        # Ring 1 (innermost): Sample of origin
        ring1_bottom = max_dist * 1.12
        ring1_height = max_dist * 0.06
        for theta, sample in zip(theta_leaves, leaf_samples):
            ax_circ.bar(theta, ring1_height, width=2 * np.pi / n_leaves, bottom=ring1_bottom,
                        color=sample_colors.get(sample, '#999999'), alpha=0.9, edgecolor='none')

        # Ring 2 (middle): Cluster number - with hatching for many clusters
        ring2_bottom = max_dist * 1.20
        ring2_height = max_dist * 0.06
        for theta, cluster in zip(theta_leaves, leaf_clusters):
            color = cluster_color_map.get(cluster, '#999999')
            hatch = cluster_hatch_map.get(cluster, '')
            ax_circ.bar(theta, ring2_height, width=2 * np.pi / n_leaves, bottom=ring2_bottom,
                        color=color, hatch=hatch, alpha=0.9, edgecolor='white', linewidth=0.1)

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

        # Add legends using figure legends (more reliable for multiple legends)
        # Sample legend (inner ring) - top right
        sample_patches = [Patch(facecolor=sample_colors[s], label=s) for s in sorted(sample_colors.keys())]
        leg1 = fig_circ.legend(handles=sample_patches, loc='upper left', bbox_to_anchor=(0.85, 0.95),
                               title='Sample (inner)', framealpha=0.9)

        # Enrichment legend (outer ring) - below sample legend
        enrich_patches = [Patch(facecolor=enrichment_colors[e], label=e) for e in sorted(enrichment_colors.keys())]
        n_samples = len(sample_patches)
        enrich_y = 0.95 - (n_samples + 2) * 0.035
        leg3 = fig_circ.legend(handles=enrich_patches, loc='upper left', bbox_to_anchor=(0.85, enrich_y),
                               title='Enrichment (outer)', framealpha=0.9)

        # Cluster legend (middle ring) - use hatches for uniqueness
        # Use more columns to reduce rows, larger patches to show hatches clearly
        n_enrichments = len(enrich_patches)
        cluster_y = enrich_y - (n_enrichments + 2) * 0.035

        cluster_patches = [Patch(facecolor=cluster_color_map[c], hatch=cluster_hatch_map[c],
                                 edgecolor='gray', label=f'C{c}')
                           for c in sorted(cluster_color_map.keys())]
        # Target ~4-5 rows max, so cols = ceil(n_clusters / 5)
        n_cols = max(6, (len(cluster_patches) + 4) // 5)
        leg2 = fig_circ.legend(handles=cluster_patches, loc='upper left', bbox_to_anchor=(0.85, cluster_y),
                               title='Cluster (middle)', ncol=n_cols, fontsize=7,
                               handlelength=2.5, handleheight=1.8, columnspacing=0.5,
                               labelspacing=0.8, framealpha=0.9)

        circ_dend_file = f"{args.output_prefix}{suffix}.circular_dendrogram.pdf"
        plt.savefig(circ_dend_file, dpi=150, bbox_inches='tight', facecolor=style['bg_color'])
        plt.close()
        print(f"  Saved circular dendrogram to: {circ_dend_file}")

# --- UMAP Visualization ---
if args.plot_umap and len(seq_names) > 10:
    print(f"\n--- Generating UMAP visualization ---")
    try:
        import umap

        # Fit UMAP
        reducer = umap.UMAP(
            n_neighbors=min(args.umap_neighbors, len(seq_names) - 1),
            min_dist=args.umap_min_dist,
            metric='euclidean',
            random_state=42
        )
        embedding = reducer.fit_transform(adj_matrix)
        print(f"  UMAP parameters: n_neighbors={args.umap_neighbors}, min_dist={args.umap_min_dist}")

        # Create read-to-cluster mapping
        seq_to_cluster = dict(zip(seq_names, cluster_labels))
        cluster_to_enrichment = dict(zip(cluster_df['cluster_id'], cluster_df['enrichment']))

        # Prepare data for plotting
        sample_list = [seq_to_sample[r] for r in seq_names]
        group_list = [sample_to_group.get(s, s) for s in sample_list]
        cluster_list = [seq_to_cluster[r] for r in seq_names]

        # Color mappings
        if args.comparison_mode == 'per-sample':
            enrichment_colors = {f'{s}-enriched': c for s, c in sample_colors.items()}
        else:
            enrichment_colors = {f'{g}-enriched': c for g, c in group_colors.items()}
        enrichment_colors['mixed'] = '#CCCCCC'

        # Generate distinct colors for clusters
        n_unique_clusters = len(set(cluster_list))
        cluster_cmap = plt.cm.tab20(np.linspace(0, 1, max(20, n_unique_clusters)))
        cluster_color_map = {c: matplotlib.colors.rgb2hex(cluster_cmap[i % 20]) for i, c in enumerate(sorted(set(cluster_list)))}

        # Generate UMAP plot for each background mode
        for bg_mode, suffix in get_backgrounds_to_generate():
            style = apply_plot_style(bg_mode)

            # Create 2x2 UMAP plot (group, enrichment, cluster numbers, cluster colors)
            fig, axes = plt.subplots(2, 2, figsize=(16, 14))
            fig.patch.set_facecolor(style['bg_color'])
            for ax in axes.flat:
                ax.set_facecolor(style['bg_color'])

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
            for r in seq_names:
                cluster = seq_to_cluster[r]
                enrich = cluster_to_enrichment.get(cluster, 'mixed')
                point_colors_enrich.append(enrichment_colors.get(enrich, '#CCCCCC'))
            axes[0, 1].scatter(embedding[:, 0], embedding[:, 1], c=point_colors_enrich, s=15, alpha=0.6)
            axes[0, 1].set_title("UMAP - Colored by Cluster Enrichment")
            axes[0, 1].set_xlabel("UMAP 1")
            axes[0, 1].set_ylabel("UMAP 2")
            enrich_patches = [Patch(facecolor=c, label=e) for e, c in enrichment_colors.items()]
            axes[0, 1].legend(handles=enrich_patches, loc='upper right')

            # 3. Bottom-left: Colored by sample with cluster number labels
            point_colors_sample = [sample_colors.get(s, '#999999') for s in sample_list]
            axes[1, 0].scatter(embedding[:, 0], embedding[:, 1], c=point_colors_sample, s=15, alpha=0.6)
            axes[1, 0].set_title("UMAP - Colored by Sample")
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
                                   bbox=dict(boxstyle='round,pad=0.2', facecolor=style['annotation_bg'], alpha=0.7, edgecolor=style['annotation_edge']))

            # Add sample legend to bottom-left panel
            sample_patches = [Patch(facecolor=sample_colors[s], label=s) for s in sorted(sample_colors.keys())]
            axes[1, 0].legend(handles=sample_patches, loc='upper right', fontsize=7)

            # 4. Bottom-right: Colored by enrichment with cluster labels
            axes[1, 1].scatter(embedding[:, 0], embedding[:, 1], c=point_colors_enrich, s=10, alpha=0.3)
            axes[1, 1].set_title("UMAP - Cluster Labels")
            axes[1, 1].set_xlabel("UMAP 1")
            axes[1, 1].set_ylabel("UMAP 2")

            # Add cluster labels with enrichment info
            for cluster_id in sorted(set(cluster_list)):
                mask = [c == cluster_id for c in cluster_list]
                cluster_points = embedding[mask]
                centroid_x = np.mean(cluster_points[:, 0])
                centroid_y = np.mean(cluster_points[:, 1])
                enrich = cluster_to_enrichment.get(cluster_id, 'mixed')
                enrich_short = enrich.replace('-enriched', '').replace('mixed', 'M')[:4]
                label = f"C{cluster_id}\n({enrich_short})"
                axes[1, 1].annotate(label, (centroid_x, centroid_y),
                                   fontsize=7, ha='center', va='center',
                                   bbox=dict(boxstyle='round,pad=0.2', facecolor=enrichment_colors.get(enrich, '#CCCCCC'),
                                            alpha=0.8, edgecolor=style['annotation_edge']))

            # Add enrichment legend to bottom-right panel
            axes[1, 1].legend(handles=enrich_patches, loc='upper right')

            plt.tight_layout()
            umap_file = f"{args.output_prefix}{suffix}.umap.pdf"
            plt.savefig(umap_file, dpi=150, bbox_inches='tight', facecolor=style['bg_color'])
            plt.close()
            print(f"  Saved UMAP plot: {umap_file}")

        # Save UMAP coordinates
        umap_coords = pd.DataFrame({
            'read': seq_names,
            'umap_1': embedding[:, 0],
            'umap_2': embedding[:, 1],
            'sample': sample_list,
            'cluster': cluster_labels,
            'enrichment': [cluster_to_enrichment.get(c, 'mixed') for c in cluster_labels]
        })
        umap_coords_file = f"{args.output_prefix}.umap_coordinates.tsv"
        umap_coords.to_csv(umap_coords_file, sep='\t', index=False)
        print(f"  Saved UMAP coordinates: {umap_coords_file}")

        # Try to generate interactive Plotly version (if enabled)
        if args.umap_html:
            try:
                import plotly.graph_objects as go
                from plotly.subplots import make_subplots

                # Create figure with dropdown menu for different colorings
                fig = go.Figure()

                # Prepare data arrays
                embedding_x = embedding[:, 0]
                embedding_y = embedding[:, 1]
                enrichment_list = [cluster_to_enrichment.get(c, 'mixed') for c in cluster_list]

                # Hover text (same for all views)
                hover_text = [f"Read: {r[:20]}...<br>Group: {g}<br>Cluster: {c}<br>Enrichment: {e}"
                             for r, g, c, e in zip(seq_names, group_list, cluster_list, enrichment_list)]

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
                    enrich = cluster_to_enrichment.get(cluster_id, 'mixed')
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
                                           {"title": f"UMAP - Colored by Group ({len(seq_names):,} reads)"}]),
                                dict(label="Color by Cluster",
                                     method="update",
                                     args=[{"visible": vis_cluster},
                                           {"title": f"UMAP - Colored by Cluster ({len(seq_names):,} reads)"}]),
                                dict(label="Color by Enrichment",
                                     method="update",
                                     args=[{"visible": vis_enrichment},
                                           {"title": f"UMAP - Colored by Enrichment ({len(seq_names):,} reads)"}]),
                            ],
                            direction="down",
                            showactive=True,
                            x=0.0,
                            xanchor="left",
                            y=1.15,
                            yanchor="top"
                        )
                    ],
                    title=f"UMAP - Colored by Group ({len(seq_names):,} reads)",
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

# --- Bubble Plot for Enrichment (all modes) ---
print(f"\n--- Generating enrichment bubble plot ---")

cluster_ids = cluster_df['cluster_id'].values
n_clusters = len(cluster_ids)

# Determine rows based on comparison mode
if args.comparison_mode == 'per-sample':
    row_labels = sample_labels
    n_rows = len(sample_labels)

    # Extract p-values and odds ratios matrices for per-sample mode
    pval_matrix = np.zeros((n_rows, n_clusters))
    odds_matrix = np.zeros((n_rows, n_clusters))
    pct_matrix = np.zeros((n_rows, n_clusters))

    for j, cid in enumerate(cluster_ids):
        row = cluster_df[cluster_df['cluster_id'] == cid].iloc[0]
        for i, sample in enumerate(sample_labels):
            pval_matrix[i, j] = row.get(f'{sample}_pval', 1.0)
            odds_matrix[i, j] = row.get(f'{sample}_odds', 1.0)
            pct_matrix[i, j] = row.get(f'{sample}_pct', 0.0)
else:
    # For two-group and multi-group modes, use groups as rows
    row_labels = all_groups
    n_rows = len(all_groups)

    # Extract percentage matrix and compute expected percentages
    pct_matrix = np.zeros((n_rows, n_clusters))

    # Calculate expected percentage for each group (based on total reads)
    total_sequences = len(seq_names)
    expected_pcts = {}
    for g in all_groups:
        group_total = sum(1 for s in seq_to_sample.values() if sample_to_group.get(s) == g)
        expected_pcts[g] = (group_total / total_sequences * 100) if total_sequences > 0 else 0

    for j, cid in enumerate(cluster_ids):
        row = cluster_df[cluster_df['cluster_id'] == cid].iloc[0]
        for i, group in enumerate(all_groups):
            pct_matrix[i, j] = row.get(f'{group}_pct', 0.0)

    # For group modes, we use overall q-value and derive enrichment from pct vs expected
    # Create pseudo-odds based on observed/expected ratio
    odds_matrix = np.zeros((n_rows, n_clusters))
    pval_matrix = np.zeros((n_rows, n_clusters))

    for j, cid in enumerate(cluster_ids):
        row = cluster_df[cluster_df['cluster_id'] == cid].iloc[0]
        cluster_qval = row.get('q_value', 1.0)
        enrichment = row.get('enrichment', 'mixed')

        for i, group in enumerate(all_groups):
            observed_pct = pct_matrix[i, j]
            expected_pct = expected_pcts[group]

            # Compute odds ratio as observed/expected
            if expected_pct > 0:
                odds_matrix[i, j] = observed_pct / expected_pct
            else:
                odds_matrix[i, j] = 1.0

            # Assign p-value: significant only if this group is enriched/depleted
            if enrichment == f'{group}-enriched':
                pval_matrix[i, j] = cluster_qval  # Use q-value for this group
            elif enrichment != 'mixed' and observed_pct < expected_pct:
                # This group is depleted (another group is enriched)
                pval_matrix[i, j] = cluster_qval
            else:
                pval_matrix[i, j] = 1.0  # Not significant for this group

# Transform p-values to -log10 (capped at 10 for visualization)
with np.errstate(divide='ignore'):
    neg_log_p = -np.log10(pval_matrix)
neg_log_p = np.clip(neg_log_p, 0, 10)  # Cap at 10 for display

# Generate enrichment plots for each background mode
for bg_mode, suffix in get_backgrounds_to_generate():
    style = apply_plot_style(bg_mode)

    # === Plot 1: Bubble plot (enrichment significance and direction) ===
    # Scale figure size for approximately square cells
    cell_size = 0.5  # inches per cell
    fig_width = max(10, n_clusters * cell_size + 4)
    fig_height = max(4, n_rows * cell_size + 2)
    fig1, ax1 = plt.subplots(figsize=(fig_width, fig_height))
    fig1.patch.set_facecolor(style['bg_color'])
    ax1.set_facecolor(style['bg_color'])

    # Create bubble plot
    for i, label in enumerate(row_labels):
        for j, cid in enumerate(cluster_ids):
            size = neg_log_p[i, j] * 50 + 10  # Scale size
            odds = odds_matrix[i, j]

            # Color by odds ratio: red = enriched (>1), blue = depleted (<1)
            if odds > 1:
                # Enriched: red scale
                intensity = min(1.0, np.log2(odds) / 3)  # Log scale, cap at 8x enrichment
                color = plt.cm.Reds(0.3 + intensity * 0.7)
            else:
                # Depleted: blue scale
                intensity = min(1.0, -np.log2(odds + 0.001) / 3)
                color = plt.cm.Blues(0.3 + intensity * 0.7)

            # Add significance indicator
            if pval_matrix[i, j] < 0.05:
                edgecolor = style['edge_color']
                linewidth = 1.5
            else:
                edgecolor = 'gray'
                linewidth = 0.5

            ax1.scatter(j, i, s=size, c=[color], edgecolors=edgecolor, linewidths=linewidth)

    ax1.set_xticks(range(n_clusters))
    ax1.set_xticklabels([str(cid) for cid in cluster_ids], rotation=90, fontsize=9)
    ax1.set_yticks(range(n_rows))
    ax1.set_yticklabels(row_labels, fontsize=11)
    ax1.set_xlabel('Cluster ID', fontsize=12)
    y_label = 'Sample' if args.comparison_mode == 'per-sample' else 'Group'
    ax1.set_ylabel(y_label, fontsize=12)
    ax1.set_title('Enrichment Bubble Plot\n(size: -log10(p-value), color: red=enriched, blue=depleted)', fontsize=13)
    ax1.set_xlim(-0.5, n_clusters - 0.5)
    ax1.set_ylim(-0.5, n_rows - 0.5)
    ax1.grid(True, alpha=0.3, color=style['grid_color'])

    # Add legend for size
    size_legend_vals = [1.3, 2, 3]  # -log10(p) values = p of 0.05, 0.01, 0.001
    for val in size_legend_vals:
        ax1.scatter([], [], s=val * 50 + 10, c='gray', alpha=0.7,
                    label=f'p = {10**(-val):.3g}')
    ax1.legend(loc='upper left', bbox_to_anchor=(1.01, 1), title='P-value', fontsize=9)

    plt.tight_layout()
    bubble_file = f"{args.output_prefix}{suffix}.enrichment_bubble.pdf"
    plt.savefig(bubble_file, dpi=150, bbox_inches='tight', facecolor=style['bg_color'])
    plt.close()
    print(f"  Saved bubble plot: {bubble_file}")

    # === Plot 2: Heatmap of sample/group percentage per cluster ===
    # Scale figure size for approximately square cells (same as bubble plot)
    fig2, ax2 = plt.subplots(figsize=(fig_width, fig_height))
    fig2.patch.set_facecolor(style['bg_color'])
    ax2.set_facecolor(style['bg_color'])

    im = ax2.imshow(pct_matrix, aspect='equal', cmap='YlOrRd', vmin=0, vmax=100)

    # Add text annotations
    for i in range(n_rows):
        for j in range(n_clusters):
            pct = pct_matrix[i, j]
            pval = pval_matrix[i, j]
            text_color = 'white' if pct > 50 else 'black'

            # Add asterisks for significance
            sig = ''
            if pval < 0.001:
                sig = '***'
            elif pval < 0.01:
                sig = '**'
            elif pval < 0.05:
                sig = '*'

            ax2.text(j, i, f'{pct:.0f}{sig}', ha='center', va='center',
                     fontsize=8, color=text_color, fontweight='bold' if sig else 'normal')

    ax2.set_xticks(range(n_clusters))
    ax2.set_xticklabels([str(cid) for cid in cluster_ids], rotation=90, fontsize=9)
    ax2.set_yticks(range(n_rows))
    ax2.set_yticklabels(row_labels, fontsize=11)
    ax2.set_xlabel('Cluster ID', fontsize=12)
    ax2.set_ylabel(y_label, fontsize=12)
    pct_title = 'Sample' if args.comparison_mode == 'per-sample' else 'Group'
    ax2.set_title(f'{pct_title} Percentage per Cluster\n(% of cluster reads from each {pct_title.lower()}; * p<0.05, ** p<0.01, *** p<0.001)', fontsize=13)

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax2, shrink=0.8)
    cbar.set_label('Percentage (%)', fontsize=11)

    plt.tight_layout()
    pct_file = f"{args.output_prefix}{suffix}.{pct_title.lower()}_percentage.pdf"
    plt.savefig(pct_file, dpi=150, bbox_inches='tight', facecolor=style['bg_color'])
    plt.close()
    print(f"  Saved percentage heatmap: {pct_file}")

# --- Summary ---
print(f"\n" + "=" * 60)
print("Summary")
print("=" * 60)
print(f"Total sequences: {len(seq_names):,}")
print(f"Number of clusters: {n_clusters}")
print(f"Valid clusters (size >= {args.min_cluster_size}): {len(cluster_df)}")
if args.comparison_mode == 'per-sample':
    for s in sample_labels:
        enriched_label = f'{s}-enriched'
        count = sum(cluster_df['enrichment'] == enriched_label)
        print(f"  - {enriched_label}: {count}")
else:
    for g in all_groups:
        enriched_label = f'{g}-enriched'
        count = sum(cluster_df['enrichment'] == enriched_label)
        print(f"  - {enriched_label}: {count}")
print(f"  - mixed: {sum(cluster_df['enrichment'] == 'mixed')}")
print(f"\nOutput files:")
print(f"  - {analysis_file}")
print(f"  - {assignments_file}")
print(f"  - {matrix_file}")
print(f"  - {plot_file}")
if args.log_file:
    print(f"  - {log_path}")
print(f"\nNext step: Use sequence_assignments.tsv with KaryoScope_cluster_plot.py")
print(f"to visualize reads from each cluster (sorted by centroid distance).")
print(f"Use --feature-matrix with the .feature_matrix.npz file for dendrogram header.")

# Print parameters table
print(f"\n" + "=" * 60)
print(f"Parameters")
print(f"=" * 60)
print(f"{'Parameter':<25} {'Value':<35}")
print(f"{'-'*25} {'-'*35}")

def fmt_param(name, value, attr_name=None):
    """Format parameter with (default) indicator if using default value."""
    if attr_name is None:
        attr_name = name.replace('-', '_')
    default_val = parser.get_default(attr_name)
    is_default = getattr(args, attr_name) == default_val
    suffix = " (default)" if is_default else ""
    return f"{name:<25} {value}{suffix}"

print(f"{'bed':<25} {len(args.bed)} file(s)")
print(f"{'output-prefix':<25} {args.output_prefix}")
print(fmt_param('sample-metadata', args.sample_metadata))
print(fmt_param('comparison-mode', args.comparison_mode))
print(fmt_param('control-group', args.control_group))
print(fmt_param('n-clusters', args.n_clusters))
print(fmt_param('min-k', args.min_k))
print(fmt_param('max-k', args.max_k))
print(fmt_param('k-selection', args.k_selection))
print(fmt_param('min-cluster-size', args.min_cluster_size))
print(fmt_param('min-sequence-length', args.min_sequence_length))
print(fmt_param('max-sequence-length', args.max_sequence_length))
print(fmt_param('exclude-features', args.exclude_features))
print(fmt_param('linkage-method', args.linkage_method))
print(fmt_param('matrix-type', args.matrix_type))
print(fmt_param('include-edges', args.include_edges, 'include_edges'))
print(fmt_param('edge-mode', args.edge_mode, 'edge_mode'))
print(fmt_param('matrix-mode', args.matrix_mode))
print(fmt_param('abundance', args.include_abundance, 'include_abundance'))
print(fmt_param('reduce-dims', args.reduce_dims))
print(fmt_param('umap', args.plot_umap, 'plot_umap'))
print(fmt_param('umap-neighbors', args.umap_neighbors))
print(fmt_param('umap-min-dist', args.umap_min_dist))
print(fmt_param('circular-dendrogram', args.plot_circular_dendrogram, 'plot_circular_dendrogram'))
print(fmt_param('perfect-threshold', args.perfect_threshold))
print(fmt_param('strong-threshold', args.strong_threshold))
print(fmt_param('early-stopping', args.early_stopping))
print(fmt_param('sequence-list', args.sequence_list))
print(fmt_param('silhouette-sample-size', args.silhouette_sample_size))
print(fmt_param('fdr-threshold', args.fdr_threshold))
print(fmt_param('fdr-method', args.fdr_method))
print(fmt_param('analysis-mode', args.analysis_mode))
print(fmt_param('structural-threshold', args.structural_threshold))
print(fmt_param('umap-html', args.umap_html))
print(fmt_param('background', args.background))
print(fmt_param('log-file', args.log_file))

# Print command used
print(f"\n" + "=" * 60)
print(f"Command")
print(f"=" * 60)
print(_original_command)
