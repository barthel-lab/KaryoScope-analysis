#!/usr/bin/env python3
"""
KaryoScope Cluster Representative Plotting

Plots representative reads from each cluster with sample and cluster annotations.
Designed to work with outputs from KaryoScope_cluster_analysis.py.

Usage:
  python KaryoScope_cluster_plot.py \
    --cluster-analysis-prefix test_aligned \
    --bed /path/to/pre/sample.bed /path/to/post/sample.bed \
    --colors /path/to/KS_human_CHM13 \
    --featuresets chromosome,region,subtelomeric \
    --output cluster_representatives.svg

With top-clusters filtering (select top N per category):
  python KaryoScope_cluster_plot.py \
    --cluster-analysis-prefix test_aligned \
    --bed /path/to/pre/sample.bed /path/to/post/sample.bed \
    --colors /path/to/KS_human_CHM13 \
    --featuresets region,subtelomeric \
    --output cluster_representatives.svg \
    --top-clusters "post:4,pre:2,mixed:3"
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict, OrderedDict
from math import floor

import drawsvg as draw
import matplotlib
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd


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


def load_cluster_analysis(cluster_analysis_file):
    """Load cluster analysis results to get enrichment info and cluster order.

    Returns:
        tuple: (cluster_enrichments dict, cluster_order list)
            - cluster_enrichments: cluster_id -> enrichment label
            - cluster_order: list of cluster_ids sorted by enrichment tier then p-value:
                Tier 0: 100% enriched (perfect)
                Tier 1: 80%+ enriched (strong)
                Tier 2: all others
    """
    cluster_enrichments = {}
    cluster_order = []

    if cluster_analysis_file and os.path.exists(cluster_analysis_file):
        try:
            df = pd.read_csv(cluster_analysis_file, sep='\t')

            # Find percentage columns (they end with _pct)
            pct_cols = [c for c in df.columns if c.endswith('_pct')]

            # Determine enrichment tier for each cluster
            def get_enrichment_tier(row):
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
                cluster_enrichments[row['cluster_id']] = row['enrichment']
                cluster_order.append(row['cluster_id'])
            print(f"  Loaded cluster analysis: {len(df)} clusters")
        except Exception as e:
            print(f"  Warning: Could not load cluster analysis: {e}")

    return cluster_enrichments, cluster_order


def load_representative_reads(reps_file, cluster_enrichments=None, cluster_order=None, max_reps=None, top_clusters=None, max_clusters=None):
    """Load read assignments from TSV file.

    Args:
        reps_file: Path to read_assignments.tsv (all reads with cluster assignments and stats)
        cluster_enrichments: Dict of cluster_id -> enrichment label from cluster_analysis.tsv
        cluster_order: List of cluster_ids in priority order (100% enriched first, then 80%+, then by p-value)
        max_reps: Maximum representatives per cluster (selects by rank, closest to centroid first)
        top_clusters: Dict of {category: n_clusters} where category is a group name or 'mixed'
        max_clusters: Maximum total clusters to include

    Returns:
        tuple: (cluster_reads OrderedDict, unique_enrichments set)
    """
    print(f"\nLoading read assignments from: {reps_file}")
    reps_df = pd.read_csv(reps_file, sep='\t')
    print(f"  Total reads: {len(reps_df)}")

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

    # Handle --top-clusters selection
    if top_clusters:
        selected_clusters = []
        for category, n_clusters in top_clusters.items():
            # Map category to enrichment label(s)
            # 'mixed' -> 'mixed'
            # 'post' -> 'post-enriched'
            # 'pre' -> 'pre-enriched'
            if category.lower() == 'mixed':
                target_enrichments = ['mixed']
            else:
                # Try to match group name to enrichment label
                target_enrichments = []
                for enrich in unique_enrichments:
                    # Match "post" to "post-enriched", case-insensitive
                    if enrich.lower() == f"{category.lower()}-enriched":
                        target_enrichments.append(enrich)
                    # Also allow direct match (e.g., if enrichment is just "post")
                    elif enrich.lower() == category.lower():
                        target_enrichments.append(enrich)

            for target_enrich in target_enrichments:
                type_df = reps_df[reps_df['enrichment'] == target_enrich]
                # Use cluster_order to get most significant clusters first
                if cluster_order:
                    type_clusters_set = set(type_df['cluster'].unique())
                    type_clusters = [c for c in cluster_order if c in type_clusters_set][:n_clusters]
                else:
                    type_clusters = list(type_df['cluster'].unique()[:n_clusters])
                selected_clusters.extend(type_clusters)
                print(f"  Selected top {len(type_clusters)} {target_enrich} clusters")

            if not target_enrichments:
                print(f"  Warning: No enrichment category matches '{category}'")

        reps_df = reps_df[reps_df['cluster'].isin(selected_clusters)]
        print(f"  Total after --top-clusters selection: {len(reps_df)}")

    # Get unique clusters - use priority order from cluster_analysis.tsv if available
    # (sorted by: 100% enriched first, then 80%+, then by p-value)
    available_clusters = set(reps_df['cluster'].unique())
    if cluster_order:
        clusters = [c for c in cluster_order if c in available_clusters]
    else:
        clusters = list(reps_df['cluster'].unique())

    if max_clusters:
        clusters = clusters[:max_clusters]
    print(f"  Clusters to plot: {len(clusters)}")

    # Group reads by cluster
    cluster_reads = OrderedDict()
    for cluster_id in clusters:
        cluster_data = reps_df[reps_df['cluster'] == cluster_id]
        enrichment = cluster_data['enrichment'].iloc[0]

        # Apply max_reps limit using rank (rank 1 = closest to centroid)
        # with proportional sampling by sample to maintain representation
        if max_reps is not None and len(cluster_data) > max_reps:
            # Count reads per sample in this cluster
            sample_counts = cluster_data['sample'].value_counts()
            total_reads = len(cluster_data)

            # Calculate proportional allocation for each sample
            selected_reads = []
            remaining_slots = max_reps

            # Sort samples by count (descending) for consistent ordering
            for sample in sample_counts.index:
                # Proportional allocation, at least 1 if sample has reads and slots remain
                proportion = sample_counts[sample] / total_reads
                n_select = max(1, round(proportion * max_reps))
                n_select = min(n_select, remaining_slots, sample_counts[sample])

                if n_select > 0:
                    # Select reads with lowest rank (closest to centroid) for this sample
                    sample_reads = cluster_data[cluster_data['sample'] == sample].sort_values('rank')
                    selected = list(zip(sample_reads['read'].iloc[:n_select],
                                       sample_reads['sample'].iloc[:n_select]))
                    selected_reads.extend(selected)
                    remaining_slots -= n_select

                if remaining_slots <= 0:
                    break

            reads = selected_reads
        else:
            reads = list(zip(cluster_data['read'], cluster_data['sample']))

        cluster_reads[cluster_id] = {
            'enrichment': enrichment,
            'reads': reads
        }

    if max_reps is not None:
        print(f"  Limited to max {max_reps} representatives per cluster (by centroid proximity, proportional sampling)")

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
                continue

            open_func = gzip.open if bed_path.endswith(".gz") else open
            mode = "rt" if bed_path.endswith(".gz") else "r"

            with open_func(bed_path, mode) as f:
                for line in f:
                    parts = line.strip().split()[:4]
                    scaffold, start, stop, feature = parts
                    start, stop = int(start), int(stop)

                    if scaffold in reads_needed:
                        read_data[scaffold][fs].append({
                            'start': start,
                            'stop': stop,
                            'feature': feature
                        })

    print(f"  Loaded data for {len(read_data)} reads")
    return read_data


# =============================================================================
# Helper Functions: Enrichment Handling
# =============================================================================

def get_enrichment_colors(group_colors, unique_enrichments):
    """Generate enrichment colors based on group colors.

    Args:
        group_colors: Dict of group -> color from metadata
        unique_enrichments: Set of enrichment labels from data

    Returns:
        dict: enrichment_label -> color
    """
    enrichment_colors = {'Mixed': '#999999'}

    for enrich in unique_enrichments:
        if enrich == 'Mixed':
            continue

        # Extract group name from enrichment label (e.g., "post-enriched" -> "post")
        group = enrich.replace('-enriched', '')

        # Try to find matching group color (case-insensitive)
        color_found = False
        for g, c in group_colors.items():
            if g.lower() == group.lower():
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

def compute_dendrogram_order(feature_matrix_data, cluster_reads):
    """Compute dendrogram ordering for reads using original pairwise distances.

    Extracts pairwise distances from the full feature matrix and performs
    hierarchical clustering with optimal leaf ordering to minimize crossings.

    Returns:
        tuple: (reordered_cluster_reads, dendro_data, read_to_original_cluster, read_to_original_enrichment)
               dendro_data contains 'linkage' and 'dist_matrix' for manual drawing
    """
    from scipy.cluster.hierarchy import linkage, leaves_list, optimal_leaf_ordering
    from scipy.spatial.distance import pdist, squareform

    read_to_original_cluster = {}
    read_to_original_enrichment = {}

    # Get all displayed reads
    all_displayed_reads = []
    read_to_sample_map = {}

    for cluster_id, data in cluster_reads.items():
        for read, sample in data['reads']:
            all_displayed_reads.append(read)
            read_to_sample_map[read] = sample
            read_to_original_cluster[read] = cluster_id
            read_to_original_enrichment[read] = data['enrichment']

    if len(all_displayed_reads) <= 2:
        return cluster_reads, None, read_to_original_cluster, read_to_original_enrichment

    try:
        full_matrix = feature_matrix_data['adj_matrix']
        full_read_names = list(feature_matrix_data['read_names'])

        read_to_idx = {r: i for i, r in enumerate(full_read_names)}
        subset_indices = [read_to_idx[r] for r in all_displayed_reads if r in read_to_idx]
        subset_reads = [r for r in all_displayed_reads if r in read_to_idx]

        if len(subset_indices) <= 2:
            return cluster_reads, None, read_to_original_cluster, read_to_original_enrichment

        # Compute pairwise distances on FULL matrix first
        full_distances = pdist(full_matrix, metric='euclidean')
        full_dist_square = squareform(full_distances)

        # Extract the subset distance matrix (preserves original distances)
        subset_dist_square = full_dist_square[np.ix_(subset_indices, subset_indices)]
        subset_distances = squareform(subset_dist_square)

        # Compute linkage on subset using 'average' method
        subset_linkage = linkage(subset_distances, method='average')

        # Apply optimal leaf ordering to rotate branches and minimize crossings
        # This reorders leaves so that adjacent leaves are most similar
        optimized_linkage = optimal_leaf_ordering(subset_linkage, subset_distances)

        # Get leaf order from optimized linkage
        leaf_order = leaves_list(optimized_linkage)

        # Reorder reads according to optimal ordering
        reordered_reads = [subset_reads[i] for i in leaf_order]

        # Store data for manual dendrogram drawing (use optimized linkage)
        dendro_data = {
            'linkage': optimized_linkage,
            'dist_matrix': subset_dist_square,
            'reads': reordered_reads,
            'leaf_order': leaf_order
        }

        # Rebuild cluster_reads with reordered reads
        cluster_reads_reordered = OrderedDict()
        cluster_reads_reordered['all'] = {
            'enrichment': 'Mixed',
            'reads': [(r, read_to_sample_map[r]) for r in reordered_reads]
        }

        print(f"  Reordered {len(reordered_reads)} reads according to dendrogram (optimal leaf ordering)")
        return cluster_reads_reordered, dendro_data, read_to_original_cluster, read_to_original_enrichment

    except Exception as e:
        print(f"  Warning: Could not compute dendrogram order: {e}")
        import traceback
        traceback.print_exc()
        return cluster_reads, None, read_to_original_cluster, read_to_original_enrichment


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


def draw_cluster_brackets(d, cluster_reads, cluster_x_start, cluster_x_end,
                          enrichment_colors, top_margin, text_color):
    """Draw cluster brackets and labels."""
    for cluster_id, data in cluster_reads.items():
        if cluster_id == 'all':  # Skip when in dendrogram mode
            continue

        x_start = cluster_x_start[cluster_id]
        x_end = cluster_x_end[cluster_id]
        enrichment = data['enrichment']
        color = enrichment_colors.get(enrichment, '#999999')

        bracket_y = top_margin - 10
        d.append(draw.Line(x_start, bracket_y, x_end, bracket_y, stroke=color, stroke_width=3))
        d.append(draw.Line(x_start, bracket_y, x_start, bracket_y + 8, stroke=color, stroke_width=2))
        d.append(draw.Line(x_end, bracket_y, x_end, bracket_y + 8, stroke=color, stroke_width=2))

        label_x = (x_start + x_end) / 2
        d.append(draw.Text(
            f"Cluster {cluster_id}",
            font_size=10, x=label_x, y=bracket_y - 15,
            fill=color, font_family='sans-serif',
            text_anchor='middle', font_weight='bold'
        ))

        d.append(draw.Text(
            enrichment,
            font_size=8, x=label_x, y=bracket_y - 3,
            fill=color, font_family='sans-serif', text_anchor='middle'
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
            orig_enrichment = read_to_original_enrichment.get(read, 'Mixed')
            enrichment_color = enrichment_colors.get(orig_enrichment, '#999999')

            # Cluster indicator bar (top, 8px height)
            d.append(draw.Rectangle(base_x, annot_start_y, group_width, 8, fill=orig_cluster_color))

            # Enrichment indicator bar (middle, thin 3px height)
            d.append(draw.Rectangle(base_x, annot_start_y + 10, group_width, 3, fill=enrichment_color))

            # Sample indicator bar (bottom, 8px height)
            d.append(draw.Rectangle(base_x, annot_start_y + 15, group_width, 8, fill=sample_color))


def draw_feature_bars(d, drawing_data, featuresets, bar_width, read_heights, num_featuresets):
    """Draw feature rectangles with borders between featuresets and around outer edge.

    Args:
        d: Drawing object
        drawing_data: Feature data per read
        featuresets: List of feature sets
        bar_width: Width of each bar
        read_heights: Dict of read -> (min_y, max_y, x_start, total_width) for borders
        num_featuresets: Number of feature sets
    """
    stroke_width = 0.5

    # First draw all feature rectangles (no individual borders)
    for read in drawing_data:
        for fs in featuresets:
            for rect in drawing_data[read][fs]:
                if rect["height"] > 0:
                    d.append(draw.Rectangle(
                        rect["x"], rect["y"],
                        bar_width, rect["height"],
                        fill=rect["fill"],
                        fill_opacity=rect["fill_opacity"]
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


def draw_read_labels(d, cluster_reads, read_x_positions, group_width, top_margin, text_color):
    """Draw read ID labels above the annotation bars, rotated 90 degrees."""
    for cluster_id, data in cluster_reads.items():
        for read, sample in data['reads']:
            if read not in read_x_positions:
                continue

            base_x = read_x_positions[read]
            label_x = base_x + group_width / 2
            label_y = top_margin - 5  # Position above annotation bars
            short_id = read[:8]

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

    for i, (sample, color) in enumerate(sample_colors.items()):
        item_x = legend_x + 60 + i * 100
        d.append(draw.Rectangle(item_x, current_y - 8, 12, 12, fill=color))
        d.append(draw.Text(
            sample, font_size=9, x=item_x + 16, y=current_y,
            fill=text_color, font_family='sans-serif'
        ))

    # --- Cluster legend (row 2) ---
    current_y += row_height
    d.append(draw.Text(
        "Cluster:", font_size=10, x=legend_x, y=current_y,
        fill=text_color, font_family='sans-serif', font_weight='bold'
    ))

    # Group clusters by enrichment
    clusters_by_enrichment = defaultdict(list)
    for read, cid in read_to_original_cluster.items():
        enrich = read_to_original_enrichment.get(read, 'Mixed')
        if cid not in clusters_by_enrichment[enrich]:
            clusters_by_enrichment[enrich].append(cid)

    # Sort clusters within each enrichment
    for enrich in clusters_by_enrichment:
        clusters_by_enrichment[enrich] = sorted(set(clusters_by_enrichment[enrich]))

    # Draw cluster items
    cluster_legend_x = legend_x + 60

    # Sort enrichment types for consistent ordering
    enrichment_order = sorted(clusters_by_enrichment.keys(),
                             key=lambda x: (x == 'Mixed', x))  # Mixed last

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


def draw_color_legends(d, featuresets, featureset_colors, featureset_color_order,
                       fs_display_names, color_legend_y_start, left_margin, text_color):
    """Draw featureset color legends at bottom."""
    color_box_size = 10
    color_text_offset = 14
    colors_per_column = 12  # More items per column = fewer columns
    item_width = 120  # Narrower columns
    row_height = 14  # Tighter row spacing

    def get_featureset_width(fs):
        num_items = len(featureset_color_order[fs])
        num_cols = (num_items + colors_per_column - 1) // colors_per_column
        return max(num_cols * item_width, 100)

    featureset_legend_x = {}
    current_legend_x = left_margin
    for fs in featuresets:
        featureset_legend_x[fs] = current_legend_x
        current_legend_x += get_featureset_width(fs) + 15  # Reduced spacing between sections

    for fs in featuresets:
        section_x = featureset_legend_x[fs]
        display_name = fs_display_names.get(fs, fs)

        d.append(draw.Text(
            display_name, font_size=9, x=section_x, y=color_legend_y_start,
            fill=text_color, font_family='sans-serif', font_weight='bold'
        ))

        for i, feature_name in enumerate(featureset_color_order[fs]):
            color, opacity = featureset_colors[fs].get(feature_name, ("#ffffff", 1.0))

            row = i % colors_per_column
            col = i // colors_per_column

            item_x = section_x + col * item_width
            item_y = color_legend_y_start + 18 + row * row_height

            d.append(draw.Rectangle(
                item_x, item_y - 7, color_box_size, color_box_size,
                fill=color, stroke=text_color, stroke_width=0.5
            ))

            d.append(draw.Text(
                feature_name, font_size=7, x=item_x + color_text_offset, y=item_y,
                fill=text_color, font_family='sans-serif'
            ))


# =============================================================================
# Main Script
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
    parser.add_argument("--bed", dest="bed_files", required=True, nargs='+',
                        help="Full paths to BED files (same files used in cluster_analysis.py)")
    parser.add_argument("--colors", dest="colors_dir", required=True,
                        help="Full path to colors database directory (contains {database}.{featureset}.colors.txt files)")
    parser.add_argument("--featuresets", required=True,
                        help="Comma-separated list of feature sets to plot")

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
    parser.add_argument("--smoothness", default="presmoothed",
                        help="Smoothness level (default: presmoothed)")

    # Filtering options
    parser.add_argument("--max-clusters", dest="max_clusters", type=int, default=None,
                        help="Maximum number of clusters to plot")
    parser.add_argument("--top-clusters", dest="top_clusters", default=None,
                        help="Select top N clusters per category, format: 'post:4,pre:2,mixed:3'. "
                             "Categories: group names from metadata (e.g., 'post', 'pre'), "
                             "'mixed' for clusters with no clear enrichment")

    # Mode options
    parser.add_argument("--hide-brackets", dest="hide_brackets", action="store_true",
                        help="Hide cluster brackets and labels (cleaner dendrogram view)")
    parser.add_argument("--max-reps-per-cluster", dest="max_reps", type=int, default=3,
                        help="Maximum representatives per cluster (default: 3)")

    return parser.parse_args()


def main():
    args = parse_args()

    # --- Auto-discover files from prefix ---
    prefix = args.cluster_prefix
    representatives_file = f"{prefix}.read_assignments.tsv"
    feature_matrix_file = f"{prefix}.feature_matrix.npz"
    sample_metadata_file = f"{prefix}.sample_metadata.tsv"
    cluster_analysis_file = f"{prefix}.cluster_analysis.tsv"

    # Verify required files exist
    if not os.path.exists(representatives_file):
        sys.stderr.write(f"Error: Read assignments file not found: {representatives_file}\n")
        sys.exit(1)

    # --- Parse BED paths to get sample directories and database ---
    sample_bed_paths, database = parse_bed_paths(args.bed_files)

    if not database:
        sys.stderr.write("Error: Could not determine database from BED file paths\n")
        sys.exit(1)

    # --- Setup ---
    background_color = args.background_color
    text_color = "#000000" if background_color == "white" else "#FFFFFF"
    featuresets = [f.strip() for f in args.featuresets.split(",")]
    num_featuresets = len(featuresets)

    # Featureset display names
    fs_display_names = {
        "chromosome": "Chromosome",
        "subtelomeric": "Subtelomere",
        "region": "Satellite",
        "acrocentric": "Acrocentric",
        "repeat": "Interspersed repeat"
    }

    max_reps = args.max_reps

    print(f"Feature sets to plot: {featuresets}")
    print(f"Background color: {background_color}")

    # --- Load data ---
    # Load sample metadata
    sample_to_group, sample_colors, group_colors = load_sample_metadata(sample_metadata_file)

    # Load cluster analysis to get enrichment info and cluster priority order
    cluster_enrichments, cluster_order = load_cluster_analysis(cluster_analysis_file)

    # Parse top_clusters if provided
    top_clusters = None
    if args.top_clusters:
        top_clusters = {}
        for part in args.top_clusters.split(','):
            key, val = part.strip().split(':')
            top_clusters[key.strip()] = int(val)

    # Load read assignments
    cluster_reads, unique_enrichments = load_representative_reads(
        representatives_file,
        cluster_enrichments=cluster_enrichments,
        cluster_order=cluster_order,
        max_reps=max_reps,
        top_clusters=top_clusters,
        max_clusters=args.max_clusters
    )

    # Load feature matrix
    feature_matrix_data = load_feature_matrix(feature_matrix_file)

    # Load color files
    featureset_colors, featureset_color_order = load_color_files(
        args.colors_dir, database, featuresets
    )

    # Get all reads we need
    all_reads_needed = set()
    read_to_sample = {}
    for cluster_id, data in cluster_reads.items():
        for read, sample in data['reads']:
            all_reads_needed.add(read)
            read_to_sample[read] = sample

    # Load BED data
    read_data = load_bed_data(
        sample_bed_paths, database, featuresets, args.smoothness, all_reads_needed
    )

    # --- Generate colors ---
    # Get all unique samples
    all_samples = sorted(set(sample for data in cluster_reads.values() for _, sample in data['reads']))
    sample_colors = generate_sample_colors(all_samples, sample_colors)

    # Generate enrichment colors from group colors
    enrichment_colors = get_enrichment_colors(group_colors, unique_enrichments)

    # --- Compute dendrogram order if feature matrix provided ---
    dendrogram_computed = False
    dendro_result = None
    read_to_original_cluster = {}
    read_to_original_enrichment = {}

    if feature_matrix_data is not None:
        cluster_reads, dendro_result, read_to_original_cluster, read_to_original_enrichment = \
            compute_dendrogram_order(feature_matrix_data, cluster_reads)
        dendrogram_computed = dendro_result is not None

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
    dendrogram_height = 100 if feature_matrix_data is not None else 0
    bracket_height = 0 if args.hide_brackets else 50
    top_margin = 100 + dendrogram_height + bracket_height

    # Calculate x positions
    read_x_positions = {}
    cluster_x_start = {}
    cluster_x_end = {}
    current_x = left_margin

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
    ratio = args.ratio
    drawing_data = defaultdict(lambda: defaultdict(list))
    uncolored_features = defaultdict(set)

    for read in read_data:
        if read not in read_x_positions:
            continue

        base_x = read_x_positions[read]
        scaffold_min_start = scaffold_min_starts.get(read, 0)

        for fs_idx, fs in enumerate(featuresets):
            x_offset = fs_idx * (args.bar_width + args.bar_spacing)

            for feat in read_data[read][fs]:
                final_start = feat['start'] - scaffold_min_start
                final_stop = feat['stop'] - scaffold_min_start

                # Feature bars start after annotation bars (annot at +22, height ~25, so start at +50)
                start_y = top_margin + 50 + floor(final_start * ratio)
                stop_y = top_margin + 50 + floor(final_stop * ratio)

                color, fill_opacity = featureset_colors[fs].get(feat['feature'], ("#ffffff", 1.0))
                if feat['feature'] not in featureset_colors[fs]:
                    uncolored_features[fs].add(feat['feature'])

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
    legend_bottom_margin = 350
    image_height = max_stop_y + 50 + legend_bottom_margin

    print(f"\nImage dimensions: {image_width} x {image_height}")

    # --- Create drawing ---
    d = draw.Drawing(image_width, image_height)
    d.append(draw.Rectangle(0, 0, image_width, image_height, fill=background_color))

    # --- Draw components ---
    # Dendrogram header
    if dendrogram_computed and dendro_result is not None:
        displayed_reads = []
        for cluster_id, data in cluster_reads.items():
            for read, sample in data['reads']:
                if read in read_x_positions:
                    displayed_reads.append(read)

        draw_dendrogram(d, dendro_result, read_x_positions, displayed_reads,
                       group_width, top_margin, dendrogram_height, background_color)

    # Cluster brackets
    if not args.hide_brackets:
        draw_cluster_brackets(d, cluster_reads, cluster_x_start, cluster_x_end,
                             enrichment_colors, top_margin, text_color)

    # Annotation bars
    draw_annotation_bars(d, cluster_reads, read_x_positions, read_to_original_cluster,
                        read_to_original_enrichment, sample_colors, cluster_colors,
                        enrichment_colors, group_width, top_margin, left_margin, text_color)

    # Feature bars
    draw_feature_bars(d, drawing_data, featuresets, args.bar_width, read_heights, num_featuresets)

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

    # Color legends at bottom
    color_legend_y_start = max_stop_y + 130
    draw_color_legends(d, featuresets, featureset_colors, featureset_color_order,
                      fs_display_names, color_legend_y_start, left_margin, text_color)

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


if __name__ == "__main__":
    main()
