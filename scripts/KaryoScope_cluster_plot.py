# KaryoScope Cluster Representative Plotting
# Plots representative reads from each cluster with sample and cluster annotations
#
# Usage:
# python KaryoScope_cluster_plot.py \
#   --featuresets chromosome,subtelomeric,region,acrocentric,repeat \
#   --bed-dirs /path/to/pre/bed /path/to/post/bed \
#   --colors-dir /path/to/colors \
#   --database KS_human_CHM13 \
#   --representatives /path/to/representative_reads.tsv \
#   --output cluster_representatives.svg \
#   --background black
#
# Compact mode (fewer reps, with dendrogram header):
# python KaryoScope_cluster_plot.py \
#   --featuresets region,subtelomeric \
#   --bed-dirs /path/to/pre/bed /path/to/post/bed \
#   --colors-dir /path/to/colors \
#   --database KS_human_CHM13 \
#   --representatives /path/to/representative_reads.tsv \
#   --feature-matrix /path/to/feature_matrix.npz \
#   --output cluster_representatives.compact.svg \
#   --compact --max-reps-per-cluster 3

from math import floor
import drawsvg as draw
import sys
import gzip
import argparse
from collections import defaultdict, OrderedDict
import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

# --- Parse command line arguments ---
parser = argparse.ArgumentParser(
    description="Generate KaryoScope SVG for cluster representative reads.",
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument("--featuresets", required=True,
                    help="Comma-separated list of feature sets to plot")
parser.add_argument("--bed-dirs", dest="bed_dirs", required=True, nargs='+',
                    help="Directories containing BED files (one per sample)")
parser.add_argument("--colors-dir", dest="colors_dir", required=True,
                    help="Directory containing the color files")
parser.add_argument("--database", required=True,
                    help="Database name (e.g., KS_human_CHM13)")
parser.add_argument("--representatives", required=True,
                    help="TSV file with representative reads (from cluster_analysis.py)")
parser.add_argument("--output", required=True,
                    help="Output SVG file path")
parser.add_argument("--background", dest="background_color", default="black",
                    choices=["white", "black"],
                    help="Background color for the SVG (default: black)")
parser.add_argument("--bar-width", dest="bar_width", type=int, default=8,
                    help="Width of each feature bar in pixels (default: 8)")
parser.add_argument("--bar-spacing", dest="bar_spacing", type=int, default=1,
                    help="Spacing between bars within a read group (default: 1)")
parser.add_argument("--read-spacing", dest="read_spacing", type=int, default=12,
                    help="Spacing between read groups (default: 12)")
parser.add_argument("--cluster-spacing", dest="cluster_spacing", type=int, default=30,
                    help="Spacing between clusters (default: 30)")
parser.add_argument("--ratio", type=float, default=1/300,
                    help="Ratio for scaling bp to pixels (default: 1/300)")
parser.add_argument("--smoothness", default="presmoothed",
                    help="Smoothness level (default: presmoothed)")
parser.add_argument("--enrichment-filter", dest="enrichment_filter", default=None,
                    choices=["Pre-enriched", "Post-enriched", "Mixed", None],
                    help="Only plot clusters with this enrichment type")
parser.add_argument("--max-clusters", dest="max_clusters", type=int, default=None,
                    help="Maximum number of clusters to plot")
parser.add_argument("--top-clusters", dest="top_clusters", default=None,
                    help="Select top N clusters per enrichment type, format: 'post:5,pre:4,mixed:3'")
parser.add_argument("--hide-brackets", dest="hide_brackets", action="store_true",
                    help="Hide cluster brackets and labels (cleaner dendrogram view)")
parser.add_argument("--feature-matrix", dest="feature_matrix_file", default=None,
                    help="NPZ file with feature matrix (from cluster_analysis.py) for dendrogram")
parser.add_argument("--compact", action="store_true",
                    help="Compact mode: fewer representatives per cluster")
parser.add_argument("--max-reps-per-cluster", dest="max_reps", type=int, default=None,
                    help="Maximum representatives per cluster (default: all, or 3 in compact mode)")
parser.add_argument("--sample-metadata", dest="sample_metadata", default=None,
                    help="TSV file with sample metadata (columns: sample, group, color).\n"
                         "If not provided, looks for .sample_metadata.tsv from cluster_analysis.py")

args = parser.parse_args()

# Set colors
BACKGROUND_COLOR = args.background_color
text_color = "#000000" if BACKGROUND_COLOR == "white" else "#FFFFFF"

featuresets = [f.strip() for f in args.featuresets.split(",")]
num_featuresets = len(featuresets)

# Handle compact mode defaults
if args.compact:
    max_reps = args.max_reps if args.max_reps else 3
    print(f"Compact mode: max {max_reps} representatives per cluster")
else:
    max_reps = args.max_reps  # None means all

# Load feature matrix if provided (for computing subset dendrogram)
feature_matrix_data = None
if args.feature_matrix_file and os.path.exists(args.feature_matrix_file):
    try:
        feature_matrix_data = np.load(args.feature_matrix_file, allow_pickle=True)
        print(f"Loaded feature matrix from: {args.feature_matrix_file}")
    except Exception as e:
        print(f"Warning: Could not load feature matrix: {e}")

print(f"Feature sets to plot: {featuresets}")
print(f"Background color: {BACKGROUND_COLOR}")

# --- Load representative reads ---
print(f"\nLoading representative reads from: {args.representatives}")
reps_df = pd.read_csv(args.representatives, sep='\t')
print(f"  Total representative reads: {len(reps_df)}")

# Handle --top-clusters selection (e.g., "treatment:5,control:4,mixed:3" or "post:5,pre:4,mixed:3")
if args.top_clusters:
    top_config = {}
    for part in args.top_clusters.split(','):
        key, val = part.strip().split(':')
        key = key.strip().lower()
        if key == 'mixed':
            top_config['Mixed'] = int(val)
        else:
            # Try to match against existing enrichment labels in the data
            # Support both old format (post/pre) and new format (group_name-enriched)
            enrichment_labels = reps_df['enrichment'].unique()
            matched = False
            for label in enrichment_labels:
                # Check for exact match (case-insensitive)
                if label.lower() == f'{key}-enriched':
                    top_config[label] = int(val)
                    matched = True
                    break
                # Backward compatibility: map 'post' -> '*Post*-enriched' patterns
                elif key == 'post' and 'post' in label.lower():
                    top_config[label] = int(val)
                    matched = True
                    break
                elif key == 'pre' and 'pre' in label.lower():
                    top_config[label] = int(val)
                    matched = True
                    break
            # If no match found, use the key directly with -enriched suffix
            if not matched:
                top_config[f'{key}-enriched'] = int(val)

    # Select top N clusters from each enrichment type
    selected_clusters = []
    for enrich_type, n_clusters in top_config.items():
        type_df = reps_df[reps_df['enrichment'] == enrich_type]
        type_clusters = type_df['cluster_id'].unique()[:n_clusters]
        selected_clusters.extend(type_clusters)
        print(f"  Selected top {len(type_clusters)} {enrich_type} clusters")

    reps_df = reps_df[reps_df['cluster_id'].isin(selected_clusters)]
    print(f"  Total after --top-clusters selection: {len(reps_df)}")

# Filter by enrichment if specified (single type)
elif args.enrichment_filter:
    reps_df = reps_df[reps_df['enrichment'] == args.enrichment_filter]
    print(f"  After filtering for {args.enrichment_filter}: {len(reps_df)}")

# Get unique clusters in order
clusters = reps_df['cluster_id'].unique()
if args.max_clusters:
    clusters = clusters[:args.max_clusters]
print(f"  Clusters to plot: {len(clusters)}")

# Group reads by cluster (with max_reps limit)
cluster_reads = OrderedDict()
for cluster_id in clusters:
    cluster_data = reps_df[reps_df['cluster_id'] == cluster_id]
    enrichment = cluster_data['enrichment'].iloc[0]
    reads = list(zip(cluster_data['read'], cluster_data['sample']))
    # Apply max_reps limit if specified
    if max_reps is not None and len(reads) > max_reps:
        reads = reads[:max_reps]
    cluster_reads[cluster_id] = {
        'enrichment': enrichment,
        'reads': reads
    }

if max_reps is not None:
    print(f"  Limited to max {max_reps} representatives per cluster")

# --- Find BED files for each sample ---
print(f"\nSearching for BED files in: {args.bed_dirs}")
sample_bed_paths = {}

for bed_dir in args.bed_dirs:
    # Look for BED files in the directory
    for f in os.listdir(bed_dir):
        if f.endswith('.features.bed.gz'):
            # Extract sample name
            parts = f.split('.')
            if len(parts) >= 2:
                sample_name = parts[0]
                if sample_name not in sample_bed_paths:
                    sample_bed_paths[sample_name] = bed_dir
                    print(f"  Found sample: {sample_name} -> {bed_dir}")

# --- Load color mappings ---
print(f"\nLoading color files...")
colors_files = {}
featureset_colors = {}
featureset_color_order = {}

for fs in featuresets:
    colors_pattern = f"{args.database}.{fs}.colors.txt"
    colors_path = os.path.join(args.colors_dir, colors_pattern)

    if not os.path.exists(colors_path):
        sys.stderr.write(f"Error: Colors file not found: {colors_path}\n")
        sys.exit(1)

    colors_files[fs] = colors_path
    featureset_colors[fs] = {0: ("#ffffff", 1.0)}
    featureset_color_order[fs] = []

    with open(colors_path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                feature = parts[0]
                color = parts[1]
                featureset_colors[fs][feature] = (color, 1.0)
                featureset_color_order[fs].append(feature)

    print(f"  {fs}: {len(featureset_color_order[fs])} colors")

# --- Load BED data for representative reads ---
print(f"\nLoading BED data for representative reads...")

# Get all reads we need
all_reads_needed = set()
read_to_sample = {}
for cluster_id, data in cluster_reads.items():
    for read, sample in data['reads']:
        all_reads_needed.add(read)
        read_to_sample[read] = sample

print(f"  Reads to load: {len(all_reads_needed)}")

# Load BED data for each featureset
read_data = defaultdict(lambda: defaultdict(list))  # read_data[read][featureset] = list of features

for fs in featuresets:
    for sample_name, bed_dir in sample_bed_paths.items():
        bed_pattern = f"{sample_name}.telogator.1.{args.database}.{fs}.{args.smoothness}.features.bed.gz"
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

                if scaffold in all_reads_needed:
                    read_data[scaffold][fs].append({
                        'start': start,
                        'stop': stop,
                        'feature': feature
                    })

print(f"  Loaded data for {len(read_data)} reads")

# --- Calculate scaffold lengths and normalize positions ---
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

# --- Calculate positions for drawing ---
ratio = args.ratio
group_width = (args.bar_width * num_featuresets) + (args.bar_spacing * (num_featuresets - 1))
left_margin = 150  # Space for cluster labels
dendrogram_height = 100 if feature_matrix_data is not None else 0
bracket_height = 0 if args.hide_brackets else 50
top_margin = 100 + dendrogram_height + bracket_height

# If feature matrix provided, compute dendrogram order and reorder reads
dendrogram_computed = False
dendro_result = None
read_to_original_cluster = {}  # Track original cluster assignments
read_to_original_enrichment = {}  # Track original enrichment

if feature_matrix_data is not None:
    try:
        from scipy.cluster.hierarchy import linkage, dendrogram as scipy_dendro

        # Get all displayed reads and track their original cluster
        all_displayed_reads = []
        read_to_sample_map = {}
        for cluster_id, data in cluster_reads.items():
            for read, sample in data['reads']:
                all_displayed_reads.append(read)
                read_to_sample_map[read] = sample
                read_to_original_cluster[read] = cluster_id
                read_to_original_enrichment[read] = data['enrichment']

        if len(all_displayed_reads) > 2:
            full_matrix = feature_matrix_data['adj_matrix']
            full_read_names = list(feature_matrix_data['read_names'])
            linkage_method = str(feature_matrix_data.get('linkage_method', 'ward'))

            read_to_idx = {r: i for i, r in enumerate(full_read_names)}
            subset_indices = [read_to_idx[r] for r in all_displayed_reads if r in read_to_idx]
            subset_reads = [all_displayed_reads[i] for i in range(len(subset_indices))]

            if len(subset_indices) > 2:
                subset_matrix = full_matrix[subset_indices]
                subset_linkage = linkage(subset_matrix, method=linkage_method)
                dendro_result = scipy_dendro(subset_linkage, no_plot=True)
                leaf_order = dendro_result['leaves']

                # Reorder reads according to dendrogram
                reordered_reads = [subset_reads[i] for i in leaf_order]

                # Rebuild cluster_reads with reordered reads
                cluster_reads_reordered = OrderedDict()
                cluster_reads_reordered['all'] = {
                    'enrichment': 'Mixed',
                    'reads': [(r, read_to_sample_map[r]) for r in reordered_reads]
                }
                cluster_reads = cluster_reads_reordered
                dendrogram_computed = True
                print(f"  Reordered {len(reordered_reads)} reads according to dendrogram")

    except Exception as e:
        print(f"  Warning: Could not compute dendrogram order: {e}")

# If not reordered, still build the mapping
if not read_to_original_cluster:
    for cluster_id, data in cluster_reads.items():
        for read, sample in data['reads']:
            read_to_original_cluster[read] = cluster_id
            read_to_original_enrichment[read] = data['enrichment']

# Calculate x positions for each read
read_x_positions = {}
current_x = left_margin
cluster_x_start = {}
cluster_x_end = {}

for cluster_id, data in cluster_reads.items():
    cluster_x_start[cluster_id] = current_x
    for read, sample in data['reads']:
        read_x_positions[read] = current_x
        current_x += group_width + args.read_spacing
    cluster_x_end[cluster_id] = current_x - args.read_spacing
    current_x += args.cluster_spacing  # Extra space between clusters

# --- Featureset display names ---
fs_display_names = {
    "chromosome": "Chromosome",
    "subtelomeric": "Subtelomere",
    "region": "Satellite",
    "acrocentric": "Acrocentric",
    "repeat": "Interspersed repeat"
}

# --- Calculate drawing data ---
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

            start_y = top_margin + 40 + floor(final_start * ratio)
            stop_y = top_margin + 40 + floor(final_stop * ratio)

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

# --- Calculate image dimensions ---
max_stop_y = 0
for read in drawing_data:
    for fs in drawing_data[read]:
        for rect in drawing_data[read][fs]:
            max_stop_y = max(max_stop_y, rect["y"] + rect["height"])

image_width = current_x + 50
legend_bottom_margin = 350
image_height = max_stop_y + 50 + legend_bottom_margin

print(f"\nImage dimensions: {image_width} x {image_height}")

# --- Create the drawing ---
d = draw.Drawing(image_width, image_height)
d.append(draw.Rectangle(0, 0, image_width, image_height, fill=BACKGROUND_COLOR))

# --- Draw linear dendrogram header (if precomputed) ---
if dendrogram_computed and dendro_result is not None:
    try:
        # Get displayed reads in their NEW order (after reordering)
        displayed_reads = []
        for cluster_id, data in cluster_reads.items():
            for read, sample in data['reads']:
                if read in read_x_positions:
                    displayed_reads.append(read)

        # Since reads are now in dendrogram leaf order, leaf i is at position i
        # Dendrogram x for leaf i is: i * 10 + 5
        n_leaves = len(displayed_reads)

        # Map dendrogram x to pixel x (leaves are now in order)
        def get_pixel_x(dx):
            # dx in dendrogram units: 5, 15, 25, ... for leaves
            leaf_idx = (dx - 5) / 10
            if leaf_idx < 0:
                leaf_idx = 0
            if leaf_idx >= n_leaves:
                leaf_idx = n_leaves - 1

            # For exact leaf positions
            low_idx = int(leaf_idx)
            high_idx = min(low_idx + 1, n_leaves - 1)

            low_x = read_x_positions[displayed_reads[low_idx]] + group_width / 2
            high_x = read_x_positions[displayed_reads[high_idx]] + group_width / 2

            frac = leaf_idx - low_idx
            return low_x * (1 - frac) + high_x * frac

        # Dendrogram y coordinates
        dendro_base_y = top_margin + 20  # Bottom of dendrogram (leaf level)
        max_height = max([max(dc) for dc in dendro_result['dcoord']]) if dendro_result['dcoord'] else 1
        if max_height == 0:
            max_height = 1

        def get_pixel_y(dy):
            return dendro_base_y - (dy / max_height) * (dendrogram_height - 15)

        line_color = '#AAAAAA' if BACKGROUND_COLOR == 'black' else '#444444'

        # Draw each U-shaped branch
        for icoord, dcoord in zip(dendro_result['icoord'], dendro_result['dcoord']):
            x1 = get_pixel_x(icoord[0])
            x2 = get_pixel_x(icoord[3])
            y_bottom_left = get_pixel_y(dcoord[0])
            y_top = get_pixel_y(dcoord[1])
            y_bottom_right = get_pixel_y(dcoord[3])

            # Left vertical line
            d.append(draw.Line(x1, y_bottom_left, x1, y_top,
                              stroke=line_color, stroke_width=1.5))
            # Horizontal line
            d.append(draw.Line(x1, y_top, x2, y_top,
                              stroke=line_color, stroke_width=1.5))
            # Right vertical line
            d.append(draw.Line(x2, y_top, x2, y_bottom_right,
                              stroke=line_color, stroke_width=1.5))

        print(f"  Drew dendrogram for {n_leaves} reads ({len(dendro_result['icoord'])} branches)")

    except Exception as e:
        import traceback
        print(f"  Warning: Could not draw dendrogram: {e}")
        traceback.print_exc()

# --- Draw cluster brackets and labels (unless hidden) ---
enrichment_colors = {
    'Pre-enriched': '#377EB8',
    'Post-enriched': '#E41A1C',
    'Mixed': '#999999'
}

if not args.hide_brackets:
    for cluster_id, data in cluster_reads.items():
        x_start = cluster_x_start[cluster_id]
        x_end = cluster_x_end[cluster_id]
        enrichment = data['enrichment']
        color = enrichment_colors.get(enrichment, '#999999')

        # Draw bracket line at top
        bracket_y = top_margin - 10
        d.append(draw.Line(x_start, bracket_y, x_end, bracket_y,
                           stroke=color, stroke_width=3))
        d.append(draw.Line(x_start, bracket_y, x_start, bracket_y + 8,
                           stroke=color, stroke_width=2))
        d.append(draw.Line(x_end, bracket_y, x_end, bracket_y + 8,
                           stroke=color, stroke_width=2))

        # Cluster label
        label_x = (x_start + x_end) / 2
        d.append(draw.Text(
            f"Cluster {cluster_id}",
            font_size=10,
            x=label_x,
            y=bracket_y - 15,
            fill=color,
            font_family='sans-serif',
            text_anchor='middle',
            font_weight='bold'
        ))

        # Enrichment label
        d.append(draw.Text(
            enrichment,
            font_size=8,
            x=label_x,
            y=bracket_y - 3,
            fill=color,
            font_family='sans-serif',
            text_anchor='middle'
        ))

# --- Load sample metadata for colors ---
sample_colors = {}
sample_to_group = {}

# Try to load from provided metadata file or auto-discover
metadata_file = args.sample_metadata
if not metadata_file:
    # Try to auto-discover metadata file from representatives file path
    reps_base = os.path.splitext(args.representatives)[0]
    reps_base = reps_base.replace('.representative_reads', '')
    potential_meta = f"{reps_base}.sample_metadata.tsv"
    if os.path.exists(potential_meta):
        metadata_file = potential_meta
        print(f"  Auto-discovered sample metadata: {metadata_file}")

if metadata_file and os.path.exists(metadata_file):
    try:
        meta_df = pd.read_csv(metadata_file, sep='\t')
        for _, row in meta_df.iterrows():
            sample = row['sample']
            sample_to_group[sample] = row.get('group', sample)
            if 'color' in meta_df.columns and pd.notna(row.get('color')):
                sample_colors[sample] = row['color']
        print(f"  Loaded sample metadata: {len(meta_df)} samples")
    except Exception as e:
        print(f"  Warning: Could not load sample metadata: {e}")

# Auto-generate colors for samples without colors
if not sample_colors:
    # Generate colors for all samples
    unique_samples = sorted(set(sample for data in cluster_reads.values() for _, sample in data['reads']))
    n_samples = len(unique_samples)
    tab10 = matplotlib.colormaps.get_cmap('tab10')

    for i, sample in enumerate(unique_samples):
        if sample not in sample_colors:
            if n_samples == 2:
                # Two samples: blue for first (alphabetically), red for second
                sample_colors[sample] = '#377EB8' if i == 0 else '#E41A1C'
            else:
                # Multiple samples: use colormap
                sample_colors[sample] = mcolors.rgb2hex(tab10(i % 10))

# Generate cluster colors using a colormap
unique_clusters = sorted(set(read_to_original_cluster.values()))
cluster_cmap = matplotlib.colormaps.get_cmap('tab20')
cluster_colors = {}
for i, cid in enumerate(unique_clusters):
    color_idx = i % 20  # tab20 has 20 colors
    cluster_colors[cid] = mcolors.rgb2hex(cluster_cmap(color_idx))

for cluster_id, data in cluster_reads.items():
    for read, sample in data['reads']:
        if read not in read_x_positions:
            continue
        base_x = read_x_positions[read]
        sample_color = sample_colors.get(sample, '#999999')

        # Get original cluster for this read
        orig_cluster = read_to_original_cluster.get(read, 'unknown')
        orig_cluster_color = cluster_colors.get(orig_cluster, '#666666')

        # Draw cluster indicator bar (top)
        d.append(draw.Rectangle(
            base_x,
            top_margin + 12,
            group_width,
            8,
            fill=orig_cluster_color
        ))

        # Draw sample indicator bar (below cluster bar)
        d.append(draw.Rectangle(
            base_x,
            top_margin + 22,
            group_width,
            8,
            fill=sample_color
        ))

# --- Draw all rectangles ---
for read in drawing_data:
    for fs in featuresets:
        for rect in drawing_data[read][fs]:
            if rect["height"] > 0:
                d.append(draw.Rectangle(
                    rect["x"],
                    rect["y"],
                    args.bar_width,
                    rect["height"],
                    fill=rect["fill"],
                    fill_opacity=rect["fill_opacity"]
                ))

# --- Draw read labels (short IDs) ---
for cluster_id, data in cluster_reads.items():
    for read, sample in data['reads']:
        if read not in read_x_positions:
            continue
        base_x = read_x_positions[read]
        label_x = base_x + group_width / 2

        # Use first 8 chars of read ID
        short_id = read[:8]

        d.append(draw.Text(
            short_id,
            font_size=5,
            x=label_x,
            y=top_margin + 35,
            fill=text_color,
            font_family='monospace',
            transform=f"rotate(90 {label_x} {top_margin + 35})",
            text_anchor='end'
        ))

# --- Draw featureset labels below bars ---
# Find max y per read
read_max_y = {}
for read in drawing_data:
    max_y = 0
    for fs in drawing_data[read]:
        for rect in drawing_data[read][fs]:
            max_y = max(max_y, rect["y"] + rect["height"])
    read_max_y[read] = max_y

# Draw labels for first read in each cluster
for cluster_id, data in cluster_reads.items():
    first_read = data['reads'][0][0]
    if first_read not in read_x_positions:
        continue

    base_x = read_x_positions[first_read]
    label_base_y = max(read_max_y.values()) + 5

    for fs_idx, fs in enumerate(featuresets):
        x_offset = fs_idx * (args.bar_width + args.bar_spacing)
        display_name = fs_display_names.get(fs, fs)
        label_x = base_x + x_offset + args.bar_width / 2

        d.append(draw.Text(
            display_name,
            font_size=6,
            x=label_x,
            y=label_base_y,
            fill=text_color,
            font_family='sans-serif',
            text_anchor='start',
            dominant_baseline='middle',
            transform=f"rotate(90 {label_x} {label_base_y})"
        ))

# --- Draw sample legend at top ---
legend_y = 20
legend_x = left_margin

d.append(draw.Text(
    "Sample:",
    font_size=10,
    x=legend_x,
    y=legend_y,
    fill=text_color,
    font_family='sans-serif',
    font_weight='bold'
))

for i, (sample, color) in enumerate(sample_colors.items()):
    item_x = legend_x + 55 + i * 90
    d.append(draw.Rectangle(item_x, legend_y - 8, 12, 12, fill=color))
    # Shorten sample names for legend
    short_sample = sample.replace('SW26_', '')
    d.append(draw.Text(
        short_sample,
        font_size=9,
        x=item_x + 16,
        y=legend_y,
        fill=text_color,
        font_family='sans-serif'
    ))

# --- Draw cluster legend at top ---
legend_x2 = legend_x + 200
d.append(draw.Text(
    "Cluster:",
    font_size=10,
    x=legend_x2,
    y=legend_y,
    fill=text_color,
    font_family='sans-serif',
    font_weight='bold'
))

# Group clusters by enrichment for the legend
clusters_by_enrichment = {'Post-enriched': [], 'Pre-enriched': [], 'Mixed': []}
for read, cid in read_to_original_cluster.items():
    enrich = read_to_original_enrichment.get(read, 'Mixed')
    if cid not in clusters_by_enrichment[enrich]:
        clusters_by_enrichment[enrich].append(cid)

# Sort clusters within each enrichment group
for enrich in clusters_by_enrichment:
    clusters_by_enrichment[enrich] = sorted(set(clusters_by_enrichment[enrich]))

# Draw cluster legend items with enrichment indication
cluster_legend_x = legend_x2 + 55
for enrich_type in ['Post-enriched', 'Pre-enriched', 'Mixed']:
    for cid in clusters_by_enrichment[enrich_type]:
        color = cluster_colors.get(cid, '#666666')
        d.append(draw.Rectangle(cluster_legend_x, legend_y - 8, 12, 12, fill=color))
        # Add small enrichment indicator
        enrich_color = enrichment_colors.get(enrich_type, '#999999')
        d.append(draw.Rectangle(cluster_legend_x, legend_y + 5, 12, 3, fill=enrich_color))
        d.append(draw.Text(
            f"C{cid}",
            font_size=8,
            x=cluster_legend_x + 15,
            y=legend_y,
            fill=text_color,
            font_family='sans-serif'
        ))
        cluster_legend_x += 45

# --- Draw enrichment legend ---
legend_x3 = cluster_legend_x + 20
d.append(draw.Text(
    "Enrichment:",
    font_size=10,
    x=legend_x3,
    y=legend_y,
    fill=text_color,
    font_family='sans-serif',
    font_weight='bold'
))

for i, (enrich, color) in enumerate(enrichment_colors.items()):
    item_x = legend_x3 + 75 + i * 100
    d.append(draw.Rectangle(item_x, legend_y - 8, 12, 12, fill=color))
    # Shorten enrichment label
    short_enrich = enrich.replace('-enriched', '')
    d.append(draw.Text(
        short_enrich,
        font_size=9,
        x=item_x + 16,
        y=legend_y,
        fill=text_color,
        font_family='sans-serif'
    ))

# --- Draw color legends at the bottom ---
color_legend_y_start = max_stop_y + 130
color_legend_x = left_margin
color_box_size = 12
color_text_offset = 16
colors_per_column = 8
item_width = 160
row_height = 18

def get_featureset_width(fs):
    num_items = len(featureset_color_order[fs])
    num_cols = (num_items + colors_per_column - 1) // colors_per_column
    return max(num_cols * item_width, 100)

featureset_legend_x = {}
current_legend_x = color_legend_x
for fs in featuresets:
    featureset_legend_x[fs] = current_legend_x
    current_legend_x += get_featureset_width(fs) + 30

for fs_idx, fs in enumerate(featuresets):
    section_x = featureset_legend_x[fs]
    display_name = fs_display_names.get(fs, fs)

    d.append(draw.Text(
        display_name,
        font_size=10,
        x=section_x,
        y=color_legend_y_start,
        fill=text_color,
        font_family='sans-serif',
        font_weight='bold'
    ))

    for i, feature_name in enumerate(featureset_color_order[fs]):
        color, opacity = featureset_colors[fs].get(feature_name, ("#ffffff", 1.0))

        row = i % colors_per_column
        col = i // colors_per_column

        item_x = section_x + col * item_width
        item_y = color_legend_y_start + 22 + row * row_height

        d.append(draw.Rectangle(
            item_x,
            item_y - 9,
            color_box_size,
            color_box_size,
            fill=color,
            stroke=text_color,
            stroke_width=0.5
        ))

        d.append(draw.Text(
            feature_name,
            font_size=8,
            x=item_x + color_text_offset,
            y=item_y,
            fill=text_color,
            font_family='sans-serif'
        ))

# --- Save the SVG ---
d.save_svg(args.output)

# --- Report ---
for fs in featuresets:
    if uncolored_features[fs]:
        sys.stderr.write(f"Warning: {fs} - features not in colors file:\n")
        for feature in sorted(list(uncolored_features[fs])):
            sys.stderr.write(f"  - {feature}\n")

print(f"\n--- Summary ---")
print(f"Clusters plotted: {len(cluster_reads)}")
total_reads = sum(len(data['reads']) for data in cluster_reads.values())
print(f"Total reads plotted: {total_reads}")
print(f"\n✅ Saved to {args.output}")
