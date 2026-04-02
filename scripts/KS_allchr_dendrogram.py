#!/usr/bin/env python3
"""
KS_allchr_dendrogram.py — All-chromosome dendrogram with feature bars

Builds a single dendrogram from ALL non-excluded chromosomes, using one
representative per cluster. Each row shows the centromere composition bar
for that representative. Chromosomes cluster together naturally because
their feature compositions are distinct.

The feature matrix is built globally (shared feature vocabulary, directional
edges, no abundance) and distances are computed globally — NOT reused from
per-chromosome clustering.

Usage:
  python3 scripts/KS_allchr_dendrogram.py \
    --assignments agent_results/allchr_structure.sequence_assignments.tsv \
    --bed /path/to/pangenome.ALLchr.centromere.KS_human_CHM13.presmoothed.region.pass.bed \
    --colors /path/to/KS_human_CHM13 \
    --output agent_results/allchr_dendrogram.svg
"""

import argparse
import gzip
import os
import subprocess
import sys
from collections import defaultdict, Counter

import drawsvg as draw
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import pdist
from sklearn.preprocessing import StandardScaler

# ─── Constants ────────────────────────────────────────────────────────────────

EXCLUDED_CHROMS = {'chr13', 'chr14', 'chr15', 'chr21', 'chr22', 'chrY'}

CHROM_ORDER = [
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
    'chr10', 'chr11', 'chr12', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chrX'
]

# Soft palette for chromosome background bands
CHROM_BAND_COLORS = [
    '#F0F4FF', '#FFF5F0', '#F0FFF4', '#FFF9F0', '#F5F0FF',
    '#F0FFFF', '#FFF0F5', '#F4FFF0', '#FFF0FF',
    '#F0F8FF', '#FFFFF0', '#FFF0F0', '#F0FFF8', '#F8F0FF',
    '#FAFFF0', '#F0FAFF', '#FFF5F5', '#F5FFF0'
]

# ─── Argument parsing ─────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="All-chromosome dendrogram with centromere feature bars")
    p.add_argument("--assignments", required=True,
                   help="sequence_assignments.tsv from structure mode")
    p.add_argument("--bed", required=True, nargs='+',
                   help="BED file(s) with feature annotations")
    p.add_argument("--colors", required=True,
                   help="Colors database directory (e.g. resources/databases/KS_human_CHM13)")
    p.add_argument("--output", required=True,
                   help="Output SVG path")
    p.add_argument("--featureset", default="region",
                   help="Featureset for color lookup (default: region)")
    p.add_argument("--exclude-features", dest="exclude_features",
                   default="novel",
                   help="Comma-separated features to exclude from matrix (default: novel)")
    p.add_argument("--edge-mode", dest="edge_mode", default="directional",
                   choices=["directional", "symmetric"],
                   help="Edge counting mode (default: directional)")
    p.add_argument("--matrix-type", dest="matrix_type",
                   default="count_log1p_zscore",
                   choices=["count", "count_log1p", "count_log1p_zscore"],
                   help="Matrix normalization (default: count_log1p_zscore)")
    p.add_argument("--linkage-method", dest="linkage_method", default="ward",
                   help="Linkage method (default: ward)")
    p.add_argument("--row-height", dest="row_height", type=int, default=10,
                   help="Height per haplotype row in pixels (default: 10)")
    p.add_argument("--bar-height", dest="bar_height", type=int, default=8,
                   help="Feature bar height within each row (default: 8)")
    p.add_argument("--dendro-width", dest="dendro_width", type=int, default=300,
                   help="Dendrogram panel width in pixels (default: 300)")
    p.add_argument("--bar-panel-width", dest="bar_panel_width", type=int, default=600,
                   help="Feature bar panel width in pixels (default: 600)")
    p.add_argument("--chrom-gap", dest="chrom_gap", type=int, default=8,
                   help="Extra gap between chromosome groups (default: 8)")
    p.add_argument("--png", action="store_true",
                   help="Also export PNG (requires rsvg-convert)")
    p.add_argument("--all-haplotypes", dest="all_haplotypes", action="store_true",
                   help="Show ALL haplotypes instead of one representative per cluster")
    p.add_argument("--dendro-order", dest="dendro_order", action="store_true",
                   help="Use raw dendrogram leaf order instead of chr1→chrX ordering")
    p.add_argument("--centroid-scan", dest="centroid_scan",
                   action=argparse.BooleanOptionalAction, default=True,
                   help="Stage 2: scan Major clusters for high-divergence outliers "
                        "missed by clustering (default: on)")
    p.add_argument("--centroid-sd", dest="centroid_sd", type=float, default=3.0,
                   help="SD threshold for centroid scan outlier detection (default: 3.0)")
    p.add_argument("--background", default="white",
                   choices=["white", "black"],
                   help="Background color (default: white)")
    return p.parse_args()


# ─── Data Loading ─────────────────────────────────────────────────────────────

def load_colors(colors_dir, featureset):
    """Load feature -> (color, opacity) mapping from colors file."""
    # Detect database name from directory
    db = os.path.basename(colors_dir)
    colors_path = os.path.join(colors_dir, f"{db}.{featureset}.colors.txt")
    if not os.path.exists(colors_path):
        print(f"Warning: colors file not found: {colors_path}")
        return {}

    color_map = {}
    with open(colors_path) as f:
        for i, line in enumerate(f):
            parts = line.strip().split()
            if len(parts) >= 2:
                feat, color = parts[0], parts[1]
                if i == 0 and feat.lower() == 'feature':
                    continue
                color_map[feat] = color
                if feat.endswith('_specific'):
                    color_map[feat[:-9]] = color
    return color_map


def load_bed_features(bed_paths, sequences_needed, exclude_features=None):
    """Load BED data for specified sequences.

    Returns:
        seq_feature_data: dict of seq -> [(feature, length), ...]  (ordered by position)
        seq_bed_data: dict of seq -> [{'start', 'stop', 'feature'}, ...]
        seq_chrom: dict of seq -> chromosome
    """
    import fnmatch
    exclude_patterns = []
    if exclude_features:
        exclude_patterns = [p.strip() for p in exclude_features.split(',')]

    def is_excluded(feat):
        for pat in exclude_patterns:
            if fnmatch.fnmatch(feat, pat):
                return True
        return False

    seq_feature_data = {}
    seq_bed_data = defaultdict(list)
    seq_chrom = {}

    for bed_path in bed_paths:
        open_func = gzip.open if bed_path.endswith('.gz') else open
        mode = 'rt' if bed_path.endswith('.gz') else 'r'

        with open_func(bed_path, mode) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 4:
                    continue
                seq, start, end, feature = parts[0], int(parts[1]), int(parts[2]), parts[3]
                chrom = parts[4] if len(parts) >= 5 else 'unknown'

                if seq not in sequences_needed:
                    continue

                seq_chrom[seq] = chrom
                length = end - start
                seq_bed_data[seq].append({
                    'start': start, 'stop': end, 'feature': feature
                })

    # Build ordered feature data and sort bed entries
    for seq in seq_bed_data:
        seq_bed_data[seq].sort(key=lambda x: x['start'])
        feats = []
        for entry in seq_bed_data[seq]:
            feat = entry['feature']
            if not is_excluded(feat):
                feats.append((feat, entry['stop'] - entry['start']))
        seq_feature_data[seq] = feats

    return seq_feature_data, dict(seq_bed_data), seq_chrom


# ─── Matrix Building ─────────────────────────────────────────────────────────

def get_edges(features, edge_mode="directional"):
    """Get transition edges from ordered feature list."""
    if len(features) <= 1:
        return []
    edges = []
    for i in range(len(features) - 1):
        f1, f2 = features[i], features[i + 1]
        if edge_mode == "symmetric":
            edges.append(tuple(sorted([f1, f2])))
        else:
            edges.append((f1, f2))
    return edges


def build_global_matrix(seq_names, seq_feature_data, all_features,
                        edge_mode, matrix_type):
    """Build a global edge-only feature matrix for all sequences.

    Returns:
        matrix: (n_seqs, n_edge_pairs) numpy array
        edge_names: list of edge pair strings
    """
    # Build edge vocabulary
    if edge_mode == "symmetric":
        all_pairs = []
        for i, f1 in enumerate(all_features):
            for f2 in all_features[i + 1:]:
                all_pairs.append(f"{f1}->{f2}")
    else:
        all_pairs = []
        for f1 in all_features:
            for f2 in all_features:
                if f1 != f2:
                    all_pairs.append(f"{f1}->{f2}")

    pair_to_idx = {pair: i for i, pair in enumerate(all_pairs)}
    matrix = np.zeros((len(seq_names), len(all_pairs)), dtype=np.float32)

    for i, seq in enumerate(seq_names):
        feat_data = seq_feature_data.get(seq, [])
        feat_names = [f for f, _ in feat_data]
        edges = get_edges(feat_names, edge_mode)
        for f1, f2 in edges:
            pair_name = f"{f1}->{f2}"
            if pair_name in pair_to_idx:
                matrix[i, pair_to_idx[pair_name]] += 1

    # Apply normalization
    if matrix_type in ("count_log1p", "count_log1p_zscore"):
        matrix = np.log1p(matrix)

    if matrix_type == "count_log1p_zscore":
        scaler = StandardScaler()
        matrix = scaler.fit_transform(matrix)

    return matrix, all_pairs


# ─── Representative Selection ────────────────────────────────────────────────

def select_representatives(assignments_df):
    """Select one representative per cluster.

    For Major clusters: pick the haplotype closest to centroid (lowest raw_divergence).
    For Outlier clusters: pick the haplotype with highest raw_divergence.
    """
    reps = []
    for cluster_id, group in assignments_df.groupby('cluster'):
        ctype = group.iloc[0]['cluster_type']
        if ctype == 'Major':
            rep = group.loc[group['raw_divergence'].idxmin()]
        else:
            rep = group.loc[group['raw_divergence'].idxmax()]
        reps.append({
            'sequence': rep['sequence'],
            'chromosome': rep['chromosome'],
            'cluster': cluster_id,
            'cluster_type': ctype,
            'sample': rep['sample'],
            'raw_divergence': rep['raw_divergence'],
            'norm_divergence': rep['norm_divergence'],
            'cluster_size': len(group),
        })
    return pd.DataFrame(reps)


# ─── Drawing ──────────────────────────────────────────────────────────────────

def draw_dendrogram_panel(d, Z, leaf_order, row_y_centers, x_base, width,
                          line_color='#555555'):
    """Draw the dendrogram to the left of feature bars.

    The dendrogram is drawn so leaf tips are at x_base + width (flush with bars)
    and root extends toward x_base.
    """
    n_leaves = len(leaf_order)
    if n_leaves < 2:
        return

    max_dist = Z[-1, 2] if len(Z) > 0 else 1.0
    if max_dist == 0:
        max_dist = 1.0

    # Map leaf index in dendrogram order -> y position
    leaf_y = {}
    for dendro_pos, orig_idx in enumerate(leaf_order):
        leaf_y[orig_idx] = row_y_centers[dendro_pos]

    # Process internal nodes
    node_x = {}
    node_y = {}

    for i in range(n_leaves):
        node_x[i] = x_base + width  # leaves at right edge
        node_y[i] = leaf_y[i]

    for i, (idx1, idx2, dist, count) in enumerate(Z):
        idx1, idx2 = int(idx1), int(idx2)
        node_idx = n_leaves + i
        y1, y2 = node_y[idx1], node_y[idx2]
        x_merge = x_base + width - (dist / max_dist) * width

        # Horizontal lines to merge point
        d.append(draw.Line(node_x[idx1], y1, x_merge, y1,
                           stroke=line_color, stroke_width=0.8))
        d.append(draw.Line(node_x[idx2], y2, x_merge, y2,
                           stroke=line_color, stroke_width=0.8))
        # Vertical connector
        d.append(draw.Line(x_merge, y1, x_merge, y2,
                           stroke=line_color, stroke_width=0.8))

        node_x[node_idx] = x_merge
        node_y[node_idx] = (y1 + y2) / 2


def draw_feature_bar(d, bed_entries, x_start, y_center, bar_height,
                     bar_width, color_map, max_bp):
    """Draw a single feature bar (centromere composition) for one haplotype."""
    if not bed_entries:
        return

    ratio = bar_width / max_bp if max_bp > 0 else 1.0

    # Draw background
    d.append(draw.Rectangle(x_start, y_center - bar_height / 2,
                            bar_width, bar_height,
                            fill='#E8E8E8', fill_opacity=0.3))

    # Find the min start to normalize positions
    min_start = min(e['start'] for e in bed_entries)

    for entry in bed_entries:
        feat = entry['feature']
        s = (entry['start'] - min_start) * ratio
        w = max((entry['stop'] - entry['start']) * ratio, 0.5)
        color = color_map.get(feat, '#CCCCCC')

        d.append(draw.Rectangle(x_start + s, y_center - bar_height / 2,
                                w, bar_height,
                                fill=color, fill_opacity=1.0))


def draw_legend(d, color_map, x_start, y_start, text_color='black'):
    """Draw a color legend for features."""
    # Collect unique features that have known colors
    # Use the region colors order
    swatch_size = 12
    col_width = 180
    items_per_col = 12
    x = x_start
    y = y_start

    items = [(feat, color) for feat, color in color_map.items()
             if not feat.endswith('_specific')]  # avoid duplicates

    # Deduplicate keeping order
    seen = set()
    unique_items = []
    for feat, color in items:
        if feat not in seen:
            seen.add(feat)
            unique_items.append((feat, color))

    for i, (feat, color) in enumerate(unique_items):
        col = i // items_per_col
        row = i % items_per_col
        cx = x + col * col_width
        cy = y + row * (swatch_size + 4)

        d.append(draw.Rectangle(cx, cy, swatch_size, swatch_size,
                                fill=color, stroke='#999', stroke_width=0.5))
        d.append(draw.Text(feat, 10, cx + swatch_size + 4, cy + swatch_size - 2,
                           fill=text_color, font_family='sans-serif'))

    return len(unique_items) // items_per_col + 1  # number of columns used


def svg_to_png(svg_path):
    """Convert SVG to PNG using rsvg-convert."""
    png_path = svg_path.rsplit('.svg', 1)[0] + '.png'
    try:
        subprocess.run(['rsvg-convert', '-z', '2', '-f', 'png',
                        '-o', png_path, svg_path],
                       check=True, capture_output=True)
        print(f"  Exported PNG: {png_path}")
    except FileNotFoundError:
        print(f"  Warning: rsvg-convert not found, skipping PNG export")
    except subprocess.CalledProcessError as e:
        print(f"  Warning: PNG export failed: {e.stderr.decode().strip()}")


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    print("=" * 60)
    print("KS All-Chromosome Dendrogram")
    print("=" * 60)

    # ── Load assignments ──────────────────────────────────────────────────
    asgn = pd.read_csv(args.assignments, sep='\t')
    asgn = asgn[~asgn['chromosome'].isin(EXCLUDED_CHROMS)]
    print(f"Loaded {len(asgn)} haplotypes across "
          f"{asgn['chromosome'].nunique()} chromosomes (after exclusion)")

    # ── Stage 2: centroid scan for missed outliers in Major clusters ─────
    if args.centroid_scan:
        print(f"\nStage 2: centroid scan (threshold={args.centroid_sd} SD)...")
        # Need BED features for Major members to build per-chrom matrices
        major_seqs_all = set(asgn[asgn['cluster_type'] == 'Major']['sequence'])
        major_feature_data, _, _ = load_bed_features(
            args.bed, major_seqs_all, args.exclude_features)

        stage2_total = 0
        for chrom in CHROM_ORDER:
            chrom_df = asgn[asgn['chromosome'] == chrom]
            major_df = chrom_df[chrom_df['cluster_type'] == 'Major']
            major_seqs = [s for s in major_df['sequence'] if s in major_feature_data]

            if len(major_seqs) < 5:
                continue

            # Build per-chromosome edge matrix for Major members
            chrom_feats = sorted(set(
                f for s in major_seqs for f, _ in major_feature_data.get(s, [])))
            if not chrom_feats:
                continue

            chrom_matrix, _ = build_global_matrix(
                major_seqs, major_feature_data, chrom_feats,
                args.edge_mode, args.matrix_type)

            centroid = chrom_matrix.mean(axis=0)
            distances = np.linalg.norm(chrom_matrix - centroid, axis=1)
            mean_d, std_d = distances.mean(), distances.std()

            if std_d == 0:
                continue

            cutoff = mean_d + args.centroid_sd * std_d
            flagged = []
            for i, seq in enumerate(major_seqs):
                if distances[i] > cutoff:
                    sd_val = (distances[i] - mean_d) / std_d
                    flagged.append((seq, distances[i], sd_val))

            if flagged:
                flagged.sort(key=lambda x: -x[1])
                stage2_total += len(flagged)
                print(f"  {chrom}: {len(flagged)} rescued from Major "
                      f"(n={len(major_seqs)})")
                for seq, dist, sd in flagged:
                    sample = seq.split('#')[0] if '#' in seq else seq
                    print(f"    {sample} ({sd:.1f} SD)")

                    # Reclassify in asgn: create new outlier cluster
                    existing_outliers = chrom_df[
                        chrom_df['cluster_type'] != 'Major']['cluster'].unique()
                    next_outlier_num = len(existing_outliers) + 1
                    # Check if we already added stage2 outliers for this chrom
                    stage2_mask = asgn['cluster'].str.contains(
                        f'{chrom}_Outlier_S2_', na=False)
                    n_existing_s2 = asgn[stage2_mask]['cluster'].nunique()
                    new_cluster = f"{chrom}_Outlier_S2_{n_existing_s2 + 1}"

                    mask = asgn['sequence'] == seq
                    asgn.loc[mask, 'cluster'] = new_cluster
                    asgn.loc[mask, 'cluster_type'] = 'Outlier'
                    asgn.loc[mask, 'enrichment'] = 'Outlier'

        print(f"  Stage 2 total: {stage2_total} outliers rescued")

    # ── Select representatives ────────────────────────────────────────────
    if args.all_haplotypes:
        reps_df = asgn.copy()
        reps_df['cluster_size'] = reps_df.groupby('cluster')['cluster'].transform('count')
        print(f"Using ALL {len(reps_df)} haplotypes")
    else:
        reps_df = select_representatives(asgn)
        print(f"Selected {len(reps_df)} cluster representatives")

    # Sort by chromosome order, then Major first, then by divergence
    chrom_rank = {c: i for i, c in enumerate(CHROM_ORDER)}
    reps_df['_chrom_rank'] = reps_df['chromosome'].map(chrom_rank).fillna(99)
    reps_df['_type_rank'] = (reps_df['cluster_type'] != 'Major').astype(int)
    reps_df = reps_df.sort_values(['_chrom_rank', '_type_rank', 'raw_divergence'])

    all_sequences = set(reps_df['sequence'].tolist())
    print(f"Loading BED data for {len(all_sequences)} sequences...")

    # ── Load BED data ─────────────────────────────────────────────────────
    seq_feature_data, seq_bed_data, seq_chrom = load_bed_features(
        args.bed, all_sequences, args.exclude_features)
    print(f"  Loaded features for {len(seq_feature_data)} sequences")

    # Filter reps to those with BED data
    reps_df = reps_df[reps_df['sequence'].isin(seq_feature_data)]
    print(f"  {len(reps_df)} representatives have BED data")

    # ── Build GLOBAL feature matrix ───────────────────────────────────────
    # Collect all unique features across all sequences
    all_features_set = set()
    for seq, feat_list in seq_feature_data.items():
        for feat, _ in feat_list:
            all_features_set.add(feat)
    all_features = sorted(all_features_set)
    print(f"  Feature vocabulary: {len(all_features)} features")

    seq_names = reps_df['sequence'].tolist()
    matrix, edge_names = build_global_matrix(
        seq_names, seq_feature_data, all_features,
        args.edge_mode, args.matrix_type)
    print(f"  Global matrix: {matrix.shape}")

    # ── Compute GLOBAL linkage ────────────────────────────────────────────
    print("Computing global pairwise distances...")
    dist = pdist(matrix, metric='euclidean')
    Z = linkage(dist, method=args.linkage_method)
    leaf_order = leaves_list(Z)
    print(f"  Linkage computed ({args.linkage_method})")

    # ── Compute max bp for bar scaling ────────────────────────────────────
    # Use 95th percentile to avoid letting huge arm-spanning outliers
    # compress all normal centromere bars into slivers.
    all_spans = []
    for seq in seq_names:
        entries = seq_bed_data.get(seq, [])
        if entries:
            min_s = min(e['start'] for e in entries)
            max_e = max(e['stop'] for e in entries)
            all_spans.append(max_e - min_s)
    if all_spans:
        max_bp = int(np.percentile(all_spans, 95))
        abs_max = max(all_spans)
        print(f"  Centromere spans: median {int(np.median(all_spans))/1e6:.1f} Mb, "
              f"95th pct {max_bp/1e6:.1f} Mb, max {abs_max/1e6:.1f} Mb")
        if abs_max > max_bp * 3:
            print(f"  Warning: {sum(1 for s in all_spans if s > max_bp)} haplotypes "
                  f"exceed 95th pct — bars will be clipped")
    else:
        max_bp = 10_000_000

    # ── Load colors ───────────────────────────────────────────────────────
    color_map = load_colors(args.colors, args.featureset)
    print(f"  Loaded {len(color_map)} color mappings")

    # ── Layout computation ────────────────────────────────────────────────
    is_dark = args.background == 'black'
    bg_color = '#1a1a1a' if is_dark else 'white'
    text_color = 'white' if is_dark else 'black'
    line_color = '#888888' if is_dark else '#555555'

    n_rows = len(seq_names)
    row_height = args.row_height
    bar_height = args.bar_height
    dendro_width = args.dendro_width
    bar_panel_width = args.bar_panel_width
    chrom_gap = args.chrom_gap
    label_width = 120  # chromosome labels
    margin_top = 60
    margin_bottom = 220  # legend space
    margin_left = 20
    margin_right = 30

    # ── Leaf ordering ────────────────────────────────────────────────────
    raw_leaf_order = leaf_order.tolist()

    if args.dendro_order:
        # Use raw dendrogram leaf order (clusters by composition similarity)
        ordered_indices = raw_leaf_order
        print("  Leaf order: raw dendrogram")
    else:
        # Constrained: chr1→chrX with within-chrom dendrogram ordering
        from collections import OrderedDict
        chrom_leaves = OrderedDict()
        seq_to_chrom_tmp = dict(zip(reps_df['sequence'], reps_df['chromosome']))
        for pos, orig_idx in enumerate(raw_leaf_order):
            seq = seq_names[orig_idx]
            chrom = seq_to_chrom_tmp[seq]
            chrom_leaves.setdefault(chrom, []).append(orig_idx)

        ordered_indices = []
        for chrom in CHROM_ORDER:
            if chrom in chrom_leaves:
                ordered_indices.extend(chrom_leaves[chrom])
        for chrom, indices in chrom_leaves.items():
            if chrom not in CHROM_ORDER:
                ordered_indices.extend(indices)
        print("  Leaf order: chr1→chrX (constrained)")

    ordered_seqs = [seq_names[i] for i in ordered_indices]

    # Build seq -> chrom lookup
    seq_to_chrom = dict(zip(reps_df['sequence'], reps_df['chromosome']))
    seq_to_cluster = dict(zip(reps_df['sequence'], reps_df['cluster']))
    seq_to_type = dict(zip(reps_df['sequence'], reps_df['cluster_type']))
    seq_to_sample = dict(zip(reps_df['sequence'], reps_df['sample']))
    seq_to_size = dict(zip(reps_df['sequence'],
                           reps_df.get('cluster_size', pd.Series(1, index=reps_df.index))))

    row_y_centers = []
    y = margin_top
    for i, seq in enumerate(ordered_seqs):
        y += row_height / 2
        row_y_centers.append(y)
        y += row_height / 2

    total_height = y + margin_bottom
    total_width = margin_left + dendro_width + label_width + bar_panel_width + margin_right

    print(f"\n  Canvas: {total_width} x {total_height} px")
    print(f"  Rows: {n_rows}, row height: {row_height}px, bar height: {bar_height}px")

    # ── Draw SVG ──────────────────────────────────────────────────────────
    print("Drawing SVG...")
    d = draw.Drawing(total_width, total_height, displayInline=False)
    d.append(draw.Rectangle(0, 0, total_width, total_height, fill=bg_color))

    # Title
    d.append(draw.Text("KaryoScope: All-Chromosome Structural Dendrogram", 20,
                        total_width / 2, 30, fill=text_color,
                        font_weight='bold', text_anchor='middle',
                        font_family='sans-serif'))

    # ── Per-row chromosome coloring ──────────────────────────────────────
    # Color each row's background by its chromosome (handles interleaved rows)
    chrom_color_map = {}
    for ci, chrom in enumerate(CHROM_ORDER):
        chrom_color_map[chrom] = CHROM_BAND_COLORS[ci % len(CHROM_BAND_COLORS)]

    band_x = margin_left + dendro_width
    band_width = label_width + bar_panel_width + margin_right
    if not is_dark:
        for i, seq in enumerate(ordered_seqs):
            chrom = seq_to_chrom[seq]
            yc = row_y_centers[i]
            band_color = chrom_color_map.get(chrom, '#F5F5F5')
            d.append(draw.Rectangle(band_x, yc - row_height / 2, band_width,
                                    row_height, fill=band_color, fill_opacity=0.4))

    # ── Labels between dendrogram and bars ─────────────────────────────────
    label_x = margin_left + dendro_width + 4

    # Chromosome label for each row (between dendrogram and bars)
    # Format: "chrN Major n=X" or "chrN [sample] n=X"
    # where n is the cluster size (so all n values per chrom sum to total haplotypes)
    for i, seq in enumerate(ordered_seqs):
        yc = row_y_centers[i]
        chrom = seq_to_chrom[seq]
        ctype = seq_to_type.get(seq, 'Major')
        sample = seq_to_sample.get(seq, 'unknown')
        if sample == 'pangenome' and '#' in seq:
            sample = seq.split('#')[0]
        size = int(seq_to_size.get(seq, 1))

        if ctype == 'Major':
            label = f"{chrom} Major n={size}"
            d.append(draw.Text(label, 8, label_x, yc + 3,
                               fill='#888888', font_family='sans-serif'))
        else:
            label = f"{chrom} [{sample}] n={size}"
            d.append(draw.Text(label, 7, label_x, yc + 3,
                               fill='#FF4444', font_family='monospace'))

    # ── Draw dendrogram ──────────────────────────────────────────────────
    draw_dendrogram_panel(d, Z, leaf_order, row_y_centers,
                          x_base=margin_left, width=dendro_width,
                          line_color=line_color)

    # ── Draw feature bars ────────────────────────────────────────────────
    bars_x = margin_left + dendro_width + label_width
    for i, seq in enumerate(ordered_seqs):
        yc = row_y_centers[i]
        bed_entries = seq_bed_data.get(seq, [])
        draw_feature_bar(d, bed_entries, bars_x, yc, bar_height,
                         bar_panel_width, color_map, max_bp)

        # Cluster type indicator (small dot)
        ctype = seq_to_type.get(seq, 'Major')
        dot_color = '#FF4444' if ctype == 'Outlier' else '#888888'
        d.append(draw.Circle(bars_x - 6, yc, 2.5, fill=dot_color))

    # (Labels are now drawn between dendrogram and bars above)

    # ── Scale bar ─────────────────────────────────────────────────────────
    scale_y = total_height - margin_bottom + 20
    scale_bp = 1_000_000  # 1 Mb
    scale_px = scale_bp * (bar_panel_width / max_bp)
    d.append(draw.Line(bars_x, scale_y, bars_x + scale_px, scale_y,
                        stroke=text_color, stroke_width=2))
    d.append(draw.Text("1 Mb", 11, bars_x + scale_px + 5, scale_y + 4,
                        fill=text_color, font_family='sans-serif'))

    # ── Legend ────────────────────────────────────────────────────────────
    legend_y = scale_y + 25
    d.append(draw.Text("Legend:", 12, margin_left, legend_y,
                        fill=text_color, font_weight='bold',
                        font_family='sans-serif'))

    # Major / Outlier indicators
    d.append(draw.Circle(margin_left + 10, legend_y + 20, 4, fill='#888888'))
    d.append(draw.Text("Major", 10, margin_left + 20, legend_y + 24,
                        fill='#888888', font_family='sans-serif'))
    d.append(draw.Circle(margin_left + 80, legend_y + 20, 4, fill='#FF4444'))
    d.append(draw.Text("Outlier", 10, margin_left + 90, legend_y + 24,
                        fill='#FF4444', font_family='sans-serif'))

    # Feature color legend
    draw_legend(d, color_map, margin_left + 180, legend_y + 5, text_color)

    # ── Save ──────────────────────────────────────────────────────────────
    d.save_svg(args.output)
    print(f"\nSaved: {args.output}")
    file_size = os.path.getsize(args.output)
    print(f"  Size: {file_size / 1024:.0f} KB")

    if args.png:
        svg_to_png(args.output)

    # ── Validation report ─────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("Validation Check")
    print("=" * 60)

    val_samples = {
        'chr3': ('NA21144#1#CM094092.1', 'Hsat1/aSat deletion'),
        'chr5': ('HG00558#1#CM088494.1', 'HSat3 deletion'),
        'chr9': ('HG02630#1#CM091811.1', 'Inversion'),
    }

    for chrom, (seq, desc) in val_samples.items():
        row = asgn[asgn['sequence'] == seq]
        if not row.empty:
            r = row.iloc[0]
            status = "OUTLIER" if r['cluster_type'] == 'Outlier' else "MAJOR (FAIL)"
            print(f"  {chrom} {seq}: {r['cluster']} — {status} — {desc}")
        else:
            print(f"  {chrom} {seq}: NOT FOUND")

    print("=" * 60)


if __name__ == '__main__':
    main()
