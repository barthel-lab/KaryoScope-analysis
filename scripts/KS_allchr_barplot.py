#!/usr/bin/env python3
"""
KS_allchr_barplot.py — Stacked bar chart of Major/Outlier haplotype counts

Shows total haplotypes per chromosome as stacked bars (Major + Outlier).
Uses the same sil-threshold + centroid-scan logic as KS_allchr_dendrogram.py
to determine final Major/Outlier assignments.

Usage:
  python3 scripts/KS_allchr_barplot.py \
    --assignments agent_results/allchr_structure.sequence_assignments.tsv \
    --bed /path/to/pangenome.BED \
    --output agent_results/allchr_barplot.svg \
    --sil-threshold 0.5 --centroid-sd 5
"""

import argparse
import os
import subprocess
import sys

import drawsvg as draw
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage as _linkage, fcluster as _fcluster
from scipy.spatial.distance import pdist
from sklearn.metrics import silhouette_score as _sil_score
from sklearn.preprocessing import StandardScaler

# Import shared functions from dendrogram script
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from KS_allchr_dendrogram import (
    EXCLUDED_CHROMS, CHROM_ORDER, CHROM_BLOCK_COLORS,
    load_bed_features, build_global_matrix
)

def parse_args():
    p = argparse.ArgumentParser(
        description="Stacked bar chart of Major/Outlier counts per chromosome")
    p.add_argument("--assignments", required=True,
                   help="sequence_assignments.tsv from structure mode")
    p.add_argument("--bed", required=True, nargs='+',
                   help="BED file(s) with feature annotations")
    p.add_argument("--output", required=True,
                   help="Output SVG path")
    p.add_argument("--exclude-features", dest="exclude_features",
                   default="novel",
                   help="Comma-separated features to exclude (default: novel)")
    p.add_argument("--edge-mode", dest="edge_mode", default="directional",
                   help="Edge counting mode (default: directional)")
    p.add_argument("--matrix-type", dest="matrix_type",
                   default="count_log1p_zscore",
                   help="Matrix normalization (default: count_log1p_zscore)")
    p.add_argument("--sil-threshold", dest="sil_threshold", type=float, default=0.0,
                   help="Silhouette threshold for collapsing weak splits (default: 0.0)")
    p.add_argument("--centroid-sd", dest="centroid_sd", type=float, default=3.0,
                   help="SD threshold for centroid scan (default: 3.0)")
    p.add_argument("--centroid-scan", dest="centroid_scan",
                   action=argparse.BooleanOptionalAction, default=True,
                   help="Enable stage 2 centroid scan (default: on)")
    p.add_argument("--png", action="store_true",
                   help="Also export PNG")
    return p.parse_args()


def apply_sil_and_centroid(asgn, args):
    """Apply silhouette threshold + centroid scan (same logic as dendrogram)."""
    if args.sil_threshold > 0:
        all_seqs_set = set(asgn['sequence'])
        all_feat_data, _, _ = load_bed_features(
            args.bed, all_seqs_set, args.exclude_features)

        for chrom in CHROM_ORDER:
            chrom_df = asgn[asgn['chromosome'] == chrom]
            if chrom_df['cluster'].nunique() <= 1:
                continue

            chrom_seqs = [s for s in sorted(chrom_df['sequence']) if s in all_feat_data]
            if len(chrom_seqs) < 5:
                continue

            chrom_feats = sorted(set(
                f for s in chrom_seqs for f, _ in all_feat_data.get(s, [])))
            if not chrom_feats:
                continue

            chrom_matrix, _ = build_global_matrix(
                chrom_seqs, all_feat_data, chrom_feats,
                args.edge_mode, args.matrix_type)

            _dist = pdist(chrom_matrix, metric='euclidean')
            _Z = _linkage(_dist, method='ward')
            _labels = _fcluster(_Z, 2, criterion='maxclust')
            sil = _sil_score(chrom_matrix, _labels) if len(set(_labels)) >= 2 else 0.0

            if sil < args.sil_threshold:
                major_cluster = f"{chrom}_Major"
                mask = asgn['chromosome'] == chrom
                asgn.loc[mask, 'cluster'] = major_cluster
                asgn.loc[mask, 'cluster_type'] = 'Major'
                asgn.loc[mask, 'enrichment'] = 'Major'

    if args.centroid_scan:
        major_seqs_all = set(asgn[asgn['cluster_type'] == 'Major']['sequence'])
        major_feature_data, _, _ = load_bed_features(
            args.bed, major_seqs_all, args.exclude_features)

        for chrom in CHROM_ORDER:
            chrom_df = asgn[asgn['chromosome'] == chrom]
            major_df = chrom_df[chrom_df['cluster_type'] == 'Major']
            major_seqs = [s for s in major_df['sequence'] if s in major_feature_data]

            if len(major_seqs) < 5:
                continue

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
            for i, seq in enumerate(major_seqs):
                if distances[i] > cutoff:
                    stage2_mask = asgn['cluster'].str.contains(
                        f'{chrom}_Outlier_S2_', na=False)
                    n_existing_s2 = asgn[stage2_mask]['cluster'].nunique()
                    new_cluster = f"{chrom}_Outlier_S2_{n_existing_s2 + 1}"
                    mask = asgn['sequence'] == seq
                    asgn.loc[mask, 'cluster'] = new_cluster
                    asgn.loc[mask, 'cluster_type'] = 'Outlier'
                    asgn.loc[mask, 'enrichment'] = 'Outlier'

    return asgn


def main():
    args = parse_args()

    print("=" * 60)
    print("KS All-Chromosome Barplot")
    print("=" * 60)

    asgn = pd.read_csv(args.assignments, sep='\t')
    asgn = asgn[~asgn['chromosome'].isin(EXCLUDED_CHROMS)]
    print(f"Loaded {len(asgn)} haplotypes across "
          f"{asgn['chromosome'].nunique()} chromosomes")

    asgn = apply_sil_and_centroid(asgn, args)

    # Count Major/Outlier per chromosome
    chrom_block_map = {c: CHROM_BLOCK_COLORS[i % len(CHROM_BLOCK_COLORS)]
                       for i, c in enumerate(CHROM_ORDER)}

    counts = []
    for chrom in CHROM_ORDER:
        cdf = asgn[asgn['chromosome'] == chrom]
        n_major = (cdf['cluster_type'] == 'Major').sum()
        n_outlier = (cdf['cluster_type'] != 'Major').sum()
        counts.append({
            'chrom': chrom, 'major': n_major, 'outlier': n_outlier,
            'total': n_major + n_outlier
        })

    counts_df = pd.DataFrame(counts)
    print("\nPer-chromosome counts:")
    for _, r in counts_df.iterrows():
        print(f"  {r['chrom']}: {r['major']} Major, {r['outlier']} Outlier "
              f"(total {r['total']})")
    print(f"  Grand total: {counts_df['total'].sum()} "
          f"({counts_df['major'].sum()} Major, {counts_df['outlier'].sum()} Outlier)")

    # ── Draw SVG ──────────────────────────────────────────────────────────
    n_chroms = len(counts_df)
    bar_width = 40
    bar_gap = 8
    margin_left = 70
    margin_right = 50
    margin_top = 60
    margin_bottom = 80
    max_count = counts_df['total'].max()
    plot_height = 500
    scale = plot_height / max_count

    total_width = margin_left + n_chroms * (bar_width + bar_gap) + margin_right
    total_height = margin_top + plot_height + margin_bottom
    # Ensure width >= height so square-crop renderers (qlmanage) don't clip
    if total_width > total_height:
        margin_bottom = total_width - margin_top - plot_height
        total_height = total_width

    d = draw.Drawing(total_width, total_height, displayInline=False)
    d.append(draw.Rectangle(0, 0, total_width, total_height, fill='white'))

    # Title
    d.append(draw.Text("Haplotypes per Chromosome: Major vs Outlier", 16,
                        total_width / 2, 25, fill='black',
                        font_weight='bold', text_anchor='middle',
                        font_family='sans-serif'))

    # Y-axis
    y_base = margin_top + plot_height
    d.append(draw.Line(margin_left, margin_top, margin_left, y_base,
                       stroke='#333', stroke_width=1))

    # Y-axis ticks
    tick_step = 50
    y_ticks = range(0, int(max_count) + tick_step, tick_step)
    for val in y_ticks:
        y = y_base - val * scale
        if y < margin_top - 5:
            break
        d.append(draw.Line(margin_left - 4, y, margin_left, y,
                           stroke='#333', stroke_width=1))
        d.append(draw.Text(str(val), 9, margin_left - 8, y + 3,
                           fill='#333', text_anchor='end',
                           font_family='sans-serif'))
        # Grid line
        if val > 0:
            d.append(draw.Line(margin_left, y,
                               margin_left + n_chroms * (bar_width + bar_gap), y,
                               stroke='#E0E0E0', stroke_width=0.5))

    # Y-axis label
    d.append(draw.Text("Number of haplotypes", 11,
                        15, margin_top + plot_height / 2,
                        fill='black', font_family='sans-serif',
                        text_anchor='middle',
                        transform=f"rotate(-90, 15, {margin_top + plot_height / 2})"))

    # X-axis
    d.append(draw.Line(margin_left, y_base,
                       margin_left + n_chroms * (bar_width + bar_gap), y_base,
                       stroke='#333', stroke_width=1))

    # Bars
    major_color = '#4A90D9'
    outlier_color = '#E63946'

    for i, (_, r) in enumerate(counts_df.iterrows()):
        x = margin_left + i * (bar_width + bar_gap) + bar_gap / 2

        # Major (bottom)
        h_major = r['major'] * scale
        d.append(draw.Rectangle(x, y_base - h_major, bar_width, h_major,
                                fill=major_color))

        # Outlier (top, stacked on major)
        h_outlier = r['outlier'] * scale
        if h_outlier > 0:
            d.append(draw.Rectangle(x, y_base - h_major - h_outlier,
                                    bar_width, h_outlier,
                                    fill=outlier_color))

        # Total count label above bar
        total_y = y_base - (h_major + h_outlier) - 5
        d.append(draw.Text(str(r['total']), 8, x + bar_width / 2, total_y,
                           fill='#333', text_anchor='middle',
                           font_family='sans-serif'))

        # Major count label (inside blue bar)
        major_label_y = y_base - h_major / 2 + 3
        d.append(draw.Text(str(r['major']), 7,
                           x + bar_width / 2, major_label_y,
                           fill='white', text_anchor='middle',
                           font_weight='bold', font_family='sans-serif'))

        # Outlier count label (if >0)
        if r['outlier'] > 0:
            outlier_label_y = y_base - h_major - h_outlier / 2 + 3
            d.append(draw.Text(str(r['outlier']), 7,
                               x + bar_width / 2, outlier_label_y,
                               fill='white', text_anchor='middle',
                               font_weight='bold', font_family='sans-serif'))

        # Chromosome label (x-axis)
        chrom_label = r['chrom'].replace('chr', '')
        d.append(draw.Text(chrom_label, 9,
                           x + bar_width / 2, y_base + 14,
                           fill='#333', text_anchor='middle',
                           font_family='sans-serif'))

        # Chromosome color dot under label
        block_color = chrom_block_map.get(r['chrom'], '#CCC')
        d.append(draw.Circle(x + bar_width / 2, y_base + 24, 4,
                             fill=block_color))

    # Legend
    legend_x = total_width - 180
    legend_y = margin_top + 10
    d.append(draw.Rectangle(legend_x, legend_y, 14, 14, fill=major_color))
    d.append(draw.Text("Major", 11, legend_x + 20, legend_y + 12,
                        fill='#333', font_family='sans-serif'))
    d.append(draw.Rectangle(legend_x, legend_y + 22, 14, 14, fill=outlier_color))
    d.append(draw.Text("Outlier", 11, legend_x + 20, legend_y + 34,
                        fill='#333', font_family='sans-serif'))

    # Save
    d.save_svg(args.output)
    print(f"\nSaved: {args.output}")
    print(f"  Size: {os.path.getsize(args.output) / 1024:.0f} KB")

    if args.png:
        png_path = args.output.rsplit('.svg', 1)[0] + '.png'
        try:
            subprocess.run(['rsvg-convert', '-z', '2', '-f', 'png',
                           '-o', png_path, args.output],
                          check=True, capture_output=True)
            print(f"  Exported PNG: {png_path}")
        except (FileNotFoundError, subprocess.CalledProcessError):
            print("  Warning: PNG export skipped (rsvg-convert not available)")


if __name__ == '__main__':
    main()
