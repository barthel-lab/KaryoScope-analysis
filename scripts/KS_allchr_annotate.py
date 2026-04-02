#!/usr/bin/env python3
"""
KS_allchr_annotate.py — Annotate outlier structural differences vs Major

Reads the dendrogram SVG to identify outlier representatives, then compares
each outlier's BED features against the Major representative for the same
chromosome. Reports block order changes (rearrangements, inversions, gains,
losses) and major abundance shifts (>15% change).

Usage:
  python3 scripts/KS_allchr_annotate.py \
    --svg agent_results/allchr_dendrogram_sil0.5_sd5.svg \
    --assignments agent_results/allchr_structure_count_log1p_zscore_blockweight.sequence_assignments.tsv \
    --bed /path/to/pangenome.ALLchr.centromere.KS_human_CHM13.presmoothed.region.pass.bed \
    --output agent_results/allchr_outlier_annotations.tsv
"""

import argparse
import re
import sys
from collections import defaultdict

import pandas as pd


EXCLUDED_CHROMS = {'chr13', 'chr14', 'chr15', 'chr21', 'chr22', 'chrY'}

KEY_FEATS = [
    'active_specific', 'inactive_specific', 'hor_multigroup1',
    'hsat1A_specific', 'hsat2_specific', 'hsat3_specific',
    'bsat_specific', 'ct_specific', 'monomeric_specific',
    'divergent_specific', 'gsat_specific', 'censat_specific',
]


def parse_args():
    p = argparse.ArgumentParser(description="Annotate outlier structural differences")
    p.add_argument("--svg", required=True, help="Dendrogram SVG to parse labels from")
    p.add_argument("--assignments", required=True, help="sequence_assignments.tsv")
    p.add_argument("--bed", required=True, help="BED file with feature annotations")
    p.add_argument("--output", required=True, help="Output TSV path")
    p.add_argument("--abundance-threshold", dest="abundance_threshold", type=float, default=15.0,
                   help="Min percent change to report abundance difference (default: 15)")
    p.add_argument("--exclude-features", dest="exclude_features", default="novel",
                   help="Comma-separated features to exclude (default: novel)")
    return p.parse_args()


def parse_svg_labels(svg_path):
    """Extract Major/Outlier labels from dendrogram SVG."""
    with open(svg_path) as f:
        content = f.read()

    labels = re.findall(r'<text\s+x="324"\s+y="([^"]+)"[^>]*>([^<]+)</text>', content)
    rows = []
    for y, label in labels:
        label = label.strip()
        if not label.startswith('chr'):
            continue
        m_major = re.match(r'(chr\w+) Major n=(\d+)', label)
        m_outlier = re.match(r'(chr\w+) \[(\w+)\] n=(\d+)', label)
        if m_major:
            rows.append({'chrom': m_major.group(1), 'type': 'Major',
                         'sample': '', 'n': int(m_major.group(2))})
        elif m_outlier:
            rows.append({'chrom': m_outlier.group(1), 'type': 'Outlier',
                         'sample': m_outlier.group(2), 'n': int(m_outlier.group(3))})
    return rows


def load_bed_features(bed_path, seqs_needed, exclude_features):
    """Load BED features for specified sequences."""
    exclude = set(f.strip() for f in exclude_features.split(','))
    feat_bp = defaultdict(lambda: defaultdict(int))
    total_bp = defaultdict(int)
    feat_blocks = defaultdict(list)

    with open(bed_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            seq, start, end, feat, chrom = (
                parts[0], int(parts[1]), int(parts[2]), parts[3], parts[4])
            if seq not in seqs_needed:
                continue
            if feat in exclude:
                continue
            bp = end - start
            feat_bp[seq][feat] += bp
            total_bp[seq] += bp
            feat_blocks[seq].append((start, end, feat))

    for seq in feat_blocks:
        feat_blocks[seq].sort()

    return feat_bp, total_bp, feat_blocks


def get_block_order(blocks):
    """Get simplified block order (collapse consecutive same-feature, skip filler)."""
    order = []
    for _, _, feat in blocks:
        if feat in ('ct_specific', 'monomeric_specific'):
            continue
        if not order or order[-1] != feat:
            order.append(feat)
    return order


def describe_order_change(major_order, outlier_order):
    """Describe how the block order changed."""
    if major_order == outlier_order:
        return None

    if major_order == list(reversed(outlier_order)):
        return "inverted block order"

    major_set = set(major_order)
    outlier_set = set(outlier_order)
    missing = major_set - outlier_set
    gained = outlier_set - major_set

    descs = []
    for f in sorted(missing):
        fname = f.replace('_specific', '').replace('_multigroup1', '')
        descs.append(f"lost {fname} block")
    for f in sorted(gained):
        fname = f.replace('_specific', '').replace('_multigroup1', '')
        descs.append(f"gained {fname} block")

    shared = sorted(major_set & outlier_set)
    if len(shared) >= 2:
        major_shared = [f for f in major_order if f in shared]
        outlier_shared = [f for f in outlier_order if f in shared]
        if major_shared != outlier_shared:
            for mf, of in zip(major_shared, outlier_shared):
                if mf != of:
                    mname = mf.replace('_specific', '').replace('_multigroup1', '')
                    oname = of.replace('_specific', '').replace('_multigroup1', '')
                    descs.append(
                        f"rearranged: {mname} and {oname} swapped positions")
                    break

    if not descs and major_order != outlier_order:
        descs.append("different block arrangement")

    return '; '.join(descs)


def annotate_outlier(major_seq, outlier_seq, feat_bp, total_bp, feat_blocks,
                     abundance_threshold):
    """Compare outlier to Major and return description string."""
    major_feats = feat_bp[major_seq]
    major_total = total_bp[major_seq]
    of = feat_bp[outlier_seq]
    ot = total_bp[outlier_seq]

    if major_total == 0 or ot == 0:
        return 'no BED data'

    diffs = []

    # Size difference
    size_ratio = ot / major_total
    if size_ratio < 0.7:
        diffs.append(f"smaller ({ot / 1e6:.1f} vs {major_total / 1e6:.1f} Mb)")
    elif size_ratio > 1.4:
        diffs.append(f"larger ({ot / 1e6:.1f} vs {major_total / 1e6:.1f} Mb)")

    # Abundance differences (only large changes)
    for feat in KEY_FEATS:
        major_pct = major_feats.get(feat, 0) / major_total * 100
        out_pct = of.get(feat, 0) / ot * 100
        fname = feat.replace('_specific', '').replace('_multigroup1', '')

        if major_pct > 1 and out_pct < 0.1:
            diffs.append(f"no {fname} ({major_pct:.0f}% to 0%)")
        elif out_pct > 1 and major_pct < 0.1:
            diffs.append(f"gained {fname} (0% to {out_pct:.0f}%)")
        elif abs(out_pct - major_pct) >= abundance_threshold:
            direction = "up" if out_pct > major_pct else "down"
            diffs.append(
                f"{fname} {direction} ({major_pct:.0f}% to {out_pct:.0f}%)")

    # Block order comparison
    major_order = get_block_order(feat_blocks.get(major_seq, []))
    outlier_order = get_block_order(feat_blocks.get(outlier_seq, []))
    order_desc = describe_order_change(major_order, outlier_order)
    if order_desc:
        diffs.append(order_desc)

    # Transition count difference
    major_trans = feat_blocks.get(major_seq, [])
    out_trans = feat_blocks.get(outlier_seq, [])
    mt = len(set(
        (major_trans[i][2], major_trans[i + 1][2])
        for i in range(len(major_trans) - 1)
        if major_trans[i][2] != major_trans[i + 1][2]
    )) if len(major_trans) > 1 else 0
    ot_t = len(set(
        (out_trans[i][2], out_trans[i + 1][2])
        for i in range(len(out_trans) - 1)
        if out_trans[i][2] != out_trans[i + 1][2]
    )) if len(out_trans) > 1 else 0
    if abs(ot_t - mt) > 3:
        diffs.append(f"transitions: {mt} to {ot_t}")

    # Deduplicate: remove "no X" when "lost X block" exists
    lost_feats = set()
    for item in diffs:
        m = re.match(r'lost (\w+) block', item)
        if m:
            lost_feats.add(m.group(1))

    final = []
    for item in diffs:
        m = re.match(r'no (\w+) \(', item)
        if m and m.group(1) in lost_feats:
            continue
        final.append(item)

    return '; '.join(final) if final else 'edge pattern difference'


def main():
    args = parse_args()

    # Parse SVG labels
    rows = parse_svg_labels(args.svg)
    outlier_rows = [r for r in rows if r['type'] == 'Outlier']
    print(f"Found {len(outlier_rows)} outliers in SVG")

    # Load assignments to map sample -> sequence
    df = pd.read_csv(args.assignments, sep='\t')
    df = df[~df['chromosome'].isin(EXCLUDED_CHROMS)]

    sample_seq_map = {}
    for _, r in df.iterrows():
        sample = r['sequence'].split('#')[0] if '#' in r['sequence'] else r['sequence']
        key = (r['chromosome'], sample)
        if key not in sample_seq_map:
            sample_seq_map[key] = r['sequence']

    # Major representatives (lowest divergence)
    major_reps = {}
    for chrom, grp in df[df['cluster_type'] == 'Major'].groupby('chromosome'):
        rep = grp.loc[grp['raw_divergence'].idxmin()]
        major_reps[chrom] = rep['sequence']

    # Collect all sequences needed
    all_seqs = set(major_reps.values())
    for r in outlier_rows:
        key = (r['chrom'], r['sample'])
        if key in sample_seq_map:
            all_seqs.add(sample_seq_map[key])

    # Load BED
    print(f"Loading BED features for {len(all_seqs)} sequences...")
    feat_bp, total_bp, feat_blocks = load_bed_features(
        args.bed, all_seqs, args.exclude_features)

    # Annotate each outlier
    output_rows = []
    for r in outlier_rows:
        chrom = r['chrom']
        sample = r['sample']
        n = r['n']

        major_seq = major_reps.get(chrom)
        outlier_seq = sample_seq_map.get((chrom, sample), '')

        if not major_seq or not outlier_seq or outlier_seq not in feat_bp:
            output_rows.append([chrom, sample, n, 'no BED data'])
            continue

        desc = annotate_outlier(major_seq, outlier_seq, feat_bp, total_bp,
                                feat_blocks, args.abundance_threshold)
        output_rows.append([chrom, sample, n, desc])

    # Write TSV
    with open(args.output, 'w') as f:
        f.write('chrom\tsample\tn\tstructural_difference\n')
        for row in output_rows:
            f.write('\t'.join(str(x) for x in row) + '\n')

    print(f"Written {len(output_rows)} outlier annotations to {args.output}")
    for row in output_rows:
        print(f"  {row[0]:>5} [{row[1]:>10}] n={row[2]:<4} {row[3]}")


if __name__ == '__main__':
    main()
