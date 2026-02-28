#!/usr/bin/env python3
"""
KaryoScope Cluster Annotation

Summarizes the dominant features for each cluster based on pre-computed per-read
sequence annotations (from KaryoScope_sequence_annotate.py).

Two-step pipeline:
  1. KaryoScope_sequence_annotate.py: per-read BED → sequence annotations TSV
  2. KaryoScope_cluster_annotate.py (this script): sequence annotations → cluster annotations TSV

Per-feature columns (for featureset 'fs' and feature 'feat'):
  {fs}_readpct__{feat}   % of cluster reads where feature > adaptive threshold
  {fs}_bppct__{feat}     Feature bp / total bp across cluster (0-100)
  {fs}_dmax__{feat}      Max density in any 1kb window (median across reads, 0-100)
  {fs}_dmin__{feat}      Min density in any 1kb window (median across reads, 0-100)
  {fs}_dmedian__{feat}   Median density across all 1kb windows (median across reads, 0-100)
  {fs}_dfirst__{feat}    Density in first 1kb of read (median across reads, 0-100)
  {fs}_dlast__{feat}     Density in last 1kb of read (median across reads, 0-100)

Reads shorter than 1kb use overall coverage fraction for all density stats.

Usage:
  python KaryoScope_cluster_annotate.py \\
    --prefix analysis_output_prefix \\
    --sequence-annotations sequence_annotations.tsv \\
    --output cluster_annotations.tsv

  # With auto-labeling:
  python KaryoScope_cluster_annotate.py \\
    --prefix analysis_output_prefix \\
    --sequence-annotations sequence_annotations.tsv \\
    --auto-label \\
    --output cluster_annotations.tsv

The script automatically finds these files from the prefix:
  - {prefix}.read_assignments.tsv
  - {prefix}.cluster_analysis.tsv
"""

import argparse
import fnmatch
import math
import os
import sys

# Capture original command line for logging
_original_command = ' '.join(sys.argv)

import pandas as pd


def load_sequence_annotations(filepath):
    """Load per-read annotations TSV from sequence_annotate.

    Returns:
        tuple: (DataFrame, detected_featuresets list, feature_names_per_fs dict)
    """
    df = pd.read_csv(filepath, sep='\t')

    # Auto-detect featuresets and feature names from column prefixes
    # Columns are like: {fs}_frac__{feat}, {fs}_bp__{feat}, {fs}_dmax__{feat}, etc.
    featuresets = []
    feature_names_per_fs = {}

    # Look for {fs}_frac__{feat} columns to detect featuresets and features
    for col in df.columns:
        if '_frac__' in col:
            parts = col.split('_frac__', 1)
            fs = parts[0]
            feat = parts[1]
            if fs not in featuresets:
                featuresets.append(fs)
                feature_names_per_fs[fs] = []
            feature_names_per_fs[fs].append(feat)

    # Sort feature names for consistency
    for fs in feature_names_per_fs:
        feature_names_per_fs[fs] = sorted(feature_names_per_fs[fs])

    print(f"  Loaded {len(df)} reads, {len(featuresets)} featuresets")
    for fs in featuresets:
        print(f"    {fs}: {len(feature_names_per_fs[fs])} features")

    return df, featuresets, feature_names_per_fs


def load_adaptive_thresholds(filepath):
    """Load thresholds TSV.

    Returns:
        dict {featureset: {feature: threshold}}
    """
    df = pd.read_csv(filepath, sep='\t')
    thresholds = {}
    for _, row in df.iterrows():
        fs = row['featureset']
        feat = row['feature']
        if fs not in thresholds:
            thresholds[fs] = {}
        thresholds[fs][feat] = row['threshold']
    return thresholds


def matches_any_pattern(feature, patterns):
    """Check if feature matches any of the exclude patterns (supports wildcards)."""
    for pattern in patterns:
        if fnmatch.fnmatch(feature, pattern):
            return True
    return False


def summarize_featureset(cluster_reads, seq_annotations, fs, top_n=3, exclude_patterns=None):
    """Summarize features for a cluster from pre-computed sequence annotations."""
    cluster_ann = seq_annotations[seq_annotations['sequence'].isin(cluster_reads)]

    if len(cluster_ann) == 0:
        return ''

    # Sum bp columns for this featureset across cluster reads
    bp_cols = [c for c in cluster_ann.columns if c.startswith(f'{fs}_bp__')]
    if not bp_cols:
        return ''

    feature_bp = {}
    for col in bp_cols:
        feat = col.split('__', 1)[1]
        feature_bp[feat] = cluster_ann[col].sum()

    feature_bp = pd.Series(feature_bp)
    total_bp = feature_bp.sum()

    if total_bp == 0:
        return ''

    # Filter out excluded features for display only
    if exclude_patterns:
        feature_bp = feature_bp[~feature_bp.index.map(lambda f: matches_any_pattern(f, exclude_patterns))]

    if len(feature_bp) == 0:
        return ''

    # Top features by coverage
    top_features = feature_bp.nlargest(top_n)
    top_str = '; '.join([f"{f}({100*v/total_bp:.1f}%)" for f, v in top_features.items()])

    return top_str


def score_cluster_features(cluster_reads, seq_annotations, fs, thresholds):
    """For each feature, compute % of cluster reads exceeding the adaptive threshold."""
    cluster_ann = seq_annotations[seq_annotations['sequence'].isin(cluster_reads)]
    n_reads = len(cluster_reads)
    scores = {}
    for feature, thresh in thresholds.items():
        col = f'{fs}_frac__{feature}'
        if col in cluster_ann.columns:
            # Zero-fill for reads in cluster_reads but not in cluster_ann
            vals = cluster_ann[col].fillna(0)
            n_sig = (vals > thresh).sum()
            # Account for reads missing from annotations (they contribute 0, which is <= thresh)
            scores[feature] = round(100 * n_sig / n_reads, 1) if n_reads > 0 else 0
        else:
            scores[feature] = 0
    return scores


def compute_cluster_bp_scores(cluster_reads, seq_annotations, fs):
    """Compute bp-weighted feature proportions for a cluster."""
    cluster_ann = seq_annotations[seq_annotations['sequence'].isin(cluster_reads)]
    if len(cluster_ann) == 0:
        return {}

    bp_cols = [c for c in cluster_ann.columns if c.startswith(f'{fs}_bp__')]
    total_bp_col = f'{fs}_total_bp'

    total_bp = cluster_ann[total_bp_col].sum() if total_bp_col in cluster_ann.columns else 0
    if total_bp == 0:
        return {}

    result = {}
    for col in bp_cols:
        feat = col.split('__', 1)[1]
        feat_bp = cluster_ann[col].sum()
        result[feat] = round(100 * feat_bp / total_bp, 2)

    return result


def compute_cluster_window_densities(cluster_reads, seq_annotations, fs, feature_names):
    """Compute cluster-level window density stats from pre-computed per-read annotations.

    Reads pre-computed density columns ({fs}_dmax__{feat} etc.) from sequence annotations,
    zero-pads for missing reads, and takes the median across cluster reads.
    Values are scaled 0-100 for cluster output (input values are 0-1 raw fractions).

    Returns dict of {feature: {'max': v, 'min': v, 'median': v, 'first': v, 'last': v,
                                'terminal': v, 'terminal_min': v, 'max_block': v}}.
    """
    n_total = len(cluster_reads)
    if n_total == 0:
        return {}

    cluster_ann = seq_annotations[seq_annotations['sequence'].isin(cluster_reads)]

    result = {}
    for feat in feature_names:
        feat_result = {}
        for stat_key, col_key in [('max', 'dmax'), ('min', 'dmin'), ('median', 'dmedian'),
                                   ('first', 'dfirst'), ('last', 'dlast'),
                                   ('terminal', 'dterminal'), ('terminal_min', 'dterminal_min'),
                                   ('max_block', 'max_block_bp')]:
            col = f'{fs}_{col_key}__{feat}'
            if col in cluster_ann.columns:
                vals = cluster_ann[col].tolist()
                # Zero-pad for reads in cluster but not in annotations
                n_pad = n_total - len(vals)
                if n_pad > 0:
                    vals.extend([0] * n_pad)
                median_val = pd.Series(vals).median()
                if stat_key == 'max_block':
                    feat_result[stat_key] = round(median_val, 0)
                else:
                    feat_result[stat_key] = round(100 * median_val, 2)
            else:
                feat_result[stat_key] = 0
        result[feat] = feat_result

    return result


def compute_cluster_interspersion(cluster_reads, seq_annotations):
    """Compute cluster-level interspersion from pre-computed per-read annotations."""
    n_total = len(cluster_reads)
    if n_total == 0:
        return {'total': 0.0, 'can_ncan': 0.0, 'tel_sat': 0.0, 'arm_tel': 0.0}

    cluster_ann = seq_annotations[seq_annotations['sequence'].isin(cluster_reads)]

    result = {}
    for key in ['total', 'can_ncan', 'tel_sat', 'arm_tel']:
        col = f'interspersion_{key}'
        if col in cluster_ann.columns:
            vals = cluster_ann[col].tolist()
            # Zero-pad for missing reads
            n_pad = n_total - len(vals)
            if n_pad > 0:
                vals.extend([0] * n_pad)
            result[key] = round(pd.Series(vals).median(), 2)
        else:
            result[key] = 0.0

    return result


def score_read_against_annotation(read_densities, read_interspersion,
                                   cluster_row, featureset_prefix):
    """Score how well a single read matches its cluster's annotation profile.

    Uses the same metrics that auto_label_cluster examines to determine how
    representative a read is of its cluster.

    Args:
        read_densities: dict {feature: {max, min, median, first, last, terminal, ...}}
                        Values are raw fractions (0-1).
        read_interspersion: dict {total, can_ncan, tel_sat, arm_tel}
        cluster_row: dict with all annotation columns for this cluster
        featureset_prefix: e.g. 'telomere_region'

    Returns:
        float: score from 0 to 1, higher = more representative
    """
    pfx = featureset_prefix
    eps = 1e-6

    def _closeness(read_val, cluster_val):
        """How close read_val is to cluster_val, normalized 0-1."""
        denom = max(abs(cluster_val), eps)
        return max(0.0, 1.0 - abs(read_val - cluster_val) / denom)

    # --- Determine cluster label type for weight adaptation ---
    cluster_name = cluster_row.get('cluster_name', '')
    is_ectr = cluster_name.startswith('ECTR') if cluster_name else False
    is_satellite = any(s in cluster_name for s in ['aSat', 'bSat', 'CenSat', 'HSat', 'GSat']) if cluster_name else False

    # --- Adaptive weights ---
    w_terminal = 0.30
    w_dmax = 0.20
    w_block = 0.15
    w_interspersion = 0.15
    w_bppct = 0.10
    w_enrichment = 0.10

    if is_ectr:
        w_terminal = 0.35
        w_dmax = 0.15
        w_bppct = 0.10
    elif is_satellite:
        w_bppct = 0.35
        w_terminal = 0.10
        w_dmax = 0.15
        w_block = 0.10
        w_interspersion = 0.10
        w_enrichment = 0.20

    # --- 1. Terminal telomere density (0-1) ---
    terminal_scores = []
    for tel_feat in ['canonical_telomere', 'noncanonical_telomere']:
        for stat in ['first', 'last', 'terminal']:
            col = f'{pfx}_d{stat}__{tel_feat}'
            cluster_val = cluster_row.get(col, 0) / 100.0  # cluster values are 0-100
            read_val = read_densities.get(tel_feat, {}).get(stat, 0)
            terminal_scores.append(_closeness(read_val, cluster_val))
    s_terminal = sum(terminal_scores) / max(len(terminal_scores), 1)

    # --- 2. Feature dmax (0-1) ---
    dmax_scores = []
    for feat, feat_stats in read_densities.items():
        col = f'{pfx}_dmax__{feat}'
        cluster_dmax = cluster_row.get(col, 0) / 100.0
        if cluster_dmax > 0.05:  # only score features cluster actually has
            read_dmax = feat_stats.get('max', 0)
            dmax_scores.append(_closeness(read_dmax, cluster_dmax))
    s_dmax = sum(dmax_scores) / max(len(dmax_scores), 1) if dmax_scores else 0.5

    # --- 3. Max contiguous block (0-1) ---
    block_scores = []
    for feat in read_densities:
        col = f'{pfx}_max_block_bp__{feat}'
        cluster_block = cluster_row.get(col, 0)
        if cluster_block > 100:  # only score meaningful blocks
            read_block = read_densities[feat].get('max_block', 0)
            block_scores.append(_closeness(read_block, cluster_block))
    s_block = sum(block_scores) / max(len(block_scores), 1) if block_scores else 0.5

    # --- 4. Interspersion rates (0-1) ---
    interspersion_scores = []
    for key in ['total', 'can_ncan', 'tel_sat', 'arm_tel']:
        col = f'interspersion_{key}'
        cluster_val = cluster_row.get(col, 0)
        read_val = read_interspersion.get(key, 0)
        interspersion_scores.append(_closeness(read_val, cluster_val))
    s_interspersion = sum(interspersion_scores) / max(len(interspersion_scores), 1)

    # --- 5. BP coverage fractions (0-1) ---
    bppct_scores = []
    for feat in read_densities:
        col = f'{pfx}_bppct__{feat}'
        cluster_bppct = cluster_row.get(col, 0) / 100.0
        if cluster_bppct > 0.01:
            # Approximate read bppct from median density
            read_median = read_densities[feat].get('median', 0)
            bppct_scores.append(_closeness(read_median, cluster_bppct))
    s_bppct = sum(bppct_scores) / max(len(bppct_scores), 1) if bppct_scores else 0.5

    # --- 6. Enrichment-specific features (0-1) ---
    enrichment_features = {
        'active': 'active aSat', 'monomeric': 'monomeric aSat',
        'bsat': 'bSat', 'censat': 'CenSat',
        'hsat1A': 'HSat1A', 'hsat2': 'HSat2', 'hsat3': 'HSat3', 'gsat': 'GSat',
        'ITS': 'ITS', 'TAR1': 'TAR1',
    }
    enrichment_scores = []
    for bed_feat, display_name in enrichment_features.items():
        col = f'{pfx}_dmax__{bed_feat}'
        cluster_enrich = cluster_row.get(col, 0) / 100.0
        if cluster_enrich > 0.20:  # cluster has this enrichment
            read_dmax = read_densities.get(bed_feat, {}).get('max', 0)
            enrichment_scores.append(_closeness(read_dmax, cluster_enrich))
    s_enrichment = sum(enrichment_scores) / max(len(enrichment_scores), 1) if enrichment_scores else 0.5

    # --- Weighted sum ---
    score = (w_terminal * s_terminal +
             w_dmax * s_dmax +
             w_block * s_block +
             w_interspersion * s_interspersion +
             w_bppct * s_bppct +
             w_enrichment * s_enrichment)

    return round(score, 4)


def compute_length_score(read_span, target_span, max_log2_dev=3.0):
    """Score how close a read's length is to the target.

    Returns 1.0 when read_span == target, decaying to 0.0 at 2^max_log2_dev × deviation.
    Uses log2 scale so 2x and 0.5x are penalized equally.
    """
    if read_span <= 0 or target_span <= 0:
        return 0.0
    log2_ratio = abs(math.log2(read_span / target_span))
    return max(0.0, 1.0 - log2_ratio / max_log2_dev)


def normalize_representatives_by_length(cluster_candidates, n_per_cluster, target_length=None):
    """Reorder candidates so that index N reads have similar lengths across clusters.

    For each rank position (1..N), pick the candidate from each cluster that
    best optimizes the combined score (feature + length) while staying close
    to the target length.

    Args:
        cluster_candidates: dict {cluster_id: [{'sequence': ..., 'read_span': ..., 'score': ...,
                            'length_score': ..., 'combined': ...}, ...]}
                            Candidates sorted by combined score descending.
        n_per_cluster: number of representatives to select per cluster
        target_length: global target span (e.g. median across all reads). If None,
                       derives from candidate pool.

    Returns:
        dict {cluster_id: [{'sequence': ..., 'read_span': ..., 'score': ..., 'rank': int}, ...]}
    """
    # Use global target, else median of all candidates
    if target_length:
        target_lengths = [target_length] * n_per_cluster
    else:
        all_spans = [c['read_span'] for cands in cluster_candidates.values()
                     for c in cands if c['read_span'] > 0]
        median_span = sorted(all_spans)[len(all_spans) // 2] if all_spans else 0
        target_lengths = [median_span] * n_per_cluster

    print(f"  Representative target lengths by rank: {[f'{l:,}bp' for l in target_lengths]}")

    # For each cluster, assign candidates to ranks by combined score + length proximity
    normalized = {}
    for cluster_id, candidates in cluster_candidates.items():
        available = list(candidates)
        assigned = []

        for rank in range(n_per_cluster):
            if not available:
                break
            target = target_lengths[rank]
            best_idx = max(range(len(available)),
                           key=lambda i: (available[i].get('combined', available[i]['score']),
                                          -abs(available[i]['read_span'] - target)))
            pick = available.pop(best_idx)
            pick['rank'] = rank + 1
            assigned.append(pick)

        normalized[cluster_id] = assigned

    return normalized


def select_annotation_representatives(results, assignments, seq_annotations,
                                       featureset_prefix, n_reps):
    """Select annotation-aware representative reads for each cluster.

    Uses pre-computed per-read density and interspersion data from sequence annotations
    instead of recomputing from BED files.
    """
    pfx = featureset_prefix

    # Detect read_span column
    span_col = 'read_span' if 'read_span' in assignments.columns else 'read_length'
    if span_col not in assignments.columns:
        span_col = None

    # Compute global target length
    if 'read_span' in assignments.columns:
        global_median_span = float(assignments['read_span'].median())
    elif 'read_length' in assignments.columns:
        global_median_span = float(assignments['read_length'].median())
    else:
        global_median_span = None

    print(f"\n  Selecting {n_reps} annotation-aware representatives per cluster...")
    print(f"  Using featureset '{pfx}' for density scoring")
    if global_median_span:
        print(f"  Target read span (global median): {global_median_span:,.0f}bp")

    # Detect feature names from columns
    density_features = []
    for col in seq_annotations.columns:
        if col.startswith(f'{pfx}_dmax__'):
            feat = col.split('__', 1)[1]
            density_features.append(feat)

    cluster_candidates = {}

    for row in results:
        cluster_id = row['cluster_id']
        cluster_reads = set(
            assignments[assignments['cluster'] == cluster_id]['sequence'].tolist()
        )

        if len(cluster_reads) == 0:
            continue

        # Build read_span lookup
        cluster_assignments = assignments[assignments['cluster'] == cluster_id]
        read_spans = {}
        if span_col:
            for _, r in cluster_assignments.iterrows():
                read_spans[r['sequence']] = r.get(span_col, 0)

        # Extract per-read densities from pre-computed columns
        cluster_ann = seq_annotations[seq_annotations['sequence'].isin(cluster_reads)]

        per_read_densities = {}
        for _, ann_row in cluster_ann.iterrows():
            read_id = ann_row['sequence']
            read_result = {}
            for feat in density_features:
                feat_stats = {}
                for stat_key, col_key in [('max', 'dmax'), ('min', 'dmin'), ('median', 'dmedian'),
                                           ('first', 'dfirst'), ('last', 'dlast'),
                                           ('terminal', 'dterminal'), ('terminal_min', 'dterminal_min'),
                                           ('max_block', 'max_block_bp')]:
                    col = f'{pfx}_{col_key}__{feat}'
                    if col in ann_row.index:
                        feat_stats[stat_key] = ann_row[col]
                    else:
                        feat_stats[stat_key] = 0
                if any(v != 0 for v in feat_stats.values()):
                    read_result[feat] = feat_stats
            per_read_densities[read_id] = read_result

        per_read_interspersion = {}
        for _, ann_row in cluster_ann.iterrows():
            read_id = ann_row['sequence']
            inter = {}
            for key in ['total', 'can_ncan', 'tel_sat', 'arm_tel']:
                col = f'interspersion_{key}'
                inter[key] = ann_row.get(col, 0) if col in ann_row.index else 0
            per_read_interspersion[read_id] = inter

        # Score each read
        scored = []
        for read_id in cluster_reads:
            read_dens = per_read_densities.get(read_id, {})
            read_inter = per_read_interspersion.get(read_id,
                         {'total': 0, 'can_ncan': 0, 'tel_sat': 0, 'arm_tel': 0})
            score = score_read_against_annotation(
                read_dens, read_inter, row, featureset_prefix)
            span = read_spans.get(read_id, 0)
            length_score = compute_length_score(span, global_median_span) if global_median_span else 1.0

            combined = math.sqrt(score * length_score) if (score > 0 and length_score > 0) else 0.0

            scored.append({
                'sequence': read_id,
                'read_span': span,
                'score': score,
                'length_score': length_score,
                'combined': combined,
            })

        scored.sort(key=lambda x: x['combined'], reverse=True)
        cluster_candidates[cluster_id] = scored[:max(20, 5 * n_reps)]

        n_reads = len(cluster_reads)
        top = scored[0] if scored else {}
        print(f"    Cluster {cluster_id} ({row.get('cluster_name', '')}): "
              f"{n_reads} reads, top score={top.get('score', 0):.3f}, "
              f"length_score={top.get('length_score', 0):.3f}, "
              f"combined={top.get('combined', 0):.3f}")

    # --- Pass 1: median-targeted selection ---
    normalized = normalize_representatives_by_length(cluster_candidates, n_reps,
                                                     target_length=global_median_span)

    # --- Pass 2: re-select targeting the longest Pass 1 representative ---
    all_rep_spans = [c['read_span'] for cands in normalized.values()
                     for c in cands if c.get('read_span', 0) > 0]
    if all_rep_spans:
        max_rep_span = max(all_rep_spans)
        print(f"\n  Pass 2: re-targeting to longest representative ({max_rep_span:,.0f}bp)")

        for cluster_id, candidates in cluster_candidates.items():
            for c in candidates:
                c['length_score'] = compute_length_score(c['read_span'], max_rep_span)
                c['combined'] = (math.sqrt(c['score'] * c['length_score'])
                                 if (c['score'] > 0 and c['length_score'] > 0) else 0.0)
            candidates.sort(key=lambda x: x['combined'], reverse=True)

        normalized = normalize_representatives_by_length(cluster_candidates, n_reps,
                                                         target_length=max_rep_span)

        pass2_spans = [c['read_span'] for cands in normalized.values()
                       for c in cands if c.get('read_span', 0) > 0]
        if pass2_spans:
            print(f"  Pass 2 representative spans: min={min(pass2_spans):,.0f}bp, "
                  f"median={sorted(pass2_spans)[len(pass2_spans)//2]:,.0f}bp, "
                  f"max={max(pass2_spans):,.0f}bp")

    return normalized


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


# Cache argparse defaults (before parse_args modifies them)
_argparse_defaults = None


def _print_params_and_command(args):
    """Print comprehensive parameters table and original command."""
    defaults = _argparse_defaults or {}

    def _fmt(value, attr_name):
        if value is None:
            s = "None"
        else:
            s = str(value)
        default_val = defaults.get(attr_name)
        if default_val is not None and value == default_val:
            s += " (default)"
        elif value is None and attr_name not in defaults:
            s += " (default)"
        return s

    params = [
        ("prefix", str(args.prefix), None),
        ("output", str(args.output), None),
        ("sequence-annotations", str(args.sequence_annotations), None),
        ("adaptive-thresholds", _fmt(args.adaptive_thresholds, "adaptive_thresholds"), None),
        ("featuresets", _fmt(args.featuresets, "featuresets"), None),
        ("top-n", _fmt(args.top_n, "top_n"), None),
        ("clusters", _fmt(args.clusters, "clusters"), None),
        ("min-size", _fmt(args.min_size, "min_size"), None),
        ("exclude-features", _fmt(args.exclude_features, "exclude_features"), None),
        ("log-file", _fmt(args.log_file, "log_file"), None),
        ("auto-label", _fmt(args.auto_label, "auto_label"), None),
        ("alt-samples", _fmt(args.alt_samples, "alt_samples"), None),
        ("alt-threshold", _fmt(args.alt_threshold, "alt_threshold"), None),
        ("select-reps", _fmt(args.select_representatives, "select_representatives"), None),
    ]

    print("\n" + "=" * 60)
    print("Parameters")
    print("=" * 60)
    print(f"{'Parameter':<25} {'Value':<35}")
    print(f"{'-' * 25} {'-' * 35}")
    for param, value, _ in params:
        print(f"{param:<25} {str(value):<35}")

    print("\n" + "=" * 60)
    print("Command")
    print("=" * 60)
    print(_original_command)


def auto_label_cluster(row, featureset_prefix):
    """Auto-label a cluster using a structural decision tree based on terminal telomere density.

    Primary classifier: terminal telomere density (dfirst/dlast) determines whether reads
    have telomere at both ends (ECTR), one end (Subtelomere), or only internally
    (Interstitial). Subtelomeres with long contiguous canonical telomere blocks
    (max_block_bp >= 6 kb) are labeled "Type II ALT subtelomere".
    Enrichment qualifiers (satellite, TAR1/ITS, rDNA, SegDup) are appended.

    Args:
        row: dict with all annotation columns for one cluster
        featureset_prefix: e.g. 'telomere_region' or 'region'

    Returns:
        Label string, or '' for unlabeled clusters
    """
    pfx = featureset_prefix

    # --- Thresholds ---
    CAN_ECTR = 70         # canonical terminal density for ECTR (each end)
    NCAN_ECTR = 10        # noncanonical terminal density for ECTR (each end)
    CAN_SUB = 15          # canonical terminal density for Subtelomere (one end)
    NCAN_SUB = 5          # noncanonical terminal density for Subtelomere (one end)
    DMAX_HIGH = 35        # dmax to count as "feature present in reads"
    ENRICH_DMAX = 35      # dmax for enrichment qualifiers (satellite, TAR1, ITS, rDNA)
    NCAN_VARIANT = 25     # noncanonical dmax for variant-enriched qualifier
    SAT_DOMINANT = 80     # readpct for satellite-dominant (rule 5)
    CT_ENRICH = 20        # {pfx}_bppct__ct for SegDup enrichment
    ALT_BLOCK_BP = 6000   # max_block_bp threshold for Type II ALT subtelomere
    ARM_PRESENT = 30      # p_arm or q_arm dmax to count as having arm sequence

    # --- Terminal telomere metrics (orientation-independent) ---
    can_dterm      = row.get(f'{pfx}_dterminal__canonical_telomere', 0)
    ncan_dterm     = row.get(f'{pfx}_dterminal__noncanonical_telomere', 0)
    can_dterm_min  = row.get(f'{pfx}_dterminal_min__canonical_telomere', 0)
    ncan_dterm_min = row.get(f'{pfx}_dterminal_min__noncanonical_telomere', 0)
    can_dmax       = row.get(f'{pfx}_dmax__canonical_telomere', 0)
    ncan_dmax      = row.get(f'{pfx}_dmax__noncanonical_telomere', 0)
    tel_dmax       = max(can_dmax, ncan_dmax)
    can_max_block  = row.get(f'{pfx}_max_block_bp__canonical_telomere', 0)
    can_dfirst     = row.get(f'{pfx}_dfirst__canonical_telomere', 0)
    can_dlast      = row.get(f'{pfx}_dlast__canonical_telomere', 0)
    ncan_dfirst    = row.get(f'{pfx}_dfirst__noncanonical_telomere', 0)
    ncan_dlast     = row.get(f'{pfx}_dlast__noncanonical_telomere', 0)

    # End-check helper: is this end telomeric?
    def _end_is_tel(can_d, ncan_d, can_thresh, ncan_thresh):
        return can_d >= can_thresh or ncan_d >= ncan_thresh

    # ECTR: both ends telomeric → check dfirst and dlast independently
    ectr = (_end_is_tel(can_dfirst, ncan_dfirst, CAN_ECTR, NCAN_ECTR) and
            _end_is_tel(can_dlast, ncan_dlast, CAN_ECTR, NCAN_ECTR))
    # Subtelomere: at least one end telomeric → use terminal (stronger end must pass)
    sub = _end_is_tel(can_dterm, ncan_dterm, CAN_SUB, NCAN_SUB)

    # --- Feature metrics ---
    its_dmax  = row.get(f'{pfx}_dmax__ITS', 0)
    tar1_dmax = row.get(f'{pfx}_dmax__TAR1', 0)
    rdna_dmax = row.get('acrocentric_dmax__rDNA', 0)
    ct        = row.get(f'{pfx}_bppct__ct', 0)

    # Arm metrics: detect whether reads extend into a chromosome arm (or SegDup)
    p_arm_dmax = row.get(f'{pfx}_dmax__p_arm', 0)
    q_arm_dmax = row.get(f'{pfx}_dmax__q_arm', 0)
    has_arm = max(p_arm_dmax, q_arm_dmax) >= ARM_PRESENT or ct >= CT_ENRICH

    its_readpct  = row.get(f'{pfx}_readpct__ITS', 0)
    tar1_readpct = row.get(f'{pfx}_readpct__TAR1', 0)
    ncan_dmax_val = row.get(f'{pfx}_dmax__noncanonical_telomere', 0)

    # Satellite: dmax for enrichment, readpct for dominance (rule 5)
    sat_dmax = {
        'active aSat':    row.get(f'{pfx}_dmax__active', 0),
        'monomeric aSat': row.get(f'{pfx}_dmax__monomeric', 0),
        'bSat':           row.get(f'{pfx}_dmax__bsat', 0),
        'CenSat':         row.get(f'{pfx}_dmax__censat', 0),
        'HSat1A':         row.get(f'{pfx}_dmax__hsat1A', 0),
        'HSat2':          row.get(f'{pfx}_dmax__hsat2', 0),
        'HSat3':          row.get(f'{pfx}_dmax__hsat3', 0),
        'GSat':           row.get(f'{pfx}_dmax__gsat', 0),
    }
    max_sat_dmax_name  = max(sat_dmax, key=sat_dmax.get)
    max_sat_dmax_score = sat_dmax[max_sat_dmax_name]

    sat_readpct = {
        'active aSat':    row.get(f'{pfx}_readpct__active', 0),
        'monomeric aSat': row.get(f'{pfx}_readpct__monomeric', 0),
        'bSat':           row.get(f'{pfx}_readpct__bsat', 0),
        'CenSat':         row.get(f'{pfx}_readpct__censat', 0),
        'HSat1A':         row.get(f'{pfx}_readpct__hsat1A', 0),
        'HSat2':          row.get(f'{pfx}_readpct__hsat2', 0),
        'HSat3':          row.get(f'{pfx}_readpct__hsat3', 0),
        'GSat':           row.get(f'{pfx}_readpct__gsat', 0),
    }
    max_sat_readpct_name  = max(sat_readpct, key=sat_readpct.get)
    max_sat_readpct_score = sat_readpct[max_sat_readpct_name]

    # --- Enrichment qualifier formatter ---
    def _format_quals(quals):
        if not quals:
            return ""
        if len(quals) == 1:
            return f" ({quals[0]}-enriched)"
        return " (" + "-, ".join(quals[:-1]) + "-, " + quals[-1] + "-enriched)"

    # --- Enrichment qualifier builder (uses dmax) ---
    def _enrichment_qualifiers():
        quals = []
        if ncan_dmax_val > NCAN_VARIANT:
            quals.append('variant')
        if max_sat_dmax_score >= ENRICH_DMAX:
            quals.append(max_sat_dmax_name)
        if tar1_dmax >= ENRICH_DMAX and its_dmax >= ENRICH_DMAX:
            quals.append('TAR1/ITS')
        elif tar1_dmax >= ENRICH_DMAX:
            quals.append('TAR1')
        elif its_dmax >= ENRICH_DMAX:
            quals.append('ITS')
        if rdna_dmax >= ENRICH_DMAX:
            quals.append('rDNA')
        if ct >= CT_ENRICH:
            quals.append('SegDup')
        return _format_quals(quals)

    # --- Decision tree ---

    # 1. ECTR: telomere at both ends, OR telomere at one end without arm/SegDup
    if ectr or (sub and not has_arm):
        return f"ECTR{_enrichment_qualifiers()}"

    # 2. Subtelomere: telomere at one end with arm or SegDup at the other
    #    Type II ALT subtelomere: long contiguous canonical telomere block (≥6 kb)
    if sub:
        if can_max_block >= ALT_BLOCK_BP:
            return f"Type II ALT subtelomere{_enrichment_qualifiers()}"
        return f"Subtelomere{_enrichment_qualifiers()}"

    # 3. Interstitial telomere: telomere dmax high but not at ends
    if tel_dmax >= DMAX_HIGH:
        return f"Interstitial telomere{_enrichment_qualifiers()}"

    # 4. Interstitial ITS/TAR1: no telomere structure, but ITS/TAR1 enriched (dmax)
    if tar1_dmax >= ENRICH_DMAX or its_dmax >= ENRICH_DMAX:
        if tar1_dmax >= ENRICH_DMAX and its_dmax >= ENRICH_DMAX:
            base = "Interstitial ITS/TAR1"
        elif tar1_dmax >= its_dmax:
            base = "Interstitial TAR1"
        else:
            base = "Interstitial ITS"
        quals = []
        if max_sat_dmax_score >= ENRICH_DMAX:
            quals.append(max_sat_dmax_name)
        if rdna_dmax >= ENRICH_DMAX:
            quals.append('rDNA')
        if ct >= CT_ENRICH:
            quals.append('SegDup')
        suffix = _format_quals(quals)
        return f"{base}{suffix}"

    # 5. Satellite dominant (uses readpct)
    if max_sat_readpct_score >= SAT_DOMINANT:
        return max_sat_readpct_name

    # 6. Unlabeled
    return ""


# --- Feature Importance Analysis ---

def analyze_annotation_importance(result_df):
    """Analysis A: Which annotation columns vary meaningfully across clusters.

    Returns a DataFrame with per-column stats (CV, sparsity, mean, std, IQR)
    and a Spearman correlation matrix for non-zero-variance columns.
    """
    import numpy as np
    from scipy import stats as scipy_stats

    # Extract numeric score columns: all __ columns + interspersion_*
    score_cols = [c for c in result_df.columns
                  if ('__' in c or c.startswith('interspersion_'))]
    if not score_cols:
        print("  WARNING: No score columns found for annotation importance analysis")
        return None, None

    numeric_df = result_df[score_cols].apply(pd.to_numeric, errors='coerce').fillna(0)

    rows = []
    for col in numeric_df.columns:
        vals = numeric_df[col].values
        mean_val = np.mean(vals)
        std_val = np.std(vals)
        cv = std_val / mean_val if mean_val != 0 else 0
        sparsity = np.mean(vals == 0) * 100
        q25, q75 = np.percentile(vals, [25, 75])
        iqr = q75 - q25

        # Determine featureset from column name
        if '__' in col:
            featureset = col.split('__')[0]
            # Strip score type prefix (e.g., telomere_region_bp -> telomere_region)
            for suffix in ('_readpct', '_bppct', '_dmax', '_dmin', '_dmedian', '_dfirst', '_dlast', '_score'):
                if featureset.endswith(suffix):
                    featureset = featureset[:-len(suffix)]
                    break
        elif col.startswith('interspersion_'):
            featureset = 'interspersion'
        else:
            featureset = 'other'

        rows.append({
            'column': col,
            'featureset': featureset,
            'mean': mean_val,
            'std': std_val,
            'cv': cv,
            'sparsity_pct': sparsity,
            'min': np.min(vals),
            'q25': q25,
            'q75': q75,
            'max': np.max(vals),
            'iqr': iqr,
        })

    importance_df = pd.DataFrame(rows).sort_values('cv', ascending=False).reset_index(drop=True)

    # Spearman correlation on non-zero-variance columns
    nonzero_cols = [c for c in numeric_df.columns if numeric_df[c].std() > 0]
    if len(nonzero_cols) >= 2:
        corr_matrix = numeric_df[nonzero_cols].corr(method='spearman')
    else:
        corr_matrix = None

    return importance_df, corr_matrix


def analyze_svd_loadings(npz_path):
    """Analysis B: Which raw BED features drive the top SVD components.

    Returns:
        feature_importance_df: All features ranked by weighted SVD importance
        layers_df: Importance aggregated by layer-1 and layer-2 components
        top_loadings: Top features x top components submatrix for heatmap
        top_feature_names: Names for heatmap rows
        top_component_labels: Labels for heatmap columns
    """
    import numpy as np

    data = np.load(npz_path, allow_pickle=True)

    # Check for required SVD data
    if 'svd_components' not in data or 'svd_feature_names' not in data:
        return None, None, None, None, None

    components = data['svd_components']          # (n_components, n_features)
    var_ratio = data['svd_explained_variance_ratio']  # (n_components,)
    feature_names = data['svd_feature_names']    # (n_features,)
    if hasattr(feature_names, 'tolist'):
        feature_names = feature_names.tolist()

    # Per-feature importance = sum(|loading| * variance_ratio) across components
    importance = np.sum(np.abs(components) * var_ratio[:, np.newaxis], axis=0)

    rows = []
    for i, fname in enumerate(feature_names):
        ftype = fname.split(':')[0] if ':' in fname else 'unknown'
        rows.append({
            'feature': fname,
            'type': ftype,
            'importance': importance[i],
        })

    feature_importance_df = pd.DataFrame(rows).sort_values('importance', ascending=False).reset_index(drop=True)

    # Decompose 2-layer features and aggregate by layer components
    layer1_importance = {}
    layer2_importance = {}
    for i, fname in enumerate(feature_names):
        imp = importance[i]
        # Parse feature name: "edge:ct:nonsubtelomeric->p_arm:canonical_telomere" or "abundance:ct:nonsubtelomeric"
        if fname.startswith('edge:'):
            body = fname[5:]  # strip "edge:"
            parts = body.split('->')
            if len(parts) == 2:
                left_layers = parts[0].split(':')
                right_layers = parts[1].split(':')
                all_l1 = []
                all_l2 = []
                for layers in [left_layers, right_layers]:
                    if len(layers) >= 1:
                        all_l1.append(layers[0])
                    if len(layers) >= 2:
                        all_l2.append(layers[1])
                share = imp / max(len(all_l1), 1)
                for l in all_l1:
                    layer1_importance[l] = layer1_importance.get(l, 0) + share
                share = imp / max(len(all_l2), 1)
                for l in all_l2:
                    layer2_importance[l] = layer2_importance.get(l, 0) + share
        elif fname.startswith('abundance:'):
            body = fname[10:]  # strip "abundance:"
            layers = body.split(':')
            if len(layers) >= 1:
                layer1_importance[layers[0]] = layer1_importance.get(layers[0], 0) + imp
            if len(layers) >= 2:
                layer2_importance[layers[1]] = layer2_importance.get(layers[1], 0) + imp

    layer_rows = []
    for name, imp in sorted(layer1_importance.items(), key=lambda x: -x[1]):
        layer_rows.append({'layer': 'layer1', 'component': name, 'importance': imp})
    for name, imp in sorted(layer2_importance.items(), key=lambda x: -x[1]):
        layer_rows.append({'layer': 'layer2', 'component': name, 'importance': imp})
    layers_df = pd.DataFrame(layer_rows)

    # Top-20 features x top-10 components submatrix for heatmap
    n_top_features = min(20, len(feature_names))
    n_top_components = min(10, components.shape[0])
    top_idx = np.argsort(-importance)[:n_top_features]
    top_loadings = components[:n_top_components, :][:, top_idx].T  # (n_top_features, n_top_components)
    top_feature_names = [feature_names[i] for i in top_idx]
    top_component_labels = [f"PC{j+1} ({var_ratio[j]*100:.1f}%)" for j in range(n_top_components)]

    return feature_importance_df, layers_df, top_loadings, top_feature_names, top_component_labels


def plot_feature_importance(annotation_importance_df, corr_matrix,
                            svd_importance_df, svd_layers_df,
                            top_loadings, top_feature_names, top_component_labels,
                            output_path):
    """Generate 6-panel feature importance PDF."""
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    from matplotlib.backends.backend_pdf import PdfPages
    from pathlib import Path

    # Register Barthel brand fonts
    FONT_DIR = Path.home() / "Documents" / "Barthel-Custom-Powerpoint-Theme" / "fonts"
    if FONT_DIR.exists():
        for font_file in FONT_DIR.glob("BasicSans-*.otf"):
            fm.fontManager.addfont(str(font_file))

    plt.rcParams.update({
        'font.family': 'Basic Sans',
        'pdf.fonttype': 42,
        'font.size': 9,
    })

    # Barthel palette
    COLORS = {
        'green': '#40D392', 'blue': '#60A5FA', 'coral': '#F07167',
        'yellow': '#FBBF24', 'emerald': '#10B981', 'royal': '#3B82F6',
        'lavender': '#C4A9E8', 'gray': '#545454',
    }
    FEATURESET_COLORS = {
        'telomere_region': COLORS['green'],
        'subtelomeric': COLORS['blue'],
        'chromosome': COLORS['coral'],
        'acrocentric': COLORS['yellow'],
        'repeat': COLORS['emerald'],
        'gene': COLORS['royal'],
        'ct': COLORS['lavender'],
        'interspersion': COLORS['gray'],
    }

    has_svd = svd_importance_df is not None

    if has_svd:
        fig, axes = plt.subplots(3, 2, figsize=(16, 20))
    else:
        fig, axes = plt.subplots(2, 2, figsize=(16, 14))
        print("  Note: SVD data not available — producing annotation-only panels (A, B, C)")

    fig.patch.set_facecolor('white')

    # --- Panel A: Annotation CV ranking (top-left) ---
    ax = axes[0, 0]
    df_a = annotation_importance_df.head(30).iloc[::-1]  # reverse for horizontal bar
    colors_a = [FEATURESET_COLORS.get(fs, '#999999') for fs in df_a['featureset']]
    ax.barh(range(len(df_a)), df_a['cv'].values, color=colors_a, edgecolor='none')
    ax.set_yticks(range(len(df_a)))
    # Truncate long labels
    labels_a = [c[:40] + '...' if len(c) > 40 else c for c in df_a['column']]
    ax.set_yticklabels(labels_a, fontsize=7)
    ax.set_xlabel('Coefficient of Variation (CV)')
    ax.set_title('A. Annotation Column Variability (top 30 by CV)', fontweight='bold', loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Legend for featuresets
    from matplotlib.patches import Patch
    legend_items = []
    seen_fs = set()
    for fs in annotation_importance_df['featureset'].unique():
        if fs not in seen_fs:
            seen_fs.add(fs)
            legend_items.append(Patch(facecolor=FEATURESET_COLORS.get(fs, '#999999'), label=fs))
    ax.legend(handles=legend_items, loc='lower right', fontsize=7, framealpha=0.8)

    # --- Panel B: Sparsity vs CV scatter (top-right) ---
    ax = axes[0, 1]
    df_b = annotation_importance_df.copy()
    colors_b = [FEATURESET_COLORS.get(fs, '#999999') for fs in df_b['featureset']]
    ax.scatter(df_b['sparsity_pct'], df_b['cv'], c=colors_b, s=30, alpha=0.7, edgecolors='none')
    # Label top 10 by CV
    for _, row in df_b.head(10).iterrows():
        label = row['column']
        if len(label) > 30:
            label = label[:27] + '...'
        ax.annotate(label, (row['sparsity_pct'], row['cv']),
                     fontsize=6, alpha=0.8, ha='left',
                     xytext=(3, 3), textcoords='offset points')
    ax.set_xlabel('Sparsity (%)')
    ax.set_ylabel('Coefficient of Variation (CV)')
    ax.set_title('B. Sparsity vs Variability', fontweight='bold', loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # --- Panel C: Annotation correlation heatmap (mid-left) ---
    ax = axes[1, 0]
    if corr_matrix is not None and len(corr_matrix) >= 2:
        from scipy.cluster.hierarchy import linkage as _linkage, leaves_list
        from scipy.spatial.distance import squareform

        # Cluster the correlation matrix
        dist = 1 - corr_matrix.abs().values
        np.fill_diagonal(dist, 0)
        dist = np.clip(dist, 0, None)
        # Make symmetric
        dist = (dist + dist.T) / 2
        condensed = squareform(dist, checks=False)
        link = _linkage(condensed, method='average')
        order = leaves_list(link)

        ordered_corr = corr_matrix.iloc[order, order]
        im = ax.imshow(ordered_corr.values, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title('C. Annotation Correlation (clustered Spearman)', fontweight='bold', loc='left')
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='Spearman r')

        # Label a few key features along the axes
        n_corr = len(ordered_corr)
        if n_corr <= 40:
            labels_c = [c[:25] if len(c) > 25 else c for c in ordered_corr.columns]
            ax.set_xticks(range(n_corr))
            ax.set_xticklabels(labels_c, rotation=90, fontsize=5)
            ax.set_yticks(range(n_corr))
            ax.set_yticklabels(labels_c, fontsize=5)
    else:
        ax.text(0.5, 0.5, 'Insufficient non-zero-variance\ncolumns for correlation',
                ha='center', va='center', transform=ax.transAxes, fontsize=10)
        ax.set_title('C. Annotation Correlation', fontweight='bold', loc='left')

    if not has_svd:
        # Panel D placeholder for annotation-only mode
        ax = axes[1, 1]
        ax.text(0.5, 0.5, 'SVD data not available\n(re-run clustering with --reduce-dims)',
                ha='center', va='center', transform=ax.transAxes, fontsize=10, color='#545454')
        ax.set_title('D. SVD Loadings (not available)', fontweight='bold', loc='left')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    else:
        # --- Panel D: SVD loadings heatmap (mid-right) ---
        ax = axes[1, 1]
        vmax = np.max(np.abs(top_loadings)) if top_loadings.size > 0 else 1
        im = ax.imshow(top_loadings, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')
        ax.set_yticks(range(len(top_feature_names)))
        labels_d = [f[:35] + '...' if len(f) > 35 else f for f in top_feature_names]
        ax.set_yticklabels(labels_d, fontsize=6)
        ax.set_xticks(range(len(top_component_labels)))
        ax.set_xticklabels(top_component_labels, rotation=45, ha='right', fontsize=7)
        ax.set_title('D. SVD Loadings (top 20 features x top 10 components)', fontweight='bold', loc='left')
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='Loading')

        # --- Panel E: SVD feature importance bar chart (bottom-left) ---
        ax = axes[2, 0]
        df_e = svd_importance_df.head(30).iloc[::-1]
        colors_e = [COLORS['blue'] if t == 'edge' else COLORS['green'] for t in df_e['type']]
        ax.barh(range(len(df_e)), df_e['importance'].values, color=colors_e, edgecolor='none')
        ax.set_yticks(range(len(df_e)))
        labels_e = [f[:40] + '...' if len(f) > 40 else f for f in df_e['feature']]
        ax.set_yticklabels(labels_e, fontsize=6)
        ax.set_xlabel('Weighted SVD Importance')
        ax.set_title('E. Raw Feature Importance (top 30)', fontweight='bold', loc='left')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.legend(handles=[
            Patch(facecolor=COLORS['blue'], label='edge'),
            Patch(facecolor=COLORS['green'], label='abundance'),
        ], loc='lower right', fontsize=7)

        # --- Panel F: Layer-aggregated importance (bottom-right) ---
        ax = axes[2, 1]
        l1 = svd_layers_df[svd_layers_df['layer'] == 'layer1'].sort_values('importance', ascending=False)
        l2 = svd_layers_df[svd_layers_df['layer'] == 'layer2'].sort_values('importance', ascending=False)

        # Two side-by-side grouped bar sections
        n_l1 = min(len(l1), 15)
        n_l2 = min(len(l2), 15)
        total = n_l1 + n_l2
        if total > 0:
            gap = 1.5
            positions_l1 = list(range(n_l1))
            positions_l2 = [n_l1 + gap + i for i in range(n_l2)]

            if n_l1 > 0:
                ax.barh(positions_l1[::-1], l1['importance'].values[:n_l1],
                        color=COLORS['royal'], edgecolor='none', label='Layer 1 (region type)')
            if n_l2 > 0:
                ax.barh([p for p in reversed(positions_l2)], l2['importance'].values[:n_l2],
                        color=COLORS['coral'], edgecolor='none', label='Layer 2 (subtype)')

            all_positions = positions_l1[::-1] + list(reversed(positions_l2))
            all_labels = list(l1['component'].values[:n_l1]) + list(l2['component'].values[:n_l2])
            ax.set_yticks(all_positions)
            ax.set_yticklabels(all_labels, fontsize=7)
            ax.set_xlabel('Aggregated SVD Importance')
            ax.legend(loc='lower right', fontsize=7)
        else:
            ax.text(0.5, 0.5, 'No layer data available', ha='center', va='center',
                    transform=ax.transAxes)
        ax.set_title('F. Layer-Aggregated Importance', fontweight='bold', loc='left')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.tight_layout()
    with PdfPages(output_path) as pdf:
        pdf.savefig(fig, dpi=150, bbox_inches='tight')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Annotate clusters with dominant features per featureset",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--prefix", required=True,
                        help="Analysis prefix (auto-finds {prefix}.read_assignments.tsv, {prefix}.cluster_analysis.tsv)")
    parser.add_argument("--sequence-annotations", dest="sequence_annotations", required=True,
                        help="Per-read sequence annotations TSV from KaryoScope_sequence_annotate.py")
    parser.add_argument("--adaptive-thresholds", dest="adaptive_thresholds", default=None,
                        help="Adaptive thresholds TSV (default: auto-derived from --sequence-annotations path)")
    parser.add_argument("--featuresets", default="region,subtelomeric,chromosome,acrocentric,repeat,gene",
                        help="Comma-separated featuresets to annotate (default: region,subtelomeric,chromosome,acrocentric,repeat,gene)")
    parser.add_argument("--output", "-o", required=True,
                        help="Output TSV file")
    parser.add_argument("--top-n", dest="top_n", type=int, default=3,
                        help="Number of top features per featureset (default: 3)")
    parser.add_argument("--clusters",
                        help="Comma-separated cluster IDs to analyze (default: all)")
    parser.add_argument("--min-size", dest="min_size", type=int, default=1,
                        help="Minimum cluster size to include (default: 1)")
    parser.add_argument("--exclude-features", dest="exclude_features",
                        default="*multigroup*,*_arm,nonsubtelomeric,nonacrocentric,nonrepeat,categorized,canonical_telomere*",
                        help="Comma-separated features to exclude, supports wildcards (default: '*multigroup*,*_arm,nonsubtelomeric,nonacrocentric,nonrepeat,categorized,canonical_telomere*')")
    parser.add_argument("--log-file", dest="log_file",
                        action=argparse.BooleanOptionalAction, default=True,
                        help="Save console output to {output}.log (default: True)")
    parser.add_argument("--auto-label", dest="auto_label",
                        action="store_true", default=False,
                        help="Auto-label clusters using a decision tree based on feature scores and interspersion metrics")
    parser.add_argument("--feature-importance", dest="feature_importance",
                        action="store_true", default=False,
                        help="Analyze which annotation columns and raw BED features matter most for cluster structure")
    parser.add_argument("--alt-samples", dest="alt_samples", default=None,
                        help="Comma-separated sample prefixes for Type I ALT relabeling (e.g. 'U2OsTrf1Fok1MUT,U2OsTrf1Fok1_WT')")
    parser.add_argument("--alt-threshold", dest="alt_threshold", type=float, default=80,
                        help="ALT sample percentage threshold for Type I ALT relabeling (default: 80)")
    parser.add_argument("--select-representatives", dest="select_representatives",
                        type=int, default=None, metavar="N",
                        help="Select N annotation-aware representative reads per cluster. "
                             "Populates 'representative_read' column(s) in the output TSV.")

    global _argparse_defaults
    _argparse_defaults = {}
    for action in parser._actions:
        if action.dest != 'help' and action.default is not None:
            _argparse_defaults[action.dest] = action.default

    args = parser.parse_args()

    # --- Set up logging ---
    if args.log_file:
        if args.output.endswith('.tsv'):
            log_path = args.output[:-4] + '.log'
        else:
            log_path = args.output + '.log'
        sys.stdout = TeeLogger(log_path)

    print("=" * 60)
    print("KaryoScope Cluster Annotation")
    print("=" * 60)

    _print_params_and_command(args)

    # Derive file paths from prefix
    read_assignments_file = f"{args.prefix}.read_assignments.tsv"
    if not os.path.exists(read_assignments_file):
        read_assignments_file = f"{args.prefix}.sequence_assignments.tsv"
    cluster_analysis_file = f"{args.prefix}.cluster_analysis.tsv"

    print(f"\nPrefix: {args.prefix}")

    # Load read assignments
    if not os.path.exists(read_assignments_file):
        print(f"ERROR: No read/sequence assignments file found for prefix: {args.prefix}")
        sys.exit(1)

    print(f"\nLoading read assignments: {read_assignments_file}")
    assignments = pd.read_csv(read_assignments_file, sep='\t')

    # Column compatibility: normalize to 'sequence'
    if 'sequence' not in assignments.columns and 'read' in assignments.columns:
        assignments.rename(columns={'read': 'sequence'}, inplace=True)

    print(f"  Total reads: {len(assignments)}")
    print(f"  Total clusters: {assignments['cluster'].nunique()}")

    # Get samples
    samples = assignments['sample'].unique().tolist()
    print(f"  Samples: {len(samples)}")

    # Load cluster analysis
    cluster_info = {}
    entity_columns = []  # list of (entity_name, suffix, column_name)
    has_odds_ratio = False
    if os.path.exists(cluster_analysis_file):
        print(f"\nLoading cluster analysis: {cluster_analysis_file}")
        ca = pd.read_csv(cluster_analysis_file, sep='\t')

        # Auto-detect entity stat columns (samples or groups)
        core_columns = {
            'cluster_id', 'size', 'odds_ratio', 'p_value', 'enrichment',
            'centroid_read', 'centroid_sample', 'centroid_group',
            'q_value', 'enrichment_raw',
        }
        stat_suffixes = ['_count', '_pct', '_pval', '_odds']
        seen_entities = []
        for col in ca.columns:
            for suffix in stat_suffixes:
                if col.endswith(suffix):
                    entity = col[:-len(suffix)]
                    if col not in core_columns:
                        entity_columns.append((entity, suffix, col))
                        if entity not in seen_entities:
                            seen_entities.append(entity)
                        break

        if entity_columns:
            print(f"  Detected entities: {seen_entities}")
            print(f"  Entity stat columns: {len(entity_columns)}")

        # Determine which entities were statistically tested (have _pval columns)
        tested_entities = set(e for e, s, _ in entity_columns if s == '_pval')
        if tested_entities:
            # Per-sample mode: keep only tested entities, drop group summaries
            before = len(entity_columns)
            entity_columns = [(e, s, c) for e, s, c in entity_columns if e in tested_entities]
            seen_entities = [e for e in seen_entities if e in tested_entities]
            if len(entity_columns) < before:
                print(f"  Filtered to tested entities: {seen_entities} ({len(entity_columns)} columns, dropped {before - len(entity_columns)} summary columns)")

        # Check if odds_ratio has any non-empty values
        if 'odds_ratio' in ca.columns:
            has_odds_ratio = ca['odds_ratio'].notna().any()

        for _, row in ca.iterrows():
            info = {
                'enrichment': row.get('enrichment', 'unknown'),
                'p_value': row.get('p_value', None),
                'q_value': row.get('q_value', None),
                'odds_ratio': row.get('odds_ratio', None),
            }
            # Store all detected entity columns
            for _entity, _suffix, col in entity_columns:
                val = row.get(col, None)
                if pd.notna(val):
                    info[col] = val
                else:
                    info[col] = None
            cluster_info[row['cluster_id']] = info
    else:
        print(f"\nWARNING: Cluster analysis file not found: {cluster_analysis_file}")

    # Parse exclude patterns
    exclude_patterns = []
    if args.exclude_features:
        exclude_patterns = [p.strip() for p in args.exclude_features.split(',') if p.strip()]
        print(f"\nExcluding features matching: {exclude_patterns}")

    # Load sequence annotations
    print(f"\nLoading sequence annotations: {args.sequence_annotations}")
    seq_annotations, detected_featuresets, feature_names_per_fs = load_sequence_annotations(args.sequence_annotations)

    # Determine featuresets to use (intersection of requested and detected)
    featuresets = [fs.strip() for fs in args.featuresets.split(',')]
    available_featuresets = [fs for fs in featuresets if fs in detected_featuresets]
    missing_featuresets = [fs for fs in featuresets if fs not in detected_featuresets]
    if missing_featuresets:
        print(f"  WARNING: Requested featuresets not in annotations: {missing_featuresets}")
    featuresets = available_featuresets
    print(f"  Featuresets to annotate: {featuresets}")

    feature_names = {fs: feature_names_per_fs[fs] for fs in featuresets}

    # Load adaptive thresholds
    if args.adaptive_thresholds:
        thresh_path = args.adaptive_thresholds
    else:
        # Auto-derive from sequence-annotations path
        sa_prefix = args.sequence_annotations[:-4] if args.sequence_annotations.endswith('.tsv') else args.sequence_annotations
        thresh_path = f"{sa_prefix}.adaptive_thresholds.tsv"

    feature_thresholds = {}
    if os.path.exists(thresh_path):
        print(f"\nLoading adaptive thresholds: {thresh_path}")
        feature_thresholds = load_adaptive_thresholds(thresh_path)
        for fs in featuresets:
            if fs in feature_thresholds:
                print(f"  {fs}: {len(feature_thresholds[fs])} features")
    else:
        print(f"\nWARNING: Adaptive thresholds not found: {thresh_path}")
        print("  Will compute thresholds from sequence annotations")
        # Fallback: compute from frac columns
        for fs in featuresets:
            feature_thresholds[fs] = {}
            for feat in feature_names.get(fs, []):
                col = f'{fs}_frac__{feat}'
                if col in seq_annotations.columns:
                    nonzero = seq_annotations[col][seq_annotations[col] > 0]
                    if len(nonzero) == 0:
                        feature_thresholds[fs][feat] = 0.001
                    else:
                        med = nonzero.median()
                        feature_thresholds[fs][feat] = max(0.001, min(0.05, med / 3))

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
        cluster_reads = set(assignments[assignments['cluster'] == cluster_id]['sequence'].tolist())

        if len(cluster_reads) < args.min_size:
            continue

        # Basic info
        row = {
            'cluster_id': cluster_id,
            'size': len(cluster_reads),
        }

        # Curation columns (empty for user to fill in)
        row['cluster_name'] = ''
        row['curated_rep_i'] = ''

        # Add info from cluster analysis
        if cluster_id in cluster_info:
            info = cluster_info[cluster_id]
            row['enrichment'] = info.get('enrichment', 'unknown')
            q = info.get('q_value')
            row['q_value'] = f"{q:.4e}" if q is not None else None

            # log2_fc only when odds_ratio is available
            if has_odds_ratio:
                odds = info.get('odds_ratio')
                row['log2_fc'] = round(math.log2(odds), 2) if odds is not None and odds > 0 else None

            # Dynamic entity stat columns
            for _entity, suffix, col in entity_columns:
                val = info.get(col)
                if suffix == '_pct' and val is not None:
                    row[col] = round(val, 1)
                elif suffix == '_pval' and val is not None:
                    row[col] = f"{val:.4e}"
                elif suffix == '_odds' and val is not None:
                    row[col] = round(val, 2)
                else:
                    row[col] = val

        # Interspersion metrics
        has_interspersion = 'interspersion_total' in seq_annotations.columns
        if has_interspersion:
            interspersion = compute_cluster_interspersion(cluster_reads, seq_annotations)
            row['interspersion_total'] = interspersion['total']
            row['interspersion_can_ncan'] = interspersion['can_ncan']
            row['interspersion_tel_sat'] = interspersion['tel_sat']
            row['interspersion_arm_tel'] = interspersion['arm_tel']

        # Annotate each featureset
        for fs in featuresets:
            row[f'{fs}_top'] = summarize_featureset(cluster_reads, seq_annotations, fs, args.top_n, exclude_patterns)
            # Per-feature columns: read-presence, bp-level, window densities
            scores = {}
            bp_scores = {}
            window_densities = {}
            if fs in feature_thresholds:
                scores = score_cluster_features(cluster_reads, seq_annotations, fs, feature_thresholds[fs])
            bp_scores = compute_cluster_bp_scores(cluster_reads, seq_annotations, fs)
            window_densities = compute_cluster_window_densities(cluster_reads, seq_annotations, fs, feature_names.get(fs, []))
            for feat in feature_names.get(fs, []):
                row[f'{fs}_readpct__{feat}'] = scores.get(feat, 0)
                row[f'{fs}_bppct__{feat}'] = bp_scores.get(feat, 0)
                feat_wd = window_densities.get(feat, {})
                row[f'{fs}_dmax__{feat}'] = feat_wd.get('max', 0)
                row[f'{fs}_dmin__{feat}'] = feat_wd.get('min', 0)
                row[f'{fs}_dmedian__{feat}'] = feat_wd.get('median', 0)
                row[f'{fs}_dfirst__{feat}'] = feat_wd.get('first', 0)
                row[f'{fs}_dlast__{feat}'] = feat_wd.get('last', 0)
                row[f'{fs}_dterminal__{feat}'] = feat_wd.get('terminal', 0)
                row[f'{fs}_dterminal_min__{feat}'] = feat_wd.get('terminal_min', 0)
                row[f'{fs}_max_block_bp__{feat}'] = feat_wd.get('max_block', 0)

        if args.auto_label:
            auto_label_fs = ['telomere_region', 'region_subtelomere_flat']
            pfx = next((fs for fs in auto_label_fs if fs in featuresets), None)
            if pfx is None:
                print("ERROR: --auto-label requires 'telomere_region' or 'region_subtelomere_flat' in --featuresets")
                sys.exit(1)
            row['cluster_name'] = auto_label_cluster(row, pfx)

            # Type I ALT relabeling: prepend "Type I ALT" for ALT-enriched clusters
            if args.alt_samples:
                alt_sample_list = [s.strip() for s in args.alt_samples.split(',')]
                alt_pct = sum(row.get(f'{s}_pct', 0) or 0 for s in alt_sample_list)
                label = row['cluster_name']
                if (label
                        and alt_pct > args.alt_threshold
                        and not label.startswith('ECTR')
                        and not label.startswith('Type II ALT')):
                    row['cluster_name'] = 'Type I ALT ' + label[0].lower() + label[1:]

        results.append(row)

    # Create output DataFrame
    result_df = pd.DataFrame(results)

    # Sort by cluster_id (ascending)
    result_df = result_df.sort_values('cluster_id', ascending=True)

    # --- Annotation-aware representative selection ---
    if args.select_representatives:
        n_reps = args.select_representatives
        pfx_for_reps = 'telomere_region' if 'telomere_region' in featuresets else featuresets[0] if featuresets else 'region'
        normalized_reps = select_annotation_representatives(
            results, assignments, seq_annotations, pfx_for_reps, n_reps)

        # Add columns: representative_read_1, representative_read_2, ..., representative_read_N
        for rank in range(1, n_reps + 1):
            col = f'representative_read_{rank}'
            rep_map = {}
            for cid, reps in normalized_reps.items():
                if len(reps) >= rank:
                    rep_map[cid] = reps[rank - 1]['sequence']
                else:
                    rep_map[cid] = ''
            result_df[col] = result_df['cluster_id'].map(rep_map).fillna('')

        # Set curated_rep_i = 1 for all clusters (default to rank 1)
        result_df['curated_rep_i'] = 1

        print(f"\n  Added {n_reps} representative read columns to output")

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

    if args.auto_label and 'cluster_name' in result_df.columns:
        print("\nBy auto-label:")
        for label in sorted(result_df['cluster_name'].unique()):
            count = (result_df['cluster_name'] == label).sum()
            display = label if label else '(unlabeled)'
            print(f"  {display}: {count}")

    if entity_columns:
        seen = []
        for entity, _suffix, _col in entity_columns:
            if entity not in seen:
                seen.append(entity)
        print(f"\nEntity stat columns included for: {', '.join(seen)}")

    # --- Feature Importance Analysis ---
    if args.feature_importance:
        import numpy as np

        # Derive output prefix from output path
        out_prefix = args.output[:-4] if args.output.endswith('.tsv') else args.output

        print(f"\n{'=' * 60}")
        print("Feature Importance Analysis")
        print("=" * 60)

        # Analysis A: Annotation column importance
        print("\n--- Analysis A: Annotation column importance ---")
        ann_importance_df, corr_matrix = analyze_annotation_importance(result_df)
        if ann_importance_df is not None:
            ann_tsv = f"{out_prefix}.feature_importance_annotations.tsv"
            ann_importance_df.to_csv(ann_tsv, sep='\t', index=False)
            print(f"  Saved annotation importance to: {ann_tsv}")
            print(f"  Total score columns analyzed: {len(ann_importance_df)}")
            top5 = ann_importance_df.head(5)
            print(f"  Top 5 by CV: {', '.join(top5['column'])}")

        # Analysis B: SVD loading importance
        npz_path = f"{args.prefix}.feature_matrix.npz"
        svd_importance_df = None
        svd_layers_df = None
        top_loadings = None
        top_feature_names_svd = None
        top_component_labels = None

        if os.path.exists(npz_path):
            print(f"\n--- Analysis B: SVD loading importance ---")
            svd_importance_df, svd_layers_df, top_loadings, top_feature_names_svd, top_component_labels = \
                analyze_svd_loadings(npz_path)

            if svd_importance_df is not None:
                svd_tsv = f"{out_prefix}.feature_importance_svd.tsv"
                svd_importance_df.to_csv(svd_tsv, sep='\t', index=False)
                print(f"  Saved SVD feature importance to: {svd_tsv}")
                print(f"  Total raw features: {len(svd_importance_df)}")
                top5_svd = svd_importance_df.head(5)
                print(f"  Top 5 by importance: {', '.join(top5_svd['feature'])}")

                layers_tsv = f"{out_prefix}.feature_importance_svd_layers.tsv"
                svd_layers_df.to_csv(layers_tsv, sep='\t', index=False)
                print(f"  Saved layer-aggregated importance to: {layers_tsv}")
            else:
                print(f"  WARNING: NPZ lacks SVD data — producing annotation-only panels")
        else:
            print(f"\n  WARNING: NPZ not found at {npz_path} — skipping SVD analysis")

        # Visualization: 6-panel (or 4-panel) PDF
        if ann_importance_df is not None:
            print(f"\n--- Generating feature importance visualization ---")
            pdf_path = f"{out_prefix}.feature_importance.pdf"
            plot_feature_importance(
                ann_importance_df, corr_matrix,
                svd_importance_df, svd_layers_df,
                top_loadings, top_feature_names_svd, top_component_labels,
                pdf_path
            )
            print(f"  Saved feature importance PDF to: {pdf_path}")


if __name__ == "__main__":
    main()
