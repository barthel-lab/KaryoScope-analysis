#!/usr/bin/env python3
"""
KaryoScope Sequence Annotation

Computes per-read (sequence-level) feature annotations from BED files and
outputs a wide-format TSV with one row per read.

Per-feature columns (for featureset 'fs' and feature 'feat'):
  {fs}_frac__{feat}           Feature fraction (0-1)
  {fs}_bp__{feat}             Feature bp (integer)
  {fs}_dmax__{feat}           Max window density (0-1 raw fraction)
  {fs}_dmin__{feat}           Min window density (0-1)
  {fs}_dmedian__{feat}        Median window density (0-1)
  {fs}_dfirst__{feat}         First 1kb density (0-1)
  {fs}_dlast__{feat}          Last 1kb density (0-1)
  {fs}_dterminal__{feat}      max(first, last) (0-1)
  {fs}_dterminal_min__{feat}  min(first, last) (0-1)
  {fs}_max_block_bp__{feat}   Longest contiguous block (bp)

Per featureset:
  {fs}_total_bp               Total annotated bp for this read in this featureset

Interspersion (from telomere_region or region featureset):
  interspersion_total, interspersion_can_ncan, interspersion_tel_sat, interspersion_arm_tel

Optional alignment stats (when --readnames-dir provided):
  sequencing_approach, n_alignments, n_secondary, n_supplementary,
  primary_mapq, primary_de, primary_align_len, primary_align_fraction,
  total_align_len, total_align_fraction

Usage:
  python KaryoScope_sequence_annotate.py \\
    --prefix analysis_output_prefix \\
    --bed-dir results \\
    --output sequence_annotations.tsv

  # With alignment stats:
  python KaryoScope_sequence_annotate.py \\
    --prefix analysis_output_prefix \\
    --bed-dir results \\
    --readnames-dir /path/to/samples \\
    --output sequence_annotations.tsv
"""

import argparse
import fnmatch
import gzip
import math
import os
import sys

# Capture original command line for logging
_original_command = ' '.join(sys.argv)

import pandas as pd

# Cache argparse defaults (before parse_args modifies them)
_argparse_defaults = None


# ---------------------------------------------------------------------------
# BED loading / discovery (from KaryoScope_cluster_annotate.py)
# ---------------------------------------------------------------------------

def load_bed_file(filepath):
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
                length = end - start
                records.append({
                    'read': read,
                    'start': start,
                    'end': end,
                    'feature': feature,
                    'length': length
                })
    return pd.DataFrame(records)


def find_featureset_beds(bed_dirs, samples, featuresets, database="KS_human_CHM13", smoothness="smoothed"):
    """Find BED files for each featureset for each sample, searching multiple directories."""
    beds_by_featureset = {fs: [] for fs in featuresets}

    for sample in samples:
        for fs in featuresets:
            found = False
            for bed_dir in bed_dirs:
                base_path = f"{bed_dir}/{sample}/telogator/1/KaryoScope/{database}"

                # Try different naming patterns (nested directory structure)
                patterns = [
                    f"{base_path}/{sample}.telogator.1.{database}.{fs}.{smoothness}.KaryoScope.bed",
                    f"{base_path}/{sample}.telogator.1.{database}.{fs}.{smoothness}.bed",
                    f"{base_path}/{sample}.telogator.1.{database}.{fs}.{smoothness}.features.bed",
                    f"{base_path}/{sample}.telogator.1.{database}.{fs}.{smoothness}.merged.bed",
                ]
                # Also try flat directory (files directly in bed_dir)
                patterns += [
                    f"{bed_dir}/{sample}.telogator.1.{database}.{fs}.{smoothness}.KaryoScope.bed",
                    f"{bed_dir}/{sample}.telogator.1.{database}.{fs}.{smoothness}.bed",
                    f"{bed_dir}/{sample}.telogator.1.{database}.{fs}.{smoothness}.features.bed",
                    f"{bed_dir}/{sample}.telogator.1.{database}.{fs}.{smoothness}.merged.bed",
                ]

                for pattern in patterns:
                    if os.path.exists(pattern):
                        beds_by_featureset[fs].append(pattern)
                        found = True
                        break
                    elif os.path.exists(pattern + '.gz'):
                        beds_by_featureset[fs].append(pattern + '.gz')
                        found = True
                        break
                if found:
                    break

    return beds_by_featureset


# ---------------------------------------------------------------------------
# Per-read feature fractions and thresholds (from KaryoScope_cluster_annotate.py)
# ---------------------------------------------------------------------------

def compute_read_feature_fractions(bed_df):
    """Compute per-read feature coverage fractions.

    Returns a DataFrame: read x feature matrix of coverage fractions.
    """
    # Total bp per read
    read_totals = bed_df.groupby('read')['length'].sum()
    # Per read x feature bp
    read_feature_bp = bed_df.groupby(['read', 'feature'])['length'].sum().unstack(fill_value=0)
    # Divide by total
    fractions = read_feature_bp.div(read_totals, axis=0)
    return fractions


def compute_adaptive_thresholds(fractions, min_thresh=0.001, max_thresh=0.05):
    """Compute adaptive significance thresholds per feature.

    For each feature: threshold = clamp(median_nonzero / 3, min_thresh, max_thresh)
    """
    thresholds = {}
    for feature in fractions.columns:
        nonzero = fractions[feature][fractions[feature] > 0]
        if len(nonzero) == 0:
            thresholds[feature] = min_thresh
        else:
            med = nonzero.median()
            thresholds[feature] = max(min_thresh, min(max_thresh, med / 3))
    return thresholds


# ---------------------------------------------------------------------------
# Block length helper (from KaryoScope_cluster_annotate.py)
# ---------------------------------------------------------------------------

BLOCK_GAP_TOL = 100  # merge blocks separated by <=100 bp gaps


def _max_block_length(coverage, gap_tol=BLOCK_GAP_TOL):
    """Longest contiguous block of 1s in coverage, merging gaps <= gap_tol bp."""
    import numpy as np
    diffs = np.diff(np.concatenate([[0], coverage, [0]]))
    run_starts = np.where(diffs == 1)[0]
    run_ends = np.where(diffs == -1)[0]
    if len(run_starts) == 0:
        return 0
    # Merge runs separated by small gaps
    merged_starts = [run_starts[0]]
    merged_ends = [run_ends[0]]
    for i in range(1, len(run_starts)):
        if run_starts[i] - merged_ends[-1] <= gap_tol:
            merged_ends[-1] = run_ends[i]
        else:
            merged_starts.append(run_starts[i])
            merged_ends.append(run_ends[i])
    return int(max(me - ms for ms, me in zip(merged_starts, merged_ends)))


# ---------------------------------------------------------------------------
# Feature classification constants for interspersion (from KaryoScope_cluster_annotate.py)
# ---------------------------------------------------------------------------

_SATELLITE_LAYER1 = frozenset({
    'hsat3', 'hsat1A', 'hsat2', 'hsat1B', 'active', 'censat', 'bsat',
    'monomeric', 'gsat', 'hor_multigroup1', 'hsat_multigroup1',
    'hsat1_multigroup1', 'asat_multigroup1',
})
_LAYER2_CANONICAL = frozenset({'canonical_telomere'})
_LAYER2_NONCANONICAL = frozenset({'noncanonical_telomere'})
_LAYER2_ITS_TAR1 = frozenset({'ITS', 'TAR1'})
_CT_LAYER1 = frozenset({'ct'})
_ARM_LAYER1 = frozenset({'p_arm', 'q_arm', 'arm_multigroup1'})


def classify_bed_feature(feature):
    """Classify a BED feature (single-layer or 2-layer) into a category.

    Handles both formats:
      - Single-layer: 'canonical_telomere', 'ct', 'active'
      - Two-layer:    'ct:nonsubtelomeric', 'arm_multigroup1:canonical_telomere'

    Priority:
    1. Satellite layer-1 always wins
    2. Layer-2 matches override ct/arm layer-1
    3. Layer-1 matches for canonical/noncanonical/ITS_TAR1 (single-layer BEDs)
    4. ct/arm layer-1 only when no layer-2 match
    5. Everything else -> other
    """
    parts = feature.split(':', 1)
    layer1 = parts[0]
    layer2 = parts[1] if len(parts) > 1 else ''

    if layer1 in _SATELLITE_LAYER1:
        return 'satellite'
    if layer2 in _LAYER2_CANONICAL:
        return 'canonical'
    if layer2 in _LAYER2_NONCANONICAL:
        return 'noncanonical'
    if layer2 in _LAYER2_ITS_TAR1:
        return 'ITS_TAR1'
    # Single-layer: check layer-1 against telomere-type sets
    if layer1 in _LAYER2_CANONICAL:
        return 'canonical'
    if layer1 in _LAYER2_NONCANONICAL:
        return 'noncanonical'
    if layer1 in _LAYER2_ITS_TAR1:
        return 'ITS_TAR1'
    if layer1 in _CT_LAYER1:
        return 'ct'
    if layer1 in _ARM_LAYER1:
        return 'arm'
    return 'other'


# ---------------------------------------------------------------------------
# Per-read window densities (bulk) (from KaryoScope_cluster_annotate.py)
# ---------------------------------------------------------------------------

def compute_per_read_window_densities_bulk(read_ids, bed_df, window_size=1000):
    """Compute per-read 1kb window density statistics for each feature.

    Same computation as compute_cluster_window_densities, but returns per-read
    stats instead of cluster medians.

    Returns:
        dict of {read_id: {feature: {'max': v, 'min': v, 'median': v,
                 'first': v, 'last': v, 'terminal': v, 'terminal_min': v,
                 'max_block': v}}}
        Values are raw fractions (0-1) except max_block which is in bp.
    """
    import numpy as np

    read_bed = bed_df[bed_df['read'].isin(read_ids)]
    if len(read_bed) == 0:
        return {}

    result = {}

    for read_id, read_records in read_bed.groupby('read'):
        read_start = read_records['start'].min()
        read_end = read_records['end'].max()
        span = read_end - read_start

        if span <= 0:
            continue

        read_result = {}

        if span < window_size:
            for feat, group in read_records.groupby('feature'):
                events = np.zeros(span + 1, dtype=np.int32)
                starts = group['start'].values - read_start
                ends = np.minimum(group['end'].values - read_start, span)
                np.add.at(events, starts, 1)
                np.add.at(events, ends, -1)
                coverage_arr = (np.cumsum(events[:span]) > 0).astype(np.int32)
                frac = coverage_arr.sum() / span
                read_result[feat] = {
                    'max': frac, 'min': frac, 'median': frac,
                    'first': frac, 'last': frac,
                    'terminal': frac, 'terminal_min': frac,
                    'max_block': _max_block_length(coverage_arr),
                }
        else:
            for feat, group in read_records.groupby('feature'):
                events = np.zeros(span + 1, dtype=np.int32)
                starts = group['start'].values - read_start
                ends = np.minimum(group['end'].values - read_start, span)
                np.add.at(events, starts, 1)
                np.add.at(events, ends, -1)
                coverage = (np.cumsum(events[:span]) > 0).astype(np.int32)

                cumsum = np.empty(len(coverage) + 1, dtype=np.int64)
                cumsum[0] = 0
                np.cumsum(coverage, out=cumsum[1:])
                window_sums = cumsum[window_size:] - cumsum[:-window_size]

                first_val = coverage[:window_size].sum() / window_size
                last_val = coverage[-window_size:].sum() / window_size
                read_result[feat] = {
                    'max': window_sums.max() / window_size,
                    'min': window_sums.min() / window_size,
                    'median': float(np.median(window_sums)) / window_size,
                    'first': first_val,
                    'last': last_val,
                    'terminal': max(first_val, last_val),
                    'terminal_min': min(first_val, last_val),
                    'max_block': _max_block_length(coverage),
                }

        result[read_id] = read_result

    return result


# ---------------------------------------------------------------------------
# Per-read interspersion (bulk) (from KaryoScope_cluster_annotate.py)
# ---------------------------------------------------------------------------

def compute_per_read_interspersion_bulk(read_ids, bed_df):
    """Compute per-read typed interspersion rates (transitions per kb).

    Same computation as compute_cluster_interspersion, but returns per-read
    stats instead of cluster medians.

    Returns:
        dict of {read_id: {'total': v, 'can_ncan': v, 'tel_sat': v, 'arm_tel': v}}
    """
    read_bed = bed_df[bed_df['read'].isin(read_ids)]
    if len(read_bed) == 0:
        return {}

    result = {}

    for read_id, read_records in read_bed.groupby('read'):
        records = read_records.sort_values('start')
        span_kb = (records['end'].max() - records['start'].min()) / 1000
        if span_kb <= 0:
            result[read_id] = {'total': 0.0, 'can_ncan': 0.0, 'tel_sat': 0.0, 'arm_tel': 0.0}
            continue

        categories = [classify_bed_feature(f) for f in records['feature']]

        total = sum(
            1 for i in range(1, len(categories))
            if categories[i] != categories[i - 1]
        )

        filtered = [c for c in categories if c != 'other']
        can_ncan = tel_sat = arm_tel = 0
        for i in range(1, len(filtered)):
            if filtered[i] == filtered[i - 1]:
                continue
            pair = frozenset({filtered[i - 1], filtered[i]})
            if pair == frozenset({'canonical', 'noncanonical'}):
                can_ncan += 1
            if 'satellite' in pair and pair & {'canonical', 'noncanonical'}:
                tel_sat += 1
            if 'arm' in pair and pair & {'canonical', 'noncanonical', 'ITS_TAR1'}:
                arm_tel += 1

        result[read_id] = {
            'total': round(total / span_kb, 2),
            'can_ncan': round(can_ncan / span_kb, 2),
            'tel_sat': round(tel_sat / span_kb, 2),
            'arm_tel': round(arm_tel / span_kb, 2),
        }

    return result


# ---------------------------------------------------------------------------
# New: per-read feature bp counts
# ---------------------------------------------------------------------------

def compute_per_read_feature_bp(bed_df):
    """Compute per-read per-feature bp counts and total bp.

    Returns:
        tuple: (read_feature_bp DataFrame with read x feature bp counts,
                read_total_bp Series with total bp per read)
    """
    read_feature_bp = bed_df.groupby(['read', 'feature'])['length'].sum().unstack(fill_value=0)
    read_total_bp = bed_df.groupby('read')['length'].sum()
    return read_feature_bp, read_total_bp


# ---------------------------------------------------------------------------
# Readnames / stats loading (from KaryoScope_annotate_sequences.py)
# ---------------------------------------------------------------------------

def load_readnames(readnames_dir, samples):
    """Load readnames.txt files for all samples.

    Args:
        readnames_dir: Base directory containing sample folders
        samples: List of sample names

    Returns:
        DataFrame with columns: read, sequencing_approach
    """
    all_readnames = []

    for sample in samples:
        readnames_file = os.path.join(readnames_dir, sample, "telogator", f"{sample}.readnames.txt")

        if not os.path.exists(readnames_file):
            raise FileNotFoundError(f"Readnames file not found: {readnames_file}")

        df = pd.read_csv(readnames_file, sep='\t', header=None, names=['read', 'sequencing_approach'])
        print(f"  {sample}: {len(df)} reads from readnames.txt")
        all_readnames.append(df)

    combined = pd.concat(all_readnames, ignore_index=True)

    # Check for duplicates
    duplicates = combined[combined.duplicated(subset=['read'], keep=False)]
    if len(duplicates) > 0:
        raise ValueError(f"Found {len(duplicates)} duplicate read names in readnames files")

    return combined


def load_stats(readnames_dir, samples, reference="CHM13"):
    """Load stats.tsv files for all samples.

    Computes per-read statistics including:
    - Primary alignment stats (mapq, de, align_len, align_fraction)
    - Alignment counts (total, secondary, supplementary)
    - Total aligned bases and fraction across all alignments

    Args:
        readnames_dir: Base directory containing sample folders
        samples: List of sample names
        reference: Reference genome name (default: CHM13)

    Returns:
        DataFrame with mapping statistics (one row per read)
    """
    all_stats = []

    for sample in samples:
        stats_file = os.path.join(readnames_dir, sample, "telogator", "aligned", f"{sample}.{reference}.stats.tsv")

        if not os.path.exists(stats_file):
            raise FileNotFoundError(f"Stats file not found: {stats_file}")

        df = pd.read_csv(stats_file, sep='\t')

        # Rename readname to read for consistency
        if 'readname' in df.columns:
            df = df.rename(columns={'readname': 'read'})

        total_rows = len(df)

        # Compute per-read aggregate statistics from ALL alignments
        if all(col in df.columns for col in ['is_primary', 'is_not_supplementary', 'align_len', 'read_len']):
            # Count alignments by type for each read
            agg_stats = df.groupby('read').agg(
                n_alignments=('read', 'count'),
                n_secondary=('is_primary', lambda x: (~x).sum()),
                n_supplementary=('is_not_supplementary', lambda x: (~x).sum()),
                read_len=('read_len', 'first'),  # Same for all alignments of a read
                max_mapq=('mapq', 'max'),
                mean_de=('de', 'mean'),
            ).reset_index()

            # For total aligned bases, only count non-secondary alignments
            # (primary + supplementary cover non-overlapping portions of the read)
            # Secondary alignments are alternative mappings for the same read portion
            # Non-secondary = primary OR supplementary
            # is_primary=True for primary, is_not_supplementary=False for supplementary
            non_secondary = df[(df['is_primary'] == True) | (df['is_not_supplementary'] == False)]
            non_secondary_align = non_secondary.groupby('read').agg(
                total_align_len=('align_len', 'sum'),
            ).reset_index()

            agg_stats = agg_stats.merge(non_secondary_align, on='read', how='left')
            agg_stats['total_align_len'] = agg_stats['total_align_len'].fillna(0)

            # Calculate total alignment fraction (primary + supplementary only)
            agg_stats['total_align_fraction'] = agg_stats['total_align_len'] / agg_stats['read_len']

            # Get primary non-supplementary alignment stats (the "main" alignment)
            primary_df = df[
                (df['is_primary'] == True) &
                (df['is_not_supplementary'] == True)
            ][['read', 'mapq', 'de', 'align_len', 'align_fraction', 'is_mapped']].copy()
            primary_df = primary_df.rename(columns={
                'mapq': 'primary_mapq',
                'de': 'primary_de',
                'align_len': 'primary_align_len',
                'align_fraction': 'primary_align_fraction',
                'is_mapped': 'is_mapped'
            })

            # Merge aggregate stats with primary alignment stats
            result_df = agg_stats.merge(primary_df, on='read', how='left')

            print(f"  {sample}: {len(result_df)} reads with alignment stats (from {total_rows} total alignments)")
            all_stats.append(result_df)
        else:
            # Fallback: just filter to primary non-supplementary
            filtered_df = df[
                (df['is_primary'] == True) &
                (df['is_not_supplementary'] == True)
            ].copy()

            # Add placeholder columns so merges work consistently
            filtered_df['n_alignments'] = 1
            filtered_df['n_secondary'] = 0
            filtered_df['n_supplementary'] = 0
            if 'align_len' in filtered_df.columns and 'read_len' in filtered_df.columns:
                filtered_df['total_align_len'] = filtered_df['align_len']
                filtered_df['total_align_fraction'] = filtered_df['align_len'] / filtered_df['read_len']
            if 'mapq' in filtered_df.columns:
                filtered_df['primary_mapq'] = filtered_df['mapq']
            if 'de' in filtered_df.columns:
                filtered_df['primary_de'] = filtered_df['de']
            if 'align_len' in filtered_df.columns:
                filtered_df['primary_align_len'] = filtered_df['align_len']
            if 'align_fraction' in filtered_df.columns:
                filtered_df['primary_align_fraction'] = filtered_df['align_fraction']

            print(f"  {sample}: {len(filtered_df)} primary alignments (columns for aggregate stats not found)")
            all_stats.append(filtered_df)

    combined = pd.concat(all_stats, ignore_index=True)

    # Verify no duplicates - should be exactly one row per read
    duplicates = combined[combined.duplicated(subset=['read'], keep=False)]
    if len(duplicates) > 0:
        dup_reads = duplicates['read'].unique()[:5]
        raise ValueError(
            f"Found {len(duplicates)} duplicate reads in stats after filtering. "
            f"Expected exactly one row per read. Examples: {list(dup_reads)}"
        )

    return combined


# ---------------------------------------------------------------------------
# TeeLogger (from KaryoScope_cluster_annotate.py)
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Parameter printing
# ---------------------------------------------------------------------------

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
        ("bed-dir", str(args.bed_dir), None),
        ("featuresets", _fmt(args.featuresets, "featuresets"), None),
        ("database", _fmt(args.database, "database"), None),
        ("smoothness", _fmt(args.smoothness, "smoothness"), None),
        ("window-size", _fmt(args.window_size, "window_size"), None),
        ("readnames-dir", _fmt(args.readnames_dir, "readnames_dir"), None),
        ("reference", _fmt(args.reference, "reference"), None),
        ("log-file", _fmt(args.log_file, "log_file"), None),
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


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compute per-read feature annotations from BED files",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--prefix", required=True,
                        help="Analysis prefix (finds {prefix}.read_assignments.tsv)")
    parser.add_argument("--bed-dir", dest="bed_dir", required=True,
                        help="Comma-separated BED directories")
    parser.add_argument("--output", "-o", required=True,
                        help="Output TSV path")
    parser.add_argument("--featuresets", default="region,subtelomeric,chromosome,acrocentric,repeat,gene",
                        help="Comma-separated featuresets (default: region,subtelomeric,chromosome,acrocentric,repeat,gene)")
    parser.add_argument("--database", default="KS_human_CHM13",
                        help="Database name (default: KS_human_CHM13)")
    parser.add_argument("--smoothness", default="smoothed",
                        help="BED smoothness (default: smoothed)")
    parser.add_argument("--window-size", dest="window_size", type=int, default=1000,
                        help="Window size in bp (default: 1000)")
    parser.add_argument("--readnames-dir", dest="readnames_dir", default=None,
                        help="Base directory for readnames.txt + stats.tsv (enables alignment stats)")
    parser.add_argument("--reference", default="CHM13",
                        help="Reference genome name for stats.tsv files (default: CHM13)")
    parser.add_argument("--log-file", dest="log_file", default=True,
                        type=lambda x: x.lower() not in ('false', '0', 'no'),
                        help="Save .log file (default: True)")

    global _argparse_defaults
    _argparse_defaults = {}
    for action in parser._actions:
        if action.dest != 'help' and action.default is not None:
            _argparse_defaults[action.dest] = action.default

    args = parser.parse_args()

    # 1. Set up logging
    if args.log_file:
        if args.output.endswith('.tsv'):
            log_path = args.output[:-4] + '.log'
        else:
            log_path = args.output + '.log'
        sys.stdout = TeeLogger(log_path)

    print("=" * 60)
    print("KaryoScope Sequence Annotation")
    print("=" * 60)
    _print_params_and_command(args)

    # 2. Load read assignments -> get read IDs and sample names
    read_assignments_file = f"{args.prefix}.read_assignments.tsv"
    if not os.path.exists(read_assignments_file):
        read_assignments_file = f"{args.prefix}.sequence_assignments.tsv"

    if not os.path.exists(read_assignments_file):
        print(f"ERROR: No read/sequence assignments file found for prefix: {args.prefix}")
        sys.exit(1)

    print(f"\nLoading read assignments: {read_assignments_file}")
    assignments = pd.read_csv(read_assignments_file, sep='\t')

    if 'sequence' not in assignments.columns and 'read' in assignments.columns:
        assignments.rename(columns={'read': 'sequence'}, inplace=True)

    all_read_ids = set(assignments['sequence'].tolist())
    samples = assignments['sample'].unique().tolist()
    print(f"  Total reads: {len(all_read_ids)}")
    print(f"  Samples: {len(samples)}")

    # 3. Find and load BED files per featureset (same logic as cluster_annotate)
    featuresets = [fs.strip() for fs in args.featuresets.split(',')]
    bed_dirs = [d.strip() for d in args.bed_dir.split(',')]

    print(f"\nFeaturesets: {featuresets}")
    print(f"Finding BED files in: {bed_dirs}")
    beds_by_featureset = find_featureset_beds(bed_dirs, samples, featuresets, args.database, args.smoothness)

    bed_data = {}
    for fs in featuresets:
        bed_files = beds_by_featureset[fs]
        if not bed_files:
            print(f"  WARNING: No BED files found for featureset '{fs}'")
            continue
        print(f"\n  Loading {fs}: {len(bed_files)} files")
        all_dfs = []
        for bf in bed_files:
            df = load_bed_file(bf)
            all_dfs.append(df)
        bed_data[fs] = pd.concat(all_dfs, ignore_index=True)
        print(f"    Total records: {len(bed_data[fs])}")

    if not bed_data:
        print("ERROR: No BED data loaded")
        sys.exit(1)

    # 4. For each featureset, compute per-read metrics
    # Build columns per featureset, then concat all at the end to avoid fragmentation
    sorted_read_ids = sorted(all_read_ids)
    master = pd.DataFrame({'sequence': sorted_read_ids})
    featureset_frames = []  # collect DataFrames to concat at end

    threshold_rows = []  # For thresholds TSV

    window_size = args.window_size

    for fs in featuresets:
        if fs not in bed_data:
            continue

        print(f"\n{'=' * 40}")
        print(f"Processing featureset: {fs}")
        print(f"{'=' * 40}")

        bed_df = bed_data[fs]

        # Compute per-read feature bp and fractions
        read_feature_bp, read_total_bp = compute_per_read_feature_bp(bed_df)
        fractions = compute_read_feature_fractions(bed_df)
        thresholds = compute_adaptive_thresholds(fractions)
        feature_names_sorted = sorted(thresholds.keys())

        # Print threshold summary
        print(f"\n  Feature thresholds for {fs}:")
        for feat in feature_names_sorted:
            nonzero = fractions[feat][fractions[feat] > 0]
            med = nonzero.median() if len(nonzero) > 0 else 0
            n_nz = len(nonzero)
            n_tot = len(fractions)
            print(f"    {feat}: {thresholds[feat]*100:.2f}% (median {med*100:.1f}%, n_nonzero={n_nz}/{n_tot})")
            threshold_rows.append({
                'featureset': fs,
                'feature': feat,
                'threshold': thresholds[feat],
                'median_nonzero': med if n_nz > 0 else 0,
                'n_nonzero': n_nz,
                'n_total': n_tot,
            })

        # Compute per-read window densities for ALL reads at once
        print(f"  Computing window densities (window_size={window_size})...")
        per_read_densities = compute_per_read_window_densities_bulk(all_read_ids, bed_df, window_size)

        # Build all columns for this featureset into a dict, then create DataFrame at once
        fs_cols = {}

        # {fs}_total_bp
        total_bp_reindexed = read_total_bp.reindex(sorted_read_ids).fillna(0).astype(int)
        fs_cols[f'{fs}_total_bp'] = total_bp_reindexed.values

        for feat in feature_names_sorted:
            # bp count
            if feat in read_feature_bp.columns:
                bp_reindexed = read_feature_bp[feat].reindex(sorted_read_ids).fillna(0).astype(int)
            else:
                bp_reindexed = pd.Series(0, index=sorted_read_ids)
            fs_cols[f'{fs}_bp__{feat}'] = bp_reindexed.values

            # fraction
            if feat in fractions.columns:
                frac_reindexed = fractions[feat].reindex(sorted_read_ids).fillna(0)
            else:
                frac_reindexed = pd.Series(0.0, index=sorted_read_ids)
            fs_cols[f'{fs}_frac__{feat}'] = frac_reindexed.values

            # Window density stats from per_read_densities
            dmax_vals = {}
            dmin_vals = {}
            dmedian_vals = {}
            dfirst_vals = {}
            dlast_vals = {}
            dterminal_vals = {}
            dterminal_min_vals = {}
            max_block_vals = {}

            for read_id, feat_dict in per_read_densities.items():
                if feat in feat_dict:
                    stats = feat_dict[feat]
                    dmax_vals[read_id] = stats['max']
                    dmin_vals[read_id] = stats['min']
                    dmedian_vals[read_id] = stats['median']
                    dfirst_vals[read_id] = stats['first']
                    dlast_vals[read_id] = stats['last']
                    dterminal_vals[read_id] = stats['terminal']
                    dterminal_min_vals[read_id] = stats['terminal_min']
                    max_block_vals[read_id] = stats['max_block']

            seq_series = master['sequence']
            fs_cols[f'{fs}_dmax__{feat}'] = seq_series.map(dmax_vals).fillna(0).values
            fs_cols[f'{fs}_dmin__{feat}'] = seq_series.map(dmin_vals).fillna(0).values
            fs_cols[f'{fs}_dmedian__{feat}'] = seq_series.map(dmedian_vals).fillna(0).values
            fs_cols[f'{fs}_dfirst__{feat}'] = seq_series.map(dfirst_vals).fillna(0).values
            fs_cols[f'{fs}_dlast__{feat}'] = seq_series.map(dlast_vals).fillna(0).values
            fs_cols[f'{fs}_dterminal__{feat}'] = seq_series.map(dterminal_vals).fillna(0).values
            fs_cols[f'{fs}_dterminal_min__{feat}'] = seq_series.map(dterminal_min_vals).fillna(0).values
            fs_cols[f'{fs}_max_block_bp__{feat}'] = seq_series.map(max_block_vals).fillna(0).values

        featureset_frames.append(pd.DataFrame(fs_cols, index=master.index))

    # Concat all featureset columns at once to avoid fragmentation
    if featureset_frames:
        master = pd.concat([master] + featureset_frames, axis=1)

    # 5. Compute interspersion using telomere_region or region BED
    interspersion_fs = 'telomere_region' if 'telomere_region' in bed_data else 'region' if 'region' in bed_data else None
    if interspersion_fs:
        print(f"\nComputing interspersion from '{interspersion_fs}'...")
        per_read_interspersion = compute_per_read_interspersion_bulk(all_read_ids, bed_data[interspersion_fs])

        inter_cols = {}
        for key in ['total', 'can_ncan', 'tel_sat', 'arm_tel']:
            vals = {}
            for read_id, inter_dict in per_read_interspersion.items():
                vals[read_id] = inter_dict[key]
            inter_cols[f'interspersion_{key}'] = master['sequence'].map(vals).fillna(0).values
        master = pd.concat([master, pd.DataFrame(inter_cols, index=master.index)], axis=1)
    else:
        print("\n  WARNING: No telomere_region or region featureset found — skipping interspersion")

    # 6. Optional: alignment stats from readnames + stats
    if args.readnames_dir:
        print(f"\n{'=' * 40}")
        print("Loading alignment statistics")
        print(f"{'=' * 40}")

        print(f"\nLoading readnames files from: {args.readnames_dir}")
        readnames = load_readnames(args.readnames_dir, samples)
        print(f"  Total reads in readnames: {len(readnames)}")

        print(f"\nLoading stats files from: {args.readnames_dir}")
        stats = load_stats(args.readnames_dir, samples, args.reference)
        print(f"  Total reads in stats: {len(stats)}")

        # Join readnames onto master (left join to keep all reads)
        master = master.merge(readnames[['read', 'sequencing_approach']],
                              left_on='sequence', right_on='read', how='left')
        if 'read' in master.columns and 'sequence' in master.columns:
            master.drop(columns=['read'], inplace=True)

        # Join stats onto master
        stats_cols = ['read', 'n_alignments', 'n_secondary', 'n_supplementary',
                      'primary_mapq', 'primary_de', 'primary_align_len', 'primary_align_fraction',
                      'total_align_len', 'total_align_fraction']
        stats_subset = stats[[c for c in stats_cols if c in stats.columns]].copy()
        master = master.merge(stats_subset, left_on='sequence', right_on='read', how='left')
        if 'read' in master.columns and 'sequence' in master.columns:
            master.drop(columns=['read'], inplace=True)

        # Report join results
        n_with_approach = master['sequencing_approach'].notna().sum() if 'sequencing_approach' in master.columns else 0
        print(f"\n  Reads with alignment stats: {n_with_approach}/{len(master)}")

    # 7. Save outputs
    output_prefix = args.output[:-4] if args.output.endswith('.tsv') else args.output

    # Main TSV
    master.to_csv(args.output, sep='\t', index=False)
    print(f"\nSaved sequence annotations to: {args.output}")
    print(f"  Rows: {len(master)}")
    print(f"  Columns: {len(master.columns)}")

    # Thresholds TSV
    if threshold_rows:
        thresh_df = pd.DataFrame(threshold_rows)
        thresh_path = f"{output_prefix}.adaptive_thresholds.tsv"
        thresh_df.to_csv(thresh_path, sep='\t', index=False)
        print(f"  Saved adaptive thresholds to: {thresh_path}")

    # Summary
    print(f"\n{'=' * 60}")
    print("Summary")
    print("=" * 60)
    print(f"Total reads annotated: {len(master)}")
    print(f"Featuresets processed: {[fs for fs in featuresets if fs in bed_data]}")

    if 'sequencing_approach' in master.columns:
        approaches = master['sequencing_approach'].value_counts()
        print(f"Sequencing approaches: {approaches.to_dict()}")


if __name__ == "__main__":
    main()
