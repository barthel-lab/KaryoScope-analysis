# KaryoScope BED Merger
# Merges multiple BED featuresets by position overlay, creating combined feature labels
#
# Usage:
# python KaryoScope_merge_beds.py \
#   --bed sample.region.bed.gz sample.chromosome.bed.gz \
#   --output sample.merged.bed.gz
#
# Three or more featuresets:
# python KaryoScope_merge_beds.py \
#   --bed sample.region.bed.gz sample.chromosome.bed.gz sample.repeat.bed.gz \
#   --output sample.merged.bed.gz
#
# The output can then be used with KaryoScope_cluster_analysis.py:
# python KaryoScope_cluster_analysis.py \
#   --bed pre.merged.bed.gz post.merged.bed.gz ...

import argparse
import gzip
import pandas as pd
import sys

parser = argparse.ArgumentParser(
    description="Merge multiple BED featuresets by position overlay.",
    formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument("--bed", required=True, nargs='+',
                    help="BED files to merge (2 or more). Features are combined in order.")
parser.add_argument("--output", "-o", required=True,
                    help="Output merged BED file (use .gz extension for compression)")
parser.add_argument("--separator", "-s", default=":",
                    help="Separator for combined feature labels (default: ':')")
parser.add_argument("--reduce-features", dest="reduce_features", type=int, default=None,
                    help="Reduce to top N most frequent combined features.\n"
                         "Less frequent features are collapsed to 'other'. (default: no reduction)")
parser.add_argument("--feature-filter", dest="feature_filters", nargs='*', default=[],
                    metavar="INDEX:KEEP:COLLAPSE",
                    help="Filter features for specific BED files before merging.\n"
                         "Format: BED_INDEX:KEEP_PATTERNS:COLLAPSE_LABEL\n"
                         "  BED_INDEX: 1-based index of BED file in --bed list\n"
                         "  KEEP_PATTERNS: Comma-separated prefixes to keep (e.g., 'p_arm,q_arm')\n"
                         "  COLLAPSE_LABEL: Label for collapsed features (e.g., 'satellite')\n"
                         "Example: --feature-filter 3:p_arm,q_arm:satellite")
parser.add_argument("--telomere-satellite-merge", dest="telomere_satellite_merge", action="store_true",
                    help="Priority merge mode: keep telomere features (canonical_telomere, noncanonical_telomere,\n"
                         "TAR1, ITS) from BED1 (subtelomeric), fill gaps with features from BED2 (satellite/region).\n"
                         "Requires exactly 2 BED files.")

args = parser.parse_args()

if len(args.bed) < 2:
    print("Error: At least 2 BED files required for merging", file=sys.stderr)
    sys.exit(1)


def load_bed_file(filepath):
    """Load a BED file into a DataFrame."""
    open_func = gzip.open if filepath.endswith('.gz') else open
    mode = 'rt' if filepath.endswith('.gz') else 'r'

    records = []
    malformed_count = 0
    with open_func(filepath, mode) as f:
        for line_num, line in enumerate(f, 1):
            parts = line.strip().split('\t')
            if len(parts) < 4:
                malformed_count += 1
                continue
            try:
                records.append({
                    'read': parts[0],
                    'start': int(parts[1]),
                    'end': int(parts[2]),
                    'feature': parts[3]
                })
            except ValueError:
                malformed_count += 1
                continue

    if malformed_count > 0:
        print(f"  Warning: Skipped {malformed_count} malformed lines in {filepath}")

    df = pd.DataFrame(records)
    if df.empty:
        print(f"  Warning: No valid records in {filepath}")
    return df


def apply_feature_filter(df, keep_patterns, collapse_label):
    """
    Filter features: keep those matching patterns, collapse others.

    Args:
        df: DataFrame with 'feature' column
        keep_patterns: List of prefixes to keep (e.g., ['p_arm', 'q_arm'])
        collapse_label: Label for non-matching features (e.g., 'satellite')

    Returns:
        DataFrame with filtered features
    """
    def filter_feature(feat):
        for pattern in keep_patterns:
            if feat.startswith(pattern):
                return feat
        return collapse_label

    df = df.copy()
    original_unique = df['feature'].nunique()
    df['feature'] = df['feature'].apply(filter_feature)
    new_unique = df['feature'].nunique()
    print(f"    Filter applied: {original_unique} -> {new_unique} unique features")
    return df


# Priority features for telomere-satellite merge mode
TELOMERE_PRIORITY_FEATURES = {'canonical_telomere', 'noncanonical_telomere', 'TAR1', 'ITS'}


def telomere_satellite_merge(df_subtelo, df_satellite):
    """
    Priority merge: keep telomere features from subtelomeric BED,
    fill remaining positions with satellite features.

    Args:
        df_subtelo: DataFrame from subtelomeric BED file
        df_satellite: DataFrame from region/satellite BED file

    Returns:
        DataFrame with merged features (single feature per interval)
    """
    common_reads = set(df_subtelo['read'].unique()) & set(df_satellite['read'].unique())
    if not common_reads:
        print("  Warning: No common reads between BED files")
        return pd.DataFrame(columns=['read', 'start', 'end', 'feature'])

    print(f"  Common reads: {len(common_reads):,}")
    print(f"  Priority features: {', '.join(sorted(TELOMERE_PRIORITY_FEATURES))}")

    results = []
    for read in common_reads:
        # Get intervals for this read
        subtelo_intervals = df_subtelo[df_subtelo['read'] == read][['start', 'end', 'feature']].values
        satellite_intervals = df_satellite[df_satellite['read'] == read][['start', 'end', 'feature']].values

        # Extract priority intervals from subtelomeric
        priority_intervals = []
        for start, end, feature in subtelo_intervals:
            if feature in TELOMERE_PRIORITY_FEATURES:
                priority_intervals.append((int(start), int(end), feature))

        # Build a coverage map of priority positions
        priority_coverage = set()
        for start, end, _ in priority_intervals:
            for pos in range(start, end):
                priority_coverage.add(pos)

        # Add priority intervals to results
        for start, end, feature in priority_intervals:
            results.append({
                'read': read,
                'start': start,
                'end': end,
                'feature': feature
            })

        # Add satellite intervals for non-priority positions
        for start, end, feature in satellite_intervals:
            start, end = int(start), int(end)
            # Find segments not covered by priority
            seg_start = None
            for pos in range(start, end + 1):
                in_priority = pos in priority_coverage
                if pos == end or in_priority:
                    # End of non-priority segment
                    if seg_start is not None:
                        seg_end = pos
                        if seg_end > seg_start:
                            results.append({
                                'read': read,
                                'start': seg_start,
                                'end': seg_end,
                                'feature': feature
                            })
                        seg_start = None
                elif seg_start is None:
                    # Start of non-priority segment
                    seg_start = pos

    return pd.DataFrame(results)


def merge_two_beds(df1, df2, sep=":"):
    """Merge two BED DataFrames by position overlay."""
    try:
        import pyranges as pr
        return _merge_pyranges(df1, df2, sep)
    except ImportError:
        return _merge_pandas(df1, df2, sep)


def _merge_pyranges(df1, df2, sep):
    """Fast merge using pyranges join/intersect."""
    import pyranges as pr

    common_reads = set(df1['read'].unique()) & set(df2['read'].unique())
    if not common_reads:
        print("  Warning: No common reads between BED files")
        return pd.DataFrame(columns=['read', 'start', 'end', 'feature'])

    df1_f = df1[df1['read'].isin(common_reads)][['read', 'start', 'end', 'feature']].copy()
    df2_f = df2[df2['read'].isin(common_reads)][['read', 'start', 'end', 'feature']].copy()

    df1_f.columns = ['Chromosome', 'Start', 'End', 'Feature1']
    df2_f.columns = ['Chromosome', 'Start', 'End', 'Feature2']

    pr1 = pr.PyRanges(df1_f)
    pr2 = pr.PyRanges(df2_f)

    joined = pr1.join(pr2)

    if len(joined) == 0:
        return pd.DataFrame(columns=['read', 'start', 'end', 'feature'])

    result_df = joined.df
    result_df['overlap_start'] = result_df[['Start', 'Start_b']].max(axis=1)
    result_df['overlap_end'] = result_df[['End', 'End_b']].min(axis=1)
    result_df = result_df[result_df['overlap_end'] > result_df['overlap_start']]
    result_df['feature'] = result_df['Feature1'] + sep + result_df['Feature2']

    return pd.DataFrame({
        'read': result_df['Chromosome'],
        'start': result_df['overlap_start'].astype(int),
        'end': result_df['overlap_end'].astype(int),
        'feature': result_df['feature']
    })


def _merge_pandas(df1, df2, sep):
    """Fallback merge using pandas (slower)."""
    common_reads = set(df1['read'].unique()) & set(df2['read'].unique())
    if not common_reads:
        return pd.DataFrame(columns=['read', 'start', 'end', 'feature'])

    df1_f = df1[df1['read'].isin(common_reads)].copy()
    df2_f = df2[df2['read'].isin(common_reads)].copy()

    grouped1 = {read: grp[['start', 'end', 'feature']].values for read, grp in df1_f.groupby('read')}
    grouped2 = {read: grp[['start', 'end', 'feature']].values for read, grp in df2_f.groupby('read')}

    results = []
    for read in common_reads:
        intervals1 = grouped1[read]
        intervals2 = grouped2[read]

        for s1, e1, f1 in intervals1:
            for s2, e2, f2 in intervals2:
                overlap_start = max(s1, s2)
                overlap_end = min(e1, e2)
                if overlap_end > overlap_start:
                    results.append({
                        'read': read,
                        'start': overlap_start,
                        'end': overlap_end,
                        'feature': f"{f1}{sep}{f2}"
                    })

    return pd.DataFrame(results)


# --- Main processing ---

# Parse feature filters
feature_filter_map = {}
for ff in args.feature_filters:
    parts = ff.split(':')
    if len(parts) != 3:
        print(f"Error: Invalid feature filter format: {ff}", file=sys.stderr)
        print("  Expected format: BED_INDEX:KEEP_PATTERNS:COLLAPSE_LABEL", file=sys.stderr)
        sys.exit(1)
    try:
        bed_idx = int(parts[0])
    except ValueError:
        print(f"Error: BED_INDEX must be an integer: {parts[0]}", file=sys.stderr)
        sys.exit(1)
    keep_patterns = [p.strip() for p in parts[1].split(',')]
    collapse_label = parts[2]
    feature_filter_map[bed_idx] = (keep_patterns, collapse_label)
    print(f"Feature filter for BED{bed_idx}: keep [{', '.join(keep_patterns)}], collapse to '{collapse_label}'")

print(f"\nMerging {len(args.bed)} BED files...")
for i, bed_file in enumerate(args.bed, 1):
    print(f"  BED{i}: {bed_file}")

# Handle telomere-satellite merge mode
if args.telomere_satellite_merge:
    if len(args.bed) != 2:
        print("Error: --telomere-satellite-merge requires exactly 2 BED files", file=sys.stderr)
        print("  BED1: subtelomeric features", file=sys.stderr)
        print("  BED2: region/satellite features", file=sys.stderr)
        sys.exit(1)

    print("\n--- Telomere-Satellite Priority Merge Mode ---")
    df_subtelo = load_bed_file(args.bed[0])
    print(f"  Subtelomeric intervals: {len(df_subtelo):,}")

    df_satellite = load_bed_file(args.bed[1])
    print(f"  Satellite intervals: {len(df_satellite):,}")

    merged_df = telomere_satellite_merge(df_subtelo, df_satellite)
    print(f"\n  Merged intervals: {len(merged_df):,}")

    if merged_df.empty:
        print("Error: Merge resulted in no intervals", file=sys.stderr)
        sys.exit(1)

    # Skip to output (no feature reduction in this mode)
    merged_df = merged_df.sort_values(['read', 'start'])

    print(f"\nWriting to: {args.output}")
    open_func = gzip.open if args.output.endswith('.gz') else open
    mode = 'wt' if args.output.endswith('.gz') else 'w'

    with open_func(args.output, mode) as f:
        for _, row in merged_df.iterrows():
            f.write(f"{row['read']}\t{row['start']}\t{row['end']}\t{row['feature']}\n")

    print("Done!")
    sys.exit(0)

# Load first BED file
merged_df = load_bed_file(args.bed[0])
print(f"\n  BED1 intervals: {len(merged_df):,}")

# Apply feature filter if specified for BED1
if 1 in feature_filter_map:
    keep, collapse = feature_filter_map[1]
    merged_df = apply_feature_filter(merged_df, keep, collapse)

if merged_df.empty:
    print("Error: First BED file has no valid records", file=sys.stderr)
    sys.exit(1)

# Iteratively merge with remaining BED files
for i, bed_file in enumerate(args.bed[1:], 2):
    df_next = load_bed_file(bed_file)
    print(f"  BED{i} intervals: {len(df_next):,}")

    # Apply feature filter if specified for this BED
    if i in feature_filter_map:
        keep, collapse = feature_filter_map[i]
        df_next = apply_feature_filter(df_next, keep, collapse)

    if df_next.empty:
        print(f"  Warning: BED{i} has no valid records, skipping")
        continue

    merged_df = merge_two_beds(merged_df, df_next, args.separator)
    print(f"  After merge {i-1}: {len(merged_df):,} intervals")

    if merged_df.empty:
        print(f"  Warning: Merge resulted in no overlapping intervals")
        break

print(f"\nFinal merged intervals: {len(merged_df):,}")

# Count unique features
feature_counts = merged_df['feature'].value_counts()
print(f"Unique combined features: {len(feature_counts):,}")

# Optional feature reduction
if args.reduce_features is not None and len(feature_counts) > args.reduce_features:
    print(f"\n--- Reducing to top {args.reduce_features} features ---")
    top_features = set(feature_counts.head(args.reduce_features).index)

    merged_df['feature'] = merged_df['feature'].apply(
        lambda x: x if x in top_features else 'other'
    )

    new_feature_counts = merged_df['feature'].value_counts()
    print(f"  Features after reduction: {len(new_feature_counts):,}")
    other_count = new_feature_counts.get('other', 0)
    if other_count > 0:
        other_pct = other_count / len(merged_df) * 100
        print(f"  Intervals collapsed to 'other': {other_count:,} ({other_pct:.1f}%)")

# Sort by read, then start position
merged_df = merged_df.sort_values(['read', 'start'])

# Write output
print(f"\nWriting to: {args.output}")
open_func = gzip.open if args.output.endswith('.gz') else open
mode = 'wt' if args.output.endswith('.gz') else 'w'

with open_func(args.output, mode) as f:
    for _, row in merged_df.iterrows():
        f.write(f"{row['read']}\t{row['start']}\t{row['end']}\t{row['feature']}\n")

print("Done!")
