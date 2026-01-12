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


def subtract_intervals(interval, blockers):
    """
    Subtract blocker intervals from a single interval.
    Returns list of non-blocked sub-intervals.

    Args:
        interval: (start, end) tuple
        blockers: list of (start, end) tuples, sorted by start

    Returns:
        list of (start, end) tuples representing uncovered portions
    """
    result = []
    current_start, current_end = interval

    for block_start, block_end in blockers:
        if block_end <= current_start:
            # Blocker is entirely before current segment
            continue
        if block_start >= current_end:
            # Blocker is entirely after current segment
            break

        # There's overlap
        if block_start > current_start:
            # Gap before blocker
            result.append((current_start, block_start))

        # Move current_start past blocker
        current_start = max(current_start, block_end)

        if current_start >= current_end:
            break

    # Add remaining segment if any
    if current_start < current_end:
        result.append((current_start, current_end))

    return result


def telomere_satellite_merge(df_subtelo, df_satellite):
    """
    Priority merge: keep telomere features from subtelomeric BED,
    fill remaining positions with satellite features.

    Uses pyranges for fast interval subtraction when available.

    Args:
        df_subtelo: DataFrame from subtelomeric BED file
        df_satellite: DataFrame from region/satellite BED file

    Returns:
        DataFrame with merged features (single feature per interval)
    """
    try:
        import pyranges as pr
        return _telomere_satellite_merge_pyranges(df_subtelo, df_satellite)
    except ImportError:
        return _telomere_satellite_merge_pandas(df_subtelo, df_satellite)


def _telomere_satellite_merge_pyranges(df_subtelo, df_satellite):
    """Fast telomere-satellite merge using pyranges subtract."""
    import pyranges as pr

    common_reads = set(df_subtelo['read'].unique()) & set(df_satellite['read'].unique())
    if not common_reads:
        print("  Warning: No common reads between BED files")
        return pd.DataFrame(columns=['read', 'start', 'end', 'feature'])

    print(f"  Common reads: {len(common_reads):,}")
    print(f"  Priority features: {', '.join(sorted(TELOMERE_PRIORITY_FEATURES))}")

    # Filter to common reads
    df_subtelo_f = df_subtelo[df_subtelo['read'].isin(common_reads)].copy()
    df_satellite_f = df_satellite[df_satellite['read'].isin(common_reads)].copy()

    # Extract priority intervals
    priority_mask = df_subtelo_f['feature'].isin(TELOMERE_PRIORITY_FEATURES)
    df_priority = df_subtelo_f[priority_mask][['read', 'start', 'end', 'feature']].copy()
    df_priority['start'] = df_priority['start'].astype(int)
    df_priority['end'] = df_priority['end'].astype(int)

    print(f"  Priority intervals: {len(df_priority):,}")

    # Prepare satellite intervals
    df_satellite_f = df_satellite_f[['read', 'start', 'end', 'feature']].copy()
    df_satellite_f['start'] = df_satellite_f['start'].astype(int)
    df_satellite_f['end'] = df_satellite_f['end'].astype(int)

    # Convert to pyranges format (Chromosome = read)
    pr_priority = pr.PyRanges(df_priority.rename(columns={'read': 'Chromosome', 'start': 'Start', 'end': 'End', 'feature': 'Feature'}))
    pr_satellite = pr.PyRanges(df_satellite_f.rename(columns={'read': 'Chromosome', 'start': 'Start', 'end': 'End', 'feature': 'Feature'}))

    # Subtract priority regions from satellite
    pr_subtracted = pr_satellite.subtract(pr_priority)

    # Convert results back to DataFrame
    if len(pr_subtracted) > 0:
        df_subtracted = pr_subtracted.df.rename(columns={'Chromosome': 'read', 'Start': 'start', 'End': 'end', 'Feature': 'feature'})
        df_subtracted = df_subtracted[['read', 'start', 'end', 'feature']]
    else:
        df_subtracted = pd.DataFrame(columns=['read', 'start', 'end', 'feature'])

    # Combine priority intervals with subtracted satellite intervals
    df_priority_out = df_priority[['read', 'start', 'end', 'feature']]
    result = pd.concat([df_priority_out, df_subtracted], ignore_index=True)

    return result


def _telomere_satellite_merge_pandas(df_subtelo, df_satellite):
    """Fallback telomere-satellite merge using pandas (slower)."""
    common_reads = set(df_subtelo['read'].unique()) & set(df_satellite['read'].unique())
    if not common_reads:
        print("  Warning: No common reads between BED files")
        return pd.DataFrame(columns=['read', 'start', 'end', 'feature'])

    print(f"  Common reads: {len(common_reads):,}")
    print(f"  Priority features: {', '.join(sorted(TELOMERE_PRIORITY_FEATURES))}")
    print("  Note: Install pyranges for faster processing")

    # Filter to common reads and convert to int once
    df_subtelo_f = df_subtelo[df_subtelo['read'].isin(common_reads)].copy()
    df_satellite_f = df_satellite[df_satellite['read'].isin(common_reads)].copy()
    df_subtelo_f['start'] = df_subtelo_f['start'].astype(int)
    df_subtelo_f['end'] = df_subtelo_f['end'].astype(int)
    df_satellite_f['start'] = df_satellite_f['start'].astype(int)
    df_satellite_f['end'] = df_satellite_f['end'].astype(int)

    # Group by read for efficient iteration
    subtelo_grouped = {read: grp[['start', 'end', 'feature']].values
                       for read, grp in df_subtelo_f.groupby('read')}
    satellite_grouped = {read: grp[['start', 'end', 'feature']].values
                         for read, grp in df_satellite_f.groupby('read')}

    results = []
    for read in common_reads:
        subtelo_intervals = subtelo_grouped.get(read, [])
        satellite_intervals = satellite_grouped.get(read, [])

        # Extract priority intervals from subtelomeric
        priority_intervals = [(s, e, f) for s, e, f in subtelo_intervals
                              if f in TELOMERE_PRIORITY_FEATURES]

        # Sort priority intervals by start position for efficient subtraction
        priority_coords = sorted([(s, e) for s, e, _ in priority_intervals])

        # Add priority intervals to results
        for start, end, feature in priority_intervals:
            results.append((read, start, end, feature))

        # Add satellite intervals for non-priority positions using interval subtraction
        for start, end, feature in satellite_intervals:
            # Subtract priority intervals from this satellite interval
            uncovered = subtract_intervals((start, end), priority_coords)
            for seg_start, seg_end in uncovered:
                if seg_end > seg_start:
                    results.append((read, seg_start, seg_end, feature))

    return pd.DataFrame(results, columns=['read', 'start', 'end', 'feature'])


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
    if args.output.endswith('.gz'):
        merged_df.to_csv(args.output, sep='\t', index=False, header=False, compression='gzip')
    else:
        merged_df.to_csv(args.output, sep='\t', index=False, header=False)

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
