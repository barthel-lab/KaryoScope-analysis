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

args = parser.parse_args()

if len(args.bed) < 2:
    print("Error: At least 2 BED files required for merging", file=sys.stderr)
    sys.exit(1)


def load_bed_file(filepath):
    """Load a BED file into a DataFrame."""
    open_func = gzip.open if filepath.endswith('.gz') else open
    mode = 'rt' if filepath.endswith('.gz') else 'r'

    records = []
    with open_func(filepath, mode) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                records.append({
                    'read': parts[0],
                    'start': int(parts[1]),
                    'end': int(parts[2]),
                    'feature': parts[3]
                })
    return pd.DataFrame(records)


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
print(f"Merging {len(args.bed)} BED files...")
for i, bed_file in enumerate(args.bed, 1):
    print(f"  BED{i}: {bed_file}")

# Load first BED file
merged_df = load_bed_file(args.bed[0])
print(f"\n  BED1 intervals: {len(merged_df):,}")

# Iteratively merge with remaining BED files
for i, bed_file in enumerate(args.bed[1:], 2):
    df_next = load_bed_file(bed_file)
    print(f"  BED{i} intervals: {len(df_next):,}")

    merged_df = merge_two_beds(merged_df, df_next, args.separator)
    print(f"  After merge {i-1}: {len(merged_df):,} intervals")

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
