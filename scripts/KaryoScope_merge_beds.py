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
import atexit
import datetime
import gzip
import pandas as pd
import sys

_original_command = ' '.join(sys.argv)

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
parser.add_argument("--priority-merge", dest="priority_merge", action="store_true",
                    help="3-way priority merge mode: subtel > region > repeat.\n"
                         "Requires exactly 3 BED files in order: subtelomeric, region, repeat.\n"
                         "Priority subtel features: canonical_telomere, noncanonical_telomere, ITS, TAR1, telomere_like_multigroup1\n"
                         "Region rules:\n"
                         "  - ct + nonrepeat -> ct; ct + other -> use repeat\n"
                         "  - noncentromeric + rRNA -> rRNA; noncentromeric + other -> rDNA\n"
                         "  - arm/p_arm/q_arm -> use repeat (background)\n"
                         "  - satellite features (censat, hsat*, etc.) -> keep region")
parser.add_argument("--chromosome-acrocentric-merge", dest="chromosome_acrocentric_merge", action="store_true",
                    help="Priority merge mode: acrocentric detail features have priority over chromosome labels.\n"
                         "Requires exactly 2 BED files: BED1=chromosome, BED2=acrocentric.\n"
                         "Acrocentric features (DJ, PHR, rDNA, SST1, PJ, array_multigroup1, acrocentric_multigroup1)\n"
                         "take priority; remaining positions filled with chromosome labels.")
parser.add_argument("--telomere-acrocentric-merge", dest="telomere_acrocentric_merge", action="store_true",
                    help="Priority merge mode: telomeric features from BED1 (subtelomeric) take priority,\n"
                         "acrocentric features from BED2 fill gaps. Composite labels (e.g., DJ_TAR1,\n"
                         "PHR_ITS) are created where TAR1/ITS overlap DJ/PHR/rDNA.\n"
                         "Requires exactly 2 BED files: BED1=subtelomeric, BED2=acrocentric.")
parser.add_argument("--collapse-non-acrocentric", dest="collapse_non_acrocentric",
                    action="store_true",
                    help="Remap specific non-acrocentric chromosome labels (chr1-12, chr16-20, chrX, chrY)\n"
                         "to 'non_acrocentric' in BED1 before merging.\n"
                         "Ambiguous labels (autosome_multigroup1, categorized, etc.) are preserved.\n"
                         "E.g., chr7 -> non_acrocentric, so merge produces non_acrocentric:rDNA.\n"
                         "Acrocentric chromosomes (chr13-15, chr21-22) are preserved.")
parser.add_argument("--log", default=None,
                    help="Path to log file. When provided, all console output is also written to this file.")

args = parser.parse_args()


# --- TeeLogger: duplicate stdout to log file ---
class TeeLogger:
    """Write to both stdout and a log file."""
    def __init__(self, log_path):
        self.terminal = sys.stdout
        self.log = open(log_path, 'w')

    def write(self, message):
        self.terminal.write(message)
        if not self.log.closed:
            self.log.write(message)
            self.log.flush()

    def flush(self):
        self.terminal.flush()
        if not self.log.closed:
            self.log.flush()

    def close(self):
        self.log.close()

if args.log:
    sys.stdout = TeeLogger(args.log)
    atexit.register(lambda: sys.stdout.close() if isinstance(sys.stdout, TeeLogger) else None)
    print(f"KaryoScope_merge_beds.py")
    print(f"  Timestamp: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Command: {_original_command}")
    print()

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

# Priority features for 3-way priority merge mode (includes telomere_like_multigroup1)
SUBTEL_PRIORITY_FEATURES = {'canonical_telomere', 'noncanonical_telomere', 'TAR1', 'ITS', 'telomere_like_multigroup1'}

# Region features that are considered "background" (can be overwritten by repeat)
REGION_BACKGROUND_FEATURES = {'p_arm', 'q_arm', 'arm_multigroup1'}

# Acrocentric features that have priority over chromosome labels
ACROCENTRIC_PRIORITY_FEATURES = {
    'DJ', 'PHR', 'rDNA', 'SST1', 'PJ',
    'array_multigroup1', 'acrocentric_multigroup1'
}

# Top-priority telomeric features for telomere-acrocentric merge mode
TELOMERE_ACRO_TOP_PRIORITY = {'canonical_telomere', 'noncanonical_telomere'}

# Subtelomeric features that can form composites with acrocentric features
TELOMERE_ACRO_COMPOSITE_SUBTEL = {'TAR1', 'ITS'}

# Acrocentric features that produce composite labels when overlapping TAR1/ITS
TELOMERE_ACRO_COMPOSITE_ACRO = {'DJ', 'PHR', 'rDNA'}

# Specific non-acrocentric chromosomes to collapse (strict mode)
STRICT_NON_ACROCENTRIC_CHROMOSOMES = {
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
    'chr10', 'chr11', 'chr12', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
    'chrX', 'chrY'
}


def remap_non_acrocentric_chromosomes(df):
    """
    Remap specific non-acrocentric chromosome labels to 'non_acrocentric'.
    Applied to chromosome BED before merging.
    """
    df = df.copy()
    original_features = set(df['feature'].unique())
    df['feature'] = df['feature'].apply(
        lambda f: 'non_acrocentric' if f in STRICT_NON_ACROCENTRIC_CHROMOSOMES else f
    )
    new_features = set(df['feature'].unique())
    collapsed = original_features - new_features
    if collapsed:
        print(f"\n  --- Non-acrocentric remap ---")
        print(f"  Chromosome features before: {len(original_features)}")
        print(f"  Chromosome features after:  {len(new_features)}")
        print(f"  Remapped {len(collapsed)} non-acrocentric chromosomes to 'non_acrocentric'")
    return df


def apply_conditional_region_repeat_rules(region_feature, repeat_feature):
    """
    Apply conditional rules for region + repeat feature combination.

    Rules:
    - ct + nonrepeat -> ct
    - ct + other repeat -> use repeat
    - noncentromeric + rRNA -> rRNA
    - noncentromeric + other -> rDNA
    - arm/p_arm/q_arm -> use repeat (background)
    - satellite features (censat, hsat*, etc.) -> keep region

    Args:
        region_feature: Feature from region BED
        repeat_feature: Feature from repeat BED

    Returns:
        Final feature label
    """
    # Background features: use repeat
    if region_feature in REGION_BACKGROUND_FEATURES:
        return repeat_feature

    # CT rules
    if region_feature == 'ct':
        if repeat_feature == 'nonrepeat':
            return 'ct'
        else:
            return repeat_feature

    # Noncentromeric rules
    if region_feature == 'noncentromeric':
        if repeat_feature == 'rRNA':
            return 'rRNA'
        else:
            return 'rDNA'

    # Satellite features (censat, hsat*, asat, bsat, gsat, etc.): keep region
    return region_feature


def priority_merge_three_way(df_subtel, df_region, df_repeat):
    """
    3-way priority merge: subtel > region (with conditional rules) > repeat.

    Args:
        df_subtel: DataFrame from subtelomeric BED file
        df_region: DataFrame from region BED file
        df_repeat: DataFrame from repeat BED file

    Returns:
        DataFrame with merged features (single feature per interval)
    """
    try:
        import pyranges as pr
        return _priority_merge_pyranges(df_subtel, df_region, df_repeat)
    except ImportError:
        return _priority_merge_pandas(df_subtel, df_region, df_repeat)


def _priority_merge_pyranges(df_subtel, df_region, df_repeat):
    """Fast 3-way priority merge using pyranges."""
    import pyranges as pr

    # Find common reads across all three
    common_reads = (set(df_subtel['read'].unique()) &
                   set(df_region['read'].unique()) &
                   set(df_repeat['read'].unique()))

    if not common_reads:
        print("  Warning: No common reads between all three BED files")
        return pd.DataFrame(columns=['read', 'start', 'end', 'feature'])

    print(f"  Common reads: {len(common_reads):,}")
    print(f"  Subtel priority features: {', '.join(sorted(SUBTEL_PRIORITY_FEATURES))}")

    # Filter to common reads
    df_subtel_f = df_subtel[df_subtel['read'].isin(common_reads)].copy()
    df_region_f = df_region[df_region['read'].isin(common_reads)].copy()
    df_repeat_f = df_repeat[df_repeat['read'].isin(common_reads)].copy()

    # Ensure integer types
    for df in [df_subtel_f, df_region_f, df_repeat_f]:
        df['start'] = df['start'].astype(int)
        df['end'] = df['end'].astype(int)

    # Step 1: Extract subtel priority intervals
    priority_mask = df_subtel_f['feature'].isin(SUBTEL_PRIORITY_FEATURES)
    df_subtel_priority = df_subtel_f[priority_mask][['read', 'start', 'end', 'feature']].copy()
    print(f"  Subtel priority intervals: {len(df_subtel_priority):,}")

    # Step 2: Join region and repeat to get conditional features
    # First, convert to pyranges
    pr_region = pr.PyRanges(df_region_f.rename(columns={
        'read': 'Chromosome', 'start': 'Start', 'end': 'End', 'feature': 'RegionFeature'
    }))
    pr_repeat = pr.PyRanges(df_repeat_f.rename(columns={
        'read': 'Chromosome', 'start': 'Start', 'end': 'End', 'feature': 'RepeatFeature'
    }))

    # Join region and repeat
    pr_region_repeat = pr_region.join(pr_repeat)

    if len(pr_region_repeat) == 0:
        print("  Warning: No overlap between region and repeat")
        # Fall back to just subtel priority features
        return df_subtel_priority.sort_values(['read', 'start'])

    # Convert to DataFrame and compute overlap coordinates
    df_rr = pr_region_repeat.df.copy()
    df_rr['overlap_start'] = df_rr[['Start', 'Start_b']].max(axis=1)
    df_rr['overlap_end'] = df_rr[['End', 'End_b']].min(axis=1)
    df_rr = df_rr[df_rr['overlap_end'] > df_rr['overlap_start']]

    # Apply conditional rules
    df_rr['feature'] = df_rr.apply(
        lambda row: apply_conditional_region_repeat_rules(row['RegionFeature'], row['RepeatFeature']),
        axis=1
    )

    # Prepare region+repeat merged DataFrame
    df_base = pd.DataFrame({
        'read': df_rr['Chromosome'],
        'start': df_rr['overlap_start'].astype(int),
        'end': df_rr['overlap_end'].astype(int),
        'feature': df_rr['feature']
    })
    print(f"  Region+repeat merged intervals: {len(df_base):,}")

    # Step 3: Subtract subtel priority regions from region+repeat base
    pr_subtel_priority = pr.PyRanges(df_subtel_priority.rename(columns={
        'read': 'Chromosome', 'start': 'Start', 'end': 'End', 'feature': 'Feature'
    }))
    pr_base = pr.PyRanges(df_base.rename(columns={
        'read': 'Chromosome', 'start': 'Start', 'end': 'End', 'feature': 'Feature'
    }))

    # Subtract subtel priority from base
    pr_subtracted = pr_base.subtract(pr_subtel_priority)

    # Convert results back
    if len(pr_subtracted) > 0:
        df_subtracted = pr_subtracted.df.rename(columns={
            'Chromosome': 'read', 'Start': 'start', 'End': 'end', 'Feature': 'feature'
        })[['read', 'start', 'end', 'feature']]
    else:
        df_subtracted = pd.DataFrame(columns=['read', 'start', 'end', 'feature'])

    # Step 4: Combine subtel priority + subtracted base
    result = pd.concat([df_subtel_priority, df_subtracted], ignore_index=True)

    return result.sort_values(['read', 'start'])


def _priority_merge_pandas(df_subtel, df_region, df_repeat):
    """Fallback 3-way priority merge using pandas (slower)."""
    common_reads = (set(df_subtel['read'].unique()) &
                   set(df_region['read'].unique()) &
                   set(df_repeat['read'].unique()))

    if not common_reads:
        print("  Warning: No common reads between all three BED files")
        return pd.DataFrame(columns=['read', 'start', 'end', 'feature'])

    print(f"  Common reads: {len(common_reads):,}")
    print(f"  Subtel priority features: {', '.join(sorted(SUBTEL_PRIORITY_FEATURES))}")
    print("  Note: Install pyranges for faster processing")

    # Filter to common reads
    df_subtel_f = df_subtel[df_subtel['read'].isin(common_reads)].copy()
    df_region_f = df_region[df_region['read'].isin(common_reads)].copy()
    df_repeat_f = df_repeat[df_repeat['read'].isin(common_reads)].copy()

    for df in [df_subtel_f, df_region_f, df_repeat_f]:
        df['start'] = df['start'].astype(int)
        df['end'] = df['end'].astype(int)

    # Group by read
    subtel_grouped = {read: grp[['start', 'end', 'feature']].values
                      for read, grp in df_subtel_f.groupby('read')}
    region_grouped = {read: grp[['start', 'end', 'feature']].values
                      for read, grp in df_region_f.groupby('read')}
    repeat_grouped = {read: grp[['start', 'end', 'feature']].values
                      for read, grp in df_repeat_f.groupby('read')}

    results = []
    for read in common_reads:
        subtel_intervals = subtel_grouped.get(read, [])
        region_intervals = region_grouped.get(read, [])
        repeat_intervals = repeat_grouped.get(read, [])

        # Extract subtel priority intervals
        priority_intervals = [(s, e, f) for s, e, f in subtel_intervals
                              if f in SUBTEL_PRIORITY_FEATURES]
        priority_coords = sorted([(s, e) for s, e, _ in priority_intervals])

        # Add subtel priority intervals
        for start, end, feature in priority_intervals:
            results.append((read, start, end, feature))

        # Compute region+repeat merged intervals
        for r_s, r_e, r_f in region_intervals:
            for rep_s, rep_e, rep_f in repeat_intervals:
                overlap_start = max(r_s, rep_s)
                overlap_end = min(r_e, rep_e)
                if overlap_end > overlap_start:
                    # Apply conditional rules
                    merged_feature = apply_conditional_region_repeat_rules(r_f, rep_f)

                    # Subtract subtel priority intervals
                    uncovered = subtract_intervals((overlap_start, overlap_end), priority_coords)
                    for seg_start, seg_end in uncovered:
                        if seg_end > seg_start:
                            results.append((read, seg_start, seg_end, merged_feature))

    return pd.DataFrame(results, columns=['read', 'start', 'end', 'feature']).sort_values(['read', 'start'])


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


def chromosome_acrocentric_merge(df_chrom, df_acro):
    """
    Priority merge: acrocentric detail features have priority over chromosome labels.
    Positions with meaningful acrocentric features (DJ, PHR, rDNA, etc.) use the
    acrocentric label; remaining positions use the chromosome label.

    Args:
        df_chrom: DataFrame from chromosome BED file
        df_acro: DataFrame from acrocentric BED file

    Returns:
        DataFrame with merged features (single feature per interval)
    """
    try:
        import pyranges as pr
        return _chromosome_acrocentric_merge_pyranges(df_chrom, df_acro)
    except ImportError:
        return _chromosome_acrocentric_merge_pandas(df_chrom, df_acro)


def _chromosome_acrocentric_merge_pyranges(df_chrom, df_acro):
    """Fast chromosome-acrocentric merge using pyranges subtract."""
    import pyranges as pr

    common_reads = set(df_chrom['read'].unique()) & set(df_acro['read'].unique())
    if not common_reads:
        print("  Warning: No common reads between BED files")
        return pd.DataFrame(columns=['read', 'start', 'end', 'feature'])

    print(f"  Common reads: {len(common_reads):,}")
    print(f"  Acrocentric priority features: {', '.join(sorted(ACROCENTRIC_PRIORITY_FEATURES))}")

    # Filter to common reads
    df_chrom_f = df_chrom[df_chrom['read'].isin(common_reads)].copy()
    df_acro_f = df_acro[df_acro['read'].isin(common_reads)].copy()

    # Extract acrocentric priority intervals
    priority_mask = df_acro_f['feature'].isin(ACROCENTRIC_PRIORITY_FEATURES)
    df_priority = df_acro_f[priority_mask][['read', 'start', 'end', 'feature']].copy()
    df_priority['start'] = df_priority['start'].astype(int)
    df_priority['end'] = df_priority['end'].astype(int)

    print(f"  Acrocentric priority intervals: {len(df_priority):,}")

    # Prepare chromosome intervals
    df_chrom_f = df_chrom_f[['read', 'start', 'end', 'feature']].copy()
    df_chrom_f['start'] = df_chrom_f['start'].astype(int)
    df_chrom_f['end'] = df_chrom_f['end'].astype(int)

    # Convert to pyranges format
    pr_priority = pr.PyRanges(df_priority.rename(columns={
        'read': 'Chromosome', 'start': 'Start', 'end': 'End', 'feature': 'Feature'
    }))
    pr_chrom = pr.PyRanges(df_chrom_f.rename(columns={
        'read': 'Chromosome', 'start': 'Start', 'end': 'End', 'feature': 'Feature'
    }))

    # Subtract acrocentric priority regions from chromosome
    pr_subtracted = pr_chrom.subtract(pr_priority)

    # Convert results back to DataFrame
    if len(pr_subtracted) > 0:
        df_subtracted = pr_subtracted.df.rename(columns={
            'Chromosome': 'read', 'Start': 'start', 'End': 'end', 'Feature': 'feature'
        })[['read', 'start', 'end', 'feature']]
    else:
        df_subtracted = pd.DataFrame(columns=['read', 'start', 'end', 'feature'])

    # Combine acrocentric priority + remaining chromosome
    result = pd.concat([df_priority, df_subtracted], ignore_index=True)

    return result.sort_values(['read', 'start'])


def _chromosome_acrocentric_merge_pandas(df_chrom, df_acro):
    """Fallback chromosome-acrocentric merge using pandas (slower)."""
    common_reads = set(df_chrom['read'].unique()) & set(df_acro['read'].unique())
    if not common_reads:
        print("  Warning: No common reads between BED files")
        return pd.DataFrame(columns=['read', 'start', 'end', 'feature'])

    print(f"  Common reads: {len(common_reads):,}")
    print(f"  Acrocentric priority features: {', '.join(sorted(ACROCENTRIC_PRIORITY_FEATURES))}")
    print("  Note: Install pyranges for faster processing")

    # Filter to common reads and convert to int
    df_chrom_f = df_chrom[df_chrom['read'].isin(common_reads)].copy()
    df_acro_f = df_acro[df_acro['read'].isin(common_reads)].copy()
    df_chrom_f['start'] = df_chrom_f['start'].astype(int)
    df_chrom_f['end'] = df_chrom_f['end'].astype(int)
    df_acro_f['start'] = df_acro_f['start'].astype(int)
    df_acro_f['end'] = df_acro_f['end'].astype(int)

    # Group by read
    chrom_grouped = {read: grp[['start', 'end', 'feature']].values
                     for read, grp in df_chrom_f.groupby('read')}
    acro_grouped = {read: grp[['start', 'end', 'feature']].values
                    for read, grp in df_acro_f.groupby('read')}

    results = []
    for read in common_reads:
        chrom_intervals = chrom_grouped.get(read, [])
        acro_intervals = acro_grouped.get(read, [])

        # Extract acrocentric priority intervals
        priority_intervals = [(s, e, f) for s, e, f in acro_intervals
                              if f in ACROCENTRIC_PRIORITY_FEATURES]
        priority_coords = sorted([(s, e) for s, e, _ in priority_intervals])

        # Add acrocentric priority intervals
        for start, end, feature in priority_intervals:
            results.append((read, start, end, feature))

        # Add chromosome intervals for non-priority positions
        for start, end, feature in chrom_intervals:
            uncovered = subtract_intervals((start, end), priority_coords)
            for seg_start, seg_end in uncovered:
                if seg_end > seg_start:
                    results.append((read, seg_start, seg_end, feature))

    return pd.DataFrame(results, columns=['read', 'start', 'end', 'feature']).sort_values(['read', 'start'])


def _resolve_telomere_acro_feature(subtel_feat, acro_feat):
    """
    Resolve feature label for a position with both subtelomeric and acrocentric annotations.

    Priority:
    1. canonical_telomere / noncanonical_telomere → always win
    2. TAR1/ITS + DJ/PHR/rDNA → composite label (e.g., DJ_TAR1)
    3. TAR1/ITS + other acro → keep subtel
    4. Other subtel + informative acro → keep acro
    5. Neither informative → nonacrocentric
    """
    subtel_is_top = subtel_feat in TELOMERE_ACRO_TOP_PRIORITY
    subtel_is_composite = subtel_feat in TELOMERE_ACRO_COMPOSITE_SUBTEL
    acro_is_informative = acro_feat in ACROCENTRIC_PRIORITY_FEATURES
    acro_is_composite = acro_feat in TELOMERE_ACRO_COMPOSITE_ACRO

    # Top-priority telomeric features always win
    if subtel_is_top:
        return subtel_feat

    # TAR1/ITS overlapping composite-eligible acrocentric → composite label
    if subtel_is_composite and acro_is_composite:
        return f"{acro_feat}_{subtel_feat}"

    # TAR1/ITS overlapping non-composite acrocentric → keep subtel
    if subtel_is_composite:
        return subtel_feat

    # telomere_like_multigroup1 or other subtel + informative acro → acro wins
    if acro_is_informative:
        return acro_feat

    # Both are background (nonsubtelomeric + nonacrocentric) → nonacrocentric
    return 'nonacrocentric'


def telomere_acrocentric_merge(df_subtel, df_acro):
    """
    Priority merge: telomeric features from subtelomeric BED take priority,
    acrocentric features fill gaps. Composite labels are created where
    TAR1/ITS overlap DJ/PHR/rDNA.

    Args:
        df_subtel: DataFrame from subtelomeric BED file
        df_acro: DataFrame from acrocentric BED file

    Returns:
        DataFrame with merged features (single feature per interval)
    """
    try:
        import pyranges as pr
        return _telomere_acrocentric_merge_pyranges(df_subtel, df_acro)
    except ImportError:
        return _telomere_acrocentric_merge_pandas(df_subtel, df_acro)


def _merge_adjacent_intervals(df):
    """Merge adjacent intervals with the same read and feature label."""
    if df.empty:
        return df

    df = df.sort_values(['read', 'start']).reset_index(drop=True)
    merged = []
    prev_read, prev_start, prev_end, prev_feat = df.iloc[0]

    for _, row in df.iloc[1:].iterrows():
        if row['read'] == prev_read and row['feature'] == prev_feat and row['start'] <= prev_end:
            prev_end = max(prev_end, row['end'])
        else:
            merged.append((prev_read, prev_start, prev_end, prev_feat))
            prev_read, prev_start, prev_end, prev_feat = row['read'], row['start'], row['end'], row['feature']

    merged.append((prev_read, prev_start, prev_end, prev_feat))
    return pd.DataFrame(merged, columns=['read', 'start', 'end', 'feature'])


def _telomere_acrocentric_merge_pyranges(df_subtel, df_acro):
    """Fast telomere-acrocentric merge using pyranges join."""
    import pyranges as pr

    common_reads = set(df_subtel['read'].unique()) & set(df_acro['read'].unique())
    if not common_reads:
        print("  Warning: No common reads between BED files")
        return pd.DataFrame(columns=['read', 'start', 'end', 'feature'])

    print(f"  Common reads: {len(common_reads):,}")
    print(f"  Top-priority telomere features: {', '.join(sorted(TELOMERE_ACRO_TOP_PRIORITY))}")
    print(f"  Composite-eligible pairs: {', '.join(f'{a}+{s}' for a in sorted(TELOMERE_ACRO_COMPOSITE_ACRO) for s in sorted(TELOMERE_ACRO_COMPOSITE_SUBTEL))}")

    # Filter to common reads and ensure int types
    df_subtel_f = df_subtel[df_subtel['read'].isin(common_reads)][['read', 'start', 'end', 'feature']].copy()
    df_acro_f = df_acro[df_acro['read'].isin(common_reads)][['read', 'start', 'end', 'feature']].copy()
    df_subtel_f['start'] = df_subtel_f['start'].astype(int)
    df_subtel_f['end'] = df_subtel_f['end'].astype(int)
    df_acro_f['start'] = df_acro_f['start'].astype(int)
    df_acro_f['end'] = df_acro_f['end'].astype(int)

    # Convert to pyranges
    pr_subtel = pr.PyRanges(df_subtel_f.rename(columns={
        'read': 'Chromosome', 'start': 'Start', 'end': 'End', 'feature': 'SubtelFeature'
    }))
    pr_acro = pr.PyRanges(df_acro_f.rename(columns={
        'read': 'Chromosome', 'start': 'Start', 'end': 'End', 'feature': 'AcroFeature'
    }))

    # Join the two interval sets
    pr_joined = pr_subtel.join(pr_acro)

    if len(pr_joined) == 0:
        print("  Warning: No overlap between subtelomeric and acrocentric BED files")
        return pd.DataFrame(columns=['read', 'start', 'end', 'feature'])

    # Convert to DataFrame and compute overlap coordinates
    df_j = pr_joined.df.copy()
    df_j['overlap_start'] = df_j[['Start', 'Start_b']].max(axis=1)
    df_j['overlap_end'] = df_j[['End', 'End_b']].min(axis=1)
    df_j = df_j[df_j['overlap_end'] > df_j['overlap_start']]

    # Resolve feature labels
    df_j['feature'] = df_j.apply(
        lambda row: _resolve_telomere_acro_feature(row['SubtelFeature'], row['AcroFeature']),
        axis=1
    )

    result = pd.DataFrame({
        'read': df_j['Chromosome'],
        'start': df_j['overlap_start'].astype(int),
        'end': df_j['overlap_end'].astype(int),
        'feature': df_j['feature']
    })

    # Merge adjacent intervals with same label
    result = _merge_adjacent_intervals(result)

    # Report feature counts
    feature_counts = result['feature'].value_counts()
    print(f"\n  Merged intervals: {len(result):,}")
    print(f"  Unique features: {len(feature_counts):,}")
    composite_feats = [f for f in feature_counts.index
                       if '_TAR1' in f or '_ITS' in f]
    if composite_feats:
        print(f"  Composite features:")
        for feat in composite_feats:
            bp = (result[result['feature'] == feat]['end'] - result[result['feature'] == feat]['start']).sum()
            print(f"    {feat}: {feature_counts[feat]:,} intervals, {bp:,} bp")

    return result


def _telomere_acrocentric_merge_pandas(df_subtel, df_acro):
    """Fallback telomere-acrocentric merge using pandas (slower)."""
    common_reads = set(df_subtel['read'].unique()) & set(df_acro['read'].unique())
    if not common_reads:
        print("  Warning: No common reads between BED files")
        return pd.DataFrame(columns=['read', 'start', 'end', 'feature'])

    print(f"  Common reads: {len(common_reads):,}")
    print(f"  Top-priority telomere features: {', '.join(sorted(TELOMERE_ACRO_TOP_PRIORITY))}")
    print("  Note: Install pyranges for faster processing")

    # Filter to common reads and convert to int
    df_subtel_f = df_subtel[df_subtel['read'].isin(common_reads)].copy()
    df_acro_f = df_acro[df_acro['read'].isin(common_reads)].copy()
    df_subtel_f['start'] = df_subtel_f['start'].astype(int)
    df_subtel_f['end'] = df_subtel_f['end'].astype(int)
    df_acro_f['start'] = df_acro_f['start'].astype(int)
    df_acro_f['end'] = df_acro_f['end'].astype(int)

    # Group by read
    subtel_grouped = {read: grp[['start', 'end', 'feature']].values
                      for read, grp in df_subtel_f.groupby('read')}
    acro_grouped = {read: grp[['start', 'end', 'feature']].values
                    for read, grp in df_acro_f.groupby('read')}

    results = []
    for read in common_reads:
        subtel_intervals = subtel_grouped.get(read, [])
        acro_intervals = acro_grouped.get(read, [])

        for s_s, s_e, s_f in subtel_intervals:
            for a_s, a_e, a_f in acro_intervals:
                overlap_start = max(int(s_s), int(a_s))
                overlap_end = min(int(s_e), int(a_e))
                if overlap_end > overlap_start:
                    feat = _resolve_telomere_acro_feature(str(s_f), str(a_f))
                    results.append((read, overlap_start, overlap_end, feat))

    result = pd.DataFrame(results, columns=['read', 'start', 'end', 'feature'])
    result = _merge_adjacent_intervals(result)

    # Report feature counts
    feature_counts = result['feature'].value_counts()
    print(f"\n  Merged intervals: {len(result):,}")
    print(f"  Unique features: {len(feature_counts):,}")
    composite_feats = [f for f in feature_counts.index
                       if '_TAR1' in f or '_ITS' in f]
    if composite_feats:
        print(f"  Composite features:")
        for feat in composite_feats:
            bp = (result[result['feature'] == feat]['end'] - result[result['feature'] == feat]['start']).sum()
            print(f"    {feat}: {feature_counts[feat]:,} intervals, {bp:,} bp")

    return result


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

    merged_df = merged_df.sort_values(['read', 'start'])

    print(f"\nWriting to: {args.output}")
    if args.output.endswith('.gz'):
        merged_df.to_csv(args.output, sep='\t', index=False, header=False, compression='gzip')
    else:
        merged_df.to_csv(args.output, sep='\t', index=False, header=False)

    print("Done!")
    sys.exit(0)

# Handle 3-way priority merge mode
if args.priority_merge:
    if len(args.bed) != 3:
        print("Error: --priority-merge requires exactly 3 BED files", file=sys.stderr)
        print("  BED1: subtelomeric features", file=sys.stderr)
        print("  BED2: region features", file=sys.stderr)
        print("  BED3: repeat features", file=sys.stderr)
        sys.exit(1)

    print("\n--- 3-Way Priority Merge Mode (subtel > region > repeat) ---")
    print("\nConditional rules:")
    print("  - ct + nonrepeat -> ct; ct + other -> use repeat")
    print("  - noncentromeric + rRNA -> rRNA; noncentromeric + other -> rDNA")
    print("  - arm/p_arm/q_arm -> use repeat (background)")
    print("  - satellite features (censat, hsat*, etc.) -> keep region")

    df_subtel = load_bed_file(args.bed[0])
    print(f"\n  Subtelomeric intervals: {len(df_subtel):,}")

    df_region = load_bed_file(args.bed[1])
    print(f"  Region intervals: {len(df_region):,}")

    df_repeat = load_bed_file(args.bed[2])
    print(f"  Repeat intervals: {len(df_repeat):,}")

    merged_df = priority_merge_three_way(df_subtel, df_region, df_repeat)
    print(f"\n  Final merged intervals: {len(merged_df):,}")

    if merged_df.empty:
        print("Error: Merge resulted in no intervals", file=sys.stderr)
        sys.exit(1)

    # Count unique features
    feature_counts = merged_df['feature'].value_counts()
    print(f"  Unique features: {len(feature_counts):,}")
    print("\n  Top 10 features:")
    for feat, count in feature_counts.head(10).items():
        pct = count / len(merged_df) * 100
        print(f"    {feat}: {count:,} ({pct:.1f}%)")

    print(f"\nWriting to: {args.output}")
    if args.output.endswith('.gz'):
        merged_df.to_csv(args.output, sep='\t', index=False, header=False, compression='gzip')
    else:
        merged_df.to_csv(args.output, sep='\t', index=False, header=False)

    print("Done!")
    sys.exit(0)

# Handle chromosome-acrocentric merge mode
if args.chromosome_acrocentric_merge:
    if len(args.bed) != 2:
        print("Error: --chromosome-acrocentric-merge requires exactly 2 BED files", file=sys.stderr)
        print("  BED1: chromosome features", file=sys.stderr)
        print("  BED2: acrocentric features", file=sys.stderr)
        sys.exit(1)

    print("\n--- Chromosome-Acrocentric Priority Merge Mode ---")
    print("  Acrocentric detail features take priority over chromosome labels")

    df_chrom = load_bed_file(args.bed[0])
    if args.collapse_non_acrocentric:
        df_chrom = remap_non_acrocentric_chromosomes(df_chrom)
    print(f"  Chromosome intervals: {len(df_chrom):,}")

    df_acro = load_bed_file(args.bed[1])
    print(f"  Acrocentric intervals: {len(df_acro):,}")

    merged_df = chromosome_acrocentric_merge(df_chrom, df_acro)
    print(f"\n  Merged intervals: {len(merged_df):,}")

    if merged_df.empty:
        print("Error: Merge resulted in no intervals", file=sys.stderr)
        sys.exit(1)

    # Count unique features
    feature_counts = merged_df['feature'].value_counts()
    print(f"  Unique features: {len(feature_counts):,}")
    print("\n  Top 10 features:")
    for feat, count in feature_counts.head(10).items():
        pct = count / len(merged_df) * 100
        print(f"    {feat}: {count:,} ({pct:.1f}%)")

    merged_df = merged_df.sort_values(['read', 'start'])

    print(f"\nWriting to: {args.output}")
    if args.output.endswith('.gz'):
        merged_df.to_csv(args.output, sep='\t', index=False, header=False, compression='gzip')
    else:
        merged_df.to_csv(args.output, sep='\t', index=False, header=False)

    print("Done!")
    sys.exit(0)

# Handle telomere-acrocentric merge mode
if args.telomere_acrocentric_merge:
    if len(args.bed) != 2:
        print("Error: --telomere-acrocentric-merge requires exactly 2 BED files", file=sys.stderr)
        print("  BED1: subtelomeric features", file=sys.stderr)
        print("  BED2: acrocentric features", file=sys.stderr)
        sys.exit(1)

    print("\n--- Telomere-Acrocentric Priority Merge Mode ---")
    print("  Telomeric features take priority; composite labels at DJ/PHR/rDNA boundaries")

    df_subtel = load_bed_file(args.bed[0])
    print(f"  Subtelomeric intervals: {len(df_subtel):,}")

    df_acro = load_bed_file(args.bed[1])
    print(f"  Acrocentric intervals: {len(df_acro):,}")

    merged_df = telomere_acrocentric_merge(df_subtel, df_acro)

    if merged_df.empty:
        print("Error: Merge resulted in no intervals", file=sys.stderr)
        sys.exit(1)

    # Count unique features
    feature_counts = merged_df['feature'].value_counts()
    print(f"\n  Top 10 features:")
    for feat, count in feature_counts.head(10).items():
        bp = (merged_df[merged_df['feature'] == feat]['end'] - merged_df[merged_df['feature'] == feat]['start']).sum()
        pct = count / len(merged_df) * 100
        print(f"    {feat}: {count:,} intervals, {bp:,} bp ({pct:.1f}%)")

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
if args.collapse_non_acrocentric:
    merged_df = remap_non_acrocentric_chromosomes(merged_df)
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
