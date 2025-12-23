#!/usr/bin/env python3
"""
KaryoScope Sequence Annotation

Joins read_assignments.tsv with readnames.txt and stats.tsv files to add
sequencing approach and mapping statistics to each sequence.

Usage:
  python KaryoScope_annotate_sequences.py \
    --read-assignments analysis.read_assignments.tsv \
    --readnames-dir /path/to/samples \
    --output analysis.read_assignments.annotated.tsv

The script expects:
  - readnames.txt files at: {readnames_dir}/{sample}/telogator/{sample}.readnames.txt
  - stats.tsv files at: {readnames_dir}/{sample}/telogator/aligned/{sample}.CHM13.stats.tsv

Added columns:
  - sequencing_approach: From readnames.txt (e.g., "hifi", "ont")
  - n_alignments: Total number of alignments for the read
  - n_secondary: Number of secondary alignments (is_primary=False)
  - n_supplementary: Number of supplementary alignments (split reads)
  - primary_mapq: Mapping quality of primary non-supplementary alignment
  - primary_de: Divergence of primary alignment
  - primary_align_len: Aligned bases in primary alignment
  - primary_align_fraction: Fraction of read aligned in primary alignment
  - total_align_len: Sum of aligned bases across ALL alignments
  - total_align_fraction: Total aligned / read length (can exceed 1.0 if overlapping)

Join behavior:
  - Inner joins are performed: every read in read_assignments.tsv must have a
    matching entry in both readnames.txt and stats.tsv
  - The script will ERROR if any read from read_assignments.tsv is missing from
    readnames.txt or stats.tsv (this would indicate data corruption or mismatch)
  - Reads present in readnames.txt/stats.tsv but NOT in read_assignments.tsv are
    expected and ignored (these are reads filtered out by cluster analysis, e.g.,
    due to minimum read length requirements)
  - The output row count will always equal the input read_assignments.tsv row count
"""

import argparse
import os
import sys
import pandas as pd


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
            non_secondary = df[df['is_primary'] == True]
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


def main():
    parser = argparse.ArgumentParser(
        description="Annotate read assignments with sequencing approach and mapping stats",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--read-assignments", dest="read_assignments", required=True,
                        help="Input read_assignments.tsv file from cluster analysis")
    parser.add_argument("--readnames-dir", dest="readnames_dir", required=True,
                        help="Base directory containing sample folders with readnames.txt files")
    parser.add_argument("--reference", default="CHM13",
                        help="Reference genome name for stats files (default: CHM13)")
    parser.add_argument("--output", "-o", required=True,
                        help="Output annotated TSV file")
    parser.add_argument("--stats-columns", dest="stats_columns",
                        default="n_alignments,n_secondary,n_supplementary,primary_mapq,primary_de,primary_align_len,primary_align_fraction,total_align_len,total_align_fraction",
                        help="Comma-separated list of stats columns to include")

    args = parser.parse_args()

    # Parse stats columns
    stats_cols = [c.strip() for c in args.stats_columns.split(',')]

    print("=" * 60)
    print("KaryoScope Read Annotation")
    print("=" * 60)

    # Load read assignments
    print(f"\nLoading read assignments: {args.read_assignments}")
    assignments = pd.read_csv(args.read_assignments, sep='\t')
    print(f"  Total reads: {len(assignments)}")

    # Get unique samples
    samples = sorted(assignments['sample'].unique())
    print(f"  Samples: {', '.join(samples)}")

    # Load readnames
    print(f"\nLoading readnames files from: {args.readnames_dir}")
    readnames = load_readnames(args.readnames_dir, samples)
    print(f"  Total reads in readnames: {len(readnames)}")

    # Load stats
    print(f"\nLoading stats files from: {args.readnames_dir}")
    stats = load_stats(args.readnames_dir, samples, args.reference)
    print(f"  Total reads in stats: {len(stats)}")

    # === Join 1: assignments + readnames ===
    print("\n--- Joining with readnames ---")
    n_before = len(assignments)
    reads_in_assignments = set(assignments['read'])
    reads_in_readnames = set(readnames['read'])

    # Check that all assignment reads are in readnames (required)
    only_in_assignments = reads_in_assignments - reads_in_readnames
    # Reads only in readnames are expected (filtered out by cluster analysis)
    only_in_readnames = reads_in_readnames - reads_in_assignments

    if only_in_assignments:
        print(f"  ERROR: {len(only_in_assignments)} reads in assignments but not in readnames")
        print(f"  Examples: {list(only_in_assignments)[:5]}")
        sys.exit(1)

    if only_in_readnames:
        print(f"  Note: {len(only_in_readnames)} reads in readnames not in assignments (filtered by cluster analysis)")

    # Inner join
    merged = assignments.merge(readnames, on='read', how='inner')
    n_after = len(merged)

    if n_after != n_before:
        print(f"  ERROR: Row count changed from {n_before} to {n_after} after readnames join")
        sys.exit(1)

    print(f"  Joined successfully: {n_after} reads (no data loss from assignments)")

    # === Join 2: merged + stats ===
    print("\n--- Joining with stats ---")
    n_before = len(merged)
    reads_in_merged = set(merged['read'])
    reads_in_stats = set(stats['read'])

    # Check that all assignment reads are in stats (required)
    only_in_merged = reads_in_merged - reads_in_stats
    # Reads only in stats are expected (filtered out by cluster analysis)
    only_in_stats = reads_in_stats - reads_in_merged

    if only_in_merged:
        print(f"  ERROR: {len(only_in_merged)} reads in assignments but not in stats")
        print(f"  Examples: {list(only_in_merged)[:5]}")
        sys.exit(1)

    if only_in_stats:
        print(f"  Note: {len(only_in_stats)} reads in stats not in assignments (filtered by cluster analysis)")

    # Select only the columns we want from stats
    stats_subset = stats[['read'] + [c for c in stats_cols if c in stats.columns]].copy()
    missing_cols = [c for c in stats_cols if c not in stats.columns]
    if missing_cols:
        print(f"  Warning: Stats columns not found: {missing_cols}")

    # Inner join
    final = merged.merge(stats_subset, on='read', how='inner')
    n_after = len(final)

    if n_after != n_before:
        print(f"  ERROR: Row count changed from {n_before} to {n_after} after stats join")
        sys.exit(1)

    print(f"  Joined successfully: {n_after} reads (no data loss from assignments)")

    # === Output ===
    print(f"\n--- Writing output ---")
    print(f"  Output file: {args.output}")
    print(f"  Columns: {', '.join(final.columns)}")

    final.to_csv(args.output, sep='\t', index=False)

    print(f"\n{'=' * 60}")
    print("Summary")
    print("=" * 60)
    print(f"Total reads annotated: {len(final)}")
    print(f"Sequencing approaches: {final['sequencing_approach'].value_counts().to_dict()}")

    # Alignment statistics summary
    if 'n_alignments' in final.columns:
        print(f"\nAlignment statistics:")
        print(f"  Total alignments: {final['n_alignments'].sum():,}")
        print(f"  Reads with secondary alignments: {(final['n_secondary'] > 0).sum():,} ({(final['n_secondary'] > 0).mean()*100:.1f}%)")
        print(f"  Reads with supplementary alignments: {(final['n_supplementary'] > 0).sum():,} ({(final['n_supplementary'] > 0).mean()*100:.1f}%)")
        print(f"  Mean alignments per read: {final['n_alignments'].mean():.2f}")

    if 'primary_align_fraction' in final.columns and 'total_align_fraction' in final.columns:
        print(f"\nAlignment coverage:")
        print(f"  Primary alignment fraction: {final['primary_align_fraction'].mean()*100:.1f}% (mean)")
        print(f"  Total alignment fraction: {final['total_align_fraction'].mean()*100:.1f}% (mean)")
        multi_align = final[final['n_alignments'] > 1]
        if len(multi_align) > 0:
            print(f"  Multi-aligned reads coverage: {multi_align['total_align_fraction'].mean()*100:.1f}% (mean of {len(multi_align)} reads)")

    print(f"\nOutput saved to: {args.output}")


if __name__ == "__main__":
    main()
