#!/usr/bin/env python3
"""
KaryoScope_ITS_analysis.py - Analyze ITS element loss in ALT tumors.

Investigates the relationship between ITS (Interstitial Telomeric Sequences) and
TAR1 elements across Normal vs Tumor samples, supporting the hypothesis that
Break-Induced Replication (BIR) causes ITS loss in ALT tumors.

Two modes of operation:

1. CLUSTER MODE (recommended): Analyze only reads used in clustering
   - Uses read_assignments.tsv to filter to clustered reads
   - Correlates ITS/TAR1 with cluster enrichment direction
   - Avoids read length bias between samples

2. SAMPLE MODE: Analyze all reads from samples
   - Original mode, may have read length bias
   - Useful for exploratory analysis

Usage (cluster mode):
  python scripts/KaryoScope_ITS_analysis.py \\
    --cluster-analysis-prefix tmp/IDH_astro \\
    --bed-prefix results \\
    --output-prefix tmp/ITS_analysis \\
    --dark-mode

Usage (sample mode):
  python scripts/KaryoScope_ITS_analysis.py \\
    --bed-prefix results \\
    --samples BJ,IMR90,2436A,... \\
    --sample-metadata samples.tsv \\
    --output-prefix tmp/ITS_analysis \\
    --dark-mode
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use('Agg')
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, linregress
from pathlib import Path

# Register Basic Sans font
FONT_DIR = Path.home() / "Documents" / "Barthel-Custom-Powerpoint-Theme" / "fonts"
if FONT_DIR.exists():
    for font_file in FONT_DIR.glob("BasicSans-*.otf"):
        fm.fontManager.addfont(str(font_file))
    plt.rcParams['font.family'] = 'Basic Sans'
    plt.rcParams['font.size'] = 10
    # Use DejaVu Sans for math text (Greek letters, etc.)
    plt.rcParams['mathtext.fontset'] = 'dejavusans'
plt.rcParams['svg.fonttype'] = 'none'

# Capture original command line for logging
_original_command = ' '.join(sys.argv)


# --- TeeLogger for console + file logging ---
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


# --- Data classes ---
@dataclass
class FeatureInterval:
    """Single feature annotation on a read."""
    start: int
    end: int
    feature: str

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class ReadFeatures:
    """All features for a single read."""
    read_id: str
    sample: str
    read_length: int
    subtelomeric_features: List[FeatureInterval]
    region_features: List[FeatureInterval]
    chromosome_features: List[FeatureInterval]

    def get_its_intervals(self) -> List[FeatureInterval]:
        return [f for f in self.subtelomeric_features if f.feature == 'ITS']

    def get_tar1_intervals(self) -> List[FeatureInterval]:
        return [f for f in self.subtelomeric_features if f.feature == 'TAR1']

    def count_solitary_tar1(self, max_distance: int = 500) -> Tuple[int, int]:
        """
        Count TAR1 elements that are NOT adjacent to ITS elements.

        Returns: (solitary_tar1_count, total_tar1_count)

        A TAR1 is considered "solitary" if no ITS element is within max_distance bp.
        """
        tar1_intervals = self.get_tar1_intervals()
        its_intervals = self.get_its_intervals()

        if not tar1_intervals:
            return (0, 0)

        solitary_count = 0
        for tar1 in tar1_intervals:
            has_adjacent_its = False
            for its in its_intervals:
                # Check if ITS is within max_distance of TAR1
                # Distance is measured from nearest edges
                if tar1.end <= its.start:
                    distance = its.start - tar1.end
                elif its.end <= tar1.start:
                    distance = tar1.start - its.end
                else:
                    distance = 0  # Overlapping

                if distance <= max_distance:
                    has_adjacent_its = True
                    break

            if not has_adjacent_its:
                solitary_count += 1

        return (solitary_count, len(tar1_intervals))

    def total_its_bp(self) -> int:
        return sum(f.length for f in self.get_its_intervals())

    def total_tar1_bp(self) -> int:
        return sum(f.length for f in self.get_tar1_intervals())

    def dominant_chromosome(self) -> Optional[str]:
        """Return the chromosome with most coverage on this read."""
        chr_bp = {}
        for f in self.chromosome_features:
            if f.feature.startswith('chr'):
                chr_bp[f.feature] = chr_bp.get(f.feature, 0) + f.length
        if chr_bp:
            return max(chr_bp.keys(), key=lambda k: chr_bp[k])
        return None

    def dominant_arm(self) -> Optional[str]:
        """Return the chromosome arm with most coverage (p or q)."""
        arm_bp = {'p_arm': 0, 'q_arm': 0}
        for f in self.region_features:
            if f.feature in arm_bp:
                arm_bp[f.feature] += f.length
        if arm_bp['p_arm'] > arm_bp['q_arm']:
            return 'p'
        elif arm_bp['q_arm'] > arm_bp['p_arm']:
            return 'q'
        return None


@dataclass
class ITSTar1Pair:
    """A pair of ITS and TAR1 features on the same read."""
    read_id: str
    sample: str
    its_start: int
    its_end: int
    tar1_start: int
    tar1_end: int
    distance: int
    its_before_tar1: bool
    chromosome: Optional[str]
    arm: Optional[str]


@dataclass
class SampleStats:
    """Per-sample ITS/TAR1 statistics."""
    sample: str
    group: str
    n_reads: int
    total_read_length: int
    its_count: int
    its_total_bp: int
    tar1_count: int
    tar1_total_bp: int
    reads_with_its: int
    reads_with_tar1: int
    reads_with_both: int

    @property
    def its_per_mb(self) -> float:
        if self.total_read_length == 0:
            return 0.0
        return self.its_count / (self.total_read_length / 1e6)

    @property
    def tar1_per_mb(self) -> float:
        if self.total_read_length == 0:
            return 0.0
        return self.tar1_count / (self.total_read_length / 1e6)

    @property
    def its_bp_per_kb(self) -> float:
        if self.total_read_length == 0:
            return 0.0
        return self.its_total_bp / (self.total_read_length / 1e3)

    @property
    def tar1_bp_per_kb(self) -> float:
        if self.total_read_length == 0:
            return 0.0
        return self.tar1_total_bp / (self.total_read_length / 1e3)


# --- Argument parsing ---
def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze ITS element loss in ALT tumors",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Cluster mode arguments (recommended)
    parser.add_argument("--cluster-analysis-prefix", dest="cluster_prefix",
                        help="Prefix for cluster analysis files (e.g., 'tmp/IDH_astro').\n"
                             "Will load {prefix}.read_assignments.tsv and {prefix}.cluster_analysis.tsv")

    # Sample mode arguments (alternative)
    parser.add_argument("--samples",
                        help="Comma-separated list of sample names, or path to file with one per line.\n"
                             "Required if --cluster-analysis-prefix not provided.")
    parser.add_argument("--sample-metadata", dest="sample_metadata",
                        help="TSV file with sample metadata (columns: sample, group).\n"
                             "Required if --cluster-analysis-prefix not provided.")

    # Common input arguments
    parser.add_argument("--bed-prefix", dest="bed_prefix", required=True,
                        help="Base directory for BED files (e.g., 'results')")
    parser.add_argument("--database", default="KS_human_CHM13",
                        help="Database name (default: KS_human_CHM13)")
    parser.add_argument("--smoothness", default="smoothed",
                        help="Smoothness level (default: smoothed)")

    # Output arguments
    parser.add_argument("--output-prefix", dest="output_prefix", required=True,
                        help="Prefix for all output files")

    # Analysis parameters
    parser.add_argument("--min-read-length", dest="min_read_length", type=int, default=2000,
                        help="Minimum read length in bp (default: 2000)")
    parser.add_argument("--max-distance", dest="max_distance", type=int, default=5000,
                        help="Maximum distance (bp) to consider ITS-TAR1 as adjacent (default: 5000)")
    parser.add_argument("--min-feature-length", dest="min_feature_length", type=int, default=10,
                        help="Minimum feature length in bp to count (default: 10)")

    # Visualization options
    parser.add_argument("--dark-mode", dest="dark_mode", action="store_true",
                        help="Use dark background for plots")
    parser.add_argument("--figsize", default="12,8",
                        help="Figure size as 'width,height' in inches (default: 12,8)")

    # Enrichment mapping
    parser.add_argument("--enrichment-mapping", dest="enrichment_mapping", default=None,
                        help="JSON file mapping enrichment labels to standard Normal/Tumor terminology.\n"
                             "Format: {\"primary-enriched\": \"Normal-enriched\", \"E6E7-enriched\": \"Tumor-enriched\", ...}\n"
                             "When not provided, auto-detects from data using built-in heuristics.")

    # Logging
    parser.add_argument("--log-file", dest="log_file",
                        action=argparse.BooleanOptionalAction, default=True,
                        help="Save console output to {output_prefix}.log (default: True)")

    return parser.parse_args()


# --- Data loading functions ---
def load_sample_metadata(metadata_file: str) -> Dict[str, str]:
    """Load sample -> group mapping from metadata TSV."""
    df = pd.read_csv(metadata_file, sep='\t')
    return dict(zip(df['sample'], df['group']))


def _load_bed_features(filepath: str, min_length: int = 10) -> Dict[str, List[FeatureInterval]]:
    """Load BED file into dict mapping read_id to list of FeatureInterval."""
    if not os.path.exists(filepath):
        return {}

    features = defaultdict(list)
    open_func = gzip.open if filepath.endswith('.gz') else open
    mode = 'rt' if filepath.endswith('.gz') else 'r'

    with open_func(filepath, mode) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                read_id = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                feature = parts[3]
                length = end - start
                if length >= min_length:
                    features[read_id].append(FeatureInterval(start, end, feature))

    # Sort features by start position for each read
    for read_id in features:
        features[read_id].sort(key=lambda x: x.start)

    return dict(features)


def load_sample_features(sample: str, bed_prefix: str, database: str,
                         smoothness: str, min_read_length: int,
                         min_feature_length: int) -> Dict[str, ReadFeatures]:
    """Load all features for a sample from BED files."""
    base_path = f"{bed_prefix}/{sample}/telogator/1/KaryoScope/{database}/{sample}.telogator.1.{database}"

    # Try different path patterns
    subtelomeric_paths = [
        f"{base_path}.subtelomeric.{smoothness}.features.bed.gz",
        f"{base_path}.subtelomeric.{smoothness}.features.bed",
    ]
    region_paths = [
        f"{base_path}.region.{smoothness}.features.bed.gz",
        f"{base_path}.region.{smoothness}.features.bed",
    ]
    chromosome_paths = [
        f"{base_path}.chromosome.{smoothness}.features.bed.gz",
        f"{base_path}.chromosome.{smoothness}.features.bed",
    ]

    def find_existing(paths):
        for p in paths:
            if os.path.exists(p):
                return p
        return None

    subtelomeric_path = find_existing(subtelomeric_paths)
    region_path = find_existing(region_paths)
    chromosome_path = find_existing(chromosome_paths)

    subtelomeric = _load_bed_features(subtelomeric_path, min_feature_length) if subtelomeric_path else {}
    region = _load_bed_features(region_path, min_feature_length) if region_path else {}
    chromosome = _load_bed_features(chromosome_path, min_feature_length) if chromosome_path else {}

    # Combine into ReadFeatures objects
    all_reads = set(subtelomeric.keys()) | set(region.keys()) | set(chromosome.keys())

    results = {}
    for read_id in all_reads:
        # Calculate read length from max end position across all features
        all_features = (subtelomeric.get(read_id, []) +
                        region.get(read_id, []) +
                        chromosome.get(read_id, []))
        if not all_features:
            continue
        read_length = max(f.end for f in all_features)

        if read_length < min_read_length:
            continue

        results[read_id] = ReadFeatures(
            read_id=read_id,
            sample=sample,
            read_length=read_length,
            subtelomeric_features=subtelomeric.get(read_id, []),
            region_features=region.get(read_id, []),
            chromosome_features=chromosome.get(read_id, [])
        )

    return results


# --- Cluster mode data loading ---
def load_read_assignments(prefix: str) -> pd.DataFrame:
    """Load read assignments from cluster analysis."""
    path = f"{prefix}.read_assignments.tsv"
    df = pd.read_csv(path, sep='\t')
    print(f"  Loaded {len(df)} read assignments from {path}")
    return df


def load_cluster_analysis(prefix: str, enrichment_mapping: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    """Load cluster analysis results.

    Args:
        prefix: Path prefix for cluster analysis files.
        enrichment_mapping: Optional dict mapping enrichment labels to standard
            Normal/Tumor terminology. When None, uses built-in heuristics to
            auto-detect from data.
    """
    path = f"{prefix}.cluster_analysis.tsv"
    df = pd.read_csv(path, sep='\t')
    print(f"  Loaded {len(df)} clusters from {path}")

    if enrichment_mapping is None:
        # Auto-detect: map common group naming conventions to Normal/Tumor
        enrichment_mapping = {}
        if 'enrichment' in df.columns:
            unique_labels = df['enrichment'].unique()
            for label in unique_labels:
                lower = label.lower()
                if lower in ('normal-enriched', 'primary-enriched', 'control-enriched'):
                    enrichment_mapping[label] = 'Normal-enriched'
                elif lower in ('tumor-enriched', 'tumour-enriched', 'case-enriched'):
                    enrichment_mapping[label] = 'Tumor-enriched'
                elif lower.endswith('-enriched'):
                    # Heuristic: if only two groups, first alphabetically is Normal
                    enriched_labels = [l for l in unique_labels if l.lower().endswith('-enriched')]
                    if len(enriched_labels) == 2:
                        sorted_labels = sorted(enriched_labels)
                        enrichment_mapping[sorted_labels[0]] = 'Normal-enriched'
                        enrichment_mapping[sorted_labels[1]] = 'Tumor-enriched'
                    break

    if 'enrichment' in df.columns and enrichment_mapping:
        original_enrichments = df['enrichment'].value_counts().to_dict()
        df['enrichment'] = df['enrichment'].map(lambda x: enrichment_mapping.get(x, x))
        new_enrichments = df['enrichment'].value_counts().to_dict()
        if original_enrichments != new_enrichments:
            print(f"  Mapped enrichment labels: {original_enrichments} -> {new_enrichments}")

    return df


def load_features_for_clustered_reads(
    read_assignments: pd.DataFrame,
    bed_prefix: str,
    database: str,
    smoothness: str,
    min_feature_length: int
) -> Dict[str, Tuple[int, int, int, int, int]]:
    """Load ITS/TAR1 features for clustered reads only.

    Returns dict: read_id -> (its_count, its_bp, tar1_count, tar1_bp, read_length)
    """
    # Get unique samples and reads needed
    samples = read_assignments['sample'].unique()
    reads_needed = set(read_assignments['read'])
    read_lengths = dict(zip(read_assignments['read'], read_assignments['read_length']))

    print(f"  Loading features for {len(reads_needed)} clustered reads from {len(samples)} samples...")

    # Result: read_id -> (its_count, its_bp, tar1_count, tar1_bp, read_length)
    read_features = {}

    for sample in samples:
        base_path = f"{bed_prefix}/{sample}/telogator/1/KaryoScope/{database}/{sample}.telogator.1.{database}"

        # Find subtelomeric BED file
        subtelomeric_paths = [
            f"{base_path}.subtelomeric.{smoothness}.features.bed.gz",
            f"{base_path}.subtelomeric.{smoothness}.features.bed",
        ]

        subtelomeric_path = None
        for p in subtelomeric_paths:
            if os.path.exists(p):
                subtelomeric_path = p
                break

        if not subtelomeric_path:
            continue

        # Load features for reads from this sample
        open_func = gzip.open if subtelomeric_path.endswith('.gz') else open
        mode = 'rt' if subtelomeric_path.endswith('.gz') else 'r'

        # Temporary storage for this sample
        sample_its = defaultdict(lambda: {'count': 0, 'bp': 0})
        sample_tar1 = defaultdict(lambda: {'count': 0, 'bp': 0})

        with open_func(subtelomeric_path, mode) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    read_id = parts[0]
                    if read_id not in reads_needed:
                        continue
                    start, end = int(parts[1]), int(parts[2])
                    feature = parts[3]
                    length = end - start

                    if length < min_feature_length:
                        continue

                    if feature == 'ITS':
                        sample_its[read_id]['count'] += 1
                        sample_its[read_id]['bp'] += length
                    elif feature == 'TAR1':
                        sample_tar1[read_id]['count'] += 1
                        sample_tar1[read_id]['bp'] += length

        # Store results for reads from this sample
        sample_reads = read_assignments[read_assignments['sample'] == sample]['read']
        for read_id in sample_reads:
            its_data = sample_its.get(read_id, {'count': 0, 'bp': 0})
            tar1_data = sample_tar1.get(read_id, {'count': 0, 'bp': 0})
            read_length = read_lengths.get(read_id, 0)
            read_features[read_id] = (
                its_data['count'],
                its_data['bp'],
                tar1_data['count'],
                tar1_data['bp'],
                read_length
            )

    print(f"  Loaded features for {len(read_features)} reads")
    return read_features


def load_adjacency_data_for_clustered_reads(
    read_assignments: pd.DataFrame,
    bed_prefix: str,
    database: str,
    smoothness: str,
    min_feature_length: int,
    max_distance: int = 500
) -> Dict[str, Tuple[int, int]]:
    """Load ITS-TAR1 adjacency data for clustered reads.

    Returns dict: read_id -> (solitary_tar1_count, total_tar1_count)

    A TAR1 is "solitary" if no ITS is within max_distance bp.
    """
    samples = read_assignments['sample'].unique()
    reads_needed = set(read_assignments['read'])

    print(f"  Loading adjacency data for {len(reads_needed)} reads (max_distance={max_distance}bp)...")

    # Result: read_id -> (solitary_tar1, total_tar1)
    adjacency_data = {}

    for sample in samples:
        base_path = f"{bed_prefix}/{sample}/telogator/1/KaryoScope/{database}/{sample}.telogator.1.{database}"

        subtelomeric_paths = [
            f"{base_path}.subtelomeric.{smoothness}.features.bed.gz",
            f"{base_path}.subtelomeric.{smoothness}.features.bed",
        ]

        subtelomeric_path = None
        for p in subtelomeric_paths:
            if os.path.exists(p):
                subtelomeric_path = p
                break

        if not subtelomeric_path:
            continue

        open_func = gzip.open if subtelomeric_path.endswith('.gz') else open
        mode = 'rt' if subtelomeric_path.endswith('.gz') else 'r'

        # Store intervals per read: read_id -> {'its': [(start, end), ...], 'tar1': [(start, end), ...]}
        read_intervals = defaultdict(lambda: {'its': [], 'tar1': []})

        with open_func(subtelomeric_path, mode) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    read_id = parts[0]
                    if read_id not in reads_needed:
                        continue
                    start, end = int(parts[1]), int(parts[2])
                    feature = parts[3]
                    length = end - start

                    if length < min_feature_length:
                        continue

                    if feature == 'ITS':
                        read_intervals[read_id]['its'].append((start, end))
                    elif feature == 'TAR1':
                        read_intervals[read_id]['tar1'].append((start, end))

        # Compute solitary TAR1 for each read
        for read_id, intervals in read_intervals.items():
            tar1_list = intervals['tar1']
            its_list = intervals['its']

            if not tar1_list:
                adjacency_data[read_id] = (0, 0)
                continue

            solitary_count = 0
            for tar1_start, tar1_end in tar1_list:
                has_adjacent_its = False
                for its_start, its_end in its_list:
                    # Distance between nearest edges
                    if tar1_end <= its_start:
                        distance = its_start - tar1_end
                    elif its_end <= tar1_start:
                        distance = tar1_start - its_end
                    else:
                        distance = 0  # Overlapping

                    if distance <= max_distance:
                        has_adjacent_its = True
                        break

                if not has_adjacent_its:
                    solitary_count += 1

            adjacency_data[read_id] = (solitary_count, len(tar1_list))

    # Fill in reads with no TAR1 features found
    for read_id in reads_needed:
        if read_id not in adjacency_data:
            adjacency_data[read_id] = (0, 0)

    print(f"  Loaded adjacency data for {len(adjacency_data)} reads")
    return adjacency_data


def load_tar1_telomere_adjacency(
    read_assignments: pd.DataFrame,
    bed_prefix: str,
    database: str,
    smoothness: str,
    min_feature_length: int,
    max_distance: int = 2000
) -> Dict[str, Tuple[int, int, int, int]]:
    """Load TAR1-telomere adjacency data for clustered reads.

    Returns dict: read_id -> (telomeric_tar1, its_associated_tar1, total_tar1, has_both)

    A TAR1 is "telomeric" if canonical/noncanonical telomere is within max_distance bp.
    A TAR1 is "ITS-associated" if ITS is within max_distance bp.
    """
    samples = read_assignments['sample'].unique()
    reads_needed = set(read_assignments['read'])

    print(f"  Loading TAR1-telomere/ITS adjacency for {len(reads_needed)} reads...")

    # Result: read_id -> (telomeric_tar1, its_tar1, total_tar1, both_count)
    adjacency_data = {}

    for sample in samples:
        base_path = f"{bed_prefix}/{sample}/telogator/1/KaryoScope/{database}/{sample}.telogator.1.{database}"

        subtelomeric_paths = [
            f"{base_path}.subtelomeric.{smoothness}.features.bed.gz",
            f"{base_path}.subtelomeric.{smoothness}.features.bed",
        ]

        subtelomeric_path = None
        for p in subtelomeric_paths:
            if os.path.exists(p):
                subtelomeric_path = p
                break

        if not subtelomeric_path:
            continue

        open_func = gzip.open if subtelomeric_path.endswith('.gz') else open
        mode = 'rt' if subtelomeric_path.endswith('.gz') else 'r'

        # Store intervals: read_id -> {'tar1': [], 'its': [], 'telomere': []}
        read_intervals = defaultdict(lambda: {'tar1': [], 'its': [], 'telomere': []})

        with open_func(subtelomeric_path, mode) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    read_id = parts[0]
                    if read_id not in reads_needed:
                        continue
                    start, end = int(parts[1]), int(parts[2])
                    feature = parts[3]
                    length = end - start

                    if length < min_feature_length:
                        continue

                    if feature == 'TAR1':
                        read_intervals[read_id]['tar1'].append((start, end))
                    elif feature == 'ITS':
                        read_intervals[read_id]['its'].append((start, end))
                    elif feature in ('canonical_telomere', 'noncanonical_telomere'):
                        read_intervals[read_id]['telomere'].append((start, end))

        # Compute adjacency for each read
        for read_id, intervals in read_intervals.items():
            tar1_list = intervals['tar1']
            its_list = intervals['its']
            telomere_list = intervals['telomere']

            if not tar1_list:
                adjacency_data[read_id] = (0, 0, 0, 0)
                continue

            telomeric_count = 0
            its_count = 0
            both_count = 0

            for tar1_start, tar1_end in tar1_list:
                has_telomere = False
                has_its = False

                # Check telomere adjacency
                for telo_start, telo_end in telomere_list:
                    if tar1_end <= telo_start:
                        distance = telo_start - tar1_end
                    elif telo_end <= tar1_start:
                        distance = tar1_start - telo_end
                    else:
                        distance = 0
                    if distance <= max_distance:
                        has_telomere = True
                        break

                # Check ITS adjacency
                for its_start, its_end in its_list:
                    if tar1_end <= its_start:
                        distance = its_start - tar1_end
                    elif its_end <= tar1_start:
                        distance = tar1_start - its_end
                    else:
                        distance = 0
                    if distance <= max_distance:
                        has_its = True
                        break

                if has_telomere:
                    telomeric_count += 1
                if has_its:
                    its_count += 1
                if has_telomere and has_its:
                    both_count += 1

            adjacency_data[read_id] = (telomeric_count, its_count, len(tar1_list), both_count)

    # Fill in reads with no data
    for read_id in reads_needed:
        if read_id not in adjacency_data:
            adjacency_data[read_id] = (0, 0, 0, 0)

    print(f"  Loaded TAR1-telomere adjacency for {len(adjacency_data)} reads")
    return adjacency_data


def calculate_cluster_tar1_context(
    read_assignments: pd.DataFrame,
    cluster_analysis: pd.DataFrame,
    adjacency_data: Dict[str, Tuple[int, int, int, int]]
) -> pd.DataFrame:
    """Calculate TAR1 context (telomeric vs ITS-associated) per cluster."""
    results = []

    for _, cluster_row in cluster_analysis.iterrows():
        cluster_id = cluster_row['cluster_id']
        cluster_reads = read_assignments[read_assignments['cluster'] == cluster_id]['read']

        total_telomeric = 0
        total_its_associated = 0
        total_tar1 = 0
        total_both = 0

        for read_id in cluster_reads:
            telo, its, total, both = adjacency_data.get(read_id, (0, 0, 0, 0))
            total_telomeric += telo
            total_its_associated += its
            total_tar1 += total
            total_both += both

        pct_telomeric = (total_telomeric / total_tar1 * 100) if total_tar1 > 0 else 0
        pct_its_associated = (total_its_associated / total_tar1 * 100) if total_tar1 > 0 else 0
        pct_both = (total_both / total_tar1 * 100) if total_tar1 > 0 else 0
        # "Solitary" = neither telomere nor ITS nearby
        pct_solitary = 100 - pct_telomeric - pct_its_associated + pct_both  # Avoid double-counting

        results.append({
            'cluster_id': cluster_id,
            'size': cluster_row['size'],
            'enrichment': cluster_row['enrichment'],
            'odds_ratio': cluster_row['odds_ratio'],
            'log2_or': np.log2(cluster_row['odds_ratio']) if cluster_row['odds_ratio'] > 0 else 0,
            'total_tar1': total_tar1,
            'telomeric_tar1': total_telomeric,
            'its_associated_tar1': total_its_associated,
            'both_tar1': total_both,
            'pct_telomeric': pct_telomeric,
            'pct_its_associated': pct_its_associated,
            'pct_both': pct_both,
        })

    return pd.DataFrame(results)


def load_fragment_counts(
    read_assignments: pd.DataFrame,
    bed_prefix: str,
    database: str,
    smoothness: str,  # 'smoothed' or 'presmoothed'
    min_feature_length: int
) -> pd.DataFrame:
    """Count number of each feature type per read.

    Returns DataFrame with columns:
        read, sample, cluster, read_length,
        its_count, tar1_count, canonical_count, noncanonical_count,
        its_bp, tar1_bp, canonical_bp, noncanonical_bp
    """
    samples = read_assignments['sample'].unique()
    reads_needed = set(read_assignments['read'])
    read_info = read_assignments.set_index('read')[['sample', 'cluster', 'read_length']].to_dict('index')

    print(f"  Loading {smoothness} feature counts for {len(reads_needed)} reads...")

    # Initialize counts
    read_counts = defaultdict(lambda: {
        'its_count': 0, 'tar1_count': 0, 'canonical_count': 0, 'noncanonical_count': 0,
        'its_bp': 0, 'tar1_bp': 0, 'canonical_bp': 0, 'noncanonical_bp': 0
    })

    # Handle both smoothed (ITS) and presmoothed (ITS_specific) naming
    feature_map = {
        'ITS': 'its',
        'ITS_specific': 'its',
        'TAR1': 'tar1',
        'TAR1_specific': 'tar1',
        'canonical_telomere': 'canonical',
        'canonical_telomere_specific': 'canonical',
        'noncanonical_telomere': 'noncanonical',
        'noncanonical_telomere_specific': 'noncanonical'
    }

    for sample in samples:
        base_path = f"{bed_prefix}/{sample}/telogator/1/KaryoScope/{database}/{sample}.telogator.1.{database}"

        subtelomeric_paths = [
            f"{base_path}.subtelomeric.{smoothness}.features.bed.gz",
            f"{base_path}.subtelomeric.{smoothness}.features.bed",
        ]

        subtelomeric_path = None
        for p in subtelomeric_paths:
            if os.path.exists(p):
                subtelomeric_path = p
                break

        if not subtelomeric_path:
            continue

        open_func = gzip.open if subtelomeric_path.endswith('.gz') else open
        mode = 'rt' if subtelomeric_path.endswith('.gz') else 'r'

        with open_func(subtelomeric_path, mode) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    read_id = parts[0]
                    if read_id not in reads_needed:
                        continue
                    start, end = int(parts[1]), int(parts[2])
                    feature = parts[3]
                    length = end - start

                    if length < min_feature_length:
                        continue

                    if feature in feature_map:
                        key = feature_map[feature]
                        read_counts[read_id][f'{key}_count'] += 1
                        read_counts[read_id][f'{key}_bp'] += length

    # Build DataFrame
    rows = []
    for read_id in reads_needed:
        if read_id not in read_info:
            continue
        info = read_info[read_id]
        counts = read_counts.get(read_id, {
            'its_count': 0, 'tar1_count': 0, 'canonical_count': 0, 'noncanonical_count': 0,
            'its_bp': 0, 'tar1_bp': 0, 'canonical_bp': 0, 'noncanonical_bp': 0
        })

        rows.append({
            'read': read_id,
            'sample': info['sample'],
            'cluster': info['cluster'],
            'read_length': info['read_length'],
            **counts
        })

    df = pd.DataFrame(rows)
    print(f"  Loaded counts for {len(df)} reads")
    return df


def analyze_fragment_counts(
    smoothed_counts: pd.DataFrame,
    presmoothed_counts: pd.DataFrame,
    cluster_analysis: pd.DataFrame,
    output_path: str
) -> pd.DataFrame:
    """Compare fragment counts between smoothed and presmoothed, Normal vs Tumor."""
    results = []

    cluster_enrichment = cluster_analysis.set_index('cluster_id')['enrichment'].to_dict()

    for label, df in [('smoothed', smoothed_counts), ('presmoothed', presmoothed_counts)]:
        df = df.copy()
        df['enrichment'] = df['cluster'].map(cluster_enrichment)

        for feature in ['its', 'tar1', 'canonical', 'noncanonical']:
            count_col = f'{feature}_count'
            bp_col = f'{feature}_bp'

            for enrichment in ['Normal-enriched', 'Tumor-enriched']:
                subset = df[df['enrichment'] == enrichment]
                # Only include reads that have at least one of this feature
                with_feature = subset[subset[count_col] > 0]

                if len(with_feature) > 0:
                    results.append({
                        'smoothness': label,
                        'feature': feature.upper() if feature != 'canonical' and feature != 'noncanonical' else feature,
                        'enrichment': enrichment,
                        'n_reads_total': len(subset),
                        'n_reads_with_feature': len(with_feature),
                        'pct_reads_with_feature': 100 * len(with_feature) / len(subset),
                        'mean_count_per_read': subset[count_col].mean(),
                        'mean_count_if_present': with_feature[count_col].mean(),
                        'total_elements': subset[count_col].sum(),
                        'mean_element_size': subset[bp_col].sum() / subset[count_col].sum() if subset[count_col].sum() > 0 else 0
                    })

    results_df = pd.DataFrame(results)
    results_df.to_csv(output_path, sep='\t', index=False)
    return results_df


def plot_fragment_counts(
    smoothed_counts: pd.DataFrame,
    presmoothed_counts: pd.DataFrame,
    cluster_analysis: pd.DataFrame,
    output_path: str,
    dark_mode: bool = False
):
    """Plot fragment count comparison: smoothed vs presmoothed, Normal vs Tumor."""
    plt.rcParams['svg.fonttype'] = 'none'

    if dark_mode:
        plt.style.use('dark_background')
        text_color = 'white'
        bbox_color = '#333333'
    else:
        text_color = 'black'
        bbox_color = 'white'

    colors = {'Normal-enriched': '#3b82f6', 'Tumor-enriched': '#ef4444'}
    cluster_enrichment = cluster_analysis.set_index('cluster_id')['enrichment'].to_dict()

    fig, axes = plt.subplots(2, 4, figsize=(18, 10))

    features = ['its', 'tar1', 'canonical', 'noncanonical']
    feature_labels = ['ITS', 'TAR1', 'Canonical Telo', 'Non-canonical Telo']

    for col, (feature, label) in enumerate(zip(features, feature_labels)):
        count_col = f'{feature}_count'

        # Top row: smoothed
        ax = axes[0, col]
        smoothed = smoothed_counts.copy()
        smoothed['enrichment'] = smoothed['cluster'].map(cluster_enrichment)
        smoothed = smoothed[smoothed['enrichment'].isin(['Normal-enriched', 'Tumor-enriched'])]

        # Only reads with this feature
        smoothed_with = smoothed[smoothed[count_col] > 0]

        normal_counts = smoothed_with[smoothed_with['enrichment'] == 'Normal-enriched'][count_col]
        tumor_counts = smoothed_with[smoothed_with['enrichment'] == 'Tumor-enriched'][count_col]

        if len(normal_counts) > 0 and len(tumor_counts) > 0:
            bp = ax.boxplot([normal_counts, tumor_counts], labels=['Normal', 'Tumor'], patch_artist=True)
            bp['boxes'][0].set_facecolor(colors['Normal-enriched'])
            bp['boxes'][1].set_facecolor(colors['Tumor-enriched'])
            for box in bp['boxes']:
                box.set_alpha(0.7)

            _, pval = mannwhitneyu(normal_counts, tumor_counts)
            ax.text(0.5, 0.95, f'N: {normal_counts.mean():.2f}\nT: {tumor_counts.mean():.2f}\np={pval:.2e}',
                    transform=ax.transAxes, ha='center', va='top', fontsize=8,
                    bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.8))

        ax.set_ylabel(f'{label} fragments/read')
        ax.set_title(f'{label}\n(smoothed)')

        # Bottom row: presmoothed
        ax = axes[1, col]
        presmooth = presmoothed_counts.copy()
        presmooth['enrichment'] = presmooth['cluster'].map(cluster_enrichment)
        presmooth = presmooth[presmooth['enrichment'].isin(['Normal-enriched', 'Tumor-enriched'])]

        presmooth_with = presmooth[presmooth[count_col] > 0]

        normal_counts_pre = presmooth_with[presmooth_with['enrichment'] == 'Normal-enriched'][count_col]
        tumor_counts_pre = presmooth_with[presmooth_with['enrichment'] == 'Tumor-enriched'][count_col]

        if len(normal_counts_pre) > 0 and len(tumor_counts_pre) > 0:
            bp = ax.boxplot([normal_counts_pre, tumor_counts_pre], labels=['Normal', 'Tumor'], patch_artist=True)
            bp['boxes'][0].set_facecolor(colors['Normal-enriched'])
            bp['boxes'][1].set_facecolor(colors['Tumor-enriched'])
            for box in bp['boxes']:
                box.set_alpha(0.7)

            _, pval = mannwhitneyu(normal_counts_pre, tumor_counts_pre)
            ax.text(0.5, 0.95, f'N: {normal_counts_pre.mean():.2f}\nT: {tumor_counts_pre.mean():.2f}\np={pval:.2e}',
                    transform=ax.transAxes, ha='center', va='top', fontsize=8,
                    bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.8))

        ax.set_ylabel(f'{label} fragments/read')
        ax.set_title(f'{label}\n(presmoothed)')

    plt.tight_layout()
    fig.savefig(output_path, bbox_inches='tight')
    png_path = output_path.replace('.svg', '.png')
    fig.savefig(png_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved fragment count plot: {output_path} and {png_path}")


def analyze_normalized_fragmentation(
    smoothed_counts: pd.DataFrame,
    presmoothed_counts: pd.DataFrame,
    cluster_analysis: pd.DataFrame,
    output_path: str
) -> pd.DataFrame:
    """Analyze fragmentation normalized by feature size (fragments per kb of feature).

    This controls for telomere/feature length differences between samples.
    A higher value indicates more fragmentation per unit of sequence.

    Returns DataFrame with normalized fragmentation metrics.
    """
    results = []

    cluster_enrichment = cluster_analysis.set_index('cluster_id')['enrichment'].to_dict()

    for label, df in [('smoothed', smoothed_counts), ('presmoothed', presmoothed_counts)]:
        df = df.copy()
        df['enrichment'] = df['cluster'].map(cluster_enrichment)

        for feature in ['its', 'tar1', 'canonical', 'noncanonical']:
            count_col = f'{feature}_count'
            bp_col = f'{feature}_bp'

            # Calculate fragments per kb of feature for each read
            # Only for reads that have this feature (bp > 0)
            df[f'{feature}_frag_per_kb'] = np.where(
                df[bp_col] > 0,
                df[count_col] / (df[bp_col] / 1000),
                np.nan
            )

            # Also calculate fragments per kb of read (normalized by read length)
            df[f'{feature}_frag_per_kb_read'] = df[count_col] / (df['read_length'] / 1000)

            for enrichment in ['Normal-enriched', 'Tumor-enriched']:
                subset = df[df['enrichment'] == enrichment]

                # Reads with this feature
                with_feature = subset[subset[bp_col] > 0].copy()

                if len(with_feature) > 0:
                    frag_per_kb = with_feature[f'{feature}_frag_per_kb']
                    frag_per_kb_read = with_feature[f'{feature}_frag_per_kb_read']

                    results.append({
                        'smoothness': label,
                        'feature': feature.upper() if feature not in ['canonical', 'noncanonical'] else feature,
                        'enrichment': enrichment,
                        'n_reads': len(with_feature),
                        'total_feature_bp': with_feature[bp_col].sum(),
                        'total_fragments': with_feature[count_col].sum(),
                        'mean_frag_per_kb_feature': frag_per_kb.mean(),
                        'median_frag_per_kb_feature': frag_per_kb.median(),
                        'std_frag_per_kb_feature': frag_per_kb.std(),
                        'mean_frag_per_kb_read': frag_per_kb_read.mean(),
                        'median_frag_per_kb_read': frag_per_kb_read.median(),
                        'mean_bp_per_fragment': with_feature[bp_col].sum() / with_feature[count_col].sum()
                    })

    results_df = pd.DataFrame(results)
    results_df.to_csv(output_path, sep='\t', index=False)
    return results_df


def plot_normalized_fragmentation(
    smoothed_counts: pd.DataFrame,
    presmoothed_counts: pd.DataFrame,
    cluster_analysis: pd.DataFrame,
    output_path: str,
    dark_mode: bool = False
):
    """Plot normalized fragmentation analysis: fragments per kb of feature.

    Creates 2x4 grid:
    - Top row: smoothed data
    - Bottom row: presmoothed data
    - Columns: ITS, TAR1, Canonical, Non-canonical
    """
    plt.rcParams['svg.fonttype'] = 'none'

    if dark_mode:
        plt.style.use('dark_background')
        text_color = 'white'
        bbox_color = '#333333'
    else:
        text_color = 'black'
        bbox_color = 'white'

    colors = {'Normal-enriched': '#3b82f6', 'Tumor-enriched': '#ef4444'}
    cluster_enrichment = cluster_analysis.set_index('cluster_id')['enrichment'].to_dict()

    fig, axes = plt.subplots(2, 4, figsize=(18, 10))

    features = ['its', 'tar1', 'canonical', 'noncanonical']
    feature_labels = ['ITS', 'TAR1', 'Canonical Telo', 'Non-canonical Telo']

    for col, (feature, label) in enumerate(zip(features, feature_labels)):
        count_col = f'{feature}_count'
        bp_col = f'{feature}_bp'

        # Top row: smoothed
        ax = axes[0, col]
        smoothed = smoothed_counts.copy()
        smoothed['enrichment'] = smoothed['cluster'].map(cluster_enrichment)
        smoothed = smoothed[smoothed['enrichment'].isin(['Normal-enriched', 'Tumor-enriched'])]

        # Calculate fragments per kb of feature
        smoothed['frag_per_kb'] = np.where(
            smoothed[bp_col] > 0,
            smoothed[count_col] / (smoothed[bp_col] / 1000),
            np.nan
        )

        # Only reads with this feature
        smoothed_with = smoothed[smoothed[bp_col] > 0]

        normal_frag = smoothed_with[smoothed_with['enrichment'] == 'Normal-enriched']['frag_per_kb'].dropna()
        tumor_frag = smoothed_with[smoothed_with['enrichment'] == 'Tumor-enriched']['frag_per_kb'].dropna()

        if len(normal_frag) > 0 and len(tumor_frag) > 0:
            bp = ax.boxplot([normal_frag, tumor_frag], labels=['Normal', 'Tumor'], patch_artist=True)
            bp['boxes'][0].set_facecolor(colors['Normal-enriched'])
            bp['boxes'][1].set_facecolor(colors['Tumor-enriched'])
            for box in bp['boxes']:
                box.set_alpha(0.7)

            # Set reasonable y-limits (exclude extreme outliers for visualization)
            combined = pd.concat([normal_frag, tumor_frag])
            q99 = combined.quantile(0.99)
            ax.set_ylim(0, min(q99 * 1.5, combined.max()))

            _, pval = mannwhitneyu(normal_frag, tumor_frag)
            fold_change = tumor_frag.mean() / normal_frag.mean() if normal_frag.mean() > 0 else np.inf
            ax.text(0.5, 0.95, f'N: {normal_frag.mean():.2f}\nT: {tumor_frag.mean():.2f}\nFC: {fold_change:.2f}x\np={pval:.2e}',
                    transform=ax.transAxes, ha='center', va='top', fontsize=8,
                    bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.8))

        ax.set_ylabel(f'{label}\nfragments/kb of {label}')
        ax.set_title(f'{label}\n(smoothed)')

        # Bottom row: presmoothed
        ax = axes[1, col]
        presmooth = presmoothed_counts.copy()
        presmooth['enrichment'] = presmooth['cluster'].map(cluster_enrichment)
        presmooth = presmooth[presmooth['enrichment'].isin(['Normal-enriched', 'Tumor-enriched'])]

        # Calculate fragments per kb of feature
        presmooth['frag_per_kb'] = np.where(
            presmooth[bp_col] > 0,
            presmooth[count_col] / (presmooth[bp_col] / 1000),
            np.nan
        )

        presmooth_with = presmooth[presmooth[bp_col] > 0]

        normal_frag_pre = presmooth_with[presmooth_with['enrichment'] == 'Normal-enriched']['frag_per_kb'].dropna()
        tumor_frag_pre = presmooth_with[presmooth_with['enrichment'] == 'Tumor-enriched']['frag_per_kb'].dropna()

        if len(normal_frag_pre) > 0 and len(tumor_frag_pre) > 0:
            bp = ax.boxplot([normal_frag_pre, tumor_frag_pre], labels=['Normal', 'Tumor'], patch_artist=True)
            bp['boxes'][0].set_facecolor(colors['Normal-enriched'])
            bp['boxes'][1].set_facecolor(colors['Tumor-enriched'])
            for box in bp['boxes']:
                box.set_alpha(0.7)

            # Set reasonable y-limits
            combined = pd.concat([normal_frag_pre, tumor_frag_pre])
            q99 = combined.quantile(0.99)
            ax.set_ylim(0, min(q99 * 1.5, combined.max()))

            _, pval = mannwhitneyu(normal_frag_pre, tumor_frag_pre)
            fold_change = tumor_frag_pre.mean() / normal_frag_pre.mean() if normal_frag_pre.mean() > 0 else np.inf
            ax.text(0.5, 0.95, f'N: {normal_frag_pre.mean():.2f}\nT: {tumor_frag_pre.mean():.2f}\nFC: {fold_change:.2f}x\np={pval:.2e}',
                    transform=ax.transAxes, ha='center', va='top', fontsize=8,
                    bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.8))

        ax.set_ylabel(f'{label}\nfragments/kb of {label}')
        ax.set_title(f'{label}\n(presmoothed)')

    fig.suptitle('Normalized Fragmentation: Fragments per kb of Feature\n(Higher = more fragmented per unit sequence)',
                 fontsize=12, fontweight='bold', y=1.02)

    plt.tight_layout()
    fig.savefig(output_path, bbox_inches='tight')
    png_path = output_path.replace('.svg', '.png')
    fig.savefig(png_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved normalized fragmentation plot: {output_path}")


def load_feature_sequence_data(
    read_assignments: pd.DataFrame,
    bed_prefix: str,
    database: str,
    smoothness: str,
    min_feature_length: int
) -> pd.DataFrame:
    """Load all subtelomeric features with positions to analyze interleaving.

    Returns DataFrame with one row per feature element:
        read, sample, cluster, start, end, feature, size
    """
    samples = read_assignments['sample'].unique()
    reads_needed = set(read_assignments['read'])
    read_info = read_assignments.set_index('read')[['sample', 'cluster', 'read_length']].to_dict('index')

    print(f"  Loading feature sequences for {len(reads_needed)} reads...")

    elements = []
    target_features = {'ITS', 'canonical_telomere', 'noncanonical_telomere', 'TAR1'}

    for sample in samples:
        base_path = f"{bed_prefix}/{sample}/telogator/1/KaryoScope/{database}/{sample}.telogator.1.{database}"

        subtelomeric_paths = [
            f"{base_path}.subtelomeric.{smoothness}.features.bed.gz",
            f"{base_path}.subtelomeric.{smoothness}.features.bed",
        ]

        subtelomeric_path = None
        for p in subtelomeric_paths:
            if os.path.exists(p):
                subtelomeric_path = p
                break

        if not subtelomeric_path:
            continue

        open_func = gzip.open if subtelomeric_path.endswith('.gz') else open
        mode = 'rt' if subtelomeric_path.endswith('.gz') else 'r'

        with open_func(subtelomeric_path, mode) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    read_id = parts[0]
                    if read_id not in reads_needed:
                        continue
                    start, end = int(parts[1]), int(parts[2])
                    feature = parts[3]
                    length = end - start

                    if length < min_feature_length:
                        continue

                    if feature in target_features and read_id in read_info:
                        info = read_info[read_id]
                        elements.append({
                            'read': read_id,
                            'sample': info['sample'],
                            'cluster': info['cluster'],
                            'read_length': info['read_length'],
                            'start': start,
                            'end': end,
                            'feature': feature,
                            'size': length
                        })

    df = pd.DataFrame(elements)
    print(f"  Loaded {len(df)} feature elements from {df['read'].nunique() if len(df) > 0 else 0} reads")
    return df


def analyze_feature_interleaving(
    feature_elements: pd.DataFrame,
    cluster_analysis: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Analyze interleaving of ITS with telomeric features.

    Looks at transitions between adjacent features along reads.
    More ITS<->telomere transitions = more mixing/recombination.

    Returns:
        read_transitions: per-read transition counts
        cluster_transitions: per-cluster transition statistics
    """
    if feature_elements.empty:
        return pd.DataFrame(), pd.DataFrame()

    # Add enrichment info
    cluster_enrichment = cluster_analysis.set_index('cluster_id')['enrichment'].to_dict()
    df = feature_elements.copy()
    df['enrichment'] = df['cluster'].map(cluster_enrichment)

    # Simplify feature names for transition analysis
    feature_map = {
        'ITS': 'ITS',
        'canonical_telomere': 'Telo',
        'noncanonical_telomere': 'Telo',
        'TAR1': 'TAR1'
    }
    df['feature_class'] = df['feature'].map(feature_map)

    # Sort features by position within each read
    df = df.sort_values(['read', 'start'])

    # Calculate transitions per read
    read_transitions = []

    for read_id, read_df in df.groupby('read'):
        if len(read_df) < 2:
            continue

        features = read_df['feature_class'].values
        enrichment = read_df['enrichment'].iloc[0]
        cluster = read_df['cluster'].iloc[0]

        # Count transitions
        n_its_telo = 0  # ITS -> Telo or Telo -> ITS
        n_its_tar1 = 0  # ITS -> TAR1 or TAR1 -> ITS
        n_telo_tar1 = 0  # Telo -> TAR1 or TAR1 -> Telo
        n_same = 0  # Same feature type adjacent
        total_transitions = len(features) - 1

        for i in range(len(features) - 1):
            f1, f2 = features[i], features[i + 1]
            if f1 == f2:
                n_same += 1
            elif (f1 == 'ITS' and f2 == 'Telo') or (f1 == 'Telo' and f2 == 'ITS'):
                n_its_telo += 1
            elif (f1 == 'ITS' and f2 == 'TAR1') or (f1 == 'TAR1' and f2 == 'ITS'):
                n_its_tar1 += 1
            elif (f1 == 'Telo' and f2 == 'TAR1') or (f1 == 'TAR1' and f2 == 'Telo'):
                n_telo_tar1 += 1

        # Count distinct features
        has_its = 'ITS' in features
        has_telo = 'Telo' in features
        has_tar1 = 'TAR1' in features

        read_transitions.append({
            'read': read_id,
            'cluster': cluster,
            'enrichment': enrichment,
            'n_features': len(features),
            'n_transitions': total_transitions,
            'n_its_telo_transitions': n_its_telo,
            'n_its_tar1_transitions': n_its_tar1,
            'n_telo_tar1_transitions': n_telo_tar1,
            'n_same_transitions': n_same,
            'pct_its_telo': 100 * n_its_telo / total_transitions if total_transitions > 0 else 0,
            'has_its': has_its,
            'has_telo': has_telo,
            'has_tar1': has_tar1,
            'has_its_and_telo': has_its and has_telo,
            'interleaving_score': n_its_telo + n_its_tar1  # Higher = more mixing
        })

    read_trans_df = pd.DataFrame(read_transitions)

    # Aggregate to cluster level
    if read_trans_df.empty:
        return pd.DataFrame(), pd.DataFrame()

    cluster_trans = read_trans_df.groupby('cluster').agg({
        'n_its_telo_transitions': ['sum', 'mean'],
        'n_its_tar1_transitions': ['sum', 'mean'],
        'interleaving_score': ['sum', 'mean'],
        'n_transitions': 'sum',
        'has_its_and_telo': 'sum',
        'enrichment': 'first'
    })
    cluster_trans.columns = ['its_telo_total', 'its_telo_mean', 'its_tar1_total', 'its_tar1_mean',
                             'interleaving_total', 'interleaving_mean', 'total_transitions',
                             'reads_with_its_and_telo', 'enrichment']
    cluster_trans = cluster_trans.reset_index()

    # Calculate interleaving rate
    cluster_trans['its_telo_rate'] = cluster_trans['its_telo_total'] / cluster_trans['total_transitions']
    cluster_trans['interleaving_rate'] = cluster_trans['interleaving_total'] / cluster_trans['total_transitions']

    # Add cluster info
    cluster_trans = cluster_trans.merge(
        cluster_analysis[['cluster_id', 'size', 'odds_ratio']],
        left_on='cluster', right_on='cluster_id', how='left'
    )

    return read_trans_df, cluster_trans


def plot_feature_interleaving(
    read_transitions: pd.DataFrame,
    cluster_transitions: pd.DataFrame,
    output_path: str,
    dark_mode: bool = False
):
    """Plot feature interleaving analysis."""
    plt.rcParams['svg.fonttype'] = 'none'

    if dark_mode:
        plt.style.use('dark_background')
        text_color = 'white'
        bbox_color = '#333333'
    else:
        text_color = 'black'
        bbox_color = 'white'

    colors = {'Normal-enriched': '#3b82f6', 'Tumor-enriched': '#ef4444'}

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Filter to Normal and Tumor
    read_df = read_transitions[read_transitions['enrichment'].isin(['Normal-enriched', 'Tumor-enriched'])]
    cluster_df = cluster_transitions[cluster_transitions['enrichment'].isin(['Normal-enriched', 'Tumor-enriched'])]

    # Panel A: ITS-Telomere transitions per read (only reads with both)
    ax = axes[0, 0]
    reads_with_both = read_df[read_df['has_its_and_telo']]
    normal_trans = reads_with_both[reads_with_both['enrichment'] == 'Normal-enriched']['n_its_telo_transitions']
    tumor_trans = reads_with_both[reads_with_both['enrichment'] == 'Tumor-enriched']['n_its_telo_transitions']

    bp = ax.boxplot([normal_trans, tumor_trans], labels=['Normal', 'Tumor'], patch_artist=True)
    bp['boxes'][0].set_facecolor(colors['Normal-enriched'])
    bp['boxes'][1].set_facecolor(colors['Tumor-enriched'])
    for box in bp['boxes']:
        box.set_alpha(0.7)

    ax.set_ylabel('ITS<->Telomere transitions')
    ax.set_title('A. ITS-Telomere Mixing\n(reads with both features)')

    if len(normal_trans) > 1 and len(tumor_trans) > 1:
        _, pval = mannwhitneyu(normal_trans, tumor_trans)
        ax.text(0.5, 0.95, f'Normal: {normal_trans.mean():.2f}\nTumor: {tumor_trans.mean():.2f}\np={pval:.2e}',
                transform=ax.transAxes, ha='center', va='top', fontsize=9,
                bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.8))

    # Panel B: Distribution of interleaving scores
    ax = axes[0, 1]
    for enrichment in ['Normal-enriched', 'Tumor-enriched']:
        subset = reads_with_both[reads_with_both['enrichment'] == enrichment]['interleaving_score']
        label = 'Normal' if 'Normal' in enrichment else 'Tumor'
        ax.hist(subset, bins=range(0, 15), alpha=0.6, label=f'{label} (n={len(subset)})',
                color=colors[enrichment], density=True)

    ax.set_xlabel('Interleaving score (ITS<->Telo + ITS<->TAR1)')
    ax.set_ylabel('Density')
    ax.set_title('B. Interleaving Score Distribution')
    ax.legend(fontsize=8)

    # Panel C: % reads with ITS-Telomere transitions
    ax = axes[0, 2]
    # Among reads with both ITS and Telo, what % have at least one transition?
    for enrichment in ['Normal-enriched', 'Tumor-enriched']:
        subset = reads_with_both[reads_with_both['enrichment'] == enrichment]
        has_transition = (subset['n_its_telo_transitions'] > 0).sum()
        total = len(subset)
        pct = 100 * has_transition / total if total > 0 else 0
        label = 'Normal' if 'Normal' in enrichment else 'Tumor'
        ax.bar(label, pct, color=colors[enrichment], alpha=0.7)
        ax.text(label, pct + 2, f'{pct:.1f}%\n(n={total})', ha='center', fontsize=9)

    ax.set_ylabel('% reads with ITS<->Telo transition')
    ax.set_title('C. Reads with Mixed ITS-Telomere')
    ax.set_ylim(0, 100)

    # Panel D: Cluster-level interleaving rate vs enrichment
    ax = axes[1, 0]
    for enrichment in ['Normal-enriched', 'Tumor-enriched']:
        subset = cluster_df[cluster_df['enrichment'] == enrichment]
        label = 'Normal' if 'Normal' in enrichment else 'Tumor'
        ax.scatter(np.log2(subset['odds_ratio'].replace(0, np.nan)),
                  subset['interleaving_rate'] * 100,
                  s=subset['size'] / 3, c=colors[enrichment], alpha=0.6,
                  label=f'{label} (n={len(subset)})')

    ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('log2(Odds Ratio)')
    ax.set_ylabel('Interleaving rate (%)')
    ax.set_title('D. Interleaving vs Enrichment')
    ax.legend(fontsize=8)

    # Correlation
    valid = cluster_df.dropna(subset=['odds_ratio', 'interleaving_rate'])
    valid = valid[valid['odds_ratio'] > 0]
    if len(valid) > 5:
        rho, pval = spearmanr(np.log2(valid['odds_ratio']), valid['interleaving_rate'])
        ax.text(0.02, 0.98, f'r = {rho:.2f}, p = {pval:.2e}',
                transform=ax.transAxes, va='top', fontsize=9,
                bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.8))

    # Panel E: Cluster-level boxplot
    ax = axes[1, 1]
    normal_rate = cluster_df[cluster_df['enrichment'] == 'Normal-enriched']['interleaving_rate'] * 100
    tumor_rate = cluster_df[cluster_df['enrichment'] == 'Tumor-enriched']['interleaving_rate'] * 100

    bp = ax.boxplot([normal_rate, tumor_rate], labels=['Normal', 'Tumor'], patch_artist=True)
    bp['boxes'][0].set_facecolor(colors['Normal-enriched'])
    bp['boxes'][1].set_facecolor(colors['Tumor-enriched'])
    for box in bp['boxes']:
        box.set_alpha(0.7)

    ax.set_ylabel('Interleaving rate (%)')
    ax.set_title('E. Cluster-Level Interleaving')

    if len(normal_rate) > 1 and len(tumor_rate) > 1:
        _, pval = mannwhitneyu(normal_rate, tumor_rate)
        ax.text(0.5, 0.95, f'Normal: {normal_rate.mean():.2f}%\nTumor: {tumor_rate.mean():.2f}%\np={pval:.2e}',
                transform=ax.transAxes, ha='center', va='top', fontsize=9,
                bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.8))

    # Panel F: Summary diagram
    ax = axes[1, 2]
    ax.axis('off')

    summary_text = """Feature Interleaving Analysis
═══════════════════════════════

Normal pattern (contiguous):
[Telo][Telo]---[ITS][ITS]---[TAR1]
  Few transitions between types

Tumor pattern (interleaved):
[Telo][ITS][Telo][ITS][TAR1][ITS]
  Many transitions = recombination

Higher interleaving suggests
BIR-mediated sequence mixing.
"""
    ax.text(0.1, 0.9, summary_text, transform=ax.transAxes, fontsize=10,
            va='top', ha='left', family='monospace',
            bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.5))

    plt.tight_layout()
    fig.savefig(output_path, bbox_inches='tight')
    png_path = output_path.replace('.svg', '.png')
    fig.savefig(png_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved interleaving plot: {output_path} and {png_path}")


def load_its_contiguity_data(
    read_assignments: pd.DataFrame,
    bed_prefix: str,
    database: str,
    smoothness: str,
    min_feature_length: int
) -> pd.DataFrame:
    """Load individual ITS element sizes to analyze contiguity/fragmentation.

    Returns DataFrame with one row per ITS element:
        read, sample, cluster, its_start, its_end, its_size
    """
    samples = read_assignments['sample'].unique()
    reads_needed = set(read_assignments['read'])
    read_info = read_assignments.set_index('read')[['sample', 'cluster', 'read_length']].to_dict('index')

    print(f"  Loading ITS element sizes for {len(reads_needed)} reads...")

    elements = []

    for sample in samples:
        base_path = f"{bed_prefix}/{sample}/telogator/1/KaryoScope/{database}/{sample}.telogator.1.{database}"

        subtelomeric_paths = [
            f"{base_path}.subtelomeric.{smoothness}.features.bed.gz",
            f"{base_path}.subtelomeric.{smoothness}.features.bed",
        ]

        subtelomeric_path = None
        for p in subtelomeric_paths:
            if os.path.exists(p):
                subtelomeric_path = p
                break

        if not subtelomeric_path:
            continue

        open_func = gzip.open if subtelomeric_path.endswith('.gz') else open
        mode = 'rt' if subtelomeric_path.endswith('.gz') else 'r'

        with open_func(subtelomeric_path, mode) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    read_id = parts[0]
                    if read_id not in reads_needed:
                        continue
                    start, end = int(parts[1]), int(parts[2])
                    feature = parts[3]
                    length = end - start

                    if length < min_feature_length:
                        continue

                    if feature == 'ITS' and read_id in read_info:
                        info = read_info[read_id]
                        elements.append({
                            'read': read_id,
                            'sample': info['sample'],
                            'cluster': info['cluster'],
                            'read_length': info['read_length'],
                            'its_start': start,
                            'its_end': end,
                            'its_size': length
                        })

    df = pd.DataFrame(elements)
    print(f"  Loaded {len(df)} ITS elements from {df['read'].nunique() if len(df) > 0 else 0} reads")
    return df


def analyze_its_contiguity(
    its_elements: pd.DataFrame,
    cluster_analysis: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Analyze ITS element contiguity at read and cluster level.

    Returns:
        read_contiguity: per-read ITS fragmentation metrics
        cluster_contiguity: per-cluster ITS fragmentation metrics
    """
    if its_elements.empty:
        return pd.DataFrame(), pd.DataFrame()

    # Add enrichment info
    cluster_enrichment = cluster_analysis.set_index('cluster_id')['enrichment'].to_dict()
    its_elements = its_elements.copy()
    its_elements['enrichment'] = its_elements['cluster'].map(cluster_enrichment)

    # Per-read metrics
    read_stats = its_elements.groupby('read').agg({
        'its_size': ['count', 'sum', 'mean', 'max', 'min', 'std'],
        'cluster': 'first',
        'sample': 'first',
        'read_length': 'first',
        'enrichment': 'first'
    })
    read_stats.columns = ['its_count', 'its_total_bp', 'its_mean_size', 'its_max_size',
                          'its_min_size', 'its_size_std', 'cluster', 'sample', 'read_length', 'enrichment']
    read_stats = read_stats.reset_index()

    # Fragmentation index: lower = more fragmented (smaller elements)
    read_stats['fragmentation_index'] = read_stats['its_mean_size']

    # Per-cluster metrics
    cluster_stats = read_stats.groupby('cluster').agg({
        'its_count': ['sum', 'mean'],
        'its_total_bp': 'sum',
        'its_mean_size': 'mean',
        'its_max_size': 'mean',
        'fragmentation_index': 'mean',
        'read_length': 'sum',
        'enrichment': 'first'
    })
    cluster_stats.columns = ['its_elements_total', 'its_elements_per_read',
                             'its_bp_total', 'mean_element_size', 'mean_max_element_size',
                             'fragmentation_index', 'total_read_length', 'enrichment']
    cluster_stats = cluster_stats.reset_index()

    # Add cluster size
    cluster_stats = cluster_stats.merge(
        cluster_analysis[['cluster_id', 'size', 'odds_ratio']],
        left_on='cluster', right_on='cluster_id', how='left'
    )

    return read_stats, cluster_stats


def plot_its_contiguity(
    read_contiguity: pd.DataFrame,
    cluster_contiguity: pd.DataFrame,
    output_path: str,
    dark_mode: bool = False
):
    """Plot ITS contiguity/fragmentation analysis."""
    plt.rcParams['svg.fonttype'] = 'none'

    if dark_mode:
        plt.style.use('dark_background')
        text_color = 'white'
        bbox_color = '#333333'
    else:
        text_color = 'black'
        bbox_color = 'white'

    colors = {'Normal-enriched': '#3b82f6', 'Tumor-enriched': '#ef4444', 'mixed': '#888888'}

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Filter to Normal and Tumor only
    read_df = read_contiguity[read_contiguity['enrichment'].isin(['Normal-enriched', 'Tumor-enriched'])]
    cluster_df = cluster_contiguity[cluster_contiguity['enrichment'].isin(['Normal-enriched', 'Tumor-enriched'])]

    # Panel A: Distribution of ITS element sizes (element-level)
    ax = axes[0, 0]
    for enrichment in ['Normal-enriched', 'Tumor-enriched']:
        # Get all element sizes for reads in this enrichment
        reads_in_group = read_contiguity[read_contiguity['enrichment'] == enrichment]['read']
        # We need to go back to the element-level data - use read_contiguity's its_mean_size as proxy
        subset = read_df[read_df['enrichment'] == enrichment]['its_mean_size']
        label = 'Normal' if 'Normal' in enrichment else 'Tumor'
        ax.hist(subset, bins=50, alpha=0.6, label=f'{label} (n={len(subset)} reads)',
                color=colors[enrichment], density=True)

    ax.set_xlabel('Mean ITS element size (bp)')
    ax.set_ylabel('Density')
    ax.set_title('A. ITS Element Size Distribution')
    ax.legend(fontsize=8)

    # Panel B: Boxplot of mean ITS element size by enrichment (read-level)
    ax = axes[0, 1]
    normal_sizes = read_df[read_df['enrichment'] == 'Normal-enriched']['its_mean_size']
    tumor_sizes = read_df[read_df['enrichment'] == 'Tumor-enriched']['its_mean_size']

    bp = ax.boxplot([normal_sizes, tumor_sizes], labels=['Normal', 'Tumor'], patch_artist=True)
    bp['boxes'][0].set_facecolor(colors['Normal-enriched'])
    bp['boxes'][1].set_facecolor(colors['Tumor-enriched'])
    for box in bp['boxes']:
        box.set_alpha(0.7)

    ax.set_ylabel('Mean ITS element size (bp)')
    ax.set_title('B. ITS Element Size by Group')

    # Add stats
    if len(normal_sizes) > 1 and len(tumor_sizes) > 1:
        _, pval = mannwhitneyu(normal_sizes, tumor_sizes)
        fc = normal_sizes.mean() / tumor_sizes.mean() if tumor_sizes.mean() > 0 else np.nan
        ax.text(0.5, 0.95, f'Normal: {normal_sizes.mean():.1f} bp\nTumor: {tumor_sizes.mean():.1f} bp\n'
                          f'FC={fc:.2f}x, p={pval:.2e}',
                transform=ax.transAxes, ha='center', va='top', fontsize=9,
                bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.8))

    # Panel C: Number of ITS elements per read
    ax = axes[0, 2]
    normal_counts = read_df[read_df['enrichment'] == 'Normal-enriched']['its_count']
    tumor_counts = read_df[read_df['enrichment'] == 'Tumor-enriched']['its_count']

    bp = ax.boxplot([normal_counts, tumor_counts], labels=['Normal', 'Tumor'], patch_artist=True)
    bp['boxes'][0].set_facecolor(colors['Normal-enriched'])
    bp['boxes'][1].set_facecolor(colors['Tumor-enriched'])
    for box in bp['boxes']:
        box.set_alpha(0.7)

    ax.set_ylabel('ITS elements per read')
    ax.set_title('C. ITS Element Count per Read')

    if len(normal_counts) > 1 and len(tumor_counts) > 1:
        _, pval = mannwhitneyu(normal_counts, tumor_counts)
        ax.text(0.5, 0.95, f'Normal: {normal_counts.mean():.2f}\nTumor: {tumor_counts.mean():.2f}\np={pval:.2e}',
                transform=ax.transAxes, ha='center', va='top', fontsize=9,
                bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.8))

    # Panel D: Cluster-level mean element size vs enrichment
    ax = axes[1, 0]
    for enrichment in ['Normal-enriched', 'Tumor-enriched']:
        subset = cluster_df[cluster_df['enrichment'] == enrichment]
        label = 'Normal' if 'Normal' in enrichment else 'Tumor'
        ax.scatter(np.log2(subset['odds_ratio'].replace(0, np.nan)),
                  subset['mean_element_size'],
                  s=subset['size'] / 3, c=colors[enrichment], alpha=0.6,
                  label=f'{label} (n={len(subset)})')

    ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('log2(Odds Ratio)')
    ax.set_ylabel('Mean ITS element size (bp)')
    ax.set_title('D. Element Size vs Enrichment')
    ax.legend(fontsize=8)

    # Correlation
    valid = cluster_df.dropna(subset=['odds_ratio', 'mean_element_size'])
    valid = valid[valid['odds_ratio'] > 0]
    if len(valid) > 5:
        rho, pval = spearmanr(np.log2(valid['odds_ratio']), valid['mean_element_size'])
        ax.text(0.02, 0.98, f'r = {rho:.2f}, p = {pval:.2e}',
                transform=ax.transAxes, va='top', fontsize=9,
                bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.8))

    # Panel E: Cluster-level boxplot
    ax = axes[1, 1]
    normal_cluster = cluster_df[cluster_df['enrichment'] == 'Normal-enriched']['mean_element_size']
    tumor_cluster = cluster_df[cluster_df['enrichment'] == 'Tumor-enriched']['mean_element_size']

    bp = ax.boxplot([normal_cluster, tumor_cluster], labels=['Normal', 'Tumor'], patch_artist=True)
    bp['boxes'][0].set_facecolor(colors['Normal-enriched'])
    bp['boxes'][1].set_facecolor(colors['Tumor-enriched'])
    for box in bp['boxes']:
        box.set_alpha(0.7)

    ax.set_ylabel('Mean ITS element size (bp)')
    ax.set_title('E. Cluster-Level Element Size')

    if len(normal_cluster) > 1 and len(tumor_cluster) > 1:
        _, pval = mannwhitneyu(normal_cluster, tumor_cluster)
        fc = normal_cluster.mean() / tumor_cluster.mean() if tumor_cluster.mean() > 0 else np.nan
        ax.text(0.5, 0.95, f'Normal: {normal_cluster.mean():.1f} bp\nTumor: {tumor_cluster.mean():.1f} bp\n'
                          f'FC={fc:.2f}x, p={pval:.2e}',
                transform=ax.transAxes, ha='center', va='top', fontsize=9,
                bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.8))

    # Panel F: Max element size comparison (largest ITS block)
    ax = axes[1, 2]
    normal_max = read_df[read_df['enrichment'] == 'Normal-enriched']['its_max_size']
    tumor_max = read_df[read_df['enrichment'] == 'Tumor-enriched']['its_max_size']

    bp = ax.boxplot([normal_max, tumor_max], labels=['Normal', 'Tumor'], patch_artist=True)
    bp['boxes'][0].set_facecolor(colors['Normal-enriched'])
    bp['boxes'][1].set_facecolor(colors['Tumor-enriched'])
    for box in bp['boxes']:
        box.set_alpha(0.7)

    ax.set_ylabel('Max ITS element size (bp)')
    ax.set_title('F. Largest ITS Block per Read')

    if len(normal_max) > 1 and len(tumor_max) > 1:
        _, pval = mannwhitneyu(normal_max, tumor_max)
        fc = normal_max.mean() / tumor_max.mean() if tumor_max.mean() > 0 else np.nan
        ax.text(0.5, 0.95, f'Normal: {normal_max.mean():.1f} bp\nTumor: {tumor_max.mean():.1f} bp\n'
                          f'FC={fc:.2f}x, p={pval:.2e}',
                transform=ax.transAxes, ha='center', va='top', fontsize=9,
                bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.8))

    plt.tight_layout()
    fig.savefig(output_path, bbox_inches='tight')
    png_path = output_path.replace('.svg', '.png')
    fig.savefig(png_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved ITS contiguity plot: {output_path} and {png_path}")


def load_read_level_features(
    read_assignments: pd.DataFrame,
    bed_prefix: str,
    database: str,
    smoothness: str,
    min_feature_length: int
) -> pd.DataFrame:
    """Load all subtelomeric features at read level for pairwise analysis.

    Returns DataFrame with columns: read, sample, cluster, enrichment, read_length,
                                    its_bp, tar1_bp, canonical_bp, noncanonical_bp,
                                    its_bp_per_kb, tar1_bp_per_kb, canonical_bp_per_kb, noncanonical_bp_per_kb
    """
    samples = read_assignments['sample'].unique()
    reads_needed = set(read_assignments['read'])
    read_info = read_assignments.set_index('read')[['sample', 'cluster', 'read_length']].to_dict('index')

    print(f"  Loading all subtelomeric features for {len(reads_needed)} reads...")

    # Feature totals per read
    read_features = defaultdict(lambda: {'its': 0, 'tar1': 0, 'canonical': 0, 'noncanonical': 0})

    for sample in samples:
        base_path = f"{bed_prefix}/{sample}/telogator/1/KaryoScope/{database}/{sample}.telogator.1.{database}"

        subtelomeric_paths = [
            f"{base_path}.subtelomeric.{smoothness}.features.bed.gz",
            f"{base_path}.subtelomeric.{smoothness}.features.bed",
        ]

        subtelomeric_path = None
        for p in subtelomeric_paths:
            if os.path.exists(p):
                subtelomeric_path = p
                break

        if not subtelomeric_path:
            continue

        open_func = gzip.open if subtelomeric_path.endswith('.gz') else open
        mode = 'rt' if subtelomeric_path.endswith('.gz') else 'r'

        with open_func(subtelomeric_path, mode) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    read_id = parts[0]
                    if read_id not in reads_needed:
                        continue
                    start, end = int(parts[1]), int(parts[2])
                    feature = parts[3]
                    length = end - start

                    if length < min_feature_length:
                        continue

                    if feature == 'ITS':
                        read_features[read_id]['its'] += length
                    elif feature == 'TAR1':
                        read_features[read_id]['tar1'] += length
                    elif feature == 'canonical_telomere':
                        read_features[read_id]['canonical'] += length
                    elif feature == 'noncanonical_telomere':
                        read_features[read_id]['noncanonical'] += length

    # Build DataFrame
    rows = []
    for read_id in reads_needed:
        if read_id not in read_info:
            continue
        info = read_info[read_id]
        feats = read_features.get(read_id, {'its': 0, 'tar1': 0, 'canonical': 0, 'noncanonical': 0})
        read_len = info['read_length']
        read_len_kb = read_len / 1000 if read_len > 0 else 1

        rows.append({
            'read': read_id,
            'sample': info['sample'],
            'cluster': info['cluster'],
            'read_length': read_len,
            'its_bp': feats['its'],
            'tar1_bp': feats['tar1'],
            'canonical_bp': feats['canonical'],
            'noncanonical_bp': feats['noncanonical'],
            'its_bp_per_kb': feats['its'] / read_len_kb,
            'tar1_bp_per_kb': feats['tar1'] / read_len_kb,
            'canonical_bp_per_kb': feats['canonical'] / read_len_kb,
            'noncanonical_bp_per_kb': feats['noncanonical'] / read_len_kb,
        })

    df = pd.DataFrame(rows)
    print(f"  Loaded features for {len(df)} reads")
    return df


def plot_pairwise_features(
    read_features_df: pd.DataFrame,
    cluster_analysis: pd.DataFrame,
    output_path: str,
    dark_mode: bool = False
):
    """Create pairwise scatter plots of feature bp/kb at read level."""
    if read_features_df.empty:
        print("  No data to plot")
        return

    plt.rcParams['svg.fonttype'] = 'none'

    if dark_mode:
        plt.style.use('dark_background')
        text_color = 'white'
        bbox_color = '#333333'
    else:
        text_color = 'black'
        bbox_color = 'white'

    # Add enrichment to read features
    cluster_enrichment = cluster_analysis.set_index('cluster_id')['enrichment'].to_dict()
    read_features_df = read_features_df.copy()
    read_features_df['enrichment'] = read_features_df['cluster'].map(cluster_enrichment)

    # Filter to Normal and Tumor only
    df = read_features_df[read_features_df['enrichment'].isin(['Normal-enriched', 'Tumor-enriched'])]

    features = [
        ('its_bp_per_kb', 'ITS'),
        ('tar1_bp_per_kb', 'TAR1'),
        ('canonical_bp_per_kb', 'Canonical'),
        ('noncanonical_bp_per_kb', 'Non-canonical')
    ]

    colors = {'Normal-enriched': '#3b82f6', 'Tumor-enriched': '#ef4444'}

    # 6 pairwise combinations in a 2x3 grid
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    pairs = [
        (0, 1),  # ITS vs TAR1
        (0, 2),  # ITS vs Canonical
        (0, 3),  # ITS vs Non-canonical
        (1, 2),  # TAR1 vs Canonical
        (1, 3),  # TAR1 vs Non-canonical
        (2, 3),  # Canonical vs Non-canonical
    ]

    for idx, (i, j) in enumerate(pairs):
        ax = axes[idx]
        feat_x, label_x = features[i]
        feat_y, label_y = features[j]

        # Sample for plotting (too many points otherwise)
        sample_size = min(5000, len(df))
        plot_df = df.sample(n=sample_size, random_state=42) if len(df) > sample_size else df

        for enrichment in ['Normal-enriched', 'Tumor-enriched']:
            subset = plot_df[plot_df['enrichment'] == enrichment]
            ax.scatter(
                subset[feat_x],
                subset[feat_y],
                c=colors[enrichment],
                alpha=0.3,
                s=5,
                label=f"{enrichment.split('-')[0]} (n={len(df[df['enrichment']==enrichment])})"
            )

        ax.set_xlabel(f'{label_x} bp/kb')
        ax.set_ylabel(f'{label_y} bp/kb')
        ax.set_title(f'{label_x} vs {label_y}')

        if idx == 0:
            ax.legend(loc='upper right', fontsize=8, markerscale=3)

        # Calculate correlations for each group
        corr_text = []
        for enrichment in ['Normal-enriched', 'Tumor-enriched']:
            subset = df[df['enrichment'] == enrichment]
            valid = subset[[feat_x, feat_y]].dropna()
            valid = valid[(valid[feat_x] > 0) | (valid[feat_y] > 0)]  # At least one non-zero
            if len(valid) > 10:
                rho, pval = spearmanr(valid[feat_x], valid[feat_y])
                label = 'N' if 'Normal' in enrichment else 'T'
                corr_text.append(f'{label}: r={rho:.2f}')

        if corr_text:
            ax.text(0.02, 0.98, '\n'.join(corr_text), transform=ax.transAxes,
                    va='top', fontsize=8, bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.8))

    plt.tight_layout()
    fig.savefig(output_path, bbox_inches='tight')
    png_path = output_path.replace('.svg', '.png')
    fig.savefig(png_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved pairwise feature plot: {output_path} and {png_path}")


def analyze_pairwise_correlations(
    read_features_df: pd.DataFrame,
    cluster_analysis: pd.DataFrame,
    output_path: str
):
    """Calculate and save pairwise correlation statistics."""
    # Add enrichment
    cluster_enrichment = cluster_analysis.set_index('cluster_id')['enrichment'].to_dict()
    df = read_features_df.copy()
    df['enrichment'] = df['cluster'].map(cluster_enrichment)

    features = ['its_bp_per_kb', 'tar1_bp_per_kb', 'canonical_bp_per_kb', 'noncanonical_bp_per_kb']
    feature_names = {'its_bp_per_kb': 'ITS', 'tar1_bp_per_kb': 'TAR1',
                     'canonical_bp_per_kb': 'Canonical', 'noncanonical_bp_per_kb': 'Non-canonical'}

    results = []
    for i, feat1 in enumerate(features):
        for feat2 in features[i+1:]:
            for enrichment in ['Normal-enriched', 'Tumor-enriched', 'All']:
                if enrichment == 'All':
                    subset = df
                else:
                    subset = df[df['enrichment'] == enrichment]

                valid = subset[[feat1, feat2]].dropna()
                valid = valid[(valid[feat1] > 0) | (valid[feat2] > 0)]

                if len(valid) > 10:
                    rho, pval = spearmanr(valid[feat1], valid[feat2])
                else:
                    rho, pval = np.nan, np.nan

                results.append({
                    'feature_1': feature_names[feat1],
                    'feature_2': feature_names[feat2],
                    'enrichment': enrichment,
                    'n_reads': len(valid),
                    'spearman_rho': rho,
                    'p_value': pval
                })

    results_df = pd.DataFrame(results)
    results_df.to_csv(output_path, sep='\t', index=False)
    return results_df


def aggregate_cluster_all_features(
    read_level_features: pd.DataFrame,
    cluster_analysis: pd.DataFrame
) -> pd.DataFrame:
    """Aggregate read-level features to cluster level for all feature types."""
    # Add enrichment to read features
    cluster_enrichment = cluster_analysis.set_index('cluster_id')['enrichment'].to_dict()
    df = read_level_features.copy()
    df['enrichment'] = df['cluster'].map(cluster_enrichment)

    # Aggregate by cluster
    cluster_stats = df.groupby('cluster').agg({
        'read_length': ['sum', 'mean', 'count'],
        'its_bp': 'sum',
        'tar1_bp': 'sum',
        'canonical_bp': 'sum',
        'noncanonical_bp': 'sum',
        'enrichment': 'first'
    })

    # Flatten column names
    cluster_stats.columns = ['_'.join(col).strip() if col[1] else col[0]
                             for col in cluster_stats.columns.values]
    cluster_stats = cluster_stats.reset_index()

    # Calculate bp per kb for each feature
    cluster_stats['its_bp_per_kb'] = cluster_stats['its_bp_sum'] / (cluster_stats['read_length_sum'] / 1000)
    cluster_stats['tar1_bp_per_kb'] = cluster_stats['tar1_bp_sum'] / (cluster_stats['read_length_sum'] / 1000)
    cluster_stats['canonical_bp_per_kb'] = cluster_stats['canonical_bp_sum'] / (cluster_stats['read_length_sum'] / 1000)
    cluster_stats['noncanonical_bp_per_kb'] = cluster_stats['noncanonical_bp_sum'] / (cluster_stats['read_length_sum'] / 1000)

    # Merge with cluster_analysis to get odds_ratio etc
    cluster_stats = cluster_stats.merge(
        cluster_analysis[['cluster_id', 'odds_ratio', 'size']],
        left_on='cluster', right_on='cluster_id', how='left'
    )
    cluster_stats['log2_or'] = np.log2(cluster_stats['odds_ratio'].replace(0, np.nan))

    return cluster_stats


def plot_cluster_pairwise_grid(
    cluster_stats: pd.DataFrame,
    output_path: str,
    dark_mode: bool = False
):
    """Create pairwise scatter plot grid at cluster level for all feature combinations."""
    plt.rcParams['svg.fonttype'] = 'none'

    if dark_mode:
        plt.style.use('dark_background')
        text_color = 'white'
        bbox_color = '#333333'
    else:
        text_color = 'black'
        bbox_color = 'white'

    # Filter to Normal and Tumor only
    df = cluster_stats[cluster_stats['enrichment_first'].isin(['Normal-enriched', 'Tumor-enriched'])].copy()

    features = [
        ('its_bp_per_kb', 'ITS'),
        ('tar1_bp_per_kb', 'TAR1'),
        ('canonical_bp_per_kb', 'Canonical'),
        ('noncanonical_bp_per_kb', 'Non-canonical')
    ]

    colors = {'Normal-enriched': '#3b82f6', 'Tumor-enriched': '#ef4444'}

    # 6 pairwise combinations in a 2x3 grid
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    pairs = [
        (0, 1),  # ITS vs TAR1
        (0, 2),  # ITS vs Canonical
        (0, 3),  # ITS vs Non-canonical
        (1, 2),  # TAR1 vs Canonical
        (1, 3),  # TAR1 vs Non-canonical
        (2, 3),  # Canonical vs Non-canonical
    ]

    for idx, (i, j) in enumerate(pairs):
        ax = axes[idx]
        feat_x, label_x = features[i]
        feat_y, label_y = features[j]

        for enrichment in ['Normal-enriched', 'Tumor-enriched']:
            subset = df[df['enrichment_first'] == enrichment]
            label = 'Normal' if 'Normal' in enrichment else 'Tumor'
            ax.scatter(
                subset[feat_x],
                subset[feat_y],
                c=colors[enrichment],
                alpha=0.6,
                s=subset['size'] / 3,  # Size by cluster size
                label=f"{label} (n={len(subset)})"
            )

        ax.set_xlabel(f'{label_x} bp/kb')
        ax.set_ylabel(f'{label_y} bp/kb')
        ax.set_title(f'{label_x} vs {label_y}')

        if idx == 0:
            ax.legend(loc='upper right', fontsize=8)

        # Calculate correlations for each group
        corr_text = []
        for enrichment in ['Normal-enriched', 'Tumor-enriched']:
            subset = df[df['enrichment_first'] == enrichment]
            valid = subset[[feat_x, feat_y]].dropna()
            if len(valid) > 5:
                rho, pval = spearmanr(valid[feat_x], valid[feat_y])
                label = 'N' if 'Normal' in enrichment else 'T'
                sig = '*' if pval < 0.05 else ''
                corr_text.append(f'{label}: r={rho:.2f}{sig}')

        if corr_text:
            ax.text(0.02, 0.98, '\n'.join(corr_text), transform=ax.transAxes,
                    va='top', fontsize=9, bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.8))

    plt.tight_layout()
    fig.savefig(output_path, bbox_inches='tight')
    png_path = output_path.replace('.svg', '.png')
    fig.savefig(png_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved cluster pairwise grid: {output_path} and {png_path}")


def analyze_cluster_pairwise_correlations(
    cluster_stats: pd.DataFrame,
    output_path: str
) -> pd.DataFrame:
    """Calculate pairwise correlations at cluster level."""
    df = cluster_stats.copy()

    features = ['its_bp_per_kb', 'tar1_bp_per_kb', 'canonical_bp_per_kb', 'noncanonical_bp_per_kb']
    feature_names = {
        'its_bp_per_kb': 'ITS',
        'tar1_bp_per_kb': 'TAR1',
        'canonical_bp_per_kb': 'Canonical',
        'noncanonical_bp_per_kb': 'Non-canonical'
    }

    results = []
    for i, feat1 in enumerate(features):
        for feat2 in features[i+1:]:
            for enrichment in ['Normal-enriched', 'Tumor-enriched', 'All']:
                if enrichment == 'All':
                    subset = df
                else:
                    subset = df[df['enrichment_first'] == enrichment]

                valid = subset[[feat1, feat2]].dropna()
                if len(valid) > 5:
                    rho, pval = spearmanr(valid[feat1], valid[feat2])
                else:
                    rho, pval = np.nan, np.nan

                results.append({
                    'feature_1': feature_names[feat1],
                    'feature_2': feature_names[feat2],
                    'enrichment': enrichment,
                    'n_clusters': len(valid),
                    'spearman_rho': rho,
                    'p_value': pval
                })

    results_df = pd.DataFrame(results)
    results_df.to_csv(output_path, sep='\t', index=False)
    return results_df


def calculate_cluster_solitary_tar1(
    read_assignments: pd.DataFrame,
    cluster_analysis: pd.DataFrame,
    adjacency_data: Dict[str, Tuple[int, int]]
) -> pd.DataFrame:
    """Calculate solitary TAR1 statistics per cluster."""
    results = []

    for _, cluster_row in cluster_analysis.iterrows():
        cluster_id = cluster_row['cluster_id']
        cluster_reads = read_assignments[read_assignments['cluster'] == cluster_id]['read']

        total_solitary = 0
        total_tar1 = 0
        reads_with_tar1 = 0
        reads_with_solitary_tar1 = 0

        for read_id in cluster_reads:
            solitary, total = adjacency_data.get(read_id, (0, 0))
            total_solitary += solitary
            total_tar1 += total
            if total > 0:
                reads_with_tar1 += 1
            if solitary > 0:
                reads_with_solitary_tar1 += 1

        pct_solitary_tar1 = (total_solitary / total_tar1 * 100) if total_tar1 > 0 else 0
        pct_reads_with_solitary = (reads_with_solitary_tar1 / reads_with_tar1 * 100) if reads_with_tar1 > 0 else 0

        results.append({
            'cluster_id': cluster_id,
            'size': cluster_row['size'],
            'enrichment': cluster_row['enrichment'],
            'odds_ratio': cluster_row['odds_ratio'],
            'log2_or': np.log2(cluster_row['odds_ratio']) if cluster_row['odds_ratio'] > 0 else 0,
            'total_tar1_elements': total_tar1,
            'solitary_tar1_elements': total_solitary,
            'pct_solitary_tar1': pct_solitary_tar1,
            'reads_with_tar1': reads_with_tar1,
            'reads_with_solitary_tar1': reads_with_solitary_tar1,
            'pct_reads_with_solitary': pct_reads_with_solitary
        })

    return pd.DataFrame(results)


def plot_solitary_tar1(solitary_df: pd.DataFrame, output_path: str,
                       dark_mode: bool = False):
    """Plot solitary TAR1 analysis results."""
    if solitary_df.empty:
        print("  No data to plot")
        return

    plt.rcParams['svg.fonttype'] = 'none'

    if dark_mode:
        plt.style.use('dark_background')
        text_color = 'white'
        bbox_color = '#333333'
    else:
        text_color = 'black'
        bbox_color = 'white'

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    colors = {'Normal-enriched': '#3b82f6', 'Tumor-enriched': '#ef4444', 'mixed': '#888888'}

    # Filter to Normal and Tumor only
    plot_data = solitary_df[solitary_df['enrichment'].isin(['Normal-enriched', 'Tumor-enriched'])]
    order = ['Normal-enriched', 'Tumor-enriched']

    # Panel 1: Boxplot of % solitary TAR1
    ax = axes[0]
    data_by_group = [plot_data[plot_data['enrichment'] == e]['pct_solitary_tar1'].values for e in order]

    bp = ax.boxplot(data_by_group, labels=['Normal', 'Tumor'], patch_artist=True)
    for patch, enrichment in zip(bp['boxes'], order):
        patch.set_facecolor(colors.get(enrichment))
        patch.set_alpha(0.7)

    ax.set_ylabel('% Solitary TAR1 elements')
    ax.set_title('Solitary TAR1 by Enrichment')

    # Add p-value
    normal_vals = plot_data[plot_data['enrichment'] == 'Normal-enriched']['pct_solitary_tar1'].values
    tumor_vals = plot_data[plot_data['enrichment'] == 'Tumor-enriched']['pct_solitary_tar1'].values
    if len(normal_vals) > 1 and len(tumor_vals) > 1:
        _, pval = mannwhitneyu(normal_vals, tumor_vals, alternative='two-sided')
        ax.text(0.5, 0.95, f'p = {pval:.2e}', transform=ax.transAxes, ha='center', va='top', fontsize=10)

    # Panel 2: Scatter plot vs enrichment
    ax = axes[1]
    for enrichment in ['Normal-enriched', 'Tumor-enriched']:
        subset = solitary_df[solitary_df['enrichment'] == enrichment]
        if len(subset) > 0:
            ax.scatter(
                subset['log2_or'],
                subset['pct_solitary_tar1'],
                s=subset['size'] / 5,
                c=colors.get(enrichment),
                alpha=0.6,
                label=f"{enrichment} (n={len(subset)})"
            )
            for _, row in subset.iterrows():
                ax.annotate(
                    str(int(row['cluster_id'])),
                    (row['log2_or'], row['pct_solitary_tar1']),
                    fontsize=6, alpha=0.7, color=text_color,
                    ha='center', va='bottom',
                    xytext=(0, 3), textcoords='offset points'
                )

    ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('log2(Odds Ratio)')
    ax.set_ylabel('% Solitary TAR1')
    ax.set_title('Solitary TAR1 vs Enrichment')
    ax.legend(loc='upper right', fontsize=8)

    # Add correlation
    valid = solitary_df.dropna(subset=['log2_or', 'pct_solitary_tar1'])
    valid = valid[valid['total_tar1_elements'] > 0]  # Only clusters with TAR1
    if len(valid) > 2:
        rho, pval = spearmanr(valid['log2_or'], valid['pct_solitary_tar1'])
        ax.text(0.02, 0.02, f'Spearman r = {rho:.3f}\np = {pval:.2e}',
                transform=ax.transAxes, va='bottom', fontsize=10,
                bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.8))

    # Panel 3: Summary schematic
    ax = axes[2]
    ax.axis('off')

    summary_text = """
    Solitary TAR1 Analysis
    ══════════════════════════════════

    Hypothesis: BIR causes ITS loss while
    TAR1 remains, creating "orphan" TAR1.

    Normal cells: TAR1-ITS pairs (adjacent)
    Tumor cells: Solitary TAR1 (ITS lost)

    Results:
      Normal-enriched: {:.1f}% solitary TAR1
      Tumor-enriched:  {:.1f}% solitary TAR1

    Interpretation:
      Higher % solitary TAR1 in tumors
      supports BIR-mediated ITS loss.
    """.format(
        normal_vals.mean() if len(normal_vals) > 0 else 0,
        tumor_vals.mean() if len(tumor_vals) > 0 else 0
    )

    ax.text(0.1, 0.95, summary_text, transform=ax.transAxes, fontsize=10,
            va='top', ha='left', family='monospace',
            bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.5))

    plt.tight_layout()
    fig.savefig(output_path, bbox_inches='tight')
    png_path = output_path.replace('.svg', '.png')
    fig.savefig(png_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved solitary TAR1 plot: {output_path} and {png_path}")


def calculate_cluster_its_tar1(
    read_assignments: pd.DataFrame,
    cluster_analysis: pd.DataFrame,
    read_features: Dict[str, Tuple[int, int, int, int, int]]
) -> pd.DataFrame:
    """Calculate ITS/TAR1 statistics for each cluster."""
    results = []

    for _, cluster_row in cluster_analysis.iterrows():
        cluster_id = cluster_row['cluster_id']
        cluster_reads = read_assignments[read_assignments['cluster'] == cluster_id]['read']

        # Aggregate ITS/TAR1 stats for this cluster
        its_counts = []
        its_bps = []
        tar1_counts = []
        tar1_bps = []
        read_lengths = []
        reads_with_its = 0
        reads_with_tar1 = 0

        for read_id in cluster_reads:
            if read_id in read_features:
                its_count, its_bp, tar1_count, tar1_bp, read_length = read_features[read_id]
                its_counts.append(its_count)
                its_bps.append(its_bp)
                tar1_counts.append(tar1_count)
                tar1_bps.append(tar1_bp)
                read_lengths.append(read_length)

                if its_count > 0:
                    reads_with_its += 1
                if tar1_count > 0:
                    reads_with_tar1 += 1

        n_reads = len(its_counts)
        if n_reads == 0:
            continue

        total_read_length = sum(read_lengths)

        results.append({
            'cluster_id': cluster_id,
            'size': cluster_row['size'],
            'enrichment': cluster_row['enrichment'],
            'odds_ratio': cluster_row['odds_ratio'],
            'q_value': cluster_row['q_value'],
            'log2_or': np.log2(cluster_row['odds_ratio']) if cluster_row['odds_ratio'] > 0 else np.nan,
            # ITS stats
            'its_count_total': sum(its_counts),
            'its_count_mean': np.mean(its_counts),
            'its_bp_total': sum(its_bps),
            'its_bp_mean': np.mean(its_bps),
            'its_bp_per_kb': sum(its_bps) / (total_read_length / 1e3) if total_read_length > 0 else 0,
            'pct_reads_with_its': 100 * reads_with_its / n_reads,
            # TAR1 stats
            'tar1_count_total': sum(tar1_counts),
            'tar1_count_mean': np.mean(tar1_counts),
            'tar1_bp_total': sum(tar1_bps),
            'tar1_bp_mean': np.mean(tar1_bps),
            'tar1_bp_per_kb': sum(tar1_bps) / (total_read_length / 1e3) if total_read_length > 0 else 0,
            'pct_reads_with_tar1': 100 * reads_with_tar1 / n_reads,
            # Read stats
            'n_reads_with_features': n_reads,
            'total_read_length': total_read_length,
            'mean_read_length': np.mean(read_lengths)
        })

    return pd.DataFrame(results)


def compare_by_enrichment(cluster_its_tar1: pd.DataFrame) -> pd.DataFrame:
    """Compare ITS/TAR1 metrics between Normal-enriched and Tumor-enriched clusters."""
    results = []

    enrichment_groups = cluster_its_tar1.groupby('enrichment')
    metrics = ['its_bp_mean', 'tar1_bp_mean', 'its_bp_per_kb', 'tar1_bp_per_kb',
               'pct_reads_with_its', 'pct_reads_with_tar1', 'its_count_mean', 'tar1_count_mean']

    for metric in metrics:
        row = {'metric': metric}

        for enrichment in ['Normal-enriched', 'Tumor-enriched', 'mixed']:
            if enrichment in enrichment_groups.groups:
                group_data = enrichment_groups.get_group(enrichment)[metric]
                row[f'{enrichment}_mean'] = group_data.mean()
                row[f'{enrichment}_median'] = group_data.median()
                row[f'{enrichment}_n'] = len(group_data)

        # Statistical test between Normal-enriched and Tumor-enriched
        if 'Normal-enriched' in enrichment_groups.groups and 'Tumor-enriched' in enrichment_groups.groups:
            normal_vals = enrichment_groups.get_group('Normal-enriched')[metric].values
            tumor_vals = enrichment_groups.get_group('Tumor-enriched')[metric].values
            if len(normal_vals) > 1 and len(tumor_vals) > 1:
                try:
                    _, pval = mannwhitneyu(normal_vals, tumor_vals, alternative='two-sided')
                    row['p_value'] = pval
                    row['fold_change'] = normal_vals.mean() / tumor_vals.mean() if tumor_vals.mean() > 0 else np.nan
                except ValueError:
                    row['p_value'] = np.nan
                    row['fold_change'] = np.nan

        results.append(row)

    return pd.DataFrame(results)


# --- Analysis 1: Positional relationships ---
def analyze_positional_relationships(read_features: Dict[str, ReadFeatures]) -> List[ITSTar1Pair]:
    """Analyze ITS-TAR1 distances and adjacency patterns."""
    pairs = []

    for read_id, rf in read_features.items():
        its_intervals = rf.get_its_intervals()
        tar1_intervals = rf.get_tar1_intervals()

        if not its_intervals or not tar1_intervals:
            continue

        chromosome = rf.dominant_chromosome()
        arm = rf.dominant_arm()

        # Find all ITS-TAR1 pairs on this read
        for its in its_intervals:
            for tar1 in tar1_intervals:
                # Calculate distance (0 if overlapping)
                if its.end <= tar1.start:
                    distance = tar1.start - its.end
                    its_before = True
                elif tar1.end <= its.start:
                    distance = its.start - tar1.end
                    its_before = False
                else:
                    distance = 0  # Overlapping
                    its_before = its.start < tar1.start

                pairs.append(ITSTar1Pair(
                    read_id=read_id,
                    sample=rf.sample,
                    its_start=its.start,
                    its_end=its.end,
                    tar1_start=tar1.start,
                    tar1_end=tar1.end,
                    distance=distance,
                    its_before_tar1=its_before,
                    chromosome=chromosome,
                    arm=arm
                ))

    return pairs


def summarize_positional_analysis(pairs: List[ITSTar1Pair],
                                  sample_to_group: Dict[str, str],
                                  max_distance: int) -> pd.DataFrame:
    """Summarize positional relationships by group."""
    if not pairs:
        return pd.DataFrame()

    df = pd.DataFrame([{
        'read_id': p.read_id,
        'sample': p.sample,
        'group': sample_to_group.get(p.sample, 'Unknown'),
        'distance': p.distance,
        'its_before_tar1': p.its_before_tar1,
        'adjacent': p.distance <= max_distance,
        'chromosome': p.chromosome,
        'arm': p.arm
    } for p in pairs])

    # Summary by group
    summary = df.groupby('group').agg({
        'read_id': 'nunique',
        'distance': ['mean', 'median', 'std', 'min', 'max'],
        'its_before_tar1': 'mean',
        'adjacent': 'mean'
    }).round(3)

    return summary


# --- Analysis 2: Sample-level quantification ---
def calculate_sample_stats(all_features: Dict[str, Dict[str, ReadFeatures]],
                           sample_to_group: Dict[str, str]) -> List[SampleStats]:
    """Calculate ITS/TAR1 statistics for each sample."""
    results = []

    for sample, read_features in all_features.items():
        group = sample_to_group.get(sample, 'Unknown')

        n_reads = len(read_features)
        total_read_length = sum(rf.read_length for rf in read_features.values())

        its_count = 0
        its_total_bp = 0
        tar1_count = 0
        tar1_total_bp = 0
        reads_with_its = 0
        reads_with_tar1 = 0
        reads_with_both = 0

        for rf in read_features.values():
            its_intervals = rf.get_its_intervals()
            tar1_intervals = rf.get_tar1_intervals()

            its_bp = rf.total_its_bp()
            tar1_bp = rf.total_tar1_bp()

            its_count += len(its_intervals)
            its_total_bp += its_bp
            tar1_count += len(tar1_intervals)
            tar1_total_bp += tar1_bp

            has_its = len(its_intervals) > 0
            has_tar1 = len(tar1_intervals) > 0

            if has_its:
                reads_with_its += 1
            if has_tar1:
                reads_with_tar1 += 1
            if has_its and has_tar1:
                reads_with_both += 1

        results.append(SampleStats(
            sample=sample,
            group=group,
            n_reads=n_reads,
            total_read_length=total_read_length,
            its_count=its_count,
            its_total_bp=its_total_bp,
            tar1_count=tar1_count,
            tar1_total_bp=tar1_total_bp,
            reads_with_its=reads_with_its,
            reads_with_tar1=reads_with_tar1,
            reads_with_both=reads_with_both
        ))

    return results


def perform_sample_comparison(sample_stats: List[SampleStats]) -> Tuple[pd.DataFrame, Dict]:
    """Perform statistical comparison between groups."""
    df = pd.DataFrame([{
        'sample': s.sample,
        'group': s.group,
        'n_reads': s.n_reads,
        'total_mb': s.total_read_length / 1e6,
        'its_count': s.its_count,
        'its_per_mb': s.its_per_mb,
        'its_bp_per_kb': s.its_bp_per_kb,
        'its_total_bp': s.its_total_bp,
        'tar1_count': s.tar1_count,
        'tar1_per_mb': s.tar1_per_mb,
        'tar1_bp_per_kb': s.tar1_bp_per_kb,
        'tar1_total_bp': s.tar1_total_bp,
        'pct_reads_with_its': 100 * s.reads_with_its / s.n_reads if s.n_reads > 0 else 0,
        'pct_reads_with_tar1': 100 * s.reads_with_tar1 / s.n_reads if s.n_reads > 0 else 0,
        'pct_reads_with_both': 100 * s.reads_with_both / s.n_reads if s.n_reads > 0 else 0,
        'its_tar1_ratio': s.its_total_bp / s.tar1_total_bp if s.tar1_total_bp > 0 else np.nan
    } for s in sample_stats])

    # Statistical tests between groups
    stats_results = {}
    groups = sorted(df['group'].unique())
    if len(groups) == 2:
        g1, g2 = groups
        metrics = ['its_per_mb', 'tar1_per_mb', 'its_bp_per_kb', 'pct_reads_with_its', 'pct_reads_with_tar1']

        print(f"\n  Statistical comparison: {g1} vs {g2}")
        for metric in metrics:
            v1 = df[df['group'] == g1][metric].values
            v2 = df[df['group'] == g2][metric].values
            if len(v1) >= 1 and len(v2) >= 1:
                try:
                    stat, pval = mannwhitneyu(v1, v2, alternative='two-sided')
                    stats_results[metric] = {
                        f'{g1}_mean': v1.mean(),
                        f'{g2}_mean': v2.mean(),
                        'p_value': pval,
                        'fold_change': v1.mean() / v2.mean() if v2.mean() > 0 else np.nan
                    }
                    print(f"    {metric}: {g1}={v1.mean():.2f}, {g2}={v2.mean():.2f}, p={pval:.4f}")
                except ValueError:
                    pass

    return df, stats_results


# --- Analysis 3: Chromosome arm specificity ---
def analyze_chromosome_specificity(all_features: Dict[str, Dict[str, ReadFeatures]],
                                   sample_to_group: Dict[str, str]) -> pd.DataFrame:
    """Calculate ITS enrichment per chromosome arm."""
    # Collect stats per (chromosome, arm, group)
    arm_data = defaultdict(lambda: {
        'n_reads': 0, 'its_count': 0, 'its_bp': 0, 'tar1_count': 0, 'tar1_bp': 0,
        'total_read_length': 0
    })

    for sample, read_features in all_features.items():
        group = sample_to_group.get(sample, 'Unknown')

        for rf in read_features.values():
            chromosome = rf.dominant_chromosome()
            arm = rf.dominant_arm()

            if chromosome is None or arm is None:
                continue

            key = (chromosome, arm, group)
            arm_data[key]['n_reads'] += 1
            arm_data[key]['its_count'] += len(rf.get_its_intervals())
            arm_data[key]['its_bp'] += rf.total_its_bp()
            arm_data[key]['tar1_count'] += len(rf.get_tar1_intervals())
            arm_data[key]['tar1_bp'] += rf.total_tar1_bp()
            arm_data[key]['total_read_length'] += rf.read_length

    # Convert to DataFrame
    rows = []
    for (chrom, arm, group), data in arm_data.items():
        rows.append({
            'chromosome': chrom,
            'arm': arm,
            'chr_arm': f"{chrom}{arm}",
            'group': group,
            'n_reads': data['n_reads'],
            'its_count': data['its_count'],
            'its_per_read': data['its_count'] / data['n_reads'] if data['n_reads'] > 0 else 0,
            'its_bp': data['its_bp'],
            'its_bp_per_kb': data['its_bp'] / (data['total_read_length'] / 1e3) if data['total_read_length'] > 0 else 0,
            'tar1_count': data['tar1_count'],
            'tar1_per_read': data['tar1_count'] / data['n_reads'] if data['n_reads'] > 0 else 0,
            'tar1_bp': data['tar1_bp'],
            'tar1_bp_per_kb': data['tar1_bp'] / (data['total_read_length'] / 1e3) if data['total_read_length'] > 0 else 0
        })

    df = pd.DataFrame(rows)
    return df


def calculate_chromosome_fold_changes(chr_df: pd.DataFrame) -> pd.DataFrame:
    """Calculate log2 fold changes per chromosome arm."""
    groups = sorted(chr_df['group'].unique())
    if len(groups) != 2:
        return pd.DataFrame()

    control, treatment = groups[0], groups[1]

    # Pivot to compare groups
    pivot = chr_df.pivot_table(
        index='chr_arm',
        columns='group',
        values=['its_per_read', 'its_bp_per_kb', 'tar1_per_read', 'n_reads'],
        aggfunc='sum'
    )

    # Flatten column names
    pivot.columns = ['_'.join(col).strip() for col in pivot.columns.values]

    # Calculate fold changes
    its_ctrl = pivot[f'its_per_read_{control}'].replace(0, 0.001)
    its_treat = pivot[f'its_per_read_{treatment}'].replace(0, 0.001)
    pivot['its_log2_fc'] = np.log2(its_treat / its_ctrl)

    tar1_ctrl = pivot[f'tar1_per_read_{control}'].replace(0, 0.001)
    tar1_treat = pivot[f'tar1_per_read_{treatment}'].replace(0, 0.001)
    pivot['tar1_log2_fc'] = np.log2(tar1_treat / tar1_ctrl)

    # Reset index
    pivot = pivot.reset_index()

    # Sort by ITS log2 fold change
    pivot = pivot.sort_values('its_log2_fc')

    print(f"\n  Chromosome arms with greatest ITS loss ({treatment} vs {control}):")
    for _, row in pivot.head(10).iterrows():
        print(f"    {row['chr_arm']}: ITS log2FC={row['its_log2_fc']:.2f}")

    return pivot


# --- Visualization functions ---
def plot_distance_distribution(pairs: List[ITSTar1Pair],
                               sample_to_group: Dict[str, str],
                               output_path: str, dark_mode: bool = False):
    """Plot ITS-TAR1 distance distributions by group."""
    if not pairs:
        print("  No ITS-TAR1 pairs to plot")
        return

    if dark_mode:
        plt.style.use('dark_background')

    df = pd.DataFrame([{
        'distance': p.distance,
        'group': sample_to_group.get(p.sample, 'Unknown')
    } for p in pairs])

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    groups = sorted(df['group'].unique())
    colors = {'Normal': '#3b82f6', 'Tumor': '#ef4444'}

    # Histogram by group
    ax1 = axes[0]
    for group in groups:
        subset = df[df['group'] == group]['distance']
        ax1.hist(subset, bins=50, alpha=0.7, label=f"{group} (n={len(subset)})",
                 color=colors.get(group, '#999999'))
    ax1.set_xlabel('ITS-TAR1 Distance (bp)')
    ax1.set_ylabel('Count')
    ax1.set_title('ITS-TAR1 Distance Distribution')
    ax1.legend()

    # Boxplot comparison
    ax2 = axes[1]
    data_by_group = [df[df['group'] == g]['distance'].values for g in groups]
    bp = ax2.boxplot(data_by_group, labels=groups, patch_artist=True)
    for patch, group in zip(bp['boxes'], groups):
        patch.set_facecolor(colors.get(group, '#999999'))
        patch.set_alpha(0.7)
    ax2.set_ylabel('ITS-TAR1 Distance (bp)')
    ax2.set_title('Distance Comparison by Group')

    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved distance plot: {output_path}")


def plot_sample_abundance(sample_stats: List[SampleStats],
                          output_path: str, dark_mode: bool = False):
    """Barplot of ITS/TAR1 abundance per sample."""
    if dark_mode:
        plt.style.use('dark_background')

    df = pd.DataFrame([{
        'sample': s.sample,
        'group': s.group,
        'its_per_mb': s.its_per_mb,
        'tar1_per_mb': s.tar1_per_mb,
        'its_bp_per_kb': s.its_bp_per_kb,
        'tar1_bp_per_kb': s.tar1_bp_per_kb
    } for s in sample_stats])

    # Sort by group then ITS abundance
    df = df.sort_values(['group', 'its_per_mb'], ascending=[True, False])
    df = df.reset_index(drop=True)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    colors = {'Normal': '#3b82f6', 'Tumor': '#ef4444'}
    bar_colors = [colors.get(g, '#999999') for g in df['group']]

    # ITS intervals per Mb
    ax = axes[0, 0]
    ax.bar(range(len(df)), df['its_per_mb'], color=bar_colors, alpha=0.8)
    ax.set_ylabel('ITS intervals per Mb')
    ax.set_title('ITS Interval Density')
    for group in df['group'].unique():
        group_mean = df[df['group'] == group]['its_per_mb'].mean()
        ax.axhline(group_mean, color=colors.get(group, '#999'),
                   linestyle='--', alpha=0.7, label=f'{group} mean')
    ax.legend()

    # TAR1 intervals per Mb
    ax = axes[0, 1]
    ax.bar(range(len(df)), df['tar1_per_mb'], color=bar_colors, alpha=0.8)
    ax.set_ylabel('TAR1 intervals per Mb')
    ax.set_title('TAR1 Interval Density')
    for group in df['group'].unique():
        group_mean = df[df['group'] == group]['tar1_per_mb'].mean()
        ax.axhline(group_mean, color=colors.get(group, '#999'),
                   linestyle='--', alpha=0.7, label=f'{group} mean')
    ax.legend()

    # ITS bp per kb
    ax = axes[1, 0]
    ax.bar(range(len(df)), df['its_bp_per_kb'], color=bar_colors, alpha=0.8)
    ax.set_ylabel('ITS bp per kb')
    ax.set_title('ITS Coverage')
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df['sample'], rotation=45, ha='right')

    # TAR1 bp per kb
    ax = axes[1, 1]
    ax.bar(range(len(df)), df['tar1_bp_per_kb'], color=bar_colors, alpha=0.8)
    ax.set_ylabel('TAR1 bp per kb')
    ax.set_title('TAR1 Coverage')
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df['sample'], rotation=45, ha='right')

    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved abundance plot: {output_path}")


def plot_chromosome_heatmap(chr_fc_df: pd.DataFrame, output_path: str,
                            dark_mode: bool = False):
    """Heatmap of ITS depletion per chromosome arm."""
    if chr_fc_df.empty:
        print("  No chromosome data to plot")
        return

    if dark_mode:
        plt.style.use('dark_background')

    # Natural chromosome sort
    def chr_sort_key(x):
        chrom = x.rstrip('pq')
        arm = x[-1]
        try:
            num = int(chrom.replace('chr', ''))
            return (0, num, arm)
        except ValueError:
            return (1, chrom, arm)

    chr_fc_df = chr_fc_df.copy()
    chr_fc_df['sort_key'] = chr_fc_df['chr_arm'].apply(chr_sort_key)
    chr_fc_df = chr_fc_df.sort_values('sort_key')

    fig, axes = plt.subplots(1, 2, figsize=(10, 12))

    cmap = 'RdBu_r'

    # ITS heatmap
    ax = axes[0]
    its_data = chr_fc_df['its_log2_fc'].values.reshape(-1, 1)
    vmax = max(abs(np.nanmin(its_data)), abs(np.nanmax(its_data)), 1)
    im1 = ax.imshow(its_data, cmap=cmap, aspect='auto', vmin=-vmax, vmax=vmax)
    ax.set_yticks(range(len(chr_fc_df)))
    ax.set_yticklabels(chr_fc_df['chr_arm'])
    ax.set_xticks([0])
    ax.set_xticklabels(['ITS log2FC'])
    ax.set_title('ITS Enrichment\n(red = loss in Tumor)')
    plt.colorbar(im1, ax=ax, label='log2 fold change')

    # TAR1 heatmap
    ax = axes[1]
    tar1_data = chr_fc_df['tar1_log2_fc'].values.reshape(-1, 1)
    vmax = max(abs(np.nanmin(tar1_data)), abs(np.nanmax(tar1_data)), 1)
    im2 = ax.imshow(tar1_data, cmap=cmap, aspect='auto', vmin=-vmax, vmax=vmax)
    ax.set_yticks(range(len(chr_fc_df)))
    ax.set_yticklabels(chr_fc_df['chr_arm'])
    ax.set_xticks([0])
    ax.set_xticklabels(['TAR1 log2FC'])
    ax.set_title('TAR1 Enrichment\n(blue = gain in Tumor)')
    plt.colorbar(im2, ax=ax, label='log2 fold change')

    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved chromosome heatmap: {output_path}")


# --- Cluster mode visualization functions ---
def plot_its_vs_enrichment(cluster_its_tar1: pd.DataFrame, output_path: str,
                           dark_mode: bool = False):
    """Scatter plot of ITS/TAR1 vs cluster enrichment (log2 odds ratio)."""
    if cluster_its_tar1.empty:
        print("  No cluster data to plot")
        return

    # SVG with editable text
    plt.rcParams['svg.fonttype'] = 'none'

    if dark_mode:
        plt.style.use('dark_background')
        text_color = 'white'
        bbox_color = '#333333'
    else:
        text_color = 'black'
        bbox_color = 'white'

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    colors = {'Normal-enriched': '#3b82f6', 'Tumor-enriched': '#ef4444', 'mixed': '#888888'}

    # ITS vs log2(OR) - using per-kb normalization to control for read length
    ax = axes[0]
    for enrichment in ['Normal-enriched', 'Tumor-enriched', 'mixed']:
        subset = cluster_its_tar1[cluster_its_tar1['enrichment'] == enrichment]
        if len(subset) > 0:
            ax.scatter(
                subset['log2_or'],
                subset['its_bp_per_kb'],
                s=subset['size'] / 5,
                c=colors.get(enrichment, '#999'),
                alpha=0.6,
                label=f"{enrichment} (n={len(subset)})"
            )
            # Add cluster labels
            for _, row in subset.iterrows():
                ax.annotate(
                    str(int(row['cluster_id'])),
                    (row['log2_or'], row['its_bp_per_kb']),
                    fontsize=6, alpha=0.7, color=text_color,
                    ha='center', va='bottom',
                    xytext=(0, 3), textcoords='offset points'
                )

    ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('log2(Odds Ratio)')
    ax.set_ylabel('ITS bp per kb read')
    ax.set_title('ITS vs Cluster Enrichment (length-normalized)')
    ax.legend(loc='upper right', fontsize=8)

    # Add correlation in bottom left
    valid = cluster_its_tar1.dropna(subset=['log2_or', 'its_bp_per_kb'])
    if len(valid) > 2:
        rho, pval = spearmanr(valid['log2_or'], valid['its_bp_per_kb'])
        ax.text(0.02, 0.02, f'Spearman r = {rho:.3f}\np = {pval:.2e}',
                transform=ax.transAxes, va='bottom', fontsize=10,
                bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.8))

    # TAR1 vs log2(OR) - using per-kb normalization to control for read length
    ax = axes[1]
    for enrichment in ['Normal-enriched', 'Tumor-enriched', 'mixed']:
        subset = cluster_its_tar1[cluster_its_tar1['enrichment'] == enrichment]
        if len(subset) > 0:
            ax.scatter(
                subset['log2_or'],
                subset['tar1_bp_per_kb'],
                s=subset['size'] / 5,
                c=colors.get(enrichment, '#999'),
                alpha=0.6,
                label=f"{enrichment} (n={len(subset)})"
            )
            # Add cluster labels
            for _, row in subset.iterrows():
                ax.annotate(
                    str(int(row['cluster_id'])),
                    (row['log2_or'], row['tar1_bp_per_kb']),
                    fontsize=6, alpha=0.7, color=text_color,
                    ha='center', va='bottom',
                    xytext=(0, 3), textcoords='offset points'
                )

    ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('log2(Odds Ratio)')
    ax.set_ylabel('TAR1 bp per kb read')
    ax.set_title('TAR1 vs Cluster Enrichment (length-normalized)')
    ax.legend(loc='upper right', fontsize=8)

    # Add correlation in bottom left
    valid = cluster_its_tar1.dropna(subset=['log2_or', 'tar1_bp_per_kb'])
    if len(valid) > 2:
        rho, pval = spearmanr(valid['log2_or'], valid['tar1_bp_per_kb'])
        ax.text(0.02, 0.02, f'Spearman r = {rho:.3f}\np = {pval:.2e}',
                transform=ax.transAxes, va='bottom', fontsize=10,
                bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.8))

    plt.tight_layout()
    # Save both SVG (editable) and PNG (viewable)
    fig.savefig(output_path, bbox_inches='tight')
    png_path = output_path.replace('.svg', '.png')
    fig.savefig(png_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved enrichment scatter: {output_path} and {png_path}")


def plot_cluster_its_distribution(cluster_its_tar1: pd.DataFrame, output_path: str,
                                  dark_mode: bool = False):
    """Box/violin plot of ITS/TAR1 by enrichment direction."""
    if cluster_its_tar1.empty:
        print("  No cluster data to plot")
        return

    # SVG with editable text
    plt.rcParams['svg.fonttype'] = 'none'

    if dark_mode:
        plt.style.use('dark_background')

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    colors = {'Normal-enriched': '#3b82f6', 'Tumor-enriched': '#ef4444'}
    order = ['Normal-enriched', 'Tumor-enriched']
    order = [e for e in order if e in cluster_its_tar1['enrichment'].unique()]

    # Filter out mixed clusters
    plot_data = cluster_its_tar1[cluster_its_tar1['enrichment'].isin(['Normal-enriched', 'Tumor-enriched'])]

    metrics = [
        ('its_bp_per_kb', 'ITS bp/kb read'),
        ('tar1_bp_per_kb', 'TAR1 bp/kb read'),
        ('pct_reads_with_its', '% reads with ITS'),
        ('pct_reads_with_tar1', '% reads with TAR1')
    ]

    for idx, (metric, label) in enumerate(metrics):
        ax = axes[idx // 2, idx % 2]

        data_by_group = [plot_data[plot_data['enrichment'] == e][metric].values
                         for e in order]

        bp = ax.boxplot(data_by_group, labels=order, patch_artist=True)
        for patch, enrichment in zip(bp['boxes'], order):
            patch.set_facecolor(colors.get(enrichment, '#999'))
            patch.set_alpha(0.7)

        ax.set_ylabel(label)
        ax.set_title(label)

        # Add p-value annotation
        if 'Normal-enriched' in order and 'Tumor-enriched' in order:
            normal_vals = plot_data[plot_data['enrichment'] == 'Normal-enriched'][metric].values
            tumor_vals = plot_data[plot_data['enrichment'] == 'Tumor-enriched'][metric].values
            if len(normal_vals) > 1 and len(tumor_vals) > 1:
                try:
                    _, pval = mannwhitneyu(normal_vals, tumor_vals, alternative='two-sided')
                    ax.text(0.5, 0.95, f'p = {pval:.2e}', transform=ax.transAxes,
                            ha='center', va='top', fontsize=10)
                except ValueError:
                    pass

    plt.tight_layout()
    # Save both SVG (editable) and PNG (viewable)
    fig.savefig(output_path, bbox_inches='tight')
    png_path = output_path.replace('.svg', '.png')
    fig.savefig(png_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved distribution plot: {output_path} and {png_path}")


def plot_its_tar1_ratio(cluster_its_tar1: pd.DataFrame, output_path: str,
                        dark_mode: bool = False):
    """Plot ITS:TAR1 ratio to show selective ITS depletion."""
    if cluster_its_tar1.empty:
        print("  No cluster data to plot")
        return

    plt.rcParams['svg.fonttype'] = 'none'

    if dark_mode:
        plt.style.use('dark_background')
        text_color = 'white'
        bbox_color = '#333333'
    else:
        text_color = 'black'
        bbox_color = 'white'

    # Calculate ITS:TAR1 ratio (ITS / (ITS + TAR1))
    df = cluster_its_tar1.copy()
    df['its_ratio'] = df['its_bp_per_kb'] / (df['its_bp_per_kb'] + df['tar1_bp_per_kb'] + 0.1)  # +0.1 to avoid div by 0

    # Filter to Normal and Tumor enriched only
    plot_data = df[df['enrichment'].isin(['Normal-enriched', 'Tumor-enriched'])]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    colors = {'Normal-enriched': '#3b82f6', 'Tumor-enriched': '#ef4444'}

    # Panel 1: Boxplot of ratio
    ax = axes[0]
    order = ['Normal-enriched', 'Tumor-enriched']
    data_by_group = [plot_data[plot_data['enrichment'] == e]['its_ratio'].values for e in order]

    bp = ax.boxplot(data_by_group, labels=order, patch_artist=True)
    for patch, enrichment in zip(bp['boxes'], order):
        patch.set_facecolor(colors.get(enrichment, '#999'))
        patch.set_alpha(0.7)

    ax.set_ylabel('ITS / (ITS + TAR1)')
    ax.set_title('ITS:TAR1 Ratio by Enrichment')

    # Add p-value
    normal_vals = plot_data[plot_data['enrichment'] == 'Normal-enriched']['its_ratio'].values
    tumor_vals = plot_data[plot_data['enrichment'] == 'Tumor-enriched']['its_ratio'].values
    if len(normal_vals) > 1 and len(tumor_vals) > 1:
        _, pval = mannwhitneyu(normal_vals, tumor_vals, alternative='two-sided')
        fc = normal_vals.mean() / tumor_vals.mean() if tumor_vals.mean() > 0 else np.nan
        ax.text(0.5, 0.95, f'p = {pval:.2e}\nFC = {fc:.1f}x', transform=ax.transAxes,
                ha='center', va='top', fontsize=10)

    # Panel 2: Scatter of ratio vs log2(OR)
    ax = axes[1]
    for enrichment in ['Normal-enriched', 'Tumor-enriched']:
        subset = df[df['enrichment'] == enrichment]
        if len(subset) > 0:
            ax.scatter(
                subset['log2_or'],
                subset['its_ratio'],
                s=subset['size'] / 5,
                c=colors.get(enrichment, '#999'),
                alpha=0.6,
                label=f"{enrichment} (n={len(subset)})"
            )
            for _, row in subset.iterrows():
                ax.annotate(
                    str(int(row['cluster_id'])),
                    (row['log2_or'], row['its_ratio']),
                    fontsize=6, alpha=0.7, color=text_color,
                    ha='center', va='bottom',
                    xytext=(0, 3), textcoords='offset points'
                )

    ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('log2(Odds Ratio)')
    ax.set_ylabel('ITS / (ITS + TAR1)')
    ax.set_title('ITS:TAR1 Ratio vs Cluster Enrichment')
    ax.legend(loc='upper right', fontsize=8)

    # Add correlation
    valid = df.dropna(subset=['log2_or', 'its_ratio'])
    if len(valid) > 2:
        rho, pval = spearmanr(valid['log2_or'], valid['its_ratio'])
        ax.text(0.02, 0.02, f'Spearman r = {rho:.3f}\np = {pval:.2e}',
                transform=ax.transAxes, va='bottom', fontsize=10,
                bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.8))

    plt.tight_layout()
    fig.savefig(output_path, bbox_inches='tight')
    png_path = output_path.replace('.svg', '.png')
    fig.savefig(png_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved ITS:TAR1 ratio plot: {output_path} and {png_path}")


def plot_its_vs_tar1_scatter(cluster_its_tar1: pd.DataFrame, output_path: str,
                             dark_mode: bool = False):
    """Scatter plot of ITS vs TAR1 to show decoupling in tumor cells."""
    if cluster_its_tar1.empty:
        print("  No cluster data to plot")
        return

    plt.rcParams['svg.fonttype'] = 'none'

    if dark_mode:
        plt.style.use('dark_background')
        text_color = 'white'
    else:
        text_color = 'black'

    fig, ax = plt.subplots(figsize=(8, 8))
    colors = {'Normal-enriched': '#3b82f6', 'Tumor-enriched': '#ef4444', 'mixed': '#888888'}

    for enrichment in ['Normal-enriched', 'Tumor-enriched', 'mixed']:
        subset = cluster_its_tar1[cluster_its_tar1['enrichment'] == enrichment]
        if len(subset) > 0:
            ax.scatter(
                subset['tar1_bp_per_kb'],
                subset['its_bp_per_kb'],
                s=subset['size'] / 3,
                c=colors.get(enrichment, '#999'),
                alpha=0.6,
                label=f"{enrichment} (n={len(subset)})"
            )
            for _, row in subset.iterrows():
                ax.annotate(
                    str(int(row['cluster_id'])),
                    (row['tar1_bp_per_kb'], row['its_bp_per_kb']),
                    fontsize=6, alpha=0.7, color=text_color,
                    ha='center', va='bottom',
                    xytext=(0, 3), textcoords='offset points'
                )

    ax.set_xlabel('TAR1 bp/kb read')
    ax.set_ylabel('ITS bp/kb read')
    ax.set_title('ITS vs TAR1: Selective ITS Depletion in Tumor Clusters')
    ax.legend(loc='upper right', fontsize=9)

    # Add annotation explaining the pattern
    ax.annotate('Normal clusters:\nHigh ITS, variable TAR1',
                xy=(0.05, 0.95), xycoords='axes fraction',
                fontsize=9, va='top', ha='left',
                bbox=dict(boxstyle='round', facecolor='#3b82f6', alpha=0.3))
    ax.annotate('Tumor clusters:\nLow ITS, TAR1 maintained',
                xy=(0.95, 0.05), xycoords='axes fraction',
                fontsize=9, va='bottom', ha='right',
                bbox=dict(boxstyle='round', facecolor='#ef4444', alpha=0.3))

    plt.tight_layout()
    fig.savefig(output_path, bbox_inches='tight')
    png_path = output_path.replace('.svg', '.png')
    fig.savefig(png_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved ITS vs TAR1 scatter: {output_path} and {png_path}")


def plot_waterfall(cluster_its_tar1: pd.DataFrame, output_path: str,
                   dark_mode: bool = False):
    """Waterfall plot showing ITS content ranked by cluster enrichment."""
    if cluster_its_tar1.empty:
        print("  No cluster data to plot")
        return

    plt.rcParams['svg.fonttype'] = 'none'

    if dark_mode:
        plt.style.use('dark_background')
        text_color = 'white'
    else:
        text_color = 'black'

    # Sort by log2(OR) - most Normal-enriched on left
    df = cluster_its_tar1.sort_values('log2_or', ascending=False).reset_index(drop=True)

    fig, axes = plt.subplots(2, 1, figsize=(14, 8), gridspec_kw={'height_ratios': [3, 1]})

    colors = {'Normal-enriched': '#3b82f6', 'Tumor-enriched': '#ef4444', 'mixed': '#888888'}
    bar_colors = [colors.get(e, '#999') for e in df['enrichment']]

    # Top panel: ITS content
    ax = axes[0]
    bars = ax.bar(range(len(df)), df['its_bp_per_kb'], color=bar_colors, alpha=0.8, edgecolor='none')

    # Add cluster labels on top of bars
    for i, (_, row) in enumerate(df.iterrows()):
        ax.text(i, row['its_bp_per_kb'] + 1, str(int(row['cluster_id'])),
                ha='center', va='bottom', fontsize=5, rotation=90, color=text_color)

    ax.set_ylabel('ITS bp/kb read')
    ax.set_title('ITS Content by Cluster (ranked by Normal->Tumor enrichment)')
    ax.set_xlim(-0.5, len(df) - 0.5)

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=colors['Normal-enriched'], label='Normal-enriched'),
                       Patch(facecolor=colors['mixed'], label='Mixed'),
                       Patch(facecolor=colors['Tumor-enriched'], label='Tumor-enriched')]
    ax.legend(handles=legend_elements, loc='upper right')

    # Bottom panel: log2(OR) to show the enrichment gradient
    ax = axes[1]
    ax.bar(range(len(df)), df['log2_or'], color=bar_colors, alpha=0.8, edgecolor='none')
    ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
    ax.set_ylabel('log2(OR)')
    ax.set_xlabel('Clusters (ranked by enrichment)')
    ax.set_xlim(-0.5, len(df) - 0.5)

    plt.tight_layout()
    fig.savefig(output_path, bbox_inches='tight')
    png_path = output_path.replace('.svg', '.png')
    fig.savefig(png_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved waterfall plot: {output_path} and {png_path}")


def plot_combined_figure(cluster_its_tar1: pd.DataFrame, output_path: str,
                         dark_mode: bool = False, solitary_tar1_df: pd.DataFrame = None,
                         tar1_context_df: pd.DataFrame = None, pairwise_df: pd.DataFrame = None,
                         read_level_features: pd.DataFrame = None, cluster_analysis: pd.DataFrame = None):
    """Publication-ready combined figure with multiple panels."""
    if cluster_its_tar1.empty:
        print("  No cluster data to plot")
        return

    plt.rcParams['svg.fonttype'] = 'none'

    if dark_mode:
        plt.style.use('dark_background')
        text_color = 'white'
        bbox_color = '#333333'
    else:
        text_color = 'black'
        bbox_color = 'white'

    # Calculate ITS:TAR1 ratio
    df = cluster_its_tar1.copy()
    df['its_ratio'] = df['its_bp_per_kb'] / (df['its_bp_per_kb'] + df['tar1_bp_per_kb'] + 0.1)

    # Determine grid size based on available data
    has_solitary = solitary_tar1_df is not None and not solitary_tar1_df.empty
    has_context = tar1_context_df is not None and not tar1_context_df.empty
    has_pairwise = pairwise_df is not None and not pairwise_df.empty

    if has_solitary and has_context and has_pairwise:
        # 3x4 grid for all panels including pairwise
        fig = plt.figure(figsize=(22, 16))
        gs = fig.add_gridspec(3, 4, hspace=0.35, wspace=0.3)
    elif has_solitary and has_context:
        # 3x3 grid for all panels
        fig = plt.figure(figsize=(18, 16))
        gs = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.3)
    elif has_solitary:
        fig = plt.figure(figsize=(20, 12))
        gs = fig.add_gridspec(2, 4, hspace=0.3, wspace=0.3)
    else:
        fig = plt.figure(figsize=(16, 12))
        gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)

    colors = {'Normal-enriched': '#3b82f6', 'Tumor-enriched': '#ef4444', 'mixed': '#888888'}
    plot_data = df[df['enrichment'].isin(['Normal-enriched', 'Tumor-enriched'])]

    # Panel A: ITS vs enrichment scatter
    ax = fig.add_subplot(gs[0, 0])
    for enrichment in ['Normal-enriched', 'Tumor-enriched', 'mixed']:
        subset = df[df['enrichment'] == enrichment]
        if len(subset) > 0:
            ax.scatter(subset['log2_or'], subset['its_bp_per_kb'],
                      s=subset['size'] / 5, c=colors.get(enrichment), alpha=0.6,
                      label=f"{enrichment} (n={len(subset)})")

    ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('log2(Odds Ratio)')
    ax.set_ylabel('ITS bp/kb read')
    ax.set_title('A. ITS vs Cluster Enrichment')
    ax.legend(loc='upper right', fontsize=7)

    valid = df.dropna(subset=['log2_or', 'its_bp_per_kb'])
    if len(valid) > 2:
        rho, pval = spearmanr(valid['log2_or'], valid['its_bp_per_kb'])
        ax.text(0.02, 0.98, f'r = {rho:.2f}, p = {pval:.1e}',
                transform=ax.transAxes, va='top', fontsize=9,
                bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.8))

    # Panel B: TAR1 vs enrichment scatter
    ax = fig.add_subplot(gs[0, 1])
    for enrichment in ['Normal-enriched', 'Tumor-enriched', 'mixed']:
        subset = df[df['enrichment'] == enrichment]
        if len(subset) > 0:
            ax.scatter(subset['log2_or'], subset['tar1_bp_per_kb'],
                      s=subset['size'] / 5, c=colors.get(enrichment), alpha=0.6)

    ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('log2(Odds Ratio)')
    ax.set_ylabel('TAR1 bp/kb read')
    ax.set_title('B. TAR1 vs Cluster Enrichment')

    valid = df.dropna(subset=['log2_or', 'tar1_bp_per_kb'])
    if len(valid) > 2:
        rho, pval = spearmanr(valid['log2_or'], valid['tar1_bp_per_kb'])
        ax.text(0.02, 0.98, f'r = {rho:.2f}, p = {pval:.1e}',
                transform=ax.transAxes, va='top', fontsize=9,
                bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.8))

    # Panel C: ITS:TAR1 ratio boxplot
    ax = fig.add_subplot(gs[0, 2])
    order = ['Normal-enriched', 'Tumor-enriched']
    data_by_group = [plot_data[plot_data['enrichment'] == e]['its_ratio'].values for e in order]

    bp = ax.boxplot(data_by_group, labels=['Normal', 'Tumor'], patch_artist=True)
    for patch, enrichment in zip(bp['boxes'], order):
        patch.set_facecolor(colors.get(enrichment))
        patch.set_alpha(0.7)

    ax.set_ylabel('ITS / (ITS + TAR1)')
    ax.set_title('C. ITS:TAR1 Ratio')

    normal_vals = plot_data[plot_data['enrichment'] == 'Normal-enriched']['its_ratio'].values
    tumor_vals = plot_data[plot_data['enrichment'] == 'Tumor-enriched']['its_ratio'].values
    if len(normal_vals) > 1 and len(tumor_vals) > 1:
        _, pval = mannwhitneyu(normal_vals, tumor_vals)
        fc = normal_vals.mean() / tumor_vals.mean() if tumor_vals.mean() > 0 else np.nan
        ax.text(0.5, 0.95, f'p = {pval:.1e}\n{fc:.1f}x higher\nin Normal',
                transform=ax.transAxes, ha='center', va='top', fontsize=9)

    # Panel D: ITS vs TAR1 scatter (decoupling)
    ax = fig.add_subplot(gs[1, 0])
    for enrichment in ['Normal-enriched', 'Tumor-enriched', 'mixed']:
        subset = df[df['enrichment'] == enrichment]
        if len(subset) > 0:
            ax.scatter(subset['tar1_bp_per_kb'], subset['its_bp_per_kb'],
                      s=subset['size'] / 3, c=colors.get(enrichment), alpha=0.6,
                      label=f"{enrichment}")

    ax.set_xlabel('TAR1 bp/kb read')
    ax.set_ylabel('ITS bp/kb read')
    ax.set_title('D. ITS vs TAR1 (Selective ITS Loss)')
    ax.legend(loc='upper right', fontsize=8)

    # Panel E: Boxplots of ITS and TAR1
    ax = fig.add_subplot(gs[1, 1])
    positions = [1, 2, 4, 5]
    labels = ['Normal\nITS', 'Tumor\nITS', 'Normal\nTAR1', 'Tumor\nTAR1']

    normal_its = plot_data[plot_data['enrichment'] == 'Normal-enriched']['its_bp_per_kb'].values
    tumor_its = plot_data[plot_data['enrichment'] == 'Tumor-enriched']['its_bp_per_kb'].values
    normal_tar1 = plot_data[plot_data['enrichment'] == 'Normal-enriched']['tar1_bp_per_kb'].values
    tumor_tar1 = plot_data[plot_data['enrichment'] == 'Tumor-enriched']['tar1_bp_per_kb'].values

    bp = ax.boxplot([normal_its, tumor_its, normal_tar1, tumor_tar1],
                    positions=positions, patch_artist=True)
    box_colors = [colors['Normal-enriched'], colors['Tumor-enriched'],
                  colors['Normal-enriched'], colors['Tumor-enriched']]
    for patch, color in zip(bp['boxes'], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_xticks(positions)
    ax.set_xticklabels(labels)
    ax.set_ylabel('bp/kb read')
    ax.set_title('E. ITS vs TAR1 Comparison')

    # Add significance annotations
    _, p_its = mannwhitneyu(normal_its, tumor_its)
    _, p_tar1 = mannwhitneyu(normal_tar1, tumor_tar1)
    ax.text(1.5, ax.get_ylim()[1] * 0.95, f'p={p_its:.1e}', ha='center', fontsize=8)
    ax.text(4.5, ax.get_ylim()[1] * 0.95, f'p={p_tar1:.1e}', ha='center', fontsize=8)

    # Panel F: Solitary TAR1 bar chart (if data available)
    solitary_stats = None
    context_stats = None
    from scipy.stats import chi2_contingency

    if has_solitary:
        if has_context:
            ax = fig.add_subplot(gs[2, 0])  # 3x3 grid: row 2, col 0
        else:
            ax = fig.add_subplot(gs[1, 2])  # 2x4 grid

        normal_sol_df = solitary_tar1_df[solitary_tar1_df['enrichment'] == 'Normal-enriched']
        tumor_sol_df = solitary_tar1_df[solitary_tar1_df['enrichment'] == 'Tumor-enriched']

        normal_solitary = normal_sol_df['solitary_tar1_elements'].sum()
        normal_total = normal_sol_df['total_tar1_elements'].sum()
        normal_paired = normal_total - normal_solitary

        tumor_solitary = tumor_sol_df['solitary_tar1_elements'].sum()
        tumor_total = tumor_sol_df['total_tar1_elements'].sum()
        tumor_paired = tumor_total - tumor_solitary

        normal_pct = normal_solitary / normal_total * 100 if normal_total > 0 else 0
        tumor_pct = tumor_solitary / tumor_total * 100 if tumor_total > 0 else 0

        # Chi-square test
        contingency = [[normal_solitary, normal_paired], [tumor_solitary, tumor_paired]]
        chi2, pval_chi2, _, _ = chi2_contingency(contingency)

        solitary_stats = {
            'normal_pct': normal_pct, 'tumor_pct': tumor_pct,
            'normal_n': normal_total, 'tumor_n': tumor_total,
            'chi2': chi2, 'pval': pval_chi2
        }

        # Bar chart
        x = [0, 1]
        heights = [normal_pct, tumor_pct]
        bar_colors = [colors['Normal-enriched'], colors['Tumor-enriched']]

        ax.bar(x, heights, color=bar_colors, alpha=0.7, width=0.6)
        ax.set_xticks(x)
        ax.set_xticklabels(['Normal', 'Tumor'])
        ax.set_ylabel('% Solitary TAR1')
        ax.set_title('F. Solitary TAR1 (no adjacent ITS)')
        ax.set_ylim(0, 100)
        ax.text(0, heights[0] + 2, f'n={normal_total}', ha='center', fontsize=8)
        ax.text(1, heights[1] + 2, f'n={tumor_total}', ha='center', fontsize=8)
        ax.text(0.5, 0.95, f'X2 = {chi2:.1f}\np = {pval_chi2:.1e}',
                transform=ax.transAxes, ha='center', va='top', fontsize=9,
                bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.8))

    # Panel G: TAR1 Context - Telomeric vs ITS-associated (if data available)
    if has_context:
        ax = fig.add_subplot(gs[2, 1])  # 3x3 grid: row 2, col 1

        normal_ctx = tar1_context_df[tar1_context_df['enrichment'] == 'Normal-enriched']
        tumor_ctx = tar1_context_df[tar1_context_df['enrichment'] == 'Tumor-enriched']

        # Calculate percentages
        normal_telo = normal_ctx['telomeric_tar1'].sum()
        normal_its = normal_ctx['its_associated_tar1'].sum()
        normal_total_ctx = normal_ctx['total_tar1'].sum()

        tumor_telo = tumor_ctx['telomeric_tar1'].sum()
        tumor_its = tumor_ctx['its_associated_tar1'].sum()
        tumor_total_ctx = tumor_ctx['total_tar1'].sum()

        normal_telo_pct = normal_telo / normal_total_ctx * 100 if normal_total_ctx > 0 else 0
        normal_its_pct = normal_its / normal_total_ctx * 100 if normal_total_ctx > 0 else 0
        tumor_telo_pct = tumor_telo / tumor_total_ctx * 100 if tumor_total_ctx > 0 else 0
        tumor_its_pct = tumor_its / tumor_total_ctx * 100 if tumor_total_ctx > 0 else 0

        # Chi-square for telomeric
        cont_telo = [[normal_telo, normal_total_ctx - normal_telo],
                     [tumor_telo, tumor_total_ctx - tumor_telo]]
        chi2_telo, pval_telo, _, _ = chi2_contingency(cont_telo)

        # Chi-square for ITS-associated
        cont_its = [[normal_its, normal_total_ctx - normal_its],
                    [tumor_its, tumor_total_ctx - tumor_its]]
        chi2_its_ctx, pval_its_ctx, _, _ = chi2_contingency(cont_its)

        context_stats = {
            'normal_telo_pct': normal_telo_pct, 'tumor_telo_pct': tumor_telo_pct,
            'normal_its_pct': normal_its_pct, 'tumor_its_pct': tumor_its_pct,
            'chi2_telo': chi2_telo, 'pval_telo': pval_telo,
            'chi2_its': chi2_its_ctx, 'pval_its': pval_its_ctx
        }

        # Grouped bar chart
        x = np.array([0, 1])
        width = 0.35

        telo_heights = [normal_telo_pct, tumor_telo_pct]
        its_heights = [normal_its_pct, tumor_its_pct]

        ax.bar(x - width/2, telo_heights, width, label='Telomeric', color='#22c55e', alpha=0.7)
        ax.bar(x + width/2, its_heights, width, label='ITS-associated', color='#f59e0b', alpha=0.7)

        ax.set_xticks(x)
        ax.set_xticklabels(['Normal', 'Tumor'])
        ax.set_ylabel('% of TAR1 elements')
        ax.set_title('G. TAR1 Context')
        ax.legend(loc='upper right', fontsize=8)
        ax.set_ylim(0, 100)

        # Add p-values
        ax.text(0.5, 0.95, f'Telomeric: p={pval_telo:.1e}\nITS-assoc: p={pval_its_ctx:.1e}',
                transform=ax.transAxes, ha='center', va='top', fontsize=8,
                bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.8))

    # Panel H: Pairwise Feature Correlations (if data available)
    pairwise_stats = None
    if has_pairwise and read_level_features is not None and cluster_analysis is not None:
        ax = fig.add_subplot(gs[2, 2])  # 3x4 grid: row 2, col 2

        # Add enrichment to read features
        cluster_enrichment = cluster_analysis.set_index('cluster_id')['enrichment'].to_dict()
        read_df = read_level_features.copy()
        read_df['enrichment'] = read_df['cluster'].map(cluster_enrichment)

        # Filter to Normal and Tumor
        read_df = read_df[read_df['enrichment'].isin(['Normal-enriched', 'Tumor-enriched'])]

        # Sample for plotting
        sample_size = min(3000, len(read_df))
        plot_read_df = read_df.sample(n=sample_size, random_state=42) if len(read_df) > sample_size else read_df

        # Plot ITS vs TAR1 at read level
        for enrichment in ['Normal-enriched', 'Tumor-enriched']:
            subset = plot_read_df[plot_read_df['enrichment'] == enrichment]
            label = 'Normal' if 'Normal' in enrichment else 'Tumor'
            ax.scatter(subset['its_bp_per_kb'], subset['tar1_bp_per_kb'],
                      c=colors[enrichment], alpha=0.3, s=8, label=label)

        ax.set_xlabel('ITS bp/kb')
        ax.set_ylabel('TAR1 bp/kb')
        ax.set_title('H. Read-Level ITS vs TAR1')
        ax.legend(loc='upper right', fontsize=8, markerscale=2)

        # Get correlations from pairwise_df
        normal_corr = pairwise_df[(pairwise_df['feature_1'] == 'ITS') &
                                   (pairwise_df['feature_2'] == 'TAR1') &
                                   (pairwise_df['enrichment'] == 'Normal-enriched')]['spearman_rho'].values
        tumor_corr = pairwise_df[(pairwise_df['feature_1'] == 'ITS') &
                                  (pairwise_df['feature_2'] == 'TAR1') &
                                  (pairwise_df['enrichment'] == 'Tumor-enriched')]['spearman_rho'].values

        if len(normal_corr) > 0 and len(tumor_corr) > 0:
            pairwise_stats = {'normal_rho': normal_corr[0], 'tumor_rho': tumor_corr[0]}
            ax.text(0.02, 0.98, f'N: r={normal_corr[0]:.2f}\nT: r={tumor_corr[0]:.2f}',
                    transform=ax.transAxes, va='top', fontsize=9,
                    bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.8))

        # Panel I: Correlation heatmap comparison
        ax2 = fig.add_subplot(gs[2, 3])  # 3x4 grid: row 2, col 3

        # Extract correlations for heatmap
        features = ['ITS', 'TAR1', 'Canonical', 'Non-canonical']
        normal_matrix = np.zeros((4, 4))
        tumor_matrix = np.zeros((4, 4))

        for i, f1 in enumerate(features):
            normal_matrix[i, i] = 1.0
            tumor_matrix[i, i] = 1.0
            for j, f2 in enumerate(features):
                if i < j:
                    row = pairwise_df[(pairwise_df['feature_1'] == f1) &
                                       (pairwise_df['feature_2'] == f2) &
                                       (pairwise_df['enrichment'] == 'Normal-enriched')]
                    if len(row) > 0:
                        normal_matrix[i, j] = row['spearman_rho'].values[0]
                        normal_matrix[j, i] = row['spearman_rho'].values[0]

                    row = pairwise_df[(pairwise_df['feature_1'] == f1) &
                                       (pairwise_df['feature_2'] == f2) &
                                       (pairwise_df['enrichment'] == 'Tumor-enriched')]
                    if len(row) > 0:
                        tumor_matrix[i, j] = row['spearman_rho'].values[0]
                        tumor_matrix[j, i] = row['spearman_rho'].values[0]

        # Show difference: Tumor - Normal
        diff_matrix = tumor_matrix - normal_matrix
        im = ax2.imshow(diff_matrix, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')

        ax2.set_xticks(range(4))
        ax2.set_yticks(range(4))
        short_labels = ['ITS', 'TAR1', 'Can.', 'Non-can.']
        ax2.set_xticklabels(short_labels, fontsize=8, rotation=45, ha='right')
        ax2.set_yticklabels(short_labels, fontsize=8)
        ax2.set_title('I. Correlation d (T-N)')

        # Add values to cells
        for i in range(4):
            for j in range(4):
                if i != j:
                    val = diff_matrix[i, j]
                    color = 'white' if abs(val) > 0.4 else text_color
                    ax2.text(j, i, f'{val:.2f}', ha='center', va='center', fontsize=7, color=color)

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax2, shrink=0.8)
        cbar.set_label('dr', fontsize=8)

    # Summary panel position depends on layout
    if has_solitary and has_context and has_pairwise:
        ax = fig.add_subplot(gs[0, 3])  # 3x4 grid: row 0, col 3
        panel_label = 'J'
    elif has_solitary and has_context:
        ax = fig.add_subplot(gs[2, 2])  # 3x3 grid: row 2, col 2
        panel_label = 'H'
    elif has_solitary:
        ax = fig.add_subplot(gs[1, 3])  # 2x4 grid
        panel_label = 'G'
    else:
        ax = fig.add_subplot(gs[1, 2])  # 2x3 grid
        panel_label = 'F'

    ax.axis('off')

    # Create summary text based on available data
    if solitary_stats and context_stats and pairwise_stats:
        summary_text = """Summary Statistics
═══════════════════════════
ITS (bp/kb): N={:.1f}, T={:.1f}
  p = {:.1e}
TAR1 (bp/kb): N={:.1f}, T={:.1f}
  p = {:.1e}

Solitary TAR1:
  N={:.1f}%, T={:.1f}%
  p = {:.1e}

ITS-TAR1 correlation:
  N: r={:.2f}, T: r={:.2f}

═══════════════════════════
BIR causes selective loss of
ITS. ITS-TAR1 decoupling is
stronger in Normal reads.
        """.format(
            normal_its.mean(), tumor_its.mean(), p_its,
            normal_tar1.mean(), tumor_tar1.mean(), p_tar1,
            solitary_stats['normal_pct'], solitary_stats['tumor_pct'], solitary_stats['pval'],
            pairwise_stats['normal_rho'], pairwise_stats['tumor_rho']
        )
    elif solitary_stats and context_stats:
        summary_text = """Summary Statistics
═════════════════════════════
ITS (bp/kb): N={:.1f}, T={:.1f}
  p = {:.1e}
TAR1 (bp/kb): N={:.1f}, T={:.1f}
  p = {:.1e}

Solitary TAR1:
  N={:.1f}%, T={:.1f}%
  p = {:.1e}

TAR1 Context:
  Telomeric: N={:.1f}%, T={:.1f}%
    p = {:.1e}
  ITS-assoc: N={:.1f}%, T={:.1f}%
    p = {:.1e}

═════════════════════════════
BIR causes selective loss of
interstitial TAR1-ITS pairs.
Remaining TAR1 is telomeric.
        """.format(
            normal_its.mean(), tumor_its.mean(), p_its,
            normal_tar1.mean(), tumor_tar1.mean(), p_tar1,
            solitary_stats['normal_pct'], solitary_stats['tumor_pct'], solitary_stats['pval'],
            context_stats['normal_telo_pct'], context_stats['tumor_telo_pct'], context_stats['pval_telo'],
            context_stats['normal_its_pct'], context_stats['tumor_its_pct'], context_stats['pval_its']
        )
    elif solitary_stats:
        summary_text = """Summary Statistics
═════════════════════════════
ITS (bp/kb): N={:.1f}, T={:.1f}
  p = {:.1e}
TAR1 (bp/kb): N={:.1f}, T={:.1f}
  p = {:.1e}

Solitary TAR1:
  Normal: {:.1f}%
  Tumor:  {:.1f}%
  p = {:.1e}

═════════════════════════════
ITS is selectively depleted
in tumor clusters. TAR1
without adjacent ITS is more
common in tumors, consistent
with BIR-mediated ITS loss.
        """.format(
            normal_its.mean(), tumor_its.mean(), p_its,
            normal_tar1.mean(), tumor_tar1.mean(), p_tar1,
            solitary_stats['normal_pct'], solitary_stats['tumor_pct'], solitary_stats['pval']
        )
    else:
        summary_text = """Summary Statistics
═════════════════════════════
ITS (bp/kb): N={:.1f}, T={:.1f}
  p = {:.1e}
TAR1 (bp/kb): N={:.1f}, T={:.1f}
  p = {:.1e}
ITS:TAR1 Ratio:
  Normal: {:.2f}
  Tumor:  {:.2f}
  p = {:.1e}

═════════════════════════════
ITS is selectively depleted
in tumor clusters while TAR1
is maintained, consistent
with BIR-mediated ITS loss.
        """.format(
            normal_its.mean(), tumor_its.mean(), p_its,
            normal_tar1.mean(), tumor_tar1.mean(), p_tar1,
            normal_vals.mean(), tumor_vals.mean(),
            mannwhitneyu(normal_vals, tumor_vals)[1]
        )

    ax.set_title(f'{panel_label}. Summary')
    ax.text(0.02, 0.95, summary_text, transform=ax.transAxes, fontsize=8,
            va='top', ha='left', family='monospace',
            bbox=dict(boxstyle='round', facecolor=bbox_color, alpha=0.5))

    plt.savefig(output_path, bbox_inches='tight')
    png_path = output_path.replace('.svg', '.png')
    plt.savefig(png_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved combined figure: {output_path} and {png_path}")


# --- Main function ---
def main():
    args = parse_args()

    # Validate arguments
    cluster_mode = args.cluster_prefix is not None
    if not cluster_mode:
        if not args.samples or not args.sample_metadata:
            print("ERROR: Either --cluster-analysis-prefix OR both --samples and --sample-metadata are required")
            sys.exit(1)

    # Set up logging
    if args.log_file:
        log_path = f"{args.output_prefix}.log"
        sys.stdout = TeeLogger(log_path)

    print("=" * 60)
    print("KaryoScope ITS Analysis")
    print("=" * 60)
    print(f"\nCommand: {_original_command}")
    print(f"Output prefix: {args.output_prefix}")
    print(f"Mode: {'CLUSTER' if cluster_mode else 'SAMPLE'}")

    # =========================================================================
    # CLUSTER MODE - Analyze only clustered reads
    # =========================================================================
    # Load enrichment mapping if provided
    enrichment_mapping = None
    if args.enrichment_mapping:
        import json
        with open(args.enrichment_mapping) as f:
            enrichment_mapping = json.load(f)
        print(f"\nLoaded enrichment mapping from: {args.enrichment_mapping}")
        for k, v in enrichment_mapping.items():
            print(f"  {k} -> {v}")

    if cluster_mode:
        print("\n--- Loading cluster analysis data ---")
        read_assignments = load_read_assignments(args.cluster_prefix)
        cluster_analysis = load_cluster_analysis(args.cluster_prefix, enrichment_mapping=enrichment_mapping)

        # Build sample to group mapping from read assignments
        sample_to_group = dict(zip(read_assignments['sample'], read_assignments['group']))

        print("\n--- Loading ITS/TAR1 features for clustered reads ---")
        read_features = load_features_for_clustered_reads(
            read_assignments,
            args.bed_prefix,
            args.database,
            args.smoothness,
            args.min_feature_length
        )

        # === Analysis: Cluster-level ITS/TAR1 ===
        print("\n" + "=" * 60)
        print("Cluster-Level ITS/TAR1 Analysis")
        print("=" * 60)

        cluster_its_tar1 = calculate_cluster_its_tar1(
            read_assignments, cluster_analysis, read_features
        )

        # Save cluster ITS/TAR1 stats
        cluster_path = f"{args.output_prefix}.cluster_its_tar1.tsv"
        cluster_its_tar1.to_csv(cluster_path, sep='\t', index=False)
        print(f"  Saved: {cluster_path}")

        # Summary by enrichment direction
        print("\n  Summary by enrichment direction:")
        for enrichment in ['Normal-enriched', 'Tumor-enriched', 'mixed']:
            subset = cluster_its_tar1[cluster_its_tar1['enrichment'] == enrichment]
            if len(subset) > 0:
                print(f"\n  {enrichment} (n={len(subset)} clusters):")
                print(f"    Mean ITS bp/read: {subset['its_bp_mean'].mean():.1f}")
                print(f"    Mean TAR1 bp/read: {subset['tar1_bp_mean'].mean():.1f}")
                print(f"    Mean % reads with ITS: {subset['pct_reads_with_its'].mean():.1f}%")
                print(f"    Mean % reads with TAR1: {subset['pct_reads_with_tar1'].mean():.1f}%")

        # === Statistical comparison ===
        print("\n" + "=" * 60)
        print("Enrichment Direction Comparison")
        print("=" * 60)

        comparison_df = compare_by_enrichment(cluster_its_tar1)
        comparison_path = f"{args.output_prefix}.enrichment_comparison.tsv"
        comparison_df.to_csv(comparison_path, sep='\t', index=False)
        print(f"  Saved: {comparison_path}")

        print("\n  Normal-enriched vs Tumor-enriched clusters:")
        for _, row in comparison_df.iterrows():
            if pd.notna(row.get('p_value')):
                metric = row['metric']
                normal_mean = row.get('Normal-enriched_mean', 'N/A')
                tumor_mean = row.get('Tumor-enriched_mean', 'N/A')
                pval = row['p_value']
                fc = row.get('fold_change', np.nan)
                print(f"    {metric}: Normal={normal_mean:.2f}, Tumor={tumor_mean:.2f}, "
                      f"FC={fc:.2f}, p={pval:.2e}")

        # === Correlation analysis ===
        print("\n" + "=" * 60)
        print("Correlation Analysis")
        print("=" * 60)

        valid = cluster_its_tar1.dropna(subset=['log2_or', 'its_bp_mean'])
        if len(valid) > 2:
            rho_its, pval_its = spearmanr(valid['log2_or'], valid['its_bp_mean'])
            print(f"  ITS bp/read vs log2(OR): Spearman r = {rho_its:.3f}, p = {pval_its:.2e}")

        valid = cluster_its_tar1.dropna(subset=['log2_or', 'tar1_bp_mean'])
        if len(valid) > 2:
            rho_tar1, pval_tar1 = spearmanr(valid['log2_or'], valid['tar1_bp_mean'])
            print(f"  TAR1 bp/read vs log2(OR): Spearman r = {rho_tar1:.3f}, p = {pval_tar1:.2e}")

        # === Read length control analysis ===
        print("\n" + "=" * 60)
        print("Read Length Control Analysis")
        print("=" * 60)

        # Check read length by enrichment
        print("\n  Mean read length by enrichment direction:")
        for enrichment in ['Normal-enriched', 'Tumor-enriched', 'mixed']:
            subset = cluster_its_tar1[cluster_its_tar1['enrichment'] == enrichment]
            if len(subset) > 0:
                print(f"    {enrichment}: {subset['mean_read_length'].mean():.0f} bp")

        # Check if ITS bp/kb correlates with read length
        valid = cluster_its_tar1.dropna(subset=['mean_read_length', 'its_bp_per_kb'])
        if len(valid) > 2:
            rho_len, pval_len = spearmanr(valid['mean_read_length'], valid['its_bp_per_kb'])
            print(f"\n  ITS bp/kb vs read length: r = {rho_len:.3f}, p = {pval_len:.2e}")
            if pval_len > 0.05:
                print("    (Not significant - per-kb normalization removes read length bias)")

        # Partial correlation controlling for read length
        valid = cluster_its_tar1.dropna(subset=['log2_or', 'its_bp_per_kb', 'mean_read_length'])
        if len(valid) > 5:
            # Residualize ITS bp/kb against read length
            slope, intercept, _, _, _ = linregress(valid['mean_read_length'], valid['its_bp_per_kb'])
            its_residual = valid['its_bp_per_kb'] - (slope * valid['mean_read_length'] + intercept)

            # Residualize log2_or against read length
            slope, intercept, _, _, _ = linregress(valid['mean_read_length'], valid['log2_or'])
            log2or_residual = valid['log2_or'] - (slope * valid['mean_read_length'] + intercept)

            # Partial correlation
            rho_partial, pval_partial = spearmanr(log2or_residual, its_residual)
            print(f"\n  Partial correlation (ITS vs enrichment, controlling for read length):")
            print(f"    ITS bp/kb vs log2(OR): r = {rho_partial:.3f}, p = {pval_partial:.2e}")

        # Stratified analysis
        median_length = cluster_its_tar1['mean_read_length'].median()
        short_clusters = cluster_its_tar1[cluster_its_tar1['mean_read_length'] < median_length]
        long_clusters = cluster_its_tar1[cluster_its_tar1['mean_read_length'] >= median_length]

        print(f"\n  Stratified analysis (median read length: {median_length:.0f} bp):")
        for name, subset in [("Short reads (<median)", short_clusters), ("Long reads (>=median)", long_clusters)]:
            normal = subset[subset['enrichment'] == 'Normal-enriched']['its_bp_per_kb']
            tumor = subset[subset['enrichment'] == 'Tumor-enriched']['its_bp_per_kb']
            if len(normal) > 0 and len(tumor) > 0:
                fc = normal.mean() / tumor.mean() if tumor.mean() > 0 else np.nan
                print(f"    {name}: Normal={normal.mean():.1f}, Tumor={tumor.mean():.1f} ITS bp/kb (FC={fc:.1f}x)")

        # === Solitary TAR1 Analysis ===
        print("\n" + "=" * 60)
        print("Solitary TAR1 Analysis (ITS-TAR1 Adjacency)")
        print("=" * 60)

        adjacency_data = load_adjacency_data_for_clustered_reads(
            read_assignments,
            args.bed_prefix,
            args.database,
            args.smoothness,
            args.min_feature_length,
            max_distance=2000  # TAR1 is "solitary" if no ITS within 2kb
        )

        solitary_tar1_df = calculate_cluster_solitary_tar1(
            read_assignments, cluster_analysis, adjacency_data
        )

        # Save solitary TAR1 stats
        solitary_path = f"{args.output_prefix}.solitary_tar1.tsv"
        solitary_tar1_df.to_csv(solitary_path, sep='\t', index=False)
        print(f"  Saved: {solitary_path}")

        # Summary
        print("\n  Solitary TAR1 by enrichment direction:")
        for enrichment in ['Normal-enriched', 'Tumor-enriched', 'mixed']:
            subset = solitary_tar1_df[solitary_tar1_df['enrichment'] == enrichment]
            if len(subset) > 0:
                print(f"\n  {enrichment} (n={len(subset)} clusters):")
                print(f"    Mean % solitary TAR1: {subset['pct_solitary_tar1'].mean():.1f}%")
                print(f"    Total TAR1 elements: {subset['total_tar1_elements'].sum()}")
                print(f"    Solitary TAR1 elements: {subset['solitary_tar1_elements'].sum()}")

        # Statistical comparison - cluster level
        normal_solitary = solitary_tar1_df[solitary_tar1_df['enrichment'] == 'Normal-enriched']['pct_solitary_tar1']
        tumor_solitary = solitary_tar1_df[solitary_tar1_df['enrichment'] == 'Tumor-enriched']['pct_solitary_tar1']
        if len(normal_solitary) > 1 and len(tumor_solitary) > 1:
            _, pval = mannwhitneyu(normal_solitary, tumor_solitary, alternative='two-sided')
            print(f"\n  Cluster-level comparison (Mann-Whitney U):")
            print(f"    Normal: {normal_solitary.mean():.1f}%, Tumor: {tumor_solitary.mean():.1f}%")
            print(f"    p-value: {pval:.2e}")

        # Element-level chi-square test (more powerful)
        from scipy.stats import chi2_contingency
        normal_df = solitary_tar1_df[solitary_tar1_df['enrichment'] == 'Normal-enriched']
        tumor_df = solitary_tar1_df[solitary_tar1_df['enrichment'] == 'Tumor-enriched']

        normal_solitary_count = normal_df['solitary_tar1_elements'].sum()
        normal_total = normal_df['total_tar1_elements'].sum()
        normal_paired = normal_total - normal_solitary_count

        tumor_solitary_count = tumor_df['solitary_tar1_elements'].sum()
        tumor_total = tumor_df['total_tar1_elements'].sum()
        tumor_paired = tumor_total - tumor_solitary_count

        # Contingency table: [[normal_solitary, normal_paired], [tumor_solitary, tumor_paired]]
        contingency = [[normal_solitary_count, normal_paired],
                       [tumor_solitary_count, tumor_paired]]
        chi2, pval_chi2, dof, expected = chi2_contingency(contingency)

        normal_pct = normal_solitary_count / normal_total * 100 if normal_total > 0 else 0
        tumor_pct = tumor_solitary_count / tumor_total * 100 if tumor_total > 0 else 0

        print(f"\n  Element-level comparison (Chi-square test):")
        print(f"    Normal: {normal_solitary_count}/{normal_total} = {normal_pct:.1f}% solitary TAR1")
        print(f"    Tumor:  {tumor_solitary_count}/{tumor_total} = {tumor_pct:.1f}% solitary TAR1")
        print(f"    Chi-square = {chi2:.2f}, p-value = {pval_chi2:.2e}")

        # === TAR1 Context Analysis (Telomeric vs ITS-associated) ===
        print("\n" + "=" * 60)
        print("TAR1 Context Analysis (Telomeric vs ITS-associated)")
        print("=" * 60)

        tar1_context_data = load_tar1_telomere_adjacency(
            read_assignments,
            args.bed_prefix,
            args.database,
            args.smoothness,
            args.min_feature_length,
            max_distance=2000
        )

        tar1_context_df = calculate_cluster_tar1_context(
            read_assignments, cluster_analysis, tar1_context_data
        )

        # Save TAR1 context stats
        context_path = f"{args.output_prefix}.tar1_context.tsv"
        tar1_context_df.to_csv(context_path, sep='\t', index=False)
        print(f"  Saved: {context_path}")

        # Summary by enrichment
        print("\n  TAR1 context by enrichment direction:")
        for enrichment in ['Normal-enriched', 'Tumor-enriched', 'mixed']:
            subset = tar1_context_df[tar1_context_df['enrichment'] == enrichment]
            if len(subset) > 0:
                print(f"\n  {enrichment} (n={len(subset)} clusters):")
                print(f"    Total TAR1 elements: {subset['total_tar1'].sum()}")
                print(f"    Telomeric TAR1: {subset['telomeric_tar1'].sum()} ({subset['pct_telomeric'].mean():.1f}%)")
                print(f"    ITS-associated TAR1: {subset['its_associated_tar1'].sum()} ({subset['pct_its_associated'].mean():.1f}%)")

        # Element-level chi-square for telomeric TAR1
        normal_ctx = tar1_context_df[tar1_context_df['enrichment'] == 'Normal-enriched']
        tumor_ctx = tar1_context_df[tar1_context_df['enrichment'] == 'Tumor-enriched']

        normal_telomeric = normal_ctx['telomeric_tar1'].sum()
        normal_non_telomeric = normal_ctx['total_tar1'].sum() - normal_telomeric
        tumor_telomeric = tumor_ctx['telomeric_tar1'].sum()
        tumor_non_telomeric = tumor_ctx['total_tar1'].sum() - tumor_telomeric

        contingency_telo = [[normal_telomeric, normal_non_telomeric],
                           [tumor_telomeric, tumor_non_telomeric]]
        chi2_telo, pval_telo, _, _ = chi2_contingency(contingency_telo)

        normal_telo_pct = normal_telomeric / normal_ctx['total_tar1'].sum() * 100 if normal_ctx['total_tar1'].sum() > 0 else 0
        tumor_telo_pct = tumor_telomeric / tumor_ctx['total_tar1'].sum() * 100 if tumor_ctx['total_tar1'].sum() > 0 else 0

        print(f"\n  Element-level comparison - Telomeric TAR1 (Chi-square):")
        print(f"    Normal: {normal_telomeric}/{normal_ctx['total_tar1'].sum()} = {normal_telo_pct:.1f}% telomeric")
        print(f"    Tumor:  {tumor_telomeric}/{tumor_ctx['total_tar1'].sum()} = {tumor_telo_pct:.1f}% telomeric")
        print(f"    Chi-square = {chi2_telo:.2f}, p-value = {pval_telo:.2e}")

        # Element-level chi-square for ITS-associated TAR1
        normal_its_assoc = normal_ctx['its_associated_tar1'].sum()
        normal_non_its = normal_ctx['total_tar1'].sum() - normal_its_assoc
        tumor_its_assoc = tumor_ctx['its_associated_tar1'].sum()
        tumor_non_its = tumor_ctx['total_tar1'].sum() - tumor_its_assoc

        contingency_its = [[normal_its_assoc, normal_non_its],
                          [tumor_its_assoc, tumor_non_its]]
        chi2_its, pval_its, _, _ = chi2_contingency(contingency_its)

        normal_its_pct = normal_its_assoc / normal_ctx['total_tar1'].sum() * 100 if normal_ctx['total_tar1'].sum() > 0 else 0
        tumor_its_pct = tumor_its_assoc / tumor_ctx['total_tar1'].sum() * 100 if tumor_ctx['total_tar1'].sum() > 0 else 0

        print(f"\n  Element-level comparison - ITS-associated TAR1 (Chi-square):")
        print(f"    Normal: {normal_its_assoc}/{normal_ctx['total_tar1'].sum()} = {normal_its_pct:.1f}% ITS-associated")
        print(f"    Tumor:  {tumor_its_assoc}/{tumor_ctx['total_tar1'].sum()} = {tumor_its_pct:.1f}% ITS-associated")
        print(f"    Chi-square = {chi2_its:.2f}, p-value = {pval_its:.2e}")

        # === Pairwise Feature Analysis (Read-Level) ===
        print("\n" + "=" * 60)
        print("Pairwise Feature Analysis (Read-Level)")
        print("=" * 60)

        read_level_features = load_read_level_features(
            read_assignments,
            args.bed_prefix,
            args.database,
            args.smoothness,
            args.min_feature_length
        )

        # Calculate pairwise correlations
        pairwise_path = f"{args.output_prefix}.pairwise_correlations.tsv"
        pairwise_df = analyze_pairwise_correlations(
            read_level_features,
            cluster_analysis,
            pairwise_path
        )
        print(f"  Saved: {pairwise_path}")

        # Print summary of correlations
        print("\n  Pairwise correlations (Spearman r) by enrichment:")
        for enrichment in ['Normal-enriched', 'Tumor-enriched']:
            subset = pairwise_df[pairwise_df['enrichment'] == enrichment]
            print(f"\n  {enrichment}:")
            for _, row in subset.iterrows():
                if pd.notna(row['spearman_rho']):
                    sig = "*" if row['p_value'] < 0.05 else ""
                    print(f"    {row['feature_1']} vs {row['feature_2']}: r={row['spearman_rho']:.3f}{sig}")

        # Generate pairwise scatter plots (read-level)
        plot_pairwise_features(
            read_level_features,
            cluster_analysis,
            f"{args.output_prefix}.pairwise_features.svg",
            args.dark_mode
        )

        # === Cluster-Level Pairwise Analysis ===
        print("\n" + "=" * 60)
        print("Cluster-Level Pairwise Feature Analysis")
        print("=" * 60)

        # Aggregate read-level features to cluster level
        cluster_all_features = aggregate_cluster_all_features(
            read_level_features,
            cluster_analysis
        )

        # Calculate cluster-level correlations
        cluster_pairwise_path = f"{args.output_prefix}.cluster_pairwise_correlations.tsv"
        cluster_pairwise_df = analyze_cluster_pairwise_correlations(
            cluster_all_features,
            cluster_pairwise_path
        )
        print(f"  Saved: {cluster_pairwise_path}")

        # Print summary
        print("\n  Cluster-level pairwise correlations (Spearman r):")
        for enrichment in ['Normal-enriched', 'Tumor-enriched']:
            subset = cluster_pairwise_df[cluster_pairwise_df['enrichment'] == enrichment]
            print(f"\n  {enrichment}:")
            for _, row in subset.iterrows():
                if pd.notna(row['spearman_rho']):
                    sig = "*" if row['p_value'] < 0.05 else ""
                    print(f"    {row['feature_1']} vs {row['feature_2']}: r={row['spearman_rho']:.3f}{sig} (n={row['n_clusters']})")

        # Generate cluster-level pairwise grid plot
        plot_cluster_pairwise_grid(
            cluster_all_features,
            f"{args.output_prefix}.cluster_pairwise_grid.svg",
            args.dark_mode
        )

        # === ITS Contiguity Analysis ===
        print("\n" + "=" * 60)
        print("ITS Contiguity/Fragmentation Analysis")
        print("=" * 60)

        its_elements = load_its_contiguity_data(
            read_assignments,
            args.bed_prefix,
            args.database,
            args.smoothness,
            args.min_feature_length
        )

        read_contiguity, cluster_contiguity = analyze_its_contiguity(
            its_elements,
            cluster_analysis
        )

        # Save contiguity data
        if not read_contiguity.empty:
            read_contiguity_path = f"{args.output_prefix}.its_contiguity_reads.tsv"
            read_contiguity.to_csv(read_contiguity_path, sep='\t', index=False)
            print(f"  Saved: {read_contiguity_path}")

            cluster_contiguity_path = f"{args.output_prefix}.its_contiguity_clusters.tsv"
            cluster_contiguity.to_csv(cluster_contiguity_path, sep='\t', index=False)
            print(f"  Saved: {cluster_contiguity_path}")

            # Summary statistics
            print("\n  ITS element size by enrichment direction (read-level):")
            for enrichment in ['Normal-enriched', 'Tumor-enriched']:
                subset = read_contiguity[read_contiguity['enrichment'] == enrichment]
                if len(subset) > 0:
                    print(f"\n  {enrichment} (n={len(subset)} reads with ITS):")
                    print(f"    Mean ITS element size: {subset['its_mean_size'].mean():.1f} bp")
                    print(f"    Mean max ITS element: {subset['its_max_size'].mean():.1f} bp")
                    print(f"    Mean ITS elements/read: {subset['its_count'].mean():.2f}")

            # Statistical comparison
            normal_sizes = read_contiguity[read_contiguity['enrichment'] == 'Normal-enriched']['its_mean_size']
            tumor_sizes = read_contiguity[read_contiguity['enrichment'] == 'Tumor-enriched']['its_mean_size']
            if len(normal_sizes) > 1 and len(tumor_sizes) > 1:
                _, pval = mannwhitneyu(normal_sizes, tumor_sizes)
                fc = normal_sizes.mean() / tumor_sizes.mean() if tumor_sizes.mean() > 0 else np.nan
                print(f"\n  Mean ITS element size comparison:")
                print(f"    Normal: {normal_sizes.mean():.1f} bp, Tumor: {tumor_sizes.mean():.1f} bp")
                print(f"    Fold change: {fc:.2f}x, p-value: {pval:.2e}")

            # Element-level analysis (all ITS elements)
            its_elements_enriched = its_elements.copy()
            cluster_enrichment = cluster_analysis.set_index('cluster_id')['enrichment'].to_dict()
            its_elements_enriched['enrichment'] = its_elements_enriched['cluster'].map(cluster_enrichment)

            normal_elements = its_elements_enriched[its_elements_enriched['enrichment'] == 'Normal-enriched']['its_size']
            tumor_elements = its_elements_enriched[its_elements_enriched['enrichment'] == 'Tumor-enriched']['its_size']
            if len(normal_elements) > 1 and len(tumor_elements) > 1:
                _, pval_elem = mannwhitneyu(normal_elements, tumor_elements)
                fc_elem = normal_elements.mean() / tumor_elements.mean() if tumor_elements.mean() > 0 else np.nan
                print(f"\n  Element-level ITS size comparison:")
                print(f"    Normal: {normal_elements.mean():.1f} bp (n={len(normal_elements)} elements)")
                print(f"    Tumor: {tumor_elements.mean():.1f} bp (n={len(tumor_elements)} elements)")
                print(f"    Fold change: {fc_elem:.2f}x, p-value: {pval_elem:.2e}")

            # Generate plot
            plot_its_contiguity(
                read_contiguity,
                cluster_contiguity,
                f"{args.output_prefix}.its_contiguity.svg",
                args.dark_mode
            )

        # === Feature Interleaving Analysis ===
        print("\n" + "=" * 60)
        print("Feature Interleaving Analysis (ITS-Telomere Mixing)")
        print("=" * 60)

        feature_elements = load_feature_sequence_data(
            read_assignments,
            args.bed_prefix,
            args.database,
            args.smoothness,
            args.min_feature_length
        )

        read_interleaving, cluster_interleaving = analyze_feature_interleaving(
            feature_elements,
            cluster_analysis
        )

        if not read_interleaving.empty:
            # Save interleaving data
            read_interleaving_path = f"{args.output_prefix}.interleaving_reads.tsv"
            read_interleaving.to_csv(read_interleaving_path, sep='\t', index=False)
            print(f"  Saved: {read_interleaving_path}")

            cluster_interleaving_path = f"{args.output_prefix}.interleaving_clusters.tsv"
            cluster_interleaving.to_csv(cluster_interleaving_path, sep='\t', index=False)
            print(f"  Saved: {cluster_interleaving_path}")

            # Summary - focus on reads that have both ITS and telomere
            reads_with_both = read_interleaving[read_interleaving['has_its_and_telo']]

            print("\n  Feature interleaving (reads with both ITS and telomere):")
            for enrichment in ['Normal-enriched', 'Tumor-enriched']:
                subset = reads_with_both[reads_with_both['enrichment'] == enrichment]
                if len(subset) > 0:
                    print(f"\n  {enrichment} (n={len(subset)} reads):")
                    print(f"    Mean ITS<->Telo transitions: {subset['n_its_telo_transitions'].mean():.2f}")
                    print(f"    Mean interleaving score: {subset['interleaving_score'].mean():.2f}")
                    has_mixing = (subset['n_its_telo_transitions'] > 0).sum()
                    print(f"    Reads with ITS-Telo mixing: {has_mixing}/{len(subset)} ({100*has_mixing/len(subset):.1f}%)")

            # Statistical comparison
            normal_trans = reads_with_both[reads_with_both['enrichment'] == 'Normal-enriched']['n_its_telo_transitions']
            tumor_trans = reads_with_both[reads_with_both['enrichment'] == 'Tumor-enriched']['n_its_telo_transitions']
            if len(normal_trans) > 1 and len(tumor_trans) > 1:
                _, pval = mannwhitneyu(normal_trans, tumor_trans)
                print(f"\n  ITS<->Telomere transitions comparison:")
                print(f"    Normal: {normal_trans.mean():.2f}, Tumor: {tumor_trans.mean():.2f}")
                print(f"    p-value: {pval:.2e}")

            # Chi-square for proportion with mixing
            normal_with_mixing = (reads_with_both[reads_with_both['enrichment'] == 'Normal-enriched']['n_its_telo_transitions'] > 0).sum()
            normal_total = len(reads_with_both[reads_with_both['enrichment'] == 'Normal-enriched'])
            tumor_with_mixing = (reads_with_both[reads_with_both['enrichment'] == 'Tumor-enriched']['n_its_telo_transitions'] > 0).sum()
            tumor_total = len(reads_with_both[reads_with_both['enrichment'] == 'Tumor-enriched'])

            if normal_total > 0 and tumor_total > 0:
                from scipy.stats import chi2_contingency
                contingency = [[normal_with_mixing, normal_total - normal_with_mixing],
                              [tumor_with_mixing, tumor_total - tumor_with_mixing]]
                chi2, pval_chi2, _, _ = chi2_contingency(contingency)
                print(f"\n  Proportion with ITS-Telo mixing (Chi-square):")
                print(f"    Normal: {normal_with_mixing}/{normal_total} ({100*normal_with_mixing/normal_total:.1f}%)")
                print(f"    Tumor: {tumor_with_mixing}/{tumor_total} ({100*tumor_with_mixing/tumor_total:.1f}%)")
                print(f"    Chi-square = {chi2:.2f}, p-value = {pval_chi2:.2e}")

            # Generate plot
            plot_feature_interleaving(
                read_interleaving,
                cluster_interleaving,
                f"{args.output_prefix}.interleaving.svg",
                args.dark_mode
            )

        # === Fragment Count Analysis (Smoothed vs Presmoothed) ===
        print("\n" + "=" * 60)
        print("Fragment Count Analysis (Smoothed vs Presmoothed)")
        print("=" * 60)

        # Load smoothed counts
        smoothed_counts = load_fragment_counts(
            read_assignments,
            args.bed_prefix,
            args.database,
            'smoothed',
            args.min_feature_length
        )

        # Load presmoothed counts
        presmoothed_counts = load_fragment_counts(
            read_assignments,
            args.bed_prefix,
            args.database,
            'presmoothed',
            args.min_feature_length
        )

        # Analyze and save
        fragment_stats = analyze_fragment_counts(
            smoothed_counts,
            presmoothed_counts,
            cluster_analysis,
            f"{args.output_prefix}.fragment_counts.tsv"
        )
        print(f"  Saved: {args.output_prefix}.fragment_counts.tsv")

        # Print summary
        cluster_enrichment = cluster_analysis.set_index('cluster_id')['enrichment'].to_dict()

        print("\n  Fragment counts per read (reads with feature present):")
        print("\n  SMOOTHED:")
        for feature in ['its', 'tar1', 'canonical', 'noncanonical']:
            smoothed_with_enrich = smoothed_counts.copy()
            smoothed_with_enrich['enrichment'] = smoothed_with_enrich['cluster'].map(cluster_enrichment)

            count_col = f'{feature}_count'
            normal = smoothed_with_enrich[(smoothed_with_enrich['enrichment'] == 'Normal-enriched') &
                                          (smoothed_with_enrich[count_col] > 0)][count_col]
            tumor = smoothed_with_enrich[(smoothed_with_enrich['enrichment'] == 'Tumor-enriched') &
                                         (smoothed_with_enrich[count_col] > 0)][count_col]
            if len(normal) > 1 and len(tumor) > 1:
                _, pval = mannwhitneyu(normal, tumor)
                print(f"    {feature.upper():12s}: N={normal.mean():.2f}, T={tumor.mean():.2f}, p={pval:.2e}")

        print("\n  PRESMOOTHED:")
        for feature in ['its', 'tar1', 'canonical', 'noncanonical']:
            presmooth_with_enrich = presmoothed_counts.copy()
            presmooth_with_enrich['enrichment'] = presmooth_with_enrich['cluster'].map(cluster_enrichment)

            count_col = f'{feature}_count'
            normal = presmooth_with_enrich[(presmooth_with_enrich['enrichment'] == 'Normal-enriched') &
                                           (presmooth_with_enrich[count_col] > 0)][count_col]
            tumor = presmooth_with_enrich[(presmooth_with_enrich['enrichment'] == 'Tumor-enriched') &
                                          (presmooth_with_enrich[count_col] > 0)][count_col]
            if len(normal) > 1 and len(tumor) > 1:
                _, pval = mannwhitneyu(normal, tumor)
                print(f"    {feature.upper():12s}: N={normal.mean():.2f}, T={tumor.mean():.2f}, p={pval:.2e}")

        # Generate plot
        plot_fragment_counts(
            smoothed_counts,
            presmoothed_counts,
            cluster_analysis,
            f"{args.output_prefix}.fragment_counts.svg",
            args.dark_mode
        )

        # === Normalized Fragmentation Analysis ===
        print("\n" + "=" * 60)
        print("Normalized Fragmentation Analysis (Fragments per kb of Feature)")
        print("=" * 60)

        # Analyze and save normalized fragmentation
        norm_frag_stats = analyze_normalized_fragmentation(
            smoothed_counts,
            presmoothed_counts,
            cluster_analysis,
            f"{args.output_prefix}.normalized_fragmentation.tsv"
        )
        print(f"  Saved: {args.output_prefix}.normalized_fragmentation.tsv")

        # Print summary
        print("\n  Fragments per kb of feature (reads with feature present):")
        print("\n  SMOOTHED:")
        for feature in ['its', 'tar1', 'canonical', 'noncanonical']:
            smoothed_with_enrich = smoothed_counts.copy()
            smoothed_with_enrich['enrichment'] = smoothed_with_enrich['cluster'].map(cluster_enrichment)
            count_col = f'{feature}_count'
            bp_col = f'{feature}_bp'

            # Calculate frag/kb for reads with this feature
            smoothed_with_enrich['frag_per_kb'] = np.where(
                smoothed_with_enrich[bp_col] > 0,
                smoothed_with_enrich[count_col] / (smoothed_with_enrich[bp_col] / 1000),
                np.nan
            )
            smoothed_with_feat = smoothed_with_enrich[smoothed_with_enrich[bp_col] > 0]

            normal = smoothed_with_feat[smoothed_with_feat['enrichment'] == 'Normal-enriched']['frag_per_kb'].dropna()
            tumor = smoothed_with_feat[smoothed_with_feat['enrichment'] == 'Tumor-enriched']['frag_per_kb'].dropna()

            if len(normal) > 1 and len(tumor) > 1:
                _, pval = mannwhitneyu(normal, tumor)
                fc = tumor.mean() / normal.mean() if normal.mean() > 0 else np.inf
                print(f"    {feature.upper():12s}: N={normal.mean():.2f}, T={tumor.mean():.2f}, FC={fc:.2f}x, p={pval:.2e}")

        print("\n  PRESMOOTHED:")
        for feature in ['its', 'tar1', 'canonical', 'noncanonical']:
            presmooth_with_enrich = presmoothed_counts.copy()
            presmooth_with_enrich['enrichment'] = presmooth_with_enrich['cluster'].map(cluster_enrichment)
            count_col = f'{feature}_count'
            bp_col = f'{feature}_bp'

            # Calculate frag/kb for reads with this feature
            presmooth_with_enrich['frag_per_kb'] = np.where(
                presmooth_with_enrich[bp_col] > 0,
                presmooth_with_enrich[count_col] / (presmooth_with_enrich[bp_col] / 1000),
                np.nan
            )
            presmooth_with_feat = presmooth_with_enrich[presmooth_with_enrich[bp_col] > 0]

            normal = presmooth_with_feat[presmooth_with_feat['enrichment'] == 'Normal-enriched']['frag_per_kb'].dropna()
            tumor = presmooth_with_feat[presmooth_with_feat['enrichment'] == 'Tumor-enriched']['frag_per_kb'].dropna()

            if len(normal) > 1 and len(tumor) > 1:
                _, pval = mannwhitneyu(normal, tumor)
                fc = tumor.mean() / normal.mean() if normal.mean() > 0 else np.inf
                print(f"    {feature.upper():12s}: N={normal.mean():.2f}, T={tumor.mean():.2f}, FC={fc:.2f}x, p={pval:.2e}")

        # Generate normalized fragmentation plot
        plot_normalized_fragmentation(
            smoothed_counts,
            presmoothed_counts,
            cluster_analysis,
            f"{args.output_prefix}.normalized_fragmentation.svg",
            args.dark_mode
        )

        # === Visualizations ===
        print("\n--- Generating visualizations ---")

        plot_its_vs_enrichment(
            cluster_its_tar1,
            f"{args.output_prefix}.its_vs_enrichment.svg",
            args.dark_mode
        )

        plot_cluster_its_distribution(
            cluster_its_tar1,
            f"{args.output_prefix}.cluster_its_distribution.svg",
            args.dark_mode
        )

        plot_its_tar1_ratio(
            cluster_its_tar1,
            f"{args.output_prefix}.its_tar1_ratio.svg",
            args.dark_mode
        )

        plot_its_vs_tar1_scatter(
            cluster_its_tar1,
            f"{args.output_prefix}.its_vs_tar1.svg",
            args.dark_mode
        )

        plot_waterfall(
            cluster_its_tar1,
            f"{args.output_prefix}.waterfall.svg",
            args.dark_mode
        )

        plot_combined_figure(
            cluster_its_tar1,
            f"{args.output_prefix}.combined_figure.svg",
            args.dark_mode,
            solitary_tar1_df=solitary_tar1_df,
            tar1_context_df=tar1_context_df,
            pairwise_df=pairwise_df,
            read_level_features=read_level_features,
            cluster_analysis=cluster_analysis
        )

        plot_solitary_tar1(
            solitary_tar1_df,
            f"{args.output_prefix}.solitary_tar1.svg",
            args.dark_mode
        )

        # === Summary ===
        print("\n" + "=" * 60)
        print("Analysis Complete")
        print("=" * 60)
        print("\nOutput files:")
        print(f"  {args.output_prefix}.cluster_its_tar1.tsv")
        print(f"  {args.output_prefix}.enrichment_comparison.tsv")
        print(f"  {args.output_prefix}.solitary_tar1.tsv")
        print(f"  {args.output_prefix}.tar1_context.tsv")
        print(f"  {args.output_prefix}.pairwise_correlations.tsv")
        print(f"  {args.output_prefix}.cluster_pairwise_correlations.tsv")
        print(f"  {args.output_prefix}.its_contiguity_reads.tsv")
        print(f"  {args.output_prefix}.its_contiguity_clusters.tsv")
        print(f"  {args.output_prefix}.interleaving_reads.tsv")
        print(f"  {args.output_prefix}.interleaving_clusters.tsv")
        print(f"  {args.output_prefix}.fragment_counts.tsv")
        print(f"  {args.output_prefix}.normalized_fragmentation.tsv")
        print(f"  {args.output_prefix}.its_vs_enrichment.svg")
        print(f"  {args.output_prefix}.cluster_its_distribution.svg")
        print(f"  {args.output_prefix}.its_tar1_ratio.svg")
        print(f"  {args.output_prefix}.its_vs_tar1.svg")
        print(f"  {args.output_prefix}.waterfall.svg")
        print(f"  {args.output_prefix}.combined_figure.svg")
        print(f"  {args.output_prefix}.solitary_tar1.svg")
        print(f"  {args.output_prefix}.pairwise_features.svg")
        print(f"  {args.output_prefix}.cluster_pairwise_grid.svg")
        print(f"  {args.output_prefix}.its_contiguity.svg")
        print(f"  {args.output_prefix}.interleaving.svg")
        print(f"  {args.output_prefix}.fragment_counts.svg")
        print(f"  {args.output_prefix}.normalized_fragmentation.svg")
        if args.log_file:
            print(f"  {args.output_prefix}.log")

        return  # End cluster mode

    # =========================================================================
    # SAMPLE MODE - Original analysis (all reads from samples)
    # =========================================================================
    print("\n--- Loading sample metadata ---")
    sample_to_group = load_sample_metadata(args.sample_metadata)
    print(f"  Loaded {len(sample_to_group)} samples from metadata")

    # Parse sample list
    if os.path.exists(args.samples):
        with open(args.samples) as f:
            samples = [line.strip() for line in f if line.strip()]
    else:
        samples = [s.strip() for s in args.samples.split(',')]
    print(f"  Processing {len(samples)} samples: {', '.join(samples)}")

    # --- Load features for all samples ---
    print("\n--- Loading feature data ---")
    all_features = {}
    for sample in samples:
        print(f"  Loading {sample}...", end=' ')
        features = load_sample_features(
            sample, args.bed_prefix, args.database,
            args.smoothness, args.min_read_length, args.min_feature_length
        )
        all_features[sample] = features
        n_with_its = sum(1 for rf in features.values() if rf.get_its_intervals())
        n_with_tar1 = sum(1 for rf in features.values() if rf.get_tar1_intervals())
        print(f"{len(features)} reads ({n_with_its} with ITS, {n_with_tar1} with TAR1)")

    total_reads = sum(len(f) for f in all_features.values())
    print(f"\n  Total: {total_reads} reads")

    # --- Analysis 1: Positional relationships ---
    print("\n" + "=" * 60)
    print("Analysis 1: ITS-TAR1 Positional Relationships")
    print("=" * 60)

    all_pairs = []
    for sample, features in all_features.items():
        pairs = analyze_positional_relationships(features)
        all_pairs.extend(pairs)

    print(f"  Found {len(all_pairs)} ITS-TAR1 pairs on {len(set(p.read_id for p in all_pairs))} reads")

    if all_pairs:
        # Save pairs TSV
        pairs_df = pd.DataFrame([{
            'read_id': p.read_id,
            'sample': p.sample,
            'group': sample_to_group.get(p.sample, 'Unknown'),
            'its_start': p.its_start,
            'its_end': p.its_end,
            'tar1_start': p.tar1_start,
            'tar1_end': p.tar1_end,
            'distance': p.distance,
            'its_before_tar1': p.its_before_tar1,
            'chromosome': p.chromosome,
            'arm': p.arm
        } for p in all_pairs])
        pairs_path = f"{args.output_prefix}.its_tar1_pairs.tsv"
        pairs_df.to_csv(pairs_path, sep='\t', index=False)
        print(f"  Saved: {pairs_path}")

        # Summary statistics
        summary = summarize_positional_analysis(all_pairs, sample_to_group, args.max_distance)
        print(f"\n  Summary by group:")
        print(summary.to_string())

        summary_path = f"{args.output_prefix}.positional_summary.tsv"
        summary.to_csv(summary_path, sep='\t')
        print(f"\n  Saved: {summary_path}")

        # Visualization
        plot_distance_distribution(
            all_pairs, sample_to_group,
            f"{args.output_prefix}.distance_distribution.png",
            args.dark_mode
        )

    # --- Analysis 2: Sample-level quantification ---
    print("\n" + "=" * 60)
    print("Analysis 2: Sample-Level ITS/TAR1 Quantification")
    print("=" * 60)

    sample_stats = calculate_sample_stats(all_features, sample_to_group)
    stats_df, comparison_stats = perform_sample_comparison(sample_stats)

    stats_path = f"{args.output_prefix}.sample_stats.tsv"
    stats_df.to_csv(stats_path, sep='\t', index=False)
    print(f"\n  Saved: {stats_path}")

    plot_sample_abundance(
        sample_stats,
        f"{args.output_prefix}.sample_abundance.png",
        args.dark_mode
    )

    # --- Analysis 3: Chromosome specificity ---
    print("\n" + "=" * 60)
    print("Analysis 3: Chromosome Arm Specificity")
    print("=" * 60)

    chr_df = analyze_chromosome_specificity(all_features, sample_to_group)

    chr_path = f"{args.output_prefix}.chromosome_arm_stats.tsv"
    chr_df.to_csv(chr_path, sep='\t', index=False)
    print(f"  Saved: {chr_path}")

    chr_fc_df = calculate_chromosome_fold_changes(chr_df)
    if not chr_fc_df.empty:
        chr_fc_path = f"{args.output_prefix}.chromosome_fold_changes.tsv"
        chr_fc_df.to_csv(chr_fc_path, sep='\t', index=False)
        print(f"  Saved: {chr_fc_path}")

        plot_chromosome_heatmap(
            chr_fc_df,
            f"{args.output_prefix}.chromosome_heatmap.png",
            args.dark_mode
        )

    # --- Summary ---
    print("\n" + "=" * 60)
    print("Analysis Complete")
    print("=" * 60)
    print("\nOutput files:")
    print(f"  {args.output_prefix}.its_tar1_pairs.tsv")
    print(f"  {args.output_prefix}.positional_summary.tsv")
    print(f"  {args.output_prefix}.sample_stats.tsv")
    print(f"  {args.output_prefix}.chromosome_arm_stats.tsv")
    print(f"  {args.output_prefix}.chromosome_fold_changes.tsv")
    print(f"  {args.output_prefix}.distance_distribution.png")
    print(f"  {args.output_prefix}.sample_abundance.png")
    print(f"  {args.output_prefix}.chromosome_heatmap.png")
    if args.log_file:
        print(f"  {args.output_prefix}.log")


if __name__ == '__main__':
    main()
