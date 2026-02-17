#!/usr/bin/env python3
"""
KaryoScope Cluster Annotation

Summarizes the dominant features for each cluster based on BED file annotations.
Annotates each featureset layer separately (e.g., region, subtelomeric, chromosome).

Usage:
  python KaryoScope_cluster_annotate.py \
    --prefix analysis_output_prefix \
    --bed-dir results \
    --output cluster_annotations.tsv

  # With specific featuresets:
  python KaryoScope_cluster_annotate.py \
    --prefix analysis_output_prefix \
    --bed-dir results \
    --featuresets region,subtelomeric,chromosome \
    --output cluster_annotations.tsv

The script automatically finds these files from the prefix:
  - {prefix}.read_assignments.tsv
  - {prefix}.cluster_analysis.tsv
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


def matches_any_pattern(feature, patterns):
    """Check if feature matches any of the exclude patterns (supports wildcards)."""
    for pattern in patterns:
        if fnmatch.fnmatch(feature, pattern):
            return True
    return False


def summarize_featureset(cluster_reads, bed_df, top_n=3, exclude_patterns=None):
    """Summarize features for a cluster from a single featureset."""
    cluster_bed = bed_df[bed_df['read'].isin(cluster_reads)]

    if len(cluster_bed) == 0:
        return ''

    # Count by feature (weighted by length)
    feature_bp = cluster_bed.groupby('feature')['length'].sum()

    if len(feature_bp) == 0:
        return ''

    # Calculate total BEFORE filtering (so percentages reflect true proportion)
    total_bp = feature_bp.sum()

    # Filter out excluded features for display only
    if exclude_patterns:
        feature_bp = feature_bp[~feature_bp.index.map(lambda f: matches_any_pattern(f, exclude_patterns))]

    if len(feature_bp) == 0:
        return ''

    # Top features by coverage (percentages relative to total including excluded)
    top_features = feature_bp.nlargest(top_n)
    top_str = '; '.join([f"{f}({100*v/total_bp:.1f}%)" for f, v in top_features.items()])

    return top_str


def compute_read_feature_fractions(bed_df):
    """Compute per-read feature coverage fractions.

    Returns a DataFrame: read × feature matrix of coverage fractions.
    """
    # Total bp per read
    read_totals = bed_df.groupby('read')['length'].sum()
    # Per read × feature bp
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


def score_cluster_features(cluster_reads, fractions, thresholds):
    """For each feature, compute % of cluster reads exceeding the adaptive threshold."""
    cluster_frac = fractions.reindex(cluster_reads).fillna(0)
    n_reads = len(cluster_reads)
    scores = {}
    for feature, thresh in thresholds.items():
        if feature in cluster_frac.columns:
            n_sig = (cluster_frac[feature] > thresh).sum()
            scores[feature] = round(100 * n_sig / n_reads, 1) if n_reads > 0 else 0
    return scores


def compute_cluster_bp_scores(cluster_reads, bed_df):
    """Compute bp-weighted feature proportions for a cluster.

    Returns dict of {feature: bp_pct} where bp_pct = 100 * feature_bp / total_bp.
    """
    cluster_bed = bed_df[bed_df['read'].isin(cluster_reads)]
    if len(cluster_bed) == 0:
        return {}
    total_bp = cluster_bed['length'].sum()
    if total_bp == 0:
        return {}
    feature_bp = cluster_bed.groupby('feature')['length'].sum()
    return {feat: round(100 * bp / total_bp, 2) for feat, bp in feature_bp.items()}


# --- Feature classification constants for interspersion metrics ---
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


def compute_cluster_ct_bp_pct(cluster_reads, bed_df):
    """Compute bp-weighted ct (centric transition) proportion for a cluster."""
    cluster_bed = bed_df[bed_df['read'].isin(cluster_reads)]
    if len(cluster_bed) == 0:
        return 0.0
    total_bp = cluster_bed['length'].sum()
    if total_bp == 0:
        return 0.0
    ct_bp = cluster_bed.loc[
        cluster_bed['feature'].str.split(':', n=1).str[0] == 'ct', 'length'
    ].sum()
    return round(100 * ct_bp / total_bp, 2)


def compute_cluster_interspersion(cluster_reads, bed_df):
    """Compute typed interspersion rates (transitions per kb) for a cluster.

    Returns dict with median rates across all reads:
      total, can_ncan, tel_sat, arm_tel
    """
    cluster_bed = bed_df[bed_df['read'].isin(cluster_reads)]
    if len(cluster_bed) == 0:
        return {'total': 0.0, 'can_ncan': 0.0, 'tel_sat': 0.0, 'arm_tel': 0.0}

    read_rates = {'total': [], 'can_ncan': [], 'tel_sat': [], 'arm_tel': []}

    for _read_id, read_records in cluster_bed.groupby('read'):
        records = read_records.sort_values('start')
        span_kb = (records['end'].max() - records['start'].min()) / 1000
        if span_kb <= 0:
            for key in read_rates:
                read_rates[key].append(0.0)
            continue

        categories = [classify_bed_feature(f) for f in records['feature']]

        # Total transitions: all category changes
        total = sum(
            1 for i in range(1, len(categories))
            if categories[i] != categories[i - 1]
        )

        # Typed transitions: filter out 'other' (noise features like
        # telomere_like_multigroup1, novel) so they don't break adjacency
        # between meaningful categories
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

        read_rates['total'].append(total / span_kb)
        read_rates['can_ncan'].append(can_ncan / span_kb)
        read_rates['tel_sat'].append(tel_sat / span_kb)
        read_rates['arm_tel'].append(arm_tel / span_kb)

    result = {}
    for key, values in read_rates.items():
        result[key] = round(pd.Series(values).median(), 2) if values else 0.0
    return result


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
        ("bed-dir", str(args.bed_dir), None),
        ("featuresets", _fmt(args.featuresets, "featuresets"), None),
        ("database", _fmt(args.database, "database"), None),
        ("smoothness", _fmt(args.smoothness, "smoothness"), None),
        ("top-n", _fmt(args.top_n, "top_n"), None),
        ("clusters", _fmt(args.clusters, "clusters"), None),
        ("min-size", _fmt(args.min_size, "min_size"), None),
        ("exclude-features", _fmt(args.exclude_features, "exclude_features"), None),
        ("log-file", _fmt(args.log_file, "log_file"), None),
        ("auto-label", _fmt(args.auto_label, "auto_label"), None),
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
    """Auto-label a cluster using a decision tree based on feature scores and interspersion metrics.

    Args:
        row: dict with all annotation columns for one cluster
        featureset_prefix: e.g. 'telomere_region' or 'region'

    Returns:
        Label string, or '' for unlabeled clusters
    """
    pfx = featureset_prefix

    sat_features = {
        'active':    row.get(f'{pfx}__active', 0),
        'monomeric': row.get(f'{pfx}__monomeric', 0),
        'bsat':      row.get(f'{pfx}__bsat', 0),
        'censat':    row.get(f'{pfx}__censat', 0),
        'hsat1A':    row.get(f'{pfx}__hsat1A', 0),
        'hsat2':     row.get(f'{pfx}__hsat2', 0),
        'hsat3':     row.get(f'{pfx}__hsat3', 0),
        'gsat':      row.get(f'{pfx}__gsat', 0),
    }
    max_sat_name = max(sat_features, key=sat_features.get)
    max_sat_score = sat_features[max_sat_name]

    ncan = row.get(f'{pfx}__noncanonical_telomere', 0)
    can  = row.get(f'{pfx}__canonical_telomere', 0)
    its  = row.get(f'{pfx}__ITS', 0)
    tar1 = row.get(f'{pfx}__TAR1', 0)
    ct   = row.get('ct_bp_pct', 0)
    icnc = row.get('interspersion_can_ncan', 0)
    rdna = row.get('acrocentric__rDNA', 0)
    acro = row.get('chromosome__acrocentric_multigroup1', 0)

    sat_label_map = {
        'active': 'active aSat', 'monomeric': 'monomeric aSat',
        'bsat': 'bSat', 'censat': 'CenSat', 'hsat1A': 'HSat1A',
        'hsat2': 'HSat2', 'hsat3': 'HSat3', 'gsat': 'GSat',
    }

    # 1. Satellite-dominant
    if max_sat_score >= 80:
        return sat_label_map[max_sat_name]

    # 2. ECTR (tightened with subtelomeric gate)
    if ct >= 25 and (ncan >= 50 or can >= 50) and tar1 < 50 and its < 50:
        return "ECTR"

    # 3. Type I ALT telomere (tightened: require low TAR1)
    if icnc >= 0.5 and ncan >= 80 and ct < 15 and tar1 < 50:
        return "Type I ALT telomere"

    # 4. ECTR (variant-rich) (tightened: require low subtelomeric signal)
    if icnc >= 0.5 and ncan >= 80 and ct >= 15 and tar1 < 50 and its < 50:
        return "ECTR (variant-rich)"

    # 5. Interspersed telomere
    if icnc >= 0.5 and (ncan >= 50 or can >= 50):
        return "Interspersed telomere"

    # 6. rDNA (require high rDNA + acrocentric chromosomes)
    if rdna >= 90 and acro >= 50:
        return "rDNA"

    # 6b. Interstitial TAR1/ITS (rDNA-associated clusters)
    if rdna >= 50 and acro < 50:
        if tar1 >= 50 or its >= 50:
            if tar1 >= its:
                return "Interstitial TAR1"
            else:
                return "Interstitial ITS"

    # 7. Type II ALT telomere
    if can >= 95 and ncan < 50 and ct < 2 and max_sat_score < 10 and icnc < 0.2:
        return "Type II ALT telomere"

    # 8. Normal telomere
    if can >= 80 and ncan < 50 and ct < 10:
        return "Normal telomere"

    # 9. Subtelomeric (moved up from #11)
    if tar1 >= 50 or its >= 50:
        if tar1 >= its:
            return "Subtelomeric (TAR1-rich)"
        else:
            return "Subtelomeric (ITS-rich)"

    # 10. Normal telomere (ct-rich) (moved down from #10)
    if can >= 50 and ct >= 10:
        return "Normal telomere (ct-rich)"

    # 11. Pericentromeric
    if ct >= 25 and max_sat_score < 30:
        return "Pericentromeric"

    # 12. Telomeric
    if ncan >= 50 or can >= 50:
        return "Telomeric"

    # 13. Fallback
    return ""


def main():
    parser = argparse.ArgumentParser(
        description="Annotate clusters with dominant features per featureset",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--prefix", required=True,
                        help="Analysis prefix (auto-finds {prefix}.read_assignments.tsv, {prefix}.cluster_analysis.tsv)")
    parser.add_argument("--bed-dir", dest="bed_dir", required=True,
                        help="Comma-separated base directories containing sample BED files")
    parser.add_argument("--featuresets", default="region,subtelomeric,chromosome,acrocentric,repeat,gene",
                        help="Comma-separated featuresets to annotate (default: region,subtelomeric,chromosome,acrocentric,repeat,gene)")
    parser.add_argument("--database", default="KS_human_CHM13",
                        help="Database name (default: KS_human_CHM13)")
    parser.add_argument("--smoothness", default="smoothed",
                        choices=["smoothed", "presmoothed"],
                        help="BED file smoothness (default: smoothed)")
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

    # Parse featuresets
    featuresets = [fs.strip() for fs in args.featuresets.split(',')]
    print(f"\nFeaturesets to annotate: {featuresets}")

    # Parse bed directories (comma-separated)
    bed_dirs = [d.strip() for d in args.bed_dir.split(',')]

    # Find BED files for each featureset
    print(f"\nFinding BED files in: {bed_dirs}")
    beds_by_featureset = find_featureset_beds(
        bed_dirs, samples, featuresets, args.database, args.smoothness
    )

    # Load BED data for each featureset
    bed_data = {}
    for fs in featuresets:
        bed_files = beds_by_featureset[fs]
        if not bed_files:
            print(f"  WARNING: No BED files found for featureset '{fs}'")
            continue

        print(f"\n  Loading {fs}: {len(bed_files)} files")
        all_data = []
        for bf in bed_files:
            df = load_bed_file(bf)
            all_data.append(df)

        bed_data[fs] = pd.concat(all_data, ignore_index=True)
        print(f"    Total records: {len(bed_data[fs])}")

    if not bed_data:
        print("ERROR: No BED data loaded")
        sys.exit(1)

    # Compute per-read feature fractions and adaptive thresholds
    feature_fractions = {}
    feature_thresholds = {}
    feature_names = {}  # sorted feature names per featureset
    for fs in featuresets:
        if fs not in bed_data:
            continue
        fractions = compute_read_feature_fractions(bed_data[fs])
        thresholds = compute_adaptive_thresholds(fractions)
        feature_fractions[fs] = fractions
        feature_thresholds[fs] = thresholds
        feature_names[fs] = sorted(thresholds.keys())

        # Print threshold summary
        print(f"\n  Feature thresholds for {fs}:")
        for feat in feature_names[fs]:
            nonzero = fractions[feat][fractions[feat] > 0]
            med = nonzero.median() if len(nonzero) > 0 else 0
            print(f"    {feat}: {thresholds[feat]*100:.2f}% (median {med*100:.1f}%)")

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

        # Interspersion metrics (telomere_region featureset)
        if 'telomere_region' in bed_data:
            interspersion = compute_cluster_interspersion(cluster_reads, bed_data['telomere_region'])
            row['interspersion_total'] = interspersion['total']
            row['interspersion_can_ncan'] = interspersion['can_ncan']
            row['interspersion_tel_sat'] = interspersion['tel_sat']
            row['interspersion_arm_tel'] = interspersion['arm_tel']

        # Annotate each featureset
        for fs in featuresets:
            if fs in bed_data:
                row[f'{fs}_top'] = summarize_featureset(cluster_reads, bed_data[fs], args.top_n, exclude_patterns)
            # Per-feature columns: read-presence then bp-level, interleaved per feature
            scores = {}
            bp_scores = {}
            if fs in feature_fractions:
                scores = score_cluster_features(cluster_reads, feature_fractions[fs], feature_thresholds[fs])
            if fs in bed_data:
                bp_scores = compute_cluster_bp_scores(cluster_reads, bed_data[fs])
            for feat in feature_names.get(fs, []):
                row[f'{fs}__{feat}'] = scores.get(feat, 0)
                row[f'{fs}_bp__{feat}'] = bp_scores.get(feat, 0)

        # ct_bp_pct from generalized bp scores (backward compat with auto_label)
        if 'telomere_region' in bed_data:
            row['ct_bp_pct'] = row.get('telomere_region_bp__ct', 0)

        if args.auto_label:
            pfx = 'telomere_region' if 'telomere_region' in bed_data else 'region'
            row['cluster_name'] = auto_label_cluster(row, pfx)

        results.append(row)

    # Create output DataFrame
    result_df = pd.DataFrame(results)

    # Sort by cluster_id (ascending)
    result_df = result_df.sort_values('cluster_id', ascending=True)

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



if __name__ == "__main__":
    main()
