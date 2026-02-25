#!/usr/bin/env python3
"""
KaryoScope Cluster Translocation Reads

Auto-discovers translocation BED files, groups by data type, merges featuresets,
and runs cluster analysis per (translocation_type, data_type_group, featureset).

Usage:
    python KaryoScope_cluster_translocation_reads.py \\
        --results-dir results \\
        --output-dir outputs/translocation_cluster_analysis \\
        --scripts-dir /path/to/KaryoScope-analysis/scripts \\
        --colors-dir /path/to/KS_human_CHM13 \\
        --sample-metadata samples.tsv
"""

import argparse
import gzip
import json
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

# Capture original command line for logging
_original_command = ' '.join(sys.argv)

import pandas as pd


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


# Default data type -> group mapping
DEFAULT_DATA_TYPE_TO_GROUP = {
    "ONT_frag": "ONT_frag",
    "ONT_nonfrag": "ONT_long",
    "ONT_UL": "ONT_long",
    "hifi_fiber": "HiFi",
    "hifi_notfiber": "HiFi",
}

# Default per-group clustering parameters
DEFAULT_GROUP_PARAMS = {
    "ONT_frag": {"min_read_length": 3000, "max_read_length": 30000, "min_k": 20},
    "ONT_long": {"min_read_length": 10000, "max_read_length": 100000, "min_k": 5},
    "HiFi": {"min_read_length": 10000, "max_read_length": 30000, "min_k": 5},
}


def discover_translocation_beds(results_dir, database, trans_types, featuresets,
                                data_type_to_group, sample_prefix=None):
    """Discover translocation BED files grouped by (sample, data_type, replicate, trans_type).

    Returns dict: (sample, data_type, replicate, trans_type) -> {featureset: Path}
    """
    pattern = re.compile(
        rf"^(.+?)\.(.+?)\.(\d+)\.{re.escape(database)}"
        rf"\.(chr\d+_chr\d+)\.(\w+)\.presmoothed\.translocations\.bed\.gz$"
    )

    combos = defaultdict(dict)

    for sample_dir in sorted(results_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        if sample_prefix is not None and not sample_dir.name.startswith(sample_prefix):
            continue

        glob_pattern = f"*.{database}.*.*.presmoothed.translocations.bed.gz"
        for bed_path in sample_dir.rglob(glob_pattern):
            m = pattern.match(bed_path.name)
            if not m:
                continue

            sample, data_type, replicate, trans_type, featureset = m.groups()
            if trans_type not in trans_types or featureset not in featuresets:
                continue
            if data_type not in data_type_to_group:
                continue

            key = (sample, data_type, replicate, trans_type)
            combos[key][featureset] = bed_path

    return dict(combos)


def merge_region_subtelomeric(region_bed, subtelomeric_bed, output_path, scripts_dir):
    """Merge region + subtelomeric BEDs using KaryoScope_merge_beds.py."""
    cmd = [
        sys.executable,
        str(scripts_dir / "KaryoScope_merge_beds.py"),
        "--bed", str(region_bed), str(subtelomeric_bed),
        "--output", str(output_path),
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)


def concatenate_beds(input_paths, output_path):
    """Concatenate gzipped BED files into a single gzipped output."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output_path, "wt") as out:
        for p in input_paths:
            with gzip.open(p, "rt") as inp:
                for line in inp:
                    out.write(line)


def create_samples_tsv(samples, output_path, metadata_df):
    """Create sample metadata TSV filtered to only samples found in the data.

    Uses the provided metadata DataFrame (with sample, group, color columns).
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Filter metadata to only samples present in data
    filtered = metadata_df[metadata_df["sample"].isin(samples)].copy()
    filtered = filtered.sort_values("sample")

    # Warn about samples missing from metadata
    missing = sorted(set(samples) - set(filtered["sample"]))
    if missing:
        print(f"    WARNING: samples not in metadata (assigned group='unknown', "
              f"color='#999999'): {missing}")
        missing_rows = pd.DataFrame({
            "sample": missing,
            "group": "unknown",
            "color": "#999999",
        })
        filtered = pd.concat([filtered, missing_rows], ignore_index=True)
        filtered = filtered.sort_values("sample")

    filtered[["sample", "group", "color"]].to_csv(output_path, sep="\t", index=False)


def run_cluster_analysis(bed_files, samples_tsv, output_prefix, params, label,
                         scripts_dir):
    """Run KaryoScope_cluster_analysis.py."""
    cmd = [
        sys.executable,
        str(scripts_dir / "KaryoScope_cluster_analysis.py"),
        "--bed", *[str(b) for b in sorted(bed_files)],
        "--sample-metadata", str(samples_tsv),
        "--output-prefix", str(output_prefix),
        "--comparison-mode", "two-group",
        "--control-group", "primary",
        "--min-read-length", str(params["min_read_length"]),
        "--max-read-length", str(params["max_read_length"]),
        "--exclude-features", "canonical_telomere*,novel,unknown",
        "--k-selection", "composite-knee",
        "--min-cluster-size", "3",
        "--min-k", str(params["min_k"]),
        "--reduce-dims", "500",
        "--umap",
        "--circular-dendrogram",
        "--background", "both",
    ]
    print(f"\n{'=' * 80}")
    print(f"CLUSTER ANALYSIS: {output_prefix.name} ({label})")
    print(f"  BED files: {len(bed_files)}")
    print(f"  Params: min_len={params['min_read_length']}, "
          f"max_len={params['max_read_length']}, min_k={params['min_k']}")
    print(f"{'=' * 80}\n")

    result = subprocess.run(cmd, text=True)
    if result.returncode != 0:
        print(f"\n  WARNING: cluster analysis exited with code {result.returncode}")
    return result.returncode


def run_cluster_plot(output_prefix, bed_files, featuresets_str, plot_output,
                     scripts_dir, colors_dir, database):
    """Run KaryoScope_cluster_plot.py."""
    cmd = [
        sys.executable,
        str(scripts_dir / "KaryoScope_cluster_plot.py"),
        "--cluster-analysis-prefix", str(output_prefix),
        "--bed", *[str(b) for b in sorted(bed_files)],
        "--colors", str(colors_dir),
        "--database", database,
        "--featuresets", featuresets_str,
        "--output", str(plot_output),
        "--n-per-cluster", "5",
        "--background", "both",
    ]
    print(f"\n  Generating cluster plot: {plot_output.name}")

    result = subprocess.run(cmd, text=True)
    if result.returncode != 0:
        print(f"  WARNING: cluster plot exited with code {result.returncode}")
    return result.returncode


def main():
    parser = argparse.ArgumentParser(
        description="Auto-discover translocation BED files, group by data type, "
                    "merge featuresets, and run cluster analysis per "
                    "(translocation_type, data_type_group, featureset).",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--results-dir", type=Path, required=True,
                        help="Directory containing per-sample results with translocation BED files")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Output directory for merged BEDs and cluster analysis results")
    parser.add_argument("--scripts-dir", type=Path, required=True,
                        help="Directory containing KaryoScope_cluster_analysis.py, "
                             "KaryoScope_cluster_plot.py, and KaryoScope_merge_beds.py")
    parser.add_argument("--colors-dir", type=Path, required=True,
                        help="Directory containing color files for cluster plotting")
    parser.add_argument("--database", default="KS_human_CHM13",
                        help="KaryoScope database name (default: KS_human_CHM13)")
    parser.add_argument("--translocation-types", nargs="+", dest="translocation_types",
                        default=["chr1_chr21", "chr2_chr13"],
                        help="Translocation types to analyze (default: chr1_chr21 chr2_chr13)")
    parser.add_argument("--featuresets", nargs="+",
                        default=["region", "subtelomeric", "chromosome"],
                        help="Featuresets to discover (default: region subtelomeric chromosome)")
    parser.add_argument("--sample-metadata", type=Path, required=True, dest="sample_metadata",
                        help="TSV file with columns: sample, group, color")
    parser.add_argument("--group-config", type=Path, default=None, dest="group_config",
                        help="Optional JSON file with 'data_type_to_group' and 'group_params' keys.\n"
                             "If not provided, built-in defaults are used.")
    parser.add_argument("--sample-prefix", default=None, dest="sample_prefix",
                        help="Only include sample dirs starting with this prefix (default: all)")
    parser.add_argument("--output-prefix", default="translocation", dest="output_prefix",
                        help="Prefix for output file names (default: translocation)")
    parser.add_argument("--log-file", dest="log_file",
                        action=argparse.BooleanOptionalAction, default=True,
                        help="Save console output to log file (default: True)")

    args = parser.parse_args()

    results_dir = args.results_dir.resolve()
    output_dir = args.output_dir.resolve()
    scripts_dir = args.scripts_dir.resolve()
    colors_dir = args.colors_dir.resolve()
    metadata_path = args.sample_metadata.resolve()

    # Validate required directories and files
    errors = []
    if not results_dir.is_dir():
        errors.append(f"Results directory does not exist: {results_dir}")
    if not scripts_dir.is_dir():
        errors.append(f"Scripts directory does not exist: {scripts_dir}")
    if not colors_dir.is_dir():
        errors.append(f"Colors directory does not exist: {colors_dir}")
    if not metadata_path.is_file():
        errors.append(f"Sample metadata file does not exist: {metadata_path}")
    if errors:
        for e in errors:
            print(e, file=sys.stderr)
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    # Set up logging
    if args.log_file:
        log_path = output_dir / f"{args.output_prefix}_cluster_analysis.log"
        sys.stdout = TeeLogger(str(log_path))

    # Load group config
    if args.group_config is not None:
        config_path = args.group_config.resolve()
        if not config_path.is_file():
            print(f"Group config file does not exist: {config_path}", file=sys.stderr)
            sys.exit(1)
        with open(config_path) as f:
            config = json.load(f)
        data_type_to_group = config.get("data_type_to_group", DEFAULT_DATA_TYPE_TO_GROUP)
        group_params = config.get("group_params", DEFAULT_GROUP_PARAMS)
    else:
        data_type_to_group = DEFAULT_DATA_TYPE_TO_GROUP
        group_params = DEFAULT_GROUP_PARAMS

    # Load sample metadata
    metadata_df = pd.read_csv(metadata_path, sep="\t")
    required_cols = {"sample", "group", "color"}
    if not required_cols.issubset(metadata_df.columns):
        missing_cols = required_cols - set(metadata_df.columns)
        print(f"Sample metadata missing columns: {missing_cols}", file=sys.stderr)
        sys.exit(1)

    # Banner
    print("=" * 60)
    print("KaryoScope Cluster Translocation Reads")
    print("=" * 60)

    print(f"\n{'Parameter':<25} {'Value':<35}")
    print(f"{'-' * 25} {'-' * 35}")
    print(f"{'results-dir':<25} {results_dir}")
    print(f"{'output-dir':<25} {output_dir}")
    print(f"{'scripts-dir':<25} {scripts_dir}")
    print(f"{'colors-dir':<25} {colors_dir}")
    print(f"{'database':<25} {args.database}")
    print(f"{'translocation-types':<25} {' '.join(args.translocation_types)}")
    print(f"{'featuresets':<25} {' '.join(args.featuresets)}")
    print(f"{'sample-metadata':<25} {metadata_path}")
    print(f"{'group-config':<25} {args.group_config or '(defaults)'}")
    print(f"{'sample-prefix':<25} {args.sample_prefix or '(all)'}")
    print(f"{'output-prefix':<25} {args.output_prefix}")
    print(f"{'log-file':<25} {args.log_file}")

    print(f"\n{'=' * 60}")
    print("Command")
    print(f"{'=' * 60}")
    print(_original_command)

    trans_types = args.translocation_types

    # Step 1: Discover translocation BED files
    print("\nStep 1: Discovering translocation BED files...")
    all_beds = discover_translocation_beds(
        results_dir, args.database, trans_types, args.featuresets,
        data_type_to_group, args.sample_prefix,
    )
    print(f"  Found {len(all_beds)} (sample, data_type, replicate, trans_type) combos")

    if not all_beds:
        print("\n  WARNING: No translocation BED files found. Nothing to do.")
        return

    # Group by (data_type_group, trans_type, sample)
    groups = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    for (sample, data_type, replicate, trans_type), featureset_beds in all_beds.items():
        group = data_type_to_group[data_type]
        groups[group][trans_type][sample].append({
            "data_type": data_type,
            "replicate": replicate,
            "beds": featureset_beds,
        })

    print("\n  Data type group summary:")
    for group in sorted(groups):
        for trans_type in trans_types:
            if trans_type in groups[group]:
                samples = sorted(groups[group][trans_type].keys())
                n_combos = sum(len(v) for v in groups[group][trans_type].values())
                print(f"    {group:10s} / {trans_type}: "
                      f"{len(samples)} samples, {n_combos} data_type/replicate combos")

    # Step 2: Merge featuresets and create per-sample BED files
    print("\nStep 2: Merging featuresets and creating per-sample BED files...")

    merged_beds = defaultdict(lambda: defaultdict(dict))
    chromosome_beds = defaultdict(lambda: defaultdict(dict))
    plot_beds = defaultdict(lambda: defaultdict(list))

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        for group in sorted(groups):
            group_dir = output_dir / group

            for trans_type in trans_types:
                if trans_type not in groups[group]:
                    continue

                merged_dir = group_dir / "merged_beds" / trans_type
                chrom_dir = group_dir / "chromosome_beds" / trans_type
                plot_dir = group_dir / "plot_beds" / trans_type
                merged_dir.mkdir(parents=True, exist_ok=True)
                chrom_dir.mkdir(parents=True, exist_ok=True)
                plot_dir.mkdir(parents=True, exist_ok=True)

                for sample in sorted(groups[group][trans_type]):
                    entries = groups[group][trans_type][sample]

                    temp_merged_parts = []
                    temp_region_parts = []
                    temp_subtel_parts = []
                    temp_chrom_parts = []

                    for entry in entries:
                        beds = entry["beds"]
                        dt = entry["data_type"]
                        rep = entry["replicate"]

                        if "region" in beds and "subtelomeric" in beds:
                            tmp_out = tmpdir / f"{sample}.{dt}.{rep}.{trans_type}.merged.bed.gz"
                            merge_region_subtelomeric(
                                beds["region"], beds["subtelomeric"], tmp_out,
                                scripts_dir,
                            )
                            temp_merged_parts.append(tmp_out)
                            print(f"    Merged: {sample}.{dt}.{rep} / {trans_type}")

                        if "region" in beds:
                            temp_region_parts.append(beds["region"])
                        if "subtelomeric" in beds:
                            temp_subtel_parts.append(beds["subtelomeric"])
                        if "chromosome" in beds:
                            temp_chrom_parts.append(beds["chromosome"])

                    if temp_merged_parts:
                        final_merged = merged_dir / f"{sample}.{trans_type}.region_subtelomeric.merged.bed.gz"
                        concatenate_beds(temp_merged_parts, final_merged)
                        merged_beds[group][trans_type][sample] = final_merged

                    if temp_region_parts:
                        final_region = plot_dir / f"{sample}.{trans_type}.region.bed.gz"
                        concatenate_beds(temp_region_parts, final_region)
                        plot_beds[group][trans_type].append(final_region)
                    if temp_subtel_parts:
                        final_subtel = plot_dir / f"{sample}.{trans_type}.subtelomeric.bed.gz"
                        concatenate_beds(temp_subtel_parts, final_subtel)
                        plot_beds[group][trans_type].append(final_subtel)

                    if temp_chrom_parts:
                        final_chrom = chrom_dir / f"{sample}.{trans_type}.chromosome.bed.gz"
                        concatenate_beds(temp_chrom_parts, final_chrom)
                        chromosome_beds[group][trans_type][sample] = final_chrom

    # Step 3: Create sample metadata per data type group
    print("\nStep 3: Creating sample metadata files...")

    for group in sorted(groups):
        group_dir = output_dir / group
        all_samples = set()
        for trans_type in trans_types:
            all_samples.update(merged_beds[group].get(trans_type, {}).keys())
            all_samples.update(chromosome_beds[group].get(trans_type, {}).keys())

        if all_samples:
            samples_tsv = group_dir / "samples.tsv"
            create_samples_tsv(all_samples, samples_tsv, metadata_df)
            print(f"    {group}: {sorted(all_samples)}")

    # Step 4: Run cluster analyses
    print("\nStep 4: Running cluster analyses...")

    analyses_run = 0
    analyses_skipped = 0

    for group in sorted(groups):
        group_dir = output_dir / group
        samples_tsv = group_dir / "samples.tsv"
        params = group_params.get(group)
        if params is None:
            print(f"\n  WARNING: No group_params for group '{group}', skipping.")
            continue

        for trans_type in trans_types:
            # Region + subtelomeric analysis
            beds = merged_beds[group].get(trans_type, {})
            if len(beds) >= 2:
                out_prefix = group_dir / f"{args.output_prefix}_{trans_type}"
                bed_files = list(beds.values())

                rc = run_cluster_analysis(
                    bed_files, samples_tsv, out_prefix, params,
                    "region_subtelomeric", scripts_dir,
                )
                if rc == 0:
                    analyses_run += 1
                    plot_bed_files = plot_beds[group].get(trans_type, [])
                    plot_output = group_dir / f"{args.output_prefix}_{trans_type}.cluster_plot.svg"
                    run_cluster_plot(
                        out_prefix, plot_bed_files,
                        "region,subtelomeric", plot_output,
                        scripts_dir, colors_dir, args.database,
                    )
            else:
                print(f"\n  SKIP: {group}/{trans_type} region_subtelomeric "
                      f"({len(beds)} sample(s), need >= 2)")
                analyses_skipped += 1

            # Chromosome-only analysis
            chrom = chromosome_beds[group].get(trans_type, {})
            if len(chrom) >= 2:
                out_prefix = group_dir / f"{args.output_prefix}_{trans_type}_chromosome"
                chrom_bed_files = list(chrom.values())

                rc = run_cluster_analysis(
                    chrom_bed_files, samples_tsv, out_prefix, params,
                    "chromosome", scripts_dir,
                )
                if rc == 0:
                    analyses_run += 1
                    plot_output = group_dir / f"{args.output_prefix}_{trans_type}_chromosome.cluster_plot.svg"
                    run_cluster_plot(
                        out_prefix, chrom_bed_files,
                        "chromosome", plot_output,
                        scripts_dir, colors_dir, args.database,
                    )
            else:
                print(f"\n  SKIP: {group}/{trans_type} chromosome "
                      f"({len(chrom)} sample(s), need >= 2)")
                analyses_skipped += 1

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print(f"  Analyses completed: {analyses_run}")
    print(f"  Analyses skipped:   {analyses_skipped}")
    print(f"  Output directory:   {output_dir}")
    print("=" * 60)


if __name__ == "__main__":
    main()
