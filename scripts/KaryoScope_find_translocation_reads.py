#!/usr/bin/env python3
"""
KaryoScope Find Translocation Reads

Auto-discovers translocation BED files from KaryoScope results directories,
extracts read IDs and per-chromosome coverage stats.

Usage:
    python KaryoScope_find_translocation_reads.py \\
        --results-dir results \\
        --output-dir outputs \\
        --translocation-types chr1_chr21 chr2_chr13 \\
        --target-chromosomes chr1 chr2 chr13 chr21
"""

import argparse
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

# Capture original command line for logging
_original_command = ' '.join(sys.argv)


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


def discover_combos(results_dir, database, trans_types, sample_prefix=None):
    """
    Auto-discover sample/data_type/replicate/trans_type combos by globbing
    for chromosome translocation BED files under sample dirs.

    Args:
        results_dir: Path to results directory.
        database: Database name (e.g., KS_human_CHM13).
        trans_types: List of translocation types to look for.
        sample_prefix: Optional prefix to filter sample dirs. If None, all dirs included.
    """
    combos = []
    pattern = re.compile(
        rf"^(.+?)\.(.+?)\.(\d+)\.{re.escape(database)}\.(chr\d+_chr\d+)\.chromosome\.presmoothed\.translocations\.bed\.gz$"
    )

    bed_glob = f"*.*.*.{database}.*.chromosome.presmoothed.translocations.bed.gz"

    for sample_dir in sorted(results_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        if sample_prefix and not sample_dir.name.startswith(sample_prefix):
            continue

        for bed_path in sample_dir.rglob(bed_glob):
            m = pattern.match(bed_path.name)
            if not m:
                continue

            sample, data_type, replicate, trans_type = m.groups()
            if trans_type not in trans_types:
                continue

            combos.append({
                "sample": sample,
                "data_type": data_type,
                "replicate": replicate,
                "trans_type": trans_type,
                "bed_path": bed_path,
            })

    combos.sort(key=lambda x: (x["sample"], x["data_type"], x["replicate"], x["trans_type"]))
    return combos


def parse_chromosome_bed(bed_path, target_chromosomes):
    """
    Parse chromosome translocation BED file and calculate per-read coverage.

    Args:
        bed_path: Path to gzipped BED file.
        target_chromosomes: List of chromosomes to track coverage for.

    Returns dict: read_id -> {length: int, chr_coverage: {chr: bp}}
    """
    reads = defaultdict(lambda: {"length": 0, "chr_coverage": defaultdict(int)})
    target_set = set(target_chromosomes)

    with gzip.open(bed_path, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue

            read_id = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            feature = parts[3]

            if end > reads[read_id]["length"]:
                reads[read_id]["length"] = end

            if feature in target_set:
                reads[read_id]["chr_coverage"][feature] += end - start

    return dict(reads)


def build_results(reads, sample, data_type, replicate, trans_type, target_chromosomes):
    """Build result rows from parsed reads with dynamic chromosome columns."""
    results = []
    for read_id, data in reads.items():
        read_length = data["length"]
        if read_length == 0:
            continue

        cc = data["chr_coverage"]
        row = {
            "read_id": read_id,
            "sample": sample,
            "data_type": data_type,
            "replicate": replicate,
            "read_length": read_length,
        }

        for chrom in target_chromosomes:
            bp = cc.get(chrom, 0)
            row[f"{chrom}_bp"] = bp
            row[f"{chrom}_pct"] = 100.0 * bp / read_length

        row["translocation_type"] = trans_type
        results.append(row)

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Find translocation reads from KaryoScope translocation BED files",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--results-dir", required=True,
                        help="Base results directory containing sample subdirectories")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for TSV results")
    parser.add_argument("--output-prefix", default="translocation_reads",
                        help="Output file prefix (default: translocation_reads)")
    parser.add_argument("--database", default="KS_human_CHM13",
                        help="Database name (default: KS_human_CHM13)")
    parser.add_argument("--translocation-types", nargs="+", dest="translocation_types",
                        default=["chr1_chr21", "chr2_chr13"],
                        help="Translocation types to search for (default: chr1_chr21 chr2_chr13)")
    parser.add_argument("--target-chromosomes", nargs="+", dest="target_chromosomes",
                        default=["chr1", "chr2", "chr13", "chr21"],
                        help="Target chromosomes to track coverage for (default: chr1 chr2 chr13 chr21)")
    parser.add_argument("--sample-prefix", dest="sample_prefix", default=None,
                        help="Only include sample dirs starting with this prefix (default: all)")
    parser.add_argument("--log-file", dest="log_file",
                        action=argparse.BooleanOptionalAction, default=True,
                        help="Save console output to {output_dir}/{output_prefix}.log (default: True)")

    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    output_dir = Path(args.output_dir)

    # Validate
    if not results_dir.is_dir():
        print(f"ERROR: Results directory does not exist: {results_dir}", file=sys.stderr)
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    # Set up logging
    if args.log_file:
        log_path = output_dir / f"{args.output_prefix}.log"
        sys.stdout = TeeLogger(str(log_path))

    print("=" * 60)
    print("KaryoScope Find Translocation Reads")
    print("=" * 60)

    print(f"\n{'Parameter':<25} {'Value':<35}")
    print(f"{'-' * 25} {'-' * 35}")
    print(f"{'results-dir':<25} {args.results_dir}")
    print(f"{'output-dir':<25} {args.output_dir}")
    print(f"{'output-prefix':<25} {args.output_prefix}")
    print(f"{'database':<25} {args.database}")
    print(f"{'translocation-types':<25} {' '.join(args.translocation_types)}")
    print(f"{'target-chromosomes':<25} {' '.join(args.target_chromosomes)}")
    print(f"{'sample-prefix':<25} {args.sample_prefix or '(all)'}")
    print(f"{'log-file':<25} {args.log_file}")

    print(f"\n{'=' * 60}")
    print("Command")
    print(f"{'=' * 60}")
    print(_original_command)

    # Discover BED files
    print("\nDiscovering translocation BED files...\n")

    combos = discover_combos(results_dir, args.database, args.translocation_types, args.sample_prefix)
    print(f"  Found {len(combos)} BED files\n")

    if not combos:
        print("  WARNING: No translocation BED files found.")
        print(f"  Searched in: {results_dir}")
        if args.sample_prefix:
            print(f"  Sample prefix filter: {args.sample_prefix}")
        return

    all_results = []
    summary = defaultdict(lambda: defaultdict(int))

    for combo in combos:
        sample = combo["sample"]
        data_type = combo["data_type"]
        replicate = combo["replicate"]
        trans_type = combo["trans_type"]
        bed_path = combo["bed_path"]

        key = f"{sample}.{data_type}.{replicate}"

        if not bed_path.exists():
            print(f"  WARNING: {bed_path} not found")
            continue

        print(f"  {key} / {trans_type}...", end=" ", flush=True)

        reads = parse_chromosome_bed(bed_path, args.target_chromosomes)
        results = build_results(reads, sample, data_type, replicate, trans_type, args.target_chromosomes)
        all_results.extend(results)

        summary[key][trans_type] = len(results)
        summary[key]["total"] = summary[key].get("total", 0) + len(results)

        print(f"{len(reads)} reads")

    # Write combined TSV
    output_path = output_dir / f"{args.output_prefix}.tsv"

    # Build dynamic header
    header = ["read_id", "sample", "data_type", "replicate", "read_length"]
    for chrom in args.target_chromosomes:
        header.extend([f"{chrom}_bp", f"{chrom}_pct"])
    header.append("translocation_type")

    with open(output_path, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in all_results:
            row_vals = [
                r["read_id"],
                r["sample"],
                r["data_type"],
                r["replicate"],
                str(r["read_length"]),
            ]
            for chrom in args.target_chromosomes:
                row_vals.append(str(r[f"{chrom}_bp"]))
                row_vals.append(f"{r[f'{chrom}_pct']:.2f}")
            row_vals.append(r["translocation_type"])
            f.write("\t".join(row_vals) + "\n")

    # Print summary
    col_width = 12
    trans_types = args.translocation_types

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    header_line = f"{'Combo':<45}"
    for tt in trans_types:
        header_line += f" {tt:>{col_width}}"
    header_line += f" {'Total':>8}"
    print(header_line)
    print("-" * 80)

    grand_totals = defaultdict(int)
    grand_total = 0

    for key in sorted(summary):
        s = summary[key]
        line = f"{key:<45}"
        row_total = 0
        for tt in trans_types:
            count = s.get(tt, 0)
            line += f" {count:>{col_width}}"
            grand_totals[tt] += count
            row_total += count
        line += f" {row_total:>8}"
        print(line)
        grand_total += row_total

    print("-" * 80)
    total_line = f"{'TOTAL':<45}"
    for tt in trans_types:
        total_line += f" {grand_totals[tt]:>{col_width}}"
    total_line += f" {grand_total:>8}"
    print(total_line)
    print("=" * 80)

    print(f"\nResults written to: {output_path}")
    print(f"Total translocation reads: {len(all_results)}")

    print("\n" + "=" * 60)
    print("Done")
    print("=" * 60)


if __name__ == "__main__":
    main()
