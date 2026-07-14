"""``karyoscope-analysis pool-samples`` — pool per-sample read BEDs for joint clustering.

Concatenates several per-sample annotation BEDs into one, namespacing each read id as
``{sample}|{read_id}`` so ids stay unique across samples, and writes a read-list TSV
(``read_id``, ``sample``) mapping each pooled read back to its sample. The pooled BED feeds
``cluster`` (one joint clustering of all samples); the read-list feeds ``test-enrichment``.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import click

#: Separator between the sample namespace and the original read id in a pooled read id.
SAMPLE_SEP = "|"


def _open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


@click.command(
    name="pool-samples", help="Pool per-sample read BEDs (namespaced) for joint clustering."
)
@click.option(
    "--bed",
    "bed_specs",
    multiple=True,
    required=True,
    help="Per-sample annotation BED as 'SAMPLE:PATH'. Repeatable.",
)
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Pooled BED output (read ids namespaced as SAMPLE|read_id).",
)
@click.option(
    "--read-list-out",
    "read_list_out",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Read-list TSV (read_id, sample) output. Default: alongside --output as *.samples.tsv.",
)
def cmd(bed_specs: tuple[str, ...], output: Path, read_list_out: Path | None) -> None:
    """Concatenate per-sample BEDs with namespaced read ids; emit the read->sample read-list."""
    specs: list[tuple[str, Path]] = []
    for spec in bed_specs:
        if ":" not in spec:
            raise click.UsageError(f"--bed must be 'SAMPLE:PATH' (got {spec!r})")
        sample, path = spec.split(":", 1)
        specs.append((sample, Path(path)))
        if not Path(path).exists():
            raise click.ClickException(f"BED file not found: {path}")

    read_list_out = read_list_out or output.with_suffix(".samples.tsv")
    # Flush (and free) each sample's pooled-id set once its LAST --bed is processed, so peak
    # memory is bounded by the samples still open, not every read id across every sample.
    last_spec = {sample: i for i, (sample, _) in enumerate(specs)}
    seen: dict[str, set[str]] = {}  # open samples only: sample -> pooled read ids so far
    n_lines = 0
    n_reads = 0
    with output.open("w") as out, read_list_out.open("w") as rl:
        rl.write("read_id\tsample\n")
        for i, (sample, path) in enumerate(specs):
            ids = seen.setdefault(sample, set())
            with _open_text(path) as fh:
                for line in fh:
                    if not line.strip():
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 4:
                        continue
                    pooled_id = f"{sample}{SAMPLE_SEP}{fields[0]}"
                    ids.add(pooled_id)
                    out.write("\t".join([pooled_id, *fields[1:]]) + "\n")
                    n_lines += 1
            if last_spec[sample] == i:  # sample complete: emit its read-list rows, free the set
                for pooled_id in sorted(ids):
                    rl.write(f"{pooled_id}\t{sample}\n")
                n_reads += len(ids)
                del seen[sample]

    click.echo(
        f"Pooled {len(specs)} sample(s), {n_reads} reads, {n_lines} BED rows -> {output}; "
        f"read-list -> {read_list_out}"
    )
