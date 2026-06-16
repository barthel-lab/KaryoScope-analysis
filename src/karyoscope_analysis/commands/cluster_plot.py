"""``karyoscope-analysis cluster-plot`` — render Engine B cluster(s) as SVG.

Reads `cluster`'s ``layout.tsv`` (per-segment placements in consensus coordinates) and
``consensus.bed`` (the union-spanning consensus) and draws each cluster as feature-colored read
tracks that stack on the consensus frame (see ``core/cluster_plot.py``). Renders through
``karyoplot`` (drawsvg); the single read-renderer for the package.
"""

from __future__ import annotations

from pathlib import Path

import click

from karyoscope_analysis.core import cluster_plot as render
from karyoscope_analysis.core.io.colors import load_colors


def _read_tsv(path: Path) -> list[dict[str, str]]:
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]


def _major_chromosomes(consensus: list[render.Interval], min_frac: float = 0.10) -> str:
    """All specific chromosomes covering ≥ ``min_frac`` of a consensus, e.g. ``chr4+chr22``.

    Labels translocation clusters with every chromosome they span (not just the dominant one);
    ambiguous chromosome layers (``autosome``/…) and tiny slivers are dropped. '' if none.
    """
    by_chrom: dict[str, int] = {}
    total = 0
    for s, e, feature in consensus:
        total += e - s
        chrom, sep, _ = feature.partition(":")
        if sep and chrom.startswith("chr"):
            by_chrom[chrom] = by_chrom.get(chrom, 0) + (e - s)
    if total == 0:
        return ""
    major = sorted(
        (c for c, v in by_chrom.items() if v >= min_frac * total), key=lambda c: -by_chrom[c]
    )
    return "+".join(major)


@click.command(name="cluster-plot", help="Render Engine B cluster(s) as an SVG of read tracks.")
@click.option(
    "--layout",
    "layout_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="`cluster` layout TSV (cluster_id, read_id, is_seed, reversed, start, end, feature).",
)
@click.option(
    "--consensus",
    "consensus_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="`cluster` consensus BED (cluster_id, start, end, feature, support, coverage).",
)
@click.option(
    "--colors",
    "colors_path",
    default=None,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Database colors.tsv (feature_set, feature, color). Falls back to an auto-palette.",
)
@click.option(
    "--cluster-id",
    default=None,
    help="Render just this cluster (e.g. cluster_0). Omit to render all clusters in one SVG.",
)
@click.option(
    "--min-cluster-size",
    default=2,
    show_default=True,
    type=int,
    help="When rendering all clusters, include only clusters with at least this many reads.",
)
@click.option(
    "--min-segment-bp",
    default=0,
    show_default=True,
    type=int,
    help="Drop feature segments shorter than this when drawing (visual denoise; off by default "
    "so gaps mean 'didn't align').",
)
@click.option(
    "--chromosome-track/--no-chromosome-track",
    default=True,
    show_default=True,
    help="Draw a chromosome-colored track directly under each read's structural track, so "
    "structure and chromosome identity line up (translocations show two chromosome colors).",
)
@click.option(
    "--consensus-track/--no-consensus-track",
    default=True,
    show_default=True,
    help="Draw the union consensus as the top row of each cluster. --no-consensus-track shows "
    "only the read tracks (clearer when judging how the reads align to each other).",
)
@click.option("--width", default=1200, show_default=True, type=int)
@click.option("--row-height", default=12, show_default=True, type=int)
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Output SVG path.",
)
def cmd(
    layout_path: Path,
    consensus_path: Path,
    colors_path: Path | None,
    cluster_id: str | None,
    min_cluster_size: int,
    min_segment_bp: int,
    chromosome_track: bool,
    consensus_track: bool,
    width: int,
    row_height: int,
    output: Path,
) -> None:
    """Render one cluster, or all clusters stacked in one SVG."""
    colors = load_colors(colors_path) if colors_path else {}

    def keep(segs: list[render.Interval]) -> list[render.Interval]:
        return [(s, e, f) for s, e, f in segs if e - s >= min_segment_bp]

    # Group layout rows by cluster, then by read (preserving first-seen order).
    layout: dict[str, dict[str, dict]] = {}
    for r in _read_tsv(layout_path):
        cid, rid = r["cluster_id"], r["read_id"]
        reads = layout.setdefault(cid, {})
        read = reads.setdefault(
            rid, {"is_seed": r["is_seed"] == "1", "reversed": r["reversed"] == "1", "segs": []}
        )
        read["segs"].append((int(r["start"]), int(r["end"]), r["feature"]))

    consensus_by_cluster: dict[str, list[render.Interval]] = {}
    for r in _read_tsv(consensus_path):
        consensus_by_cluster.setdefault(r["cluster_id"], []).append(
            (int(r["start"]), int(r["end"]), r["feature"])
        )

    if cluster_id is not None:
        if cluster_id not in layout:
            raise click.ClickException(f"cluster {cluster_id!r} not found in {layout_path}")
        selected = [cluster_id]
    else:
        selected = [c for c, reads in layout.items() if len(reads) >= min_cluster_size]
        if not selected:
            raise click.ClickException(f"no clusters with >= {min_cluster_size} reads")

    panels: list[render.ClusterPanel] = []
    for cid in selected:
        placed = [
            render.PlacedRead(rid, r["is_seed"], r["reversed"], keep(r["segs"]))
            for rid, r in layout[cid].items()
        ]
        consensus = keep(consensus_by_cluster.get(cid, []))
        span = max(
            [e for r in placed for _s, e, _f in r.segments] + [e for _s, e, _f in consensus] + [1]
        )
        chrom = _major_chromosomes(consensus_by_cluster.get(cid, []))
        title = f"{cid}  n={len(placed)}" + (f"  {chrom}" if chrom else "")
        panels.append(render.ClusterPanel(title=title, width=span, placed=placed, consensus=consensus))

    svg = render.render_clusters_svg(
        panels, colors, width=width, row_height=row_height, chromosome_track=chromosome_track,
        consensus_track=consensus_track,
    )
    output.write_text(svg)
    n_reads = sum(len(p.placed) for p in panels)
    click.echo(f"Rendered {len(panels)} cluster(s), {n_reads} reads, to {output}")
