"""Enrichment heatmap: clusters x samples, colored by log2 fold-enrichment.

The capstone figure of the enrichment analysis — one glance shows which structural haplotypes
(rows, labeled by ``cluster-annotate``) concentrate in which sample/group (columns), via the
per-cluster log2 fold-enrichment from ``test-enrichment``. Rows are the enriched clusters, sorted
by their strongest enrichment.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass


@dataclass(frozen=True)
class EnrichmentRow:
    """One cluster's per-group log2 fold-enrichment, plus its label, for the heatmap."""

    cluster_id: str
    n_total: int
    label: str
    log2fc: dict[str, float]  # group -> log2 fold-enrichment (-inf if the group is absent)


def _parse_float(text: str) -> float:
    if text in ("-inf", "-Inf", "-Infinity"):
        return float("-inf")
    try:
        return float(text)
    except ValueError:
        return float("-inf")


def select_rows(
    enrichment: Sequence[Mapping[str, str]],
    groups: Sequence[str],
    labels: Mapping[str, str],
    *,
    enriched_only: bool = True,
    max_clusters: int | None = None,
) -> list[EnrichmentRow]:
    """Pick + order clusters for the heatmap (enriched first, by strongest enrichment desc)."""
    rows: list[EnrichmentRow] = []
    for r in enrichment:
        if enriched_only and r.get("enriched") != "1":
            continue
        log2fc = {g: _parse_float(r.get(f"log2fc_{g}", "-inf")) for g in groups}
        rows.append(
            EnrichmentRow(
                cluster_id=r["cluster_id"],
                n_total=int(r["n_total"]),
                label=labels.get(r["cluster_id"], ""),
                log2fc=log2fc,
            )
        )
    rows.sort(key=lambda row: max(row.log2fc.values()), reverse=True)
    return rows[:max_clusters] if max_clusters else rows


Segment = tuple[int, int, str]  # (start, end, feature)


def _segment_color(feature: str, colors: Mapping[str, str]) -> str:
    """DB color for a (possibly composite ``chrom:structural``) consensus feature.

    Colors by the structural layer; ``novel`` -> white; an unknown feature -> light gray (this is
    a visualization, so it degrades gracefully rather than erroring).
    """
    struct = feature.split(":")[-1]
    if struct in colors:
        return colors[struct]
    return "#ffffff" if struct == "novel" else "#dddddd"


def _draw_consensus_panel(ax, rows, consensus, colors, fg) -> set[str]:
    """Draw each row's consensus as a normalized feature-colored bar; return features shown."""
    shown: set[str] = set()
    for i, row in enumerate(rows):
        segs = sorted(consensus.get(row.cluster_id, []), key=lambda s: s[0])
        if not segs:
            continue
        span_start = segs[0][0]
        span = max(1, max(e for _s, e, _f in segs) - span_start)
        xranges = [((s - span_start) / span, (e - s) / span) for s, e, _f in segs]
        facecolors = [_segment_color(f, colors) for _s, _e, f in segs]
        shown.update(f.split(":")[-1] for _s, _e, f in segs)
        ax.broken_barh(xranges, (i - 0.4, 0.8), facecolors=facecolors, edgecolors="none")
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.set_title("consensus structure (normalized)", color=fg, fontsize=9)
    return shown


def render_heatmap(
    rows: Sequence[EnrichmentRow],
    groups: Sequence[str],
    output_path: str,
    *,
    clamp: float = 4.0,
    dark_mode: bool = False,
    consensus: Mapping[str, Sequence[Segment]] | None = None,
    colors: Mapping[str, str] | None = None,
    sort_key=None,
) -> None:
    """Render the clusters x groups log2-fold-enrichment heatmap to ``output_path``.

    Values are clamped to ``[-clamp, clamp]`` for the diverging color scale (a group a cluster is
    absent from, log2fc = -inf, shows at the floor). The format follows the output extension. When
    ``consensus`` (cluster -> segments) and ``colors`` (feature -> hex) are given, a feature-colored
    consensus-structure panel is drawn beside the heatmap (so the label can be checked against the
    actual structure), with a feature color legend ordered by ``sort_key`` if provided.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from karyoplot.mpl.style import apply_default_style, fg_color
    from matplotlib.patches import Patch

    if not rows:
        raise ValueError("no clusters to plot")

    apply_default_style(dark_mode=dark_mode)
    fg = fg_color(dark_mode)
    show_consensus = bool(consensus and colors)

    def _clamp(v: float) -> float:
        if v == float("-inf"):
            return -clamp
        return max(-clamp, min(clamp, v))

    matrix = np.array([[_clamp(row.log2fc[g]) for g in groups] for row in rows])
    row_labels = [
        f"{r.cluster_id} (n={r.n_total})" + (f"  {r.label}" if r.label else "") for r in rows
    ]

    fig_h = max(2.5, 0.3 * len(rows) + 1.6)
    heat_w = 1.1 * len(groups) + 3.0
    if show_consensus:
        fig_w = heat_w + 6.0
        fig, (ax_cons, ax) = plt.subplots(
            1, 2, figsize=(fig_w, fig_h), sharey=True,
            gridspec_kw={"width_ratios": [6.0, max(1.5, 1.1 * len(groups))], "wspace": 0.04},
        )
        shown = _draw_consensus_panel(ax_cons, rows, consensus, colors, fg)
        label_ax = ax_cons
    else:
        fig_w = max(4.0, heat_w + 1.0)
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        shown = set()
        label_ax = ax

    im = ax.imshow(matrix, cmap="RdBu_r", vmin=-clamp, vmax=clamp, aspect="auto")
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, color=fg)
    ax.set_title("log2 fold vs pool", color=fg, fontsize=9)

    label_ax.set_yticks(range(len(rows)))
    label_ax.set_yticklabels(row_labels, fontsize=7, color=fg)
    if show_consensus:
        ax.tick_params(axis="y", labelleft=False)

    for i, row in enumerate(rows):  # annotate each heatmap cell with its (unclamped) value
        for j, g in enumerate(groups):
            v = row.log2fc[g]
            txt = "·" if v == float("-inf") else f"{v:.1f}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=6,
                    color="black" if abs(_clamp(v)) < clamp * 0.6 else "white")

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("log2 fold-enrichment", color=fg)
    cbar.ax.yaxis.set_tick_params(color=fg)
    for t in cbar.ax.get_yticklabels():
        t.set_color(fg)

    if shown:  # feature color legend for the consensus panel
        feats = sorted(shown, key=sort_key) if sort_key else sorted(shown)
        handles = [Patch(facecolor=_segment_color(f, colors), edgecolor="none", label=f)
                   for f in feats]
        ncol = max(1, min(8, len(feats)))
        fig.legend(handles=handles, loc="lower center", ncol=ncol, fontsize=6,
                   frameon=False, labelcolor=fg, bbox_to_anchor=(0.5, -0.02))

    fig.suptitle("Cluster enrichment by sample", color=fg)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
