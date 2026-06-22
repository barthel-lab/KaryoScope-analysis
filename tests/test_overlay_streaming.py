"""Tests for the streaming internals of overlay-annotations.

Cover the properties that make the tool a true single-pass, O(featuresets)-memory
overlay: each input is consumed exactly once as a one-shot iterator, the lockstep
sweep detects sequence-order / coverage / span disagreements, the file and in-memory
entry points agree, and the CLI writes its output atomically (no partial file on error).
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterator
from pathlib import Path

import pytest

from karyoscope_analysis.cli import main
from karyoscope_analysis.core import overlay_annotations as core
from karyoscope_analysis.core.annotation_resolution import load_builtin_preset, load_spec
from karyoscope_analysis.core.feature_vocab import FeatureHierarchy
from karyoscope_analysis.core.io.bed import BedRow, iter_annotation_rows

HIERARCHY_TSV = Path(__file__).resolve().parent / "data" / "hierarchy.tsv"


@pytest.fixture(scope="module")
def h() -> FeatureHierarchy:
    return FeatureHierarchy.from_tsv(HIERARCHY_TSV)


def _composite_spec(h: FeatureHierarchy):
    return load_spec(
        {
            "name": "o",
            "precedence": ["region", "repeat"],
            "rules": [{"emit": {"composite": "all"}}],
        },
        h,
    )


def _counting(rows: list[BedRow], counter: Counter, key: str) -> Iterator[BedRow]:
    """A one-shot iterator that tallies how many rows it actually yields."""
    for row in rows:
        counter[key] += 1
        yield row


def test_each_input_consumed_exactly_once(h):
    """The sweep pulls each track once, fully — and works on one-shot iterators."""
    region = [("r1", 0, 5, "arm"), ("r1", 5, 10, "bSat"), ("r2", 0, 10, "arm")]
    repeat = [("r1", 0, 10, "LINE"), ("r2", 0, 4, "LINE"), ("r2", 4, 10, "nonrepeat")]
    pulled: Counter = Counter()
    streams = {
        "region": _counting(region, pulled, "region"),
        "repeat": _counting(repeat, pulled, "repeat"),
    }
    rows = list(core.overlay_streams(streams, _composite_spec(h), h))
    # Every input row was read exactly once (single full pass, no re-iteration).
    assert pulled["region"] == len(region)
    assert pulled["repeat"] == len(repeat)
    # And the overlay produced sensible composite rows for both sequences.
    assert rows == [
        ("r1", 0, 5, "arm:LINE"),
        ("r1", 5, 10, "bSat:LINE"),
        ("r2", 0, 4, "arm:LINE"),
        ("r2", 4, 10, "arm:nonrepeat"),
    ]


def test_detects_seq_id_order_mismatch(h):
    """Tracks pointing at different sequences at a boundary is an error."""
    streams = {
        "region": iter([("r1", 0, 10, "arm")]),
        "repeat": iter([("r2", 0, 10, "LINE")]),
    }
    with pytest.raises(ValueError, match="same order"):
        list(core.overlay_streams(streams, _composite_spec(h), h))


def test_detects_seq_id_reordering(h):
    """Same sequence set, different order across tracks -> error (order is required)."""
    # Both tracks have exactly {A, B}, but listed in opposite order.
    streams = {
        "region": iter([("A", 0, 10, "arm"), ("B", 0, 10, "arm")]),
        "repeat": iter([("B", 0, 10, "LINE"), ("A", 0, 10, "LINE")]),
    }
    with pytest.raises(ValueError, match="same order"):
        list(core.overlay_streams(streams, _composite_spec(h), h))


def test_detects_early_exhaustion(h):
    """A track that ends before the others (fewer sequences) is an error."""
    streams = {
        "region": iter([("r1", 0, 10, "arm"), ("r2", 0, 10, "arm")]),
        "repeat": iter([("r1", 0, 10, "LINE")]),  # missing r2
    }
    with pytest.raises(ValueError, match="ran out of sequences"):
        list(core.overlay_streams(streams, _composite_spec(h), h))


def test_detects_span_mismatch(h):
    """Tracks of a sequence that don't share a span is an error."""
    streams = {
        "region": iter([("r1", 0, 10, "arm")]),
        "repeat": iter([("r1", 0, 8, "LINE")]),  # span 0-8 != 0-10
    }
    with pytest.raises(ValueError, match="different spans"):
        list(core.overlay_streams(streams, _composite_spec(h), h))


def test_streaming_matches_in_memory_on_fixtures(h, v2_subset_beds):
    """The file-streaming path and the in-memory dict path produce identical rows."""
    from karyoscope_analysis.core.io.bed import read_annotation_bed

    fs_used = ("region", "repeat", "subtelomeric")
    spec = load_builtin_preset("priority", h)

    streamed = list(
        core.overlay_streams(
            {fs: iter_annotation_rows(v2_subset_beds[fs]) for fs in fs_used}, spec, h
        )
    )
    in_memory = list(
        core.overlay_annotations(
            {fs: read_annotation_bed(v2_subset_beds[fs]) for fs in fs_used}, spec, h
        )
    )
    assert streamed == in_memory


def test_cli_no_partial_output_on_error(cli_runner, tmp_path: Path):
    """A mid-stream error leaves no output file (atomic temp-file write)."""
    # r1 overlays fine; r2 has a span mismatch -> error after r1 was already swept.
    region = tmp_path / "region.bed"
    subtel = tmp_path / "subtelomeric.bed"
    region.write_text("r1\t0\t10\tarm\nr2\t0\t10\tarm\n")
    subtel.write_text("r1\t0\t10\tnonsubtelomeric\nr2\t0\t8\tnonsubtelomeric\n")
    out = tmp_path / "out.bed.gz"
    result = cli_runner.invoke(
        main,
        [
            "overlay-annotations",
            "--bed",
            f"region={region}",
            "--bed",
            f"subtelomeric={subtel}",
            "--hierarchy",
            str(HIERARCHY_TSV),
            "-o",
            str(out),
        ],
    )
    assert result.exit_code != 0
    assert "different spans" in result.output
    assert not out.exists(), "a partial output file was left behind"
    # No stray temp files either.
    assert not list(tmp_path.glob("out.bed.gz.*.tmp"))
