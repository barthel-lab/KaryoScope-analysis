"""Tests for annotation-BED reading/writing and the C4 invariant."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from karyoscope_analysis.core.io import bed

# A valid, C4-conformant 2-sequence annotation BED (gapless partition per seq).
VALID = (
    "read1\t0\t5\tarm\n"
    "read1\t5\t10\tcentromeric\n"
    "read2\t0\t3\tcanonical_telomere\n"
    "read2\t3\t8\ttelomere_like\n"
)


def _write(path: Path, text: str) -> Path:
    if path.suffix == ".gz":
        with gzip.open(path, "wt", newline="") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return path


def test_read_valid(tmp_path: Path):
    groups = bed.read_annotation_bed(_write(tmp_path / "a.bed", VALID))
    assert list(groups) == ["read1", "read2"]
    assert groups["read1"] == [(0, 5, "arm"), (5, 10, "centromeric")]
    assert groups["read2"] == [(0, 3, "canonical_telomere"), (3, 8, "telomere_like")]


def test_read_gzip(tmp_path: Path):
    groups = bed.read_annotation_bed(_write(tmp_path / "a.bed.gz", VALID))
    assert groups["read1"] == [(0, 5, "arm"), (5, 10, "centromeric")]


def test_skips_blank_and_comment_lines(tmp_path: Path):
    text = "# header comment\n\nread1\t0\t5\tarm\nread1\t5\t10\tct\n"
    groups = bed.read_annotation_bed(_write(tmp_path / "a.bed", text))
    assert groups["read1"] == [(0, 5, "arm"), (5, 10, "ct")]


def test_too_few_fields(tmp_path: Path):
    p = _write(tmp_path / "a.bed", "read1\t0\t5\n")
    with pytest.raises(ValueError, match=r"a\.bed:1: expected >=4"):
        bed.read_annotation_bed(p)


def test_non_integer_coords(tmp_path: Path):
    p = _write(tmp_path / "a.bed", "read1\t0\tfive\tarm\n")
    with pytest.raises(ValueError, match="non-integer coordinates"):
        bed.read_annotation_bed(p)


def test_non_positive_interval(tmp_path: Path):
    p = _write(tmp_path / "a.bed", "read1\t5\t5\tarm\n")
    with pytest.raises(ValueError, match="non-positive interval"):
        bed.read_annotation_bed(p)


def test_gap_is_error(tmp_path: Path):
    p = _write(tmp_path / "a.bed", "read1\t0\t5\tarm\nread1\t6\t10\tct\n")
    with pytest.raises(ValueError, match="gap in seq_id"):
        bed.read_annotation_bed(p)


def test_overlap_is_error(tmp_path: Path):
    p = _write(tmp_path / "a.bed", "read1\t0\t6\tarm\nread1\t5\t10\tct\n")
    with pytest.raises(ValueError, match="overlap in seq_id"):
        bed.read_annotation_bed(p)


def test_non_contiguous_seq_id_is_error(tmp_path: Path):
    text = "read1\t0\t5\tarm\nread2\t0\t5\tarm\nread1\t5\t10\tct\n"
    p = _write(tmp_path / "a.bed", text)
    with pytest.raises(ValueError, match="not contiguous"):
        bed.read_annotation_bed(p)


def test_validate_false_allows_gaps(tmp_path: Path):
    p = _write(tmp_path / "a.bed", "read1\t0\t5\tarm\nread1\t6\t10\tct\n")
    groups = bed.read_annotation_bed(p, validate=False)
    assert groups["read1"] == [(0, 5, "arm"), (6, 10, "ct")]


def test_write_roundtrip_plain_and_gz(tmp_path: Path):
    rows = [("read1", 0, 5, "arm"), ("read1", 5, 10, "centromeric")]
    for name in ("out.bed", "out.bed.gz"):
        path = tmp_path / name
        bed.write_annotation_bed(path, rows)
        assert bed.read_annotation_bed(path)["read1"] == [
            (0, 5, "arm"),
            (5, 10, "centromeric"),
        ]


def test_iter_bed_rows(tmp_path: Path):
    p = _write(tmp_path / "a.bed", VALID)
    assert next(iter(bed.iter_bed_rows(p))) == ("read1", 0, 5, "arm")
