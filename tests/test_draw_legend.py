"""Tests for the ``draw-legend`` CLI command (DB-driven standalone legend)."""

from __future__ import annotations

from pathlib import Path

from karyoscope_analysis.cli import main

COLORS_TSV = Path(__file__).resolve().parent / "data" / "colors.tsv"


def test_draw_legend_grouped_by_featureset(cli_runner, tmp_path):
    out = tmp_path / "legend.svg"
    result = cli_runner.invoke(main, ["draw-legend", "--colors", str(COLORS_TSV), "-o", str(out)])
    assert result.exit_code == 0, result.output
    assert out.exists()
    svg = out.read_text()
    assert svg.lstrip().startswith("<")
    # Feature-set headers and at least one DB color are present.
    assert "region" in svg and "chromosome" in svg
    assert "#fa0000" in svg.lower()  # region/aSat from the DB

    # Regression: every feature must render (grouped layout must not truncate groups
    # past a column boundary). One swatch <rect> per feature, plus the background rect.
    n_features = sum(
        1
        for line in COLORS_TSV.read_text().splitlines()
        if line.strip() and not line.startswith("feature_set")
    )
    assert svg.count("<rect") == n_features + 1


def test_draw_legend_feature_set_filter(cli_runner, tmp_path):
    out = tmp_path / "legend.svg"
    result = cli_runner.invoke(
        main,
        ["draw-legend", "--colors", str(COLORS_TSV), "--feature-set", "subtelomeric", "-o", str(out)],
    )
    assert result.exit_code == 0, result.output
    svg = out.read_text()
    assert "subtelomeric" in svg
    assert "chromosome" not in svg  # other sets excluded


def test_draw_legend_empty_after_filter_errors(cli_runner, tmp_path):
    out = tmp_path / "legend.svg"
    result = cli_runner.invoke(
        main,
        ["draw-legend", "--colors", str(COLORS_TSV), "--feature-set", "nonexistent", "-o", str(out)],
    )
    assert result.exit_code != 0
    assert "no legend entries" in result.output
    assert not out.exists()
