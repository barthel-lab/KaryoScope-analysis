"""Core analysis logic for KaryoScope-analysis.

Pure, importable, side-effect-free logic lives here (clustering, enrichment,
interval algebra, feature-matrix construction, representative selection, etc.),
separated from the thin click command wrappers in
:mod:`karyoscope_analysis.commands`. Readers/writers live in
:mod:`karyoscope_analysis.core.io`.
"""
