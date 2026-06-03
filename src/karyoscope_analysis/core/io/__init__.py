"""I/O readers and writers for KaryoScope-analysis.

Annotation-BED readers/writers (bgzip-aware), result-layout path resolution,
colors/hierarchy loading, and TSV/NPZ readers for clustering artifacts. Keeping
I/O here (mirroring the core ``karyoscope`` engine's ``core/io`` layout) keeps the
logic in :mod:`karyoscope_analysis.core` pure and testable.
"""
