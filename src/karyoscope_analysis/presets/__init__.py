"""Built-in overlay-annotations resolution presets (YAML).

These ship the four legacy `KaryoScope_merge_beds` priority modes as explicit,
hierarchy-validated resolution specs (decision M2). They are loaded via
:func:`karyoscope_analysis.core.annotation_resolution.load_builtin_preset`. Users
may also supply their own spec file. The basic "overlay" mode (join all featuresets
with ``:``) is the tool's built-in default and is not a preset.
"""
