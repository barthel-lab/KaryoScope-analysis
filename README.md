# KaryoScope-analysis

Clustering, annotation, cross-sample enrichment, and visualization tools for
[KaryoScope](https://github.com/barthel-lab/KaryoScope) sequence annotations.

## Overview

KaryoScope-analysis is an installable Python package exposing a single
`karyoscope-analysis` command-line interface. It consumes the per-read feature-annotation BEDs
produced by the KaryoScope engine's Snakemake pipeline and turns them into clustered, annotated,
and visualized structural-haplotype calls — a chromosome-end "karyotype" view of long-read data.

It sits on top of two sibling packages: the [KaryoScope](https://github.com/barthel-lab/KaryoScope)
engine (`karyoscope`, which provides the canonical database parsers) and the shared plotting
library [KaryoScope-plotlib](https://github.com/barthel-lab/KaryoScope-plotlib) (`karyoplot`). See
`KaryoScope-plotlib/docs/ECOSYSTEM.md` for how the repositories fit together.

> **Note:** this repository was reorganized from a flat collection of `KaryoScope_*.py` scripts into
> this package + CLI. The legacy scripts and their per-script audit live under `docs/audit/`.

## Install

This package depends on two sibling KaryoScope-ecosystem packages that are **not on PyPI**
(`karyoplot`, `karyoscope`). Install them editable from their sibling checkouts first, then this
package. From the parent directory holding all three repos:

```bash
python3 -m venv KaryoScope-analysis/.venv
source KaryoScope-analysis/.venv/bin/activate
python -m pip install --upgrade pip

pip install -e KaryoScope-plotlib          # provides `karyoplot`
pip install -e KaryoScope                   # provides `karyoscope`
pip install -e 'KaryoScope-analysis[dev]'   # this package + dev tools
pip install ruff                            # pinned in .pre-commit-config.yaml
```

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for the full development setup.

## Commands

Run `karyoscope-analysis --help`, or `karyoscope-analysis <command> --help` for any command.

### Data foundation
| Command | Description |
|---|---|
| `bin-annotations` | Hierarchy-aware rolling-window mode filter that denoises a featureset BED before overlay |
| `overlay-annotations` | Combine per-featureset BEDs into one composite (`chrom:structural`) annotation |
| `build-feature-matrix` | Build a per-read feature matrix from annotations |

### Rearrangements (Engine A)
| Command | Description |
|---|---|
| `detect-rearrangements` | Differential-colocalization rearrangement detection |
| `genome-weights` | Reference-genome information-content feature weights (input to clustering) |

### Clustering & enrichment (Engine B)
| Command | Description |
|---|---|
| `cluster` | Overlap-layout-consensus clustering → clusters, consensus, layout |
| `pool-samples` | Namespace + pool per-sample overlays for one joint clustering |
| `test-enrichment` | Per-cluster cross-sample enrichment (descriptive log2 fold vs pool) |
| `cluster-annotate` | Label each cluster from its consensus structure (hierarchy-derived classes) |
| `select-representatives` | Consensus-as-representative catalog |
| `compare-clusterings` | Compare two clusterings (ARI / NMI + overlap) |

### Visualization
| Command | Description |
|---|---|
| `cluster-plot` | Render clustered reads + consensus per cluster (SVG) |
| `plot-enrichment` | Clusters × samples enrichment heatmap (+ optional consensus-structure panel) |
| `plot-reads` | Per-read feature bars (SVG/PNG; heatmap tracks, grouping, telogator preset) |
| `draw-legend` | Standalone SVG legend from the DB color palette |
| `version` | Print the package version |

## Typical pipeline

Per sample, denoise each featureset and overlay them; then pool, cluster once, and
annotate/enrich/visualize:

```
KaryoScope engine output (per-read featureset BEDs)
        │   bin-annotations ×3  (region, subtelomeric, chromosome)
        ▼
   overlay-annotations  ──►  per-sample composite BED
        │   pool-samples (all samples)
        ▼
   cluster  ──►  clusters.tsv, consensus.bed, layout.tsv
        ├──► test-enrichment ──►  enrichment.tsv ──┐
        ├──► cluster-annotate ──► annot.tsv ───────┼──► plot-enrichment (heatmap + consensus)
        └──► cluster-plot (per-group structure SVGs)┘
```

The whole-sample clustering pipeline is wrapped as `scripts/run_cluster_pipeline.sh`
(`--sample S --prefix P --db DB`). The clustering step needs `--block-min-bp` (blocking index;
without it clustering is O(N²)). See `docs/audit/rearrangement_detection.md` §13 for the runbook.

## Documentation

Built with MkDocs Material and hosted on GitHub Pages:
[https://barthel-lab.github.io/KaryoScope-analysis/](https://barthel-lab.github.io/KaryoScope-analysis/).
The durable design record (conventions, decisions, per-script audit) lives under `docs/audit/`.

```bash
pip install -r requirements.txt   # mkdocs + theme
mkdocs serve                       # preview at http://127.0.0.1:8000/
```

## Fonts

Plot outputs default to the generic `sans-serif` family. The optional **Basic Sans** brand font is
auto-registered if available (`karyoplot.core.fonts.register_fonts()`); missing fonts silently fall
back to `sans-serif`.

## Contributing

Actively maintained by the Barthel Lab. See [`CONTRIBUTING.md`](CONTRIBUTING.md); for questions or
issues, open an issue on GitHub.

## License

GPL-3.0-or-later. See [`LICENSE`](LICENSE).

## Citation

If you use KaryoScope in your research, please cite:

[Citation information to be added]

## Contact

[Barthel Lab](https://www.barthel-lab.com/) ·
[Translational Genomics Research Institute (TGen)](https://www.tgen.org/)
