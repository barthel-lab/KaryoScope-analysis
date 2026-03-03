# KaryoScope Analysis

Analysis scripts and documentation for KaryoScope clustering methods.

## Overview

This repository contains the downstream analysis pipeline for [KaryoScope](https://github.com/barthel-lab/KaryoScope), a graph-based method for analyzing long-read sequencing data that identifies structural patterns based on genomic feature composition and transition patterns. The scripts here consume KaryoScope Snakemake pipeline outputs (BED files) and perform clustering, enrichment analysis, and visualization.

## Documentation

The documentation is built using MkDocs Material and hosted on GitHub Pages.

**Live Documentation:** [https://barthel-lab.github.io/KaryoScope-analysis/](https://barthel-lab.github.io/KaryoScope-analysis/)

## Scripts

### Core Analysis

| Script | Description |
|--------|-------------|
| `KaryoScope_cluster_analysis.py` | Hierarchical clustering with enrichment testing and sequence assignments |
| `KaryoScope_merge_beds.py` | Merge multiple BED featuresets by position overlay |
| `KaryoScope_annotate_sequences.py` | Join sequence assignments with read names and mapping stats |

### Visualization

| Script | Description |
|--------|-------------|
| `KaryoScope_cluster_plot.py` | Plot representative reads per cluster with sample and cluster annotations |
| `KaryoScope_select_representatives.py` | Select best representative reads per cluster (feeds cluster_plot) |
| `KaryoScope_cluster_diagnostics.py` | Diagnostic stats and plots comparing clusters across metrics |
| `KaryoScope_enrichment_bubbles.py` | Bubble plot of cluster enrichment with curated labels |
| `KaryoScope_volcano_plot.py` | Volcano plot of odds ratio vs q-value per cluster |

### Comparison & Annotation

| Script | Description |
|--------|-------------|
| `KaryoScope_cluster_annotate.py` | Summarize dominant features per cluster from BED annotations |
| `KaryoScope_compare_clusterings.py` | Compare two clustering runs (ARI, NMI, Sankey diagram) |

### Utilities

| Script | Description |
|--------|-------------|
| `KaryoScope_plot_reads.py` | Visualize telomeric reads as vertical bars with region features |
| `KaryoScope_draw_legend.py` | Generate SVG legends from KaryoScope color mapping files |
| `create_panning_animation.py` | Create panning animations from wide or tall PNG images |

### Typical Workflow

```
KaryoScope Snakemake Output (BED files)
        │
        ▼
KaryoScope_merge_beds.py ──► merged BED
        │
        ▼
KaryoScope_cluster_analysis.py ──► clusters, enrichment, assignments
        │
        ├──► KaryoScope_cluster_annotate.py ──► feature summaries per cluster
        ├──► KaryoScope_cluster_diagnostics.py ──► diagnostic plots
        ├──► KaryoScope_enrichment_bubbles.py ──► enrichment bubble plot
        ├──► KaryoScope_volcano_plot.py ──► volcano plot
        ├──► KaryoScope_compare_clusterings.py ──► cross-run comparison
        │
        ├──► KaryoScope_select_representatives.py ──► representative reads
        │         │
        │         ▼
        └──► KaryoScope_cluster_plot.py ──► cluster visualization
```

## Local Development

To build and preview the documentation locally:

```bash
# Install dependencies
pip install -r requirements.txt

# Serve documentation locally
mkdocs serve

# Build static site
mkdocs build
```

The documentation will be available at `http://127.0.0.1:8000/`

### Fonts

Plot outputs use the **Basic Sans** font family, bundled in the `fonts/` directory. The font is registered automatically by scripts that generate figures. No manual installation is required.

## Contributing

This repository is actively maintained by the Barthel Lab. For questions or issues, please open an issue on GitHub.

## License

Documentation and code in this repository are available under the MIT License.

## Citation

If you use KaryoScope in your research, please cite:

[Citation information to be added]

## Contact

[Barthel Lab](https://www.barthel-lab.com/)

[Translational Genomics Research Institute (TGen)](https://www.tgen.org/)
