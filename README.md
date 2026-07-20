# CXCL12 regional Hi-C figure generator

This repository contains an R script for making a publication-style regional Hi-C figure centered on **CXCL12**. The figure combines a triangular Hi-C heatmap with loop calls, TAD/domain outlines, CTCF and RAD21 peak tracks, H3K27ac and BRD4 signal tracks, and a filtered gene annotation track.

The code is written for a Juicer `.hic` file and hg38-style genomic coordinates. By default, the script uses a 5 kb resolution and an asymmetric CXCL12 window with an 800 kb left flank and a 560 kb right flank. These values are easy to change in the **USER SETTINGS** section of the script.

## Repository structure

```text
hic-cxcl12-figure/
├── README.md
├── CITATION.cff
├── LICENSE
├── .gitignore
├── data/
│   └── README.md
├── docs/
│   └── github_setup_checklist.md
├── results/
│   └── .gitkeep
└── scripts/
    ├── install_dependencies.R
    └── plot_cxcl12_hic_figure.R
```

## Inputs

For privacy, the public script uses placeholder paths such as `/path/to/inter_30.hic` instead of private lab or cluster paths. Before running the workflow, either replace `/path/to` in your private local copy or pass real paths at runtime with command-line options.

| Input | Default filename | Description |
|---|---|---|
| Hi-C matrix | `inter_30.hic` | Juicer `.hic` file |
| Loops | `merged_loops.bedpe` | Loop calls in BEDPE format |
| Domains/TADs | `10000_blocks.bedpe` | Contact-domain/TAD calls in BEDPE format |
| CTCF peaks | `CTCF.bed` | BED or narrowPeak-style CTCF peaks |
| RAD21 peaks | `RAD21.liver.bed` | BED or narrowPeak-style RAD21 peaks |
| H3K27ac signal | `CtrlVeh_H3K27ac_log2.bw` | BigWig signal track |
| BRD4 signal | `CtrlVeh_BRD4_log2.bw` | BigWig signal track |
| Gene annotation | `genes.gtf` | GTF file containing `CXCL12` |

The included `.gitignore` intentionally excludes large genomics files and generated outputs. For a public GitHub repository, upload code and documentation, not protected or unpublished raw data. Avoid committing absolute paths that reveal private lab, cluster, username, or institution information.

## Install R dependencies

From the repository root:

```r
source("scripts/install_dependencies.R")
```

Or install manually:

```r
install.packages(c("strawr", "ggplot2", "patchwork", "ragg"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("rtracklayer", "GenomicRanges", "IRanges"))
```

## Run the analysis

### Option 1: pass file paths from the command line, recommended

This keeps the public script unchanged while letting you use private local or cluster paths on your own computer.


```bash
Rscript scripts/plot_cxcl12_hic_figure.R \
  --hic=/path/to/inter_30.hic \
  --loops=/path/to/merged_loops.bedpe \
  --domains=/path/to/10000_blocks.bedpe \
  --ctcf=/path/to/CTCF.bed \
  --rad21=/path/to/RAD21.liver.bed \
  --h3k27ac=/path/to/CtrlVeh_H3K27ac_log2.bw \
  --brd4=/path/to/CtrlVeh_BRD4_log2.bw \
  --gtf=/path/to/genes.gtf \
  --out_dir=results \
  --sample_label="Sample 1" \
  --out_prefix="cxcl12_sample1"
```

### Option 2: edit the placeholder paths locally

In your private working copy, replace `/path/to` with the real folder that contains your files. Do not commit your edited private paths back to the public repository.

## Main outputs

The script writes results to `results/` by default:

- `*_regional_hic_hg38_5kb.tiff`: final multi-track Hi-C figure.
- `*_genes_plotted.tsv`: genes included in the gene track.
- `*_CTCF_peaks_in_window.tsv`: CTCF peaks overlapping the plotted window.
- `*_RAD21_peaks_in_window.tsv`: RAD21 peaks overlapping the plotted window.
- `*_CXCL12_loop_connected_genes.tsv`: CXCL12 loop-connected genes used for the gene-track filter.
- `*_HiC_balance_QC_bins.tsv`, `*_HiC_balance_QC_summary.tsv`, and `*_HiC_balance_QC.tiff`: optional balance/QC summaries.
- `*_session_info.txt`: R package versions used for the run.

## Methods summary

1. Import the GTF and identify the target gene coordinates.
2. Define a regional window around the target gene.
3. Extract observed Hi-C contacts from a Juicer `.hic` file with `strawr`.
4. Convert the Hi-C matrix to triangular polygons for plotting with `ggplot2`.
5. Overlay loops, TAD/domain outlines, CTCF/RAD21 peaks, H3K27ac/BRD4 signal, and gene models.
6. Export a high-resolution TIFF figure and TSV tables for reproducibility.

## Notes for adapting this project

- Change `gene_name`, `left_flank_bp`, and `right_flank_bp` in `scripts/plot_cxcl12_hic_figure.R` to make figures for other loci.
- Use `norm_method <- "NONE"` for unbalanced observed counts, or change it to an available normalization such as `"KR"` if your `.hic` file supports it.
- Keep the output TSV files with your manuscript or report so reviewers can see exactly which loops, peaks, domains, and genes were plotted.
- Add a small example dataset only if you have permission to share it publicly. Otherwise, provide a note explaining where public users can obtain equivalent reference data.
- Keep public code paths generic, for example `/path/to/inter_30.hic`, instead of private cluster paths.

## License

This template uses the MIT License by default. Edit `LICENSE` before publishing if your lab or institution requires a different license.
