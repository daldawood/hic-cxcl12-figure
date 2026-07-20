# CXCL12 regional Hi-C figure generator
# Triangular Hi-C heatmap with loops, TAD/domain outlines, CTCF/RAD21 peaks, H3K27ac/BRD4 signal, and gene annotations.
# Default use case: human hg38, Juicer .hic, 5 kb resolution.

# ============================================================
# 1) USER SETTINGS
# ============================================================
# Public GitHub privacy note:
#   The default input paths below intentionally use /path/to/<file name> placeholders
#   instead of private lab, cluster, or institutional paths.
#
# Before running the script, either:
#   1) replace /path/to with your real local/cluster folder in your private copy, or
#   2) keep this script unchanged and pass real paths at runtime with --key=value.
#
# Examples:
#   Rscript scripts/plot_cxcl12_hic_figure.R --data_dir=/path/to --out_dir=results
#   Rscript scripts/plot_cxcl12_hic_figure.R --hic=/path/to/inter_30.hic --loops=/path/to/merged_loops.bedpe
#
# Do not commit private, protected, or very large genomics files to a public GitHub repository.

get_cli_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  key <- paste0("--", name, "=")
  hit <- args[startsWith(args, key)]
  if (length(hit) == 0) return(default)
  sub(key, "", hit[[length(hit)]], fixed = TRUE)
}

detect_project_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_file <- cmd[startsWith(cmd, file_arg)]
  if (length(script_file) > 0) {
    script_file <- sub(file_arg, "", script_file[[1]], fixed = TRUE)
    return(normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

project_dir <- detect_project_dir()
# Default input paths are placeholders for privacy in the public GitHub version.
# Replace /path/to locally or override with command-line options.
data_dir <- get_cli_arg("data_dir", "/path/to")
out_dir <- get_cli_arg("out_dir", file.path(project_dir, "results"))

gene_name <- get_cli_arg("gene", "CXCL12")
resolution <- as.integer(get_cli_arg("resolution", "5000"))
sample_label <- get_cli_arg("sample_label", "Regional Hi-C")

hic_file <- get_cli_arg("hic", file.path(data_dir, "inter_30.hic"))
loops_file <- get_cli_arg("loops", file.path(data_dir, "merged_loops.bedpe"))
tad_file <- get_cli_arg("domains", file.path(data_dir, "10000_blocks.bedpe"))
h3k27ac_bw_file <- get_cli_arg("h3k27ac", file.path(data_dir, "CtrlVeh_H3K27ac_log2.bw"))
brd4_bw_file <- get_cli_arg("brd4", file.path(data_dir, "CtrlVeh_BRD4_log2.bw"))
ctcf_bed_file <- get_cli_arg("ctcf", file.path(data_dir, "CTCF.bed"))
rad21_bed_file <- get_cli_arg("rad21", file.path(data_dir, "RAD21.liver.bed"))
gtf_file <- get_cli_arg("gtf", file.path(data_dir, "genes.gtf"))

out_prefix <- get_cli_arg("out_prefix", paste0(tolower(gene_name), "_regional_hic_hg38_5kb"))
out_tiff <- file.path(out_dir, paste0(out_prefix, ".tiff"))
out_gene_table <- file.path(out_dir, paste0(out_prefix, "_genes_plotted.tsv"))
out_ctcf_table <- file.path(out_dir, paste0(out_prefix, "_CTCF_peaks_in_window.tsv"))
out_rad21_table <- file.path(out_dir, paste0(out_prefix, "_RAD21_peaks_in_window.tsv"))
out_cxcl12_connected_gene_table <- file.path(out_dir, paste0(out_prefix, "_CXCL12_loop_connected_genes.tsv"))
out_cxcl12_loop_gene_pairs_table <- file.path(out_dir, paste0(out_prefix, "_CXCL12_loop_connected_loop_gene_pairs.tsv"))
out_selected_tad_table <- file.path(out_dir, paste0(out_prefix, "_selected_CXCL12_domain.tsv"))
out_domain_candidates_table <- file.path(out_dir, paste0(out_prefix, "_CXCL12_domain_candidates.tsv"))
out_balance_table <- file.path(out_dir, paste0(out_prefix, "_HiC_balance_QC_bins.tsv"))
out_balance_summary <- file.path(out_dir, paste0(out_prefix, "_HiC_balance_QC_summary.tsv"))
out_balance_tiff <- file.path(out_dir, paste0(out_prefix, "_HiC_balance_QC.tiff"))
out_session_info <- file.path(out_dir, paste0(out_prefix, "_session_info.txt"))

# Window setting.
# Asymmetric CXCL12 window requested for this version:
#   - keep the left side expanded to show additional upstream loops/interactions: CXCL12 - 800 kb
#   - expand the right side slightly so the TMEM72/TMEM72-AS1 gene body is visible,
#     while preserving the same overall figure style.
# Change these only if you want a different genomic window.
left_flank_bp <- 800000
right_flank_bp <- 560000

# Do NOT use automatic TAD-based cropping for this version.
# We still plot TAD/domain outlines that overlap this asymmetric fixed window.
use_cxcl12_tad_window <- FALSE
plot_only_selected_tad_outline <- FALSE

# These TAD-selection settings are kept only as a backup if you later set
# use_cxcl12_tad_window <- TRUE again.
cxcl12_tad_choice <- "large_centered"  # options: "large_centered", "smallest", "largest"
large_domain_min_width_bp <- 700000
large_domain_max_width_bp <- 1500000
large_domain_target_width_bp <- 1100000
tad_padding_bp <- 0
fallback_flank_bp <- max(left_flank_bp, right_flank_bp)

# Optional manual override. Leave FALSE unless you want to force exact coordinates.
use_manual_window <- FALSE
manual_region_start <- 43700000
manual_region_end <- 44850000

# No Hi-C normalization: unbalanced observed contacts from the .hic file.
norm_method <- "NONE"
matrix_type <- "observed"

# Figure export settings.
fig_width_in <- 9
fig_height_in <- 7
fig_dpi <- 600
fig_background <- "white"

# Track display settings. These affect only the plotted range, not the underlying data.
h3k27ac_cap_quantile <- NA_real_
brd4_cap_quantile <- NA_real_
ctcf_cap_quantile <- 0.99
rad21_cap_quantile <- 0.99
h3k27ac_floor_zero <- FALSE
brd4_floor_zero <- FALSE
ctcf_floor_zero <- TRUE
rad21_floor_zero <- TRUE

# Hi-C heatmap color-display settings.
# These make the red triangular interaction domains/pyramids richer without changing the matrix itself.
# The Hi-C data are still norm_method = "NONE" and matrix_type = "observed".
heatmap_transform <- "log1p"
heatmap_cap_quantile <- 0.985
heatmap_cap_ignore_near_diagonal_bins <- 2
heatmap_palette <- c("white", "#fff5f0", "#fcbba1", "#fb6a4a", "#de2d26", "#67000d")

# Gene-track settings.
# Removes RP-style clone/predicted genes such as RP11-... and LINC RNA genes,
# while keeping the target gene CXCL12 even if a filter would otherwise match it.
remove_gene_name_regex <- "^(RP[0-9]+-|LINC)"
remove_gene_biotypes <- c("lincRNA")
max_gene_lanes <- 3
gene_label_padding_bp <- 30000
gene_label_size <- 1.55

# In addition to CXCL12-loop-connected genes, keep the TMEM72/TMEM72-AS1 bodies if
# they are present in the GTF and overlap the plotted window. The labels are still
# suppressed; this only ensures the gene body is drawn when visible in the window.
extra_gene_bodies_to_keep <- c("TMEM72", "TMEM72-AS1")

# Make small/right-edge gene bodies visible without changing which genes are shown.
# These display-only widths prevent tiny terminal gene fragments from disappearing
# at publication scale. The exact coordinates are still written to the gene table.
gene_body_min_plot_width_bp <- 18000
gene_exon_min_plot_width_bp <- 5500
tmem_gene_body_min_plot_width_bp <- 65000
gene_body_line_size <- 0.34
gene_body_rect_half_height <- 0.026
gene_exon_rect_half_height <- 0.055

# Keep all loop arcs in the figure, but only label genes with loop evidence to CXCL12.
# A gene is called CXCL12-loop-connected when one loop anchor overlaps the CXCL12 gene body
# plus cxcl12_loop_anchor_padding_bp, and the other anchor overlaps that gene body plus
# partner_gene_anchor_padding_bp. CXCL12 itself is always kept.
show_only_cxcl12_loop_connected_genes <- TRUE
cxcl12_loop_anchor_padding_bp <- 25000
partner_gene_anchor_padding_bp <- 30000

# Final visual tweaks requested:
# - keep the same CTCF/RAD21 peak scaling as the previous figure
# - make CTCF/RAD21 peaks a little wider only for display so they are easier to see
# - make peak colors richer/darker
# - make the gene panel much shorter by labeling only CXCL12-loop-connected genes
ctcf_peak_color <- "#2A0055"
rad21_peak_color <- "#8F0000"
ctcf_min_plot_width_bp <- 9000
rad21_min_plot_width_bp <- 9000
gene_track_top_y <- 0.62
gene_track_bottom_y <- 0.22
gene_label_y_offset <- 0.070
gene_panel_layout_height <- 0.30
gene_h3k27ac_spacer_height <- 0.070
loop_panel_layout_height <- 0.52
loop_to_ctcf_spacer_height <- 0.085
regulatory_track_gap_height <- 0.055

# Publication-display options.
# Peak heights are drawn on a relative 0-1 scale for visual clarity; the raw/capped BED score
# is retained in the output peak tables as score. This avoids unreadable scientific-notation
# axes for BED peak tracks and makes CTCF/RAD21 easier to compare visually.
peak_display_transform <- "sqrt"   # options: "linear", "sqrt", "log1p"
peak_display_rescale <- TRUE
hide_peak_y_axis_numbers <- TRUE
cxcl12_highlight_fill <- "grey75"
cxcl12_highlight_alpha <- 0.28
track_label_size <- 6.8
track_tick_label_size <- 5.4
loop_line_size <- 0.28
cxcl12_loop_line_size <- 0.36
tad_line_size <- 0.24

# Fine-tuning for the gene panel.
# These settings do not affect which genes are plotted; they only move the text labels.
cxcl12_label_below_offset <- 0.145       # put the CXCL12 label below the CXCL12 gene body
znf_label_base_offset <- 0.090           # base offset above the ZNF gene line
znf_label_vertical_step <- 0.145         # vertical separation among clustered ZNF/ZNF-AS labels
znf_label_x_nudge_bp <- 36000            # tiny horizontal separation among clustered ZNF/ZNF-AS labels
right_edge_label_margin_bp <- 125000     # labels this close to the right edge are right-aligned inward
right_edge_label_inset_bp <- 18000       # keep right-edge labels just inside the plot boundary

# Publication gene-label cleanup. Labels are drawn in small white boxes with thin leader lines
# so dense genes such as ZNF32 / ZNF32-AS genes do not overlap the gene bodies.
gene_label_box_fill <- "white"
gene_label_box_alpha <- 0.96
gene_label_leader_color <- "grey70"
gene_label_leader_size <- 0.16

show_loop_anchor_guides <- TRUE

# Balance QC settings.
# The main heatmap stays NONE/observed. If KR is available, it is extracted only as a reference
# so you can visually see what a balanced matrix's bin marginals look like.
make_balance_qc <- TRUE
include_KR_reference_in_balance_QC <- TRUE
balance_reference_norm <- "KR"

# ============================================================
# 2) LOAD PACKAGES
# ============================================================

cran_pkgs <- c("strawr", "ggplot2", "patchwork")
bioc_pkgs <- c("rtracklayer", "GenomicRanges", "IRanges")

missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_cran) > 0) {
  stop("Install these CRAN packages first: ", paste(missing_cran, collapse = ", "), call. = FALSE)
}
if (length(missing_bioc) > 0) {
  stop("Install these Bioconductor packages first: ", paste(missing_bioc, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(strawr)
  library(ggplot2)
  library(patchwork)
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
})

# ============================================================
# 3) HELPER FUNCTIONS
# ============================================================

need_file <- function(path, label) {
  if (!file.exists(path)) stop(label, " not found: ", path, call. = FALSE)
}

strip_chr <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)
add_chr <- function(x) paste0("chr", strip_chr(x))

choose_hic_chromosome <- function(hic_file, chrom_from_gtf) {
  hic_chroms <- strawr::readHicChroms(hic_file)
  hic_names <- if (is.data.frame(hic_chroms)) as.character(hic_chroms[[1]]) else as.character(hic_chroms)
  candidates <- unique(c(as.character(chrom_from_gtf), strip_chr(chrom_from_gtf), add_chr(chrom_from_gtf)))
  hit <- candidates[candidates %in% hic_names][1]
  if (is.na(hit)) {
    stop(
      "Could not match chromosome ", chrom_from_gtf, " to the .hic file. .hic chromosomes start with: ",
      paste(head(hic_names), collapse = ", "),
      call. = FALSE
    )
  }
  hit
}

read_bedpe <- function(path) {
  x <- utils::read.table(
    path,
    sep = "\t",
    comment.char = "#",
    header = FALSE,
    fill = TRUE,
    quote = "",
    stringsAsFactors = FALSE
  )
  
  if (nrow(x) == 0) return(x)
  
  bedpe_names <- c(
    "chr1", "x1", "x2", "chr2", "y1", "y2", "name", "score",
    "strand1", "strand2", "color", "observed", "expectedBL",
    "expectedDonut", "expectedH", "expectedV", "fdrBL", "fdrDonut",
    "fdrH", "fdrV", "numCollapsed", "centroid1", "centroid2", "radius"
  )
  
  names(x)[seq_len(min(ncol(x), length(bedpe_names)))] <- bedpe_names[seq_len(min(ncol(x), length(bedpe_names)))]
  
  numeric_cols <- intersect(
    c(
      "x1", "x2", "y1", "y2", "score", "observed", "expectedBL",
      "expectedDonut", "expectedH", "expectedV", "fdrBL", "fdrDonut",
      "fdrH", "fdrV", "numCollapsed", "centroid1", "centroid2", "radius"
    ),
    names(x)
  )
  
  for (nm in numeric_cols) x[[nm]] <- suppressWarnings(as.numeric(x[[nm]]))
  x
}

read_peak_bed <- function(path, chrom, start_bp, end_bp, floor_zero = TRUE,
                          cap_quantile = 0.99, peak_prefix = "peak", min_plot_width_bp = 0) {
  x <- utils::read.table(
    path,
    sep = "\t",
    comment.char = "#",
    header = FALSE,
    fill = TRUE,
    quote = "",
    stringsAsFactors = FALSE
  )
  
  if (nrow(x) == 0 || ncol(x) < 3) {
    return(data.frame(chrom = character(), start = numeric(), end = numeric(), score = numeric(), name = character()))
  }
  
  bed_names <- c("chrom", "start0", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
  names(x)[seq_len(min(ncol(x), length(bed_names)))] <- bed_names[seq_len(min(ncol(x), length(bed_names)))]
  
  x$start0 <- suppressWarnings(as.numeric(x$start0))
  x$end <- suppressWarnings(as.numeric(x$end))
  
  # BED start is 0-based. Convert to 1-based genomic position for plotting with GTF/Hi-C coordinates.
  x$start <- x$start0 + 1
  
  x <- x[
    strip_chr(x$chrom) == strip_chr(chrom) &
      is.finite(x$start) & is.finite(x$end) &
      x$end >= start_bp & x$start <= end_bp,
  ]
  
  if (nrow(x) == 0) {
    return(data.frame(chrom = character(), start = numeric(), end = numeric(), score = numeric(), name = character()))
  }
  
  # Choose the best available height column for ENCODE BED/narrowPeak-style files.
  if ("signalValue" %in% names(x)) {
    score <- suppressWarnings(as.numeric(x$signalValue))
  } else if ("score" %in% names(x)) {
    score <- suppressWarnings(as.numeric(x$score))
  } else if ("qValue" %in% names(x)) {
    score <- suppressWarnings(as.numeric(x$qValue))
  } else {
    score <- rep(1, nrow(x))
  }
  
  if (all(!is.finite(score) | is.na(score))) {
    score <- rep(1, nrow(x))
  }
  score[!is.finite(score) | is.na(score)] <- 0
  if (floor_zero) score <- pmax(score, 0)
  
  if (!is.na(cap_quantile)) {
    positive_scores <- score[is.finite(score) & score > 0]
    if (length(positive_scores) > 0) {
      cap <- as.numeric(stats::quantile(positive_scores, cap_quantile, na.rm = TRUE, names = FALSE))
      if (is.finite(cap) && cap > 0) score <- pmin(score, cap)
    }
  }
  
  name <- if ("name" %in% names(x)) as.character(x$name) else paste0(peak_prefix, "_", seq_len(nrow(x)))
  
  out <- data.frame(
    chrom = as.character(x$chrom),
    start = pmax(x$start, start_bp),
    end = pmin(x$end, end_bp),
    score = score,
    name = name,
    stringsAsFactors = FALSE
  )
  out$mid <- round((out$start + out$end) / 2)
  
  # Preserve the actual BED coordinates in start/end, but optionally widen only the
  # displayed bar so narrow BED peaks remain visible in a 1.3 Mb regional figure.
  out$plot_start <- out$start
  out$plot_end <- out$end
  if (is.finite(min_plot_width_bp) && min_plot_width_bp > 0) {
    actual_half_width <- pmax((out$end - out$start) / 2, 1)
    display_half_width <- pmax(actual_half_width, min_plot_width_bp / 2)
    out$plot_start <- pmax(start_bp, round(out$mid - display_half_width))
    out$plot_end <- pmin(end_bp, round(out$mid + display_half_width))
  }
  
  out[order(out$start, out$end), ]
}

make_peak_display_scores <- function(peaks, transform = "sqrt", rescale = TRUE) {
  if (nrow(peaks) == 0) {
    peaks$display_score <- numeric()
    return(peaks)
  }
  
  x <- peaks$score
  x[!is.finite(x) | is.na(x)] <- 0
  x <- pmax(x, 0)
  
  display <- switch(
    transform,
    linear = x,
    sqrt = sqrt(x),
    log1p = log1p(x),
    stop("peak_display_transform must be one of: linear, sqrt, log1p", call. = FALSE)
  )
  
  if (isTRUE(rescale)) {
    mx <- max(display, na.rm = TRUE)
    if (is.finite(mx) && mx > 0) display <- display / mx
  }
  
  peaks$display_score <- display
  peaks
}

hic_to_triangles <- function(hic, bin_size) {
  ids <- seq_len(nrow(hic))
  display_value <- if ("display_count" %in% names(hic)) hic$display_count else hic$count
  
  p <- rbind(
    data.frame(id = ids, corner = 1, x0 = hic$x,            y0 = hic$y,            count = hic$count, display_count = display_value),
    data.frame(id = ids, corner = 2, x0 = hic$x + bin_size, y0 = hic$y,            count = hic$count, display_count = display_value),
    data.frame(id = ids, corner = 3, x0 = hic$x + bin_size, y0 = hic$y + bin_size, count = hic$count, display_count = display_value),
    data.frame(id = ids, corner = 4, x0 = hic$x,            y0 = hic$y + bin_size, count = hic$count, display_count = display_value)
  )
  
  p$xplot <- (p$x0 + p$y0) / 2
  p$yplot <- (p$y0 - p$x0) / 2
  p[order(p$id, p$corner), ]
}

make_hic_display_values <- function(hic, bin_size, transform = "log1p", cap_quantile = 0.985,
                                    ignore_near_diagonal_bins = 2) {
  raw <- hic$count
  raw[!is.finite(raw) | is.na(raw)] <- 0
  raw <- pmax(raw, 0)
  
  display <- switch(
    transform,
    linear = raw,
    sqrt = sqrt(raw),
    log1p = log1p(raw),
    stop("heatmap_transform must be one of: linear, sqrt, log1p", call. = FALSE)
  )
  
  # Estimate the display cap mainly from off-diagonal pixels so the very strong diagonal
  # does not wash out the domain/pyramid structure.
  diag_distance_bins <- abs(hic$y - hic$x) / bin_size
  cap_pool <- display[
    is.finite(display) & display > 0 &
      diag_distance_bins > ignore_near_diagonal_bins
  ]
  
  if (length(cap_pool) < 20) {
    cap_pool <- display[is.finite(display) & display > 0]
  }
  
  cap <- NA_real_
  if (!is.na(cap_quantile) && length(cap_pool) > 0) {
    cap <- as.numeric(stats::quantile(cap_pool, cap_quantile, na.rm = TRUE, names = FALSE))
    if (is.finite(cap) && cap > 0) display <- pmin(display, cap)
  }
  
  attr(display, "display_cap") <- cap
  display
}

compute_hic_balance_qc <- function(hic_sparse, region_start, region_end, bin_size, norm_label, chrom_label) {
  # straw returns sparse upper-triangular contacts for chrA:region x chrA:region.
  # For a symmetric matrix marginal, off-diagonal contacts contribute to both bins;
  # diagonal contacts contribute once.
  if (nrow(hic_sparse) == 0) {
    return(data.frame())
  }
  
  first_bin <- floor(region_start / bin_size) * bin_size
  last_bin <- floor((region_end - 1) / bin_size) * bin_size
  bins <- seq(first_bin, last_bin, by = bin_size)
  bins <- bins[bins >= region_start & bins <= region_end]
  bins <- sort(unique(c(bins, hic_sparse$x, hic_sparse$y)))
  bins <- bins[bins >= region_start & bins <= region_end]
  
  marginal_sum <- numeric(length(bins))
  pixels_with_contacts <- integer(length(bins))
  
  ix <- match(hic_sparse$x, bins)
  iy <- match(hic_sparse$y, bins)
  keep <- is.finite(hic_sparse$count) & !is.na(ix) & !is.na(iy)
  ix <- ix[keep]
  iy <- iy[keep]
  counts <- hic_sparse$count[keep]
  
  if (length(counts) > 0) {
    for (k in seq_along(counts)) {
      marginal_sum[ix[k]] <- marginal_sum[ix[k]] + counts[k]
      pixels_with_contacts[ix[k]] <- pixels_with_contacts[ix[k]] + 1L
      if (iy[k] != ix[k]) {
        marginal_sum[iy[k]] <- marginal_sum[iy[k]] + counts[k]
        pixels_with_contacts[iy[k]] <- pixels_with_contacts[iy[k]] + 1L
      }
    }
  }
  
  positive_median <- stats::median(marginal_sum[marginal_sum > 0], na.rm = TRUE)
  if (!is.finite(positive_median) || positive_median <= 0) positive_median <- NA_real_
  
  data.frame(
    chrom = chrom_label,
    bin_start = bins,
    bin_end = bins + bin_size,
    bin_mid = bins + bin_size / 2,
    norm_label = norm_label,
    marginal_sum = marginal_sum,
    marginal_ratio_to_median = marginal_sum / positive_median,
    pixels_with_contacts = pixels_with_contacts,
    stringsAsFactors = FALSE
  )
}

summarize_hic_balance_qc <- function(balance_df) {
  if (nrow(balance_df) == 0) return(data.frame())
  
  do.call(rbind, lapply(split(balance_df, balance_df$norm_label), function(z) {
    positive <- z$marginal_sum[z$marginal_sum > 0 & is.finite(z$marginal_sum)]
    ratio <- z$marginal_ratio_to_median[is.finite(z$marginal_ratio_to_median)]
    
    data.frame(
      norm_label = z$norm_label[1],
      n_bins = nrow(z),
      n_zero_bins = sum(z$marginal_sum <= 0 | !is.finite(z$marginal_sum)),
      median_marginal_sum = if (length(positive) > 0) stats::median(positive) else NA_real_,
      mean_marginal_sum = if (length(positive) > 0) mean(positive) else NA_real_,
      sd_marginal_sum = if (length(positive) > 1) stats::sd(positive) else NA_real_,
      cv_marginal_sum = if (length(positive) > 1 && mean(positive) > 0) stats::sd(positive) / mean(positive) else NA_real_,
      min_ratio_to_median = if (length(ratio) > 0) min(ratio) else NA_real_,
      max_ratio_to_median = if (length(ratio) > 0) max(ratio) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}

save_hic_balance_qc_plot <- function(balance_df, balance_summary, out_file, dpi = 600) {
  if (nrow(balance_df) == 0) return(invisible(NULL))
  
  p <- ggplot(balance_df, aes(x = bin_mid / 1e6, y = marginal_ratio_to_median)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.25) +
    geom_line(linewidth = 0.25, na.rm = TRUE) +
    facet_wrap(~ norm_label, ncol = 1, scales = "free_y") +
    labs(
      title = "Hi-C balance QC across the plotted CXCL12 window",
      subtitle = "Each point is a 5 kb bin. A well-balanced matrix has bin marginal sums closer to the dashed line at 1.",
      x = "Genomic position (Mb)",
      y = "Bin marginal sum / median marginal sum"
    ) +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 7),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(size = 8),
      plot.margin = margin(4, 4, 4, 4)
    )
  
  message("Writing Hi-C balance QC plot: ", out_file)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_tiff(filename = out_file, width = 9, height = 4.5, units = "in", res = dpi,
                   compression = "lzw", background = "white")
  } else {
    grDevices::tiff(filename = out_file, width = 9, height = 4.5, units = "in", res = dpi,
                    compression = "lzw", bg = "white", type = "cairo")
  }
  print(p)
  grDevices::dev.off()
  
  invisible(p)
}

loops_to_dotted_segments <- function(loops) {
  if (nrow(loops) == 0) {
    return(data.frame(
      x = numeric(), y = numeric(), xend = numeric(), yend = numeric(),
      loop_id = integer(), is_cxcl12_loop = logical()
    ))
  }
  
  a1 <- (loops$x1 + loops$x2) / 2
  a2 <- (loops$y1 + loops$y2) / 2
  left <- pmin(a1, a2)
  right <- pmax(a1, a2)
  top_x <- (left + right) / 2
  top_y <- (right - left) / 2
  loop_ids <- if ("plot_loop_id" %in% names(loops)) loops$plot_loop_id else seq_along(left)
  is_cx <- if ("is_cxcl12_loop" %in% names(loops)) as.logical(loops$is_cxcl12_loop) else rep(FALSE, nrow(loops))
  
  data.frame(
    x = c(left, right),
    y = 0,
    xend = c(top_x, top_x),
    yend = c(top_y, top_y),
    loop_id = rep(loop_ids, 2),
    is_cxcl12_loop = rep(is_cx, 2)
  )
}

loops_to_arches <- function(loops, n = 80) {
  if (nrow(loops) == 0) {
    return(data.frame(x = numeric(), y = numeric(), loop_id = integer(), is_cxcl12_loop = logical()))
  }
  
  out <- vector("list", nrow(loops))
  loop_ids <- if ("plot_loop_id" %in% names(loops)) loops$plot_loop_id else seq_len(nrow(loops))
  is_cx <- if ("is_cxcl12_loop" %in% names(loops)) as.logical(loops$is_cxcl12_loop) else rep(FALSE, nrow(loops))
  for (i in seq_len(nrow(loops))) {
    a1 <- mean(c(loops$x1[i], loops$x2[i]))
    a2 <- mean(c(loops$y1[i], loops$y2[i]))
    left <- min(a1, a2)
    right <- max(a1, a2)
    t <- seq(0, 1, length.out = n)
    out[[i]] <- data.frame(
      x = left + (right - left) * t,
      y = 4 * t * (1 - t),
      loop_id = loop_ids[i],
      is_cxcl12_loop = is_cx[i]
    )
  }
  do.call(rbind, out)
}

loops_to_anchor_guides <- function(loops, region_start, region_end) {
  if (nrow(loops) == 0) return(data.frame(x = numeric()))
  anchors <- c((loops$x1 + loops$x2) / 2, (loops$y1 + loops$y2) / 2)
  anchors <- anchors[anchors >= region_start & anchors <= region_end]
  data.frame(x = sort(unique(round(anchors))))
}

tads_to_triangles <- function(tads, region_start, region_end) {
  if (nrow(tads) == 0) return(data.frame(x = numeric(), y = numeric(), tad_id = integer()))
  
  tads <- tads[tads$x2 > region_start & tads$x1 < region_end, ]
  if (nrow(tads) == 0) return(data.frame(x = numeric(), y = numeric(), tad_id = integer()))
  
  tads$start_clip <- pmax(tads$x1, region_start)
  tads$end_clip <- pmin(tads$x2, region_end)
  tads <- tads[tads$end_clip > tads$start_clip, ]
  if (nrow(tads) == 0) return(data.frame(x = numeric(), y = numeric(), tad_id = integer()))
  
  tads$mid <- (tads$start_clip + tads$end_clip) / 2
  tads$height <- (tads$end_clip - tads$start_clip) / 2
  
  tri <- rbind(
    data.frame(x = tads$start_clip, y = 0,           tad_id = seq_len(nrow(tads)), order = 1),
    data.frame(x = tads$mid,        y = tads$height, tad_id = seq_len(nrow(tads)), order = 2),
    data.frame(x = tads$end_clip,   y = 0,           tad_id = seq_len(nrow(tads)), order = 3)
  )
  tri[order(tri$tad_id, tri$order), ]
}

select_cxcl12_tad <- function(tads, chrom, gene_start, gene_end,
                              choice = "large_centered",
                              min_width_bp = 700000,
                              max_width_bp = 1500000,
                              target_width_bp = 1100000) {
  if (nrow(tads) == 0) return(NULL)
  
  z <- tads[
    strip_chr(tads$chr1) == strip_chr(chrom) &
      strip_chr(tads$chr2) == strip_chr(chrom) &
      is.finite(tads$x1) & is.finite(tads$x2),
  ]
  if (nrow(z) == 0) return(NULL)
  
  gene_center <- round((gene_start + gene_end) / 2)
  
  # Prefer domains that contain the CXCL12 center.
  contains_center <- z$x1 <= gene_center & z$x2 >= gene_center
  if (any(contains_center, na.rm = TRUE)) {
    candidates <- z[contains_center, ]
    candidates$width <- candidates$x2 - candidates$x1
    candidates$mid <- (candidates$x1 + candidates$x2) / 2
    candidates$center_distance_bp <- abs(candidates$mid - gene_center)
    candidates$target_width_distance_bp <- abs(candidates$width - target_width_bp)
    
    if (choice == "smallest") {
      candidates <- candidates[order(candidates$width, candidates$center_distance_bp, candidates$x1), ]
      candidates$selection_reason <- "smallest_domain_containing_CXCL12_center"
      return(candidates[1, , drop = FALSE])
    }
    
    if (choice == "largest") {
      candidates <- candidates[order(-candidates$width, candidates$center_distance_bp, candidates$x1), ]
      candidates$selection_reason <- "largest_domain_containing_CXCL12_center"
      return(candidates[1, , drop = FALSE])
    }
    
    if (choice == "large_centered") {
      # First restrict to the red-box-like larger domains.
      large_candidates <- candidates[
        candidates$width >= min_width_bp & candidates$width <= max_width_bp,
      ]
      
      if (nrow(large_candidates) > 0) {
        # Choose the large domain where CXCL12 is closest to the center. Use target width only as a tie-breaker.
        large_candidates <- large_candidates[
          order(large_candidates$center_distance_bp, large_candidates$target_width_distance_bp, -large_candidates$width, large_candidates$x1),
        ]
        large_candidates$selection_reason <- "large_centered_domain_containing_CXCL12_center"
        return(large_candidates[1, , drop = FALSE])
      }
      
      # Fallback if no domain falls within the requested width range:
      # use a combined rank that favors centered domains and domains close to the target width.
      candidates$rank_center <- rank(candidates$center_distance_bp, ties.method = "first")
      candidates$rank_width <- rank(candidates$target_width_distance_bp, ties.method = "first")
      candidates$selection_score <- candidates$rank_center + 0.5 * candidates$rank_width
      candidates <- candidates[order(candidates$selection_score, candidates$center_distance_bp, candidates$target_width_distance_bp), ]
      candidates$selection_reason <- "fallback_best_centered_domain_containing_CXCL12_center"
      return(candidates[1, , drop = FALSE])
    }
    
    stop("cxcl12_tad_choice must be one of: large_centered, smallest, largest", call. = FALSE)
  }
  
  # Fallback: use the domain with the largest overlap with the CXCL12 gene body.
  overlap_bp <- pmax(0, pmin(z$x2, gene_end) - pmax(z$x1, gene_start))
  if (any(overlap_bp > 0, na.rm = TRUE)) {
    candidates <- z[overlap_bp > 0, ]
    candidates$overlap_bp <- overlap_bp[overlap_bp > 0]
    candidates$width <- candidates$x2 - candidates$x1
    candidates$mid <- (candidates$x1 + candidates$x2) / 2
    candidates$center_distance_bp <- abs(candidates$mid - gene_center)
    candidates <- candidates[order(-candidates$overlap_bp, candidates$center_distance_bp, -candidates$width, candidates$x1), ]
    candidates$selection_reason <- "fallback_largest_gene_overlap"
    return(candidates[1, , drop = FALSE])
  }
  
  NULL
}

import_bigwig_region <- function(path, chrom, start_bp, end_bp) {
  choices <- unique(c(as.character(chrom), add_chr(chrom), strip_chr(chrom)))
  for (ch in choices) {
    gr <- GenomicRanges::GRanges(seqnames = ch, ranges = IRanges::IRanges(start_bp, end_bp))
    ans <- tryCatch(
      rtracklayer::import(path, which = gr, as = "GRanges"),
      error = function(e) NULL
    )
    if (!is.null(ans)) return(ans)
  }
  stop("Could not import BigWig region from ", path, ". Check chromosome naming in the BigWig file.", call. = FALSE)
}

bigwig_to_track <- function(path, chrom, start_bp, end_bp, floor_zero = FALSE, cap_quantile = NA_real_) {
  gr <- import_bigwig_region(path, chrom, start_bp, end_bp)
  
  if (length(gr) == 0) {
    signal <- data.frame(start = start_bp, end = end_bp, score = 0)
  } else {
    signal <- data.frame(
      start = start(gr),
      end = end(gr),
      score = as.numeric(mcols(gr)$score)
    )
    signal$score[is.na(signal$score)] <- 0
  }
  
  if (floor_zero) signal$score <- pmax(signal$score, 0)
  
  if (!is.na(cap_quantile)) {
    positive_scores <- signal$score[is.finite(signal$score) & signal$score > 0]
    if (length(positive_scores) > 0) {
      cap <- as.numeric(stats::quantile(positive_scores, cap_quantile, na.rm = TRUE, names = FALSE))
      if (is.finite(cap) && cap > 0) signal$score <- pmin(signal$score, cap)
    }
  }
  
  signal$ymin <- pmin(0, signal$score)
  signal$ymax <- pmax(0, signal$score)
  signal
}

make_signal_plot <- function(signal, label, fill_color, gene_start, gene_end, x_limits, anchor_guides = NULL) {
  p <- ggplot(signal) +
    annotate("rect", xmin = gene_start, xmax = gene_end, ymin = -Inf, ymax = Inf,
             fill = cxcl12_highlight_fill, alpha = cxcl12_highlight_alpha) +
    geom_rect(aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax), fill = fill_color, color = NA)
  
  if (!is.null(anchor_guides) && nrow(anchor_guides) > 0 && isTRUE(show_loop_anchor_guides)) {
    p <- p + geom_vline(
      data = anchor_guides,
      aes(xintercept = x),
      color = "grey55",
      size = 0.14,
      linetype = "dotted",
      alpha = 0.35
    )
  }
  
  p +
    scale_x_continuous(limits = x_limits, labels = function(x) sprintf("%.2f", x / 1e6), expand = expansion(mult = 0)) +
    labs(y = label) +
    theme_classic(base_size = 8) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = track_label_size),
      axis.text.y = if (isTRUE(hide_peak_y_axis_numbers)) element_blank() else element_text(size = track_tick_label_size),
      axis.ticks.y = if (isTRUE(hide_peak_y_axis_numbers)) element_blank() else element_line(size = 0.2),
      plot.margin = margin(0.4, 2, 0.4, 2)
    )
}

make_peak_plot <- function(peaks, label, fill_color, gene_start, gene_end, x_limits, anchor_guides = NULL) {
  p <- ggplot() +
    annotate("rect", xmin = gene_start, xmax = gene_end, ymin = -Inf, ymax = Inf,
             fill = cxcl12_highlight_fill, alpha = cxcl12_highlight_alpha)
  
  if (nrow(peaks) > 0) {
    if (!("plot_start" %in% names(peaks))) peaks$plot_start <- peaks$start
    if (!("plot_end" %in% names(peaks))) peaks$plot_end <- peaks$end
    peaks$plot_score <- if ("display_score" %in% names(peaks)) peaks$display_score else peaks$score
    p <- p + geom_rect(
      data = peaks,
      aes(xmin = plot_start, xmax = plot_end, ymin = 0, ymax = plot_score),
      fill = fill_color,
      color = NA
    )
  } else {
    p <- p + annotate("text", x = mean(x_limits), y = 0.5, label = paste("No", label, "peaks in window"), size = 2.3)
  }
  
  if (!is.null(anchor_guides) && nrow(anchor_guides) > 0 && isTRUE(show_loop_anchor_guides)) {
    p <- p + geom_vline(
      data = anchor_guides,
      aes(xintercept = x),
      color = "grey55",
      size = 0.14,
      linetype = "dotted",
      alpha = 0.35
    )
  }
  
  y_max <- if (nrow(peaks) > 0 && "display_score" %in% names(peaks) && any(is.finite(peaks$display_score))) max(peaks$display_score, na.rm = TRUE) else if (nrow(peaks) > 0 && any(is.finite(peaks$score))) max(peaks$score, na.rm = TRUE) else 1
  if (!is.finite(y_max) || y_max <= 0) y_max <- 1
  
  p +
    scale_x_continuous(limits = x_limits, labels = function(x) sprintf("%.2f", x / 1e6), expand = expansion(mult = 0)) +
    coord_cartesian(ylim = c(0, y_max * 1.05), expand = FALSE, clip = "on") +
    labs(y = label) +
    theme_classic(base_size = 8) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = track_label_size),
      axis.text.y = if (isTRUE(hide_peak_y_axis_numbers)) element_blank() else element_text(size = track_tick_label_size),
      axis.ticks.y = if (isTRUE(hide_peak_y_axis_numbers)) element_blank() else element_line(size = 0.2),
      plot.margin = margin(0.4, 2, 0.4, 2)
    )
}

make_gene_track_data <- function(gtf, gene_name, plot_chrom, region_start, region_end, gene_names_to_keep = NULL) {
  biotype_vec <- if ("plot_gene_biotype" %in% names(mcols(gtf))) as.character(gtf$plot_gene_biotype) else NA_character_
  
  all_genes <- data.frame(
    chrom = as.character(seqnames(gtf)),
    start = start(gtf),
    end = end(gtf),
    strand = as.character(strand(gtf)),
    type = gtf$plot_type,
    gene_name = gtf$plot_gene_name,
    gene_biotype = biotype_vec,
    stringsAsFactors = FALSE
  )
  
  all_genes <- all_genes[
    strip_chr(all_genes$chrom) == strip_chr(plot_chrom) &
      all_genes$end >= region_start & all_genes$start <= region_end &
      !is.na(all_genes$gene_name) & all_genes$gene_name != "",
  ]
  
  gene_rows <- all_genes[all_genes$type %in% "gene", ]
  if (nrow(gene_rows) == 0) {
    gene_rows <- do.call(rbind, lapply(split(all_genes, all_genes$gene_name), function(z) {
      data.frame(
        chrom = z$chrom[1],
        start = min(z$start),
        end = max(z$end),
        strand = z$strand[1],
        type = "gene",
        gene_name = z$gene_name[1],
        gene_biotype = z$gene_biotype[which(!is.na(z$gene_biotype))[1]],
        stringsAsFactors = FALSE
      )
    }))
  }
  
  # Collapse duplicate gene rows by gene name.
  gene_rows <- do.call(rbind, lapply(split(gene_rows, gene_rows$gene_name), function(z) {
    biotype <- unique(z$gene_biotype[!is.na(z$gene_biotype) & z$gene_biotype != ""])
    if (length(biotype) == 0) biotype <- NA_character_
    data.frame(
      chrom = z$chrom[1],
      start = min(z$start),
      end = max(z$end),
      strand = z$strand[1],
      type = "gene",
      gene_name = z$gene_name[1],
      gene_biotype = biotype[1],
      stringsAsFactors = FALSE
    )
  }))
  
  gene_rows$is_target <- gene_rows$gene_name == gene_name
  keep_by_requested_name <- rep(FALSE, nrow(gene_rows))
  if (!is.null(gene_names_to_keep)) {
    gene_names_to_keep <- unique(c(gene_name, as.character(gene_names_to_keep)))
    keep_by_requested_name <- gene_rows$gene_name %in% gene_names_to_keep
  }
  
  remove_by_name <- grepl(remove_gene_name_regex, gene_rows$gene_name, ignore.case = TRUE)
  remove_by_biotype <- !is.na(gene_rows$gene_biotype) & gene_rows$gene_biotype %in% remove_gene_biotypes
  
  # Requested genes, including CXCL12-loop-connected genes and forced TMEM bodies,
  # are allowed through even if their biotype/name would otherwise be filtered.
  gene_rows <- gene_rows[gene_rows$is_target | keep_by_requested_name | !(remove_by_name | remove_by_biotype), ]
  
  if (!is.null(gene_names_to_keep)) {
    gene_rows <- gene_rows[gene_rows$is_target | gene_rows$gene_name %in% gene_names_to_keep, ]
  }
  
  gene_rows <- gene_rows[order(gene_rows$start, gene_rows$end), ]
  
  # Assign genes to lanes so nearby gene labels are more readable.
  if (nrow(gene_rows) > 0) {
    lane_end <- rep(-Inf, max_gene_lanes)
    gene_rows$lane <- NA_integer_
    for (i in seq_len(nrow(gene_rows))) {
      possible <- which(gene_rows$start[i] > lane_end + gene_label_padding_bp)
      lane <- if (length(possible) > 0) possible[1] else which.min(lane_end)
      gene_rows$lane[i] <- lane
      lane_end[lane] <- max(lane_end[lane], gene_rows$end[i])
    }
    used_lanes <- max(gene_rows$lane, na.rm = TRUE)
    gene_rows$y <- gene_track_top_y -
      (gene_rows$lane - 1) * ((gene_track_top_y - gene_track_bottom_y) / max(1, used_lanes - 1))
    if (used_lanes == 1) gene_rows$y <- mean(c(gene_track_top_y, gene_track_bottom_y))
    gene_rows$label_x <- pmin(pmax((gene_rows$start + gene_rows$end) / 2, region_start + gene_label_padding_bp), region_end - gene_label_padding_bp)
  }
  
  exon_rows <- all_genes[all_genes$type %in% "exon" & all_genes$gene_name %in% gene_rows$gene_name, ]
  if (nrow(exon_rows) > 0) {
    exon_rows <- merge(exon_rows, gene_rows[, c("gene_name", "y", "is_target")], by = "gene_name")
  } else {
    exon_rows <- data.frame(
      gene_name = character(), chrom = character(), start = numeric(), end = numeric(),
      strand = character(), type = character(), gene_biotype = character(), y = numeric(), is_target = logical()
    )
  }
  
  list(gene_rows = gene_rows, exon_rows = exon_rows)
}

make_collapsed_gene_rows <- function(gtf, gene_name, plot_chrom, region_start, region_end) {
  biotype_vec <- if ("plot_gene_biotype" %in% names(mcols(gtf))) as.character(gtf$plot_gene_biotype) else NA_character_
  
  all_genes <- data.frame(
    chrom = as.character(seqnames(gtf)),
    start = start(gtf),
    end = end(gtf),
    strand = as.character(strand(gtf)),
    type = gtf$plot_type,
    gene_name = gtf$plot_gene_name,
    gene_biotype = biotype_vec,
    stringsAsFactors = FALSE
  )
  
  all_genes <- all_genes[
    strip_chr(all_genes$chrom) == strip_chr(plot_chrom) &
      all_genes$end >= region_start & all_genes$start <= region_end &
      !is.na(all_genes$gene_name) & all_genes$gene_name != "",
  ]
  
  gene_rows <- all_genes[all_genes$type %in% "gene", ]
  if (nrow(gene_rows) == 0) {
    gene_rows <- do.call(rbind, lapply(split(all_genes, all_genes$gene_name), function(z) {
      data.frame(
        chrom = z$chrom[1],
        start = min(z$start),
        end = max(z$end),
        strand = z$strand[1],
        type = "gene",
        gene_name = z$gene_name[1],
        gene_biotype = z$gene_biotype[which(!is.na(z$gene_biotype))[1]],
        stringsAsFactors = FALSE
      )
    }))
  }
  
  if (nrow(gene_rows) == 0) return(gene_rows)
  
  gene_rows <- do.call(rbind, lapply(split(gene_rows, gene_rows$gene_name), function(z) {
    biotype <- unique(z$gene_biotype[!is.na(z$gene_biotype) & z$gene_biotype != ""])
    if (length(biotype) == 0) biotype <- NA_character_
    data.frame(
      chrom = z$chrom[1],
      start = min(z$start),
      end = max(z$end),
      strand = z$strand[1],
      type = "gene",
      gene_name = z$gene_name[1],
      gene_biotype = biotype[1],
      stringsAsFactors = FALSE
    )
  }))
  
  gene_rows$is_target <- gene_rows$gene_name == gene_name
  remove_by_name <- grepl(remove_gene_name_regex, gene_rows$gene_name, ignore.case = TRUE)
  remove_by_biotype <- !is.na(gene_rows$gene_biotype) & gene_rows$gene_biotype %in% remove_gene_biotypes
  gene_rows <- gene_rows[gene_rows$is_target | !(remove_by_name | remove_by_biotype), ]
  gene_rows[order(gene_rows$start, gene_rows$end), ]
}

find_cxcl12_loop_connected_genes <- function(gtf, loops, gene_name, plot_chrom,
                                             gene_start, gene_end, region_start, region_end,
                                             cxcl12_padding_bp = 25000,
                                             partner_padding_bp = 30000) {
  genes <- make_collapsed_gene_rows(gtf, gene_name, plot_chrom, region_start, region_end)
  if (nrow(genes) == 0) {
    return(list(gene_names = gene_name, loop_gene_pairs = data.frame()))
  }
  
  genes$is_target <- genes$gene_name == gene_name
  partner_genes <- genes[!genes$is_target, ]
  
  if (nrow(loops) == 0 || nrow(partner_genes) == 0) {
    return(list(gene_names = gene_name, loop_gene_pairs = data.frame()))
  }
  
  overlap_one <- function(a_start, a_end, b_start, b_end) {
    is.finite(a_start) & is.finite(a_end) & a_end >= b_start & a_start <= b_end
  }
  
  cx_start <- gene_start - cxcl12_padding_bp
  cx_end <- gene_end + cxcl12_padding_bp
  out <- list()
  idx <- 1L
  
  for (i in seq_len(nrow(loops))) {
    loop_id <- if ("plot_loop_id" %in% names(loops)) loops$plot_loop_id[i] else i
    
    x_start <- loops$x1[i]
    x_end <- loops$x2[i]
    y_start <- loops$y1[i]
    y_end <- loops$y2[i]
    
    x_hits_cxcl12 <- overlap_one(x_start, x_end, cx_start, cx_end)
    y_hits_cxcl12 <- overlap_one(y_start, y_end, cx_start, cx_end)
    
    partner_intervals <- data.frame(
      anchor = character(),
      start = numeric(),
      end = numeric(),
      stringsAsFactors = FALSE
    )
    
    if (isTRUE(x_hits_cxcl12) && !isTRUE(y_hits_cxcl12)) {
      partner_intervals <- rbind(partner_intervals, data.frame(anchor = "anchor2", start = y_start, end = y_end))
    }
    if (isTRUE(y_hits_cxcl12) && !isTRUE(x_hits_cxcl12)) {
      partner_intervals <- rbind(partner_intervals, data.frame(anchor = "anchor1", start = x_start, end = x_end))
    }
    if (isTRUE(x_hits_cxcl12) && isTRUE(y_hits_cxcl12)) {
      partner_intervals <- rbind(
        partner_intervals,
        data.frame(anchor = "anchor1", start = x_start, end = x_end),
        data.frame(anchor = "anchor2", start = y_start, end = y_end)
      )
    }
    
    if (nrow(partner_intervals) == 0) next
    
    for (j in seq_len(nrow(partner_intervals))) {
      p_start <- partner_intervals$start[j] - partner_padding_bp
      p_end <- partner_intervals$end[j] + partner_padding_bp
      
      hit <- partner_genes$end >= p_start & partner_genes$start <= p_end
      if (!any(hit, na.rm = TRUE)) next
      
      z <- partner_genes[hit, ]
      for (k in seq_len(nrow(z))) {
        out[[idx]] <- data.frame(
          loop_id = loop_id,
          cxcl12_anchor_padding_bp = cxcl12_padding_bp,
          partner_gene_anchor_padding_bp = partner_padding_bp,
          cxcl12_start = gene_start,
          cxcl12_end = gene_end,
          partner_anchor = partner_intervals$anchor[j],
          partner_anchor_start = partner_intervals$start[j],
          partner_anchor_end = partner_intervals$end[j],
          gene_name = z$gene_name[k],
          gene_start = z$start[k],
          gene_end = z$end[k],
          gene_biotype = z$gene_biotype[k],
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  
  loop_gene_pairs <- if (length(out) > 0) do.call(rbind, out) else data.frame()
  connected_gene_names <- if (nrow(loop_gene_pairs) > 0) sort(unique(loop_gene_pairs$gene_name)) else character()
  
  list(
    gene_names = unique(c(gene_name, connected_gene_names)),
    loop_gene_pairs = loop_gene_pairs
  )
}


# ============================================================
# 4) CHECK INPUTS, FIND CXCL12, AND SET A FIXED CENTERED WINDOW
# ============================================================

# If BED files were moved into the output folder, use those copies automatically.
if (!file.exists(ctcf_bed_file)) {
  ctcf_fallback <- file.path(out_dir, basename(ctcf_bed_file))
  if (file.exists(ctcf_fallback)) ctcf_bed_file <- ctcf_fallback
}
if (!file.exists(rad21_bed_file)) {
  rad21_fallback <- file.path(out_dir, basename(rad21_bed_file))
  if (file.exists(rad21_fallback)) rad21_bed_file <- rad21_fallback
}
if (!file.exists(brd4_bw_file)) {
  brd4_fallback <- file.path(out_dir, basename(brd4_bw_file))
  if (file.exists(brd4_fallback)) brd4_bw_file <- brd4_fallback
}

need_file(hic_file, "Hi-C .hic file")
need_file(loops_file, "Loop BEDPE file")
need_file(tad_file, "TAD/domain BEDPE file")
need_file(h3k27ac_bw_file, "H3K27ac BigWig file")
need_file(brd4_bw_file, "BRD4 BigWig file")
need_file(ctcf_bed_file, "CTCF BED file")
need_file(rad21_bed_file, "RAD21 BED file")
need_file(gtf_file, "GTF file")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading GTF: ", gtf_file)
gtf <- rtracklayer::import(gtf_file)
meta <- as.data.frame(mcols(gtf))

name_col <- intersect(c("gene_name", "gene", "Name", "gene_id", "ID"), names(meta))[1]
if (is.na(name_col)) stop("Could not find a gene-name column in the GTF metadata.", call. = FALSE)

biotype_col <- intersect(c("gene_biotype", "gene_type", "biotype", "transcript_biotype", "transcript_type"), names(meta))[1]

gtf$plot_gene_name <- as.character(meta[[name_col]])
gtf$plot_type <- if ("type" %in% names(meta)) as.character(meta$type) else NA_character_
gtf$plot_gene_biotype <- if (!is.na(biotype_col)) as.character(meta[[biotype_col]]) else NA_character_

cxcl12_idx <- which(gtf$plot_gene_name == gene_name & (gtf$plot_type == "gene" | is.na(gtf$plot_type)))
if (length(cxcl12_idx) == 0) cxcl12_idx <- which(gtf$plot_gene_name == gene_name)
if (length(cxcl12_idx) == 0) stop("Could not find ", gene_name, " in the GTF file.", call. = FALSE)

cxcl12_gr <- gtf[cxcl12_idx]
plot_chrom <- as.character(seqnames(cxcl12_gr)[1])
gene_start <- min(start(cxcl12_gr))
gene_end <- max(end(cxcl12_gr))
gene_center <- round((gene_start + gene_end) / 2)

message(gene_name, " from GTF: ", plot_chrom, ":", gene_start, "-", gene_end)

message("Reading TAD/domain file: ", tad_file)
tads_all <- read_bedpe(tad_file)
selected_tad <- select_cxcl12_tad(
  tads_all,
  plot_chrom,
  gene_start,
  gene_end,
  choice = cxcl12_tad_choice,
  min_width_bp = large_domain_min_width_bp,
  max_width_bp = large_domain_max_width_bp,
  target_width_bp = large_domain_target_width_bp
)

# Save all candidate domains containing the CXCL12 center, so you can check which domain was selected.
cxcl12_domain_candidates <- tads_all[
  strip_chr(tads_all$chr1) == strip_chr(plot_chrom) &
    strip_chr(tads_all$chr2) == strip_chr(plot_chrom) &
    is.finite(tads_all$x1) & is.finite(tads_all$x2) &
    tads_all$x1 <= gene_center & tads_all$x2 >= gene_center,
]
if (nrow(cxcl12_domain_candidates) > 0) {
  cxcl12_domain_candidates$width <- cxcl12_domain_candidates$x2 - cxcl12_domain_candidates$x1
  cxcl12_domain_candidates$mid <- (cxcl12_domain_candidates$x1 + cxcl12_domain_candidates$x2) / 2
  cxcl12_domain_candidates$center_distance_bp <- abs(cxcl12_domain_candidates$mid - gene_center)
  cxcl12_domain_candidates$target_width_distance_bp <- abs(cxcl12_domain_candidates$width - large_domain_target_width_bp)
  cxcl12_domain_candidates <- cxcl12_domain_candidates[
    order(cxcl12_domain_candidates$center_distance_bp, cxcl12_domain_candidates$target_width_distance_bp, -cxcl12_domain_candidates$width),
  ]
  utils::write.table(cxcl12_domain_candidates, out_domain_candidates_table, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Candidate CXCL12 domains written: ", out_domain_candidates_table)
}

if (isTRUE(use_manual_window)) {
  region_start <- floor(manual_region_start / resolution) * resolution
  region_end <- ceiling(manual_region_end / resolution) * resolution
  message("Using manual window: ", plot_chrom, ":", region_start, "-", region_end)
} else if (isTRUE(use_cxcl12_tad_window) && !is.null(selected_tad) && nrow(selected_tad) > 0) {
  region_start <- max(0, floor((selected_tad$x1[1] - tad_padding_bp) / resolution) * resolution)
  region_end <- ceiling((selected_tad$x2[1] + tad_padding_bp) / resolution) * resolution
  message("Using CXCL12 TAD/domain window from ", basename(tad_file), " (choice=", cxcl12_tad_choice, "): ", plot_chrom, ":", region_start, "-", region_end)
} else {
  region_start <- round(gene_center - left_flank_bp)
  region_end <- round(gene_center + right_flank_bp)
  message(
    "Using asymmetric CXCL12 window: ", plot_chrom, ":", region_start, "-", region_end,
    " (left flank=", left_flank_bp / 1000, " kb; right flank=", right_flank_bp / 1000, " kb; total=",
    (left_flank_bp + right_flank_bp) / 1000, " kb)"
  )
}

# Save the selected TAD/domain for record-keeping.
if (!is.null(selected_tad) && nrow(selected_tad) > 0) {
  selected_tad_out <- selected_tad
  selected_tad_out$selected_for_window <- TRUE
  utils::write.table(selected_tad_out, out_selected_tad_table, sep = "\t", quote = FALSE, row.names = FALSE)
}

# ============================================================
# 5) EXTRACT 5 KB Hi-C CONTACTS FROM THE .hic FILE
# ============================================================

if (!(resolution %in% strawr::readHicBpResolutions(hic_file))) {
  stop("Resolution ", resolution, " is not available in this .hic file.", call. = FALSE)
}

allowed_matrix_types <- c("observed", "oe", "expected")
if (!(matrix_type %in% allowed_matrix_types)) {
  stop(
    "matrix_type must be one of observed/oe/expected for strawr. ",
    "Use norm_method = 'KR' for a balanced observed map, or norm_method = 'NONE' for unbalanced observed counts.",
    call. = FALSE
  )
}

available_norms <- unique(c("NONE", strawr::readHicNormTypes(hic_file)))
if (!(norm_method %in% available_norms)) {
  stop(
    "Normalization '", norm_method, "' is not available. Available normalizations: ",
    paste(available_norms, collapse = ", "),
    call. = FALSE
  )
}

hic_chrom <- choose_hic_chromosome(hic_file, plot_chrom)
hic_region <- paste(hic_chrom, region_start, region_end, sep = ":")

message("Extracting Hi-C contacts: ", hic_region, ", resolution=", resolution, ", norm=", norm_method)
hic <- strawr::straw(norm_method, hic_file, hic_region, hic_region, "BP", resolution, matrix = matrix_type)
names(hic)[1:3] <- c("x", "y", "count")
hic <- hic[hic$x <= hic$y & hic$x >= region_start & hic$y <= region_end & is.finite(hic$count), ]
if (nrow(hic) == 0) stop("No Hi-C contacts returned for this region.", call. = FALSE)

# Balance QC: the main heatmap remains NONE/observed. This only checks whether bin marginals
# are flat and, when available, adds KR as a reference balanced matrix.
if (isTRUE(make_balance_qc)) {
  balance_qc_list <- list(
    compute_hic_balance_qc(
      hic_sparse = hic,
      region_start = region_start,
      region_end = region_end,
      bin_size = resolution,
      norm_label = paste(norm_method, matrix_type),
      chrom_label = plot_chrom
    )
  )
  
  if (isTRUE(include_KR_reference_in_balance_QC) &&
      balance_reference_norm %in% available_norms &&
      balance_reference_norm != norm_method) {
    message("Extracting ", balance_reference_norm, " observed contacts only for balance QC/reference. Main heatmap remains ", norm_method, ".")
    
    hic_ref <- tryCatch(
      strawr::straw(balance_reference_norm, hic_file, hic_region, hic_region, "BP", resolution, matrix = "observed"),
      error = function(e) {
        warning("Could not extract ", balance_reference_norm, " reference for balance QC: ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    
    if (!is.null(hic_ref) && nrow(hic_ref) > 0) {
      names(hic_ref)[1:3] <- c("x", "y", "count")
      hic_ref <- hic_ref[hic_ref$x <= hic_ref$y & hic_ref$x >= region_start & hic_ref$y <= region_end & is.finite(hic_ref$count), ]
      if (nrow(hic_ref) > 0) {
        balance_qc_list[[length(balance_qc_list) + 1]] <- compute_hic_balance_qc(
          hic_sparse = hic_ref,
          region_start = region_start,
          region_end = region_end,
          bin_size = resolution,
          norm_label = paste(balance_reference_norm, "observed reference"),
          chrom_label = plot_chrom
        )
      }
    }
  } else if (isTRUE(include_KR_reference_in_balance_QC) && !(balance_reference_norm %in% available_norms)) {
    message(balance_reference_norm, " is not available in this .hic file, so only NONE balance QC will be written.")
  }
  
  balance_qc <- do.call(rbind, balance_qc_list)
  balance_summary <- summarize_hic_balance_qc(balance_qc)
  
  utils::write.table(balance_qc, out_balance_table, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(balance_summary, out_balance_summary, sep = "\t", quote = FALSE, row.names = FALSE)
  save_hic_balance_qc_plot(balance_qc, balance_summary, out_balance_tiff, dpi = fig_dpi)
  
  message("Hi-C balance QC table: ", out_balance_table)
  message("Hi-C balance QC summary: ", out_balance_summary)
}

hic$display_count <- make_hic_display_values(
  hic,
  bin_size = resolution,
  transform = heatmap_transform,
  cap_quantile = heatmap_cap_quantile,
  ignore_near_diagonal_bins = heatmap_cap_ignore_near_diagonal_bins
)

hic_display_cap <- attr(hic$display_count, "display_cap")
hic_display_max <- max(hic$display_count, na.rm = TRUE)
# Keep the colorbar but remove the verbose scale title from the figure.
# The displayed values are still log1p-transformed NONE observed counts with the configured display cap.
hic_fill_legend_title <- NULL
message(
  "Hi-C color display: transform=", heatmap_transform,
  ", cap_quantile=", heatmap_cap_quantile,
  ", display_cap=", signif(hic_display_cap, 4)
)

hic_poly <- hic_to_triangles(hic, resolution)
heatmap_ymax <- (region_end - region_start) / 2

# ============================================================
# 6) READ LOOPS AND TADs/DOMAINS
# ============================================================

loops <- read_bedpe(loops_file)
loops <- loops[
  strip_chr(loops$chr1) == strip_chr(hic_chrom) &
    strip_chr(loops$chr2) == strip_chr(hic_chrom) &
    loops$x2 >= region_start & loops$x1 <= region_end &
    loops$y2 >= region_start & loops$y1 <= region_end,
]
if (nrow(loops) > 0) {
  loops$plot_loop_id <- seq_len(nrow(loops))
  cx_loop_start <- gene_start - cxcl12_loop_anchor_padding_bp
  cx_loop_end <- gene_end + cxcl12_loop_anchor_padding_bp
  loops$is_cxcl12_loop <- (
    loops$x2 >= cx_loop_start & loops$x1 <= cx_loop_end
  ) | (
    loops$y2 >= cx_loop_start & loops$y1 <= cx_loop_end
  )
} else {
  loops$is_cxcl12_loop <- logical()
}
message("Loops in window: ", nrow(loops))
message("Loops with an anchor overlapping CXCL12 +/- ", cxcl12_loop_anchor_padding_bp, " bp: ", sum(loops$is_cxcl12_loop, na.rm = TRUE))

loop_segments <- loops_to_dotted_segments(loops)
loop_arches <- loops_to_arches(loops)
loop_anchor_guides <- loops_to_anchor_guides(loops, region_start, region_end)

if (isTRUE(plot_only_selected_tad_outline) && !is.null(selected_tad) && nrow(selected_tad) > 0) {
  tads <- selected_tad
} else {
  tads <- tads_all[
    strip_chr(tads_all$chr1) == strip_chr(hic_chrom) &
      strip_chr(tads_all$chr2) == strip_chr(hic_chrom) &
      tads_all$x2 >= region_start & tads_all$x1 <= region_end,
  ]
}
message("TAD/domain outlines plotted: ", nrow(tads))

tad_triangles <- tads_to_triangles(tads, region_start, region_end)

# ============================================================
# 7) IMPORT CTCF BED, RAD21 BED, H3K27ac/BRD4 BIGWIGS, AND BUILD GENE TRACK
# ============================================================

message("Importing CTCF BED: ", ctcf_bed_file)
ctcf_peaks <- read_peak_bed(
  ctcf_bed_file,
  plot_chrom,
  region_start,
  region_end,
  floor_zero = ctcf_floor_zero,
  cap_quantile = ctcf_cap_quantile,
  peak_prefix = "CTCF_peak",
  min_plot_width_bp = ctcf_min_plot_width_bp
)
ctcf_peaks <- make_peak_display_scores(ctcf_peaks, peak_display_transform, peak_display_rescale)
message("CTCF peaks in window: ", nrow(ctcf_peaks))
if (nrow(ctcf_peaks) == 0) {
  warning("No CTCF BED peaks were found in the plotted window. Check that the CTCF BED file is hg38 and uses the expected chromosome naming.")
}

message("Importing RAD21 BED: ", rad21_bed_file)
rad21_peaks <- read_peak_bed(
  rad21_bed_file,
  plot_chrom,
  region_start,
  region_end,
  floor_zero = rad21_floor_zero,
  cap_quantile = rad21_cap_quantile,
  peak_prefix = "RAD21_peak",
  min_plot_width_bp = rad21_min_plot_width_bp
)
rad21_peaks <- make_peak_display_scores(rad21_peaks, peak_display_transform, peak_display_rescale)
message("RAD21 peaks in window: ", nrow(rad21_peaks))
if (nrow(rad21_peaks) == 0) {
  warning("No RAD21 BED peaks were found in the plotted window. Check that RAD21.liver.bed is hg38 and uses the expected chromosome naming.")
}

message("Importing H3K27ac BigWig: ", h3k27ac_bw_file)
h3k27ac_signal <- bigwig_to_track(
  h3k27ac_bw_file,
  plot_chrom,
  region_start,
  region_end,
  floor_zero = h3k27ac_floor_zero,
  cap_quantile = h3k27ac_cap_quantile
)

message("Importing BRD4 BigWig: ", brd4_bw_file)
brd4_signal <- bigwig_to_track(
  brd4_bw_file,
  plot_chrom,
  region_start,
  region_end,
  floor_zero = brd4_floor_zero,
  cap_quantile = brd4_cap_quantile
)

cxcl12_gene_names_to_plot <- NULL
cxcl12_loop_gene_pairs <- data.frame()

if (isTRUE(show_only_cxcl12_loop_connected_genes)) {
  cxcl12_connected <- find_cxcl12_loop_connected_genes(
    gtf = gtf,
    loops = loops,
    gene_name = gene_name,
    plot_chrom = plot_chrom,
    gene_start = gene_start,
    gene_end = gene_end,
    region_start = region_start,
    region_end = region_end,
    cxcl12_padding_bp = cxcl12_loop_anchor_padding_bp,
    partner_padding_bp = partner_gene_anchor_padding_bp
  )
  
  cxcl12_gene_names_to_plot <- unique(c(cxcl12_connected$gene_names, extra_gene_bodies_to_keep))
  cxcl12_loop_gene_pairs <- cxcl12_connected$loop_gene_pairs
  
  connected_gene_summary <- data.frame(
    gene_name = cxcl12_gene_names_to_plot,
    is_CXCL12 = cxcl12_gene_names_to_plot == gene_name,
    stringsAsFactors = FALSE
  )
  
  utils::write.table(connected_gene_summary, out_cxcl12_connected_gene_table, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(cxcl12_loop_gene_pairs, out_cxcl12_loop_gene_pairs_table, sep = "\t", quote = FALSE, row.names = FALSE)
  
  message("Gene bodies restricted to CXCL12-loop-connected genes plus forced TMEM bodies. Genes kept: ", paste(cxcl12_gene_names_to_plot, collapse = ", "))
  message("CXCL12 loop-gene pair table: ", out_cxcl12_loop_gene_pairs_table)
}

gene_track <- make_gene_track_data(
  gtf,
  gene_name,
  plot_chrom,
  region_start,
  region_end,
  gene_names_to_keep = cxcl12_gene_names_to_plot
)
gene_rows <- gene_track$gene_rows
exon_rows <- gene_track$exon_rows

# ------------------------------------------------------------
# Gene-body display intervals
# ------------------------------------------------------------
# Small genes or genes clipped by the plot edge can disappear at journal scale.
# The *_plot columns are used only for drawing the schematic gene track.
# Original genomic coordinates remain unchanged in gene_rows/exon_rows and in the TSV output.
expand_interval_for_display <- function(start, end, min_width, region_start, region_end) {
  if (length(start) == 0) {
    return(data.frame(start_plot = numeric(), end_plot = numeric()))
  }
  if (length(min_width) == 1) {
    min_width <- rep(min_width, length(start))
  }
  
  start_vis <- pmax(start, region_start)
  end_vis <- pmin(end, region_end)
  mid <- (start_vis + end_vis) / 2
  current_width <- pmax(0, end_vis - start_vis)
  needs_expand <- is.finite(current_width) & current_width < min_width
  
  start_plot <- start_vis
  end_plot <- end_vis
  
  start_plot[needs_expand] <- mid[needs_expand] - min_width[needs_expand] / 2
  end_plot[needs_expand] <- mid[needs_expand] + min_width[needs_expand] / 2
  
  # If the expanded interval falls outside the plotting window, shift it back in.
  left_over <- region_start - start_plot
  idx <- which(left_over > 0)
  if (length(idx) > 0) {
    start_plot[idx] <- start_plot[idx] + left_over[idx]
    end_plot[idx] <- end_plot[idx] + left_over[idx]
  }
  
  right_over <- end_plot - region_end
  idx <- which(right_over > 0)
  if (length(idx) > 0) {
    start_plot[idx] <- start_plot[idx] - right_over[idx]
    end_plot[idx] <- end_plot[idx] - right_over[idx]
  }
  
  start_plot <- pmax(start_plot, region_start)
  end_plot <- pmin(end_plot, region_end)
  data.frame(start_plot = start_plot, end_plot = end_plot)
}

if (nrow(gene_rows) > 0) {
  gene_rows$display_min_width_bp <- ifelse(
    grepl("^TMEM72", gene_rows$gene_name, ignore.case = TRUE),
    tmem_gene_body_min_plot_width_bp,
    gene_body_min_plot_width_bp
  )
  gene_disp <- expand_interval_for_display(
    gene_rows$start, gene_rows$end,
    gene_rows$display_min_width_bp,
    region_start, region_end
  )
  gene_rows$start_plot <- gene_disp$start_plot
  gene_rows$end_plot <- gene_disp$end_plot
} else {
  gene_rows$start_plot <- numeric()
  gene_rows$end_plot <- numeric()
}

if (nrow(exon_rows) > 0) {
  exon_rows$display_min_width_bp <- ifelse(
    grepl("^TMEM72", exon_rows$gene_name, ignore.case = TRUE),
    max(gene_exon_min_plot_width_bp, 8000),
    gene_exon_min_plot_width_bp
  )
  exon_disp <- expand_interval_for_display(
    exon_rows$start, exon_rows$end,
    exon_rows$display_min_width_bp,
    region_start, region_end
  )
  exon_rows$start_plot <- exon_disp$start_plot
  exon_rows$end_plot <- exon_disp$end_plot
} else {
  exon_rows$start_plot <- numeric()
  exon_rows$end_plot <- numeric()
}

# ------------------------------------------------------------
# Gene-label placement
# ------------------------------------------------------------
# For the final publication figure, gene bodies are kept but only CXCL12 is labeled.
# This removes clutter from dense neighboring genes while preserving their genomic positions.
gene_label_rows <- data.frame(
  gene_name = character(), start = numeric(), end = numeric(), y = numeric(),
  is_target = logical(), label_x_plot = numeric(), label_y = numeric(),
  label_hjust = numeric(), stringsAsFactors = FALSE
)
gene_label_leaders <- data.frame(
  leader_x = numeric(), leader_xend = numeric(), leader_y = numeric(), leader_yend = numeric()
)

if (nrow(gene_rows) > 0) {
  gene_rows$gene_mid <- (gene_rows$start + gene_rows$end) / 2
  
  cx_idx <- which(gene_rows$is_target | gene_rows$gene_name == gene_name)
  if (length(cx_idx) > 0) {
    gene_rows$label_x_plot <- pmin(
      pmax(gene_rows$gene_mid, region_start + gene_label_padding_bp),
      region_end - gene_label_padding_bp
    )
    gene_rows$label_y <- pmax(gene_rows$y - cxcl12_label_below_offset, 0.055)
    gene_rows$label_hjust <- 0.5
    
    gene_label_rows <- gene_rows[cx_idx, c(
      "gene_name", "start", "end", "y", "is_target",
      "label_x_plot", "label_y", "label_hjust"
    ), drop = FALSE]
  }
}

message("Genes plotted after RP/LINC filtering and CXCL12-loop filtering: ", nrow(gene_rows))

utils::write.table(gene_rows, out_gene_table, sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(ctcf_peaks, out_ctcf_table, sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(rad21_peaks, out_rad21_table, sep = "\t", quote = FALSE, row.names = FALSE)

# ============================================================
# 8) MAKE THE FIGURE
# ============================================================

mb_labels <- function(x) sprintf("%.2f", x / 1e6)
x_limits <- c(region_start, region_end)
base_x <- scale_x_continuous(limits = x_limits, labels = mb_labels, expand = expansion(mult = 0))

p_hic <- ggplot() +
  annotate("rect", xmin = gene_start, xmax = gene_end, ymin = 0, ymax = heatmap_ymax,
           fill = cxcl12_highlight_fill, alpha = cxcl12_highlight_alpha) +
  geom_polygon(data = hic_poly, aes(x = xplot, y = yplot, group = id, fill = display_count), color = NA) +
  geom_path(data = tad_triangles, aes(x = x, y = y, group = tad_id),
            color = "grey20", size = tad_line_size) +
  geom_segment(data = loop_segments[!loop_segments$is_cxcl12_loop, , drop = FALSE],
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "black", size = loop_line_size, linetype = "dotted") +
  geom_segment(data = loop_segments[loop_segments$is_cxcl12_loop, , drop = FALSE],
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "black", size = cxcl12_loop_line_size, linetype = "solid") +
  scale_fill_gradientn(
    colors = heatmap_palette,
    limits = c(0, hic_display_max),
    name = hic_fill_legend_title
  ) +
  coord_cartesian(xlim = x_limits, ylim = c(0, heatmap_ymax), expand = FALSE, clip = "on") +
  labs(title = paste0(sample_label, ": ", gene_name)) +
  theme_void(base_size = 8) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold", hjust = 0),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 5.8),
    legend.key.height = grid::unit(0.75, "cm"),
    legend.key.width = grid::unit(0.20, "cm"),
    plot.margin = margin(2, 2, 1, 2)
  )

p_loop <- ggplot() +
  annotate("rect", xmin = gene_start, xmax = gene_end, ymin = -Inf, ymax = Inf,
           fill = cxcl12_highlight_fill, alpha = cxcl12_highlight_alpha) +
  geom_path(data = loop_arches[!loop_arches$is_cxcl12_loop, , drop = FALSE],
            aes(x = x, y = y, group = loop_id),
            color = "black", size = loop_line_size, linetype = "dotted") +
  geom_path(data = loop_arches[loop_arches$is_cxcl12_loop, , drop = FALSE],
            aes(x = x, y = y, group = loop_id),
            color = "black", size = cxcl12_loop_line_size, linetype = "solid") +
  base_x +
  coord_cartesian(ylim = c(0, 1.05), expand = FALSE, clip = "on") +
  labs(y = "Loops") +
  theme_classic(base_size = 8) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_text(size = track_label_size),
    plot.margin = margin(0, 2, 0, 2)
  )

p_ctcf <- make_peak_plot(
  peaks = ctcf_peaks,
  label = "CTCF",
  fill_color = ctcf_peak_color,
  gene_start = gene_start,
  gene_end = gene_end,
  x_limits = x_limits,
  anchor_guides = loop_anchor_guides
)

p_rad21 <- make_peak_plot(
  peaks = rad21_peaks,
  label = "RAD21",
  fill_color = rad21_peak_color,
  gene_start = gene_start,
  gene_end = gene_end,
  x_limits = x_limits,
  anchor_guides = loop_anchor_guides
)

p_h3k27ac <- make_signal_plot(
  signal = h3k27ac_signal,
  label = "H3K27ac",
  fill_color = "#006D2C",
  gene_start = gene_start,
  gene_end = gene_end,
  x_limits = x_limits,
  anchor_guides = loop_anchor_guides
)

p_brd4 <- make_signal_plot(
  signal = brd4_signal,
  label = "BRD4",
  fill_color = "#08519C",
  gene_start = gene_start,
  gene_end = gene_end,
  x_limits = x_limits,
  anchor_guides = loop_anchor_guides
)

p_genes <- ggplot() +
  annotate("rect", xmin = gene_start, xmax = gene_end, ymin = -Inf, ymax = Inf,
           fill = cxcl12_highlight_fill, alpha = cxcl12_highlight_alpha) +
  geom_segment(data = gene_rows, aes(x = start_plot, xend = end_plot, y = y, yend = y),
               size = gene_body_line_size, color = "grey28") +
  geom_rect(data = gene_rows, aes(xmin = start_plot, xmax = end_plot,
                                  ymin = y - gene_body_rect_half_height, ymax = y + gene_body_rect_half_height,
                                  fill = is_target),
            color = NA) +
  geom_rect(data = exon_rows, aes(xmin = start_plot, xmax = end_plot,
                                  ymin = y - gene_exon_rect_half_height, ymax = y + gene_exon_rect_half_height,
                                  fill = is_target),
            color = NA) +
  geom_segment(data = gene_label_leaders,
               aes(x = leader_x, xend = leader_xend, y = leader_y, yend = leader_yend),
               color = gene_label_leader_color, size = gene_label_leader_size) +
  geom_label(data = gene_label_rows, aes(x = label_x_plot, y = label_y, label = gene_name,
                                         fontface = ifelse(is_target, "bold", "plain"),
                                         hjust = label_hjust),
             size = gene_label_size, label.size = 0, fill = gene_label_box_fill,
             alpha = gene_label_box_alpha, label.padding = grid::unit(0.035, "lines"),
             check_overlap = FALSE) +
  scale_fill_manual(values = c("FALSE" = "grey65", "TRUE" = "black"), guide = "none") +
  base_x +
  coord_cartesian(ylim = c(0, 1), expand = FALSE, clip = "off") +
  labs(x = paste0(plot_chrom, " position (Mb)"), y = "Genes") +
  theme_classic(base_size = 8) +
  theme(
    axis.title.y = element_text(size = track_label_size),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_line(color = "black", size = 0.30),
    plot.margin = margin(1.5, 8, 2, 2)
  )

track_gap_loop_ctcf <- patchwork::plot_spacer()
track_gap1 <- patchwork::plot_spacer()
track_gap2 <- patchwork::plot_spacer()
track_gap3 <- patchwork::plot_spacer()
track_gap4 <- patchwork::plot_spacer()

final_plot <- p_hic / p_loop / track_gap_loop_ctcf / p_ctcf / track_gap1 / p_rad21 / track_gap2 / p_h3k27ac / track_gap3 / p_brd4 / track_gap4 / p_genes +
  patchwork::plot_layout(
    heights = c(
      3.88,
      loop_panel_layout_height,
      loop_to_ctcf_spacer_height,
      0.54,
      regulatory_track_gap_height,
      0.54,
      regulatory_track_gap_height,
      0.56,
      regulatory_track_gap_height,
      0.56,
      gene_h3k27ac_spacer_height,
      gene_panel_layout_height
    ),
    guides = "collect"
  )

# ============================================================
# 9) EXPORT FIGURE
# ============================================================

message("Writing: ", out_tiff)

if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_tiff(
    filename = out_tiff,
    width = fig_width_in,
    height = fig_height_in,
    units = "in",
    res = fig_dpi,
    compression = "lzw",
    background = fig_background
  )
} else {
  grDevices::tiff(
    filename = out_tiff,
    width = fig_width_in,
    height = fig_height_in,
    units = "in",
    res = fig_dpi,
    compression = "lzw",
    bg = fig_background,
    type = "cairo"
  )
}

print(final_plot)
grDevices::dev.off()

utils::capture.output(sessionInfo(), file = out_session_info)

message("Done: ", normalizePath(out_tiff))
message("Gene table: ", normalizePath(out_gene_table))
message("CTCF peak table: ", normalizePath(out_ctcf_table))
message("RAD21 peak table: ", normalizePath(out_rad21_table))
if (file.exists(out_cxcl12_connected_gene_table)) message("CXCL12-connected gene table: ", normalizePath(out_cxcl12_connected_gene_table))
if (file.exists(out_cxcl12_loop_gene_pairs_table)) message("CXCL12 loop-gene pair table: ", normalizePath(out_cxcl12_loop_gene_pairs_table))
if (file.exists(out_selected_tad_table)) message("Selected CXCL12 domain table: ", normalizePath(out_selected_tad_table))
if (file.exists(out_domain_candidates_table)) message("Candidate CXCL12 domain table: ", normalizePath(out_domain_candidates_table))
message("Session info: ", normalizePath(out_session_info))
if (isTRUE(make_balance_qc)) {
  message("Hi-C balance QC bins: ", normalizePath(out_balance_table))
  message("Hi-C balance QC summary: ", normalizePath(out_balance_summary))
  message("Hi-C balance QC TIFF: ", normalizePath(out_balance_tiff))
}
