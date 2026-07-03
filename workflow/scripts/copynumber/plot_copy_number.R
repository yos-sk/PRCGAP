#!/usr/bin/env Rscript
# Copy Number Visualization Script (Genome-wide view)
# This script creates copy number plots for two haplotypes with CBS segmentation
# Using a single-row genome-wide layout for efficient space usage

library(scales)
library(dplyr)
library(optparse)
library(tidyverse)
library(patchwork)
library(slider)
library(viridis)
library(ggnewscale)

# Peak-comb calibration helpers (detect_peaks, peak_comb_fit). estimate_ploidy.R
# lives alongside this script; resolve this script's own directory so it can be
# sourced regardless of the current working directory (the workflow invokes
# plot_copy_number.R by absolute path, not from the scripts/ directory).
get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  return(getwd())
}
source(file.path(get_script_dir(), "estimate_ploidy.R"))

# Constants
MIN_COVERAGE_THRESHOLD <- 100000
MIN_TUMOR_COVERAGE <- 10000
MIN_SEGMENT_SIZE <- 100
DEFAULT_PLOIDY <- 1
DEFAULT_MAX_COPY <- 6
DEFAULT_BINWIDTH <- 0.05
SCALING_FACTOR <- 1e6  # Scale to Mb for genome-wide view
SEGMENT_COLOR <- "#B63679"  # Bright red for CBS segments
ANNOTATION_COLOR <- "#D3D3D3"
POINT_ALPHA <- 0.3
POINT_SIZE <- 0.15

# CHM13 sex chromosome lengths (bp)
CHM13_CHRX_LENGTH_BP <- 154259566
CHM13_CHRY_LENGTH_BP <- 62460029

# Parse command-line arguments
parse_arguments <- function() {
  option_list <- list(
    make_option(c("-i", "--input1"),
                type = "character",
                help = "Input haplotype 1 copy-number text file"),
    make_option(c("-j", "--input2"),
                type = "character",
                help = "Input haplotype 2 copy-number text file"),
    make_option(c("-k", "--cbs1"),
                type = "character",
                help = "Input haplotype 1 CBS text file"),
    make_option(c("-l", "--cbs2"),
                type = "character",
                help = "Input haplotype 2 CBS text file"),
    make_option(c("-o", "--output"),
                type = "character",
                help = "Output PDF file"),
    make_option(c("-s", "--sample"),
                type = "character",
                help = "Sample name"),
    make_option(c("-p", "--ploidy1"),
                type = "integer",
                default = DEFAULT_PLOIDY,
                help = "Tumor hap1 ploidy [default: %default]"),
    make_option(c("-t", "--ploidy2"),
                type = "integer",
                default = DEFAULT_PLOIDY,
                help = "Tumor hap2 ploidy [default: %default]"),
    make_option(c("-q", "--table1"),
                type = "character",
                help = "Input reference-contig table (hap1)"),
    make_option(c("-r", "--table2"),
                type = "character",
                help = "Input reference-contig table (hap2)"),
    make_option(c("-a", "--annotation"),
                type = "character",
                help = "Annotation file"),
    make_option(c("-c", "--max_copy"),
                type = "numeric",
                default = DEFAULT_MAX_COPY,
                help = "Max copy number for visualization [default: %default]"),
    make_option(c("-w", "--binwidth"),
                type = "double",
                default = DEFAULT_BINWIDTH,
                help = "Bin width for getting mode value of depth ratio [default: %default]"),

    make_option(c("--hap1_label"),
                type = "character",
                default = "Haplotype1",
                help = "haplotype1 label"),

    make_option(c("--hap2_label"),
                type = "character",
                default = "Haplotype2",
                help = "haplotype2 label"),

    make_option(c("--plot-sex-chrom"),
                type = "character",
                default = "TRUE",
                help = paste("Force chrX/chrY onto the plot even when absent from",
                             "the data (estimating length from CHM13). Set to",
                             "FALSE to show sex chromosomes only when present",
                             "[default: %default]"))
  )

  opt <- parse_args(OptionParser(option_list = option_list))

  # Validate required arguments
  required_args <- c("input1", "input2", "cbs1", "cbs2", "output",
                     "sample", "table1", "table2", "annotation")
  missing_args <- required_args[sapply(required_args, function(x) is.null(opt[[x]]))]

  if (length(missing_args) > 0) {
    stop(paste("Error: Missing required arguments:", paste(missing_args, collapse = ", ")))
  }

  return(opt)
}

# Calculate chromosome offsets for genome-wide coordinates
calculate_chr_offsets <- function(chr_lengths) {
  chr_order <- c(paste0("chr", 1:22), "chrX", "chrY")
  # Only use chromosomes that exist in chr_lengths
  chr_order <- chr_order[chr_order %in% names(chr_lengths)]
  chr_lengths_ordered <- as.numeric(chr_lengths[chr_order])
  names(chr_lengths_ordered) <- chr_order

  offsets <- c(0, cumsum(as.numeric(chr_lengths_ordered[-length(chr_lengths_ordered)])))
  names(offsets) <- chr_order

  # Calculate midpoints for x-axis labels
  midpoints <- offsets + chr_lengths_ordered / 2

  # Calculate boundaries for chromosome separator lines
  boundaries <- cumsum(as.numeric(chr_lengths_ordered))
  names(boundaries) <- chr_order

  list(
    offsets = offsets,
    midpoints = midpoints,
    boundaries = boundaries,
    total_length = sum(as.numeric(chr_lengths_ordered)),
    chr_order = chr_order,
    chr_lengths = chr_lengths_ordered
  )
}

# Calculate max chromosome lengths from copynumber data
calculate_chr_lengths_from_copynumber <- function(copynumber1, copynumber2) {
  # Get max position for each chromosome from both haplotypes
  max_per_chr1 <- copynumber1 %>%
    group_by(chr) %>%
    summarize(max_pos = max(pos, na.rm = TRUE), .groups = "drop")

  max_per_chr2 <- copynumber2 %>%
    group_by(chr) %>%
    summarize(max_pos = max(pos, na.rm = TRUE), .groups = "drop")

  # Combine and take max from both haplotypes
  combined <- bind_rows(max_per_chr1, max_per_chr2) %>%
    group_by(chr) %>%
    summarize(max_pos = max(max_pos, na.rm = TRUE), .groups = "drop")

  # Convert to named vector
  chr_lengths <- setNames(combined$max_pos, combined$chr)

  return(chr_lengths)
}

# Ensure sex chromosomes are included in chr_lengths.
# If a sex chromosome is missing from both haplotypes, estimate its bin length
# using CHM13 reference lengths and the median bp_per_bin from autosomes.
ensure_sex_chromosomes <- function(chr_lengths, ref_table1, ref_table2) {
  # Get bp lengths from ref_tables
  ref_max1 <- ref_table1 %>% group_by(chr) %>% summarize(max_bp = max(ref_end), .groups = "drop")
  ref_max2 <- ref_table2 %>% group_by(chr) %>% summarize(max_bp = max(ref_end), .groups = "drop")
  ref_max <- bind_rows(ref_max1, ref_max2) %>%
    group_by(chr) %>%
    summarize(max_bp = max(max_bp), .groups = "drop")

  # Compute bp_per_bin for autosomes that have both ref_table and copynumber data
  autosomes <- names(chr_lengths)[grepl("^chr[0-9]+$", names(chr_lengths))]
  bp_per_bin_values <- c()
  for (chr_name in autosomes) {
    if (chr_name %in% ref_max$chr) {
      bp <- ref_max$max_bp[ref_max$chr == chr_name]
      bins <- chr_lengths[chr_name]
      bp_per_bin_values <- c(bp_per_bin_values, bp / bins)
    }
  }
  median_bp_per_bin <- median(bp_per_bin_values)
  cat("Median bp_per_bin:", median_bp_per_bin, "\n")

  # Add chrX if missing
  if (!("chrX" %in% names(chr_lengths))) {
    chr_lengths["chrX"] <- CHM13_CHRX_LENGTH_BP / median_bp_per_bin
    cat("chrX not in data, using CHM13 estimate:", chr_lengths["chrX"], "bins\n")
  }
  # Add chrY if missing
  if (!("chrY" %in% names(chr_lengths))) {
    chr_lengths["chrY"] <- CHM13_CHRY_LENGTH_BP / median_bp_per_bin
    cat("chrY not in data, using CHM13 estimate:", chr_lengths["chrY"], "bins\n")
  }

  return(chr_lengths)
}

# Define custom theme for genome-wide plots
create_plot_theme <- function() {
  theme_bw(base_family = "Helvetica") %+replace%
    theme(
      title = element_text(size = 7),
      plot.title = element_text(hjust = 0.5, size = 7, margin = margin(t = -3, b = 0)),
      panel.border = element_rect(colour = "grey50", fill = NA, linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray80", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      axis.text = element_text(size = 6),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(linewidth = 0.3),
      axis.title = element_text(size = 8),
      axis.title.x = element_blank(),
      legend.position = "none",
      plot.margin = margin(1, 5, 1, 5),
      complete = TRUE
    )
}

# Load and normalize copy number data
# File format: chr, contig, pos_bin, tumor_depth, normal_depth
# pos_bin is a bin index (not actual Mb), where each bin is ~50 kb
load_and_normalize_copynumber <- function(input_file, tumor_ploidy, binwidth) {
  # Load data
  data <- read.table(input_file, sep = "\t", header = FALSE, comment.char = "")
  colnames(data) <- c("chr", "contig", "pos_bin", "tumor_depth", "normal_depth")

  # Use bin position directly (same coordinate system as CBS file)
  data$pos <- data$pos_bin

  # Calculate normalization factors
  normal_ploidy <- 1
  normal_count <- sum(data[data$normal_depth > MIN_COVERAGE_THRESHOLD, ]$normal_depth)
  tumor_count <- sum(data[data$normal_depth > MIN_COVERAGE_THRESHOLD, ]$tumor_depth)
  normal_nrow <- nrow(data)
  tumor_nrow <- nrow(data[data$tumor_depth > MIN_TUMOR_COVERAGE, ])

  # Ploidy-agnostic depth ratio (R0, tumor_ploidy = 1)
  data$depth_ratio <- data$tumor_depth / data$normal_depth *
    normal_count / tumor_count *
    tumor_nrow / (normal_nrow * normal_ploidy)

  # Filter data
  filtered_data <- data %>%
    filter(!is.na(depth_ratio) & is.finite(depth_ratio),
           normal_depth > MIN_COVERAGE_THRESHOLD)

  # Calibrate so the dominant copy-number state sits at the integer ploidy,
  # using the peak comb (consistent with cbs.R / estimate_ploidy).
  filtered_data <- calibrate_depth_ratio(filtered_data, tumor_ploidy, binwidth)

  return(filtered_data)
}

# Calculate mode value from histogram
calculate_mode_value <- function(data, binwidth) {
  positive_data <- data[data$depth_ratio > 0, ]

  hist_plot <- ggplot(positive_data, aes(x = depth_ratio)) +
    geom_histogram(binwidth = binwidth)
  hist_data <- ggplot_build(hist_plot)$data[[1]]

  # Get mode value only from bins where x > 0
  hist_data_positive <- hist_data[hist_data$x > 0, ]
  mode_value <- hist_data_positive$x[which.max(hist_data_positive$count)]

  return(mode_value)
}

# Calibrate the ploidy-agnostic depth ratio (R0) so the dominant copy-number
# state sits at the integer ploidy, matching cbs.R: the dominant state is the
# tallest KDE peak (detect_peaks), and mu = tallest_peak / ploidy places it
# exactly at the integer ploidy. Falls back to the histogram mode when no peak
# is found.
calibrate_depth_ratio <- function(filtered_data, ploidy, binwidth) {
  peaks <- detect_peaks(filtered_data$depth_ratio)
  if (nrow(peaks) > 0) {
    mu <- peaks$x[which.max(peaks$y)] / ploidy
  } else {
    mu <- calculate_mode_value(filtered_data, binwidth) / ploidy
  }
  filtered_data$depth_ratio <- filtered_data$depth_ratio / mu
  return(filtered_data)
}

# Load CBS segmentation data with genome-wide coordinates
# CBS file format: sample, chr_num, start_bin, end_bin, num_markers, seg_mean
# start_bin and end_bin use the same coordinate system as copynumber data
load_cbs_data <- function(cbs_file, chr_info) {
  # Define chromosome order
  chr_levels <- c(paste0("chr", 1:22), "chrX", "chrY")

  CNA_data <- read.table(cbs_file, sep = "\t", header = FALSE, comment.char = "")
  colnames(CNA_data) <- c("sample", "chr_num", "start_bin", "end_bin", "num_markers", "seg_mean")
  # Add "chr" prefix (positions are already in bin coordinates)
  CNA_data <- CNA_data %>%
    mutate(
      chr = paste0("chr", chr_num),
      start = start_bin,
      end = end_bin
    )

  # Convert chr to factor with proper order and add genome-wide coordinates
  CNA_data <- CNA_data %>%
    mutate(chr = factor(chr, levels = chr_levels)) %>%
    filter(chr %in% chr_info$chr_order) %>%
    mutate(
      genome_start = start + chr_info$offsets[as.character(chr)],
      genome_end = end + chr_info$offsets[as.character(chr)]
    )

  return(CNA_data)
}

# Transform annotation data to chromosome coordinates using reference table
# Then scale to bin coordinates to match copynumber data
transform_annotation_to_chr_coords <- function(annotation_file, ref_table, chr_info, copynumber_data) {
  # Define chromosome order
  chr_levels <- c(paste0("chr", 1:22), "chrX", "chrY")

  contigs <- ref_table$name

  # Calculate scaling factors: for each chromosome, find bp_per_bin
  # Using copynumber data max position and reference table max position
  copynumber_max <- copynumber_data %>%
    group_by(chr) %>%
    summarize(max_bin = max(pos, na.rm = TRUE), .groups = "drop")

  ref_max <- ref_table %>%
    group_by(chr) %>%
    summarize(max_bp = max(ref_end, na.rm = TRUE), .groups = "drop")

  scale_factors <- copynumber_max %>%
    inner_join(ref_max, by = "chr") %>%
    mutate(bp_per_bin = max_bp / max_bin)

  # Read annotation file and select only first 3 columns
  annotation_raw <- read.table(gzfile(annotation_file), sep = "\t", header = FALSE, comment.char = "")
  annotation_data <- annotation_raw[, 1:3]
  colnames(annotation_data) <- c("contig", "start", "end")

  # Filter to contigs in reference table and join
  transformed <- annotation_data %>%
    filter(contig %in% contigs) %>%
    inner_join(ref_table, by = c("contig" = "name"), relationship = "many-to-many") %>%
    filter(start >= contig_start & end <= contig_end) %>%
    mutate(
      # Transform coordinates based on strand (to chromosome bp)
      chr_start_bp = ifelse(
        strand == "+",
        ref_start + (start - contig_start),
        ref_end - (end - contig_start)
      ),
      chr_end_bp = ifelse(
        strand == "+",
        ref_start + (end - contig_start),
        ref_end - (start - contig_start)
      ),
      chr = factor(chr, levels = chr_levels)
    ) %>%
    filter(chr %in% chr_info$chr_order) %>%
    # Join with scale factors to convert bp to bin coordinates
    left_join(scale_factors, by = "chr") %>%
    mutate(
      chr_start = chr_start_bp / bp_per_bin,
      chr_end = chr_end_bp / bp_per_bin,
      genome_start = chr_start + chr_info$offsets[as.character(chr)],
      genome_end = chr_end + chr_info$offsets[as.character(chr)]
    )

  cat("  Final annotation entries:", nrow(transformed), "\n")

  return(transformed)
}

# Load and process reference table
load_reference_table <- function(table_file) {
  ref_table <- read.table(table_file, sep = "\t", header = FALSE, comment.char = "")
  colnames(ref_table) <- c("name", "contig_start", "contig_end",
                           "strand", "chr", "ref_start", "ref_end")

  return(ref_table)
}

# Add genome-wide coordinates and rolling mean to copynumber data
# Copynumber data already has chromosome coordinates (chr, pos)
filter_and_transform_copynumber <- function(copynumber_data, chr_info) {
  # Define chromosome order
  chr_levels <- c(paste0("chr", 1:22), "chrX", "chrY")

  transformed <- copynumber_data %>%
    mutate(chr = factor(chr, levels = chr_levels)) %>%
    filter(chr %in% chr_info$chr_order) %>%
    mutate(
      genome_pos = pos + chr_info$offsets[as.character(chr)]
    ) %>%
    # Apply rolling mean grouped by chromosome
    group_by(chr) %>%
    arrange(pos) %>%
    mutate(rollmean = slide_dbl(depth_ratio, mean,
                                .before = 3, .after = 3,
                                .complete = TRUE)) %>%
    ungroup() %>%
    filter(!is.na(rollmean))

  return(transformed)
}


# Create chromosome background rectangles for alternating colors
create_chr_background <- function(chr_info, max_copy) {
  chr_order <- chr_info$chr_order
  n_chr <- length(chr_order)

  data.frame(
    chr = chr_order,
    xmin = c(0, chr_info$boundaries[-n_chr]),
    xmax = chr_info$boundaries,
    ymin = 0,
    ymax = max_copy,
    odd = seq_along(chr_order) %% 2 == 1
  )
}

# Create genome-wide copy number plot for one haplotype
create_copynumber_plot <- function(filtered_data, cbs_data, annotation_data,
                                   chr_info, haplotype_name, max_copy, plot_theme,
                                   show_legend = FALSE) {
  # Create chromosome labels (remove "chr" prefix for cleaner display)
  chr_labels <- gsub("chr", "", chr_info$chr_order)

  # Get unique contigs for color mapping
  n_contigs <- length(unique(filtered_data$contig))

  p <- ggplot() +
    # Chromosome boundary lines (vertical lines at chromosome boundaries)
    geom_vline(xintercept = chr_info$boundaries[-length(chr_info$boundaries)],
               color = "gray40", linewidth = 0.2, linetype = "solid") +
    # Reference line at CN = 1
    geom_hline(yintercept = 1, linetype = "dashed",
               color = "gray50", linewidth = 0.3) +
    # Data points (colored by contig using viridis)
    geom_point(data = filtered_data,
               aes(x = genome_pos, y = pmin(rollmean, max_copy), color = contig),
               size = POINT_SIZE, alpha = POINT_ALPHA) +
    scale_color_viridis(discrete = TRUE, option = "D", begin = 0, end = 1,
                        guide = "none") +
    # CBS segments (bright red, thick)
    geom_segment(data = cbs_data,
                 aes(x = genome_start, xend = genome_end,
                     y = pmin(seg_mean, max_copy),
                     yend = pmin(seg_mean, max_copy)),
                 colour = SEGMENT_COLOR, linewidth = 0.8) +
    # Axis settings
    scale_x_continuous(
      breaks = chr_info$midpoints,
      labels = chr_labels,
      limits = c(0, chr_info$total_length),
      expand = c(0.005, 0)
    ) +
    scale_y_continuous(
      breaks = seq(0, max_copy, by = 1),
      limits = c(0, max_copy),
      expand = c(0.02, 0)
    ) +
    labs(title = haplotype_name, y = "Copy number") +
    plot_theme

  # Add annotation rectangles if available
  if (nrow(annotation_data) > 0) {
    p <- p + geom_rect(data = annotation_data,
                       aes(xmin = genome_start, xmax = genome_end,
                           ymin = 0, ymax = max_copy,
                           fill = "Centromere/Satellite"),
                       alpha = 0.4,
                       inherit.aes = FALSE, linewidth = 0) +
      scale_fill_manual(values = c("Centromere/Satellite" = ANNOTATION_COLOR),
                        name = NULL)
    if (show_legend) {
      p <- p + theme(legend.position = "inside",
                     legend.position.inside = c(1, 1.1),
                     legend.justification = c(1, 0.5),
                     legend.background = element_blank(),
                     legend.key.size = unit(3, "mm"),
                     legend.text = element_text(size = 6),
                     legend.margin = margin(0, 2, 0, 0),
                     plot.margin = margin(8, 5, 1, 5))
    }
  }

  return(p)
}

# Process one haplotype (load data and create plot)
process_haplotype <- function(copynumber_data, cbs_file, ref_table,
                              annotation_file, haplotype_name, max_copy,
                              chr_info, plot_theme, show_legend = FALSE) {
  cat("Processing", haplotype_name, "...\n")

  # Load CBS data with genome-wide coordinates
  cbs_data <- load_cbs_data(cbs_file, chr_info)

  # Add genome-wide coordinates and rolling mean
  filtered_data <- filter_and_transform_copynumber(copynumber_data, chr_info)

  # Transform annotation data to chromosome coordinates (scaled to bin coordinates)
  annotation_data <- transform_annotation_to_chr_coords(annotation_file, ref_table, chr_info, copynumber_data)

  cat("Annotation data loaded:", nrow(annotation_data), "regions\n")
  cat("Filtered data points:", nrow(filtered_data), "\n")
  cat("CBS segments:", nrow(cbs_data), "\n")

  # Create plot
  plot <- create_copynumber_plot(filtered_data, cbs_data, annotation_data,
                                 chr_info, haplotype_name, max_copy, plot_theme,
                                 show_legend = show_legend)

  return(plot)
}

# Save combined plot (vertical stacking for genome-wide view)
save_combined_plot <- function(plot1, plot2, output_file, sample_name) {
  # Set output dimensions for genome-wide view (wider, shorter)
  width_mm <- 175
  height_mm <- 48

  # Combine plots vertically using patchwork with sample title
  combined_plot <- plot1 / plot2 +
    plot_annotation(title = sample_name,
                    theme = theme(plot.title = element_text(
                      family = "Helvetica", size = 8, face = "bold",
                      hjust = 0.5, margin = margin(b = 0))))

  # Save the combined plot
  ggsave(
    filename = output_file,
    plot = combined_plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 300
  )

  cat("Plot saved to:", output_file, "\n")
}

# Main function
main <- function() {
  # Parse arguments
  opt <- parse_arguments()

  cat("Creating copy number plots for sample:", opt$sample, "\n")

  # Load and normalize copy number data first
  cat("Loading copy number data...\n")
  copynumber1 <- load_and_normalize_copynumber(opt$input1, opt$ploidy1, opt$binwidth)
  copynumber2 <- load_and_normalize_copynumber(opt$input2, opt$ploidy2, opt$binwidth)

  cat("Haplotype 1 data points:", nrow(copynumber1), "\n")
  cat("Haplotype 2 data points:", nrow(copynumber2), "\n")

  # Load reference tables for annotation coordinate transformation
  cat("Loading reference tables...\n")
  ref_table1 <- load_reference_table(opt$table1)
  ref_table2 <- load_reference_table(opt$table2)

  # Calculate chromosome lengths from copynumber data (max of both haplotypes)
  chr_lengths <- calculate_chr_lengths_from_copynumber(copynumber1, copynumber2)

  # Force chrX/chrY onto the plot (estimate length from CHM13 if absent) unless
  # --plot-sex-chrom is FALSE, in which case they appear only when present in the
  # data (calculate_chr_offsets drops chromosomes missing from chr_lengths).
  plot_sex_chrom <- !tolower(opt[["plot-sex-chrom"]]) %in% c("false", "f", "no", "n", "0", "none")
  if (plot_sex_chrom) {
    chr_lengths <- ensure_sex_chromosomes(chr_lengths, ref_table1, ref_table2)
  } else {
    cat("plot-sex-chrom=FALSE: not forcing chrX/chrY (shown only if present in data)\n")
  }
  cat("Chromosome lengths calculated\n")

  # Calculate chromosome offsets for genome-wide coordinates
  chr_info <- calculate_chr_offsets(chr_lengths)
  cat("Total genome length:", chr_info$total_length, "\n")
  cat("Chromosomes:", paste(chr_info$chr_order, collapse = ", "), "\n")

  # Create plot theme
  plot_theme <- create_plot_theme()

  # Process haplotype 1
  plot1 <- process_haplotype(
    copynumber_data = copynumber1,
    cbs_file = opt$cbs1,
    ref_table = ref_table1,
    annotation_file = opt$annotation,
    haplotype_name = opt$hap1_label,
    max_copy = opt$max_copy,
    chr_info = chr_info,
    plot_theme = plot_theme,
    show_legend = TRUE
  )

  # Process haplotype 2
  plot2 <- process_haplotype(
    copynumber_data = copynumber2,
    cbs_file = opt$cbs2,
    ref_table = ref_table2,
    annotation_file = opt$annotation,
    haplotype_name = opt$hap2_label,
    max_copy = opt$max_copy,
    chr_info = chr_info,
    plot_theme = plot_theme
  )

  # Save combined plot
  display_name <- ifelse(opt$sample == "HG008T", "HG008", opt$sample)
  save_combined_plot(plot1, plot2, opt$output, display_name)

  cat("Analysis complete!\n")
}

# Run main function
if (!interactive()) {
  main()
}
