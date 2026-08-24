#!/usr/bin/env Rscript
# Circular Binary Segmentation (CBS) analysis for copy number variation detection
# This script performs CBS on copy number data from PRCGAP output

library(DNAcopy)
library(tidyverse)
library(optparse)

# CBS-segment ploidy estimator (estimate_ploidy_halfwin(), segment_dominant_level()).
# estimate_ploidy.R lives alongside this script; resolve this script's own
# directory so it can be sourced regardless of the current working directory
# (the workflow invokes cbs.R by absolute path, not from the scripts/ directory).
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
DEFAULT_PLOIDY <- 1
DEFAULT_BINWIDTH <- 0.05
UNDO_SD_THRESHOLD <- 1.5
# segment() defaults to p.method = "hybrid", whose permutation test draws from
# the RNG, so an unseeded run does not reproduce: the same input segmented twice
# differed at a median of 16 breakpoints on H2009 hap1.
DEFAULT_SEED <- 42
# half_win auto-ploidy (see estimate_ploidy_halfwin): a segment cluster at
# [0.4,0.6]*baseline carrying > HALFWIN_CUTOFF of the genome length marks a
# whole-genome doubling; recurse up to MAX_DOUBLINGS times (ploidy 1,2,4,8).
HALFWIN_CUTOFF <- 0.09
MAX_DOUBLINGS <- 3

# Command-line argument parsing
parse_arguments <- function() {
  option_list <- list(
    make_option(c("-i", "--input"),
                type = "character",
                help = "Input copy number file from PRCGAP"),
    make_option(c("-s", "--sample"),
                type = "character",
                help = "Sample name"),
    make_option(c("-o", "--output"),
                type = "character",
                help = "Output file of circular binary segmentation result"),
    make_option(c("-p", "--ploidy"),
                type = "integer",
                default = DEFAULT_PLOIDY,
                help = "Tumor ploidy [default: %default]"),
    make_option(c("-a", "--auto-ploidy"),
                action = "store_true",
                default = FALSE,
                help = "Estimate tumor ploidy automatically from depth ratio (overrides --ploidy) [default: %default]"),
    make_option(c("-w", "--binwidth"),
                type = "double",
                default = DEFAULT_BINWIDTH,
                help = "Bin width for getting mode value of depth ratio [default: %default]"),
    make_option(c("-y", "--ploidy-out"),
                type = "character",
                default = NULL,
                help = "Write the resolved tumor ploidy (integer) to this file, for downstream use"),
    make_option(c("-S", "--seed"),
                type = "integer",
                default = DEFAULT_SEED,
                help = "RNG seed for the CBS permutation test; 0 leaves the RNG unseeded [default: %default]")
  )

  opt <- parse_args(OptionParser(option_list = option_list))

  # Validate required arguments
  if (is.null(opt$input) || is.null(opt$sample) || is.null(opt$output)) {
    stop("Error: --input, --sample, and --output are required arguments")
  }

  if (!file.exists(opt$input)) {
    stop(paste("Error: Input file does not exist:", opt$input))
  }

  return(opt)
}

# Load copy number data
load_copynumber_data <- function(input_file) {
  data <- read.table(input_file, sep = "\t", header = FALSE, comment.char = "")

  # Assign meaningful column names
  colnames(data) <- c("chromosome", "start", "end", "tumor_depth", "normal_depth")

  if (nrow(data) == 0) {
    stop("Error: Input file is empty")
  }

  return(data)
}

# Calculate normalization factors
calculate_normalization_factors <- function(data, tumor_ploidy) {
  normal_ploidy <- 1

  # Calculate total counts for high-coverage regions
  normal_count <- sum(data[data$normal_depth > MIN_COVERAGE_THRESHOLD, ]$normal_depth)
  tumor_count <- sum(data[data$normal_depth > MIN_COVERAGE_THRESHOLD, ]$tumor_depth)

  # Calculate number of rows
  normal_nrow <- nrow(data)
  tumor_nrow <- nrow(data[data$tumor_depth > MIN_TUMOR_COVERAGE, ])

  list(
    normal_count = normal_count,
    tumor_count = tumor_count,
    normal_nrow = normal_nrow,
    tumor_nrow = tumor_nrow,
    normal_ploidy = normal_ploidy,
    tumor_ploidy = tumor_ploidy
  )
}

# Calculate normalized depth ratio
calculate_depth_ratio <- function(data, norm_factors) {
  # Calculate raw depth ratio with normalization
  data$depth_ratio <- data$tumor_depth / data$normal_depth *
    norm_factors$normal_count / norm_factors$tumor_count *
    norm_factors$tumor_nrow * norm_factors$tumor_ploidy /
    (norm_factors$normal_nrow * norm_factors$normal_ploidy)

  return(data)
}

# Filter data for CBS analysis
filter_data_for_cbs <- function(data) {
  filtered_data <- data %>%
    filter(!is.na(depth_ratio) & is.finite(depth_ratio),
           normal_depth > MIN_COVERAGE_THRESHOLD)

  if (nrow(filtered_data) == 0) {
    stop("Error: No data remaining after filtering")
  }

  return(filtered_data)
}

# Prepare data for CBS analysis
prepare_cbs_data <- function(filtered_data) {
  cna_data <- data.frame(
    chromosome = gsub("chr", "", filtered_data$chromosome),
    maploc = filtered_data$end,
    log_ratio = filtered_data$depth_ratio
  )

  return(cna_data)
}

# Perform CBS segmentation
perform_cbs <- function(cna_data, sample_name, seed = DEFAULT_SEED) {
  # Create CNA object
  CNA_object <- CNA(
    genomdat = cna_data$log_ratio,
    chrom = cna_data$chromosome,
    maploc = cna_data$maploc,
    data.type = "logratio",
    sampleid = sample_name
  )

  # Smooth the data
  smoothed_CNA_object <- smooth.CNA(CNA_object)

  # Seed immediately before the only RNG consumer in the script.
  if (!is.null(seed) && !is.na(seed) && seed != 0) {
    set.seed(seed)
  }

  # Segment the data
  segment_CNA_object <- segment(
    smoothed_CNA_object,
    undo.splits = "sdundo",
    undo.SD = UNDO_SD_THRESHOLD,
    verbose = 1
  )

  return(segment_CNA_object)
}

# Write segmentation results to file
write_results <- function(segment_result, output_file) {
  write.table(
    segment_result$output,
    file = output_file,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t",
    quote = FALSE
  )

  cat("Results written to:", output_file, "\n")
}

# Main function
main <- function() {
  # Parse command-line arguments
  opt <- parse_arguments()

  cat("Loading copy number data from:", opt$input, "\n")
  data <- load_copynumber_data(opt$input)

  # Ploidy-agnostic depth ratio R0 (tumor_ploidy = 1). CBS is run on the RAW
  # ratio: ploidy is decided AFTER segmentation, from the segment levels, rather
  # than calibrating the per-window ratio beforehand.
  norm_factors <- calculate_normalization_factors(data, tumor_ploidy = 1)
  data <- calculate_depth_ratio(data, norm_factors)
  filtered_data <- filter_data_for_cbs(data)

  cna_data <- prepare_cbs_data(filtered_data)
  cat("Performing circular binary segmentation (on raw depth ratio)...\n")
  if (opt$seed != 0) {
    cat("CBS permutation-test seed:", opt$seed, "\n")
  } else {
    cat("CBS permutation-test RNG left unseeded; segments will not reproduce\n")
  }
  segment_result <- perform_cbs(cna_data, opt$sample, opt$seed)
  segs <- segment_result$output

  # Resolve ploidy and the per-copy unit mu from the CBS segment levels.
  if (isTRUE(opt[["auto-ploidy"]])) {
    est <- estimate_ploidy_halfwin(segs$seg.mean, segs$num.mark,
                                   cutoff = HALFWIN_CUTOFF, max_doublings = MAX_DOUBLINGS)
    ploidy <- est$ploidy; mu <- est$mu
    cat(sprintf("[auto-ploidy] half_win: ploidy = %d (mu = %.3f, L* = %.3f; half_win chain = %s)\n",
                ploidy, mu, est$Lstar,
                paste(sprintf("%.3f", est$half_win_chain), collapse = ";")))
  } else {
    ploidy <- opt$ploidy
    Lstar <- segment_dominant_level(segs$seg.mean, segs$num.mark)
    mu <- Lstar / ploidy
    cat(sprintf("[manual] ploidy = %d (mu = L*/ploidy = %.3f, L* = %.3f)\n", ploidy, mu, Lstar))
  }

  # Optionally write the resolved ploidy for downstream steps.
  if (!is.null(opt[["ploidy-out"]])) {
    writeLines(as.character(ploidy), opt[["ploidy-out"]])
    cat("Wrote resolved ploidy to:", opt[["ploidy-out"]], "\n")
  }

  # Calibrate segment levels to copy-number units (divide by mu). Segmenting the
  # raw ratio then scaling seg.mean is equivalent to segmenting the calibrated
  # ratio (CBS boundaries are scale-invariant), but lets mu be derived from the
  # segments themselves.
  if (is.finite(mu) && mu > 0) {
    segment_result$output$seg.mean <- segment_result$output$seg.mean / mu
  } else {
    warning("Non-finite mu; writing uncalibrated segment means.")
  }

  write_results(segment_result, opt$output)
  cat("Analysis complete!\n")
}

# Run main function
if (!interactive()) {
  main()
}
