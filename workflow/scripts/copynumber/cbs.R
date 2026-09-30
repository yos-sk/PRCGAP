#!/usr/bin/env Rscript
# Circular Binary Segmentation (CBS) analysis for copy number variation detection
# This script performs CBS on copy number data from PRCGAP output
#
# Two passes over the ploidy-agnostic depth ratio R0 (tumor_ploidy = 1):
#   Pass 1  CBS on the raw R0 with plain sdundo, used only to resolve ploidy from
#           the segment levels (estimate_ploidy_halfwin(); L* = dominant level).
#   mu      depth ratio per copy = L_bin / ploidy, where L_bin is the tallest KDE
#           peak of the per-window R0 (the calibration plot_copy_number.R uses
#           for the points); L* / ploidy only when no peak is found.
#   Pass 2  CBS on R0 / mu, i.e. in copy-number units, so CBS_INPUT_CAP and
#           MIN_LONG_STEP apply directly: winsorized input, length-aware undo,
#           then segment means recomputed from the uncapped values.

library(DNAcopy)
library(tidyverse)
library(optparse)

# CBS-segment ploidy estimator (estimate_ploidy_halfwin(), segment_dominant_level(),
# detect_peaks()). estimate_ploidy.R lives alongside this script; resolve this
# script's own directory so it can be sourced regardless of the current working
# directory (the workflow invokes cbs.R by absolute path, not from the scripts/
# directory).
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
# the RNG, so an unseeded run does not reproduce -- the same input segmented
# twice returns different breakpoints.
DEFAULT_SEED <- 42
# half_win auto-ploidy (see estimate_ploidy_halfwin): a segment cluster at
# [0.4,0.6]*baseline carrying > HALFWIN_CUTOFF of the genome length marks a
# whole-genome doubling; recurse up to MAX_DOUBLINGS times (ploidy 1,2,4,8).
HALFWIN_CUTOFF <- 0.09
MAX_DOUBLINGS <- 3

# Length-aware undo (see length_aware_undo): a changepoint between two segments
# of >= MIN_LONG_MARKS marks is merged only if its step < min(UNDO_SD_THRESHOLD*SD,
# max(UNDO_SD_LONG*SD, MIN_LONG_STEP)); all other changepoints use plain sdundo
# (step < UNDO_SD_THRESHOLD*SD). MIN_LONG_STEP is in copy-number units; pass 2
# runs on the mu-calibrated ratio, so it is used directly, with no conversion.
UNDO_SD_LONG <- 1.0
MIN_LONG_STEP <- 0.3
MIN_LONG_MARKS <- 100

# Winsorization cap (copy-number units) for the CBS input only: extreme focal
# amplifications (HCC1954 ERBB2, ~250 copies) inflate the within-chromosome
# variance and collapse the segmentation. Segment means are restored from
# uncapped values afterwards. Applied directly to the pass-2 mu-calibrated
# ratio (mu is only known after the first, unwinsorized, pass).
CBS_INPUT_CAP <- 20

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
                help = "RNG seed for the CBS permutation test; 0 leaves the RNG unseeded [default: %default]"),
    make_option(c("--undo-sd"),
                type = "double",
                default = UNDO_SD_THRESHOLD,
                help = "SD threshold (in trimmed SDs) for the short-pair undo rule [default: %default]"),
    make_option(c("--undo-sd-long"),
                type = "double",
                default = UNDO_SD_LONG,
                help = "SD threshold (in trimmed SDs) for the long-pair undo rule [default: %default]"),
    make_option(c("--min-long-marks"),
                type = "double",
                default = MIN_LONG_MARKS,
                help = "Minimum number of marks for a segment to count as 'long' (use Inf to disable the long-pair rule) [default: %default]"),
    make_option(c("--min-long-step"),
                type = "double",
                default = MIN_LONG_STEP,
                help = "Floor, in copy-number units, on the long-pair undo threshold [default: %default]")
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

# Build the smoothed CNA object DNAcopy works on. Same bins throughout; only the
# values differ between the (possibly winsorized) changepoint-detection pass and
# the uncapped pass used to recompute segment means.
smooth_ratios <- function(cna_data, values, sample_name) {
  CNA_object <- CNA(
    genomdat = values,
    chrom = cna_data$chromosome,
    maploc = cna_data$maploc,
    data.type = "logratio",
    sampleid = sample_name
  )

  smooth.CNA(CNA_object)
}

# Find changepoints with no undo step (undo.splits = "none"); the caller applies
# length_aware_undo() afterwards (min_long = Inf reproduces plain sdundo). Returns
# the segmentation plus the smoothed CNA object, which length_aware_undo() needs
# to recompute segment medians. The seed is set before smoothing; smooth.CNA
# draws no RNG, so only segment()'s permutation test depends on it.
find_changepoints <- function(cna_data, values, sample_name, seed) {
  if (!is.null(seed) && !is.na(seed) && seed != 0) {
    set.seed(seed)
  }

  smoothed <- smooth_ratios(cna_data, values, sample_name)

  segment_result <- segment(
    smoothed,
    undo.splits = "none",
    verbose = 1
  )

  list(segment_result = segment_result, smoothed = smoothed)
}

# Length-aware analogue of DNAcopy's changepoints.sdundo (same loop; min_long = Inf
# reproduces plain segment(undo.splits = "sdundo")). Merges the weakest adjacent
# pair below thr_long (both flanking segments >= min_long marks) or thr_short
# otherwise. undo_sd and undo_sd_long are multiples of the trimmed SD;
# min_long_step is in the units of `smoothed` (copy-number units in pass 2). When
# `smoothed` is in raw depth-ratio units, pass `mu` (raw units per copy) to also
# log the thresholds in copy-number units.
length_aware_undo <- function(segment_result, smoothed, undo_sd,
                               undo_sd_long, min_long, min_long_step,
                               mu = NA_real_) {
  all_values <- smoothed[[3]]
  trimmed_sd <- sqrt(DNAcopy:::trimmed.variance(all_values[!is.na(all_values)], 0.025))
  thr_short <- undo_sd * trimmed_sd
  thr_long <- min(thr_short, max(undo_sd_long * trimmed_sd, min_long_step))

  out <- segment_result$output
  chrom_order <- unique(as.character(out$chrom))

  smoothed_chrom <- as.character(smoothed$chrom)
  smoothed_maploc <- smoothed$maploc
  smoothed_value <- smoothed[[3]]

  new_rows <- vector("list", length(chrom_order))
  n_cp_before <- 0L
  n_cp_after <- 0L
  n_removed_short <- 0L
  n_removed_long <- 0L

  for (ci in seq_along(chrom_order)) {
    ch <- chrom_order[ci]
    idx <- which(as.character(out$chrom) == ch)
    lseg <- out$num.mark[idx]
    sample_id <- out$ID[idx][1]

    chrom_idx <- which(smoothed_chrom == ch)
    ord <- chrom_idx[order(smoothed_maploc[chrom_idx])]
    genomdat <- smoothed_value[ord]
    maploc_sorted <- smoothed_maploc[ord]

    n_cp_before <- n_cp_before + (length(lseg) - 1L)

    cpt.loc <- cumsum(lseg)

    repeat {
      k <- length(cpt.loc)
      if (k <= 1) break

      segments0 <- cbind(c(1, 1 + cpt.loc[-k]), cpt.loc)
      segmed <- apply(segments0, 1, function(i, x) median(x[i[1]:i[2]]), genomdat)
      adsegmed <- abs(diff(segmed))

      seg_lengths <- diff(c(0, cpt.loc))
      n_left <- seg_lengths[-length(seg_lengths)]
      n_right <- seg_lengths[-1]

      is_long_pair <- pmin(n_left, n_right) >= min_long
      removable <- ifelse(is_long_pair, adsegmed < thr_long, adsegmed < thr_short)

      if (!any(removable)) break

      cand <- which(removable)
      best <- cand[which.min(adsegmed[cand])]

      if (is_long_pair[best]) {
        n_removed_long <- n_removed_long + 1L
      } else {
        n_removed_short <- n_removed_short + 1L
      }

      cpt.loc <- cpt.loc[-best]
    }

    n_cp_after <- n_cp_after + (length(cpt.loc) - 1L)

    seg_lengths_final <- diff(c(0, cpt.loc))
    starts <- c(1, 1 + cpt.loc[-length(cpt.loc)])
    ends <- cpt.loc

    new_rows[[ci]] <- data.frame(
      ID = sample_id,
      chrom = ch,
      loc.start = maploc_sorted[starts],
      loc.end = maploc_sorted[ends],
      num.mark = seg_lengths_final,
      seg.mean = vapply(seq_along(starts), function(s) {
        round(mean(genomdat[starts[s]:ends[s]]), 4)
      }, numeric(1)),
      stringsAsFactors = FALSE
    )
  }

  new_out <- do.call(rbind, new_rows)
  rownames(new_out) <- NULL

  cat(sprintf(
    paste0("[length-aware-undo] trimmed SD = %.6f; short-pair threshold = %.6f; ",
           "long-pair threshold = %.6f; changepoints before = %d, after = %d; ",
           "removed by short-pair rule = %d, removed by long-pair rule = %d\n"),
    trimmed_sd, thr_short, thr_long, n_cp_before, n_cp_after, n_removed_short, n_removed_long))

  if (is.finite(mu) && mu > 0) {
    cat(sprintf(
      "[length-aware-undo] (copy units) short-pair threshold = %.6f copies; long-pair threshold = %.6f copies\n",
      thr_short / mu, thr_long / mu))
  }

  segment_result$output <- new_out
  segment_result
}

# Recompute every segment mean from the uncapped smoothed values (DNAcopy reports
# the mean of the possibly-capped input) so amplicon segments show their true
# magnitude. Always recomputes: length_aware_undo() may leave a placeholder mean.
recompute_segment_means <- function(segment_result, cna_data, sample_name) {
  smoothed <- smooth_ratios(cna_data, cna_data$log_ratio, sample_name)
  values <- smoothed[[3]]

  out <- segment_result$output

  out$seg.mean <- vapply(seq_len(nrow(out)), function(i) {
    in_segment <- as.character(smoothed$chrom) == as.character(out$chrom[i]) &
      smoothed$maploc >= out$loc.start[i] &
      smoothed$maploc <= out$loc.end[i]
    round(mean(values[in_segment]), 4)
  }, numeric(1))

  segment_result$output <- out
  return(segment_result)
}

# Tallest peak (by KDE height) of the raw depth-ratio distribution, via the
# detect_peaks() KDE from estimate_ploidy.R. This gives an independent,
# per-window estimate of the dominant copy-number state's raw depth ratio
# (L_bin), used to calibrate mu alongside the CBS-segment-level L* from
# estimate_ploidy_halfwin()/segment_dominant_level().
bin_dominant_level <- function(ratio) {
  vals <- ratio[is.finite(ratio) & ratio > 0.05]
  peaks <- detect_peaks(vals)
  if (nrow(peaks) == 0) return(NA_real_)
  peaks$x[which.max(peaks$y)]
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

  # Ploidy-agnostic depth ratio R0 (tumor_ploidy = 1). Pass 1 segments the raw
  # ratio to decide ploidy from the segment levels (half_win); pass 2 re-segments
  # the mu-calibrated ratio in copy-number units.
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

  # --- Pass 1: unwinsorized segmentation, plain sdundo, to resolve ploidy/mu ---
  cp0 <- find_changepoints(cna_data, cna_data$log_ratio, opt$sample, opt$seed)
  segment_result1 <- length_aware_undo(
    cp0$segment_result, cp0$smoothed,
    undo_sd = opt[["undo-sd"]], undo_sd_long = opt[["undo-sd"]],
    min_long = Inf, min_long_step = 0)
  segs1 <- segment_result1$output

  # Resolve ploidy from the pass-1 segment levels.
  if (isTRUE(opt[["auto-ploidy"]])) {
    est <- estimate_ploidy_halfwin(segs1$seg.mean, segs1$num.mark,
                                   cutoff = HALFWIN_CUTOFF, max_doublings = MAX_DOUBLINGS)
    ploidy <- est$ploidy; mu <- est$mu; Lstar <- est$Lstar
    cat(sprintf("[auto-ploidy] half_win: ploidy = %d (mu = %.3f, L* = %.3f; half_win chain = %s)\n",
                ploidy, mu, est$Lstar,
                paste(sprintf("%.3f", est$half_win_chain), collapse = ";")))
  } else {
    ploidy <- opt$ploidy
    Lstar <- segment_dominant_level(segs1$seg.mean, segs1$num.mark)
    mu <- Lstar / ploidy
    cat(sprintf("[manual] ploidy = %d (mu = L*/ploidy = %.3f, L* = %.3f)\n", ploidy, mu, Lstar))
  }

  # Refine mu using the tallest peak of the raw per-window depth-ratio KDE
  # (bin_dominant_level), falling back to the CBS-segment-level L*/ploidy when
  # no peak is found.
  L_bin <- bin_dominant_level(cna_data$log_ratio)
  if (is.finite(L_bin)) {
    mu <- L_bin / ploidy
    cat(sprintf("mu = L_bin / ploidy = %.4f (L_bin = %.4f, half_win L* = %.4f)\n",
                mu, L_bin, Lstar))
  } else {
    mu <- Lstar / ploidy
    cat(sprintf("L_bin not available (no peak detected); falling back to mu = L*/ploidy = %.4f\n", mu))
  }

  # Optionally write the resolved ploidy for downstream steps.
  if (!is.null(opt[["ploidy-out"]])) {
    writeLines(as.character(ploidy), opt[["ploidy-out"]])
    cat("Wrote resolved ploidy to:", opt[["ploidy-out"]], "\n")
  }

  # --- Pass 2: re-segment the mu-calibrated ratio in copy-number units (cap,
  # length-aware undo and mean recomputation all in copy units). segment() is
  # not numerically scale-invariant, so calibrating first is required for
  # reproducible copy-unit thresholds. ---
  cal <- cna_data
  cal$log_ratio <- cna_data$log_ratio / mu

  n_capped <- sum(cal$log_ratio > CBS_INPUT_CAP)
  if (n_capped > 0) {
    cat("Winsorizing", n_capped, "bin(s) above", CBS_INPUT_CAP,
        "for changepoint detection (max =",
        round(max(cal$log_ratio), 1), ")\n")
  }

  cp <- find_changepoints(cal, pmin(cal$log_ratio, CBS_INPUT_CAP), opt$sample, opt$seed)

  segment_result <- length_aware_undo(
    cp$segment_result, cp$smoothed,
    undo_sd = opt[["undo-sd"]], undo_sd_long = opt[["undo-sd-long"]],
    min_long = opt[["min-long-marks"]], min_long_step = opt[["min-long-step"]])
  segment_result <- recompute_segment_means(segment_result, cal, opt$sample)

  write_results(segment_result, opt$output)
  cat("Analysis complete!\n")
}

# Run main function
if (!interactive()) {
  main()
}
