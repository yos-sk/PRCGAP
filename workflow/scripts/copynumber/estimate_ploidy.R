# !/usr/bin/env Rscript

# Automatic ploidy estimation (CN=1 vs CN=2 baseline) from per-haplotype

detect_peaks <- function(depth_ratio, min_dr = 0.30, max_dr = 5,
                         bw_factor = 1.0,
                         min_prominence_ratio = 0.20,
                         n_peaks_max = 8) {
  # Keep only depth ratios in a plausible range; drop near-zero and outlier bins.
  positive <- depth_ratio[depth_ratio > min_dr & depth_ratio < max_dr]
  # Need enough bins for a stable density estimate; otherwise return no peaks.
  if (length(positive) < 100) return(data.frame())
  bw <- bw.SJ(positive) * bw_factor # Sheather & Jones method for automatic bandwidth selection
  # Kernel density estimate of the depth-ratio distribution on a fixed grid.
  dens <- density(positive, bw = bw, n = 2048, from = min_dr, to = max_dr)
  # Local maxima of the KDE: a point is a peak where the sign of the slope
  # (diff(dens$y)) flips from + to - (diff(sign(...)) == -2). The FALSE pads keep
  # the logical vector aligned with dens$x and exclude the two endpoints.
  is_max <- c(FALSE, diff(sign(diff(dens$y))) == -2, FALSE)
  cand <- data.frame(x = dens$x[is_max], y = dens$y[is_max])
  # Prominence filter: keep peaks at least min_prominence_ratio of the tallest
  # peak's height, dropping shallow noise bumps.
  cand <- cand[cand$y >= max(dens$y) * min_prominence_ratio, ]
  # Keep the n_peaks_max tallest peaks, then return them sorted by position x.
  cand <- cand[order(cand$y, decreasing = TRUE), ]
  cand <- head(cand, n_peaks_max)
  cand[order(cand$x), ]
}

# ---------------------------------------------------------------------------
# CBS-segment ploidy (half_win): the ploidy / per-copy unit are decided from the
# (length-weighted, denoised) CBS segment levels rather than from the raw
# per-window depth-ratio distribution. (detect_peaks() above is retained because
# plot_copy_number.R still uses it for its own peak-comb calibration.)
# ---------------------------------------------------------------------------

# Length-weighted dominant segment level L*, with a robust fallback when the
# kernel density cannot be estimated (too few / too sparse segments).
segment_dominant_level <- function(seg_mean, num_mark) {
  keep <- is.finite(seg_mean) & seg_mean > 0.05
  lv <- seg_mean[keep]; w <- num_mark[keep]
  if (length(lv) == 0) return(NA_real_)
  wf <- w / sum(w)
  tryCatch({
    d <- density(lv, weights = wf, bw = "SJ", n = 2048, from = 0.1, to = max(lv) * 1.05)
    d$x[which.max(d$y)]
  }, error = function(e) lv[which.max(w)])
}

# Recursive half_win ploidy. While a distinct segment cluster sits at HALF the
# current baseline (the whole-genome-doubling signature), halve the baseline and
# double k. Yields ploidy in {1, 2, 4, 8, ...} and mu = L*/ploidy (depth ratio
# per single copy). Diploid scattered deletions do NOT concentrate at half, so
# `half_win` separates diploid from doubled with a wide margin (cutoff ~0.09).
estimate_ploidy_halfwin <- function(seg_mean, num_mark, cutoff = 0.09, max_doublings = 3) {
  keep <- is.finite(seg_mean) & seg_mean > 0.05
  lv <- seg_mean[keep]; w <- num_mark[keep]
  if (length(lv) == 0) {
    return(list(ploidy = 1L, mu = NA_real_, Lstar = NA_real_, half_win_chain = numeric(0)))
  }
  wf <- w / sum(w)
  L <- segment_dominant_level(seg_mean, num_mark)
  k <- 1L; B <- L; chain <- numeric(0)
  for (i in seq_len(max_doublings)) {
    hw <- sum(wf[lv > 0.40 * B & lv < 0.60 * B])   # weight at [0.4,0.6]*baseline
    chain <- c(chain, hw)
    if (hw > cutoff) { B <- B / 2; k <- k * 2L } else break
  }
  list(ploidy = k, mu = L / k, Lstar = L, half_win_chain = chain)
}
