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

peak_comb_fit <- function(peaks, k_tallest = 1L, max_cn = 6, n_grid = 400) {
  if (nrow(peaks) == 0) return(list(mu = NA, max_rel_residual = NA))
  tallest <- peaks[which.max(peaks$y), ]
  mu0 <- tallest$x / k_tallest
  # mu = depth per single copy. Grid search over candidate mu values within
  # +/-15% of the initial guess mu0, spaced evenly on a log scale (multiplicative
  # steps are natural for a depth-per-copy quantity).
  grid <- exp(seq(log(mu0 * 0.85), log(mu0 * 1.15), length.out = n_grid))
  # vapply evaluates the function once per candidate mu in grid, so scores is a
  # vector aligned 1:1 with grid (scores[i] is the score for grid[i]).
  # For each mu: assign each peak its integer copy number k = round(x / mu),
  # clipped to [1, max_cn], then score = weighted sum of squared residuals
  # between observed peak position x and its theoretical position k * mu.
  # Smaller score = peaks sit more cleanly on integer multiples of mu.
  scores <- vapply(grid, function(mu) {
    k <- pmin(pmax(round(peaks$x / mu), 1), max_cn)
    sum(peaks$y * (peaks$x - k * mu)^2)
  }, numeric(1))
  # Pick the mu with the smallest residual score.
  best_mu <- grid[which.min(scores)]
  # Finalize integer copy numbers under best_mu and report the worst relative
  # residual (normalized by mu) as a goodness-of-fit / confidence indicator.
  k_assign <- pmin(pmax(round(peaks$x / best_mu), 1), max_cn)
  rel_residual <- (peaks$x - k_assign * best_mu) / best_mu
  list(mu = best_mu, max_rel_residual = max(abs(rel_residual)))
}

# Coverage score. For each peak find its integer tooth k. If within
# dist_threshold*mu, it rewards tooth k (best peak per tooth). Otherwise it is
# a rogue. A sub-baseline rogue (peak below the first tooth, x < mu) is the
# LOH / single-copy signature => the dominant peak is CN>=2, so it is weighted
# by sub_penalty (> rogue_penalty) to penalize the CN=1 prior. Empty teeth in
# the active range are penalized too.
# Tuning parameters (penalty ordering: sub_penalty > rogue_penalty > empty_penalty):
#   dist_threshold : max relative distance |x - k*mu| / mu for a peak to count as
#                    a hit on tooth k. Reward = y * (1 - rel_dist/dist_threshold),
#                    so max on the tooth and ~0 at the threshold. Beyond it the
#                    peak is a rogue.
#   empty_penalty  : penalty per empty tooth within the used range k_min:k_max
#                    (discourages skipped teeth, e.g. CN=3 present but CN=2 missing).
#   rogue_penalty  : weight for a normal rogue peak (no tooth hit) sitting at or
#                    above the baseline (x >= mu).
#   sub_penalty    : weight for a sub-baseline rogue (x < mu), the LOH/single-copy
#                    signature. Heavier than rogue_penalty to penalize the CN=1 prior.
comb_coverage_score <- function(peaks, mu,
                                dist_threshold = 0.25,
                                empty_penalty  = 0.30,
                                rogue_penalty  = 1.0,
                                sub_penalty    = 2.5,
                                max_cn         = 8) {
  if (is.na(mu) || nrow(peaks) == 0) return(NA_real_)
  k_int <- pmin(pmax(round(peaks$x / mu), 1), max_cn)
  rel_dist <- abs(peaks$x - k_int * mu) / mu
  reward_per_tooth <- numeric(max_cn)
  rogue <- 0
  for (i in seq_along(k_int)) {
    if (rel_dist[i] < dist_threshold) {
      r <- peaks$y[i] * (1 - rel_dist[i] / dist_threshold)
      reward_per_tooth[k_int[i]] <- max(reward_per_tooth[k_int[i]], r)
    } else {
      w <- if (peaks$x[i] < mu) sub_penalty else rogue_penalty
      rogue <- rogue + w * peaks$y[i]
    }
  }
  k_min <- min(k_int); k_max <- max(k_int)
  empty <- sum(reward_per_tooth[k_min:k_max] == 0)
  sum(reward_per_tooth) - empty_penalty * empty - rogue
}

# Main entry point. Returns a list with the integer ploidy, the baseline mu,
# a confidence label, and the underlying scores/peaks for logging.
estimate_ploidy <- function(depth_ratio, score_margin = 0.05, sub_penalty = 2.5) {
  peaks <- detect_peaks(depth_ratio)
  if (nrow(peaks) == 0) {
    return(list(ploidy = 1L, mu = NA_real_, confidence = "low (no peaks)",
                cn1_score = NA_real_, cn2_score = NA_real_, peaks = peaks))
  }
  f1 <- peak_comb_fit(peaks, k_tallest = 1)
  f2 <- peak_comb_fit(peaks, k_tallest = 2)
  s1 <- comb_coverage_score(peaks, f1$mu, sub_penalty = sub_penalty)
  s2 <- comb_coverage_score(peaks, f2$mu, sub_penalty = sub_penalty)

  ploidy <- if (!is.na(s1) && !is.na(s2) && s2 > s1 + score_margin) 2L else 1L
  mu     <- if (ploidy == 2L) f2$mu else f1$mu
  diff   <- if (is.na(s1) || is.na(s2)) NA_real_ else abs(s2 - s1)

  confidence <- if (nrow(peaks) <= 1) "high (single peak -> CN=1)"
                else if (is.na(diff))         "low (manual review)"
                else if (diff > 0.5)          "high"
                else if (diff > score_margin) "medium"
                else                          "low (manual review)"

  list(ploidy = ploidy, mu = mu, confidence = confidence,
       cn1_score = s1, cn2_score = s2, peaks = peaks)
}
