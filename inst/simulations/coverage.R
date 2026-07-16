# Monte Carlo coverage of boot_pgt() subsampling intervals under a
# known data-generating process (envelope model).
#
# DGP: true frontier b_min(y) = 2 + 0.5 y (convex, attainable);
# y ~ U(1, 10); observed b = b_min(y) * exp(u) with u ~ Exp(mean 0.3),
# so the true efficiency is theta = b_min(y) / b = exp(-u) in (0, 1].
# The convex lower (y, b) envelope of the sample converges to the true
# frontier, so theta is the estimand of pgt(model = "envelope").
#
# For each sample size L, R Monte Carlo replicates draw a sample, fit
# boot_pgt(model = "envelope", B = 100), and record whether each DMU's
# nominal 95 percent interval covers its true theta, separately for
# near-frontier units (theta > 0.9) and interior units, plus the
# group-mean intervals against the group's mean true theta.
#
# Run from the package root:  Rscript inst/simulations/coverage.R
# Results: inst/simulations/coverage-results.csv (shipped with the
# package; the productivity vignette reports this table).

library(pgt)

run_one <- function(L, seed) {
  set.seed(seed)
  y <- runif(L, 1, 10)
  u <- rexp(L, rate = 1 / 0.3)
  b_min <- 2 + 0.5 * y
  b <- b_min * exp(u)
  theta <- b_min / b
  grp <- factor(rep(c("A", "B"), length.out = L))
  tech <- pgt_tech(x = matrix(1, L, 1), y = y, b = b, group = grp)
  bt <- suppressWarnings(
    boot_pgt(tech, model = "envelope", B = 100, seed = seed)
  )
  pd <- bt$per_dmu
  covered <- pd$lower <= theta & theta <= pd$upper
  gm <- bt$group_means
  theta_g <- tapply(theta, grp, mean)[gm$group]
  data.frame(
    L = L,
    coverage_all = mean(covered, na.rm = TRUE),
    coverage_frontier = mean(covered[theta > 0.9], na.rm = TRUE),
    coverage_interior = mean(covered[theta <= 0.9], na.rm = TRUE),
    width_median = stats::median(pd$upper - pd$lower, na.rm = TRUE),
    coverage_groups = mean(gm$lower <= theta_g & theta_g <= gm$upper,
                           na.rm = TRUE)
  )
}

R_REPS <- 200
grid <- c(50, 100)
res <- do.call(rbind, lapply(grid, function(L) {
  reps <- do.call(rbind, lapply(seq_len(R_REPS), function(r) {
    run_one(L, seed = 1000 * L + r)
  }))
  data.frame(
    L = L,
    reps = R_REPS,
    coverage_all = mean(reps$coverage_all, na.rm = TRUE),
    coverage_frontier = mean(reps$coverage_frontier, na.rm = TRUE),
    coverage_interior = mean(reps$coverage_interior, na.rm = TRUE),
    median_width = mean(reps$width_median, na.rm = TRUE),
    coverage_group_means = mean(reps$coverage_groups, na.rm = TRUE)
  )
}))

print(res, row.names = FALSE)
write.csv(res, file.path("inst", "simulations", "coverage-results.csv"),
          row.names = FALSE)
