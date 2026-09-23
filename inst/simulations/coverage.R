# Monte Carlo coverage of boot_pgt() subsampling intervals under known
# data-generating processes (envelope model).
#
# DGPs: y ~ U(1, 10); observed b = b_min(y) * exp(u) with u ~ Exp(mean
# 0.3), so the true efficiency is theta = b_min(y) / b = exp(-u) in
# (0, 1]. Two frontiers:
#   affine  b_min(y) = 2 + 0.5 y
#   curved  b_min(y) = 1 + 0.05 y^2   (convex and increasing)
# The VRS lower (y, b) envelope of the sample converges to b_min in both
# cases, so theta is the estimand of pgt(model = "envelope"). The affine
# frontier is interpolated without error by any two frontier units; the
# curved one adds the chord bias of a piecewise-linear envelope.
#
# Cells: both DGPs at L = 50 and L = 100 under boot_pgt()'s defaults
# (m = round(L^0.7), kappa = 2/3 for this two-dimensional frontier), and
# at L = 100 the subsample-size sensitivity m = round(L^0.6),
# round(L^0.8) and the rate sensitivity kappa = 1/2, 1. The L = 100
# cells of one DGP share their seeds, so they see the same samples.
#
# boot_pgt() returns NA bounds for a unit whose subsample LP is feasible
# in fewer than B/2 replicates. Under VRS a unit is feasible only in
# subsamples containing a unit with at least its output, so the units
# with the largest y can go without an interval. Each replication r
# therefore records, over its L units:
#   available   share of units with an interval (A_r)
#   cond        coverage among units with an interval (C_r)
#   delivered   share of all L units whose interval exists and covers
#               the true score (D_r)
#   avail_top   availability among the 10 per cent of units with the
#               largest y, and avail_rest among the others
#   frontier, interior   C_r restricted to theta > 0.9 and theta <= 0.9
#   width       median interval width among units with an interval
#   groups      share of the two group-mean intervals covering the
#               group's mean true score
# and each cell reports the means over R replications with Monte Carlo
# standard errors sd / sqrt(R) across replications.
#
# Run from the package root:  Rscript inst/simulations/coverage.R
# (about 20 minutes on 12 worker processes). Results:
# inst/simulations/coverage-reps.csv (one row per replication and cell)
# and inst/simulations/coverage-results.csv (cell summaries, read by the
# R Journal article's Table 3).

library(pgt)

frontiers <- list(
  affine = function(y) 2 + 0.5 * y,
  curved = function(y) 1 + 0.05 * y^2
)

run_one <- function(dgp, L, m_exp, kappa, seed) {
  set.seed(seed)
  y <- runif(L, 1, 10)
  u <- rexp(L, rate = 1 / 0.3)
  b_min <- frontiers[[dgp]](y)
  b <- b_min * exp(u)
  theta <- b_min / b
  grp <- factor(rep(c("A", "B"), length.out = L))
  tech <- pgt_tech(x = matrix(1, L, 1), y = y, b = b, group = grp)
  bt <- suppressWarnings(
    boot_pgt(tech, model = "envelope", B = 100, m = round(L^m_exp),
             kappa = kappa, seed = seed)
  )
  pd <- bt$per_dmu
  avail <- !is.na(pd$lower)
  covered <- avail & pd$lower <= theta & theta <= pd$upper
  top <- y >= stats::quantile(y, 0.9)
  gm <- bt$group_means
  theta_g <- tapply(theta, grp, mean)[gm$group]
  cond <- function(keep) {
    k <- keep & avail
    if (any(k)) mean(covered[k]) else NA_real_
  }
  data.frame(
    dgp = dgp, L = L, m = round(L^m_exp), kappa = kappa, seed = seed,
    available = mean(avail),
    cond = cond(rep(TRUE, L)),
    delivered = mean(covered),
    avail_top = mean(avail[top]),
    avail_rest = mean(avail[!top]),
    frontier = cond(theta > 0.9),
    interior = cond(theta <= 0.9),
    width = stats::median((pd$upper - pd$lower)[avail]),
    groups = mean(gm$lower <= theta_g & theta_g <= gm$upper),
    stringsAsFactors = FALSE
  )
}

R_REPS <- 200
cells <- rbind(
  expand.grid(dgp = names(frontiers), L = c(50, 100), m_exp = 0.7,
              kappa = 2 / 3, stringsAsFactors = FALSE),
  expand.grid(dgp = names(frontiers), L = 100, m_exp = c(0.6, 0.8),
              kappa = 2 / 3, stringsAsFactors = FALSE),
  expand.grid(dgp = names(frontiers), L = 100, m_exp = 0.7,
              kappa = c(1 / 2, 1), stringsAsFactors = FALSE)
)
# the affine seeds are those of the study shipped with pgt 0.7.0
seed_of <- function(dgp, L, r) {
  (dgp == "curved") * 1e6 + 1000 * L + r
}
jobs <- do.call(c, lapply(seq_len(nrow(cells)), function(k) {
  lapply(seq_len(R_REPS), function(r) {
    c(as.list(cells[k, ]), seed = seed_of(cells$dgp[k], cells$L[k], r))
  })
}))

cl <- parallel::makeCluster(max(1L, parallel::detectCores() - 2L))
invisible(parallel::clusterEvalQ(cl, library(pgt)))
parallel::clusterExport(cl, c("frontiers", "run_one"))
reps <- do.call(rbind, parallel::parLapplyLB(cl, jobs, function(j) {
  run_one(j$dgp, j$L, j$m_exp, j$kappa, j$seed)
}))
parallel::stopCluster(cl)

se <- function(z) stats::sd(z, na.rm = TRUE) / sqrt(sum(!is.na(z)))
res <- do.call(rbind, lapply(
  split(reps, list(reps$dgp, reps$L, reps$m, reps$kappa), drop = TRUE),
  function(d) data.frame(
    dgp = d$dgp[1], L = d$L[1], m = d$m[1], kappa = d$kappa[1],
    reps = nrow(d),
    available = mean(d$available), available_se = se(d$available),
    avail_top = mean(d$avail_top), avail_rest = mean(d$avail_rest),
    coverage_cond = mean(d$cond), coverage_cond_se = se(d$cond),
    coverage_delivered = mean(d$delivered),
    coverage_delivered_se = se(d$delivered),
    coverage_frontier = mean(d$frontier, na.rm = TRUE),
    coverage_interior = mean(d$interior, na.rm = TRUE),
    mean_median_width = mean(d$width),
    coverage_group_means = mean(d$groups),
    stringsAsFactors = FALSE
  )))
res <- res[order(res$dgp, res$L, res$m, res$kappa), ]
rownames(res) <- NULL

print(res, row.names = FALSE, digits = 3)
# written into the package sources when run from the source root, and
# next to a temporary copy otherwise (a sourced installed copy has no
# inst/ tree to write into)
out_dir <- file.path("inst", "simulations")
if (!dir.exists(out_dir)) out_dir <- tempdir()
write.csv(reps, file.path(out_dir, "coverage-reps.csv"), row.names = FALSE)
write.csv(res, file.path(out_dir, "coverage-results.csv"), row.names = FALSE)
message("coverage results written to ", out_dir)
