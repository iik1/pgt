#' Subsampling inference for pollution-generating efficiency scores
#'
#' Constructs subsampling confidence intervals for the environmental
#' efficiency scores of a pgt model and for group means, by re-solving
#' each DMU against random subsamples of the reference set. Subsampling
#' (m out of n, without replacement) is the resampling scheme with the
#' firmest footing for the boundary-estimator, slow-convergence setting
#' of nonparametric frontiers (Kneip, Simar and Wilson 2008; Simar and
#' Wilson 2011).
#'
#' The intervals are heuristic. Bootstrap and subsampling theory for
#' data envelopment analysis has been developed for the radial and
#' directional-distance technologies, not for the materials-balance-
#' constrained programs implemented here, so the coverage of these
#' intervals is not backed by a convergence-rate result. Use the
#' subsample-size sensitivity of [boot_pgt_sensitivity()] to check that
#' the intervals are stable in \code{m} before relying on them.
#'
#' @param tech A [pgt_tech()] object.
#' @param model One of the environmental-efficiency models
#'   (\code{"wgd"}, \code{"envelope"}, \code{"byprod"}, \code{"mb_cost"},
#'   \code{"wd"}); the directional model \code{"fdmo"} is not supported.
#' @param B Number of subsampling replicates.
#' @param m Subsample size. Defaults to \code{round(L^0.7)}, a
#'   conventional choice between \eqn{\sqrt{L}} and \eqn{L}.
#' @param level Confidence level for the percentile intervals.
#' @param returns,peers,pollutant Passed to the underlying model, as in
#'   [pgt()].
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An object of class \code{"pgt_boot"}: a list with
#'   \code{per_dmu} (a data frame of point estimate, lower and upper
#'   bound, and bootstrap standard error per DMU), \code{group_means}
#'   (the same for group means, when a group is present), and the call
#'   settings.
#'
#' @references
#' Kneip, A., Simar, L., & Wilson, P. W. (2008). Asymptotics and
#' consistent bootstraps for DEA estimators in nonparametric frontier
#' models. \emph{Econometric Theory}, 24(6), 1663--1697.
#'
#' @seealso [pgt()], [boot_pgt_sensitivity()]
#' @examples
#' data(steeldemo)
#' tech <- pgt_tech(
#'   x = steeldemo[, c("coal_coke", "other_fuel", "raw_material", "flux")],
#'   y = steeldemo$production, b = steeldemo$emissions, v = 0.01467,
#'   group = steeldemo$route, id = steeldemo$plant
#' )
#' \donttest{
#' bt <- boot_pgt(tech, model = "wgd", B = 100, seed = 1)
#' head(bt$per_dmu)
#' }
#' @export
boot_pgt <- function(tech, model = c("wgd", "envelope", "byprod",
                                     "mb_cost", "wd"),
                     B = 200, m = NULL, level = 0.95,
                     returns = c("vrs", "crs"), peers = c("all", "group"),
                     pollutant = 1L, seed = NULL) {
  stopifnot(inherits(tech, "pgt_tech"))
  model <- match.arg(model)
  returns <- match.arg(returns)
  peers <- match.arg(peers)
  vrs <- returns == "vrs"
  p <- .pollutant_index(tech, pollutant)
  L <- tech$L
  if (is.null(m)) m <- max(2L, round(L^0.7))
  if (m < 2L || m > L) {
    stop("'m' must lie in 2..L.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  point <- .model_efficiency(tech, model, vrs, p,
                             .peer_sets(tech, peers))
  boot <- matrix(NA_real_, L, B)
  grp <- tech$group
  do_group <- peers == "group" && !is.null(grp)
  for (rrep in seq_len(B)) {
    if (do_group) {
      # subsample within each group so every DMU keeps in-group peers,
      # allocating the overall size m across groups in proportion to
      # group size so that m still drives the resampling
      sub <- unlist(lapply(split(seq_len(L), grp), function(ix) {
        mm <- max(2L, min(length(ix), round(length(ix) * m / L)))
        if (length(ix) <= mm) ix else sample(ix, mm)
      }))
    } else {
      sub <- sample(seq_len(L), m)
    }
    peer_sets <- lapply(seq_len(L), function(i) {
      base <- if (do_group) sub[grp[sub] == grp[i]] else sub
      unique(c(base, i))
    })
    boot[, rrep] <- .model_efficiency(tech, model, vrs, p, peer_sets)
  }

  a <- (1 - level) / 2
  per_dmu <- data.frame(
    id = tech$id,
    estimate = point,
    lower = apply(boot, 1, stats::quantile, probs = a, na.rm = TRUE),
    upper = apply(boot, 1, stats::quantile, probs = 1 - a, na.rm = TRUE),
    se = apply(boot, 1, stats::sd, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  if (!is.null(grp)) per_dmu <- cbind(per_dmu[1], group = grp, per_dmu[-1])
  rownames(per_dmu) <- NULL

  group_means <- NULL
  if (!is.null(grp)) {
    glev <- levels(grp)
    gm_boot <- vapply(glev, function(g) {
      apply(boot[grp == g, , drop = FALSE], 2, mean, na.rm = TRUE)
    }, numeric(B))
    group_means <- data.frame(
      group = glev,
      estimate = vapply(glev, function(g) mean(point[grp == g], na.rm = TRUE),
                        numeric(1)),
      lower = apply(gm_boot, 2, stats::quantile, probs = a, na.rm = TRUE),
      upper = apply(gm_boot, 2, stats::quantile, probs = 1 - a, na.rm = TRUE),
      se = apply(gm_boot, 2, stats::sd, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    rownames(group_means) <- NULL
  }

  structure(
    list(per_dmu = per_dmu, group_means = group_means, model = model,
         B = B, m = m, level = level, returns = returns, peers = peers),
    class = "pgt_boot"
  )
}

# Environmental efficiency of every DMU under `model`, evaluated against
# the supplied per-DMU peer sets.
.model_efficiency <- function(tech, model, vrs, p, peer_sets) {
  b_i <- tech$b[, p]
  vapply(seq_len(tech$L), function(i) {
    sol <- .lp_solve_one(model, i, tech, peer_sets[[i]], vrs, p = p)
    if (!is.null(sol$status) && sol$status != 0) return(NA_real_)
    switch(model,
      byprod = sol$emission_eff,
      mb_cost = sol$mbe,
      sol$b_star / b_i[i])
  }, numeric(1))
}

#' Subsample-size sensitivity of subsampling intervals
#'
#' Re-runs [boot_pgt()] over a grid of subsample sizes \code{m} and
#' reports how the median interval width responds. Stable widths across
#' the grid support the choice of \code{m}; a strong trend warns that the
#' intervals are driven by the resampling size rather than the data.
#'
#' @param tech A [pgt_tech()] object.
#' @param model,B,level,returns,peers,pollutant,seed As in [boot_pgt()].
#' @param m_grid Integer vector of subsample sizes. Defaults to five
#'   sizes spanning \eqn{L^{0.5}} to \eqn{L^{0.9}}.
#'
#' @return A data frame with one row per \code{m}: the size, the median
#'   per-DMU interval width, and the mean bootstrap standard error.
#'
#' @seealso [boot_pgt()]
#' @examples
#' data(steeldemo)
#' tech <- pgt_tech(
#'   x = steeldemo[, c("coal_coke", "other_fuel", "raw_material", "flux")],
#'   y = steeldemo$production, b = steeldemo$emissions, v = 0.01467,
#'   id = steeldemo$plant
#' )
#' \donttest{
#' boot_pgt_sensitivity(tech, model = "wgd", B = 50, seed = 1)
#' }
#' @export
boot_pgt_sensitivity <- function(tech, model = "wgd", B = 100,
                                 m_grid = NULL, level = 0.95,
                                 returns = "vrs", peers = "all",
                                 pollutant = 1L, seed = NULL) {
  L <- tech$L
  if (is.null(m_grid)) {
    m_grid <- unique(round(L^seq(0.5, 0.9, length.out = 5)))
    m_grid <- pmin(pmax(m_grid, 2L), L)
  }
  do.call(rbind, lapply(m_grid, function(mm) {
    bt <- boot_pgt(tech, model = model, B = B, m = mm, level = level,
                   returns = returns, peers = peers, pollutant = pollutant,
                   seed = seed)
    w <- bt$per_dmu$upper - bt$per_dmu$lower
    data.frame(m = mm, median_width = stats::median(w, na.rm = TRUE),
               mean_se = mean(bt$per_dmu$se, na.rm = TRUE))
  }))
}

#' @export
print.pgt_boot <- function(x, ...) {
  cat(sprintf("pgt subsampling inference: model = %s, B = %d, m = %d, level = %.2f\n",
              x$model, x$B, x$m, x$level))
  cat(sprintf("  DMUs: %d\n", nrow(x$per_dmu)))
  w <- x$per_dmu$upper - x$per_dmu$lower
  cat(sprintf("  median interval width: %.4f   mean SE: %.4f\n",
              stats::median(w, na.rm = TRUE),
              mean(x$per_dmu$se, na.rm = TRUE)))
  if (!is.null(x$group_means)) {
    cat("\nGroup means:\n")
    print(x$group_means, row.names = FALSE, digits = 4)
  }
  cat("\n  intervals are heuristic; see ?boot_pgt.\n")
  invisible(x)
}

#' @export
as.data.frame.pgt_boot <- function(x, ...) {
  x$per_dmu
}
