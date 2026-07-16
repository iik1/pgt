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
#' The intervals follow the subsampling recipe of Politis and Romano
#' (1994): the law of \eqn{L^\kappa(\hat\theta - \theta)} is
#' approximated by the recentred, rate-scaled subsample distribution
#' \eqn{m^\kappa(\theta^*_m - \hat\theta)}, giving
#' \deqn{[\hat\theta - (m/L)^\kappa q_{1-\alpha/2}(\theta^*_m -
#'   \hat\theta),\; \hat\theta - (m/L)^\kappa q_{\alpha/2}(\theta^*_m -
#'   \hat\theta)].}
#' Each DMU is evaluated against the subsample frontier alone: the
#' evaluated unit is not added to the reference set, so frontier units
#' carry genuine resampling variation instead of a degenerate
#' zero-width interval. Replicates in which a DMU's program is
#' infeasible are dropped as \code{NA} and counted: \code{per_dmu$n_ok}
#' reports the number of feasible replicates behind each interval, and
#' a warning is issued when it falls below half of \code{B}. Because a
#' subsample frontier lies weakly inside the full-sample frontier, the
#' deviations are non-negative and the interval extends downward from
#' the point estimate, matching the direction of the frontier bias; the
#' point estimate itself typically lies at or above its own upper
#' endpoint, since the interval targets the true score, which the
#' estimate overestimates. The reported \code{se} is the subsample
#' standard deviation rescaled by \eqn{(m/L)^\kappa}.
#'
#' With \code{peers = "group"} the subsample is stratified: each group
#' contributes in proportion to its size, with a floor of two units, so
#' groups smaller than the floor are included whole and effectively not
#' resampled.
#'
#' Group-mean intervals are a descriptive companion, not a supported
#' inference: the replicate statistic averages all group members
#' against a common subsample frontier, so replicate variation reflects
#' frontier movement only, and the per-DMU rate exponent is reused
#' although means of frontier estimators converge at their own rates
#' (see Kneip, Simar and Wilson 2015 for the aggregate theory). Read
#' them as a stability band around the group mean, not as a confidence
#' interval with nominal coverage.
#'
#' The intervals remain heuristic. Subsampling theory for data
#' envelopment analysis backs the radial and directional-distance
#' technologies, not the materials-balance-constrained programs
#' implemented here, so \code{kappa} borrows the DEA convergence rate
#' as an approximation. Use the subsample-size sensitivity of
#' [boot_pgt_sensitivity()] to check that the intervals are stable in
#' \code{m} (its implied \code{kappa_hat} makes the borrowed rate
#' checkable) before relying on them.
#'
#' @param tech A [pgt_tech()] object.
#' @param model One of the environmental-efficiency models
#'   (\code{"wgd"}, \code{"envelope"}, \code{"byprod"}, \code{"mb_cost"},
#'   \code{"wd"}); the directional model \code{"fdmo"} is not supported.
#' @param B Number of subsampling replicates.
#' @param m Subsample size. Defaults to \code{round(L^0.7)}, a point in
#'   the interior of the admissible range; Simar and Wilson (2011) show
#'   that sensitivity to \code{m} is the central practical issue, so
#'   check the default with [boot_pgt_sensitivity()].
#' @param level Confidence level for the intervals.
#' @param kappa Convergence-rate exponent used to rescale the subsample
#'   distribution. Defaults to the DEA rate \eqn{2/(d+1)} under VRS and
#'   \eqn{2/d} under CRS (Kneip, Simar and Wilson 2008) for the model's
#'   effective frontier dimension \eqn{d}: \eqn{N+2} for \code{"wgd"}
#'   and \code{"wd"} (inputs, good and bad output), \eqn{N+1} for
#'   \code{"mb_cost"} and \code{"byprod"} (whose programs live in the
#'   \eqn{(x, y)} and sub-technology spaces), and \eqn{2} for the
#'   input-free \code{"envelope"} model, whose frontier lives in the
#'   \eqn{(y, b)} plane. A heuristic; supply your own exponent to
#'   override.
#' @param returns,peers,pollutant Passed to the underlying model, as in
#'   [pgt()].
#' @param seed Optional integer seed for reproducibility. The caller's
#'   RNG state is restored on exit.
#'
#' @return An object of class \code{"pgt_boot"}: a list with
#'   \code{per_dmu} (a data frame of point estimate, lower and upper
#'   bound, rescaled subsampling standard error, and \code{n_ok}, the
#'   number of feasible replicates, per DMU), \code{group_means} (the
#'   same for group means, when a group is present), and the call
#'   settings.
#'
#' @references
#' Kneip, A., Simar, L., & Wilson, P. W. (2008). Asymptotics and
#' consistent bootstraps for DEA estimators in nonparametric frontier
#' models. \emph{Econometric Theory}, 24(6), 1663--1697.
#' \doi{10.1017/S0266466608080651}
#'
#' Kneip, A., Simar, L., & Wilson, P. W. (2015). When bias kills the
#' variance: Central limit theorems for DEA and FDH efficiency scores.
#' \emph{Econometric Theory}, 31(2), 394--422.
#' \doi{10.1017/S0266466614000413}
#'
#' Politis, D. N., & Romano, J. P. (1994). Large sample confidence
#' regions based on subsamples under minimal assumptions. \emph{The
#' Annals of Statistics}, 22(4), 2031--2050.
#' \doi{10.1214/aos/1176325770}
#'
#' Politis, D. N., Romano, J. P., & Wolf, M. (1999).
#' \emph{Subsampling}. Springer, New York.
#' \doi{10.1007/978-1-4612-1554-7}
#'
#' Simar, L., & Wilson, P. W. (2011). Inference by the m out of n
#' bootstrap in nonparametric frontier models. \emph{Journal of
#' Productivity Analysis}, 36(1), 33--53.
#' \doi{10.1007/s11123-010-0200-4}
#'
#' @seealso [pgt()], [boot_pgt_sensitivity()]
#' @examples
#' data(steeldemo)
#' steel60 <- steeldemo[1:60, ]
#' tech <- pgt_tech(
#'   x = steel60[, c("coal_coke", "other_fuel", "raw_material", "flux")],
#'   y = steel60$production, b = steel60$emissions, v = 0.01467,
#'   group = steel60$route, id = steel60$plant
#' )
#' \donttest{
#' bt <- boot_pgt(tech, model = "wgd", B = 40, seed = 1)
#' head(bt$per_dmu)
#' }
#' @export
boot_pgt <- function(tech, model = c("wgd", "envelope", "byprod",
                                     "mb_cost", "wd"),
                     B = 200, m = NULL, level = 0.95, kappa = NULL,
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
  if (is.null(kappa)) {
    # effective frontier dimension: wgd/wd use (x, y, b); mb_cost lives
    # in (x, y); byprod's headline score comes from the T2 sub-LP over
    # the polluting inputs and b; the input-free envelope lives in the
    # (y, b) plane
    dims <- switch(model,
      envelope = 2L,
      mb_cost = tech$N + 1L,
      byprod = length(.polluting_inputs(tech, p)) + 1L,
      tech$N + 2L)
    kappa <- if (vrs) 2 / (dims + 1) else 2 / dims
  }
  stopifnot(is.numeric(kappa), length(kappa) == 1L, kappa > 0)
  if (!is.null(tech$period)) {
    warning("the technology carries a 'period': boot_pgt() subsamples ",
            "rows as independent cross-sectional draws and ignores ",
            "within-unit dependence across periods; subset to one ",
            "period or interpret with care.", call. = FALSE)
  }
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = globalenv(),
                           inherits = FALSE)) {
      get(".Random.seed", envir = globalenv(), inherits = FALSE)
    } else NULL
    on.exit(
      if (is.null(old_seed)) {
        rm(".Random.seed", envir = globalenv())
      } else {
        assign(".Random.seed", old_seed, envir = globalenv())
      },
      add = TRUE
    )
    set.seed(seed)
  }

  ctx <- .solve_ctx(tech, model, p)
  point <- .model_efficiency(tech, model, vrs, p,
                             .peer_sets(tech, peers), ctx)
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
    # each DMU is evaluated against the subsample frontier alone; the
    # evaluated unit is NOT appended, so frontier units carry genuine
    # resampling variation and infeasible replicates surface as NA
    peer_sets <- if (do_group) {
      lapply(seq_len(L), function(i) sub[grp[sub] == grp[i]])
    } else {
      rep(list(sub), L)
    }
    boot[, rrep] <- .model_efficiency(tech, model, vrs, p, peer_sets, ctx)
  }

  # Politis-Romano subsampling: recentre the subsample scores at the
  # point estimate and rescale by the relative convergence rate. Both
  # endpoints are truncated to [0, max(1, estimate)] AFTER the quantile
  # flip so the truncation preserves lower <= upper and does not force
  # the upper bound below a legitimate above-1 estimate (wgd scores of
  # materials-balance violators).
  rate <- (m / L)^kappa
  a <- (1 - level) / 2
  qdev <- function(z, p) {
    z <- z[!is.na(z)]
    if (!length(z)) return(NA_real_)
    stats::quantile(z, p, names = FALSE)
  }
  dev <- boot - point
  n_ok <- rowSums(!is.na(boot))
  hi <- pmax(1, point)
  clamp <- function(z, top) pmin(pmax(z, 0), top)
  per_dmu <- data.frame(
    id = tech$id,
    estimate = point,
    lower = clamp(point - rate *
                    apply(dev, 1, qdev, p = 1 - a), hi),
    upper = clamp(point - rate *
                    apply(dev, 1, qdev, p = a), hi),
    se = rate * apply(boot, 1, stats::sd, na.rm = TRUE),
    n_ok = n_ok,
    stringsAsFactors = FALSE
  )
  per_dmu <- .insert_group(per_dmu, grp)
  rownames(per_dmu) <- NULL
  n_thin <- sum(n_ok < B / 2)
  if (n_thin > 0) {
    warning(sprintf(paste0(
      "%d of %d DMUs have infeasible LPs in more than half of their ",
      "subsample replicates; their intervals rest on few replicates ",
      "(see per_dmu$n_ok)."), n_thin, L), call. = FALSE)
  }

  group_means <- NULL
  if (!is.null(grp)) {
    glev <- levels(grp)
    est_g <- vapply(glev, function(g) mean(point[grp == g], na.rm = TRUE),
                    numeric(1))
    gm_boot <- vapply(glev, function(g) {
      apply(boot[grp == g, , drop = FALSE], 2,
            function(z) if (all(is.na(z))) NA_real_ else
              mean(z, na.rm = TRUE))
    }, numeric(B))
    gm_dev <- sweep(gm_boot, 2, est_g)
    hi_g <- pmax(1, est_g)
    group_means <- data.frame(
      group = glev,
      estimate = est_g,
      lower = clamp(est_g - rate *
                      apply(gm_dev, 2, qdev, p = 1 - a), hi_g),
      upper = clamp(est_g - rate *
                      apply(gm_dev, 2, qdev, p = a), hi_g),
      se = rate * apply(gm_boot, 2, stats::sd, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    rownames(group_means) <- NULL
  }

  structure(
    list(per_dmu = per_dmu, group_means = group_means, model = model,
         B = B, m = m, level = level, kappa = kappa, returns = returns,
         peers = peers),
    class = "pgt_boot"
  )
}

# Environmental efficiency of every DMU under `model`, evaluated against
# the supplied per-DMU peer sets. The headline-score definition is
# shared with pgt() through .headline_score().
.model_efficiency <- function(tech, model, vrs, p, peer_sets, ctx = NULL) {
  if (is.null(ctx)) ctx <- .solve_ctx(tech, model, p)
  b_i <- ctx$b_p
  vapply(seq_len(tech$L), function(i) {
    sol <- .lp_solve_one(model, i, tech, peer_sets[[i]], vrs, p = p,
                         ctx = ctx)
    .headline_score(model, sol, b_i[i])
  }, numeric(1))
}

#' Subsample-size sensitivity of subsampling intervals
#'
#' Re-runs [boot_pgt()] over a grid of subsample sizes \code{m} and
#' reports how the median interval width responds. Stable widths across
#' the grid support the choice of \code{m}; a strong trend warns that the
#' intervals are driven by the resampling size rather than the data.
#'
#' When at least three grid points have positive widths, the attribute
#' \code{"kappa_hat"} reports the convergence-rate exponent implied by
#' the log-spread regression of Politis, Romano and Wolf (1999): the
#' rescaled width scales as \eqn{m^{\kappa - \kappa_0}} for true rate
#' \eqn{\kappa_0} and working rate \eqn{\kappa}, so
#' \eqn{\hat\kappa_0 = \kappa - } slope of \eqn{\log} width on
#' \eqn{\log m}. Agreement between \code{kappa_hat} and the working
#' \code{kappa} (attribute \code{"kappa_used"}) supports the borrowed
#' DEA rate; disagreement suggests supplying \code{kappa_hat} to
#' [boot_pgt()]. The same \code{seed} is reset before every grid point,
#' so replicate draws are comparable across \code{m}.
#'
#' @param tech A [pgt_tech()] object.
#' @param model,B,level,kappa,returns,peers,pollutant,seed As in
#'   [boot_pgt()].
#' @param m_grid Integer vector of subsample sizes. Defaults to five
#'   sizes spanning \eqn{L^{0.5}} to \eqn{L^{0.9}}.
#'
#' @return A data frame with one row per \code{m}: the size, the median
#'   per-DMU interval width, the median over strictly positive widths,
#'   and the mean rescaled standard error, with attributes
#'   \code{"kappa_used"} and (when estimable) \code{"kappa_hat"}.
#'
#' @seealso [boot_pgt()]
#' @examples
#' data(steeldemo)
#' steel60 <- steeldemo[1:60, ]
#' tech <- pgt_tech(
#'   x = steel60[, c("coal_coke", "other_fuel", "raw_material", "flux")],
#'   y = steel60$production, b = steel60$emissions, v = 0.01467,
#'   id = steel60$plant
#' )
#' \donttest{
#' boot_pgt_sensitivity(tech, model = "wgd", B = 25,
#'                      m_grid = c(15, 25, 40), seed = 1)
#' }
#' @export
boot_pgt_sensitivity <- function(tech, model = "wgd", B = 100,
                                 m_grid = NULL, level = 0.95,
                                 kappa = NULL,
                                 returns = "vrs", peers = "all",
                                 pollutant = 1L, seed = NULL) {
  L <- tech$L
  if (is.null(m_grid)) {
    m_grid <- unique(round(L^seq(0.5, 0.9, length.out = 5)))
    m_grid <- pmin(pmax(m_grid, 2L), L)
  }
  kappa_used <- NA_real_
  out <- do.call(rbind, lapply(m_grid, function(mm) {
    bt <- boot_pgt(tech, model = model, B = B, m = mm, level = level,
                   kappa = kappa, returns = returns, peers = peers,
                   pollutant = pollutant, seed = seed)
    kappa_used <<- bt$kappa
    w <- bt$per_dmu$upper - bt$per_dmu$lower
    data.frame(m = mm, median_width = stats::median(w, na.rm = TRUE),
               median_width_pos = stats::median(w[w > 0], na.rm = TRUE),
               mean_se = mean(bt$per_dmu$se, na.rm = TRUE))
  }))
  attr(out, "kappa_used") <- kappa_used
  ok <- is.finite(out$median_width_pos) & out$median_width_pos > 0
  if (sum(ok) >= 3L) {
    slope <- stats::coef(stats::lm(log(out$median_width_pos[ok]) ~
                                     log(out$m[ok])))[[2L]]
    attr(out, "kappa_hat") <- kappa_used - slope
  }
  out
}

#' @export
print.pgt_boot <- function(x, ...) {
  cat(sprintf("pgt subsampling inference: model = %s, B = %d, m = %d, level = %.2f, kappa = %.3g\n",
              x$model, x$B, x$m, x$level, x$kappa))
  cat(sprintf("  DMUs: %d\n", nrow(x$per_dmu)))
  w <- x$per_dmu$upper - x$per_dmu$lower
  cat(sprintf("  median interval width: %.4f   mean SE: %.4f\n",
              stats::median(w, na.rm = TRUE),
              mean(x$per_dmu$se, na.rm = TRUE)))
  n_thin <- sum(x$per_dmu$n_ok < x$B / 2)
  if (n_thin > 0) {
    cat(sprintf("  DMUs with fewer than B/2 feasible replicates: %d (see per_dmu$n_ok)\n",
                n_thin))
  }
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
