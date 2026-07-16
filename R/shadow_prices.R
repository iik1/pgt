#' Extract shadow values from a pgt fit
#'
#' Returns the informative constraint diagnostics of the per-DMU linear
#' programs. \code{dual_output} is the shadow value of the output
#' constraint, \eqn{\partial b^* / \partial y \ge 0}: the marginal
#' bad-output content of the good output along the frontier (for CO2,
#' the marginal emission intensity in tonnes of CO2 per tonne of
#' output). \code{mb_headroom} (weak-G-disposability model only) is the
#' slack \eqn{u'x_l - v y_l - b^*_l} between the projection and the
#' DMU's materials-balance ceiling.
#'
#' The dual of the peer emission envelope
#' \eqn{\sum_l \lambda_l b_l \le b} is not reported: in this reduced
#' form it equals \eqn{-1} whenever the LP solves (dual feasibility on
#' the emission variable, with the materials-balance row slack at any
#' optimum), so it carries no cross-DMU information. In the boundary
#' case \eqn{b^* = u'x_l - v y_l} exactly, the split of this dual
#' between the envelope and cap rows is basis-dependent. A
#' materials-balance cap that would bind manifests as an infeasible LP,
#' not as a capped projection.
#'
#' Values are reported in the units of the linear program (quantities,
#' not money). To monetise the output dual, combine it with an output
#' price via [mac_curve()].
#'
#' @param fit A [pgt()] fit with \code{model = "wgd"} or
#'   \code{"envelope"}, the models that report output duals.
#'
#' @return A data frame with one row per DMU: \code{id}, \code{group}
#'   (if present), \code{b}, \code{b_star}, \code{dual_output} and
#'   \code{mb_headroom} (\code{NA} for the envelope model).
#'
#' @references
#' Fare, R., Grosskopf, S., Lovell, C. A. K., & Yaisawarng, S. (1993).
#' Derivation of shadow prices for undesirable outputs: A distance
#' function approach. \emph{The Review of Economics and Statistics},
#' 75(2), 374--380. \doi{10.2307/2109448}
#'
#' Rodseth, K. L. (2025). On the development of a unified, nonparametric
#' materials balance-based efficiency analysis model and its
#' applications. \emph{Journal of Productivity Analysis}, 64(3),
#' 305--319. \doi{10.1007/s11123-025-00768-0}
#'
#' @seealso [pgt()], [mac_curve()]
#' @examples
#' data(steeldemo)
#' tech <- pgt_tech(
#'   x = steeldemo[, c("coal_coke", "other_fuel", "raw_material", "flux")],
#'   y = steeldemo$production,
#'   b = steeldemo$emissions,
#'   v = 0.01467,
#'   group = steeldemo$route,
#'   id = steeldemo$plant
#' )
#' fit <- pgt(tech, model = "wgd")
#' head(shadow_prices(fit))
#' @export
shadow_prices <- function(fit) {
  stopifnot(inherits(fit, "pgt"))
  if (!fit$model %in% c("wgd", "envelope")) {
    stop("shadow_prices() requires a fit with output duals ",
         "(model = \"wgd\" or \"envelope\"); got model = \"",
         fit$model, "\".", call. = FALSE)
  }
  keep <- intersect(c("id", "group", "b", "b_star",
                      "dual_output", "mb_headroom"),
                    names(fit$results))
  fit$results[keep]
}

#' Marginal abatement cost curve
#'
#' Builds a marginal abatement cost (MAC) curve from a [pgt()] fit.
#' Along the estimated frontier, reducing the bad output by one unit
#' requires forgoing \eqn{1 / \eta_l} units of the good output, where
#' \eqn{\eta_l = \partial b^* / \partial y} is the output-constraint
#' dual. The DMU's marginal abatement cost is therefore
#' \eqn{p / \eta_l} for an output price \eqn{p}; with no price supplied
#' the curve is expressed in output units per unit of bad output
#' (\eqn{1 / \eta_l}). Each DMU's abatement potential is its distance to
#' the frontier, \eqn{b_l - b^*_l}. DMUs are ordered by increasing
#' marginal cost and the potential is accumulated.
#'
#' The curve therefore combines two margins. It orders DMUs by the
#' shadow price at their frontier projection and plots against it the
#' abatement available from eliminating inefficiency (the gap
#' \eqn{b_l - b^*_l}), which the model itself prices at zero output
#' loss; the \code{mac} value is the marginal cost of abatement beyond
#' the frontier point. The area under the curve is consequently not a
#' total-cost estimate.
#'
#' DMUs with \eqn{\eta_l = 0} are excluded: a zero output dual arises
#' when the output constraint is slack, i.e. the DMU projects onto the
#' flat segment of the \eqn{(y, b)} frontier, where abatement via
#' output contraction has locally infinite marginal cost. Only
#' exact-zero duals are excluded. Unsolved LPs are excluded as well.
#' The number of excluded DMUs and the summed abatement potential of
#' the solved-but-excluded ones (which is NOT part of the curve) are
#' attached as attributes \code{"n_excluded"} and
#' \code{"excluded_abatement"} and reported by \code{print}.
#'
#' A DMU that projects onto a vertex of the piecewise-linear frontier
#' has a range of optimal duals, and the solver reports the dual of
#' whichever optimal basis it terminates in, so the MAC of such a DMU
#' is one point of an interval and can differ across solver versions or
#' platforms. Frontier-interior projections have unique duals.
#'
#' @param fit A [pgt()] fit with \code{model = "wgd"} or
#'   \code{"envelope"}, the models that report output duals.
#' @param price Optional scalar price of the good output. If supplied,
#'   \code{mac} is in money per unit of bad output.
#'
#' @return A data frame of class \code{"pgt_mac"}, ordered by
#'   \code{mac}: \code{id}, \code{group} (if present), \code{mac},
#'   \code{abatement} (\eqn{b - b^*}) and \code{cum_abatement}.
#'
#' @references
#' Fare, R., Grosskopf, S., Lovell, C. A. K., & Yaisawarng, S. (1993).
#' Derivation of shadow prices for undesirable outputs: A distance
#' function approach. \emph{The Review of Economics and Statistics},
#' 75(2), 374--380. \doi{10.2307/2109448}
#'
#' Rodseth, K. L. (2025). On the development of a unified, nonparametric
#' materials balance-based efficiency analysis model and its
#' applications. \emph{Journal of Productivity Analysis}, 64(3),
#' 305--319. \doi{10.1007/s11123-025-00768-0}
#'
#' @seealso [pgt()], [shadow_prices()]
#' @examples
#' data(steeldemo)
#' tech <- pgt_tech(
#'   x = steeldemo[, c("coal_coke", "other_fuel", "raw_material", "flux")],
#'   y = steeldemo$production,
#'   b = steeldemo$emissions,
#'   v = 0.01467,
#'   group = steeldemo$route,
#'   id = steeldemo$plant
#' )
#' fit <- pgt(tech, model = "wgd")
#' mac <- mac_curve(fit, price = 550)
#' plot(mac)
#' @export
mac_curve <- function(fit, price = NULL) {
  stopifnot(inherits(fit, "pgt"))
  if (!fit$model %in% c("wgd", "envelope")) {
    stop("mac_curve() requires a fit with output duals ",
         "(model = \"wgd\" or \"envelope\"); got model = \"",
         fit$model, "\".", call. = FALSE)
  }
  if (!is.null(price)) {
    stopifnot(is.numeric(price), length(price) == 1L, price > 0)
  }
  r <- fit$results
  ok <- r$status == 0 & !is.na(r$dual_output) & r$dual_output > 0
  n_excluded <- sum(!ok)
  solved_excluded <- !ok & r$status == 0
  excluded_abatement <- sum(pmax(r$b[solved_excluded] -
                                   r$b_star[solved_excluded], 0))
  d <- r[ok, , drop = FALSE]

  unit_cost <- 1 / d$dual_output
  mac <- if (is.null(price)) unit_cost else price * unit_cost

  out <- data.frame(
    id = d$id,
    mac = mac,
    # solver noise can leave b* a few ulps above b; potential is >= 0
    abatement = pmax(d$b - d$b_star, 0),
    stringsAsFactors = FALSE
  )
  out <- .insert_group(out, d$group)
  out <- out[order(out$mac), , drop = FALSE]
  out$cum_abatement <- cumsum(out$abatement)
  rownames(out) <- NULL

  class(out) <- c("pgt_mac", "data.frame")
  attr(out, "price") <- price
  attr(out, "n_excluded") <- n_excluded
  attr(out, "excluded_abatement") <- excluded_abatement
  out
}

#' @export
print.pgt_mac <- function(x, ...) {
  price <- attr(x, "price")
  cat("Marginal abatement cost curve\n")
  if (nrow(x) == 0) {
    cat(sprintf("  no DMUs with a positive output dual (excluded: %d)\n",
                attr(x, "n_excluded")))
    return(invisible(x))
  }
  cat(sprintf("  DMUs: %d (excluded: %d)   total abatement potential: %.4g\n",
              nrow(x), attr(x, "n_excluded"), max(x$cum_abatement)))
  ex_ab <- attr(x, "excluded_abatement")
  if (!is.null(ex_ab) && ex_ab > 0) {
    cat(sprintf("  abatement potential of excluded solved DMUs (not in curve): %.4g\n",
                ex_ab))
  }
  cat(sprintf("  mac units: %s\n",
              if (is.null(price)) "good output per unit of bad output"
              else sprintf("money per unit of bad output (price = %g)", price)))
  cat(sprintf("  mac quartiles: %.4g / %.4g / %.4g\n",
              stats::quantile(x$mac, 0.25), stats::median(x$mac),
              stats::quantile(x$mac, 0.75)))
  invisible(x)
}

#' @export
plot.pgt_mac <- function(x, ...) {
  if (nrow(x) == 0) {
    stop("empty MAC curve: no DMUs with a positive output dual.",
         call. = FALSE)
  }
  graphics::plot(x$cum_abatement, x$mac, type = "s",
                 xlab = "Cumulative abatement potential (units of b)",
                 ylab = if (is.null(attr(x, "price")))
                   "Marginal abatement cost (units of y per b)"
                 else "Marginal abatement cost (money per unit of b)",
                 ...)
  invisible(x)
}
