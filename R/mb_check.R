#' Audit the materials-balance identity
#'
#' Pre-estimation feasibility audit of the materials-balance condition
#' \eqn{u'x_l - v y_l \ge b_l} for every DMU (and every pollutant). DMUs
#' that violate the identity carry inconsistent material accounts (more
#' pollutant leaves than enters). A violation makes self-reference
#' infeasible in the weak-G-disposability program, so that DMU's LP
#' solves only if some peer mix meets every constraint within its
#' materials-balance cap; conversely, every infeasible
#' weak-G-disposability LP belongs to a DMU with an exact violation,
#' \code{gap < 0} (satisfying DMUs always solve via self-reference).
#' Audit before estimation.
#'
#' @param tech A [pgt_tech()] object.
#' @param tol Relative tolerance: a DMU-pollutant account is flagged
#'   (\code{violated}) when its closure gap is below \code{-tol} times
#'   its pollutant potential \eqn{u'x_l}. The flag is a data-quality
#'   screen; LP infeasibility tracks the exact sign of the gap, so the
#'   flagged set is a subset of the accounts with \code{gap < 0}
#'   (counted separately in attribute \code{"n_exact"}).
#'
#' @return A data frame of class \code{"pgt_mb"} with one row per DMU
#'   (per pollutant when several are present): \code{id}, \code{group}
#'   (if present), \code{pollutant} (when several), \code{potential}
#'   (\eqn{u'x_l}), \code{retained} (\eqn{v y_l}), \code{b}, \code{gap}
#'   (\eqn{u'x_l - v y_l - b_l}), \code{rel_gap} (\code{gap / potential})
#'   and \code{violated}. The number of violations is attached as
#'   attribute \code{"n_violations"}.
#'
#' @seealso [pgt_tech()], [pgt()]
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
#' mb <- mb_check(tech)
#' mb
#' @export
mb_check <- function(tech, tol = 1e-8) {
  stopifnot(inherits(tech, "pgt_tech"))
  P <- tech$P
  parts <- lapply(seq_len(P), function(p) {
    potential <- .mb_potential(tech, p)
    retained <- tech$v[, p] * tech$y
    gap <- potential - retained - tech$b[, p]
    rel_gap <- gap / pmax(potential, .Machine$double.eps)
    d <- data.frame(
      id = tech$id, potential = potential, retained = retained,
      b = tech$b[, p], gap = gap, rel_gap = rel_gap,
      violated = rel_gap < -tol, stringsAsFactors = FALSE
    )
    if (!is.null(tech$group)) {
      d <- cbind(d[1], group = tech$group, d[-1])
    }
    if (P > 1L) {
      d <- cbind(d[if (is.null(tech$group)) 1 else 1:2],
                 pollutant = tech$pollutants[p],
                 d[-(if (is.null(tech$group)) 1 else 1:2)])
    }
    d
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  class(out) <- c("pgt_mb", "data.frame")
  attr(out, "tol") <- tol
  attr(out, "n_violations") <- sum(out$violated)
  attr(out, "n_exact") <- sum(out$gap < 0)
  attr(out, "P") <- P
  out
}

#' @export
print.pgt_mb <- function(x, ...) {
  nv <- attr(x, "n_violations")
  L <- nrow(x)
  cat("Materials-balance audit (u'x - v y >= b)\n")
  cat(sprintf("  accounts: %d   violations: %d (%.1f%%)   tolerance: %g\n",
              L, nv, 100 * nv / L, attr(x, "tol")))
  cat(sprintf("  closure gap (gap / potential): min %.4f, median %.4f, max %.4f\n",
              min(x$rel_gap), stats::median(x$rel_gap), max(x$rel_gap)))
  ne <- attr(x, "n_exact")
  if (!is.null(ne)) {
    cat(sprintf("  accounts with gap < 0 (exact): %d; any infeasible weak-G LPs are confined to these\n",
                ne))
  }
  if (nv > 0) {
    cat("\n  worst violations:\n")
    bad <- x[x$violated, , drop = FALSE]
    bad <- bad[order(bad$rel_gap), , drop = FALSE]
    print.data.frame(utils::head(bad, 5), row.names = FALSE, digits = 4)
    cat("\n  inspect the material accounts of the flagged DMUs.\n")
  }
  invisible(x)
}
