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
#'   and \code{violated}. When the technology records an abatement
#'   output, also \code{a} and \code{closure} (\eqn{gap - a_l}): the
#'   equality residual \eqn{u'x_l - v y_l - b_l - a_l}, zero when the
#'   account closes exactly, which \code{model = "fdmo"} requires. The
#'   number of violations is attached as attribute
#'   \code{"n_violations"}.
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
    if (!is.null(tech$a)) {
      d$a <- tech$a[, p]
      d$closure <- gap - tech$a[, p]
    }
    d <- .insert_group(d, tech$group)
    if (P > 1L) {
      lead <- if (is.null(tech$group)) "id" else c("id", "group")
      d <- cbind(d[lead], pollutant = tech$pollutants[p],
                 d[setdiff(names(d), lead)])
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
  # counts are recomputed from the printed rows, so a subset of the
  # audit prints counts that describe the subset, not the full table
  L <- nrow(x)
  cat("Materials-balance audit (u'x - v y >= b)\n")
  if (L == 0) {
    cat("  no accounts\n")
    return(invisible(x))
  }
  nv <- sum(x$violated)
  tol <- attr(x, "tol")
  cat(sprintf("  accounts: %d   violations: %d (%.1f%%)   tolerance: %s\n",
              L, nv, 100 * nv / L,
              if (is.null(tol)) "unknown" else format(tol)))
  cat(sprintf("  closure gap (gap / potential): min %.4f, median %.4f, max %.4f\n",
              min(x$rel_gap), stats::median(x$rel_gap), max(x$rel_gap)))
  ne <- sum(x$gap < 0)
  cat(sprintf("  accounts with gap < 0 (exact): %d; any infeasible weak-G LPs are confined to these\n",
              ne))
  if (!is.null(x$closure)) {
    cat(sprintf("  equality closure (gap - a): largest |closure| / potential = %.2g (model = \"fdmo\" requires exact closure)\n",
                max(abs(x$closure) / pmax(x$potential,
                                          .Machine$double.eps))))
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
