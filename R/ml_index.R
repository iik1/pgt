#' Global Malmquist-Luenberger productivity index
#'
#' Computes the global Malmquist-Luenberger (GML) productivity index of
#' Oh (2010) for a panel of pollution-generating technologies, decomposed
#' into efficiency change and best-practice (frontier) change. The index
#' is built on the directional distance function of the pollution-
#' generating technology, which expands the good output and contracts the
#' bad output jointly, so productivity change accounts for emissions
#' rather than ignoring them.
#'
#' The directional distance to a reference technology \eqn{R} is
#' \deqn{D_R(x, y, b) = \max\{\beta : (y + \beta g_y, b - \beta g_b)
#'   \textrm{ feasible in } R \textrm{ given } x\},}
#' with direction \eqn{g = (y, b)} (each observation is scaled by its own
#' level). The global reference pools every period into one technology
#' (Oh 2010), which removes the cross-period infeasibility of the
#' adjacent-period Malmquist-Luenberger index (Aparicio, Pastor and Zofio
#' 2013) and makes the index circular. For a DMU observed in periods
#' \eqn{t} and \eqn{t+1},
#' \deqn{GML = \frac{1 + D_G(t)}{1 + D_G(t+1)}
#'   = \underbrace{\frac{1 + D_t(t)}{1 + D_{t+1}(t+1)}}_{EC}
#'   \times \underbrace{GML / EC}_{BPC},}
#' where \eqn{D_G} is the distance to the global frontier and \eqn{D_t}
#' the distance to period \eqn{t}'s own frontier. \eqn{GML > 1} is
#' productivity growth; \eqn{EC} is catch-up to the frontier and
#' \eqn{BPC} is frontier movement.
#'
#' Three reference technologies are available. \code{"wd"} imposes weak
#' disposability of the bad output (the Kuosmanen 2005 VRS form, with
#' the intensity weights split into an active and an abatement part):
#' this is the technology under which Chung, Fare and Grosskopf (1997)
#' and Oh (2010) define the (global) Malmquist-Luenberger index, so
#' \code{technology = "wd"} is the faithful Oh (2010) comparator.
#' \code{"wgd"} (the default) and \code{"envelope"} instead treat the
#' bad output as reducible down to the peer emission envelope with no
#' proportional output sacrifice, the same lower-envelope treatment as
#' the corresponding [pgt()] programs, with (\code{"wgd"}) and without
#' (\code{"envelope"}) the input constraints. Under constant returns
#' the free-disposal technologies contain the weak-disposability set,
#' so their distances are weakly larger than under \code{"wd"}; under
#' variable returns the sets are not nested in general, because the
#' Kuosmanen form lets the active intensity weights sum below one.
#'
#' This estimator is experimental: the global Malmquist-Luenberger
#' index under a materials-balance technology is not yet settled in the
#' literature. Under \code{"wgd"} and \code{"envelope"} the index does
#' not re-impose the per-DMU materials-balance cap on the projected
#' point, and the bad output is freely disposable toward the envelope
#' rather than weakly disposable, so these variants are exploratory
#' companions to the \code{"wd"} index, not implementations of Oh
#' (2010).
#'
#' @param tech A [pgt_tech()] object with a non-\code{NULL} \code{period}
#'   and \code{id} identifying the same DMU across periods.
#' @param technology Reference technology: \code{"wgd"} keeps the input
#'   constraints and the lower emission envelope (default);
#'   \code{"envelope"} frees the inputs, so the frontier becomes the
#'   \eqn{(y, b)} envelope; \code{"wd"} imposes weak disposability of
#'   the bad output, the technology of Oh (2010). See Details.
#' @param returns Returns to scale: \code{"vrs"} (default) or
#'   \code{"crs"}.
#' @param pollutant For a multi-pollutant technology, the pollutant to
#'   index. Defaults to the first.
#'
#' @return An object of class \code{"pgt_ml"}: a list with
#'   \code{results}, a data frame with one row per DMU per consecutive
#'   period transition (\code{id}, \code{from}, \code{to}, \code{gml},
#'   \code{ec}, \code{bpc}), and the call settings. Transitions are
#'   formed only for DMUs present in both adjacent periods.
#'
#' @references
#' Chambers, R. G., Chung, Y., & Fare, R. (1996). Benefit and distance
#' functions. \emph{Journal of Economic Theory}, 70(2), 407--419.
#' \doi{10.1006/jeth.1996.0096}
#'
#' Chung, Y. H., Fare, R., & Grosskopf, S. (1997). Productivity and
#' undesirable outputs: A directional distance function approach.
#' \emph{Journal of Environmental Management}, 51(3), 229--240.
#' \doi{10.1006/jema.1997.0146}
#'
#' Oh, D.-h. (2010). A global Malmquist-Luenberger productivity index.
#' \emph{Journal of Productivity Analysis}, 34(3), 183--197.
#' \doi{10.1007/s11123-010-0178-y}
#'
#' Aparicio, J., Pastor, J. T., & Zofio, J. L. (2013). On the
#' inconsistency of the Malmquist-Luenberger index. \emph{European
#' Journal of Operational Research}, 229(3), 738--742.
#' \doi{10.1016/j.ejor.2013.03.031}
#'
#' @seealso [pgt()], [pgt_tech()]
#' @examples
#' # Two-period panel: period 2's frontier emits less for the same output.
#' set.seed(1)
#' n <- 12
#' d1 <- data.frame(id = 1:n, y = runif(n, 5, 10), x = runif(n, 8, 12),
#'                  b = runif(n, 2, 6))
#' d2 <- data.frame(id = 1:n, y = d1$y * 1.05, x = d1$x, b = d1$b * 0.85)
#' d <- rbind(cbind(d1, period = 1), cbind(d2, period = 2))
#' tech <- pgt_tech(x = d[, "x", drop = FALSE], y = d$y, b = d$b,
#'                  period = d$period, id = d$id)
#' ml <- pgt_ml(tech)
#' summary(ml)
#' @export
pgt_ml <- function(tech, technology = c("wgd", "envelope", "wd"),
                   returns = c("vrs", "crs"), pollutant = 1L) {
  stopifnot(inherits(tech, "pgt_tech"))
  technology <- match.arg(technology)
  returns <- match.arg(returns)
  vrs <- returns == "vrs"
  p <- .pollutant_index(tech, pollutant)
  if (is.null(tech$period)) {
    stop("pgt_ml() requires a 'period' in pgt_tech().", call. = FALSE)
  }
  use_inputs <- technology == "wgd"

  b_p <- tech$b[, p]
  periods <- levels(tech$period)
  global <- seq_len(tech$L)
  by_period <- split(global, tech$period)

  dd <- function(k, ref) {
    if (technology == "wd") {
      .ddf_ml_wd(k, ref, tech$x, tech$y, b_p, vrs)
    } else {
      .ddf_ml(k, ref, tech$x, tech$y, b_p, use_inputs, vrs)
    }
  }
  # global and contemporaneous directional distances for every row
  dg <- vapply(global, function(k) dd(k, global), numeric(1))
  dc <- vapply(global, function(k) {
    dd(k, by_period[[as.character(tech$period[k])]])
  }, numeric(1))
  n_failed <- sum(is.na(dg)) + sum(is.na(dc))
  if (n_failed > 0) {
    warning(sprintf(paste0(
      "%d distance LPs failed; the transitions of the affected ",
      "observations are NA."), n_failed), call. = FALSE)
  }

  rows <- list()
  for (j in seq_len(length(periods) - 1L)) {
    from <- periods[j]; to <- periods[j + 1L]
    it <- by_period[[from]]
    it1 <- by_period[[to]]
    common <- intersect(tech$id[it], tech$id[it1])
    for (idv in common) {
      k <- it[match(idv, tech$id[it])]
      k1 <- it1[match(idv, tech$id[it1])]
      gml <- (1 + dg[k]) / (1 + dg[k1])
      ec <- (1 + dc[k]) / (1 + dc[k1])
      rows[[length(rows) + 1L]] <- data.frame(
        id = idv, from = from, to = to,
        gml = gml, ec = ec, bpc = gml / ec,
        stringsAsFactors = FALSE
      )
    }
  }
  results <- if (length(rows)) do.call(rbind, rows) else
    data.frame(id = character(0), from = character(0), to = character(0),
               gml = numeric(0), ec = numeric(0), bpc = numeric(0))
  rownames(results) <- NULL

  structure(
    list(results = results, technology = technology, returns = returns,
         pollutant = tech$pollutants[p]),
    class = "pgt_ml"
  )
}

# Directional distance of observation k to reference set `ref`, direction
# g = (y_k, b_k): max beta s.t. sum_ref lam y >= y_k(1+beta),
# [sum_ref lam x <= x_k], sum_ref lam b <= b_k(1-beta), sum lam = 1 [VRS].
.ddf_ml <- function(k, ref, X, y, b, use_inputs, vrs) {
  L <- length(ref)
  N <- ncol(X)
  nc <- L + 1L
  ib <- nc
  lp <- lpSolveAPI::make.lp(nrow = 0, ncol = nc)
  invisible(lpSolveAPI::lp.control(lp, sense = "max"))
  o <- numeric(nc); o[ib] <- 1; lpSolveAPI::set.objfn(lp, o)
  r <- numeric(nc); r[seq_len(L)] <- y[ref]; r[ib] <- -y[k]
  lpSolveAPI::add.constraint(lp, r, ">=", y[k])
  if (use_inputs) {
    for (n in seq_len(N)) {
      r <- numeric(nc); r[seq_len(L)] <- X[ref, n]
      lpSolveAPI::add.constraint(lp, r, "<=", X[k, n])
    }
  }
  r <- numeric(nc); r[seq_len(L)] <- b[ref]; r[ib] <- b[k]
  lpSolveAPI::add.constraint(lp, r, "<=", b[k])
  if (vrs) {
    r <- numeric(nc); r[seq_len(L)] <- 1
    lpSolveAPI::add.constraint(lp, r, "=", 1)
  }
  lpSolveAPI::set.bounds(lp, lower = c(rep(0, L), -Inf))
  if (lpSolveAPI::solve.lpExtPtr(lp) != 0) return(NA_real_)
  lpSolveAPI::get.variables(lp)[ib]
}

# Weak-disposability directional distance (Kuosmanen 2005 VRS form):
# the technology under which Chung, Fare and Grosskopf (1997) and Oh
# (2010) define the (global) Malmquist-Luenberger index. Intensity
# weights split into an active part z (scales goods and bads together)
# and an abatement part w (scales inputs only); direction g = (y_k, b_k).
#
#   max beta
#   s.t. sum_ref z y          >= y_k (1 + beta)
#        sum_ref z b           = b_k (1 - beta)   [weak disposability]
#        sum_ref (z + w) x_n  <= x_kn  (n = 1..N)
#        sum_ref (z + w)       = 1               [VRS]
#        z, w >= 0, beta free
.ddf_ml_wd <- function(k, ref, X, y, b, vrs) {
  L <- length(ref)
  N <- ncol(X)
  nc <- 2L * L + 1L
  ib <- nc
  lp <- lpSolveAPI::make.lp(nrow = 0, ncol = nc)
  invisible(lpSolveAPI::lp.control(lp, sense = "max"))
  o <- numeric(nc); o[ib] <- 1; lpSolveAPI::set.objfn(lp, o)
  r <- numeric(nc); r[seq_len(L)] <- y[ref]; r[ib] <- -y[k]
  lpSolveAPI::add.constraint(lp, r, ">=", y[k])
  r <- numeric(nc); r[seq_len(L)] <- b[ref]; r[ib] <- b[k]
  lpSolveAPI::add.constraint(lp, r, "=", b[k])
  for (n in seq_len(N)) {
    r <- numeric(nc); r[seq_len(L)] <- X[ref, n]
    r[(L + 1):(2L * L)] <- X[ref, n]
    lpSolveAPI::add.constraint(lp, r, "<=", X[k, n])
  }
  if (vrs) {
    r <- numeric(nc); r[seq_len(2L * L)] <- 1
    lpSolveAPI::add.constraint(lp, r, "=", 1)
  }
  lpSolveAPI::set.bounds(lp, lower = c(rep(0, 2L * L), -Inf))
  if (lpSolveAPI::solve.lpExtPtr(lp) != 0) return(NA_real_)
  lpSolveAPI::get.variables(lp)[ib]
}

# Header line shared by print.pgt_ml and print.summary.pgt_ml.
.print_ml_header <- function(x) {
  cat(sprintf("pgt global Malmquist-Luenberger index (technology = %s, returns = %s)\n",
              x$technology, x$returns))
}

#' @export
print.pgt_ml <- function(x, ...) {
  r <- x$results
  .print_ml_header(x)
  cat(sprintf("  transitions: %d   DMUs: %d\n",
              nrow(r), length(unique(r$id))))
  if (nrow(r)) {
    n_na <- sum(is.na(r$gml))
    if (n_na > 0) {
      cat(sprintf("  transitions with NA index (failed LPs): %d\n", n_na))
    }
    cat(sprintf("  GML median %.4f   EC median %.4f   BPC median %.4f\n",
                stats::median(r$gml, na.rm = TRUE),
                stats::median(r$ec, na.rm = TRUE),
                stats::median(r$bpc, na.rm = TRUE)))
  }
  invisible(x)
}

#' @export
summary.pgt_ml <- function(object, ...) {
  r <- object$results
  qs <- c(0, 0.25, 0.5, 0.75, 1)
  tab <- rbind(
    GML = stats::quantile(r$gml, qs, na.rm = TRUE),
    EC = stats::quantile(r$ec, qs, na.rm = TRUE),
    BPC = stats::quantile(r$bpc, qs, na.rm = TRUE)
  )
  by_period <- NULL
  if (nrow(r)) {
    by_period <- do.call(rbind, lapply(split(r, r$to), function(d) {
      data.frame(to = d$to[1], n = nrow(d),
                 gml = exp(mean(log(d$gml), na.rm = TRUE)),
                 ec = exp(mean(log(d$ec), na.rm = TRUE)),
                 bpc = exp(mean(log(d$bpc), na.rm = TRUE)))
    }))
    rownames(by_period) <- NULL
  }
  structure(list(quantiles = tab, by_period = by_period,
                 technology = object$technology, returns = object$returns),
            class = "summary.pgt_ml")
}

#' @export
print.summary.pgt_ml <- function(x, ...) {
  .print_ml_header(x)
  cat("\nDistribution (GML > 1 is productivity growth):\n")
  print(round(x$quantiles, 4))
  if (!is.null(x$by_period)) {
    cat("\nGeometric mean by period transition:\n")
    print(x$by_period, row.names = FALSE, digits = 4)
  }
  invisible(x)
}

#' @export
as.data.frame.pgt_ml <- function(x, ...) {
  x$results
}
