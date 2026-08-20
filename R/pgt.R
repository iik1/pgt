#' Estimate a pollution-generating technology model
#'
#' Fits a nonparametric materials-balance efficiency model by solving one
#' linear program per DMU. The minimising models return the least bad
#' output \eqn{b^*_l} attainable at the DMU's activity level, with
#' environmental efficiency \eqn{b^*_l / b_l}; the directional
#' model returns the joint good-output expansion and bad-output
#' contraction.
#'
#' The models fall into two axiom families. The materials-balance models
#' (\code{"wgd"}, \code{"wgd_input_fixed"}, \code{"envelope"},
#' \code{"fdmo"}, \code{"mb_cost"}) enforce the identity
#' \eqn{u'x_l - v'y_l \ge b_l}; the reference models (\code{"byprod"},
#' \code{"wd"}) implement competing systems for cross-comparison (see
#' [compare_models()]).
#' \describe{
#'   \item{\code{"wgd"}}{The weak-G-disposability model of Rodseth
#'     (2025), Eq. 6 in reduced form. Equation 6 fixes only the
#'     intended outputs at the evaluated DMU's levels; the inputs are
#'     decision variables. Substituting the weak-G summing-up condition
#'     on the disposal slacks therefore collapses the programme to
#'     minimising the peer emission envelope plus the retained content
#'     of the good output produced beyond the DMU's own levels,
#'     \deqn{b^*_i = \min_\lambda \sum_l \lambda_l b_l +
#'       v_i'\left(\sum_l \lambda_l y_l - y_i\right),}
#'     over the output rows \eqn{\sum_l \lambda_l y_{ml} \ge y_{mi}}
#'     and the returns-to-scale row, valued at the DMU's own
#'     coefficients \eqn{v_i}. Self-reference is always feasible, so
#'     every score lies in \eqn{(0, 1]}. When the abatement output
#'     \code{a} is observed, the optimal peer mix also delivers the
#'     implied uncontrolled emission \code{z_star}
#'     \eqn{= \sum_l \lambda_l z_l} and abatement \code{a_star}
#'     \eqn{= \sum_l \lambda_l a_l}, the Table 2 quantities of Rodseth
#'     (2025). \code{model = "wgd_rodseth"} is an alias.}
#'   \item{\code{"wgd_input_fixed"}}{The input-fixed companion benchmark
#'     (the model estimated as \code{"wgd"} before version 0.6.0):
#'     minimise the bad output subject to output feasibility, the input
#'     rows \eqn{\sum_l \lambda_l x_{nl} \le x_{ni}}, the peer emission
#'     envelope, and the evaluated DMU's own materials-balance cap
#'     \eqn{b \le u_i'x_i - v_i'y_i}. With several pollutants the caps
#'     of the remaining pollutants constrain the peer mix as well, so
#'     the projection respects every pollutant's materials balance.
#'     DMUs whose data violate a materials-balance identity lose the
#'     self-reference guarantee: their LP may be infeasible (returning
#'     \code{NA}), and every infeasible LP belongs to such a DMU. A DMU
#'     violating another pollutant's identity can also score above 1,
#'     since its MB-consistent projection may emit more of the selected
#'     pollutant than the DMU reports. Audit with [mb_check()] first.}
#'   \item{\code{"envelope"}}{The \code{"wgd"} programme without the
#'     retained-content term (its \eqn{v = 0} case): the convex lower
#'     envelope of the \eqn{(y, b)} scatter. Self-reference is always
#'     feasible, so all scores lie in \eqn{(0, 1]}.}
#'   \item{\code{"fdmo"}}{The factorially determined multi-output
#'     (directional) representation of Rodseth (2025), Eq. 13 with the
#'     abatement output fixed at each DMU's own level. It jointly
#'     maximises the expansion of the good output (\eqn{\theta_y}) and
#'     the contraction of the bad output (\eqn{\theta_b}) along the
#'     materials-balance frontier; gross inefficiency is
#'     \eqn{\theta_y + \theta_b}, with 0 for a frontier DMU. The gross
#'     score adds good-output units to bad-output units (direction
#'     \eqn{(1, 1)}), so it is unit-dependent and comparable only
#'     within one unit system, as in Rodseth's Table 3; rescaling
#'     \code{y} or \code{b} changes the direction's implicit weighting
#'     and can change the optimiser, not merely the reported number. Requires an
#'     abatement output \code{a} in [pgt_tech()] and a single intended
#'     output. The model imposes the
#'     materials-balance identity as an exact equality
#'     \eqn{u'x_l - v y_l = b_l + a_l}: accounts that do not close
#'     exactly either make the LP infeasible or shift the closure gap
#'     into the scores, and \code{pgt()} warns when it detects open
#'     accounts (absolute closure gap \eqn{|u'x_l - v y_l - b_l - a_l|}
#'     exceeding \eqn{10^{-6}} times the pollutant potential
#'     \eqn{u'x_l}). \code{model = "ddf"} is an alias.}
#'   \item{\code{"mb_cost"}}{The materials-balance cost model of Coelli,
#'     Lauwers and Van Huylenbroeck (2007): environmental efficiency is
#'     the ratio of minimal to observed aggregate material inflow
#'     \eqn{u'x}, decomposed as \eqn{EE = TE \times EAE} into technical
#'     and environmental-allocative efficiency (returned as \code{te}
#'     and \code{eae}). A unit that uses none of the pollutant-bearing
#'     inputs has an undefined inflow ratio: its scores are \code{NA}
#'     with \code{status} 2, reported without a solver call.}
#'   \item{\code{"byprod"}}{The by-production intersection technology of
#'     Murty, Russell and Levkoff (2012): the good output is produced by
#'     one sub-technology and the bad output generated by another over
#'     the emission-causing inputs. The implemented model is the
#'     original Murty-Russell-Levkoff intersection technology without
#'     integrated abatement; Hampf (2014) develops the
#'     abatement-integrated extension. Returns emission efficiency
#'     (\eqn{b^*/b}, reported as \code{efficiency}), output
#'     efficiency (\code{output_eff})
#'     and \code{mean_eff}, the arithmetic mean of the two sub-efficiencies:
#'     a package summary in the spirit of, but not identical to, the
#'     Fare-Grosskopf-Lovell graph measure, and not a measure defined
#'     in the source. Set the
#'     emission-causing inputs with \code{polluting} in [pgt_tech()].}
#'   \item{\code{"wd"}}{The weak-disposability model (Kuosmanen 2005
#'     correct VRS formulation; see the exchange settled in Kuosmanen
#'     and Podinovski 2009), a reference axiom system that imposes no
#'     materials balance. Returns emission efficiency \eqn{b^*/b}.}
#' }
#'
#' Several intended outputs are supported by every model except
#' \code{"fdmo"}: supply \code{y} as an \eqn{L \times M} matrix in
#' [pgt_tech()]. The output rows then hold all \eqn{M} outputs at the
#' evaluated DMU's levels, and the retained content generalises to
#' \eqn{v_i'y_i} with one coefficient per output and pollutant.
#'
#' @param tech A [pgt_tech()] object.
#' @param model Character: \code{"wgd"} (alias \code{"wgd_rodseth"}),
#'   \code{"wgd_input_fixed"}, \code{"envelope"}, \code{"fdmo"} (alias
#'   \code{"ddf"}), \code{"mb_cost"}, \code{"byprod"} or \code{"wd"}.
#'   See Details.
#' @param returns Returns to scale: \code{"vrs"} (default) or
#'   \code{"crs"}.
#' @param peers Reference set: \code{"all"} pools every DMU;
#'   \code{"group"} restricts each DMU's peers to its own technology
#'   group (requires \code{group} in [pgt_tech()]).
#' @param pollutant For a multi-pollutant technology, the pollutant to
#'   estimate: a column name or index of \code{b}. Defaults to the first
#'   pollutant. This selects which bad output the objective acts on. For
#'   \code{model = "wgd_input_fixed"} the materials-balance caps of every
#'   pollutant constrain the projection; the other models use only the
#'   selected pollutant's data.
#'
#' @return An object of class \code{"pgt"}: a list with
#'   \describe{
#'     \item{\code{results}}{Data frame with one row per DMU. The good
#'       outputs appear as a single column \code{y} when \eqn{M = 1}
#'       and as one named column per output otherwise, with the output
#'       duals named accordingly (\code{dual_output}, or
#'       \code{dual_<name>} per output). For \code{"wgd"}: \code{id},
#'       \code{group} (if present), the output columns, \code{b},
#'       \code{b_star}, \code{efficiency} (\eqn{b^*/b}), the output
#'       duals (total derivatives
#'       \eqn{\partial b^*/\partial y_{mi} = \mu_m - v_{mi}}, which can
#'       be negative when the retained-content coefficient is large),
#'       \code{z_star} and \code{a_star} (only when \code{a} is
#'       observed) and \code{status}, the \code{lp_solve} solver code
#'       (0 = solved; 2 = infeasible; 5 = numerically failed). For
#'       \code{"wgd_input_fixed"}, \code{"envelope"} and \code{"wd"}: the
#'       same score columns with non-negative output duals (\code{NA}
#'       for \code{"wd"}, which reports no duals) and
#'       \code{mb_headroom}, the cap slack
#'       \eqn{u_i'x_i - v_i'y_i - b^*_i} (\code{NA} except for
#'       \code{"wgd_input_fixed"}). For \code{"fdmo"}: \code{id},
#'       \code{group} (if present), \code{y}, \code{b}, \code{gross}
#'       (\eqn{\theta_y + \theta_b}), \code{good_eff} (\eqn{\theta_y}),
#'       \code{bad_eff} (\eqn{\theta_b}), \code{maximal_y}
#'       (\eqn{y_l + \theta_y}) and \code{status}.}
#'     \item{\code{weights}}{List of named vectors of non-zero peer
#'       weights per DMU, keyed by \code{make.unique(id)}.}
#'     \item{\code{model}, \code{returns}, \code{peers},
#'       \code{pollutant}}{The call settings.}
#'   }
#'
#' @references
#' Coelli, T., Lauwers, L., & Van Huylenbroeck, G. (2007).
#' Environmental efficiency measurement and the materials balance
#' condition. \emph{Journal of Productivity Analysis}, 28(1--2), 3--12.
#' \doi{10.1007/s11123-007-0052-8}
#'
#' Fare, R., Grosskopf, S., & Lovell, C. A. K. (1985). \emph{The
#' Measurement of Efficiency of Production}. Kluwer-Nijhoff, Boston.
#' \doi{10.1007/978-94-015-7721-2}
#'
#' Hampf, B. (2014). Separating environmental efficiency into production
#' and abatement efficiency: A nonparametric model with application to
#' US power plants. \emph{Journal of Productivity Analysis}, 41(3),
#' 457--473. \doi{10.1007/s11123-013-0357-8}
#'
#' Kuosmanen, T. (2005). Weak disposability in nonparametric production
#' analysis with undesirable outputs. \emph{American Journal of
#' Agricultural Economics}, 87(4), 1077--1082.
#' \doi{10.1111/j.1467-8276.2005.00788.x}
#'
#' Kuosmanen, T., & Podinovski, V. V. (2009). Weak disposability in
#' nonparametric production analysis: Reply to Fare and Grosskopf.
#' \emph{American Journal of Agricultural Economics}, 91(2), 539--545.
#' \doi{10.1111/j.1467-8276.2008.01238.x}
#'
#' Murty, S., Russell, R. R., & Levkoff, S. B. (2012). On modeling
#' pollution-generating technologies. \emph{Journal of Environmental
#' Economics and Management}, 64(1), 117--135.
#' \doi{10.1016/j.jeem.2012.02.005}
#'
#' Rodseth, K. L. (2025). On the development of a unified, nonparametric
#' materials balance-based efficiency analysis model and its
#' applications. \emph{Journal of Productivity Analysis}, 64(3),
#' 305--319. \doi{10.1007/s11123-025-00768-0}
#'
#' @seealso [pgt_tech()], [pgt_decompose()], [shadow_prices()],
#'   [mac_curve()]
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
#' summary(fit)
#' @export
pgt <- function(tech, model = c("wgd", "wgd_rodseth", "wgd_input_fixed",
                                "envelope", "fdmo", "ddf", "byprod",
                                "mb_cost", "wd"),
                returns = c("vrs", "crs"), peers = c("all", "group"),
                pollutant = 1L) {
  stopifnot(inherits(tech, "pgt_tech"))
  model <- match.arg(model)
  if (model == "wgd_rodseth") model <- "wgd"
  if (model == "ddf") model <- "fdmo"
  returns <- match.arg(returns)
  peers <- match.arg(peers)
  vrs <- returns == "vrs"
  p <- .pollutant_index(tech, pollutant)

  if (model == "fdmo") {
    if (tech$M > 1L) {
      stop("model = \"fdmo\" is defined for a single intended output; ",
           "the direction of Rodseth (2025, Eq. 13) has no ",
           "multi-output form.", call. = FALSE)
    }
    if (is.null(tech$a)) {
      stop("model = \"fdmo\" requires an abatement output 'a' in ",
           "pgt_tech().", call. = FALSE)
    }
    closure <- .mb_cap(tech, p) - tech$b[, p] - tech$a[, p]
    rel <- abs(closure) / pmax(.mb_potential(tech, p),
                               .Machine$double.eps)
    n_open <- sum(rel > 1e-6)
    if (n_open > 0) {
      warning(sprintf(paste0(
        "%d of %d DMUs do not close the materials-balance identity ",
        "exactly (u'x - v y != b + a; largest relative gap %.2g). ",
        "model = \"fdmo\" imposes the identity as an equality, so open ",
        "accounts either make the LP infeasible or shift the closure ",
        "gap into the scores."), n_open, tech$L, max(rel)), call. = FALSE)
    }
  }

  peer_sets <- .peer_sets(tech, peers)
  L <- tech$L
  ctx <- .solve_ctx(tech, model, p)
  weights <- vector("list", L)
  names(weights) <- make.unique(tech$id)
  status <- rep(NA_integer_, L)
  sols <- vector("list", L)

  for (i in seq_len(L)) {
    ps <- peer_sets[[i]]
    sol <- .lp_solve_one(model, i, tech, ps, vrs, p = p, ctx = ctx)
    sols[[i]] <- sol
    status[i] <- sol$status
    w <- if (model == "byprod") sol$mu else sol$lambda
    if (!is.null(w)) {
      weights[[i]] <- .named_weights(w, tech, ps)
    }
  }
  num <- function(field) {
    vapply(sols, function(s) {
      val <- s[[field]]
      if (is.null(val)) NA_real_ else as.numeric(val)
    }, numeric(1))
  }
  num_m <- function(field, m) {
    vapply(sols, function(s) {
      val <- s[[field]]
      if (is.null(val) || length(val) < m) NA_real_ else
        as.numeric(val[m])
    }, numeric(1))
  }
  b_p <- ctx$b_p

  # one intended-output column when M = 1, one per named output else;
  # duals follow the same naming
  y_cols <- if (tech$M == 1L) {
    stats::setNames(data.frame(tech$y[, 1L]), "y")
  } else {
    stats::setNames(as.data.frame(tech$y), tech$outputs)
  }
  dual_cols <- function() {
    if (tech$M == 1L) {
      stats::setNames(data.frame(num_m("dual_output", 1L)),
                      "dual_output")
    } else {
      stats::setNames(
        as.data.frame(lapply(seq_len(tech$M), function(m)
          num_m("dual_output", m))),
        paste0("dual_", tech$outputs))
    }
  }

  results <- switch(model,
    fdmo = cbind(
      data.frame(id = tech$id, stringsAsFactors = FALSE), y_cols,
      data.frame(
        b = b_p, gross = num("gross"), good_eff = num("theta_y"),
        bad_eff = num("theta_b"),
        maximal_y = tech$y[, 1L] + num("theta_y"),
        status = status, stringsAsFactors = FALSE
      )
    ),
    byprod = cbind(
      data.frame(id = tech$id, stringsAsFactors = FALSE), y_cols,
      data.frame(
        b = b_p, b_star = num("b_star"),
        efficiency = num("emission_eff"),
        output_eff = num("output_eff"), mean_eff = num("mean_eff"),
        status = status, stringsAsFactors = FALSE
      )
    ),
    mb_cost = cbind(
      data.frame(id = tech$id, stringsAsFactors = FALSE), y_cols,
      data.frame(
        b = b_p, b_star = num("b_star"), efficiency = num("mbe"),
        te = num("te"), eae = num("eae"), status = status,
        stringsAsFactors = FALSE
      )
    ),
    wgd = {
      out <- cbind(
        data.frame(id = tech$id, stringsAsFactors = FALSE), y_cols,
        data.frame(
          b = b_p, b_star = num("b_star"),
          efficiency = num("b_star") / b_p, stringsAsFactors = FALSE
        ),
        dual_cols()
      )
      if (!is.null(tech$a)) {
        out$z_star <- num("z_star")
        out$a_star <- num("a_star")
      }
      out$status <- status
      out
    },
    cbind(
      data.frame(id = tech$id, stringsAsFactors = FALSE), y_cols,
      data.frame(
        b = b_p, b_star = num("b_star"),
        efficiency = num("b_star") / b_p, stringsAsFactors = FALSE
      ),
      dual_cols(),
      data.frame(
        mb_headroom = num("mb_rhs") - num("b_star"),
        status = status, stringsAsFactors = FALSE
      )
    )
  )

  n_failed <- sum(status != 0)
  if (n_failed > 0) {
    hint <- switch(model,
      wgd_input_fixed = paste0(" For model = \"wgd_input_fixed\" this flags ",
                            "DMUs violating a materials-balance ",
                            "identity; run mb_check()."),
      fdmo = paste0(" For model = \"fdmo\" this flags DMUs whose ",
                    "accounts do not close exactly; run mb_check()."),
      "")
    warning(sprintf(
      "%d of %d LPs infeasible or failed; their scores are NA.%s",
      n_failed, L, hint), call. = FALSE)
  }
  if (model == "mb_cost") {
    n_neg <- sum(results$b_star < 0, na.rm = TRUE)
    if (n_neg > 0) {
      warning(sprintf(paste0(
        "%d of %d DMUs have a negative implied minimal emission b_star ",
        "(the minimised material inflow lies below the DMU's retained ",
        "content, which positive retained-content coefficients permit); ",
        "interpret b_star with care."), n_neg, L), call. = FALSE)
    }
  }
  results <- .insert_group(results, tech$group)

  structure(
    list(results = results, weights = weights, model = model,
         returns = returns, peers = peers,
         pollutant = tech$pollutants[p]),
    class = "pgt"
  )
}

# Resolve a pollutant selector (name or index) to a column position.
.pollutant_index <- function(tech, pollutant) {
  if (is.character(pollutant)) {
    if (length(pollutant) != 1L) {
      stop("'pollutant' must be a single column name or integer index.",
           call. = FALSE)
    }
    p <- match(pollutant, tech$pollutants)
    if (is.na(p)) {
      stop("pollutant '", pollutant, "' not found; available: ",
           paste(tech$pollutants, collapse = ", "), ".", call. = FALSE)
    }
  } else {
    if (!is.numeric(pollutant) || length(pollutant) != 1L ||
        is.na(pollutant) || pollutant != round(pollutant)) {
      stop("'pollutant' must be a single column name or integer index.",
           call. = FALSE)
    }
    p <- as.integer(pollutant)
    if (p < 1L || p > tech$P) {
      stop("'pollutant' must index one column of 'b' (1..", tech$P, ").",
           call. = FALSE)
    }
  }
  p
}

# Non-zero peer weights, named by DMU id.
.named_weights <- function(lambda, tech, peers) {
  nz <- which(lambda > 1e-8)
  if (!length(nz)) return(NULL)
  w <- lambda[nz]
  names(w) <- tech$id[peers[nz]]
  w
}

# Insert the group column directly after the id column (shared by the
# result frames of pgt(), boot_pgt(), mac_curve() and mb_check()).
.insert_group <- function(d, group) {
  if (is.null(group)) return(d)
  cbind(d[1], group = group, d[-1])
}

# TRUE for the directional model, whose score columns differ.
.is_directional <- function(x) identical(x$model, "fdmo")

# Header lines shared by print.pgt and print.summary.pgt.
.print_pgt_header <- function(x) {
  cat(sprintf("pgt fit: model = %s, returns = %s, peers = %s\n",
              x$model, x$returns, x$peers))
  if (length(x$pollutant)) {
    cat(sprintf("  pollutant: %s\n", x$pollutant))
  }
}

#' @export
print.pgt <- function(x, ...) {
  r <- x$results
  .print_pgt_header(x)
  if (.is_directional(x)) {
    cat(sprintf("  DMUs: %d   solved: %d   gross inefficiency: median %.3f\n",
                nrow(r), sum(r$status == 0),
                stats::median(r$gross, na.rm = TRUE)))
  } else {
    cat(sprintf("  DMUs: %d   solved: %d   efficiency (b*/b): median %.3f\n",
                nrow(r), sum(r$status == 0),
                stats::median(r$efficiency, na.rm = TRUE)))
  }
  invisible(x)
}

#' @export
summary.pgt <- function(object, ...) {
  r <- object$results
  qs <- c(0, 0.25, 0.5, 0.75, 1)
  directional <- .is_directional(object)
  score <- if (directional) r$gross else r$efficiency
  overall <- stats::quantile(score, qs, na.rm = TRUE)
  by_group <- NULL
  if (!is.null(r$group)) {
    by_group <- do.call(rbind, lapply(split(r, r$group), function(d) {
      s <- if (directional) d$gross else d$efficiency
      out <- data.frame(
        group = d$group[1], n = nrow(d), n_na = sum(is.na(s)),
        median = stats::median(s, na.rm = TRUE),
        mean = mean(s, na.rm = TRUE)
      )
      if (directional) {
        out$median_good <- stats::median(d$good_eff, na.rm = TRUE)
        out$median_bad <- stats::median(d$bad_eff, na.rm = TRUE)
      } else if (!is.null(d$dual_output)) {
        out$median_dual_output <- stats::median(d$dual_output, na.rm = TRUE)
      }
      out
    }))
    rownames(by_group) <- NULL
  }
  structure(
    list(model = object$model, returns = object$returns,
         peers = object$peers, pollutant = object$pollutant,
         directional = directional, n = nrow(r),
         n_failed = sum(r$status != 0),
         overall = overall, by_group = by_group),
    class = "summary.pgt"
  )
}

#' @export
print.summary.pgt <- function(x, ...) {
  .print_pgt_header(x)
  cat(sprintf("  DMUs: %d   failed LPs: %d\n", x$n, x$n_failed))
  cat(sprintf("\n%s:\n",
              if (x$directional) "Gross inefficiency (theta_y + theta_b)"
              else "Efficiency (b*/b)"))
  print(round(x$overall, 4))
  if (!is.null(x$by_group)) {
    cat("\nBy group:\n")
    print(x$by_group, row.names = FALSE, digits = 4)
  }
  invisible(x)
}

#' @export
plot.pgt <- function(x, ...) {
  r <- x$results
  directional <- .is_directional(x)
  score <- if (directional) r$gross else r$efficiency
  lab <- if (directional) "Gross inefficiency (theta_y + theta_b)"
         else "Environmental efficiency (b*/b)"
  if (all(is.na(score))) {
    stop("all scores are NA (all LPs failed); run mb_check().",
         call. = FALSE)
  }
  if (!is.null(r$group)) {
    graphics::boxplot(score ~ r$group, ylab = lab, xlab = "", ...)
  } else {
    graphics::hist(score, main = "", xlab = lab, ...)
  }
  invisible(x)
}

#' @export
as.data.frame.pgt <- function(x, ...) {
  x$results
}
