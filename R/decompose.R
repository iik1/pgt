#' Decompose environmental efficiency across technology groups
#'
#' Metafrontier decompositions of the environmental efficiency score
#' \eqn{b^*/b} into within-group and technology-gap components, built on
#' Rodseth's (2025) weak-G-disposability technology and the DEA
#' metafrontier and technology-gap ratio of O'Donnell, Rao and Battese
#' (2008); the metafrontier idea originates with Battese, Rao and
#' O'Donnell (2004).
#'
#' \describe{
#'   \item{\code{"envelope"}}{The \eqn{(y, b)}-envelope decomposition
#'     built on Rodseth's (2025) Eq. 6 with inputs free. Stages differ by
#'     peer set only, so the weak-G-disposability slack equality is
#'     respected throughout:
#'     \deqn{Total = \frac{b^*_{all}}{b} = WR \times TGR, \quad
#'           WR = \frac{b^*_{group}}{b}, \quad
#'           TGR = \frac{b^*_{all}}{b^*_{group}}.}
#'     \eqn{WR} is within-group reallocation efficiency; \eqn{TGR} is the
#'     technology-gap ratio (O'Donnell, Rao and Battese 2008).
#'     Self-reference is always feasible, so every component lies in
#'     \eqn{(0, 1]} with no feasibility screen.}
#'   \item{\code{"rodseth"}}{A three-stage decomposition of the full
#'     weak-G-disposability model (Rodseth 2025, Eq. 11, collapsed to
#'     three components when there are no dedicated abatement inputs):
#'     \deqn{EnvEff = \frac{b^*_3}{b}
#'       = \underbrace{\frac{b^*_1}{b}}_{TE}
#'         \times \underbrace{\frac{b^*_2}{b^*_1}}_{Technology}
#'         \times \underbrace{\frac{b^*_3}{b^*_2}}_{AE},}
#'     where stage 1 solves the WGD program against own-group peers
#'     (technical efficiency), stage 2 against all peers (technology
#'     gap), and stage 3 additionally drops the input constraints
#'     (input-mix / allocative component). Each stage is a strict
#'     relaxation of the previous one, so every component lies in
#'     \eqn{(0, 1]} where feasible for DMUs satisfying every
#'     pollutant's materials-balance identity; as in [pgt()], a DMU
#'     violating another pollutant's identity can have \code{te} and
#'     \code{total} above 1, since its MB-consistent projection may
#'     emit more of the selected pollutant than the DMU reports. When
#'     an early stage is infeasible
#'     (materials-balance violators, see [mb_check()]) its component is
#'     \code{NA} while later-stage ratios and \code{total} can remain
#'     defined, matching the reference implementation; the
#'     multiplicative identity holds for rows with all components
#'     present.}
#' }
#'
#' @param tech A [pgt_tech()] object with a non-\code{NULL} \code{group}.
#' @param type \code{"envelope"} (default) or \code{"rodseth"}. See
#'   Details.
#' @param returns Returns to scale: \code{"vrs"} (default) or
#'   \code{"crs"}.
#' @param pollutant For a multi-pollutant technology, the pollutant to
#'   decompose: a column name or index of \code{b}. Defaults to the
#'   first.
#'
#' @return An object of class \code{"pgt_decomp"}: a list with
#'   \code{results} (one row per DMU with the stage minima and the
#'   multiplicative components), \code{type} and \code{returns}. The
#'   component columns are \code{WR}, \code{TGR}, \code{total} for
#'   \code{type = "envelope"} and \code{te}, \code{technology},
#'   \code{ae}, \code{total} for \code{type = "rodseth"}.
#'
#' @references
#' Battese, G. E., Rao, D. S. P., & O'Donnell, C. J. (2004). A
#' metafrontier production function for estimation of technical
#' efficiencies and technology gaps for firms operating under different
#' technologies. \emph{Journal of Productivity Analysis}, 21(1), 91--103.
#' \doi{10.1023/B:PROD.0000012454.06094.29}
#'
#' O'Donnell, C. J., Rao, D. S. P., & Battese, G. E. (2008).
#' Metafrontier frameworks for the study of firm-level efficiencies and
#' technology ratios. \emph{Empirical Economics}, 34(2), 231--255.
#' \doi{10.1007/s00181-007-0119-4}
#'
#' Rodseth, K. L. (2025). On the development of a unified, nonparametric
#' materials balance-based efficiency analysis model and its
#' applications. \emph{Journal of Productivity Analysis}, 64(3),
#' 305--319. \doi{10.1007/s11123-025-00768-0}
#'
#' @seealso [pgt()], [pgt_tech()]
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
#' dec <- pgt_decompose(tech, type = "envelope")
#' summary(dec)
#' @export
pgt_decompose <- function(tech, type = c("envelope", "rodseth"),
                          returns = c("vrs", "crs"), pollutant = 1L) {
  stopifnot(inherits(tech, "pgt_tech"))
  type <- match.arg(type)
  returns <- match.arg(returns)
  if (type == "envelope" && is.null(tech$group)) {
    stop("type = \"envelope\" requires a 'group' in pgt_tech().",
         call. = FALSE)
  }
  vrs <- returns == "vrs"
  L <- tech$L
  p <- .pollutant_index(tech, pollutant)
  b_p <- tech$b[, p]
  group_sets <- if (!is.null(tech$group)) .peer_sets(tech, "group")
  all_peers <- seq_len(L)
  if (type == "envelope") {
    ctx <- .solve_ctx(tech, "envelope", p)
    b_group <- b_all <- rep(NA_real_, L)
    for (i in seq_len(L)) {
      b_group[i] <- .lp_solve_one("envelope", i, tech, group_sets[[i]],
                                  vrs, p = p, ctx = ctx)$b_star
      b_all[i] <- .lp_solve_one("envelope", i, tech, all_peers,
                                vrs, p = p, ctx = ctx)$b_star
    }
    results <- data.frame(
      id = tech$id,
      group = tech$group,
      b = b_p,
      b_star_group = b_group,
      b_star_all = b_all,
      WR = b_group / b_p,
      TGR = b_all / b_group,
      total = b_all / b_p,
      stringsAsFactors = FALSE
    )
  } else {
    # Rodseth (2025, Eq. 11): stage programmes of the extended
    # representation (Eq. 9), successively freeing quality, the
    # production inputs, the abatement output and the abatement
    # inputs. Components collapse to 1 when the data cannot separate
    # them: quality requires producer-specific coefficients, the
    # abatement stages require an observed abatement output, and the
    # last stage requires dedicated pollution-control inputs.
    has_quality <- !(.u_is_uniform(tech, p) &&
                       all(apply(matrix(tech$v[, , p], L, tech$M), 2L,
                                 function(col) diff(range(col)) == 0)))
    has_a <- !is.null(tech$a)
    has_xa <- length(tech$x_abate) > 0L
    stage <- function(i, ps, hold_xp, hold_xa, hold_a, hold_q) {
      .lp_wgd_stage(i, tech, ps, vrs, p = p, hold_xp = hold_xp,
                    hold_xa = hold_xa, hold_a = hold_a,
                    hold_quality = hold_q)$b_star
    }
    # the abatement output is held at its observed level in the early
    # stages only when it is observed; without abatement data the
    # programme is Eq. 6, whose abatement is endogenous throughout
    b1 <- b2 <- b3 <- b4 <- b5 <- rep(NA_real_, L)
    for (i in seq_len(L)) {
      ps <- all_peers
      b1[i] <- stage(i, ps, TRUE, TRUE, has_a, TRUE)
      b2[i] <- if (has_quality) stage(i, ps, TRUE, TRUE, has_a, FALSE)
               else b1[i]
      b3[i] <- stage(i, ps, FALSE, TRUE, has_a, FALSE)
      b4[i] <- if (has_a) stage(i, ps, FALSE, TRUE, FALSE, FALSE)
               else b3[i]
      b5[i] <- if (has_xa) stage(i, ps, FALSE, FALSE, FALSE, FALSE)
               else b4[i]
    }
    n_na <- sum(is.na(b1) | is.na(b5))
    if (n_na > 0) {
      warning(sprintf(paste0(
        "%d of %d DMUs have infeasible or failed stage LPs; the ",
        "affected components are NA."), n_na, L), call. = FALSE)
    }
    results <- data.frame(
      id = tech$id,
      b = b_p,
      b_star = b5,
      te_production = b1 / b_p,
      quality = b2 / b1,
      ae_production = b3 / b2,
      te_abatement = b4 / b3,
      ae_abatement = b5 / b4,
      total = b5 / b_p,
      stringsAsFactors = FALSE
    )
    results <- .insert_group(results, tech$group)
  }

  structure(
    list(results = results, type = type, returns = returns),
    class = "pgt_decomp"
  )
}

# Component column names per decomposition type.
.decomp_components <- function(type) {
  if (type == "envelope") c("WR", "TGR") else
    c("te_production", "quality", "ae_production", "te_abatement",
      "ae_abatement")
}

# Per-DMU attribution of the total reduction potential 1 - total to the
# sequential components: with components c1, c2, ..., cK (product =
# total), the reduction attributed to component k is
# c1 ... c_{k-1} (1 - c_k), and the attributions sum to 1 - total
# exactly.
.decomp_attribution <- function(results, comps) {
  total_red <- 1 - results$total
  cum <- rep(1, nrow(results))
  shares <- matrix(NA_real_, nrow(results), length(comps),
                   dimnames = list(NULL, paste0("share_", comps)))
  for (k in seq_along(comps)) {
    red_k <- cum * (1 - results[[comps[k]]])
    shares[, k] <- ifelse(total_red > 1e-9, red_k / total_red, NA_real_)
    cum <- cum * results[[comps[k]]]
  }
  shares
}

#' @export
print.pgt_decomp <- function(x, ...) {
  comps <- .decomp_components(x$type)
  cat(sprintf("pgt decomposition: type = %s, returns = %s\n",
              x$type, x$returns))
  cat(sprintf("  DMUs: %d   components: %s (product = total)\n",
              nrow(x$results), paste(comps, collapse = " x ")))
  meds <- vapply(x$results[c(comps, "total")],
                 stats::median, numeric(1), na.rm = TRUE)
  cat("  medians: ",
      paste(sprintf("%s = %.3f", names(meds), meds), collapse = ", "),
      "\n", sep = "")
  invisible(x)
}

#' @export
summary.pgt_decomp <- function(object, ...) {
  r <- object$results
  comps <- .decomp_components(object$type)
  shares <- .decomp_attribution(r, comps)
  r_sh <- cbind(r, shares)
  if (is.null(r_sh$group)) r_sh$group <- factor("all")

  by_group <- do.call(rbind, lapply(split(r_sh, r_sh$group), function(d) {
    out <- data.frame(group = d$group[1], n = nrow(d),
                      n_na = sum(!stats::complete.cases(
                        d[c(comps, "total")])))
    for (cc in c(comps, "total")) {
      out[[cc]] <- stats::median(d[[cc]], na.rm = TRUE)
    }
    for (sc in colnames(shares)) {
      out[[sc]] <- stats::median(d[[sc]], na.rm = TRUE)
    }
    out
  }))
  rownames(by_group) <- NULL

  structure(
    list(type = object$type, returns = object$returns,
         n = nrow(r), by_group = by_group, components = comps),
    class = "summary.pgt_decomp"
  )
}

#' @export
print.summary.pgt_decomp <- function(x, ...) {
  cat(sprintf("pgt decomposition: type = %s, returns = %s, DMUs = %d\n",
              x$type, x$returns, x$n))
  cat(sprintf("  components: %s (product = total)\n",
              paste(x$components, collapse = " x ")))
  cat("\nGroup medians (components, total, attribution shares):\n")
  print(x$by_group, row.names = FALSE, digits = 3)
  cat("\n(share_* are medians of DMU-level shares; they need not sum to 1.\n")
  cat(" n_na counts DMUs with any NA component.)\n")
  invisible(x)
}

#' @export
plot.pgt_decomp <- function(x, ...) {
  comps <- c(.decomp_components(x$type), "total")
  r <- x$results
  if (is.null(r$group)) r$group <- factor("all")
  med <- t(vapply(split(r, r$group), function(d) {
    vapply(comps, function(cc) stats::median(d[[cc]], na.rm = TRUE),
           numeric(1))
  }, numeric(length(comps))))
  graphics::barplot(t(med), beside = TRUE, ylim = c(0, 1),
                    legend.text = comps,
                    ylab = "Median efficiency component", ...)
  invisible(x)
}
