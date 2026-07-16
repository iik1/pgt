#' Compare competing axiom systems on identical data
#'
#' Fits several pollution-generating technology models to the same
#' [pgt_tech()] object and reports how their environmental efficiency
#' scores and rankings agree. The materials-balance debate turns on which
#' disposability axiom is appropriate for bad outputs (Forsund 2021 and
#' the accompanying Journal of Productivity Analysis symposium); this
#' function puts the competing systems side by side so the sensitivity of
#' conclusions to that choice is visible.
#'
#' All models report an environmental efficiency \eqn{b^*/b \in (0, 1]}
#' (the directional model \code{"fdmo"} is excluded, since its score is a
#' gross inefficiency on a different scale). Rank agreement is measured
#' by Spearman correlation over the DMUs solved by both members of each
#' pair.
#'
#' @param tech A [pgt_tech()] object.
#' @param models Character vector of models to compare. Defaults to
#'   \code{c("wgd", "byprod", "mb_cost", "wd")}: the weak-G-disposability
#'   (materials balance), by-production, materials-balance cost and
#'   weak-disposability systems.
#' @param returns Returns to scale passed to every model: \code{"vrs"}
#'   (default) or \code{"crs"}.
#' @param peers Reference set passed to every model: \code{"all"}
#'   (default) or \code{"group"}.
#' @param pollutant For a multi-pollutant technology, the pollutant to
#'   compare on. Defaults to the first.
#'
#' @return An object of class \code{"pgt_compare"}: a list with
#'   \describe{
#'     \item{\code{scores}}{Data frame of environmental efficiency, one
#'       column per model, one row per DMU (with \code{id} and
#'       \code{group}).}
#'     \item{\code{spearman}}{Matrix of pairwise Spearman rank
#'       correlations.}
#'     \item{\code{agreement}}{Data frame with, per model, the number of
#'       DMUs solved, the median score, and the share of DMUs it places
#'       in the bottom efficiency quartile that the weak-G-disposability
#'       model also places there (top-of-agenda overlap).}
#'     \item{\code{models}, \code{returns}, \code{peers}}{Call settings.}
#'   }
#'
#' @references
#' Forsund, F. R. (2021). Performance measurement and joint production of
#' intended and unintended outputs. \emph{Journal of Productivity
#' Analysis}, 55(3), 157--175.
#'
#' Murty, S., Russell, R. R., & Levkoff, S. B. (2012). On modeling
#' pollution-generating technologies. \emph{Journal of Environmental
#' Economics and Management}, 64(1), 117--135.
#'
#' @seealso [pgt()]
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
#' cmp <- compare_models(tech, models = c("wgd", "byprod", "mb_cost"))
#' cmp
#' @export
compare_models <- function(tech,
                           models = c("wgd", "byprod", "mb_cost", "wd"),
                           returns = c("vrs", "crs"),
                           peers = c("all", "group"), pollutant = 1L) {
  stopifnot(inherits(tech, "pgt_tech"))
  returns <- match.arg(returns)
  peers <- match.arg(peers)
  valid <- c("wgd", "envelope", "byprod", "mb_cost", "wd")
  bad <- setdiff(models, valid)
  if (length(bad)) {
    stop("compare_models() compares efficiency-scored models; drop: ",
         paste(bad, collapse = ", "),
         ". Available: ", paste(valid, collapse = ", "), ".",
         call. = FALSE)
  }
  if (length(models) < 2L) {
    stop("supply at least two models to compare.", call. = FALSE)
  }

  fits <- lapply(models, function(m) {
    suppressWarnings(pgt(tech, model = m, returns = returns,
                         peers = peers, pollutant = pollutant))
  })
  names(fits) <- models

  scores <- data.frame(id = tech$id, stringsAsFactors = FALSE)
  if (!is.null(tech$group)) scores$group <- tech$group
  for (m in models) scores[[m]] <- fits[[m]]$results$efficiency

  sp <- matrix(NA_real_, length(models), length(models),
               dimnames = list(models, models))
  for (a in models) for (bm in models) {
    ok <- stats::complete.cases(scores[[a]], scores[[bm]])
    sp[a, bm] <- if (sum(ok) >= 3L)
      stats::cor(scores[[a]][ok], scores[[bm]][ok], method = "spearman")
    else NA_real_
  }

  ref <- models[1]
  ref_bottom <- .bottom_quartile(scores[[ref]])
  agreement <- do.call(rbind, lapply(models, function(m) {
    s <- scores[[m]]
    mb <- .bottom_quartile(s)
    both <- sum(mb & ref_bottom, na.rm = TRUE)
    # share of this model's own worst-quartile DMUs that the reference
    # model also places in its worst quartile
    denom <- sum(mb, na.rm = TRUE)
    data.frame(
      model = m, n_solved = sum(!is.na(s)),
      median = stats::median(s, na.rm = TRUE),
      bottom_q_overlap = if (denom > 0) both / denom else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  rownames(agreement) <- NULL

  structure(
    list(scores = scores, spearman = sp, agreement = agreement,
         models = models, reference = ref, returns = returns,
         peers = peers),
    class = "pgt_compare"
  )
}

# Logical flag for the worst (bottom) efficiency quartile.
.bottom_quartile <- function(s) {
  q <- stats::quantile(s, 0.25, na.rm = TRUE)
  s <= q
}

#' @export
print.pgt_compare <- function(x, ...) {
  cat(sprintf("pgt model comparison: %d models, returns = %s, peers = %s\n",
              length(x$models), x$returns, x$peers))
  cat(sprintf("  models: %s\n", paste(x$models, collapse = ", ")))
  cat("\nEnvironmental efficiency (b*/b) summary:\n")
  print(x$agreement, row.names = FALSE, digits = 4)
  cat("\nSpearman rank correlation:\n")
  print(round(x$spearman, 3))
  worst <- which(x$spearman == min(x$spearman, na.rm = TRUE), arr.ind = TRUE)
  if (nrow(worst)) {
    w <- worst[1, ]
    cat(sprintf("\n  largest ranking disagreement: %s vs %s (rho = %.3f)\n",
                x$models[w[1]], x$models[w[2]],
                x$spearman[w[1], w[2]]))
  }
  invisible(x)
}

#' @export
as.data.frame.pgt_compare <- function(x, ...) {
  x$scores
}

#' @export
plot.pgt_compare <- function(x, ...) {
  score_cols <- x$models
  m <- as.matrix(x$scores[score_cols])
  graphics::matplot(
    seq_len(nrow(m)), m[order(m[, 1]), , drop = FALSE],
    type = "l", lty = 1, col = seq_along(score_cols),
    xlab = "DMU (ordered by first model's efficiency)",
    ylab = "Environmental efficiency (b*/b)", ...
  )
  graphics::legend("bottomright", legend = score_cols,
                   col = seq_along(score_cols), lty = 1, bty = "n")
  invisible(x)
}
