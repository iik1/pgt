#' Construct a pollution-generating technology
#'
#' Bundles the data of a pollution-generating technology, one or more
#' good outputs, one or more bad outputs, and material inputs, together
#' with the material flow coefficients that tie them to the
#' materials-balance identity. The returned object is the input to
#' [pgt()], [pgt_decompose()] and [mb_check()].
#'
#' The materials-balance principle requires that the pollutant content
#' of the inputs is either embodied in the good output, removed by
#' end-of-pipe abatement, or emitted:
#' \deqn{u'x_l - v'y_l \ge b_l,}
#' where \eqn{u} holds the material flow coefficients of the inputs (for
#' example tonnes of CO2 potential per unit of input) and \eqn{v} the
#' coefficients of the good outputs (pollutant retained in the
#' products). Data already expressed in pollutant-potential units use
#' the default \code{u = 1} for every input. When the abatement output
#' \code{a} is observed, the identity closes as an equality,
#' \eqn{u'x_l - v'y_l = b_l + a_l}. Use [mb_check()] to audit the
#' identity before estimation.
#'
#' Material flow coefficients may vary across DMUs (heterogeneous input
#' or output quality; Eder 2022, Rodseth 2025). Supply \code{u} as an
#' \eqn{L \times N} matrix and \code{v} as a length-\eqn{L} vector
#' (single output) or \eqn{L \times M} matrix to attach DMU-specific
#' coefficients. With homogeneous coefficients a length-\eqn{N} vector
#' \code{u} and scalar \code{v} suffice, and all estimators reduce to
#' the homogeneous-quality case.
#'
#' With several intended outputs, supply \code{y} as an
#' \eqn{L \times M} matrix (one named column per output). Every
#' estimator supports \eqn{M > 1} except \code{model = "fdmo"} in
#' [pgt()] and [pgt_ml()]; the retained content of pollutant \eqn{p}
#' generalises to \eqn{v_p'y_l} with one coefficient per output.
#'
#' With several pollutants, supply \code{b} as an \eqn{L \times P}
#' matrix (one named column per pollutant), \code{u} as a named list
#' with one element per pollutant (each a length-\eqn{N} vector or
#' \eqn{L \times N} matrix), and \code{v} as a named list with one
#' element per pollutant (each a scalar, length-\eqn{M} vector or
#' \eqn{L \times M} matrix), an \eqn{M \times P} matrix, or an
#' \eqn{L \times M \times P} array; the single-output shapes of
#' earlier versions (length-\eqn{P} vector, \eqn{L \times P} matrix)
#' remain accepted. Each pollutant carries its own materials-balance
#' identity.
#'
#' @param x Numeric matrix or data frame of material inputs, one row per
#'   decision-making unit (DMU), one column per input.
#' @param y Good outputs, strictly positive: a numeric vector for a
#'   single intended output, or an \eqn{L \times M} matrix or data
#'   frame with one named column per output.
#' @param b Bad outputs (controlled emissions), strictly positive: a
#'   numeric vector for a single pollutant, or an \eqn{L \times P}
#'   matrix or data frame with one column per pollutant.
#' @param u Material flow coefficients of the inputs. A length-\eqn{N}
#'   vector (one coefficient per input column, shared by all DMUs), an
#'   \eqn{L \times N} matrix (DMU-specific coefficients), or, with
#'   several pollutants, a named list with one such element per
#'   pollutant. Defaults to 1 for every input and pollutant.
#' @param v Material flow coefficients of the good outputs (pollutant
#'   embodied in the products). A scalar, a length-\eqn{M} vector (one
#'   coefficient per output), a length-\eqn{L} vector (single output,
#'   DMU-specific), or an \eqn{L \times M} matrix; with several
#'   pollutants, a named list with one such element per pollutant, an
#'   \eqn{M \times P} matrix, an \eqn{L \times M \times P} array, or
#'   the single-output shapes of earlier versions (length-\eqn{P}
#'   vector, \eqn{L \times P} matrix). Defaults to 0.
#' @param a Optional abatement output (pollutant removed by end-of-pipe
#'   control), non-negative: a numeric vector for a single pollutant or
#'   an \eqn{L \times P} matrix. Required by the models that treat
#'   pollution control explicitly; see [pgt()].
#' @param x_abate Optional marker of the pollution-control input columns
#'   of \code{x} (Rodseth 2025, Eq. 7): a character vector of column
#'   names or an integer vector of column positions. The remaining
#'   columns are production inputs. Used by the Rodseth decomposition
#'   of [pgt_decompose()], whose final stage frees the
#'   pollution-control inputs to separate abatement-input allocative
#'   efficiency; the [pgt()] estimators treat all inputs alike.
#' @param polluting Optional marker of the emission-generating input
#'   columns of \code{x} (the residual-generating sub-technology of the
#'   by-production model; Murty, Russell and Levkoff 2012): a character
#'   vector of column names or an integer vector of column positions.
#'   Defaults, when the by-production model is fitted, to the inputs
#'   with a positive material flow coefficient.
#' @param group Optional factor (or vector coercible to factor) assigning
#'   each DMU to a technology group, for example a production route.
#'   Required by [pgt_decompose()] and by group-referenced estimation.
#' @param period Optional factor (or vector coercible to factor) giving
#'   the time period of each row, for panel data. Required by [pgt_ml()].
#'   The pair \code{(id, period)} identifies an observation.
#' @param id Optional vector of DMU identifiers. Defaults to the row
#'   names of \code{x}, or a running index.
#'
#' @return An object of class \code{"pgt_tech"}: a list with elements
#'   \code{x} (\eqn{L \times N} matrix), \code{y} (\eqn{L \times M}
#'   matrix), \code{b} (\eqn{L \times P} matrix), \code{u}
#'   (\eqn{L \times N \times P} array), \code{v}
#'   (\eqn{L \times M \times P} array), \code{a} (\eqn{L \times P}
#'   matrix or \code{NULL}), \code{x_abate} (integer vector of
#'   pollution-control input columns, possibly empty), \code{group},
#'   \code{id}, \code{L}, \code{N}, \code{M}, \code{P},
#'   \code{outputs} and \code{pollutants} (character vectors of output
#'   and pollutant names). Single-output, single-pollutant input is
#'   stored in the same canonical shapes with \eqn{M = P = 1}.
#'
#' @references
#' Rodseth, K. L. (2025). On the development of a unified, nonparametric
#' materials balance-based efficiency analysis model and its
#' applications. \emph{Journal of Productivity Analysis}, 64(3),
#' 305--319. \doi{10.1007/s11123-025-00768-0}
#'
#' Coelli, T., Lauwers, L., & Van Huylenbroeck, G. (2007). Environmental
#' efficiency measurement and the materials balance condition.
#' \emph{Journal of Productivity Analysis}, 28(1--2), 3--12.
#' \doi{10.1007/s11123-007-0052-8}
#'
#' Eder, A. (2022). Environmental efficiency measurement when producers
#' control pollutants under heterogeneous conditions: a generalization
#' of the materials balance approach. \emph{Journal of Productivity
#' Analysis}, 57(2), 157--176.
#' \doi{10.1007/s11123-021-00623-y}
#'
#' @seealso [mb_check()], [pgt()], [pgt_decompose()]
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
#' tech
#'
#' # Explicit abatement: Rodseth (2025) pig-finishing example, with the
#' # uncontrolled-emission aggregate as the material carrier
#' data(pigfarms)
#' tech2 <- pgt_tech(
#'   x = pigfarms[, c("uncontrolled", "labor", "capital")],
#'   y = pigfarms$meat,
#'   b = pigfarms$controlled,
#'   u = c(1, 0, 0),
#'   a = pigfarms$abatement,
#'   id = pigfarms$farm
#' )
#' @export
pgt_tech <- function(x, y, b, u = NULL, v = 0, a = NULL, x_abate = NULL,
                     polluting = NULL, group = NULL, period = NULL,
                     id = NULL) {
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop("'x' must be numeric.", call. = FALSE)
  }
  storage.mode(x) <- "double"
  y <- .canon_good(y)
  L <- nrow(y)
  M <- ncol(y)
  outputs <- colnames(y)
  N <- ncol(x)

  if (nrow(x) != L) {
    stop("'x' and 'y' must describe the same number of DMUs.",
         call. = FALSE)
  }
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("input", seq_len(N))
  }

  b <- .canon_bad(b, L, "b")
  P <- ncol(b)
  pollutants <- colnames(b)

  if (anyNA(x) || anyNA(y) || anyNA(b)) {
    stop("missing values in 'x', 'y' or 'b'; impute or drop incomplete ",
         "DMUs before calling pgt_tech().", call. = FALSE)
  }
  if (any(is.infinite(x)) || any(is.infinite(y)) || any(is.infinite(b))) {
    stop("'x', 'y' and 'b' must be finite.", call. = FALSE)
  }
  if (any(x < 0)) {
    stop("'x' must be non-negative.", call. = FALSE)
  }
  if (any(y <= 0)) {
    stop("'y' must be strictly positive.", call. = FALSE)
  }
  if (any(b <= 0)) {
    stop("'b' must be strictly positive (efficiency scores are ratios ",
         "b*/b).", call. = FALSE)
  }

  u <- .canon_u(u, L, N, P, pollutants, colnames(x))
  v <- .canon_v(v, L, M, P, outputs, pollutants)

  if (!is.null(a)) {
    a <- .canon_bad(a, L, "a", pollutants = pollutants)
    if (ncol(a) != P) {
      stop("'a' must have one column per pollutant in 'b'.",
           call. = FALSE)
    }
    if (anyNA(a) || any(!is.finite(a)) || any(a < 0)) {
      stop("'a' must be non-negative and finite.", call. = FALSE)
    }
  }

  x_abate <- .canon_x_abate(x_abate, colnames(x), "x_abate",
                            allow_all = FALSE)
  polluting <- .canon_x_abate(polluting, colnames(x), "polluting",
                              allow_all = TRUE)

  if (!is.null(group)) {
    group <- as.factor(group)
    if (length(group) != L || anyNA(group)) {
      stop("'group' must have one non-missing value per DMU.",
           call. = FALSE)
    }
    group <- droplevels(group)
  }

  if (!is.null(period)) {
    period <- as.factor(period)
    if (length(period) != L || anyNA(period)) {
      stop("'period' must have one non-missing value per DMU.",
           call. = FALSE)
    }
    period <- droplevels(period)
  }

  if (is.null(id)) {
    id <- if (!is.null(rownames(x))) rownames(x) else as.character(seq_len(L))
  }
  id <- as.character(id)
  if (length(id) != L) {
    stop("'id' must have one value per DMU.", call. = FALSE)
  }
  if (!is.null(period)) {
    key <- paste(id, as.character(period), sep = "\r")
    if (anyDuplicated(key)) {
      dup <- unique(id[duplicated(key)])
      stop("the pair (id, period) must identify each observation ",
           "uniquely; duplicated id(s): ",
           paste(utils::head(dup, 5L), collapse = ", "), ".",
           call. = FALSE)
    }
  }

  structure(
    list(x = x, y = y, b = b, u = u, v = v, a = a, x_abate = x_abate,
         polluting = polluting, group = group, period = period,
         id = id, L = L, N = N, M = M, P = P, outputs = outputs,
         pollutants = pollutants),
    class = "pgt_tech"
  )
}

# Canonicalise a bad-output (or abatement) specification to an L x P
# matrix with pollutant column names.
.canon_bad <- function(b, L, what, pollutants = NULL) {
  if (is.data.frame(b)) b <- as.matrix(b)
  if (is.matrix(b)) {
    if (!is.numeric(b)) {
      stop("'", what, "' must be numeric.", call. = FALSE)
    }
    if (nrow(b) != L) {
      stop("'", what, "' must have one row per DMU.", call. = FALSE)
    }
    storage.mode(b) <- "double"
    cn <- colnames(b)
    named <- !is.null(cn) && !any(is.na(cn) | cn == "")
    if (!is.null(pollutants)) {
      # 'a' aligned to the pollutant set already fixed by 'b': match by
      # name when the columns are named (as u and v do), else by
      # position.
      if (named) {
        if (!setequal(cn, pollutants)) {
          stop("'", what, "' column names must match the pollutant ",
               "names: ", paste(pollutants, collapse = ", "), ".",
               call. = FALSE)
        }
        b <- b[, pollutants, drop = FALSE]
      } else if (ncol(b) == length(pollutants)) {
        colnames(b) <- pollutants
      }
    } else if (!named) {
      colnames(b) <- if (ncol(b) == 1L) what else
        paste0(what, seq_len(ncol(b)))
    }
  } else {
    if (is.factor(b) || !is.numeric(b)) {
      stop("'", what, "' must be numeric (not a factor).", call. = FALSE)
    }
    if (length(b) != L) {
      stop("'", what, "' must have one value per DMU.", call. = FALSE)
    }
    b <- matrix(as.numeric(b), ncol = 1L,
                dimnames = list(NULL, if (is.null(pollutants)) what else
                  pollutants[1L]))
  }
  rownames(b) <- NULL
  b
}

# Canonicalise one pollutant's input coefficients to an L x N matrix.
.canon_u_one <- function(spec, L, N, input_names, label) {
  if (is.null(spec)) {
    spec <- matrix(1, L, N)
  } else if (is.matrix(spec)) {
    if (!is.numeric(spec) || nrow(spec) != L || ncol(spec) != N) {
      stop("matrix 'u'", label, " must be L x N (one row per DMU, one ",
           "column per input).", call. = FALSE)
    }
    storage.mode(spec) <- "double"
  } else {
    if (is.factor(spec) || !is.numeric(spec) || length(spec) != N) {
      stop("'u'", label, " must be a non-negative numeric vector with ",
           "one coefficient per input column, or an L x N matrix.",
           call. = FALSE)
    }
    spec <- matrix(as.numeric(spec), L, N, byrow = TRUE)
  }
  if (any(!is.finite(spec)) || any(spec < 0)) {
    stop("'u'", label, " must be non-negative and finite.", call. = FALSE)
  }
  dimnames(spec) <- list(NULL, input_names)
  spec
}

# Canonicalise the full input-coefficient specification to an
# L x N x P array.
.canon_u <- function(u, L, N, P, pollutants, input_names) {
  if (is.array(u) && length(dim(u)) == 3L) {
    if (!all(dim(u) == c(L, N, P))) {
      stop("array 'u' must be L x N x P.", call. = FALSE)
    }
    u <- lapply(seq_len(P), function(p) u[, , p, drop = TRUE])
    u <- lapply(u, function(m) matrix(m, L, N))
  } else if (is.list(u)) {
    if (length(u) != P) {
      stop("list 'u' must have one element per pollutant (", P, ").",
           call. = FALSE)
    }
    if (!is.null(names(u))) {
      if (!setequal(names(u), pollutants)) {
        stop("names of list 'u' must match the pollutant names: ",
             paste(pollutants, collapse = ", "), ".", call. = FALSE)
      }
      u <- u[pollutants]
    }
  } else if (P > 1L) {
    if (!is.null(u)) {
      stop("with several pollutants, supply 'u' as a named list with ",
           "one element per column of 'b' (or leave NULL for 1s).",
           call. = FALSE)
    }
    u <- vector("list", P)
  } else {
    u <- list(u)
  }
  if (!is.list(u)) u <- list(u)
  out <- array(NA_real_, dim = c(L, N, P),
               dimnames = list(NULL, input_names, pollutants))
  for (p in seq_len(P)) {
    label <- if (P > 1L) paste0("[['", pollutants[p], "']]") else ""
    out[, , p] <- .canon_u_one(u[[p]], L, N, input_names, label)
  }
  out
}

# Canonicalise the good-output specification to an L x M matrix with
# output column names.
.canon_good <- function(y) {
  if (is.data.frame(y)) y <- as.matrix(y)
  if (is.matrix(y)) {
    if (!is.numeric(y)) {
      stop("'y' must be numeric.", call. = FALSE)
    }
    storage.mode(y) <- "double"
    cn <- colnames(y)
    if (is.null(cn) || any(is.na(cn) | cn == "")) {
      colnames(y) <- if (ncol(y) == 1L) "y" else
        paste0("y", seq_len(ncol(y)))
    }
  } else {
    if (is.factor(y) || !is.numeric(y)) {
      stop("'y' must be numeric (not a factor).", call. = FALSE)
    }
    y <- matrix(as.numeric(y), ncol = 1L, dimnames = list(NULL, "y"))
  }
  rownames(y) <- NULL
  y
}

# Canonicalise one pollutant's good-output coefficients to an
# L x M matrix.
.canon_v_one <- function(spec, L, M, output_names, label) {
  if (is.null(spec)) {
    spec <- matrix(0, L, M)
  } else if (is.matrix(spec)) {
    if (!is.numeric(spec) || nrow(spec) != L || ncol(spec) != M) {
      stop("matrix 'v'", label, " must be L x M (one row per DMU, one ",
           "column per good output).", call. = FALSE)
    }
    storage.mode(spec) <- "double"
  } else {
    if (is.factor(spec) || !is.numeric(spec)) {
      stop("'v'", label, " must be numeric.", call. = FALSE)
    }
    spec <- as.numeric(spec)
    if (length(spec) == 1L) {
      spec <- matrix(spec, L, M)
    } else if (M == 1L && length(spec) == L) {
      spec <- matrix(spec, L, 1L)
    } else if (M > 1L && length(spec) == M) {
      if (L == M) {
        stop("ambiguous 'v'", label, ": its length equals both the ",
             "number of DMUs and the number of good outputs; supply an ",
             "L x M matrix.", call. = FALSE)
      }
      spec <- matrix(spec, L, M, byrow = TRUE)
    } else {
      stop("'v'", label, " must be a scalar, a length-M vector (one ",
           "coefficient per good output), a length-L vector (single ",
           "output), or an L x M matrix.", call. = FALSE)
    }
  }
  spec
}

# Canonicalise the full good-output coefficient specification to an
# L x M x P array (DMU x output x pollutant).
.canon_v <- function(v, L, M, P, outputs, pollutants) {
  if (is.array(v) && length(dim(v)) == 3L) {
    if (!all(dim(v) == c(L, M, P))) {
      stop("array 'v' must be L x M x P.", call. = FALSE)
    }
    v <- lapply(seq_len(P), function(p) matrix(v[, , p], L, M))
  } else if (is.list(v)) {
    if (length(v) != P) {
      stop("list 'v' must have one element per pollutant (", P, ").",
           call. = FALSE)
    }
    if (!is.null(names(v))) {
      if (!setequal(names(v), pollutants)) {
        stop("names of list 'v' must match the pollutant names: ",
             paste(pollutants, collapse = ", "), ".", call. = FALSE)
      }
      v <- v[pollutants]
    }
  } else if (is.matrix(v) && M == 1L && P > 1L && nrow(v) == L &&
             ncol(v) == P) {
    # historical L x P shape for a single good output
    v <- lapply(seq_len(P), function(p) matrix(v[, p], L, 1L))
  } else if (!is.matrix(v) && is.numeric(v) && M == 1L && P > 1L &&
             length(v) == P && P != L) {
    # historical length-P vector for a single good output
    v <- lapply(as.numeric(v), function(el) matrix(el, L, 1L))
  } else if (!is.matrix(v) && is.numeric(v) && M == 1L && P > 1L &&
             length(v) == P && P == L) {
    stop("ambiguous 'v': its length equals both the number of DMUs ",
         "and the number of pollutants; supply an L x P matrix or ",
         "a named list.", call. = FALSE)
  } else if (P > 1L) {
    if (is.matrix(v) && nrow(v) == M && ncol(v) == P) {
      # per-output-per-pollutant coefficients shared by all DMUs
      v <- lapply(seq_len(P), function(p)
        matrix(v[, p], L, M, byrow = TRUE))
    } else if (!is.matrix(v) && is.numeric(v) && length(v) == 1L) {
      v <- rep(list(matrix(as.numeric(v), L, M)), P)
    } else {
      stop("with several pollutants, supply 'v' as a named list with ",
           "one element per column of 'b' (each a scalar, length-M ",
           "vector, or L x M matrix), an M x P matrix, or an ",
           "L x M x P array.", call. = FALSE)
    }
  } else {
    v <- list(v)
  }
  out <- array(NA_real_, dim = c(L, M, P),
               dimnames = list(NULL, outputs, pollutants))
  for (p in seq_len(P)) {
    label <- if (P > 1L) paste0("[['", pollutants[p], "']]") else ""
    out[, , p] <- .canon_v_one(v[[p]], L, M, outputs, label)
  }
  if (any(!is.finite(out)) || any(out < 0)) {
    stop("'v' must be non-negative and finite.", call. = FALSE)
  }
  out
}

# Pollutant p's retained content v_p'y per DMU (one value per DMU).
.retained <- function(tech, p = 1L) {
  rowSums(matrix(tech$v[, , p], tech$L, tech$M) * tech$y)
}

# Canonicalise an input-column marker (pollution-control or
# emission-generating) to sorted column positions.
.canon_x_abate <- function(marker, input_names, what = "x_abate",
                           allow_all = FALSE) {
  if (is.null(marker)) {
    return(integer(0))
  }
  if (is.character(marker)) {
    idx <- match(marker, input_names)
    if (anyNA(idx)) {
      stop("'", what, "' names not found in the columns of 'x': ",
           paste(marker[is.na(idx)], collapse = ", "), ".",
           call. = FALSE)
    }
  } else if (is.numeric(marker)) {
    idx <- as.integer(marker)
    if (any(idx < 1L) || any(idx > length(input_names)) ||
        anyNA(idx)) {
      stop("'", what, "' positions must lie in 1..N.", call. = FALSE)
    }
  } else {
    stop("'", what, "' must be a character or integer vector.",
         call. = FALSE)
  }
  idx <- sort(unique(idx))
  if (!allow_all && length(idx) == length(input_names)) {
    stop("'", what, "' cannot mark every input.", call. = FALSE)
  }
  idx
}

# TRUE when the coefficients of pollutant p are shared by all DMUs.
.u_is_uniform <- function(tech, p = 1L) {
  um <- tech$u[, , p, drop = FALSE]
  all(apply(matrix(um, tech$L, tech$N), 2L,
            function(col) diff(range(col)) == 0))
}

#' @export
print.pgt_tech <- function(x, ...) {
  cat("Pollution-generating technology\n")
  cat(sprintf("  DMUs: %d   inputs: %d   good outputs: %d   bad outputs: %d\n",
              x$L, x$N, x$M, x$P))
  cat(sprintf("  inputs: %s\n", paste(colnames(x$x), collapse = ", ")))
  if (length(x$x_abate)) {
    cat(sprintf("  pollution-control inputs: %s\n",
                paste(colnames(x$x)[x$x_abate], collapse = ", ")))
  }
  for (p in seq_len(x$P)) {
    tag <- if (x$P > 1L) paste0(" [", x$pollutants[p], "]") else ""
    vm <- matrix(x$v[, , p], x$L, x$M)
    v_uniform <- all(apply(vm, 2L, function(col) diff(range(col)) == 0))
    if (.u_is_uniform(x, p) && v_uniform) {
      cat(sprintf("  material flow coefficients%s: u = [%s], v = [%s]\n",
                  tag,
                  paste(format(x$u[1L, , p], digits = 4), collapse = ", "),
                  paste(format(vm[1L, ], digits = 4), collapse = ", ")))
    } else {
      cat(sprintf(
        "  material flow coefficients%s: DMU-specific (u in [%s], v in [%s])\n",
        tag,
        paste(format(range(x$u[, , p]), digits = 4), collapse = ", "),
        paste(format(range(vm), digits = 4), collapse = ", ")))
    }
  }
  if (!is.null(x$a)) {
    cat("  abatement output: observed\n")
  }
  if (!is.null(x$group)) {
    tab <- table(x$group)
    cat(sprintf("  groups: %s\n",
                paste(sprintf("%s (%d)", names(tab), tab), collapse = ", ")))
  }
  mb <- mb_check(x)
  nv <- attr(mb, "n_violations")
  if (nv > 0) {
    cat(sprintf("  materials balance: %d of %d DMU-pollutant accounts violate the identity (see mb_check())\n",
                nv, x$L * x$P))
  } else {
    cat("  materials balance: satisfied for all DMUs\n")
  }
  invisible(x)
}
