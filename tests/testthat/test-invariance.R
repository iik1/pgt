# Generic estimator properties: DMU-order invariance, duplicate
# neutrality, and the retained-content dichotomy example.

make_perm_data <- function() {
  set.seed(7)
  L <- 20
  x <- matrix(runif(2 * L, 10, 100), L, 2)
  u <- c(1, 1)
  pot <- as.vector(x %*% u)
  v <- 0.02
  y <- runif(L, 1, 40)  # 0.02 * 40 << min(pot), so every account holds
  b <- (pot - v * y) * runif(L, 0.6, 0.98)
  list(x = x, y = y, b = b, u = u, v = v,
       id = sprintf("D%02d", seq_len(L)))
}

test_that("scores are invariant to DMU ordering", {
  d <- make_perm_data()
  perm <- sample(length(d$id))
  for (m in c("wgd", "wgd_input_fixed", "envelope")) {
    f1 <- pgt(pgt_tech(d$x, d$y, d$b, u = d$u, v = d$v, id = d$id),
              model = m)$results
    f2 <- pgt(pgt_tech(d$x[perm, , drop = FALSE], d$y[perm], d$b[perm],
                       u = d$u, v = d$v, id = d$id[perm]),
              model = m)$results
    expect_equal(f2$efficiency[match(d$id, f2$id)], f1$efficiency,
                 tolerance = 1e-9)
    expect_equal(f2$status[match(d$id, f2$id)], f1$status)
  }
})

test_that("duplicating a DMU changes no score", {
  # A duplicate leaves the convex hull unchanged, so every original
  # keeps its score and the copy scores as its source row.
  d <- make_perm_data()
  k <- 3L
  for (m in c("wgd", "wgd_input_fixed", "envelope")) {
    f1 <- pgt(pgt_tech(d$x, d$y, d$b, u = d$u, v = d$v, id = d$id),
              model = m)$results
    f2 <- pgt(pgt_tech(rbind(d$x, d$x[k, ]), c(d$y, d$y[k]),
                       c(d$b, d$b[k]), u = d$u, v = d$v,
                       id = c(d$id, "DUP")),
              model = m)$results
    expect_equal(f2$efficiency[match(d$id, f2$id)], f1$efficiency,
                 tolerance = 1e-9)
    expect_equal(f2$efficiency[f2$id == "DUP"], f1$efficiency[k],
                 tolerance = 1e-9)
  }
})

test_that("the retained-content term moves wgd exactly when the emission boundary falls in y", {
  # Two units, (y, b) = (1, 2) and (2, 1): expanding output buys lower
  # peer emission, so the equality boundary falls at rate 1. A's
  # envelope score is 1/2 through B; the wgd score rises in v and
  # switches to self-reference once v exceeds the descent rate.
  # x = 5 with u = 1 keeps both accounts satisfied at every v used
  # (potential 5 >= v*y + b); x and u never enter the wgd/envelope LPs.
  ab <- function(v) pgt_tech(x = matrix(5, 2, 1), y = c(1, 2),
                             b = c(2, 1), u = 1, v = v,
                             id = c("A", "B"))
  effA <- vapply(c(0, 0.5, 2), function(v)
    pgt(ab(v), model = "wgd")$results$efficiency[1], numeric(1))
  expect_equal(effA, c(0.5, 0.75, 1), tolerance = 1e-9)
  expect_equal(pgt(ab(0), model = "envelope")$results$efficiency[1],
               0.5, tolerance = 1e-9)
  expect_equal(names(pgt(ab(0.5), model = "wgd")$weights[[1]]), "B")
  expect_equal(names(pgt(ab(2), model = "wgd")$weights[[1]]), "A")
})
