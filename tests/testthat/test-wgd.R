test_that("wgd model solves a hand-computable example", {
  # Two DMUs, identical inputs and output, different emissions.
  # The dirty DMU is benchmarked against the clean one: b* = 5.
  x <- matrix(c(10, 10), 2, 1)
  tech <- pgt_tech(x, y = c(5, 5), b = c(10, 5))
  fit <- pgt(tech, model = "wgd")

  expect_equal(fit$results$status, c(0L, 0L))
  expect_equal(fit$results$b_star, c(5, 5), tolerance = 1e-8)
  expect_equal(fit$results$efficiency, c(0.5, 1), tolerance = 1e-8)
  # peer weights: DMU1 projects onto DMU2
  expect_equal(fit$weights[[1]], c(`2` = 1), tolerance = 1e-8)
})

test_that("materials-balance violations surface as infeasible LPs", {
  # DMU1 violates closure (b = 2 > u'x = 1) and, having the highest
  # output, can only reference itself: its LP is infeasible.
  x <- matrix(c(1, 1), 2, 1)
  tech <- pgt_tech(x, y = c(10, 1), b = c(2, 0.5))
  expect_warning(fit <- pgt(tech, model = "wgd"), "infeasible")
  expect_true(is.na(fit$results$b_star[1]))
  expect_true(is.na(fit$results$efficiency[1]))
  expect_equal(fit$results$efficiency[2], 1, tolerance = 1e-8)
})

test_that("envelope model matches the convex lower (y,b) envelope", {
  # y = (1, 3, 2), b = (1, 3, 4): DMU3 projects onto the midpoint of
  # DMUs 1 and 2 (y = 2, b = 2), so b* = 2 and efficiency = 0.5.
  x <- matrix(10, 3, 1)
  tech <- pgt_tech(x, y = c(1, 3, 2), b = c(1, 3, 4))
  fit <- pgt(tech, model = "envelope")

  expect_equal(fit$results$b_star, c(1, 3, 2), tolerance = 1e-8)
  expect_equal(fit$results$efficiency, c(1, 1, 0.5), tolerance = 1e-8)
  # dual of the output constraint = slope of the (y,b) envelope = 1
  expect_equal(fit$results$dual_output[3], 1, tolerance = 1e-6)
})

test_that("vrs and crs differ as expected for the envelope", {
  # y = (1, 2), b = (2, 2). VRS: DMU1 cannot scale peers down, b* = 2.
  # CRS: DMU2 scaled by 0.5 gives (y, b) = (1, 1), so b* = 1.
  x <- matrix(10, 2, 1)
  tech <- pgt_tech(x, y = c(1, 2), b = c(2, 2))
  vrs <- pgt(tech, model = "envelope", returns = "vrs")
  crs <- pgt(tech, model = "envelope", returns = "crs")

  expect_equal(vrs$results$b_star[1], 2, tolerance = 1e-8)
  expect_equal(crs$results$b_star[1], 1, tolerance = 1e-8)
})

test_that("crs envelope equals y_i * min(b/y)", {
  tech <- make_random_tech(L = 25, seed = 42)
  crs <- pgt(tech, model = "envelope", returns = "crs")
  expected <- tech$y * min(tech$b / tech$y)
  expect_equal(crs$results$b_star, expected, tolerance = 1e-6)
})

test_that("self-reference bounds envelope scores in (0, 1]", {
  tech <- make_random_tech(L = 50, N = 4, seed = 7)
  fit <- pgt(tech, model = "envelope")
  expect_true(all(fit$results$status == 0))
  expect_true(all(fit$results$efficiency > 0))
  expect_true(all(fit$results$efficiency <= 1 + 1e-8))
})

test_that("group peers never beat pooled peers", {
  tech <- make_random_tech(L = 40, seed = 11)
  grp <- pgt(tech, model = "envelope", peers = "group")
  all <- pgt(tech, model = "envelope", peers = "all")
  expect_true(all(all$results$b_star <= grp$results$b_star + 1e-8))
})

test_that("group peers require a group", {
  x <- matrix(1:4, 2, 2)
  tech <- pgt_tech(x, y = c(1, 2), b = c(1, 1))
  expect_error(pgt(tech, peers = "group"), "requires a 'group'")
})

test_that("wgd_rodseth is an alias for wgd", {
  tech <- make_random_tech(seed = 3)
  a <- pgt(tech, model = "wgd")
  b <- pgt(tech, model = "wgd_rodseth")
  expect_equal(a$results, b$results)
  expect_equal(a$model, "wgd")
})

test_that("envelope efficiency is invariant to scaling b", {
  tech <- make_random_tech(L = 30, seed = 5)
  fit1 <- pgt(tech, model = "envelope")
  tech2 <- pgt_tech(tech$x, tech$y, 3.7 * tech$b, u = tech$u,
                    v = tech$v, group = tech$group)
  fit2 <- pgt(tech2, model = "envelope")
  expect_equal(fit1$results$efficiency, fit2$results$efficiency,
               tolerance = 1e-8)
})

test_that("crs relaxes vrs for the wgd model", {
  tech <- make_random_tech(L = 25, N = 3, seed = 55)
  vrs <- pgt(tech, model = "wgd", returns = "vrs")
  crs <- pgt(tech, model = "wgd", returns = "crs")
  ok <- vrs$results$status == 0 & crs$results$status == 0
  expect_true(any(ok))
  expect_true(all(crs$results$b_star[ok] <= vrs$results$b_star[ok] + 1e-8))
})

test_that("wgd diagnostics carry the documented signs", {
  tech <- make_random_tech(L = 30, N = 3, seed = 77)
  fit <- pgt(tech, model = "wgd")
  ok <- fit$results$status == 0
  # output row is a >= row whose relaxation cannot lower b*: dual >= 0
  expect_true(all(fit$results$dual_output[ok] >= -1e-10))
  # the projection can never exceed the DMU's materials-balance ceiling
  expect_true(all(fit$results$mb_headroom[ok] >= -1e-8))
})

test_that("the v-term of the materials-balance cap binds feasibility", {
  # u'x = 10, v*y = 5, so the cap is 5. With b = (6, 5.5) no peer mix
  # fits under any DMU's cap: every LP must be infeasible. Any loosening
  # of the v-term (cap 10) would make both solvable.
  x <- matrix(c(10, 10), 2, 1)
  tech <- pgt_tech(x, y = c(10, 10), b = c(6, 5.5), v = 0.5)
  expect_warning(fit <- pgt(tech, model = "wgd"), "infeasible")
  expect_true(all(fit$results$status != 0))
  expect_true(all(is.na(fit$results$b_star)))

  # Companion: DMU1 violates its cap (b = 6 > 5) but solves through the
  # peer mix lambda = DMU2 (b = 4 <= 5): a feasible violator.
  tech2 <- pgt_tech(x, y = c(10, 10), b = c(6, 4), v = 0.5)
  fit2 <- pgt(tech2, model = "wgd")
  expect_equal(fit2$results$status, c(0L, 0L))
  expect_equal(fit2$results$b_star, c(4, 4), tolerance = 1e-8)
  expect_equal(fit2$results$efficiency[1], 4 / 6, tolerance = 1e-8)
  expect_equal(fit2$results$mb_headroom, c(1, 1), tolerance = 1e-8)
})

test_that("infeasible wgd LPs are confined to DMUs with exact gap < 0", {
  x <- matrix(c(10, 10), 2, 1)
  tech <- pgt_tech(x, y = c(10, 10), b = c(6, 5.5), v = 0.5)
  expect_warning(fit <- pgt(tech, model = "wgd"), "infeasible")
  mb <- mb_check(tech)
  infeasible <- fit$results$status != 0
  expect_true(all(mb$gap[infeasible] < 0))
  expect_equal(attr(mb, "n_exact"), 2L)
})

test_that("wgd with group peers restricts the reference set", {
  x <- matrix(10, 3, 1)
  tech <- pgt_tech(x, y = c(5, 5, 5), b = c(10, 5, 1),
                   group = c("A", "A", "B"), id = c("d1", "d2", "d3"))
  fit <- pgt(tech, model = "wgd", peers = "group")
  # d1 must not benchmark against group B's b = 1
  expect_equal(fit$results$b_star, c(5, 5, 1), tolerance = 1e-8)
  expect_equal(names(fit$weights[["d1"]]), "d2")
  expect_equal(names(fit$weights[["d3"]]), "d3")
})

test_that("wgd group peers equal the rodseth stage-1 minima", {
  tech <- make_random_tech(L = 25, N = 3, seed = 88)
  fit <- pgt(tech, model = "wgd", peers = "group")
  dec <- pgt_decompose(tech, type = "rodseth")
  expect_equal(fit$results$b_star, dec$results$b_star_te,
               tolerance = 1e-10)
})

test_that("wgd agrees with an independent reference implementation", {
  # Independent, self-contained statement of the Rodseth (2025) Eq. 6
  # reduced-form LP, against which the package must agree to the solver
  # tolerance. Guards the packaged kernel against refactoring drift.
  # u_row: the evaluated DMU's own material flow coefficients (one row of
  # the L x N x P array, pollutant 1); v_l: its good-output coefficient.
  ref_wgd <- function(l_prime, X, y, b, u_row, v_l, vrs = TRUE) {
    L <- length(y); N <- ncol(X)
    lp <- lpSolveAPI::make.lp(nrow = 0, ncol = L + 1)
    invisible(lpSolveAPI::lp.control(lp, sense = "min"))
    lpSolveAPI::set.objfn(lp, c(rep(0, L), 1))
    lpSolveAPI::add.constraint(lp, c(y, 0), ">=", y[l_prime])
    for (n in 1:N) {
      lpSolveAPI::add.constraint(lp, c(X[, n], 0), "<=", X[l_prime, n])
    }
    lpSolveAPI::add.constraint(lp, c(b, -1), "<=", 0)
    mb_rhs <- sum(u_row * X[l_prime, ]) - v_l * y[l_prime]
    lpSolveAPI::add.constraint(lp, c(rep(0, L), 1), "<=", mb_rhs)
    if (vrs) lpSolveAPI::add.constraint(lp, c(rep(1, L), 0), "=", 1)
    lpSolveAPI::set.bounds(lp, lower = rep(0, L + 1))
    if (lpSolveAPI::solve.lpExtPtr(lp) != 0) return(NA_real_)
    lpSolveAPI::get.variables(lp)[L + 1]
  }

  tech <- make_random_tech(L = 20, N = 3, seed = 99)
  fit <- pgt(tech, model = "wgd")
  ref <- vapply(seq_len(tech$L), function(i) {
    ref_wgd(i, tech$x, tech$y, tech$b[, 1], tech$u[i, , 1], tech$v[i, 1])
  }, numeric(1))

  expect_equal(fit$results$b_star, ref, tolerance = 1e-8)
})
