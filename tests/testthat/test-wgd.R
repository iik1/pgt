test_that("wgd model solves a hand-computable example", {
  # Two DMUs, identical inputs and output, different emissions.
  # The dirty DMU is benchmarked against the clean one: b* = 5.
  x <- matrix(c(10, 10), 2, 1)
  tech <- pgt_tech(x, y = c(5, 5), b = c(10, 5))
  fit <- pgt(tech, model = "wgd")

  expect_equal(fit$results$status, c(0L, 0L))
  expect_equal(fit$results$b_star, c(5, 5), tolerance = 1e-8)
  expect_equal(fit$results$efficiency, c(0.5, 1), tolerance = 1e-8)
  expect_equal(fit$weights[[1]], c(`2` = 1), tolerance = 1e-8)
})

test_that("wgd is always self-feasible, also for violators", {
  # Equation 6 fixes only the intended output, so every unit's
  # programme is feasible; a materials-balance violation does not
  # affect the faithful wgd score (audit it with mb_check()).
  x <- matrix(c(1, 1), 2, 1)
  tech <- pgt_tech(x, y = c(10, 1), b = c(2, 0.5))
  fit <- pgt(tech, model = "wgd")
  expect_equal(fit$results$status, c(0L, 0L))
  expect_equal(fit$results$efficiency, c(1, 1), tolerance = 1e-8)
  expect_true(any(mb_check(tech)$violated))
})

test_that("the v-penalty distinguishes wgd from the plain envelope", {
  # Evaluating DMU A (y = 1, b = 1) with v = 0.5 against B
  # (y = 10, b = 0.5): the envelope ignores retained content and picks
  # B (b* = 0.5); Eq. 6 prices the output overshoot at v and keeps A
  # (b* = 1 + 0.5*1 - 0.5*1 = 1).
  x <- matrix(c(10, 10), 2, 1)
  tech <- pgt_tech(x, y = c(1, 10), b = c(1, 0.5), v = 0.5)
  env <- pgt(tech, model = "envelope")
  wgd <- pgt(tech, model = "wgd")
  expect_equal(env$results$b_star[1], 0.5, tolerance = 1e-8)
  expect_equal(wgd$results$b_star[1], 1, tolerance = 1e-8)
  expect_equal(wgd$results$efficiency[1], 1, tolerance = 1e-8)
})

test_that("wgd equals the envelope when v = 0", {
  tech <- make_random_tech(L = 30, N = 3, seed = 21)
  tech0 <- pgt_tech(tech$x, tech$y[, 1], tech$b[, 1], u = tech$u[1, , 1],
                    v = 0, group = tech$group)
  expect_equal(pgt(tech0, model = "wgd")$results$b_star,
               pgt(tech0, model = "envelope")$results$b_star,
               tolerance = 1e-8)
})

test_that("wgd handles several intended outputs", {
  # C (y = (2,2), b = 4) projects onto the A-B midpoint
  # (y = (2.5, 2.5), b = 2): both output rows hold and b* = 2.
  x <- matrix(10, 3, 1)
  Y <- cbind(elec = c(4, 1, 2), heat = c(1, 4, 2))
  tech <- pgt_tech(x, y = Y, b = c(2, 2, 4))
  expect_equal(tech$M, 2L)
  fit <- pgt(tech, model = "wgd")
  expect_equal(fit$results$b_star, c(2, 2, 2), tolerance = 1e-8)
  expect_equal(fit$results$efficiency[3], 0.5, tolerance = 1e-8)
  expect_true(all(c("elec", "heat", "dual_elec", "dual_heat") %in%
                    names(fit$results)))
})

test_that("envelope model matches the convex lower (y,b) envelope", {
  x <- matrix(10, 3, 1)
  tech <- pgt_tech(x, y = c(1, 3, 2), b = c(1, 3, 4))
  fit <- pgt(tech, model = "envelope")

  expect_equal(fit$results$b_star, c(1, 3, 2), tolerance = 1e-8)
  expect_equal(fit$results$efficiency, c(1, 1, 0.5), tolerance = 1e-8)
  expect_equal(fit$results$dual_output[3], 1, tolerance = 1e-6)
})

test_that("vrs and crs differ as expected for the envelope", {
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
  y1 <- tech$y[, 1]
  expected <- y1 * min(tech$b[, 1] / y1)
  expect_equal(crs$results$b_star, expected, tolerance = 1e-6)
})

test_that("self-reference bounds envelope scores in (0, 1]", {
  tech <- make_random_tech(L = 50, N = 4, seed = 7)
  fit <- pgt(tech, model = "envelope")
  expect_true(all(fit$results$status == 0))
  expect_true(all(fit$results$efficiency > 0))
  expect_true(all(fit$results$efficiency <= 1 + 1e-8))
})

test_that("wgd scores lie in (0, 1] on consistent data", {
  tech <- make_random_tech(L = 40, N = 3, seed = 13)
  fit <- pgt(tech, model = "wgd")
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
  expect_true(all(crs$results$b_star <= vrs$results$b_star + 1e-8))
})

test_that("wgd_anchored keeps the cap and infeasibility semantics", {
  # u'x = 10, v*y = 5, so the cap is 5. With b = (6, 5.5) no peer mix
  # fits under any DMU's cap: every wgd_anchored LP must be infeasible.
  x <- matrix(c(10, 10), 2, 1)
  tech <- pgt_tech(x, y = c(10, 10), b = c(6, 5.5), v = 0.5)
  expect_warning(fit <- pgt(tech, model = "wgd_anchored"), "infeasible")
  expect_true(all(fit$results$status != 0))
  expect_true(all(is.na(fit$results$b_star)))
  mb <- mb_check(tech)
  expect_true(all(mb$gap[fit$results$status != 0] < 0))
  expect_equal(attr(mb, "n_exact"), 2L)

  # Companion: DMU1 violates its cap (b = 6 > 5) but solves through the
  # peer mix lambda = DMU2 (b = 4 <= 5): a feasible violator.
  tech2 <- pgt_tech(x, y = c(10, 10), b = c(6, 4), v = 0.5)
  fit2 <- pgt(tech2, model = "wgd_anchored")
  expect_equal(fit2$results$status, c(0L, 0L))
  expect_equal(fit2$results$b_star, c(4, 4), tolerance = 1e-8)
  expect_equal(fit2$results$efficiency[1], 4 / 6, tolerance = 1e-8)
  expect_equal(fit2$results$mb_headroom, c(1, 1), tolerance = 1e-8)
})

test_that("wgd_anchored diagnostics carry the documented signs", {
  tech <- make_random_tech(L = 30, N = 3, seed = 77)
  fit <- pgt(tech, model = "wgd_anchored")
  ok <- fit$results$status == 0
  expect_true(all(fit$results$dual_output[ok] >= -1e-10))
  expect_true(all(fit$results$mb_headroom[ok] >= -1e-8))
})

test_that("wgd_anchored never beats the faithful wgd", {
  # The wgd_anchored programme adds input rows and the cap to a
  # programme whose remaining rows coincide with Eq. 6 only at v = 0,
  # so the comparison is made there: extra constraints cannot lower b*.
  tech <- make_random_tech(L = 30, N = 3, seed = 31)
  tech0 <- pgt_tech(tech$x, tech$y[, 1], tech$b[, 1], v = 0,
                    group = tech$group)
  wgd <- pgt(tech0, model = "wgd")
  anch <- pgt(tech0, model = "wgd_anchored")
  ok <- anch$results$status == 0
  expect_true(all(wgd$results$b_star[ok] <=
                    anch$results$b_star[ok] + 1e-8))
})

test_that("wgd with group peers restricts the reference set", {
  x <- matrix(10, 3, 1)
  tech <- pgt_tech(x, y = c(5, 5, 5), b = c(10, 5, 1),
                   group = c("A", "A", "B"), id = c("d1", "d2", "d3"))
  fit <- pgt(tech, model = "wgd", peers = "group")
  expect_equal(fit$results$b_star, c(5, 5, 1), tolerance = 1e-8)
  expect_equal(names(fit$weights[["d1"]]), "d2")
  expect_equal(names(fit$weights[["d3"]]), "d3")
})

test_that("the stage kernel reproduces the reduced form (Proposition 1)", {
  # The fully freed stage programme of the extended representation
  # (Eq. 9) must return the same minimal emissions as the reduced-form
  # wgd kernel, with and without producer-specific coefficients.
  tech <- make_random_tech(L = 20, N = 3, seed = 61)
  fit <- pgt(tech, model = "wgd")
  b5 <- vapply(seq_len(tech$L), function(i) {
    pgt:::.lp_wgd_stage(i, tech, seq_len(tech$L), vrs = TRUE, p = 1L,
                        hold_xp = FALSE, hold_xa = FALSE,
                        hold_a = FALSE, hold_quality = FALSE)$b_star
  }, numeric(1))
  expect_equal(b5, fit$results$b_star, tolerance = 1e-7)

  U <- matrix(runif(tech$L * tech$N, 0.5, 2), tech$L, tech$N)
  vv <- runif(tech$L, 0, 0.3)
  het <- pgt_tech(tech$x, tech$y[, 1], tech$b[, 1], u = U, v = vv)
  fit_h <- pgt(het, model = "wgd")
  b5_h <- vapply(seq_len(het$L), function(i) {
    pgt:::.lp_wgd_stage(i, het, seq_len(het$L), vrs = TRUE, p = 1L,
                        hold_xp = FALSE, hold_xa = FALSE,
                        hold_a = FALSE, hold_quality = FALSE)$b_star
  }, numeric(1))
  expect_equal(b5_h, fit_h$results$b_star, tolerance = 1e-7)
})

test_that("wgd agrees with an independent reference implementation", {
  # Independent statement of the Eq. 6 reduced form: minimise
  # sum lambda (b_l + v_i' y_l) - v_i' y_i subject to the output rows
  # and VRS. Guards the packaged kernel against refactoring drift.
  ref_wgd <- function(l_prime, Y, b, v_row, vrs = TRUE) {
    L <- nrow(Y)
    lp <- lpSolveAPI::make.lp(nrow = 0, ncol = L)
    invisible(lpSolveAPI::lp.control(lp, sense = "min"))
    lpSolveAPI::set.objfn(lp, b + as.vector(Y %*% v_row))
    for (m in seq_len(ncol(Y))) {
      lpSolveAPI::add.constraint(lp, Y[, m], ">=", Y[l_prime, m])
    }
    if (vrs) lpSolveAPI::add.constraint(lp, rep(1, L), "=", 1)
    lpSolveAPI::set.bounds(lp, lower = rep(0, L))
    if (lpSolveAPI::solve.lpExtPtr(lp) != 0) return(NA_real_)
    lpSolveAPI::get.objective(lp) - sum(v_row * Y[l_prime, ])
  }

  tech <- make_random_tech(L = 20, N = 3, seed = 99)
  fit <- pgt(tech, model = "wgd")
  ref <- vapply(seq_len(tech$L), function(i) {
    ref_wgd(i, tech$y, tech$b[, 1],
            as.vector(tech$v[i, , 1]))
  }, numeric(1))

  expect_equal(fit$results$b_star, ref, tolerance = 1e-8)
})
