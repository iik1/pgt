# Solver robustness at tonne magnitudes. steeldemo is in tonnes
# (quantities of order 1e5-1e7); multiplied by 1e6 its quantities reach
# 1e13. Solved as given, lp_solve stopped the decomposition's first
# stage at suboptimal vertices with status 0 (96 of 180 plants at 1e6,
# two with production technical efficiency above 1) and reported
# self-feasible wgd_input_fixed programmes as infeasible. Every
# programme is homogeneous of degree one in the quantities, so a common
# rescaling must leave every score unchanged.

steel_k <- function(k, b = pgt::steeldemo$emissions) {
  d <- pgt::steeldemo
  pgt_tech(
    x = as.matrix(d[, c("coal_coke", "other_fuel", "raw_material", "flux")]) * k,
    y = d$production * k, b = b * k, v = 0.01467,
    group = d$route, id = d$plant
  )
}

# Stage 1 of pgt_decompose(type = "rodseth") (every input held, uniform
# coefficients, no abatement output) reduces by the summing-up condition
# to
#   b1_i = min sum_l lambda_l (b_l - u'x_l + v y_l) + u'x_i - v y_i
#   s.t. sum_l lambda_l x_l <= x_i, sum_l lambda_l y_l >= y_i,
#        sum_l lambda_l = 1, lambda >= 0,
# solved here on a copy in million tonnes (where it agrees with GLPK to
# 5e-9 of b) and returned at magnitude k.
stage1_reduced <- function(k) {
  d <- pgt::steeldemo
  x <- as.matrix(d[, c("coal_coke", "other_fuel", "raw_material", "flux")]) / 1e6
  y <- d$production / 1e6
  b <- d$emissions / 1e6
  v <- 0.01467
  L <- nrow(x)
  k * 1e6 * vapply(seq_len(L), function(i) {
    lp <- lpSolveAPI::make.lp(0, L)
    lpSolveAPI::set.objfn(lp, b - rowSums(x) + v * y)
    for (n in seq_len(ncol(x))) {
      lpSolveAPI::add.constraint(lp, x[, n], "<=", x[i, n])
    }
    lpSolveAPI::add.constraint(lp, y, ">=", y[i])
    lpSolveAPI::add.constraint(lp, rep(1, L), "=", 1)
    stopifnot(lpSolveAPI::solve.lpExtPtr(lp) == 0)
    lpSolveAPI::get.objective(lp) + sum(x[i, ]) - v * y[i]
  }, numeric(1))
}

test_that("decomposition stage 1 attains its reduced-form optimum at 1e6 magnitude", {
  for (k in c(1, 1e6)) {
    r <- pgt_decompose(steel_k(k), type = "rodseth")$results
    b1 <- r$te_production * r$b
    expect_false(anyNA(b1))
    expect_true(all(r$te_production <= 1 + 1e-8))
    expect_lt(max(abs(b1 - stage1_reduced(k)) / r$b), 1e-7)
  }
})

test_that("scores are invariant to a common rescaling of the quantities", {
  small <- steel_k(1e-6)
  big <- steel_k(1e6)
  for (m in c("wgd", "wgd_input_fixed", "envelope")) {
    f_small <- pgt(small, model = m)$results
    f_big <- pgt(big, model = m)$results
    expect_equal(f_big$status, f_small$status, label = m)
    expect_equal(f_big$efficiency, f_small$efficiency, tolerance = 1e-8,
                 label = m)
  }
  d_small <- pgt_decompose(small, type = "rodseth")$results
  d_big <- pgt_decompose(big, type = "rodseth")$results
  for (cc in c(pgt:::.decomp_components("rodseth"), "total")) {
    expect_equal(d_big[[cc]], d_small[[cc]], tolerance = 1e-8, label = cc)
  }
})

test_that("wgd_input_fixed solves every unit that meets its cap at 1e6 magnitude", {
  # steeldemo's accounts close, so self-reference is feasible for every
  # plant: no programme may be infeasible or score above 1
  fit <- pgt(steel_k(1e6), model = "wgd_input_fixed")$results
  expect_equal(fit$status, rep(0L, nrow(fit)))
  expect_true(all(fit$efficiency <= 1 + 1e-8))

  # 40 plants pushed 30% above their cap: the others must still solve,
  # and a violator's programme ends optimal or infeasible, never failed
  d <- pgt::steeldemo
  cap <- rowSums(d[, c("coal_coke", "other_fuel", "raw_material", "flux")]) -
    0.01467 * d$production
  set.seed(3)
  viol <- sample(nrow(d), 40)
  b <- d$emissions
  b[viol] <- 1.3 * cap[viol]
  st <- suppressWarnings(
    pgt(steel_k(1e6, b), model = "wgd_input_fixed", peers = "group")
  )$results$status
  expect_equal(st[-viol], rep(0L, nrow(d) - 40L))
  expect_true(all(st[viol] %in% c(0L, 2L)))
  st_small <- suppressWarnings(
    pgt(steel_k(1e-6, b), model = "wgd_input_fixed", peers = "group")
  )$results$status
  expect_equal(st, st_small)
})

test_that("a solve above the self-reference bound is retried, then flagged", {
  ok <- pgt:::.lp_ok(b_i = 1, self_feasible = TRUE)
  # suboptimal under the default settings, optimal under the first
  # alternative: the retry returns the alternative's solution
  calls <- 0L
  sol <- pgt:::.lp_retry(function(control) {
    calls <<- calls + 1L
    list(status = 0L, b_star = if (is.null(control)) 1.2 else 0.8,
         lambda = c(1, 0))
  }, ok)
  expect_equal(sol$b_star, 0.8)
  expect_equal(calls, 2L)
  # suboptimal under every setting: NA values under lp_solve's
  # SUBOPTIMAL code 1
  sol <- pgt:::.lp_retry(function(control) {
    list(status = 0L, b_star = 1.2, lambda = c(1, 0), dual_output = 0.5)
  }, ok)
  expect_equal(sol$status, 1L)
  expect_true(is.na(sol$b_star))
  expect_true(is.na(sol$dual_output))
  expect_null(sol$lambda)
  # self-feasible: infeasibility is retried and, if every attempt
  # agrees, the first status is kept
  sol <- pgt:::.lp_retry(function(control) {
    list(status = 2L, b_star = NA_real_)
  }, ok)
  expect_equal(sol$status, 2L)
  # not self-feasible: a numerical failure is retried until a definitive
  # status, and infeasibility is accepted
  st <- c(5L, 2L)
  calls <- 0L
  sol <- pgt:::.lp_retry(function(control) {
    calls <<- calls + 1L
    list(status = st[min(calls, 2L)], b_star = NA_real_)
  }, pgt:::.lp_ok(b_i = 1, self_feasible = FALSE))
  expect_equal(sol$status, 2L)
  expect_equal(calls, 2L)
})

test_that("pgt() retries and reports solves that exceed b", {
  x <- matrix(c(10, 10, 10), 3, 1)
  tech <- pgt_tech(x, y = c(5, 5, 5), b = c(10, 5, 8))
  # an envelope kernel that returns 2b under the default settings
  local_mocked_bindings(.lp_envelope_one = function(i, Y, b, peers,
                                                    vrs = TRUE,
                                                    control = NULL) {
    list(status = 0L, b_star = if (is.null(control)) 2 * b[i] else 0.5 * b[i],
         lambda = as.numeric(peers == i), dual_output = 0, mb_rhs = NA_real_)
  })
  fit <- pgt(tech, model = "envelope")
  expect_equal(fit$results$efficiency, rep(0.5, 3))
  # ... and under every setting: status 1 with NA scores and a warning
  local_mocked_bindings(.lp_envelope_one = function(i, Y, b, peers,
                                                    vrs = TRUE,
                                                    control = NULL) {
    list(status = 0L, b_star = 2 * b[i], lambda = as.numeric(peers == i),
         dual_output = 0, mb_rhs = NA_real_)
  })
  expect_warning(fit <- pgt(tech, model = "envelope"), "status 1")
  expect_equal(fit$results$status, rep(1L, 3))
  expect_true(all(is.na(fit$results$efficiency)))
})

# uscoal with the sulfur coefficients scaled by 0.9. Painesville (plant
# 2936, 38 MWh) is the smallest generator and the lowest emitter, so no
# peer mix emits less than its own 6.543 t, above its cap of 6.153 t:
# its input-fixed programme is infeasible. On the unit-magnitude copy its
# cap row is of order 1e-7. 0.7.1 took the plant for one meeting its cap
# (an absolute tolerance of 1e-8 there is 0.78 t), retried lp_solve's
# infeasibility verdict and accepted the first alternative scaling's
# point 6 per cent above the cap with status 0.
uscoal_u <- function(ucoal) {
  d <- pgt::uscoal
  pgt_tech(x = d[, c("coal", "other_heat", "capacity", "sorbent")],
           y = d$gen, b = d$so2,
           u = cbind(coal = ucoal, other_heat = 0, capacity = 0,
                     sorbent = 0),
           v = 0, x_abate = "sorbent", id = d$plant)
}

test_that("wgd_input_fixed reports no point above the cap as optimal", {
  d <- pgt::uscoal
  u <- 0.9 * 2 * d$sulfur / 100
  r <- suppressWarnings(pgt(uscoal_u(u), model = "wgd_input_fixed"))$results
  i <- which(d$plant == "2936")
  expect_equal(r$status[i], 2L)
  expect_true(is.na(r$b_star[i]) && is.na(r$mb_headroom[i]))
  cap <- u * d$coal
  ok <- r$status == 0
  expect_true(all(r$mb_headroom[ok] >= -1e-5 * cap[ok]))
  # the programme is feasible exactly when its optimum without the cap
  # row lies within the cap, and then keeps that optimum; a cap out of
  # reach gives the uncapped programme
  nocap <- pgt(uscoal_u(rep(1e6, nrow(d))), model = "wgd_input_fixed")$results
  expect_equal(r$status != 0, nocap$b_star > cap * (1 + 1e-5))
  expect_equal(r$b_star[ok], nocap$b_star[ok], tolerance = 1e-8)
})

test_that("the input-fixed kernel rejects lp_solve's point above the cap", {
  # the solve that produced the point in 0.7.1: Painesville's programme
  # on the unit-magnitude copy under the first alternative scaling
  d <- pgt::uscoal
  sc <- pgt:::.scale_tech(uscoal_u(0.9 * 2 * d$sulfur / 100))$tech
  ctx <- pgt:::.solve_ctx(sc, "wgd_input_fixed", 1L)
  i <- which(d$plant == "2936")
  expect_false(pgt:::.meets_cap(ctx$b_p[i], ctx$mb_cap[i]))
  sol <- pgt:::.lp_wgd_input_fixed_one(
    i, sc$x, sc$y, ctx$b_p, ctx$mb_cap, seq_len(sc$L),
    control = pgt:::.lp_alt_controls[[1]])
  expect_equal(sol$status, 2L)
  expect_true(is.na(sol$b_star))
})

test_that("the row check measures each row against its own magnitude", {
  lp <- lpSolveAPI::make.lp(0, 2)
  lpSolveAPI::set.objfn(lp, c(1, 1))
  lpSolveAPI::add.constraint(lp, c(1, 0), ">=", 1e-7)
  lpSolveAPI::add.constraint(lp, c(0, 1), ">=", 1)
  expect_equal(lpSolveAPI::solve.lpExtPtr(lp), 0L)
  expect_true(pgt:::.lp_verified(lp, n_w = 0))
  # lpSolveAPI keeps the point when a right-hand side changes: the first
  # row, of order 1e-7, is now missed by 5 per cent, then by 1e-7
  lpSolveAPI::set.rhs(lp, 1.05e-7, 1)
  expect_false(pgt:::.lp_verified(lp, n_w = 0))
  lpSolveAPI::set.rhs(lp, 1e-7 * (1 + 1e-7), 1)
  expect_true(pgt:::.lp_verified(lp, n_w = 0))

  # two peer weights summing to one, the second forced to 5e-9 and
  # carrying an input the unit does not use: once that input's row
  # allows none, it is broken by the small weight alone, which passes
  # when the columns are peer weights and fails when they are not
  lp <- lpSolveAPI::make.lp(0, 2)
  lpSolveAPI::set.objfn(lp, c(0, 1))
  lpSolveAPI::add.constraint(lp, c(1, 1), "=", 1)
  lpSolveAPI::add.constraint(lp, c(0, 1), ">=", 5e-9)
  lpSolveAPI::add.constraint(lp, c(0, 0.5), "<=", 1e-8)
  expect_equal(lpSolveAPI::solve.lpExtPtr(lp), 0L)
  lpSolveAPI::set.rhs(lp, 0, 3)
  expect_true(pgt:::.lp_verified(lp, n_w = 2))
  expect_false(pgt:::.lp_verified(lp, n_w = 0))
})
