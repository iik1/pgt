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
