# Materials-balance cost model (Coelli, Lauwers and Van Huylenbroeck
# 2007) and the weak-disposability reference model (Kuosmanen 2005).

test_that("mb_cost satisfies the EE = TE x EAE decomposition identity", {
  tech <- make_random_tech(L = 20, N = 3, seed = 11)
  r <- pgt(tech, model = "mb_cost")$results
  expect_equal(r$efficiency, r$te * r$eae, tolerance = 1e-9)
  expect_true(all(r$efficiency <= r$te + 1e-9))          # EE <= TE
  expect_true(all(r$eae > 0 & r$eae <= 1 + 1e-9))        # EAE in (0, 1]
  expect_true(all(r$te > 0 & r$te <= 1 + 1e-9))
})

test_that("mb_cost matches a hand-computed two-input example", {
  # Two DMUs, one input mix each; CRS. DMU2 uses the high-nitrogen input
  # heavily and can cut material inflow by moving toward DMU1's mix.
  x <- rbind(c(4, 1), c(1, 4))
  colnames(x) <- c("clean", "dirty")
  y <- c(1, 1)
  b <- c(3, 6)
  u <- c(clean = 1, dirty = 2)     # dirty input carries twice the pollutant
  tech <- pgt_tech(x, y, b, u = u, v = 0)
  r <- pgt(tech, model = "mb_cost", returns = "crs")$results
  # Observed material inflow: DMU1 = 4*1 + 1*2 = 6; DMU2 = 1*1 + 4*2 = 9.
  # Both produce y = 1; the least-material way to make y = 1 (CRS) is
  # DMU1's mix at unit scale, inflow 6. So EE1 = 6/6 = 1, EE2 = 6/9 = 2/3.
  expect_equal(r$efficiency[1], 1, tolerance = 1e-7)
  expect_equal(r$efficiency[2], 2 / 3, tolerance = 1e-7)
})

test_that("mb_cost uses Coelli et al. (2007) phosphorus coefficients", {
  # Coelli, Lauwers & Van Huylenbroeck (2007), footnote 16: phosphorus
  # contents (kg P per kg) of feed 0.0124, piglet 0.0117; output meat
  # 0.0117. A minimal three-farm check that the coefficients flow through
  # the material-inflow objective.
  x <- rbind(c(200, 20), c(210, 22), c(190, 18))
  colnames(x) <- c("feed", "piglet")
  y <- c(100, 100, 95)
  b <- c(1.5, 1.6, 1.4)
  tech <- pgt_tech(x, y, b, u = c(feed = 0.0124, piglet = 0.0117),
                   v = 0.0117)
  r <- pgt(tech, model = "mb_cost", returns = "crs")$results
  expect_true(all(r$efficiency > 0 & r$efficiency <= 1 + 1e-9))
  expect_true(all(r$efficiency <= r$te + 1e-9))
})

test_that("mb_cost returns NA (not NaN) for a zero-potential DMU", {
  # DMU2 uses none of the only pollutant-bearing input (coal), so its
  # aggregate material inflow is zero and the ratio is undefined.
  x <- rbind(c(2, 3), c(0, 4), c(3, 2))
  colnames(x) <- c("coal", "scrap")
  tech <- pgt_tech(x, y = c(5, 5, 5), b = c(4, 2, 3),
                   u = c(coal = 1, scrap = 0), v = 0)
  suppressWarnings(r <- pgt(tech, model = "mb_cost", returns = "crs")$results)
  expect_true(is.na(r$efficiency[2]))
  expect_false(is.nan(r$efficiency[2]))
  expect_true(r$status[2] != 0)                 # flagged, not silently solved
  expect_true(all(is.finite(r$efficiency[c(1, 3)])))
})

test_that("wd gives emission efficiency in (0, 1]", {
  tech <- make_random_tech(L = 18, N = 2, seed = 5)
  r <- pgt(tech, model = "wd")$results
  expect_true(all(r$efficiency > 0 & r$efficiency <= 1 + 1e-9))
  expect_true(all(r$status == 0))
})

test_that("wd frontier point scores 1 on its own emission", {
  # The minimum-emission-per-output DMU is weakly efficient.
  x <- matrix(c(10, 10, 10), ncol = 1)
  y <- c(5, 5, 5)
  b <- c(2, 4, 6)
  tech <- pgt_tech(x, y, b)
  r <- pgt(tech, model = "wd", returns = "vrs")$results
  expect_equal(r$efficiency[1], 1, tolerance = 1e-7)   # lowest emission
  expect_true(r$efficiency[3] < 1)                     # highest emission
})

test_that("mb_cost warns when the implied minimal emission is negative", {
  # The peer's low material inflow undercuts DMU 1's retained content
  # v * y, so pot_star - v * y < 0.
  x <- matrix(c(2, 0.5), 2, 1)
  tech <- pgt_tech(x, y = c(1, 2), b = c(0.1, 0.1), v = 0.9)
  expect_warning(fit <- pgt(tech, model = "mb_cost"),
                 "negative implied minimal emission")
  expect_true(fit$results$b_star[1] < 0)
})

test_that("a failed sub-LP yields all-NA scores under nonzero status", {
  # T2 of the by-production model is infeasible for DMU 2 against peer 1
  # alone (peer's polluting input is below DMU 2's), while T1 solves;
  # the whole intersection measure must then be NA.
  X <- matrix(c(1, 5), 2, 1)
  sol <- pgt:::.lp_byprod_one(2, X, y = c(1, 1), b = c(1, 1), pol = 1L,
                              peers = 1L, vrs = TRUE)
  expect_true(sol$status != 0)
  expect_true(is.na(sol$output_eff))
  expect_true(is.na(sol$emission_eff))
  expect_true(is.na(sol$fgl))
  expect_true(is.na(sol$b_star))
})
