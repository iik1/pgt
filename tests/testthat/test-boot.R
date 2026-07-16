# Subsampling inference for pgt efficiency scores.

test_that("boot_pgt intervals extend downward from the point estimate", {
  tech <- make_random_tech(L = 30, N = 3, seed = 2)
  bt <- boot_pgt(tech, model = "wgd", B = 80, seed = 1)
  expect_s3_class(bt, "pgt_boot")
  pd <- bt$per_dmu
  expect_equal(nrow(pd), tech$L)
  # Point estimate matches a plain pgt fit.
  point <- suppressWarnings(pgt(tech, model = "wgd"))$results$efficiency
  expect_equal(pd$estimate, point, tolerance = 1e-9)
  # Recentred, rate-scaled subsampling intervals: a subsample frontier
  # can only lie inside the full-sample frontier, so the interval covers
  # values at or below the (upward-biased) point estimate.
  expect_true(all(pd$lower <= pd$upper + 1e-9, na.rm = TRUE))
  expect_true(all(pd$upper <= pd$estimate + 1e-6, na.rm = TRUE))
  expect_true(all(pd$se >= 0, na.rm = TRUE))
})

test_that("intervals stay ordered when the rescaled deviation exceeds the estimate", {
  # One dominant low-emission peer that rarely enters a size-2
  # subsample: the raw upper endpoint point - rate * q_a(dev) goes
  # negative for the dominated DMUs, so both endpoints must be
  # truncated to [0, 1] together to preserve lower <= upper.
  L <- 200
  x <- matrix(1, L, 1)
  tech <- pgt_tech(x, y = rep(1, L), b = c(0.1, rep(100, L - 1)))
  bt <- boot_pgt(tech, model = "envelope", B = 100, m = 2, seed = 1)
  pd <- bt$per_dmu
  expect_true(all(pd$lower <= pd$upper + 1e-9, na.rm = TRUE))
  expect_true(all(pd$lower >= 0 & pd$upper >= 0, na.rm = TRUE))
})

test_that("boot_pgt rescales by the convergence rate through kappa", {
  tech <- make_random_tech(L = 25, N = 2, seed = 5)
  # (m/L)^kappa shrinks with kappa, so a near-zero exponent gives wider
  # intervals and larger rescaled standard errors than a large one.
  bt_wide <- boot_pgt(tech, model = "wgd", B = 40, seed = 1, kappa = 1e-6)
  bt_narrow <- boot_pgt(tech, model = "wgd", B = 40, seed = 1, kappa = 1)
  w_wide <- bt_wide$per_dmu$upper - bt_wide$per_dmu$lower
  w_narrow <- bt_narrow$per_dmu$upper - bt_narrow$per_dmu$lower
  expect_true(median(w_wide, na.rm = TRUE) >
                median(w_narrow, na.rm = TRUE))
  expect_true(mean(bt_wide$per_dmu$se, na.rm = TRUE) >
                mean(bt_narrow$per_dmu$se, na.rm = TRUE))
})

test_that("boot_pgt restores the caller's RNG state", {
  tech <- make_random_tech(L = 15, N = 2, seed = 3)
  set.seed(123)
  r_direct <- rnorm(3)
  set.seed(123)
  invisible(boot_pgt(tech, model = "envelope", B = 5, seed = 99))
  r_after <- rnorm(3)
  expect_identical(r_direct, r_after)
})

test_that("boot_pgt is reproducible under a seed and produces group means", {
  tech <- make_random_tech(L = 24, N = 2, seed = 6)
  bt1 <- boot_pgt(tech, model = "wgd", B = 60, seed = 42)
  bt2 <- boot_pgt(tech, model = "wgd", B = 60, seed = 42)
  expect_equal(bt1$per_dmu$lower, bt2$per_dmu$lower)
  expect_false(is.null(bt1$group_means))
  expect_true(all(bt1$group_means$lower <= bt1$group_means$upper + 1e-9))
})

test_that("boot_pgt works for by-production and mb_cost models", {
  tech <- make_random_tech(L = 20, N = 2, seed = 9)
  for (m in c("byprod", "mb_cost", "wd")) {
    bt <- boot_pgt(tech, model = m, B = 40, seed = 3)
    expect_equal(nrow(bt$per_dmu), tech$L)
    expect_true(all(bt$per_dmu$estimate > 0 &
                      bt$per_dmu$estimate <= 1 + 1e-9, na.rm = TRUE))
  }
})

test_that("boot_pgt validates m", {
  tech <- make_random_tech(L = 15, seed = 1)
  expect_error(boot_pgt(tech, m = 1), "2..L")
  expect_error(boot_pgt(tech, m = 99), "2..L")
})

test_that("boot_pgt_sensitivity reports width across a grid of m", {
  tech <- make_random_tech(L = 40, N = 2, seed = 4)
  s <- boot_pgt_sensitivity(tech, model = "wgd", B = 30, seed = 1)
  expect_true(all(c("m", "median_width", "mean_se") %in% names(s)))
  expect_true(nrow(s) >= 3)
  expect_true(all(s$median_width >= 0))
})

test_that("group-mode subsampling responds to m", {
  # The subsample size must drive the group-mode resampling, so a
  # sensitivity sweep is not flat.
  tech <- make_random_tech(L = 40, N = 2, groups = c("A", "B"), seed = 4)
  s <- boot_pgt_sensitivity(tech, model = "wgd", peers = "group",
                            m_grid = c(6, 12, 25, 38), B = 40, seed = 1)
  expect_true(length(unique(round(s$median_width, 6))) > 1)
})
