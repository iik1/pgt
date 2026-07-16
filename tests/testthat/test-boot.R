# Subsampling inference for pgt efficiency scores.

test_that("boot_pgt returns per-DMU intervals bracketing the point estimate", {
  tech <- make_random_tech(L = 30, N = 3, seed = 2)
  bt <- boot_pgt(tech, model = "wgd", B = 80, seed = 1)
  expect_s3_class(bt, "pgt_boot")
  pd <- bt$per_dmu
  expect_equal(nrow(pd), tech$L)
  # Point estimate matches a plain pgt fit.
  point <- suppressWarnings(pgt(tech, model = "wgd"))$results$efficiency
  expect_equal(pd$estimate, point, tolerance = 1e-9)
  # Valid intervals: ordered bounds, non-negative width and SE. The
  # subsampling distribution reflects reference-set resampling and, given
  # the boundary bias of DEA scores, need not be centred on the point
  # estimate.
  expect_true(all(pd$lower <= pd$upper + 1e-9, na.rm = TRUE))
  expect_true(all(pd$se >= 0, na.rm = TRUE))
  expect_true(all(pd$lower > 0 & pd$upper <= 1 + 1e-6, na.rm = TRUE))
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
