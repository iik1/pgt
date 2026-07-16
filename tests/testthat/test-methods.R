test_that("fit methods run and return the right shapes", {
  tech <- make_random_tech(L = 30, seed = 17)
  fit <- pgt(tech, model = "wgd")

  expect_output(print(fit), "pgt fit")
  s <- summary(fit)
  expect_s3_class(s, "summary.pgt")
  expect_output(print(s), "Efficiency")
  expect_equal(nrow(s$by_group), nlevels(tech$group))

  df <- as.data.frame(fit)
  expect_equal(nrow(df), tech$L)
  expect_true(all(c("id", "group", "b_star", "efficiency") %in% names(df)))

  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  plot(fit)
  grDevices::dev.off()
  unlink(f)
})

test_that("shadow_prices returns the diagnostics", {
  tech <- make_random_tech(L = 25, seed = 23)
  fit <- pgt(tech, model = "wgd")
  sp <- shadow_prices(fit)
  expect_equal(nrow(sp), tech$L)
  expect_true(all(c("id", "group", "dual_output", "mb_headroom") %in%
                    names(sp)))

  env <- pgt(tech, model = "envelope")
  expect_true(all(is.na(shadow_prices(env)$mb_headroom)))
})

test_that("empty MAC curves print and refuse to plot gracefully", {
  # both DMUs infeasible: the curve is empty, everything excluded
  x <- matrix(c(10, 10), 2, 1)
  tech <- pgt_tech(x, y = c(10, 10), b = c(6, 5.5), v = 0.5)
  suppressWarnings(fit <- pgt(tech, model = "wgd"))
  mac <- mac_curve(fit)
  expect_equal(nrow(mac), 0L)
  expect_equal(attr(mac, "n_excluded"), 2L)
  expect_output(print(mac), "no DMUs with a positive output dual")
  expect_error(plot(mac), "empty MAC curve")
})

test_that("plot.pgt refuses an all-NA fit", {
  x <- matrix(c(10, 10), 2, 1)
  tech <- pgt_tech(x, y = c(10, 10), b = c(6, 5.5), v = 0.5)
  suppressWarnings(fit <- pgt(tech, model = "wgd"))
  expect_error(plot(fit), "all scores are NA")
})

test_that("mac_curve orders and accumulates correctly", {
  tech <- make_random_tech(L = 40, seed = 29)
  fit <- pgt(tech, model = "envelope")

  mac <- mac_curve(fit)
  expect_s3_class(mac, "pgt_mac")
  expect_true(!is.unsorted(mac$mac))
  expect_true(!is.unsorted(mac$cum_abatement))
  expect_equal(max(mac$cum_abatement), sum(mac$abatement), tolerance = 1e-10)
  expect_equal(nrow(mac) + attr(mac, "n_excluded"), tech$L)

  mac_p <- mac_curve(fit, price = 100)
  expect_equal(mac_p$mac, 100 * mac$mac, tolerance = 1e-10)

  expect_output(print(mac), "Marginal abatement cost")
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  plot(mac)
  grDevices::dev.off()
  unlink(f)
})

test_that("mac_curve validates price", {
  tech <- make_random_tech(seed = 2)
  fit <- pgt(tech, model = "envelope")
  expect_error(mac_curve(fit, price = -1))
  expect_error(mac_curve(fit, price = c(1, 2)))
})

test_that("mac_curve and shadow_prices reject models without output duals", {
  tech <- make_random_tech(L = 12, N = 2, seed = 37)
  fit <- pgt(tech, model = "mb_cost")
  expect_error(mac_curve(fit), "output duals")
  expect_error(shadow_prices(fit), "output duals")
})

test_that("the pollutant line is printed for every model and by summary", {
  tech <- make_random_tech(L = 12, N = 2, seed = 41)
  env <- pgt(tech, model = "envelope")
  expect_output(print(env), "pollutant")
  expect_output(print(summary(env)), "pollutant")
})
