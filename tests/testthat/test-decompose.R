test_that("envelope decomposition satisfies its exact identity", {
  tech <- make_random_tech(L = 60, N = 4, seed = 13)
  dec <- pgt_decompose(tech, type = "envelope")
  r <- dec$results

  expect_false(anyNA(r$total))
  expect_equal(r$total, r$WR * r$TGR, tolerance = 1e-10)
  expect_true(all(r$WR > 0 & r$WR <= 1 + 1e-8))
  expect_true(all(r$TGR > 0 & r$TGR <= 1 + 1e-8))
  # pooled envelope cannot lie above the group envelope
  expect_true(all(r$b_star_all <= r$b_star_group + 1e-8))
  # decomposition agrees with the corresponding pgt() fits
  fit_grp <- pgt(tech, model = "envelope", peers = "group")
  fit_all <- pgt(tech, model = "envelope", peers = "all")
  expect_equal(r$b_star_group, fit_grp$results$b_star, tolerance = 1e-10)
  expect_equal(r$b_star_all, fit_all$results$b_star, tolerance = 1e-10)
})

test_that("Eq. 11 components telescope to the wgd efficiency", {
  tech <- make_random_tech(L = 40, N = 3, seed = 21)
  dec <- pgt_decompose(tech, type = "rodseth")
  r <- dec$results
  comps <- r$te_production * r$quality * r$ae_production *
    r$te_abatement * r$ae_abatement
  expect_false(anyNA(r$total))
  expect_equal(r$total, comps, tolerance = 1e-8)
  fit <- pgt(tech, model = "wgd")
  expect_equal(r$total, fit$results$efficiency, tolerance = 1e-7)
  for (cc in pgt:::.decomp_components("rodseth")) {
    expect_true(all(r[[cc]] > 0 & r[[cc]] <= 1 + 1e-6), label = cc)
  }
  # homogeneous coefficients, no abatement data: the quality and
  # abatement components collapse to exactly 1
  expect_equal(r$quality, rep(1, nrow(r)))
  expect_equal(r$te_abatement, rep(1, nrow(r)))
  expect_equal(r$ae_abatement, rep(1, nrow(r)))
})

test_that("attribution shares sum to one where inefficiency exists", {
  tech <- make_random_tech(L = 50, seed = 31)
  dec <- pgt_decompose(tech, type = "envelope")
  shares <- pgt:::.decomp_attribution(dec$results,
                                      pgt:::.decomp_components("envelope"))
  tot <- rowSums(shares)
  ineff <- (1 - dec$results$total) > 1e-9
  expect_equal(tot[ineff], rep(1, sum(ineff)), tolerance = 1e-8)
})

test_that("quality and abatement components activate with the data", {
  set.seed(9)
  L <- 25
  x <- cbind(fuel = runif(L, 20, 60), sorbent = runif(L, 1, 8),
             labour = runif(L, 30, 70))
  U <- cbind(fuel = runif(L, 0.8, 1.6), sorbent = 0, labour = 0)
  y <- runif(L, 5, 20)
  pot <- rowSums(U * x)
  a <- runif(L, 0.05, 0.25) * pot
  b <- (pot - a) * runif(L, 0.7, 0.98)
  z <- b + a
  tech <- pgt_tech(x, y, b, u = U, a = a, x_abate = "sorbent",
                   id = seq_len(L))
  dec <- pgt_decompose(tech, type = "rodseth")
  r <- dec$results
  comps <- r$te_production * r$quality * r$ae_production *
    r$te_abatement * r$ae_abatement
  expect_equal(r$total, comps, tolerance = 1e-7)
  fit <- pgt(tech, model = "wgd")
  expect_equal(r$total, fit$results$efficiency, tolerance = 1e-7)
  # heterogeneous u activates the quality stage; a and x_abate activate
  # the abatement stages: at least one unit separates each component
  expect_true(any(r$quality < 1 - 1e-8))
  expect_true(any(r$te_abatement < 1 - 1e-8) ||
                any(r$ae_abatement < 1 - 1e-8))
  # stage minima are monotone under the successive relaxations
  expect_true(all(r$te_production >= r$total - 1e-8))
})

test_that("only the envelope decomposition requires a group", {
  x <- matrix(c(10, 10), 2, 1)
  tech <- pgt_tech(x, y = c(1, 2), b = c(1, 1))
  expect_error(pgt_decompose(tech, type = "envelope"),
               "requires a 'group'")
  dec <- pgt_decompose(tech, type = "rodseth")
  expect_s3_class(dec, "pgt_decomp")
  expect_output(print(summary(dec)), "Group medians")
})

test_that("decomposition methods run", {
  tech <- make_random_tech(seed = 8)
  dec <- pgt_decompose(tech, type = "envelope")
  expect_output(print(dec), "envelope")
  s <- summary(dec)
  expect_s3_class(s, "summary.pgt_decomp")
  expect_output(print(s), "Group medians")
  expect_true(all(c("WR", "TGR", "total") %in% names(s$by_group)))

  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  plot(dec)
  grDevices::dev.off()
  unlink(f)
})
