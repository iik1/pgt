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

test_that("rodseth decomposition stages are strict relaxations", {
  tech <- make_random_tech(L = 40, N = 3, seed = 21)
  dec <- pgt_decompose(tech, type = "rodseth")
  r <- dec$results

  ok <- !is.na(r$total)
  expect_true(all(ok))
  expect_true(all(r$b_star_te[ok] >= r$b_star_technology[ok] - 1e-8))
  expect_true(all(r$b_star_technology[ok] >= r$b_star_ae[ok] - 1e-8))
  expect_equal(r$total[ok], (r$te * r$technology * r$ae)[ok],
               tolerance = 1e-10)
  expect_true(all(r$te[ok] > 0 & r$te[ok] <= 1 + 1e-8))
  expect_true(all(r$technology[ok] > 0 & r$technology[ok] <= 1 + 1e-8))
  expect_true(all(r$ae[ok] > 0 & r$ae[ok] <= 1 + 1e-8))
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

test_that("rodseth partial-stage infeasibility keeps reference semantics", {
  # DMU1 (group A) violates its cap (b = 12 > u'x = 10): stage 1 (group
  # peers) is infeasible, but stages 2-3 solve through group B's DMU3.
  # Matching the reference implementation, later-stage ratios and total
  # stay defined while te/technology are NA.
  x <- matrix(c(10, 10, 10), 3, 1)
  tech <- pgt_tech(x, y = c(10, 5, 10), b = c(12, 9, 8),
                   group = c("A", "A", "B"))
  expect_warning(dec <- pgt_decompose(tech, type = "rodseth"),
                 "infeasible")
  r <- dec$results
  expect_true(is.na(r$te[1]))
  expect_true(is.na(r$technology[1]))
  expect_equal(r$ae[1], 1, tolerance = 1e-8)
  expect_equal(r$total[1], 8 / 12, tolerance = 1e-8)
  # identity holds for complete rows
  cc <- stats::complete.cases(r[c("te", "technology", "ae")])
  expect_equal(r$total[cc], (r$te * r$technology * r$ae)[cc],
               tolerance = 1e-10)
  # summary counts the partially-NA DMU
  s <- summary(dec)
  expect_equal(s$by_group$n_na[s$by_group$group == "A"], 1L)
})

test_that("decomposition requires a group", {
  x <- matrix(1:4, 2, 2)
  tech <- pgt_tech(x, y = c(1, 2), b = c(1, 1))
  expect_error(pgt_decompose(tech), "requires a 'group'")
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
