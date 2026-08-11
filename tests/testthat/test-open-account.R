test_that("z_star uses the accounting-cap convention on open accounts", {
  # Deliberately open accounts: with u = 1 and v = 0 every peer's cap
  # u'x - v'y = x differs from b + a (D1: 30 vs 10, ..., D6: 18 vs 5),
  # so the two candidate conventions for the implied uncontrolled
  # emission disagree and the fit pins down the accounting-cap
  # convention against the closed-account formula sum_l lambda_l
  # (b_l + a_l).
  x <- matrix(c(30, 26, 24, 22, 20, 18), ncol = 1)
  y <- c(4, 5, 6, 7, 8, 9)
  b <- c(9, 8, 7, 6, 5, 4)
  a <- rep(1, 6)
  tech <- pgt_tech(x, y = y, b = b, u = 1, v = 0, a = a,
                   id = paste0("D", 1:6))
  expect_true(all(abs(mb_check(tech)$closure) > 1))

  fit <- pgt(tech, model = "wgd")
  r <- fit$results
  expect_equal(r$status, rep(0L, 6))

  # b falls one-for-one in y, so every DMU's optimal peer mix is the
  # unique vertex D6 (y = 9, b = 4): b* = 4, z* = x_6 = 18, a* = 1.
  expect_equal(r$b_star, rep(4, 6), tolerance = 1e-8)
  expect_equal(r$z_star, rep(18, 6), tolerance = 1e-8)
  expect_equal(r$a_star, rep(1, 6), tolerance = 1e-8)

  # Recompute both candidates from the returned peer weights.
  cap <- as.numeric(x)  # u'x - v'y with u = 1, v = 0
  peer_sum <- function(q) vapply(fit$weights, function(w)
    sum(w * q[match(names(w), tech$id)]), numeric(1), USE.NAMES = FALSE)
  expect_equal(r$z_star, peer_sum(cap), tolerance = 1e-8)
  expect_equal(r$a_star, peer_sum(a), tolerance = 1e-8)
  # The closed-account formula sum lambda (b + a) gives 5 here: on an
  # open account it is not what z_star reports, and z* - a* != b*.
  expect_equal(peer_sum(b + a), rep(5, 6), tolerance = 1e-8)
  expect_true(all(abs(r$z_star - peer_sum(b + a)) > 1))
  expect_true(all(abs(r$z_star - r$a_star - r$b_star) > 1))
})
