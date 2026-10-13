# Replication of the numerical example of Rodseth (2025, JPA 64(3),
# 305-319): Table 1 data, Table 2 results.

pig_tech <- function() {
  pgt_tech(
    x = pigfarms[, c("feed", "piglet", "labor", "capital", "uncontrolled")],
    y = pigfarms$meat,
    b = pigfarms$controlled,
    u = c(0, 0, 0, 0, 1),
    v = 0,
    id = pigfarms$farm
  )
}

test_that("pigfarms data are internally consistent (Table 1)", {
  expect_equal(nrow(pigfarms), 5L)
  expect_equal(pigfarms$farm, c("A", "B", "C", "D", "E"))
  expect_equal(pigfarms$uncontrolled - pigfarms$controlled,
               pigfarms$abatement, tolerance = 1e-12)
})

test_that("mb_check recovers the paper's abatement column", {
  mb <- mb_check(pig_tech())
  expect_equal(mb$gap, pigfarms$abatement, tolerance = 1e-12)
  expect_equal(attr(mb, "n_violations"), 0L)
})

test_that("wgd reproduces minimal controlled emissions of Table 2", {
  fit <- pgt(pig_tech(), model = "wgd")
  expect_equal(fit$results$status, rep(0L, 5))
  # Rodseth (2025), Table 2: Eqs. 6, 7 and 10 all give these minima
  expect_equal(fit$results$b_star, c(16, 16, 16, 20, 16),
               tolerance = 1e-8)
  expect_equal(fit$results$efficiency,
               c(1, 1, 0.8, 1, 16 / 21), tolerance = 1e-8)
})

test_that("the (y,b) envelope agrees on this dataset", {
  fit <- pgt(pig_tech(), model = "envelope")
  expect_equal(fit$results$b_star, c(16, 16, 16, 20, 16),
               tolerance = 1e-8)
})

test_that("the two documented pigfarms mappings give identical minima", {
  # Vignette mapping: uncontrolled aggregate as the only carrier input.
  data(pigfarms, package = "pgt", envir = environment())
  t3 <- pgt_tech(
    x = pigfarms[, c("uncontrolled", "labor", "capital")],
    y = pigfarms$meat, b = pigfarms$controlled,
    u = c(1, 0, 0), a = pigfarms$abatement, id = pigfarms$farm
  )
  # Help-page mapping: all five inputs, u = 0 on the non-carriers.
  t5 <- pgt_tech(
    x = pigfarms[, c("feed", "piglet", "labor", "capital", "uncontrolled")],
    y = pigfarms$meat, b = pigfarms$controlled,
    u = c(0, 0, 0, 0, 1), id = pigfarms$farm
  )
  expect_equal(pgt(t3, model = "wgd")$results$b_star,
               pgt(t5, model = "wgd")$results$b_star, tolerance = 1e-9)
})

test_that("the extended solution reproduces Table 2 of Rodseth (2025)", {
  # Table 2 reports, alongside the minimal controlled emissions, the
  # optimal uncontrolled emissions and abatement of the projection. The
  # paper prints no coefficients; v = 7/6 follows from farms A and B,
  # which share inputs (meat differs by 3, uncontrolled by 3.5). The
  # retained content of farm A's extra meat makes farm B the unique peer
  # of farms B and C; at v = 0, A and B would tie and the targets of B
  # and C would be solver-dependent.
  data(pigfarms, package = "pgt", envir = environment())
  v <- 7 / 6
  # Printed data: each farm's feed coefficient closes its account
  # exactly to the printed uncontrolled column.
  u_feed <- (pigfarms$uncontrolled + v * pigfarms$meat) / pigfarms$feed
  tech <- pgt_tech(
    x = pigfarms[, c("feed", "piglet", "labor", "capital")],
    y = pigfarms$meat, b = pigfarms$controlled,
    u = cbind(u_feed, 0, 0, 0), v = v,
    a = pigfarms$abatement, id = pigfarms$farm
  )
  r <- pgt(tech, model = "wgd")$results
  expect_equal(r$b_star, c(16, 16, 16, 20, 16), tolerance = 1e-8)
  expect_equal(r$z_star, c(16.7, 20.2, 20.2, 20.4, 16.7), tolerance = 1e-8)
  expect_equal(r$a_star, c(0.7, 4.2, 4.2, 0.4, 0.7), tolerance = 1e-8)
  expect_equal(r$z_star - r$a_star, r$b_star, tolerance = 1e-8)
  # Common coefficients consistent with the printed rounding (the
  # mapping of the replication vignette): the targets round to Table 2.
  u_c <- c(feed = 31 / 25, piglet = 349 / 300, labor = 0, capital = 0)
  x_c <- as.matrix(pigfarms[, names(u_c)])
  z_c <- as.vector(x_c %*% u_c) - v * pigfarms$meat
  a_c <- z_c - pigfarms$controlled
  rc <- pgt(pgt_tech(x = x_c, y = pigfarms$meat, b = pigfarms$controlled,
                     u = u_c, v = v, a = a_c, id = pigfarms$farm),
            model = "wgd")$results
  expect_equal(rc$b_star, c(16, 16, 16, 20, 16), tolerance = 1e-8)
  expect_equal(round(rc$z_star, 1), c(16.7, 20.2, 20.2, 20.4, 16.7),
               tolerance = 1e-12)
  expect_equal(round(rc$a_star, 1), c(0.7, 4.2, 4.2, 0.4, 0.7),
               tolerance = 1e-12)
})

test_that("the five-component decomposition telescopes on the pig data", {
  data(pigfarms, package = "pgt", envir = environment())
  tech <- pgt_tech(
    x = pigfarms[, c("uncontrolled", "labor", "capital")],
    y = pigfarms$meat, b = pigfarms$controlled, u = c(1, 0, 0),
    a = pigfarms$abatement, id = pigfarms$farm
  )
  dec <- pgt_decompose(tech, type = "rodseth")
  r <- dec$results
  fit <- pgt(tech, model = "wgd")
  expect_equal(r$total, fit$results$efficiency, tolerance = 1e-8)
  expect_equal(r$total,
               r$te_production * r$quality * r$ae_production *
                 r$te_abatement * r$ae_abatement, tolerance = 1e-10)
  # homogeneous carrier coefficients and no dedicated abatement input:
  # quality and ae_abatement collapse to exactly 1
  expect_equal(r$quality, rep(1, 5))
  expect_equal(r$ae_abatement, rep(1, 5))

  # The same holds under the common coefficients of the replication
  # vignette (v = 7/6), whose decomposition the vignette prints.
  v <- 7 / 6
  u_c <- c(feed = 31 / 25, piglet = 349 / 300, labor = 0, capital = 0)
  x_c <- as.matrix(pigfarms[, names(u_c)])
  a_c <- as.vector(x_c %*% u_c) - v * pigfarms$meat - pigfarms$controlled
  tech_c <- pgt_tech(x = x_c, y = pigfarms$meat, b = pigfarms$controlled,
                     u = u_c, v = v, a = a_c, id = pigfarms$farm)
  rc <- pgt_decompose(tech_c, type = "rodseth")$results
  expect_equal(rc$total, pgt(tech_c, model = "wgd")$results$efficiency,
               tolerance = 1e-8)
  expect_equal(rc$total,
               rc$te_production * rc$quality * rc$ae_production *
                 rc$te_abatement * rc$ae_abatement, tolerance = 1e-10)
  expect_equal(rc$quality, rep(1, 5))
  expect_equal(rc$ae_abatement, rep(1, 5))
})
