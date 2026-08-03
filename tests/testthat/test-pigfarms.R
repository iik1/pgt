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
  # optimal uncontrolled emissions and abatement of the projection.
  # For farms A, D and E the optimum is unique (a peer vertex), so the
  # values are coefficient-robust replication targets; B and C have
  # degenerate optima (A and B tie), so only b* = z* - a* is asserted.
  data(pigfarms, package = "pgt", envir = environment())
  tech <- pgt_tech(
    x = pigfarms[, c("uncontrolled", "labor", "capital")],
    y = pigfarms$meat, b = pigfarms$controlled, u = c(1, 0, 0),
    a = pigfarms$abatement, id = pigfarms$farm
  )
  fit <- pgt(tech, model = "wgd")
  r <- fit$results
  expect_equal(r$b_star, c(16, 16, 16, 20, 16), tolerance = 1e-8)
  expect_equal(r$z_star - r$a_star, r$b_star, tolerance = 1e-8)
  expect_equal(r$z_star[c(1, 4, 5)], c(16.7, 20.4, 16.7),
               tolerance = 1e-8)
  expect_equal(r$a_star[c(1, 4, 5)], c(0.7, 0.4, 0.7),
               tolerance = 1e-8)
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
})
