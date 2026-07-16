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
