test_that("mb_check computes closure gaps by hand", {
  # potential = u'x = (10, 5); retained = 0; b = (8, 6)
  # gaps = (2, -1): second DMU violates
  x <- matrix(c(10, 5), 2, 1)
  tech <- pgt_tech(x, y = c(1, 1), b = c(8, 6))
  mb <- mb_check(tech)

  expect_equal(mb$potential, c(10, 5))
  expect_equal(mb$gap, c(2, -1))
  expect_equal(mb$violated, c(FALSE, TRUE))
  expect_equal(attr(mb, "n_violations"), 1L)
  expect_output(print(mb), "violations: 1")
})

test_that("mb_check respects v and u", {
  x <- matrix(c(10, 10), 2, 1)
  tech <- pgt_tech(x, y = c(100, 200), b = c(4, 4), u = 2, v = 0.1)
  mb <- mb_check(tech)
  # potential = 20; retained = (10, 20); gap = 20 - retained - 4
  expect_equal(mb$gap, c(6, -4))
  expect_equal(attr(mb, "n_violations"), 1L)
})

test_that("clean technologies pass", {
  tech <- make_random_tech()
  mb <- mb_check(tech)
  expect_equal(attr(mb, "n_violations"), 0L)
  expect_true(all(mb$gap >= 0))
})

test_that("printing a subset of the audit reports the subset's counts", {
  x <- matrix(c(10, 5), 2, 1)
  tech <- pgt_tech(x, y = c(1, 1), b = c(8, 6))
  mb <- mb_check(tech)
  sub <- mb[!mb$violated, , drop = FALSE]
  expect_output(print(sub), "violations: 0")
})

test_that("mb_check reports the equality closure when abatement is observed", {
  data(pigfarms, package = "pgt", envir = environment())
  tech <- pgt_tech(
    x = pigfarms[, c("uncontrolled", "labor", "capital")],
    y = pigfarms$meat, b = pigfarms$controlled,
    u = c(1, 0, 0), a = pigfarms$abatement, id = pigfarms$farm
  )
  mb <- mb_check(tech)
  expect_true(all(c("a", "closure") %in% names(mb)))
  # the pig-farm accounts close exactly: uncontrolled = controlled + a
  expect_equal(mb$closure, rep(0, nrow(mb)), tolerance = 1e-9)
  expect_output(print(mb), "equality closure")
})
