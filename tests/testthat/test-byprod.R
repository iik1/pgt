# Murty, Russell and Levkoff (2012) by-production intersection
# technology. Example 1 (their Eq. 4.1) is a CRS technology with printed,
# analytically solved efficiency scores.

mrl_example1 <- function() {
  # Columns (x, y, z): one input, one good output, one bad output.
  x <- c(1, 1, 1, 2, 2)
  y <- c(2, 3 / 2, 2 / 3, 3, 2)
  z <- c(4, 1, 2, 5, 3)
  pgt_tech(
    x = matrix(x, ncol = 1, dimnames = list(NULL, "x")),
    y = y, b = z, u = 1, polluting = "x",
    id = paste0("DMU", 1:5)
  )
}

test_that("byprod reproduces MRL (2012) Example 1 analytic scores", {
  fit <- pgt(mrl_example1(), model = "byprod", returns = "crs")
  r <- fit$results
  rownames(r) <- r$id

  # DMU3: output efficiency E1 = 1/3, emission efficiency E2 = 1/2,
  # graph measure E_FGL = 5/12.
  expect_equal(r["DMU3", "output_eff"], 1 / 3, tolerance = 1e-7)
  expect_equal(r["DMU3", "efficiency"], 1 / 2, tolerance = 1e-7)
  expect_equal(r["DMU3", "fgl"], 5 / 12, tolerance = 1e-7)

  # DMU2: E1 = 3/4, E2 = 1, E_FGL = 7/8.
  expect_equal(r["DMU2", "output_eff"], 3 / 4, tolerance = 1e-7)
  expect_equal(r["DMU2", "efficiency"], 1, tolerance = 1e-7)
  expect_equal(r["DMU2", "fgl"], 7 / 8, tolerance = 1e-7)
})

test_that("byprod scores are ratios in (0, 1] with FGL their average", {
  fit <- pgt(mrl_example1(), model = "byprod", returns = "crs")
  r <- fit$results
  expect_true(all(r$efficiency > 0 & r$efficiency <= 1 + 1e-9))
  expect_true(all(r$output_eff > 0 & r$output_eff <= 1 + 1e-9))
  expect_equal(r$fgl, (r$output_eff + r$efficiency) / 2, tolerance = 1e-9)
})

test_that("byprod defaults polluting inputs to those with u > 0", {
  x <- cbind(fuel = c(10, 12, 14), labour = c(5, 4, 6))
  tech <- pgt_tech(x, y = c(3, 4, 5), b = c(8, 9, 11),
                   u = c(fuel = 1, labour = 0))
  expect_equal(unname(pgt:::.polluting_inputs(tech)), 1L)
  fit <- pgt(tech, model = "byprod")
  expect_true(all(fit$results$efficiency > 0 &
                    fit$results$efficiency <= 1 + 1e-9))
})
