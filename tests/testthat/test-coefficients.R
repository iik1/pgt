# DMU-specific material flow coefficients (Eder 2022, Rodseth 2025) and
# multi-pollutant technologies.

test_that("constant-matrix u matches vector u (back-compatibility)", {
  tech_vec <- make_random_tech(L = 25, N = 3, seed = 7)
  U <- matrix(tech_vec$u[1, , 1], tech_vec$L, tech_vec$N, byrow = TRUE)
  tech_mat <- pgt_tech(
    x = tech_vec$x, y = tech_vec$y, b = tech_vec$b[, 1],
    u = U, v = tech_vec$v[1, 1], group = tech_vec$group
  )
  expect_equal(pgt(tech_vec, model = "wgd")$results$b_star,
               pgt(tech_mat, model = "wgd")$results$b_star)
  expect_equal(mb_check(tech_vec)$gap, mb_check(tech_mat)$gap)
})

test_that("DMU-specific u enters the materials-balance potential", {
  x <- matrix(c(10, 10, 10, 20, 20, 20), 3, 2)
  y <- c(5, 5, 5)
  b <- c(2, 2, 2)
  # First DMU carries double the pollutant content of the others.
  U <- rbind(c(2, 2), c(1, 1), c(1, 1))
  tech <- pgt_tech(x, y, b, u = U, v = 0)
  mb <- mb_check(tech)
  expect_equal(mb$potential, c(2 * 30, 30, 30))
  # Heterogeneous coefficients are reported by print.
  expect_output(print(tech), "DMU-specific")
})

test_that("heterogeneous u is not equivalent to its column means", {
  # Columns differ within each DMU, so a DMU-specific coefficient matrix
  # is not reproduced by applying the column means to every DMU.
  x <- matrix(c(10, 40, 40, 10), 2, 2)
  y <- c(3, 3)
  b <- c(1, 1)
  U_het <- rbind(c(2, 0.5), c(0.5, 2))
  U_avg <- matrix(c(1.25, 1.25), 2, 2, byrow = TRUE)
  expect_false(isTRUE(all.equal(
    mb_check(pgt_tech(x, y, b, u = U_het))$potential,
    mb_check(pgt_tech(x, y, b, u = U_avg))$potential
  )))
})

test_that("single pollutant stored in canonical matrix/array shape", {
  tech <- pgt_tech(matrix(1:6, 3, 2), c(1, 2, 3), c(1, 1, 1))
  expect_equal(dim(tech$b), c(3L, 1L))
  expect_equal(dim(tech$v), c(3L, 1L))
  expect_equal(dim(tech$u), c(3L, 2L, 1L))
  expect_equal(tech$pollutants, "b")
})

test_that("two pollutants each carry a materials-balance identity", {
  x <- matrix(c(10, 12, 14, 8, 9, 10), 3, 2,
              dimnames = list(NULL, c("coal", "gas")))
  y <- c(4, 5, 6)
  B <- cbind(co2 = c(20, 24, 28), so2 = c(1.0, 1.2, 1.4))
  U <- list(
    co2 = c(coal = 2, gas = 1),
    so2 = c(coal = 0.1, gas = 0)
  )
  tech <- pgt_tech(x, y, B, u = U, v = c(co2 = 0.5, so2 = 0))
  expect_equal(tech$P, 2L)
  expect_equal(tech$pollutants, c("co2", "so2"))
  mb <- mb_check(tech)
  expect_equal(nrow(mb), 6L)
  expect_true(all(c("pollutant") %in% names(mb)))

  # CO2 potential of DMU 1: 2*10 + 1*8 = 28; SO2 potential: 0.1*10 = 1.
  co2 <- mb[mb$pollutant == "co2", ]
  so2 <- mb[mb$pollutant == "so2", ]
  expect_equal(co2$potential[1], 28)
  expect_equal(so2$potential[1], 1)
})

test_that("pollutant selector drives the estimated bad output", {
  x <- matrix(c(10, 12, 14, 8, 9, 10), 3, 2)
  y <- c(4, 5, 6)
  B <- cbind(co2 = c(20, 24, 28), so2 = c(3, 2, 4))
  tech <- pgt_tech(x, y, B, u = list(co2 = c(2, 1), so2 = c(1, 1)),
                   v = c(0, 0))
  fit_co2 <- pgt(tech, model = "envelope", pollutant = "co2")
  fit_so2 <- pgt(tech, model = "envelope", pollutant = 2L)
  expect_equal(fit_co2$results$b, B[, "co2"])
  expect_equal(fit_so2$results$b, B[, "so2"])
  expect_equal(fit_so2$pollutant, "so2")
  expect_error(pgt(tech, pollutant = "nox"), "not found")
})

test_that("abatement column count must match pollutants", {
  x <- matrix(c(10, 12, 14, 8, 9, 10), 3, 2)
  y <- c(4, 5, 6)
  B <- cbind(co2 = c(20, 24, 28), so2 = c(3, 2, 4))
  expect_error(
    pgt_tech(x, y, B, a = c(1, 1, 1)),
    "one column per pollutant"
  )
})

test_that("named abatement columns align to pollutants by name, not position", {
  x <- matrix(c(10, 12, 14, 8, 9, 10), 3, 2,
              dimnames = list(NULL, c("coal", "gas")))
  y <- c(4, 5, 6)
  B <- cbind(co2 = c(20, 24, 28), so2 = c(3, 2, 4))
  A_ord <- cbind(co2 = c(1, 2, 3), so2 = c(0.1, 0.2, 0.3))
  A_rev <- A_ord[, c("so2", "co2")]        # user supplies reversed order
  tech_ord <- pgt_tech(x, y, B, a = A_ord, v = c(0, 0))
  tech_rev <- pgt_tech(x, y, B, a = A_rev, v = c(0, 0))
  # both must store abatement aligned to the pollutant order of b
  expect_equal(tech_ord$a, tech_rev$a)
  expect_equal(tech_rev$a[, "co2"], c(1, 2, 3))
  expect_equal(tech_rev$a[, "so2"], c(0.1, 0.2, 0.3))
  expect_error(
    pgt_tech(x, y, B, a = cbind(nox = c(1, 2, 3), so2 = c(1, 1, 1))),
    "column names must match"
  )
})

test_that("fractional pollutant index is rejected", {
  x <- matrix(c(10, 12, 14, 8, 9, 10), 3, 2)
  B <- cbind(co2 = c(20, 24, 28), so2 = c(3, 2, 4))
  tech <- pgt_tech(x, y = c(4, 5, 6), b = B, u = list(co2 = c(1, 1), so2 = c(1, 1)))
  expect_error(pgt(tech, pollutant = 1.9), "integer index")
  expect_silent(pgt(tech, model = "envelope", pollutant = 2))  # integer double ok
})
