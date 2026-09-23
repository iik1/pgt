# Rodseth (2025) Table 3: the factorially determined multi-output
# (directional) representation, Eq. 14 (Eq. 13 with abatement fixed at
# each DMU's own level), applied to the Table 1 pig-finishing example.
#
# The evaluated DMU's uncontrolled emission is z_i = u_i'x_i - v y_i, and
# equals the data's `uncontrolled` column by definition. The good-output
# coefficient v = 7/6 follows exactly from farms A and B, which share
# inputs (feed 21, piglet 2, labor 50, capital 100): their uncontrolled
# emissions differ by 3.5 and their meat by 3, so v = 3.5 / 3 = 7/6. Each
# farm's material coefficient on feed is then fixed so that
# u_i * feed_i - v * meat_i equals its uncontrolled emission (the two
# output sets carry different input quality, i.e. DMU-specific u).

pigfarms_fdmo_tech <- function() {
  data(pigfarms, package = "pgt", envir = environment())
  v <- 7 / 6
  # DMU-specific feed coefficient closing each farm's account to its
  # uncontrolled-emission column.
  u_feed <- (pigfarms$uncontrolled + v * pigfarms$meat) / pigfarms$feed
  U <- cbind(feed = u_feed, piglet = 0, labor = 0, capital = 0)
  pgt_tech(
    x = pigfarms[, c("feed", "piglet", "labor", "capital")],
    y = pigfarms$meat,
    b = pigfarms$controlled,
    u = U,
    v = v,
    a = pigfarms$abatement,
    id = pigfarms$farm
  )
}

test_that("fdmo reproduces Rodseth (2025) Table 3 for farms A-D exactly", {
  fit <- pgt(pigfarms_fdmo_tech(), model = "fdmo")
  r <- fit$results
  rownames(r) <- r$id

  # Table 3: gross, good (theta_y), bad (theta_b), maximal good output.
  expect_equal(r["A", "gross"], 0, tolerance = 1e-7)
  expect_equal(r["B", "gross"], 0, tolerance = 1e-7)
  expect_equal(r["D", "gross"], 0, tolerance = 1e-7)

  expect_equal(r["C", "good_eff"], 3, tolerance = 1e-7)
  expect_equal(r["C", "bad_eff"], 3.5, tolerance = 1e-7)
  expect_equal(r["C", "gross"], 6.5, tolerance = 1e-7)
  expect_equal(r["C", "maximal_y"], 10, tolerance = 1e-7)
})

test_that("fdmo bad-output efficiency equals v * good-output efficiency", {
  # Materials-balance identity: with z_i = b_i + a_i (uncontrolled =
  # controlled + abatement) the MB row forces theta_b = v * theta_y for
  # every DMU. This is the exact ratio the paper's rounded Table 3 row
  # for farm E (0.8, 1.0) approximates.
  fit <- pgt(pigfarms_fdmo_tech(), model = "fdmo")
  r <- fit$results
  pos <- r$good_eff > 1e-8
  expect_equal(r$bad_eff[pos] / r$good_eff[pos],
               rep(7 / 6, sum(pos)), tolerance = 1e-6)

  # Farm E: maximal good output rounds to the paper's 10.8.
  e <- r[r$id == "E", ]
  expect_equal(round(e$maximal_y, 1), 10.8)
})

test_that("common coefficients within the paper's rounding reproduce all of Table 3", {
  # The paper prints one decimal, so the farm E divergence on the printed
  # data is conditional. Common coefficients v = 7/6, u_feed = 31/25 and
  # u_piglet = 349/300 with the printed inputs, meat and controlled
  # emissions give uncontrolled emissions and abatement that round to
  # Table 1, and every fitted row, farm E included, rounds to Table 3.
  data(pigfarms, package = "pgt", envir = environment())
  v <- 7 / 6
  u <- c(feed = 31 / 25, piglet = 349 / 300, labor = 0, capital = 0)
  x <- as.matrix(pigfarms[, names(u)])
  z <- as.vector(x %*% u) - v * pigfarms$meat
  a <- z - pigfarms$controlled
  expect_equal(round(z, 1), pigfarms$uncontrolled)
  expect_equal(round(a, 1), pigfarms$abatement)
  r <- pgt(pgt_tech(x = x, y = pigfarms$meat, b = pigfarms$controlled,
                    u = u, v = v, a = a, id = pigfarms$farm),
           model = "fdmo")$results
  expect_equal(round(r$gross, 1), c(0, 0, 6.5, 0, 1.8))
  expect_equal(round(r$good_eff, 1), c(0, 0, 3, 0, 0.8))
  expect_equal(round(r$bad_eff, 1), c(0, 0, 3.5, 0, 1.0))
  expect_equal(round(r$maximal_y, 1), c(10, 7, 10, 11, 10.8))
})

test_that("fdmo requires an abatement output", {
  data(pigfarms, package = "pgt", envir = environment())
  tech <- pgt_tech(
    x = pigfarms[, c("feed", "piglet", "labor", "capital")],
    y = pigfarms$meat, b = pigfarms$controlled, v = 7 / 6,
    id = pigfarms$farm
  )
  expect_error(pgt(tech, model = "fdmo"), "abatement")
})

test_that("fdmo is a frontier at the good-output maximum", {
  # A DMU that attains the maximum good output in its comparison set and
  # sits on its own materials-balance identity has zero gross
  # inefficiency.
  fit <- pgt(pigfarms_fdmo_tech(), model = "fdmo")
  r <- fit$results
  expect_true(all(r$gross >= -1e-8))
  expect_true(all(r$maximal_y >= r$y - 1e-8))
  # ddf is an alias for fdmo.
  fit2 <- pgt(pigfarms_fdmo_tech(), model = "ddf")
  expect_equal(fit$results$gross, fit2$results$gross)
})

test_that("fdmo warns when material accounts do not close exactly", {
  data(pigfarms, package = "pgt", envir = environment())
  v <- 7 / 6
  u_feed <- (pigfarms$uncontrolled + v * pigfarms$meat) / pigfarms$feed
  U <- cbind(feed = u_feed, piglet = 0, labor = 0, capital = 0)
  b_open <- pigfarms$controlled + c(1, rep(0, nrow(pigfarms) - 1L))
  tech <- pgt_tech(
    x = pigfarms[, c("feed", "piglet", "labor", "capital")],
    y = pigfarms$meat, b = b_open, u = U, v = v,
    a = pigfarms$abatement, id = pigfarms$farm
  )
  w <- capture_warnings(pgt(tech, model = "fdmo"))
  expect_true(any(grepl("close the materials-balance identity", w)))
})

test_that("fdmo solves every steeldemo row with capture_energy as control input", {
  data(steeldemo)
  tech <- pgt_tech(
    x = steeldemo[, c("coal_coke", "other_fuel", "raw_material", "flux",
                      "capture_energy")],
    y = steeldemo$production, b = steeldemo$emissions,
    a = steeldemo$captured, v = 0.01467, x_abate = "capture_energy",
    group = steeldemo$route, id = steeldemo$plant
  )
  # the identity closes exactly, so no row may be infeasible; tonne
  # magnitudes are solved on the unit-magnitude copy of the technology
  r <- pgt(tech, model = "fdmo")$results
  expect_true(all(r$status == 0))
  expect_equal(r$bad_eff, 0.01467 * r$good_eff, tolerance = 1e-6)
})

test_that("fdmo keeps projected emissions non-negative under heterogeneous coefficients", {
  # Two closed accounts (u x - v y = b + a) with different input
  # coefficients. Without the admissibility row unit A (u = 0.2) borrows
  # B's output 9 and contracts its emission by 8, to -7; with it the
  # contraction stops at A's own emission: theta_y = theta_b = 1.
  tech <- pgt_tech(x = matrix(c(10, 10), 2, 1), y = c(1, 9), b = c(1, 1),
                   a = c(0, 0), u = matrix(c(0.2, 1), 2, 1), v = c(1, 1),
                   id = c("A", "B"))
  r <- pgt(tech, model = "fdmo")$results
  expect_equal(r$status, c(0L, 0L))
  expect_true(all(r$b - r$bad_eff >= -1e-8))
  expect_equal(r$good_eff, c(1, 0), tolerance = 1e-8)
  expect_equal(r$bad_eff, c(1, 0), tolerance = 1e-8)
  expect_equal(r$gross, c(2, 0), tolerance = 1e-8)
})
