# Rodseth (2025) Table 3: the factorially determined multi-output
# (directional) representation, Eq. 13 with abatement fixed at each
# DMU's own level, applied to the Table 1 pig-finishing example.
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
