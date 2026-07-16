# Global Malmquist-Luenberger index (Oh 2010) under the pgt technology.

two_period_panel <- function(shrink = 0.85, grow = 1.05, seed = 1, n = 12) {
  set.seed(seed)
  d1 <- data.frame(id = seq_len(n), y = runif(n, 5, 10),
                   x = runif(n, 8, 12), b = runif(n, 2, 6))
  d2 <- data.frame(id = seq_len(n), y = d1$y * grow, x = d1$x,
                   b = d1$b * shrink)
  d <- rbind(cbind(d1, period = 1), cbind(d2, period = 2))
  pgt_tech(x = d[, "x", drop = FALSE], y = d$y, b = d$b,
           period = d$period, id = d$id)
}

test_that("ML index decomposition GML = EC x BPC holds exactly", {
  ml <- pgt_ml(two_period_panel())
  r <- ml$results
  expect_equal(r$gml, r$ec * r$bpc, tolerance = 1e-9)
  expect_true(all(is.finite(r$gml)))         # global frontier: no infeasibility
})

test_that("a uniform frontier improvement is recovered as growth via BPC", {
  # Period 2 emits 15% less for the same output: a pure frontier shift.
  ml <- pgt_ml(two_period_panel(shrink = 0.85))
  r <- ml$results
  expect_true(median(r$gml) > 1)             # productivity growth
  expect_true(median(r$bpc) > 1)             # driven by best-practice change
  expect_equal(median(r$ec), 1, tolerance = 1e-6)  # each DMU keeps its rank
})

test_that("no productivity change when the frontier is unchanged", {
  # Period 2 identical to period 1: no output growth, no emission change.
  ml <- pgt_ml(two_period_panel(shrink = 1, grow = 1))
  r <- ml$results
  expect_equal(r$gml, rep(1, nrow(r)), tolerance = 1e-6)
})

test_that("ML directional distance matches an independent recompute", {
  tech <- two_period_panel(seed = 3)
  ml <- pgt_ml(tech, technology = "envelope")
  # Independent global directional distance for one DMU-period, envelope
  # (inputs free), direction g = (y, b).
  ddf_ref <- function(k) {
    L <- tech$L
    lp <- lpSolveAPI::make.lp(0, L + 1)
    invisible(lpSolveAPI::lp.control(lp, sense = "max"))
    o <- numeric(L + 1); o[L + 1] <- 1; lpSolveAPI::set.objfn(lp, o)
    r <- c(tech$y, -tech$y[k]); lpSolveAPI::add.constraint(lp, r, ">=", tech$y[k])
    r <- c(tech$b[, 1], tech$b[k, 1]); lpSolveAPI::add.constraint(lp, r, "<=", tech$b[k, 1])
    lpSolveAPI::add.constraint(lp, c(rep(1, L), 0), "=", 1)
    lpSolveAPI::set.bounds(lp, lower = c(rep(0, L), -Inf))
    lpSolveAPI::solve.lpExtPtr(lp)
    lpSolveAPI::get.variables(lp)[L + 1]
  }
  # Reconstruct DMU 1's GML from independent global distances.
  k1 <- which(tech$id == 1 & tech$period == 1)
  k2 <- which(tech$id == 1 & tech$period == 2)
  gml_ref <- (1 + ddf_ref(k1)) / (1 + ddf_ref(k2))
  got <- r <- ml$results[ml$results$id == "1", "gml"]
  expect_equal(got, gml_ref, tolerance = 1e-7)
})

test_that("pgt_ml requires a period", {
  data(steeldemo)
  tech <- pgt_tech(x = steeldemo[, c("coal_coke", "other_fuel")],
                   y = steeldemo$production, b = steeldemo$emissions,
                   id = steeldemo$plant)
  expect_error(pgt_ml(tech), "period")
})

test_that("summary.pgt_ml tolerates NA transitions", {
  ml <- pgt_ml(two_period_panel())
  ml$results$gml[1] <- NA_real_
  s <- summary(ml)
  expect_s3_class(s, "summary.pgt_ml")
  expect_true(all(is.finite(s$quantiles["EC", ])))
  expect_output(print(s), "Distribution")
  expect_output(print(ml), "transitions with NA index")
})

test_that("technology = 'wd' gives the weak-disposability GML, nested in free disposal", {
  tech <- two_period_panel(seed = 5)
  ml_wd <- pgt_ml(tech, technology = "wd")
  r <- ml_wd$results
  expect_equal(r$gml, r$ec * r$bpc, tolerance = 1e-9)
  expect_true(all(is.finite(r$gml)))
  # Under CRS the free-disposal technology contains the
  # weak-disposability set (take lambda = z), so directional distances
  # are weakly larger there, observation by observation. Under VRS the
  # sets are not nested (the Kuosmanen split lets active intensity sum
  # below one), so the comparison is only made under CRS.
  X <- tech$x; y <- tech$y; b <- tech$b[, 1]
  gl <- seq_len(tech$L)
  d_wd <- vapply(gl, function(k)
    pgt:::.ddf_ml_wd(k, gl, X, y, b, vrs = FALSE), numeric(1))
  d_free <- vapply(gl, function(k)
    pgt:::.ddf_ml(k, gl, X, y, b, use_inputs = TRUE, vrs = FALSE),
    numeric(1))
  expect_true(all(d_wd <= d_free + 1e-8))
})
