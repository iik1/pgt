test_that("uscoal ships the documented sample", {
  data(uscoal)
  expect_equal(nrow(uscoal), 212L)
  expect_equal(sum(uscoal$fgd), 180L)
  expect_true(all(uscoal$coal > 0 & uscoal$gen > 0 & uscoal$so2 > 0 &
                    uscoal$capacity > 0 & uscoal$sulfur > 0))
  expect_true(all(uscoal$sorbent >= 0 & uscoal$other_heat >= 0))
  expect_false(anyDuplicated(uscoal$plant) > 0)
})

test_that("uscoal carries a genuinely measured account", {
  data(uscoal)
  tech <- pgt_tech(
    x = uscoal[, c("coal", "other_heat", "capacity", "sorbent")],
    y = uscoal$gen,
    b = uscoal$so2,
    u = cbind(coal = 2 * uscoal$sulfur / 100, other_heat = 0,
              capacity = 0, sorbent = 0),
    v = 0,
    x_abate = "sorbent",
    group = factor(ifelse(uscoal$fgd, "FGD", "No FGD")),
    id = uscoal$plant
  )
  mb <- mb_check(tech)
  # measured accounts violate and stay open in the data, not by design
  expect_equal(attr(mb, "n_negative"), 8L)
  # the faithful programme solves every plant, violators included
  fit <- pgt(tech, model = "wgd")
  expect_true(all(fit$results$status == 0))
  # input-fixed infeasibility is confined to violating plants
  fitx <- suppressWarnings(pgt(tech, model = "wgd_input_fixed"))
  viol <- unique(as.data.frame(mb)$id[as.data.frame(mb)$gap < 0])
  bad <- fitx$results$id[fitx$results$status != 0]
  expect_true(all(bad %in% viol))
  # the five components telescope on real data with sorbent active
  dec <- suppressWarnings(pgt_decompose(tech, type = "rodseth"))
  expect_equal(dec$results$total, fit$results$efficiency,
               tolerance = 1e-6)
  expect_true(any(dec$results$ae_abatement < 1 - 1e-8, na.rm = TRUE))
  expect_true(any(dec$results$quality < 1 - 1e-8, na.rm = TRUE))
})
