test_that("usfgd closes its account and solves every model", {
  data(usfgd)
  expect_equal(nrow(usfgd), 154L)
  expect_true(all(usfgd$so2_removed > 0))
  pot <- 2 * usfgd$sulfur / 100 * usfgd$coal
  expect_equal(pot, usfgd$so2 + usfgd$so2_removed, tolerance = 1e-10)
  tech <- pgt_tech(
    x = usfgd[, c("coal", "other_heat", "capacity", "sorbent", "fgd_mwh")],
    y = usfgd$gen, b = usfgd$so2, a = usfgd$so2_removed,
    u = cbind(coal = 2 * usfgd$sulfur / 100, other_heat = 0,
              capacity = 0, sorbent = 0, fgd_mwh = 0),
    v = 0, x_abate = c("sorbent", "fgd_mwh"), id = usfgd$plant
  )
  expect_equal(attr(mb_check(tech), "n_violations"), 0L)
  for (m in c("wgd", "wgd_input_fixed", "fdmo")) {
    r <- pgt(tech, model = m)$results
    expect_true(all(r$status == 0), info = m)
  }
  dec <- pgt_decompose(tech, type = "rodseth")$results
  expect_true(all(stats::complete.cases(
    dec[c("te_production", "quality", "ae_production",
          "te_abatement", "ae_abatement", "total")])))
})
