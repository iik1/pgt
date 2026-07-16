# compare_models(): competing axiom systems on identical data.

test_that("compare_models returns aligned scores and a correlation matrix", {
  tech <- make_random_tech(L = 30, N = 3, seed = 3)
  cmp <- compare_models(tech, models = c("wgd", "byprod", "mb_cost", "wd"))
  expect_s3_class(cmp, "pgt_compare")
  expect_equal(nrow(cmp$scores), tech$L)
  expect_true(all(c("wgd", "byprod", "mb_cost", "wd") %in% names(cmp$scores)))
  expect_equal(dim(cmp$spearman), c(4L, 4L))
  expect_equal(unname(diag(cmp$spearman)), rep(1, 4), tolerance = 1e-9)
  # Symmetric correlation matrix.
  expect_equal(cmp$spearman, t(cmp$spearman), tolerance = 1e-9)
  expect_equal(nrow(cmp$agreement), 4L)
})

test_that("compare_models rejects the directional model and short lists", {
  tech <- make_random_tech(L = 12, seed = 8)
  expect_error(compare_models(tech, models = c("wgd", "fdmo")), "fdmo")
  expect_error(compare_models(tech, models = "wgd"), "at least two")
})

test_that("compare_models scores match standalone pgt fits", {
  tech <- make_random_tech(L = 15, N = 2, seed = 21)
  cmp <- compare_models(tech, models = c("wgd", "mb_cost"))
  wgd <- suppressWarnings(pgt(tech, model = "wgd"))$results$efficiency
  mbc <- pgt(tech, model = "mb_cost")$results$efficiency
  expect_equal(cmp$scores$wgd, wgd, tolerance = 1e-9)
  expect_equal(cmp$scores$mb_cost, mbc, tolerance = 1e-9)
})

test_that("bottom_q_overlap is the model's own worst-quartile share", {
  tech <- make_random_tech(L = 40, N = 3, seed = 12)
  cmp <- compare_models(tech, models = c("wgd", "byprod", "mb_cost", "wd"))
  ag <- cmp$agreement
  # The reference model (first) agrees with itself completely.
  expect_equal(ag$bottom_q_overlap[ag$model == "wgd"], 1)
  # Overlap is a share in [0, 1].
  expect_true(all(ag$bottom_q_overlap >= -1e-9 &
                    ag$bottom_q_overlap <= 1 + 1e-9, na.rm = TRUE))
  # Definition: |bottom_m intersect bottom_wgd| / |bottom_m|.
  s <- cmp$scores
  ref_bottom <- s$wgd <= stats::quantile(s$wgd, 0.25, na.rm = TRUE)
  mb <- s$mb_cost <= stats::quantile(s$mb_cost, 0.25, na.rm = TRUE)
  expect_equal(ag$bottom_q_overlap[ag$model == "mb_cost"],
               sum(mb & ref_bottom, na.rm = TRUE) / sum(mb, na.rm = TRUE))
})

test_that("compare_models print and as.data.frame work", {
  tech <- make_random_tech(L = 20, seed = 4)
  cmp <- compare_models(tech, models = c("wgd", "wd", "mb_cost"))
  expect_output(print(cmp), "model comparison")
  expect_s3_class(as.data.frame(cmp), "data.frame")
})
