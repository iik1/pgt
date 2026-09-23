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
  s <- lapply(cmp$scores[cmp$models], pgt:::.tie_scores)
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

test_that("compare_models accepts documented aliases and rejects duplicates", {
  tech <- make_random_tech(L = 15, N = 2, seed = 31)
  cmp <- compare_models(tech, models = c("wgd_rodseth", "wd"))
  expect_true("wgd" %in% cmp$models)
  expect_error(compare_models(tech, models = c("ddf", "wd")), "fdmo")
  expect_warning(
    cmp2 <- compare_models(tech, models = c("wgd", "wd", "wgd")),
    "duplicate"
  )
  expect_equal(cmp2$models, c("wgd", "wd"))
  expect_true(all(is.finite(cmp2$spearman)))
})

test_that("print.pgt_compare never reports a self-comparison", {
  # Perfectly correlated rankings put the matrix minimum on the
  # diagonal; the disagreement line must still name two distinct models.
  x <- matrix(c(10, 10, 10, 10), 4, 1)
  tech <- pgt_tech(x, y = c(5, 5, 5, 5), b = c(2, 3, 4, 5))
  cmp <- compare_models(tech, models = c("wgd", "wd"))
  out <- paste(capture.output(print(cmp)), collapse = "\n")
  expect_false(grepl("wgd vs wgd", out))
  expect_false(grepl("wd vs wd", out))
})

test_that("near-tied scores rank as ties, so solver noise cannot reorder them", {
  s <- c(1, 1 - 5e-10, 0.5, 1 - 3e-10, NA, 0.5 + 1e-6)
  expect_equal(pgt:::.tie_scores(s),
               c(rep(1 - 5e-10, 2), 0.5, 1 - 5e-10, NA, 0.5 + 1e-6))
  # two copies differing only by noise below the tolerance rank
  # identically; ranked raw, the noise reorders the efficient block
  set.seed(1)
  base <- c(rep(1, 20), seq(0.3, 0.9, length.out = 20))
  s1 <- base - runif(40, 0, 1e-9) * (base == 1)
  s2 <- base - runif(40, 0, 1e-9) * (base == 1)
  expect_lt(stats::cor(s1, s2, method = "spearman"), 1)
  expect_equal(stats::cor(pgt:::.tie_scores(s1), pgt:::.tie_scores(s2),
                          method = "spearman"), 1)
})
