# Generate the synthetic steeldemo panel shipped in data/steeldemo.rda.
# Two routes with different carbon intensities; inputs in CO2-potential
# units (u = 1); emissions satisfy the materials balance
# u'x - v*y >= b with v = 0.01467.

set.seed(20260713)

n_int <- 35  # integrated (BF-BOF style) plants
n_min <- 25  # mini-mill (EAF style) plants
years <- 2021:2023
v <- 0.01467

input_names <- c("coal_coke", "other_fuel", "raw_material", "flux")
share_base <- list(
  Integrated = c(0.75, 0.08, 0.12, 0.05),
  Minimill   = c(0.25, 0.35, 0.30, 0.10)
)

plants <- data.frame(
  plant = sprintf("P%03d", seq_len(n_int + n_min)),
  route = rep(c("Integrated", "Minimill"), c(n_int, n_min)),
  stringsAsFactors = FALSE
)
# Plant-level scale and carbon intensity (tCO2 potential per tonne steel)
plants$scale <- ifelse(plants$route == "Integrated",
                       rlnorm(n_int + n_min, log(3e6), 0.5),
                       rlnorm(n_int + n_min, log(1e6), 0.5))
plants$intensity <- ifelse(plants$route == "Integrated",
                           runif(n_int + n_min, 1.6, 2.4),
                           runif(n_int + n_min, 0.30, 0.65))

steeldemo <- do.call(rbind, lapply(years, function(yr) {
  d <- plants
  d$year <- yr
  d$production <- d$scale * runif(nrow(d), 0.85, 1.10)
  potential <- d$intensity * d$production * runif(nrow(d), 0.98, 1.08)

  # Input shares by route, with noise, rescaled to the potential
  shares <- t(vapply(seq_len(nrow(d)), function(j) {
    s <- share_base[[d$route[j]]] * runif(4, 0.8, 1.2)
    s / sum(s)
  }, numeric(4)))
  X <- shares * potential
  colnames(X) <- input_names

  # Emissions: potential minus carbon retained in product, minus a small
  # non-emitted remainder (process recovery), so the MB identity holds.
  d$emissions <- (potential - v * d$production) * runif(nrow(d), 0.90, 0.995)

  cbind(d[, c("plant", "year", "route", "production")],
        as.data.frame(X),
        d[, "emissions", drop = FALSE])
}))
rownames(steeldemo) <- NULL
steeldemo$route <- factor(steeldemo$route)
steeldemo <- steeldemo[order(steeldemo$plant, steeldemo$year), ]
rownames(steeldemo) <- NULL

# sanity: materials balance holds everywhere
stopifnot(all(rowSums(steeldemo[, input_names]) -
              v * steeldemo$production - steeldemo$emissions >= 0))

save(steeldemo, file = file.path("data", "steeldemo.rda"),
     compress = "bzip2")
cat("steeldemo:", nrow(steeldemo), "rows written to data/steeldemo.rda\n")
