# Generate the synthetic steeldemo panel shipped in data/steeldemo.rda.
# Two routes with different carbon intensities; inputs in CO2-potential
# units (u = 1). The materials-balance identity closes exactly on every
# row:
#     u'x - v * y = emissions + captured,     v = 0.01467,
# where `captured` is the CO2 that leaves the plant as an abatement
# output (post-combustion CCS, capture for utilisation, or slag
# mineralisation) and `capture_energy` is the CO2 potential of the energy
# the capture process itself consumes, a pollution-control input that
# enters the account like any other input.

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

# Abatement technology, assigned at plant level.
#   CCS: post-combustion capture on integrated plants only, adopted in
#        2022 or 2023 so the panel records uptake; captures 60-90 per
#        cent of gross emissions with an energy penalty of 0.10-0.20 t
#        of CO2 potential per tonne captured.
#   CCU: capture for utilisation, either route, 5-15 per cent, same
#        penalty range, all years.
#   mineralisation: slag carbonation, either route, 1-3 per cent, small
#        penalty (0.01-0.03), all years.
plants$tech <- "none"
plants$rate <- 0
plants$penalty <- 0
plants$ccs_year <- NA_integer_

ccs <- sample(which(plants$route == "Integrated"), 8)
plants$tech[ccs] <- "CCS"
plants$rate[ccs] <- runif(8, 0.60, 0.90)
plants$penalty[ccs] <- runif(8, 0.10, 0.20)
plants$ccs_year[ccs] <- rep(c(2022L, 2023L), c(3, 5))

rest <- setdiff(seq_len(nrow(plants)), ccs)
ccu <- sample(rest, 6)
plants$tech[ccu] <- "CCU"
plants$rate[ccu] <- runif(6, 0.05, 0.15)
plants$penalty[ccu] <- runif(6, 0.10, 0.20)

rest <- setdiff(rest, ccu)
mnr <- sample(rest, 10)
plants$tech[mnr] <- "mineralisation"
plants$rate[mnr] <- runif(10, 0.01, 0.03)
plants$penalty[mnr] <- runif(10, 0.01, 0.03)

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

  # Abatement active this year? CCS only from its adoption year on.
  active <- d$tech != "none" & !(d$tech == "CCS" & yr < d$ccs_year)
  r   <- ifelse(active, d$rate, 0)
  pen <- ifelse(active, d$penalty, 0)

  # Gross emissions E0 include the CO2 of the capture energy, which is
  # itself proportional to the captured tonnage:
  #   E0 = (potential - v y) + pen * r * E0.
  net <- potential - v * d$production
  E0 <- net / (1 - pen * r)
  d$captured <- r * E0
  d$capture_energy <- pen * d$captured
  d$emissions <- E0 - d$captured
  d$abatement_tech <- ifelse(active, d$tech, "none")

  cbind(d[, c("plant", "year", "route", "production")],
        as.data.frame(X),
        d[, c("emissions", "captured", "capture_energy",
              "abatement_tech")])
}))
rownames(steeldemo) <- NULL
steeldemo$route <- factor(steeldemo$route)
steeldemo$abatement_tech <- factor(steeldemo$abatement_tech,
  levels = c("none", "CCS", "CCU", "mineralisation"))
steeldemo <- steeldemo[order(steeldemo$plant, steeldemo$year), ]
rownames(steeldemo) <- NULL

# sanity: the identity closes exactly on every row
resid <- rowSums(steeldemo[, input_names]) + steeldemo$capture_energy -
  v * steeldemo$production - steeldemo$emissions - steeldemo$captured
stopifnot(all(abs(resid) < 1e-6 * steeldemo$emissions),
          all(steeldemo$captured >= 0),
          all(steeldemo$emissions > 0))

save(steeldemo, file = file.path("data", "steeldemo.rda"),
     compress = "bzip2")
cat("steeldemo:", nrow(steeldemo), "rows written to data/steeldemo.rda\n")
