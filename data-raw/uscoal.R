# Construction of the uscoal dataset: US coal-fired power plants, 2022,
# with a measured sulfur-dioxide materials-balance account.
#
# Sources (public, accessed 2026-08-04):
#   EIA-923 (2022 final revision): fuel consumption, receipts with
#     sulfur content, and Schedule 8C air-emissions-control information.
#     https://www.eia.gov/electricity/data/eia923/archive/xls/f923_2022.zip
#   EPA eGRID2022: plant-level measured SO2 emissions, net generation,
#     nameplate capacity, primary fuel category.
#     https://www.epa.gov/system/files/documents/2024-01/egrid2022_data.xlsx
#
# Plants are matched on the ORIS plant code (EIA "Plant Id" = eGRID
# "ORISPL"). Sample: plants with eGRID primary fuel category COAL,
# positive coal consumption, positive net generation, positive measured
# SO2, positive capacity, and at least one 2022 coal receipt with
# reported sulfur content. Plants burning only stockpiled coal (no 2022
# receipts) drop out because their sulfur content is unobserved.
#
# Conventions: coal consumption in short tons (EIA-923 Page 1 year
# total over fuel codes ANT, BIT, SUB, LIG, WC, RC); sulfur content is
# the receipt-tonnage-weighted average percent by weight (Page 5,
# FUEL_GROUP == "Coal"); the SO2 potential of one short ton of coal is
# 2 x sulfur/100 short tons (molar mass ratio SO2/S = 2, full
# conversion, no ash retention); other_heat is the year-total MMBtu of
# all non-coal fuels; sorbent is the Schedule 8C FGD sorbent quantity
# (converted from thousand short tons to short tons; 0 both for plants
# without SO2 controls and for FGD plants that list a control but
# report no sorbent quantity); fgd marks plants with at least one SO2
# control listed in Schedule 8C.
#
# Run from the package root: Rscript data-raw/uscoal.R
# Requires: readxl. Downloads about 38 MB into tempdir().

library(readxl)

dir <- file.path(tempdir(), "uscoal-src")
dir.create(dir, showWarnings = FALSE)
ua <- "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"
eia_zip <- file.path(dir, "f923_2022.zip")
egrid_xlsx <- file.path(dir, "egrid2022_data.xlsx")
if (!file.exists(eia_zip)) {
  download.file(
    "https://www.eia.gov/electricity/data/eia923/archive/xls/f923_2022.zip",
    eia_zip, mode = "wb", headers = c(`User-Agent` = ua))
  unzip(eia_zip, exdir = dir)
}
if (!file.exists(egrid_xlsx)) {
  download.file(
    "https://www.epa.gov/system/files/documents/2024-01/egrid2022_data.xlsx",
    egrid_xlsx, mode = "wb", headers = c(`User-Agent` = ua))
}
main <- file.path(dir,
  "EIA923_Schedules_2_3_4_5_M_12_2022_Final_Revision.xlsx")
env8 <- file.path(dir,
  "EIA923_Schedule_8_Annual_Environmental_Information_2022_Final.xlsx")

COAL <- c("ANT", "BIT", "SUB", "LIG", "WC", "RC")

# EIA-923 Page 1: plant-year fuel consumption
p1 <- read_excel(main, sheet = "Page 1 Generation and Fuel Data",
                 skip = 5, guess_max = 20000)
names(p1) <- gsub("[\r\n]+", " ", names(p1))
p1$pid <- suppressWarnings(as.integer(p1[["Plant Id"]]))
p1$fuel <- p1[["Reported Fuel Type Code"]]
p1$qty <- suppressWarnings(as.numeric(p1[["Total Fuel Consumption Quantity"]]))
p1$mmbtu <- suppressWarnings(as.numeric(p1[["Total Fuel Consumption MMBtu"]]))
p1 <- p1[!is.na(p1$pid), ]
coal_cons <- aggregate(cbind(coal = qty) ~ pid,
                       data = p1[p1$fuel %in% COAL, ], sum)
oth_cons <- aggregate(cbind(other_heat = mmbtu) ~ pid,
                      data = p1[!(p1$fuel %in% COAL), ], sum)

# EIA-923 Page 5: receipt-tonnage-weighted sulfur content
p5 <- read_excel(main, sheet = "Page 5 Fuel Receipts and Costs",
                 skip = 4, guess_max = 200000)
names(p5) <- gsub("[\r\n]+", " ", names(p5))
p5$pid <- suppressWarnings(as.integer(p5[["Plant Id"]]))
p5$q <- suppressWarnings(as.numeric(p5[["QUANTITY"]]))
p5$s <- suppressWarnings(as.numeric(p5[["Average Sulfur Content"]]))
p5c <- p5[p5[["FUEL_GROUP"]] == "Coal" & !is.na(p5$q) & !is.na(p5$s) &
            p5$q > 0, ]
sulf <- do.call(rbind, lapply(split(p5c, p5c$pid), function(d)
  data.frame(pid = d$pid[1], sulfur = sum(d$q * d$s) / sum(d$q))))

# EIA-923 Schedule 8C: SO2 controls and FGD sorbent
c8 <- read_excel(env8, sheet = "8C Air Emissions Control Info",
                 skip = 4, guess_max = 20000)
names(c8) <- gsub("[\r\n]+", " ", names(c8))
c8$pid <- suppressWarnings(as.integer(c8[["Plant ID"]]))
c8$so2id <- c8[["SO2  Control ID"]]
c8$sorb <- suppressWarnings(
  as.numeric(c8[["FGD Sorbent Quantity  (thousand tons)"]]))
c8s <- c8[!is.na(c8$pid) & !is.na(c8$so2id) & c8$so2id != "", ]
fgd <- do.call(rbind, lapply(split(c8s, c8s$pid), function(d)
  data.frame(pid = d$pid[1],
             sorbent = 1000 * sum(d$sorb, na.rm = TRUE))))

# eGRID2022 plant sheet: measured SO2, generation, capacity
pl <- read_excel(egrid_xlsx, sheet = "PLNT22", skip = 1,
                 guess_max = 20000)
pl2 <- data.frame(pid = suppressWarnings(as.integer(pl$ORISPL)),
                  name = pl$PNAME, state = pl$PSTATABB,
                  fuelct = pl$PLFUELCT,
                  capacity = suppressWarnings(as.numeric(pl$NAMEPCAP)),
                  gen = suppressWarnings(as.numeric(pl$PLNGENAN)),
                  so2 = suppressWarnings(as.numeric(pl$PLSO2AN)))
pl2 <- pl2[!is.na(pl2$pid) & !is.na(pl2$fuelct) &
             pl2$fuelct == "COAL", ]

# merge and filter
d <- Reduce(function(a, b) merge(a, b, by = "pid"),
            list(pl2, coal_cons, sulf))
d <- merge(d, oth_cons, by = "pid", all.x = TRUE)
d$other_heat[is.na(d$other_heat)] <- 0
d <- merge(d, fgd, by = "pid", all.x = TRUE)
d$fgd <- !is.na(d$sorbent)
d$sorbent[is.na(d$sorbent)] <- 0
d <- d[!is.na(d$so2) & d$so2 > 0 &
         !is.na(d$gen) & d$gen > 0 &
         !is.na(d$coal) & d$coal > 0 &
         !is.na(d$capacity) & d$capacity > 0, ]
d <- d[order(d$pid), ]

uscoal <- data.frame(
  plant = as.character(d$pid),
  name = d$name,
  state = d$state,
  fgd = d$fgd,
  coal = d$coal,
  other_heat = d$other_heat,
  capacity = d$capacity,
  sorbent = d$sorbent,
  gen = d$gen,
  so2 = d$so2,
  sulfur = d$sulfur,
  stringsAsFactors = FALSE
)
rownames(uscoal) <- NULL

stopifnot(nrow(uscoal) == 212, sum(uscoal$fgd) == 180,
          sum(2 * uscoal$sulfur / 100 * uscoal$coal - uscoal$so2 < 0) == 8)

save(uscoal, file = "data/uscoal.rda", compress = "xz")
cat("uscoal:", nrow(uscoal), "plants written to data/uscoal.rda\n")
