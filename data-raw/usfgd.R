# Construction of the usfgd dataset: US coal-fired power plants with
# flue-gas desulphurisation (FGD), 2023, with a derived sulfur-dioxide
# abatement output and the FGD sorbent and electricity as
# pollution-control inputs.
#
# Sources (public, accessed 2026-09-09):
#   EIA-923 (2023 final revision): fuel consumption and net generation
#     (Page 1), coal receipts with sulfur content (Page 5), and
#     Schedule 8C air-emissions-control information (SO2 control units,
#     SO2 removal efficiency, FGD sorbent quantity, FGD electricity).
#     https://www.eia.gov/electricity/data/eia923/archive/xls/f923_2023.zip
#   EPA Clean Air Markets Division (CAMD) Power Sector Emissions Data:
#     unit-level annual SO2 mass and heat input from continuous emission
#     monitoring, primary fuel and SO2 control equipment, through the
#     CAMPD streaming API (annual apportioned emissions); and the
#     facility attributes bulk file for generator nameplate capacity.
#     https://api.epa.gov/easey/streaming-services/emissions/apportioned/annual
#     https://api.epa.gov/easey/bulk-files/facility/facility-2023.csv
#
# Plants are matched on the ORIS plant code (EIA "Plant Id" = CAMD
# "Facility ID"). Sample: plants with at least one CAMD unit whose
# primary fuel is coal and which lists SO2 control equipment, at least
# one operating SO2 control unit in Schedule 8C, positive coal
# consumption, net generation, receipt sulfur content and measured SO2,
# and a positive derived abatement.
#
# Conventions: coal consumption in short tons (Page 1 year total over
# fuel codes ANT, BIT, SUB, LIG, WC, RC, total-fuel basis as in uscoal);
# sulfur is the receipt-tonnage-weighted average percent by weight
# (Page 5, FUEL_GROUP == "Coal"); the SO2 potential of one short ton of
# coal is 2 x sulfur/100 short tons; other_heat is the year-total MMBtu
# of all non-coal fuels; so2 is the CAMD SO2 mass summed over all the
# plant's monitored units; so2_removed is the DERIVED abatement,
# potential minus so2, which also absorbs sulfur retained in ash and
# any error in the sulfur share (it is not a metered quantity);
# efficiency is the mean "SO2 removal efficiency rate at annual
# operating factor" over the plant's operating SO2 control units, which
# EIA instructs respondents to base on monitoring data where available
# and on design specifications otherwise; sorbent is the Schedule 8C
# FGD sorbent quantity (thousand short tons converted to short tons);
# fgd_mwh the Schedule 8C FGD electricity consumption; capacity the sum
# of the nameplate capacities of the generators associated with the
# plant's CAMD units; coal_units and scrubbed_units count the plant's
# monitored coal-fired units and those listing SO2 controls, so a
# partially scrubbed plant (scrubbed_units < coal_units) can be told
# apart from a fully scrubbed one when reading so2_removed.
#
# Run from the package root: Rscript data-raw/usfgd.R
# Requires: readxl. Downloads about 25 MB into tempdir(). The CAMPD API
# accepts the public DEMO_KEY (rate-limited); a free personal key from
# https://api.data.gov/signup can be supplied in the environment
# variable CAMD_API_KEY.

library(readxl)

dir <- file.path(tempdir(), "usfgd-src")
dir.create(dir, showWarnings = FALSE)
ua <- "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"
key <- Sys.getenv("CAMD_API_KEY", "DEMO_KEY")

eia_zip <- file.path(dir, "f923_2023.zip")
if (!file.exists(eia_zip)) {
  for (u in c("https://www.eia.gov/electricity/data/eia923/archive/xls/f923_2023.zip",
              "https://www.eia.gov/electricity/data/eia923/xls/f923_2023.zip")) {
    ok <- tryCatch({
      download.file(u, eia_zip, mode = "wb", headers = c(`User-Agent` = ua))
      TRUE
    }, error = function(e) FALSE)
    if (ok) break
  }
  unzip(eia_zip, exdir = dir)
}
camd_csv <- file.path(dir, "camd_annual_2023.csv")
if (!file.exists(camd_csv)) {
  download.file(
    paste0("https://api.epa.gov/easey/streaming-services/emissions/",
           "apportioned/annual?year=2023&API_KEY=", key),
    camd_csv, mode = "wb",
    headers = c(Accept = "text/csv", `User-Agent` = ua))
}
fac_csv <- file.path(dir, "facility-2023.csv")
if (!file.exists(fac_csv)) {
  download.file("https://api.epa.gov/easey/bulk-files/facility/facility-2023.csv",
                fac_csv, mode = "wb", headers = c(`User-Agent` = ua))
}
main <- file.path(dir,
  "EIA923_Schedules_2_3_4_5_M_12_2023_Final_Revision.xlsx")
env8 <- file.path(dir, "EIA923_Schedule_8_Annual_Envir_Infor_2023_Final.xlsx")

COAL <- c("ANT", "BIT", "SUB", "LIG", "WC", "RC")
num <- function(x) suppressWarnings(as.numeric(x))

# EIA-923 Page 1: plant-year fuel consumption and net generation
p1 <- read_excel(main, sheet = "Page 1 Generation and Fuel Data",
                 skip = 5, guess_max = 20000)
names(p1) <- gsub("[\r\n]+", " ", names(p1))
p1$pid <- suppressWarnings(as.integer(p1[["Plant Id"]]))
p1$fuel <- p1[["Reported Fuel Type Code"]]
p1$qty <- num(p1[["Total Fuel Consumption Quantity"]])
p1$mmbtu <- num(p1[["Total Fuel Consumption MMBtu"]])
p1$gen <- num(p1[["Net Generation (Megawatthours)"]])
p1 <- p1[!is.na(p1$pid), ]
coal_cons <- aggregate(cbind(coal = qty) ~ pid,
                       data = p1[p1$fuel %in% COAL, ], sum)
oth_cons <- aggregate(cbind(other_heat = mmbtu) ~ pid,
                      data = p1[!(p1$fuel %in% COAL), ], sum)
gen_tot <- aggregate(cbind(gen = gen) ~ pid, data = p1, sum)

# EIA-923 Page 5: receipt-tonnage-weighted sulfur content
p5 <- read_excel(main, sheet = "Page 5 Fuel Receipts and Costs",
                 skip = 4, guess_max = 200000)
names(p5) <- gsub("[\r\n]+", " ", names(p5))
p5$pid <- suppressWarnings(as.integer(p5[["Plant Id"]]))
p5$q <- num(p5[["QUANTITY"]])
p5$s <- num(p5[["Average Sulfur Content"]])
p5c <- p5[p5[["FUEL_GROUP"]] == "Coal" & !is.na(p5$q) & !is.na(p5$s) &
            p5$q > 0, ]
sulf <- do.call(rbind, lapply(split(p5c, p5c$pid), function(d)
  data.frame(pid = d$pid[1], sulfur = sum(d$q * d$s) / sum(d$q))))

# EIA-923 Schedule 8C: operating SO2 control units, their removal
# efficiency, sorbent and electricity
c8 <- read_excel(env8, sheet = "8C Air Emissions Control Info",
                 skip = 4, guess_max = 20000)
names(c8) <- gsub("[\r\n]+", " ", names(c8))
c8$pid <- suppressWarnings(as.integer(c8[["Plant ID"]]))
c8$so2id <- c8[["SO2  Control ID"]]
c8$status <- c8[["Status"]]
c8$eff <- num(c8[["SO2 Removal Efficiency Rate  at Annual Operating Factor"]])
c8$sorb <- num(c8[["FGD Sorbent Quantity  (thousand tons)"]])
c8$mwh <- num(c8[["FGD Electricity  Consumption  (Megawatthours)"]])
c8s <- c8[!is.na(c8$pid) & !is.na(c8$so2id) & c8$so2id != "" &
            !is.na(c8$status) & c8$status == "OP", ]
fgd <- do.call(rbind, lapply(split(c8s, c8s$pid), function(d)
  data.frame(pid = d$pid[1],
             n_fgd = length(unique(d$so2id)),
             efficiency = if (any(!is.na(d$eff))) mean(d$eff, na.rm = TRUE)
                          else NA_real_,
             sorbent = 1000 * sum(d$sorb, na.rm = TRUE),
             fgd_mwh = sum(d$mwh, na.rm = TRUE))))

# CAMD annual apportioned emissions: measured SO2 and heat input per
# unit, coal units and their SO2 controls
ca <- read.csv(camd_csv, check.names = FALSE, stringsAsFactors = FALSE)
ca$pid <- suppressWarnings(as.integer(ca[["Facility ID"]]))
ca$so2 <- num(ca[["SO2 Mass (short tons)"]])
ca$hi <- num(ca[["Heat Input (mmBtu)"]])
ca$coal_unit <- !is.na(ca[["Primary Fuel Type"]]) &
  ca[["Primary Fuel Type"]] == "Coal"
ca$scrubbed <- ca$coal_unit & !is.na(ca[["SO2 Controls"]]) &
  ca[["SO2 Controls"]] != ""
ca <- ca[!is.na(ca$pid), ]
camd <- do.call(rbind, lapply(split(ca, ca$pid), function(d)
  data.frame(pid = d$pid[1],
             so2 = sum(d$so2, na.rm = TRUE),
             heat_input = sum(d$hi, na.rm = TRUE),
             n_coal_units = sum(d$coal_unit),
             n_scrubbed = sum(d$scrubbed))))
camd <- camd[camd$n_coal_units > 0 & camd$n_scrubbed > 0, ]

# CAMD facility attributes: name, state, nameplate capacity of the
# generators associated with the plant's units ("1 (153.1)|2 (200)")
fa <- read.csv(fac_csv, check.names = FALSE, stringsAsFactors = FALSE)
fa$pid <- suppressWarnings(as.integer(fa[["Facility ID"]]))
fa <- fa[!is.na(fa$pid), ]
fac <- do.call(rbind, lapply(split(fa, fa$pid), function(d) {
  gens <- unlist(strsplit(d[["Associated Generators & Nameplate Capacity (MWe)"]],
                          "|", fixed = TRUE))
  gens <- trimws(gens[!is.na(gens) & gens != ""])
  gid <- sub(" *\\(.*$", "", gens)
  cap <- num(sub("^.*\\(([^)]*)\\).*$", "\\1", gens))
  keep <- !duplicated(gid) & !is.na(cap)
  data.frame(pid = d$pid[1], name = d[["Facility Name"]][1],
             state = d[["State"]][1], capacity = sum(cap[keep]))
}))

# merge and filter
d <- Reduce(function(a, b) merge(a, b, by = "pid"),
            list(camd, fac, coal_cons, gen_tot, sulf, fgd))
d <- merge(d, oth_cons, by = "pid", all.x = TRUE)
d$other_heat[is.na(d$other_heat)] <- 0
n_matched <- nrow(d)
d <- d[!is.na(d$coal) & d$coal > 0 & !is.na(d$gen) & d$gen > 0 &
         !is.na(d$so2) & d$so2 > 0 & !is.na(d$sulfur) & d$sulfur > 0 &
         !is.na(d$capacity) & d$capacity > 0, ]
n_positive <- nrow(d)
d$potential <- 2 * d$sulfur / 100 * d$coal
d$so2_removed <- d$potential - d$so2
n_violating <- sum(d$so2_removed <= 0)
d <- d[d$so2_removed > 0, ]
d <- d[order(d$pid), ]

usfgd <- data.frame(
  plant = as.character(d$pid),
  name = d$name,
  state = d$state,
  coal = d$coal,
  sulfur = d$sulfur,
  other_heat = d$other_heat,
  capacity = d$capacity,
  gen = d$gen,
  heat_input = d$heat_input,
  coal_units = d$n_coal_units,
  scrubbed_units = d$n_scrubbed,
  so2 = d$so2,
  so2_removed = d$so2_removed,
  efficiency = d$efficiency,
  sorbent = d$sorbent,
  fgd_mwh = d$fgd_mwh,
  n_fgd = d$n_fgd,
  stringsAsFactors = FALSE
)
rownames(usfgd) <- NULL

cat(sprintf("matched %d plants; %d with positive quantities; %d dropped with non-positive derived abatement; %d kept\n",
            n_matched, n_positive, n_violating, nrow(usfgd)))
share <- usfgd$so2_removed / (usfgd$so2_removed + usfgd$so2)
cat(sprintf("derived removal share: median %.3f (IQR %.3f-%.3f); reported efficiency: median %.3f, missing %d; Spearman %.3f\n",
            median(share), quantile(share, .25), quantile(share, .75),
            median(usfgd$efficiency, na.rm = TRUE), sum(is.na(usfgd$efficiency)),
            cor(share, usfgd$efficiency, method = "spearman", use = "complete.obs")))
cat(sprintf("sorbent > 0: %d; fgd_mwh > 0: %d\n", sum(usfgd$sorbent > 0), sum(usfgd$fgd_mwh > 0)))

cat(sprintf("partially scrubbed plants (scrubbed_units < coal_units): %d\n",
            sum(usfgd$scrubbed_units < usfgd$coal_units)))
stopifnot(nrow(usfgd) == 154, n_violating == 1,
          all(usfgd$so2_removed > 0))

save(usfgd, file = "data/usfgd.rda", compress = "xz")
cat("usfgd:", nrow(usfgd), "plants written to data/usfgd.rda\n")
