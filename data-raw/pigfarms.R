# pigfarms: Table 1 of Rodseth (2025, Journal of Productivity Analysis,
# 64(3), 305-319, doi 10.1007/s11123-025-00768-0): synthetic data on
# manure transport in pig finishing, extending the Rodseth (2016, EJOR)
# reworking of the Coelli, Lauwers and Van Huylenbroeck (2007) example
# with end-of-pipe abatement. Five farms; nitrogen is the material.
#
# The data satisfy uncontrolled - controlled = abatement exactly, so the
# materials-balance closure gap of mb_check() recovers the abatement
# column when the uncontrolled-emission aggregate is used as the
# material carrier (u = 1 on that column, v = 0).

pigfarms <- data.frame(
  farm = c("A", "B", "C", "D", "E"),
  feed = c(21, 21, 21, 24, 24),
  piglet = c(2, 2, 2, 3, 3),
  labor = c(50, 50, 50, 50, 50),
  capital = c(100, 100, 100, 100, 100),
  meat = c(10, 7, 7, 11, 10),
  controlled = c(16, 16, 20, 20, 21),
  uncontrolled = c(16.7, 20.2, 20.2, 20.4, 21.6),
  abatement = c(0.7, 4.2, 0.2, 0.4, 0.6),
  stringsAsFactors = FALSE
)

stopifnot(isTRUE(all.equal(pigfarms$uncontrolled - pigfarms$controlled,
                           pigfarms$abatement)))

save(pigfarms, file = file.path("data", "pigfarms.rda"), compress = "bzip2")
cat("pigfarms:", nrow(pigfarms), "rows written to data/pigfarms.rda\n")
