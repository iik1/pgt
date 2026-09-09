#' Pig-finishing farms from Rodseth (2025), Table 1
#'
#' The synthetic manure-transport example printed as Table 1 of Rodseth
#' (2025): five pig-finishing farms with piglets and feed as material
#' inputs, labour and capital as non-material inputs, saleable meat as
#' the good output, and nitrogen emissions as the bad output, extended
#' with end-of-pipe abatement. The example descends from Coelli, Lauwers
#' and Van Huylenbroeck (2007) via Rodseth (2016). The columns satisfy
#' \code{uncontrolled - controlled = abatement} exactly.
#'
#' The paper does not print the nitrogen coefficients of feed and
#' piglets, so the technology is constructed in pollutant-potential
#' units: the uncontrolled-emission aggregate carries the material with
#' \code{u = 1} on that column and \code{v = 0}, which reproduces the
#' paper's materials-balance cap exactly. With this mapping,
#' \code{pgt(tech, model = "wgd")} reproduces the minimal controlled
#' emissions of the paper's Table 2 (16, 16, 16, 20, 16), and
#' [mb_check()]'s closure gap recovers the abatement column.
#'
#' @format A data frame with 5 rows and 9 columns:
#' \describe{
#'   \item{farm}{Farm identifier, \code{"A"} to \code{"E"}.}
#'   \item{feed}{Feed input.}
#'   \item{piglet}{Piglet input.}
#'   \item{labor}{Labour input.}
#'   \item{capital}{Capital input.}
#'   \item{meat}{Saleable meat (good output).}
#'   \item{controlled}{Controlled nitrogen emissions (bad output).}
#'   \item{uncontrolled}{Uncontrolled (ex ante) nitrogen emissions.}
#'   \item{abatement}{End-of-pipe abatement,
#'     \code{uncontrolled - controlled}.}
#' }
#' @source Rodseth, K. L. (2025). On the development of a unified,
#'   nonparametric materials balance-based efficiency analysis model and
#'   its applications. \emph{Journal of Productivity Analysis}, 64(3),
#'   305--319, Table 1. \doi{10.1007/s11123-025-00768-0}
#' @references
#' Rodseth, K. L. (2016). Environmental efficiency measurement and the
#' materials balance condition reconsidered. \emph{European Journal of
#' Operational Research}, 250(1), 342--346.
#' \doi{10.1016/j.ejor.2015.10.061}
#' @examples
#' data(pigfarms)
#' tech <- pgt_tech(
#'   x = pigfarms[, c("feed", "piglet", "labor", "capital", "uncontrolled")],
#'   y = pigfarms$meat,
#'   b = pigfarms$controlled,
#'   u = c(0, 0, 0, 0, 1),
#'   v = 0,
#'   id = pigfarms$farm
#' )
#' # closure gap recovers the paper's abatement column
#' mb_check(tech)$gap
#' # minimal controlled emissions of Rodseth (2025), Table 2
#' pgt(tech, model = "wgd")$results$b_star
"pigfarms"

#' Synthetic steel plant panel
#'
#' A simulated panel of steel plants loosely calibrated to the structure
#' of global plant-level data: two production routes with different
#' carbon intensities, four material inputs expressed in CO2-potential
#' units (so the material flow coefficients are \eqn{u = 1}), crude
#' steel output, Scope 1 CO2 emissions, and an observed abatement
#' output with its pollution-control input. The materials-balance
#' identity closes exactly in every row,
#' \eqn{u'x - v y = b + a}, with \eqn{v = 0.01467} (carbon retained in
#' the product), \eqn{b} the emissions and \eqn{a} the captured CO2;
#' the five input columns, \code{capture_energy} included, make up
#' \eqn{x}. The synthetic generator (see \code{data-raw/steeldemo.R})
#' assumes roughly 0.4 per cent retained carbon by mass, converted to
#' CO2 units: \eqn{0.004 \times 44/12 = 0.01467} tonnes of CO2 per
#' tonne of steel. Abatement is assigned at plant level: eight
#' integrated plants adopt post-combustion carbon capture and storage
#' (CCS) in 2022 or 2023 and capture 60 to 90 per cent of their gross
#' emissions from the adoption year on, six plants of either route
#' capture 5 to 15 per cent for utilisation (CCU), and ten plants
#' carbonate slag (mineralisation, 1 to 3 per cent). Capturing costs
#' energy: \code{capture_energy} is the CO2 potential of that energy
#' (0.10 to 0.20 tonnes per tonne captured for CCS and CCU, 0.01 to
#' 0.03 for mineralisation), and its CO2 is part of the gross emissions
#' the capture process treats. The data are synthetic; they mimic
#' magnitudes, not any real plant. Versions before 0.6.2 shipped this
#' panel without abatement and with a small unexplained closure gap in
#' every row; the values of all columns changed with the regeneration.
#'
#' @format A data frame with 180 rows (60 plants over 3 years) and 12
#'   columns:
#' \describe{
#'   \item{plant}{Plant identifier.}
#'   \item{year}{Observation year.}
#'   \item{route}{Production route: \code{"Integrated"} or
#'     \code{"Minimill"}.}
#'   \item{production}{Crude steel output (tonnes).}
#'   \item{coal_coke}{Coal and coke input (tonnes CO2 potential).}
#'   \item{other_fuel}{Other fuel input (tonnes CO2 potential).}
#'   \item{raw_material}{Raw material input (tonnes CO2 potential).}
#'   \item{flux}{Flux and alloy input (tonnes CO2 potential).}
#'   \item{emissions}{Scope 1 CO2 emissions (tonnes).}
#'   \item{captured}{CO2 leaving as an abatement output (tonnes): the
#'     \code{a} slot of \code{\link{pgt_tech}}; zero where no
#'     abatement is active.}
#'   \item{capture_energy}{CO2 potential of the energy consumed by the
#'     capture process (tonnes): a pollution-control input, the
#'     natural \code{x_abate} marker.}
#'   \item{abatement_tech}{Abatement active in the plant-year:
#'     \code{"none"}, \code{"CCS"}, \code{"CCU"} or
#'     \code{"mineralisation"}.}
#' }
#' @source Simulated; see \code{data-raw/steeldemo.R} in the package
#'   sources.
#' @examples
#' data(steeldemo)
#' tech <- pgt_tech(
#'   x = steeldemo[, c("coal_coke", "other_fuel", "raw_material", "flux",
#'                     "capture_energy")],
#'   y = steeldemo$production, b = steeldemo$emissions,
#'   a = steeldemo$captured, v = 0.01467, x_abate = "capture_energy",
#'   group = steeldemo$route, id = steeldemo$plant)
#' attr(mb_check(tech), "n_violations")
#' table(steeldemo$abatement_tech, steeldemo$year)
"steeldemo"

#' US coal-fired power plants with a measured SO2 account, 2022
#'
#' A cross-section of 212 US coal-fired power plants for 2022 with a
#' measured sulfur-dioxide materials-balance account, built from two
#' public sources matched on the ORIS plant code: EIA-923 (fuel
#' consumption, coal receipts with sulfur content, and Schedule 8C
#' air-emissions-control information) and EPA eGRID2022 (measured SO2
#' emissions, net generation, nameplate capacity). The SO2 potential of
#' a plant's coal is \code{2 * sulfur/100} short tons per short ton of
#' coal (molar mass ratio, full conversion), giving a producer-specific
#' material flow coefficient; electricity retains no sulfur, so
#' \code{v = 0}. Unlike the synthetic \code{steeldemo} panel, the two
#' sides of the account here come from separate data collections, and
#' it does not close by construction: 8 of the 212 accounts violate the materials-balance
#' condition (measured SO2 exceeds the sulfur-implied potential,
#' a data inconsistency the audit exists to catch), and the remaining
#' accounts are open by the sulfur retained in ash and, at the 180
#' plants with flue-gas desulphurisation (FGD), by the sulfur removed
#' by the scrubbers. The FGD sorbent quantity is a dedicated
#' pollution-control input for \code{x_abate}. Note that
#' \code{sorbent} is 0 both at plants without SO2 controls and at the
#' 36 FGD plants that list a control but report no sorbent quantity,
#' so a zero does not by itself mark an unscrubbed plant; \code{fgd}
#' is the scrubbing indicator. \code{coal} is the EIA-923 total fuel
#' consumption, which at combined-heat-and-power plants includes coal
#' burned for useful thermal output, while \code{gen} is net electric
#' generation only and eGRID allocates a CHP plant's SO2 to
#' electricity with a factor built from EIA heat-input data; at such
#' plants the sulfur-implied potential is overstated relative to the
#' emissions side, and their accounts lean toward looser closure
#' (fewer violations, larger relative gaps).
#'
#' eGRID's SO2 value is a monitored stack measurement only for units
#' reporting to the EPA Clean Air Markets Division's (CAMD) Power
#' Sector Emissions Data, which generally covers fossil units serving
#' generators above 25 MW; for other units eGRID estimates SO2 from
#' EIA heat input and fuel-specific emission factors, that is, from
#' the same EIA-923 fuel data that build the potential (see the
#' eGRID2022 Technical Guide). No filter on this measurement
#' provenance is applied, so the two sides of the account are fully
#' independent only at CAMD-monitored plants.
#'
#' Construction conventions, filters and download URLs are documented
#' in \code{data-raw/uscoal.R} in the package sources (sources accessed
#' 2026-08-04). Plants burning only stockpiled coal drop out because
#' their 2022 receipts carry no sulfur observation.
#'
#' @format A data frame with 212 rows and 11 columns:
#' \describe{
#'   \item{plant}{ORIS plant code (character).}
#'   \item{name}{Plant name.}
#'   \item{state}{State abbreviation.}
#'   \item{fgd}{\code{TRUE} if the plant reports at least one SO2
#'     control in EIA-923 Schedule 8C.}
#'   \item{coal}{Coal consumption (short tons).}
#'   \item{other_heat}{Non-coal fuel consumption (MMBtu).}
#'   \item{capacity}{Nameplate capacity (MW).}
#'   \item{sorbent}{FGD sorbent quantity (short tons; 0 both without
#'     SO2 controls and where an FGD plant reports no quantity).}
#'   \item{gen}{Net generation (MWh).}
#'   \item{so2}{Annual SO2 emissions from eGRID (short tons; see
#'     Details for measurement provenance).}
#'   \item{sulfur}{Receipt-tonnage-weighted sulfur content of the
#'     coal (percent by weight).}
#' }
#' @source US Energy Information Administration, Form EIA-923 (2022
#'   final revision), \url{https://www.eia.gov/electricity/data/eia923/};
#'   US Environmental Protection Agency, eGRID2022,
#'   \url{https://www.epa.gov/egrid}. Both public; matched and filtered
#'   as documented in \code{data-raw/uscoal.R}.
#' @examples
#' data(uscoal)
#' tech <- pgt_tech(
#'   x = uscoal[, c("coal", "other_heat", "capacity", "sorbent")],
#'   y = uscoal$gen,
#'   b = uscoal$so2,
#'   u = cbind(coal = 2 * uscoal$sulfur / 100, other_heat = 0,
#'             capacity = 0, sorbent = 0),
#'   v = 0,
#'   x_abate = "sorbent",
#'   group = factor(ifelse(uscoal$fgd, "FGD", "No FGD")),
#'   id = uscoal$plant
#' )
#' mb <- mb_check(tech)
#' attr(mb, "n_violations")
"uscoal"
