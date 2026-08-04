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
#' steel output, and Scope 1 CO2 emissions satisfying the
#' materials-balance identity \eqn{u'x - v y \ge b} with
#' \eqn{v = 0.01467} (carbon retained in the product). The synthetic
#' generator (see \code{data-raw/steeldemo.R}) assumes roughly 0.4 per
#' cent retained carbon by mass, converted to CO2 units:
#' \eqn{0.004 \times 44/12 = 0.01467} tonnes of CO2 per tonne of steel;
#' emissions are drawn as the CO2 potential minus this retained content,
#' minus a small non-emitted remainder, so the identity holds in every
#' row. The data are synthetic; they mimic magnitudes, not any real
#' plant.
#'
#' @format A data frame with 180 rows (60 plants over 3 years) and 9
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
#' }
#' @source Simulated; see \code{data-raw/steeldemo.R} in the package
#'   sources.
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
#' \code{v = 0}. Unlike the synthetic \code{steeldemo} panel, the
#' account here is measured on both sides, and it does not close by
#' construction: 8 of the 212 accounts violate the materials-balance
#' condition (measured SO2 exceeds the sulfur-implied potential,
#' a data inconsistency the audit exists to catch), and the remaining
#' accounts are open by the sulfur retained in ash and, at the 180
#' plants with flue-gas desulphurisation (FGD), by the sulfur removed
#' by the scrubbers. The FGD sorbent quantity is a dedicated
#' pollution-control input for \code{x_abate}.
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
#'   \item{sorbent}{FGD sorbent quantity (short tons; 0 without SO2
#'     controls).}
#'   \item{gen}{Net generation (MWh).}
#'   \item{so2}{Measured SO2 emissions (short tons).}
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
