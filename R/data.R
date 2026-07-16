#' Pig-finishing farms from Rodseth (2025), Table 1
#'
#' The synthetic manure-transport example printed as Table 1 of Rodseth
#' (2025): five pig-finishing farms with piglets and feed as material
#' inputs, labor and capital as non-material inputs, saleable meat as
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
#'   \item{labor}{Labor input.}
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
#' \eqn{v = 0.01467} (carbon retained in the product). The data are
#' synthetic; they mimic magnitudes, not any real plant.
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
