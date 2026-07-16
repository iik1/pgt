#' pgt: Pollution-generating technologies under the materials balance
#'
#' Nonparametric efficiency analysis of pollution-generating
#' technologies built on the materials-balance condition
#' \eqn{u'x - v y \ge b}: the pollutant bound in a DMU's material
#' inputs (\eqn{u'x}, with material flow coefficients \eqn{u}) is
#' either retained in the good output (\eqn{v y}) or leaves as the bad
#' output \eqn{b}. How bad outputs should enter the technology is
#' contested, and the package implements the competing axiom systems
#' (weak-G-disposability, by-production, materials-balance cost, weak
#' disposability and the directional representation) behind one
#' interface, so results can be compared across systems on identical
#' data.
#'
#' The workflow: construct the technology with [pgt_tech()], audit the
#' materials-balance accounts with [mb_check()], estimate with [pgt()],
#' and post-process with [pgt_decompose()] (metafrontier
#' decompositions), [shadow_prices()] and [mac_curve()] (output duals
#' and marginal abatement costs) and [compare_models()] (cross-system
#' agreement). [pgt_ml()] computes global Malmquist-Luenberger
#' productivity change and [boot_pgt()] bootstrap inference.
#'
#' Scope: each technology carries one good output (several pollutants
#' and DMU-specific material flow coefficients are supported). The
#' efficiency estimators treat the rows as an independent
#' cross-section; panel structure enters only through [pgt_ml()]'s
#' pooled global frontier, and [boot_pgt()] warns on panel
#' technologies.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom lpSolveAPI make.lp lp.control set.objfn add.constraint
#'   set.bounds get.variables get.objective get.sensitivity.rhs
#' @importFrom stats median quantile aggregate complete.cases
#' @importFrom utils head
#' @importFrom graphics hist boxplot barplot lines legend axis abline mtext
NULL
