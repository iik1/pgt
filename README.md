# pgt: Data Envelopment Analysis for Pollution-Generating Technologies

`pgt` implements nonparametric efficiency analysis for
pollution-generating technologies under the materials-balance
principle. It packages the competing axiom systems for modelling bad
outputs behind one interface, with an enforced materials-balance
identity, a pre-estimation feasibility audit, metafrontier
decompositions, bad-output shadow prices, marginal abatement cost
curves, a productivity index and subsampling inference.

The materials-balance principle goes back to Ayres and Kneese (1969);
Lauwers (2009) makes the case for building it into frontier models, and
Dakpo, Jeanneaux and Latruffe (2016) survey the modelling landscape.
Existing R packages handle undesirable outputs by data translation
(`deaR`'s Seiford-Zhu undesirable-output models) or by weak
disposability and directional distances (`Benchmarking::dea.direct`,
`nonparaeff`'s directional-distance routines), but, to our knowledge,
no maintained package in R or any other major statistical ecosystem
ships materials-balance-constrained DEA; applied studies implement
these estimators in ad-hoc optimisation code. `pgt` packages them
behind one coherent API.

## Scope

Each technology has one good output; several pollutants (with their own
materials-balance accounts) and DMU-specific material flow coefficients
are supported. The efficiency estimators treat the rows as an
independent cross-section: panel structure enters only through
`pgt_ml()`'s pooled global frontier, and `boot_pgt()` warns when a
panel technology is subsampled.

## Installation

``` r
# development version
remotes::install_github("iik1/pgt")
# the peer-reviewed snapshot
remotes::install_github("iik1/pgt@v0.4.1")
```

Cite the package with `citation("pgt")`.

## Quick start

``` r
library(pgt)
data(steeldemo)   # synthetic steel-plant panel shipped with the package

# 1. Build the technology: inputs in CO2-potential units (u = 1),
#    v = carbon retained in the product.
tech <- pgt_tech(
  x = steeldemo[, c("coal_coke", "other_fuel", "raw_material", "flux")],
  y = steeldemo$production,
  b = steeldemo$emissions,
  v = 0.01467,
  group = steeldemo$route,
  id = steeldemo$plant
)

# 2. Audit the materials-balance identity before estimating.
mb_check(tech)

# 3. Fit the weak-G-disposability model (Rodseth 2025, Eq. 6).
fit <- pgt(tech, model = "wgd")
summary(fit)

# 4. Decompose environmental efficiency across production routes.
summary(pgt_decompose(tech, type = "envelope"))

# 5. Shadow prices and the marginal abatement cost curve.
head(shadow_prices(fit))
plot(mac_curve(fit, price = 550))

# 6. Compare competing axiom systems on the same data.
compare_models(tech, models = c("wgd", "byprod", "mb_cost", "wd"))
```

## Models

| Function | What it does |
|---|---|
| `pgt_tech()` | Technology constructor: inputs, good/bad outputs, material flow coefficients `u`, `v`, abatement `a`, technology groups, panel `period` |
| `mb_check()` | Audit of `u'x - v y >= b` per DMU and pollutant |
| `pgt(model = "wgd")` | Rodseth (2025) weak-G-disposability LP |
| `pgt(model = "envelope")` | Eq. 6 with inputs free: convex lower (y, b) envelope |
| `pgt(model = "fdmo")` | Rodseth (2025) directional representation, Eq. 13 (alias `"ddf"`) |
| `pgt(model = "mb_cost")` | Coelli et al. (2007) materials-balance cost model, `EE = TE x EAE` |
| `pgt(model = "byprod")` | Murty-Russell-Levkoff (2012) by-production intersection technology |
| `pgt(model = "wd")` | Kuosmanen (2005) weak-disposability reference model |
| `pgt_decompose()` | Metafrontier decompositions (envelope WR x TGR; staged Rodseth) |
| `compare_models()` | Competing axiom systems on identical data, with rank agreement |
| `pgt_ml()` | Global Malmquist-Luenberger productivity index (Oh 2010; experimental) |
| `boot_pgt()` | Subsampling inference for scores and group means (heuristic intervals) |
| `shadow_prices()`, `mac_curve()` | Constraint duals and marginal abatement cost curves |

DMU-specific material flow coefficients and multiple pollutants are
supported: pass `u` as a matrix or list and `b` as a matrix.

## Validation

Validation runs as unit tests on every check, at three levels of
strength.

- **Published-table replication.** `pgt(model = "wgd")` reproduces the
  Rødseth (2025) Table 2 minimal controlled emissions
  (16, 16, 16, 20, 16), `pgt(model = "fdmo")` the Table 3 directional
  scores for farms A to D exactly (farm E is discussed in the
  replication vignette), and `pgt(model = "byprod")` the analytic
  efficiency scores of Murty, Russell and Levkoff's (2012) Example 1.
  The package ships the pig-finishing example as `data(pigfarms)`.
- **Analytic hand-computed checks.** The `mb_cost`, `wd` and
  `envelope` models are verified against small problems solved by
  hand, and the `wgd` kernel against an independent reference
  implementation in the test suite. The `mb_cost` check uses the
  phosphorus material-flow coefficients of Coelli, Lauwers and Van
  Huylenbroeck (2007) and confirms the `EE = TE x EAE` decomposition
  (an internal-consistency check, not a replication of a printed
  table).
- **Identity and property checks.** The decompositions' multiplicative
  identities, `GML = EC x BPC`, score ranges and infeasibility
  semantics are asserted across randomised technologies.

The vignettes reproduce the published-table replications in the open:
`vignette("replication", "pgt")`.

## Vignettes

- `introduction`: the core workflow on the synthetic steel panel and a
  rice nitrogen-balance example.
- `models`: the estimating linear programs, stated in full.
- `replication`: the published-result replications above.
- `comparing-axioms`: `compare_models()` across the axiom systems.
- `multiple-pollutants`: multi-pollutant technologies, DMU-specific
  coefficients and the staged decomposition.
- `productivity`: `pgt_ml()` productivity change and `boot_pgt()`
  inference.

## References

- Ayres, R. U., & Kneese, A. V. (1969). Production, consumption, and
  externalities. *American Economic Review*, 59(3), 282-297.
- Battese, G. E., Rao, D. S. P., & O'Donnell, C. J. (2004). A
  metafrontier production function for estimation of technical
  efficiencies and technology gaps for firms operating under different
  technologies. *Journal of Productivity Analysis*, 21(1), 91-103.
  doi:10.1023/B:PROD.0000012454.06094.29
- Coelli, T., Lauwers, L., & Van Huylenbroeck, G. (2007). Environmental
  efficiency measurement and the materials balance condition. *Journal
  of Productivity Analysis*, 28(1-2), 3-12. doi:10.1007/s11123-007-0052-8
- Dakpo, K. H., Jeanneaux, P., & Latruffe, L. (2016). Modelling
  pollution-generating technologies in performance benchmarking: Recent
  developments, limits and future prospects in the nonparametric
  framework. *European Journal of Operational Research*, 250(2),
  347-359. doi:10.1016/j.ejor.2015.07.024
- Kuosmanen, T. (2005). Weak disposability in nonparametric production
  analysis with undesirable outputs. *American Journal of Agricultural
  Economics*, 87(4), 1077-1082. doi:10.1111/j.1467-8276.2005.00788.x
- Lauwers, L. (2009). Justifying the incorporation of the materials
  balance principle into frontier-based eco-efficiency models.
  *Ecological Economics*, 68(6), 1605-1614.
  doi:10.1016/j.ecolecon.2008.08.022
- Murty, S., Russell, R. R., & Levkoff, S. B. (2012). On modeling
  pollution-generating technologies. *Journal of Environmental
  Economics and Management*, 64(1), 117-135.
  doi:10.1016/j.jeem.2012.02.005
- Oh, D.-h. (2010). A global Malmquist-Luenberger productivity index.
  *Journal of Productivity Analysis*, 34(3), 183-197.
  doi:10.1007/s11123-010-0178-y
- Rødseth, K. L. (2025). On the development of a unified, nonparametric
  materials balance-based efficiency analysis model and its
  applications. *Journal of Productivity Analysis*, 64(3), 305-319.
  doi:10.1007/s11123-025-00768-0
