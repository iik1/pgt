# pgt: Data Envelopment Analysis for Pollution-Generating Technologies

`pgt` implements nonparametric efficiency analysis for
pollution-generating technologies under the materials-balance
principle. It packages the competing axiom systems for modelling bad
outputs behind one interface, with an enforced materials-balance
identity, a pre-estimation feasibility audit, metafrontier
decompositions, bad-output shadow prices, marginal abatement cost
curves, a productivity index and subsampling inference.

Existing R packages cover directional-distance and weak-disposability
baselines (`Benchmarking`, `deaR`, `nonparaeff`, `eat`), but, to our
knowledge, no maintained package in R or any other major statistical
ecosystem ships materials-balance-constrained DEA: applied energy and
steel papers re-code these estimators in GAMS or MATLAB. `pgt` packages
them behind one coherent API.

## Installation

``` r
# development version
remotes::install_github("iik1/pgt")
```

## Quick start

``` r
library(pgt)
data(steeldemo)

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
| `pgt_ml()` | Global Malmquist-Luenberger productivity index (Oh 2010) |
| `boot_pgt()` | Subsampling inference for scores and group means |
| `shadow_prices()`, `mac_curve()` | Constraint duals and marginal abatement cost curves |

DMU-specific material flow coefficients and multiple pollutants are
supported: pass `u` as a matrix or list and `b` as a matrix.

## Validation

Every estimator is checked against the numerical example printed in its
source paper, and these run as unit tests on every check.

- **Rodseth (2025).** The package ships the paper's pig-finishing example
  (`data(pigfarms)`). `pgt(model = "wgd")` reproduces the Table 2 minimal
  controlled emissions (16, 16, 16, 20, 16), and `pgt(model = "fdmo")`
  reproduces the Table 3 directional scores.
- **Murty, Russell and Levkoff (2012).** `pgt(model = "byprod")`
  reproduces the analytic efficiency scores of the paper's Example 1.
- **Coelli, Lauwers and Van Huylenbroeck (2007).** The `mb_cost` model
  uses the paper's phosphorus material-flow coefficients and satisfies
  its `EE = TE x EAE` decomposition.

The vignettes reproduce these in the open:
`vignette("replication", "pgt")`.

## Vignettes

- `introduction`: the core workflow on the steel panel and a rice
  nitrogen-balance example.
- `replication`: the published-result replications above.
- `comparing-axioms`: `compare_models()` across the axiom systems.
- `productivity`: `pgt_ml()` productivity change and `boot_pgt()`
  inference.

## References

- Rodseth, K. L. (2025). On the development of a unified, nonparametric
  materials balance-based efficiency analysis model and its
  applications. *Journal of Productivity Analysis*, 64(3), 305-319.
  doi:10.1007/s11123-025-00768-0
- Coelli, T., Lauwers, L., & Van Huylenbroeck, G. (2007). Environmental
  efficiency measurement and the materials balance condition. *Journal
  of Productivity Analysis*, 28(1-2), 3-12.
- Murty, S., Russell, R. R., & Levkoff, S. B. (2012). On modeling
  pollution-generating technologies. *Journal of Environmental
  Economics and Management*, 64(1), 117-135.
- Kuosmanen, T. (2005). Weak disposability in nonparametric production
  analysis with undesirable outputs. *American Journal of Agricultural
  Economics*, 87(4), 1077-1082.
- Oh, D.-h. (2010). A global Malmquist-Luenberger productivity index.
  *Journal of Productivity Analysis*, 34(3), 183-197.
- Battese, G. E., Rao, D. S. P., & O'Donnell, C. J. (2004). A
  metafrontier production function for estimation of technical
  efficiencies and technology gaps for firms operating under different
  technologies. *Journal of Productivity Analysis*, 21(1), 91-103.
