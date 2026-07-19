# pgt 0.5.1

* Terminology: the score reported in the `efficiency` column is called
  the principal environmental-efficiency score throughout the
  documentation and printed output; the internal helper was renamed
  accordingly.

# pgt 0.5.0

Revisions from a full peer review of the package as a software
submission.

## Inference

* `boot_pgt()` no longer forces the evaluated DMU into every subsample:
  each unit is scored against the subsample frontier alone, so frontier
  units carry genuine resampling variation instead of a degenerate
  zero-width interval. Replicates with an infeasible programme surface
  as `NA`; the new `per_dmu$n_ok` column counts feasible replicates and
  a warning fires when they run thin. Interval endpoints truncate to
  `[0, max(1, estimate)]`, respecting above-1 violator scores.
* The default `kappa` for `model = "byprod"` now uses the
  polluting-input dimension of the T2 sub-programme; `boot_pgt()` warns
  when a panel technology is subsampled as a cross-section; group-mean
  intervals are documented as a descriptive stability band (the
  aggregate theory of Kneip, Simar and Wilson 2015 is not implemented).
* `boot_pgt_sensitivity()` reports the implied rate exponent
  `kappa_hat` from the Politis-Romano-Wolf log-spread regression, so
  the borrowed DEA rate is checkable, plus the median over strictly
  positive widths.
* A Monte Carlo coverage study under a known envelope-model DGP ships
  as `inst/simulations/coverage.R`; its results are reported in the
  productivity vignette.

## Productivity index

* `pgt_ml()` gains `technology = "wd"`: the weak-disposability
  directional distance in the Kuosmanen form, the technology under
  which Chung, Fare and Grosskopf (1997) and Oh (2010) define the
  index, making it the faithful Oh (2010) comparator. The free-disposal
  semantics of the default `"wgd"` and `"envelope"` options are now
  stated plainly and the claim that they coincide with the
  weak-disposability index when the cap does not bind has been removed
  (it was not what the kernel computed).

## Documentation and validation

* The validation claims are scoped honestly: published-table
  replications (wgd, fdmo, byprod) are distinguished from analytic
  hand-computed checks (mb_cost, wd, envelope) and property checks; the
  replication vignette shows the farm E divergence from Rodseth (2025)
  Table 3 side by side and explains why the printed row cannot satisfy
  the materials-balance-forced ratio; the two documented pigfarms
  mappings are stated and tested to agree.
* Two new vignettes: `models` states every estimating programme in
  full, and `multiple-pollutants` demonstrates multi-pollutant
  technologies, DMU-specific coefficients and the staged decomposition.
* `compare_models()` documents that the models' principal scores are
  different
  quantities (EE for mb_cost) and that only rank-based statistics are
  strictly comparable; `mac_curve()` documents which margin the curve
  prices; foundational and origin literature is cited (Ayres and
  Kneese 1969; Lauwers 2009; Dakpo, Jeanneaux and Latruffe 2016; Chung,
  Fare and Grosskopf 1997; Chambers, Chung and Fare 1996; Kuosmanen and
  Podinovski 2009; O'Donnell, Rao and Battese 2008; Hampf 2014; Fare,
  Grosskopf, Lovell and Yaisawarng 1993), every vignette carries a
  reference list, and DOIs are attached throughout.
* `inst/CITATION` added; a package overview with a scope statement
  replaces the empty `pgt-package` help page; the README states the
  package's scope and labels the experimental and heuristic modules.

# pgt 0.4.1

Fixes from a whole-package review.

## Inference

* `boot_pgt()` now builds proper subsampling intervals (Politis and
  Romano 1994): the subsample distribution is recentred at the point
  estimate and rescaled by the relative convergence rate `(m/L)^kappa`,
  so intervals extend downward from the estimate, matching the
  direction of the frontier bias. Previous versions built percentile
  intervals from the raw subsample distribution, which lies entirely at
  or above the point estimate. The reported `se` is rescaled the same
  way. A new `kappa` argument controls the rate exponent, defaulting to
  the DEA rate for the model's effective dimension.
* `boot_pgt()` restores the caller's RNG state on exit instead of
  permanently reseeding the session.

## Model correctness

* Multi-pollutant `pgt(model = "wgd")` now enforces every pollutant's
  materials-balance cap on the peer mix, as documented; earlier
  versions constrained only the selected pollutant.
* `pgt(model = "fdmo")` warns when material accounts do not close
  exactly (`u'x - v y != b + a`): the model imposes the identity as an
  equality, so open accounts either make the LP infeasible or shift the
  closure gap into the scores. `mb_check()` reports the equality
  residual in a new `closure` column when abatement is observed.
* A failed sub-LP in `pgt(model = "byprod")` or `"mb_cost"` now sets
  every score of that DMU to `NA`, so `status != 0` always means the
  scores are `NA`, matching the warning text and `boot_pgt()`.
* `pgt(model = "mb_cost")` warns when the implied minimal emission
  `b_star` is negative (possible under DMU-specific coefficients).
* `pgt()` and `pgt_decompose()` warn that `x_abate` is recorded but not
  yet used by any estimator.
* `pgt_tech()` rejects an ambiguous `v` whose length equals both the
  number of DMUs and the number of pollutants, and rejects duplicated
  `(id, period)` pairs, which `pgt_ml()` silently mismatched before.

## Interface

* `compare_models()` recognises the documented aliases: `wgd_rodseth`
  is accepted, and `ddf` now yields the clear directional-model
  rejection instead of a misleading unknown-model error. Duplicated
  entries are ignored with a warning, and the largest-disagreement line
  no longer reports a model's self-correlation.
* `mac_curve()` and `shadow_prices()` now error clearly for models
  without output duals (`fdmo`, `byprod`, `mb_cost`, `wd`); earlier
  versions returned an empty curve reporting zero exclusions.
* `print()` on a subset of an `mb_check()` audit reports counts for the
  printed rows, not the full table.
* `summary()` of a `pgt_ml` fit no longer errors when transitions
  contain `NA`; failed distance LPs are counted and warned about.
* The pollutant line is printed for every model (previously omitted for
  `envelope`) and by `summary()` as well.
* `pgt(pollutant = <character vector>)` gives a clear error instead of
  an internal condition-length failure.

## Housekeeping

* DMU-invariant quantities are computed once per fit instead of once
  per DMU (and once per LP solve in `boot_pgt()`).
* `Benchmarking` dropped from `Suggests` (never used).
* Example and vignette runtimes reduced.

# pgt 0.4.0

Panel productivity measurement and inference.

* `pgt_ml()`: the global Malmquist-Luenberger productivity index (Oh
  2010) under the pollution-generating technology, decomposed into
  efficiency change and best-practice change (`GML = EC x BPC`). The
  global reference frontier removes the cross-period infeasibility of the
  adjacent-period index. Marked experimental: the index under a
  materials-balance technology is not yet settled in the literature.
* `boot_pgt()` and `boot_pgt_sensitivity()`: subsampling (m out of n)
  confidence intervals for per-DMU scores and group means, with a
  subsample-size sensitivity check. The intervals are heuristic; the
  documentation states why.
* `pgt_tech()` gains a `period` argument for panel data.
* S3 methods for the new classes `pgt_ml` and `pgt_boot`.

# pgt 0.3.0

Competing axiom systems and a comparison harness.

* `pgt(model = "byprod")`: the by-production intersection technology of
  Murty, Russell and Levkoff (2012). Reports output efficiency, emission
  efficiency and their graph average. Unit tests reproduce the analytic
  scores of the paper's Example 1.
* `pgt(model = "mb_cost")`: the materials-balance cost model of Coelli,
  Lauwers and Van Huylenbroeck (2007), with the `EE = TE x EAE`
  decomposition into technical and environmental-allocative efficiency.
* `pgt(model = "wd")`: the weak-disposability model (Kuosmanen 2005
  correct VRS formulation) as a reference axiom system.
* `compare_models()`: fits the efficiency-scored models on one technology
  and reports Spearman rank agreement, worst-quartile overlap and the
  largest ranking disagreement.
* `pgt_tech()` gains a `polluting` argument marking the emission-causing
  inputs of the by-production model.
* Four package vignettes: the core workflow, the published-result
  replications, the cross-axiom comparison, and productivity with
  inference.

# pgt 0.2.0

The directional representation and richer material accounting.

* `pgt(model = "fdmo")`: the factorially determined multi-output
  (directional) representation of Rodseth (2025), Eq. 13 with abatement
  fixed at each DMU's own level (alias `"ddf"`). Jointly expands the good
  output and contracts the bad output; returns good-output, bad-output
  and gross efficiency. Unit tests reproduce the paper's Table 3.
* DMU-specific material flow coefficients: `u` may be an `L x N` matrix
  and `v` a length-`L` vector, so material or output quality can vary
  across units (Eder 2022, Rodseth 2025).
* Multi-pollutant technologies: `b` may be an `L x P` matrix, with `u` a
  named list and `v` a vector or matrix; each pollutant carries its own
  materials-balance identity, and a `pollutant` argument selects the one
  to estimate.
* `pgt_tech()` gains an `a` (abatement output) argument, required by the
  directional model.

# pgt 0.1.0

Initial release: the Rodseth (2025) weak-G-disposability slice.

* `pgt_tech()`: pollution-generating technology constructor with
  materials-balance metadata (material flow coefficients `u`, `v`).
* `mb_check()`: pre-estimation audit of the materials-balance identity
  `u'x - v y >= b`.
* `pgt()`: per-DMU minimum-bad-output estimation; models `"wgd"`
  (Rodseth 2025, Eq. 6 reduced form, alias `"wgd_rodseth"`) and
  `"envelope"` (Eq. 6 with inputs free: the convex lower (y, b)
  envelope); VRS/CRS; pooled or group reference sets.
* `pgt_decompose()`: metafrontier decompositions of environmental
  efficiency; `"envelope"` (Total = WR x TGR, exact identity, no
  feasibility screen) and `"rodseth"` (TE x Technology x AE from staged
  relaxations of the full WGD program).
* `shadow_prices()` and `mac_curve()`: constraint duals and marginal
  abatement cost curves.
* S3 `print`, `summary`, `plot`, `as.data.frame` methods; synthetic
  `steeldemo` panel; testthat suite with hand-computed LP checks and a
  cross-check against an independent reference implementation of the
  reduced-form program.
* `pigfarms` data: the numerical example of Rodseth (2025), Table 1.
  Unit tests replicate the paper's Table 2 minimal controlled
  emissions (16, 16, 16, 20, 16) and recover its abatement column from
  the `mb_check()` closure gap.
