# cran-comments

## Submission

Update release: pgt 0.7.1 (previous CRAN version: 0.7.0, published
2026-09-22).

This update follows 0.7.0 closely because it corrects results that
0.7.0 can return wrongly:

* `pgt(model = "fdmo")` could report a negative projected emission
  when material flow coefficients differ across producers: a unit
  whose inputs carry little pollutant could take a pollutant-rich
  peer's output. The programme now bounds the bad-output contraction
  by the unit's own emission. The bound cannot bind under common
  coefficients with closed accounts, so the published replication of
  Rodseth (2025, Table 3) is unchanged.
* `pgt(model = "wgd")` reported an implied uncontrolled emission
  `z_star` that omitted the retained content of any output produced
  beyond the unit's own level, so `z_star - a_star` did not equal
  `b_star` whenever such an overshoot occurred with a positive
  retained-content coefficient. `b_star` and the scores were correct
  and are unchanged.
* `pgt(model = "wgd_input_fixed")` could report an infeasible
  programme as a numerical failure (status 5) rather than status 2.
  Every solve of the `"wgd"`, `"envelope"`, `"wgd_input_fixed"` and
  decomposition-stage programmes is now checked against the bound
  that self-reference gives, and a failed or out-of-bound solve is
  retried under alternative lp_solve settings. Results that passed at
  the first attempt are unchanged.
* `compare_models()` ranked efficient units' scores (1 up to solver
  noise) by that noise, so its Spearman matrix depended on the solver
  build; scores within 1e-8 now rank as ties.

NEWS.md lists the remaining changes: a new test file of large-magnitude
cases, extended precomputed Monte Carlo results in `inst/simulations/`
(read by a vignette; the study script is not run during checks), and
documentation.

## Test environments

* local Windows 11 Enterprise, R 4.5.1 (x86_64-w64-mingw32)
* GitHub Actions: Windows (release), macOS (release), Ubuntu (devel,
  release, oldrel-1)

## R CMD check results

0 errors | 0 warnings | 1 note

The note is the incoming-feasibility check's "Days since last update",
explained under Submission above.

A second local note, "unable to verify current time", is an artifact of
the checking machine's restricted network access and does not concern
the package.

## Notes for the reviewers

* Estimators are validated against the numerical examples printed in
  their source papers (Rodseth 2025, Journal of Productivity Analysis;
  Murty, Russell and Levkoff 2012, JEEM); the replications run as unit
  tests and are reproduced in vignette("replication", "pgt").
* The two shipped US data sets are built from public EIA and EPA files
  by scripts in data-raw/ (not part of the package); the package itself
  needs no internet access.
* boot_pgt() saves and restores the caller's .Random.seed when a seed
  is supplied, so the user's RNG state is never left altered.
* All examples run in well under 5 seconds each on the test machine.
