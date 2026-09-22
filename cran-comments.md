# cran-comments

## Submission

Update release: pgt 0.7.0 (previous CRAN version: 0.6.1).

This release adds a second measured data set, `usfgd` (154 US
coal-fired plants with flue-gas desulphurisation in 2023, with an
abatement output derived from EIA-923 and EPA CAMD records), and
regenerates the synthetic `steeldemo` panel with observed carbon
capture: three new columns (`captured`, `capture_energy`,
`abatement_tech`), an account that closes exactly in every row, and
new values in every column. Every result computed on `steeldemo`
therefore differs from 0.6.1; the head of NEWS.md says so. No estimator
changed its programme. Internally, every linear programme is now solved
on a copy of the technology rescaled to unit magnitude and the level
fields scaled back, which removes the solver failures and suboptimal
vertices that lp_solve produced at tonne magnitudes; scores are ratios
and are unaffected where the previous solve succeeded. The
documentation now cites the directional model as Rodseth (2025) Eq. 14
(previously Eq. 13) and no longer describes the materials-balance
account as "enforced", stating instead how it enters each programme.

## Test environments

* local Windows 11 Enterprise, R 4.5.1 (x86_64-w64-mingw32)
* GitHub Actions: Windows (release), macOS (release), Ubuntu (devel,
  release, oldrel-1)

## R CMD check results

0 errors | 0 warnings | 0 notes

A local note, "unable to verify current time", is an artifact of the
checking machine's restricted network access and does not concern the
package.

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
