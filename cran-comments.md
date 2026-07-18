# cran-comments

## Submission

First submission of pgt 0.5.0 to CRAN.

The package implements nonparametric efficiency analysis for
pollution-generating technologies under the materials-balance principle
(weak-G-disposability, by-production, materials-balance cost and
weak-disposability estimators, with a feasibility audit, decompositions,
shadow prices, a productivity index and subsampling inference).

## Test environments

* local Windows 11 Enterprise, R 4.5.1 (x86_64-w64-mingw32)
* GitHub Actions (R-CMD-check workflow)

## R CMD check results

0 errors | 0 warnings | 1 note

* "New submission": this is the package's first CRAN submission.

A second local note, "unable to verify current time", is an artifact of
the checking machine's restricted network access and does not concern
the package.

## Notes for the reviewers

* Estimators are validated against the numerical examples printed in
  their source papers (Rodseth 2025, Journal of Productivity Analysis;
  Murty, Russell and Levkoff 2012, JEEM); the replications run as unit
  tests and are reproduced in vignette("replication", "pgt").
* boot_pgt() saves and restores the caller's .Random.seed when a seed
  is supplied, so the user's RNG state is never left altered.
* All \donttest examples run in under 5 seconds each on the test
  machine.
