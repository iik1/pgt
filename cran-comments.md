# cran-comments

## Submission

Update release: pgt 0.6.0 (previous CRAN version: 0.5.0).

This release redefines the package's principal estimator to follow its
source paper exactly: `pgt(model = "wgd")` now implements Equation 6 of
Rodseth (2025) as printed (inputs are decision variables), and the
previous programme survives unchanged as `model = "wgd_anchored"`. The
release also adds support for several intended outputs, the
five-component efficiency decomposition of Rodseth (2025, Eq. 11), and
a numerically robust retry ladder for the decomposition's stage
programmes. The breaking change is documented at the head of NEWS.md;
scores on data with a positive retained-content coefficient differ
from 0.5.x by design.

## Test environments

* local Windows 11 Enterprise, R 4.5.1 (x86_64-w64-mingw32)
* GitHub Actions (R-CMD-check workflow)

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
* boot_pgt() saves and restores the caller's .Random.seed when a seed
  is supplied, so the user's RNG state is never left altered.
* All \donttest examples run in under 5 seconds each on the test
  machine.
