## Resubmission

This is a resubmission. In this version I have:

* Replaced `print()` with `stop()` in `R/fit_development_pattern.R` for invalid
  input handling, as console output should use `message()`, `warning()` or
  `stop()` rather than `print()` or `cat()`
* Removed `set.seed()` calls with hardcoded values from `R/create_simulations.R`
  so that users control the seed themselves
* Replaced `\dontrun{}` with `\donttest{}` for examples that are long-running
  but executable. Unwrapped examples in `R/layer_loss.R` which run in under 5
  seconds using the bundled `losses` dataset

## R CMD check results

0 errors | 0 warnings | 3 notes

The 3 notes are all local machine artefacts and will not appear on CRAN servers:

* `New submission` — expected for a resubmission
* `unable to verify current time` — network issue on local machine
* `README.md or NEWS.md cannot be checked without 'pandoc' being installed` —
  pandoc is not installed locally but is available on CRAN servers
