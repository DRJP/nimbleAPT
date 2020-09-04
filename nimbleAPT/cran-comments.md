## Test environments
* local OS X install, R 3.6.2
* ubuntu 14.04 (on travis-ci), R 3.6.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.




## Here are the results from (rhub::check_for_cran())$cran_summary(). There are two notes: 
## 1) the check prints the Maintainer information from DESCRIPTION - these details are correct.
## 2) The check fails to verify the timestep - this is not a package issue, but a broken URL used by the check program.

> results$cran_summary()
For a CRAN submission we recommend that you fix all NOTEs, WARNINGs and ERRORs.
## Test environments
- R-hub windows-x86_64-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)

## R CMD check results
❯ On windows-x86_64-devel (r-devel), ubuntu-gcc-release (r-release), fedora-clang-devel (r-devel)
  checking CRAN incoming feasibility ... NOTE
  
  Maintainer: 'David Pleydell <david.pleydell@inrae.fr>'
  New submission

❯ On windows-x86_64-devel (r-devel), ubuntu-gcc-release (r-release)
  checking for future file timestamps ... NOTE
  unable to verify current time

0 errors ✔ | 0 warnings ✔ | 2 notes ✖
