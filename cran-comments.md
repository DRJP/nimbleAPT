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


## In response to Gregor Seyer email, "CRAN submission nimbleAPT 1.0.3", 11 September 2020 05:28

_Please omit the space within the doi specification to make it clickable._

This has been done.

_Please add small executable examples in your Rd-files to illustrate the
use of the exported function but also enable automatic testing._

More examples have been added.


_Please always add all authors, contributors and copyright holders in the
Authors@R field with the appropriate roles.
 From CRAN policies you agreed to:
"... Where code is copied (or derived) from the work of others
(including from R itself), care must be taken that any copyright/license
statements are preserved and authorship is not misrepresented.
Preferably, an ‘Authors@R’ would be used with ‘ctb’ roles for the
authors of such code. Alternatively, the ‘Author’ field should list
these authors as contributors.
Where copyrights are held by an entity other than the package authors,
this should preferably be indicated via ‘cph’ roles in the ‘Authors@R’
field, or using a ‘Copyright’ field (if necessary referring to an
inst/COPYRIGHTS file)."
e.g.: Daniel Turek
Please explain in the submission comments what you did about this issue._

I have contacted the NIMBLE Development Team regarding this issue and they have listed four authors who contibuted to the original NIMBLE code which I modified to make this package.
So I have added an Authors@R section in the DESCRIPTION file.
This now lists myself (as the package developer) and the four members of the NIMBLE development team.
These contributions have been cited using "ctb" (as suggested by Gregor Seyer of CRAN), however, the NIMBLE development team have suggested that "cph" would be an appropriate alternative since they do not consider themselves authors of this package.
I have included a inst/COPYRIGHTS file which cites the appropriate lines from the inst/COPYRIGHTS file of the NIMBLE pacakge.
In NIMBLE, only C++ code is licenced under GPL (>=2), therefore, since I have only adapted R code (and not C++ code) from NIMBLE I only need to cite the BSD_3_clause licence in inst/COPYRIGHTS.
I have also added a NOTE to my DESCRIPTION file that clarifies the contributions of each contributor.
I hope you find these modifications accepatable.
