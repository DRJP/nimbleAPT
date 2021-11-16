## How to write your own R package and publish it on CRAN
## https://www.mzes.uni-mannheim.de/socialsciencedatalab/article/r-package/
## Work through the checks listed in this blog post

## Tests the package can build and install on various platforms
library(rhub)
library(here)
setwd(here::here("nimbleAPT"))
getwd()

rhub::check()
## 11h53 ran command - 12h03 preperror email
## 15: Ubuntu Linux 20.04.1 LTS, R-devel, GCC (ubuntu-gcc-devel)

## Info on last submission to rhub
rhub::last_check()
rhub::last_check()$cran_summary()


# From https://www.mzes.uni-mannheim.de/socialsciencedatalab/article/r-package/#subsection6-1
# Check for CRAN specific requirements using rhub and
# save it in the results objects
#
# results <- (rhub::check_for_cran())$cran_summary()
results <- rhub::check_for_cran()
# Get the summary of your results
results$cran_summary()

## 18h25 launched
## 18h email 1
## 18h email 2
## 18h email 3


# Alternatively
devtools::check_rhub()
devtools::check()

# Or... alternatively, perhaps this can help...
# https://stackoverflow.com/questions/63586750/error-with-ubuntu-linux-when-running-check-rhub-due-to-xml-library
# platforms()
res  = rhub::check(platform="ubuntu-gcc-release")
#res2 = rhub::check(platform="ubuntu-gcc-devel")
## SUGGESTS IT'S AN EIGEN or RHUB PROBLEM
##
## 13h48ish ran command - 13h4
res$cran_summary()



# Generate your cran-comments.md, then you copy-paste the output from the function above
usethis::use_cran_comments()

# Spell checkl
devtools::spell_check()

# Checks for windows version
devtools::check_win_devel()

## Submit to CRAN with the following
devtools::release()
