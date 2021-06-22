## How to write your own R package and publish it on CRAN
## https://www.mzes.uni-mannheim.de/socialsciencedatalab/article/r-package/
## Work through the checks listed in this blog post

## Tests the package can build and install on various platforms
library(rhub)
library(here)
setwd(here())
setwd("nimbleAPT")

rhub::check()

# From https://www.mzes.uni-mannheim.de/socialsciencedatalab/article/r-package/#subsection6-1
# Check for CRAN specific requirements using rhub and
# save it in the results objects
results <- (rhub::check_for_cran())$cran_summary()
# Get the summary of your results
results$cran_summary()

# Generate your cran-comments.md, then you copy-paste the output from the function above
usethis::use_cran_comments()

## Submit to CRAN with the following
devtools::release()
