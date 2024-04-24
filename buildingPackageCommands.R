sessionInfo()

detach("package:nimbleAPT")
unloadNamespace("nimbleAPT")

setwd(here::here('nimbleAPT'))
roxygen2::roxygenise()

setwd(here::here())
devtools::build("nimbleAPT", vignettes=FALSE) # TRUE # Can save time with FALSE

# devtools::check("nimbleAPT")

devtools::install("nimbleAPT", build_vignettes = FALSE) # TRUE

library("nimbleAPT")

help.start()

setwd("nimbleAPT")
devtools::build_vignettes()
