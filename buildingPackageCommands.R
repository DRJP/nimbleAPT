sessionInfo()

detach("package:nimbleAPT")
unloadNamespace("nimbleAPT")

setwd(here::here('nimbleAPT'))
roxygen2::roxygenise()

setwd(here::here())
devtools::build("nimbleAPT", vignettes=TRUE) # FALSE

devtools::check("nimbleAPT")

devtools::install("nimbleAPT", build_vignettes = TRUE) # FALSE

help.start()

setwd("nimbleAPT")
devtools::build_vignettes()
