sessionInfo()

detach("package:nimbleAPT")
unloadNamespace("nimbleAPT")

setwd('/home/pleydell/nimbleProjects/nimbleAPT/nimbleAPT/')
roxygen2::roxygenise()

setwd('/home/pleydell/nimbleProjects/nimbleAPT')
devtools::build("nimbleAPT", vignettes=TRUE) # FALSE

devtools::check("nimbleAPT")

devtools::install("nimbleAPT", build_vignettes = TRUE) # FALSE

help.start()
