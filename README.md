nimbleAPT
===============

A collection of functions for sampling the posterior distributions of Bayesian models (written in Nimble's implementation of the BUGS language) via an **adaptive parallel tempering (APT)** algorithm based on Miasojedow's [2013](https://www.tandfonline.com/doi/abs/10.1080/10618600.2013.778779) and [2016](https://link.springer.com/article/10.1007/s11222-015-9579-0) papers.

The package was previously released under the name *nimble-snippets*. 

Installation
------------
The package can be installed from github using the R package 'remotes' as follows.
remotes::install_git(url="https://github.com/DRJP/nimble-snippets.git", subdir="nimbleAPT", build_vignettes = TRUE)

Switching to 'build_vignettes = FALSE' makes for a faster installation, and thus may be a prefered option for users already familiar with the package. 


Citation
--------
If you use this package in you analyses please cite it using the following DOI.
[![DOI](https://zenodo.org/badge/106051349.svg)](https://zenodo.org/badge/latestdoi/106051349)

