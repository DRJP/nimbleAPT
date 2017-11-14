######################################################################################################################
## This is a DEBUG version of APT_demo.R with everything, except lines required to reproduce the bug, commented out ##
######################################################################################################################

library(coda)
## library(nimble)    
## library(NimbleSnippets) ## The bug is related to the namespace of functions in this package - it reports a class is not found
                           ## Hack: if we source the R files of this package, so the functions are in .GlobalEnv, not error is returned

####################
## Some Functions ##
####################
k2theta <- nimbleFunction(
    run = function(k = double(),              ## A categorical random variable in continuous form.
                   K = double(),              ## Upper limit of k.
                   randomEffects = double(1), ## A set of "randomEffects" which will be fixed during MCMC in order to generate poor mixing
                   theta0=double()            ## An intercept term
                   ) {
        ## Uses a categorical rv to apply a shift in the value of a parameter theta
        k <- ceiling(k)
        if (k < 1 | k > K)
            return(-Inf) 
        theta <- theta0 + randomEffects[k]
        return(theta)
        returnType(double())
    }
)

k2sigma <- nimbleFunction(
    run = function(k = double(), K = double(), randomEffects = double(1), sigma0=double()) {
        ## Uses a categorical rv to determine a change in a standard deviation parameter theta
        k <- ceiling(k)
        if (k < 1 | k > K)
            return(-Inf) 
        sigma <- abs(sigma0 + randomEffects[k])
        return(sigma)
        returnType(double())
    }
)

##################################
## BUGS code for the 'right' model
##################################
rightCode <- nimbleCode({
    theta ~ dnorm(0, sd=1E6)
    for (i in 1:N) 
        y[i] ~ dnorm(mean=theta, sd=sigma0)
})
##################################
## BUGS code for the 'wrong' model
##################################
wrongCode <- nimbleCode({
    k      ~ dcat(pk[1:K])
    theta0 ~ dnorm(0, sd=1E6)
    theta <- k2theta(k, K, rfx_theta[1:K], theta0) 
    sigma0 ~ dexp(1)
    sigma <- abs(k2sigma(k, K, rfx_sigma[1:K], sigma0)) 
    for (i in 1:N)
        y[i] ~ dnorm(mean=theta, sd=sigma)
})


##############################################
## Parameters, constants and initial values ##
##############################################
N         <- 50  ## low N => MCMC jumps between modes easily, high N => MCMC never jumps between modes, intermediate N => MCMC occasionally jumps between modes.
theta0    <- 0
sigma0    <- 1
K         <- 10
pk        <- rep(1/K, K)
rfx_theta <- rnorm(K)
rfx_sigma <- rnorm(K)
##
rConstants <- list(N=N, pk=pk, sigma0=sigma0)
rInits     <- list(theta=theta0) 
##
wConstants <- list(N=N, K=K)
wInits     <- list(theta0=theta0, sigma0=sigma0, rfx_theta=rfx_theta, rfx_sigma=rfx_sigma, k=1, pk=pk) 

###########################
## Initialise the models ##
###########################
right <- nimbleModel(rightCode, constants=rConstants, inits=rInits)
wrong <- nimbleModel(wrongCode, constants=wConstants, inits=wInits)

#############################################
## Simulate data but not the stochastic nodes
#############################################
set.seed(11)
right$simulate('y')
wrong$simulate(wrong$getDependencies(c('theta0','sigma0'), self=FALSE), includeData = TRUE)

#####################################################
## Transfer data simulated with 'right' to 'wrong' ##
#####################################################
Data    <- list(y=right$y)
wrong   <- nimbleModel(wrongCode, constants=wConstants, inits=wInits, data=Data) ## For standard MCMC samplers
wrong2  <- nimbleModel(wrongCode, constants=wConstants, inits=wInits, data=Data) ## For APT samplers
wrong$calculate()  ## Total logProb of wrong
wrong2$calculate() ## Total logProb of wrong2

#########################
## Compile wrong model ##
#########################
cWrong2 <- compileNimble(wrong2) 

###################################
## Setup APT samplers for wrong2 ##
###################################
nimCopy(wrong, cWrong2, logProb=TRUE)
apt <- configureMCMC(wrong2, print=TRUE)
apt$removeSamplers()
apt$addSampler("theta0", type="sampler_RW_tempered", control=list(temperPriors=TRUE))
apt$addSampler("sigma0", type="sampler_RW_tempered", control=list(temperPriors=TRUE))
apt$addSampler("k",      type="sampler_slice_tempered", control=list(temperPriors=TRUE))
apt$addSampler(c("theta0","sigma0"), type="sampler_RW_block_tempered", control=list(temperPriors=TRUE))
apt$printSamplers()
## Add monitors for logProbs
(allLogProbs <- wrong2$getVarNames(TRUE)[grep("logProb", wrong2$getVarNames(TRUE))])
apt$addMonitors2(allLogProbs)
apt$getMonitors()
apt$getMonitors2()

#################################################
## HERE COMES THE BUG                          ##
## Build and compile the MCMC & APT algorithms ##
#################################################
APT <- buildAPT(apt, Temps=1:7, monitorTmax=TRUE)  
#################################################################
## Error in getClass(Class, where = topenv(parent.frame())) :  ##
##   "nfRefClass_R_GlobalEnv45" is not a defined class         ##
#################################################################





################################################
## THE HACK: source the R files and try again ##
################################################
(packDir <- "~/nimbleProject/nimble-snippets/nimble-snippets") ##find.package("NimbleSnippets"))
source(paste0(packDir, "/R/APT_samplers.R"))
source(paste0(packDir, "/R/APT_build.R"))


###################################
## Setup APT samplers for wrong2 ##
###################################
nimCopy(wrong, cWrong2, logProb=TRUE)
apt <- configureMCMC(wrong2, print=TRUE)
apt$removeSamplers()
apt$addSampler("theta0", type="sampler_RW_tempered", control=list(temperPriors=TRUE))
apt$addSampler("sigma0", type="sampler_RW_tempered", control=list(temperPriors=TRUE))
apt$addSampler("k",      type="sampler_slice_tempered", control=list(temperPriors=TRUE))
apt$addSampler(c("theta0","sigma0"), type="sampler_RW_block_tempered", control=list(temperPriors=TRUE))
apt$printSamplers()
## Add monitors for logProbs
(allLogProbs <- wrong2$getVarNames(TRUE)[grep("logProb", wrong2$getVarNames(TRUE))])
apt$addMonitors2(allLogProbs)
apt$getMonitors()
apt$getMonitors2()

#################################################
## Build and compile the MCMC & APT algorithms ##
#################################################
APT  <- buildAPT(apt, Temps=1:7, monitorTmax=TRUE)  ## THIS NOW BUILDS THE APT WITHOU COMPLAINT ##
cAPT <- compileNimble(APT)

###############
## APT in C ##
###############
nIter   <- 1E4
cAPT$thinPrintTemps <- floor(nIter/10)       ## Limits printing temperature ladder to 10 lines
cAPT$run(nIter, printTemps=TRUE, progressBar=FALSE)
samples <- as.matrix(cAPT$mvSamples)
codaAPT <- as.mcmc(samples)

## Plot 1: trajectories
dev.set()
plot(codaAPT)                              
## Plot 2: autocorrelation
dev.set()
autocorr.plot(codaAPT)
## Plot 3: sigma0~theta0
dev.set()
par(mfrow=c(1,1))
plot(samples[,"theta0"],samples[,"sigma0"], xlab="theta0", ylab="sigma0")
## Plot 4: temperature ladder trajectories
dev.set()
plotTempTraj(cAPT)

