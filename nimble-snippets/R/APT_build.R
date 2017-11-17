######################################################################################
## buildAPT was adapted from Nimble's buildMCMC to provide adaptive parallel tempering
######################################################################################

##' Create an APT function, from an MCMCconf object
##' 
##' Adapted from buildMCMC. Accepts a single required argument, which
##' may be of class MCMCconf, or inherit from class modelBaseClass (a
##' NIMBLE model object).  Returns an APT function; see details
##' section.
##' 
##' Calling buildAPT(conf,Temps,monitorTmax,ULT,thinPrintTemps) will
##' produce an uncompiled (R) APT function object, say 'myAPT'.
##' 
##' The uncompiled MCMC function will have arguments:
##' 
##' \code{niter} The number of iterations to run the MCMC.
##' 
##' \code{reset} Boolean specifying whether to reset the internal MCMC
##' sampling algorithms to their initial state (in terms of self-adapting
##' tuning parameters), and begin recording posterior sample chains anew.
##' Specifying \code{reset=FALSE} allows the MCMC algorithm to continue running
##' from where it left off, appending additional posterior samples to the
##' already existing sample chains. Generally, \code{reset=FALSE} should only
##' be used when the MCMC has already been run (default = TRUE).
##' 
##' \code{resetTempering} Boolean specifying whether to reset the
##' flexibility of the temperature ladder's adaptation process.
##' 
##' \code{simulateAll} Boolean specifying whether to simulate into all
##' stochastic nodes.  This will overwrite the current values in all stochastic
##' nodes (default = FALSE).
##' 
##' \code{time} Boolean specifying whether to record runtimes of the
##' individual internal MCMC samplers.  When \code{time=TRUE}, a vector of
##' runtimes (measured in seconds) can be extracted from the MCMC using the
##' method \code{mcmc$getTimes()} (default = FALSE).
##'
##' \code{adaptTemps} Boolean specifying whether the temperature
##' ladder will be adapted or not.
##'
##' \code{printTemps} Boolean specifying whether the temperature
##' ladder will be printed during the MCMC. THe print frequency is
##' controlled by thinPrintTemps.
##'
##' \code{tuneTemper1} Numeric tuning parameter of the adaptation
##' process of the temperature ladder. See source code for
##' buildAPT. Defaults to 10.
##' 
##' \code{tuneTemper2} Numeric tuning parameter of the adaptation
##' process of the temperature ladder. See source code for
##' buildAPT. Defaults to 1.
##' 
##' \code{progressBar} Boolean specifying whether to display a progress bar
##' during MCMC execution (default = TRUE).  The progress bar can be
##' permanently disabled by setting the system option
##' \code{nimbleOptions(MCMCprogressBar = FALSE)}.
##' 
##' Samples corresponding to the \code{monitors} and \code{monitors2} from the
##' MCMCconf are stored into the interval variables \code{mvSamples} and
##' \code{mvSamples2}, respectively. These may be accessed and converted into R
##' matrix objects via: \code{as.matrix(mcmc$mvSamples)}
##' \code{as.matrix(mcmc$mvSamples2)}
##' 
##' The uncompiled (R) MCMC function may be compiled to a compiled MCMC object,
##' taking care to compile in the same project as the R model object, using:
##' \code{Cmcmc <- compileNimble(Rmcmc, project=Rmodel)}
##' 
##' The compiled function will function identically to the uncompiled object,
##' except acting on the compiled model object.
##' 
##' @param conf An object of class MCMCconf that specifies the model, samplers,
##' monitors, and thinning intervals for the resulting MCMC function.  See
##' \code{configureMCMC} for details of creating MCMCconf objects.
##' Alternatively, \code{MCMCconf} may a NIMBLE model object, in which case an
##' MCMC function corresponding to the default MCMC configuration for this
##' model is returned.
##' @param Temps A numeric vector giving the initial temperature ladder. 
##' @param monitorTmax A logical indicator (default = TRUE) controlling if MCMC output should be stored at the hottest rung of the temperature ladder. Useful when monitorring the behaviour of APT. When TRUE mvSamples and mvSamples2 monitor T=1 and mvSamplesTmax and mvSamples2Tmax provide identically defined monitors (i.e. for exaclty the same nodes) for T=Tmax.
##' @param ULT A numberic value (default = 1E6) that provides an upper limit to temperature during APT. 
##' @param thinPrintTemps A numeric value controlling how often temperatures of the temperature ladder should be printed to screen when runtime parameter printTemps is TRUE. The default value of 1 is often too verbose. A good value to use is niter/10. 
##' 
##' @param ... Additional arguments to be passed to \code{configureMCMC} if
##' \code{conf} is a NIMBLE model object
## ##' @section Calculating WAIC:
## ##' 
## ##' After the MCMC has been run, calling the \code{calculateWAIC()} method of
## ##' the MCMC object will return the WAIC for the model, calculated using the
## ##' posterior samples from the MCMC run.
## ##' 
## ##' \code{calculateWAIC()} has a single arugment:
## ##' 
## ##' \code{nburnin} The number of iterations to subtract from the beginning of
## ##' the posterior samples of the MCMC object for WAIC calculation.  Defaults to
## ##' 0.
## ##' 
## ##' The \code{calculateWAIC} method calculates the WAIC of the model that the
## ##' MCMC was performed on. The WAIC (Watanabe, 2010) is calculated from
## ##' Equations 5, 12, and 13 in Gelman (2014).  Note that the set of all
## ##' parameters monitored by the mcmc object will be treated as \eqn{theta} for
## ##' the purposes of e.g. Equation 5 from Gelman (2014).  All parameters
## ##' downstream of the monitored parameters that are necessary to calculate
## ##' \eqn{p(y|theta)} will be simulated from the posterior samples of
## ##' \eqn{theta}.
##' @author David Pleydell (adapted from code by Daniel Turek).
####################################################################################################################
## ##' @examples                                                                                                  ##
## ##'                                                                                                            ##
## ##' \dontrun{                                                                                                  ##
## ##' code <- nimbleCode({                                                                                       ##
## ##'     mu ~ dnorm(0, 1)                                                                                       ##
## ##'     x ~ dnorm(mu, 1)                                                                                       ##
## ##' })                                                                                                         ##
## ##' Rmodel <- nimbleModel(code)                                                                                ##
## ##' conf <- configureMCMC(Rmodel)                                                                              ##
## ##' nIter <- 1E5                                                                                               ##
## ##' Rmcmc <- buildAPT(conf, Temps=1:10, monitorTmax=TRUE, thinPrintTemps=nIter/10)                             ##
## ##' Cmodel <- compileNimble(Rmodel)                                                                            ##
## ##' Cmcmc <- compileNimble(Rmcmc, project=Rmodel)                                                              ##
## ##' Cmcmc$run(nIter, reset=TRUE, resetTempering=TRUE, adaptTempts=TRUE, printTemps=TRUE, progressBar=FALSE)    ##
## ##' plot.tempTraj(Cmcmc)  ## Plots the trajectories of the temperature ladder                                  ##
## ##' Cmcmc$run(nIter, reset=FALSE, resetTempering=FALSE, adaptTempts=FALSE, printTemps=FALSE, progressBar=TRUE) ##
## ##' samples <- tail(as.matrix(Cmcmc$mvSamples), n=nIter)                                                       ##
## ##' summary(samples)                                                                                           ##
## ##' samplesTM <- tail(as.matrix(Cmcmc$mvSamplesTmax), n=nIter)                                                 ##
## ##' summary(samplesTM)                                                                                         ##
## ## ##' WAIC <- Cmcmc$calculateWAIC(nburnin = 1000)                                                             ##
## ##' }                                                                                                          ##
####################################################################################################################
#' @import nimble
#' @export 
buildAPT <- nimbleFunction(
    setup = function(conf,                ## As for buildMCMC 
                     Temps,               ## Vector of temperatures. Typically, lowest temperature should be 1.
                     monitorTmax=TRUE,    ## Logical, save MCMC output for Tmax too.
                     ULT=1E6,             ## Scalar, Upper Limit on Temperatures
                     thinPrintTemps=1,    ## Nb. iterations between prints of temps when adaptTemps==TRUE
                     ...) {
        if(inherits(conf, 'modelBaseClass')) {
            conf <- configureMCMC(conf, ...)
        } else if(!inherits(conf, 'MCMCconf')) {
            stop('conf must either be a nimbleModel or a MCMCconf object (created by configureMCMC(...) )')
        }
        if (missing(ULT)) {
            ULT <- 1E6
            nimPrint("ULT set to ", ULT)
        }
        model              <- conf$model
        my_initializeModel <- initializeModel(model)
        mvSaved            <- modelValues(model)                   ## Used to restore model following rejection in MCMC
        nSamplersPerT      <- length(seq_along(conf$samplerConfs)) ## Nb. samplers / temp'
        samplerFunctions   <- nimbleFunctionList(sampler_APT)
        ## Tempering related objects
        Temps        <- as.numeric(sort(Temps)) ## Ensures compiler knows Temps is not integer
        nTemps       <- length(Temps)           ## Number of temperatures for tempering
        invTemps     <- 1/Temps[1:nTemps]
        TempsCurrent <- Temps[1:nTemps]
        TempsOri     <- Temps[1:nTemps]
        nTemps_1     <- nTemps - 1        
        mvTemps      <- modelValues(model, nTemps) ## Stores state of MCMC at each temperature
        tempTraj     <- nimMatrix(0, nrow=1000, ncol=nTemps) ## Stores trajectories of temperatures
        logProbTemps <- numeric(nTemps)
        pSwapMatrix  <- nimMatrix(0, nTemps, nTemps)
        temporary    <- numeric(nTemps)
        accCountSwap <- nimMatrix(0, nTemps, nTemps)
        nimPrint("Initial temperatures:", Temps,"\n")
        ## 
        for (tt in 1:nTemps) {
            for(ss in 1:nSamplersPerT) { 
                ist <- ss + (tt-1) * nSamplersPerT ## index for Sampler & Temperature
                samplerFunctions[[ist]] <- conf$samplerConfs[[ss]]$buildSampler(model=model, mvSaved=mvSaved)
                samplerFunctions[[ist]]$setTemp(temp=Temps[tt])
            }
        }
        ##
        thin              <- conf$thin
        thin2             <- conf$thin2
        mvSamplesConf     <- conf$getMvSamplesConf(1)
        mvSamples2Conf    <- conf$getMvSamplesConf(2)
        monitors          <- processMonitorNames(model, conf$monitors)
        monitors2         <- processMonitorNames(model, conf$monitors2)
        mvSamples         <- modelValues(mvSamplesConf)   ## For storing MCMC output (T=1)
        mvSamples2        <- modelValues(mvSamples2Conf)  ## For storing MCMC output (T=Tmax)
        if (monitorTmax==TRUE) {
            mvSamplesTmax  <- modelValues(mvSamplesConf)  ## For MCMC output (T=Tmax) 
            mvSamples2Tmax <- modelValues(mvSamples2Conf) ## For MCMC output (T=Tmax) too
        }
        samplerTimes      <- c(0,0) ## Establish as a vector
        progressBarLength <- 52     ## Multiples of 4 only
        resetAnyway       <- TRUE
        totalIters        <- 1
    },
    #################################################
    run = function(niter          = integer(),
                   reset          = logical(default=TRUE),
                   resetTempering = logical(default=FALSE),
                   simulateAll    = logical(default=FALSE),
                   time           = logical(default=FALSE),
                   adaptTemps     = logical(default=TRUE),
                   printTemps     = double(default=FALSE),
                   tuneTemper1    = double(default=10), 
                   tuneTemper2    = double(default=1), 
                   progressBar    = logical(default=TRUE)) {
        if(simulateAll) {
            simulate(model)     ## Default behavior excludes data nodes
        }
        if(resetAnyway==TRUE) { ## Force reset to be TRUE on first usage to avoid segfault.
            reset          <- resetAnyway
            resetTempering <- resetAnyway
            resetAnyway   <<- !reset
            totalIters    <<- 1
        }
        my_initializeModel$run()
        if(resetTempering) {
            nimPrint("Resetting adaptation rate for tempering") 
            totalIters <<- 1
        }
        accCountSwap[,] <<- 0 * accCountSwap[,] 
        if(reset) {
            nimPrint("Resetting MCMC samplers and initial values set from mvTemps row 1") 
            for (tt in 1:nTemps) {
                nimCopy(from = model, to = mvTemps, row = tt, logProb = TRUE)
            }
            for(i in seq_along(samplerFunctions)) {
                samplerFunctions[[i]]$reset()
            }
            mvSamples_offset  <- 0
            mvSamples2_offset <- 0
            setSize(tempTraj,  niter/thin, nTemps)
            resize(mvSamples,  niter/thin) 
            resize(mvSamples2, niter/thin2)
            if (monitorTmax==TRUE) {
                resize(mvSamplesTmax,  niter/thin)
                resize(mvSamples2Tmax, niter/thin2)
            }
        } else {
            mvSamples_offset  <- getsize(mvSamples)
            mvSamples2_offset <- getsize(mvSamples2)
            setSize(tempTraj,  niter/thin, nTemps)
            resize(mvSamples,  mvSamples_offset  + niter/thin)
            resize(mvSamples2, mvSamples2_offset + niter/thin2)
            if (monitorTmax==TRUE) {
                resize(mvSamplesTmax,  mvSamples_offset  + niter/thin)
                resize(mvSamples2Tmax, mvSamples2_offset + niter/thin2)
            } 
        }
        ##### Monitors & Progress Bar #####
        if(dim(samplerTimes)[1] != length(samplerFunctions))            
            setSize(samplerTimes, length(samplerFunctions)) 
        if(niter < progressBarLength+3)
            progressBar <- progressBar & 0 ## avoids compiler warning
        if(progressBar) {
            for(iPB1 in 1:4) {
                cat('|'); for(iPB2 in 1:(progressBarLength/4)) cat('-')
            }
            print('|'); cat('|')
        }
        progressBarIncrement <- niter/(progressBarLength+3)
        progressBarNext      <- progressBarIncrement
        progressBarNextFloor <- floor(progressBarNext)
        ##########################
        ######## SAMPLING ########
        for(iter in 1:niter) {
            ## nimPrint(iter)
            checkInterrupt()
            ## ###############################################################
            ## Random Walk & Adaptation Phase:                              ##
            ## i.e. MCMC at each temperature, with or without time tracking ##
            if(time) { ## time == TRUE
                for(tt in 1:nTemps) {
                    ## Copy state of MCMC at temperature tt to model & mvSaved and continue
                    nimCopy(from=mvTemps, to=model,   row=tt,  logProb = TRUE) 
                    nimCopy(from=model,   to=mvSaved, rowTo=1, logProb = TRUE) 
                    for(ss in 1:nSamplersPerT) {
                        iST <- ss + (tt-1) * nSamplersPerT
                        samplerFunctions[[iST]]$setTemp(temp=Temps[tt])
                        samplerTimes[iST] <<- samplerTimes[iST] +
                            run.time(samplerFunctions[[iST]]$run())
                    }
                    ## Copy state of MCMC at temperature tt back to mvTemps
                    nimCopy(from=model, to=mvTemps, rowTo=tt, logProb = TRUE)
                    logProbTemps[tt] <<- model$getLogProb()
                } 
            } else {   ## time == FALSE 
                for(tt in 1:nTemps) {
                    ## browser()
                    ## Copy state of MCMC at temperature tt to model & mvSaved and continue
                    nimCopy(from=mvTemps, to=model,   row=tt,  logProb = TRUE) 
                    nimCopy(from=model,   to=mvSaved, rowTo=1, logProb = TRUE) 
                    for(ss in 1:nSamplersPerT) {
                        iST <- ss + (tt-1) * nSamplersPerT
                        samplerFunctions[[iST]]$setTemp(temp=Temps[tt])
                        samplerFunctions[[iST]]$run()
                    }
                    ## Copy state of MCMC at temperature tt back to mvTemps
                    nimCopy(from=model, to=mvTemps, rowTo=tt, logProb = TRUE)
                    logProbTemps[tt] <<- model$getLogProb()
                }   
            }       
            ## ################################################
            ## Random Swap Phase: Jumps between temperatures ##
            pSwapCalc() ## Updates pSwapMatrix            
            for (ii in 1:nTemps) { 
                if (iter %% 2 == 0) {
                    iFrom <- ii
                } else {
                    iFrom <- nTemps + 1 - ii
                }
                iTo <- rcat(n=1, prob=pSwapMatrix[iFrom, 1:nTemps]) 
                pProp <- pSwapMatrix[iFrom,iTo]
                pRev  <- pSwapMatrix[iTo,iFrom]
                lMHR  <- (invTemps[iTo]-invTemps[iFrom]) * 
                    (logProbTemps[iFrom]-logProbTemps[iTo]) + log(pRev) - log(pProp)
                lu    <- log(runif(1))
                if (lu < lMHR) { 
                    ## Swap model values
                    nimCopy(from=mvTemps, to=mvSaved, row=iTo,   rowTo=1,     logProb = TRUE)
                    nimCopy(from=mvTemps, to=mvTemps, row=iFrom, rowTo=iTo,   logProb = TRUE)
                    nimCopy(from=mvSaved, to=mvTemps, row=1,     rowTo=iFrom, logProb = TRUE)
                    ## Swap logProbTemps elements
                    temporary[1]        <<- logProbTemps[iTo]
                    logProbTemps[iTo]   <<- logProbTemps[iFrom]
                    logProbTemps[iFrom] <<- temporary[1]
                    ## Swap rows of pSwapMat
                    temporary[1:nTemps]         <<- pSwapMatrix[iTo,1:nTemps]
                    pSwapMatrix[iTo,1:nTemps]   <<- pSwapMatrix[iFrom,1:nTemps]                    
                    pSwapMatrix[iFrom,1:nTemps] <<- temporary[1:nTemps]
                    ## Swap cols of pSwapMat
                    temporary[1:nTemps]         <<- pSwapMatrix[1:nTemps,iTo]
                    pSwapMatrix[1:nTemps,iTo]   <<- pSwapMatrix[1:nTemps,iFrom]
                    pSwapMatrix[1:nTemps,iFrom] <<- temporary[1:nTemps]
                    ## Acceptance counter
                    accCountSwap[iFrom, iTo]    <<- accCountSwap[iFrom, iTo] + 1 
                } 
            }
            ## nimPrint(pSwapMatrix)
            ## nimPrint(Temps, "     ", invTemps)
            ###################################
            ## Temperature adaptation scheme ##
            totalIters <<- totalIters + 1
            if (adaptTemps) {
                ## gammaTSA <- 1 / ((totalIters/thin + 3) ^ 0.8) 
                gammaTSA    <- 1 / ((totalIters/tuneTemper1 + 3) ^ tuneTemper2)  
                TempsCurrent[1:nTemps] <<- Temps
                for (ii in 1:nTemps_1) {
                    accProb          <- min(1, exp( (invTemps[ii+1]-invTemps[ii]) * (logProbTemps[ii]-logProbTemps[ii+1]) ))
                    Temps[ii+1]     <<- Temps[ii] + (TempsCurrent[ii+1] - TempsCurrent[ii]) * exp(gammaTSA*(accProb-0.234))
                    if (Temps[ii+1] > ULT)
                        Temps[ii+1] <<- ULT + ii ## This prevents the temperature ladder from exploding 
                } 
                for (ii in 1:nTemps_1) 
                    invTemps[ii+1]  <<- 1 / Temps[ii+1]
                ## nimPrint("iter: ", iter)
                ## nimPrint("logProbTemps: ", logProbTemps)
                ## nimPrint(accCountSwap)
            }            
            ##############################################
            ## Automatically drop unneeded temperatures ##
            ## Not implimented                          ##                         
            ##                                         
            #################
            ## MCMC Output ##
            if(iter %% thin  == 0) {
                nimCopy(from = mvTemps, to = mvSamples,
                        row = 1, rowTo = mvSamples_offset + iter/thin,  nodes = monitors)
                tempTraj[iter/thin, 1:nTemps] <<- Temps[1:nTemps]
                if(printTemps == 1.0) {
                    if(totalIters %% thinPrintTemps == 0) {
                        nimPrint(iter, asRow(Temps))  
                        if(!adaptTemps)
                            printTemps <- 0.0
                    }
                }
            }
            if(iter %% thin2 == 0) 
                nimCopy(from = mvTemps, to = mvSamples2, 
                        row = 1, rowTo = mvSamples2_offset + iter/thin2, nodes = monitors2) 
            if (monitorTmax==TRUE) {
                if(iter %% thin  == 0)
                    nimCopy(from = mvTemps, to = mvSamplesTmax, row = nTemps,
                            rowTo = mvSamples_offset  + iter/thin,  nodes = monitors) 
                if(iter %% thin2 == 0)
                    nimCopy(from = mvTemps, to = mvSamples2Tmax, row = nTemps,
                            rowTo = mvSamples2_offset + iter/thin2, nodes = monitors2)
            } 
            ## Progress Bar
            if(progressBar & (iter == progressBarNextFloor)) {
                cat('-')
                progressBarNext <- progressBarNext + progressBarIncrement
                progressBarNextFloor <- floor(progressBarNext)
            }
        }
        if(progressBar) print('|')
        nimCopy(from = mvTemps, to = model, row = 1, logProb = TRUE) ## Ensures model is parameterised with T=1 samples prior to exiting.
    },
    ## #########################################################
    methods = list(
        getTimes = function() {
            returnType(double(1))
            return(samplerTimes)
        }, 
        pSwapCalc = function() {
            ## Fill matrix
            for (ii in 2:nTemps) {
                for (jj in 1:(ii-1)) {
                    pSwapMatrix[ii,jj] <<- exp(-abs(logProbTemps[ii] - logProbTemps[jj]))
                    pSwapMatrix[jj,ii] <<- pSwapMatrix[ii,jj]
                }
            }
            ## Normalise rows            
            for (ii in 1:nTemps) {
                rowSum <- sum(pSwapMatrix[ii,1:nTemps])
                if(rowSum>0) {
                    for (jj in 1:nTemps)
                        pSwapMatrix[ii,jj] <<- pSwapMatrix[ii,jj] / rowSum
                } else {
                    for (jj in 1:nTemps) {
                        pSwapMatrix[ii,jj] <<- (ii!=jj) / nTemps_1
                    }
                }
            }
        },
        mvTemps2model = function(row = double()) {
            ## Useful in R for exploring node values at each temp
            nimCopy(mvTemps, model, row=row, logProb=TRUE)
        }
    )    ## ,
    ## where = getLoadingNamespace()
    ## where = getNamespace("NimbleSnippets")    
)


## This is a function that will weed out missing indices from the monitors
## processMonitorNames <- function(model, nodes){
## 	isLogProbName <- grepl('logProb', nodes)
## 	expandedNodeNames <- model$expandNodeNames(nodes[!isLogProbName])
## 	origLogProbNames <- nodes[isLogProbName]
## 	expandedLogProbNames <- character()
## 	if(length(origLogProbNames) > 0){
## 		nodeName_fromLogProbName <- gsub('logProb_', '', origLogProbNames)
## 		expandedLogProbNames <- model$modelDef$nodeName2LogProbName(nodeName_fromLogProbName)
## 	}
## 	return( c(expandedNodeNames, expandedLogProbNames) )
## }



