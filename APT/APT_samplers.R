####################################################################
### virtual nimbleFunction template, included for ALL samplers #####
####################################################################

sampler_APT <- nimbleFunctionVirtual(
    ## run = function(temperture=double(0, default=1)) {}, 
    methods = list(
        setTemp           = function(temp = double()) {},
        ## turnOffAdaptation = function() {},
        reset             = function() {}
    )
)


####################################################################
### scalar RW sampler with normal proposal distribution ############
####################################################################

sampler_RW_tempered <- nimbleFunction(
    contains = sampler_APT,
    setup = function(model, mvSaved, target, control) {
        ## Control list extraction
        logScale      <- control$log
        reflective    <- control$reflective
        adaptive      <- control$adaptive
        adaptInterval <- control$adaptInterval
        scale         <- control$scale
        temperPriors  <- control$temperPriors
        if (is.null(temperPriors))
            stop("control$temperPriors unspecified in sampler_RW_tempered")
        ## Node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes      <- model$getDependencies(target) 
        targetNode     <- model$expandNodeNames(target)
        dependantNodes <- calcNodes[!is.element(calcNodes, targetNode)]
        ## Numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        optimalAR     <- 0.44
        gamma1        <- 0
        range         <- unlist(getDistributionInfo(model$getDistribution(target))$range) ## getDistribution(model$getNodeDistribution(target))$range
        ## Checks
        if(length(targetAsScalar) > 1) stop('cannot use RW sampler on more than one target; try RW_block sampler')
        if(model$isDiscrete(target))   stop('cannot use RW sampler on discrete-valued target; try slice sampler')
        if(logScale & reflective)      stop('cannot use reflective RW sampler on a log scale (i.e. with options log=TRUE and reflective=TRUE')
        ## Initialise temperature
        temperature <- nimNumeric(length=2, value=1) ## Length 2 is just a hack. Only 1st element is used.
    },  
    run = function() { 
        ## nimPrint(temperature[1])
        currentValue <- model[[target]]
        propLogScale <- 0
        if(logScale) {
            propLogScale <- rnorm(1, mean = 0, sd = scale)
            propValue    <- currentValue * exp(propLogScale)
        } else {
            propValue <- rnorm(1, mean = currentValue,  sd = scale)
        }
        if(reflective) { 
            while(propValue < range[1] | propValue > range[2]) {
                if(propValue < range[1]) propValue <- 2*range[1] - propValue
                if(propValue > range[2]) propValue <- 2*range[2] - propValue
            } 
        } 
        if (temperPriors) {
            model[[target]] <<- propValue 
            logMHR <- calculateDiff(model, calcNodes) / temperature[1] + propLogScale ## Original. Tempers everything.
        } else {
            logPriorWeightOri   <- getLogProb(model, target)
            model[[target]]    <<- propValue
            diffLogProb         <- calculateDiff(model, calcNodes)
            logPriorWeightProp  <- getLogProb(model, target)
            diffLogPriorWeights <- logPriorWeightProp - logPriorWeightOri            
            logMHR <- (diffLogProb - diffLogPriorWeights) / temperature[1] + diffLogPriorWeights + propLogScale ## POSSIBLY BUGGY ???            
        }
        jump <- decide(logMHR)
        if (jump) {
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        } else {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        }
        if(adaptive) {
            adaptiveProcedure(jump)
        }
    },                
    methods = list(   
        setTemp = function(temp=double()) {
            temperature[1] <<- temp
        },            
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                gamma1 <<- 1/((timesAdapted + 3)^0.8)
                gamma2 <- 10 * gamma1
                adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
                scale <<- scale * adaptFactor
                timesRan <<- 0
                timesAccepted <<- 0
            }
        },
        ## turnOffAdaptation = function() {
        ##     adaptive <<- FALSE
        ## },
        reset = function() {
            scale         <<- scaleOriginal
            timesRan      <<- 0
            timesAccepted <<- 0
            timesAdapted  <<- 0
            gamma1        <<- 0
        }
    ), where = getLoadingNamespace()
)


########################################################################
### block RW sampler with multi-variate normal proposal distribution ###
########################################################################

sampler_RW_block_tempered <- nimbleFunction(
    contains = sampler_APT,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive       <- control$adaptive
        adaptScaleOnly <- control$adaptScaleOnly
        adaptInterval  <- control$adaptInterval
        scale          <- control$scale
        propCov        <- control$propCov
        temperPriors   <- control$temperPriors
        if (is.null(temperPriors)) {
            stop("control$temperPriors unspecified in sampler_RW_block_tempered")
        }
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes      <- model$getDependencies(target)
        targetNodes    <- model$expandNodeNames(target)
        ## numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        d <- length(targetAsScalar)
        if(is.character(propCov) && propCov == 'identity') {
            propCov <- diag(d)
        }
        propCovOriginal <- propCov
        chol_propCov <- chol(propCov)
        chol_propCov_scale <- scale * chol_propCov
        empirSamp <- matrix(0, nrow=adaptInterval, ncol=d)
        ## nested function and function list definitions
        my_setAndCalculateDiff  <- setAndCalculateDiff(model, target)
        my_decideAndJump        <- decideAndJump(model, mvSaved, calcNodes)
        my_calcAdaptationFactor <- calcAdaptationFactor(d)
        ## checks
        if(class(propCov)      != 'matrix')  stop('propCov must be a matrix\n')
        if(class(propCov[1,1]) != 'numeric') stop('propCov matrix must be numeric\n')
        if(!all(dim(propCov)   == d))        stop('propCov matrix must have dimension ', d, 'x', d, '\n')
        if(!isSymmetric(propCov))            stop('propCov matrix must be symmetric')
        ## Initialise temperature
        temperature <- nimNumeric(length=2, value=1) ## Length 2 is a hack. Only temperature[1] is used.
    },
    run = function() {
        ## browser()
        propValueVector <- generateProposalVector()
        if (temperPriors) {
            lpMHR <- my_setAndCalculateDiff$run(propValueVector) / temperature[1]
        } else {
            logPriorWeightOri   <- getLogProb(model, targetNodes)
            lpMHR               <- my_setAndCalculateDiff$run(propValueVector)
            diffLogPriorWeights <- getLogProb(model, targetNodes) - logPriorWeightOri
            lpMHR               <- (lpMHR - diffLogPriorWeights) / temperature[1] + diffLogPriorWeights
        }        
        jump <- my_decideAndJump$run(lpMHR, 0, 0, 0) ## will use lpMHR - 0
        ## nimPrint(temperature[1], " ", jump)
        if (adaptive) {
            adaptiveProcedure(jump)
        }
    },
    methods = list(
        setTemp = function(temp=double()) {
            temperature[1] <<- temp
        },
        generateProposalVector = function() {
            propValueVector <- rmnorm_chol(1, values(model,target), chol_propCov_scale, 0)  ## last argument specifies prec_param = FALSE
            returnType(double(1))
            return(propValueVector)
        },
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump) {
                timesAccepted <<- timesAccepted + 1
            }
            if(!adaptScaleOnly) {
                empirSamp[timesRan, 1:d] <<- values(model, target)
            }
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                adaptFactor <- my_calcAdaptationFactor$run(acceptanceRate)
                scale <<- scale * adaptFactor
                ## calculate empirical covariance, and adapt proposal covariance
                if(!adaptScaleOnly) {
                    gamma1 <- my_calcAdaptationFactor$gamma1
                    for(i in 1:d)     empirSamp[, i] <<- empirSamp[, i] - mean(empirSamp[, i])
                    empirCov <- (t(empirSamp) %*% empirSamp) / (timesRan-1)
                    propCov <<- propCov + gamma1 * (empirCov - propCov)
                    chol_propCov <<- chol(propCov)
                }
                chol_propCov_scale <<- chol_propCov * scale
                timesRan <<- 0
                timesAccepted <<- 0
            }
        },
        ## turnOffAdaptation = function() {
        ##     adaptive <<- FALSE
        ## },
        reset = function() {
            scale              <<- scaleOriginal
            propCov            <<- propCovOriginal
            chol_propCov       <<- chol(propCov)
            chol_propCov_scale <<- chol_propCov * scale
            timesRan           <<- 0
            timesAccepted      <<- 0
            timesAdapted       <<- 0
            my_calcAdaptationFactor$reset()
        }
    ), where = getLoadingNamespace()
)


####################################################################
### slice sampler (discrete or continuous) #########################
####################################################################

sampler_slice_tempered <- nimbleFunction(
    contains = sampler_APT,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive       <- control$adaptive
        adaptInterval  <- control$adaptInterval
        width          <- control$sliceWidth
        maxSteps       <- control$sliceMaxSteps
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        calcNodes      <- model$getDependencies(target)
        ## numeric value generation
        widthOriginal <- width
        timesRan      <- 0
        timesAdapted  <- 0
        sumJumps      <- 0
        discrete      <- model$isDiscrete(target)
        ## checks     
        if(length(targetAsScalar) > 1)     stop('cannot use slice sampler on more than one target node')
        ## initialise temperature
        temperature <- nimNumeric(length=2, value=1) ## Length 2 is a hack. Only first element is used.
    },
    run = function() { 
        ## nimPrint(temperature[1])
        u  <- getLogProb(model, calcNodes) / temperature[1] - rexp(1, 1)    # generate (log)-auxiliary variable: exp(u) ~ uniform(0, exp(lp))
        x0 <- model[[target]]    # create random interval (L,R), of width 'width', around current value of target
        L  <- x0 - runif(1, 0, 1) * width
        R  <- L + width
        maxStepsL <- floor(runif(1, 0, 1) * maxSteps)    # randomly allot (maxSteps-1) into maxStepsL and maxStepsR
        maxStepsR <- maxSteps - 1 - maxStepsL
        lp <- setAndCalculateTarget(L) / temperature[1]
        while(maxStepsL > 0 & !is.nan(lp) & lp >= u) {   # step L left until outside of slice (max maxStepsL steps)
            L <- L - width
            lp <- setAndCalculateTarget(L) / temperature[1]
            maxStepsL <- maxStepsL - 1
        }
        lp <- setAndCalculateTarget(R) / temperature[1]
        while(maxStepsR > 0 & !is.nan(lp) & lp >= u) {   # step R right until outside of slice (max maxStepsR steps)
            R <- R + width
            lp <- setAndCalculateTarget(R) / temperature[1]
            maxStepsR <- maxStepsR - 1
        }
        x1 <- L + runif(1, 0, 1) * (R - L)
        lp <- setAndCalculateTarget(x1) / temperature[1]
        while(is.nan(lp) | lp < u) {   # must be is.nan()
            if(x1 < x0) {
                L <- x1
            } else {
                R <- x1
            }
            x1 <- L + runif(1, 0, 1) * (R - L)           # sample uniformly from (L,R) until sample is inside of slice (with shrinkage)
            lp <- setAndCalculateTarget(x1) / temperature[1]
        }
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        jumpDist <- abs(x1 - x0)
        if(adaptive)     adaptiveProcedure(jumpDist)
    },
    methods = list(
        setTemp = function(temp = double()) {
            temperature[1] <<- temp
        },
        setAndCalculateTarget = function(value = double()) {
            if(discrete)     value <- floor(value)
            model[[target]] <<- value
            lp <- calculate(model, calcNodes)
            returnType(double())
            return(lp)
        },
        adaptiveProcedure = function(jumpDist = double()) {
            timesRan <<- timesRan + 1
            sumJumps <<- sumJumps + jumpDist   # cumulative (absolute) distance between consecutive values
            if(timesRan %% adaptInterval == 0) {
                adaptFactor <- (3/4) ^ timesAdapted
                meanJump <- sumJumps / timesRan
                width <<- width + (2*meanJump - width) * adaptFactor   # exponentially decaying adaptation of 'width' -> 2 * (avg. jump distance)
                timesAdapted <<- timesAdapted + 1
                timesRan <<- 0
                sumJumps <<- 0
            }
        },
        ## turnOffAdaptation = function() {
        ##     adaptive <<- FALSE
        ## },
        reset = function() {
            width        <<- widthOriginal
            timesRan     <<- 0
            timesAdapted <<- 0
            sumJumps     <<- 0
        }
    ), where = getLoadingNamespace()
)


#######################################################################################
### RW_multinomial sampler for multinomial distributions ##############################
#######################################################################################

sampler_RW_multinomial_tempered <- nimbleFunction( 
    contains = sampler_APT,
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        adaptive      <- control$adaptive
        adaptInterval <- control$adaptInterval
        useTempering  <- control$useTempering
        ## node list generation
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        targetAllNodes <- unique(model$expandNodeNames(target))
        calcNodes      <- model$getDependencies(target) 
        lTarget        <- length(targetAsScalar)
        Ntotal         <- sum(values(model,target))
        NOverL         <- Ntotal / lTarget
        ## numeric value generation
        Zeros             <- matrix(0, lTarget, lTarget)
        Ones              <- matrix(1, lTarget, lTarget)
        timesRan          <- Zeros
        AcceptRates       <- Zeros
        ScaleShifts       <- Zeros
        totalAdapted      <- Zeros
        timesAccepted     <- Zeros
        ENSwapMatrix      <- Ones
        ENSwapDeltaMatrix <- Ones
        RescaleThreshold  <- 0.2 * Ones
        lpProp  <- 0
        lpRev   <- 0
        Pi      <- pi 
        PiOver2 <- Pi / 2 ## Irrational number prevents recycling becoming degenerate
        u       <- runif(1, 0, Pi)
        ## nested function and function list definitions
        my_setAndCalculateDiff <- setAndCalculateDiff(model, target)
        my_decideAndJump       <- decideAndJump(model, mvSaved, calcNodes)
        ## checks
        ## if(model$getNodeDistribution(target) != 'dmulti')   stop('can only use RW_multinomial sampler for multinomial distributions')
        if(model$getDistribution(target) != 'dmulti') stop('can only use RW_multinomial sampler for multinomial distributions')
        if(length(targetAllNodes) > 1)                stop('cannot use RW_multinomial sampler on more than one target')
        if(adaptive & adaptInterval < 100)            stop('adaptInterval < 100 is not recommended for RW_multinomial sampler')
        ## initialise temperature
        temperature <- nimNumeric(length=2, value=1) ## Length 2 is a hack. Only first element is used.
    },
    run = function() {
        for(iFROM in 1:lTarget) {            
            for(iTO in 1:(lTarget-1)) {
                if(u > PiOver2) {                
                    iFrom <- iFROM
                    iTo   <- iTO
                    if (iFrom == iTo)
                        iTo <- lTarget
                    u <<- 2 * (u - PiOver2)   # recycle u
                } else {
                    iFrom <- iTO
                    iTo   <- iFROM
                    if (iFrom == iTo)
                        iFrom <- lTarget
                    u <<- 2 * (PiOver2 - u)   # recycle u
                }
                propValueVector <- generateProposalVector(iFrom, iTo)
                if (useTempering) {
                    lpMHR <- my_setAndCalculateDiff$run(propValueVector) / temperature[1] + lpRev - lpProp
                } else {
                    lpMHR <- my_setAndCalculateDiff$run(propValueVector) + lpRev - lpProp
                }
                jump <- my_decideAndJump$run(lpMHR, 0, 0, 0) ## returns lpMHR + 0 - 0 + 0
                if(adaptive)   adaptiveProcedure(jump=jump, iFrom=iFrom, iTo=iTo)
            }
        }
    },
    methods = list(
        setTemp = function(temp=double()) {
            temperature[1] <<- temp
        },
        generateProposalVector = function(iFrom = integer(), iTo = integer()) { 
            propVector <- values(model,target) 
            pSwap      <- min(1, max(1, ENSwapMatrix[iFrom,iTo]) / propVector[iFrom]) 
            nSwap      <- rbinom(n=1,   size=propVector[iFrom], prob=pSwap) 
            lpProp    <<- dbinom(nSwap, size=propVector[iFrom], prob=pSwap, log=TRUE) 
            propVector[iFrom] <- propVector[iFrom] - nSwap 
            propVector[iTo]   <- propVector[iTo]   + nSwap 
            pRevSwap   <- min(1, max(1, ENSwapMatrix[iTo,iFrom]) / (propVector[iTo] + nSwap)) 
            lpRev     <<- dbinom(nSwap, size=propVector[iTo], prob=pRevSwap, log=TRUE) 
            returnType(double(1)) 
            return(propVector) 
        },
        adaptiveProcedure = function(jump=logical(), iFrom=integer(), iTo=integer()) {
            NVector <- values(model,target) 
            timesRan[iFrom, iTo] <<- timesRan[iFrom, iTo] + 1
            if(jump)
                timesAccepted[iFrom, iTo] <<- timesAccepted[iFrom, iTo] + 1
            if (timesRan[iFrom, iTo] %% adaptInterval == 0) {
                totalAdapted[iFrom, iTo] <<- totalAdapted[iFrom, iTo] + 1
                accRate                   <- timesAccepted[iFrom, iTo] / timesRan[iFrom, iTo]
                AcceptRates[iFrom, iTo]  <<- accRate
                if (accRate > 0.5) {
                    ENSwapMatrix[iFrom, iTo] <<-
                        min(Ntotal,
                            ENSwapMatrix[iFrom,iTo] + ENSwapDeltaMatrix[iFrom, iTo] / totalAdapted[iFrom,iTo])
                } else {
                    ENSwapMatrix[iFrom, iTo] <<-
                        max(1,
                            ENSwapMatrix[iFrom,iTo] - ENSwapDeltaMatrix[iFrom,iTo] / totalAdapted[iFrom,iTo])
                } 
                if(accRate<RescaleThreshold[iFrom,iTo] | accRate>(1-RescaleThreshold[iFrom,iTo])) {
                    ## rescale iff ENSwapMatrix[iFrom, iTo] is not set to an upper or lower bound 
                    if (ENSwapMatrix[iFrom, iTo] > 1 & ENSwapMatrix[iFrom, iTo] < Ntotal) {
                        ScaleShifts[iFrom, iTo]       <<- ScaleShifts[iFrom, iTo] + 1 
                        ENSwapDeltaMatrix[iFrom, iTo] <<- min(NOverL, ENSwapDeltaMatrix[iFrom, iTo] * totalAdapted[iFrom,iTo] / 10)
                        ENSwapDeltaMatrix[iTo, iFrom] <<- ENSwapDeltaMatrix[iFrom, iTo] 
                        RescaleThreshold[iFrom,iTo]   <<- 0.2 * 0.95^ScaleShifts[iFrom, iTo]
                    }
                }
                ## lower Bound 
                if(ENSwapMatrix[iFrom, iTo] < 1)
                    ENSwapMatrix[iFrom, iTo] <<- 1                
                ## symmetry in ENSwapMatrix helps maintain good acceptance rates
                ENSwapMatrix[iTo,iFrom]   <<- ENSwapMatrix[iFrom,iTo]
                timesRan[iFrom, iTo]      <<- 0
                timesAccepted[iFrom, iTo] <<- 0
            }
        },
        ## turnOffAdaptation = function() {
        ##     adaptive <<- FALSE
        ## },
        reset = function() {
            timesRan          <<- Zeros
            AcceptRates       <<- Zeros
            ScaleShifts       <<- Zeros
            totalAdapted      <<- Zeros
            timesAccepted     <<- Zeros
            ENSwapMatrix      <<- Ones
            ENSwapDeltaMatrix <<- Ones
            RescaleThreshold  <<- 0.2 * Ones
        }
    ), where = getLoadingNamespace()
)


#' APT Sampling Algorithms
#'
#' Details of the APT sampling algorithms adapted from NIMBLE's MCMC samplers.
#'
#' @param model (uncompiled) model on which the APT algorithm is to be run
#' 
#' @param mvSaved \code{modelValues} object to be used to store MCMC samples
#' @param target node(s) on which the sampler will be used
#' @param control named list that controls the precise behavior of the sampler, with elements specific to \code{samplertype}.  The default values for control list are specified in the setup code of each sampling algorithm.  Descriptions of each sampling algorithm, and the possible customizations for each sampler (using the \code{control} argument) appear below.
#'
#' @section \code{sampler_APT}: base class for APT samplers
#'
#' When you write a new sampler for use in a NIMBLE MCMC with APT, you must include \code{contains = sampler_APT}.
#'
#' @section RW sampler:
#'
#' The RW sampler executes adaptive Metropolis-Hastings sampling with a normal proposal distribution (Metropolis, 1953), implementing the adaptation routine given in Shaby and Wells, 2011.  This sampler can be applied to any scalar continuous-valued stochastic node, and can optionally sample on a log scale.
#'
#' The RW sampler accepts the following control list elements:
#' \itemize{
#' \item logScale. A logical argument, specifying whether the sampler should operate on the log scale. (default = FALSE)
#' \item reflective. A logical argument, specifying whether the normal proposal distribution should reflect to stay within the range of the target distribution. (default = FALSE)
#' \item adaptive. A logical argument, specifying whether the sampler should adapt the scale (proposal standard deviation) throughout the course of MCMC execution to achieve a theoretically desirable acceptance rate. (default = TRUE)
#' \item adaptInterval. The interval on which to perform adaptation.  Every adaptInterval MCMC iterations (prior to thinning), the RW sampler will perform its adaptation procedure.  This updates the scale variable, based upon the sampler's achieved acceptance rate over the past adaptInterval iterations. (default = 200)
#' \item scale. The initial value of the normal proposal standard deviation.  If adaptive = FALSE, scale will never change. (default = 1)
#' \item temperPriors. Logical indicator determining if temperring should apply to prior likelihoods. Usually can be set to TRUE. But setting to FALSE can help avoid degeneracy issues for complex problems where bounded uniform priors have been transformed to other (e.g. logit) scales. 
#' }
#'
#' The RW sampler cannot be used with options log=TRUE and reflective=TRUE, i.e. it cannot do reflective sampling on a log scale.
#'
#' @section RW_block sampler:
#'
#' The RW_block sampler performs a simultaneous update of one or more model nodes, using an adaptive Metropolis-Hastings algorithm with a multivariate normal proposal distribution (Roberts and Sahu, 1997), implementing the adaptation routine given in Shaby and Wells, 2011.  This sampler may be applied to any set of continuous-valued model nodes, to any single continuous-valued multivariate model node, or to any combination thereof. \cr
#'
#' The RW_block sampler accepts the following control list elements:
#' \itemize{
#' \item adaptive. A logical argument, specifying whether the sampler should adapt the scale (a coefficient for the entire proposal covariance matrix) and propCov (the multivariate normal proposal covariance matrix) throughout the course of MCMC execution.  If only the scale should undergo adaptation, this argument should be specified as TRUE. (default = TRUE)
#' \item adaptScaleOnly. A logical argument, specifying whether adaption should be done only for scale (TRUE) or also for provCov (FALSE).  This argument is only relevant when adaptive = TRUE.  When adaptScaleOnly = FALSE, both scale and propCov undergo adaptation; the sampler tunes the scaling to achieve a theoretically good acceptance rate, and the proposal covariance to mimic that of the empirical samples.  When adaptScaleOnly = TRUE, only the proposal scale is adapted. (default = FALSE)
#' \item adaptInterval. The interval on which to perform adaptation.  Every adaptInterval MCMC iterations (prior to thinning), the RW_block sampler will perform its adaptation procedure, based on the past adaptInterval iterations. (default = 200)
#' \item scale. The initial value of the scalar multiplier for propCov.  If adaptive = FALSE, scale will never change. (default = 1)
#' \item propCov. The initial covariance matrix for the multivariate normal proposal distribution.  This element may be equal to the character string 'identity', in which case the identity matrix of the appropriate dimension will be used for the initial proposal covariance matrix. (default = 'identity')
#' \item temperPriors. Logical indicator determining if temperring should apply to prior likelihoods. Usually can be set to TRUE. But setting to FALSE can help avoid degeneracy issues for complex problems where bounded uniform priors have been transformed to other (e.g. logit) scales. 
#' }
#'
#' @section slice sampler:
#'
#' The slice sampler performs slice sampling of the scalar node to which it is applied (Neal, 2003).  This sampler can operate on either continuous-valued or discrete-valued scalar nodes.  The slice sampler performs a 'stepping out' procedure, in which the slice is iteratively expanded to the left or right by an amount sliceWidth.  This sampler is optionally adaptive, governed by a control list element, whereby the value of sliceWidth is adapted towards the observed absolute difference between successive samples.
#'
#' The slice sampler accepts the following control list elements:
#' \itemize{
#' \item adaptive. A logical argument, specifying whether the sampler will adapt the value of sliceWidth throughout the course of MCMC execution. (default = TRUE)
#' \item adaptInterval. The interval on which to perform adaptation. (default = 200)
#' \item width. The initial value of the width of each slice, and also the width of the expansion during the iterative 'stepping out' procedure. (default = 1)
#' \item maxSteps. The maximum number of expansions which may occur during the 'stepping out' procedure. (default = 100)
#' }
#'
#' @section RW_multinomial sampler:
#'
#' This sampler is designed for sampling multinomial target distributions.  The sampler performs a series of Metropolis-Hastings steps between pairs of groups.  Proposals are generated via a draw from a binomial distribution, whereafter the proposed number density is moved from one group to another group.  The acceptance or rejection of these proposals follows a standard Metropolis-Hastings procedure.  Probabilities for the random binomial proposals are adapted to a target acceptance rate of 0.5.
#'
#' The \code{RW_multinomial} sampler accepts the following control list elements:
#' \itemize{
#' \item adaptive.  A logical argument, specifying whether the sampler should adapt the binomial proposal probabilities throughout the course of MCMC execution. (default = TRUE)
#' \item adaptInterval.  The interval on which to perform adaptation.  A minimum value of 100 is required. (default = 200)
#' \item useTempering. A logical argument to optionally turn temporing off (i.e. assume all temperatures are 1) for this sampler. 
#' }
#'
#' @name samplers
#'
#' @aliases sampler sampler_RW_tempered sampler_RW_block_tempered sampler_RW_multinomial_tempered sampler_slice_tempered
#'
#' @examples
#' 
#' @seealso \code{\link{configureMCMC}} \code{\link{addSampler}} \code{\link{buildMCMC}} \code{\link{runMCMC}}
#'
#' @author David Pleydell (based on code from Daniel Turek)
#'
#' @references
#'
#' Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., Teller, A. H., and Teller, E. (1953). Equation of State Calculations by Fast Computing Machines. \emph{The Journal of Chemical Physics}, 21(6), 1087-1092.
#'
#' Neal, Radford M. (2003). Slice Sampling. \emph{The Annals of Statistics}, 31(3), 705-741.
#'
#' Roberts, G. O. and S. K. Sahu (1997). Updating Schemes, Correlation Structure, Blocking and Parameterization for the Gibbs Sampler. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 59(2), 291-317.
#'
#' Shaby, B. and M. Wells (2011). \emph{Exploring an Adaptive Metropolis Algorithm}. 2011-14. Department of Statistics, Duke University.
NULL



