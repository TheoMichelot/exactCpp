
#' Main loop of MCMC algorithm
#' 
#' @param allArgs List of all arguments, as returned by \code{\link{setupMCMC}}:
#' \item{obs}{Matrix of observations}
#' \item{map}{Map of habitats}
#' \item{nbState}{Number of states}
#' \item{nbIter}{Number of iterations}
#' \item{adapt}{TRUE if adaptative model, FALSE otherwise}
#' \item{par0}{Initial movement parameters}
#' \item{rates0}{Initial switching rates}
#' \item{priorMean}{Mean of priors for movement parameters}
#' \item{priorSD}{SD of priors for movement parameters}
#' \item{priorShape}{Shape of prior for switching rates}
#' \item{proposalSD}{SD of proposal distribution for movement parameters}
#' \item{controls}{See args of \code{\link{setupMCMC}}}
#' \item{homog}{See args of \code{\link{setupMCMC}}}
#' \item{fileparams}{Output file for sampled movement parameters}
#' \item{filerates}{Output file for sampled switching rates}
#' \item{aSwitches}{Matrix of actual switches (initial value)}
#' \item{nbActual}{Number of actual switches (initial value)}
#' \item{whichActual}{Indices of actual switches (initial value)}
MCMCloop <- function(allArgs)
{
    # unpack arguments
    obs <- allArgs$obs
    nbObs <- nrow(obs)
    map <- allArgs$map
    nbState <- allArgs$nbState
    nbIter <- allArgs$nbIter
    adapt <- allArgs$adapt
    par0 <- allArgs$par0
    rates0 <- allArgs$rates0
    priorMean <- allArgs$priorMean
    priorSD <- allArgs$priorSD
    priorShape <- allArgs$priorShape
    proposalSD <- allArgs$proposalSD
    controls <- allArgs$controls
    homog <- allArgs$homog
    fileparams <- allArgs$fileparams
    filerates <- allArgs$filerates
    aSwitches <- allArgs$aSwitches
    nbActual <- allArgs$nbActual
    whichActual <- allArgs$whichActual
    
    t <- Sys.time()
    accTraj <- 0
    accTrajAll <- rep(NA,nbIter/controls$thin)
    accPar <- 0
    accParAll <- rep(NA,nbIter/controls$thin)
    
    # enable references by "name" 
    colX <- 1; colY <- 2; colTime <- 3; colState <- 4; colHabitat <- 5; colJump <- 6; colBehav <- 7
    
    # initialise movement parameters and rates
    par <- allArgs$par0
    rates <- allArgs$rates0
 
    bk <- FALSE
       
    # MCMC loop
    for(iter in 1:nbIter) {
        # length of considered interval of observations
        len <- sample(controls$lenmin:controls$lenmax,size=1)
        
        # indices of first and last selected obs
        point1 <- sample(1:(nbObs-len+1),size=1)
        point2 <- point1+len-1
        
        # selected observations
        subObs <- obs[point1:point2,]
        
        # times of beginning and end
        Tbeg <- subObs[1,colTime]
        Tend <- subObs[len,colTime]
        
        ####################################
        ## Simulate movement and switches ##
        ####################################
        sim <- simMove_rcpp(subObs,par,controls$kappa,rates,nbState,map,adapt)
        subData <- sim[[1]]
        indObs <- sim[[2]]
        indObs <- indObs[order(indObs)]
        indSwitch <- sim[[3]]
        indSwitch <- indSwitch[order(indSwitch)]
        bk <- sim[[4]]
        
        
        if(!bk) {
            # compute the likelihood of the trajectory
            HR <- moveLike_rcpp(subData,indObs-1,indSwitch-1,par,aSwitches,nbState)
            
            if(runif(1)<HR) {
                
                accTraj <- accTraj + 1
                
                #######################
                ## Accept trajectory ##
                #######################
                # Update proposals - states etc have changed
                switches <- subData[c(0,indSwitch),,drop=FALSE]
                
                # Update Data states
                obs[obs[,colTime]>=Tbeg & obs[,colTime]<=Tend,colState] <- subData[indObs,colState]
                
                # data's jump do not need updating it is always 0
                
                # there is nothing outside the interval of Tbeg:Tend so that we can use the new Tswitch over old ATswitch
                if(!any(aSwitches[,colTime]>Tend | aSwitches[,colTime]<Tbeg)) { # Rare!
                    aSwitches <- switches
                } else {
                    # there IS something outside the interval we have to keep those 
                    # and replace the stuff in the interval and then resort
                    
                    whichOut <- which(aSwitches[,colTime]>Tend | aSwitches[,colTime]<Tbeg)
                    aSwitches <- rbind(aSwitches[whichOut,],switches)     
                    
                    # sort switches by time
                    orderA <- order(aSwitches[,colTime])
                    aSwitches <- aSwitches[orderA,,drop=FALSE]
                }
            }
        }
        
        # update nbActual for new aSwitches
        nbActual <- nrow(aSwitches)
        
        ################################
        ## Update movement parameters ##
        ################################
        # actual switches and observations
        allData <- rbind(aSwitches,obs)
        
        # ranks of data by time
        ranksAll <- rank(allData[,colTime])
        # indices of actual switches among observations
        indSwitchAll <- ranksAll[1:nbActual]
        # indices of observations among actual switches
        indObsAll <- ranksAll[(nbActual+1):nrow(allData)]
        
        # order data in time
        ord <- order(allData[,colTime])
        allData <- allData[ord,]
        
        parCopy <- par
        
        # parameter update
        if(runif(1)<controls$prUpdateMove)
            par <- updatePar_rcpp(allData,par,priorMean,priorSD,proposalSD,nbState,
                                  homog$mHomog,homog$bHomog,homog$vHomog,obs)
        
        if(!all(par==parCopy))
            accPar <- accPar + 1
        
        # print parameters to file
        if(iter%%controls$thin==0)
            cat(file=fileparams, round(par,6), "\n", append = TRUE)
        
        # update jump rate k
        if(!bk & nbActual>0) {
            if(adapt)
                rates <- updateRate(allData, indSwitchAll, controls$kappa, priorShape[1], priorShape[2], nbState)
            else
                rates <- updateRate_unconstr(allData, indSwitchAll, controls$kappa, priorShape[1], priorShape[2], nbState)
        }
        
        # print switching rates to file
        if(iter%%controls$thin==0)
            cat(file=filerates,rates,"\n", append = TRUE)
        
        bk <- FALSE
        
        #########################
        ## Insert local update ##
        #########################
        # Local update to predecessor to obs j
        j <- sample(2:nbObs,size=1)
        jorder <- indObsAll[j]
        
        if((jorder-1)%in%indSwitchAll) # should be moveable, else can't do anything locally
            aSwitches <- localUpdate_rcpp(allData,aSwitches,jorder-1,par,rates,controls$kappa,nbState,
                                          controls$SDP,map,adapt)
        
        if(iter%%controls$thin==0) {
            cat("End iteration ",iter,"/",nbIter," -- ",Sys.time()-t,"\n",sep="")
            accTrajAll[iter/controls$thin] <- accTraj/iter*100
            accParAll[iter/controls$thin] <- accPar/iter*100
        }
        
        rm(allData)
    }
    
    return(list(fileparams=fileparams,
                filerates=filerates,
                accTrajAll=accTrajAll,
                accParAll=accParAll))
}
