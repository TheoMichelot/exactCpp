
#' Movement likelihood
#' 
#' @param subData Observations and potential switches during considered interval
#' @param indObs Vector of indices of observations in subData
#' @param indSwitch Vector of indices of switches in subData
#' @param aSwitches Actual switches
#' @param nbStates Number of states
moveLike <- function(subData,indObs,indSwitch,par,aSwitches,nbState)
{
    # enable references by "name" 
    colX <- 1; colY <- 2; colTime <- 3; colState <- 4; colHabitat <- 5; colJump <- 6; colBehav <- 7
 
    Tbeg <- subData[1,colTime]
    Tend <- subData[nrow(subData),colTime]
       
    ################
    ## Likelihood ##
    ################
    ## 1. Likelihood with potential switches
    
    # intervals between data points (except 1st one) and previous switches
    deltaT <- subData[indObs[-1],colTime] - subData[indObs[-1]-1,colTime]
    
    # states at each switch preceding a data point
    states <- subData[indObs[-1]-1,colState] 
    # index of each switch preceding a data point
    which <- indObs[-1]-1
    
    # compute distributions of next locations
    # (i.e. distribution of actual data points, to deduce likelihood)
    move <- rawMove(par,state=states,deltaT=deltaT,
                    x=subData[which,colX],y=subData[which,colY],
                    nbState=nbState)
    
    RWX <- dnorm(subData[indObs[-1],colX],
                 mean = subData[which,colX]+move$emx,
                 sd = move$sdx,
                 log = TRUE)
    
    RWY <- dnorm(subData[indObs[-1],colY],
                 mean = subData[which,colY]+move$emy,
                 sd = move$sdy,
                 log = TRUE)
    
    ## 2. Likelihood with actual switches
    
    # indices of "actual switches" happening between Tbeg and Tend
    whichActual <- which(aSwitches[,colTime]<Tend & aSwitches[,colTime]>Tbeg)
    
    # number of those switches
    nbActual <- length(whichActual)
    
    # Actual switches and observations between Tbeg and Tend
    aSubData <- rbind(aSwitches[whichActual,],subObs)
    
    # ranks of data by time
    aRanks <- rank(aSubData[,colTime])
    # indices of actual switches among observations
    aIndSwitch <- aRanks[1:nbActual]
    # indices of observations among actual switches
    aIndObs <- aRanks[(nbActual+1):nrow(aSubData)]
    
    # order data in time
    ord <- order(aSubData[,colTime])
    aSubData <- aSubData[ord,]
    
    # intervals between obs and predecessors
    tmpDeltaT <- aSubData[aIndObs[-1],colTime] - aSubData[aIndObs[-1]-1,colTime]
    # states of predecessors
    tmpStates <- aSubData[aIndObs[-1]-1,colState]
    
    # compute distributions of next data point
    aMove <- rawMove(par,state=tmpStates,deltaT=tmpDeltaT,
                     x=aSubData[aIndObs[-1]-1,colX],y=aSubData[aIndObs[-1]-1,colY],
                     nbState=nbState)
    
    oldLikeX <- dnorm(aSubData[aIndObs[-1],colX],
                      mean = aSubData[aIndObs[-1]-1,colX]+aMove$emx,
                      sd = aMove$sdx,
                      log=TRUE)
    
    oldLikeY <- dnorm(aSubData[aIndObs[-1],colY],
                      mean = aSubData[aIndObs[-1]-1,colY]+aMove$emy,
                      sd = aMove$sdy,
                      log=TRUE)
    
    # Hastings ratio
    HR <- exp(sum(RWX+RWY)-sum(oldLikeX+oldLikeY))
    
    return(HR=HR)
}