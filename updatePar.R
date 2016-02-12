
#' Update movement parameters
updatePar <- function(allData,par,priorMean,priorSD,proposalSD,nbState,mHomog,bHomog,vHomog,obs)
{
    # enable references by "name" 
    colX <- 1; colY <- 2; colTime <- 3; colState <- 4; colHabitat <- 5; colJump <- 6; colBehav <- 7
    
    # pick out points that are informative about movement, and the states during movement
    whichInfo <- which(allData[,colTime]>min(obs[,colTime]))
    states <- allData[whichInfo-1,colState]
    
    # calculate changes in position and time
    preX <- allData[whichInfo-1,colX]
    dX <- allData[whichInfo,colX]-preX
    preY <- allData[whichInfo-1,colY]
    dY <- allData[whichInfo,colY]-preY
    deltaT <-  allData[whichInfo,colTime]-allData[whichInfo-1,colTime]
    
    moveStep <- updateMove(par,priorMean,priorSD,proposalSD,nbState,mHomog,bHomog,vHomog)
    
    # old & new likelihoods
    oldMove <- rawMove(par,state=states,deltaT=deltaT,x=preX,y=preY,nbState=nbState)
    newMove <- rawMove(moveStep$newPar,state=states,deltaT=deltaT,x=preX,y=preY,nbState=nbState)
    
    oldLogLX <- sum(dnorm(dX,mean=oldMove$emx,sd=oldMove$sdx,log=TRUE))
    oldLogLY <- sum(dnorm(dY,mean=oldMove$emy,sd=oldMove$sdy,log=TRUE))
    
    newLogLX <- sum(dnorm(dX,mean=newMove$emx,sd=newMove$sdx,log=TRUE))
    newLogLY <- sum(dnorm(dY,mean=newMove$emy,sd=newMove$sdy,log=TRUE))
    
    logHR <- moveStep$newLogPrior-moveStep$oldLogPrior+newLogLX+newLogLY-oldLogLX-oldLogLY
    
    if(runif(1)<exp(logHR)) # accept new movement parameters
        par <- moveStep$newPar
    
    return(par)
}