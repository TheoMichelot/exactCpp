
#' Local update
#' 
#' We want to perturb (T,X,Y)[jorder-1]
localUpdate <- function(allData,aSwitches,jorder,par,lambdapar,kappa,nbState,SDP,map)
{
    # enable references by "name" 
    colX <- 1; colY <- 2; colTime <- 3; colState <- 4; colHabitat <- 5; colJump <- 6; colBehav <- 7
    
    T2 <- allData[jorder-2,colTime]
    X2 <- allData[jorder-2,colX]
    Y2 <- allData[jorder-2,colY]
    S2 <- allData[jorder-2,colState]
    
    T1 <- allData[jorder-1,colTime]
    X1 <- allData[jorder-1,colX]
    Y1 <- allData[jorder-1,colY]
    S1 <- allData[jorder-1,colState]
    
    H1 <- findRegion(X1,Y1,map)
    
    Tj <- allData[jorder,colTime]
    Xj <- allData[jorder,colX]
    Yj <- allData[jorder,colY]
    
    # Old log-likelihood
    
    oldMove2 <- rawMove(par,state=S2,deltaT=T1-T2,x=X2,y=Y2,nbState=nbState)
    
    oldLogLX2 <- dnorm(X1,mean=X2+oldMove2$emx,sd=oldMove2$sdx,log=TRUE)
    oldLogLY2 <- dnorm(Y1,mean=Y2+oldMove2$emy,sd=oldMove2$sdy,log=TRUE)
    
    oldMove1 <- rawMove(par,state=S1,deltaT=Tj-T1,x=X1,y=Y1,nbState=nbState) 
    
    oldLogLX1 <- dnorm(Xj,mean=X1+oldMove1$emx,sd=oldMove1$sdx,log=TRUE)
    oldLogLY1 <- dnorm(Yj,mean=Y1+oldMove1$emy,sd=oldMove1$sdy,log=TRUE)
    
    oldLL <- oldLogLX1+oldLogLY1+oldLogLX2+oldLogLY2  
    
    # Propose new location
    Xnew <- rnorm(1,mean=X1,sd=SDP)
    Ynew <- rnorm(1,mean=Y1,sd=SDP)
    Tnew <- runif(1,min=T2,max=Tj)
    Hnew <- findRegion(Xnew,Ynew,map)
    
    # Effect of habitat
    rate1 <- switchRate(S2,H1,lambdapar)
    rate1[S2] <- kappa-sum(rate1)
    rateNew <- switchRate(S2,Hnew,lambdapar)
    rateNew[S2] <- kappa-sum(rateNew)
    newHfactor <- rateNew[S1]/rate1[S1]
    
    # New log-likelihood
    newMove2 <- rawMove(par,state=S2,deltaT=Tnew-T2,x=X2,y=Y2,nbState=nbState)
    
    newLogLX2 <- dnorm(Xnew,mean=X2+newMove2$emx,sd=newMove2$sdx,log=TRUE)
    newLogLY2 <- dnorm(Ynew,mean=Y2+newMove2$emy,sd=newMove2$sdy,log=TRUE)
    
    newMove1 <- rawMove(par,state=S1,deltaT=Tj-Tnew,x=Xnew,y=Ynew,nbState=nbState) 
    
    newLogLX1 <- dnorm(Xj,mean=Xnew+newMove1$emx,sd=newMove1$sdx,log=TRUE)
    newLogLY1 <- dnorm(Yj,mean=Ynew+newMove1$emy,sd=newMove1$sdy,log=TRUE)
    
    newLL <- newLogLX1+newLogLY1+newLogLX2+newLogLY2  
    
    logHR <- newLL-oldLL+log(newHfactor)
    # Calc needs to be done in this order, to avoid special case: exp(newLL-oldLL)=Inf, Hfactor=0
    
    if(runif(1)<exp(logHR)) { # accept local update (T2 to Tj)
        # index of switch "jorder-1" in aSwitches
        prev <- which(aSwitches[,colTime]==allData[jorder-1,colTime])
        
        aSwitches[prev,colTime] <- Tnew
        aSwitches[prev,colX] <- Xnew
        aSwitches[prev,colY] <- Ynew
    } # end accept
    
    return(aSwitches)
}