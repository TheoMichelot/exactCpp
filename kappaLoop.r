
source("basicSetup.r")
source("sharedSetup.r")

for (iii in 1:nbIter)
{
    # length of considered interval of observations
    len <- sample(lenmin:lenmax,size=1)
    
    # indices of first and last selected obs
    point1 <- sample(1:(ndata-len+1),size=1)
    point2 <- point1+len-1
    
    # selected observations
    subObs <- obs[point1:point2,] 
    
    # times and states of beginning and end
    Tbeg <- subObs[1,colTime]
    Tend <- subObs[len,colTime]
    Sbeg <- subObs[1,colState]
    Send <- subObs[len,colState]
    
    #################################
    ## Simulate potential switches ##
    #################################
    # number of potential switches
    nbSwitch <- rpois(1,(Tend-Tbeg)*kappa)
    
    # initialize simulated switches
    switches <- cbind(X=rep(NA,nbSwitch),
                      Y=rep(NA,nbSwitch),
                      Time=sort(runif(nbSwitch,Tbeg,Tend)), # Poisson process -> uniformly distributed
                      State=rep(NA,nbSwitch),
                      Habitat=rep(NA,nbSwitch),
                      Jump=rep(NA,nbSwitch),
                      Behav=rep(0,nbSwitch))
    
    # all simulated switches are in the first part, data are in second part
    subData <- rbind(switches,subObs)
    
    # ranks of data by time
    ranks <- rank(subData[,colTime])
    # indices of potential switches among observations
    indSwitch <- ranks[1:nbSwitch]
    # indices of observations among potential switches
    indObs <- ranks[(nbSwitch+1):nrow(subData)]
    
    # order data in time
    ord <- order(subData[,colTime])
    subData <- subData[ord,]
    
    if(!bk) {
        
        ####################################
        ## Simulate movement and switches ##
        ####################################  
        for (t in 2:nrow(subData)) {
            
            # data point or potential jump? potential jump has Habitat NA
            if( is.na(subData[t,colHabitat]) ) {
                state <- subData[t-1,colState] # state prior to t
                deltaT <- subData[t,colTime]-subData[t-1,colTime] # time interval between t and its predecessor
                
                Xfrom <- subData[t-1,colX]
                Yfrom <- subData[t-1,colY]
                
                # compute distribution of next location (i.e. location j)
                move <- rawMove(mpar,bpar,vpar,state=state,deltaT=deltaT,xx=Xfrom,yy=Yfrom)
                
                # simulate location j
                subData[t,colX] <- rnorm(1,mean=Xfrom+move$emx,sd=move$sdx)
                subData[t,colY] <- rnorm(1,mean=Yfrom+move$emy,sd=move$sdy)
                
                # Mapping back on to known habitat
                subData[t,colHabitat] <- findRegion(subData[t,colX],subData[t,colY],map)
                
                # Prob of actual switch depends on rates and on region
                probs <- switchRate(subData[t-1,colState],subData[t,colHabitat],subData[t,colTime],lambdapar)/kappa
                
                jumpNow <- runif(1)<sum(probs) # is there a jump?
                
                if(jumpNow)
                    subData[t,colState] <- successor(probs) # pick state of successor in case of jump
                else
                    subData[t,colState] <- subData[t-1,colState] # keep previous state if no jump
                
                subData[t,colJump] <- as.integer(jumpNow)
                
            } else {
                # we have a data point
                
                if( (t==nrow(subData) | subData[t,colBehav]==1) & 
                    (subData[t-1,colState] != subData[t,colState]) ) {
                    # state at Tend cannot be changed
                    
                    bk <- TRUE
                    break            
                }
                
                # data point state is always same as previous  and no jump
                subData[t,colState] <- subData[t-1,colState]  
            }
            
        } # End of loop over t>=2
    }
    
    if(!bk) {
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
        move <- rawMove(mpar,bpar,vpar,state=states,deltaT=deltaT,
                         xx=subData[which,colX],yy=subData[which,colY])
        
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
        aMove <- rawMove(mpar,bpar,vpar,state=tmpStates,deltaT=tmpDeltaT,
                         xx=aSubData[aIndObs[-1]-1,colX],yy=aSubData[aIndObs[-1]-1,colY])
        
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
        
        if(runif(1)<HR) {
            #######################
            ## Accept trajectory ##
            #######################
            
            # Update proposals - states etc have changed
            switches <- subData[c(0,indSwitch),,drop=FALSE]
            
            # Update Data states
            obs[obs[,colTime]>(Tbeg-0.1) & obs[,colTime]<(Tend+0.1),colState] <- subData[indObs,colState]
            
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
            
            acctraj <- acctraj+1         
        } # end of "if accept trajectory"      
    } # end of "if(!bk)" 
    
    nbActual <- nrow(aSwitches) # this nbActual has been updated to match new aSwitches
    
    ################################
    ## Update movement parameters ##
    ################################
    
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
    
    if(runif(1)<prUpdateMove) # start of "if" on updating movement
    {   
        # Update movement parameters
        # Pick out points that are informative about movement,
        # their predecessors, and the states during movement
        
        whichInfo <- which(allData[,colTime]>min(obs[,colTime]))
        states <- allData[whichInfo-1,colState]
        
        # Calculate changes in position and time
        
        preX <- allData[whichInfo-1,colX]
        dX <- allData[whichInfo,colX]-preX
        preY <- allData[whichInfo-1,colY]
        dY <- allData[whichInfo,colY]-preY
        deltaT <-  allData[whichInfo,colTime]-allData[whichInfo-1,colTime]
        
        moveStep <- updateMove(bpar, vpar, bHomog, bProposalSD, nbState, vHomog, vProposalSD, 
                                mpar, mProposalSD, mPriorMean, mPriorSD, bPriorMean, bPriorSD, 
                                vPriorMean, vPriorSD)
        
        # Old & new likelihoods
        
        oldMove <- rawMove(mpar,bpar,vpar,state=states,deltaT=deltaT,xx=preX,yy=preY)
        newMove <- rawMove(moveStep$mprime,moveStep$bprime,moveStep$vprime,state=states,
                            deltaT=deltaT,xx=preX,yy=preY)
        
        oldLogLX <- sum(dnorm(dX,mean=oldMove$emx,sd=oldMove$sdx,log=TRUE))
        oldLogLY <- sum(dnorm(dY,mean=oldMove$emy,sd=oldMove$sdy,log=TRUE))
        
        newLogLX <- sum(dnorm(dX,mean=newMove$emx,sd=newMove$sdx,log=TRUE))
        newLogLY <- sum(dnorm(dY,mean=newMove$emy,sd=newMove$sdy,log=TRUE))
        
        logHR <- moveStep$newLogPrior-moveStep$oldLogPrior+newLogLX+newLogLY-oldLogLX-oldLogLY
        
        if(runif(1)<exp(logHR)) {
            # Accept movement parameters
            accmove <- accmove+1
            mpar <- moveStep$mprime
            bpar <- moveStep$bprime
            vpar <- moveStep$vprime
        }
    } # end of "if" on updating movement
    
    if(iii%%thin==0)
        cat(file=fileparams, round(mpar,6), round(bpar,6), round(vpar,6), "\n", append = TRUE)
    
    # update jump rate k
    if(!bk & nbActual>0)
        lambdapar <- updateRate(allData, indSwitchAll, kappa, shape1, shape2)
    
    bk <- FALSE  
    
    #########################
    ## Insert local update ##
    #########################
    
    # Local update to predecessor to obs j
    j <- sample(2:nrow(obs),size=1)
    jorder <- indObsAll[j]
    
    if((jorder-1)%in%indSwitchAll) # Should be moveable
    {
        # Want to perturb (T,X,Y)[jorder-1]
        
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
        
        oldMove2 <- rawMove(mpar,bpar,vpar,state=S2,deltaT=T1-T2,xx=X2,yy=Y2)
        
        oldLogLX2 <- dnorm(X1,mean=X2+oldMove2$emx,sd=oldMove2$sdx,log=TRUE)
        oldLogLY2 <- dnorm(Y1,mean=Y2+oldMove2$emy,sd=oldMove2$sdy,log=TRUE)
        
        oldMove1 <- rawMove(mpar,bpar,vpar,state=S1,deltaT=Tj-T1,xx=X1,yy=Y1) 
        
        oldLogLX1 <- dnorm(Xj,mean=X1+oldMove1$emx,sd=oldMove1$sdx,log=TRUE)
        oldLogLY1 <- dnorm(Yj,mean=Y1+oldMove1$emy,sd=oldMove1$sdy,log=TRUE)
        
        oldLL <- oldLogLX1+oldLogLY1+oldLogLX2+oldLogLY2  
        
        # Propose new location
        
        Xnew <- rnorm(1,mean=X1,sd=SDP)
        Ynew <- rnorm(1,mean=Y1,sd=SDP)
        Tnew <- runif(1,min=T2,max=Tj)
        Hnew <- findRegion(Xnew,Ynew,map)
        
        # Effect of habitat
        
        rate1 <- switchRate(S2,H1,T1,lambdapar)
        rate1[S2] <- kappa-sum(rate1)
        rate.new <- switchRate(S2,Hnew,Tnew,lambdapar)
        rate.new[S2] <- kappa-sum(rate.new)
        newHfactor <- rate.new[S1]/rate1[S1]
        
        # New log-likelihood
        
        newMove2 <- rawMove(mpar,bpar,vpar,state=S2,deltaT=Tnew-T2,xx=X2,yy=Y2)
        
        newLogLX2 <- dnorm(Xnew,mean=X2+newMove2$emx,sd=newMove2$sdx,log=TRUE)
        newLogLY2 <- dnorm(Ynew,mean=Y2+newMove2$emy,sd=newMove2$sdy,log=TRUE)
        
        newMove1 <- rawMove(mpar,bpar,vpar,state=S1,deltaT=Tj-Tnew,xx=Xnew,yy=Ynew) 
        
        newLogLX1 <- dnorm(Xj,mean=Xnew+newMove1$emx,sd=newMove1$sdx,log=TRUE)
        newLogLY1 <- dnorm(Yj,mean=Ynew+newMove1$emy,sd=newMove1$sdy,log=TRUE)
        
        newLL <- newLogLX1+newLogLY1+newLogLX2+newLogLY2  
        
        logHR <- newLL-oldLL+log(newHfactor)
        # Calc needs to be done in this order, to avoid special case: exp(newLL-oldLL)=Inf, Hfactor=0
        
        if(runif(1)<exp(logHR)) { # accept local update (T2 to Tj)
            accloc <- accloc+1
        
            # index of switch "jorder-1" in aSwitches
            prev <- which(aSwitches[,colTime]==allData[jorder-1,colTime])
            
            aSwitches[prev,colTime] <- Tnew
            aSwitches[prev,colX] <- Xnew
            aSwitches[prev,colY] <- Ynew
        } # end accept
    } # end previous point imputed
    # else can't do anything locally
    
    if(iii%%thin==0) {
        cat(file=filekappa,lambdapar,"\n", append = TRUE)
        cat(file=filedata, obs[1:ndata,4], "\n", append = TRUE)
    }
    
    # Diagnostics
    if(iii%%thin==0) {
        cat("\nEnd iteration",iii,"\n")
        cat(file=filediag, acctraj, accmove, accloc, "\n", append = TRUE)
    }
    
    rm(allData)
    
} # end of for loop niter
