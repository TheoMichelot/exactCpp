
source("setup.R")

for (iii in 1:nbIter)
{
    # length of considered interval of observations
    len <- sample(lenmin:lenmax,size=1)
    
    # indices of first and last selected obs
    point1 <- sample(1:(nbObs-len+1),size=1)
    point2 <- point1+len-1
    
    # selected observations
    subObs <- obs[point1:point2,] 
    
    # times and states of beginning and end
    Tbeg <- subObs[1,colTime]
    Tend <- subObs[len,colTime]
    
    ####################################
    ## Simulate movement and switches ##
    ####################################
    sim <- simMove_rcpp(subObs,par,kappa,lambdapar,nbState,map)
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
    
    nbActual <- nrow(aSwitches) # this nbActual has been updated to match new aSwitches
    
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
    
    # parameter update
    if(runif(1)<prUpdateMove) {
        par <- updatePar(allData,par,priorMean,priorSD,proposalSD,nbState,mHomog,bHomog,vHomog,obs)
        par2 <- updatePar_rcpp(allData,par,priorMean,priorSD,proposalSD,nbState,mHomog,bHomog,vHomog,obs)
        
        print(par)
        print(par2)
        stop()
    }
    
    # print parameters to file
    if(iii%%thin==0)
        cat(file=fileparams, round(par,6), "\n", append = TRUE)
    
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
    
    if((jorder-1)%in%indSwitchAll) # should be moveable, else can't do anything locally
        aSwitches <- localUpdate(allData,aSwitches,jorder,par,lambdapar,kappa,nbState,SDP,map)
    
    # print switching rates to file
    if(iii%%thin==0)
        cat(file=filekappa,lambdapar,"\n", append = TRUE)
    
    if(iii%%thin==0)
        cat("\nEnd iteration",iii,"\n")
    
    rm(allData)
}
