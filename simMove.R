
#' Simulation movement and switches
#' 
#' @param subObs Observations between Tbeg and Tend
#' @param par Vector of parameters (m,b,v)
#' @param kappa Maximum switching rate
#' @param lambdapar Vector of non-diagonal switching rates, row-wise
#' @param nbState Number of states
#' @param map Map of habitat (matrix)
simMove <- function(subObs,par,kappa,lambdapar,nbState,map)
{
    # enable references by "name" 
    colX <- 1; colY <- 2; colTime <- 3; colState <- 4; colHabitat <- 5; colJump <- 6; colBehav <- 7
    
    Tbeg <- subObs[1,colTime] # time of beginning
    Tend <- subObs[nrow(subObs),colTime] # time of end
    
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
    
    ####################################
    ## Simulate movement and switches ##
    ####################################  
    t <- 2
    bk <- FALSE
    while (t<=nrow(subData) & !bk) {
        
        # data point or potential jump? potential jump has Habitat NA
        if( is.na(subData[t,colHabitat]) ) {
            state <- subData[t-1,colState] # state prior to t
            deltaT <- subData[t,colTime]-subData[t-1,colTime] # time interval between t and its predecessor
            
            Xfrom <- subData[t-1,colX]
            Yfrom <- subData[t-1,colY]
            
            # compute distribution of next location (i.e. location t)
            move <- rawMove(par,state=state,deltaT=deltaT,x=Xfrom,y=Yfrom,nbState=nbState)
            
            # simulate location t
            subData[t,colX] <- rnorm(1,mean=Xfrom+move$emx,sd=move$sdx)
            subData[t,colY] <- rnorm(1,mean=Yfrom+move$emy,sd=move$sdy)
            
            # map back on to known habitat
            subData[t,colHabitat] <- findRegion(subData[t,colX],subData[t,colY],map)
            
            # prob of actual switch depends on rates and on region
            probs <- switchRate(subData[t-1,colState],subData[t,colHabitat],lambdapar)/kappa
            
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
            } else {
                # data point state is always same as previous  and no jump
                subData[t,colState] <- subData[t-1,colState]  
            }
        }
        
        t <- t+1
    } # End of loop over t>=2
    
    return(list(subData=subData,indObs=indObs,indSwitch=indSwitch,bk=bk))
}
