
#' Simulate data from a switching OU process
#' 
#' @param mu Matrix of coordinates of centres of attraction (one row for each state)
#' @param b Vector of b parameters.
#' @param v Vector of v parameters.
#' @param rates Vector of switching rates, e.g (lambda_12, lambda_13, lambda_21, lambda_23, lambda_31, lambda_32).
#' @param interval Interval between observations.
#' @param duration Duration between first and last observations.
#' @param write If TRUE, the simulated data are written into the file obs.csv (default: FALSE).
simDataOU <- function(mu, b, v, rates, map, interval=1, duration=500, write=FALSE)
{
    nbState <- length(b)
    
    ############################
    ## Prepare the parameters ##
    ############################
    B <- array(NA,c(2,2,nbState))
    for(state in 1:nbState)
        B[,,state] <- b[state]*diag(2)
    
    V <- array(NA,c(2,2,nbState))
    for(state in 1:nbState)
        V[,,state] <- v[state]*diag(2)
    
    lambda <- diag(nbState)
    lambda[!lambda] <- rates
    lambda <- t(lambda) # filled column-wise
    diag(lambda) <- 0
    
    kappa <- max(rowSums(lambda))
    
    # number of observations to simulate
    nbObs <- duration/interval
    
    # times of observations to simulate
    obsTimes <- cumsum(rep(interval,nbObs))
    
    # number of switches
    nbSwitch <- rpois(1,duration*kappa)
    # times of potential switches
    switchTimes <- runif(nbSwitch,0,duration)
    
    # times of all data (obs + switches)
    times <- c(obsTimes,switchTimes)
    times <- times[order(times)]
    nbData <- length(times)
    
    data <- data.frame(x=rep(NA,nbData),
                       y=rep(NA,nbData),
                       state=rep(NA,nbData),
                       time=times,
                       habitat=rep(NA,nbData))
    
    # initial location
    data$x[1] <- 4
    data$y[1] <- 6
    data$habitat[1] <- findRegion(data$x[1],data$y[1],map)
    data$state[1] <- data$habitat[1]
    
    # is the data point an observation? (or a potential switch)
    isObs <- (times %in% obsTimes)
    
    # loop over simulated observations
    for(t in 2:nbData) {
        if(t%%500==0) 
            cat(t/nbData*100,"%\n",sep="")
        
        # time since last observation
        dt <- times[t] - times[t-1]
        
        meanx <- mu[data$state[t-1],1] + exp(b[data$state[t-1]]*dt)*(data$x[t-1]-mu[data$state[t-1],1])
        meany <- mu[data$state[t-1],2] + exp(b[data$state[t-1]]*dt)*(data$y[t-1]-mu[data$state[t-1],2])
        sd <- sqrt(v[data$state[t-1]]*(1-exp(2*b[data$state[t-1]]*dt)))
        data$x[t] <- rnorm(1,meanx,sd)
        data$y[t] <- rnorm(1,meany,sd)
        
        data$habitat[t] <- findRegion(data$x[t],data$y[t],map)
        
        # if potential switch
        if(!isObs[t]) {
            prJumpNow <- lambda[data$state[t-1],data$habitat[t]]/kappa
            
            if(runif(1)<prJumpNow)
                data$state[t] <- data$habitat[t]
            else
                data$state[t] <- data$state[t-1]
        } else
            data$state[t] <- data$state[t-1]
    }
    
    # only keep obervations (remove switches)
    obs <- data[isObs,]
    
    # plot observations on map (fisher-type)
    par(mfrow=c(1,1))
    pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")
    # image(seq(0.5,16.5,by=1),seq(0.5,10.5,by=1),map,xlab="x",ylab="y",col=pal[1:nbState]) # fisher map
    image(seq(0.5,11.5,by=1),seq(0.5,11.5,by=1),map,xlab="x",ylab="y",col=pal[1:nbState]) # 2-state map
    points(obs[,1],obs[,2],pch=21,bg=pal[obs[,3]],cex=0.8,type="o")

    if(write)
        write.csv(obs[,c(1,2,4)],file="obs.csv",row.names=FALSE)        

    return(obs)
}


