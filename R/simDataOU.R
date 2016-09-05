
#' Simulate data from a switching OU process
#' 
#' @param mu Matrix of coordinates of centres of attraction (one row for each state)
#' @param b Vector of b parameters.
#' @param v Vector of v parameters.
#' @param rates Vector of switching rates, 
#' e.g (lambda_12, lambda_13, lambda_21, lambda_23, lambda_31, lambda_32).
#' @param mty Vector indicating the type of model used in each state (1 for Brownian motion,
#' 2 for location OU).
#' @param map Map of habitats.
#' @param interval Interval between observations.
#' @param duration Duration between first and last observations.
#' @param write If TRUE, the simulated data are written into the file obs.csv (default: FALSE).
#' @param showPlot If TRUE, the simulated data are plotted (default).
#' @param points If TRUE, the simulated data are plotted as dots (default), 
#' otherwise only segments are plotted.
simDataOU <- function(mu, b, v, rates, mty, map=NULL, interval=1, duration=500, 
                      write=FALSE, showPlot=TRUE, points=TRUE)
{
    nbState <- length(b)
    
    adapt <- TRUE # adaptative case?
    if(is.null(map)) {
        map <- matrix(1,1,1)
        adapt <- FALSE
    }
    
    # matrix of switching rates
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
    data$state[1] <- sample(1:nbState,size=1)
    if(adapt)
        data$habitat[1] <- data$state[1]
    data$x[1] <- mu[data$state[1],1]
    data$y[1] <- mu[data$state[1],2]
    
    # is the data point an observation? (or a potential switch)
    isObs <- (times %in% obsTimes)
    
    # loop over simulated observations
    for(t in 2:nbData) {
        if(t%%500==0) 
            cat(t/nbData*100,"%\n",sep="")
        
        # time since last observation
        dt <- times[t] - times[t-1]
        
        # state
        state <- data$state[t-1]
        
        if(mty[state]==1) { # Brownian motion
            meanx <- data$x[t-1]
            meany <- data$y[t-1]
            sd <- v[state]
        } else { # Ornstein-Uhlenbeck
            meanx <- mu[state,1] + exp(b[state]*dt)*(data$x[t-1]-mu[state,1])
            meany <- mu[state,2] + exp(b[state]*dt)*(data$y[t-1]-mu[state,2])
            sd <- sqrt(v[state]*(1-exp(2*b[state]*dt)))
        }
        
        data$x[t] <- rnorm(1,meanx,sd)
        data$y[t] <- rnorm(1,meany,sd)
        
        if(adapt)
            data$habitat[t] <- findRegion(data$x[t],data$y[t],map)
        
        # if potential switch
        if(!isObs[t]) {
            if(adapt)
                prJumpNow <- lambda[data$state[t-1],data$habitat[t]]/kappa
            else
                prJumpNow <- sum(lambda[data$state[t-1],])/kappa
            
            if(runif(1)<prJumpNow) {
                if(adapt)
                    data$state[t] <- data$habitat[t]
                else
                    data$state[t] <- sample(1:nbState,size=1,prob=lambda[data$state[t-1],])
            } else
                data$state[t] <- data$state[t-1]
        } else
            data$state[t] <- data$state[t-1]
    }
    
    # only keep obervations (remove switches)
    obs <- data[isObs,]
    
    if(showPlot) {
        par(mfrow=c(1,1))
        pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")
        
        if(adapt)
            image(seq(0.5,nrow(map)-0.5,by=1),seq(0.5,ncol(map)-0.5,by=1),map,xlab="x",ylab="y",
                  col=pal[1:nbState])
        else {
            xmin <- min(obs[,1])
            xmax <- max(obs[,1])
            xmid <- (xmin+xmax)/2
            ymin <- min(obs[,2])
            ymax <- max(obs[,2])
            ymid <- (ymin+ymax)/2
            l <- max(xmax-xmin,ymax-ymin)/2
            
            if(points) {
                plot(obs[1,1],obs[1,2],pch=21,bg=pal[obs[1,3]],cex=0.8,
                     xlim=c(xmid-l,xmid+l),ylim=c(ymid-l,ymid+l),xlab="x",ylab="y")                
            } else {
                plot(NA,xlim=c(xmid-l,xmid+l),ylim=c(ymid-l,ymid+l),xlab="x",ylab="y")    
            }
        }
        
        if(points) {
            points(obs[,1],obs[,2],pch=21,bg=pal[obs[,3]],cex=0.8,type="o")
        } else {
            for(t in 2:nrow(obs))
                segments(obs[t-1,1],obs[t-1,2],obs[t,1],obs[t,2],col=pal[obs[t-1,3]])
        }
    }
    
    if(write)
        write.csv(obs[,c(1,2,4)],file="obs.csv",row.names=FALSE)        
    
    return(obs)
}


