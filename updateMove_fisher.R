
#' Update movement parameters
#' 
#' @param par Vector of movement parameters (m,b,v)
#' @param priorMean Vector of prior means
#' @param priorSD Vector of prior standard deviations
#' @param proposalSD Vector of standard deviations of the proposal distributions
#' @param nbState Number of states
#' @param mHomog TRUE if m is the same in all states, FALSE otherwise (default)
#' @param bHomog TRUE if b is the same in all states, FALSE otherwise (default)
#' @param vHomog TRUE if v is the same in all states, FALSE otherwise (default)
#' 
#' Suggest new movement parameters from the proposal distributions.
updateMove <- function(par,priorMean,priorSD,proposalSD,nbState,mHomog=FALSE,bHomog=FALSE,vHomog=FALSE)
{
    par <- n2w(par,nbState)
    # unpack the vector of parameters
    m <- par[1:(2*nbState)]
    b <- par[(2*nbState+1):(3*nbState)]
    v <- par[(3*nbState+1):(4*nbState)]
    
    # unpack prior
    mPriorMean <- priorMean[1:(2*nbState)]
    bPriorMean <- priorMean[(2*nbState+1):(3*nbState)]
    vPriorMean <- priorMean[(3*nbState+1):(4*nbState)]
    mPriorSD <- priorSD[1:(2*nbState)]
    bPriorSD <- priorSD[(2*nbState+1):(3*nbState)]
    vPriorSD <- priorSD[(3*nbState+1):(4*nbState)]
    
    # unpack proposal
    mProposalSD <- proposalSD[1:(2*nbState)]
    bProposalSD <- proposalSD[(2*nbState+1):(3*nbState)]
    vProposalSD <- proposalSD[(3*nbState+1):(4*nbState)]
    
    # allow for homogeneous and non-homogeneous cases
    if(mHomog) {
        mprime <- m + rep(rnorm(2,0,mProposalSD),nbState)
        m1 <- m[1:2]
        m2 <- mprime[1:2]
    } else {
        mprime <- m + rnorm(2*nbState,0,mProposalSD)
        m1 <- m
        m2 <- mprime
    }
    
    if(bHomog) {
        bprime <- b + rnorm(1,0,bProposalSD)
        b1 <- b[1]
        b2 <- bprime[1]
    } else {
        bprime <- b + rnorm(nbState,0,bProposalSD)
        b1 <- b
        b2 <- bprime
    }
    
    if(vHomog) {
        vprime <- v + rnorm(1,0,vProposalSD) 
        v1 <- v[1]
        v2 <- vprime[1]
    } else {
        vprime <- v + rnorm(nbState,0,vProposalSD)
        v1 <- v
        v2 <- vprime
    }

    # Old log prior
    oldLogPrior <- sum(dnorm(m1,mPriorMean,mPriorSD,log=TRUE)) + 
        sum(dnorm(b1,bPriorMean,bPriorSD,log=TRUE)) + 
        sum(dnorm(v1,vPriorMean,vPriorSD,log=TRUE))
    
    # New log prior
    newLogPrior <- sum(dnorm(m2,mPriorMean,mPriorSD,log=TRUE)) + 
        sum(dnorm(b2,bPriorMean,bPriorSD,log=TRUE)) + 
        sum(dnorm(v2,vPriorMean,vPriorSD,log=TRUE))
    
    newPar <- w2n(c(mprime,bprime,vprime),nbState)
    
    return(list(newPar=newPar,oldLogPrior=oldLogPrior,newLogPrior=newLogPrior))  
}