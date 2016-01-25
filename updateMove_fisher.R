
#' Update movement parameters
#' 
#' Suggest new movement parameters from the proposal distributions.
updateMove <- function(bpar, vpar, bHomog, bProposalSD, nbState, vHomog, vProposalSD, 
                       mpar, mProposalSD, mPriorMean, mPriorSD, bPriorMean, bPriorSD, 
                       vPriorMean, vPriorSD)
{
    breal <- log(bpar)
    vreal <- log(vpar)
    
    # allow for homogeneous and non-homogeneous cases
    if(bHomog) {
        bprime <- breal+rnorm(1,0,bProposalSD)
        b1 <- breal[1]
        b2 <- bprime[1]
    } else {
        bprime <- breal+rnorm(nbState,0,bProposalSD)
        b1 <- breal
        b2 <- bprime
    }
    
    if(vHomog) {
        vprime <- vreal+rnorm(1,0,vProposalSD) 
        v1 <- vreal[1]
        v2 <- vprime[1]
    } else {
        vprime <- rnorm(nbState,vreal,vProposalSD)
        v1 <- vreal
        v2 <- vprime
    } 
    
    mprime <- rnorm(2,mpar,mProposalSD)
    
    # Old log prior
    OLP <- sum(dnorm(mpar,mPriorMean,mPriorSD,log=TRUE))+
        sum(dnorm(b1,bPriorMean,bPriorSD,log=TRUE))+
        sum(dnorm(v1,vPriorMean,vPriorSD,log=TRUE))
    
    # New log prior
    NLP=sum(dnorm(mprime,mPriorMean,mPriorSD,log=TRUE))+
        sum(dnorm(b2,bPriorMean,bPriorSD,log=TRUE))+
        sum(dnorm(v2,vPriorMean,vPriorSD,log=TRUE))
    
    return(list(mprime=mprime,bprime=exp(bprime),vprime=exp(vprime),oldLogPrior=OLP,newLogPrior=NLP))  
}