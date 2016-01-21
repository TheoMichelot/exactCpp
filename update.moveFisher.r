
#' Update movement
updateMove <- function(b.paras, v.paras, homog.b, b.proposal.sd, nstate, homog.v, v.proposal.sd, m.paras,
                       m.proposal.sd, m.prior.mean, m.prior.sd, b.prior.mean, b.prior.sd, v.prior.mean,
                       v.prior.sd)
{
    b.real <- log(b.paras)
    v.real <- log(v.paras)
    
    if(homogb) {
        b.prime <- b.real+rnorm(1,0,b.proposal.sd)
        b1 <- b.real[1]
        b2 <- b.prime[1]
    } else {
        b.prime <- b.real+rnorm(nstate,0,b.proposal.sd)
        b1 <- b.real
        b2 <- b.prime
    }
    
    if(homogv) {
        v.prime <- v.real+rnorm(1,0,v.proposal.sd) 
        v1 <- v.real[1]
        v2 <- v.prime[1]
    } else {
        v.prime <- rnorm(nstate,v.real,v.proposal.sd)
        v1 <- v.real
        v2 <- v.prime
    } 
    
    m.prime <- rnorm(2,m.paras,m.proposal.sd)
    
    # Old log prior
    OLP <- sum(dnorm(m.paras,m.prior.mean,m.prior.sd,log=TRUE))+
        sum(dnorm(b1,b.prior.mean,b.prior.sd,log=TRUE))+
        sum(dnorm(v1,v.prior.mean,v.prior.sd,log=TRUE))
    
    # New log prior
    NLP=sum(dnorm(m.prime,m.prior.mean,m.prior.sd,log=TRUE))+
        sum(dnorm(b2,b.prior.mean,b.prior.sd,log=TRUE))+
        sum(dnorm(v2,v.prior.mean,v.prior.sd,log=TRUE))
    
    return(list(m.prime=m.prime,b.prime=exp(b.prime),v.prime=exp(v.prime),old.logprior=OLP,new.logprior=NLP))  
}