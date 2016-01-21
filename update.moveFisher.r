update.move.fn <- function (b.paras, v.paras, homog.b, b.proposal.sd, nstate, homog.v, v.proposal.sd, m.paras, m.proposal.sd, m.prior.mean, m.prior.sd, b.prior.mean, b.prior.sd, v.prior.mean, v.prior.sd) 
{
  b.real=log(b.paras)
  v.real=log(v.paras)
  
  if(homog.b) 
    b.prime=b.real+rnorm(1,0,b.proposal.sd) else
      b.prime=b.real+rnorm(nstate,0,b.proposal.sd)
  if(homog.v)
    v.prime=v.real+rnorm(1,0,v.proposal.sd) else 
      v.prime=rnorm(nstate,v.real,v.proposal.sd)
  
  m.prime=rnorm(2,m.paras,m.proposal.sd)
  
  # Old & new priors
  
  OLP=sum(dnorm(m.paras,m.prior.mean,m.prior.sd,log=TRUE))+
    ifelse(homog.b,sum(dnorm(b.real[1],b.prior.mean,b.prior.sd,log=TRUE)),
           sum(dnorm(b.real,b.prior.mean,b.prior.sd,log=TRUE)))+
    ifelse(homog.v,sum(dnorm(v.real[1],v.prior.mean,v.prior.sd,log=TRUE)), 
           sum(dnorm(v.real,v.prior.mean,v.prior.sd,log=TRUE)))
  NLP=sum(dnorm(m.prime,m.prior.mean,m.prior.sd,log=TRUE))+
    ifelse(homog.b,sum(dnorm(b.prime[1],b.prior.mean,b.prior.sd,log=TRUE)),
           sum(dnorm(b.prime,b.prior.mean,b.prior.sd,log=TRUE)))+
    ifelse(homog.v,sum(dnorm(v.prime[1],v.prior.mean,v.prior.sd,log=TRUE)), 
           sum(dnorm(v.prime,v.prior.mean,v.prior.sd,log=TRUE)))
  
  list(m.prime=m.prime,b.prime=exp(b.prime),v.prime=exp(v.prime),old.logprior=OLP,new.logprior=NLP)
}

# Some calculations rather inelegant, permitting state-homogeneous special case for comparison