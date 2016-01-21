
switchRate <- function(state,habitat,par,ns=3,nr=3)
{
    rate <- array(0,c(ns,ns,nr))
    z <- matrix(0,ns,ns)
    off <- lower.tri(z)|upper.tri(z)
    cz <- col(z)[off]
    rz <- row(z)[off]
    rate[cbind(rz,cz,cz)] <- par
    return(rate[state,,habitat])
}