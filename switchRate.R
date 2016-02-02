
#' Switching rates
#' 
#' Vector of switching rates when in state 'state' and in habitat 'habitat'.
#' All zeros if state==habitat; otherwise, all zeros except for habitat
#' (only possible to switch to state corresponding to current habitat).
#' 
#' @param state State
#' @param habitat Habitat
#' @param par Lambda parameters
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