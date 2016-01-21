# Adaptive i.e. can switch only to one state matching current region
# Diagonals are zero, unlike actual generator

struct.fn <- function(par,ns,nr)
{
  rate <- array(0,c(ns,ns,nr))
  z <- matrix(0,ns,ns)
  off <- lower.tri(z)|upper.tri(z)
  cz <- col(z)[off]
  rz <- row(z)[off]
  rate[cbind(rz,cz,cz)] <- par
  rate
}