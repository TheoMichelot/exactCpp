
setwd("~/git/moveMCMC")
source("minitabGeneric.r")
source("updateRate_fisher.R")
source("updateRate_unconstr.R")
source("utilities.R")
source("allPlots.R")

library(Rcpp)
library(RcppArmadillo)
library(gtools) # for rdirichlet in homogeneous case
sourceCpp("simMove.cpp")
sourceCpp("moveLike.cpp")
sourceCpp("updatePar.cpp")
sourceCpp("localUpdate.cpp")

source("simDataOU.R")
source("setupMCMC.R")
source("MCMCloop.R")

# define map
map <- matrix(1,nrow=12,ncol=12)
map[7:12,] <- 2

# movement parameters
mu <- matrix(c(5,6,5,6), ncol=2, byrow=TRUE)
b <- c(-2,-0.5)
v <- c(3,0.5)

rates <- rep(0.5,2)

set.seed(1)

# simulate observations
sim <- simDataOU(mu=mu, b=b, v=v, rates=rates, map=map)
obs <- as.matrix(sim[,c(1,2,4)])

# initial parameters
nbState <- length(b)
mpar <- rep(c(0,0),nbState) # centers of attraction
bpar <- rep(1,nbState) # coefficients for matrix B
vpar <- rep(2,nbState) # coefficients for matrix Lambda
par0 <- list(m=mpar,b=bpar,v=vpar)
rates0 <- rep(1,nbState*(nbState-1))

homog=list(mHomog=TRUE, bHomog=FALSE, vHomog=FALSE)

allArgs <- setupMCMC(obs, par0, rates0, map=map, homog=homog)
mod <- MCMCloop(allArgs)

allPlots(nbState=nbState, fileparams=mod$fileparams, filerates=mod$filerates,
         truePar=c(as.vector(t(mu)),b,v))