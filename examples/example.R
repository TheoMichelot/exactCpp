
setwd("~/git/moveMCMC")
source("R/minitabGeneric.r")
source("R/updateRate_fisher.R")
source("R/updateRate_unconstr.R")
source("R/utilities.R")
source("R/allPlots.R")

library(Rcpp)
library(RcppArmadillo)
library(gtools) # for rdirichlet in homogeneous case
sourceCpp("source/simMove.cpp")
sourceCpp("source/moveLike.cpp")
sourceCpp("source/updatePar.cpp")
sourceCpp("source/localUpdate.cpp")

source("R/simDataOU.R")
source("R/setupMCMC.R")
source("R/MCMCloop.R")

# define map
map <- matrix(1,nrow=12,ncol=12)
map[7:12,] <- 2

# movement parameters
mu <- matrix(c(0,0,10,0), ncol=2, byrow=TRUE)
b <- c(-2,-0.5)
v <- c(3,1)

rates <- rep(0.1,2)

set.seed(1)

# simulate observations
sim <- simDataOU(mu=mu, b=b, v=v, rates=rates, map=NULL)
obs <- as.matrix(sim[,c(1,2,4)])

# initial parameters
nbState <- length(b)
mpar <- rep(c(0,0),nbState) # centers of attraction
bpar <- rep(1,nbState) # coefficients for matrix B
vpar <- rep(2,nbState) # coefficients for matrix Lambda
par0 <- list(m=mpar,b=bpar,v=vpar)
rates0 <- rep(0.1,nbState*(nbState-1))

homog=list(mHomog=FALSE, bHomog=FALSE, vHomog=FALSE)

allArgs <- setupMCMC(obs, par0, rates0, map=NULL, homog=homog, nbState=nbState)
mod <- MCMCloop(allArgs)

estim <- read.table(mod$fileparams,header=TRUE)

allPlots(nbState=nbState, fileparams=mod$fileparams, filerates=mod$filerates,
         truePar=c(as.vector(t(mu)),b,v))