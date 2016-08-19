
setwd("~/git/moveMCMC")
source("minitabGeneric.r")
source("updateRate_fisher.R")
source("updateRate_unconstr.R")
source("utilities.R")

library(Rcpp)
library(RcppArmadillo)
library(gtools) # for rdirichlet in homogeneous case
sourceCpp("simMove.cpp")
sourceCpp("moveLike.cpp")
sourceCpp("updatePar.cpp")
sourceCpp("localUpdate.cpp")

source("setupMCMC.R")
source("MCMCloop.R")

obs <- as.matrix(read.csv("~/Grive/software/simulation/obs.csv"))
map <- matrix(1,nrow=12,ncol=12)
map[7:12,] <- 2

mpar <- rep(c(0,0),nbState) # centers of attraction
bpar <- rep(1,nbState) # coefficients for matrix B
vpar <- rep(2,nbState) # coefficients for matrix Lambda
par0 <- list(m=mpar,b=bpar,v=vpar)

rates0 <- rep(1,nbState*(nbState-1))

set.seed(1)
allArgs <- setupMCMC(obs, par0, rates0, map=map)
mod <- MCMCloop(allArgs)
