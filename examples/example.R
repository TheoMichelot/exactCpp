
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
b <- c(-0.5,-1)
v <- c(1,0.5)

rates <- rep(0.1,2)

set.seed(1)

# simulate observations
sim <- simDataOU(mu=mu, b=b, v=v, rates=rates, map=NULL, duration=500)
obs <- as.matrix(sim[,c(1,2,4)])

# initial parameters
nbState <- length(b)
par0 <- list(m=c(t(mu)),b=-b,v=v)
rates0 <- rates

homog=list(mHomog=FALSE, bHomog=FALSE, vHomog=FALSE)
controls=list(kappa=0.2,lenmin=3,lenmax=6,thin=100,prUpdateMove=1,SDP=0.15)

allArgs <- setupMCMC(obs, par0, rates0, map=NULL, nbState=2, homog=homog, controls=controls, 
                     priorShape=c(12,4))
mod <- MCMCloop(allArgs)

estim <- read.table(mod$fileparams,header=TRUE)
rates <- read.table(mod$filerates,header=TRUE)

allPlots(nbState=nbState, fileparams=mod$fileparams, filerates=mod$filerates,
         truePar=c(as.vector(t(mu)),b,v))