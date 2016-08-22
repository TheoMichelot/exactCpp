
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

obs <- as.matrix(read.csv("~/Grive/work_other/HMMworkshop_MosselBay/reedbuck/reedbuck_data.csv")[,c(2,3)])
obs <- cbind(obs,1:nrow(obs))

set.seed(1)

# initial parameters
nbState <- 2
mpar <- rep(c(30,-30),nbState) # centers of attraction
bpar <- rep(0.5,nbState) # coefficients for matrix B
vpar <- rep(0.5,nbState) # coefficients for matrix Lambda
par0 <- list(m=mpar,b=bpar,v=vpar)
rates0 <- rep(1,nbState*(nbState-1))

homog=list(mHomog=TRUE, bHomog=FALSE, vHomog=FALSE)

controls=list(kappa=1,lenmin=2,lenmax=4,thin=100,prUpdateMove=0.5,SDP=0.15)

allArgs <- setupMCMC(obs, par0, rates0, nbState=nbState, homog=homog, controls=controls, nbIter=2e6)
mod <- MCMCloop(allArgs)

estim <- read.table(mod$fileparams,header=TRUE)
allPlots(nbState=nbState, fileparams=mod$fileparams, filerates=mod$filerates)
