
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

# Are we working with the adaptative model?
adaptative <- FALSE

###############
## Read data ##
###############
# observations
# obs <- cbind(as.matrix(read.csv("fisherXYT.csv",header=FALSE,skip=3)),NA,NA,NA,NA)
# obs <- cbind(as.matrix(read.csv("~/Grive/software/simulateOU/obs.csv")),NA,NA,NA,NA)

obs <- cbind(as.matrix(read.table("~/Grive/data/simdata_myopic_review.txt",header=TRUE)),NA,NA,NA,NA)
colnames(obs) <- c("X","Y","Time","State","Habitat","Jump","Behav")
nbObs <- nrow(obs)

# enable references by "name" 
colX <- 1; colY <- 2; colTime <- 3; colState <- 4; colHabitat <- 5; colJump <- 6; colBehav <- 7

# map fisher
# map <- as.matrix(read.csv("habitat.csv",header=FALSE,skip=3))

# 1-state model map
map <- matrix(1,nrow=1,ncol=1)

# 2-state model map
# map <- matrix(1,nrow=12,ncol=12)
# map[7:12,] <- 2

nbHabitat <- length(unique(c(map))) # count habitat types
if(adaptative) {
    nbState <- nbHabitat
} else
    nbState <- 3

# Set up rest of data matrix
obs[,colHabitat] <- findRegion(obs[,colX],obs[,colY],map)
obs[,colState] <- obs[,colHabitat]
obs[,colJump] <- 0 # jump for data point is always 0
obs[,colBehav] <- 0 # behavioural states not known

########################
## Initial parameters ##
########################
mpar <- rep(c(0,0),nbState) # centers of attraction
# mpar <- c(4,5,8,7) # centers of attraction
bpar <- rep(1,nbState) # coefficients for matrix B
vpar <- rep(2,nbState) # coefficients for matrix Lambda
par <- c(mpar,bpar,vpar)

# Priors (on log scale for b and v)
priorMean <- n2w(par,nbState)

mPriorSD <- rep(c(10,10),nbState)
bPriorSD <- rep(10,nbState)
vPriorSD <- rep(10,nbState)
priorSD <- c(mPriorSD,bPriorSD,vPriorSD)

# MH proposals (on log scale for b and v)
mProposalSD <- rep(c(0.03,0.03),nbState)
bProposalSD <- rep(0.03,nbState)
vProposalSD <- rep(0.03,nbState)
proposalSD <- c(mProposalSD,bProposalSD,vProposalSD)

# for local update
SDP <- 0.15

# Initial lambda (non-diagonal elements, filled row-wise)
lambdapar <- rep(1,nbState*(nbState-1))

# prior beta parameters for lambdas
shape1 <- 4
shape2 <- 4
kappa <- 2

# homogeneity between states
mHomog <- TRUE
bHomog <- FALSE
vHomog <- FALSE

# probability of movement parameter update
prUpdateMove <- 1

####################
## Prepare output ##
####################
d <- format(Sys.time(), "%Y-%m-%d-%H%M")
fileparams <- paste("params", d, ".txt", sep = "")

cat(file=fileparams, mpar,bpar,vpar, "\n", append = TRUE)

filekappa <- paste("rates", d, ".txt", sep = "")
cat(file=filekappa, lambdapar, "\n", append = TRUE)

set.seed(1)

####################################
## Prepare set of actual swtiches ##
####################################
# indices of state switches
whichActual <- which(obs[-1,colState]!=obs[-nbObs,colState])+1
nbActual <- length(whichActual)

dt <- 0.1*min(diff(obs[,colTime]))

aSwitches <- cbind(X=obs[whichActual,colX],
                   Y=obs[whichActual,colY],
                   Time=obs[whichActual,colTime]-dt,
                   State=obs[whichActual,colState],
                   Habitat=obs[whichActual,colHabitat],
                   Jump=rep(1,nbActual),
                   Behav=rep(0,nbActual))

# Initialise
bk <- FALSE

# Controls
lenmin <- 3
lenmax <- 6
nbIter <- 5e5
thin <- 100
