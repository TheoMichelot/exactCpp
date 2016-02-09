
source("switchRate.R")
source("minitabGeneric.r")
source("updateMove_fisher.R")
source("updateRate_fisher.R")
source("rawMove.R")
source("utilities.R")

###############
## Read data ##
###############
# observations
obs <- cbind(as.matrix(read.table("fisherXYT.csv",sep=",",header=FALSE,skip=3)),NA,NA,NA,NA)
colnames(obs) <- c("X","Y","Time","State","Habitat","Jump","Behav")
nbObs <- nrow(obs)

# enable references by "name" 
colX <- 1
colY <- 2
colTime <- 3
colState <- 4
colHabitat <- 5
colJump <- 6
colBehav <- 7

# map
map <- as.matrix(read.table("habitat.csv",sep=",",header=FALSE,skip=3))

nbHabitat <- length(unique(c(map))) # count habitat types
nbState <- nbHabitat

# Set up rest of data matrix
obs[,colHabitat] <- findRegion(obs[,colX],obs[,colY],map)
obs[,colState] <- findRegion(obs[,colX],obs[,colY],map)
obs[,colJump] <- 0 # jump for data point is always 0
obs[,colBehav] <- 0 # behavioural states not known

########################
## Initial parameters ##
########################
# mpar <- rep(c(9,3.5),nbState) # centers of attraction
mpar <- rep(c(9,3.5),1) # centers of attraction
bpar <- rep(1,nbState) # coefficients for matrix B
vpar <- rep(10,nbState) # coefficients for matrix Lambda
par <- c(mpar,bpar,vpar)

# Priors (on log scale for b, v)
# priorMean <- n2w(par,nbState)
mPriorMean <- mpar
bPriorMean <- log(bpar)
vPriorMean <- log(vpar)

# mPriorSD <- rep(c(0.7,0.7),nbState)
mPriorSD <- rep(c(0.7,0.7),1)
bPriorSD <- rep(2.0,nbState)
vPriorSD <- rep(2.0,nbState)
priorSD <- c(mPriorSD,bPriorSD,vPriorSD)

# MH proposals (on log scale for b, v)
# mProposalSD <- rep(c(0.01,0.01),nbState)
mProposalSD <- rep(c(0.01,0.01),1)
bProposalSD <- rep(0.2,nbState)
vProposalSD <- rep(0.1,nbState)
proposalSD <- c(mProposalSD,bProposalSD,vProposalSD)

# for "local update"
SDP <- 0.15

# Initial lambda (non-diagonal elements, filled row-wise)
lambdapar <- rep(3,6)

# prior beta parameters for lambdas
shape1 <- 12
shape2 <- 4
kappa <- 4

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

cat(file=fileparams, "mux","muy","b1","b2","b3","v1","v2","v3", "\n")
cat(file=fileparams, mpar,bpar,vpar, "\n", append = TRUE)

filekappa <- paste("rates", d, ".txt", sep = "")
cat(file=filekappa, "L12","L13","L21","L23","L31","L32", "\n")
cat(file=filekappa, lambdapar, "\n", append = TRUE)

set.seed(1)

####################################
## Prepare set of actual swtiches ##
####################################
# indices of state switches
whichActual <- which(obs[-1,colState]!=obs[-nbObs,colState])+1
nbActual <- length(whichActual)

dt <- 0.1*min(diff(obs[,colTime])) ## ?

aSwitches <- cbind(X=obs[whichActual,colX],
                   Y=obs[whichActual,colY],
                   Time=obs[whichActual,colTime]-dt,
                   State=obs[whichActual,colState],
                   Habitat=obs[whichActual,colHabitat],
                   Jump=rep(1,nbActual),
                   Behav=rep(0,nbActual))

# Initialise
bk <- FALSE
acctraj <- 0
accmove <- 0
accloc <- 0

# Controls
lenmin <- 3
lenmax <- 6
nbIter <- 1000
thin <- 100
