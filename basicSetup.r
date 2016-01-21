# Read data

obs <- cbind(as.matrix(read.table("fisherXYT.csv",sep=",",header=FALSE,skip=3)),NA,NA,NA,NA)
colnames(obs) <- c("X","Y","Time","State","Habitat","Jump","Behav")
ndata <- nrow(obs)

# Enable references by "name" 
colX <- 1 # "X"
colY <- 2 # "Y"
colTime <- 3 # "Time"
colState <- 4 # "State"
colHabitat <- 5 # "Habitat"
colJump <- 6 # "Jump"
colBehav <- 7 # "Behav"

# Read map of habitat
map <- as.matrix(read.table("habitat.csv",sep=",",header=FALSE,skip=3))
nbHabitat <- 3
nbState <- nbHabitat

#' Find region
#' 
#' Identify the region of (x,y) on 'map'.
#' 
#' @param x 'x' coordinate
#' @param y 'y' coordinate
#' @param map A matrix of habitats
findRegion <- function(x,y,map)
{
    nx <- nrow(map) # height of the map
    ny <- ncol(map) # width of the map
    
    xi <- ceiling(x)
    yi <- ceiling(y)
    
    # keep coordinates within the bounds of the map
    xi[xi<=0] <- 1
    xi[xi>nx] <- nx
    yi[yi<=0] <- 1
    yi[yi>ny] <- ny
    
    return(map[cbind(xi,yi)])
}

# Set up rest of data matrix
obs[,colHabitat] <- findRegion(obs[,colX],obs[,colY],map)
obs[,colState] <- findRegion(obs[,colX],obs[,colY],map)
obs[,colJump] <- 0 # jump for data point is always 0
obs[,colBehav] <- 0 # behavioural states not known

# Initial OU parameters
m <- c(9,3.5) # center of attraction
b <- 1.0 # factor for matrix B
v <- 10 # factor for matrix Lambda
mpar <- m
bpar <- rep(b,nbState) # all same to be uninformative AND for correct results with simpler models
vpar <- rep(v,nbState) #  " 

# Priors (on log scale for b, v)
mPriorMean <- mpar
bPriorMean <- log(bpar)
vPriorMean <- log(vpar)

mPriorSD <- rep(0.7,2)
bPriorSD <- rep(2.0,nbState)
vPriorSD <- rep(2.0,nbState)

# MH proposals (on log scale for b, v)
mProposalSD <- rep(0.01,2)
bProposalSD <- rep(0.2,nbState)
vProposalSD <- rep(0.1,nbState)

SDP <- 0.15 # for "local update"

# Other model-specific functions
source("adaptive.r") # to create struct.fn
source("setup.lambda.adaptive.r") 
source("update.moveFisher.r")
source("update.rateFisher.r")
source("mv.fnOUcommon.r")

#' Truncated beta
#' 
#' Random generation for the truncated Beta distribution
#' 
#' @param n Number of values to generate
#' @param shape1 Non-negative parameter
#' @param shape2 Non-negative parameter
#' @param low Lower limit of the truncated distribution
rtruncbeta <- function(n=1,shape1=1,shape2=1,low=0.5)
{
    p.min <- pbeta(low,shape1,shape2)
    p <- p.min+runif(n)*(1-p.min)
    
    return(qbeta(p,shape1,shape2))
}

# Initial lambda as a function of current behaviour and location
initL=3.0
lambda21 <- initL
lambda31 <- initL
lambda12 <- initL
lambda32 <- initL
lambda13 <- initL
lambda23 <- initL
lambdapar <- c(lambda12,lambda13,lambda21,lambda23,lambda31,lambda32)

# prior beta parameters for lambdas
shape1 <- 12
shape2 <- 4
kappa <- 4

# Output
d <- format(Sys.time(), "%Y-%m-%d-%H%M")
fileparams <- paste("params", d, ".txt", sep = "")

cat(file=fileparams, "mux","muy","b1","b2","b3","v1","v2","v3", "\n")
cat(file=fileparams, mpar,bpar,vpar, "\n", append = TRUE)

filekappa <- paste("rates", d, ".txt", sep = "")
cat(file=filekappa, "L12","L13","L21","L23","L31","L32", "\n")
cat(file=filekappa, lambdapar, "\n", append = TRUE)

# Homogeneity between states - for model variation inc DIC calculation
bHomog <- FALSE
vHomog <- FALSE

prUpdateMove <- 1

set.seed(1)
