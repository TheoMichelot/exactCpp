# Create initial set of actual switches

nbActual <- 0
ATswitch <- c()
ASswitch <- c()
AXswitch <- c()
AYswitch <- c()
AHswitch <- c()

dt <- 0.1*min(diff(obs[,colTime]))

for(kkk in 2:ndata) { 
    if( (obs[kkk,4]!=obs[kkk-1,4])) {
        ATswitch <- c(ATswitch,obs[kkk,3]-dt)
        ASswitch <- c(ASswitch,obs[kkk,4])
        AXswitch <- c(AXswitch,obs[kkk,1])
        AYswitch <- c(AYswitch,obs[kkk,2])
        AHswitch <- c(AHswitch,obs[kkk,5])
        nbActual <- nbActual+1
    }
}

Ajump <- matrix(c(1),1,nbActual)
Ajump <- c(Ajump)
Abehav <- rep(0,length(Ajump))

aSwitches <- cbind(X=AXswitch,
                   Y=AYswitch,
                   Time=ATswitch,
                   State=ASswitch,
                   Habitat=AHswitch,
                   Jump=Ajump,
                   Behav=Abehav)

# Define functions

successor <- function(jumpLambda)
    sample(1:nbState,size=1,prob=jumpLambda)

source("minitabGeneric.r") 


# Initialise

oldLikeX <- -Inf
oldLikeY <- -Inf
bk <- FALSE
acctraj <- 0
accmove <- 0
accloc <- 0

# Specific output files

filedata <- paste("datastate", d, ".txt", sep = "")
cat(file=filedata,paste("data",1:ndata,sep=""), "\n")
cat(file=filedata, obs[,4], "\n", append = TRUE)

filediag <- paste("diagnostic", d, ".txt", sep = "")
cat(file=filediag, "acctraj", "accmove", "accloc","\n")
cat(file=filediag, acctraj, accmove, accloc,"\n", append = TRUE)

# Controls

lenmin <- 3
lenmax <- 6
nbIter <- 1000
thin <- 100
