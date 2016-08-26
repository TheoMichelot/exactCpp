
#' Setup data and parameters for MCMC algorithm
#' 
#' @param obs Matrix or data frame of observations, with columns x, y, and time.
#' @param par0 List of initial values for the movement parameters. Must have three elements: m, b, and v,
#' each of length the number of states.
#' @param rates0 Vector of initial values for the switching rates. Must be of length 
#' number of states * (number of states - 1).
#' @param homog List of mHomog, bHomog, and vHomog, which each indicate whether the corresponding movement
#' parameter is homogeneous across states (TRUE) or not (FALSE -- default).
#' @param priorMean List of means of priors for the movement parameters. Must have three elements: m, b, and v,
#' each of length the number of states.
#' @param priorSD List of standard deviations of priors for the movement parameters. Must have three elements: 
#' m, b, and v, each of length the number of states.
#' @param proposalSD List of standard deviations of proposals for the movement parameters. Must have three elements: 
#' m, b, and v, each of length the number of states.
#' @param shape Prior shapes for beta distribution of switching rates.
#' @param nbIter Number of iterations of the MCMC algorithm.
#' @param map Map of habitats, if adaptative model.
#' @param nbState Number of states.
#' @param controls List of control parameters: 
#' \item{kappa}{Kappa parameter from Blackwell et al (2015): upper bound for the switching rate out of any state;}
#' \item{lenmin}{Minimum length of updated interval;}
#' \item{lenmax}{Maximum length of updated interval;}
#' \item{thin}{Thinning factor;}
#' \item{prUpdateMove}{Probability of updating movement parameters at each iteration.}
setupMCMC <- function(obs, par0, rates0, homog=list(mHomog=FALSE,bHomog=FALSE,vHomog=FALSE), 
                      priorMean=NULL, priorSD=NULL, proposalSD=NULL, priorShape=c(4,4), nbIter=5e5, 
                      map=NULL, nbState=NULL, 
                      controls=list(kappa=2,lenmin=3,lenmax=6,thin=100,prUpdateMove=1,SDP=0.15))
{
    par <- c(par0$m,par0$b,par0$v)
    
    # Are we working with the adaptative model?
    if(is.null(map)) {
        if(is.null(nbState))
            stop("'nbState' needs to be specified if no map is given.")
        adapt <- FALSE
        
        map <- matrix(1,nrow=1,ncol=1)
    } else
        adapt <- TRUE
    
    if(any(par0$b<0))
        stop("Initial values for b should be positive (we really estimate -b).")
    
    ###############
    ## Read data ##
    ###############
    # observations
    obs <- cbind(obs,NA,NA,NA,NA)
    colnames(obs) <- c("X","Y","Time","State","Habitat","Jump","Behav")
    nbObs <- nrow(obs)
    
    # enable references by "name" 
    colX <- 1; colY <- 2; colTime <- 3; colState <- 4; colHabitat <- 5; colJump <- 6; colBehav <- 7
    
    nbHabitat <- length(unique(c(map))) # count habitat types
    if(adapt)
        nbState <- nbHabitat
    
    # Set up rest of data matrix
    obs[,colHabitat] <- findRegion(obs[,colX],obs[,colY],map)
    obs[,colState] <- obs[,colHabitat]
    obs[,colJump] <- 0 # jump for data point is always 0
    obs[,colBehav] <- 0 # behavioural states not known
    
    #######################
    ## Set up parameters ##
    #######################
    # Priors (on log scale for b and v)
    if(is.null(priorMean))
        priorMean <- n2w(par,nbState)
    else
        priorMean <- n2w(c(priorMean$m,priorMean$b,priorMean$v),nbState)
    
    if(is.null(priorSD)) {
        mPriorSD <- rep(c(10,10),nbState)
        bPriorSD <- rep(10,nbState)
        vPriorSD <- rep(10,nbState)
        priorSD <- c(mPriorSD,bPriorSD,vPriorSD) 
    } else
        priorSD <- c(priorSD$m,priorSD$b,priorSD$v) 
    
    # MH proposals (on log scale for b and v)
    if(is.null(proposalSD)) {
        mProposalSD <- rep(c(0.03,0.03),nbState)
        bProposalSD <- rep(0.1,nbState)
        vProposalSD <- rep(0.1,nbState)
        proposalSD <- c(mProposalSD,bProposalSD,vProposalSD)
    } else
        proposalSD <- c(proposalSD$m,proposalSD$b,proposalSD$v)
    
    # Initial lambda (non-diagonal elements, filled row-wise)
    if(is.null(rates0))
        rates0 <- rep(1,nbState*(nbState-1))
    
    ####################
    ## Prepare output ##
    ####################
    d <- format(Sys.time(), "%Y-%m-%d-%H%M")
    
    fileparams <- paste("params", d, ".txt", sep = "")
    
    # initialize file (to make sure it's empty)
    cat(file=fileparams, "", sep="")
    # Header
    for(state in 1:nbState)
        cat(file=fileparams, "mux", state, " muy", state, " ", append=TRUE, sep = "")
    for(state in 1:nbState)
        cat(file=fileparams, "b", state, " ", append=TRUE, sep = "")
    for(state in 1:nbState)
        cat(file=fileparams, "v", state, " ", append=TRUE, sep = "")
    # First row
    cat(file=fileparams, "\n", par, "\n", append=TRUE)
    
    filerates <- paste("rates", d, ".txt", sep = "")
    
    # initialize file (to make sure it's empty)
    cat(file=filerates, "", sep="")
    # Header
    for(state1 in 1:nbState)
        for(state2 in 1:nbState)
            if(state1!=state2)
                cat(file=filerates, "lambda", state1, state2, " ", append=TRUE, sep = "")
    # First row
    cat(file=filerates, "\n", rates0, "\n", append = TRUE)
    
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
    
    # Controls
    if(is.null(controls))
        controls <- list()
    if(is.null(controls$kappa))
        controls$kappa <- 3
    if(is.null(controls$lenmin))
        controls$lenmin <- 3
    if(is.null(controls$lenmax))
        controls$lenmax <- 6
    if(is.null(controls$thin))
        controls$thin <- 100
    if(is.null(controls$prUpdateMove))
        controls$prUpdateMove <- 1
    if(is.null(controls$SDP))
        controls$SDP <- 0.15
    
    return(list(obs=obs,
                map=map,
                nbState=nbState,
                nbIter=nbIter,
                adapt=adapt,
                par0=par,
                rates0=rates0,
                priorMean=priorMean,
                priorSD=priorSD,
                priorShape=priorShape,
                proposalSD=proposalSD,
                controls=controls,
                homog=homog,
                fileparams=fileparams,
                filerates=filerates,
                aSwitches=aSwitches,
                nbActual=nbActual,
                whichActual=whichActual))
}
