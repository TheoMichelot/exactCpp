
#' Update rate (unconstrained case)
#' 
#' @param allData Matrix of all observations and switches, sorted by time.
#' @param indSwitch Indices of switches among observations.
#' @param kappa Upper limit for switching rates
#' @param shape1 First shape parameter
#' @param shape2 Second shape parameter
#' @param nbState Number of states 
updateRate_unconstr <- function (allData, indSwitch, kappa, shape1, shape2, nbState)
{
    # counts[i,j,k] = number of switches from state i to state j while in habitat k
    counts <- minitabGeneric(allData[indSwitch-1,colState],allData[indSwitch,colState],
                             allData[indSwitch,colHabitat],nbState=nbState,nbHabitat=1)
       
    alpha0 <- rep(shape1,nbState)
    LL <- matrix(NA,nbState,nbState-1)
    
    for (kk in 1:nbState) {
        alpha <- alpha0
        alpha[kk] <- shape2
        alpha <- alpha+counts[kk,,1]
        LL[kk,] <- kappa*rdirichlet(n=1,alpha)[-kk]
    }
    
    return(c(t(LL)))
}
