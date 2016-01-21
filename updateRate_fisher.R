
#' Update rate
#' 
#' @param allData Matrix of all observations and actual switches, sorted by time.
#' @param indSwitch Indices of actual switches among observations.
#' @param kappa Upper limit for switching rates
#' @param shape1 First shape parameter
#' @param shape2 Second shape parameter
updateRate <- function (allData, indSwitch, kappa, shape1, shape2)
{
    # counts[i,j,k] = number of switches from state i to state j while in habitat k
    counts <- minitabGeneric(allData[indSwitch-1,colState],allData[indSwitch,colState],
                             allData[indSwitch,colHabitat],ns=3,nr=3)
    
    lambda21 <- kappa*rtruncbeta(1,shape1+counts[2,1,1],shape2+counts[2,2,1])
    lambda31 <- kappa*rtruncbeta(1,shape1+counts[3,1,1],shape2+counts[3,3,1])
    lambda12 <- kappa*rtruncbeta(1,shape1+counts[1,2,2],shape2+counts[1,1,2])
    lambda32 <- kappa*rtruncbeta(1,shape1+counts[3,2,2],shape2+counts[3,3,2])
    lambda13 <- kappa*rtruncbeta(1,shape1+counts[1,3,3],shape2+counts[1,1,3])
    lambda23 <- kappa*rtruncbeta(1,shape1+counts[2,3,3],shape2+counts[2,2,3])
    
    return(c(lambda12,lambda13,lambda21,lambda23,lambda31,lambda32))
}
