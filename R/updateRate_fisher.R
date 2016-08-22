
#' Update rate
#' 
#' @param allData Matrix of all observations and actual switches, sorted by time.
#' @param indSwitch Indices of actual switches among observations.
#' @param kappa Upper limit for switching rates
#' @param shape1 First shape parameter
#' @param shape2 Second shape parameter
#' @param nbState Number of states
#' 
#' @details This function requires nbState=nbHabitat
updateRate <- function (allData, indSwitch, kappa, shape1, shape2, nbState)
{
    # enable references by "name" 
    colX <- 1; colY <- 2; colTime <- 3; colState <- 4; colHabitat <- 5; colJump <- 6; colBehav <- 7
    
    # counts[i,j,k] = number of switches from state i to state j while in habitat k
    counts <- minitabGeneric(allData[indSwitch-1,colState],allData[indSwitch,colState],
                             allData[indSwitch,colHabitat],nbState=nbState,nbHabitat=nbState)

    # rates contains the non-diagonal elements of the switching rate matrix, ordered row-wise
    # (e.g. lambda12, lambda13, lambda21, lambda23, lambda31, lambda32 for 3-state model)
    rates <- rep(NA,nbState*(nbState-1))
    k <- 1
    for(i in 1:nbState) {
        for(j in 1:nbState) {
            if(i!=j) {
                rates[k] <- kappa*rtruncbeta(1,shape1+counts[i,j,j],shape2+counts[i,i,j])
                k <- k+1
            }    
        }
    }
    
    return(rates)
}
