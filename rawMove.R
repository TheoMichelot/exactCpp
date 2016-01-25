
#' Raw movement
#' 
#' Returns mean and standard deviation of the animal's next locations, from the movement
#' parameters and the previous location.
#' 
#' @param PAR1 PAR2 PAR3 Movement parameters
#' @param state State of underlying process
#' @param deltaT Time interval between current and next location
#' @param xx X-coordinate of current location
#' @param yy Y-coordinate of current location
#' 
#' @details Assume for now exactly 3 parameter objects!
#' Argumetns 'state', 'deltaT', 'xx', 'yy' may be scalars or vectors.
rawMove <- function(PAR1,PAR2,PAR3,state,deltaT,xx,yy)
{
    #
    # All states OU, with common centre
    #
    mux <- PAR1[1]
    muy <- PAR1[2]
    b <- -PAR2[state]
    V <- PAR3[state]
    phi <- exp(b*deltaT)
    c1 <- 1-phi
    SD <- sqrt(V*(1-phi^2))
    return(list(emx=c1*(mux-xx),emy=c1*(muy-yy),sdx=SD,sdy=SD))
}