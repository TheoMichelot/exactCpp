
#' Raw movement
#' 
#' Returns mean and standard deviation of the animal's next locations, from the movement
#' parameters and the previous location.
#' 
#' @param par Vector of movement parameters (m, b, v)
#' @param state State of underlying process
#' @param deltaT Time interval between current and next location
#' @param x X-coordinate of current location
#' @param y Y-coordinate of current location
#' @param nbState Number of states
#' 
#' @details Arguments 'state', 'deltaT', 'x', 'y' may be scalars or vectors.
rawMove <- function(par,state,deltaT,x,y,nbState)
{
    # unpack the vector of parameters
    m <- par[1:(2*nbState)]
    b <- par[(2*nbState+1):(3*nbState)]
    v <- par[(3*nbState+1):(4*nbState)]
    
    # select the right parameters
    mux <- m[2*state-1]
    muy <- m[2*state]
    b <- -b[state]
    v <- v[state]
    
    phi <- exp(b*deltaT)
    sd <- sqrt(v*(1-phi^2))
    
    return(list(emx=(1-phi)*(mux-x), emy=(1-phi)*(muy-y), sdx=sd, sdy=sd))
}