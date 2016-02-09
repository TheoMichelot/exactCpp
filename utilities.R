
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

#' Natural to working
#' 
#' Transforms the parameters from the natural scale to the working scale.
#' 
#' @param par Vector of parameters, in the order: m, b, v.
#' @param nbState Number of states.
n2w <- function(par,nbState)
{
    m <- par[1:(2*nbState)]
    b <- par[(2*nbState+1):(3*nbState)]
    v <- par[(3*nbState+1):(4*nbState)]
    
    wpar <- c(m,log(b),log(v))
    return(wpar)
}

#' Working to natural
#' 
#' Transforms the parameters from the working scale to the natural scale.
#' 
#' @param wpar Vector of parameters, in the order: m, b, v.
#' @param nbState Number of states.
n2w <- function(wpar,nbState)
{
    m <- par[1:(2*nbState)]
    b <- par[(2*nbState+1):(3*nbState)]
    v <- par[(3*nbState+1):(4*nbState)]
    
    par <- c(m,exp(b),exp(v))
    return(par)
}

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

# pick successor state
successor <- function(jumpLambda)
    sample(1:nbState,size=1,prob=jumpLambda)
