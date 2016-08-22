
#' Returns an array of dimensions (nbState,nbState,nbHabitat) filled with number of occurrences of each element in
#' (a1,a2,a3)
#' 
#' @param a1 Vector of values for 1st dimension
#' @param a2 Vector of values for 2nd dimension
#' @param a3 Vector of values for 3rd dimension
#' @param nbState Number of states (number of rows and columns of the array)
#' @param nbHabitat Number of habitats (number of layers of the array)
minitabGeneric <- function(a1,a2,a3,nbState,nbHabitat)
{
    dims <- c(nbState,nbState,nbHabitat)
    pd <- prod(dims)
    bin <- a1+nbState*(a2-1)+nbState*nbState*(a3-1)
    return(array(tabulate(bin, pd), dims))
}
