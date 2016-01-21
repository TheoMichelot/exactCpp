
#' Returns an array of dimensions (ns,ns,nr) filled with number of occurrences of each element in
#' (a1,a2,a3)
#' 
#' @param a1 Vector of values for 1st dimension
#' @param a2 Vector of values for 2nd dimension
#' @param a3 Vector of values for 3rd dimension
#' @param ns Number of states (number of rows and columns of the array)
#' @param nr Number of observations (number of layers of the array)
minitabGeneric <- function(a1,a2,a3,ns,nr)
{
  dims <- c(ns,ns,nr)
  pd <- prod(dims)
  bin <- a1+ns*(a2-1)+ns*ns*(a3-1)
  return(array(tabulate(bin, pd), dims))
}
