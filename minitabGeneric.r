minitabGeneric<-
function(a1,a2,a3,ns,nr)
{
dims=c(ns,ns,nr)
pd=prod(dims)
bin=a1+ns*(a2-1)+ns*ns*(a3-1)
array(tabulate(bin, pd), dims)
}
