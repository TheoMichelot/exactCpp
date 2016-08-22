
library(scales) # for nice colors
par(mfrow=c(1,1))

# movement parameters
par <- c(3,6,9,6, # mu
         -2,-2, # b
         2,2) # v

lambdapar <- rep(0.5,2)

estim <- read.table(mod$fileparams,header=TRUE)
start <- nrow(estim)/2
end <- nrow(estim)

# plot bounds for mu
muxmax <- max(estim[(start:end),(2*(1:nbState)-1)])
muxmin <- min(estim[(start:end),(2*(1:nbState)-1)])
muymax <- max(estim[(start:end),(2*(1:nbState))])
muymin <- min(estim[(start:end),(2*(1:nbState))])

# plot bounds for b vs v
bmax <- max(log(estim[(start:end),(2*nbState+1):(3*nbState)]))
bmin <- min(log(estim[(start:end),(2*nbState+1):(3*nbState)]))
vmax <- max(log(estim[(start:end),(3*nbState+1):(4*nbState)]))
vmin <- min(log(estim[(start:end),(3*nbState+1):(4*nbState)]))

# plot mu
plot(estim[start:end,1],estim[start:end,2],pch=19,cex=0.2,col=alpha(2,0.3),
     xlab="mu_x",ylab="mu_y",xlim=c(muxmin,muxmax),ylim=c(muymin,muymax))
points(par[1],par[2],pch=19)

if(nbState>1) {
    for(i in 2:nbState) {
        points(estim[start:end,2*i-1],estim[start:end,2*i],pch=19,cex=0.2,col=alpha(i+1,0.3))
        points(par[2*i-1],par[2*i],pch=19)
    }
}

# plot b vs v
plot(log(estim[start:end,2*nbState+1]),log(estim[start:end,3*nbState+1]),pch=19,cex=0.2,
     col=alpha(2,0.3),xlim=c(bmin,bmax),ylim=c(vmin,vmax),xlab="log(b)",ylab="log(v)")
points(log(-par[2*nbState+1]),log(par[3*nbState+1]),pch=19)

if(nbState>1) {
    for(i in 2:nbState) {
        points(log(estim[start:end,2*nbState+i]),log(estim[start:end,3*nbState+i]),pch=19,cex=0.2,
               col=alpha(i+1,0.3))
        points(log(-par[2*nbState+i]),log(par[3*nbState+i]),pch=19)
    }
}

# trace plots of rates
if(nbState>1) {
    rates <- read.table(filekappa,header=TRUE)
    par(mfrow=c(nbState,nbState-1))
    for(i in 1:(nbState*(nbState-1)))
        plot(rates[start:end,i],type="l",main=paste("True value:",lambdapar[i]))
    par(mfrow=c(1,1))
}
