
# estim <- read.table("params2016-04-07-1121.txt",header=TRUE)

logb1 <- log(estim[20000:50000,7])
logb2 <- log(estim[20000:50000,8])
logb3 <- log(estim[20000:50000,9])
logv1 <- log(estim[20000:50000,10])
logv2 <- log(estim[20000:50000,11])
logv3 <- log(estim[20000:50000,12])

bmax <- max(c(logb1,logb2,logb3))
bmin <- min(c(logb1,logb2,logb3))
vmax <- max(c(logv1,logv2,logv3))
vmin <- min(c(logv1,logv2,logv3))

plot(logb1,logv1,pch=19,cex=0.2,col=rgb(1,0,0,alpha=0.3),xlim=c(bmin,bmax),ylim=c(vmin,vmax),
     xlab="log(b)",ylab="log(v)")
points(logb2,logv2,pch=19,cex=0.2,col=rgb(0,1,0,alpha=0.3))
points(logb3,logv3,pch=19,cex=0.2,col=rgb(0,0,1,alpha=0.3))

# values of parameters used in simulation
points(log(1.5),log(5),pch=19)
points(log(0.8),log(3),pch=19)
points(log(0.2),log(1),pch=19)

plot(estim[20000:50000,1],estim[20000:50000,2],pch=19,cex=0.3,col=rgb(0,0,0,alpha=0.5),
     xlab="mu_x",ylab="mu_y")

points(6,6,pch=19,col="firebrick3")

# rates <- read.table("rates2016-04-07-1551.txt",header=TRUE)
par(mfrow=c(3,2))
for(i in 1:6)
    plot(rates[20000:50000,i],type="l")