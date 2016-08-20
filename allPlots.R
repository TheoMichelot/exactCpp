
allPlots <- function(nbState, fileparams, filerates, truePar=NULL)
{
    library(scales)
    
    estim <- read.table(fileparams,header=TRUE)
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
    
    bmax <- log(0.5)
    vmin <- log(2)
    
    # plot mu
    plot(estim[start:end,1],estim[start:end,2],pch=19,cex=0.2,col=alpha(2,0.3),
         xlab="mu_x",ylab="mu_y",xlim=c(muxmin,muxmax),ylim=c(muymin,muymax))
    
    if(!is.null(truePar))
        points(truePar[1],truePar[2],pch=19)
    
    if(nbState>1) {
        for(i in 2:nbState) {
            points(estim[start:end,2*i-1],estim[start:end,2*i],pch=19,cex=0.2,col=alpha(i+1,0.3))
            
            if(!is.null(truePar))
                points(truePar[2*i-1],par[2*i],pch=19)
        }
    }
    
    # plot log(b) vs log(v)
    plot(log(estim[start:end,2*nbState+1]),log(estim[start:end,3*nbState+1]),pch=19,cex=0.2,
         col=alpha(2,0.3),xlim=c(bmin,bmax),ylim=c(vmin,vmax),xlab="log(b)",ylab="log(v)")
    
    if(!is.null(truePar))
        points(log(-truePar[2*nbState+1]),log(truePar[3*nbState+1]),pch=19)
    
    if(nbState>1) {
        for(i in 2:nbState) {
            points(log(estim[start:end,2*nbState+i]),log(estim[start:end,3*nbState+i]),pch=19,cex=0.2,
                   col=alpha(i+1,0.3))
            
            if(!is.null(truePar))
                points(log(-truePar[2*nbState+i]),log(truePar[3*nbState+i]),pch=19)
        }
    }
    
    # trace plots of rates
    if(nbState>1) {
        rates <- read.table(filekappa,header=TRUE)
        par(mfrow=c(nbState,nbState-1))
        for(i in 1:(nbState*(nbState-1)))
            plot(rates[start:end,i],type="l")
        par(mfrow=c(1,1))
    }
}