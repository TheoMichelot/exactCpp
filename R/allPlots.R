
allPlots <- function(nbState, fileparams, filerates, states, mty, truePar=NULL, trueState=NULL)
{
    library(scales)
    pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")
    par(mfrow=c(1,1))
    
    estim <- read.table(fileparams,header=TRUE)
    start <- floor(nrow(estim)/2)
    end <- nrow(estim)
    
    # which states are BM and which are OU?
    whichBM <- which(mty==1)
    whichOU <- which(mty==2)
    
    mux <- estim[(start:end),2*(1:nbState)-1]
    muy <- estim[(start:end),2*(1:nbState)]
    b <- estim[(start:end),(2*nbState+1):(3*nbState)]
    v <- estim[(start:end),(3*nbState+1):(4*nbState)]
    
    # plot bounds for mu
    muxmax <- max(mux[,whichOU])
    muxmin <- min(mux[,whichOU])
    muymax <- max(muy[,whichOU])
    muymin <- min(muy[,whichOU])
    
    # plot bounds for b vs v
    bmax <- max(log(b[,whichOU]))
    bmin <- min(log(b[,whichOU]))
    vmax <- max(log(v[,whichOU]))
    vmin <- min(log(v[,whichOU]))
    
    # make sure true values are within plot bounds, if provided
    if(!is.null(truePar)) {
        truemux <- truePar[2*(1:nbState)-1]
        truemuy <- truePar[2*(1:nbState)]
        trueb <- truePar[(2*nbState+1):(3*nbState)]
        truev <- truePar[(3*nbState+1):(4*nbState)]
        
        muxmax <- max(muxmax, truemux[whichOU])
        muxmin <- min(muxmin, truemux[whichOU])
        muymax <- max(muymax, truemuy[whichOU])
        muymin <- min(muymin, truemuy[whichOU])
        bmax <- max(bmax, log(-trueb[whichOU]))
        bmin <- min(bmin, log(-trueb[whichOU]))
        vmax <- max(vmax, log(truev[whichOU]))
        vmin <- min(vmin, log(truev[whichOU]))
    }

    # to ensure x and y scales are identical
    muxmid <- (muxmin+muxmax)/2
    muymid <- (muymin+muymax)/2
    lxy <- max(muxmax-muxmin,muymax-muymin)/2
    bmid <- (bmin+bmax)/2
    vmid <- (vmin+vmax)/2
    lbv <- max(bmax-bmin,vmax-vmin)/2
    
    # plot mu
    plot(mux[,whichOU[1]],muy[,whichOU[1]],pch=19,cex=0.2,col=alpha(pal[1],0.5),
         xlab="mu_x",ylab="mu_y",xlim=c(muxmid-lxy,muxmid+lxy),ylim=c(muymid-lxy,muymid+lxy))
    
    if(!is.null(truePar))
        points(truemux[whichOU[1]],truemuy[whichOU[1]],pch=19)
    
    if(length(whichOU)>1) {
        for(i in whichOU[-1]) {
            points(mux[,i],muy[,i],pch=19,cex=0.2,col=alpha(pal[i],0.5))
            
            if(!is.null(truePar))
                points(truemux[i],truemuy[i],pch=19)
        }
    }
    
    # plot log(b) vs log(v)
    plot(log(b[,whichOU[1]]),log(v[,whichOU[1]]),pch=19,cex=0.2,
         col=alpha(pal[1],0.5),xlim=c(bmid-lbv,bmid+lbv),ylim=c(vmid-lbv,vmid+lbv),
         xlab="log(b)",ylab="log(v)")
    
    if(!is.null(truePar))
        points(log(-trueb[whichOU[1]]),log(truev[whichOU[1]]),pch=19)
    
    if(length(whichOU)>1) {
        for(i in whichOU[-1]) {
            points(log(b[,i]),log(v[,i]),pch=19,cex=0.2,
                   col=alpha(pal[i],0.5))
            
            if(!is.null(truePar))
                points(log(-trueb[i]),log(truev[i]),pch=19)
        }
    }
    
    # trace plots of BM variance
    if(length(whichBM)>0) {
        vmax <- max(log(v[,whichBM]))
        vmin <- min(log(v[,whichBM]))
        if(!is.null(truePar)) {
            vmax <- max(vmax, log(truev[whichBM]))
            vmin <- min(vmin, log(truev[whichBM]))
        }
        
        plot(log(v[,whichBM[1]]),type="l",ylab="v",main=paste("State",whichBM[1]),
             ylim=c(vmin,vmax))
        
        if(!is.null(truePar))
            abline(h=log(truev[whichBM[1]]),col=2)
        
        for(i in whichBM[-1]) {
            plot(log(v[,i]),type="l",ylab="v",main=paste("State",i))
            
            if(!is.null(truePar))
                abline(h=log(truev[i]),col=2)
        }
    }
    
    # trace plots of rates
    if(nbState>1) {
        rates <- read.table(filerates,header=TRUE)
        par(mfrow=c(nbState,nbState-1))
        for(i in 1:(nbState*(nbState-1)))
            plot(rates[start:end,i],type="l")
        par(mfrow=c(1,1))
    }
    
    # plot most probable states
    maxState <- apply(states,1,which.max)
    if(is.null(trueState))
        plot(maxState,type="o",pch=19,cex=0.5,ylab="state",main="Decoded states")
    else {
        par(mfrow=c(2,1))
        plot(maxState,type="o",pch=19,cex=0.5,ylab="state",main="Decoded states")
        plot(trueState,type="o",pch=19,cex=0.5,ylab="state",main="True states")
    }
}