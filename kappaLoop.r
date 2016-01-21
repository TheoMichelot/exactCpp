# Version 1.5 draft for MEE
# functionally identical to version 1.4.5

source("basicSetup.r")
source("sharedSetup.r")

for (iii in 1:niter)
{
    #cat("Start iteration",iii,"\n")
    pp.length=sample(pp.min:pp.max,size=1)
    ppoint1=sample(1:(ndata-pp.length+1),size=1)
    ppoint2=ppoint1+pp.length-1
    
    fData <- fixes[ppoint1:ppoint2,]
    
    interv <- dim(fData)[1]
    Tbeg <- fData[1,q.Time]
    Tend <- fData[interv,q.Time]
    Sbeg <- fData[1,q.State]
    Send <- fData[interv,q.State]
    
    # Simulate
    #
    # Potential switches
    
    Nswitch <- rpois(1,(Tend-Tbeg)*kappa)
    
    NswitchNA <- rep(NA,Nswitch)
    F.switch <- cbind(X=NswitchNA,
                      Y=NswitchNA,
                      Time=sort(runif(Nswitch,Tbeg,Tend)),
                      State=NswitchNA,
                      Habitat=NswitchNA,
                      Jump=NswitchNA,
                      Behav=rep(0,Nswitch))
    
    #all simulation switch are in the first part, data are in second part
    
    F.all <- rbind(F.switch,fData)
    
    ranks <- rank(F.all[,q.Time])
    orders <- order(F.all[,q.Time])
    
    # First switch - switch from Tbeg
    if(bk==0)
        for (iij in 2:(Nswitch+interv))
        {
            j <- orders[iij]
            if( is.na(F.all[j,q.Habitat]) )    #data point or potential jump? potential jump has Habitat NA
            {
                orj1 <- orders[ranks[j]-1]
                S <- F.all[orj1,q.State]
                T <- F.switch[j,q.Time]-F.all[orj1,q.Time]
                
                Xfrom <- F.all[orj1,q.X]
                Yfrom <- F.all[orj1,q.Y]      
                
                move <- mv.fn.raw(m.paras,b.paras,v.paras,state=S,deltat=T,xx=Xfrom,yy=Yfrom)
                
                F.all[j,q.X] <- rnorm(1,mean=Xfrom+move$emx,sd=move$sdx)
                F.all[j,q.Y] <- rnorm(1,mean=Yfrom+move$emy,sd=move$sdy)
                
                # Mapping back on to known habitat
                
                F.all[j,q.Habitat] <- region.fn(F.all[j,q.X],F.all[j,q.Y],map)
                
                # Prob of actual switch depends on rates and on region
                probs <- lambda.fn(F.all[orders[ranks[j]-1],q.State],F.all[j,q.Habitat],F.switch[j,q.Time],para.la)/kappa
                jump.now <- runif(1)<sum(probs)
                F.all[j,q.State] <- ifelse(jump.now,successor(probs),F.all[orders[ ranks[j]-1],q.State])
                F.all[j,q.Jump] <- as.integer(jump.now)           
                
            } else   # we have a data point           
            {
                if((iij==(Nswitch+interv)|F.all[j,q.Behav]==1)&(F.all[orders[ranks[j]-1],q.State]!=F.all[j,q.State]))  
                    # state at Tend cannot be changed
                {
                    #cat("Wrong final/known state (",iij==(Nswitch+interv),Ball[j],")\n")
                    bk <- 1
                    break            
                }                
                F.all[j,q.State] <- F.all[orders[ranks[j]-1],q.State]  
                #data point state is always same as previous  and no jump 
            }                 
        }  # End of loop over j>=2
    
    # likelihood
    if(bk==0)
    {       
        gdeltaT <- F.all[(Nswitch+2):(Nswitch+interv),q.Time]-F.all[orders[ranks[(Nswitch+2):(Nswitch+interv)]-1],q.Time]
        # ranks[...-1]  find what is the rank of the currnt point and what rank of the point before the current
        #orders [ ranks[...-1] ]  for the rank we looking for, what index of that point   order show you the index
        # of the points with rank 1 2 3 4 5...
        
        gStates <- F.all[orders[ranks[(Nswitch+2):(Nswitch+interv)]-1],q.State] 
        gwhich <- orders[ranks[(Nswitch+2):(Nswitch+interv)]-1]
        
        gmove <- mv.fn.raw(m.paras,b.paras,v.paras,state=gStates,deltat=gdeltaT,xx=F.all[gwhich,q.X],yy=F.all[gwhich,q.Y])
        RWX <- (dnorm(F.all[(Nswitch+2):(Nswitch+interv),q.X],mean=F.all[gwhich,q.X]+gmove$emx,sd=gmove$sdx,log=TRUE))
        RWY <- (dnorm(F.all[(Nswitch+2):(Nswitch+interv),q.Y],mean=F.all[gwhich,q.Y]+gmove$emy,sd=gmove$sdy,log=TRUE))
        
        #oldLike
        #only select the proposed trajectories fall in the interval    if nothing, use the data, for interval cover
        #jump from different habitat   we make some artificial jump point to compare with the proposed 
        
        Owhich <- (F.Aswitch[,q.Time]<Tend&F.Aswitch[,q.Time]>Tbeg)
        ANswitch <- sum(Owhich)
        F.Aall <- rbind(F.Aswitch[Owhich,],fData)     
        Aranks <- rank(F.Aall[,q.Time])
        Aorders <- order(F.Aall[,q.Time])
        
        A.precede <- Aorders[Aranks[(ANswitch+2):(ANswitch+interv)]-1]
        temp.AdeltaT <- (F.Aall[(ANswitch+2):(ANswitch+interv),q.Time]-F.Aall[A.precede,q.Time])
        temp.AStates <- F.Aall[A.precede,q.State]
        
        Amove <- mv.fn.raw(m.paras,b.paras,v.paras,state=temp.AStates,deltat=temp.AdeltaT,
                           xx=F.Aall[A.precede,q.X],yy=F.Aall[A.precede,q.Y])
        
        oldLikeX <- (dnorm(F.Aall[(ANswitch+2):(ANswitch+interv),q.X],
                           mean=F.Aall[A.precede,q.X]+Amove$emx,sd=Amove$sdx,log=TRUE))
        oldLikeY <- (dnorm(F.Aall[(ANswitch+2):(ANswitch+interv),q.Y],
                           mean=F.Aall[A.precede,q.Y]+Amove$emy,sd=Amove$sdy,log=TRUE))   
        
        HR <- exp(sum(RWX+RWY)-sum(oldLikeX+oldLikeY))
        
        if(runif(1)<HR) 
        {
            #cat("Accept trajectory:",Tbeg,"to",Tend,"\n")
            
            # Update proposals - states etc have changed
            F.switch <- F.all[0:Nswitch,,drop=FALSE]
            
            # Update Data states
            fixes[fixes[,q.Time]>(Tbeg-0.1)&fixes[,q.Time]<(Tend+0.1),q.State]=F.all[(Nswitch+1):(Nswitch+interv),q.State]
            
            #data's jump do not need updating it is always 0
            
            #there is nothing outside the interval of Tbeg:Tend so that we can use the new Tswitch over old ATswitch
            if(!any(F.Aswitch[,q.Time]>Tend|F.Aswitch[,q.Time]<Tbeg))# Rare!      
            {  
                #cat("Rare!\n")
                F.Aswitch <- F.switch
            } else  # there IS something outside the interval we have to keep those 
                # and replace the stuff in the interval and then resort
            {                   
                out.which <- F.Aswitch[,q.Time]>Tend|F.Aswitch[,q.Time]<Tbeg
                
                F.Aswitch <- rbind(F.Aswitch[out.which,],F.switch)     
                
                orderA <- order(F.Aswitch[,q.Time])
                
                F.Aswitch <- F.Aswitch[orderA,,drop=FALSE]
                ANswitch <- nrow(F.Aswitch)            
            }
            
            acc.traj <- acc.traj+1         
        }#end of  "if accept trajectory"      
    } 
    
    #cat("F.Aswitch\n")
    #print(F.Aswitch)
    ANswitch <- nrow(F.Aswitch) # this ANswitch has been updated to match new F.Aswitch
    #cat("a) ANswitch",ANswitch,"\n")
    
    # All data need to be added to the saved trajectories to update lambda and kappa. because Ajump is not only for
    # proposed location, animal could have jump from data to the proposed location.
    
    F.Kall <- rbind(F.Aswitch,fixes)
    Kranks <- rank(F.Kall[,q.Time])
    Korders <- order(F.Kall[,q.Time])
    
    if(runif(1)<p.update.mv) # start of "if" on updating movement
    {   
        # Update movement parameters
        # Pick out points that are informative about movement,
        # their predecesors, and the states during movement
        
        whichi <- F.Kall[,q.Time]>min(fixes[,q.Time])
        predec <- Korders[Kranks[whichi]-1]
        states <- F.Kall[predec,q.State]
        
        # Calculate changes in position and time
        
        x.pre <- F.Kall[predec,q.X]
        dx.any <- F.Kall[whichi,q.X]-x.pre
        y.pre <- F.Kall[predec,q.Y]
        dy.any <- F.Kall[whichi,q.Y]-y.pre
        deltaT.any <-  F.Kall[whichi,q.Time]-F.Kall[predec,q.Time]
        
        move.step <- update.move.fn(b.paras, v.paras, homog.b, b.proposal.sd, nstate, homog.v, v.proposal.sd, m.paras, m.proposal.sd, m.prior.mean, m.prior.sd, b.prior.mean, b.prior.sd, v.prior.mean, v.prior.sd)

        # Old & new likelihoods
        
        old.move <- mv.fn.raw(m.paras,b.paras,v.paras,state=states,deltat=deltaT.any,xx=x.pre,yy=y.pre)
        new.move <- mv.fn.raw(move.step$m.prime,move.step$b.prime,move.step$v.prime,state=states,deltat=deltaT.any,xx=x.pre,yy=y.pre)
        
        oldLogLX <- sum(dnorm(dx.any,mean=old.move$emx,sd=old.move$sdx,log=TRUE))
        oldLogLY <- sum(dnorm(dy.any,mean=old.move$emy,sd=old.move$sdy,log=TRUE))
        
        newLogLX <- sum(dnorm(dx.any,mean=new.move$emx,sd=new.move$sdx,log=TRUE))
        newLogLY <- sum(dnorm(dy.any,mean=new.move$emy,sd=new.move$sdy,log=TRUE))
        
        logHR <- move.step$new.logprior-move.step$old.logprior+newLogLX+newLogLY-oldLogLX-oldLogLY
        
        if(runif(1)<exp(logHR))
        {
            # Accept movement parameters
            # cat("Accept movement parameters\n")
            acc.move <- acc.move+1
            m.paras <- move.step$m.prime
            b.paras <- move.step$b.prime
            v.paras <- move.step$v.prime
        }
    } # end of "if" on updating movement
    
    if(iii%%thin==0)
    {
        cat(file=fileparams, round(m.paras,6), round(b.paras,6), round(v.paras,6), "\n", append = TRUE)
    }
    
    #update jump rate k
    
    #cat("b) ANswitch",ANswitch,"\n")
    if(bk==0&ANswitch>0)
    {
        para.la <- update.rate.fn(ANswitch, F.Kall[,q.State], Korders, Kranks, F.Aswitch[,q.Habitat], F.Aswitch[,q.Jump], F.Aswitch[,q.Time], kappa, shape1, shape2, para.la)     
    }
    bk <- 0                      
    
    ### Insert local update
    
    # Local update to predecessor to obs j
    
    j <- sample(2:(dim(fixes)[1]),size=1)
    
    ANswitch <- nrow(F.Aswitch)
    F.Kall <- rbind(F.Aswitch,fixes)
    Kranks <- rank(F.Kall[,q.Time])
    Korders <- order(F.Kall[,q.Time])
    
    jorder <- ANswitch+j
    jrank <- Kranks[ANswitch+j]
    j1order <- Korders[jrank-1]
    
    if(j1order<=ANswitch) # Should be moveable
    {
        j2order <- Korders[jrank-2]
        
        # Want to perturb (T,X,Y)[j1order]
        
        T2 <- F.Kall[j2order,q.Time]
        X2 <- F.Kall[j2order,q.X]
        Y2 <- F.Kall[j2order,q.Y]
        S2 <- F.Kall[j2order,q.State]
        
        T1 <- F.Kall[j1order,q.Time]
        X1 <- F.Kall[j1order,q.X]
        Y1 <- F.Kall[j1order,q.Y]
        S1 <- F.Kall[j1order,q.State]
        
        H1 <- region.fn(X1,Y1,map)
        
        Tj <- F.Kall[jorder,q.Time]
        Xj <- F.Kall[jorder,q.X]
        Yj <- F.Kall[jorder,q.Y]
        
        # Old log-likelihood
        
        oldMove2 <- mv.fn.raw(m.paras,b.paras,v.paras,state=S2,deltat=T1-T2,xx=X2,yy=Y2)
        
        oldLogLX2 <- dnorm(X1,mean=X2+oldMove2$emx,sd=oldMove2$sdx,log=TRUE)
        oldLogLY2 <- dnorm(Y1,mean=Y2+oldMove2$emy,sd=oldMove2$sdy,log=TRUE)
        
        oldMove1 <- mv.fn.raw(m.paras,b.paras,v.paras,state=S1,deltat=Tj-T1,xx=X1,yy=Y1) 
        
        oldLogLX1 <- dnorm(Xj,mean=X1+oldMove1$emx,sd=oldMove1$sdx,log=TRUE)
        oldLogLY1 <- dnorm(Yj,mean=Y1+oldMove1$emy,sd=oldMove1$sdy,log=TRUE)
        
        oldLL <- oldLogLX1+oldLogLY1+oldLogLX2+oldLogLY2  
        
        # Propose new location
        
        Xnew <- rnorm(1,mean=X1,sd=SDP)
        Ynew <- rnorm(1,mean=Y1,sd=SDP)
        Tnew <- runif(1,min=T2,max=Tj)
        Hnew <- region.fn(Xnew,Ynew,map)
        
        # Effect of habitat
        
        rate1 <- lambda.fn(S2,H1,T1,para.la);rate1[S2]=kappa-sum(rate1)
        rate.new <- lambda.fn(S2,Hnew,Tnew,para.la);rate.new[S2]=kappa-sum(rate.new)
        newHfactor <- rate.new[S1]/rate1[S1]
        
        # New log-likelihood
        
        newMove2 <- mv.fn.raw(m.paras,b.paras,v.paras,state=S2,deltat=Tnew-T2,xx=X2,yy=Y2)
        
        newLogLX2 <- dnorm(Xnew,mean=X2+newMove2$emx,sd=newMove2$sdx,log=TRUE)
        newLogLY2 <- dnorm(Ynew,mean=Y2+newMove2$emy,sd=newMove2$sdy,log=TRUE)
        
        newMove1 <- mv.fn.raw(m.paras,b.paras,v.paras,state=S1,deltat=Tj-Tnew,xx=Xnew,yy=Ynew) 
        
        newLogLX1 <- dnorm(Xj,mean=Xnew+newMove1$emx,sd=newMove1$sdx,log=TRUE)
        newLogLY1 <- dnorm(Yj,mean=Ynew+newMove1$emy,sd=newMove1$sdy,log=TRUE)
        
        newLL <- newLogLX1+newLogLY1+newLogLX2+newLogLY2  
        
        logHR <- newLL-oldLL+log(newHfactor)
        # Calc needs to be done in this order, to avoid special case: exp(newLL-oldLL)=Inf, Hfactor=0
        
        if(runif(1)<exp(logHR))
        {
            #cat("Accept local update",T2,"to",Tj,"\n")
            acc.loc <- acc.loc+1
            F.Aswitch[j1order,q.Time] <- Tnew
            F.Aswitch[j1order,q.X] <- Xnew
            F.Aswitch[j1order,q.Y] <- Ynew
        } # end accept
    } # end previous point imputed
    # else can't do anything locally
    
    ### End of local update
    
    if(iii%%thin==0)
    {
        cat(file=filekappa,para.la,"\n", append = TRUE)
        cat(file=filedata, fixes[1:ndata,4], "\n", append = TRUE)
    }
    
    ### Diagnostics
    
    if(iii%%thin==0)
    {
        cat("\nEnd iteration",iii,"\n")
        cat(file=filediag, acc.traj, acc.move, acc.loc, "\n", append = TRUE)
    }
    
    rm(F.Kall)

} #end of for loop niter