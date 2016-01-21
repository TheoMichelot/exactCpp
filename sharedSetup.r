# Create initial set of actual switches

ANswitch=0
ATswitch<-c()
ASswitch<-c()
AXswitch<-c()
AYswitch<-c()
AHswitch<-c()

dt=0.1*min(diff(fixes[,q.Time]))
for(kkk in 2:ndata)
{ 
   if( (fixes[kkk,4]!=fixes[kkk-1,4]))
   {
      ATswitch<-c(ATswitch,fixes[kkk,3]-dt)
      ASswitch<-c(ASswitch,fixes[kkk,4])
      AXswitch<-c(AXswitch,fixes[kkk,1])
      AYswitch<-c(AYswitch,fixes[kkk,2])
      AHswitch<-c(AHswitch,fixes[kkk,5])
      ANswitch=ANswitch+1
   }
}

Ajump<-matrix(c(1),1,ANswitch)
Ajump=c(Ajump)
Abehav=rep(0,length(Ajump))
F.Aswitch=cbind(X=AXswitch,
                Y=AYswitch,
                Time=ATswitch,
                State=ASswitch,
                Habitat=AHswitch,
                Jump=Ajump,
                Behav=Abehav)


# Define functions

successor=function(jumpLambda)
   sample(1:nstate,size=1,prob=jumpLambda)

source("minitabGeneric.r") 


# Initialise

oldLikeX=-Inf
oldLikeY=-Inf
bk=0
acc.traj=0
acc.move=0
acc.loc=0

# Specific output files

filedata <- paste("datastate", d, ".txt", sep = "")
cat(file=filedata,paste("data",1:ndata,sep=""), "\n")
cat(file=filedata, fixes[,4], "\n", append = TRUE)

filediag <- paste("diagnostic", d, ".txt", sep = "")
cat(file=filediag, "acc.traj", "acc.move", "acc.loc","\n")
cat(file=filediag, acc.traj, acc.move, acc.loc,"\n", append = TRUE)

# Controls

pp.min=3
pp.max=6
niter=1000
thin=100
