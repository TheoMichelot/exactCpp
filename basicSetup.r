# Read data

fixes=cbind(as.matrix(read.table("fisherXYT.csv",sep=",",header=FALSE,skip=3)),NA,NA,NA,NA)
dimnames(fixes) <- list(NULL, c("X","Y","Time","State","Habitat","Jump","Behav"))
ndata=nrow(fixes)

# Enable references by "name" 

q.X=1#"X"
q.Y=2#"Y"
q.Time=3#"Time"
q.State=4#"State"
q.Habitat=5#"Habitat"
q.Jump=6#"Jump"
q.Behav=7#"Behav"

# Read map of habitat

map=as.matrix(read.table("habitat.csv",sep=",",header=FALSE,skip=3))
nhabitat=3
nstate=nhabitat

# Define function to interrogate habitat

region.fn=function(xx,yy,map)
{
  nx=nrow(map)
  ny=ncol(map)
  xi=ceiling(xx)
  yi=ceiling(yy)
  xi[xi<=0]=1
  xi[xi>nx]=nx
  yi[yi<=0]=1
  yi[yi>ny]=ny
  map[cbind(xi,yi)]
}

# Set up rest of data matrix

fixes[,q.State]=fixes[,q.Habitat]=region.fn(fixes[,q.X],fixes[,q.Y],map) 
fixes[,q.Jump]=0 #jump for data point is always 0
fixes[,q.Behav]=0 #behavioural states not known

# Initial OU parameters

mu.x=9
mu.y=3.5
b=1.0
V=10
m.paras=c(mu.x,mu.y)
b.paras=rep(b,nstate) # all same to be uninformative AND for correct results with simpler models
v.paras=rep(V,nstate) #  " 

# Priors (on log scale for b, v)

m.prior.mean=m.paras
b.prior.mean=log(b.paras)
v.prior.mean=log(v.paras)

m.prior.sd=rep(0.7,2)
b.prior.sd=rep(2.0,nstate)
v.prior.sd=rep(2.0,nstate)

# MH proposals (on log scale for b, v)

m.proposal.sd=rep(0.01,2)
b.proposal.sd=rep(0.2,nstate)
v.proposal.sd=rep(0.1,nstate)

SDP=0.15 # for "local update"

# Other model-specific functions

source("adaptive.r") # to create struct.fn
source("setup.lambda.adaptive.r") 

source("update.moveFisher.r")
source("update.rateFisher.r")
source("mv.fnOUcommon.r")

lower.limit=0.5 # on lambda/kappa
rtruncbeta=function(n=1,a=1,b=1,low=lower.limit)
{
  p.min=pbeta(low,a,b)
  qbeta(p.min+runif(n)*(1-p.min),a,b)
}


# Initial lambda as a function of current behaviour and location

initL=3.0
lambda21=initL
lambda31=initL
lambda12=initL
lambda32=initL
lambda13=initL
lambda23=initL
para.la<-c(lambda12,lambda13,lambda21,lambda23,lambda31,lambda32)

# prior beta parameters for lambdas

shape1=12
shape2=4

kappa=4


# Output

d <- format(Sys.time(), "%Y-%m-%d-%H%M")
fileparams <- paste("params", d, ".txt", sep = "")

cat(file=fileparams, "mux","muy","b1","b2","b3","v1","v2","v3", "\n")
cat(file=fileparams, m.paras,b.paras,v.paras, "\n", append = TRUE)

filekappa <- paste("rates", d, ".txt", sep = "")
cat(file=filekappa, "L12","L13","L21","L23","L31","L32", "\n")
cat(file=filekappa, para.la, "\n", append = TRUE)

# Homogeneity between states - for model variation inc DIC calculation

homog.b=FALSE
homog.v=FALSE

p.update.mv=1

set.seed(1)