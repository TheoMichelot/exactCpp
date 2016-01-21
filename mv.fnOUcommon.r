# To avoid having to name arguments and use "...", assume for now exactly 3 parameter objects!
# state, deltat, xx, yy may be scalars or vectors

mv.fn.raw=function(PAR1,PAR2,PAR3,state,deltat,xx,yy)
{
  #
  # All states OU, with common centre
  #
  mux=PAR1[1]
  muy=PAR1[2]
  b=-PAR2[state]
  V=PAR3[state]
  phi=exp(b*deltat);c1=1-phi;SD=sqrt(V*(1-phi^2))
  return(list(emx=c1*(mux-xx),emy=c1*(muy-yy),sdx=SD,sdy=SD))
}