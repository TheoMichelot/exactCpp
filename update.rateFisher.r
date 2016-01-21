update.rate.fn <- function (ANswitch, KSall, Korders, Kranks, AHswitch, Ajump, ATswitch, kappa, shape1, shape2, para.la)
{      
  jj=1:ANswitch
  
  ns=3;nr=3
  counts=minitabGeneric(KSall[Korders[Kranks[jj]-1]],KSall[Korders[Kranks[jj]]],AHswitch[jj],
                        ns,nr)
  lambda21=kappa*rtruncbeta(1,shape1+counts[2,1,1],shape2+counts[2,2,1])
  lambda31=kappa*rtruncbeta(1,shape1+counts[3,1,1],shape2+counts[3,3,1])
  lambda12=kappa*rtruncbeta(1,shape1+counts[1,2,2],shape2+counts[1,1,2])
  lambda32=kappa*rtruncbeta(1,shape1+counts[3,2,2],shape2+counts[3,3,2])
  lambda13=kappa*rtruncbeta(1,shape1+counts[1,3,3],shape2+counts[1,1,3])
  lambda23=kappa*rtruncbeta(1,shape1+counts[2,3,3],shape2+counts[2,2,3])
  c(lambda12,lambda13,lambda21,lambda23,lambda31,lambda32)
}
