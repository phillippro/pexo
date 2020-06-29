t0 <- proc.time()
Ntry <- 1000
for(i in 1:Ntry){
  for(j in 1:Nset){
                instr <- sets[j]
                inds <- which(Data[index,'instrument']==instr)
#                rv <- Data[index[inds],2]
#                erv <- Data[index[inds],3]
#                jitter <- 0
#                n <- paste0('logjitterRv','.',star,'.',instr)
#                n1 <- paste0('jitterRv','.',star,'.',instr)
#                if(any(names(ParFit)==n)){
#                    jitter <- exp(ParFit[n])
#                }else if(any(names(ParFit)==n1)){
#                    jitter <- ParFit[n1]
#                }
#                ll <- sum(dnorm(rv,mean=RvHat[inds],sd=sqrt(erv^2+jitter^2),log=T))
#                llike <- llike+ll
            }
###calculate red noise component
}

t1 <- fit_TimeCount(t0,'dur1',ofac=Ntry)
