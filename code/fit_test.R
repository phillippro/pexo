ParAll <- Par
  Data0 <- Data
    Ncore <- ParAll$Ncore
    Niter0 <- ParAll$Niter
    ParOpt <- ParIni0 <- ParIni
    ParMin0 <- ParMin
    ParMax0 <- ParMax
    Res <- ParML <- ParMP <- McOpt <- list()
    Nsig <- 0
    for(np in 0:ParAll$Nmax){
#    for(np in 1:ParAll$Nmax){
        ParAll$Np <- np
###Parallel signal identification and constraints
        if(np>0){
            ParAll$Np <- 1
            ParIni <- fit_Add1Kep(ParIni0,KepIni,0,ParAll$KepName)
            ParMin <- fit_Add1Kep(ParMin0,KepMin,0,ParAll$KepName)
            ParMax <- fit_Add1Kep(ParMax0,KepMax,0,ParAll$KepName)
        }
        ParAll$Npar <- length(ParIni)

####hot chain
        cat('\n Find signal',np,'using hot chain!\n')
###               dlogP <- (KepMax['logP']-KepMin['logP'])/Ncore
###                ParIni['logP1'] <- KepMin['logP']+dlogP*runif(1,ncore-1,ncore)
####initial MCMC constraints of a signal; output posterior sample and par.opt
        tmp <- fit_HotChain(Data,ParIni,ParMin,ParMax,ParAll,verbose=verbose)
        ParIni <- tmp$ParOpt
#    }

####cold chain combined constraint
        cat('\n Constrian',np,'signals using cold chain!\n')
        mcmc <- foreach(ncore=1:Ncore,.errorhandling='pass') %dopar% {
            ParAll$Np <- np
            Data <- Data0
            ParAll$Niter <- Niter0
            Data <- Data0
            if(np>1){
                KepOpt <- tmp$ParOpt[paste0(ParAll$KepName,1)]
                ParIni <- fit_Add1Kep(ParOpt,KepOpt,np-1,ParAll$KepName)
                ParMin <- fit_Add1Kep(ParOpt,KepMin,np-1,ParAll$KepName)
                ParMax <- fit_Add1Kep(ParOpt,KepMax,np-1,ParAll$KepName)
            }

            tmp <- try(fit_AM(Data,RateBary,ParIni,ParMin,ParMax,ParAll,verbose=verbose),TRUE)
###burnin
            if(class(tmp)!='try-error'){
                out <- tmp[-(1:floor(nrow(tmp)/2)),]
                out
            }
        }

###remove bad parallel chains
        ind <- which(sapply(1:length(mcmc),function(k) is.null(dim(mcmc[[k]]))))
        if(length(ind)>0) mcmc <- mcmc[-ind]

#        ind.rm <- which(is.na(mcmc[,ParAll$Npar+2]) | (abs(mcmc[,ParAll$Npar+2])==Inf))
#        if(length(ind.rm)>0) mcmc <- mcmc[-ind.rm,]

###Save mcmc results
        llmax <- sapply(1:length(mcmc), function(i) max(mcmc[[i]][,ncol(mcmc[[i]])-1]))
        cat('Maximum likelihood for all chains:',llmax,'\n')
        mcopt <- mcmc[[which.max(llmax)]]
        McOpt[[paste0('sig',np)]] <- mcopt

###subtract the best-fit in order to find the next signal
        ParAll$Npar <- ncol(mcopt)-3
        ParOpt <- mcopt[which.max(mcopt[,'loglike']),1:ParAll$Npar]
        ParML[[paste0('sig',np)]] <- ParOpt
        ParMP[[paste0('sig',np)]] <- mcopt[which.max(mcopt[,'logpost']),1:ParAll$Npar]
        tmp <- fit_CalRes(Data,RateBary,ParOpt,ParAll)
        Res[[paste0('sig',np)]] <- tmp
        if(ParAll$Np>0){
            Data <- tmp$DataRes
        }

###whether to stop DRAM according to BF criterion
        ll <-  max(mcopt[,'loglike'])
        if(ParAll$Np>0){
            lnbf3 <- ll-ll0-1.5*log(ParAll$Nepoch)
            cat('lnbf3=',lnbf3,'\n')
            if(lnbf3<5){
                break()
            }else{
                Nsig <- Nsig+1
            }
        }
        ll0 <- ll
    }
    ParAll$Niter <- Niter0
