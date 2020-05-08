###Dec. 15, 2019, fixed the bug in fit_ChangePar() by adding "ParNew$logmC1 <- log(Par$mT)"
###This file contains functions for fitting.
tol0 <- 1e-20
fit_PolyTrend <- function(t,b,coef){
####################################
## Polynomical trend
##
## Input:
##   t - input time
##   CoefPoly - polynomial coefficient
##
## Output:
##   trend - the values of the trend at the input times
####################################
    trend <- rep(b,length(t))
    for(j3 in 1:length(coef)){
        trend <- trend+coef[j3]*t^j3#t[,j3]
    }
    return(trend)
}

fit_OptPar <- function(mc,llp){
####################################
## find the optimal parameter at the maximum likelihood or maximum posterior
##
## Input:
##   mc - mcmc chain
##   llp - log likelihood/poseterior
##
## Output:
##   ParOpt - Optimal parameters
##   llpmax - maximum llp
####################################
    ind <- which.max(llp)
    return(list(ParOpt=mc[ind,],llpmax=llp[ind]))
}


fit_ARMA <- function(t,res,p,q,ParFit,instr,star,ARMAtype='abs',Dtype='rv'){
####################################
## ARMA noise model
##
## Input:
##   t - input time
##   res - data-model
##   ParFit - fitable parameters
##   p - order of the AR model
##   q - order of the MA model
##
## Output:
##   trend - the values of the trend at the input times
####################################
    dv.arma <- 0

    ##AR(p) model
    if(p>0){
        dVr.ar <- 0
        ampAR <- ParFit[grep(paste0('ampAR',Dtype,1:p,'.',star,'.',instr),names(ParFit))]
        logtauAR <- ParFit[paste0('logtauAR',Dtype,'.',star,'.',instr)]
        for(j in 1:p){
            dVr.ar <- dVr.ar + c(rep(0,j),ampAR[j]*exp((t[-(length(t)+1-(1:j))]-t[-(1:j)])/exp(logtauAR))*x[-(length(t)+1-(1:j))])
        }
        dv.arma <- dv.arma + dVr.ar
    }

    ##MA(q) model
    if(q>0){
        dVr.ma <- 0
        ampMA <- ParFit[paste0('ampMA',Dtype,1:q,'.',star,'.',instr)]
        logtauMA <- ParFit[[paste0('logtauMA',Dtype,'.',star,'.',instr)]]
        for(j in 1:q){
            dVr.ma <- dVr.ma + c(rep(0,j),ampMA[j]*exp(-abs(t[-(length(t)+1-(1:j))]-t[-(1:j)])/exp(logtauMA))*res[-(length(t)+1-(1:j))])
        }
        dv.arma <- dv.arma + dVr.ma
    }
    return(dv.arma)
}

fit_AstroTrend <- function(Data,ParFit,Par,star,type='abs'){
####################################
## Trend in astrometry
##
## Input:
##   Data - a list of data objects
##   ParFit - fitable parameters
##   Par - all input parameters
##   star - star name
##   type - astrometry type
##
## AstroTrend:
##   Output - trend in astromery
####################################
###global trend
    AstroTrend <- array(0,dim=c(nrow(Data),2))
    for(coord in c('ra','dec')){
        if(coord=='ra'){
            colum <- 1
        }else{
            colum <- 2
        }
        sets <- unique(Data[,'instrument'])
        Nset <- length(sets)
        for(j in 1:Nset){
            instr <- sets[j]
            index <- which(Data[,'instrument']==instr)
            b <- 0
            if(type=='abs'){
                n <- paste0('b',toupper(coord),'.',star,'.',instr)
                if(any(names(ParFit)==n)) b <- ParFit[n]*DMAS2R
            }else{
                n <- paste0('b',toupper(coord),'.rel.',instr)
                if(any(names(ParFit)==n)) b <- ParFit[n]
            }
            AstroTrend[index,] <- AstroTrend[index,]+b
        }
    }
    AstroTrend
}

fit_AstroKep <- function(OutTime,Data,OutObsT,OutObsC,ParFit,Par,star,geometry=TRUE){
####################################
## Calculate Keplerian RV
##
## Input:
##   OutTime - Output of time_Ta2te()
##   Data - a list of data objects
##   OutObs - Output of time_Utc2tb()
##   ParFit - fitable parameters
##   Par - all input parameters
##   star - star name
##
## Output:
##   AstroKep - astrometry in RA and DEC in units of rad
####################################
    ind <- which(Data$type=='abs')
    AstroKep <- array(0,dim=c(length(ind),2))
    if(length(ind)>0){
        if(star==Par$companion){
            if(Par$geometry){
                AstroKep[ind,] <- gen_Xyz2lb(OutTime$uOC[ind,])
            }else{
                AstroNew <- astro_FullModel(OutObs,OutTime,Par,Mlens=Par$mT,component='C')
            }
        }else{
            if(Par$geometry){
                AstroKep[ind,] <- gen_Xyz2lb(OutTime$uOT[ind,])
            }else{
                AstroNew <- astro_FullModel(OutObs,OutTime,Par,Mlens=Par$mC1,component='T')
            }
        }
        if(!Par$geometry)   AstroKep[ind,] <- AstroNew$DirObs[ind,]
    }
    AstroKep
}

fit_AstroRed <- function(OutTime,AstroHat,Data,ParFit,Par,star,type='abs'){
####################################
## red noise in astrometry
##
## Input:
##   OutTime - Output of time_Ta2te()
##   AstroHat - AstroTrend+AstroKep
##   Data - a list of data objects
##   ParFit - fitable parameters
##   Par - all input parameters
##   star - star name
##   type - astrometry type
##
## Output:
##   AstroRed - astrometry red noise
####################################
    AstroRed <- array(0,dim=c(nrow(Data),2))
    ind <- which(Data$type==type)
    if(length(ind)>0){
        sets <- unique(Data[ind,'instrument'])
        Nset <- length(sets)
        for(j in 1:Nset){
            instr <- sets[j]
            index <- ind[Data[ind,'instrument']==instr]
            x <- Data[index,2*k]
            q <- Par$ObsInfo[index[1],'q']
            p <- Par$ObsInfo[index[1],'p']
            if(p>0 | q>0){
                tauE <- rowSums(OutTime$tauE)[index]
                res <- (x-AstroHat[index,k])
                if(type=='abs'){
###from rad to mas
                    res <- res/DMAS2R
                }
                astro.arma <- fit_ARMA(t=tauE,res=res,p=p,q=q,ParFit=ParFit,instr=instr,star=star,Dtype='astro')
                AstroRed[index,k] <- astro.arma
            }
        }
    }
    if(type=='abs'){
###from mas to rad to be compared with data
        AstroRed*DMAS2R
    }else{
        AstroRed
    }
}

fit_AstroLike <- function(AstroHat,Data,ParFit,Par,star,type='abs'){
####################################
## Likelihood for absolulte astrometry model
##
## Input:
##   AstroHat - AstroTrend+AstroKep (unit: rad for abs; mas for rel)
##   Data - a list of data objects
##   ParFit - fitable parameters
##   Par - all input parameters
##   star - star name
##   type - astrometry type
##
## Output:
##   loglike - log likelihood
####################################
    loglike <- 0
    sets <- unique(Data[,'instrument'])
    Nset <- length(sets)
    for(j in 1:Nset){
        instr <- sets[j]
        index <- which(Data[,'instrument']==instr)
        astro <- Data[index,c(2,4)]
        eastro <- Data[index,c(3,5)]
        jitter <- 0
        n <- paste0('logjitterAstro','.',star,'.',instr)
        n1 <- paste0('jitterAstro','.',star,'.',instr)
        if(any(names(ParFit)==n)){
            jitter <- exp(ParFit[n])
        }else if(any(names(ParFit)==n1)){
            jitter <- ParFit[n1]
        }
        res <- astro-AstroHat[index,]
        if(type=='abs') res[,1] <- res[,1]*cos(AstroHat[index,2])
        dastro <- unlist(res)#2D residual to 1D residual
#        dastro <- res
        if(type=='abs') dastro <- dastro/DMAS2R
        ll <- sum(dnorm(dastro,mean=0,sd=unlist(sqrt(eastro^2+jitter^2)),log=T))
#        ll <- ll+sum(dnorm(dastro[,2],mean=0,sd=sqrt(eastro[ind,2]^2+jitter^2),log=T))
        loglike <- loglike+ll
    }
    loglike
}

fit_RvTrend <- function(OutTime,Data,ParFit,Par,star){
####################################
## Calculate reference RV or the RV of reference spectrum
##
## Input:
##   OutTime - Output of time_Ta2te()
##   Data - a list of data objects
##   ParFit - fitable parameters
##   Par - all input parameters
##   star - star name
##
## Output:
##   RvTrend - RV trend
####################################
    RvTrend <- rep(0,nrow(Data))
    n <- paste0('aRv','.',star)
    if(any(names(ParFit)==n)){
        tt <- (OutTime$tauE[,1]-Par$T0)+OutTime$tauE[,2]
        RvTrend <- fit_PolyTrend(tt/DJY,0,ParFit[n])
    }
    sets <- unique(Data[,'instrument'])
    Nset <- length(sets)
    for(j in 1:Nset){
        instr <- sets[j]
        index <- which(Data[,'instrument']==instr)
        b <- 0
        n <- paste0('bRv','.',star,'.',instr)
        if(any(names(ParFit)==n)) b <- ParFit[n]
        RvTrend[index] <- RvTrend[index]+b*1e3
    }
    RvTrend
}

fit_RvRef <- function(OutTime,Data,ParFit,Par,star){
####################################
## Calculate reference RV or the RV of reference spectrum
##
## Input:
##   OutTime - Output of time_Ta2te()
##   Data - a list of data objects
##   ParFit - fitable parameters
##   Par - all input parameters
##   star - star name
##
## Output:
##   RvRef - RV trend
####################################
    RvRef <- rep(0,nrow(Data))
    sets <- unique(Data[,'instrument'])
    Nset <- length(sets)
    for(j in 1:Nset){
        instr <- sets[j]
        index <- which(Data[,'instrument']==instr)
        b <- 0
        n <- paste0('bRv','.',star,'.',instr)
        if(any(names(ParFit)==n)) b <- ParFit[n]
        RvRef[index] <- RvRef[index]+b*1e3
    }
###add minor signal because 1+Zref=(1+Zabs)/(1+Zmeas) and Zref*c=RvRef, so RvRef=Zabs*c-Zmeas*c-Zabs*Zmeas*c=RvKep-RvData-RvKep*RvData/c; where Zabs is the absolute redshift.
    RvRef
}

fit_RvKep <- function(OutTime,Data,OutObs,Par){
####################################
## Calculate Keplerian RV
##
## Input:
##   OutTime - Output of time_Ta2te()
##   Data - a list of data objects
##   OutObs - Output of time_Utc2tb()
##   Par - all input parameters
##   star - star name
##
## Output:
##   RvKep - Keplerian RV
####################################
    ind <- which(Data$type=='rv')
    RvKep <- array(NA,dim=c(nrow(Data),1))
    if(length(ind)>0){
#            OutRvNew <- rv_FullModel(OutObs,OutTime,Par,component='C')
        OutRvNew <- rv_FullModel(OutObs,OutTime,Par,component='T')
        if(!Par$binary){
            RvKep[ind] <- -OutRvNew$Zcomb$ZB[ind]*CMPS
###                RvHat <- -OutRvNew$Zcomb$ZBwe[ind]*CMPS
###            RvHat <- OutRvNew$RvTot[Par$IndRv]-Par$rv*1e3#remove systematic velocity
        }else{
#            RvKep[ind] <- OutRvNew$RvTot[ind]-Par$rv*1e3#remove systematic velocity
            RvKep[ind] <- OutRvNew$RvTot[ind]
#	    cat('RvKep[ind]=',head(RvKep[ind],'\n'))
        }
    }
    RvKep
}

fit_RvRed <- function(OutTime,RvHat,Data,ParFit,Par,star){
####################################
## RV red noise
##
## Input:
##   OutTime - Output of time_Ta2te()
##   RvHat - RvTrend+RvKep
##   Data - a list of data objects
##   ParFit - fitable parameters
##   Par - all input parameters
##   star - star name
##
## Output:
##   RvRed - RV red noise
####################################
    sets <- unique(Data[,'instrument'])
    Nset <- length(sets)
    RvRed <- rep(0,nrow(Data))
    for(j in 1:Nset){
        instr <- sets[j]
        index <- which(Data[,'instrument']==instr)
        rv <- Data[,2]
        q <- Par$ObsInfo[index[1],'q']
        p <- Par$ObsInfo[index[1],'p']
        if(p>0 | q>0){
            tauE <- rowSums(OutTime$tauE)[index]
            rv.arma <- fit_ARMA(t=tauE,x=rv,xhat=RvHat[index],p=p,q=q,ParFit=ParFit,instr=instr,star=star,Dtype='rv')
            RvRed[index] <- rv.arma
        }
    }
    RvRed
}

fit_RvLike <- function(RvHat,Data,ParFit,star){
####################################
## Likelihood for RV model
##
## Input:
##   RvHat - a function of RvRef and RvKep; RvHat=RvHat-RvRef-RvHat*RvHat/c
##   Data - a list of data objects
##   ParFit - fitable parameters
##   star - star name
##
## Output:
##   loglike - log likelihood
####################################
    loglike <- 0
    sets <- unique(Data[,'instrument'])
    Nset <- length(sets)
    for(j in 1:Nset){
        instr <- sets[j]
        index <- which(Data[,'instrument']==instr)
        rv <- Data[index,2]
        erv <- Data[index,3]
        jitter <- 0
        n <- paste0('logjitterRv','.',star,'.',instr)
        n1 <- paste0('jitterRv','.',star,'.',instr)
        if(any(names(ParFit)==n)){
            jitter <- exp(ParFit[n])
        }else if(any(names(ParFit)==n1)){
            jitter <- ParFit[n1]
        }
        ll <- sum(dnorm(rv,mean=RvHat[index],sd=sqrt(erv^2+jitter^2),log=T))
        loglike <- loglike+ll
    }
    loglike
}

fit_ChangePar <- function(Par,component='T'){
####################################
## change target parameters into companion parameters
##
## Input:
##   Par - input parameters
##
## Output:
##   ParNew - New parameters
####################################
    ParNew <- Par
    ParNew$mT <- exp(Par$logmC1)
    ParNew$star <- Par$stars[2]
    ParNew$mC1 <- Par$mT
    ParNew$logmC1 <- log(Par$mT)
    ParNew$omegaT1 <- (Par$omegaT1+pi)%%(2*pi)
    ParNew
}

fit_ModelPredict <- function(Data,Par,ParOpt,Nsim=1e3){
####################################
## change target parameters into companion parameters
##
## Input:
##   Data - input data
##   Par - input parameters
##   ParOpt - optimal parameters
##
## Output:
##   pred - Model prediction
##   sim - Simulated data
####################################
    stars <- unique(Data$star)
    Sim <- c()
    Obs <- c()
    for(star in stars){
        ind0 <- which(Data$star==star)
        types <- unique(Data$type[ind0])
        for(type in types){
            ind1 <- which(Data$star==star & Data$type==type)
            inss <- unique(Data[ind1,'instrument'])
            for(ins in inss){
                ind <- which(Data$star==star & Data$type==type & Data$instrument==ins)
                tmin <- min(Data[ind,1])
                tmax <- max(Data[ind,1])
#                tsim <- seq(tmin,tmax,length.out=Nsim)
                tsim <- sort(runif(Nsim,tmin,tmax))
                Sim <- rbind(Sim,cbind(tsim,0,0,0,0,Data$star[ind[1]],Data$type[ind[1]],Data$instrument[ind[1]],Data$wavelength[ind[1]]))
                Obs <- rbind(Obs,t(replicate(Nsim,unlist(Par$ObsInfo[ind[1],]))))
            }
        }
    }
    Par$ObsInfo <- fit_changeType(Obs,coltype=c('n','n','n','n','n','n','c','c','n','n','n','n','n','n','n'))
    Sim <- fit_changeType(Sim,coltype=c(rep('n',5),rep('c',3),'n'))
    colnames(Sim) <- colnames(Data)
    utc <- time_ChangeBase(cbind(Sim[,1],0))#simulated utc
    Par$Nepoch <- nrow(utc)
    ParNew <- update_par(Par,ParOpt)
    tmp <- gen_CombineModel(utc,Sim,ParNew,component=Par$component)
    OutObs <- tmp$OutObs
    OutTime0 <- OutTime <- tmp$OutTime
    OutAstro <- tmp$OutAstro
    OutRv <- tmp$OutRv
    fit <- fit_LogLike(Sim,OutObs,RateObs,ParOpt,ParNew,OutTime0=OutTime0)
    list(sim=Sim,pred=fit$model)
}

#fit_LogLike <- function(Data,OutObs,RateObs,ParFit,Par,OutTime0=NULL,verbose=TRUE){
fit_LogLike <- function(Data,OutObs,RateObs,ParFit,Par,OutTime0=NULL,verbose=FALSE){
####################################
## Calculate log Likelihood function
##
## Input:
##   Data - a list of data objects
##   OutObs - time_Utc2tb() output
##   RateObs - derivatives of OutObs
##   ParFit - fitable parameters
##   Par - all input parameters
##
## Output:
##   LogLike - log likelihood
##   ModelTrend - Trend model
##   ModelKep - Keplerian model
##   ModelRed - Red noise models
##   model - Full model
##   sigma - total residual RMS
####################################
    if(length(RateObs)>0){
        OutObsNew <- update_OutObs(ParFit,OutObs,RateObs)
    }else{
        OutObsNew <- OutObs
    }
    ParNew <- update_par(Par,ParFit)
#    ParNew <- Par
    t11 <- proc.time()
    out <- update_CombineModel(Data,ParNew,OutObsNew,geometry=Par$geometry,OutTime0=OutTime0)
    dur <- proc.time()-t11
    stars <- unique(Data$star)
    types <- unique(Data$type)
    LogLike <- 0
    ModelTrend <- ModelRed <- ModelKep <- model <- Data
    for(star in stars){
###RV model
        index <- which(Data$star==star & Data$type=='rv')
        if(length(index)>0){
            RvKep <- out$OutRv[[star]]$rv
            ParNew$ObsInfo <- Par$ObsInfo[index,]
            ModelKep[index,2] <- RvKep
            RvRef <- fit_RvRef(OutTime[[star]]$rv,Data[index,],ParFit,ParNew,star)
            ModelTrend[index,2] <- -RvRef
#            RvHat <- RvKep+RvRef
            RvHat <- ((1+RvKep/CMPS)/(1+RvRef/CMPS)-1)*CMPS#1+Zmeas = (1+Zabs)/(1+Zref); Zabs is the absolute observed redshift while Zref is the absolute redshift of the reference spectrum.
            RvRed <- fit_RvRed(OutTimes[[star]]$rv,RvHat,Data[index,],ParFit,ParNew)
            ModelRed[index,2] <- RvRed
            RvHat <- RvHat+RvRed
            model[index,1] <- rowSums(out$OutTime[[star]]$rv$tauE)
            model[index,2] <- RvHat
#cat('names(ParFit)=',names(ParFit),'\n')
#cat('ParFit=',ParFit,'\n')
            llike <- fit_RvLike(RvHat,Data[index,],ParFit,star)
            if(verbose) cat('star',star,';rv;sd(res)=',round(sd(RvHat-Data[index,2]),2),';llike=',llike,';b=',ParFit[grep('bRv',names(ParFit))],'km/s\n')
            LogLike <- LogLike+llike
        }

###absolute astrometry model
        index <- which(Data$type=='abs' & Data$star==star)
        if(length(index)>0){
            AstroKep <- out$OutAbs[[star]]$abs
            ParNew$ObsInfo <- Par$ObsInfo[index,]
            ModelKep[index,c(2,4)] <- AstroKep
            AstroTrend <- fit_AstroTrend(Data[index,],ParFit,ParNew,star=star,type='abs')#rad
            ModelTrend[index,c(2,4)] <- AstroTrend
            AstroHat <- AstroKep+AstroTrend
            AstroRed <- fit_AstroRed(OutTime[[star]]$abs,AstroHat,Data[index,],ParFit,ParNew,star=star,type='abs')#rad
            ModelRed[index,c(2,4)] <- AstroRed
            model[index,1] <- rowSums(out$OutTime[[star]]$abs$tauE)
            AstroHat <- AstroHat+AstroRed
            model[index,c(2,4)] <- AstroHat
            llike <- fit_AstroLike(AstroHat,Data[index,],ParFit,ParNew,star=star,type='abs')
            dastro <- AstroHat-Data[index,c(2,4)]
            if(verbose) cat('star',star,';abs;sd=',round(sd(unlist(dastro))/DMAS2R,2),'mas; llike=',llike,'\n')
            LogLike <- LogLike+llike
        }
    }

###relative astrometry model
    for(star in stars){
        index <- which(Data$type=='rel' & Data$star==star)
        if(length(index)>0){
            AstroKep <- out$OutRel[[star]]$rel
            ModelKep[index,c(2,4)] <- AstroKep
            AstroTrend <- fit_AstroTrend(Data[index,],ParFit,ParNew,star,type='rel')
            ModelTrend[index,c(2,4)] <- AstroTrend
            AstroHat <- AstroKep+AstroTrend
            ParNew$ObsInfo <- Par$ObsInfo[index,]
            AstroRed <- fit_AstroRed(OutTime[[star]]$rel,AstroHat,Data[index,],ParFit,ParNew,star,type='rel')
            ModelRed[index,c(2,4)] <- AstroRed
            AstroHat <- AstroHat+AstroRed
            model[index,1] <- rowSums(out$OutTime[[star]]$rel$tauE)
            model[index,c(2,4)] <- AstroHat
            llike <- fit_AstroLike(AstroHat,Data[index,],ParFit,ParNew,star,type='rel')
            dastro <- AstroHat-Data[index,c(2,4)]
            if(verbose) cat('star',star,';rel;sd=',round(sd(unlist(dastro)),2),'llike=',llike,'\n')
            LogLike <- LogLike+llike
        }
    }
    return(list(llike=LogLike,ModelTrend=ModelTrend,ModelKep=ModelKep,ModelRed=ModelRed,model=model))
}

fit_LogPrior <- function(ParFit,ParMin,ParMax,Par){
####################################
## Calculate log prior
##
## Input:
##   ParFit - values of fitable parameters
##   ParMin - lower boundary of fitable parameters
##   ParMax - upper boundary of fitable parameters
##   ParMax - upper boundary of fitable parameters
##   Par - Other input parameters
##
## Output:
##   logprior - log prior
####################################
    logprior <- 0
##assume uniform distribution
    for(n in names(ParFit)){
        ind <- which(names(ParFit)==n)
        if(length(ind)>0){
            if(Par$prior[ind]=='N'){
                lp <- dnorm(ParFit[n],Par$PriorPar1[ind],Par$PriorPar2[ind],log=TRUE)
            }else if(Par$prior[ind]=='U'){
                lp <- log(1/(ParMax[n]-ParMin[n]))
            }else if(Par$prior[ind]=='logU'){
                lp <- log(1/(log(ParMax[n])-log(ParMin[n]))/ParFit[n])
            }else{
                stop(paste('no',n,'found in Par$prior!'))
            }
            logprior <- logprior+lp
#            if(n=='bRv.alphaCenB.PFS') cat('prior type:',Par$prior[ind],';logprior=',lp,';')
#            cat('n=',n,';prior=',Par$prior[ind],';lp=',lp,'\n')
        }else{
            stop(paste('no',n,'found in Par$Ini!'))
        }
    }
    return(logprior)
}

fit_LogPost <- function(Data,OutObs,RateObs,ParFit,ParMin,ParMax,Par,tem=1,OutTime0=NULL){
####################################
## Calculate log posterior
##
## Input:
##   Data - data frame
##   OutObs - time_Utc2tb() output
##   RateObs - derivatives of OutObs
##   ParFit - fitable parameters
##   tem - tempering parameter (or temperature)
##
## Output:
##   loglike - log likelihood
##   logprior - log prior
##   logpost - log posterior
####################################
    t1 <- proc.time()
#    tmp <- fit_LogLike(Data,OutObs,RateObs,ParFit,Par,OutTime0=OutTime0)
    tmp <- try(fit_LogLike(Data,OutObs,RateObs,ParFit,Par,OutTime0=OutTime0),TRUE)
#    if(class(tmp)=='try-error') cat('ParFit=',paste(ParFit,collapse=','),'\n')
    if(class(tmp)=='try-error') save(Data,OutObs,RateObs,ParFit,Par,OutTime0,file='test.Robj')
    loglike <- tmp$llike
    logprior <- fit_LogPrior(ParFit,ParMin,ParMax,Par)
    logpost <- loglike*tem+logprior
    return(list(loglike=loglike,logprior=logprior,logpost=logpost))
}

fit_proposal <- function(ParOld, ParMin,ParMax,CovNew,PostSample,Sd,eps=1e-16,silent=TRUE){
####################################
## propose a new set of parameters
##
## Input:
##   ParOld - current or old parameter values
##   CovAdapt - adaptive covariance
##
## Output:
##   ParProp - proposed or new parameter values
####################################
    Ntt <- 1e4
    for(k in 1:Ntt){
        tol <- tol0
        ParProp <- try(mvrnorm(n=1,mu=ParOld,Sigma=CovNew,tol=tol0),TRUE)
        if(class(ParProp)=='try-error'){
            if(!silent) cat('Warning: covariance is calculated directly from the MCMC sample!\n')
            ParProp <-try(mvrnorm(n=1,mu=ParOld,Sigma=fit_CovDir(PostSample),tol=tol0),TRUE)
        }
        if(class(ParProp)=='try-error'){
            if(!silent) cat('Warning: covariance is calculated using make.positive.definite()!\n')
            ParProp <- try(mvrnorm(n=1,mu=ParOld,Sigma=make.positive.definite(CovNew)),TRUE)
        }

        if(class(ParProp)=='try-error'){
            if(!silent) cat('Warning: covariance is calculated by adding a small diagnal matrix!\n')
            Sd <- 2.4^2/length(ParOld)
            ParProp <-try(mvrnorm(n=1,mu=ParOld,Sigma=CovNew+eps*diag(length(ParOld)),tol=tol0),TRUE)
        }
        if(class(ParProp)=='try-error'){
        if(!silent) cat('Warning: covariance is calculated using nearPD()!\n')
        ParProp <- try(mvrnorm(n=1,mu=ParOld,Sigma=nearPD(CovNew)),TRUE)
	    }
        if(class(ParProp)=='try-error'){
            if(!silent) cat('Warning: covariance is calculated using make.positive.definite()!\n')
            ParProp <- try(mvrnorm(n=1,mu=ParOld,Sigma=make.positive.definite(CovNew)),TRUE)
        }
        if(class(ParProp)=='try-error'){
            if(!silent) cat('Warning: covariance is calculated using small covariance matrix!\n')
            ParProp <- mvrnorm(n=1,mu=ParOld,Sigma=eps*diag(length(ParOld)),TRUE)
        }
        if(class(ParProp)=='try-error') stop('Error: adatpive covariance is not positive definitive!\n')
        if(all(ParProp>ParMin & ParProp<ParMax)) break()
    }
    if(k==Ntt){
        ind <- which(ParProp<ParMin | ParProp>ParMax)
	cat('problematic parameter:',names(ParProp)[ind],'\n')
	stop('Error: maximum iterations reached and Par is always outside of the boundary!\n')
    }
    return(ParProp)
}

fit_CovDir <- function(mat,eps=1e-16){
####################################
## Calculation covariance directly from MCMC samples
##
## Input:
##   mat - MCMC samples
##   Nupdate - Number of steps for update
##   eps - parameter for a small matrix to be added to the resulting matrix to avoid singularity
##
## Output:
##   cov2 - Covariance of the MCMC posterior sample of parameters
####################################
    Sd <- 2.4^2/ncol(mat)
    k <- nrow(mat)-1
    cmat <- cov(mat)
if(FALSE){
cat('sd(mat[,1])=',sd(mat[,1]),'\n')
cat('sd(mat[,2])=',sd(mat[,2]),'\n')
cat('sd(mat[,11])=',sd(mat[,11]),'\n')
    cat('diag(cmat)=',diag(cmat),'\n')
}
    par.mean <- colSums(mat)/k
    xx <- 0
    for(j in 1:k){
        xx <- xx+mat[j,]%*%t(mat[j,])
    }
yy <- par.mean%*%t(par.mean)
    cmat1 <- 1/k*(xx-(k+1)*yy)
if(FALSE){
cat('dim(xx)=',dim(xx),'\n')
cat('dim(yy)=',dim(yy),'\n')
cat('dim(cmat1)=',dim(cmat1),'\n')
    cat('diag(cmat1)=',diag(cmat1),'\n')
}
    Sd*cmat+Sd*eps*diag(ncol(mat))
}

fit_CovIter <- function(ParFit,cov1,mu1,n,Sd,Nupdate=1,eps=1e-16){
####################################
## Calculation of the covariance matrix through iterations
## It is more efficient than fit_CovDirect
##
## Input:
##   ParFitSample - Posterior sample of fitable parameters
##   cov1 - old covariance
##   mu1 - old mean
##   n - iteration step
##   Nupdate - Number of steps for update
##   eps - parameter for a small matrix to be added to the resulting matrix to avoid singularity
##
## Output:
##   cov2 - Covariance of the MCMC posterior sample of parameters
####################################
    if(n%%Nupdate==0){
        Sd <- 2.4^2/length(ParFit)
        mu2 <- (mu1*(n-1)+ParFit)/n#mu1 and mu2 are the mean vectors of parameters for n-1 and n iterations
        N <- n-1
        cov2 <- (N-1)*cov1/N+Sd/N*(N*mu1%*%t(mu1)-(N+1)*mu2%*%t(mu2)+ParFit%*%t(ParFit)+eps*diag(length(ParFit)))
    }else{
        cov2 <- cov1
    }
    return(cov2)
}

fit_AM <- function(Data,OutObs,RateObs,ParIni,ParMin,ParMax,Par,tem=1, verbose=FALSE,OutTime0=NULL){
####################################
## Modified verion of the delayed rejection adaptive MCMC (DRAM; https://link.springer.com/article/10.1007/s11222-006-9438-0)
## Draw posterior samples using MCMC walkers with adaptive time steps
##
## Input:
##   Data - list of data sets
##   OutObs - time_Utc2tb() output
##   RateObs - derivatives of OutObs
##   ParIni - Initial parameter values
##   Niter - Number of iterations for each chain or chain length
##   Verbose- Print output for diagnosis
##   tem - Tempering parameter (or temperature)
##
## Output:
##   chain - MCMC posterior samples of parameters
##   logpost - Log posterior
##   loglike - Log likelihood
####################################
    Niter <- Par$Niter
    Npar <- length(ParIni)
    CovIni  <- 1e-6*diag(Npar)
#    CovIni  <- 1e-3*diag(Npar)
    chain  <-  array(dim=c(Niter+1,Npar))
    colnames(chain) <- names(ParIni)
    sig <- logpost <-  loglike  <-  rep(NA,Niter+1)
    colnames(chain) <- names(ParIni)
    chain[1,]<- unlist(ParIni)
    mu1 <- chain[1,]
    logpost.out  <- fit_LogPost(Data,OutObs,RateObs,chain[1,],ParMin,ParMax,Par,tem=tem,OutTime0=OutTime0)
    logpost[1] <- logpost.pre <- logpost.out$logpost
    loglike[1] <- loglike.pre <- logpost.out$loglike
    logprior.pre <- logpost.out$logprior
    dt0 <- 0
    CovAdapt <- array(data=0,dim=c(Npar,Npar))
    Sd <- 2.4^2/Npar
    t.start <- proc.time()
    n0 <- 2
    for(i in 1:Niter){
        if(i <= n0){
            CovAdapt <- CovIni
        }else if(i > n0){
            CovAdapt <- fit_CovIter(chain[i,],CovAdapt,mu1,i,Sd)
#            CovAdapt1 <- fit_CovDir(mat=chain[1:n0,])
#	    cat('diag(CovAdapt)=',diag(CovAdapt),'\n')
#	    cat('diag(CovAdapt1)=',diag(CovAdapt1),'\n')
        }
        proposal  <-  fit_proposal(chain[i,],ParMin,ParMax,CovAdapt,PostSample=chain[1:i,],Sd=Sd)
        proprop <- fit_LogPost(Data,OutObs,RateObs,proposal,ParMin,ParMax,Par,tem=tem,OutTime0=OutTime0)
        logpost.prop <- proprop$logpost
	logprior.prop <- proprop$logprior
        loglike.prop <- proprop$loglike

        logpost.cur <- logpost.pre
        loglike.cur <- loglike.pre
	logprior.cur <- logprior.pre
        if(is.na(logpost.prop)){
            cat('Warning: NA values in log posteriors are changed into zeros!\n')
            probab = 0
        }else{
            probab = exp(logpost.prop-logpost.cur)
        }
        if(runif(1)<probab){
            chain[i+1,]=proposal
###values for next proposal estimation
            logpost.pre <- logpost.prop
            loglike.pre <- loglike.prop
	    logprior.pre <- logprior.prop
        }else{
            chain[i+1,]=chain[i,]
            logpost.pre <- logpost.cur
            loglike.pre <- loglike.cur
            logprior.pre <- logprior.cur
        }
##save values
        logpost[i+1] <- logprior.pre + loglike.pre
        loglike[i+1] <- loglike.pre
        mu1 <- (mu1*i+chain[i+1,])/(i+1)
##moniter parameters
        if((i%%(floor(Niter/10))==0 | i==1) & verbose){
#            cat('i=',i,';diag(CovAdapt)=',diag(CovAdapt),';tem=',tem,',acceptance percentage:',fit_acceptance(chain[1:(i+1),]),'\n')
            cat('i=',i,';tem=',tem,',acceptance percentage:',fit_acceptance(chain[1:(i+1),]),'diag(CovAdapt)=',diag(CovAdapt),'\n')
            for(n in colnames(chain)){
                cat(n,'=',chain[i+1,n],';')
            }
            ind <- which.max(loglike[1:(i+1)])
            ParFit <- chain[ind,]
            cat('max(logpost)=',max(logpost[1:(i+1)]),'; maximum likelihood:',max(loglike[1:(i+1)]),'\n\n')
            if(Par$Np>0 & FALSE){
#            if(ParNew$Np>0){
#                dev.new()
                ParNew <- update_par(Par,ParFit)
                pdf('test.pdf',8,4)
                par(mfrow=c(1,2))
                ind <- which(Data$type=='rv')
                yfit <- fit_LogLike(Data,OutObs,RateObs,ParFit,ParNew)$model[ind,2]
                t <- Data[ind,1]
                rv <- Data[ind,2]
                plot(t,rv,xlab='t',ylab='rv',ylim=range(yfit,rv),main=paste('sd=',round(sd(rv-yfit,2),2),'m/s;logLmax',round(max(loglike[1:(i+1)]),2)))
                points(t,yfit,col='red')
                plot(t,rv-yfit,xlab='t',ylab='O-C [m/s]')
                dev.off()
                t.start <- proc.time()
            }
        }
    }
    cbind(chain,logpost,loglike)
}

fit_acceptance <- function(mc){
####################################
## Determine the optimal temperture of hot chains; similar to simulated annealing
##
## Input:
##   mc - mcmc chains
## Output:
##  acceptance - acceptance rate
####################################
    (1-mean(duplicated(mc)))*100
}

fit_multiChain <- function(Data,OutObs,RateObs,ParIni,ParMin,ParMax,Par,tem=1,verbose=FALSE,OutTime0=NULL){
####################################
    ## Parallal adpative MCMC
    ##
    ## Input:
    ##   Data - list of data sets
    ##   ParIni - Initial parameter values
    ##   ParMin - Minimum parameter values
    ##   ParMax - Maximum parameter values
    ##   Par - All input parameters
    ##   AccUp - upper limit of the acceptance rate to choose the optimal temperture
    ##   TemLow - Lower limit of the temperture trials
    ##
    ## Output:
    ##   ParIni - Optimal parameter value from hot chains
    ##   TemOpt - Optimal chain temperture
####################################
    if(Par$Ncore>1){
#        foreach(ncore=1:Par$Ncore,.combine='rbind',.errorhandling='pass') %dopar% {
        p <- foreach(ncore=1:Par$Ncore,.errorhandling='pass') %dopar% {
            fit_AM(Data,OutObs,RateObs,ParIni,ParMin,ParMax,Par,tem=tem,verbose=verbose,OutTime0=OutTime0)
        }
        if(any(grepl('error',class(p[[1]]))))  stop(as.character(p[[1]]))
	lmax <- c()
	for(j in 1:length(p)){
		if(!is.null(dim(p[[j]]))){
        		lmax <- c(lmax,max(p[[j]][,'loglike']))
		}else{
	        	lmax <- c(lmax,-1e10)
		}
	}
	if(length(lmax)==0) cat(p[[1]],'\n')
	ind <- which.max(lmax)
	p[[ind]]
    }else{
        fit_AM(Data,OutObs,RateObs,ParIni,ParMin,ParMax,Par,tem=tem,verbose=verbose,OutTime0=OutTime0)
    }
}

fit_HotChain <- function(Data,OutObs,ParIni,ParMin,ParMax,Par,AccUp=20, TemLow=1e-3,verbose=FALSE,OutTime0=NULL,OffUpdate=FALSE){
####################################
## Determine the optimal temperture of hot chains; similar to simulated annealing
##
## Input:
##   Data - list of data sets
##   ParIni - Initial parameter values
##   ParMin - Minimum parameter values
##   ParMax - Maximum parameter values
##   Par - All input parameters
##   AccUp - upper limit of the acceptance rate to choose the optimal temperture
##   TemLow - Lower limit of the temperture trialss
##   OffUpdate - whether or not update offsets for each run
##
## Output:
##   ParIni - Optimal parameter value from hot chains
##   TemOpt - Optimal chain temperture
####################################
    Niter0 <- Par$Niter
    Npar <- length(ParIni)
    Niter1 <- 1e3
    Par$Niter <- Niter1
####check whether the lowest temperture is appropriate and redermine TemLow
    mcmc <- fit_multiChain(Data,OutObs,RateObs,ParIni,ParMin,ParMax,Par,tem=TemLow,verbose=verbose,OutTime0=OutTime0)
    ParIni <- mcmc[which.max(mcmc[,'loglike']),1:Npar]
    acceptance <- fit_acceptance(mcmc)
    cat('initial acceptance=',acceptance,'\n')
    beta <- 10
#    beta <- 2
    Nbasic <- Niter1
    if(Niter0>=1e4) Nbasic <- 1e4
    if(acceptance<50){
        TemLow <- 1e-3*TemLow
#	beta <- 2
    }

#####Determine the optimal temperture
    Ntry <- 20
    i1 <- 0
    tem <- TemLow
    for(i in 1:Ntry){
        Par$Niter <- Nbasic
        if(acceptance>(AccUp+5) & tem<1){
            tem <- min(tem*beta,1)
        }else if(acceptance<=(AccUp+5) & tem<1){
            tem <- min(tem*2,1)
            if(Niter0>=1e4) Par$Niter <- 1e4
        }
        mcmc <- fit_multiChain(Data=Data,OutObs=OutObs,RateObs=RateObs,ParIni=ParIni,ParMin=ParMin,ParMax=ParMax,Par=Par,tem=tem,verbose=verbose,OutTime0=OutTime0)
        acceptance <-  (1-mean(duplicated(mcmc)))*100

        if((acceptance<AccUp & tem>1e-3) | tem==1) break()
        ParIni <- mcmc[which.max(mcmc[,'loglike']),1:length(ParIni)]
	if(OffUpdate) ParIni <- fit_OptIni(Data,OutObs,RateObs,ParIni,Par)#optimize offsets

#        if(Niter0<1e5){
            cat('Adaptive tempering burning: tem=',format(tem,digit=3),'\n')
            cat('maximum likelihood=',max(mcmc[,'loglike']),'\n')
            cat('acceptance:',acceptance,'\n')
            cat('names(ParIni)=',names(ParIni),'\n')
            cat('ParIni=',ParIni,'\n\n')
#        }

        ind <- which(ParIni<ParMin | ParIni>ParMax)
        if(length(ind)>0) cat('The following parameters out of boundary after using fit_OptIni:',names(ParIni)[ind],'\n')
    }
    if(Niter0>=1e4){
	Par$Niter <- max(min(Niter0/10,1e4),1e4)
    }else{
	Par$Niter <- max(min(Niter0/10,1e4),Niter1)
    }

####Run MCMC with the optimally tempered chain
    cat('\nRun hot chain with tem =',tem,'\n')
    mcmc <- fit_multiChain(Data=Data,OutObs=OutObs,RateObs=RateObs,ParIni=ParIni,ParMin=ParMin,ParMax=ParMax,Par=Par,tem=tem,verbose=verbose,OutTime0=OutTime0)
    Npar <- ncol(mcmc)-2
    ParIni <- mcmc[which.max(mcmc[,'loglike']),1:length(ParIni)]
    cat('OffUpdate=',OffUpdate,'\n')
    if(OffUpdate)  ParIni <- fit_OptIni(Data,OutObs,RateObs,ParIni,Par)#optimize offsets

####Initial constraint of a signal with cold chain; optional
    cat('\nRun short cold chain to constrain the signal\n')
    mcmc <- fit_multiChain(Data=Data,OutObs=OutObs,RateObs=RateObs,ParIni=ParIni,ParMin=ParMin,ParMax=ParMax,Par=Par,tem=1,verbose=verbose,OutTime0=OutTime0)
    ParOpt <- mcmc[which.max(mcmc[,'loglike']),1:Npar]
    LogLikeMax <- max(mcmc[,'loglike'])
    return(list(ParOpt=ParOpt,TemOpt=tem,LogLikeMax=LogLikeMax))
}

fit_CalRes <- function(Data,RateObs,ParOpt, Par){
####################################
## Calculate the residual
##
## Input:
##   Data - list of data sets
##   ParOpt - MAP parameter values
##   Par - All input parameters
##
## Output:
##   fit - output of fit_LogLike
####################################
    fit <- fit_LogLike(Data,OutObs=OutObs,RateObs,ParOpt,Par)
    res <- Data
    res[,c(2,4)] <- Data[,c(2,4)]-fit$model[,c(2,4)]
    c(fit,list(res=res))
}

fit_AddKep <- function(par,kep,Np,kep.name){
####################################
## add keplerian parameters for arbitrary number of planets/companions to a parameter array
##
## Input:
##   par - parameter vector
##   kep - Keplerian parameters
##   Np - number of planets
##   kep.name - names of Keplerian parameters
##
## Output:
##   par - new parameter vector
####################################
    for(j in 1:Np){
        names(kep) <- paste0(kep.name,j)
        par <- c(par,kep)
    }
    par
}

fit_Add1Kep <- function(par,kep,np,kep.name){
####################################
## add one more set of Keplerian parameters to a parameter array
##
## Input:
##   par - parameter vector
##   kep - Keplerian parameters
##   np - current
##   kep.name - names of Keplerian parameters
##
## Output:
##   par - new parameter vector
####################################
    names(kep) <- paste0(kep.name,np+1)
    c(par,kep)
}

fit_OptIni <- function(Data,OutObs,RateObs,ParIni,Par){
####################################
## add one more set of Keplerian parameters to a parameter array
##
## Input:
##   Data - combined data
##   OutObs - output of time_Utc2tb()
##   RateObs - Derivatives of OutObs
##   ParIni - Initial parameter values
##   Par - All input parameters
##
## Output:
##   ParIni - optimized initial values of parameters
####################################
    Par1 <- update_par(Par,ParIni)
    ParIni0 <- ParIni
    tmp <- fit_LogLike(Data,OutObs=OutObs,RateObs,ParIni,Par1)
    model <- tmp$model
#    cat('star:',star,';type:',tps,'\n')
    par.all <- parRAall <- parDECall <- list()
    for(star in Par$stars){
        ind1 <- which(Data[,'star']==star)
        if(length(ind1)>0){
            tps <- unique(Data[ind1,'type'])
            for(tp in tps){
                cat('star:',star,'type:',tp,'\n')
                ind2 <- which(Data[,'star']==star & Data[,'type']==tp)
                if(length(ind2)>0){
                    inss <- unique(Data[ind2,'instrument'])
                    for(ins in inss){
                        ind3 <- which(Data[,'star']==star & Data[,'type']==tp & Data[,'instrument']==ins)
                        if(length(ind3)>0){
                            if(tp=='rv'){
                                n <- paste0('bRv.',star,'.',ins)#reference RV; RVref=RVabs-RVmeas-RVabs*RVmeas/c
                                if(any(names(ParIni)==n)){
                                    ParIni[n] <- ParIni0[n]+mean(model[ind3,2]-Data[ind3,2]-Data[ind3,2]*model[ind3,2]/CMPS)/1e3
                                }else if(any(names(ParIni)=='rvOff')){
                                    ParIni['rvOff'] <- ParIni0['rvOff']+mean(model[ind3,2]-Data[ind3,2]-model[ind3,2]*Data[ind3,2]/CMPS)/1e3
                                }
                            }else{
                                for(coord in c('RA','DEC')){
                                    if(coord=='RA') colum <- 2
                                    if(coord=='DEC') colum <- 4
                                    if(tp=='abs'){
#                                        n <- paste0('b',toupper(coord),'.',star,'.',ins)
                                        n <- paste0(tolower(coord),'Off')
                                    }else{
                                        n <- paste0('b',toupper(coord),'.rel.',ins)
                                    }
                                    if(any(names(ParIni)==n)){
                                        if(tp=='abs'){
                                            res <- (Data[ind3,colum]-model[ind3,colum])/DMAS2R
                                            dt <- (Data[ind3,1]-Par$epoch)/DJY#
#                                            dt <- (Data[ind3,1]-median(Data[ind3,1]))/DJY#Par$epoch
                                            tmp <- lm(res~dt)
                                            cat(n,'=',tmp$coefficients[1],'\n')
                                            par.all[[star]] <- max(min(Par$Max[n],tmp$coefficients[1]),Par$Min[n])
                                            ParIni[n] <- ParIni0[n]+par.all[[star]]
                                            n1 <- paste0('pm',tolower(coord),'Off')
                                            if(any(names(ParIni)==n1)){
                                                if(coord=='RA'){
                                                    dra.dt <- tmp$coefficients[2]*median(cos(Data[ind3,4]))#This is the correct pm offset.
                                                    cat('dpmra=',dra.dt,'\n')
                                       ###  dra.dt <- tmp$coefficients[2]#This is a wrong version
                                                    parRAall[[star]] <- max(min(Par$Max[n],dra.dt),Par$Min[n])
                                                    ParIni[n1] <- ParIni0[n1]+parRAall[[star]]
                                                }else{
                                                    ddec.dt <- tmp$coefficients[2]
                                                    cat('dpmdec=',ddec.dt,'\n')
                                                    parDECall[[star]] <- max(min(Par$Max[n],ddec.dt),Par$Min[n])
                                                    ParIni[n1] <- ParIni0[n1]+parDECall[[star]]
                                                }
                                            }
###barycentric offsets
                                            if(star==Par$stars[length(Par$stars)] & length(Par$stars)>1 & length(par.all)>1){
                                                mC <- exp(Par$logmC1)
                                                mT <- Par$mT
                                                mCT <- mC+mT
                                                ParIni[n] <- (par.all[[Par$star]]*mT+par.all[[Par$companion]]*mC)/mCT
                                                if(any(names(ParIni)==n1)){
                                                    if(coord=='RA') ParIni[n1] <- (parRAall[[Par$star]]*mT+parRAall[[Par$companion]]*mC)/mCT
                                                    if(coord=='DEC') ParIni[n1] <- (parDECall[[Par$star]]*mT+parDECall[[Par$companion]]*mC)/mCT
                                                }
                                            }
                                        }else{
                                            ParIni[n] <- mean(Data[ind3,colum]-model[ind3,colum])
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    ParIni
}


fit_PTAM <- function(Data,OutObs,RateObs,ParIni,ParMin,ParMax,Par,KepIni,KepMin,KepMax,verbose=FALSE,OutTime0=NULL,OffUpdate=FALSE){
####################################
## parallal tempering adaptive MCMC
##
## Input:
##   Data - list of data sets
##   OutObs - time_Utc2tb() output
##   RateObs - derivatives of OutObs
##   ParIni - Initial parameter values
##   ParMin - Minimum parameter values
##   ParMax - Maximum parameter values
##   Par - All input parametersp
##   KepIni - Initial values of Keplerian parameters
##   KepMin - Minimum values of Keplerian parameters
##   KepMax - Maximum values of Keplerian parameters
##
## Output:
##   mcmc - output of the optimal mcmc chains
##   ParOpt - optimal parameters
##   res - output of fit_CalRes()
####################################
    Data0 <- Data
    Ncore <- Par$Ncore
    Niter0 <- Par$Niter
    ParOpt <- ParIni0 <- ParIni <- Par$Ini
    ParMin0 <- ParMin <- Par$Min
    ParMax0 <- ParMax <- Par$Max
    prior <- Par$prior
    PriorPar1 <- Par$PriorPar1
    PriorPar2 <- Par$PriorPar2
    Res <- ParML <- ParMP <- McOpt <- list()
    Np0 <- Par$Np
    Nsig <- Np0
    ll0 <- -1e6#initial baseline model likelihood in case np start from 1
    for(np in Np0:Par$Nmax){
        Par$Np <- np
###Parallel signal identification and constraints
        if(np>0){
            Par$Np <- 1
            ParIni <- fit_Add1Kep(ParIni0,Par$KepIni,0,Par$KepName)
            ParMin <- fit_Add1Kep(ParMin0,Par$KepMin,0,Par$KepName)
            ParMax <- fit_Add1Kep(ParMax0,Par$KepMax,0,Par$KepName)
            Par$prior <- fit_Add1Kep(prior,Par$KepPrior,0,Par$KepName)
            Par$PriorPar1 <- fit_Add1Kep(PriorPar1,Par$KepPar1,0,Par$KepName)
            Par$PriorPar2 <- fit_Add1Kep(PriorPar2,Par$KepPar2,0,Par$KepName)
        }
        Par$Npar <- length(ParIni)
###whether or not to update initial parameters depends on whether updating increase the likelihood
        l1 <- fit_LogLike(Data,OutObs,RateObs,ParIni,Par)$llike
        ParIni2 <- fit_OptIni(Data,OutObs,RateObs,ParIni,Par)#optimize offsets
        l2 <- fit_LogLike(Data,OutObs,RateObs,ParIni2,Par)$llike
        if(l2>l1) ParIni <- ParIni2

###hot chain
        cat('\n Find signal',np,'using hot chain!\n')
###    dlogP <- (KepMax['logP']-KepMin['logP'])/Ncore
###    ParIni['logP1'] <- KepMin['logP']+dlogP*runif(1,ncore-1,ncore)
####initial MCMC constraints of a signal; output posterior sample and par.opt
       tmp <- fit_HotChain(Data,OutObs,ParIni,ParMin,ParMax,Par,verbose=verbose,OutTime0=OutTime0,OffUpdate=OffUpdate)
       ParIni <- tmp$ParOpt

####cold chain combined constraint
        cat('\n Constrian',np,'signals using cold chain!\n')
        mcmc <- list()
        if(Ncore>1){
            mcmc <- foreach(ncore=1:Ncore,.errorhandling='pass') %dopar% {
                Par$Np <- np
                Data <- Data0
                Par$Niter <- Niter0
                Data <- Data0
                if(np>1){
                    KepOpt <- tmp$ParOpt[paste0(Par$KepName,1)]
                    ParIni <- fit_Add1Kep(ParOpt,KepOpt,np-1,Par$KepName)
                    ParMin <- fit_Add1Kep(ParOpt,KepMin,np-1,Par$KepName)
                    ParMax <- fit_Add1Kep(ParOpt,KepMax,np-1,Par$KepName)
                }

                                        #            tmp <- try(fit_AM(Data,RateObs,ParIni,ParMin,ParMax,Par,verbose=verbose),TRUE)
                tmp <- fit_AM(Data,OutObs,RateObs,ParIni,ParMin,ParMax,Par,verbose=verbose,tem=1,OutTime0=OutTime0)
###burnin
                if(class(tmp)!='try-error'){
                    out <- tmp[-(1:floor(nrow(tmp)/2)),]
                    out
                }
            }
        }else{
            Par$Np <- np
            Data <- Data0
            Par$Niter <- Niter0
            Data <- Data0
            if(np>1){
                KepOpt <- tmp$ParOpt[paste0(Par$KepName,1)]
                ParIni <- fit_Add1Kep(ParOpt,KepOpt,np-1,Par$KepName)
                ParMin <- fit_Add1Kep(ParOpt,KepMin,np-1,Par$KepName)
                ParMax <- fit_Add1Kep(ParOpt,KepMax,np-1,Par$KepName)
            }

            tmp <- fit_AM(Data,OutObs,RateObs,ParIni,ParMin,ParMax,Par,verbose=verbose,tem=1,OutTime0=OutTime0)
            mcmc[[1]] <- tmp[-(1:floor(nrow(tmp)/2)),]
        }

###remove bad parallel chains
        ind <- which(sapply(1:length(mcmc),function(k) is.null(dim(mcmc[[k]]))))
        if(length(ind)>0) mcmc <- mcmc[-ind]

###Save mcmc results
        llmax <- sapply(1:length(mcmc), function(i) max(mcmc[[i]][,'loglike']))
        cat('Maximum likelihood for all chains:',llmax,'\n')
        mcopt <- mcmc[[which.max(llmax)]]
        McOpt[[paste0('sig',np)]] <- mcopt

###subtract the best-fit in order to find the next signal
        Par$Npar <- ncol(mcopt)-2
        ParOpt <- mcopt[which.max(mcopt[,'loglike']),1:Par$Npar]
        ParML[[paste0('sig',np)]] <- ParOpt
        ParMP[[paste0('sig',np)]] <- mcopt[which.max(mcopt[,'logpost']),1:Par$Npar]
        tmp <- fit_CalRes(Data,RateObs,ParOpt,Par)
        Res[[paste0('sig',np)]] <- tmp
        if(Par$Np>0){
            Data <- tmp$res
        }

###whether to stop DRAM according to BF criterion
        ll <-  max(mcopt[,'loglike'])
        if(Par$Np>0){
            lnbf3 <- ll-ll0-1.5*log(Par$Nepoch)
            cat('lnbf3=',lnbf3,'\n')
            if(lnbf3<5){
                break()
            }else{
                Nsig <- Nsig+1
            }
        }
        ll0 <- ll
    }
    Par$Niter <- Niter0
    list(McOpt=McOpt,ParOpt=ParOpt,ParML=ParML,ParMP=ParMP,res=Res,Nsig=Nsig)
}


fit_PrepPar <- function(Par){
####################################
## parallal tempering adaptive MCMC
##
## Input:
##   Par - input parameters
##
## Output:
##   ParIni - Initial values of fitable parameters
##   ParMin - Minimum values of fitable parameters
##   ParMax - Maximum values of fitable parameters
####################################
    ParIni <- ParMin <- ParMax <- list()
    ns <- Par$FitName
    for(i in 1:length(ns)){
        var <- as.character(ns[i])
        ParIni[[var]] <- as.numeric(Par$FitIni[i])
        if(Par$FitType[i]=='U'){
            ParMin[[var]] <- Par$Fit1[i]
            ParMax[[var]] <- Par$Fit2[i]
        }else if(Par$FitType[i]=='N'){
            ParMin[[var]] <- Par$Fit1[i]-10*Par$Fit2[i]
            ParMax[[var]] <- Par$Fit1[i]+10*Par$Fit2[i]
        }
        if(grepl('omega|Omega|Mo',var)){
            ParMin[[var]] <- 0
            ParMax[[var]] <- 2*pi
        }
        if(var=='I'){
            ParMin[[var]] <- 0
            ParMax[[var]] <- pi
        }
        if(grepl('jitter',var))   ParMin[[var]] <- 0
        if(var=='e'){
            ParMin[[var]] <- 0
            ParMax[[var]] <- 1
        }
        if(var=='plxOff') ParMin[[var]] <- -Par$plx
    }
    kep.name <- c('logmC','mC','logP','P','e','I','Omega','omegaT','omegaC','Mo','Tp','Tc','Tasc')
    ind <- match(kep.name,ns)
    ind <- ind[!is.na(ind)]
    if(length(ind)>0){
        KepIni <- unlist(ParIni)[ind]
        KepMin <- unlist(ParMin)[ind]
        KepMax <- unlist(ParMax)[ind]
        ParIni <- ParIni[-ind]
        ParMax <- ParMax[-ind]
        ParMin <- ParMin[-ind]
    }else{
        KepIni <- KepMin <- KepMax <- NULL
    }
    return(list(ParIni=ParIni,ParMin=ParMin,ParMax=ParMax,KepIni=KepIni,KepMax=KepMax,KepMin=KepMin))
}

fit_modify <- function(Fit,Par){
####################################
## modify input priors
##
## Input:
##   Fit - prior table
##   Par - input parameters
##
## Output:
##   Par - new Par
####################################
    Fit <- Fit[!grepl('#',Fit[,1]),]
    Fit <- fit_changeType(Fit,c('c','n','n','c'))
    ns <- as.character(Fit[,1])
    par.name <- names(Par)
    fitname <- fittype <- fitini <- fit1 <- fit2 <- c()
    for(i in 1:length(ns)){
        var <- as.character(ns[i])
        if(grepl('logtau|jitter|bRA|bDEC|bDelay|bRv|amp',var)){
##instrument-dependent variables
            for(star in Par$stars){
                if(grepl('Rv',var)){
                    inss <- Par$ins[[star]]$rv
                }else if(grepl('Astro',var)){
                    inss <- Par$ins[[star]]$astro
                }else if(grepl('Delay',var)){
                    inss <- Par$ins[[star]]$delay
                }
                for(ins in inss){
                    ind <- which(par.name==var)
                    if(!grepl('amp',var)){
                        Nrep <- 1
                    }else if(var=='ampARrv'){
                        Nrep <- Par$p[[star]]$rv
                    }else if(var=='ampMArv'){
                        Nrep <- Par$q[[star]]$rv
                    }else if(var=='ampARastro'){
                        Nrep <- Par$p[[star]]$astro
                    }else if(var=='ampMAastro'){
                        Nrep <- Par$q[[star]]$astro
                    }
                    if(Nrep==1){
                        fitname <- c(fitname,paste0(var,'.',star,'.',ins))
                    }else{
                        fitname <- c(fitname,paste0(var,1:Nrep,'.',star,'.',ins))
                    }
                    if(length(ind)>0){
                        fitini <- c(fitini,rep(Par[[ind]],Nrep))
                    }else{
                        cat('Initial value of ',var,'is not found and zero is used!\n')
                        fitini <- c(fitini,rep(0,Nrep))
                    }
                    fit1 <- c(fit1,rep(as.numeric(as.character(Fit[i,2])),Nrep))
                    fit2 <- c(fit2,rep(as.numeric(as.character(Fit[i,3])),Nrep))
                    fittype <- c(fittype,rep(as.character(Fit[i,4]),Nrep))
                }
            }
        }else if(var=='aRv'| var=='aAstro' | var=='aDelay'){
##star-dependent variables
            for(star in Par$stars){
                fitname <- c(fitname,paste0(var,'.',star))
                ind <- which(par.name==var)
                if(length(ind)>0){
                    fitini <- c(fitini,Par[[ind]])
                }else{
                    cat('initial value of ',var,'is not found and zero is used!\n')
                    fitini <- c(fitini,0)
                }
                fit1 <- c(fit1,as.numeric(as.character(Fit[i,2])))
                fit2 <- c(fit2,as.numeric(as.character(Fit[i,3])))
                fittype <- c(fittype,as.character(Fit[i,4]))
            }
        }else{
###universal variables such as Keplerian variables
            fitname <- c(fitname,var)
            ind <- which(par.name==var)
            if(length(ind)>0){
                fitini <- c(fitini,Par[[ind]])
            }else{
                cat('initial value of ',var,'is not found and zero is used!\n')
                fitini <- c(fitini,0)
            }
            fit1 <- c(fit1,as.numeric(as.character(Fit[i,2])))
            fit2 <- c(fit2,as.numeric(as.character(Fit[i,3])))
            fittype <- c(fittype,as.character(Fit[i,4]))
        }
    }
    Par$FitName <- as.character(fitname)
    Par$FitType <- fittype
    Par$Fit1 <- fit1
    Par$Fit2 <- fit2
    Par$FitIni <- fitini
    return(Par)
}

fit_changeType <- function(df,coltype){
####################################
## change data type of dataframe columns
##
## Input:
##   df - data frame
##   coltype - data types, could be character, numeric, or factor
##
## Output:
##   df - new df
####################################
    if(class(df)!='data.frame') df <- data.frame(df)
    cn <- colnames(df)
    for(j in 1:length(coltype)){
        type <- coltype[j]
        if(type=='c'){
            df[[cn[j]]] <- as.character(df[[cn[j]]])
        }else if(type=='n'){
            df[[cn[j]]] <- as.numeric(as.character(df[[cn[j]]]))
        }else if(type=='f'){
            df[[cn[j]]] <- as.factor(df[[cn[j]]])
        }
    }
    df
}
