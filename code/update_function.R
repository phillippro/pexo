###This file contains functions used to update the outputs from timing, astrometry and RV modeling.

update_NumDerivObs <- function(Data,ParAll,ParFit,drTel=1,drGeo=1,dvGeo=0.001){
####################################
    ## Derivative of OutObs w.r.t. certain variables
    ##
    ## Input:
    ##   utc - utc time
    ##   ParAll - input parameters
    ##   ParFit - fitable parameters
    ##   drTel - Offset of the geocentric position of telescope (m; )
    ##   drGeo - Offset of the heliocentric position of the geocenter (km; default 5; recommend >=5)
    ##   dvGeo - Offset of the heliocentric velocity of the geocenter (m/s; default 0.001)
    ##
    ## Output:
    ##   tSf - first order coordinate time at SSB
    ##   tBf - first order coordinate time at TSB
####################################
    nn <- names(ParFit)
    rate <- list()
    for(n in nn){
        Par1 <- Par2 <- ParAll
        NumD <- FALSE

        if(n=='xtelOff'){
            dvar <- drTel*1e-3#km
            Par1$xtel <- Par1$xtel+dvar#km
            Par2$xtel <- Par2$xtel-dvar#km
            NumD <- TRUE
        }
        if(n=='ytelOff'){
            dvar <- drTel*1e-3#km
            Par1$ytel <- Par1$ytel+dvar#km
            Par2$ytel <- Par2$ytel-dvar#km
            NumD <- TRUE
        }
        if(n=='ztelOff'){
            dvar <- drTel*1e-3#km
            Par1$ztel <- Par1$ztel+dvar#km
            Par2$ztel <- Par2$ztel-dvar#km
            NumD <- TRUE
        }
        if(grepl('telOff',n)){
            eph <- sofa_Gc2gd(ParAll$n,xyz=c(Par1$xtel,Par1$ytel,Par1$ztel))
            Par1$elong <- eph['elong']*pi/180#rad
            Par1$phi <- eph['phi']*pi/180#rad
            Par1$height <- eph['height']/1e3#km
            eph <- sofa_Gc2gd(ParAll$n,xyz=c(Par2$xtel,Par2$ytel,Par2$ztel))
            Par2$elong <- eph['elong']*pi/180#rad
            Par2$phi <- eph['phi']*pi/180#rad
            Par2$height <- eph['height']/1e3#km
        }

        if(n=='xgeoOff'){
            dvar <- drGeo
            Par1$xgeoOff <- dvar#km
            Par2$xgeoOff <- -dvar#km
            NumD <- TRUE
        }
        if(n=='ygeoOff'){
            dvar <- drGeo
            Par1$ygeoOff <- dvar#km
            Par2$ygeoOff <- -dvar#km
            NumD <- TRUE
        }
        if(n=='zgeoOff'){
            dvar <- drGeo
            Par1$zgeoOff <- dvar#km
            Par2$zgeoOff <- -dvar#km
            NumD <- TRUE
        }
        if(n=='vxgeoOff'){
            dvar <- dvGeo*1e-6
            Par1$vxgeoOff <- dvar#km/s
            Par2$vxgeoOff <- -dvar#km/s
            NumD <- TRUE
        }
        if(n=='vygeoOff'){
            dvar <- dvGeo*1e-6#km/s
            Par1$vygeoOff <- dvar#km/s
            Par2$vygeoOff <- -dvar#km/s
            NumD <- TRUE
        }
        if(n=='vzgeoOff'){
            dvar <- dvGeo*1e-6
            Par1$vzgeoOff <- dvar#km/s
            Par2$vzgeoOff <- -dvar#km/s
            NumD <- TRUE
        }

        if(NumD){
            OutObs1 <- gen_CombineModel(utc,Data,Par1,component='')$OutObs
            OutObs2 <- gen_CombineModel(utc,Data,Par2,component='')$OutObs
            types <- unique(Data$type)
#            stars <- unique(Data$star)
            stars <- ParAll$stars
            for(star in stars){
                for(type in types){
                    obs1 <- OutObs1[[star]][[type]]
                    obs2 <- OutObs2[[star]][[type]]
                    RateObs <- obs1
                    n1 <- names(obs1)
                    for(ii in n1){
                        if(is.list(obs1[[ii]])){
                            n2 <- names(obs1[[ii]])
                            for(jj in n2){
                                RateObs[[ii]][[jj]] <- (obs1[[ii]][[jj]]-obs2[[ii]][[jj]])/(2*dvar)
                            }
                        }else{
                            RateObs[[ii]] <- (obs1[[ii]]-obs2[[ii]])/(2*dvar)
                        }
                    }
                    rate[[n]][[star]][[type]] <- RateObs
                }
            }
        }
    }
    rate
}

update_OutObs  <- function(ParFit,OutObs,RateObs){
####################################
    ## Update OutObs
    ##
    ## Input:
    ##   ParFit - Fitable parameters
    ##   OutObs - current OutObs
    ##   RateObs - derivatives of OutObs
    ##
    ## Output:
    ##   OutObsNew - new OutObs
####################################
    if(length(RateObs)>0){
        stars <- names(RateObs)
        OutObsNew <- OutObs
        for(star in stars){
            types <- names(RateObs[[star]][[type]])
            for(type in types){
                nfit <- names(ParFit)
                nrate <- names(RateObs[[star]][[type]])
                ind.fit <- match(nrate,nfit)
                rate <- list()
                for(kk in 1:length(nrate)){
                    n  <- nrate[kk]
                    rate <- RateObs[[star]][[type]][[n]]
                    n1 <- names(rate)
                    for(ii in n1){
                        if(is.list(OutObs[[star]][[type]][[ii]])){
                            n2 <- names(OutObs1[[ii]])
                            for(jj in n2){
                                OutObsNew[[star]][[type]][[ii]][[jj]] <- OutObs[[star]][[type]][[ii]][[jj]]+ParFit[ind.fit[kk]]*rate[[ii]][[jj]]
                            }
                        }else{
                            OutObsNew[[star]][[type]][[ii]] <- OutObs[[star]][[type]][[ii]]+ParFit[ind.fit[kk]]*rate[[ii]]
                        }
                    }
                }
            }
        }
        return(OutObsNew)
    }else{
        return(OutObs)
    }
}
update_par <- function(ParOld,ParFit){
####################################
## Update parameters
##
## Input:
##   ParOld - Old parameters
##   ParFit - Fitable parameters
##
## Output:
##   ParNew - new parameters
####################################
    ParNew <- ParOld
    for(n in names(ParFit)) ParNew[[n]] <- ParFit[n]
    cn <- names(ParFit)
    if(any(grepl('raOff',cn))) ParNew$ra <- ParOld$ra+ParFit['raOff']*DMAS2R
    if(any(grepl('decOff',cn))) ParNew$dec <- ParOld$dec+ParFit['decOff']*DMAS2R
    if(any(grepl('pmraOff',cn))) ParNew$pmra <- ParOld$pmra+ParFit['pmraOff']
    if(any(grepl('pmdecOff',cn))) ParNew$pmdec <- ParOld$pmdec+ParFit['pmdecOff']
    if(any(grepl('plxOff',cn))) ParNew$plx <- ParOld$plx+ParFit['plxOff']
    if(any(grepl('rvOff',cn))) ParNew$rv <- ParOld$rv+ParFit['rvOff']
    indP <- grep('logP',cn)
    indmC <- grep('logmC',cn)
    if(length(indP)>0){
        for(j in indP){
            ParNew[[gsub('log','',cn[j])]] <- exp(ParFit[cn[j]])
        }
    }
    if(length(indmC)>0){
        for(j in indmC){
            ParNew[[gsub('log','',cn[j])]] <- exp(ParFit[cn[j]])
        }
    }
    ind <- grep('vOff',cn)
    if(length(ind)>0){
        for(j in ind){
            ParNew[[cn[j]]] <- ParOld[[cn[j]]]
        }
    }
    ParNew$pqu <- astro_RaDec2pqu(ParNew$ra,ParNew$dec)
    ParNew
}

update_OutAstro  <- function(ParFit,OutObs,RateOut){
####################################
## Update OutAstro
##
## Input:
##   ParFit - Fitable parameters
##   OutObs - current OutObs
##   RateObs - derivatives of OutObs
##
## Output:
##   OutObsNew - new OutObs
####################################
    if(length(RateObs)>0){
        nfit <- names(ParFit)
        nrate <- names(RateObs)
        ind.fit <- match(nrate,nfit)
        OutObsNew <- OutObs
        rate <- list()
        for(kk in 1:length(nrate)){
            n  <- nrate[kk]
            rate <- RateObs[[n]]
            n1 <- names(rate)
            for(ii in n1){
                if(is.list(OutObs[[ii]])){
                    n2 <- names(OutObs1[[ii]])
                    for(jj in n2){
                        OutObsNew[[ii]][[jj]] <- OutObs[[ii]][[jj]]+ParFit[ind.fit[kk]]*rate[[ii]][[jj]]
                    }
                }else{
                    OutObsNew[[ii]] <- OutObs[[ii]]+ParFit[ind.fit[kk]]*rate[[ii]]
                }
            }
        }
        return(OutObsNew)
    }else{
        return(OutObs)
    }
}

update_CombineModel <- function(Data,Par,OutObs,geometry=TRUE,OutTime0=NULL,TimeUpdate=FALSE){
####################################
## combine all models for different stars and data types
##
## Input:
##   Data - Data frame
##   Par - Input parameters
##   OutObs - Output of time_Utc2tb()
##   geometry - without gravitaional/atmospheric deflection
##
## Output:
##   Out - Model Output
####################################
    ParNew <- Par
    types <- unique(Data$type)
#    stars <- unique(Data$star)
    stars <- Par$stars
    OutTime <- OutTime0
    OutAstro <- OutRel <- OutRv <- OutAbs <- list()
    for(star in stars){
        if(star==Par$companion){
            ParNew <- fit_ChangePar(Par)
        }else{
            ParNew <- Par
        }
        for(type in types){
	    ind <- which(Data$star==star & Data$type==type)
            if(type=='rel') ind <- which(Data$type=='rel')
            if(length(ind)>0){
                ParNew$ObsInfo <- Par$ObsInfo[ind,]
                if(TimeUpdate | is.null(OutTime0)){
                    if(!is.null(OutTime0) | length(OutTime0)>0){
                        OutTime[[star]][[type]] <- time_Ta2te(OutObs[[star]][[type]],ParNew,fit=TRUE,OutTime0=OutTime0[[star]][[type]])
                    }else{
                        OutTime[[star]][[type]] <- time_Ta2te(OutObs[[star]][[type]],ParNew,fit=FALSE)
                    }
                }
                if(type=='rv'){
                    if(Par$binary | Par$Np>0){
                        rvmodel <- rv_FullModel(OutObs[[star]][[type]],OutTime[[star]][[type]],ParNew)

                        if(Par$RvType=='RAW'){
                            OutRv[[star]][[type]] <- rvmodel$RvTot#-Par$rvx*1e3
                        }else{
                            OutRv[[star]][[type]] <- rvmodel$RvBT
                        }
###                       OutRv[[star]][[type]] <- -rv_FullModel(OutObs[[star]][[type]],OutTime[[star]][[type]],ParNew)$Zcomb$ZB*CMPS
                    }else{
#                        OutRv[[star]][[type]] <- rv_FullModel(OutObs[[star]][[type]],OutTime[[star]][[type]],ParNew)$RvTot#-Par$rv*1e3
                        if(Par$RvType=='RAW'){
                            OutRv[[star]][[type]] <- -rv_FullModel(OutObs[[star]][[type]],OutTime[[star]][[type]],ParNew)$Zcomb$ZB*CMPS
                        }else{
                            OutRv[[star]][[type]] <- rep(0,length(ind))
                        }
                    }
                }

                if(type=='abs'){
                    if(!geometry){
                        OutAbs[[star]][[type]] <- astro_FullModel(OutObs[[star]][[type]],OutTime[[star]][[type]],ParNew,Mlens=ParNew$mC1,component='T')$DirObs
                    }else{
                        OutAbs[[star]][[type]] <- gen_Xyz2lb(OutTime[[star]][[type]]$uOT)
                    }
                }
                if(type=='rel'){
                    if(!geometry){
                        OutAstro[[star]][[type]] <- astro_FullModel(OutObs[[star]][[type]],OutTime[[star]][[type]],ParNew,Mlens=ParNew$mC1,component='T')$DirObs
                    }else{
                        OutAstro[[star]][[type]] <- gen_Xyz2lb(OutTime[[star]][[type]]$uOT)
                    }
                }
            }else{
                OutAstro[[star]][[type]] <- OutAbs[[star]][[type]] <- OutRv[[star]][[type]] <- OutTime[[star]][[type]] <- NULL
            }
        }
    }
    for(star in stars){
        ind <- which(Data$star==star)
        types <- unique(Data[ind,]$type)
        if(length(types)>0){
            if(any(types=='rel')){
###	     cat('star=',star,';type=',type,'\n')
                if(star==Par$star){
                    astro1 <- OutAstro[[Par$star]][[type]]
                    astro2 <- OutAstro[[Par$companion]][[type]]
                }else{
                    astro1 <- OutAstro[[Par$companion]][[type]]
                    astro2 <- OutAstro[[Par$star]][[type]]
                }
                astro <- astro1-astro2
                astro[,1] <- astro[,1]*cos((astro1[,2]+astro2[,2])/2)
                OutRel[[star]][[type]] <- astro/DMAS2R#mas
            }
        }
    }
    return(list(OutObs=OutObs,OutTime=OutTime,OutAbs=OutAbs,OutRel=OutRel,OutAstro=OutAstro,OutRv=OutRv))
}
