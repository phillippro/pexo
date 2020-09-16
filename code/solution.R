##initial covariance matrix
###launch MDRAM
for(n in names(Par$Ini)) Par[[n]] <- Par$Ini[[n]]
Par$Npar <- length(Par$Ini)
###parallax tempering adaptive MCMC (PTAM)
###numerical derivation of OutObs objects w.r.t. observatory site and ephemerides
RateObs <- update_NumDerivObs(utc,Par,Par$Ini)
ParNew <- Par
if(!all(Par$Ini==0)) ParNew <- update_par(Par,Par$Ini)
mcmc <- fit_PTAM(Data=Data,OutObs=OutObs,RateObs=RateObs,ParIni=Par$Ini,ParMin=Par$Min,ParMax=Par$Max,Par=ParNew,KepIni=Par$KepIni,KepMin=Par$KepMin,KepMax=Par$KepMax,verbose=opt$verbose,OutTime0=OutTime0,OffUpdate=FALSE)
#mcmc <- fit_PTAM(Data=Data,OutObs=OutObs,RateObs=RateObs,ParIni=Par$Ini,ParMin=Par$Min,ParMax=Par$Max,Par=Par,KepIni=Par$KepIni,KepMin=Par$KepMin,KepMax=Par$KepMax,verbose=verbose,OutTime0=OutTime0,OffUpdate=TRUE)
#Nsig <- fit$Nsig-1
#Nsig <- max(mcmc$Nsig,0)
#Nsig <- min(mcmc$Nsig-1,1)
Nsig <- Par$Np
mc <- mcmc$McOpt[[paste0('sig',Nsig)]]
ll <- mc[,'loglike']
lp <- mc[,'logpost']
Npar <- ncol(mc)-2
indl <- which.max(ll)
indp <- which.max(lp)
cat('logPmax=',lp[indp],'\n')
cat('logLmax=',ll[indl],'\n')
ParOpt <- mc[indl,1:Npar]
ParStat <- c()
for(j in 1:Npar){
    ParStat <- cbind(ParStat,data.distr(mc[,j],lp=lp,plotf=FALSE))
}
colnames(ParStat) <- colnames(mc)[1:Npar]

###update outputs
OutObs <- update_OutObs(ParOpt,OutObs,RateObs)
tmp <- gen_CombineModel(utc,Data,Par,component=Par$component)
OutObs <- tmp$OutObs
OutTime <- tmp$OutTime
OutAstro <- tmp$OutAstro
OutRv <- tmp$OutRv

###output residuals
#source('get_residual.R')

###add fit results as output variables
#OutRv$RvRes <- RvRes
for(n in names(mcmc)) OutRv[[n]] <- mcmc[[n]]

if(Par$unit=='TDB'){
###change unit
    source('output_tcb2tdb.R')
}

