if(!exists('ParIni')){
    ParAll <- Par
    for(n in names(Par$Ini)) ParAll[[n]] <- ParAll$Ini[[n]]
    ParAll$Npar <- length(ParAll$Ini)
###parallax tempering adaptive MCMC (PTAM)
###numerical derivation of OutObs objects w.r.t. observatory site and ephemerides
    OutObs <- time_Utc2tb(utc,Par)
    if(!exists('RateObs')){
        RateObs <- update_NumDerivObs(utc,ParAll,Par$Ini)
    }
    if(Par$Niter>1e6){
        verbose <- FALSE
    }else{
        verbose <- TRUE
    }
    Data0 <- Data
    Ncore <- ParAll$Ncore
    Niter0 <- ParAll$Niter
    ParIni0 <- ParIni <- ParAll$Ini
    ParMin0 <- ParMin <- ParAll$Min
    ParMax0 <- ParMax <- ParAll$Max
    Res <- ParML <- ParMP <- McOpt <- list()
    Nsig <- 0
    Np0 <- ParAll$Np
    ll0 <- -1e6#initial
    ParIni <- fit_Add1Kep(ParIni0,ParAll$KepIni,0,ParAll$KepName)
    ParMin <- fit_Add1Kep(ParMin0,ParAll$KepMin,0,ParAll$KepName)
    ParMax <- fit_Add1Kep(ParMax0,ParAll$KepMax,0,ParAll$KepName)
    ParAll$Npar <- length(ParIni)
}
ParIni <- fit_OptIni(Data,ParIni,Par)#only optimize offsets
ParNew <- update_par(Par,ParIni)
tmp <- fit_LogLike(Data,OutObs,RateObs,ParIni,ParNew)
model <- tmp$model

#model <- tmp$modeln
#model <- tmp$ModelKep
fout <- '../results/initial.pdf'
cat(fout,'\n')
pdf(fout,16,16)
par(mfrow=c(4,4))
par(mfrow=c(4,4))
for(star in Par$stars){
    index <- which(Data$type=='rv' & Data$star==star)
    if(length(index)>0){
        plot(Data[index,1],Data[index,2],xlab='jd',ylab='rv',main=paste('rv for',star),xlim=range(model[index,1],Data[index,1]),ylim=range(model[index,2],Data[index,2]))
        points(model[index,1],model[index,2],col='red')
        plot(Data[index,1],Data[index,2]-model[index,2],xlab='jd',ylab='rv',main=paste('rv residual for',star,';sd(res)=',round(sd(Data[index,2]-model[index,2]))))
    }

    inds <- which(Data$type=='abs' & Data$star==star)
    if(length(inds)>0){
        inss <- unique(Data[inds,'instrument'])
        for(instr in inss){
            index <- inds[Data[inds,'instrument']==instr]
            plot(Data[index,2]*180/pi,Data[index,4]*180/pi,xlab='ra[deg]',ylab='dec[deg]',main=paste(instr,'absolute astrometry for',star),xlim=range(model[index,2]*180/pi,Data[index,2]*180/pi),ylim=range(model[index,4]*180/pi,Data[index,4]*180/pi))
            points(model[index,2]*180/pi,model[index,4]*180/pi,col='red')
            plot((Data[index,2]-model[index,2])*206264.8,(Data[index,4]-model[index,4])*206264.8,xlab='ra[as]',ylab='dec[as]',main=paste(instr,'astrometry residual for',star))
            plot(Data[index,1],(Data[index,2]-model[index,2])*206264.8,xlab='JD',ylab='ra[as]',main=paste(instr,'astrometry residual for',star))
            plot(Data[index,1],(Data[index,4]-model[index,4])*206264.8,xlab='JD',ylab='dec[as]',main=paste(instr,'astrometry residual for',star))
        }
    }
}

inds <- which(Data$type=='rel')
if(length(inds)>0){
    inss <- unique(Data[inds,'instrument'])
    for(instr in inss){
        index <- inds[Data[inds,'instrument']==instr]
        plot(Data[index,2],Data[index,4],xlab='ra*',ylab='dec',main=paste(instr,'relative astrometry'),xlim=range(model[index,2],Data[index,2]),ylim=range(model[index,4],Data[index,4]))
        points(model[index,2],model[index,4],col='red')
                                        #bt <- gen_CalOffset(OutTime$uOT+OutTime$rBT/Par$mC*Par$mTC*Par$plx*DMAS2R,OutTime$uOT)*1e3#mas
                                        #points(bt[,1],bt[,2],col='red')
        plot(Data[index,2]-model[index,2],Data[index,4]-model[index,4],xlab='ra*',ylab='dec',main=paste(instr,'relative astrometry residual'))
    }
}
dev.off()
