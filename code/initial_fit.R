source('mcmc_func.R')
#if(!exists('ParIni')){
    Data0 <- Data
    Ncore <- Par$Ncore
    Niter0 <- Par$Niter
    ParIni0 <- ParIni <- Par$Ini
    ParMin0 <- ParMin <- Par$Min
    ParMax0 <- ParMax <- Par$Max
    Res <- ParML <- ParMP <- McOpt <- list()
    Nsig <- 0
    Np0 <- Par$Np
    ll0 <- -1e6#initial
    if(Par$Np>0){
        ParIni <- fit_Add1Kep(ParIni0,Par$KepIni,0,Par$KepName)
        ParMin <- fit_Add1Kep(ParMin0,Par$KepMin,0,Par$KepName)
        ParMax <- fit_Add1Kep(ParMax0,Par$KepMax,0,Par$KepName)
    }
    Par$Npar <- length(ParIni)
#}
if(TRUE){
#if(FALSE){
    ParFit <- fit_OptIni(Data,OutObs,RateObs,ParIni,Par)#only optimize offsets
}else{
    ParFit <- ParIni
}
fit <- fit_LogLike(Data,OutObs,RateObs,ParFit,Par,verbose=TRUE,OutTime0=OutTime0)
model <- fit$model
cat('loglike=',fit$llike,'\n')
cat('RMS=',sd(Data[,2]-model[,2]),'\n')
ParNew <- update_par(Par,ParFit)
tmp <- gen_CombineModel(utc,Data,ParNew,component='TAR')
OutObs <- tmp$OutObs
OutTime <- tmp$OutTime
OutAstro <- tmp$OutAstro
OutRv <- tmp$OutRv

####plot
fout <- '../results/initial.pdf'
cat(fout,'\n')
pdf(fout,16,16)
par(mfrow=c(4,4))
for(star in Par$stars){
    index <- which(Data$type=='rv' & Data$star==star)
    if(length(index)>0){
        inss <- unique(Data[index,'instrument'])
        for(ins in inss){
            ind <- index[Data[index,'instrument']==ins]
            plot(Data[ind,1],Data[ind,2],xlab='jd',ylab='rv',main=paste('rv for',star,';',ins),xlim=range(model[ind,1],Data[ind,1]),ylim=range(model[ind,2],Data[ind,2]))
            points(model[ind,1],model[ind,2],col='red')
            res <- Data[ind,2]-model[ind,2]
            plot(Data[ind,1],res,xlab='jd',ylab='rv',main=paste('rv residual for',star,';sd(res)=',round(sd(Data[ind,2]-model[ind,2]),3)))
            if(ins=='HARPS'){
                fout <- paste0(star,'_',ins,'res.dat')
                cat(fout,'\n')
                out <- cbind(model[ind,1],res,model[ind,3])
                write.table(out,file=fout,quote=FALSE,row.names=FALSE,col.names=c('BJD','RV','eRV'))
                fout <- paste0(star,'_',ins,'bin.dat')
                cat(fout,'\n')
                write.table(wtb(out[,1],out[,2],out[,3],dt=15/24/60),file=fout,quote=FALSE,row.names=FALSE,col.names=c('BJD','RV','eRV'))
            }
        }
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
        plot(Data[index,2]-model[index,2],Data[index,4]-model[index,4],xlab='ra*',ylab='dec',main=paste(instr,'relative astrometry residual'))
    }
}

if(Par$star=='alphaCenB' & FALSE){
###reduced  RV
acb <- read.table('/Users/ffeng/Documents/projects/pexo/input/HD128621v1/HD128621_HARPS.rv')
aca <- read.table('/Users/ffeng/Documents/projects/pexo/input/HD128620v1/HD128620_HARPS.rv')
aca[,1] <- aca[,1]+2400000
acb[,1] <- acb[,1]+2400000
tmin <- min(acb[,1])
tmax <- max(acb[,1])
jd1 <- seq(tmin,tmax,length.out=1e3)
tmin <- min(aca[,1])
tmax <- max(aca[,1])
jd2 <- seq(tmin,tmax,length.out=1e3)
utc <- time_ChangeBase(cbind(c(jd1,jd2),0))
Data1 <- cbind(jd1,NA,NA,NA,NA,Par$stars[1],'rv',NA,NA)
Data2 <- cbind(jd1,NA,NA,NA,NA,Par$stars[2],'rv',NA,NA)
DataSim <- rbind(Data1,Data2)
colnames(DataSim) <- c('JD','V1','eV1','V2','eV2','star','rv')
DataSim <- fit_changeType(DataSim,c('n','n','n','n','n','c','c'))
ParNew <- update_par(Par,ParFit)
tmp <- gen_CombineModel(utc,DataSim,ParNew,component='TAR')
OutObs <- tmp$OutObs
OutTime <- tmp$OutTime
OutAstro <- tmp$OutAstro
OutRv <- tmp$OutRv

rvT <- OutTime[[Par$stars[1]]]$rv$RvST-OutTime[[Par$stars[1]]]$rv$RvST[1]+acb[1,2]
rvC <- OutTime[[Par$stars[2]]]$rv$RvST-OutTime[[Par$stars[2]]]$rv$RvST[1]+aca[1,2]

plot(acb[,1],acb[,2],xlab='bjd',ylab='rv',main='alphaCenB')
lines(jd1,rvT,col='red')

plot(aca[,1],aca[,2],xlab='bjd',ylab='rv',main='alphaCenA')
lines(jd1,rvC,col='red')
}

dev.off()


