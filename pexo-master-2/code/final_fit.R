ParNew <- update_par(Par,ParOpt)
#ParFit <- fit_OptIni(Data,OutObs,RateObs,ParOpt,ParNew)
ParFit <- ParOpt
fit <- fit_LogLike(Data,OutObs,RateObs,ParFit,ParNew,OutTime0=OutTime0)
model <- fit$model

####plot
fout <- '../results/final.pdf'
cat(fout,'\n')
pdf(fout,16,16)
par(mfrow=c(4,4))
for(star in Par$stars){
    index <- which(Data$type=='rv' & Data$star==star)
    if(length(index)>0){
        inss <- unique(Data[index,'instrument'])
        for(instr in inss){
            ind <- index[Data[index,'instrument']==instr]
            plot(Data[ind,1],Data[ind,2],xlab='jd',ylab='rv',main=paste('rv for',star),xlim=range(model[ind,1],Data[ind,1]),ylim=range(model[ind,2],Data[ind,2]))
            points(model[ind,1],model[ind,2],col='red')
            plot(Data[ind,1],Data[ind,2]-model[ind,2],xlab='jd',ylab='rv',main=paste('rv residual for',star,';',instr,';sd(res)=',round(sd(Data[ind,2]-model[ind,2]))))
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

###reduced  RV
acb <- read.table('/Users/ffeng/Documents/projects/pexo/input/HD128621v1/HD128621_HARPS.rv')
aca <- read.table('/Users/ffeng/Documents/projects/pexo/input/HD128620v1/HD128620_HARPS.rv')
aca[,1] <- aca[,1]+2400000
acb[,1] <- acb[,1]+2400000
ts <- c(aca[,1],acb[,1])
tmin <- min(ts)
tmax <- max(ts)
jd1 <- seq(tmin,tmax,length.out=1e3)
Data1 <- cbind(jd1,NA,NA,NA,NA,Par$stars[1],'rv',NA,NA)
Data2 <- cbind(jd1,NA,NA,NA,NA,Par$stars[2],'rv',NA,NA)
DataSim <- rbind(Data1,Data2)

utc1 <- time_ChangeBase(cbind(jd1,0))
Par1 <- ParNew
Par1$Nepoch <- nrow(utc1)
ind <- which(Data$type=='rv' & Data$star==Par$star)
p <- t(replicate(nrow(Data1),unlist(Par$ObsInfo[ind[1],])))
Par1$ObsInfo <- fit_changeType(p,coltype=c('n','n','n','n','n','n','c','c','n','n','n','n','n','n','n'))
OutObs1 <- time_Utc2tb(utc1,Par1)
OutTime1 <- time_Ta2te(OutObs1,Par1)
OutT <- rv_FullModel(OutObs1,OutTime1,Par1,component='T')
OutC <- rv_FullModel(OutObs1,OutTime1,Par1,component='C')
rvT <- OutT$RvST-OutT$RvST[1]+acb[1,2]
rvC <- OutC$RvST-OutC$RvST[1]+aca[1,2]

if(Par$star=='alphaCenB'){
plot(acb[,1],acb[,2],xlab='bjd',ylab='rv',main='alphaCenB')
lines(jd1,rvT,col='red')

plot(aca[,1],aca[,2],xlab='bjd',ylab='rv',main='alphaCenA')
lines(jd1,rvC,col='red')
}

dev.off()

