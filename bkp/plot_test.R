####plot
fout <- '../results/test.pdf'
cat(fout,'\n')
pdf(fout,16,16)
par(mfrow=c(4,4))
ZgsO <- c(OutRv$alphaCenB$rv$Zcomb$ZgsO,OutRv$alphaCenBC$rv$Zcomb$ZgsO)
ZgsT <- c(OutRv$alphaCenB$rv$Zcomb$ZgsT,OutRv$alphaCenBC$rv$Zcomb$ZgsT)
ZSO <- c(OutRv$alphaCenB$rv$Zcomb$ZSO,OutRv$alphaCenBC$rv$Zcomb$ZSO)
ZST <- c(OutRv$alphaCenB$rv$Zcomb$ZST,OutRv$alphaCenBC$rv$Zcomb$ZST)
ZST0 <- c(OutRv$alphaCenB$rv$Zcomb$ZST0,OutRv$alphaCenBC$rv$Zcomb$ZST0)
ZlO <- c(OutRv$alphaCenB$rv$Zcomb$ZlO,OutRv$alphaCenBC$rv$Zcomb$ZlO)
ZlT <- c(OutRv$alphaCenB$rv$Zcomb$ZlT,OutRv$alphaCenBC$rv$Zcomb$ZlT)
Ztropo <- c(OutRv$alphaCenB$rv$Zcomb$Ztropo,OutRv$alphaCenBC$rv$Zcomb$Ztropo)
#Z <- (1-ZgsO)/(1-ZgsT)*(1+ZST-ZlT)/(1+ZSO-ZlO-Ztropo)-1
Z <- (1-ZgsO)/(1-ZgsT)*(1+ZST0-ZlT)/(1+ZSO-ZlO-Ztropo)-1
#Z <- (1+ZST0)/(1+ZSO)-1
RV <- Z*CMPS
for(star in Par$stars){
    index <- which(Data$type=='rv' & Data$star==star)
    if(length(index)>0){
        rv <- RV[index]-RV[index[1]]+Data[index[1],2]
#        plot(Data[index,1],Data[index,2],xlab='jd',ylab='rv',main=paste('rv for',star),xlim=range(model[index,1],Data[index,1]),ylim=range(model[index,2],Data[index,2]))
        plot(Data[index,1],Data[index,2],xlab='jd',ylab='rv',main=paste('rv for',star),xlim=range(model[index,1],Data[index,1]),ylim=range(model[index,2],Data[index,2]))
        points(model[index,1],rv,col='red')
#        plot(Data[index,1],Data[index,2]-model[index,2],xlab='jd',ylab='rv',main=paste('rv residual for',star,';sd(res)=',round(sd(Data[index,2]-model[index,2]))))
        plot(Data[index,1],Data[index,2]-rv,xlab='jd',ylab='rv',main=paste('rv residual for',star,';sd(res)=',round(sd(Data[index,2]-rv))))
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

rvT <- OutRv1[[Par$stars[1]]]$rv$RvST-OutRv1[[Par$stars[1]]]$rv$RvST[1]+acb[1,2]
rvC <- OutRv1[[Par$stars[2]]]$rv$RvST-OutRv1[[Par$stars[2]]]$rv$RvST[1]+aca[1,2]

index <- which(Data$type=='rv' & Data$star==Par$stars[1])[1]
rv0 <- OutRv1[[Par$stars[1]]]$rv$RvTot
ind <- which.min(abs(rowSums(utc.new)-Data[index,1]))
rvOT <- rv0-rv0[ind]+Data[index,2]

index <- which(Data$type=='rv' & Data$star==Par$stars[2])[1]
if(!all(is.na(index))){
    rv0 <- OutRv1[[Par$stars[2]]]$rv$RvTot
    ind <- which.min(abs(rowSums(utc.new)-Data[index,1]))
    rvOC <- rv0-rv0[ind]+Data[index,2]
}else{
    rvOC <- NULL
}

astroT <- gen_Xyz2lb(OutTime1[[Par$stars[1]]]$rv$uOT)
astroC <- gen_Xyz2lb(OutTime1[[Par$stars[2]]]$rv$uOT)
astroC1 <- gen_Xyz2lb(OutTime1[[Par$stars[1]]]$rv$uOC)
AstroRel <- astroC-astroT
AstroRel1 <- astroC1-astroT
AstroRel[,1] <- AstroRel[,1]*cos(astroT[,2])
AstroRel1[,1] <- AstroRel1[,1]*cos(astroT[,2])

plot(AstroRel[,1]*pc2au,AstroRel[,2]*pc2au,xlab='ra',ylab='dec',main='relative astrometry')
points(AstroRel1[,1]*pc2au,AstroRel1[,2]*pc2au,col='red',pch='.')

plot(acb[,1],acb[,2],xlab='bjd',ylab='rv',main='alphaCenB')
lines(jd1,rvT,col='red')
points(acbP[,1],acbP[,2],col='blue')

plot(aca[,1],aca[,2],xlab='bjd',ylab='rv',main='alphaCenA')
lines(jd2,rvC,col='red')
points(acaP[,1],acaP[,2],col='blue')

####DATA + SIMULATIN
for(star in Par$stars){
    if(!is.null(rvOC)){
        ind1 <- which(Data$type=='rv' & Data$star==star)
        ind2 <- which(DataSim$type=='rv' & DataSim$star==star)
        if(star==Par$stars[1]){
            rv <- rvOT
        }else if(!is.null(rvOC)){
            rv <- rvOC
        }else{
            rv <- NULL
        }
        if(length(ind1)>0 & !is.null(rv)){
            plot(Data[ind1,1],Data[ind1,2],xlab='jd',ylab='rv',main=paste('rv for',star),xlim=range(2458000,max(Data[ind1,1])),ylim=range(15000,max(Data[ind1,2])))
            lines(rowSums(utc.new[ind2,]),rv,col='red')
        }
    }
}
if(!is.null(rvOC)){
    plot(rowSums(utc.new[1:1e3,]),rvOC-rvOT,xlab='JD',ylab='drv')
}
dev.off()
