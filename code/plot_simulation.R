pdf('simulation.pdf',16,16)
par(mfrow=c(4,4))
jd <- rowSums(utc)
rvTB <- rowSums(OutT$uOT*OutT$vST)*auyr2kms*1e3
rvTS <- RVT$RvTot
rvCB <- rowSums(OutC$uOT*OutC$vST)*auyr2kms*1e3
rvCS <- RVC$RvTot

rvTD <- RVTd$RvTot
rvCD <- RVCd$RvTot-300

##T simulate
indT <- which(Data$type=='rv' & Data$star==Par$star)
tT <- Data[indT,1]
rvT <- Data[indT,2]-Data[1,2]+rvTS[1]

##S simulate
indC <- which(Data$type=='rv' & Data$star==Par$secondary)
tC <- Data[indC,1]
rvC <- Data[indC,2]-Data[nrow(Data),2]+rvCS[length(rvCS)]

plot(jd,rvTS,xlab='JD',ylab='uOT',ylim=range(rvTS,rvT))
points(tT,rvT,col='red')
points(tT,rvTD[indT],col='green')

plot(jd,rvCS,xlab='JD',ylab='uOC',ylim=range(rvCS,rvC))
points(tC,rvC,col='red')
points(tC,rvCD[indC],col='green')

plot(jd,rvTB,xlab='JD',ylab='uOT',type='l')
plot(jd,rvCB,xlab='JD',ylab='uOC',type='l')

##residual
rvTR <- rvT-rvTD[indT]
rvCR <- rvC-rvCD[indC]
plot(tT,rvTR,xlab='t',ylab='RVres',main=paste0('T;RMS=',round(sd(rvTR),2)))
if(length(rvCR)>0) plot(tC,rvCR,xlab='t',ylab='RVres',main=paste0('C;RMS=',round(sd(rvCR),2)))
dev.off()
