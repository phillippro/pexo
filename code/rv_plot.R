options(scipen=10)
if((Par$star=='TC' | Par$star=='TauCeti' | Par$star=='Tau Ceti' | Par$star=='tau Ceti' | Par$star=='tau ceti') & Par$CompareT2 & nrow(utc)==1001){
    source('TC_plot.R')
}
star <- Par$star
jd.utc <- rowSums(utc)
fname <- paste0(star,'_RV_',Par$BinaryModel,'_Ntime',Par$Nepoch,'_',Par$RefType,'.pdf')
dir.out <- '../results/'
if(!file.exists(dir.out)) system(paste('mkdir',dir.out))
fout <- paste0(dir.out,'paper_',fname)
cat(fout,'\n')
#####pdf plot
pdf(fout,16,12)
size1 <- 1.5
size2 <- 1.2
tS <- OutTime$tS
tpos <- Par$tpos
par(mfrow=c(3,4),mar=c(5,5,2,1),cex.lab=size1,cex.axis=size2,cex=1)
#tt3 <- (rowSums(tS)-sum(tpos))/DJY
tt3 <- time_Jd2yr(tS)
#xlab <- expression((t[a]^{SSB}-t[pos])*' [yr]')
xlab <- 'tS [yr]'
type <- 'o'
np <- 0
if(Par$binary & Par$Np>0){
    np <- np+1
    plot(tt3,OutRv$RvgT,xlab=xlab,ylab=expression(Delta*v[r]*' [m/s]'),main=paste0('P',np,': General relativity in TS'),type=type)
}
np <- np+1
plot(tt3,OutRv$RvsT,xlab=xlab,ylab=expression(Delta*v[r]*' [m/s]'),main=paste0('P',np,': Special relativity in TS'),type=type)
if(Par$binary & Par$Np>0){
    np <- np+1
    plot(tt3,OutRv$RvBT,xlab=xlab,ylab=expression(Delta*v[r]*' [m/s]'),main=paste0('P',np,': Motion of T w.r.t. TSB'),type=type)
    np <- np+1
    plot(tt3,OutRv$RvSB,xlab=xlab,ylab=expression(Delta*v[r]*' [m/s]'),main=paste0('P',np,': Motion of TSB w.r.t. SSB'),type=type)
    np <- np+1
    plot(tt3,OutRv$RvlT,xlab=xlab,ylab=expression(Delta*v[r]*' [m/s]'),main=paste0('P',np,': Lensing in TS'),type=type)
}
np <- np+1
plot(tt3,OutRv$RvlO,xlab=xlab,ylab=expression(Delta*v[r]*' [m/s]'),main=paste0('P',np,': Lensing in SS'),type=type)
np <- np+1
plot(tt3,OutRv$RvgsO,xlab=xlab,ylab=expression(Delta*v[r]*' [m/s]'),main=paste0('P',np,': Relativistic effects in SS'),type=type)
np <- np+1
plot(tt3,OutRv$RvSG,xlab=xlab,ylab=expression(Delta*v[r]*' [m/s]'),main=paste0('P',np,': Motion of geocenter w.r.t. SSB'),type=type)
np <- np+1
plot(tt3,OutRv$RvGO,xlab=xlab,ylab=expression(Delta*v[r]*' [m/s]'),main=paste0('P',np,': Earth rotation'),type=type)#col=tcol('grey',0.2)
np <- np+1
plot(tt3,OutRv$RvSO,xlab=xlab,ylab=expression(Delta*v[r]*' [m/s]'),main=paste0('P',np,': Motion of observer w.r.t. SSB'),type=type)
ind <- which(abs(OutTime$elevation*180/pi)>10)
if(length(ind)>0){
    np <- np+1
    plot(tt3[ind],OutRv$RvTropo[ind],xlab=xlab,ylab=expression(Delta*v[r]*' [m/s]'),main=paste0('P',np,': Troposphere refraction'),type=type)
}

#plot(tt3,OutRv$local,xlab=xlab,ylab=expression(Delta*v[r]*' [m/s]'),main='P12: Total local effect',type=type)
if(Par$binary & Par$Np>0){
#    plot(tt3,OutRv$remote,xlab=xlab,ylab=expression(Delta*v[r]*' [m/s]'),main='P13: Total remote effect',type=type)
    np <- np+1
    plot(tt3,OutRv$RvTot,xlab=xlab,ylab=expression(Delta*v[r]*' [m/s]'),main=paste0('P',np,': All effects'),type=type)
}
dev.off()

pdf('../results/remoteSB.pdf',6,6)
size <- 1.2
par(mar=c(5,5,1,1),cex=size,cex.lab=size,cex.axis=size)
plot(tt3[-1],diff(OutRv$RvSB)/diff(tt3),xlab=xlab,ylab=expression(Delta*dot(v)[r]*' [m/s/year]'),type=type,main='')
dev.off()

adjust.data <- function(t1,y1,t2,y2,factor=1e-3){
    if(max(t1)<2e4) t1 <- t1+4e4
    if(max(t1)<2400000) t1 <- t1+2400000
    t1 <- time_Jd2yr(cbind(t1%/%1,t1%%1))
    orb <- approxfun(t2,y2)
    y1 <- (y1-mean(y1))*factor+mean(orb(t1))
    cbind(t1,y1)
}

