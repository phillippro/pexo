if(Par$star!='TC'){
###tempo2 results
###using freq
###barycorr.py
rv.up <- 1e6
tempo.type <- 'FB90'
if(DE==430 & tempo.type=='IF99'){
    tf <- read.table('../input/Vtempo0.txt')[,1]
}else if(DE==430 & tempo.type=='FB90'){
    tf <- read.table('../input/Vtempo5.txt')[,1]
}else if(DE==405 & tempo.type=='IF99'){
    tf <- read.table('../input/Vtempo7.txt')[,1]
}else if(DE==405 & tempo.type=='FB90'){
    tf <- read.table('../input/Vtempo3.txt')[,1]
}
#bcpy <- bcpy*IFTE.K
#tf <- read.table('tempo2/Vtempo4.txt')[,1]
#tf <- read.table('tempo2/Ztempo0.txt')[,1]*CMPS
#N <- nrow(utc)/3
#ind1 <- 3*(1:N)-1
numerical <- FALSE
if(nrow(utc)>length(tf)){
    N <- nrow(utc)/3
    ind1 <- 3*(1:N)-1
    numerical <- TRUE
}else{
    ind1 <- 1:nrow(utc)
}
ind <- which(abs(bcpy)<rv.up & abs(tf)<rv.up & RvFull$zcomb$B[ind1]<rv.up)
ind1 <- ind1[ind]
rv.bary <- RvFull$all[ind1]
rv.bary.de <- rv.bary+RvFull$zcomb$gso[ind1]*CMPS-RvFull$zcomb$gso.de[ind1]*CMPS
#plot(jd.utc[ind],tf[ind]-RvFull$zcomb$B[ind]*CMPS-2*tb$zobs[ind]*CMPS,main='T2-pexo(DE)',type='l')
####paper
}
fout <- paste0('../results/paper_RV_',Par$star,'_tempo',tempo.type,'_DE',DE,'_ttt2tdb',TtTdbMethod,'_',Par$RefType,'.pdf')
cat(fout,'\n')
pdf(fout,6,6)
size <- 1.2
par(mar=c(5,5,1,1),cex=size,cex.axis=size,cex.lab=size)
plot(jd.utc[ind1],(RvFull$zcomb$B[ind1]*CMPS-tf[ind])*1e3,xlab=expression(JD[UTC]),ylab=expression(v[bary]*' [mm/s]'),type='l')
#plot(jd.utc[ind1],tf[ind]-RvFull$zcomb$B[ind1]*CMPS,xlab=expression(JD[UTC]),ylab=expression(RV[bary]*'mm/s'),pch=20,cex=0.2)
abline(h=0,lty=2)
dev.off()

####numerical approach
if(numerical & star!='TC'){
    for(k in 1:2){
        if(k==1){
            fout <- paste0('../results/paper_pexoRV_',Par$star,'_tempo',tempo.type,'_DE',DE,'_ttt2tdb',TtTdbMethod,'_dTstep',gsub('\\.','',dTstep),'d_',Par$RefType,'.pdf')
            bjd.tdb.pexo <- cbind(OutBary$tdb[,1],OutBary$tdb[,2]+(roemer.all$tempo$all/IFTE.K)/DAYSEC-dt.shapiroS)#-TropoDelay
        }else{
            fout <- paste0('../results/paper_t2RV_',Par$star,'_tempo',tempo.type,'_DE',DE,'_ttt2tdb',TtTdbMethod,'_dTstep',gsub('\\.','',dTstep),'d_',Par$RefType,'.pdf')
            bjd.tdb.pexo <- cbind(utc[,1],utc[,2]+(tempo2$tt+tempo2$tt2tb+tempo2$roemer-tempo2$shapiro-tempo2$tropo)/DAYSEC)
        }
        cat(fout,'\n')
        pdf(fout,6,6)
        size <- 1.2
        par(mar=c(5,5,1,1),cex=size,cex.axis=size,cex.lab=size)
        if(dTstep==0.01){
            tempo2 <- read.table('../input/10700.out10',header=TRUE)
        }else if(dTstep==0.1){
            tempo2 <- read.table('../input/10700.out11',header=TRUE)
        }
        ##ztempo <- (tempo2$tt[ind1+1]+tempo2$tt2tb[ind1+1]+tempo2$roemer[ind1+1]-tempo2$shapiro[ind1+1]-tempo2$tropo[ind1+1]-(tempo2$tt[ind1-1]+tempo2$tt2tb[ind1-1]+tempo2$roemer[ind1-1]-tempo2$shapiro[ind1-1]-tempo2$tropo[ind1-1]))/DAYSEC/(2*dTstep)
        ztempo <- (tempo2$tt[ind1+1]+tempo2$tt2tb[ind1+1]+tempo2$roemer[ind1+1]-tempo2$shapiro[ind1+1]-(tempo2$tt[ind1-1]+tempo2$tt2tb[ind1-1]+tempo2$roemer[ind1-1]-tempo2$shapiro[ind1-1]))/DAYSEC/(2*dTstep)
        zpexo <- t1mt2(bjd.tdb.pexo[ind1+1,],bjd.tdb.pexo[ind1-1,])/(2*dTstep)-1
        plot(jd.utc[ind1],(zpexo-ztempo)*CMPS*1e3,xlab=expression(JD[UTC]),ylab=expression(v[bary]*' [mm/s]'),type='l')
        abline(h=0,lty=2)
        dev.off()
    }
}

####more
if(Par$star=='TC'){
fout <- paste0('../results/',Par$star,'_tempo',tempo.type,'_DE',DE,'_ttt2tdb',TtTdbMethod,'_',Par$RefType,'.pdf')
cat(fout,'\n')
pdf(fout,12,12)
par(mfrow=c(3,3))
plot(jd.utc[ind1],tf[ind]-RvFull$zcomb$Bt[ind1]*CMPS,main='T2-pexo(tempo)',type='l')
abline(h=0,lty=2)
#z <- (bjd.tcb[ind1+1,2]-bjd.tcb[ind1-1,2])/(tb$tt[ind1+1,2]-tb$tt[ind1-1,2])-1
#plot(jd.utc[ind1],z*CMPS-tf[ind],main='Numercial - T2',type='l')
#abline(h=0,lty=2)
plot(jd.utc[ind1],tf[ind]-RvFull$zcomb$B0[ind1]*CMPS,main='T2-pexo(analytic)',type='l')
abline(h=0,lty=2)
plot(jd.utc[ind1],RvFull$zcomb$lo[ind1]*CMPS,main='zlo*c',type='l')
abline(h=0,lty=2)
plot(jd.utc[ind1],RvFull$zcomb$kpo[ind1]*CMPS,main='zkpo*c',type='l')
abline(h=0,lty=2)
plot(jd.utc[ind1],RvFull$zcomb$kpt[ind1]*CMPS,main='zkpt*c',type='l')
abline(h=0,lty=2)
plot(jd.utc[ind1],RvFull$zcomb$lt[ind1]*CMPS,main='zlt*c',type='l')
abline(h=0,lty=2)
plot(jd.utc[ind1],bcpy[ind]-RvFull$zcomb$B[ind1]*CMPS,main='barycorr-pexo',type='l')
abline(h=0,lty=2)
plot(jd.utc[ind1],rv.bary[ind]-rv.bary.de[ind],main='pexo-pexo.de',type='l')
abline(h=0,lty=2)
plot(jd.utc[ind1],bcpy[ind]-tf[ind],main='barycorr',type='l')
abline(h=0,lty=2)
dev.off()
if(FALSE){
plot(jd.utc,tf*CMPS+(rv.bary-RvFull$gt),main='pexo')
plot(jd.utc,tf*CMPS+(rv.bary+RvFull$gt),main='pexo')
plot(jd.utc,tf*CMPS+(rv.bary-RvFull$ko),main='pexo')
plot(jd.utc,tf*CMPS+(rv.bary+RvFull$ko),main='pexo')
plot(jd.utc,tf*CMPS+(rv.bary-RvFull$kt),main='pexo')
plot(jd.utc,tf*CMPS+(rv.bary+RvFull$kt),main='pexo')
plot(jd.utc,tf*CMPS+(rv.bary-RvFull$lo),main='pexo')
plot(jd.utc,tf*CMPS+(rv.bary+RvFull$lo),main='pexo')
plot(jd.utc,tf*CMPS+(rv.bary-RvFull$lt),main='pexo')
plot(jd.utc,tf*CMPS+(rv.bary+RvFull$lt),main='pexo')
plot(jd.utc,bcpy-tf*CMPS,main='barycorr')

#tem <- read.table('../input/10700v1.out1',header=TRUE)
tem <- read.table('../input/10700.out4',header=TRUE)
dt <- rowSums(tem[,c('tt','tt2tb','roemer')])-rowSums(tem[,c('shapiro','tropo')])
N <- length(dt)/2
djd <- tem[(1:N)*2,'sat']-tem[(1:N)*2-1,'sat']
ztem <- (dt[(1:N)*2]-dt[(1:N)*2-1])/(djd[1]*3600*24)
dbjd <- tem[(1:N)*2,'bat']-tem[(1:N)*2-1,'bat']
ztem2 <- dbjd/djd-1

###compare tempo and pexo
ind <- which(jd.utc-DJM0>tem[1,'sat'])
plot(jd.utc[ind],RvFull$local[ind]+ztem*CMPS,xlab='JD',ylab='RV [m/s]',type='l',main='tempo1 & pexo')
plot(jd.utc[ind],RvFull$local[ind]+ztem2*CMPS,xlab='JD',ylab='RV [m/s]',type='l',main='tempo2 & pexo')


###compare tempo and barycorr
ind <- which(jd.utc-DJM0>tem[1,'sat'])
plot(jd.utc[ind],bcpy[ind]-ztem*CMPS,xlab='JD',ylab='RV [m/s]',type='l',main='tempo1 & BC')
plot(jd.utc[ind],bcpy[ind]-ztem2*CMPS,xlab='JD',ylab='RV [m/s]',type='l',main='tempo2 & BC')

###compare pexo and barycorr
plot(jd.utc,RvFull$local+bcpy,xlab='JD',ylab='RV [m/s]',type='l')
#lines(jd.utc,bcpy,col='red')
}
}
