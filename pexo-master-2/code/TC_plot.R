TtTdbMethod <- Par$TtTdbMethod
rv.up <- 1e6

if(DE==430){
    bcpy <- read.table('../input/tc_rv4.txt')[,1]#barycentric velocity in m/s from barycorr.py
}else{
    bcpy <- read.table('../input/tc_rv3.txt')[,1]#barycentric velocity in m/s from barycorr.py
}
tempo.type <- 'FB90'
if(DE==430 & tempo.type=='IF99'){
    tf <- read.table('../input/Vtempo0.txt')[,1]
}else if(DE==430 & tempo.type=='FB90'){
    tf <- read.table('../input/Vtempo5.txt')[,1]#chose by default
}else if(DE==405 & tempo.type=='IF99'){
    tf <- read.table('../input/Vtempo7.txt')[,1]
}else if(DE==405 & tempo.type=='FB90'){
    tf <- read.table('../input/Vtempo3.txt')[,1]
}

if(Par$RVmethod=='analytical'){
    rv.bary <- OutRv$RvLocal
    rv.bary.de <- rv.bary+OutRv$Zcomb$ZgsO*CMPS-OutRv$Zcomb$ZgsO.de*CMPS
                                        #plot(jd.utc[ind],tf[ind]-OutRv$zcomb$ZB[ind]*CMPS-2*P18$zobs[ind]*CMPS,main='T2-pexo(DE)',type='l')
####paper

    fout <- paste0('../results/paper_RV_',Par$star,'_tempo',tempo.type,'_DE',DE,'_ttt2tdb',TtTdbMethod,'_',Par$RefType,'.pdf')
    cat(fout,'\n')
    pdf(fout,6,6)
    size <- 1.2
    jd.utc <- rowSums(utc)
    par(mar=c(5,5,1,1),cex=size,cex.axis=size,cex.lab=size)
    drv <- (OutRv$Zcomb$ZBwe*CMPS-tf)*1e3#mm/s
#    drv <- (OutRv$RvTot-tf)*1e3#mm/s
    #rv.up <- 1#m/s
    rv.up <- 1e5#m/s
    inds <- which(abs(drv)<rv.up)
    ind.rm <- which(abs(drv)>=rv.up)
    if(length(ind.rm)>1){
        cat('warning: more outliers occur in the comparison of the RVs simulated by PEXO and TEMPO2!\n')
    }
    plot(jd.utc[inds],drv[inds],xlab=expression(JD[UTC]),ylab=expression(Delta*v[bary]*' [mm/s]'),type='l')
    ##plot(jd.utc,drv,xlab=expression(JD[UTC]),ylab=expression(v[bary]*' [mm/s]'),type='l')
    ##plot(jd.utc[ind1],tf[ind]-OutRv$zcomb$B[ind1]*CMPS,xlab=expression(JD[UTC]),ylab=expression(RV[bary]*'mm/s'),pch=20,cex=0.2)
    abline(h=0,lty=2)
    dev.off()
}

####numerical approach
if(Par$RVmethod=='numerical'){
    fin <- paste0('../input/10700_from42000to52000by10day_Tstep',gsub('\\.','',Par$Tstep),'.out')
    tempo2 <- read.table(fin,header=TRUE)
#    tempo2 <- read.table('../input/10700.out11',header=TRUE)
    N <- nrow(utc)/3
    ind1 <- 3*(1:N)-1
    ind <- 1:N
    fout <- paste0('../results/paper_pexo_vs_T2_',Par$star,'_tempo',tempo.type,'_DE',DE,'_ttt2tdb',TtTdbMethod,'_Tstep',gsub('\\.','',Par$Tstep),'d_',Par$RefType,'.pdf')
#    bjd.tdb.pexo <- OutTime$BJDtdb
#    bjd.tdb.pexo <- cbind(OutBary$JDtdb[,1],OutBary$JDtdb[,2]+(OutTime$RoemerSolar-OutTime$ShapiroPlanet$Sun)/DAYSEC)
    bjd.tdb.pexo <- cbind(OutBary$JDtdb[,1],OutBary$JDtdb[,2]+OutTime$RoemerSolar/DAYSEC)
#    bjd.tdb.t2 <- cbind(utc[,1],utc[,2]+(tempo2$tt+tempo2$tt2tb+tempo2$roemer-tempo2$shapiro-tempo2$tropo)/DAYSEC)
    bjd.tdb.t2 <- cbind(tempo2$utc[,1],utc[,2]+(tempo2$tt+tempo2$tt2tb)/DAYSEC)
    deinstein <- time_T2mT2(OutBary$JDtdb,OutBary$JDtt)-tempo2$tt2tb
#    droemer <- OutTime$RoemerSolar-tempo2$roemer
    droemer <- OutTime$RoemerT2$all+tempo2$roemer
    dshapiro <- OutTime$ShapiroPlanet$Sun-tempo2$shapiro
#    Deinstein <- time_T2mT2

    cat(fout,'\n')
    pdf(fout,6,6)
    size <- 1.2
    par(mar=c(5,5,1,1),cex=size,cex.axis=size,cex.lab=size)
    ztempo <- (tempo2$tt[ind1+1]+tempo2$tt2tb[ind1+1]+tempo2$roemer[ind1+1]-tempo2$shapiro[ind1+1]-(tempo2$tt[ind1-1]+tempo2$tt2tb[ind1-1]+tempo2$roemer[ind1-1]-tempo2$shapiro[ind1-1]))/DAYSEC/(2*Par$Tstep)
#    dtdb <- time_T2mT2(bjd.tdb.pexo[ind1+1,],bjd.tdb.pexo[ind1-1,])
#    dtropo <- OutTime$TropoDelay[ind1+1]-OutTime$TropoDelay[ind1-1]
    dtt <- time_T2mT2(OutBary$JDtt[ind1+1,],OutBary$JDtt[ind1-1,])
#    zpexo <- (dtdb+dtropo)/dtt-1
#    drv <- (zpexo[inds]-ztempo[inds])*CMPS*1e3
    drv <- (deinstein[ind1+1]+droemer[ind1+1]+dshapiro[ind1+1]-(deinstein[ind1-1]+droemer[ind1-1]+dshapiro[ind1-1]))*CMPS*1e3/dtt
    inds <- which(abs(drv)<1e6)
    plot(jd.utc[inds],drv[inds],xlab=expression(JD[UTC]),ylab=expression(Delta*v[bary]*' [mm/s]'),type='l')
    abline(h=0,lty=2)
    dev.off()
}
