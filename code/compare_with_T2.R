options(scipen=0)
Npar <- 4
tempo.type <- 'FB90'
original <- TRUE
DE <- Par$DE
if(DE==430 & tempo.type=='IF99'){
    tempo2 <- read.table('../input/10700.out1',header=TRUE)
}else if(DE==430 & tempo.type=='FB90'){
    tempo2 <- read.table('../input/10700.out6',header=TRUE)
}else if(DE==405 & tempo.type=='IF99'){
    tempo2 <- read.table('../input/10700.out8',header=TRUE)
}else if(DE==405 & tempo.type=='FB90'){
    tempo2 <- read.table('../input/10700.out2',header=TRUE)
}
if(Par$pmra!=0 & Par$pmdec!=0 & DE==405) tempo2 <- read.table('../input/10700.out4',header=TRUE)
if(Par$pmra!=0 & Par$pmdec!=0 & DE==430) tempo2 <- read.table('../input/10700.out44',header=TRUE)
tempo2.TEMPO <- read.table('../input/10700.out41',header=TRUE)
tempo2.DE405 <- read.table('../input/10700.out4',header=TRUE)
tempo2.IF99 <- read.table('../input/10700.out42',header=TRUE)
if(Par$pmra!=0 & Par$pmdec!=0 & Par$rv!=0 & DE==430) tempo2 <- read.table('../input/10700.out44',header=TRUE)
if(nrow(utc)==nrow(tempo2)){
    dt <- diff(rowSums(utc))[1]
    Ntime <- nrow(utc)
    fname <- paste0(Par$star,'_RV_',Par$BinaryModel,'_dt',gsub('\\.','',dt),'day_Ntime',Ntime,'unit',Par$Unit,'.pdf')
    dir.out <- '../results/'
    fout1 <- paste0('../results/timing_',Par$star,'_tempo',tempo.type,'_DE',DE,'_ttt2tdb',Par$TtTdbMethod,'_',Par$RefType,'.pdf')
    cat('output pdf:\n',fout1,'\n')
    pdf(fout1,16,16)
    par(mfrow=c(4,4))
    dtt <- time_T2mMjd(OutObs$JDtt,tempo2[,'Ttt'])
    ylim <- c(-1e-5,1e-5)
###initial UTC input
    plot(jd.utc,time_T2mMjd(utc,tempo2$sat),xlab='JD[UTC]',ylab='dutc',type='l')
    abline(h=0,lty=2)
###UTC to TAI
    plot(jd.utc,time_T2mMjd(OutObs$JDtai,tempo2$sat)-tempo2$clock0,xlab='JD[UTC]',ylab='d(tai-utc)',type='l')
    abline(h=0,lty=2)
###TAI to TT
    plot(jd.utc,time_T2mT2(OutObs$JDtt,utc)-tempo2$tt,xlab='JD[UTC]',ylab='d(tt-utc)',type='l')
    abline(h=0,lty=2)
###TT to TB
    plot(jd.utc,time_T2mT2(OutObs$JDtdb,OutObs$JDtt)-(tempo2$tt2tb)/IFTE.K,xlab='JD[UTC]',ylab='d(tdb-tt)',type='l',main=paste('TtTdbMethod:',Par$TtTdbMethod))
    abline(h=0,lty=2)
###Roemer delay; main difference
    plot(jd.utc,OutTime$RoemerSolar/IFTE.K+tempo2$roemer,xlab='JD[UTC]',ylab='Difference in Solar System Roemer delay',type='l',main='PEXO method')
    abline(h=0,lty=2)
###Roemer delay; tempo2 approach
    plot(jd.utc,OutTime$RoemerT2$all/IFTE.K+tempo2$roemer,xlab='JD[UTC]',ylab='d(roemer tempo)',type='l',main='TEMPO method')
    abline(h=0,lty=2)
                                        #abline(v=max(eop[,'MJD'])+DJM0,col='red')
###Shapiro delay
    plot(jd.utc,OutTime$ShapiroSolar-tempo2$shapiro,xlab='JD[UTC]',ylab='d(shapiro)',type='l')
    abline(h=0,lty=2)
###tropospheric delay; main difference
    plot(jd.utc,OutTime$elevation-tempo2$elevation*pi/180,xlab='JD[UTC]',ylab='d(elevation)',type='l',main='PEXO method')
    if(!is.null(OutTime$elevation)){
        plot(jd.utc,OutTime$elevationT2-tempo2$elevation*pi/180,xlab='JD[UTC]',ylab='d(ele tempo)',type='l',main='TEMPO method')
        ind <- which(abs(OutTime$elevationT2*180/pi)>5)
    }else{
        ind <- 1:length(jd.utc)
    }
    plot(jd.utc[ind],OutTime$TropoDelay[ind]-tempo2$tropo[ind],xlab='JD[UTC]',ylab='d(tropo)',type='l',main='PEXO method')
    plot(jd.utc[ind],OutTime$TropoDelayT2[ind]-tempo2$tropo[ind],xlab='JD[UTC]',ylab='d(tropo,tempo)',type='l',main='TEMPO method')
    abline(h=0,lty=2)
    dev.off()

######
###total BJD[TDB] difference
    if(!original){
        dt.t2 <- tempo2$tt+tempo2$tt2tb+tempo2$roemer-tempo2$shapiro-tempo2$tropo#second
    }else{
        dt.t2 <- tempo2$tt+tempo2$tt2tb+tempo2$roemer-(tempo2$shapiro+tempo2$shapiroJ+tempo2$shapiroS+tempo2$shapiroV+tempo2$shapiroU+tempo2$shapiroN)-tempo2$tropo#second
    }
    dt.t2.TEMPO <- tempo2.TEMPO$tt+tempo2.TEMPO$tt2tb+tempo2.TEMPO$roemer-tempo2.TEMPO$shapiro-tempo2.TEMPO$tropo#second
    dt.t2.IF99 <- tempo2.IF99$tt+tempo2.IF99$tt2tb+tempo2.IF99$roemer-tempo2.IF99$shapiro-tempo2.IF99$tropo#second
    dt.t2.DE405 <- tempo2.DE405$tt+tempo2.DE405$tt2tb+tempo2.DE405$roemer-tempo2.DE405$shapiro-tempo2.DE405$tropo#second
    bjd.tdb.t2 <- cbind(utc[,1],utc[,2]+dt.t2/DAYSEC)
    bjd.tdb.t2.TEMPO <- cbind(utc[,1],utc[,2]+dt.t2.TEMPO/DAYSEC)
    bjd.tdb.t2.DE405 <- cbind(utc[,1],utc[,2]+dt.t2.DE405/DAYSEC)
    bjd.tdb.t2.IF99 <- cbind(utc[,1],utc[,2]+dt.t2.IF99/DAYSEC)
    bjd.tdb.pexo1 <- cbind(OutObs$JDtdb[,1],OutObs$JDtdb[,2]-(OutTime$RoemerT2$all+OutTime$ShapiroSolar+OutTime$TropoDelay)/DAYSEC)#tempo-version
    bjd.tdb.pexo2 <- OutTime$BJDtdb#pexo-version
    dt.pexo1 <- time_T2mT2(bjd.tdb.pexo1,utc)
    dt.pexo2 <- time_T2mT2(OutTime$BJDtdb,utc)

###load E10 output
    app <- read.table('../data/TC_E10.txt')[,1]
    ddt.app <- ((app-bjd.tdb.t2[,1])-bjd.tdb.t2[,2])*DAYSEC#second
                                        #
    if(Npar==4){
        idl0 <- unlist(readLines('../input/bjd_fractionOT2.txt'))
    }else if(Npar==2){
        idl0 <- unlist(readLines('../input/bjd_fraction.txt'))
    }
    if(original){
        idl0 <- unlist(readLines('../input/bjd_fraction_original.txt'))
    }else{
        idl0 <- unlist(readLines('../input/bjd_fraction.txt'))
    }
    idl.mod <- unlist(readLines('../input/bjd_fraction.txt'))

    idl <-c()
    tmp <- c()
    j <- 1
    for(i in idl0){
        i1 <- unlist(strsplit(i,' '))
        i2 <- as.numeric(i1[i1!=''])
        if(j==2){
            tmp <- c(tmp,i2)
            idl <- rbind(idl,tmp)
            j <- 0
            tmp <- c()
        }else{
            tmp <- c(tmp,i2)
        }
        j <- j+1
    }


    idl1 <-c()
    tmp <- c()
    j <- 1
    for(i in idl.mod){
        i1 <- unlist(strsplit(i,' '))
        i2 <- as.numeric(i1[i1!=''])
        if(j==2){
            tmp <- c(tmp,i2)
            idl1 <- rbind(idl1,tmp)
            j <- 0
            tmp <- c()
        }else{
            tmp <- c(tmp,i2)
        }
        j <- j+1
    }
    droemer.idl <- idl[,2]#geo_corr
    dtdb.geo.idl <- idl[,3]#clock_corr
    dtdb.obs.idl <- idl[,4]#einstein_corr; second
    dtdb.idl <- dtdb.geo.idl+dtdb.obs.idl
    dshapiro.idl <- idl[,5]#shapiro_corr
    dbjd.idl <- droemer.idl+dtdb.idl+dshapiro.idl

    droemer.idl1 <- idl1[,2]#geo_corr
    dtdb.geo.idl1 <- idl1[,3]#clock_corr
    dtdb.obs.idl1 <- idl1[,4]#einstein_corr; second
    dtdb.idl1 <- dtdb.geo.idl1+dtdb.obs.idl1
    dshapiro.idl1 <- idl1[,5]#shapiro_corr
    dbjd.idl1 <- droemer.idl1+dtdb.idl1+dshapiro.idl1

    ddt.idl <- dbjd.idl-dt.pexo2
    if(original){
        ddt.idl1 <- dbjd.idl1-dt.pexo2
    }else{
        ddt.idl1 <- dbjd.idl1-dt.t2#second; bjd-utc
    }
    if(original & Npar==4){
        fout2 <- paste0('../results/timing_E10original_',Par$star,'_tempo',tempo.type,'_DE',DE,'_ttt2tdb',Par$TtTdbMethod,'_tempo_par',Npar,'_',Par$RefType,'.pdf')
        cat('output pdf:\n',fout2,'\n')
        pdf(fout2,6,6)
        size <- 1.2
        par(mar=c(5,5,2,1),cex=size,cex.axis=size,cex.lab=size)
        plot(jd.utc[ind],ddt.idl[ind]-ddt.idl1[ind],xlab=expression(JD[UTC]),ylab=expression(Delta*BJD[TDB]* ' [s]'),type='l',main='Original E10 - Corrected E10')
        dev.off()
    }

###plot
    ddt1 <- ((bjd.tdb.pexo1[,1]-bjd.tdb.t2[,1])+(bjd.tdb.pexo1[,2]-bjd.tdb.t2[,2]))*DAYSEC
    ddt1.IF99 <- ((bjd.tdb.pexo1[,1]-bjd.tdb.t2.IF99[,1])+(bjd.tdb.pexo1[,2]-bjd.tdb.t2.IF99[,2]))*DAYSEC
    ddt2 <- ((bjd.tdb.pexo2[,1]-bjd.tdb.t2[,1])+(bjd.tdb.pexo2[,2]-bjd.tdb.t2[,2]))*DAYSEC
    if(original){
        ddt1 <- -ddt1
        ddt2 <- -ddt2
    }
    ddt2.TEMPO <- -((bjd.tdb.pexo2[,1]-bjd.tdb.t2.TEMPO[,1])+(bjd.tdb.pexo2[,2]-bjd.tdb.t2.TEMPO[,2]))*DAYSEC
    ddt2.IF99 <- -((bjd.tdb.pexo2[,1]-bjd.tdb.t2.IF99[,1])+(bjd.tdb.pexo2[,2]-bjd.tdb.t2.IF99[,2]))*DAYSEC
    ddt2.DE405 <- -((bjd.tdb.pexo2[,1]-bjd.tdb.t2.DE405[,1])+(bjd.tdb.pexo2[,2]-bjd.tdb.t2.DE405[,2]))*DAYSEC

    for(j in 1:4){
        if(j==1){
            ddt <-ddt1
            fout2 <- paste0('../results/timing_',Par$star,'_tempo',tempo.type,'_DE',DE,'_ttt2tdb',Par$TtTdbMethod,'_tempo_par',Npar,'_',Par$RefType,'.pdf')
        }else if(j==2){
            ddt <-ddt2
            fout2 <- paste0('../results/timing_',Par$star,'_tempo',tempo.type,'_DE',DE,'_ttt2tdb',Par$TtTdbMethod,'_pexo_par',Npar,'_',Par$RefType,'.pdf')
        }else if(j==3){
            fout2 <- paste0('../results/timing_',Par$star,'_tempo',tempo.type,'_DE',DE,'_ttt2tdb',Par$TtTdbMethod,'_e10_par',Npar,'_',Par$RefType,'.pdf')
        }else if(j==4 & Npar==4){
            fout2 <- paste0('../results/pexot_',Par$star,'_tempo',tempo.type,'_DE',DE,'_ttt2tdb',Par$TtTdbMethod,'_e10_par',Npar,'_',Par$RefType,'.pdf')
        }
        fout2 <- gsub('\\.pdf',paste0('_original',original,'\\.pdf'),fout2)
        if((j==4 & Npar==4) | j<4){
            cat('output pdf:\n',fout2,'\n')
            pdf(fout2,6,6)
            size <- 1.2
            par(mar=c(5,5,2,1),cex=size,cex.axis=size,cex.lab=size)
            if(j<3){
                plot(jd.utc[ind],ddt[ind]*1e9,xlab=expression(JD[UTC]),ylab=expression(Delta*BJD[TDB]* ' [ns]'),type='l')
            }else if(j<4){
                if(Npar==4){
                    if(original){
                        plot(jd.utc[ind],ddt.idl[ind],xlab=expression(JD[UTC]),ylab=expression(Delta*BJD[TDB]* ' [s]'),type='l',main='Original E10 - PEXO',col='red')
                                        #                lines(jd.utc[ind],ddt.idl1[ind],col='darkgrey')
                        abline(h=0)
                        legend('topright',xpd=NA,bty='n',lty=1,col=c('black','red','darkgrey'),legend=c('PEXO','Original E10'))
                    }else{
                        plot(jd.utc[ind],ddt.idl[ind]*1e6,xlab=expression(JD[UTC]),ylab=expression(Delta*BJD[TDB]* ' ['*mu*'s]'),type='l',col='red')#ylim=range(ddt.app[ind]*1e3,ddt.idl[ind]*1e3))
                                        #        lines(jd.utc[ind],ddt1[ind]*1e6,col='blue')#tempo-version
                        lines(jd.utc[ind],ddt2[ind]*1e6,col='black')#pexo-version
                        legend('top',inset=-0.1,xpd=NA,horiz=TRUE,bty='n',lty=1,col=c('black','red','blue'),legend=c('PEXO','E10 IDL','TEMPO2'))
                    }
                }else{
                    plot(jd.utc[ind],ddt.app[ind]*1e3,xlab=expression(JD[UTC]),ylab=expression(Delta*BJD[TDB]* ' [ms]'),type='l')#ylim=range(ddt.app[ind]*1e3,ddt.idl[ind]*1e3))
                    lines(jd.utc[ind],ddt.idl[ind]*1e3,col='red')#pexo-version
                    legend('top',inset=-0.1,xpd=NA,horiz=TRUE,bty='n',lty=1,col=c('black','red','blue'),legend=c('E10 applet','E10 IDL','TEMPO2'))
                }
                                        #        lines(jd.utc[ind],ddt.idl[ind]*1e3,col='red')
            }else if(j==4){
                if(!original){
                    plot(jd.utc[ind],ddt2[ind]*1e9,xlab=expression(JD[UTC]),ylab=expression(Delta*BJD[TDB]* ' [ns]'),type='l',col='black')
                    lines(jd.utc[ind],ddt1[ind]*1e9,col='red')
                    legend('top',inset=-0.1,xpd=NA,horiz=TRUE,bty='n',lty=1,col=c('black','red','blue'),legend=c('PEXO','PEXOt','TEMPO2'))
                }else{
                                        #            ylim <- range(ddt2[ind]*1e9,ddt2.DE405[ind]*1e9,ddt2.IF99[ind]*1e9,ddt2.TEMPO[ind]*1e9)
                                        #            ylim <- range(ddt2[ind]*1e9,ddt2.IF99[ind]*1e9)
                    ylim <- range(ddt2[ind]*1e9)
                    plot(jd.utc[ind],ddt2[ind]*1e9,xlab=expression(JD[UTC]),ylab=expression(Delta*BJD[TDB]* ' [ns]'),type='l',col='blue',ylim=ylim,main='Original TEMPO2 - PEXO',)
                                        #            lines(jd.utc[ind],ddt2.DE405[ind]*1e9,col='red')
                                        #            lines(jd.utc[ind],ddt2.IF99[ind]*1e9,col='orange')
                                        #            lines(jd.utc[ind],ddt2.TEMPO[ind]*1e9,col='green')
                    abline(h=0,col='black')
                                        #            legend('top',inset=-0.1,xpd=NA,horiz=TRUE,bty='n',lty=1,col=c('black','blue','red','orange','green'),legend=c('PEXO','TEMPO2 Default','TEMPO2 DE405','TEMPO2 IF99','TEMPO2 TEMPO'))
                    legend('top',xpd=NA,bty='n',lty=1,col=c('black','blue'),legend=c('PEXO','Original TEMPO2'))#Original TEMPO2 with IF99'
                }
            }
            if(!original) abline(h=0,col='blue')
            dev.off()
        }
    }
    write.table(bjd.tdb.pexo2,file=paste0('../results/pexo_DE',DE,'.txt'),quote=FALSE,row.names=FALSE,col.names=FALSE)
###shapiro delays
    ns <- c('Sun','Jupiter','Saturn','Uranus')
    for(n in ns){
        fout2 <- paste0('../results/',n,'_',Par$star,'_tempo',tempo.type,'_DE',DE,'_ttt2tdb',Par$TtTdbMethod,'_par',Npar,'_',Par$RefType,'.pdf')
        cat('output pdf:\n',fout2,'\n')
        pdf(fout2,6,6)
        size <- 1.2
        par(mar=c(5,5,2,1),cex=size,cex.axis=size,cex.lab=size)
        plot(jd.utc,OutTime$ShapiroPlanet[[n]]*1e9,xlab=expression(JD[UTC]),ylab='Shapiro delay [ns]',main=n,type='l')
        abline(h=0,lty=2)
        dev.off()
    }
}
