if(file.exists('mcmc_func.R')){
source('mcmc_func.R')
}else{
source('../code/mcmc_func.R')
}
if(!exists('outf')) outf <- TRUE
#ParNew <- update_par(Par,ParOpt)
fit <- fit_LogLike(Data,OutObs,RateObs,ParOpt,Par,OutTime0=NULL)
model <- fit$model
llmax <- round(fit$llike)

####plot
if(outf){
    fpdf <- paste0('../results/',star,'_',opt$component,'_Nmax',Par$Nmax,'_llmax',llmax,'_N',Par$Niter,'_einstein',Par$Einstein,'.pdf')
    cat(fpdf,'\n')
    pdf(fpdf,16,16)
    par(mfrow=c(4,4))
}
for(star in Par$stars){
    index <- which(Data$type=='rv' & Data$star==star)
    if(length(index)>0){
        tauE <- model[index,1]
        rv <- Data[index,2]-model[index,2]#residual
        erv <- Data[index,3]
	ind <- sort(tauE,index.return=TRUE)$ix
        out <- cbind(tauE,rv,erv)[ind,]
        fout <- paste0('../../dwarfs/bary/data/',star,'_comb.rv')
        cat(fout,'\n')
        write.table(out,file=fout,quote=FALSE,row.names=FALSE,col.names=c('tauE','RV','eRV'))
        if(star=='HD128620'){
            ind <- which((out[,2]-mean(out[,2]))<3*sd(out[,2]) & out[,2]> -15)
        }else{
            ind <- which((out[,2]-mean(out[,2]))<3*sd(out[,2]))
        }
        fcons <- gsub('comb','cons',fout)
        cat(fcons,'\n')
        write.table(out[ind,],file=fcons,quote=FALSE,row.names=FALSE,col.names=c('tauE','RV','eRV'))

        out2 <- wtb(out[,1],out[,2],out[,3],dt=15/24/60)
        fout2 <- gsub('.rv','Bin.rv',fout)
        cat(fout2,'\n')
        write.table(out2,file=fout2,quote=FALSE,row.names=FALSE,col.names=c('tauE','RV','eRV'))
###save residual RVs
        inss <- unique(Data$instrument[index])
        for(instr in inss){
            ind <- index[Data[index,'instrument']==instr]
            if(outf){
                plot(Data[ind,1],Data[ind,2],xlab='jd',ylab='rv',main=paste('rv for',instr,star),xlim=range(model[ind,1],Data[ind,1]),ylim=range(model[ind,2],Data[ind,2]))
#                points(model[ind,1],model[ind,2],col='red')
                points(Data[ind,1],model[ind,2],col='red')
                plot(Data[ind,1],Data[ind,2]-model[ind,2],xlab='jd',ylab='rv',main=paste('rv residual for',instr,star,';sd(res)=',round(sd(Data[ind,2]-model[ind,2]),3)))
            }
            tauE <- model[ind,1]
            rv <- Data[ind,2]-model[ind,2]#residual
            erv <- Data[ind,3]
            out <- cbind(tauE,rv,erv)
            fout <- paste0('../../dwarfs/bary/data/',star,'_',instr,'.rv')
            cat(fout,'\n')
            write.table(out,file=fout,quote=FALSE,row.names=FALSE,col.names=c('tauE','RV','eRV'))
            ind <- which((out[,2]-median(out[,2]))<5*sd(out[,2]))
            if(instr=='HARPS' & FALSE){
                rv0 <- read.table(paste0('../input/',star,'/',star,'_HARPS.rv'),header=TRUE)
                cat('old data:',nrow(rv0),' RVs\n')
                cat('new data:',length(ind),' RVs\n')
                fnew <- paste0('../input/',star,'/',star,'_HARPSnew.rv')
                cat(fnew,'\n')
                write.table(rv0[ind,],file=fnew,quote=FALSE,row.names=FALSE)
            }
            out1 <- out[ind,]
            fout1 <- gsub('.rv','cons.rv',fout)
            cat(fout1,'\n')
            write.table(out1,file=fout1,quote=FALSE,row.names=FALSE,col.names=c('tauE','RV','eRV'))
            out2 <- wtb(out1[,1],out1[,2],out1[,3],dt=15/24/60)
            fout2 <- gsub('.rv','bin.rv',fout)
            cat(fout2,'\n')
            write.table(out2,file=fout2,quote=FALSE,row.names=FALSE,col.names=c('tauE','RV','eRV'))
	    if(star=='HD128621' & instr=='HARPS') source('linear_fit.R')
        }
    }

    inds <- which(Data$type=='abs' & Data$star==star)
    if(length(inds)>0){
        inss <- unique(Data[inds,'instrument'])
        for(instr in inss){
            index <- inds[Data[inds,'instrument']==instr]
            if(outf){
                plot(Data[index,2]*180/pi,Data[index,4]*180/pi,xlab='ra[deg]',ylab='dec[deg]',main=paste(instr,'absolute astrometry for',star),xlim=range(model[index,2]*180/pi,Data[index,2]*180/pi),ylim=range(model[index,4]*180/pi,Data[index,4]*180/pi))
                points(model[index,2]*180/pi,model[index,4]*180/pi,col='red')
                plot((Data[index,2]-model[index,2])*206264.8,(Data[index,4]-model[index,4])*206264.8,xlab='ra[as]',ylab='dec[as]',main=paste(instr,'astrometry residual for',star))
                plot(Data[index,1],(Data[index,2]-model[index,2])*206264.8,xlab='JD',ylab='ra[as]',main=paste(instr,'astrometry residual for',star))
                plot(Data[index,1],(Data[index,4]-model[index,4])*206264.8,xlab='JD',ylab='dec[as]',main=paste(instr,'astrometry residual for',star))
            }
###save residual RVs
            dra <- (Data[index,2]-model[index,2])*cos(model[index,4])/DMAS2R
            ddec <- (Data[index,4]-model[index,4])/DMAS2R
            tauE <- model[index,1]
            out <- cbind(tauE,dra,Data[index,3],ddec,Data[index,5])
            fout <- paste0('../../dwarfs/bary/data/',star,'_',instr,'.abs')
            cat(fout,'\n')
            write.table(out,file=fout,quote=FALSE,row.names=FALSE,col.names=c('tauE','dra','era','ddec','edec'))
            fout1 <- paste0('../../dwarfs/bary/data/',star,'_',instr,'_dRA.abs')
            cat(fout1,'\n')
            write.table(out[,1:3],file=fout1,quote=FALSE,row.names=FALSE,col.names=c('tauE','dra','era'))
            fout2 <- paste0('../../dwarfs/bary/data/',star,'_',instr,'_dDEC.abs')
            cat(fout2,'\n')
            write.table(out[,c(1,4:5)],file=fout2,quote=FALSE,row.names=FALSE,col.names=c('tauE','ddec','edec'))
        }
    }
}

inds <- which(Data$type=='rel')
if(length(inds)>0){
    inss <- unique(Data[inds,'instrument'])
    for(instr in inss){
        index <- inds[Data[inds,'instrument']==instr]
        if(outf){
            plot(Data[index,2],Data[index,4],xlab='ra*',ylab='dec',main=paste(instr,'relative astrometry'),xlim=range(model[index,2],Data[index,2]),ylim=range(model[index,4],Data[index,4]))
            points(model[index,2],model[index,4],col='red')
            plot(Data[index,2]-model[index,2],Data[index,4]-model[index,4],xlab='ra*',ylab='dec',main=paste(instr,'relative astrometry residual'))
        }
        dra <- Data[index,2]-model[index,2]
        ddec <- Data[index,4]-model[index,4]
        tauE <- model[index,1]
        out <- cbind(tauE,dra,Data[index,3],ddec,Data[index,5])
        fout <- paste0('../../dwarfs/bary/data/',star,'_',instr,'.rel')
        cat(fout,'\n')
        write.table(out,file=fout,quote=FALSE,row.names=FALSE,col.names=c('tauE','dra','era','ddec','edec'))
        fout1 <- paste0('../../dwarfs/bary/data/',star,'_',instr,'_dRA.rel')
        cat(fout1,'\n')
        write.table(out[,1:3],file=fout1,quote=FALSE,row.names=FALSE,col.names=c('tauE','dra','era'))
        fout2 <- paste0('../../dwarfs/bary/data/',star,'_',instr,'_dDEC.rel')
        cat(fout2,'\n')
        write.table(out[,c(1,4:5)],file=fout2,quote=FALSE,row.names=FALSE,col.names=c('tauE','ddec','edec'))
    }
}

if(outf){
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
    Par1 <- update_par(Par,ParOpt)
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
###save all variables
    fobj <- gsub('pdf','Robj',fpdf)
    cat(fobj,'\n')
    save(list=ls(all=TRUE),file=fobj)
}


