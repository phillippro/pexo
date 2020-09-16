if(!exists('fit')){
##ParNew <- update_par(Par,ParOpt)
fit <- fit_LogLike(Data,OutObs,RateObs,ParOpt,Par,,OutTime0=OutTime0)
model <- fit$model
llmax <- round(fit$llike)
}
if(!exists('pred')){
##model prediction
tmp <- fit_ModelPredict(Data,Par,ParOpt,Nsim=1e3)
pred <- tmp$pred
}

####plot
if(!exists('fname')){
    fpdf <- paste0('../results/',opt$primary,'_',opt$Companion,'companion_llmax',llmax,'_N',Par$Niter,'_einstein',Par$Einstein,'_fancy.pdf')
}else{
    fpdf <- paste0(fname,'_fancy.pdf')
}
cat('Fancy fitting plots:\n',fpdf,'\n\n')
#pdf(fpdf,16,16)
pdf(fpdf,6,6)
#par(mfrow=c(4,4))
for(star in Par$stars){
    index <- which(Data$type=='rv' & Data$star==star)
    if(length(index)>0){
###save residual RVs
#        cat('RV fit!\n')
        inss <- unique(Data$instrument[index])
	for(instr in inss){
            ind <- index[Data[index,'instrument']==instr]
            ind1 <- which(pred$type=='rv' & pred$star==star & pred$instrument==instr)
            res <- data.frame(Data[ind,1],Data[ind,2]-model[ind,2],Data[ind,3])
            plot_OC1D(raw=Data[ind,],mp=pred[ind1,],res=res,FitType='line',alpha.fit=0.4,alpha.data=0.2,bsize=20,title=instr)
            if(FALSE){
                tauE <- model[ind,1]
                rv <- Data[ind,2]-model[ind,2]#residual
                erv <- Data[ind,3]
                out <- cbind(tauE,rv,erv)
                fout <- paste0('../../dwarfs/bary/data/',star,'_',toupper(instr),'.rv')
                cat(fout,'\n')
                write.table(out,file=fout,quote=FALSE,row.names=FALSE,col.names=c('tauE','RV','eRV'))
            }
	}
    }
    inds <- which(Data$type=='abs' & Data$star==star)
    if(length(inds)>0){
#        cat('absolute astrometry fit!\n')
        inss <- unique(Data[inds,'instrument'])
        for(instr in inss){
            index <- inds[Data[inds,'instrument']==instr]
            ind1 <- which(pred$type=='abs' & pred$star==star & pred$instrument==instr)
            raw <- Data[index,]
            raw[,c(2,4)] <- raw[,c(2,4)]*180/pi
            mp <- pred[ind1,]
            mp[,c(2,4)] <- mp[,c(2,4)]*180/pi
            res <- data.frame(Data[index,1],(Data[index,2]-model[index,2])/DMAS2R,Data[index,3],(Data[index,4]-model[index,4])/DMAS2R,Data[index,5])
            colnames(res) <- colnames(mp)[1:5]
            plot_OC2D(raw=raw,mp=mp,res=res,xlab1='RA [deg]',xlab2='RA residual [mas]',ylab1='DEC [deg]',ylab2='DEC residual [mas]',alpha.fit=0.4,alpha.data=0.4,bsize=20)
            if(FALSE){
                dra <- (Data[index,2]-model[index,2])*cos(model[index,4])/DMAS2R
                ddec <- (Data[index,4]-model[index,4])/DMAS2R
                tauE <- model[index,1]
                out <- cbind(tauE,dra,Data[index,3],ddec,Data[index,5])
                fout <- paste0('../../dwarfs/bary/data/',star,'_',toupper(instr),'.abs')
                cat(fout,'\n')
                write.table(out,file=fout,quote=FALSE,row.names=FALSE,col.names=c('tauE','dra','era','ddec','edec'))
            }
        }
    }

    inds <- which(Data$type=='rel' & Data$star==star)
    if(length(inds)>0){
        inss <- unique(Data[inds,'instrument'])
#        cat('relative astrometry fit!\n')
        for(instr in inss){
            index <- inds[Data[inds,'instrument']==instr]
            ind1 <- which(pred$type=='rel' & pred$star==star & pred$instrument==instr)
            dra <- Data[index,2]-model[index,2]
            ddec <- Data[index,4]-model[index,4]
            res <- data.frame(Data[index,1],dra,Data[index,3],ddec,Data[index,5])
            plot_OC2D(raw=Data[index,],mp=pred[ind1,],res=res,xlab1='RA [deg]',xlab2='RA residual [mas]',ylab1='DEC [deg]',ylab2='DEC residual [mas]',alpha.fit=0.4,alpha.data=0.4,bsize=20)
            if(FALSE){
                dra <- Data[index,2]-model[index,2]
                ddec <- Data[index,4]-model[index,4]
                tauE <- model[index,1]
                out <- cbind(tauE,dra,Data[index,3],ddec,Data[index,5])
                fout <- paste0('../../dwarfs/bary/data/',star,'_',toupper(instr),'.rel')
                cat(fout,'\n')
                write.table(out,file=fout,quote=FALSE,row.names=FALSE,col.names=c('tauE','dra','era','ddec','edec'))
            }
        }
    }
}

if(FALSE){
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
ParNew <- update_par(Par,ParOpt)
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
}

dev.off()
###save all variables
#fobj <- gsub('pdf','Robj',fpdf)
#cat(fobj,'\n')
#save(list=ls(all=TRUE),file=fobj)


