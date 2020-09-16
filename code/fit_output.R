if(file.exists('mcmc_func.R')){
    source('mcmc_func.R')
}else{
    source('../code/mcmc_func.R')
}
if(!exists('outf')) outf <- TRUE
#if(!exists('outf')) outf <- FALSE
fit <- fit_LogLike(Data,OutObs,RateObs=RateObs,ParFit=ParOpt,Par=Par,OutTime0=OutTime0,TimeUpdate=TRUE)
model <- fit$model
llmax <- round(fit$llike)
##model prediction
tmp <- fit_ModelPredict(Data,Par,ParOpt)
pred <- tmp$pred

dir.out <- '../results/'
if(!dir.exists(dir.out)) system(paste('mkdir',dir.out))
####plot
fname <- paste0(dir.out,opt$primary,'_',opt$Companion,'companion_llmax',llmax,'_N',Par$Niter,'_einstein',Par$Einstein)
fpdf <- paste0(fname,'.pdf')
pdf(fpdf,16,16)
par(mfrow=c(4,4))
for(star in Par$stars){
    index <- which(Data$type=='rv' & Data$star==star)
    if(length(index)>0){
        tauE <- model[index,1]
        rv <- Data[index,2]-model[index,2]#residual
        erv <- Data[index,3]
        ins <- Data[index,'instrument']
	ind <- sort(tauE,index.return=TRUE)$ix
        if(outf){
            out <- data.frame(tauE,rv,erv,ins)[ind,]
            fout <- paste0(dir.out,star,'_barycorrected.rv')
            cat('Barycentrically corrected RV file:\n',fout,'\n\n')
            write.table(out,file=fout,quote=FALSE,row.names=FALSE,col.names=c('tauE','RV','eRV','Instrument'))
            if(star=='HD128620'){
                ind <- which((out[,2]-mean(out[,2]))<3*sd(out[,2]) & out[,2]> -15)
            }else{
                ind <- which((out[,2]-mean(out[,2]))<3*sd(out[,2]))
            }
        }
###save residual RVs
        inss <- unique(Data$instrument[index])
        for(instr in inss){
            ind <- index[Data[index,'instrument']==instr]
            plot(Data[ind,1],Data[ind,2],xlab='jd',ylab='RV',main=paste('RV for',instr,star),xlim=range(model[ind,1],Data[ind,1]),ylim=range(model[ind,2],Data[ind,2]))
            points(Data[ind,1],model[ind,2],col='red')
            plot(Data[ind,1],Data[ind,2]-model[ind,2],xlab='JD',ylab='RV [m/s]',main=paste('rv residual for',instr,star,';sd(res)=',round(sd(Data[ind,2]-model[ind,2]),3)))
            tauE <- model[ind,1]
            rv <- Data[ind,2]-model[ind,2]#residual
            erv <- Data[ind,3]
            out <- cbind(tauE,rv,erv)
            ind <- which((out[,2]-median(out[,2]))<5*sd(out[,2]))
            out1 <- out[ind,]
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
            if(outf){
                ##save residual RVs
                dra <- (Data[index,2]-model[index,2])*cos(model[index,4])/DMAS2R
                ddec <- (Data[index,4]-model[index,4])/DMAS2R
                tauE <- model[index,1]
                out <- cbind(tauE,dra,Data[index,3],ddec,Data[index,5])
                fout <- paste0(dir.out,star,'_',instr,'.abs')
                cat(fout,'\n\n')
                write.table(out,file=fout,quote=FALSE,row.names=FALSE,col.names=c('tauE','dra','era','ddec','edec'))
            }
        }
    }
}

inds <- which(Data$type=='rel')
if(length(inds)>0 & outf){
    inss <- unique(Data[inds,'instrument'])
    for(instr in inss){
        index <- inds[Data[inds,'instrument']==instr]
        plot(Data[index,2],Data[index,4],xlab='ra*',ylab='dec',main=paste(instr,'relative astrometry'),xlim=range(model[index,2],Data[index,2]),ylim=range(model[index,4],Data[index,4]))
        points(model[index,2],model[index,4],col='red')
        plot(Data[index,2]-model[index,2],Data[index,4]-model[index,4],xlab='ra*',ylab='dec',main=paste(instr,'relative astrometry residual'))
        dra <- Data[index,2]-model[index,2]
        ddec <- Data[index,4]-model[index,4]
        tauE <- model[index,1]
        out <- cbind(tauE,dra,Data[index,3],ddec,Data[index,5])
        fout <- paste0(dir.out,star,'_',instr,'.rel')
        cat(fout,'\n\n')
        write.table(out,file=fout,quote=FALSE,row.names=FALSE,col.names=c('tauE','dra','era','ddec','edec'))
    }
}

####trace and posterior plots
if(any(Par$types=='rv')) source('mcmc_analysis.R')
cat('Summary plots:\n',fpdf,'\n\n')
dev.off()

####more fancy plots
source('sim_fit.R')

###save all outputs as R objects
fobj <- paste0(fname,'.Robj')
if(grepl('Robj',opt$out)){
    fobj <- opt$out
}
cat('R object with all variables:\n')
cat(fobj,'\n\n')
save(list=ls(all=TRUE),file=fobj)

# Save to CSV
txt_path <- gsub('.Robj', '_ParStat.txt', fobj)
write.table(ParStat, file=txt_path, quote=FALSE)

